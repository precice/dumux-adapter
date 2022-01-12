// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \brief The free flow sub problem
 */
#ifndef DUMUX_STOKES_SUBPROBLEM_HH
#define DUMUX_STOKES_SUBPROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#if EXNUMBER >= 3
#include <dumux/io/grid/subgridgridcreator.hh>
#endif

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include <dumux/discretization/staggered/freeflow/properties.hh>
#include <dumux/freeflow/navierstokes/model.hh>
#include <dumux/freeflow/navierstokes/problem.hh>

namespace Dumux
{
template<class TypeTag>
class StokesSubProblem;

namespace Properties
{
// Create new type tags
namespace TTag
{
struct StokesOneP {
    using InheritsFrom = std::tuple<NavierStokes, StaggeredFreeFlowModel>;
};
}  // end namespace TTag

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::StokesOneP> {
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type =
        FluidSystems::OnePLiquid<Scalar, Dumux::Components::SimpleH2O<Scalar> >;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::StokesOneP> {
    static constexpr auto dim = 2;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using TensorGrid =
        Dune::YaspGrid<2, Dune::TensorProductCoordinates<Scalar, dim> >;

#if EXNUMBER < 3  // use "normal" grid
    using type = TensorGrid;
#else  // use dune-subgrid
    using HostGrid = TensorGrid;
    using type = Dune::SubGrid<dim, HostGrid>;
#endif
};

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::StokesOneP> {
    using type = Dumux::StokesSubProblem<TypeTag>;
};

template<class TypeTag>
struct EnableFVGridGeometryCache<TypeTag, TTag::StokesOneP> {
    static constexpr bool value = true;
};
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::StokesOneP> {
    static constexpr bool value = true;
};
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::StokesOneP> {
    static constexpr bool value = true;
};
}  // namespace Properties

/*!
 * \brief The free flow sub problem
 */
template<class TypeTag>
class StokesSubProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;

    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    using Indices =
        typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;

    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolumeFace =
        typename FVElementGeometry::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

   public:
    StokesSubProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                     std::shared_ptr<CouplingManager> couplingManager)
        : ParentType(fvGridGeometry, "Stokes"),
          eps_(1e-6),
          couplingManager_(couplingManager)
    {
        deltaP_ = getParamFromGroup<Scalar>(this->paramGroup(),
                                            "Problem.PressureDifference");
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Return the temperature within the domain in [K].
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const { return 273.15 + 10; }  // 10°C

    /*!
     * \brief Return the sources within the domain.
     *
     * \param globalPos The global position
     */
    NumEqVector sourceAtPos(const GlobalPosition &globalPos) const
    {
        return NumEqVector(0.0);
    }
    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param element The finite element
     * \param scvf The sub control volume face
     */
    BoundaryTypes boundaryTypes(const Element &element,
                                const SubControlVolumeFace &scvf) const
    {
        BoundaryTypes values;

        const auto &globalPos = scvf.dofPosition();

#if EXNUMBER == 0  // flow from top to bottom
        if (onUpperBoundary_(globalPos)) {
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);
        }

        if (onRightBoundary_(globalPos) || (onLeftBoundary_(globalPos))) {
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);
        }
#else  // flow flom left to right
        if (onLeftBoundary_(globalPos) || onRightBoundary_(globalPos))
            values.setDirichlet(Indices::pressureIdx);
        else {
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);
        }
#endif

        if (couplingManager().isCoupledEntity(CouplingManager::stokesIdx,
                                              scvf)) {
            values.setCouplingNeumann(Indices::conti0EqIdx);
#if EXNUMBER < 3
            values.setCouplingNeumann(Indices::momentumYBalanceIdx);
#else
            //consider orientation of face automatically
            values.setCouplingNeumann(scvf.directionIndex());
#endif

#if EXNUMBER < 2
            values.setDirichlet(
                Indices::velocityXIdx);  // assume no slip on interface
#elif EXNUMBER == 2
            // set the Beaver-Joseph-Saffman slip condition for the tangential momentum balance equation
            values.setBJS(Indices::momentumXBalanceIdx);
#else
            // set the Beaver-Joseph-Saffman slip condition for the tangential momentum balance equation,
            // consider orientation of face automatically
            values.setBJS(1 - scvf.directionIndex());
#endif
        }

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a Dirichlet control volume.
     *
     * \param globalPos The global position
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);
        values = initialAtPos(globalPos);

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a Neumann control volume.
     *
     * \param element The element for which the Neumann boundary condition is set
     * \param fvGeomentry The fvGeometry
     * \param elemVolVars The element volume variables
     * \param elemFaceVars The element face variables
     * \param scvf The boundary sub control volume face
     */
    template<class ElementVolumeVariables, class ElementFaceVariables>
    NumEqVector neumann(const Element &element,
                        const FVElementGeometry &fvGeometry,
                        const ElementVolumeVariables &elemVolVars,
                        const ElementFaceVariables &elemFaceVars,
                        const SubControlVolumeFace &scvf) const
    {
        NumEqVector values(0.0);

        if (couplingManager().isCoupledEntity(CouplingManager::stokesIdx,
                                              scvf)) {
            values[Indices::conti0EqIdx] =
                couplingManager().couplingData().massCouplingCondition(
                    element, fvGeometry, elemVolVars, elemFaceVars, scvf);
#if EXNUMBER < 3
            values[Indices::momentumYBalanceIdx] =
                couplingManager().couplingData().momentumCouplingCondition(
                    element, fvGeometry, elemVolVars, elemFaceVars, scvf);
#else
            values[scvf.directionIndex()] =
                couplingManager().couplingData().momentumCouplingCondition(
                    element, fvGeometry, elemVolVars, elemFaceVars, scvf);
#endif
        }
        return values;
    }

    // \}

    //! Set the coupling manager
    void setCouplingManager(std::shared_ptr<CouplingManager> cm)
    {
        couplingManager_ = cm;
    }

    //! Get the coupling manager
    const CouplingManager &couplingManager() const { return *couplingManager_; }

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);
#if EXNUMBER == 0
        values[Indices::velocityYIdx] =
            -1e-6 * globalPos[0] *
            (this->fvGridGeometry().bBoxMax()[0] - globalPos[0]);
#else
        // set fixed pressures on the left and right boundary
        if (onLeftBoundary_(globalPos))
            values[Indices::pressureIdx] = deltaP_;
        if (onRightBoundary_(globalPos))
            values[Indices::pressureIdx] = 0.0;
#endif

        return values;
    }

    /*!
     * \brief Returns the intrinsic permeability of required as input parameter for the Beavers-Joseph-Saffman boundary condition
     */
    Scalar permeability(const Element &element,
                        const SubControlVolumeFace &scvf) const
    {
        return couplingManager().couplingData().darcyPermeability(element,
                                                                  scvf);
    }

    /*!
     * \brief Returns the alpha value required as input parameter for the Beavers-Joseph-Saffman boundary condition
     */
    Scalar alphaBJ(const SubControlVolumeFace &scvf) const
    {
        return couplingManager()
            .problem(CouplingManager::darcyIdx)
            .spatialParams()
            .beaversJosephCoeffAtPos(scvf.center());
    }

    /*!
     * \brief calculate the analytical velocity in x direction based on Beavers & Joseph (1967)
     */
    void calculateAnalyticalVelocityX() const
    {
        analyticalVelocityX_.resize(this->fvGridGeometry().gridView().size(0));

        using std::sqrt;
        const Scalar dPdX = -deltaP_ / (this->fvGridGeometry().bBoxMax()[0] -
                                        this->fvGridGeometry().bBoxMin()[0]);
        static const Scalar mu = FluidSystem::viscosity(temperature(), 1e5);
        static const Scalar alpha =
            getParam<Scalar>("Darcy.SpatialParams.AlphaBeaversJoseph");
        static const Scalar K =
            getParam<Scalar>("Darcy.SpatialParams.Permeability");
        static const Scalar sqrtK = sqrt(K);
        const Scalar sigma = (this->fvGridGeometry().bBoxMax()[1] -
                              this->fvGridGeometry().bBoxMin()[1]) /
                             sqrtK;

        const Scalar uB =
            -K / (2.0 * mu) *
            ((sigma * sigma + 2.0 * alpha * sigma) / (1.0 + alpha * sigma)) *
            dPdX;

        for (const auto &element :
             elements(this->fvGridGeometry().gridView())) {
            const auto eIdx =
                this->fvGridGeometry().gridView().indexSet().index(element);
            const Scalar y = element.geometry().center()[1] -
                             this->fvGridGeometry().bBoxMin()[1];

            const Scalar u =
                uB * (1.0 + alpha / sqrtK * y) +
                1.0 / (2.0 * mu) * (y * y + 2 * alpha * y * sqrtK) * dPdX;
            analyticalVelocityX_[eIdx] = u;
        }
    }

    /*!
     * \brief Get the analytical velocity in x direction
     */
    const std::vector<Scalar> &getAnalyticalVelocityX() const
    {
        if (analyticalVelocityX_.empty())
            calculateAnalyticalVelocityX();
        return analyticalVelocityX_;
    }

    // \}

   private:
    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] < this->fvGridGeometry().bBoxMin()[0] + eps_;
    }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] > this->fvGridGeometry().bBoxMax()[0] - eps_;
    }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] < this->fvGridGeometry().bBoxMin()[1] + eps_;
    }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] > this->fvGridGeometry().bBoxMax()[1] - eps_;
    }

    Scalar eps_;
    Scalar deltaP_;

    std::shared_ptr<CouplingManager> couplingManager_;

    mutable std::vector<Scalar> analyticalVelocityX_;
};
}  // namespace Dumux

#endif  // DUMUX_STOKES_SUBPROBLEM_HH
