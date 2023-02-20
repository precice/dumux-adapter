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
#if DUMUX_VERSION_MAJOR >= 3 & DUMUX_VERSION_MINOR >= 4
#include <dumux/common/numeqvector.hh>
#endif

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include <dumux/discretization/staggered/freeflow/properties.hh>
#include <dumux/freeflow/navierstokes/model.hh>
#if DUMUX_VERSION_MAJOR >= 3 & DUMUX_VERSION_MINOR >= 6
#include <dumux/freeflow/navierstokes/staggered/problem.hh>
#else
#include <dumux/freeflow/navierstokes/problem.hh>
#endif

#include <dumux-precice/couplingadapter.hh>

namespace Dumux
{
template<class TypeTag>
class StokesSubProblem;

namespace Properties
{
// Create new type tags
namespace TTag
{
struct FreeFlowModel {
    using InheritsFrom = std::tuple<NavierStokes, StaggeredFreeFlowModel>;
};
}  // end namespace TTag

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::FreeFlowModel> {
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type =
        FluidSystems::OnePLiquid<Scalar, Dumux::Components::SimpleH2O<Scalar> >;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::FreeFlowModel> {
    using type = Dune::YaspGrid<2,
                                Dune::EquidistantOffsetCoordinates<
                                    GetPropType<TypeTag, Properties::Scalar>,
                                    2> >;
};

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::FreeFlowModel> {
    using type = Dumux::StokesSubProblem<TypeTag>;
};

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::FreeFlowModel> {
    static constexpr bool value = true;
};
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::FreeFlowModel> {
    static constexpr bool value = true;
};
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::FreeFlowModel> {
    static constexpr bool value = true;
};
}  // end namespace Properties

/*!
 * \brief The free flow sub problem
 */
template<class TypeTag>
#if DUMUX_VERSION_MAJOR >= 3 & DUMUX_VERSION_MINOR >= 6
class StokesSubProblem : public NavierStokesStaggeredProblem<TypeTag>
{
    using ParentType = NavierStokesStaggeredProblem<TypeTag>;
#else
class StokesSubProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;
#endif

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    using Indices =
        typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    using BoundaryTypes = Dumux::NavierStokesBoundaryTypes<
        GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;

    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace =
        typename FVElementGeometry::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;

#if DUMUX_VERSION_MAJOR >= 3 & DUMUX_VERSION_MINOR >= 4
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
#else
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
#endif

    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

public:
    StokesSubProblem(std::shared_ptr<const GridGeometry> gridGeometry)
        : ParentType(gridGeometry, "FreeFlow"),
          eps_(1e-6),
          couplingInterface_(Dumux::Precice::CouplingAdapter::getInstance()),
          pressureId_(0),
          velocityId_(0),
          dataIdsWereSet_(false)
    {
        deltaP_ = getParamFromGroup<Scalar>(this->paramGroup(),
                                            "Problem.PressureDifference");
    }

    /*!
     * \name Problem parameters
     */
    // \{

#if DUMUX_VERSION_MAJOR >= 3 & DUMUX_VERSION_MINOR < 5
    /*!
     * \brief Return the temperature within the domain in [K].
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const { return 273.15 + 10; }  // 10°C
#endif

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
        const auto faceId = scvf.index();

        // left/right wall
        if (onRightBoundary_(globalPos) || (onLeftBoundary_(globalPos))) {
            values.setDirichlet(Indices::pressureIdx);
        }
        // coupling interface
        else if (couplingInterface_.isCoupledEntity(faceId)) {
            assert(dataIdsWereSet_);

            values.setDirichlet(Indices::velocityYIdx);
            values.setBeaversJoseph(Indices::momentumXBalanceIdx);
        } else {
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);
        }

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a Dirichlet control volume.
     *
     * \param globalPos The global position
     */
    using ParentType::dirichlet;
    PrimaryVariables dirichlet(const Element &element,
                               const SubControlVolumeFace &scvf) const
    {
        PrimaryVariables values(0.0);
        values = initialAtPos(scvf.center());

        const auto faceId = scvf.index();
        if (couplingInterface_.isCoupledEntity(faceId)) {
            values[Indices::velocityYIdx] =
                couplingInterface_.getScalarQuantityOnFace(velocityId_, faceId);
        }

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

        assert(dataIdsWereSet_);
        const auto faceId = scvf.index();
        if (couplingInterface_.isCoupledEntity(faceId)) {
            const Scalar density =
                1000;  // TODO how to handle compressible fluids?
            values[Indices::conti0EqIdx] = density *
                                           elemFaceVars[scvf].velocitySelf() *
                                           scvf.directionSign();
            values[Indices::momentumYBalanceIdx] =
                scvf.directionSign() *
                (couplingInterface_.getScalarQuantityOnFace(pressureId_,
                                                            faceId) -
                 initialAtPos(scvf.center())[Indices::pressureIdx]);
        }
        return values;
    }

    // \}

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
        //values[Indices::velocityYIdx] = -1e-6 * globalPos[0] * (this->gridGeometry().bBoxMax()[0] - globalPos[0]);
        if (onLeftBoundary_(globalPos))
            values[Indices::pressureIdx] = deltaP_;
        if (onRightBoundary_(globalPos))
            values[Indices::pressureIdx] = 0.0;

        return values;
    }

    /*!
     * \brief Returns the intrinsic permeability of required as input parameter for the Beavers-Joseph-Saffman boundary condition
     */
    Scalar permeability(const Element &element,
                        const SubControlVolumeFace &scvf) const
    {
        return 1e-10;  // TODO transfer information or just use constant value
    }

    /*!
     * \brief Returns the alpha value required as input parameter for the Beavers-Joseph-Saffman boundary condition
     */
    Scalar alphaBJ(const SubControlVolumeFace &scvf) const
    {
        return 1.0;  // TODO transfer information or just use constant value
    }

    /*!
     * \brief calculate the analytical velocity in x direction based on Beavers & Joseph (1967)
     */
    void calculateAnalyticalVelocityX() const
    {
        analyticalVelocityX_.resize(this->gridGeometry().gridView().size(0));

        using std::sqrt;
        const Scalar dPdX = -deltaP_ / (this->gridGeometry().bBoxMax()[0] -
                                        this->gridGeometry().bBoxMin()[0]);
#if DUMUX_VERSION_MAJOR >= 3 & DUMUX_VERSION_MINOR > 4
        static const Scalar mu = FluidSystem::viscosity(273.15 + 10, 1e5);
#else
        static const Scalar mu = FluidSystem::viscosity(temperature(), 1e5);
#endif
        static const Scalar alpha =
            getParam<Scalar>("Darcy.SpatialParams.AlphaBeaversJoseph");
        static const Scalar K =
            getParam<Scalar>("Darcy.SpatialParams.Permeability");
        static const Scalar sqrtK = sqrt(K);
        const Scalar sigma = (this->gridGeometry().bBoxMax()[1] -
                              this->gridGeometry().bBoxMin()[1]) /
                             sqrtK;

        const Scalar uB =
            -K / (2.0 * mu) *
            ((sigma * sigma + 2.0 * alpha * sigma) / (1.0 + alpha * sigma)) *
            dPdX;

        for (const auto &element : elements(this->gridGeometry().gridView())) {
            const auto eIdx =
                this->gridGeometry().gridView().indexSet().index(element);
            const Scalar y = element.geometry().center()[1] -
                             this->gridGeometry().bBoxMin()[1];

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

    void updatePreciceDataIds()
    {
        pressureId_ = couplingInterface_.getIdFromName("Pressure");
        velocityId_ = couplingInterface_.getIdFromName("Velocity");
        dataIdsWereSet_ = true;
    }

    // \}

private:
    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_;
    }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_;
    }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] < this->gridGeometry().bBoxMin()[1] + eps_;
    }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_;
    }

    Scalar eps_;
    Scalar deltaP_;

    Dumux::Precice::CouplingAdapter &couplingInterface_;
    size_t pressureId_;
    size_t velocityId_;
    bool dataIdsWereSet_;

    mutable std::vector<Scalar> analyticalVelocityX_;
};
}  // namespace Dumux

#endif  // DUMUX_STOKES_SUBPROBLEM_HH
