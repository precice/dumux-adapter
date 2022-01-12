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
 * \ingroup NavierStokesTests
 * \brief A simple Stokes test problem for the staggered grid (Navier-)Stokes model.
 */
#ifndef DUMUX_STOKES1P2C_SUBPROBLEM_HH
#define DUMUX_STOKES1P2C_SUBPROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/material/fluidsystems/1padapter.hh>
#include <dumux/material/fluidsystems/h2oair.hh>

#include <dumux/discretization/staggered/freeflow/properties.hh>
#include <dumux/freeflow/compositional/navierstokesncmodel.hh>
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
struct StokesNC {
    using InheritsFrom = std::tuple<NavierStokesNC, StaggeredFreeFlowModel>;
};
}  // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::StokesNC> {
    using type = Dune::YaspGrid<2,
                                Dune::EquidistantOffsetCoordinates<
                                    GetPropType<TypeTag, Properties::Scalar>,
                                    2>>;
};

// The fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::StokesNC> {
    using H2OAir =
        FluidSystems::H2OAir<GetPropType<TypeTag, Properties::Scalar>>;
    using type = FluidSystems::OnePAdapter<H2OAir, H2OAir::gasPhaseIdx>;
};

// Do not replace one equation with a total mass balance
template<class TypeTag>
struct ReplaceCompEqIdx<TypeTag, TTag::StokesNC> {
    static constexpr int value = 3;
};

// Use formulation based on mass fractions
template<class TypeTag>
struct UseMoles<TypeTag, TTag::StokesNC> {
    static constexpr bool value = true;
};

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::StokesNC> {
    using type = Dumux::StokesSubProblem<TypeTag>;
};

template<class TypeTag>
struct EnableFVGridGeometryCache<TypeTag, TTag::StokesNC> {
    static constexpr bool value = true;
};
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::StokesNC> {
    static constexpr bool value = true;
};
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::StokesNC> {
    static constexpr bool value = true;
};
}  // namespace Properties

/*!
 * \ingroup NavierStokesTests
 * \brief  Test problem for the one-phase compositional (Navier-) Stokes problem.
 *
 * Horizontal flow from left to right with a parabolic velocity profile.
 */
template<class TypeTag>
class StokesSubProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;

    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Indices =
        typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;

    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolumeFace =
        typename FVElementGeometry::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementVolumeVariables =
        typename GetPropType<TypeTag,
                             Properties::GridVolumeVariables>::LocalView;
    using ElementFaceVariables =
        typename GetPropType<TypeTag, Properties::GridFaceVariables>::LocalView;
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using TimeLoopPtr = std::shared_ptr<TimeLoop<Scalar>>;

    static constexpr bool useMoles =
        GetPropType<TypeTag, Properties::ModelTraits>::useMoles();

    static constexpr auto dim =
        GetPropType<TypeTag, Properties::ModelTraits>::dim();
    static constexpr auto transportCompIdx = 1;

   public:
    StokesSubProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                     std::shared_ptr<CouplingManager> couplingManager)
        : ParentType(fvGridGeometry, "Stokes"),
          eps_(1e-6),
          couplingManager_(couplingManager)
    {
        velocity_ =
            getParamFromGroup<Scalar>(this->paramGroup(), "Problem.Velocity");
        pressure_ =
            getParamFromGroup<Scalar>(this->paramGroup(), "Problem.Pressure");
        moleFraction_ = getParamFromGroup<Scalar>(this->paramGroup(),
                                                  "Problem.MoleFraction");
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Return the temperature within the domain in [K].
     */
    Scalar temperature() const { return 293.15; }

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

        const auto &globalPos = scvf.center();

        if (onLeftBoundary_(globalPos)) {
            values.setDirichlet(Indices::conti0EqIdx + 1);
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);
        } else if (onRightBoundary_(globalPos)) {
            values.setDirichlet(Indices::pressureIdx);
            values.setOutflow(Indices::conti0EqIdx + 1);
        } else {
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);
            values.setNeumann(Indices::conti0EqIdx);
            values.setNeumann(Indices::conti0EqIdx + 1);
        }

        if (couplingManager().isCoupledEntity(CouplingManager::stokesIdx,
                                              scvf)) {
            values.setCouplingNeumann(Indices::conti0EqIdx);
            values.setCouplingNeumann(Indices::conti0EqIdx + 1);
            values.setCouplingNeumann(Indices::momentumYBalanceIdx);
            values.setBJS(Indices::momentumXBalanceIdx);
        }

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a Dirichlet control volume.
     *
     * \param element The element
     * \param scvf The subcontrolvolume face
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &pos) const
    {
        PrimaryVariables values(0.0);
        values = initialAtPos(pos);

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
    NumEqVector neumann(const Element &element,
                        const FVElementGeometry &fvGeometry,
                        const ElementVolumeVariables &elemVolVars,
                        const ElementFaceVariables &elemFaceVars,
                        const SubControlVolumeFace &scvf) const
    {
        PrimaryVariables values(0.0);

        if (couplingManager().isCoupledEntity(CouplingManager::stokesIdx,
                                              scvf)) {
            values[Indices::momentumYBalanceIdx] =
                couplingManager().couplingData().momentumCouplingCondition(
                    element, fvGeometry, elemVolVars, elemFaceVars, scvf);

            const auto massFlux =
                couplingManager().couplingData().massCouplingCondition(
                    element, fvGeometry, elemVolVars, elemFaceVars, scvf);
            values[Indices::conti0EqIdx] = massFlux[0];
            values[Indices::conti0EqIdx + 1] = massFlux[1];
        }
        return values;
    }

    // \}

    /*!
     * \brief Set the coupling manager
     */
    void setCouplingManager(std::shared_ptr<CouplingManager> cm)
    {
        couplingManager_ = cm;
    }

    /*!
     * \brief Get the coupling manager
     */
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

        // This is only an approximation of the actual hydorostatic pressure gradient.
        // Air is compressible (the density depends on pressure), thus the actual
        // vertical pressure profile is non-linear.
        // This discrepancy can lead to spurious flows at the outlet if the inlet
        // velocity is chosen too small.
        FluidState fluidState;
        updateFluidStateForBC_(fluidState, pressure_);
        const Scalar density = FluidSystem::density(fluidState, 0);
        values[Indices::pressureIdx] =
            pressure_ -
            density * this->gravity()[1] *
                (this->fvGridGeometry().bBoxMax()[1] - globalPos[1]);

        // gravity has negative sign
        values[Indices::conti0EqIdx + 1] = moleFraction_;
        values[Indices::velocityXIdx] =
            4.0 * velocity_ *
            (globalPos[1] - this->fvGridGeometry().bBoxMin()[1]) *
            (this->fvGridGeometry().bBoxMax()[1] - globalPos[1]) /
            (height_() * height_());

        return values;
    }

    void setTimeLoop(TimeLoopPtr timeLoop) { timeLoop_ = timeLoop; }

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

    //! \brief updates the fluid state to obtain required quantities for IC/BC
    void updateFluidStateForBC_(FluidState &fluidState,
                                const Scalar pressure) const
    {
        fluidState.setTemperature(temperature());
        fluidState.setPressure(0, pressure);
        fluidState.setSaturation(0, 1.0);
        fluidState.setMoleFraction(0, 1, moleFraction_);
        fluidState.setMoleFraction(0, 0, 1.0 - moleFraction_);

        typename FluidSystem::ParameterCache paramCache;
        paramCache.updatePhase(fluidState, 0);

        const Scalar density = FluidSystem::density(fluidState, paramCache, 0);
        fluidState.setDensity(0, density);

        const Scalar molarDensity =
            FluidSystem::molarDensity(fluidState, paramCache, 0);
        fluidState.setMolarDensity(0, molarDensity);

        const Scalar enthalpy =
            FluidSystem::enthalpy(fluidState, paramCache, 0);
        fluidState.setEnthalpy(0, enthalpy);
    }
    // the height of the free-flow domain
    const Scalar height_() const
    {
        return this->fvGridGeometry().bBoxMax()[1] -
               this->fvGridGeometry().bBoxMin()[1];
    }

    Scalar eps_;

    Scalar velocity_;
    Scalar pressure_;
    Scalar moleFraction_;

    TimeLoopPtr timeLoop_;

    std::shared_ptr<CouplingManager> couplingManager_;
};
}  // namespace Dumux

#endif  // DUMUX_STOKES1P2C_SUBPROBLEM_HH
