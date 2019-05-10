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
 * \brief The free-flow sub problem
 */
#ifndef DUMUX_FREEFLOW1P2C_SUBPROBLEM_HH
#define DUMUX_FREEFLOW1P2C_SUBPROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/material/fluidsystems/1padapter.hh>
#include <dumux/material/fluidsystems/h2oair.hh>
#include <dumux/discretization/staggered/freeflow/properties.hh>

#if EXNUMBER >= 1
#include <dumux/freeflow/compositional/zeroeqncmodel.hh>
#include <dumux/freeflow/rans/zeroeq/problem.hh>
#else
#include <dumux/freeflow/compositional/navierstokesncmodel.hh>
#include <dumux/freeflow/navierstokes/problem.hh>
#endif

namespace Dumux {

template <class TypeTag>
class FreeFlowSubProblem;

namespace Properties {

// Create new type tags
namespace TTag {
#if EXNUMBER >= 1
struct StokesZeroEq { using InheritsFrom = std::tuple<ZeroEqNCNI, StaggeredFreeFlowModel>; };
#else
struct StokesZeroEq { using InheritsFrom = std::tuple<NavierStokesNCNI, StaggeredFreeFlowModel>; };
#endif
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::StokesZeroEq> { using type = Dune::YaspGrid<2, Dune::TensorProductCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2> >; };

// The fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::StokesZeroEq>
{
  using H2OAir = FluidSystems::H2OAir<GetPropType<TypeTag, Properties::Scalar>>;
  static constexpr auto phaseIdx = H2OAir::gasPhaseIdx; // simulate the air phase
  using type = FluidSystems::OnePAdapter<H2OAir, phaseIdx>;
};

template<class TypeTag>
struct ReplaceCompEqIdx<TypeTag, TTag::StokesZeroEq> { static constexpr int value = 3; };

// Use formulation based on mass fractions
template<class TypeTag>
struct UseMoles<TypeTag, TTag::StokesZeroEq> { static constexpr bool value = true; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::StokesZeroEq> { using type = Dumux::FreeFlowSubProblem<TypeTag> ; };

template<class TypeTag>
struct EnableFVGridGeometryCache<TypeTag, TTag::StokesZeroEq> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::StokesZeroEq> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::StokesZeroEq> { static constexpr bool value = true; };
}

/*!
 * \brief The free-flow sub problem
 */
template <class TypeTag>
#if EXNUMBER >= 1
class FreeFlowSubProblem : public ZeroEqProblem<TypeTag>
{
    using ParentType = ZeroEqProblem<TypeTag>;
#else
class FreeFlowSubProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;
#endif

    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;

    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using ElementFaceVariables = typename GetPropType<TypeTag, Properties::GridFaceVariables>::LocalView;
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using TimeLoopPtr = std::shared_ptr<TimeLoop<Scalar>>;

    using DiffusionCoefficientAveragingType = typename StokesDarcyCouplingOptions::DiffusionCoefficientAveragingType;

    static constexpr bool useMoles = GetPropType<TypeTag, Properties::ModelTraits>::useMoles();

public:
    FreeFlowSubProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(fvGridGeometry, "Stokes"), eps_(1e-6), couplingManager_(couplingManager)
    {
        refVelocity_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.RefVelocity");
        refPressure_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.RefPressure");
        refMoleFrac_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.refMoleFrac");
        refTemperature_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.RefTemperature");

        diffCoeffAvgType_ = StokesDarcyCouplingOptions::stringToEnum(DiffusionCoefficientAveragingType{},
                                                                     getParamFromGroup<std::string>(this->paramGroup(), "Problem.InterfaceDiffusionCoefficientAvg"));
    }

   /*!
     * \name Problem parameters
     */
    // \{

   /*!
     * \brief Return the temperature within the domain in [K].
     */
    Scalar temperature() const
    { return refTemperature_; }

   /*!
     * \brief Return the sources within the domain.
     *
     * \param globalPos The global position
     */
    NumEqVector sourceAtPos(const GlobalPosition &globalPos) const
    { return NumEqVector(0.0); }

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
    BoundaryTypes boundaryTypes(const Element& element,
                                const SubControlVolumeFace& scvf) const
    {
        BoundaryTypes values;

        const auto& globalPos = scvf.center();

        if (onLeftBoundary_(globalPos))
        {
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);
            values.setDirichlet(Indices::conti0EqIdx + 1);
            values.setDirichlet(Indices::energyEqIdx);
        }

        if (onLowerBoundary_(globalPos))
        {
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);
            values.setNeumann(Indices::conti0EqIdx);
            values.setNeumann(Indices::conti0EqIdx + 1);
            values.setNeumann(Indices::energyEqIdx);
        }

        if (onUpperBoundary_(globalPos))
        {
#if EXNUMBER >=2
            values.setAllSymmetry();
#else
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);
            values.setNeumann(Indices::conti0EqIdx);
            values.setNeumann(Indices::conti0EqIdx + 1);
            values.setNeumann(Indices::energyEqIdx);
#endif
        }

        if (onRightBoundary_(globalPos))
        {
            values.setDirichlet(Indices::pressureIdx);
            values.setOutflow(Indices::conti0EqIdx + 1);
            values.setOutflow(Indices::energyEqIdx);
        }

        if (couplingManager().isCoupledEntity(CouplingManager::stokesIdx, scvf))
        {
            values.setCouplingNeumann(Indices::conti0EqIdx);
            values.setCouplingNeumann(Indices::conti0EqIdx + 1);
            values.setCouplingNeumann(Indices::energyEqIdx);
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
    PrimaryVariables dirichletAtPos(const GlobalPosition& pos) const
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
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFaceVariables& elemFaceVars,
                        const SubControlVolumeFace& scvf) const
    {
        PrimaryVariables values(0.0);
        if(couplingManager().isCoupledEntity(CouplingManager::stokesIdx, scvf))
        {
            values[Indices::momentumYBalanceIdx] = couplingManager().couplingData().momentumCouplingCondition(element, fvGeometry, elemVolVars, elemFaceVars, scvf);

            const auto massFlux = couplingManager().couplingData().massCouplingCondition(element, fvGeometry, elemVolVars, elemFaceVars, scvf, diffCoeffAvgType_);
            values[Indices::conti0EqIdx] = massFlux[0];
            values[Indices::conti0EqIdx + 1] = massFlux[1];
            values[Indices::energyEqIdx] = couplingManager().couplingData().energyCouplingCondition(element, fvGeometry, elemVolVars, elemFaceVars, scvf, diffCoeffAvgType_);
        }
        return values;
    }

    // \}

    /*!
     * \brief Set the coupling manager
     */
    void setCouplingManager(std::shared_ptr<CouplingManager> cm)
    { couplingManager_ = cm; }

    /*!
     * \brief Get the coupling manager
     */
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

#if EXNUMBER >= 2
    bool isOnWallAtPos(const GlobalPosition& globalPos) const
    {
        return (onLowerBoundary_(globalPos));
    }
#elif EXNUMBER >= 1
    bool isOnWallAtPos(const GlobalPosition& globalPos) const
    {
        return (onLowerBoundary_(globalPos) || onUpperBoundary_(globalPos));
    }
#endif

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
        FluidState fluidState;
        updateFluidStateForBC_(fluidState, refTemperature(), refPressure(), refMoleFrac());

        const Scalar density = FluidSystem::density(fluidState, 0);

        PrimaryVariables values(0.0);
        values[Indices::pressureIdx] = refPressure() + density*this->gravity()[1]*(globalPos[1] - this->fvGridGeometry().bBoxMin()[1]);
        values[Indices::conti0EqIdx + 1] = refMoleFrac();
        values[Indices::velocityXIdx] = refVelocity();
        values[Indices::temperatureIdx] = refTemperature();

#if EXNUMBER >= 2
        if(onLowerBoundary_(globalPos))
            values[Indices::velocityXIdx] = 0.0;
#else
        if(onUpperBoundary_(globalPos) || onLowerBoundary_(globalPos))
            values[Indices::velocityXIdx] = 0.0;
#endif

        return values;
    }

    //! \brief Returns the reference velocity.
    const Scalar refVelocity() const
    { return refVelocity_ ;}

    //! \brief Returns the reference pressure.
    const Scalar refPressure() const
    { return refPressure_; }

    //! \brief Returns the reference mass fraction.
    const Scalar refMoleFrac() const
    { return refMoleFrac_; }

    //! \brief Returns the reference temperature.
    const Scalar refTemperature() const
    { return refTemperature_; }


    void setTimeLoop(TimeLoopPtr timeLoop)
    { timeLoop_ = timeLoop; }

    /*!
     * \brief Returns the intrinsic permeability of required as input parameter for the Beavers-Joseph-Saffman boundary condition
     */
    Scalar permeability(const Element& element, const SubControlVolumeFace& scvf) const
    {
        return couplingManager().problem(CouplingManager::darcyIdx).spatialParams().permeabilityAtPos(scvf.center());
    }

    /*!
     * \brief Returns the alpha value required as input parameter for the Beavers-Joseph-Saffman boundary condition
     */
    Scalar alphaBJ(const SubControlVolumeFace& scvf) const
    {
        return couplingManager().problem(CouplingManager::darcyIdx).spatialParams().beaversJosephCoeffAtPos(scvf.center());
    }

    // \}

private:
    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] < this->fvGridGeometry().bBoxMin()[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] > this->fvGridGeometry().bBoxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] < this->fvGridGeometry().bBoxMin()[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] > this->fvGridGeometry().bBoxMax()[1] - eps_; }

    //! \brief updates the fluid state to obtain required quantities for IC/BC
    void updateFluidStateForBC_(FluidState& fluidState, const Scalar temperature,
                                const Scalar pressure, const Scalar moleFraction) const
    {
        fluidState.setTemperature(temperature);
        fluidState.setPressure(0, pressure);
        fluidState.setSaturation(0, 1.0);
        fluidState.setMoleFraction(0, 1, moleFraction);
        fluidState.setMoleFraction(0, 0, 1.0 - moleFraction);

        typename FluidSystem::ParameterCache paramCache;
        paramCache.updatePhase(fluidState, 0);

        const Scalar density = FluidSystem::density(fluidState, paramCache, 0);
        fluidState.setDensity(0, density);

        const Scalar molarDensity = FluidSystem::molarDensity(fluidState, paramCache, 0);
        fluidState.setMolarDensity(0, molarDensity);

        const Scalar enthalpy = FluidSystem::enthalpy(fluidState, paramCache, 0);
        fluidState.setEnthalpy(0, enthalpy);
    }

    // the height of the free-flow domain
    const Scalar height_() const
    { return this->fvGridGeometry().bBoxMax()[1] - this->fvGridGeometry().bBoxMin()[1]; }

    Scalar eps_;

    Scalar refVelocity_;
    Scalar refPressure_;
    Scalar refMoleFrac_;
    Scalar refTemperature_;

    TimeLoopPtr timeLoop_;

    std::shared_ptr<CouplingManager> couplingManager_;

    DiffusionCoefficientAveragingType diffCoeffAvgType_;
};
} //end namespace

#endif // DUMUX_STOKES1P2C_SUBPROBLEM_HH
