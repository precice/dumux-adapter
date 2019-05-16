// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup BoundaryTests
 * \brief A simple free flow test problem for the staggered grid (Navier-)Stokes model.
 */

#ifndef DUMUX_FREEFLOW_SUBPROBLEM_HH
#define DUMUX_FREEFLOW_SUBPROBLEM_HH

#ifndef ENABLEMONOLITHIC
#define ENABLEMONOLITHIC 0
#endif

#include <dune/grid/yaspgrid.hh>

#include <dumux/material/fluidsystems/1pgas.hh>
#include <dumux/material/components/air.hh>

#include <dumux/freeflow/navierstokes/problem.hh>
#include <dumux/discretization/staggered/freeflow/properties.hh>
#include <dumux/freeflow/navierstokes/model.hh>

#if !ENABLEMONOLITHIC
#include "../common/preciceadapter.hh"
#endif

namespace Dumux {
template <class TypeTag>
class FreeFlowSubProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct FreeFlowModel { using InheritsFrom = std::tuple<NavierStokesNI, StaggeredFreeFlowModel>; };
} // end namespace TTag

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::FreeFlowModel>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePGas<Scalar, Dumux::Components::Air<Scalar> > ;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::FreeFlowModel> { using type = Dune::YaspGrid<2, Dune::TensorProductCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2> >; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::FreeFlowModel> { using type = Dumux::FreeFlowSubProblem<TypeTag> ; };

template<class TypeTag>
struct EnableFVGridGeometryCache<TypeTag, TTag::FreeFlowModel> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::FreeFlowModel> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::FreeFlowModel> { static constexpr bool value = true; };
} // end namespace Properties

/*!
 * \ingroup BoundaryTests
 * \brief Test problem for the one-phase (Navier-) Stokes problem.
 *
 * Horizontal flow from left to right with a parabolic velocity profile.
 */
template <class TypeTag>
class FreeFlowSubProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;

    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;

    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;

#if ENABLEMONOLITHIC
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
#endif

public:
#if ENABLEMONOLITHIC
    FreeFlowSubProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(fvGridGeometry, "FreeFlow"), eps_(1e-6), couplingManager_(couplingManager)
#else
    FreeFlowSubProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
      : ParentType(fvGridGeometry, "FreeFlow"), eps_(1e-6), couplingInterface_(precice_adapter::PreciceAdapter::getInstance())
#endif
    {
        problemName_  =  getParam<std::string>("Vtk.OutputName") + "_" + getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");
        inletVelocity_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.InletVelocity");
        inletTemperature_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.InletTemperature");
    }

    /*!
     * \brief The problem name.
     */
    const std::string& name() const
    {
        return problemName_;
    }

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

        const auto& globalPos = scvf.dofPosition();

        // inlet
        if (onLeftBoundary_(globalPos))
        {
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);
            values.setDirichlet(Indices::temperatureIdx);
        }
        //  outlet
        else if (onRightBoundary_(globalPos))
        {
            values.setDirichlet(Indices::pressureIdx);
            values.setOutflow(Indices::energyEqIdx);
        }
        else // wall
        {
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);
            values.setNeumann(Indices::energyEqIdx);
        }

#if ENABLEMONOLITHIC
        if (couplingManager().isCoupledEntity(CouplingManager::freeFlowIdx, scvf))
        {
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);
            values.setCouplingNeumann(Indices::energyEqIdx);
        }
#else
        //TODO preCICE
        const auto faceId = scvf.index();
        if ( couplingInterface_.isCoupledEntity(faceId) )
        {
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);
            values.setNeumann(Indices::energyEqIdx);
        }
#endif

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet control volume.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        return initialAtPos(globalPos);
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann control volume.
     *
     * \param element The element for which the Neumann boundary condition is set
     * \param fvGeometry The fvGeometry
     * \param elemVolVars The element volume variables
     * \param elemFaceVars The element face variables
     * \param scvf The boundary sub control volume face
     */
    template<class ElementVolumeVariables, class ElementFaceVariables>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFaceVariables& elemFaceVars,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);

#if ENABLEMONOLITHIC
        static const auto avgType = FreeFlowHeatCouplingOptions::stringToEnum(getParamFromGroup<std::string>(this->paramGroup(), "Problem.DiffCoeffAvgType", "FreeFlowOnly"));
        if (couplingManager().isCoupledEntity(CouplingManager::freeFlowIdx, scvf))
            values[Indices::energyEqIdx] = couplingManager().couplingData().energyCouplingCondition(element, fvGeometry, elemVolVars, elemFaceVars, scvf, avgType);
#else
        // TODO preCICE
        const auto faceId = scvf.index();
        if ( couplingInterface_.isCoupledEntity( faceId ) )
        {
          // Actually gets temperature
          const auto& scv = fvGeometry.scv(scvf.insideScvIdx());
          const auto& volVars = elemVolVars[scv];
          const Scalar cellCenterTemperature = volVars.temperature();
          const Scalar distance = (scvf.center() - scv.center()).two_norm();
          const Scalar insideLambda = volVars.effectiveThermalConductivity();
          const Scalar boundaryTemperature = couplingInterface_.getTemperatureOnFace(faceId);
          // q = -lambda * (t_face - t_cc) / dx
          //const Scalar qHeat = -insideLambda * (cellCenterTemperature-boundaryTemperature) / distance;
          const Scalar qHeat = insideLambda * (cellCenterTemperature-boundaryTemperature) / distance;

          values[Indices::energyEqIdx] = qHeat;
        }
#endif

        return values;
    }


#if ENABLEMONOLITHIC
    //! Get the coupling manager
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }
#endif

   /*!
     * \name Volume terms
     */
    // \{

   /*!
     * \brief Evaluates the initial value for a control volume.
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values(0.0);

        values[Indices::pressureIdx] = 1e5;
        values[Indices::energyEqIdx] = inletTemperature_;

        // parabolic velocity profile
        values[Indices::velocityXIdx] =  inletVelocity_*(globalPos[1] - this->fvGridGeometry().bBoxMin()[1])*(this->fvGridGeometry().bBoxMax()[1] - globalPos[1])
                                      / (0.25*(this->fvGridGeometry().bBoxMax()[1] - this->fvGridGeometry().bBoxMin()[1])*(this->fvGridGeometry().bBoxMax()[1] - this->fvGridGeometry().bBoxMin()[1]));

        values[Indices::velocityYIdx] = 0.0;

        return values;
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

    Scalar eps_;
    std::string problemName_;
    Scalar inletVelocity_;
    Scalar inletTemperature_;

#if ENABLEMONOLITHIC
    std::shared_ptr<CouplingManager> couplingManager_;
#else
   precice_adapter::PreciceAdapter& couplingInterface_;
#endif

};
} // end namespace Dumux

#endif // DUMUX_FREEFLOW_SUBPROBLEM_HH
