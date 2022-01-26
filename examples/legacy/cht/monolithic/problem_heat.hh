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
 * \brief A simple heat test problem (cell-centered finite volume method).
 */

#ifndef DUMUX_HEAT_SUBPROBLEM_HH
#define DUMUX_HEAT_SUBPROBLEM_HH

#ifndef ENABLEMONOLITHIC
#define ENABLEMONOLITHIC 1
#endif

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/porousmediumflow/solidenergy/model.hh>

#include <dumux/porousmediumflow/problem.hh>

#include "spatialparams.hh"

#if !ENABLEMONOLITHIC
#include "../common/preciceadapter.hh"
#endif

namespace Dumux
{
template<class TypeTag>
class HeatSubProblem;

namespace Properties
{
// Create new type tags
namespace TTag
{
struct HeatModel {
    using InheritsFrom = std::tuple<SolidEnergy, CCTpfaModel>;
};
}  // end namespace TTag

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::HeatModel> {
    using type = Dumux::HeatSubProblem<TypeTag>;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::HeatModel> {
    using type = Dune::YaspGrid<
        2,
        Dune::TensorProductCoordinates<GetPropType<TypeTag, Properties::Scalar>,
                                       2> >;
};

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::HeatModel> {
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = OnePSpatialParams<FVGridGeometry, Scalar>;
};
}  // end namespace Properties

template<class TypeTag>
class HeatSubProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    using FVElementGeometry =
        typename GetPropType<TypeTag, Properties::FVGridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace =
        typename FVElementGeometry::SubControlVolumeFace;
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;

    using Indices =
        typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

#if ENABLEMONOLITHIC
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
#endif

public:
#if ENABLEMONOLITHIC
    HeatSubProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                   std::shared_ptr<CouplingManager> couplingManager)
        : ParentType(fvGridGeometry, "SolidEnergy"),
          eps_(1e-7),
          couplingManager_(couplingManager)
#else
    HeatSubProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
        : ParentType(fvGridGeometry, "SolidEnergy"),
          eps_(1e-7),
          couplingInterface_(precice_adapter::PreciceAdapter::getInstance())
#endif
    {
        problemName_ =
            getParam<std::string>("Vtk.OutputName") + "_" +
            getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");
        bottomTemperature_ = getParamFromGroup<Scalar>(
            this->paramGroup(), "Problem.BottomTemperature");
        initialTemperature_ = getParamFromGroup<Scalar>(
            this->paramGroup(), "Problem.InitialTemperature");
    }

    /*!
     * \brief The problem name.
     */
    const std::string &name() const { return problemName_; }

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
      * \brief Specifies which kind of boundary condition should be
      *        used for which equation on a given boundary control volume.
      *
      * \param element The element
      * \param scvf The boundary sub control volume face
      */
    BoundaryTypes boundaryTypes(const Element &element,
                                const SubControlVolumeFace &scvf) const
    {
        BoundaryTypes values;
        values.setAllNeumann();

        if (onLowerBoundary_(scvf.center()))
            values.setAllDirichlet();

#if ENABLEMONOLITHIC
        if (couplingManager().isCoupledEntity(CouplingManager::solidEnergyIdx,
                                              scvf))
            values.setAllCouplingNeumann();
#endif
        // Note: No special treatment for preCICE needed, all BCs are Neumann anyways

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet control volume.
     *
     * \param element The element for which the Dirichlet boundary condition is set
     * \param scvf The boundary subcontrolvolumeface
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichlet(const Element &element,
                               const SubControlVolumeFace &scvf) const
    {
        return initialAtPos(scvf.center());
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann control volume.
     *
     * \param element The element for which the Neumann boundary condition is set
     * \param fvGeometry The fvGeometry
     * \param elemVolVars The element volume variables
     * \param scvf The boundary sub control volume face
     *
     * For this method, the \a values variable stores primary variables.
     */
    template<class ElementVolumeVariables>
    NumEqVector neumann(const Element &element,
                        const FVElementGeometry &fvGeometry,
                        const ElementVolumeVariables &elemVolVars,
                        const SubControlVolumeFace &scvf) const
    {
        NumEqVector values(0.0);

#if ENABLEMONOLITHIC
        static const auto avgType = FreeFlowHeatCouplingOptions::stringToEnum(
            getParamFromGroup<std::string>(this->paramGroup(),
                                           "Problem.DiffCoeffAvgType",
                                           "FreeFlowOnly"));
        if (couplingManager().isCoupledEntity(CouplingManager::solidEnergyIdx,
                                              scvf))
            values = couplingManager().couplingData().energyCouplingCondition(
                element, fvGeometry, elemVolVars, scvf, avgType);
#else
        //TODO preCICE
        const auto faceId = scvf.index();
        if (couplingInterface_.isCoupledEntity(faceId)) {
            // ALEX: I made the heat flux negative since the normal vectors of both solvers are different!
            values[Indices::energyEqIdx] =
                -couplingInterface_.getHeatFluxOnFace(faceId);
        }
#endif

        return values;
    }

    // \}

    /*!
     * \brief Evaluates the initial value for a control volume.
     *
     * \param element The element
     *
     * For this method, the \a priVars parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition &pos) const
    {
        PrimaryVariables values(initialTemperature_);

        if (onLowerBoundary_(pos))
            values = bottomTemperature_;

        return values;
    }

    // \}

#if ENABLEMONOLITHIC
    //! Get the coupling manager
    const CouplingManager &couplingManager() const { return *couplingManager_; }
#endif

private:
    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] < this->fvGridGeometry().bBoxMin()[1] + eps_;
    }

    Scalar eps_;
    std::string problemName_;
    Scalar bottomTemperature_;
    Scalar initialTemperature_;

#if ENABLEMONOLITHIC
    std::shared_ptr<CouplingManager> couplingManager_;
#else
    precice_adapter::PreciceAdapter &couplingInterface_;
#endif
};
}  // end namespace Dumux

#endif  //DUMUX_HEAT_SUBPROBLEM_HH
