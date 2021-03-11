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
 * \brief A simple Stokes test problem for the staggered grid (Navier-)Stokes model.
 */

#ifndef DUMUX_STOKES_SUBPROBLEM_HH
#define DUMUX_STOKES_SUBPROBLEM_HH

#include <dune/grid/yaspgrid.hh>

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
    using type = Dune::YaspGrid<2,
                                Dune::EquidistantOffsetCoordinates<
                                    GetPropType<TypeTag, Properties::Scalar>,
                                    2> >;
};

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::StokesOneP> {
    using type = Dumux::StokesSubProblem<TypeTag>;
};

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::StokesOneP> {
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
}  // end namespace Properties

/*!
 * \ingroup BoundaryTests
 * \brief Test problem for the one-phase (Navier-) Stokes problem.
 *
 * Horizontal flow from left to right with a parabolic velocity profile.
 */
template<class TypeTag>
class StokesSubProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using BoundaryTypes = Dumux::NavierStokesBoundaryTypes<
        GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace =
        typename FVElementGeometry::SubControlVolumeFace;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

   public:
    using Indices = typename ModelTraits::Indices;

    StokesSubProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                     std::shared_ptr<CouplingManager> couplingManager)
        : ParentType(gridGeometry, "Stokes"),
          eps_(1e-6),
          time_(0.0),
          couplingManager_(couplingManager)
    {
        problemName_ =
            getParam<std::string>("Vtk.OutputName") + "_" +
            getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");
    }

    /*!
     * \brief The problem name.
     */
    const std::string &name() const { return problemName_; }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Returns the temperature within the domain in [K].
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const { return 273.15 + 10; }  // 10°C

    /*!
     * \brief Returns the sources within the domain.
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

        //const auto &globalPos = scvf.dofPosition();

        //values.setDirichlet(Indices::velocityXIdx);
        //values.setDirichlet(Indices::velocityYIdx);

        if (couplingManager().isCoupledEntity(CouplingManager::stokesIdx,
                                              scvf)) {
            assert(
                couplingManager().couplingMode() !=
                CouplingManager::CouplingMode::reconstructPorousMediumPressure);

            values.setCouplingDirichlet(Indices::velocityYIdx);
            values.setBeaversJoseph(Indices::momentumXBalanceIdx);
        } else {
            values.setDirichlet(Indices::pressureIdx);
        }

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet control volume.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        return initialAtPos(globalPos);
    }

    PrimaryVariables dirichlet(const Element &element,
                               const SubControlVolume &scv) const
    {
        return initialAtPos(scv.center());
    }

    PrimaryVariables dirichlet(const Element &element,
                               const SubControlVolumeFace &scvf) const
    {
        if (couplingManager().isCoupledEntity(CouplingManager::stokesIdx,
                                              scvf)) {
            PrimaryVariables values(0.0);
            values[Indices::velocityYIdx] =
                couplingManager().couplingData().porousMediumInterfaceVelocity(
                    element, scvf);
            return values;
        } else
            return initialAtPos(scvf.center());
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
    NumEqVector neumann(const Element &element,
                        const FVElementGeometry &fvGeometry,
                        const ElementVolumeVariables &elemVolVars,
                        const ElementFaceVariables &elemFaceVars,
                        const SubControlVolumeFace &scvf) const
    {
        NumEqVector values(0.0);

        //std::cout << "Neumann BC: " << values << std::endl;

        if (couplingManager().isCoupledEntity(CouplingManager::stokesIdx,
                                              scvf)) {
            // std::cout << "vel is " << elemFaceVars[scvf].velocitySelf() << std::endl;
            values[Indices::conti0EqIdx] =
                couplingManager().couplingData().massCouplingCondition(
                    element, fvGeometry, elemVolVars, elemFaceVars, scvf);
            values[Indices::momentumYBalanceIdx] =
                couplingManager().couplingData().momentumCouplingCondition(
                    element, fvGeometry, elemVolVars, elemFaceVars, scvf);
        }
        return values;
    }

    // \}

    //! Get the coupling manager
    const CouplingManager &couplingManager() const { return *couplingManager_; }

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluates the initial value for a control volume.
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);

        values[Indices::velocityXIdx] = xVelocityAt(time_, globalPos);
        values[Indices::velocityYIdx] = yVelocityAt(time_, globalPos);
        values[Indices::pressureIdx] = pressureAt(time_, globalPos);

        std::cout << values[Indices::velocityXIdx] << std::endl;
        std::cout << values[Indices::velocityYIdx] << std::endl;
        std::cout << values[Indices::pressureIdx] << std::endl;
        //if (onLeftBoundary_(globalPos))
        //    values[Indices::pressureIdx] = deltaP;
        return values;
    }

    /*!
     * \brief Returns the intrinsic permeability of required as input parameter
              for the Beavers-Joseph-Saffman boundary condition
     */
    Scalar permeability(const Element &element,
                        const SubControlVolumeFace &scvf) const
    {
        return couplingManager().couplingData().darcyPermeability(element,
                                                                  scvf);
    }

    /*!
     * \brief Returns the alpha value required as input parameter for the
              Beavers-Joseph-Saffman boundary condition.
     */
    Scalar alphaBJ(const SubControlVolumeFace &scvf) const
    {
        return couplingManager()
            .problem(CouplingManager::darcyIdx)
            .spatialParams()
            .beaversJosephCoeffAtPos(scvf.center());
    }

    void setTime(const Scalar time) { time_ = time; }

    PrimaryVariables getExactSolution(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);

        values[Indices::pressureIdx] = pressureAt(time_, globalPos);
        std::cout << values[Indices::pressureIdx] << std::endl;
        //values[Indices::velocityXIdx] = xVelocityAt( time_, globalPos);
        //values[Indices::velocityYIdx] = yVelocityAt( time_, globalPos);

        return values;
    }

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

    Scalar pressureAt(const Scalar time, const GlobalPosition &globalPos) const
    {
        const auto x = globalPos[0];
        const auto y = globalPos[1];
        return (-std::pow(x, 2) * y + x * y + std::pow(y, 2)) *
               std::cos(M_PI * time);
    }

    Scalar xVelocityAt(const Scalar time, const GlobalPosition &globalPos) const
    {
        const auto x = globalPos[0];
        const auto y = globalPos[1];
        return (std::pow(y, 2) - 2.0 * y + 2.0 * x - 4.0 * std::pow(y, 3) * x -
                3) *
               std::cos(M_PI * time);
    };

    Scalar yVelocityAt(const Scalar time, const GlobalPosition &globalPos) const
    {
        const auto x = globalPos[0];
        const auto y = globalPos[1];
        return (std::pow(x, 2) - x - 2.0 * y + std::pow(y, 4)) *
               std::cos(M_PI * time);
    };

    Scalar eps_;
    std::string problemName_;
    Scalar time_;
    //std::vector<Scalar> analyticalVelocityX_;
    //std::vector<Scalar> analyticalVelocityY_;
    //std::vector<Scalar> analyticalPressure_;

    std::shared_ptr<CouplingManager> couplingManager_;
};
}  // end namespace Dumux

#endif  // DUMUX_STOKES_SUBPROBLEM_HH
