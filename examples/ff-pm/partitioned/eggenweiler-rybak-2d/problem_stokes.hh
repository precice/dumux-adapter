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
 * \brief The free-flow sub-problem of coupled FreeFlow/Darcy convergence test
 */

#ifndef DUMUX_FREEFLOW_SUBPROBLEM_HH
#define DUMUX_FREEFLOW_SUBPROBLEM_HH

#include <dumux/common/numeqvector.hh>
#include <dune/common/fvector.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/staggered/freeflow/properties.hh>

#include <dumux/freeflow/navierstokes/boundarytypes.hh>
#include <dumux/freeflow/navierstokes/model.hh>
#include <dumux/freeflow/navierstokes/problem.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

namespace Dumux
{
template<class TypeTag>
class FreeFlowSubProblem;

/*!
 * \ingroup BoundaryTests
 * \brief The Stokes sub-problem of coupled Stokes-Darcy convergence test
 */
template<class TypeTag>
class FreeFlowSubProblem : public NavierStokesProblem<TypeTag>
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
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;

public:
    //! export the Indices
    using Indices = typename ModelTraits::Indices;

    FreeFlowSubProblem(std::shared_ptr<const GridGeometry> gridGeometry)
        : ParentType(gridGeometry, "FreeFlow"),
          couplingInterface_(Dumux::Precice::CouplingAdapter::getInstance()),
          pressureId_(0),
          velocityId_(0),
          dataIdsWereSet_(false)
    {
        problemName_ =
            "EggenweilerRybak_" +
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
        return rhsNewICNonSymmetrized_(globalPos);
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

        values.setDirichlet(Indices::velocityXIdx);
        values.setDirichlet(Indices::velocityYIdx);

        if (couplingInterface_.isCoupledEntity(scvf.index())) {
            values.setCouplingNeumann(Indices::conti0EqIdx);
            values.setCouplingNeumann(Indices::momentumYBalanceIdx);
            values.setSlipCondition(Indices::momentumXBalanceIdx);
        }

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet control volume.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        return analyticalSolution(globalPos);
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

        // if (couplingManager().isCoupledEntity(CouplingManager::freeFlowIdx,
        //                                       scvf)) {
        if (couplingInterface_.isCoupledEntity(scvf.index())) {
            //TODO
            // values[Indices::conti0EqIdx] =
            //     couplingManager().couplingData().massCouplingCondition(
            //         element, fvGeometry, elemVolVars, elemFaceVars, scvf);
            // values[Indices::momentumYBalanceIdx] =
            //     couplingManager().couplingData().momentumCouplingCondition(
            //         element, fvGeometry, elemVolVars, elemFaceVars, scvf);
        }
        return values;
    }

    // \}

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
        return PrimaryVariables(0.0);
    }

    /*!
    * \brief Returns the velocity in the porous medium (which is 0 by default according to Saffmann).
    */
    auto porousMediumTerm(const Element &element,
                          const SubControlVolumeFace &scvf) const
    {
        //TODO
        // Must return a vector
        //return couplingInterface_.getScalarQuantityOnFace(velocityId_, scvf.index());
        using VelocityVector = std::decay_t<decltype(scvf.unitOuterNormal())>;
        return VelocityVector(0.0);
        // return couplingManager().couplingData().porousMediumVelocity(element,
        //                                                             scvf);
    }

    // /*!
    //  * \brief Returns the intrinsic permeability of required as input parameter
    //           for the Beavers-Joseph-Saffman boundary condition
    //  */
    // auto permeability(const Element &element,
    //                   const SubControlVolumeFace &scvf) const
    // {
    //     // return couplingManager().couplingData().darcyPermeability(element,
    //     //                                                           scvf);
    // }

    /*!
     * \brief Returns the alpha value required as input parameter for the
              Beavers-Joseph-Saffman boundary condition.
     */
    Scalar alphaBJ(const SubControlVolumeFace &scvf) const
    {
        //TODO
        return Scalar{0.};
        // return couplingManager()
        //     .problem(CouplingManager::porousMediumIdx)
        //     .spatialParams()
        //     .beaversJosephCoeffAtPos(scvf.center());
    }

    /*!
     * \brief Returns the scale separation parameter epsilon required as input parameter for the
              new coupling conditions
     */
    Scalar epsInterface(const SubControlVolumeFace &scvf) const
    {
        //TODO
        return Scalar{0.};
        // return couplingManager()
        //     .problem(CouplingManager::porousMediumIdx)
        //     .spatialParams()
        //     .epsInterfaceAtPos(scvf.center());
    }

    /*!
     * \brief Returns the boundary layer constant N_1_bl required as input parameter for the
              new coupling condition for the tangential component
     */
    Scalar factorNTangential(const SubControlVolumeFace &scvf) const
    {
        //TODO
        return 0;
        // return couplingManager()
        //     .problem(CouplingManager::porousMediumIdx)
        //     .spatialParams()
        //     .factorNTangentialAtPos(scvf.center());
    }

    /*!
     * \brief Returns the boundary layer constant N_s_bl required as input parameter for the
              new coupling condition for the tangential component
     */
    Scalar factorNMomentum(const SubControlVolumeFace &scvf) const
    {
        //TODO
        return 0;
        // return couplingManager()
        //     .problem(CouplingManager::porousMediumIdx)
        //     .spatialParams()
        //     .factorNMomentumAtPos(scvf.center());
    }

    /*!
     * \brief Returns the boundary layer matrix M_bl required as input parameter for the
              new coupling condition for the tangential component
     */
    auto matrixNTangential(const SubControlVolumeFace &scvf) const
    {
        //TODO
        return 0;
        // return couplingManager()
        //     .problem(CouplingManager::porousMediumIdx)
        //     .spatialParams()
        //     .matrixNTangentialAtPos(scvf.center());
    }

    /*!
     * \brief Returns the analytical solution of the problem at a given position.
     *
     * \param globalPos The global position
     * \param time A parameter for consistent signatures. It is ignored here as this is a stationary test
     */
    PrimaryVariables analyticalSolution(const GlobalPosition &globalPos,
                                        Scalar time = 0.0) const
    {
        return analyticalSolutionNewICNonSymmetrized_(globalPos);
    }

    // \}

private:
    // exact solution for new IC with non-symmetrized stress tensor (by Elissa Eggenweiler)
    PrimaryVariables analyticalSolutionNewICNonSymmetrized_(
        const GlobalPosition &globalPos) const
    {
        PrimaryVariables sol(0.0);
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];

        using std::cos;
        using std::exp;
        using std::sin;
        sol[Indices::velocityXIdx] =
            (y - 1.0) * (y - 1.0) + x * (y - 1.0) + 3.0 * x + 1.5;
        sol[Indices::velocityYIdx] =
            0.5 * x * (x - 1.0) - 0.5 * (y - 1.0) * (y - 1.0) - 3.0 * y + 2.0;
        sol[Indices::pressureIdx] = 2.0 * x + y - 4.0;
        return sol;
    }

    // exact solution for new IC with non-symmetrized stress tensor (by Elissa Eggenweiler)
    NumEqVector rhsNewICNonSymmetrized_(const GlobalPosition &globalPos) const
    {
        using std::cos;
        using std::exp;
        using std::sin;
        NumEqVector source(0.0);
        source[Indices::momentumXBalanceIdx] = -2.0;
        source[Indices::momentumYBalanceIdx] = 1.0;
        return source;
    }

    static constexpr Scalar eps_ = 1e-7;
    Dumux::Precice::CouplingAdapter &couplingInterface_;
    size_t pressureId_;
    size_t velocityId_;
    bool dataIdsWereSet_;
    std::string problemName_;
};
}  // end namespace Dumux

#endif
