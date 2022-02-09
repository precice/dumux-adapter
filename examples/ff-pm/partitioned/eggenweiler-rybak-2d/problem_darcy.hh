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
 * \brief The Darcy sub-problem of coupled Stokes-Darcy convergence test
 */

#ifndef DUMUX_DARCY_SUBPROBLEM_HH
#define DUMUX_DARCY_SUBPROBLEM_HH

#include <dumux/common/numeqvector.hh>
#include <dune/common/fvector.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/box.hh>
#include <dumux/discretization/cctpfa.hh>

#include <dumux/common/boundarytypes.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/problem.hh>

#include "spatialparams.hh"

#include <dumux-precice/couplingadapter.hh>

namespace Dumux
{
template<class TypeTag>
class DarcySubProblem;

/*!
 * \ingroup BoundaryTests
 * \brief The Darcy sub-problem of coupled Stokes-Darcy convergence test
 */
template<class TypeTag>
class DarcySubProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView =
        typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<
        GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using FVElementGeometry =
        typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace =
        typename FVElementGeometry::SubControlVolumeFace;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    static constexpr auto velocityXIdx = 0;
    static constexpr auto velocityYIdx = 1;
    static constexpr auto pressureIdx = 2;

public:
    //! export the Indices
    using Indices =
        typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    DarcySubProblem(std::shared_ptr<const GridGeometry> gridGeometry)
        : ParentType(gridGeometry, "Darcy"),
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
     */
    Scalar temperature() const { return 273.15 + 10; }  // 10°C
    // \}

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

        if (couplingInterface_.isCoupledEntity(scvf.index()))
            values.setAllCouplingNeumann();
        else
            values.setAllDirichlet();

        return values;
    }

    /*!
      * \brief Specifies which kind of boundary condition should be
      *        used for which equation on a given boundary control volume.
      *
      * \param element The element
      * \param scv The boundary sub control volume
      */
    BoundaryTypes boundaryTypes(const Element &element,
                                const SubControlVolume &scv) const
    {
        BoundaryTypes values;

        values.setAllDirichlet();

        //corners will not be handled as coupled
        if (onLeftBoundary_(scv.dofPosition()) ||
            onRightBoundary_(scv.dofPosition())) {
            values.setAllDirichlet();
            return values;
        }

        auto fvGeometry = localView(this->gridGeometry());
        fvGeometry.bindElement(element);
        for (auto &&scvf : scvfs(fvGeometry)) {
            if (couplingInterface_.isCoupledEntity(scvf.index())) {
                values.setAllCouplingNeumann();
            }
        }
        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet control volume.
     *
     * \param globalPos The position for which the Dirichlet value is set
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        const auto p = analyticalSolution(globalPos)[pressureIdx];
        return PrimaryVariables(p);
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann control volume.
     *
     * \param element The element for which the Neumann boundary condition is set
     * \param fvGeometry The fvGeometry
     * \param elemVolVars The element volume variables
     * \param elemFluxVarsCache Flux variables caches for all faces in stencil
     * \param scvf The boundary sub control volume face
     *
     * For this method, the \a values variable stores primary variables.
     */
    template<class ElementVolumeVariables, class ElementFluxVarsCache>
    NumEqVector neumann(const Element &element,
                        const FVElementGeometry &fvGeometry,
                        const ElementVolumeVariables &elemVolVars,
                        const ElementFluxVarsCache &elemFluxVarsCache,
                        const SubControlVolumeFace &scvf) const
    {
        NumEqVector values(0.0);

        if (couplingInterface_.isCoupledEntity(scvf.index()))
            //values[Indices::conti0EqIdx] = couplingInterface_.getMassCouplingCondition( scvf );
            values[Indices::conti0EqIdx] =
                couplingInterface_.getScalarQuantityOnFace(velocityId_,
                                                           scvf.index());
        ;
        // couplingManager().couplingData().massCouplingCondition(
        //     element, fvGeometry, elemVolVars, elemFluxVarsCache, scvf);

        return values;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{
    /*!
     * \brief Evaluates the source term for all phases within a given
     *        sub control volume.
     *
     * \param globalPos The global position
     */
    NumEqVector sourceAtPos(const GlobalPosition &globalPos) const
    {
        return rhsNewICNonSymmetrized_(globalPos);
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
    PrimaryVariables initial(const Element &element) const
    {
        return PrimaryVariables(0.0);
    }

    /*!
     * \brief Returns the analytical solution of the problem at a given position.
     *
     * \param globalPos The global position
     */
    auto analyticalSolution(const GlobalPosition &globalPos) const
    {
        return analyticalSolutionNewICNonSymmetrized_(globalPos);
    }

    // \}

private:
    // see Rybak et al., 2015: "Multirate time integration for coupled saturated/unsaturated porous medium and free flow systems"
    Dune::FieldVector<Scalar, 3> analyticalSolutionRybak_(
        const GlobalPosition &globalPos) const
    {
        Dune::FieldVector<Scalar, 3> sol(0.0);
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];

        using std::cos;
        using std::exp;
        using std::sin;
        sol[velocityXIdx] = -0.5 * M_PI * y * y * cos(M_PI * x);
        sol[velocityYIdx] = -1.0 * y * sin(M_PI * x);
        sol[pressureIdx] = 0.5 * y * y * sin(M_PI * x);
        return sol;
    }

    // see Rybak et al., 2015: "Multirate time integration for coupled saturated/unsaturated porous medium and free flow systems"
    NumEqVector rhsRybak_(const GlobalPosition &globalPos) const
    {
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];
        using std::sin;
        return NumEqVector((0.5 * M_PI * y * M_PI * y - 1) * sin(M_PI * x));
    }

    // exact solution for new IC with non-symmetrized stress tensor (by Elissa Eggenweiler)
    Dune::FieldVector<Scalar, 3> analyticalSolutionNewICNonSymmetrized_(
        const GlobalPosition &globalPos) const
    {
        Dune::FieldVector<Scalar, 3> sol(0.0);
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];

        using std::cos;
        using std::exp;
        using std::sin;
        sol[velocityXIdx] = x * (y - 1.0) + (x - 1.0) * (y - 1.0) - 3.0;
        sol[velocityYIdx] = x * (x - 1.0) - (y - 1.0) * (y - 1.0) - 2.0;
        sol[pressureIdx] = x * (1.0 - x) * (y - 1.0) +
                           1.0 / 3.0 * (y - 1.0) * (y - 1.0) * (y - 1.0) +
                           3.0 * x + 2.0 * y + 1.0;
        return sol;
    }

    // exact solution for new IC with non-symmetrized stress tensor (by Elissa Eggenweiler)
    NumEqVector rhsNewICNonSymmetrized_(const GlobalPosition &globalPos) const
    {
        return NumEqVector(0.0);
    }

    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_;
    }
    bool onRightBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_;
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
