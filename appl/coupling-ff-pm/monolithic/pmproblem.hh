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
 *
 * \brief The porous medium flow sub problem
 */
#ifndef DUMUX_DARCY_SUBPROBLEM_HH
#define DUMUX_DARCY_SUBPROBLEM_HH

#ifndef ENABLEMONOLITHIC
#define ENABLEMONOLITHIC 1
#endif

#include <dune/grid/yaspgrid.hh>

//****** uncomment for the last exercise *****//
// #include <dumux/io/grid/subgridgridcreator.hh>

#include <dumux/discretization/cctpfa.hh>

#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/problem.hh>

#include "1pspatialparams.hh"

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include "../precice/preciceadapter.hh"

namespace Dumux
{
template <class TypeTag>
class DarcySubProblem;

namespace Properties
{
// Create new type tags
namespace TTag {
struct DarcyOneP { using InheritsFrom = std::tuple<OneP, CCTpfaModel>; };
} // end namespace TTag

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::DarcyOneP> { using type = Dumux::DarcySubProblem<TypeTag>; };

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::DarcyOneP>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Dumux::Components::SimpleH2O<Scalar> > ;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::DarcyOneP>
{
    static constexpr auto dim = 2;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using TensorGrid = Dune::YaspGrid<2, Dune::TensorProductCoordinates<Scalar, dim> >;

//****** comment out for the last exercise *****//
    using type = TensorGrid;

//****** uncomment for the last exercise *****//
    // using HostGrid = TensorGrid;
    // using type = Dune::SubGrid<dim, HostGrid>;
};

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::DarcyOneP> {
    using type = OnePSpatialParams<GetPropType<TypeTag, FVGridGeometry>, GetPropType<TypeTag, Scalar>>;
};

} // end namespace Properties

/*!
 * \brief The porous medium flow sub problem
 */
template <class TypeTag>
class DarcySubProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::FVGridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;

    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

#if ENABLEMONOLITHIC
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
#endif

public:
#if ENABLEMONOLITHIC
    DarcySubProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                   std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(fvGridGeometry, "Darcy"), eps_(1e-7), couplingManager_(couplingManager)
#else
DarcySubProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry, "Darcy"), eps_(1e-7),
      couplingInterface_(precice_adapter::PreciceAdapter::getInstance() ),
      pressureId_(couplingInterface_.getIdFromName( "Pressure" ) ),
      velocityId_(couplingInterface_.getIdFromName( "Velocity" ) )
#endif
    {}

    /*!
     * \name Simulation steering
     */
    // \{

    /*!
     * \brief Return the temperature within the domain in [K].
     *
     */
    Scalar temperature() const
    { return 273.15 + 10; } // 10°C
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
    BoundaryTypes boundaryTypes(const Element &element, const SubControlVolumeFace &scvf) const
    {
        BoundaryTypes values;

        // set Neumann BCs to all boundaries first
        values.setAllNeumann();

#if ENABLEMONOLITHIC
        // set the coupling boundary condition at the interface
        if (couplingManager().isCoupledEntity(CouplingManager::darcyIdx, scvf))
            values.setAllCouplingNeumann();
#endif
        return values;
    }

        /*!
     * \brief Evaluate the boundary conditions for a Dirichlet control volume.
     *
     * \param element The element for which the Dirichlet boundary condition is set
     * \param scvf The boundary subcontrolvolumeface
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichlet(const Element &element, const SubControlVolumeFace &scvf) const
    {
        // set p = 0 at the bottom
        PrimaryVariables values(0.0);
        values = initial(element);

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a Neumann control volume.
     *
     * \param element The element for which the Neumann boundary condition is set
     * \param fvGeomentry The fvGeometry
     * \param elemVolVars The element volume variables
     * \param scvf The boundary sub control volume face
     *
     * For this method, the \a values variable stores primary variables.
     */
    template<class ElementVolumeVariables>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const SubControlVolumeFace& scvf) const
    {
        // no-flow everywhere ...
        NumEqVector values(0.0);

#if ENABLEMONOLITHIC
        // ... except at the coupling interface
        if (couplingManager().isCoupledEntity(CouplingManager::darcyIdx, scvf))
            values[Indices::conti0EqIdx] = couplingManager().couplingData().massCouplingCondition(element, fvGeometry, elemVolVars, scvf);
#else
        const auto faceId = scvf.index();
        if ( couplingInterface_.isCoupledEntity(faceId) )
        {
          values[Indices::conti0EqIdx] = - couplingInterface_.getQuantityOnFace( velocityId_, faceId );
        }
        // if (/*preCICE*/)
        // {
        //     values[Indices::conti0EqIdx] = /*mass flux from stokes*/;
        // }
#endif
        return values;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{
    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * \param element The element for which the source term is set
     * \param fvGeomentry The fvGeometry
     * \param elemVolVars The element volume variables
     * \param scv The subcontrolvolume
     */
    template<class ElementVolumeVariables>
    NumEqVector source(const Element &element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume &scv) const
    { return NumEqVector(0.0); }

    // \}

    /*!
     * \brief Evaluate the initial value for a control volume.
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

    // \}

#if ENABLEMONOLITHIC
    //! Get the coupling manager
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }
#endif

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

#if ENABLEMONOLITHIC
    std::shared_ptr<CouplingManager> couplingManager_;
#else
   precice_adapter::PreciceAdapter& couplingInterface_;
   const size_t pressureId_;
   const size_t velocityId_;
#endif
};
} //end namespace

#endif //DUMUX_DARCY_SUBPROBLEM_HH
