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

#ifndef ENABLEMONOLITHIC
#define ENABLEMONOLITHIC 0
#endif

#include <dune/grid/yaspgrid.hh>

#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/simpleh2o.hh>

#include <dumux/freeflow/navierstokes/problem.hh>
#include <dumux/discretization/staggered/freeflow/properties.hh>
#include <dumux/freeflow/navierstokes/model.hh>

#include "../../precice-adapter/include/preciceadapter.hh"

namespace Dumux
{
template <class TypeTag>
class StokesSubProblem;

namespace Properties
{
// Create new type tags
namespace TTag {
struct FreeFlowModel { using InheritsFrom = std::tuple<NavierStokes, StaggeredFreeFlowModel>; };
} // end namespace TTag

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::FreeFlowModel>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Dumux::Components::SimpleH2O<Scalar> > ;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::FreeFlowModel> { using type = Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2> >; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::FreeFlowModel> { using type = Dumux::StokesSubProblem<TypeTag> ; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::FreeFlowModel> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::FreeFlowModel> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::FreeFlowModel> { static constexpr bool value = true; };
} // end namespace Properties

/*!
 * \brief The free flow sub problem
 */
template <class TypeTag>
class StokesSubProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;

    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

#if ENABLEMONOLITHIC
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
#endif

public:
#if ENABLEMONOLITHIC
    StokesSubProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, "Stokes"), eps_(1e-6), couplingManager_(couplingManager)
#else
    StokesSubProblem(std::shared_ptr<const GridGeometry> gridGeometry)
      : ParentType(gridGeometry, "FreeFlow"),
        eps_(1e-6),
        couplingInterface_(precice_adapter::PreciceAdapter::getInstance() ),
        pressureId_(0),
        velocityId_(0),
        dataIdsWereSet_(false)
#endif
    {
        deltaP_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.PressureDifference");
//        pressureId_ =  couplingInterface_.getIdFromName( "Pressure" );
//        velocityId_ = couplingInterface_.getIdFromName( "Velocity" );
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
    Scalar temperature() const
    { return 273.15 + 10; } // 10°C

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

        const auto& globalPos = scvf.dofPosition();
        const auto faceId = scvf.index();

        // left/right wall
        if (onRightBoundary_(globalPos) || (onLeftBoundary_(globalPos)))
        {
          values.setDirichlet(Indices::pressureIdx);
        }


        // coupling interface
#if ENABLEMONOLITHIC
        else if(couplingManager().isCoupledEntity(CouplingManager::stokesIdx, scvf))
        {
            values.setCouplingNeumann(Indices::conti0EqIdx);
            values.setCouplingNeumann(Indices::momentumYBalanceIdx);
            values.setBeaversJoseph(Indices::momentumXBalanceIdx);
        }
#else

        else if ( couplingInterface_.isCoupledEntity(faceId) )
        {
        // // TODO do preCICE stuff in analogy to heat transfer
            assert( dataIdsWereSet_ );
          //TODO What do I want to do here?
        //  values.setCouplingNeumann(Indices::conti0EqIdx);
        //  values.setCouplingNeumann(Indices::momentumYBalanceIdx);
          values.setDirichlet(Indices::velocityYIdx);

//          values.setNeumann(Indices::conti0EqIdx);
//          values.setNeumann(Indices::momentumYBalanceIdx);
          values.setBeaversJoseph(Indices::momentumXBalanceIdx);
        }
#endif
        else
        {
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
    PrimaryVariables dirichlet(const Element& element, const SubControlVolumeFace& scvf) const
    {
        PrimaryVariables values(0.0);
        values = initialAtPos(scvf.center());

        const auto faceId = scvf.index();
        if( couplingInterface_.isCoupledEntity( faceId ) )
        {
          values[Indices::velocityYIdx] =
              couplingInterface_.getScalarQuantityOnFace( velocityId_, faceId );
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
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFaceVariables& elemFaceVars,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);

#if ENABLEMONOLITHIC
        if(couplingManager().isCoupledEntity(CouplingManager::stokesIdx, scvf))
        {
            values[Indices::conti0EqIdx] = couplingManager().couplingData().massCouplingCondition(element, fvGeometry, elemVolVars, elemFaceVars, scvf);
            values[Indices::momentumYBalanceIdx] = couplingManager().couplingData().momentumCouplingCondition(element, fvGeometry, elemVolVars, elemFaceVars, scvf);
        }
#else
        assert( dataIdsWereSet_ );
        const auto faceId = scvf.index();
        if( couplingInterface_.isCoupledEntity( faceId ) )
        {
          const Scalar density = 1000; // TODO how to handle compressible fluids?
          const auto& volVars = elemVolVars[scvf.insideScvIdx()];
          const Scalar density_ = volVars.density();
          values[Indices::conti0EqIdx] = density * elemFaceVars[scvf].velocitySelf() * scvf.directionSign();
          values[Indices::momentumYBalanceIdx] = scvf.directionSign() * (couplingInterface_.getScalarQuantityOnFace( pressureId_, faceId ) - initialAtPos(scvf.center())[Indices::pressureIdx]) ;
        }
#endif

        return values;
    }

    // \}

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
     * \brief Evaluate the initial value for a control volume.
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);
        //values[Indices::velocityYIdx] = -1e-6 * globalPos[0] * (this->gridGeometry().bBoxMax()[0] - globalPos[0]);
        if(onLeftBoundary_(globalPos))
          values[Indices::pressureIdx] = deltaP_;
        if(onRightBoundary_(globalPos))
          values[Indices::pressureIdx] = 0.0;

        return values;
    }

    /*!
     * \brief Returns the intrinsic permeability of required as input parameter for the Beavers-Joseph-Saffman boundary condition
     */
    Scalar permeability(const Element& element, const SubControlVolumeFace& scvf) const
    {
#if ENABLEMONOLITHIC
        return couplingManager().couplingData().darcyPermeability(element, scvf);
#else
        return 1e-10; // TODO transfer information or just use constant value
#endif
    }

    /*!
     * \brief Returns the alpha value required as input parameter for the Beavers-Joseph-Saffman boundary condition
     */
    Scalar alphaBJ(const SubControlVolumeFace& scvf) const
    {
#if ENABLEMONOLITHIC
        return couplingManager().problem(CouplingManager::darcyIdx).spatialParams().beaversJosephCoeffAtPos(scvf.center());
#else
        return 1.0; // TODO transfer information or just use constant value
#endif
    }

    /*!
     * \brief calculate the analytical velocity in x direction based on Beavers & Joseph (1967)
     */
    void calculateAnalyticalVelocityX() const
    {
        analyticalVelocityX_.resize(this->gridGeometry().gridView().size(0));

        using std::sqrt;
        const Scalar dPdX = -deltaP_ / (this->gridGeometry().bBoxMax()[0] - this->gridGeometry().bBoxMin()[0]);
        static const Scalar mu = FluidSystem::viscosity(temperature(), 1e5);
        static const Scalar alpha = getParam<Scalar>("Darcy.SpatialParams.AlphaBeaversJoseph");
        static const Scalar K = getParam<Scalar>("Darcy.SpatialParams.Permeability");
        static const Scalar sqrtK = sqrt(K);
        const Scalar sigma = (this->gridGeometry().bBoxMax()[1] - this->gridGeometry().bBoxMin()[1])/sqrtK;

        const Scalar uB =  -K/(2.0*mu) * ((sigma*sigma + 2.0*alpha*sigma) / (1.0 + alpha*sigma)) * dPdX;

        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            const auto eIdx = this->gridGeometry().gridView().indexSet().index(element);
            const Scalar y = element.geometry().center()[1] - this->gridGeometry().bBoxMin()[1];

            const Scalar u = uB*(1.0 + alpha/sqrtK*y) + 1.0/(2.0*mu) * (y*y + 2*alpha*y*sqrtK) * dPdX;
            analyticalVelocityX_[eIdx] = u;
        }
    }

    /*!
     * \brief Get the analytical velocity in x direction
     */
    const std::vector<Scalar>& getAnalyticalVelocityX() const
    {
        if(analyticalVelocityX_.empty())
            calculateAnalyticalVelocityX();
        return analyticalVelocityX_;
    }

#if !ENABLEMONOLITHIC
    void updatePreciceDataIds()
    {
      pressureId_ = couplingInterface_.getIdFromName( "Pressure" );
      velocityId_ = couplingInterface_.getIdFromName( "Velocity" );
      dataIdsWereSet_ = true;
    }
#endif

    // \}

private:
    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] < this->gridGeometry().bBoxMin()[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_; }

    Scalar eps_;
    Scalar deltaP_;

#if ENABLEMONOLITHIC
    std::shared_ptr<CouplingManager> couplingManager_;
#else
   precice_adapter::PreciceAdapter& couplingInterface_;
   size_t pressureId_;
   size_t velocityId_;
   bool dataIdsWereSet_;
#endif

    mutable std::vector<Scalar> analyticalVelocityX_;
};
} //end namespace

#endif // DUMUX_STOKES_SUBPROBLEM_HH
