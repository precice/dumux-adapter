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
 * \brief The coupled exercise properties file or the interface case.
 */
#ifndef DUMUX_EXERCISE_COUPLED_INTERFACE_PROPERTIES_HH
#define DUMUX_EXERCISE_COUPLED_INTERFACE_PROPERTIES_HH

// Both domains
#include <dumux/multidomain/boundary/stokesdarcy/couplingmanager.hh>
#include <dumux/multidomain/staggeredtraits.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/io/grid/gridmanager_sub.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

// Porous medium flow domain
#include <dumux/discretization/cctpfa.hh>
#include <dumux/porousmediumflow/1p/model.hh>

#include "1pspatialparams.hh"
#include "porousmediumsubproblem.hh"

// Free-flow domain
#include <dumux/discretization/staggered/freeflow/properties.hh>
#include <dumux/freeflow/navierstokes/model.hh>

#include "freeflowsubproblem.hh"

namespace Dumux::Properties
{
// Create new type tags
namespace TTag
{
struct DarcyOneP {
    using InheritsFrom = std::tuple<OneP, CCTpfaModel>;
};
struct StokesOneP {
    using InheritsFrom = std::tuple<NavierStokes, StaggeredFreeFlowModel>;
};
}  // end namespace TTag

// Set the coupling manager
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::StokesOneP> {
    using Traits = StaggeredMultiDomainTraits<TypeTag,
                                              TypeTag,
                                              Properties::TTag::DarcyOneP>;
    using type = Dumux::StokesDarcyCouplingManager<Traits>;
};
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::DarcyOneP> {
    using Traits = StaggeredMultiDomainTraits<Properties::TTag::StokesOneP,
                                              Properties::TTag::StokesOneP,
                                              TypeTag>;
    using type = Dumux::StokesDarcyCouplingManager<Traits>;
};

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::DarcyOneP> {
    using type = Dumux::PorousMediumSubProblem<TypeTag>;
};
template<class TypeTag>
struct Problem<TypeTag, TTag::StokesOneP> {
    using type = Dumux::FreeFlowSubProblem<TypeTag>;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::DarcyOneP> {
    static constexpr auto dim = 2;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using TensorGrid =
        Dune::YaspGrid<2, Dune::TensorProductCoordinates<Scalar, dim>>;

    using HostGrid = TensorGrid;
    using type = Dune::SubGrid<dim, HostGrid>;
};
template<class TypeTag>
struct Grid<TypeTag, TTag::StokesOneP> {
    static constexpr auto dim = 2;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using TensorGrid =
        Dune::YaspGrid<2, Dune::TensorProductCoordinates<Scalar, dim>>;

    using HostGrid = TensorGrid;
    using type = Dune::SubGrid<dim, HostGrid>;
};

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::DarcyOneP> {
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type =
        FluidSystems::OnePLiquid<Scalar, Dumux::Components::SimpleH2O<Scalar>>;
};
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::StokesOneP> {
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type =
        FluidSystems::OnePLiquid<Scalar, Dumux::Components::SimpleH2O<Scalar>>;
};

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::DarcyOneP> {
    using type = OnePSpatialParams<GetPropType<TypeTag, GridGeometry>,
                                   GetPropType<TypeTag, Scalar>>;
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

}  // end namespace Dumux::Properties

#endif
