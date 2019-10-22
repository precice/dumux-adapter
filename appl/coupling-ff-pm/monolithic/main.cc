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
 * \brief A test problem for the coupled Stokes/Darcy problem (1p)
 */
#include <config.h>

#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/istl/io.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/partial.hh>
#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/discretization/method.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/staggeredvtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager.hh>

#include <dumux/multidomain/staggeredtraits.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/newtonsolver.hh>

#include <dumux/multidomain/boundary/stokesdarcy/couplingmanager.hh>

#include "pmproblem.hh"
#include "ffproblem.hh"

namespace Dumux {
namespace Properties {

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::FreeFlowModel>
{
    using Traits = StaggeredMultiDomainTraits<TypeTag, TypeTag, Properties::TTag::DarcyOneP>;
    using type = Dumux::StokesDarcyCouplingManager<Traits>;
};

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::DarcyOneP>
{
    using Traits = StaggeredMultiDomainTraits<Properties::TTag::FreeFlowModel, Properties::TTag::FreeFlowModel, TypeTag>;
    using type = Dumux::StokesDarcyCouplingManager<Traits>;
};

} // end namespace Properties
} // end namespace Dumux


/*!
 * \brief Returns the velocity at the interface using Darcy's law for reconstruction
 */
template<class FluxVariables, class Problem, class Element, class FVElementGeometry,
         class ElementVolumeVariables, class SubControlVolumeFace,
         class ElementFluxVariablesCache>
auto velocityAtInterface(const Problem& problem,
                           const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const SubControlVolumeFace& scvf,
                           const ElementFluxVariablesCache& elemFluxVarsCache)
{
  const int phaseIdx = 0;
  FluxVariables fluxVars;
  fluxVars.init(problem, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);
  auto upwindTerm = [phaseIdx](const auto& volVars) { return volVars.mobility(phaseIdx); };
  const auto scalarVelocity = fluxVars.advectiveFlux(phaseIdx, upwindTerm)/scvf.area();
  return scalarVelocity;
}

template<class FluxVariables, class CouplingManager, class Problem, class GridVariables, class SolutionVector>
void writeVelocitiesOnInterfaceToFile( const std::string& filename,
                                       const CouplingManager& couplingManager,
                                       const Problem& problem,
                                       const GridVariables& gridVars,
                                       const SolutionVector& sol)
{
  const auto& fvGridGeometry = problem.fvGridGeometry();
  auto fvGeometry = localView(fvGridGeometry);
  auto elemVolVars = localView(gridVars.curGridVolVars());
  auto elemFluxVarsCache = localView(gridVars.gridFluxVarsCache());



  std::ofstream ofs( filename+".csv", std::ofstream::out | std::ofstream::trunc);
  ofs << "x,y,";
  ofs << "velocityY" << "\n";
  for (const auto& element : elements(fvGridGeometry.gridView()))
  {
    fvGeometry.bind(element);
    elemVolVars.bind(element, fvGeometry, sol);
    elemFluxVarsCache.bind(element, fvGeometry, elemVolVars);

    for (const auto& scvf : scvfs(fvGeometry))
    {

      if (couplingManager.isCoupledEntity(CouplingManager::darcyIdx, scvf))
      {
        const auto& pos = scvf.center();
        for (int i = 0; i < 2; ++i )
        {
          ofs << pos[i] << ",";
        }
        const double v = couplingManager.couplingData().massCouplingCondition(element, fvGeometry, elemVolVars, scvf);
        ofs << v / 1e3 << "\n";
      }
    }
  }

  ofs.close();
}

template<class CouplingManager, class Problem, class SolutionVector>
void writeStokesVelocitiesOnInterfaceToFile( const std::string& filename,
                                       const CouplingManager& couplingManager,
                                       const Problem& problem,
                                       const SolutionVector& sol)
{
  const auto& fvGridGeometry = problem.fvGridGeometry();
  auto fvGeometry = localView(fvGridGeometry);
  using FVGridGeometry = std::decay_t<decltype (fvGridGeometry)>;

  std::ofstream ofs( filename+".csv", std::ofstream::out | std::ofstream::trunc);
  ofs << "x,y,";
  ofs << "velocityY" << "\n";
  for (const auto& element : elements(fvGridGeometry.gridView()))
  {
    fvGeometry.bind(element);

    for (const auto& scvf : scvfs(fvGeometry))
    {

      if (couplingManager.isCoupledEntity(CouplingManager::stokesIdx, scvf))
      {
        const auto& pos = scvf.center();
        for (int i = 0; i < 2; ++i )
        {
          ofs << pos[i] << ",";
        }
        const double v = sol[scvf.dofIndex()];
        ofs << v << "\n";
      }
    }
  }

  ofs.close();
}

int main(int argc, char** argv) try
{
    using namespace Dumux;

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // Define the sub problem type tags
    using StokesTypeTag = Properties::TTag::FreeFlowModel;
    using DarcyTypeTag = Properties::TTag::DarcyOneP;



    // ******************** comment-out this section for the last exercise **************** //

    // create two individual grids (from the given grid file or the input file)
    // for both sub-domains
    using DarcyGridManager = Dumux::GridManager<GetPropType<DarcyTypeTag, Properties::Grid>>;
    DarcyGridManager darcyGridManager;
    darcyGridManager.init("Darcy"); // pass parameter group

    using StokesGridManager = Dumux::GridManager<GetPropType<StokesTypeTag, Properties::Grid>>;
    StokesGridManager stokesGridManager;
    stokesGridManager.init("Stokes"); // pass parameter group

    // we compute on the leaf grid view
    const auto& darcyGridView = darcyGridManager.grid().leafGridView();
    const auto& stokesGridView = stokesGridManager.grid().leafGridView();

    // ************************************************************************************ //


    // ******************** uncomment this section for the last exercise ****************** //

    // // use dune-subgrid to create the individual grids
    // static constexpr int dim = 2;
    // using HostGrid = Dune::YaspGrid<2, Dune::TensorProductCoordinates<double, dim> >;
    // using HostGridManager = Dumux::GridManager<HostGrid>;
    // HostGridManager hostGridManager;
    // hostGridManager.init();
    // auto& hostGrid = hostGridManager.grid();
    //
    // struct Params
    // {
    //     double amplitude = getParam<double>("Grid.Amplitude");
    //     double baseline = getParam<double>("Grid.Baseline");
    //     double offset = getParam<double>("Grid.Offset");
    //     double scaling = getParam<double>("Grid.Scaling");
    // };
    //
    // Params params;
    //
    // auto elementSelectorStokes = [&](const auto& element)
    // {
    //     double interface = params.amplitude * std::sin(( element.geometry().center()[0] -params.offset) / params.scaling * 2.0 * M_PI) + params.baseline;
    //     return element.geometry().center()[1] > interface;
    // };
    //
    // auto elementSelectorDarcy = [&](const auto& element)
    // {
    //     double interface  =  params.amplitude * std::sin(( element.geometry().center()[0] - params.offset) / params.scaling * 2.0 * M_PI) + params.baseline;
    //     return element.geometry().center()[1] < interface;
    // };
    //
    // // subgrid Pointer
    // auto stokesGridPtr = SubgridGridCreator<HostGrid>::makeGrid(hostGrid, elementSelectorStokes, "Stokes");
    // auto darcyGridPtr = SubgridGridCreator<HostGrid>::makeGrid(hostGrid, elementSelectorDarcy, "Darcy");
    //
    // // we compute on the leaf grid view
    // const auto& darcyGridView = darcyGridPtr->leafGridView();
    // const auto& stokesGridView = stokesGridPtr->leafGridView();

    // ************************************************************************************ //


    // create the finite volume grid geometry
    using StokesFVGridGeometry = GetPropType<StokesTypeTag, Properties::FVGridGeometry>;
    auto stokesFvGridGeometry = std::make_shared<StokesFVGridGeometry>(stokesGridView);
    stokesFvGridGeometry->update();
    using DarcyFVGridGeometry = GetPropType<DarcyTypeTag, Properties::FVGridGeometry>;
    auto darcyFvGridGeometry = std::make_shared<DarcyFVGridGeometry>(darcyGridView);
    darcyFvGridGeometry->update();

    using Traits = StaggeredMultiDomainTraits<StokesTypeTag, StokesTypeTag, DarcyTypeTag>;

    // the coupling manager
    using CouplingManager = StokesDarcyCouplingManager<Traits>;
    auto couplingManager = std::make_shared<CouplingManager>(stokesFvGridGeometry, darcyFvGridGeometry);

    // the indices
    constexpr auto stokesCellCenterIdx = CouplingManager::stokesCellCenterIdx;
    constexpr auto stokesFaceIdx = CouplingManager::stokesFaceIdx;
    constexpr auto darcyIdx = CouplingManager::darcyIdx;

    // the problem (initial and boundary conditions)
    using StokesProblem = GetPropType<StokesTypeTag, Properties::Problem>;
    auto stokesProblem = std::make_shared<StokesProblem>(stokesFvGridGeometry, couplingManager);
    using DarcyProblem = GetPropType<DarcyTypeTag, Properties::Problem>;
    auto darcyProblem = std::make_shared<DarcyProblem>(darcyFvGridGeometry, couplingManager);

    // the solution vector
    Traits::SolutionVector sol;
    sol[stokesCellCenterIdx].resize(stokesFvGridGeometry->numCellCenterDofs());
    sol[stokesFaceIdx].resize(stokesFvGridGeometry->numFaceDofs());
    sol[darcyIdx].resize(darcyFvGridGeometry->numDofs());

    auto stokesSol = partial(sol, stokesCellCenterIdx, stokesFaceIdx);

    stokesProblem->applyInitialSolution(stokesSol);
    darcyProblem->applyInitialSolution(sol[darcyIdx]);

    couplingManager->init(stokesProblem, darcyProblem, sol);

    // the grid variables
    using StokesGridVariables = GetPropType<StokesTypeTag, Properties::GridVariables>;
    auto stokesGridVariables = std::make_shared<StokesGridVariables>(stokesProblem, stokesFvGridGeometry);
    stokesGridVariables->init(stokesSol);
    using DarcyGridVariables = GetPropType<DarcyTypeTag, Properties::GridVariables>;
    auto darcyGridVariables = std::make_shared<DarcyGridVariables>(darcyProblem, darcyFvGridGeometry);
    darcyGridVariables->init(sol[darcyIdx]);

    // intialize the vtk output module
    const auto stokesName = getParam<std::string>("Problem.Name") + "_" + stokesProblem->name();
    const auto darcyName = getParam<std::string>("Problem.Name") + "_" + darcyProblem->name();

    StaggeredVtkOutputModule<StokesGridVariables, decltype(stokesSol)> stokesVtkWriter(*stokesGridVariables, stokesSol, stokesName);
    GetPropType<StokesTypeTag, Properties::IOFields>::initOutputModule(stokesVtkWriter);

    //****** uncomment the add analytical solution of v_x *****//
    stokesVtkWriter.addField(stokesProblem->getAnalyticalVelocityX(), "analyticalV_x");

    stokesVtkWriter.write(0.0);

    VtkOutputModule<DarcyGridVariables, GetPropType<DarcyTypeTag, Properties::SolutionVector>> darcyVtkWriter(*darcyGridVariables, sol[darcyIdx], darcyName);
    using DarcyVelocityOutput = GetPropType<DarcyTypeTag, Properties::VelocityOutput>;
    darcyVtkWriter.addVelocityOutput(std::make_shared<DarcyVelocityOutput>(*darcyGridVariables));
    GetPropType<DarcyTypeTag, Properties::IOFields>::initOutputModule(darcyVtkWriter);
    darcyVtkWriter.write(0.0);

    // the assembler for a stationary problem
    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(std::make_tuple(stokesProblem, stokesProblem, darcyProblem),
                                                 std::make_tuple(stokesFvGridGeometry->cellCenterFVGridGeometryPtr(),
                                                                 stokesFvGridGeometry->faceFVGridGeometryPtr(),
                                                                 darcyFvGridGeometry),
                                                 std::make_tuple(stokesGridVariables->cellCenterGridVariablesPtr(),
                                                                 stokesGridVariables->faceGridVariablesPtr(),
                                                                 darcyGridVariables),
                                                 couplingManager);

    // the linear solver
    using LinearSolver = UMFPackBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);

    // solve the non-linear system
    nonLinearSolver.solve(sol);

    // write vtk output
    stokesVtkWriter.write(1.0);
    darcyVtkWriter.write(1.0);

    {
      using FluxVariables = GetPropType<DarcyTypeTag, Properties::FluxVariables>;
      const std::string filename = getParam<std::string>("Problem.Name") + "-" + darcyProblem->name() + "-interface-velocity";
      writeVelocitiesOnInterfaceToFile<FluxVariables>( filename,
                                                       *couplingManager,
                                                       *darcyProblem,
                                                       *darcyGridVariables,
                                                       sol[darcyIdx] );
    }

    //TODO make freeflow
    {
      const std::string filename = getParam<std::string>("Problem.Name") + "-" + stokesProblem->name() + "-interface-velocity";
      writeStokesVelocitiesOnInterfaceToFile( filename,
                                              *couplingManager,
                                              *stokesProblem,
                                              sol[stokesFaceIdx] );
    }


    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////

    // print dumux end message
    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }

    return 0;
} // end main
catch (Dumux::ParameterException &e)
{
    std::cerr << std::endl << e << " ---> Abort!" << std::endl;
    return 1;
}
catch (Dune::DGFException & e)
{
    std::cerr << "DGF exception thrown (" << e <<
                 "). Most likely, the DGF file name is wrong "
                 "or the DGF file is corrupted, "
                 "e.g. missing hash at end of file or wrong number (dimensions) of entries."
                 << " ---> Abort!" << std::endl;
    return 2;
}
catch (Dune::Exception &e)
{
    std::cerr << "Dune reported error: " << e << " ---> Abort!" << std::endl;
    return 3;
}
catch (...)
{
    std::cerr << "Unknown exception thrown! ---> Abort!" << std::endl;
    return 4;
}
