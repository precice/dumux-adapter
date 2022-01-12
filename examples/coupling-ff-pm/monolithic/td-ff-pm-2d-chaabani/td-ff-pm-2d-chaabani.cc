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
 * \brief A test problem for the coupled Stokes/Darcy problem (1p).
 */

#include <config.h>

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/istl/io.hh>

#include <dumux/assembly/diffmethod.hh>
#include <dumux/assembly/fvassembler.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/partial.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/method.hh>
#include <dumux/io/grid/gridmanager.hh>
#include <dumux/io/staggeredvtkoutputmodule.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/linear/seqsolverbackend.hh>

#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/newtonsolver.hh>
#include <dumux/multidomain/staggeredtraits.hh>

#include <dumux-precice/dumux-addon/multidomain/boundary/stokesdarcy/couplingmanager.hh>

#include "problem_darcy.hh"
#include "problem_stokes.hh"

namespace Dumux
{
namespace Properties
{
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

}  // end namespace Properties
}  // end namespace Dumux

int main(int argc, char **argv)
try {
    using namespace Dumux;

    // initialize MPI, finalize is done automatically on exit
    const auto &mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // Define the sub problem type tags
    using StokesTypeTag = Properties::TTag::StokesOneP;
    using DarcyTypeTag = Properties::TTag::DarcyOneP;

    // try to create a grid (from the given grid file or the input file)
    // for both sub-domains
    using DarcyGridManager =
        Dumux::GridManager<GetPropType<DarcyTypeTag, Properties::Grid>>;
    DarcyGridManager darcyGridManager;
    darcyGridManager.init("Darcy");  // pass parameter group

    using StokesGridManager =
        Dumux::GridManager<GetPropType<StokesTypeTag, Properties::Grid>>;
    StokesGridManager stokesGridManager;
    stokesGridManager.init("Stokes");  // pass parameter group

    // we compute on the leaf grid view
    const auto &darcyGridView = darcyGridManager.grid().leafGridView();
    const auto &stokesGridView = stokesGridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using StokesGridGeometry =
        GetPropType<StokesTypeTag, Properties::GridGeometry>;
    auto stokesGridGeometry =
        std::make_shared<StokesGridGeometry>(stokesGridView);
    stokesGridGeometry->update();
    using DarcyGridGeometry =
        GetPropType<DarcyTypeTag, Properties::GridGeometry>;
    auto darcyGridGeometry = std::make_shared<DarcyGridGeometry>(darcyGridView);
    darcyGridGeometry->update();

    using Traits =
        StaggeredMultiDomainTraits<StokesTypeTag, StokesTypeTag, DarcyTypeTag>;

    // the coupling manager
    using CouplingManager = StokesDarcyCouplingManager<Traits>;
    auto couplingManager = std::make_shared<CouplingManager>(stokesGridGeometry,
                                                             darcyGridGeometry);

    // the indices
    constexpr auto stokesCellCenterIdx = CouplingManager::stokesCellCenterIdx;
    constexpr auto stokesFaceIdx = CouplingManager::stokesFaceIdx;
    constexpr auto darcyIdx = CouplingManager::darcyIdx;

    // the problem (initial and boundary conditions)
    using StokesProblem = GetPropType<StokesTypeTag, Properties::Problem>;
    auto stokesProblem =
        std::make_shared<StokesProblem>(stokesGridGeometry, couplingManager);
    using DarcyProblem = GetPropType<DarcyTypeTag, Properties::Problem>;
    auto darcyProblem =
        std::make_shared<DarcyProblem>(darcyGridGeometry, couplingManager);

    // the solution vector
    Traits::SolutionVector sol;
    //sol[stokesCellCenterIdx].resize(stokesGridGeometry->numCellCenterDofs());
    //sol[stokesFaceIdx].resize(stokesGridGeometry->numFaceDofs());
    //sol[darcyIdx].resize(darcyGridGeometry->numDofs());
    sol[stokesCellCenterIdx].resize(stokesGridGeometry->numCellCenterDofs());
    sol[stokesFaceIdx].resize(stokesGridGeometry->numFaceDofs());
    sol[darcyIdx].resize(darcyGridGeometry->numDofs());

    // get a solution vector storing references to the two Stokes solution vectors
    auto stokesSol = partial(sol, stokesFaceIdx, stokesCellCenterIdx);

    // the grid variables
    using StokesGridVariables =
        GetPropType<StokesTypeTag, Properties::GridVariables>;
    auto stokesGridVariables = std::make_shared<StokesGridVariables>(
        stokesProblem, stokesGridGeometry);
    stokesGridVariables->init(stokesSol);
    using DarcyGridVariables =
        GetPropType<DarcyTypeTag, Properties::GridVariables>;
    auto darcyGridVariables =
        std::make_shared<DarcyGridVariables>(darcyProblem, darcyGridGeometry);
    darcyGridVariables->init(sol[darcyIdx]);

    /*
  const auto couplingMode = [] {
    const auto mode = getParam<std::string>("Problem.CouplingMode", "ReconstructPorousMediumPressure");
    if (mode == "ReconstructPorousMediumPressure")
      return CouplingManager::CouplingMode::reconstructPorousMediumPressure;
    else if (mode == "ReconstructFreeFlowNormalStress")
      return CouplingManager::CouplingMode::reconstructFreeFlowNormalStress;
    else
      DUNE_THROW(Dune::InvalidStateException, mode << " is not a valid coupling mode. Use ReconstructPorousMediumPressure or ReconstructFreeFlowNormalStress");
  }();
  */

    //couplingManager->init(stokesProblem, darcyProblem, sol);

    // initialize the coupling manager
    couplingManager->init(
        stokesProblem, darcyProblem,
        std::make_tuple(stokesGridVariables->faceGridVariablesPtr(),
                        stokesGridVariables->cellCenterGridVariablesPtr(),
                        darcyGridVariables),
        sol, CouplingManager::CouplingMode::reconstructFreeFlowNormalStress);

    //couplingManager->setGridVariables(std::make_tuple(stokesGridVariables->cellCenterGridVariablesPtr(),
    //                                                  stokesGridVariables->faceGridVariablesPtr(),
    //                                                  darcyGridVariables));

    // intialize the vtk output module
    StaggeredVtkOutputModule<StokesGridVariables, decltype(stokesSol)>
        stokesVtkWriter(*stokesGridVariables, stokesSol, stokesProblem->name());
    GetPropType<StokesTypeTag, Properties::IOFields>::initOutputModule(
        stokesVtkWriter);
    stokesVtkWriter.write(0.0);

    VtkOutputModule<DarcyGridVariables,
                    GetPropType<DarcyTypeTag, Properties::SolutionVector>>
        darcyVtkWriter(*darcyGridVariables, sol[darcyIdx],
                       darcyProblem->name());
    using DarcyVelocityOutput =
        GetPropType<DarcyTypeTag, Properties::VelocityOutput>;
    darcyVtkWriter.addVelocityOutput(
        std::make_shared<DarcyVelocityOutput>(*darcyGridVariables));
    GetPropType<DarcyTypeTag, Properties::IOFields>::initOutputModule(
        darcyVtkWriter);
    darcyVtkWriter.write(0.0);

    // the assembler for a stationary problem
    using Assembler =
        MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(
        std::make_tuple(stokesProblem, stokesProblem, darcyProblem),
        std::make_tuple(stokesGridGeometry->faceFVGridGeometryPtr(),
                        stokesGridGeometry->cellCenterFVGridGeometryPtr(),
                        darcyGridGeometry),
        std::make_tuple(stokesGridVariables->faceGridVariablesPtr(),
                        stokesGridVariables->cellCenterGridVariablesPtr(),
                        darcyGridVariables),
        couplingManager);

    // the linear solver
    using LinearSolver = UMFPackBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver =
        MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);

    // solve the non-linear system
    nonLinearSolver.solve(sol);

    // write vtk output
    stokesVtkWriter.write(1.0);
    darcyVtkWriter.write(1.0);

    //    for (const auto& element : elements(darcyGridView))
    //    {
    //        auto fvGeometry = localView(*darcyGridGeometry);
    //        fvGeometry.bind(element);

    //        for (const auto& scvf : scvfs(fvGeometry))
    //        {
    //            if (couplingManager->isCoupledEntity(CouplingManager::darcyIdx, scvf))
    //            {
    //                std::ostream tmp(std::cout.rdbuf());

    //                tmp << std::scientific << "Interface pressure ff at " << scvf.center() << " is " << couplingManager->couplingData().freeFlowInterfacePressure(element, scvf) << std::endl;
    //            }
    //        }
    //    }

    //    for (const auto& element : elements(stokesGridView))
    //    {
    //        auto fvGeometry = localView(*stokesGridGeometry);
    //        auto elemVolVars = localView(stokesGridVariables->cellCenterGridVariablesPtr()->curGridVolVars());
    //        auto elemFaceVars = localView(stokesGridVariables->faceGridVariablesPtr()->curGridFaceVars());
    //        // auto elemFluxVarsCache = localView(gridVars_(stokesIdx).gridFluxVarsCache());

    //        fvGeometry.bind(element);
    //        elemVolVars.bind(element, fvGeometry, sol[stokesCellCenterIdx]);
    //        elemFaceVars.bind(element, fvGeometry, sol[stokesFaceIdx]);

    //        for (const auto& scvf : scvfs(fvGeometry))
    //        {
    //            if (couplingManager->isCoupledEntity(CouplingManager::stokesIdx, scvf))
    //            {
    //                std::ostream tmp(std::cout.rdbuf());
    //                tmp << std::scientific << "Interface pressure pm at " << scvf.center() << " is " << stokesProblem->neumann(element,
    //                                    fvGeometry,
    //                                    elemVolVars,
    //                                    elemFaceVars,
    //                                    scvf)[scvf.directionIndex()] << ", vself " << elemFaceVars[scvf].velocitySelf() << ", vDarcy " << couplingManager->couplingData().darcyInterfaceVelocity(element,
    //                                                                  fvGeometry,
    //                                                                  elemVolVars,
    //                                                                  elemFaceVars,
    //                                                                  scvf)  << std::endl;
    //            }
    //        }
    //    }

    /*
  {
    double min = std::numeric_limits<double>::max();
    double max = std::numeric_limits<double>::min();
    double sum = 0.;
    using FluxVariables = GetPropType<DarcyTypeTag, Properties::FluxVariables>;
    const std::string filename = getParam<std::string>("Problem.Name") + "-" + darcyProblem->name() + "-interface-velocity";
    std::tie(min, max, sum) = writeVelocitiesOnInterfaceToFile<FluxVariables>(filename,
                                                                              *couplingManager,
                                                                              *darcyProblem,
                                                                              *darcyGridVariables,
                                                                              sol[darcyIdx]);
    const auto prec = std::cout.precision();
    std::cout << "Velocity statistics (Darcy):" << std::endl
              << std::setprecision(std::numeric_limits<double>::digits10 + 1)
              << "  min: " << min << std::endl
              << "  max: " << max << std::endl
              << "  sum: " << sum << std::endl;
    std::cout.precision(prec);

    {
      const std::string filenameDarcy = "darcy-flow-statistics.txt";
      std::ofstream ofs(filenameDarcy + ".txt", std::ofstream::out | std::ofstream::trunc);
      const auto prec = ofs.precision();
      ofs << "Velocity statistics (Darcy):" << std::endl
          << std::setprecision(std::numeric_limits<double>::digits10 + 1)
          << "  min: " << min << std::endl
          << "  max: " << max << std::endl
          << "  sum: " << sum << std::endl;
      ofs.precision(prec);
      ofs.close();
    }

    {
      double min = std::numeric_limits<double>::max();
      double max = std::numeric_limits<double>::min();
      double sum = 0.;
      const std::string filename = getParam<std::string>("Problem.Name") + "-" + darcyProblem->name() + "-interface-pressure";
      std::tie(min, max, sum) = writePressuresOnInterfaceToFile(filename,
                                                                *darcyProblem,
                                                                *darcyGridVariables,
                                                                sol[darcyIdx]);
      const auto prec = std::cout.precision();
      std::cout << "Pressure statistics (Darcy):" << std::endl
                << std::setprecision(std::numeric_limits<double>::digits10 + 1)
                << "  min: " << min << std::endl
                << "  max: " << max << std::endl
                << "  sum: " << sum << std::endl;
      std::cout.precision(prec);
    }
  }

  //TODO make freeflow
  {
    double min = std::numeric_limits<double>::max();
    double max = std::numeric_limits<double>::min();
    double sum = 0.;
    const std::string filename = getParam<std::string>("Problem.Name") + "-" + stokesProblem->name() + "-interface-velocity";
    std::tie(min, max, sum) = writeStokesVelocitiesOnInterfaceToFile(filename,
                                                                     *couplingManager,
                                                                     *stokesProblem,
                                                                     sol[stokesFaceIdx]);
    const auto prec = std::cout.precision();
    std::cout << "Velocity statistics (free flow):" << std::endl
              << std::setprecision(std::numeric_limits<double>::digits10 + 1)
              << "  min: " << min << std::endl
              << "  max: " << max << std::endl
              << "  sum: " << sum << std::endl;
    std::cout.precision(prec);
    {
      const std::string filenameFlow = "free-flow-statistics";
      std::ofstream ofs(filenameFlow + ".txt", std::ofstream::out | std::ofstream::trunc);
      const auto prec = ofs.precision();
      ofs << "Velocity statistics (free flow):" << std::endl
          << std::setprecision(std::numeric_limits<double>::digits10 + 1)
          << "  min: " << min << std::endl
          << "  max: " << max << std::endl
          << "  sum: " << sum << std::endl;
      ofs.precision(prec);
      ofs.close();
    }
  }
*/
    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////

    // print dumux end message
    if (mpiHelper.rank() == 0) {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }

    return 0;
}  // end main
catch (Dumux::ParameterException &e) {
    std::cerr << std::endl << e << " ---> Abort!" << std::endl;
    return 1;
} catch (Dune::DGFException &e) {
    std::cerr << "DGF exception thrown (" << e
              << "). Most likely, the DGF file name is wrong "
                 "or the DGF file is corrupted, "
                 "e.g. missing hash at end of file or wrong number "
                 "(dimensions) of entries."
              << " ---> Abort!" << std::endl;
    return 2;
} catch (Dune::Exception &e) {
    std::cerr << "Dune reported error: " << e << " ---> Abort!" << std::endl;
    return 3;
} catch (...) {
    std::cerr << "Unknown exception thrown! ---> Abort!" << std::endl;
    return 4;
}
