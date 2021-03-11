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

#include <dumux/multidomain/boundary/freeflowsolidenergy/couplingmanager.hh>

#include "problem_freeflow.hh"
#include "problem_heat.hh"

namespace Dumux
{
namespace Properties
{
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::FreeFlowModel> {
    using Traits = StaggeredMultiDomainTraits<TypeTag,
                                              TypeTag,
                                              Properties::TTag::HeatModel>;
    using type = Dumux::FreeFlowSolidEnergyCouplingManager<Traits>;
};

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::HeatModel> {
    using Traits = StaggeredMultiDomainTraits<Properties::TTag::FreeFlowModel,
                                              Properties::TTag::FreeFlowModel,
                                              TypeTag>;
    using type = Dumux::FreeFlowSolidEnergyCouplingManager<Traits>;
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
    using FreeFlowTypeTag = Properties::TTag::FreeFlowModel;
    using SolidEnergyTypeTag = Properties::TTag::HeatModel;

    // try to create a grid (from the given grid file or the input file)
    // for both sub-domains
    using SolidEnergyGridManager =
        Dumux::GridManager<GetPropType<SolidEnergyTypeTag, Properties::Grid>>;
    SolidEnergyGridManager solidEnergyGridManager;
    solidEnergyGridManager.init("SolidEnergy");  // pass parameter group

    using FreeFlowGridManager =
        Dumux::GridManager<GetPropType<FreeFlowTypeTag, Properties::Grid>>;
    FreeFlowGridManager freeFlowGridManager;
    freeFlowGridManager.init("FreeFlow");  // pass parameter group

    // we compute on the leaf grid view
    const auto &solidEnergyGridView =
        solidEnergyGridManager.grid().leafGridView();
    const auto &freeFlowGridView = freeFlowGridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using FreeFlowFVGridGeometry =
        GetPropType<FreeFlowTypeTag, Properties::FVGridGeometry>;
    auto freeFlowFvGridGeometry =
        std::make_shared<FreeFlowFVGridGeometry>(freeFlowGridView);
    freeFlowFvGridGeometry->update();
    using SolidEnergyFVGridGeometry =
        GetPropType<SolidEnergyTypeTag, Properties::FVGridGeometry>;
    auto solidEnergyFvGridGeometry =
        std::make_shared<SolidEnergyFVGridGeometry>(solidEnergyGridView);
    solidEnergyFvGridGeometry->update();

    using Traits = StaggeredMultiDomainTraits<FreeFlowTypeTag, FreeFlowTypeTag,
                                              SolidEnergyTypeTag>;

    // the coupling manager
    using CouplingManager = FreeFlowSolidEnergyCouplingManager<Traits>;
    auto couplingManager = std::make_shared<CouplingManager>(
        freeFlowFvGridGeometry, solidEnergyFvGridGeometry);

    // the indices
    constexpr auto freeFlowCellCenterIdx =
        CouplingManager::freeFlowCellCenterIdx;
    constexpr auto freeFlowFaceIdx = CouplingManager::freeFlowFaceIdx;
    constexpr auto solidEnergyIdx = CouplingManager::solidEnergyIdx;

    // the problem (initial and boundary conditions)
    using FreeFlowProblem = GetPropType<FreeFlowTypeTag, Properties::Problem>;
    auto freeFlowProblem = std::make_shared<FreeFlowProblem>(
        freeFlowFvGridGeometry, couplingManager);
    using SolidEnergyProblem =
        GetPropType<SolidEnergyTypeTag, Properties::Problem>;
    auto solidEnergyProblem = std::make_shared<SolidEnergyProblem>(
        solidEnergyFvGridGeometry, couplingManager);

    // the solution vector
    Traits::SolutionVector sol;
    sol[freeFlowCellCenterIdx].resize(
        freeFlowFvGridGeometry->numCellCenterDofs());
    sol[freeFlowFaceIdx].resize(freeFlowFvGridGeometry->numFaceDofs());
    sol[solidEnergyIdx].resize(solidEnergyFvGridGeometry->numDofs());

    // get a solution vector storing references to the two Stokes solution vectors
    auto stokesSol = partial(sol, freeFlowCellCenterIdx, freeFlowFaceIdx);

    // apply initial solution for instationary problems
    freeFlowProblem->applyInitialSolution(stokesSol);
    solidEnergyProblem->applyInitialSolution(sol[solidEnergyIdx]);

    auto solOld = sol;

    couplingManager->init(freeFlowProblem, solidEnergyProblem, sol);

    // the grid variables
    using FreeFlowGridVariables =
        GetPropType<FreeFlowTypeTag, Properties::GridVariables>;
    auto freeFlowGridVariables = std::make_shared<FreeFlowGridVariables>(
        freeFlowProblem, freeFlowFvGridGeometry);
    freeFlowGridVariables->init(stokesSol);
    using SolidEnergyGridVariables =
        GetPropType<SolidEnergyTypeTag, Properties::GridVariables>;
    auto solidEnergyGridVariables = std::make_shared<SolidEnergyGridVariables>(
        solidEnergyProblem, solidEnergyFvGridGeometry);
    solidEnergyGridVariables->init(sol[solidEnergyIdx]);

    // intialize the vtk output module
    StaggeredVtkOutputModule<FreeFlowGridVariables, decltype(stokesSol)>
        freeFlowVtkWriter(*freeFlowGridVariables, stokesSol,
                          freeFlowProblem->name());
    GetPropType<FreeFlowTypeTag, Properties::IOFields>::initOutputModule(
        freeFlowVtkWriter);
    freeFlowVtkWriter.write(0.0);

    VtkOutputModule<SolidEnergyGridVariables,
                    GetPropType<SolidEnergyTypeTag, Properties::SolutionVector>>
        solidEnergyVtkWriter(*solidEnergyGridVariables, sol[solidEnergyIdx],
                             solidEnergyProblem->name());
    GetPropType<SolidEnergyTypeTag, Properties::IOFields>::initOutputModule(
        solidEnergyVtkWriter);
    solidEnergyVtkWriter.write(0.0);

    // instantiate time loop
    using Scalar = GetPropType<FreeFlowTypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    // the assembler for a stationary problem
    using Assembler =
        MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(
        std::make_tuple(freeFlowProblem, freeFlowProblem, solidEnergyProblem),
        std::make_tuple(freeFlowFvGridGeometry->cellCenterFVGridGeometryPtr(),
                        freeFlowFvGridGeometry->faceFVGridGeometryPtr(),
                        solidEnergyFvGridGeometry),
        std::make_tuple(freeFlowGridVariables->cellCenterGridVariablesPtr(),
                        freeFlowGridVariables->faceGridVariablesPtr(),
                        solidEnergyGridVariables),
        couplingManager, timeLoop);

    // the linear solver
    using LinearSolver = UMFPackBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver =
        MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);

    // time loop
    timeLoop->start();
    do {
        // set previous solution for storage evaluations
        assembler->setPreviousSolution(solOld);

        // solve the non-linear system with time step control
        nonLinearSolver.solve(sol, *timeLoop);

        // make the new solution the old solution
        solOld = sol;
        freeFlowGridVariables->advanceTimeStep();
        solidEnergyGridVariables->advanceTimeStep();

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();

        // write vtk output
        freeFlowVtkWriter.write(timeLoop->time());
        solidEnergyVtkWriter.write(timeLoop->time());

        // report statistics of this time step
        timeLoop->reportTimeStep();

        // set new dt as suggested by newton solver
        timeLoop->setTimeStepSize(
            nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));

    } while (!timeLoop->finished());

    timeLoop->finalize(freeFlowGridView.comm());
    timeLoop->finalize(solidEnergyGridView.comm());

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
