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

 #include <ctime>
 #include <iostream>

 #include <dune/common/parallel/mpihelper.hh>
 #include <dune/common/timer.hh>
 #include <dune/istl/io.hh>

 #include <dumux/common/properties.hh>
 #include <dumux/common/parameters.hh>
 #include <dumux/common/partial.hh>
 #include <dumux/common/dumuxmessage.hh>
 #include <dumux/linear/seqsolverbackend.hh>
 #include <dumux/assembly/fvassembler.hh>
 #include <dumux/assembly/diffmethod.hh>
 #include <dumux/discretization/method.hh>
 #include <dumux/io/vtkoutputmodule.hh>
 #include <dumux/io/staggeredvtkoutputmodule.hh>
 #include <dumux/io/grid/gridmanager.hh>

 #include <dumux/assembly/staggeredfvassembler.hh>
 #include <dumux/nonlinear/newtonsolver.hh>

#include "../monolithic/ffproblem.hh"

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
    using FreeFlowTypeTag = Properties::TTag::FreeFlowModel;

    // try to create a grid (from the given grid file or the input file)
    using FreeFlowGridManager = Dumux::GridManager<GetPropType<FreeFlowTypeTag, Properties::Grid>>;
    FreeFlowGridManager freeFlowGridManager;
    freeFlowGridManager.init("FreeFlow"); // pass parameter group

    // we compute on the leaf grid view
    const auto& freeFlowGridView = freeFlowGridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using FreeFlowFVGridGeometry = GetPropType<FreeFlowTypeTag, Properties::FVGridGeometry>;
    auto freeFlowFvGridGeometry = std::make_shared<FreeFlowFVGridGeometry>(freeFlowGridView);
    freeFlowFvGridGeometry->update();

    // the problem (initial and boundary conditions)
    using FreeFlowProblem = GetPropType<FreeFlowTypeTag, Properties::Problem>;
    auto freeFlowProblem = std::make_shared<FreeFlowProblem>(freeFlowFvGridGeometry);

    // the solution vector
    GetPropType<FreeFlowTypeTag, Properties::SolutionVector> sol;
    sol[FreeFlowFVGridGeometry::cellCenterIdx()].resize(freeFlowFvGridGeometry->numCellCenterDofs());
    sol[FreeFlowFVGridGeometry::faceIdx()].resize(freeFlowFvGridGeometry->numFaceDofs());


    // apply initial solution for instationary problems
    freeFlowProblem->applyInitialSolution(sol);

    // the grid variables
    using FreeFlowGridVariables = GetPropType<FreeFlowTypeTag, Properties::GridVariables>;
    auto freeFlowGridVariables = std::make_shared<FreeFlowGridVariables>(freeFlowProblem, freeFlowFvGridGeometry);
    freeFlowGridVariables->init(sol);

    // intialize the vtk output module
    StaggeredVtkOutputModule<FreeFlowGridVariables, decltype(sol)> freeFlowVtkWriter(*freeFlowGridVariables, sol, freeFlowProblem->name());
    GetPropType<FreeFlowTypeTag, Properties::IOFields>::initOutputModule(freeFlowVtkWriter);
    freeFlowVtkWriter.addField(freeFlowProblem->getAnalyticalVelocityX(), "analyticalV_x");
    freeFlowVtkWriter.write(0.0);

    // the assembler for a stationary problem
    using Assembler = StaggeredFVAssembler<FreeFlowTypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(freeFlowProblem, freeFlowFvGridGeometry, freeFlowGridVariables);

    // the linear solver
    using LinearSolver = UMFPackBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    // solve the non-linear system
    nonLinearSolver.solve(sol);

    // write vtk output
    freeFlowVtkWriter.write(1.0);

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
