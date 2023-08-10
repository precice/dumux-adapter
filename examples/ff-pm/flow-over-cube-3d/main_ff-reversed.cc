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

#include <dumux/assembly/staggeredfvassembler.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include "ffproblem-reversed.hh"

//TODO
// Helper function to put pressure on interface

template<class ElementFaceVariables, class SubControlVolumeFace>
auto velocityAtInterface(const ElementFaceVariables &elemFaceVars,
                         const SubControlVolumeFace &scvf)
{
    const double scalarVelocity = elemFaceVars[scvf].velocitySelf();
    auto velocity = scvf.unitOuterNormal();
    velocity[scvf.directionIndex()] = scalarVelocity;
    return velocity;
}

template<class FluxVariables,
         class Problem,
         class Element,
         class SubControlVolumeFace,
         class FVElementGeometry,
         class ElementVolumeVariables,
         class ElementFaceVariables,
         class ElementFluxVariablesCache>
auto pressureAtInterface(const Problem &problem,
                         const Element &element,
                         const SubControlVolumeFace &scvf,
                         const FVElementGeometry &fvGeometry,
                         const ElementVolumeVariables &elemVolVars,
                         const ElementFaceVariables &elemFaceVars,
                         const ElementFluxVariablesCache &elemFluxVarsCache)
{
    FluxVariables fluxVars;
    return fluxVars.computeMomentumFlux(problem, element, scvf, fvGeometry,
                                        elemVolVars, elemFaceVars,
                                        elemFluxVarsCache.gridFluxVarsCache()) /
           scvf.area();
}

template<class FluxVariables,
         class Problem,
         class GridVariables,
         class SolutionVector>
void setInterfacePressures(const Problem &problem,
                           const GridVariables &gridVars,
                           const SolutionVector &sol,
                           const precice::string_view meshNameView,
                           const precice::string_view dataNameView)
{
    const auto &gridGeometry = problem.gridGeometry();
    auto fvGeometry = localView(gridGeometry);
    auto elemVolVars = localView(gridVars.curGridVolVars());
    auto elemFaceVars = localView(gridVars.curGridFaceVars());
    auto elemFluxVarsCache = localView(gridVars.gridFluxVarsCache());

    auto &couplingParticipant = Dumux::Precice::CouplingAdapter::getInstance();

    for (const auto &element : elements(gridGeometry.gridView())) {
        fvGeometry.bind(element);
        elemVolVars.bind(element, fvGeometry, sol);
        elemFaceVars.bind(element, fvGeometry, sol);
        elemFluxVarsCache.bind(element, fvGeometry, elemVolVars);

        for (const auto &scvf : scvfs(fvGeometry)) {
            if (couplingParticipant.isCoupledEntity(scvf.index())) {
                //TODO: What to do here?
                const auto p = pressureAtInterface<FluxVariables>(
                    problem, element, scvf, fvGeometry, elemVolVars,
                    elemFaceVars, elemFluxVarsCache);
                couplingParticipant.writeScalarQuantityOnFace(meshNameView,
                                                            dataNameView,
                                                            scvf.index(), p);
            }
        }
    }
}

template<class Problem, class GridVariables, class SolutionVector>
void setInterfaceVelocities(const Problem &problem,
                            const GridVariables &gridVars,
                            const SolutionVector &sol,
                            const precice::string_view meshNameView,
                            const precice::string_view dataNameView)
{
    const auto &gridGeometry = problem.gridGeometry();
    auto fvGeometry = localView(gridGeometry);
    auto elemVolVars = localView(gridVars.curGridVolVars());
    auto elemFaceVars = localView(gridVars.curGridFaceVars());

    auto &couplingParticipant = Dumux::Precice::CouplingAdapter::getInstance();

    for (const auto &element : elements(gridGeometry.gridView())) {
        fvGeometry.bindElement(element);
        elemVolVars.bindElement(element, fvGeometry, sol);
        elemFaceVars.bindElement(element, fvGeometry, sol);

        for (const auto &scvf : scvfs(fvGeometry)) {
            if (couplingParticipant.isCoupledEntity(scvf.index())) {
                //TODO: What to do here?
                const auto v = velocityAtInterface(elemFaceVars,
                                                   scvf)[scvf.directionIndex()];
                couplingParticipant.writeScalarQuantityOnFace(meshNameView, dataNameView,
                                                            scvf.index(), v);
            }
        }
    }
}

template<class Problem, class GridVariables, class SolutionVector>
std::tuple<double, double, double> writeVelocitiesOnInterfaceToFile(
    const precice::string_view &meshNameView,
    const std::string &filename,
    const Problem &problem,
    const GridVariables &gridVars,
    const SolutionVector &sol)
{
    const auto &gridGeometry = problem.gridGeometry();
    auto fvGeometry = localView(gridGeometry);
    auto elemVolVars = localView(gridVars.curGridVolVars());
    auto elemFaceVars = localView(gridVars.curGridFaceVars());

    const auto &couplingParticipant =
        Dumux::Precice::CouplingAdapter::getInstance();

    std::ofstream ofs(filename + ".csv",
                      std::ofstream::out | std::ofstream::trunc);
    ofs << "x,y,";
    if (couplingParticipant.getMeshDimensions(meshNameView) == 3)
        ofs << "z,";
    ofs << "velocityY"
        << "\n";

    double min = std::numeric_limits<double>::max();
    double max = std::numeric_limits<double>::min();
    double sum = 0.;
    for (const auto &element : elements(gridGeometry.gridView())) {
        fvGeometry.bind(element);
        elemVolVars.bind(element, fvGeometry, sol);
        elemFaceVars.bindElement(element, fvGeometry, sol);

        for (const auto &scvf : scvfs(fvGeometry)) {
            if (couplingParticipant.isCoupledEntity(scvf.index())) {
                const auto &pos = scvf.center();
                for (int i = 0; i < couplingParticipant.getMeshDimensions(meshNameView); ++i) {
                    ofs << pos[i] << ",";
                }
                const double v = problem.dirichlet(element, scvf)[1];
                max = std::max(v, max);
                min = std::min(v, min);
                sum += v;
                //velocityAtInterface(elemFaceVars, scvf)[scvf.directionIndex()];
                const int prec = ofs.precision();
                ofs << std::setprecision(std::numeric_limits<double>::digits10 +
                                         1)
                    << v << "\n";
                ofs.precision(prec);
            }
        }
    }

    ofs.close();

    return std::make_tuple(min, max, sum);
}

template<class FluxVariables,
         class Problem,
         class GridVariables,
         class SolutionVector>
void writePressuresOnInterfaceToFile(const precice::string_view &meshNameView,
                                     const std::string &filename,
                                     const Problem &problem,
                                     const GridVariables &gridVars,
                                     const SolutionVector &sol)
{
    const auto &gridGeometry = problem.gridGeometry();
    auto fvGeometry = localView(gridGeometry);
    auto elemVolVars = localView(gridVars.curGridVolVars());
    auto elemFaceVars = localView(gridVars.curGridFaceVars());
    auto elemFluxVarsCache = localView(gridVars.gridFluxVarsCache());

    const auto &couplingParticipant =
        Dumux::Precice::CouplingAdapter::getInstance();

    std::ofstream ofs(filename + ".csv",
                      std::ofstream::out | std::ofstream::trunc);
    ofs << "x,y,";
    if (couplingParticipant.getMeshDimensions(meshNameView) == 3)
        ofs << "z,";
    ofs << "pressure"
        << "\n";
    for (const auto &element : elements(gridGeometry.gridView())) {
        fvGeometry.bind(element);
        elemVolVars.bind(element, fvGeometry, sol);
        elemFaceVars.bind(element, fvGeometry, sol);
        elemFluxVarsCache.bind(element, fvGeometry, elemVolVars);

        for (const auto &scvf : scvfs(fvGeometry)) {
            if (couplingParticipant.isCoupledEntity(scvf.index())) {
                const auto &pos = scvf.center();
                for (int i = 0; i < couplingParticipant.getMeshDimensions(meshNameView); ++i) {
                    ofs << pos[i] << ",";
                }
                const double p = pressureAtInterface<FluxVariables>(
                    problem, element, scvf, fvGeometry, elemVolVars,
                    elemFaceVars, elemFluxVarsCache);
                ofs << p << "\n";
            }
        }
    }

    ofs.close();
} 

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

    // try to create a grid (from the given grid file or the input file)
    using FreeFlowGridManager =
        Dumux::GridManager<GetPropType<FreeFlowTypeTag, Properties::Grid>>;
    FreeFlowGridManager freeFlowGridManager;
    freeFlowGridManager.init("FreeFlow");  // pass parameter group

    // we compute on the leaf grid view
    const auto &freeFlowGridView = freeFlowGridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using FreeFlowGridGeometry =
        GetPropType<FreeFlowTypeTag, Properties::GridGeometry>;
    auto freeFlowGridGeometry =
        std::make_shared<FreeFlowGridGeometry>(freeFlowGridView);
#if DUMUX_VERSION_MAJOR >= 3 & DUMUX_VERSION_MINOR >= 5
    freeFlowGridGeometry->update(freeFlowGridManager.grid().leafGridView());
#else
    freeFlowGridGeometry->update();
#endif

    // the problem (initial and boundary conditions)
    using FreeFlowProblem = GetPropType<FreeFlowTypeTag, Properties::Problem>;
    auto freeFlowProblem =
        std::make_shared<FreeFlowProblem>(freeFlowGridGeometry);

    // the solution vector
    GetPropType<FreeFlowTypeTag, Properties::SolutionVector> sol;
    sol[FreeFlowGridGeometry::cellCenterIdx()].resize(
        freeFlowGridGeometry->numCellCenterDofs());
    sol[FreeFlowGridGeometry::faceIdx()].resize(
        freeFlowGridGeometry->numFaceDofs());

    // Initialize preCICE.Tell preCICE about:
    // - Name of solver
    // - What rank of how many ranks this instance is
    // Configure preCICE. For now the config file is hardcoded.
    //couplingParticipant.createInstance( "FreeFlow", mpiHelper.rank(), mpiHelper.size() );
    std::string preciceConfigFilename = "precice-config.xml";
    //    if (argc == 3)
    //      preciceConfigFilename = argv[2];
    if (argc > 2)
        preciceConfigFilename = argv[argc - 1];

    auto &couplingParticipant = Dumux::Precice::CouplingAdapter::getInstance();
    couplingParticipant.announceSolver("FreeFlow", preciceConfigFilename,
                                     mpiHelper.rank(), mpiHelper.size());

    const precice::string_view meshNameView = std::string("FreeFlowMesh");
    const int dim = couplingParticipant.getMeshDimensions(meshNameView);
    std::cout << dim << "  " << int(FreeFlowGridGeometry::GridView::dimension)
              << std::endl;
    if (dim != int(FreeFlowGridGeometry::GridView::dimension))
        DUNE_THROW(Dune::InvalidStateException, "Dimensions do not match"); 

    // GET mesh corodinates
    const double xMin =
        getParamFromGroup<std::vector<double>>("Darcy", "Grid.LowerLeft")[0];
    const double xMax =
        getParamFromGroup<std::vector<double>>("Darcy", "Grid.UpperRight")[0];
    std::vector<double> coords;  //( dim * vertexSize );
    std::vector<int> coupledScvfIndices;
    precice::span<double> coordsSpan(coords);

    for (const auto &element : elements(freeFlowGridView)) {
        auto fvGeometry = localView(*freeFlowGridGeometry);
        fvGeometry.bindElement(element);

        for (const auto &scvf : scvfs(fvGeometry)) {
            static constexpr auto eps = 1e-7;
            const auto &pos = scvf.center();
            if (pos[1] < freeFlowGridGeometry->bBoxMin()[1] + eps) {
                if (pos[0] > xMin - eps && pos[0] < xMax + eps) {
                    coupledScvfIndices.push_back(scvf.index());
                    for (const auto p : pos)
                        coords.push_back(p);
                }
            }
        }
    }

    const auto numberOfPoints = coords.size() / dim;
    double preciceDt = couplingParticipant.getMaxTimeStepSize();
    couplingParticipant.setMesh(meshNameView, coordsSpan);
    couplingParticipant.createIndexMapping(coupledScvfIndices);

    const precice::string_view dataNameViewV = std::string("Velocity");
    const precice::string_view dataNameViewP = std::string("Pressure");
    couplingParticipant.announceQuantity(meshNameView, dataNameViewV);
    couplingParticipant.announceQuantity(meshNameView, dataNameViewP);

    // apply initial solution for instationary problems
    freeFlowProblem->applyInitialSolution(sol);

    // the grid variables
    using FreeFlowGridVariables =
        GetPropType<FreeFlowTypeTag, Properties::GridVariables>;
    auto freeFlowGridVariables = std::make_shared<FreeFlowGridVariables>(
        freeFlowProblem, freeFlowGridGeometry);
    freeFlowGridVariables->init(sol);

    // intialize the vtk output module
    StaggeredVtkOutputModule<FreeFlowGridVariables, decltype(sol)>
        freeFlowVtkWriter(*freeFlowGridVariables, sol, freeFlowProblem->name());
    GetPropType<FreeFlowTypeTag, Properties::IOFields>::initOutputModule(
        freeFlowVtkWriter);
    freeFlowVtkWriter.addField(freeFlowProblem->getAnalyticalVelocityX(),
                               "analyticalV_x");
    freeFlowVtkWriter.write(0.0);

    using FluxVariables =
        GetPropType<FreeFlowTypeTag, Properties::FluxVariables>;

    if (couplingParticipant.hasToWriteInitialData()) {
        //TODO
        //      couplingParticipant.writeQuantityVector( pressureId );

        setInterfacePressures<FluxVariables>(*freeFlowProblem,
                                             *freeFlowGridVariables, sol, meshNameView, dataNameViewP);
        couplingParticipant.writeQuantityToOtherSolver(meshNameView, dataNameViewP);
    }
    couplingParticipant.initialize();

    // the assembler for a stationary problem
    using Assembler =
        StaggeredFVAssembler<FreeFlowTypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(
        freeFlowProblem, freeFlowGridGeometry, freeFlowGridVariables);

    // the linear solver
    using LinearSolver = UMFPackBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    auto dt = preciceDt;
    auto sol_checkpoint = sol;

    double vtkTime = 1.0;
    size_t iter = 0;

    while (couplingParticipant.isCouplingOngoing()) {
        if (couplingParticipant.hasToWriteIterationCheckpoint()) {
            //DO CHECKPOINTING
            sol_checkpoint = sol;
        }

        // TODO
        couplingParticipant.readQuantityFromOtherSolver(meshNameView, dataNameViewV, dt);
        //        // For testing
        //        {
        //          const auto v = couplingParticipant.getQuantityVector( velocityId );
        //          const double sum = std::accumulate( v.begin(), v.end(), 0. );
        //          std::cout << "Sum of velocities over boundary to pm: \n" << sum << std::endl;
        //        }

        // solve the non-linear system
        nonLinearSolver.solve(sol);

        // TODO
        setInterfacePressures<FluxVariables>(*freeFlowProblem,
                                             *freeFlowGridVariables, sol, meshNameView, dataNameViewP);
        couplingParticipant.writeQuantityToOtherSolver(meshNameView, dataNameViewP);

        //Read checkpoint
        freeFlowVtkWriter.write(vtkTime);
        vtkTime += 1.;
        preciceDt = couplingParticipant.getMaxTimeStepSize();
        dt = std::min(preciceDt, dt);

        ++iter;

        if (couplingParticipant.hasToReadIterationCheckpoint()) {
            //            //Read checkpoint
            //            freeFlowVtkWriter.write(vtkTime);
            //            vtkTime += 1.;
            sol = sol_checkpoint;
            freeFlowGridVariables->update(sol);
            freeFlowGridVariables->advanceTimeStep();
            //freeFlowGridVariables->init(sol);
        } else  // coupling successful
        {
            // write vtk output
            freeFlowVtkWriter.write(vtkTime);
        }
    }
    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////

    couplingParticipant.finalize();

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
