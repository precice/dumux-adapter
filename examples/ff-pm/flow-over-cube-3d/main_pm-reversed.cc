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
#include <fstream>
#include <iostream>
#include <string>

bool printstuff = false;

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/istl/io.hh>

#include <dumux/common/defaultusagemessage.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/linear/amgbackend.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/assembly/diffmethod.hh>
#include <dumux/assembly/fvassembler.hh>

#include <dumux/discretization/method.hh>

#include <dumux/io/grid/gridmanager.hh>
#include <dumux/io/vtkoutputmodule.hh>

#include "pmproblem-reversed.hh"

/*!
  * \brief Returns the pressure at the interface using Darcy's law for reconstruction
  */
template<class Problem,
         class Element,
         class FVElementGeometry,
         class ElementVolumeVariables,
         class SubControlVolumeFace,
         class ElementFluxVariablesCache>
auto pressureAtInterface(const Problem &problem,
                         const Element &element,
                         const FVElementGeometry &fvGeometry,
                         const ElementVolumeVariables &elemVolVars,
                         const SubControlVolumeFace &scvf,
                         const ElementFluxVariablesCache &elemFluxVarsCache)
{
    using Scalar = typename ElementVolumeVariables::VolumeVariables::
        PrimaryVariables::value_type;
    const auto &volVars = elemVolVars[scvf.insideScvIdx()];

    const Scalar boundaryFlux = problem.neumann(
        element, fvGeometry, elemVolVars, elemFluxVarsCache, scvf);

    const auto K = volVars.permeability();
    const Scalar ccPressure = volVars.pressure();
    const Scalar mobility = volVars.mobility();
    const Scalar density = volVars.density();

    // v = -K/mu * (gradP + rho*g)
    auto velocity = scvf.unitOuterNormal();
    velocity *= boundaryFlux;  // TODO check sign
    velocity /= density;

    // v = -kr/mu*K * (gradP + rho*g) = -mobility*K * (gradP + rho*g)
    const auto alpha =
        Dumux::vtmv(scvf.unitOuterNormal(), K, problem.gravity());

    auto distanceVector = scvf.center() - element.geometry().center();
    distanceVector /= distanceVector.two_norm2();
    const Scalar ti = Dumux::vtmv(distanceVector, K, scvf.unitOuterNormal());

    return (1 / mobility * (scvf.unitOuterNormal() * velocity) +
            density * alpha) /
               ti +
           ccPressure;
}

template<class Problem, class GridVariables, class SolutionVector>
void setInterfacePressures(const Problem &problem,
                           const GridVariables &gridVars,
                           const SolutionVector &sol,
                           const precice::string_view meshNameView,
                           const precice::string_view dataNameView)
{
    const auto &gridGeometry = problem.gridGeometry();
    auto fvGeometry = localView(gridGeometry);
    auto elemVolVars = localView(gridVars.curGridVolVars());
    auto elemFluxVarsCache = localView(gridVars.gridFluxVarsCache());

    auto &couplingParticipant = Dumux::Precice::CouplingAdapter::getInstance();

    for (const auto &element : elements(gridGeometry.gridView())) {
        fvGeometry.bindElement(element);
        elemVolVars.bindElement(element, fvGeometry, sol);

        //sstd::cout << "Pressure by reconstruction" << std::endl;
        for (const auto &scvf : scvfs(fvGeometry)) {
            if (couplingParticipant.isCoupledEntity(scvf.index())) {
                //TODO: What to do here?
                const double p =
                    pressureAtInterface(problem, element, gridGeometry,
                                        elemVolVars, scvf, elemFluxVarsCache);
                couplingParticipant.writeScalarQuantityOnFace(meshNameView, dataNameView,
                                                            scvf.index(), p);
            }
        }
    }
}

/*!
  * \brief Returns the velocity at the interface using Darcy's law for reconstruction
  */
template<class FluxVariables,
         class Problem,
         class Element,
         class FVElementGeometry,
         class ElementVolumeVariables,
         class SubControlVolumeFace,
         class ElementFluxVariablesCache>
auto velocityAtInterface(const Problem &problem,
                         const Element &element,
                         const FVElementGeometry &fvGeometry,
                         const ElementVolumeVariables &elemVolVars,
                         const SubControlVolumeFace &scvf,
                         const ElementFluxVariablesCache &elemFluxVarsCache)
{
    const int phaseIdx = 0;
    FluxVariables fluxVars;
    fluxVars.init(problem, element, fvGeometry, elemVolVars, scvf,
                  elemFluxVarsCache);
    auto upwindTerm = [phaseIdx](const auto &volVars) {
        return volVars.mobility(phaseIdx);
    };
    const auto scalarVelocity =
        fluxVars.advectiveFlux(phaseIdx, upwindTerm) / scvf.area();
    return scalarVelocity;
}

template<class FluxVariables,
         class Problem,
         class GridVariables,
         class SolutionVector>
void setInterfaceVelocities(const Problem &problem,
                            const GridVariables &gridVars,
                            const SolutionVector &sol,
                            const precice::string_view meshNameView,
                            const precice::string_view dataNameView)
{
    const auto &gridGeometry = problem.gridGeometry();
    auto fvGeometry = localView(gridGeometry);
    auto elemVolVars = localView(gridVars.curGridVolVars());
    auto elemFluxVarsCache = localView(gridVars.gridFluxVarsCache());

    auto &couplingParticipant = Dumux::Precice::CouplingAdapter::getInstance();

    for (const auto &element : elements(gridGeometry.gridView())) {
        fvGeometry.bind(element);
        elemVolVars.bind(element, fvGeometry, sol);
        elemFluxVarsCache.bind(element, fvGeometry, elemVolVars);

        for (const auto &scvf : scvfs(fvGeometry)) {
            if (couplingParticipant.isCoupledEntity(scvf.index())) {
                //TODO: What to do here?
                const double v = velocityAtInterface<FluxVariables>(
                    problem, element, fvGeometry, elemVolVars, scvf,
                    elemFluxVarsCache);
                couplingParticipant.writeScalarQuantityOnFace(meshNameView, dataNameView,
                                                            scvf.index(), v);
            }
        }
    }
}

template<class FluxVariables,
         class Problem,
         class GridVariables,
         class SolutionVector>
std::tuple<double, double, double> writeVelocitiesOnInterfaceToFile(
    const precice::string_view &meshName,
    const std::string &filename,
    const Problem &problem,
    const GridVariables &gridVars,
    const SolutionVector &sol)
{
    const auto &gridGeometry = problem.gridGeometry();
    auto fvGeometry = localView(gridGeometry);
    auto elemVolVars = localView(gridVars.curGridVolVars());
    auto elemFluxVarsCache = localView(gridVars.gridFluxVarsCache());

    const auto &couplingParticipant =
        Dumux::Precice::CouplingAdapter::getInstance();

    std::ofstream ofs(filename + ".csv",
                      std::ofstream::out | std::ofstream::trunc);
    ofs << "x,y,";
    if (couplingParticipant.getMeshDimensions(meshName) == 3)
        ofs << "z,";
    ofs << "velocityY"
        << "\n";

    double min = std::numeric_limits<double>::max();
    double max = std::numeric_limits<double>::min();
    double sum = 0.;
    for (const auto &element : elements(gridGeometry.gridView())) {
        fvGeometry.bind(element);
        elemVolVars.bind(element, fvGeometry, sol);
        elemFluxVarsCache.bind(element, fvGeometry, elemVolVars);

        for (const auto &scvf : scvfs(fvGeometry)) {
            if (couplingParticipant.isCoupledEntity(scvf.index())) {
                const auto &pos = scvf.center();
                for (int i = 0; i < couplingParticipant.getMeshDimensions(meshName); ++i) {
                    ofs << pos[i] << ",";
                }
                const double v = velocityAtInterface<FluxVariables>(
                    problem, element, fvGeometry, elemVolVars, scvf,
                    elemFluxVarsCache);
                max = std::max(v, max);
                min = std::min(v, min);
                sum += v;
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

template<class Problem, class GridVariables, class SolutionVector>
void writePressuresOnInterfaceToFile(const precice::string_view &meshName,
                                     std::string &filename,
                                     const Problem &problem,
                                     const GridVariables &gridVars,
                                     const SolutionVector &sol)
{
    const auto &gridGeometry = problem.gridGeometry();
    auto fvGeometry = localView(gridGeometry);
    auto elemVolVars = localView(gridVars.curGridVolVars());
    auto elemFluxVarsCache = localView(gridVars.gridFluxVarsCache());

    const auto &couplingParticipant =
        Dumux::Precice::CouplingAdapter::getInstance();

    std::ofstream ofs(filename + ".csv",
                      std::ofstream::out | std::ofstream::trunc);
    ofs << "x,y,";
    if (couplingParticipant.getMeshDimensions(meshName) == 3)
        ofs << "z,";
    ofs << "pressure"
        << "\n";
    for (const auto &element : elements(gridGeometry.gridView())) {
        fvGeometry.bind(element);
        elemVolVars.bind(element, fvGeometry, sol);
        elemFluxVarsCache.bind(element, fvGeometry, elemVolVars);

        for (const auto &scvf : scvfs(fvGeometry)) {
            if (couplingParticipant.isCoupledEntity(scvf.index())) {
                const auto &pos = scvf.center();
                for (int i = 0; i < couplingParticipant.getMeshDimensions(meshName); ++i) {
                    ofs << pos[i] << ",";
                }
                const double p =
                    pressureAtInterface(problem, element, gridGeometry,
                                        elemVolVars, scvf, elemFluxVarsCache);
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

    using DarcyTypeTag = Properties::TTag::DarcyOneP;

    using DarcyGridManager =
        Dumux::GridManager<GetPropType<DarcyTypeTag, Properties::Grid>>;
    DarcyGridManager darcyGridManager;
    darcyGridManager.init("Darcy");  // pass parameter group

    // we compute on the leaf grid view
    const auto &darcyGridView = darcyGridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using DarcyGridGeometry =
        GetPropType<DarcyTypeTag, Properties::GridGeometry>;
    auto darcyGridGeometry = std::make_shared<DarcyGridGeometry>(darcyGridView);
#if DUMUX_VERSION_MAJOR >= 3 & DUMUX_VERSION_MINOR >= 5
    darcyGridGeometry->update(darcyGridManager.grid().leafGridView());
#else
    darcyGridGeometry->update();
#endif

    using DarcyProblem = GetPropType<DarcyTypeTag, Properties::Problem>;
    auto darcyProblem = std::make_shared<DarcyProblem>(darcyGridGeometry);

    // the solution vector
    GetPropType<DarcyTypeTag, Properties::SolutionVector> sol;
    sol.resize(darcyGridGeometry->numDofs());

    // Initialize preCICE.Tell preCICE about:
    // - Name of solver
    // - What rank of how many ranks this instance is
    // Configure preCICE. For now the config file is hardcoded.
    //couplingParticipant.createInstance( "darcy", mpiHelper.rank(), mpiHelper.size() );
    std::string preciceConfigFilename = "precice-config.xml";
    //    if (argc == 3)
    //      preciceConfigFilename = argv[2];
    if (argc > 2)
        preciceConfigFilename = argv[argc - 1];

    auto &couplingParticipant = Dumux::Precice::CouplingAdapter::getInstance();
    couplingParticipant.announceSolver("Darcy", preciceConfigFilename,
                                     mpiHelper.rank(), mpiHelper.size());

    const precice::string_view meshNameView = std::string("DarcyMesh");
    const int dim = couplingParticipant.getMeshDimensions(meshNameView);
    std::cout << dim << "  " << int(DarcyGridGeometry::GridView::dimension)
              << std::endl;
    if (dim != int(DarcyGridGeometry::GridView::dimension))
        DUNE_THROW(Dune::InvalidStateException, "Dimensions do not match");

    // GET mesh corodinates
    const double xMin =
        getParamFromGroup<std::vector<double>>("Darcy", "Grid.LowerLeft")[0];
    const double xMax =
        getParamFromGroup<std::vector<double>>("Darcy", "Grid.UpperRight")[0];
    std::vector<double> coords;  //( dim * vertexSize );
    std::vector<int> coupledScvfIndices;
    precice::span<double> coordsSpan(coords);

    for (const auto &element : elements(darcyGridView)) {
        auto fvGeometry = localView(*darcyGridGeometry);
        fvGeometry.bindElement(element);

        for (const auto &scvf : scvfs(fvGeometry)) {
            static constexpr auto eps = 1e-7;
            const auto &pos = scvf.center();
            if (pos[1] > darcyGridGeometry->bBoxMax()[1] - eps) {
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
    couplingParticipant.announceQuantity(meshNameView, dataNameViewP);
    couplingParticipant.announceQuantity(meshNameView, dataNameViewV);

    darcyProblem->applyInitialSolution(sol);

    // the grid variables
    using DarcyGridVariables =
        GetPropType<DarcyTypeTag, Properties::GridVariables>;
    auto darcyGridVariables =
        std::make_shared<DarcyGridVariables>(darcyProblem, darcyGridGeometry);
    darcyGridVariables->init(sol);

    // intialize the vtk output module
    const auto darcyName =
        getParam<std::string>("Problem.Name") + "_" + darcyProblem->name();

    VtkOutputModule<DarcyGridVariables,
                    GetPropType<DarcyTypeTag, Properties::SolutionVector>>
        darcyVtkWriter(*darcyGridVariables, sol, darcyProblem->name());
    using DarcyVelocityOutput =
        GetPropType<DarcyTypeTag, Properties::VelocityOutput>;
    darcyVtkWriter.addVelocityOutput(
        std::make_shared<DarcyVelocityOutput>(*darcyGridVariables));
    GetPropType<DarcyTypeTag, Properties::IOFields>::initOutputModule(
        darcyVtkWriter);
    darcyVtkWriter.write(0.0);

    using FluxVariables = GetPropType<DarcyTypeTag, Properties::FluxVariables>;
    if (couplingParticipant.hasToWriteInitialData()) {
        //TODO
        setInterfaceVelocities<FluxVariables>(*darcyProblem,
                                              *darcyGridVariables, sol, meshNameView, dataNameViewV);
        couplingParticipant.writeQuantityToOtherSolver(meshNameView, dataNameViewV);
    }
    couplingParticipant.setMeshAndInitialize(
        "DarcyMesh", numberOfPoints, coords);// couplingParticipant.initialize();

    // the assembler for a stationary problem
    using Assembler = FVAssembler<DarcyTypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(
        darcyProblem, darcyGridGeometry, darcyGridVariables);

    // the linear solver
    using LinearSolver = UMFPackBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
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

        couplingParticipant.readQuantityFromOtherSolver(meshNameView, dataNameViewP);
    
        // solve the non-linear system
        nonLinearSolver.solve(sol);
        setInterfaceVelocities<FluxVariables>(*darcyProblem,
                                              *darcyGridVariables, sol, meshNameView, dataNameViewV);
        // For testing
        {
            const auto v = couplingParticipant.getQuantityVector(velocityId);
            for (size_t i = 0; i < v.size(); ++i) {
                for (size_t d = 0; d < dim; ++d) {
                    std::cout << coords[i * dim + d] << " ";
                }
                std::cout << "| v[" << i << "]=" << v[i] << std::endl;
            }

            const double sum = std::accumulate(v.begin(), v.end(), 0.);
            std::cout << "Velocities to be sent to ff" << std::endl;
            //          for (size_t i = 0; i < v.size(); ++i) {
            //            std::cout << "v[" << i << "]=" << v[i] << std::endl;
            //          }
            std::cout << "Sum of velocities over boundary to ff: \n"
                      << sum << std::endl;
        }
        couplingParticipant.writeQuantityToOtherSolver(meshNameView, dataNameViewV);

        couplingParticipant.advance(dt);
        preciceDt = couplingParticipant.getMaxTimeStepSize();
        dt = std::min(preciceDt, dt);

        ++iter;

        if (couplingParticipant.hasToReadIterationCheckpoint()) {
            //Read checkpoint
            darcyVtkWriter.write(vtkTime);
            vtkTime += 1.;
            sol = sol_checkpoint;
            darcyGridVariables->update(sol);
            darcyGridVariables->advanceTimeStep();
        } else  // coupling successful
        {
            // write vtk output
            darcyVtkWriter.write(vtkTime);
        }
    }
    // write vtk output
    darcyVtkWriter.write(1.0);

    couplingParticipant.finalize();

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
