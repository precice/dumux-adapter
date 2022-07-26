// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*!
 * \file TODO
 *
 * \brief TODO
 */
#include <config.h>

#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>

#include "dumux-precice/couplingadapter.hh"


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

    const std::string solverName = getParamFromGroup<std::string>"preCICE", "SolverName");
    const std::string preciceConfigFilename =    getParamFromGroup<std::string>"preCICE", "ConfigFileName");
    const std::string meshName =    getParamFromGroup<std::string>"preCICE", "MeshName");


    // Initialize preCICE.Tell preCICE about:
    // - Name of solver
    // - What rank of how many ranks this instance is
    // Configure preCICE. For now the config file is hardcoded.
    // std::string preciceConfigFilename = "precice-config.xml";
    //    if (argc == 3)
    //      preciceConfigFilename = argv[2];
    // if (argc > 2)
    //     preciceConfigFilename = argv[argc - 1];

    auto &couplingInterface = Dumux::Precice::CouplingAdapter::getInstance();
    couplingInterface.announceSolver("FreeFlow", preciceConfigFilename,
                                     mpiHelper.rank(), mpiHelper.size());

    const int dim = couplingInterface.getDimensions();
    // std::cout << dim << "  " << int(FreeFlowGridGeometry::GridView::dimension)
    //           << std::endl;

    // GET mesh corodinates
    std::vector<double> coords;  //( dim * vertexSize );
    std::vector<int> coupledScvfIndices;

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
    const double preciceDt = couplingInterface.setMeshAndInitialize(
        "FreeFlowMesh", numberOfPoints, coords);
    couplingInterface.createIndexMapping(coupledScvfIndices);

    const auto velocityId = couplingInterface.announceScalarQuantity("Velocity");
    const auto pressureId = couplingInterface.announceScalarQuantity("Pressure");


    if (couplingInterface.hasToWriteInitialData()) {
        couplingInterface.writeScalarQuantityToOtherSolver(pressureId);
        couplingInterface.announceInitialDataWritten();
    }
    couplingInterface.initializeData();

    auto dt = preciceDt;
    auto sol_checkpoint = sol;

    double vtkTime = 1.0;
    size_t iter = 0;

    while (couplingInterface.isCouplingOngoing()) {
        if (couplingInterface.hasToWriteIterationCheckpoint()) {
            //DO CHECKPOINTING
            sol_checkpoint = sol;
            couplingInterface.announceIterationCheckpointWritten();
        }

        couplingInterface.readScalarQuantityFromOtherSolver(velocityId);
        // solve the non-linear system
        // Nothing to do in dummy

        // TODO
        setInterfacePressures<FluxVariables>(*freeFlowProblem,
                                             *freeFlowGridVariables, sol);
        couplingInterface.writeScalarQuantityToOtherSolver(pressureId);

        //Read checkpoint
        freeFlowVtkWriter.write(vtkTime);
        vtkTime += 1.;
        const double preciceDt = couplingInterface.advance(dt);
        dt = std::min(preciceDt, dt);

        ++iter;

        if (couplingInterface.hasToReadIterationCheckpoint()) {
            //            //Read checkpoint
            //            freeFlowVtkWriter.write(vtkTime);
            //            vtkTime += 1.;
            sol = sol_checkpoint;
            freeFlowGridVariables->update(sol);
            freeFlowGridVariables->advanceTimeStep();
            //freeFlowGridVariables->init(sol);
            couplingInterface.announceIterationCheckpointRead();
        } else  // coupling successful
        {
            // write vtk output
            freeFlowVtkWriter.write(vtkTime);
        }
    }
    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////

    couplingInterface.finalize();

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
