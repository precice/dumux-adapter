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
#include <dune/istl/io.hh>

#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/parameters.hh>

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

    // Initialize preCICE. Tell preCICE about:
    // - Name of solver
    // - Configuration file name
    // - Solver rank
    const std::string solverName =
        getParamFromGroup<std::string>("preCICE", "SolverName");
    const std::string preciceConfigFilename =
        getParamFromGroup<std::string>("preCICE", "ConfigFileName");
    const std::string meshName =
        getParamFromGroup<std::string>("preCICE", "MeshName");

    auto &couplingParticipant = Dumux::Precice::CouplingAdapter::getInstance();
    couplingParticipant.announceSolver(solverName, preciceConfigFilename,
                                       mpiHelper.rank(), mpiHelper.size());
    std::cout << "DUMMY (" << mpiHelper.rank()
              << "): Running solver dummy with preCICE config file \""
              << preciceConfigFilename << "\", participant name \""
              << solverName << "\", and mesh name \"" << meshName << "\".\n";

    const int dimensions = couplingParticipant.getMeshDimensions(meshName);
    assert(dimensions == 3);
    const std::string scalarDataWriteName =
        (solverName == "SolverOne") ? "scalarDataOne" : "scalarDataTwo";
    const std::string scalarDataReadName =
        (solverName == "SolverOne") ? "scalarDataTwo" : "scalarDataOne";

    const int numberOfVertices = 3;

    std::vector<double> writeScalarData(numberOfVertices);
    std::vector<double> writeScalarDataInDumux(
        numberOfVertices *
        2);  // dumux IDs "i+numberOfVertices" is mapped to preCICE IDs "i"

    std::vector<double> vertices(numberOfVertices * dimensions);  // coordinates
    std::vector<int> dumuxVertexIDs(numberOfVertices);

    for (int i = 0; i < numberOfVertices; i++) {
        writeScalarData.at(i) = i + numberOfVertices;
        dumuxVertexIDs.at(i) = i + numberOfVertices;
        writeScalarDataInDumux.at(i + numberOfVertices) = i + numberOfVertices;
        for (int j = 0; j < dimensions; j++) {
            vertices.at(j + dimensions * i) = i;
        }
    }

    std::cout << "DUMMY (" << mpiHelper.rank()
              << "): Initialize preCICE and set mesh\n";
    couplingParticipant.setMesh(meshName, vertices);

    // Create index mapping between DuMuX's index numbering and preCICE's numbering
    std::cout << "DUMMY (" << mpiHelper.rank() << "): Create index mapping\n";
    couplingParticipant.createIndexMapping(dumuxVertexIDs);

    if (couplingParticipant.requiresToWriteInitialData()) {
        std::cout << "DUMMY (" << mpiHelper.rank()
                  << "): Writing initial data\n";
        couplingParticipant.writeToPreCICE(meshName, scalarDataWriteName,
                                           writeScalarDataInDumux);
    }
    std::cout << "DUMMY (" << mpiHelper.rank() << "): Exchange initial\n";
    couplingParticipant.initialize();
    double preciceDt = 0;

    // Check exchanged initial data
    if (solverName == "SolverOne") {
        std::cout << "DUMMY (" << mpiHelper.rank()
                  << "): Reading initial data\n";

        std::vector<double> readScalarQuantityInDumux(numberOfVertices * 2);

        couplingParticipant.readFromPreCICE(
            meshName, scalarDataReadName, preciceDt, readScalarQuantityInDumux);

        std::cout << "DUMMY (" << mpiHelper.rank() << "): Scalar data\n";
        for (const double &value : readScalarQuantityInDumux)
            std::cout << value << ",";
        std::cout << "\n";

        for (int i = 0; i < numberOfVertices; i++) {
            if (readScalarQuantityInDumux.at(i + numberOfVertices) !=
                writeScalarData.at(i)) {
                std::cout << "DUMMY (" << mpiHelper.rank()
                          << "): Reading initialized SCALAR data error\n"
                          << "Index: " << i << ", Expected "
                          << writeScalarData.at(i) << ", Found "
                          << readScalarQuantityInDumux.at(i + numberOfVertices)
                          << "\n";
                throw(std::runtime_error("Did not find expected SCALAR data."));
            }
        }
    }

    int iter = 0;

    while (couplingParticipant.isCouplingOngoing()) {
        if (couplingParticipant.writeCheckpointIfRequired()) {
            std::cout << "DUMMY (" << mpiHelper.rank()
                      << "): Writing iteration checkpoint\n";
        }

        //Read data
        std::cout << "DUMMY (" << mpiHelper.rank() << "): Reading data\n";
        std::vector<double> readScalarQuantityInDumux(numberOfVertices * 2);
        couplingParticipant.readFromPreCICE(
            meshName, scalarDataReadName, preciceDt, readScalarQuantityInDumux);

        // Check data
        if (iter > 0) {
            int offset = (solverName == "SolverOne") ? 0 : 1;

            for (int i = 0; i < numberOfVertices; i++) {
                if (readScalarQuantityInDumux.at(i + numberOfVertices) !=
                    writeScalarData.at(i) + offset) {
                    std::cout
                        << "DUMMY (" << mpiHelper.rank()
                        << "): Reading SCALAR data error\n"
                        << "Index " << i << ", Expected "
                        << writeScalarData.at(i) + offset << ", Found "
                        << readScalarQuantityInDumux.at(i + numberOfVertices)
                        << "\n";
                    throw(std::runtime_error(
                        "Did not find expected SCALAR data."));
                }
            }
        }

        ++iter;

        std::cout << "DUMMY (" << mpiHelper.rank() << "): Writing data\n";
        // Write scalar data via DuMuX ID <-> preCICE ID mapping
        std::vector<double> toWriteData(numberOfVertices * 2);
        for (int i = 0; i < numberOfVertices; i++) {
            const double value = i + iter;
            toWriteData.at(dumuxVertexIDs[i]) = value;
            writeScalarData.at(i) = value;
        }
        couplingParticipant.writeToPreCICE(meshName, scalarDataWriteName,
                                           toWriteData);

        preciceDt = couplingParticipant.getMaxTimeStepSize();
        couplingParticipant.advance(preciceDt);

        if (couplingParticipant.readCheckpointIfRequired()) {
            std::cout << "DUMMY (" << mpiHelper.rank()
                      << "): Reading iteration checkpoint\n";
        } else {
            std::cout << "DUMMY (" << mpiHelper.rank()
                      << "): Advancing in time\n";
        }
    }
    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////

    couplingParticipant.finalize();
    std::cout << "DUMMY (" << mpiHelper.rank()
              << "): Closing C++ solver dummy...\n";

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
} catch (Dune::Exception &e) {
    std::cerr << "Dune reported error: " << e << " ---> Abort!" << std::endl;
    return 3;
} catch (std::runtime_error &e) {
    std::cerr << std::endl << e.what() << " ---> Abort!" << std::endl;
    return 4;
} catch (...) {
    std::cerr << "Unknown exception thrown! ---> Abort!" << std::endl;
    return 5;
}
