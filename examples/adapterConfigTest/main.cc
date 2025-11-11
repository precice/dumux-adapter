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

    // Initialize preCICE by reading the configuration in .input file
    auto &couplingParticipant = Dumux::Precice::CouplingAdapter::getInstance();
    couplingParticipant.readAdapterConfig(
        argv[1]);
    const std::string solverName =
        couplingParticipant.getParticipantName();
    const std::string preciceConfigFilename =
        couplingParticipant.getPreCICEConfigFileName();

    couplingParticipant.announceSolver(rank=mpiHelper.rank(), size=mpiHelper.size());
    std::cout << "DUMMY (" << mpiHelper.rank()
              << "): Running solver dummy with preCICE config file \""
              << preciceConfigFilename << "\", participant name \""
              << solverName << "\".\n";

    const int dimensions = couplingParticipant.getMeshDimensions(meshName);
    assert(dimensions == 3);

    const std::vector<std::string> meshNames = couplingParticipant.getMeshNames();
    const std::vector<std::string> dataWriteNames = couplingParticipant.getReadDataNames();
    const std::vector<std::string> dataReadNames = couplingParticipant.getWriteDataNames();

    const int numberOfVertices = 3;

    std::vector<double> writeScalarData(numberOfVertices);
    std::vector<double> readScalarData(numberOfVertices);
    std::vector<double> writeVectorData(numberOfVertices * dimensions);
    std::vector<double> readVectorData(numberOfVertices * dimensions);

    std::vector<double> vertices(numberOfVertices * dimensions);
    std::vector<int> dumuxVertexIDs(numberOfVertices);

    for (int i = 0; i < numberOfVertices; i++) {
        writeScalarData.at(i) = i + numberOfVertices;
        dumuxVertexIDs.at(i) = i + numberOfVertices;
        for (int j = 0; j < dimensions; j++) {
            vertices.at(j + dimensions * i) = i;
            writeVectorData.at(j + dimensions * i) = i;
        }
    }

    std::cout << "adapterConfigTest (" << mpiHelper.rank()
              << "): Initialize preCICE and set mesh\n";
    couplingParticipant.setMesh(meshNames[0], vertices); // only one mesh on each participant

    // Create index mapping between DuMuX's index numbering and preCICE's numbering
    std::cout << "adapterConfigTest (" << mpiHelper.rank() << "): Create index mapping\n";
    couplingParticipant.createIndexMapping(dumuxVertexIDs);

    if (couplingParticipant.requiresToWriteInitialData()) {
        std::cout << "adapterConfigTest (" << mpiHelper.rank()
                  << "): Writing initial data\n";
        couplingParticipant.writeQuantityVector(meshNames[0], dataWriteNames[0],
                                                writeScalarData);
        couplingParticipant.writeQuantityToOtherSolver(meshNames[0],
                                                       dataWriteNames[0]);
        couplingParticipant.writeQuantityVector(meshNames[0], dataWriteNames[1],
                                                writeVectorData);
        couplingParticipant.writeQuantityToOtherSolver(meshNames[0],
                                                       dataWriteNames[1]);
    }
    std::cout << "DUMMY (" << mpiHelper.rank() << "): Exchange initial\n";
    couplingParticipant.initialize();
    double preciceDt = 0;

    // Check exchanged initial data
    if (solverName == "SolverOne") {
        std::cout << "DUMMY (" << mpiHelper.rank()
                  << "): Reading initial data\n";
        couplingParticipant.readQuantityFromOtherSolver(
            meshNames[0], scalarDataReadName, preciceDt);
        couplingParticipant.readQuantityFromOtherSolver(
            meshNames[0], vectorDataReadName, preciceDt);

        const std::vector<double> &readScalarQuantity =
            couplingParticipant.getQuantityVector(meshNames[0], scalarDataReadName);

        std::cout << "DUMMY (" << mpiHelper.rank() << "): Scalar data\n";
        for (const double &value : readScalarQuantity)
            std::cout << value << ",";
        std::cout << "\n";

        const std::vector<double> &readVectorQuantity =
            couplingParticipant.getQuantityVector(meshName, vectorDataReadName);

        std::cout << "DUMMY (" << mpiHelper.rank() << "): Vector data\n";
        for (const double &value : readVectorQuantity)
            std::cout << value << ",";
        std::cout << "\n";

        for (int i = 0; i < numberOfVertices; i++) {
            if (readScalarQuantity.at(i) != writeScalarData.at(i)) {
                std::cout << "DUMMY (" << mpiHelper.rank()
                          << "): Reading initialized SCALAR data error\n"
                          << "Index: " << i << ", Expected "
                          << writeScalarData.at(i) << ", Found "
                          << readScalarQuantity.at(i) << "\n";
                throw(std::runtime_error("Did not find expected SCALAR data."));
            }

            for (int j = 0; j < dimensions; j++) {
                if (readVectorQuantity.at(j + dimensions * i) !=
                    writeVectorData.at(j + dimensions * i)) {
                    std::cout
                        << "DUMMY (" << mpiHelper.rank()
                        << "): Reading initialized VECTOR data error\n"
                        << "Expected " << writeVectorData.at(j + dimensions * i)
                        << ", Found "
                        << readVectorQuantity.at(j + dimensions * i) << "\n";
                    throw(std::runtime_error(
                        "Did not find expected VECTOR data."));
                }
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
        couplingParticipant.readQuantityFromOtherSolver(
            meshName, scalarDataReadName, preciceDt);
        couplingParticipant.readQuantityFromOtherSolver(
            meshName, vectorDataReadName, preciceDt);

        // Check data
        if (iter > 0) {
            int offset = (solverName == "SolverOne") ? 0 : 1;

            const std::vector<double> &readScalarQuantity =
                couplingParticipant.getQuantityVector(meshName,
                                                      scalarDataReadName);
            const std::vector<double> &readVectorQuantity =
                couplingParticipant.getQuantityVector(meshName,
                                                      vectorDataReadName);

            for (int i = 0; i < numberOfVertices; i++) {
                if (readScalarQuantity.at(i) !=
                    writeScalarData.at(i) + offset) {
                    std::cout << "DUMMY (" << mpiHelper.rank()
                              << "): Reading initialized SCALAR data error\n"
                              << "Index " << i << ", Expected "
                              << writeScalarData.at(i) + offset << ", Found "
                              << readScalarQuantity.at(i) << "\n";
                    throw(std::runtime_error(
                        "Did not find expected SCALAR data."));
                }

                for (int j = 0; j < dimensions; j++) {
                    if (readVectorQuantity.at(j + dimensions * i) !=
                        writeVectorData.at(j + dimensions * i) + offset) {
                        std::cout
                            << "DUMMY (" << mpiHelper.rank()
                            << "): Reading initialized VECTOR data error\n"
                            << "Index " << j + dimensions * i << ", Expected "
                            << writeVectorData.at(j + dimensions * i) + offset
                            << ", Found "
                            << readVectorQuantity.at(j + dimensions * i)
                            << "\n";
                        throw(std::runtime_error(
                            "Did not find expected VECTOR data."));
                    }
                }
            }
        }

        ++iter;

        std::cout << "DUMMY (" << mpiHelper.rank() << "): Writing data\n";
        for (int i = 0; i < numberOfVertices; i++) {
            writeScalarData.at(i) = i + iter;
            for (int j = 0; j < dimensions; j++) {
                writeVectorData.at(j + dimensions * i) = i + iter;
            }
        }

        // Write scalar data via DuMuX ID <-> preCICE ID mapping
        for (int i = 0; i < numberOfVertices; i++) {
            const double value = i + iter;
            couplingParticipant.writeScalarQuantityOnFace(
                meshName, scalarDataWriteName, dumuxVertexIDs[i], value);
        }
        couplingParticipant.writeQuantityToOtherSolver(meshName,
                                                       scalarDataWriteName);

        // Write vector data
        couplingParticipant.writeQuantityVector(meshName, vectorDataWriteName,
                                                writeVectorData);
        couplingParticipant.writeQuantityToOtherSolver(meshName,
                                                       vectorDataWriteName);
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
