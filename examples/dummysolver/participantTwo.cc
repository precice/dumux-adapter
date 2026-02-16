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
{
    using namespace Dumux;

    // initialize MPI, finalize is done automatically on exit
    const auto &mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // Initialize preCICE.
    auto &couplingParticipant = Dumux::Precice::CouplingAdapter::getInstance();
    couplingParticipant.announceConfig(mpiHelper.rank(), mpiHelper.size());
    const std::string meshName = couplingParticipant.getMeshNames()[0];

    const int dimensions = couplingParticipant.getMeshDimensions(meshName);
    assert(dimensions == 3);

    const std::string scalarDataWriteName =
        couplingParticipant.getWriteDataNamesOnMesh(meshName)[0];
    const std::string scalarDataReadName =
        couplingParticipant.getReadDataNamesOnMesh(meshName)[0];
    const std::string vectorDataWriteName =
        couplingParticipant.getWriteDataNamesOnMesh(meshName)[1];
    const std::string vectorDataReadName =
        couplingParticipant.getReadDataNamesOnMesh(meshName)[1];

    const int numberOfVertices = 3;

    std::vector<double> writeScalarData(numberOfVertices);
    std::vector<double> readScalarData(numberOfVertices);
    std::vector<double> writeVectorData(numberOfVertices * dimensions);
    std::vector<double> readVectorData(numberOfVertices * dimensions);

    std::vector<double> vertices(numberOfVertices * dimensions);  // coordinates
    std::vector<int> dumuxVertexIDs(numberOfVertices);

    for (int i = 0; i < numberOfVertices; i++) {
        writeScalarData.at(i) = i + numberOfVertices;
        dumuxVertexIDs.at(i) = i + numberOfVertices;
        for (int j = 0; j < dimensions; j++) {
            vertices.at(j + dimensions * i) = i;
            writeVectorData.at(j + dimensions * i) = i;
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
        couplingParticipant.writeQuantityVector(meshName, scalarDataWriteName,
                                                writeScalarData);
        couplingParticipant.writeQuantityToOtherSolver(meshName,
                                                       scalarDataWriteName);
        couplingParticipant.writeQuantityVector(meshName, vectorDataWriteName,
                                                writeVectorData);
        couplingParticipant.writeQuantityToOtherSolver(meshName,
                                                       vectorDataWriteName);
    }
    std::cout << "DUMMY (" << mpiHelper.rank() << "): Exchange initial\n";
    couplingParticipant.initialize();
    double preciceDt = 0;

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
}
