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

    const std::string solverName = getParamFromGroup<std::string>("preCICE", "SolverName");
    const std::string preciceConfigFilename =    getParamFromGroup<std::string>("preCICE", "ConfigFileName");
    const std::string meshName =    getParamFromGroup<std::string>("preCICE", "MeshName");


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
    couplingInterface.announceSolver(solverName, preciceConfigFilename,
                                     mpiHelper.rank(), mpiHelper.size());


    std::cout << "DUMMY (" << mpiHelper.rank()
              << "): Running solver dummy with preCICE config file \"" << preciceConfigFilename
              << "\", participant name \"" << solverName
              << "\", and mesh name \"" << meshName
              << "\".\n";

    const int dimensions = couplingInterface.getDimensions();
    assert( dimensions == 3 );
    const std::string scalarDataWriteName = (solverName == "SolverOne") ? "scalarDataOne" : "scalarDataTwo" ;
    const std::string scalarDataReadName = (solverName == "SolverOne") ? "scalarDataOne" : "scalarDataTwo" ;

    const std::string vectorDataWriteName = (solverName == "SolverOne") ? "vectorDataOne" : "vectorDataTwo" ;
    const std::string vectorDataReadName = (solverName == "SolverOne") ? "vectorDataOne" : "vectorDataTwo" ;

    const int numberOfVertices = 3;

    std::vector<double> readScalarData(numberOfVertices);
    std::vector<double> writeScalarData(numberOfVertices);
    std::vector<double> readVectorData(numberOfVertices * dimensions);
    std::vector<double> writeVectorData(numberOfVertices * dimensions);
    std::vector<double> vertices(numberOfVertices * dimensions);
    std::vector<int>    preciceVertexIDs(numberOfVertices);
    std::vector<int>    dumuxVertexIDs(numberOfVertices);

    for (int i = 0; i < numberOfVertices; i++) {
        writeScalarData.at(i) = i + numberOfVertices;
        dumuxVertexIDs.at(i) = i + numberOfVertices;
        for (int j = 0; j < dimensions; j++) {
            vertices.at(j + dimensions * i)  = i;
            writeVectorData.at(j + dimensions * i) = i;
        }
    }

    double preciceDt = couplingInterface.setMeshAndInitialize(meshName, numberOfVertices, vertices);

    // Create index mapping between DuMuX's index numbering and preCICE's numbering
    couplingInterface.createIndexMapping(dumuxVertexIDs);

    // const int readScalarDataID  = couplingInterface.getDataID(scalarDataReadName, meshID);
    // const int writeScalarDataID = couplingInterface.getDataID(scalarDataWriteName, meshID);

    // const int readVectorDataID  = couplingInterface.getDataID(vectorDataReadName, meshID);
    // const int writeVectorDataID = couplingInterface.getDataID(vectorDataWriteName, meshID);

    const int readScalarDataID  = couplingInterface.announceScalarQuantity(scalarDataReadName);
    const int writeScalarDataID = couplingInterface.announceScalarQuantity(scalarDataWriteName);

    const int readVectorDataID  = couplingInterface.announceVectorQuantity(vectorDataReadName);
    const int writeVectorDataID = couplingInterface.announceVectorQuantity(vectorDataWriteName);

    // const auto velocityId = couplingInterface.announceScalarQuantity("Velocity");
    // const auto pressureId = couplingInterface.announceScalarQuantity("Pressure");

    if (couplingInterface.hasToWriteInitialData()) {
        // Scalar data
        couplingInterface.writeQuantityVector(writeScalarDataID, writeScalarData);
        couplingInterface.writeScalarQuantityToOtherSolver(writeScalarDataID);
        // Vector data
        couplingInterface.writeQuantityVector(writeVectorDataID, writeVectorData);
        couplingInterface.writeScalarQuantityToOtherSolver(writeVectorDataID);
        couplingInterface.announceInitialDataWritten();
    }
    couplingInterface.initializeData();

    // Check exchanged initial data
    // {
    //     for (int i = 0; i < numberOfVertices; i++) {
    //         if ( readScalarData.at(i) != writeScalarData.at(i) ) {
    //             std::cout << "DUMMY (" << mpiHelper.rank()
    //                       << "): Reading initialized SCALAR data error\n"
    //                       << "Expected " << writeScalarData.at(i)
    //                       << ", Found " << readScalarData.at(i)
    //                       << "\n";
    //         }

    //         for (int j = 0; j < dimensions; j++) {
    //             if ( readVectorData.at(j + dimensions * i) != writeVectorData.at(j + dimensions * i) ) {
    //                 std::cout << "DUMMY (" << mpiHelper.rank()
    //                           << "): Reading initialized VECTOR data error\n"
    //                           << "Expected " << writeVectorData.at(j + dimensions * i)
    //                           << ", Found " << readVectorData.at(j + dimensions * i)
    //                           << "\n";
    //             }
    //         }
    //     }
    // }
    std::fill( readScalarData.begin(), readScalarData.end(), 0 );
    std::fill( readVectorData.begin(), readVectorData.end(), 0 );

    int iter = 0;

    while (couplingInterface.isCouplingOngoing()) {
        if (couplingInterface.hasToWriteIterationCheckpoint()) {
            std::cout << "DUMMY (" << mpiHelper.rank() << "): Writing iteration checkpoint\n";
            couplingInterface.announceIterationCheckpointWritten();
        }

        //Read data
        couplingInterface.readQuantityFromOtherSolver( readScalarDataID, Dumux::Precice::QuantityType::Scalar );
        couplingInterface.readQuantityFromOtherSolver( readVectorDataID, Dumux::Precice::QuantityType::Vector );

        // Check data
        for (int i = 0; i < numberOfVertices; i++) {
            if ( readScalarData.at(i) != writeScalarData.at(i) ) {
                std::cout << "DUMMY (" << mpiHelper.rank()
                          << "): Reading initialized SCALAR data error\n"
                          << "Expected " << writeScalarData.at(i)
                          << ", Found " << readScalarData.at(i)
                          << "\n";
            }

            for (int j = 0; j < dimensions; j++) {
                if ( readVectorData.at(j + dimensions * i) != writeVectorData.at(j + dimensions * i) ) {
                    std::cout << "DUMMY (" << mpiHelper.rank()
                              << "): Reading initialized VECTOR data error\n"
                              << "Expected " << writeVectorData.at(j + dimensions * i)
                              << ", Found " << readVectorData.at(j + dimensions * i)
                              << "\n";
                }
            }
        }

        ++iter;

        for (int i = 0; i < numberOfVertices; i++) {
            writeScalarData.at(i) = i + numberOfVertices + iter;
            for (int j = 0; j < dimensions; j++) {
                writeVectorData.at(j + dimensions * i) = i + iter;
            }
        }

        // Write scalar data via DuMuX ID <-> preCICE ID mapping
        for (int i = 0; i < numberOfVertices; i++) {
            const double value = i + iter;
            couplingInterface.writeScalarQuantityOnFace(writeScalarDataID,
                                dumuxVertexIDs[i],
                                value);
        }
        couplingInterface.writeScalarQuantityToOtherSolver(writeScalarDataID);

        // Write vector data
        couplingInterface.writeQuantityVector(writeVectorDataID, writeVectorData);
        couplingInterface.writeScalarQuantityToOtherSolver(writeVectorDataID);

        preciceDt = couplingInterface.advance(preciceDt);

        if (couplingInterface.hasToReadIterationCheckpoint()) {
            std::cout << "DUMMY (" << mpiHelper.rank() << "): Reading iteration checkpoint\n";
            couplingInterface.announceIterationCheckpointRead();
        } else
        {
            std::cout << "DUMMY (" << mpiHelper.rank() << "): Advancing in time\n";
        }
    }
    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////

    couplingInterface.finalize();
      std::cout  << "DUMMY (" << mpiHelper.rank() << "): Closing C++ solver dummy...\n";

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
} catch (...) {
    std::cerr << "Unknown exception thrown! ---> Abort!" << std::endl;
    return 4;
}
