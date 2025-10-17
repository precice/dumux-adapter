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

#include <memory>
#include <numeric>
#include <type_traits>

// Mock TimeLoop class to test checkpointing functionality
class TimeLoop
{
public:
    double t{0.0};
    long idx{0};
    double time() const { return t; }
    long timeStepIndex() const { return idx; }
    void setTime(double tt, long i)
    {
        t = tt;
        idx = i;
    }
    void setTimeStepSize(double /*dt*/) {}
};

// Mock GridVariables class to test checkpointing functionality
class GridVariables
{
public:
    bool updated{false};
    bool advanced{false};
    void update(const std::vector<double> & /*x*/) { updated = true; }
    void advanceTimeStep() { advanced = true; }
};

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
    const std::string dataToWrite =
        (solverName == "SolverOne") ? "dataOne" : "dataTwo";
    const std::string dataToRead =
        (solverName == "SolverOne") ? "dataTwo" : "dataOne";

    const int numberOfVertices = 3;

    std::vector<double> writeScalarData(numberOfVertices);
    std::vector<double> readScalarData(numberOfVertices);
    std::vector<double> dataToKeep(numberOfVertices);

    std::vector<double> vertices(numberOfVertices * dimensions);  // coordinates
    std::vector<int> dumuxVertexIDs(numberOfVertices);

    // initialize writeScalarData and dumuxVertexIDs with consecutive values
    std::iota(writeScalarData.begin(), writeScalarData.end(), numberOfVertices);
    std::iota(dumuxVertexIDs.begin(), dumuxVertexIDs.end(), numberOfVertices);
    // set vertex coordinates: for each vertex i fill its `dimensions` entries with i
    for (int i = 0; i < numberOfVertices; ++i) {
        std::fill_n(vertices.begin() + i * dimensions, dimensions,
                    static_cast<double>(i));
    }

    std::cout << "DUMMY (" << mpiHelper.rank()
              << "): Initialize preCICE and set mesh\n";
    couplingParticipant.setMesh(meshName, vertices);

    // Create index mapping between DuMuX's index numbering and preCICE's numbering
    std::cout << "DUMMY (" << mpiHelper.rank() << "): Create index mapping\n";
    couplingParticipant.createIndexMapping(dumuxVertexIDs);

    couplingParticipant.announceQuantity(meshName, dataToWrite);
    couplingParticipant.announceQuantity(meshName, dataToRead);

    if (couplingParticipant.requiresToWriteInitialData()) {
        std::cout << "DUMMY (" << mpiHelper.rank()
                  << "): Writing initial data\n";
        couplingParticipant.writeQuantityVector(meshName, dataToWrite,
                                                writeScalarData);
        couplingParticipant.writeQuantityToOtherSolver(meshName, dataToWrite);
    }
    std::cout << "DUMMY (" << mpiHelper.rank() << "): Exchange initial\n";
    couplingParticipant.initialize();
    double preciceDt = 0;

    // Create instances of the file-scope mock types and register them for
    // checkpointing. Types are declared in the anonymous namespace above so
    // template instantiation will see their definitions.
    TimeLoop timeLoop;
    GridVariables gridVars;

    // Register writeScalarData as the solver state for checkpointing.
    couplingParticipant.initializeCheckpoint(writeScalarData, timeLoop,
                                             gridVars);

    // Check exchanged initial data
    if (solverName == "SolverOne") {
        std::cout << "SolverOne: Reading initial data\n";
        couplingParticipant.readQuantityFromOtherSolver(meshName, dataToRead,
                                                        preciceDt);
    }

    int iter = 0;
    dataToKeep = writeScalarData;

    while (couplingParticipant.isCouplingOngoing()) {
        if (solverName == "SolverOne") {
            couplingParticipant.writeCheckpointIfRequired();

            ++iter;

            std::iota(writeScalarData.begin(), writeScalarData.end(), iter);

            preciceDt = couplingParticipant.getMaxTimeStepSize();
            couplingParticipant.advance(preciceDt);
            couplingParticipant.readCheckpointIfRequired(preciceDt);
            if (writeScalarData != dataToKeep) {
                throw std::runtime_error(
                    "SolverOne: Checkpointing failed, data not restored "
                    "correctly");
            }
        } else {
            couplingParticipant.writeCheckpointIfRequired();
            preciceDt = couplingParticipant.getMaxTimeStepSize();
            couplingParticipant.advance(preciceDt);
            couplingParticipant.readCheckpointIfRequired(preciceDt);
        }
    }
    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////

    couplingParticipant.finalize();
    std::cout << "DUMMY (" << mpiHelper.rank()
              << "): Closing single checkpoint test.\n";

    // print dumux end message
    if (mpiHelper.rank() == 0) {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }

    return 0;
}  // end main
// namespace
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
