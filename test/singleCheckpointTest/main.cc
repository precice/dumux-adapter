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

#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/parameters.hh>

#include "dumux-precice/couplingadapter.hh"

#include <memory>
#include <numeric>
#include <type_traits>

class TimeLoop
{
private:
    double currentTime_{0.0};
    double currentDt_{0.0};
    long currentStep_{0};

public:
    double time() const { return currentTime_; }
    long timeStepIndex() const { return currentStep_; }
    double timeStepSize() const { return currentDt_; }
    void advanceTime(double dt)
    {
        currentTime_ += dt;
        currentDt_ = dt;
        ++currentStep_;
    }
    void setTime(double time, long step)
    {
        currentTime_ = time;
        currentStep_ = step;
    }
    void setTimeStepSize(double dt) { currentDt_ = dt; }
};
// Mock GridVariables class to test checkpointing functionality
class GridVariables
{
    bool updated{false};

public:
    void update(const std::vector<double> & /*x*/) { updated = true; }
};

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

    auto &couplingParticipant = Dumux::Precice::CouplingAdapter::getInstance();
    couplingParticipant.announceConfig(mpiHelper.rank(), mpiHelper.size());

    const std::string meshName =
        (couplingParticipant.getSolverName() == "SolverOne") ? "MeshOne"
                                                             : "MeshTwo";
    const int dimensions = couplingParticipant.getMeshDimensions(meshName);
    assert(dimensions == 3);

    const int numberOfVertices = 3;

    std::vector<double> writeScalarData(numberOfVertices);
    std::vector<double> readScalarData(numberOfVertices);

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

    if (couplingParticipant.requiresToWriteInitialData()) {
        std::cout << "DUMMY (" << mpiHelper.rank()
                  << "): Writing initial data\n";
        couplingParticipant.writeQuantityVector(
            meshName, couplingParticipant.getWriteDataNamesOnMesh(meshName)[0],
            writeScalarData);
        couplingParticipant.writeQuantityToOtherSolver(
            meshName, couplingParticipant.getWriteDataNamesOnMesh(meshName)[0]);
    }
    std::cout << "DUMMY (" << mpiHelper.rank() << "): Exchange initial\n";
    couplingParticipant.initialize();
    double preciceDt = 0;

    // Create instances of the file-scope mock types for checkpointing.
    GridVariables gridVars;
    TimeLoop timeLoop;

    // Register writeScalarData as the solver state for checkpointing.
    couplingParticipant.initializeCheckpoint(writeScalarData, gridVars,
                                             timeLoop);

    // Check exchanged initial data
    if (couplingParticipant.getSolverName() == "SolverOne") {
        std::cout << "SolverOne: Reading initial data\n";
        couplingParticipant.readQuantityFromOtherSolver(
            meshName, couplingParticipant.getReadDataNamesOnMesh(meshName)[0],
            preciceDt);
    }

    int iter = 0;
    std::vector<double> dataToKeep(numberOfVertices);
    double timeToKeep = 0.0;
    long timeStepIndexToKeep = 0;
    double timeStepSizeToKeep = 0.0;

    while (couplingParticipant.isCouplingOngoing()) {
        if (couplingParticipant.getSolverName() == "SolverOne") {
            if (couplingParticipant.writeCheckpointIfRequired()) {
                // Keep data for later comparison
                dataToKeep = writeScalarData;
                timeToKeep = timeLoop.time();
                timeStepIndexToKeep = timeLoop.timeStepIndex();
                timeStepSizeToKeep = timeLoop.timeStepSize();
            }

            ++iter;

            std::iota(writeScalarData.begin(), writeScalarData.end(), iter);

            preciceDt = couplingParticipant.getMaxTimeStepSize();
            couplingParticipant.advance(preciceDt);
            timeLoop.advanceTime(preciceDt);
            timeLoop.setTimeStepSize(preciceDt / 2.);
            if (couplingParticipant.readCheckpointIfRequired()) {
                if (writeScalarData != dataToKeep) {
                    throw std::runtime_error(
                        "SolverOne: Checkpointing failed, data not restored "
                        "correctly");
                }
                if ((timeStepIndexToKeep != timeLoop.timeStepIndex()) ||
                    (timeToKeep != timeLoop.time()) ||
                    (timeStepSizeToKeep != timeLoop.timeStepSize())) {
                    throw std::runtime_error(
                        "SolverOne: Checkpointing failed, time step not "
                        "restored correctly");
                }
            }
        } else {
            couplingParticipant.writeCheckpointIfRequired();
            preciceDt = couplingParticipant.getMaxTimeStepSize();
            couplingParticipant.advance(preciceDt);
            couplingParticipant.readCheckpointIfRequired();
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
