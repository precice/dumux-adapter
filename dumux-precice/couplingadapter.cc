#include "couplingadapter.hh"

#include <algorithm>
#include <cassert>
#include <exception>
#include <limits>

using namespace Dumux::Precice;

CouplingAdapter::CouplingAdapter()
    : wasCreated_(false),
      precice_(nullptr),
      meshWasCreated_(false),
      preciceWasInitialized_(false),
      hasIndexMapper_(false)
{
}

CouplingAdapter &CouplingAdapter::getInstance()
{
    static CouplingAdapter instance;
    return instance;
}

void CouplingAdapter::announceSolver(const std::string &name,
                                     const std::string &configurationFileName,
                                     const int rank,
                                     const int size)
{
    assert(precice_ == nullptr);
    precice_ = std::make_unique<precice::Participant>(
        name, configurationFileName, rank, size);
    wasCreated_ = true;
}

int CouplingAdapter::getMeshDimensions(const std::string &meshName) const
{
    assert(wasCreated_);
    return precice_->getMeshDimensions(meshName);
}

void CouplingAdapter::setMesh(const std::string &meshName,
                              const std::vector<double> &positions)
{
    assert(wasCreated_);
    vertexIDs_.resize(positions.size() / getMeshDimensions(meshName));
    precice_->setMeshVertices(meshName, positions, vertexIDs_);
    meshWasCreated_ = true;
}

void CouplingAdapter::initialize()
{
    assert(wasCreated_);
    assert(meshWasCreated_);
    assert(!preciceWasInitialized_);

    precice_->initialize();

    preciceWasInitialized_ = true;
    assert(preciceWasInitialized_);
}

double CouplingAdapter::getMaxTimeStepSize() const
{
    return precice_->getMaxTimeStepSize();
}

void CouplingAdapter::createIndexMapping(
    const std::vector<int> &dumuxFaceIndices)
{
    assert(meshWasCreated_);
    indexMapper_.createMapping(dumuxFaceIndices, vertexIDs_);
    hasIndexMapper_ = true;
}

void CouplingAdapter::finalize()
{
    assert(wasCreated_);
    if (preciceWasInitialized_)
        precice_->finalize();
}

void CouplingAdapter::advance(const double computedTimeStepLength)
{
    assert(wasCreated_);
    precice_->advance(computedTimeStepLength);
}

bool CouplingAdapter::isCouplingOngoing()
{
    assert(wasCreated_);
    return precice_->isCouplingOngoing();
}

size_t CouplingAdapter::getNumberOfVertices()
{
    assert(wasCreated_);
    return vertexIDs_.size();
}

bool CouplingAdapter::isCoupledEntity(const int faceID) const
{
    assert(wasCreated_);
    return indexMapper_.isDumuxIdMapped(faceID);
}

void CouplingAdapter::print(std::ostream &os)
{
    os << indexMapper_;
}

void CouplingAdapter::readFromPreCICE(const std::string &meshName,
                                      const std::string &dataName,
                                      double relativeReadTime,
                                      std::vector<double> &dataValues)
{
    std::vector<double> readValues(vertexIDs_.size());
    precice_->readData(meshName, dataName, vertexIDs_, relativeReadTime,
                       readValues);

    if (hasIndexMapper_) {
        for (size_t i = 0; i < vertexIDs_.size(); ++i) {
            const auto dumuxId = indexMapper_.getDumuxId(vertexIDs_[i]);
            dataValues[dumuxId] = readValues[i];
        }
    } else {
        dataValues = readValues;
    }
}

void CouplingAdapter::writeToPreCICE(const std::string &meshName,
                                     const std::string &dataName,
                                     std::vector<double> &dataValues)
{
    std::vector<double> writeValues(vertexIDs_.size());

    if (hasIndexMapper_) {
        for (size_t i = 0; i < vertexIDs_.size(); ++i) {
            const auto dumuxId = indexMapper_.getDumuxId(vertexIDs_[i]);
            writeValues[i] = dataValues[dumuxId];
        }
    } else {
        writeValues = dataValues;
    }
    precice_->writeData(meshName, dataName, vertexIDs_, writeValues);
}

bool CouplingAdapter::requiresToWriteInitialData()
{
    assert(wasCreated_);
    return precice_->requiresInitialData();
}

bool CouplingAdapter::writeCheckpointIfRequired()
{
    assert(wasCreated_);
    if (!precice_->requiresWritingCheckpoint()) {
        return false;
    }
    for (auto &state : states_) {
        state->writeState();
    }
    return true;
}

bool CouplingAdapter::readCheckpointIfRequired()
{
    assert(wasCreated_);
    if (!precice_->requiresReadingCheckpoint()) {
        return false;
    }
    for (auto &state : states_) {
        state->readState();
    }
    return true;
}

CouplingAdapter::~CouplingAdapter() {}
