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
      hasIndexMapper_(false),
      timeStepSize_(0.)
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

void CouplingAdapter::announceQuantity(const precice::string_view &meshName,
                                       const precice::string_view &dataName)
{
    assert(meshWasCreated_);
    const std::string key = createKeyFromName(meshName, dataName);
    if (dataMap_.find(key) != dataMap_.end()) {
        throw(std::runtime_error(" Error! Duplicate quantity announced! "));
    }

    int dataDimension = precice_->getDataDimensions(meshName, dataName);
    std::vector<double> dataValues(vertexIDs_.size() * dataDimension);
    dataMap_.insert(std::make_pair(key, dataValues));
}

int CouplingAdapter::getMeshDimensions(
    const precice::string_view &meshName) const
{
    assert(wasCreated_);
    return precice_->getMeshDimensions(meshName);
}

void CouplingAdapter::setMesh(const precice::string_view &meshName,
                              precice::span<const double> positions)
{
    assert(wasCreated_);
    vertexIDsSpan_ = precice::span(vertexIDs_);
    precice_->setMeshVertices(meshName, positions, vertexIDsSpan_);
    meshWasCreated_ = true;
}

void CouplingAdapter::initialize()
{
    assert(wasCreated_);
    assert(meshWasCreated_);
    assert(!preciceWasInitialized_);

    precice_->initialize();
    timeStepSize_ = precice_->getMaxTimeStepSize();
    assert(timeStepSize_ > 0);

    preciceWasInitialized_ = true;
    assert(preciceWasInitialized_);
}

double CouplingAdapter::getMaxTimeStepSize()
{
    return precice_->getMaxTimeStepSize();
}

void CouplingAdapter::createIndexMapping(
    const std::vector<int> &dumuxFaceIndices)  // TODO what does this do?
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
    return precice_->advance(computedTimeStepLength);
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

double CouplingAdapter::getScalarQuantityOnFace(
    const precice::string_view &meshName,
    const precice::string_view &dataName,
    const int faceID)
{
    assert(wasCreated_);
    assert(hasIndexMapper_);
    if (!hasIndexMapper_) {
        throw std::runtime_error(
            "Reading quantity using faceID, but index mapping was not "
            "created!");
    }
    const auto idx = indexMapper_.getPreciceId(faceID);
    std::vector<double> &dataVector = getQuantityVector(meshName, dataName);
    assert(idx < dataVector.size());
    return dataVector[idx];
}

void CouplingAdapter::writeScalarQuantityOnFace(
    const precice::string_view &meshName,
    const precice::string_view &dataName,
    const int faceID,
    const double value)
{
    assert(wasCreated_);
    assert(hasIndexMapper_);
    if (!hasIndexMapper_) {
        throw std::runtime_error(
            "Writing quantity using faceID, but index mapping was not "
            "created!");
    }
    const auto idx = indexMapper_.getPreciceId(faceID);
    std::vector<double> &dataVector = getQuantityVector(meshName, dataName);
    assert(idx < dataVector.size());
    dataVector[idx] = value;
}

std::vector<double> &CouplingAdapter::getQuantityVector(
    const precice::string_view &meshName,
    const precice::string_view &dataName)
{
    std::string key = createKeyFromName(meshName, dataName);
    assert(dataMap_.find(key) != dataMap_.end());
    return dataMap_[key];
}

void CouplingAdapter::writeQuantityVector(const precice::string_view &meshName,
                                          const precice::string_view &dataName,
                                          std::vector<double> &values)
{
    std::vector<double> &dataVector = getQuantityVector(meshName, dataName);
    assert(dataVector.size() == values.size());
    dataVector = values;
}

bool CouplingAdapter::isCoupledEntity(const int faceID) const
{
    assert(wasCreated_);
    return indexMapper_.isDumuxIdMapped(faceID);
}

std::string CouplingAdapter::createKeyFromName(
    const precice::string_view &meshName,
    const precice::string_view &dataName) const
{
    assert(wasCreated_);
    std::string combinedKey;

    for (int i = 0; i < (meshName.size() + 1 + dataName.size()); i++) {
        if (i < meshName.size())
            combinedKey += meshName[i];
        else if (i == meshName.size())
            combinedKey += ':';
        else
            combinedKey += dataName[i - meshName.size()];
    }

    return combinedKey;
}

void CouplingAdapter::print(std::ostream &os)
{
    os << indexMapper_;
}

void CouplingAdapter::readQuantityFromOtherSolver(
    const precice::string_view &meshName,
    const precice::string_view &dataName,
    double relativeReadTime)
{
    precice::span<double> dataValuesSpan(getQuantityVector(meshName, dataName));
    precice_->readData(meshName, dataName, vertexIDsSpan_, relativeReadTime,
                       dataValuesSpan);
}

void CouplingAdapter::writeQuantityToOtherSolver(
    const precice::string_view &meshName,
    const precice::string_view &dataName)
{
    precice::span<const double> dataValuesSpan(
        getQuantityVector(meshName, dataName));
    precice_->writeData(meshName, dataName, vertexIDsSpan_, dataValuesSpan);
}

bool CouplingAdapter::hasToWriteInitialData()
{
    assert(wasCreated_);
    return precice_->requiresInitialData();
}

bool CouplingAdapter::hasToReadIterationCheckpoint()
{
    assert(wasCreated_);
    return precice_->requiresReadingCheckpoint();
}

bool CouplingAdapter::hasToWriteIterationCheckpoint()
{
    assert(wasCreated_);
    return precice_->requiresWritingCheckpoint();
}
CouplingAdapter::~CouplingAdapter() {}
