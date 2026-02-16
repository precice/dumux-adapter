#include "couplingadapter.hh"
#include <dumux/common/parameters.hh>

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

void CouplingAdapter::announceConfig(const int rank, const int size)
{
    preciceConfigName_ = Dumux::getParamFromGroup<std::string>(
        "precice-adapter-config", "precice_config_file_path");
    participantName_ = Dumux::getParamFromGroup<std::string>(
        "precice-adapter-config", "participant_name");

    assert(precice_ == nullptr);
    precice_ = std::make_unique<precice::Participant>(
        participantName_, preciceConfigName_, rank, size);
    wasCreated_ = true;

    int interfaceIndex = 0;
    do {
        interfaceIndex++;

        std::string meshNameTag =
            "interfaces." + std::to_string(interfaceIndex) + ".mesh_name";

        if (!Dumux::hasParamInGroup("precice-adapter-config", meshNameTag)) {
            break;
        }

        const std::string meshName = Dumux::getParamFromGroup<std::string>(
            "precice-adapter-config", meshNameTag);

        int dataTag = 0;
        do {
            dataTag++;
            const std::string readDataTag =
                "interfaces." + std::to_string(interfaceIndex) + ".read_data." +
                std::to_string(dataTag) + ".name";

            if (!Dumux::hasParamInGroup("precice-adapter-config",
                                        readDataTag)) {
                break;
            }

            const std::string readDataName =
                Dumux::getParamFromGroup<std::string>("precice-adapter-config",
                                                      readDataTag);
            auto key = std::make_pair(meshName, readDataName);
            dataRead_.try_emplace(key);
        } while (true);

        dataTag = 0;
        do {
            dataTag++;
            const std::string writeDataTag =
                "interfaces." + std::to_string(interfaceIndex) +
                ".write_data." + std::to_string(dataTag) + ".name";

            if (!Dumux::hasParamInGroup("precice-adapter-config",
                                        writeDataTag)) {
                break;
            }

            const std::string writeDataName =
                Dumux::getParamFromGroup<std::string>("precice-adapter-config",
                                                      writeDataTag);
            auto key = std::make_pair(meshName, writeDataName);
            dataWrite_.try_emplace(key);
        } while (true);
    } while (true);
}

void CouplingAdapter::announceSolver(const std::string &name,
                                     const std::string &configurationFileName,
                                     const int rank,
                                     const int size)
{
    participantName_ = name;
    preciceConfigName_ = configurationFileName;

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

    // compute size of data vectors for coupling data on this mesh
    auto dataToReadOnMesh = getReadDataNamesOnMesh(meshName);
    auto dataToWriteOnMesh = getWriteDataNamesOnMesh(meshName);

    for (auto dataName : dataToReadOnMesh) {
        int dataDimension = precice_->getDataDimensions(meshName, dataName);
        dataRead_[std::make_pair(meshName, dataName)].resize(vertexIDs_.size() *
                                                             dataDimension);
    }

    for (auto dataName : dataToWriteOnMesh) {
        int dataDimension = precice_->getDataDimensions(meshName, dataName);
        dataWrite_[std::make_pair(meshName, dataName)].resize(
            vertexIDs_.size() * dataDimension);
    }
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

std::string CouplingAdapter::getSolverName() const
{
    assert(wasCreated_);
    return participantName_;
}

std::vector<std::string> CouplingAdapter::getMeshNames() const
{
    assert(wasCreated_);
    std::vector<std::string> meshNames;

    for (const auto &[key, value] : dataRead_) {
        auto it = std::find(meshNames.begin(), meshNames.end(), key.first);
        if (it == meshNames.end()) {
            meshNames.push_back(key.first);
        }
    }
    for (const auto &[key, value] : dataWrite_) {
        auto it = std::find(meshNames.begin(), meshNames.end(), key.first);
        if (it == meshNames.end()) {
            meshNames.push_back(key.first);
        }
    }
    return meshNames;
}

std::vector<std::string> CouplingAdapter::getReadDataNamesOnMesh(
    const std::string &meshName) const
{
    assert(wasCreated_);
    std::vector<std::string> readNames;

    for (const auto &[key, value] : dataRead_) {
        if (key.first == meshName) {
            readNames.push_back(key.second);
        }
    }
    return readNames;
}

std::vector<std::string> CouplingAdapter::getWriteDataNamesOnMesh(
    const std::string &meshName) const
{
    assert(wasCreated_);
    std::vector<std::string> writeNames;

    for (const auto &[key, value] : dataWrite_) {
        if (key.first == meshName) {
            writeNames.push_back(key.second);
        }
    }

    return writeNames;
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

bool CouplingAdapter::isTimeWindowComplete()
{
    assert(wasCreated_);
    return precice_->isTimeWindowComplete();
}

size_t CouplingAdapter::getNumberOfVertices()
{
    assert(wasCreated_);
    return vertexIDs_.size();
}

double CouplingAdapter::getScalarQuantityOnFace(const std::string &meshName,
                                                const std::string &dataName,
                                                const int faceID)
{
    assert(wasCreated_);
    assert(hasIndexMapper_);
    if (!hasIndexMapper_) {
        throw std::runtime_error(
            "Reading quantity using faceID, but index mapping was not "
            "created!");
    }
    auto key = std::make_pair(meshName, dataName);
    std::vector<double> &dataVector = dataRead_[key];

    const auto idx = indexMapper_.getPreciceId(faceID);
    assert(idx < dataVector.size());

    return dataVector[idx];
}

void CouplingAdapter::writeScalarQuantityOnFace(const std::string &meshName,
                                                const std::string &dataName,
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

    auto key = std::make_pair(meshName, dataName);
    std::vector<double> &dataVector = dataWrite_[key];

    const auto idx = indexMapper_.getPreciceId(faceID);
    assert(idx < dataVector.size());

    dataVector[idx] = value;
}

void CouplingAdapter::writeQuantityVector(const std::string &meshName,
                                          const std::string &dataName,
                                          const std::vector<double> &values)
{
    auto key = std::make_pair(meshName, dataName);
    std::vector<double> &dataVector = dataWrite_[key];
    assert(dataVector.size() == values.size());
    dataVector = values;
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

void CouplingAdapter::readQuantityFromOtherSolver(const std::string &meshName,
                                                  const std::string &dataName,
                                                  double relativeReadTime)
{
    auto key = std::make_pair(meshName, dataName);
    for (std::map<std::pair<std::string, std::string>,
                  std::vector<double>>::const_iterator it = dataRead_.begin();
         it != dataRead_.end(); ++it) {
        std::cout << it->first.first << " " << it->first.second << " "
                  << it->second.size() << "\n";
    }
    std::vector<double> &dataVector = dataRead_[key];

    precice_->readData(meshName, dataName, vertexIDs_, relativeReadTime,
                       dataVector);
}

void CouplingAdapter::writeQuantityToOtherSolver(const std::string &meshName,
                                                 const std::string &dataName)
{
    auto key = std::make_pair(meshName, dataName);
    std::vector<double> &dataVector = dataWrite_[key];
    precice_->writeData(meshName, dataName, vertexIDs_, dataVector);
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
