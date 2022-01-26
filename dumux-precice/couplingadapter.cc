#include "couplingadapter.hh"

#include <algorithm>
#include <cassert>
#include <exception>

using namespace Dumux::Precice;

CouplingAdapter::CouplingAdapter()
    : _wasCreated(false),
      _precice(nullptr),
      _meshWasCreated(false),
      _preciceWasInitialized(false),
      _hasIndexMapper(false),
      _meshID(0),
      _timeStepSize(0.)
{
    _preciceDataID.reserve(_reserveSize);
    _dataNames.reserve(_reserveSize);
    _dataVectors.reserve(_reserveSize);
}

CouplingAdapter &CouplingAdapter::getInstance()
{
    static CouplingAdapter instance;
    return instance;
}

void CouplingAdapter::announceSolver(const std::string &name,
                                     const std::string configurationFileName,
                                     const int rank,
                                     const int size)
{
    assert(_precice == nullptr);
    _precice = std::make_unique<precice::SolverInterface>(
        name, configurationFileName, rank, size);
    _wasCreated = true;
}

size_t CouplingAdapter::announceQuantity(const std::string &name)
{
    assert(_meshWasCreated);
    auto it = std::find(_dataNames.begin(), _dataNames.end(), name);
    if (it != _dataNames.end()) {
        throw(std::runtime_error(" Error! Duplicate quantity announced! "));
    }
    _dataNames.push_back(name);
    _preciceDataID.push_back(_precice->getDataID(name, _meshID));
    _dataVectors.push_back(std::vector<double>(_vertexIDs.size()));

    return getNumberOfQuantities() - 1;
}

int CouplingAdapter::getDimensions() const
{
    assert(_wasCreated);
    return _precice->getDimensions();
}
/*
void CouplingAdapter::setMeshName(const std::string& meshName)
{
  assert( _wasCreated );
  _meshID = _precice->getMeshID(meshName);
}
*/

void CouplingAdapter::setMesh(const std::string &meshName,
                              const size_t numPoints,
                              std::vector<double> &coordinates)
{
    assert(_wasCreated);
    assert(numPoints == coordinates.size() / getDimensions());
    _meshID = _precice->getMeshID(meshName);
    _vertexIDs.resize(numPoints);
    _precice->setMeshVertices(_meshID, numPoints, coordinates.data(),
                              _vertexIDs.data());
    _meshWasCreated = true;
}

double CouplingAdapter::initialize()
{
    assert(_wasCreated);
    assert(_meshWasCreated);
    assert(!_preciceWasInitialized);

    _timeStepSize = _precice->initialize();
    assert(_timeStepSize > 0);

    _preciceWasInitialized = true;
    return _timeStepSize;
}

void CouplingAdapter::createIndexMapping(const std::vector<int> &dumuxFaceIDs)
{
    assert(_meshWasCreated);
    _indexMapper.createMapping(dumuxFaceIDs, _vertexIDs);
    _hasIndexMapper = true;
}

double CouplingAdapter::setMeshAndInitialize(const std::string &meshName,
                                             const size_t numPoints,
                                             std::vector<double> &coordinates)
{
    setMesh(meshName, numPoints, coordinates);
    return initialize();
}

void CouplingAdapter::initializeData()
{
    assert(_preciceWasInitialized);
    _precice->initializeData();
}

void CouplingAdapter::finalize()
{
    assert(_wasCreated);
    if (_preciceWasInitialized)
        _precice->finalize();
}

double CouplingAdapter::advance(const double computedTimeStepLength)
{
    assert(_wasCreated);
    return _precice->advance(computedTimeStepLength);
}

bool CouplingAdapter::isCouplingOngoing()
{
    assert(_wasCreated);
    return _precice->isCouplingOngoing();
}

size_t CouplingAdapter::getNumberOfVertices()
{
    assert(_wasCreated);
    return _vertexIDs.size();
}

double CouplingAdapter::getScalarQuantityOnFace(const size_t dataID,
                                                const int faceID) const
{
    assert(_wasCreated);
    assert(_hasIndexMapper);
    if (!_hasIndexMapper) {
        throw std::runtime_error(
            "Reading quantity using faceID, but index mapping was not "
            "created!");
    }
    const auto idx = _indexMapper.getPreciceId(faceID);
    assert(dataID < _dataVectors.size());
    const std::vector<double> &quantityVector = _dataVectors[dataID];
    assert(idx < quantityVector.size());
    return quantityVector[idx];
}

void CouplingAdapter::writeScalarQuantityOnFace(const size_t dataID,
                                                const int faceID,
                                                const double value)
{
    assert(_wasCreated);
    assert(_hasIndexMapper);
    if (!_hasIndexMapper) {
        throw std::runtime_error(
            "Writing quantity using faceID, but index mapping was not "
            "created!");
    }
    const auto idx = _indexMapper.getPreciceId(faceID);
    assert(dataID < _dataVectors.size());
    std::vector<double> &quantityVector = _dataVectors[dataID];
    assert(idx < quantityVector.size());
    quantityVector[idx] = value;
}

//void CouplingAdapter::writeVectorQuantityOnFace(const size_t dataID,
//                                               const int faceID,
//                                               const double* value,
//                                               const size_t size)
//{
//  assert( _wasCreated );
//  assert( _hasIndexMapper );
//  assert( size == getDimensions() );
//  if ( !_hasIndexMapper )
//  {
//    throw std::runtime_error("Writing quantity using faceID, but index mapping was not created!");
//  }
//  const auto idx = _indexMapper.getPreciceId( faceID ) * size;
//  assert( dataID < _dataVectors.size() );
//  std::vector<double>& quantityVector = _dataVectors[ dataID ];
//  assert( idx < quantityVector.size() );
//  //quantityVector[idx] = value;
//  std::copy_n( value, size, quantityVector[idx] );
//}

std::vector<double> &CouplingAdapter::getQuantityVector(const size_t dataID)
{
    assert(_wasCreated);
    assert(dataID < _dataVectors.size());
    return _dataVectors[dataID];
}

const std::vector<double> &CouplingAdapter::getQuantityVector(
    const size_t dataID) const
{
    assert(_wasCreated);
    return getQuantityVector(dataID);
}

void CouplingAdapter::writeScalarQuantityVector(const size_t dataID,
                                                std::vector<double> &values)
{
    assert(_wasCreated);
    assert(dataID < _dataVectors.size());
    assert(_dataVectors[dataID].size() == values.size());
    _dataVectors[dataID] = values;
}

void CouplingAdapter::writeScalarQuantityToOtherSolver(const size_t dataID)
{
    assert(_wasCreated);
    assert(dataID < _dataVectors.size());
    assert(dataID < _preciceDataID.size());
    assert(dataID < std::numeric_limits<int>::max());
    writeBlockScalarDataToPrecice(_preciceDataID[dataID], _dataVectors[dataID]);
}

void CouplingAdapter::readScalarQuantityFromOtherSolver(const size_t dataID)
{
    assert(_wasCreated);
    assert(dataID < _dataVectors.size());
    assert(dataID < _preciceDataID.size());
    assert(dataID < std::numeric_limits<int>::max());
    readBlockScalarDataFromPrecice(_preciceDataID[dataID],
                                   _dataVectors[dataID]);
}

bool CouplingAdapter::isCoupledEntity(const int faceID) const
{
    assert(_wasCreated);
    return _indexMapper.isDumuxIdMapped(faceID);
}

size_t CouplingAdapter::getIdFromName(const std::string &dataName) const
{
    assert(_wasCreated);
    const auto it = std::find(_dataNames.begin(), _dataNames.end(), dataName);
    if (it == _dataNames.end()) {
        throw(std::runtime_error(" Error! Name of data not found! "));
    }
    const auto idx = std::distance(_dataNames.begin(), it);
    assert(idx > -1);
    return size_t(idx);
}

std::string CouplingAdapter::getNameFromId(const size_t dataID) const
{
    assert(_wasCreated);
    assert(dataID < _dataNames.size());
    return _dataNames[dataID];
}

void CouplingAdapter::print(std::ostream &os)
{
    os << _indexMapper;
}

bool CouplingAdapter::checkIfActionIsRequired(const std::string &condition)
{
    assert(_wasCreated);
    return _precice->isActionRequired(condition);
}

void CouplingAdapter::actionIsFulfilled(const std::string &condition)
{
    assert(_wasCreated);
    _precice->markActionFulfilled(condition);
}

void CouplingAdapter::readBlockScalarDataFromPrecice(const int dataID,
                                                     std::vector<double> &data)
{
    assert(_wasCreated);
    assert(_vertexIDs.size() == data.size());
    _precice->readBlockScalarData(dataID, _vertexIDs.size(), _vertexIDs.data(),
                                  data.data());
}

void CouplingAdapter::writeBlockScalarDataToPrecice(const int dataID,
                                                    std::vector<double> &data)
{
    assert(_wasCreated);
    assert(_vertexIDs.size() == data.size());
    _precice->writeBlockScalarData(dataID, _vertexIDs.size(), _vertexIDs.data(),
                                   data.data());
}

bool CouplingAdapter::hasToWriteInitialData()
{
    assert(_wasCreated);
    return checkIfActionIsRequired(
        precice::constants::actionWriteInitialData());
}

void CouplingAdapter::announceInitialDataWritten()
{
    assert(_wasCreated);
    _precice->markActionFulfilled(precice::constants::actionWriteInitialData());
}

bool CouplingAdapter::hasToReadIterationCheckpoint()
{
    assert(_wasCreated);
    return checkIfActionIsRequired(
        precice::constants::actionReadIterationCheckpoint());
}

void CouplingAdapter::announceIterationCheckpointRead()
{
    assert(_wasCreated);
    actionIsFulfilled(precice::constants::actionReadIterationCheckpoint());
}

bool CouplingAdapter::hasToWriteIterationCheckpoint()
{
    assert(_wasCreated);
    return checkIfActionIsRequired(
        precice::constants::actionWriteIterationCheckpoint());
}

void CouplingAdapter::announceIterationCheckpointWritten()
{
    assert(_wasCreated);
    actionIsFulfilled(precice::constants::actionWriteIterationCheckpoint());
}

CouplingAdapter::~CouplingAdapter() {}
