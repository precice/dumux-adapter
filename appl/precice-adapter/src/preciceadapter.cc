#include "../include/preciceadapter.hh"

#include <algorithm>
#include <exception>
#include <cassert>

using namespace precice_adapter;

PreciceAdapter::PreciceAdapter():
  wasCreated_(false), precice_(nullptr), meshWasCreated_(false), preciceWasInitialized_(false), hasIndexMapper_(false),
  meshID_(0), timeStepSize_(0.)
{
  interfaceData_.reserve(reserveSize_);
}

PreciceAdapter& PreciceAdapter::getInstance()
{
  static PreciceAdapter instance;
  return instance;
}

void PreciceAdapter::announceSolver(const std::string &name, const std::string configurationFileName, const int rank, const int size)
{
  assert( precice_ == nullptr );
  precice_ = std::make_unique<precice::SolverInterface>(name, rank, size);
  wasCreated_ = true;
  precice_->configure( configurationFileName );
}

size_t PreciceAdapter::announceQuantity(const std::string &name)
{
  return announceQuantity<1>( name );
}

int PreciceAdapter::getDimensions() const
{
  assert( wasCreated_ );
  return precice_->getDimensions();
}
/*
void PreciceAdapter::setMeshName(const std::string& meshName)
{
  assert( wasCreated_ );
  meshID_ = precice_->getMeshID(meshName);
}
*/


void PreciceAdapter::setMesh(const std::string& meshName,
                             const size_t numPoints,
                             std::vector<double>& coordinates)
{
  assert( wasCreated_ );
  assert( numPoints == coordinates.size() / getDimensions() );
  meshID_ = precice_->getMeshID(meshName);
  vertexIDs_.resize( numPoints );
  precice_->setMeshVertices( meshID_, numPoints, coordinates.data(), vertexIDs_.data() );
  meshWasCreated_ = true;
}

double PreciceAdapter::initialize()
{
  assert( wasCreated_ );
  assert( meshWasCreated_ );
  assert( !preciceWasInitialized_ );

  timeStepSize_ = precice_->initialize();
  assert( timeStepSize_ > 0 );

  preciceWasInitialized_ = true;
  return timeStepSize_;
}

 void PreciceAdapter::createIndexMapping( const std::vector<int>& dumuxFaceIDs )
 {
   assert( meshWasCreated_ );
   indexMapper_.createMapping( dumuxFaceIDs, vertexIDs_);
   hasIndexMapper_ = true;
 }

double PreciceAdapter::setMeshAndInitialize(const std::string& meshName,
                                            const size_t numPoints,
                                            std::vector<double>& coordinates)
{
  setMesh( meshName, numPoints, coordinates );
  return initialize();
}

void PreciceAdapter::initializeData()
{
  assert( preciceWasInitialized_ );
  precice_->initializeData();
}

void PreciceAdapter::finalize()
{
  assert( wasCreated_ );
  if (preciceWasInitialized_) precice_->finalize();
}

double PreciceAdapter::advance( const double computedTimeStepLength )
{
  assert( wasCreated_ );
  return precice_->advance( computedTimeStepLength );
}

bool PreciceAdapter::isCouplingOngoing()
{
  assert( wasCreated_ );
  return precice_->isCouplingOngoing();
}

size_t PreciceAdapter::getNumberOfVertices()
{
  assert( wasCreated_ );
  return vertexIDs_.size();
}

double PreciceAdapter::getScalarQuantityOnFace(const size_t dataID, const int faceID) const
{
  assert( wasCreated_ );
  assert( hasIndexMapper_ );
  if ( !hasIndexMapper_ )
  {
    throw std::runtime_error("Reading quantity using faceID, but index mapping was not created!");
  }
  const auto idx = indexMapper_.getPreciceId( faceID );
  assert( dataID < interfaceData_.size() );
  const std::vector<double>& quantityVector = interfaceData_[ dataID ].data;
  assert(idx < quantityVector.size() );
  return quantityVector[idx];
}

void PreciceAdapter::writeScalarQuantityOnFace(const size_t dataID, const int faceID, const double value)
{
  assert( wasCreated_ );
  assert( hasIndexMapper_ );
  if ( !hasIndexMapper_ )
  {
    throw std::runtime_error("Writing quantity using faceID, but index mapping was not created!");
  }
  const auto idx = indexMapper_.getPreciceId( faceID );
  assert( dataID < interfaceData_.size() );
  std::vector<double>& quantityVector = interfaceData_[ dataID ].data;
  assert( idx < quantityVector.size() );
  quantityVector[idx] = value;
}

//void PreciceAdapter::writeVectorQuantityOnFace(const size_t dataID,
//                                               const int faceID,
//                                               const double* value,
//                                               const size_t size)
//{
//  assert( wasCreated_ );
//  assert( hasIndexMapper_ );
//  assert( size == getDimensions() );
//  if ( !hasIndexMapper_ )
//  {
//    throw std::runtime_error("Writing quantity using faceID, but index mapping was not created!");
//  }
//  const auto idx = indexMapper_.getPreciceId( faceID ) * size;
//  assert( dataID < interfaceData_.size() );
//  std::vector<double>& quantityVector = interfaceData_[ dataID ].data;
//  assert( idx < quantityVector.size() );
//  //quantityVector[idx] = value;
//  std::copy_n( value, size, quantityVector[idx] );
//}

std::vector<double>& PreciceAdapter::getQuantityVector(const size_t dataID )
{
  assert( wasCreated_ );
  assert( dataID < interfaceData_.size() );
  return interfaceData_[dataID].data;
}

const std::vector<double> &PreciceAdapter::getQuantityVector(const size_t dataID ) const
{
  assert( wasCreated_ );
  return getQuantityVector( dataID );
}

void PreciceAdapter::writeScalarQuantityVector(const size_t dataID,
                                               std::vector<double> &values)
{
 assert( wasCreated_);
 assert( dataID < interfaceData_.size() );
 assert( interfaceData_[dataID].data.size() == values.size() );
 interfaceData_[dataID].data = values;
}

void PreciceAdapter::writeScalarQuantityToOtherSolver(const size_t dataID)
{
  assert( wasCreated_ );
  assert( dataID < interfaceData_.size() );
  assert( dataID < std::numeric_limits<int>::max() );
  writeBlockScalarDataToPrecice( interfaceData_[dataID].preciceDataID, interfaceData_[dataID].data );
}

void PreciceAdapter::readScalarQuantityFromOtherSolver(const size_t dataID)
{
  assert( wasCreated_ );
  assert( dataID < interfaceData_.size() );
  assert( dataID < std::numeric_limits<int>::max() );
  readBlockScalarDataFromPrecice( interfaceData_[dataID].preciceDataID, interfaceData_[dataID].data );
}

bool PreciceAdapter::isCoupledEntity(const int faceID) const
{
  assert( wasCreated_ );
  return indexMapper_.isDumuxIdMapped( faceID );
}

size_t PreciceAdapter::getIdFromName(const std::string &dataName) const
{
  assert( wasCreated_ );
  // Check if data is actually there
  const auto it = getDataIterator( dataName );
  if ( it == interfaceData_.end() )
  {
    throw( std::runtime_error(" Error! Name of data not found! ") );
  }
  const auto idx = std::distance( interfaceData_.begin(), it );
  assert( idx > -1 );
  return size_t(idx);
}

std::string PreciceAdapter::getNameFromId(const size_t dataID) const
{
  assert( wasCreated_ );
  assert( dataID < interfaceData_.size() );
  return interfaceData_[dataID].name;
}

void PreciceAdapter::print(std::ostream& os)
{
  os << indexMapper_;
}

bool PreciceAdapter::checkIfActionIsRequired( const std::string& condition )
{
  assert( wasCreated_ );
  return precice_->isActionRequired( condition );
}

void PreciceAdapter::actionIsFulfilled(const std::string& condition)
{
  assert( wasCreated_ );
  precice_->fulfilledAction( condition );
}

void PreciceAdapter::readBlockScalarDataFromPrecice(const int dataID, std::vector<double> &data)
{
  assert( wasCreated_ );
  assert( vertexIDs_.size() == data.size() );
  precice_->readBlockScalarData( dataID, vertexIDs_.size(), vertexIDs_.data(), data.data() );
}

void PreciceAdapter::writeBlockScalarDataToPrecice(const int dataID, std::vector<double> &data)
{
  assert( wasCreated_ );
  assert( vertexIDs_.size() == data.size() );
  precice_->writeBlockScalarData( dataID, vertexIDs_.size(), vertexIDs_.data(), data.data() );
}

bool PreciceAdapter::hasToWriteInitialData()
{
  assert( wasCreated_ );
  return checkIfActionIsRequired(precice::constants::actionWriteInitialData());
}

void PreciceAdapter::announceInitialDataWritten()
{
  assert( wasCreated_ );
  precice_->fulfilledAction( precice::constants::actionWriteInitialData() );
}

bool PreciceAdapter::hasToReadIterationCheckpoint()
{
  assert( wasCreated_ );
  return checkIfActionIsRequired(precice::constants::actionReadIterationCheckpoint());
}

void PreciceAdapter::announceIterationCheckpointRead()
{
  assert( wasCreated_ );
  actionIsFulfilled( precice::constants::actionReadIterationCheckpoint() );
}

bool PreciceAdapter::hasToWriteIterationCheckpoint()
{
  assert( wasCreated_ );
  return checkIfActionIsRequired(precice::constants::actionWriteIterationCheckpoint());
}

void PreciceAdapter::announceIterationCheckpointWritten()
{
  assert( wasCreated_ );
  actionIsFulfilled( precice::constants::actionWriteIterationCheckpoint() );
}

PreciceAdapter::~PreciceAdapter()
{
}

