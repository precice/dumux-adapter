#include "precicewrapper.hh"

#include <cassert>

using namespace precice_wrapper;

PreciceWrapper::PreciceWrapper():
  wasCreated_(false), precice_(nullptr), meshWasCreated_(false), preciceWasInitialized_(false),
  meshID_(0), freeFlowHeatFluxID_(0), solidHeatFluxID_(0), timeStepSize_(0.), writeHeatFluxType_(UNDEFINED), readHeatFluxType_(UNDEFINED)
{

}

PreciceWrapper& PreciceWrapper::getInstance()
{
  static PreciceWrapper instance;
  return instance;
}

void PreciceWrapper::configure( const std::string& configurationFileName )
{
  precice_->configure( configurationFileName );
}

void PreciceWrapper::announceHeatFluxToWrite(const HeatFluxType heatFluxType)
{
  writeHeatFluxType_ = heatFluxType;
}

void PreciceWrapper::announceHeatFluxToRead(const HeatFluxType heatFluxType)
{
  readHeatFluxType_ = heatFluxType;
}

void PreciceWrapper::announceSolver( const std::string& name, const int rank, const int size )
{
  assert( precice_ == nullptr );
  precice_ = std::make_unique<precice::SolverInterface>(name, rank, size);
  wasCreated_ = true;
}

int PreciceWrapper::getDimensions()
{
  assert( wasCreated_ );
  return precice_->getDimensions();
}
/*
void PreciceWrapper::setMeshName(const std::string& meshName)
{
  assert( wasCreated_ );
  meshID_ = precice_->getMeshID(meshName);
}
*/
void PreciceWrapper::setMesh(const std::string& meshName,
                             const size_t numPoints,
                              std::vector<double>& coordinates,
                              const std::vector<int>& dumuxFaceIDs )
{
  assert( wasCreated_ );
  assert( numPoints == dumuxFaceIDs.size() );
  meshID_ = precice_->getMeshID(meshName);
  vertexIDs_.resize( numPoints );
  precice_->setMeshVertices( meshID_, numPoints, coordinates.data(), vertexIDs_.data() );
  indexMapper_.createMapping( dumuxFaceIDs, vertexIDs_);
  meshWasCreated_ = true;
}
/*
int PreciceWrapper::getDataID( const std::string& dataName, const int meshID )
{
  assert( wasCreated_ );
  return precice_->getDataID( dataName, meshID );
}
*/
double PreciceWrapper::initialize()
{
  assert( wasCreated_ );
  assert( meshWasCreated_ );

  solidHeatFluxID_ = precice_->getDataID( "Solid-Heat-Flux", meshID_ );
  freeFlowHeatFluxID_ = precice_->getDataID( "FreeFlow-Heat-Flux", meshID_ );

  freeFlowHeatFlux_.resize( getNumberOfVertices() );
  solidHeatFlux_.resize( getNumberOfVertices() );

  timeStepSize_ = precice_->initialize();
  assert( timeStepSize_ > 0 );

  precice_->initializeData();
  preciceWasInitialized_ = true;
  return timeStepSize_;
}

void PreciceWrapper::finalize()
{
  assert( wasCreated_ );
  precice_->finalize();
}

/*
void PreciceWrapper::initializeData()
{
  assert( wasCreated_ );
  precice_->initializeData();
}
*/

double PreciceWrapper::advance( const double computedTimeStepLength )
{
  assert( wasCreated_ );
  return precice_->advance( computedTimeStepLength );
}

bool PreciceWrapper::isCouplingOngoing()
{
  assert( wasCreated_ );
  return precice_->isCouplingOngoing();
}

size_t PreciceWrapper::getNumberOfVertices()
{
  assert( wasCreated_ );
  return vertexIDs_.size();
}

double PreciceWrapper::getHeatFluxAtFace( const int faceID) const
{
  assert( wasCreated_ );
  const auto idx = indexMapper_.getPreciceId( faceID );
  assert( readHeatFluxType_ != HeatFluxType::UNDEFINED );
  if (readHeatFluxType_ == HeatFluxType::FreeFlow)
  {
    assert(idx < freeFlowHeatFlux_.size() );
    return freeFlowHeatFlux_[idx];
  }
  else
  {
    assert(idx < solidHeatFlux_.size() );
    return solidHeatFlux_[idx];
  }
}

void PreciceWrapper::writeHeatFluxOnFace(const int faceID,
                                         const double value)
{
  assert( wasCreated_ );
  const auto idx = indexMapper_.getPreciceId( faceID );
  assert( writeHeatFluxType_ != HeatFluxType::UNDEFINED );
  if ( writeHeatFluxType_ == HeatFluxType::FreeFlow )
  {
    freeFlowHeatFlux_[idx] = value;
  }
  else
  {
    solidHeatFlux_[idx] = value;
  }
}

void PreciceWrapper::writeHeatFluxToOtherSolver()
{
  assert( wasCreated_ );
  assert( writeHeatFluxType_ != HeatFluxType::UNDEFINED );

  if ( writeHeatFluxType_ == HeatFluxType::FreeFlow )
    writeBlockScalarDataToPrecice( freeFlowHeatFluxID_, freeFlowHeatFlux_ );
  else
    writeBlockScalarDataToPrecice( solidHeatFluxID_, solidHeatFlux_ );

}

void PreciceWrapper::readHeatFluxFromOtherSolver()
{
  assert( wasCreated_ );
  assert( readHeatFluxType_ != HeatFluxType::UNDEFINED );

  if ( readHeatFluxType_ == HeatFluxType::FreeFlow )
    readBlockScalarDataFromPrecice( freeFlowHeatFluxID_, freeFlowHeatFlux_ );
  else
    readBlockScalarDataFromPrecice( solidHeatFluxID_, solidHeatFlux_ );
}

bool PreciceWrapper::isCoupledEntity(const int faceID) const
{
  assert( wasCreated_ );
  return indexMapper_.isDumuxIdMapped( faceID );
}

std::vector<double>& PreciceWrapper::getHeatFluxToWrite()
{
  assert( wasCreated_ );
  assert( heatFluxToWrite_ != HeatFluxType::UNDEFINED);
  if ( writeHeatFluxType_ == HeatFluxType::FreeFlow )
    return freeFlowHeatFlux_;
  else
    return solidHeatFlux_;
}

//void PreciceWrapper::readScalarQuantitiy(const int dataID, std::vector<double> &data)
//{
//  assert( wasCreated_ );
//  precice_->readBlockScalarData( dataID, vertexIDs_.size(),
//                                       vertexIDs_.data(), data.data() );
//}
//
//void PreciceWrapper::writeScalarQuantitiy(const int dataID, std::vector<double> &data)
//{
//  assert( wasCreated_ );
//  precice_->writeBlockScalarData( dataID, vertexIDs_.size(),
//                                        vertexIDs_.data(), data.data() );
//}

void PreciceWrapper::print(std::ostream& os)
{
  os << indexMapper_;
}

/*
void PreciceWrapper::writeSolidTemperature(std::vector<double> &temperature)
{
  assert( wasCreated_ );
  precice_->writeBlockScalarData( solidTemperatureID_, vertexIDs_.size(),
                                       vertexIDs_.data(), temperature.data() );

}

void PreciceWrapper::readSolidTemperature(std::vector<double> &temperature)
{
  assert( wasCreated_ );
  precice_->readBlockScalarData( solidTemperatureID_, vertexIDs_.size(),
                                       vertexIDs_.data(), temperature.data() );
}
*/

/*
void PreciceWrapper::readBlockScalarData( const int dataID,
                                          const int size,
                                          int* const valueIndices,
                                          double* const values )
{
  assert( wasCreated_ );
  precice_->readBlockScalarData( dataID, size, valueIndices, values );
}

void PreciceWrapper::writeBlockScalarData( const int dataID,
                                           const int size,
                                           int* const valueIndices,
                                           double* const values )
{
  assert( wasCreated_ );
  return precice_->writeBlockScalarData( dataID, size, valueIndices, values );
}
*/

bool PreciceWrapper::checkIfActionIsRequired( const std::string& condition )
{
  assert( wasCreated_ );
  return precice_->isActionRequired( condition );
}

void PreciceWrapper::actionIsFulfilled(const std::string& condition)
{
  assert( wasCreated_ );
  precice_->fulfilledAction( condition );
}

void PreciceWrapper::readBlockScalarDataFromPrecice(const int dataID, std::vector<double> &data)
{
  assert( wasCreated_ );
  assert( vertexIDs.size() == data.size() );
  precice_->readBlockScalarData( dataID, vertexIDs_.size(), vertexIDs_.data(), data.data() );
}

void PreciceWrapper::writeBlockScalarDataToPrecice(const int dataID, std::vector<double> &data)
{
  assert( wasCreated_ );
  assert( vertexIDs.size() == data.size() );
  precice_->writeBlockScalarData( dataID, vertexIDs_.size(), vertexIDs_.data(), data.data() );
}

bool PreciceWrapper::hasToWriteInitialData()
{
  assert( wasCreated_ );
  return checkIfActionIsRequired(precice::constants::actionWriteInitialData());
}


bool PreciceWrapper::hasToReadIterationCheckpoint()
{
  assert( wasCreated_ );
  return checkIfActionIsRequired(precice::constants::actionReadIterationCheckpoint());
}

void PreciceWrapper::announceIterationCheckpointRead()
{
  assert( wasCreated_ );
  actionIsFulfilled( precice::constants::actionWriteIterationCheckpoint() );
}

bool PreciceWrapper::hasToWriteIterationCheckpoint()
{
  assert( wasCreated_ );
  return checkIfActionIsRequired(precice::constants::actionWriteIterationCheckpoint());
}

void PreciceWrapper::announceIterationCheckpointWritten()
{
  assert( wasCreated_ );
  actionIsFulfilled( precice::constants::actionWriteIterationCheckpoint() );
}

/*
void PreciceWrapper::writeInitialBlockScalarData( const int dataID,
                                                  const int size,
                                                  int* const valueIndices,
                                                  double* const values )
{
  assert( wasCreated_ );
  if ( hasToWriteInitialData() )
  {
    precice_->writeBlockScalarData( dataID, size, valueIndices, values );
  }
}

void PreciceWrapper::announceAllInitialDataWritten()
{
  assert( wasCreated_ );
  if ( hasToWriteInitialData() )
  {
    precice_->fulfilledAction( precice::constants::actionWriteInitialData() );
  }
}
*/
PreciceWrapper::~PreciceWrapper()
{
}

