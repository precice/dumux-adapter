#include "preciceadapter.hh"

#include <cassert>

using namespace precice_adapter;

PreciceAdapter::PreciceAdapter():
  wasCreated_(false), precice_(nullptr), meshWasCreated_(false), preciceWasInitialized_(false),
  meshID_(0), heatFluxID_(0), temperatureID_(0), timeStepSize_(0.)
{

}

PreciceAdapter& PreciceAdapter::getInstance()
{
  static PreciceAdapter instance;
  return instance;
}

void PreciceAdapter::announceSolver( const std::string& name,
                                     const std::string& configurationFileName,
                                     const int rank,
                                     const int size )
{
  assert( precice_ == nullptr );
  precice_ = std::make_unique<precice::SolverInterface>(name, rank, size);
  precice_->configure( configurationFileName );
  wasCreated_ = true;
}

int PreciceAdapter::getDimensions()
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
/*
void PreciceAdapter::setMesh(const std::string& meshName,
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
*/
/*
int PreciceAdapter::getDataID( const std::string& dataName, const int meshID )
{
  assert( wasCreated_ );
  return precice_->getDataID( dataName, meshID );
}
*/
double PreciceAdapter::initialize()
{
  assert( wasCreated_ );
  assert( meshWasCreated_ );

  heatFluxID_ = precice_->getDataID( "Heat-Flux", meshID_ );
  temperatureID_ = precice_->getDataID( "Temperature", meshID_ );

  heatFlux_.resize( getNumberOfVertices() );
  temperature_.resize( getNumberOfVertices() );

  timeStepSize_ = precice_->initialize();
  assert( timeStepSize_ > 0 );

  preciceWasInitialized_ = true;
  return timeStepSize_;
}

void PreciceAdapter::initializeData()
{
  assert( preciceWasInitialized_ );
  precice_->initializeData();
}

void PreciceAdapter::finalize()
{
  assert( wasCreated_ );
  precice_->finalize();
}

/*
void PreciceAdapter::initializeData()
{
  assert( wasCreated_ );
  precice_->initializeData();
}
*/

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

double PreciceAdapter::getHeatFluxOnFace( const int faceID) const
{
  assert( wasCreated_ );
  const auto idx = indexMapper_.getPreciceId( faceID );
  assert(idx < heatFlux_.size() );
  return heatFlux_[idx];
}

void PreciceAdapter::writeHeatFluxOnFace(const int faceID,
                                         const double value)
{
  assert( wasCreated_ );
  const auto idx = indexMapper_.getPreciceId( faceID );
  assert(idx < heatFlux_.size() );
  heatFlux_[idx] = value;
}

double PreciceAdapter::getTemperatureOnFace(const int faceID) const
{
  assert( wasCreated_ );
  const auto idx = indexMapper_.getPreciceId( faceID );
  assert(idx < temperature_.size() );
  return temperature_[idx];
}

void PreciceAdapter::writeTemperatureOnFace(const int faceID, const double value)
{
  assert( wasCreated_ );
  const auto idx = indexMapper_.getPreciceId( faceID );
  assert(idx < temperature_.size() );
  temperature_[idx] = value;
}

void PreciceAdapter::writeHeatFluxToOtherSolver()
{
  assert( wasCreated_ );
  writeBlockScalarDataToPrecice( heatFluxID_, heatFlux_ );
  for (auto v: heatFlux_ )
  {
    std::cout << "Written heat-flux is :" << v << std::endl;
  }
}

void PreciceAdapter::readHeatFluxFromOtherSolver()
{
  assert( wasCreated_ );
  readBlockScalarDataFromPrecice( heatFluxID_, heatFlux_ );
  for (auto v: heatFlux_ )
  {
    std::cout << "Read heat-flux is :" << v << std::endl;
  }
}

void PreciceAdapter::writeTemperatureToOtherSolver()
{
  assert( wasCreated_ );
  writeBlockScalarDataToPrecice( temperatureID_, temperature_ );
  for (auto v: temperature_ )
  {
    std::cout << "Written temperature is :" << v << std::endl;
  }
}

void PreciceAdapter::readTemperatureFromOtherSolver()
{
  assert( wasCreated_ );
  readBlockScalarDataFromPrecice( temperatureID_, temperature_ );
  for (auto v: temperature_ )
  {
    std::cout << "Read temperature is :" << v << std::endl;
  }
}

bool PreciceAdapter::isCoupledEntity(const int faceID) const
{
  assert( wasCreated_ );
  return indexMapper_.isDumuxIdMapped( faceID );
}
/*
std::vector<double>& PreciceAdapter::getHeatFluxToWrite()
{
  assert( wasCreated_ );
  assert( writeHeatFluxType_ != HeatFluxType::UNDEFINED);
  if ( writeHeatFluxType_ == HeatFluxType::FreeFlow )
    return freeFlowHeatFlux_;
  else
    return solidHeatFlux_;
}
*/
//void PreciceAdapter::readScalarQuantitiy(const int dataID, std::vector<double> &data)
//{
//  assert( wasCreated_ );
//  precice_->readBlockScalarData( dataID, vertexIDs_.size(),
//                                       vertexIDs_.data(), data.data() );
//}
//
//void PreciceAdapter::writeScalarQuantitiy(const int dataID, std::vector<double> &data)
//{
//  assert( wasCreated_ );
//  precice_->writeBlockScalarData( dataID, vertexIDs_.size(),
//                                        vertexIDs_.data(), data.data() );
//}

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

double PreciceAdapter::setMeshAndInitialize(const std::string &meshName, const size_t numPoints, std::vector<double> &coordinates, const std::vector<int> &dumuxFaceIDs)
{
  assert( wasCreated_ );
  assert( numPoints == dumuxFaceIDs.size() );
  meshID_ = precice_->getMeshID(meshName);
  vertexIDs_.resize( numPoints );
  int tmp = numPoints;
  precice_->setMeshVertices( meshID_, tmp, coordinates.data(), vertexIDs_.data() );
  indexMapper_.createMapping( dumuxFaceIDs, vertexIDs_);
  meshWasCreated_ = true;
  return initialize();
}


bool PreciceAdapter::isInitialDataAvailable()
{
  assert( wasCreated_ );
  return precice_->isReadDataAvailable();
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

/*
void PreciceAdapter::writeInitialBlockScalarData( const int dataID,
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

void PreciceAdapter::announceAllInitialDataWritten()
{
  assert( wasCreated_ );
  if ( hasToWriteInitialData() )
  {
    precice_->fulfilledAction( precice::constants::actionWriteInitialData() );
  }
}
*/
PreciceAdapter::~PreciceAdapter()
{
  //finalize();
}

