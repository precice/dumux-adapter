#include "precicewrapper.hh"

#include <cassert>

PreciceWrapper::PreciceWrapper():
  wasCreated_(false), precice_(nullptr), meshWasCreated_(false), preciceWasInitialized_(false),
  meshID_(0), temperatureID_(0), heatFluxID_(0), timeStepSize_(0.)
{

}

PreciceWrapper& PreciceWrapper::get()
{
  static PreciceWrapper instance;
  return instance;
}

void PreciceWrapper::configure( const std::string& configurationFileName )
{
  get().precice_->configure( configurationFileName );
}

void PreciceWrapper::createInstance( const std::string& name, const int rank, const int size )
{
  assert( get().precice_ == nullptr );
  get().precice_ = std::make_unique<precice::SolverInterface>(name, rank, size);
  get().wasCreated_ = true;
}

int PreciceWrapper::getDimensions()
{
  assert( get().wasCreated_ );
  return get().precice_->getDimensions();
}
/*
void PreciceWrapper::setMeshName(const std::string& meshName)
{
  assert( get().wasCreated_ );
  get().meshID_ = get().precice_->getMeshID(meshName);
}
*/
void PreciceWrapper::setMesh(const std::string& meshName,
                             const size_t numPoints,
                              std::vector<double>& coordinates,
                              const std::vector<int>& dumuxFaceIDs )
{
  assert( get().wasCreated_ );
  assert( numPoints == dumuxFaceIDs.size() );
  get().meshID_ = get().precice_->getMeshID(meshName);
  get().vertexIDs_.resize( numPoints );
  get().precice_->setMeshVertices( get().meshID_, numPoints, coordinates.data(), get().vertexIDs_.data() );
  get().indexMapper_.createMapping( dumuxFaceIDs, get().vertexIDs_);
  get().meshWasCreated_ = true;
}
/*
int PreciceWrapper::getDataID( const std::string& dataName, const int meshID )
{
  assert( get().wasCreated_ );
  return get().precice_->getDataID( dataName, meshID );
}
*/
double PreciceWrapper::initialize()
{
  assert( get().wasCreated_ );
  assert( get().meshWasCreated_ );

  get().temperatureID_ = get().precice_->getDataID( "Temperature", get().meshID_ );
  get().heatFluxID_ = get().precice_->getDataID( "Heat-Flux", get().meshID_ );

  get().timeStepSize_ = get().precice_->initialize();
  assert( get().timeStepSize_ > 0 );

  get().precice_->initializeData();
  get().preciceWasInitialized_ = true;
  return get().timeStepSize_;
}

/*
void PreciceWrapper::initializeData()
{
  assert( get().wasCreated_ );
  get().precice_->initializeData();
}
*/

double PreciceWrapper::advance( const double computedTimeStepLength )
{
  assert( get().wasCreated_ );
  return get().precice_->advance( computedTimeStepLength );
}

bool PreciceWrapper::isCouplingOngoing()
{
  assert( get().wasCreated_ );
  return get().precice_->isCouplingOngoing();
}

size_t PreciceWrapper::getNumberOfVertices()
{
  assert( get().wasCreated_ );
  return get().vertexIDs_.size();
}

void PreciceWrapper::readScalarQuantitiy(const int dataID, std::vector<double> &data)
{
  assert( get().wasCreated_ );
  get().precice_->readBlockScalarData( dataID, get().vertexIDs_.size(),
                                       get().vertexIDs_.data(), data.data() );
}

void PreciceWrapper::writeScalarQuantitiy(const int dataID, std::vector<double> &data)
{
  assert( get().wasCreated_ );
  get().precice_->writeBlockScalarData( dataID, get().vertexIDs_.size(),
                                        get().vertexIDs_.data(), data.data() );
}

void PreciceWrapper::print(std::ostream& os)
{
  os << get().indexMapper_;
}

/*
void PreciceWrapper::writeSolidTemperature(std::vector<double> &temperature)
{
  assert( get().wasCreated_ );
  get().precice_->writeBlockScalarData( get().solidTemperatureID_, get().vertexIDs_.size(),
                                       get().vertexIDs_.data(), temperature.data() );

}

void PreciceWrapper::readSolidTemperature(std::vector<double> &temperature)
{
  assert( get().wasCreated_ );
  get().precice_->readBlockScalarData( get().solidTemperatureID_, get().vertexIDs_.size(),
                                       get().vertexIDs_.data(), temperature.data() );
}
*/

/*
void PreciceWrapper::readBlockScalarData( const int dataID,
                                          const int size,
                                          int* const valueIndices,
                                          double* const values )
{
  assert( get().wasCreated_ );
  get().precice_->readBlockScalarData( dataID, size, valueIndices, values );
}

void PreciceWrapper::writeBlockScalarData( const int dataID,
                                           const int size,
                                           int* const valueIndices,
                                           double* const values )
{
  assert( get().wasCreated_ );
  return get().precice_->writeBlockScalarData( dataID, size, valueIndices, values );
}
*/

bool PreciceWrapper::checkIfActionIsRequired( const std::string& condition )
{
  assert( get().wasCreated_ );
  return get().precice_->isActionRequired( condition );
}

void PreciceWrapper::actionIsFulfilled(const std::string& condition)
{
  assert( get().wasCreated_ );
  get().precice_->fulfilledAction( condition );
}

bool PreciceWrapper::hasToWriteInitialData()
{
  assert( get().wasCreated_ );
  return get().checkIfActionIsRequired(precice::constants::actionWriteInitialData());
}


bool PreciceWrapper::hasToReadIterationCheckpoint()
{
  assert( get().wasCreated_ );
  return get().checkIfActionIsRequired(precice::constants::actionReadIterationCheckpoint());
}

void PreciceWrapper::announceIterationCheckpointRead()
{
  assert( get().wasCreated_ );
  get().actionIsFulfilled( precice::constants::actionWriteIterationCheckpoint() );
}

bool PreciceWrapper::hasToWriteIterationCheckpoint()
{
  assert( get().wasCreated_ );
  return get().checkIfActionIsRequired(precice::constants::actionWriteIterationCheckpoint());
}

void PreciceWrapper::announceIterationCheckpointWritten()
{
  assert( get().wasCreated_ );
  get().actionIsFulfilled( precice::constants::actionWriteIterationCheckpoint() );
}

/*
void PreciceWrapper::writeInitialBlockScalarData( const int dataID,
                                                  const int size,
                                                  int* const valueIndices,
                                                  double* const values )
{
  assert( get().wasCreated_ );
  if ( get().hasToWriteInitialData() )
  {
    get().precice_->writeBlockScalarData( dataID, size, valueIndices, values );
  }
}

void PreciceWrapper::announceAllInitialDataWritten()
{
  assert( get().wasCreated_ );
  if ( get().hasToWriteInitialData() )
  {
    get().precice_->fulfilledAction( precice::constants::actionWriteInitialData() );
  }
}
*/
PreciceWrapper::~PreciceWrapper()
{
  if ( !get().wasCreated_ )
    get().initialize();
  get().precice_->finalize();
}

