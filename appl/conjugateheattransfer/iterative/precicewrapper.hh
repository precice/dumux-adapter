#ifndef PRECICEWRAPPER_HH
#define PRECICEWRAPPER_HH

#include<string>
#include<ostream>
#include<precice/SolverInterface.hpp>

#include "dumuxpreciceindexwrapper.hh"


class PreciceWrapper
{

private:
  bool wasCreated_;
  std::unique_ptr<precice::SolverInterface> precice_;

  PreciceWrapper();

  bool checkIfActionIsRequired( const std::string& condition );
  void actionIsFulfilled( const std::string& condition );
  bool hasToWriteInitialData();

  bool meshWasCreated_;
  bool preciceWasInitialized_;
  int meshID_;
  int dimension_;
  int temperatureID_;
  int heatFluxID_;

  double timeStepSize_;

  std::vector<int> vertexIDs_; //should be size_t
  std::vector<double> temperature_;
  std::vector<double> heatFlux_;

  DumuxPreciceIndexMapper<int> indexMapper_;


  ~PreciceWrapper();
public:
  PreciceWrapper(const PreciceWrapper&) = delete;
  void operator=(const PreciceWrapper&) = delete;

  static PreciceWrapper& getInstance();

  void announceSolver( const std::string& name, const int rank, const int size );
  void configure( const std::string& configurationFileName );

  int getDimensions();
  // static int getMeshID( const std::string& meshName );
  //static void setMeshName( const std::string& meshName );
  //static int getDataID( const std::string& dataName, const int meshID );

  //static void writeInitialBlockScalarData( const int dataID,
  //                                         const int size,
  //                                         int* const valueIndices,
  //                                         double* const values );
  //static void announceAllInitialDataWritten();

  bool hasToReadIterationCheckpoint();
  void announceIterationCheckpointRead();
  bool hasToWriteIterationCheckpoint();
  void announceIterationCheckpointWritten();

  void setMesh( const std::string& meshName,
                const size_t numPoints,
                std::vector<double>& coordinates,
                const std::vector<int>& dumuxFaceIDs ) ;

  double initialize();
  void finalize();
  //static void initializeData();

  double advance( const double computedTimeStepLength );
  bool isCouplingOngoing();

  size_t getNumberOfVertices();


  double getHeatFluxAtFace( const int faceID ) const;
  double getTemperatureAtFace( const int faceID ) const;

  void writeHeatFluxOnFace( const int faceID );
  void writeTemperatureOnFace( const int faceID );

  bool isCoupledEntity( const int faceID ) const;


//  static void readScalarQuantitiy( const int dataID, std::vector<double>& data );
//  static void writeScalarQuantitiy( const int dataID, std::vector<double>& data );


  void print( std::ostream& os );

  /*
  static void writeSolidTemperature( std::vector<double>& temperature );
  static void readSolidTemperature( std::vector<double>& temperature );
  static void writeFluidTemperature( const std::vector<double>& temperature );
  static void readFluidTemperature( std::vector<double>& temperature );
  */

  /*
  static void readBlockScalarData( const int dataID,
                                   const int size,
                                   int* const valueIndices,
                                   double* const values );

  static void writeBlockScalarData( const int dataID,
                                    const int size,
                                    int* const valueIndices,
                                    double* const values );
  */

};


#endif
