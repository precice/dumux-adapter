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

  static bool checkIfActionIsRequired( const std::string& condition );
  static void actionIsFulfilled( const std::string& condition );
  static bool hasToWriteInitialData();

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

  static PreciceWrapper& get();

  static void createInstance( const std::string& name, const int rank, const int size );
  static void configure( const std::string& configurationFileName );

  static int getDimensions();
  // static int getMeshID( const std::string& meshName );
  //static void setMeshName( const std::string& meshName );
  //static int getDataID( const std::string& dataName, const int meshID );

  //static void writeInitialBlockScalarData( const int dataID,
  //                                         const int size,
  //                                         int* const valueIndices,
  //                                         double* const values );
  //static void announceAllInitialDataWritten();

  static bool hasToReadIterationCheckpoint();
  static void announceIterationCheckpointRead();
  static bool hasToWriteIterationCheckpoint();
  static void announceIterationCheckpointWritten();

  static void setMesh( const std::string& meshName,
                       const size_t numPoints,
                       std::vector<double>& coordinates,
                       const std::vector<int>& dumuxFaceIDs ) ;

  static double initialize();
  static void finalize();
  //static void initializeData();

  static double advance( const double computedTimeStepLength );
  static bool isCouplingOngoing();

  static size_t getNumberOfVertices();


  static void readScalarQuantitiy( const int dataID, std::vector<double>& data );
  static void writeScalarQuantitiy( const int dataID, std::vector<double>& data );


  static void print( std::ostream& os );

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
