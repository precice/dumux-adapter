#ifndef PRECICEWRAPPER_HH
#define PRECICEWRAPPER_HH

#include<string>
#include<ostream>
#include<precice/SolverInterface.hpp>

#include "dumuxpreciceindexwrapper.hh"

namespace precice_wrapper{

  enum HeatFluxType
  {
    UNDEFINED, FreeFlow, Solid
  };

class PreciceWrapper
{

private:
  bool wasCreated_;
  std::unique_ptr<precice::SolverInterface> precice_;

  PreciceWrapper();

  bool checkIfActionIsRequired( const std::string& condition );
  void actionIsFulfilled( const std::string& condition );

  void readBlockScalarDataFromPrecice( const int dataID, std::vector<double>& data );
  void writeBlockScalarDataToPrecice( const int dataID, std::vector<double>& data );

  bool hasToWriteInitialData();

  bool meshWasCreated_;
  bool preciceWasInitialized_;
  int meshID_;
  int dimension_;
  int freeFlowHeatFluxID_;
  int solidHeatFluxID_;

  double timeStepSize_;

  HeatFluxType writeHeatFluxType_;
  HeatFluxType readHeatFluxType_;


  std::vector<int> vertexIDs_; //should be size_t
  std::vector<double> freeFlowHeatFlux_;
  std::vector<double> solidHeatFlux_;

  DumuxPreciceIndexMapper<int> indexMapper_;


  ~PreciceWrapper();
public:
  PreciceWrapper(const PreciceWrapper&) = delete;
  void operator=(const PreciceWrapper&) = delete;

  static PreciceWrapper& getInstance();

  void announceSolver( const std::string& name, const int rank, const int size );
  void configure( const std::string& configurationFileName );

  void announceHeatFluxToWrite( const HeatFluxType heatFluxType );
  void announceHeatFluxToRead( const HeatFluxType heatFluxType );

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
  void writeHeatFluxOnFace( const int faceID, const double value );

  void writeHeatFluxToOtherSolver();
  void readHeatFluxFromOtherSolver();


  bool isCoupledEntity( const int faceID ) const;


  std::vector<double>& getHeatFluxToWrite();

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

}
#endif
