#ifndef PRECICEWRAPPER_HH
#define PRECICEWRAPPER_HH

#include<string>
#include<ostream>
#include<precice/SolverInterface.hpp>

#include <dune/common/fvector.hh>

#include "../src/dumuxpreciceindexwrapper.hh"


namespace precice_adapter{

  template<unsigned DIM>
  using VectorXd = Dune::DenseVector<Dune::FieldVector<double, DIM>>;

  struct DataContainer {
      const std::string name;
      const int preciceDataID;
      std::vector<double> data;
      const bool isVectorData;

      DataContainer( const std::string& name,
                     const int dataID,
                     const size_t dataSize,
                     const bool isVectorData ) : name(name), preciceDataID( dataID ), isVectorData( isVectorData )
      {
        data.resize( dataSize );
      }

  };

  class PreciceAdapter
  {

    private:
      bool wasCreated_;
      std::unique_ptr<precice::SolverInterface> precice_;

      PreciceAdapter();

      bool checkIfActionIsRequired( const std::string& condition );
      void actionIsFulfilled( const std::string& condition );

      void readBlockScalarDataFromPrecice( const int dataID, std::vector<double>& data );
      void writeBlockScalarDataToPrecice( const int dataID, std::vector<double>& data );

      bool meshWasCreated_;
      bool preciceWasInitialized_;
      bool hasIndexMapper_;
      int meshID_;

      double timeStepSize_;

      //      std::vector< int > preciceDataID_;
//      std::vector< std::vector< double > > dataVectors_;

      std::vector< DataContainer > interfaceData_;

      std::vector<int> vertexIDs_; //should be size_t

      DumuxPreciceIndexMapper<int> indexMapper_;

      size_t getNumberOfQuantities() const { return interfaceData_.size(); }

      auto getDataIterator( const std::string& name ) const {
        auto dataName = [name](const DataContainer& item) { return item.name == name; };
        auto it = std::find_if( interfaceData_.begin(), interfaceData_.end(), dataName );
        return it;
      }

      static constexpr size_t reserveSize_ = 4;

      template<unsigned DIM>
      size_t announceQuantity( const std::string& name );

      ~PreciceAdapter();
    public:
      PreciceAdapter(const PreciceAdapter&) = delete;
      void operator=(const PreciceAdapter&) = delete;

      static PreciceAdapter& getInstance();

      void announceSolver( const std::string& name,
                           const std::string configurationFileName,
                           const int rank,
                           const int size );

      [[deprecated("Please use annoucScalarQuantity or annountVectorQuantity<DIM>!")]]
      size_t announceQuantity( const std::string& name );

      size_t announceScalarQuantity( const std::string& name ) {
        return announceQuantity<1>( name );
      }

      template<unsigned DIM>
      size_t announceVectorQuantity( const std::string& name ) {
        return announceQuantity<DIM>( name );
      }

      int getDimensions() const;

      bool hasToReadIterationCheckpoint();
      void announceIterationCheckpointRead();
      bool hasToWriteIterationCheckpoint();
      void announceIterationCheckpointWritten();

      bool hasToWriteInitialData();
      void announceInitialDataWritten();

      void setMesh( const std::string& meshName,
                    const size_t numPoints,
                    std::vector<double>& coordinates );

      double initialize();


      void createIndexMapping( const std::vector<int>& dumuxFaceIDs );

      double setMeshAndInitialize( const std::string& meshName,
                                   const size_t numPoints,
                                   std::vector<double>& coordinates) ;

      void initializeData();
      void finalize();

      double advance( const double computedTimeStepLength );
      bool isCouplingOngoing();

      size_t getNumberOfVertices();


      double getScalarQuantityOnFace( const size_t dataID, const int faceID ) const;

      template<unsigned DIM>
      const VectorXd<DIM>& getVectorScalarQuantityOnFace( const size_t dataID, const int faceID ) const;

      void writeScalarQuantityOnFace( const size_t dataID,
                                      const int faceID,
                                      const double value );

      template<unsigned DIM>
      void writeVectorQuantityOnFace( const size_t dataID,
                                      const int faceID,
                                      const VectorXd<DIM>& value );

      std::vector<double>& getQuantityVector( const size_t dataID );

      const std::vector<double>& getQuantityVector( const size_t dataID ) const;

      void writeScalarQuantityVector( const size_t dataID,
                                      std::vector<double>& values );

      void writeScalarQuantityToOtherSolver( const size_t dataID );
      void readScalarQuantityFromOtherSolver( const size_t dataID );

      bool isCoupledEntity( const int faceID ) const;

      size_t getIdFromName( const std::string& dataName) const;

      std::string getNameFromId( const size_t dataID ) const;

      void print( std::ostream& os );
  };


  template<unsigned DIM>
  size_t PreciceAdapter::announceQuantity( const std::string& name ) {
    assert( meshWasCreated_ );
    // Check if data is already announced
    auto it = getDataIterator( name );
    if ( it != interfaceData_.end() )
    {
      throw( std::runtime_error(" Error! Duplicate quantity announced! ") );
    }

    interfaceData_.push_back( DataContainer(name,
                                            precice_->getDataID( name, meshID_ ),
                                            vertexIDs_.size() * DIM,
                                            DIM == getDimensions() )
                              );


    return getNumberOfQuantities()-1;
  }

  template<unsigned DIM>
  const VectorXd<DIM>& PreciceAdapter::getVectorScalarQuantityOnFace( const size_t dataID,
                                                                      const int faceID ) const {
    //TODO
    assert( DIM == getDimensions() );
    const auto idx = indexMapper_.getPreciceId( faceID );
    assert( dataID < interfaceData_.size() );
    assert( interfaceData_[dataID].isVectorData );
    const std::vector<double>& quantityVector = interfaceData_[ dataID ].data;
    assert(idx < quantityVector.size() );

    VectorXd<DIM> vector;
    assert( idx*DIM+DIM-1 < quantityVector.size() );
    std::copy( quantityVector[idx*DIM], quantityVector[idx*DIM+DIM-1], vector );
    return vector;
  }

  template<unsigned DIM>
  void PreciceAdapter::writeVectorQuantityOnFace( const size_t dataID,
                                                  const int faceID,
                                                  const VectorXd<DIM>& value ) {
    assert( wasCreated_ );
    assert( hasIndexMapper_ );
    if ( !hasIndexMapper_ )
    {
      throw std::runtime_error("Writing quantity using faceID, but index mapping was not created!");
    }

    const auto idx = indexMapper_.getPreciceId( faceID );
    assert( dataID < interfaceData_.size() );
    assert( interfaceData_[dataID].isVectorData );
    std::vector<double>& quantityVector = interfaceData_[ dataID ].data;
    assert( idx < quantityVector.size() );

    assert( idx*DIM + DIM - 1 < quantityVector.size() );
    std::copy( value.begin(), value.end(), quantityVector[idx*DIM] );
  }


}
#endif
