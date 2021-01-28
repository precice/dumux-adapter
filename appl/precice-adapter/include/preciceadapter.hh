#ifndef PRECICEWRAPPER_HH
#define PRECICEWRAPPER_HH

#include <ostream>
#include <precice/SolverInterface.hpp>
#include <string>

#include "../src/dumuxpreciceindexwrapper.hh"

namespace precice_adapter
{
class PreciceAdapter
{
   private:
    bool wasCreated_;
    std::unique_ptr<precice::SolverInterface> precice_;

    PreciceAdapter();

    bool checkIfActionIsRequired(const std::string &condition);
    void actionIsFulfilled(const std::string &condition);

    void readBlockScalarDataFromPrecice(const int dataID,
                                        std::vector<double> &data);
    void writeBlockScalarDataToPrecice(const int dataID,
                                       std::vector<double> &data);

    size_t numberOfQuantities() const { return dataNames_.size(); }

    bool meshWasCreated_;
    bool preciceWasInitialized_;
    bool hasIndexMapper_;
    int meshID_;

    double timeStepSize_;

    std::vector<std::string> dataNames_;
    std::vector<int> preciceDataID_;
    std::vector<std::vector<double> > dataVectors_;

    std::vector<int> vertexIDs_;  //should be size_t

    DumuxPreciceIndexMapper<int> indexMapper_;

    size_t getNumberOfQuantities() const { return dataNames_.size(); }

    static constexpr size_t reserveSize_ = 4;

    ~PreciceAdapter();

   public:
    PreciceAdapter(const PreciceAdapter &) = delete;
    void operator=(const PreciceAdapter &) = delete;

    static PreciceAdapter &getInstance();

    void announceSolver(const std::string &name,
                        const std::string configurationFileName,
                        const int rank,
                        const int size);

    size_t announceQuantity(const std::string &name);

    int getDimensions() const;

    bool hasToReadIterationCheckpoint();
    void announceIterationCheckpointRead();
    bool hasToWriteIterationCheckpoint();
    void announceIterationCheckpointWritten();

    bool hasToWriteInitialData();
    void announceInitialDataWritten();

    void setMesh(const std::string &meshName,
                 const size_t numPoints,
                 std::vector<double> &coordinates);

    double initialize();

    void createIndexMapping(const std::vector<int> &dumuxFaceIDs);

    double setMeshAndInitialize(const std::string &meshName,
                                const size_t numPoints,
                                std::vector<double> &coordinates);

    void initializeData();
    void finalize();

    double advance(const double computedTimeStepLength);
    bool isCouplingOngoing();

    size_t getNumberOfVertices();

    double getScalarQuantityOnFace(const size_t dataID, const int faceID) const;

    const std::vector<double> &getVectorScalarQuantityOnFace(
        const size_t dataID,
        const int faceID) const;

    void writeScalarQuantityOnFace(const size_t dataID,
                                   const int faceID,
                                   const double value);

    //  void writeVectorQuantityOnFace( const size_t dataID,
    //                                  const int faceID,
    //                                  const double* value,
    //                                  const size_t size );

    std::vector<double> &getQuantityVector(const size_t dataID);

    const std::vector<double> &getQuantityVector(const size_t dataID) const;

    void writeScalarQuantityVector(const size_t dataID,
                                   std::vector<double> &values);

    void writeScalarQuantityToOtherSolver(const size_t dataID);
    void readScalarQuantityFromOtherSolver(const size_t dataID);

    bool isCoupledEntity(const int faceID) const;

    size_t getIdFromName(const std::string &dataName) const;

    std::string getNameFromId(const size_t dataID) const;

    void print(std::ostream &os);
};

}  // namespace precice_adapter
#endif
