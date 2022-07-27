#ifndef PRECICEWRAPPER_HH
#define PRECICEWRAPPER_HH

#include <ostream>
#include <precice/SolverInterface.hpp>
#include <string>

#include "dumuxpreciceindexmapper.hh"

/*!
 * @brief Namespace of dumux-precice
 *
 */
namespace Dumux::Precice
{
enum class QuantityType { Scalar, Vector };

/*!
 * @brief A DuMuX-preCICE coupling adapter class
 *
 * The class provides an interface to DuMuX to couple simulations
 * via the coupling tool preCICE. The class aims to provide an
 * easy-to-use interface that is reasonably close to the coupling
 * interface for monolithic couplings that is integrated into DuMuX.
 *
 * \note The coupling adapter is currently implemented as a Singleton.
 *
 */
class CouplingAdapter
{
private:
    //! True if preCICE instance was initiated
    bool wasCreated_;
    //! Pointer to preCICE instance
    std::unique_ptr<precice::SolverInterface> precice_;
    //! Constructor
    CouplingAdapter();
    /*!
     * @brief Checks whether an action predefined by preCICE
     *        needs to be carried out.
     *
     * @param[in] condition Name of the action.
     * @return true Action must be carried out.
     * @return false Action must not be carried out.
     */
    bool checkIfActionIsRequired(const std::string &condition);
    /*!
     * @brief Announce to preCICE that an action was carried out.
     *
     * @param[in] condition Name of the action.
     */
    void actionIsFulfilled(const std::string &condition);
    /*!
     * @brief Reads full block of data from preCICE.
     *
     * @param[in] dataID Identifier of dataset to read.
     * @param[out] data Vector to store the read data to.
     */
    void readBlockDataFromPrecice(const int dataID,
                                  std::vector<double> &data,
                                  const QuantityType quantity_type);
    /*!
     * @brief Writes full block of data to preCICE.
     *
     * @param[in] dataID Identifier of dataset to read.
     * @param[in] data Vector containing data to write into preCICE's buffer.
     */
    void writeBlockDataToPrecice(const int dataID,
                                 std::vector<double> &data,
                                 const QuantityType quantity_type);
    /*!
     * @brief Gives the number of quantities/datasets defined on coupling interface.
     *
     * @return size_t Number of quantities defined on the coupling interface.
     */
    size_t numberOfQuantities() const { return dataNames_.size(); }
    //! True if the coupling mesh was created.
    bool meshWasCreated_;
    //! True if precice::SolverInterface.initialize() has been called.
    bool preciceWasInitialized_;
    //! True if instance owns an instance of DumuxPreciceIndexMapper.
    bool hasIndexMapper_;
    //! Stores identifier of the coupling mesh provided by preCICE.
    int meshID_;
    //! Time step size.
    double timeStepSize_;
    //! Vector of names of data exchanged over coupling interface.
    std::vector<std::string> dataNames_;
    //! Vector of identifiers of data exchanged over coupling interface.
    std::vector<int> preciceDataID_;
    //! Vector storing data vectors of the data exchanged over the coupling interface.
    std::vector<std::vector<double> > dataVectors_;
    //! Vector of identifiers of the vertices of the coupling mesh.
    std::vector<int> vertexIDs_;  //should be size_t
    /*!
     * @brief Instance of DumuxPreciceIndexMapper that translates between
     *        DuMuX' identifiers of vertices and preCICE's identifiers.
     *
     */
    Internal::DumuxPreciceIndexMapper<int> indexMapper_;
    /*!
     * @brief Get the of quantities exchanged.
     *
     * @return size_t Number of quantities defined on coupling interface.
     */
    size_t getNumberOfQuantities() const { return dataNames_.size(); }
    //! Number of expected quantities on the coupling interface.
    static constexpr size_t reserveSize_ = 4;
    /*!
     * @brief Destroy the CouplingAdapter object
     *
     */
    ~CouplingAdapter();

public:
    CouplingAdapter(const CouplingAdapter &) = delete;
    void operator=(const CouplingAdapter &) = delete;

    /*!
     * @brief Get the instance of the CouplingAdapter
     *
     * @return CouplingAdapter& Reference to current instance of the CouplingAdapter
     */
    static CouplingAdapter &getInstance();
    /*!
     * @brief Announces the DuMuX solver.
     *
     * @param[in] name Name of the DuMuX solver.
     * @param[in] configurationFileName  Path and file name to preCICE configuration file.
     * @param[in] rank Rank of the current process of the DuMuX solver.
     * @param[in] size Total number of processes of the DuMuX solver.
     */
    void announceSolver(const std::string &name,
                        const std::string configurationFileName,
                        const int rank,
                        const int size);
    /*!
     * @brief Announces an additional quantity on the coupling interface.
     *
     * Internally, the quantity is announced to preCICE and the corresponding
     * data structures are initilized to store information about the quantity.
     *
     * @param[in] name Name of the scalar quantity.
     * @param[in] quantity_type Type (Scalar or Vector) of the quantity
     * @return size_t Number of currently announced quantities.
     */
    size_t announceQuantity(const std::string &name,
                            const QuantityType quantity_type);

    /*!
     * @brief Announces an additional scalar quantity on the coupling interface.
     *
     * Internally, the scalar quantity is announced to preCICE and the corresponding
     * data structures are initilized to store information about the quantity.
     *
     * @param[in] name Name of the scalar quantity.
     * @return size_t Number of currently announced quantities.
     */
    size_t announceScalarQuantity(const std::string &name);
    /*!
     * @brief Announces an additional vector quantity on the coupling interface.
     *
     * Internally, the vector quantity is announced to preCICE and the corresponding
     * data structures are initilized to store information about the quantity.
     *
     * @param[in] name Name of the vector quantity.
     * @return size_t Number of currently announced quantities.
     */
    size_t announceVectorQuantity(const std::string &name);
    /*!
     * @brief Get the number of spatial dimensions
     *
     * @return int Number of space dimensions. Legal values are 2 and 3.
     */
    int getDimensions() const;
    /*!
     * @brief Checks if simulation checkpoint needs to be restored.
     *
     * @return true Simulation checkpoint has to be restored.
     * @return false No further action is needed.
     */
    bool hasToReadIterationCheckpoint();
    //! Announce that the simulation checkpoint was read.
    void announceIterationCheckpointRead();
    /*!
     * @brief Checks if simulation checkpoint needs to be saved.
     *
     * @return true Simulation checkpoints needs to be stored.
     * @return false No further action is needed.
     */
    bool hasToWriteIterationCheckpoint();
    //! Announce that the simulation checkpoint was written.
    void announceIterationCheckpointWritten();
    /*!
     * @brief Checks if initial coupling data has to be wrutten.
     *
     * @return true Initial coupling data has to be provided.
     * @return false No further action is needed.
     */
    bool hasToWriteInitialData();
    //! Announce that the initial coupling data has been written.
    void announceInitialDataWritten();
    /*!
     * @brief Adds mesh for coupling of solvers.
     *
     * @param[in] meshName Name of the mesh.
     * @param[in] numPoints Number of points/vertices.
     * @param[in] coordinates Coordinates of the points.
     *
     * \note The coordinates need to be stored consecutively
     *       according to their spatial coordinates as.\n
     *       Example 2D:\n
     *       [x_1, y_1, x_2, y_2,...x_numPoints, y_numPoints]\n
     *       Example 3D:\n
     *       [x_1, y_1, z_1, x_2, y_2, z_2,...x_numPoints, y_numPoints, z_numPoints]
     */
    void setMesh(const std::string &meshName,
                 const size_t numPoints,
                 std::vector<double> &coordinates);
    /*!
     * @brief Initializes the coupling
     *
     * The coupling needs be initialized after all quantities/datasets and coupling meshes
     * are known.
     *
     * @return double Maximum allowed time step size.
     */
    double initialize();
    /*!
     * @brief Creates mapping between DuMuX' face identifiers and preCICE's
     *        vertex identifiers.
     *
     * @param[in] dumuxFaceIDs Vector containing the face identifiers on the coupling interface.
     *
     * \note The order of the face identifiers must be correspond to the order of coordinates
     *       passed in setMesh.
     */
    void createIndexMapping(const std::vector<int> &dumuxFaceIDs);
    /*!
     * @brief Sets the coupling mesh and initializes coupling.
     *
     * This is a convenience function that sets the coupling mesh using setMesh.
     * Afterwards, the coupling is initialized via initialzie.
     *
     * @param[in] meshName Name of the mesh.
     * @param[in] numPoints Number of points/vertices.
     * @param[in] coordinates Coordinates of the points.
     * @return double Maximum allowed time step size.
     */
    double setMeshAndInitialize(const std::string &meshName,
                                const size_t numPoints,
                                std::vector<double> &coordinates);

    /*!
     * @brief Initializes the coupling data.
     *
     * If one wants to set non-zero data, one has to write data to the
     * corresponding quantities via one of the `write` functions first.
     *
     */
    void initializeData();
    /*!
     * @brief Destroys the coupling.
     *
     * This function must called at the end of a simulation.
     *
     */
    void finalize();
    /*!
     * @brief Advances coupling by the given time step length.
     *
     * @param[in] computedTimeStepLength Time step lengths of the current simulation stel.
     * @return double Maximum time step length for successive time steps.
     */
    double advance(const double computedTimeStepLength);
    /*!
     * @brief Checks whether the coupling is still ongoing.
     *
     * @return true Coupling is still ongoing.
     * @return false Coupling finished.
     */
    bool isCouplingOngoing();
    /*!
     * @brief Get the number of vertices on the coupling interface.
     *
     * @return size_t Number of vertices on the coupling interface.
     */
    size_t getNumberOfVertices();
    /*!
     * @brief Gets value of a scalar quantity.
     *
     * @param[in] dataID Identifier of the quantity.
     * @param[in] faceID Identifier of the face according to DuMuX' numbering.
     * @return double Value of scalar quantity.
     */
    double getScalarQuantityOnFace(const size_t dataID, const int faceID) const;
    // /*!
    //  * @brief Gets value of a vector quantity.
    //  *
    //  * @param[in] dataID Identifier of the quantity.
    //  * @param[in] faceID Identifier of the face according to DuMuX' numbering.
    //  * @return std::vector<double> Value of vector quantity.
    //  */
    // std::vector<double> getVectorQuantityOnFace(const size_t dataID, const int faceID) const;
    // std::vector<double> getVectorQuantity(const size_t dataID) const;
    /*!
     * @brief Gets value of a vector quantity.
     *
     * @param[in] dataID Identifier of the quantity.
     * @param[in] faceID Identifier of the face according to DuMuX' numbering.
     * @return const std::vector<double>& Value of vector quantity.
     */
    const std::vector<double> &getVectorScalarQuantityOnFace(
        const size_t dataID,
        const int faceID) const;
    /*!
     * @brief Writes value of scalar quantity on given face.
     *
     * @param[in] dataID Identifier of the quantity.
     * @param[in] faceID Identifier of the face according to DuMuX' numbering.
     * @param[in] value  Value of scalar quantity.
     */
    void writeScalarQuantityOnFace(const size_t dataID,
                                   const int faceID,
                                   const double value);

    /*!
     * @brief Returns reference to data vector of quantity with given identifier.
     *
     * @param dataID Identifier of the quantity.
     * @return[in] std::vector<double>& Reference to data vector.
     */
    std::vector<double> &getQuantityVector(const size_t dataID);
    /*!
     * @brief Returns const reference to data vector of quantity with given identifier.
     *
     * @param[in] dataID Identifier of the quantity.
     * @return std::vector<double>& Const reference to data vector.
     */
    const std::vector<double> &getQuantityVector(const size_t dataID) const;

    /*!
     * @brief Writes value of scalar or vector quantity on all vertices.
     *
     * @param[in] dataID Identifier of the quantity.
     * @param[in] values Value of the scalar or vector quantity.
     */
    void writeQuantityVector(const size_t dataID, std::vector<double> &values);
    /*!
     * @brief Writes data from adapter's buffer into preCICE's communication buffer.
     *
     * @param[in] dataID Identifier of the quantity to write into communication buffer.
     */
    void writeQuantityToOtherSolver(const size_t dataID,
                                    const QuantityType quantity_type);
    /*!
     * @brief Reads data from preCICE's communication buffer and puts it into adapter's buffer.
     *
     * @param dataID Identifier of the quantity to read into adapter buffer.
     */
    void readQuantityFromOtherSolver(const size_t dataID,
                                     const QuantityType quantity_type);
    /*!
     * @brief Writes data from adapter's buffer into preCICE's communication buffer.
     *
     * @param[in] dataID Identifier of the quantity to write into communication buffer.
     */
    void writeScalarQuantityToOtherSolver(const size_t dataID);
    /*!
     * @brief Reads data from preCICE's communication buffer and puts it into adapter's buffer.
     *
     * @param dataID Identifier of the quantity to read into adapter buffer.
     */
    void readScalarQuantityFromOtherSolver(const size_t dataID);
    /*!
     * @brief Writes data from adapter's buffer into preCICE's communication buffer.
     *
     * @param[in] dataID Identifier of the quantity to write into communication buffer.
     */
    /*!
     * @brief Checks whether face with given identifier is part of coupling interface.
     *
     * @param[in] faceID Identifier of the face according to DuMuX' numbering.
     * @return true Face is part of coupling interface.
     * @return false Face is not part of coupling interface.
     */
    bool isCoupledEntity(const int faceID) const;
    /*!
     * @brief Get a quantity's numeric identifier from its name.
     *
     * @param[in] dataName Name of the quantity.
     * @return size_t Numeric identifier of quantity.
     */
    size_t getIdFromName(const std::string &dataName) const;
    /*!
     * @brief Get a quantitiy's name from its numeric identifier.
     *
     * @param dataID Identifier of the quantity to read into adapter buffer.
     * @return std::string Name of the quantity.
     */
    std::string getNameFromId(const size_t dataID) const;
    /*!
     * @brief Prints status of coupling adapter to given output stream.
     *
     * @param os Output stream.
     */
    void print(std::ostream &os);
};

}  // namespace Dumux::Precice
#endif
