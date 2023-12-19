#ifndef PRECICEWRAPPER_HH
#define PRECICEWRAPPER_HH

#include <ostream>
#include <precice/precice.hpp>
#include <string>

#include "dumuxpreciceindexmapper.hh"

/*!
 * @brief Namespace of dumux-precice
 *
 */
namespace Dumux::Precice
{
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
    std::unique_ptr<precice::Participant> precice_;
    //! True if the coupling mesh was created.
    bool meshWasCreated_;
    //! True if precice::Participant.initialize() has been called.
    bool preciceWasInitialized_;
    //! True if instance owns an instance of DumuxPreciceIndexMapper.
    bool hasIndexMapper_;
    //! Time step size.
    double timeStepSize_;
    //! Map storing meshName:dataName and data vectors
    std::map<std::string, std::vector<double>> dataMap_;
    //! Vector of identifiers (in preCICE) of the vertices of the coupling mesh.
    std::vector<int> vertexIDs_;  //should be size_t
    //! Span of the precice vertex indices vector vertexIDs_
    precice::span<precice::VertexID> vertexIDsSpan_;
    //! Constructor
    CouplingAdapter();
    /*!
     * @brief Instance of DumuxPreciceIndexMapper that translates between
     *        DuMuX' identifiers of vertices and preCICE's identifiers.
     *
     */
    Internal::DumuxPreciceIndexMapper<int> indexMapper_;
    /*!
     * @brief Get the number of quantities exchanged.
     *
     * @return size_t Number of quantities defined on coupling interface.
     */
    size_t getNumberOfQuantities() const { return dataMap_.size(); }
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
                        const std::string &configurationFileName,
                        const int rank,
                        const int size);
    /*!
     * @brief Announces a quantity on the coupling interface.
     *
     * Internally, the quantity is announced to preCICE and the corresponding
     * data structures are initilized to store information about the quantity.
     *
     * @param[in] meshName Name of the mesh.
     * @param[in] dataName Name of the data.
     */
    void announceQuantity(const std::string &meshName,
                          const std::string &dataName);
    /*!
     * @brief Get the number of spatial dimensions
     *
     * @param[in] meshName Name of the mesh
     * @return int Number of space dimensions.
     */
    int getMeshDimensions(const std::string &meshName) const;
    /*!
     * @brief Get the maximum time step size from preCICE
     *
     * @return double time step size
     */
    double getMaxTimeStepSize() const;
    /*!
     * @brief Checks if the participant is required to read an iteration checkpoint. If true, the participant is required to read an iteration checkpoint before calling advance(). 
     *
     * @return true Simulation checkpoint has to be restored.
     * @return false No further action is needed.
     */
    bool requiresToReadCheckpoint();

    /*!
     * @brief Checks if the participant is required to write an iteration checkpoint. If true, the participant is required to write an iteration checkpoint before calling advance(). 
     *
     * @return true Simulation checkpoints needs to be stored.
     * @return false No further action is needed.
     */
    bool requiresToWriteCheckpoint();

    /*!
     * @brief Checks if the participant is required to provide initial data. If true, the participant needs to write initial data to defined vertices prior to calling initialize().
     *
     * @return true Initial coupling data has to be provided.
     * @return false No further action is needed.
     */
    bool requiresToWriteInitialData();

    /*!
     * @brief Adds mesh for coupling of solvers.
     *
     * @param[in] meshName The name of the mesh to add the vertices to.
     * @param[in] positions A span to the coordinates of the vertices.
     * 
     * \note The coordinates need to be stored consecutively
     *       according to their spatial coordinates as.\n
     *       Example 2D:\n
     *       [x_1, y_1, x_2, y_2,...x_numPoints, y_numPoints]\n
     *       Example 3D:\n
     *       [x_1, y_1, z_1, x_2, y_2, z_2,...x_numPoints, y_numPoints, z_numPoints]
     */
    void setMesh(const std::string &meshName, std::vector<double> &positions);
    /*!
     * @brief Initializes the coupling
     *
     * The coupling needs be initialized after all quantities/datasets and coupling meshes
     * are known.
     *
     */
    void initialize();
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
    void advance(const double computedTimeStepLength);
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
     * @brief Reads full block of data from preCICE.
     *
     * @param[in] meshName Name of the mesh.
     * @param[in] dataName Name of the data.
     * @param[in] relativeReadTime The relative time tagged to the data to be read.
     */
    void readQuantityFromOtherSolver(const std::string &meshName,
                                     const std::string &dataName,
                                     double relativeReadTime);
    /*!
     * @brief Writes full block of data to preCICE.
     *
     * @param[in] meshName Name of the mesh.
     * @param[in] dataName Name of the data.
     */
    void writeQuantityToOtherSolver(const std::string &meshName,
                                    const std::string &dataName);
    /*!
     * @brief Gets value of a scalar quantity on a finite volume face.
     *
     * @param[in] meshName Name of the mesh.
     * @param[in] dataName Name of the data.
     * @param[in] faceID Identifier of the face according to DuMuX' numbering.
     * @return double Value of scalar quantity.
     */
    double getScalarQuantityOnFace(const std::string &meshName,
                                   const std::string &dataName,
                                   const int faceID);
    /*!
     * @brief Writes value of scalar quantity on a given finite volume face to data map.
     *
     * @param[in] meshName Name of the mesh.
     * @param[in] dataName Name of the data.
     * @param[in] faceID Identifier of the face according to DuMuX' numbering.
     * @param[in] value  Value of scalar quantity.
     */
    void writeScalarQuantityOnFace(const std::string &meshName,
                                   const std::string &dataName,
                                   const int faceID,
                                   const double value);
    /*!
     * @brief Gets the quantity value vector from the data map according to the mesh and data name.
     *
     * @param[in] meshName Name of the mesh.
     * @param[in] dataName Name of the data.
     * @return The value vector of the quantity.
     */
    std::vector<double> &getQuantityVector(const std::string &meshName,
                                           const std::string &dataName);
    /*!
     * @brief Writes the quantity value vector into the data map.
     *
     * @param[in] meshName Name of the mesh.
     * @param[in] dataName Name of the data.
     * @param[in] values Value of the scalar or vector quantity.
     */
    void writeQuantityVector(const std::string &meshName,
                             const std::string &dataName,
                             std::vector<double> &values);
    /*! 
     * @brief Checks whether face with given identifier is part of coupling interface.
     *
     * @param[in] faceID Identifier of the face according to DuMuX' numbering.
     * @return true Face is part of coupling interface.
     * @return false Face is not part of coupling interface.
     */
    bool isCoupledEntity(const int faceID) const;
    /*!
     * @brief Get a quantity's identifier from its name.
     *
     * @param[in] meshName Name of the mesh.
     * @param[in] dataName Name of the quantity.
     * @return size_t Numeric identifier of quantity.
     */
    std::string meshAndDataKey(const std::string &meshName,
                               const std::string &dataName) const;
    /*!
     * @brief Prints status of coupling adapter to given output stream.
     *
     * @param os Output stream.
     */
    void print(std::ostream &os);
};

}  // namespace Dumux::Precice
#endif
