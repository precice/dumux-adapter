#ifndef PRECICEWRAPPER_HH
#define PRECICEWRAPPER_HH

#include <map>
#include <ostream>
#include <precice/precice.hpp>
#include <string>

#include "dumuxpreciceindexmapper.hh"
#include "solverstate.hh"

/*!
 * @brief Namespace of dumux-precice
 *
 */
namespace Dumux::Precice
{

//! Type of Dumux face IDs
using FaceID = int;

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
    //! Map storing <meshName, dataName> and data vectors for reading
    std::map<std::pair<std::string, std::string>, std::vector<double>>
        dataRead_;
    //! Map storing <meshName, dataName> and data vectors for writing
    std::map<std::pair<std::string, std::string>, std::vector<double>>
        dataWrite_;
    //! Vector of identifiers (in preCICE) of the vertices of the coupling mesh.
    std::vector<int> vertexIDs_;  //should be size_t
    //! Configuration file name and psth
    std::string preciceConfigName_;
    //! Participant or solver name
    std::string participantName_;
    //! Constructor
    CouplingAdapter();
    /*!
     * @brief Instance of DumuxPreciceIndexMapper that translates between
     *        DuMuX' identifiers of vertices and preCICE's identifiers.
     *
     */
    Internal::DumuxPreciceIndexMapper<FaceID, precice::VertexID> indexMapper_;
    /*!
     * @brief Store the states of the solver for checkpointing.
     *
     */
    std::vector<std::unique_ptr<SolverStateBase>> states_;
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
     * @brief Read in configuration and add the paired mesh and data names into respective maps, initilize with empty vector
     *
    */
    void announceConfig(const int rank, const int size);
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
     * @brief Get the participant name from config
     *
     * @return vector of datanames
     */
    std::string getSolverName() const;
    /*!
     * @brief Get the coupled meshnames on this participant
     *
     * @return vector of neshnames
     */
    std::vector<std::string> getMeshNames() const;
    /*!
     * @brief Get the datanames for reading on this mesh
     *
     * @return vector of datanames
     */
    std::vector<std::string> getReadDataNamesOnMesh(
        const std::string &meshName) const;
    /*!
     * @brief Get the datanames for writing on this mesh
     *
     * @return vector of datanames
     */
    std::vector<std::string> getWriteDataNamesOnMesh(
        const std::string &meshName) const;

    /*!
     * @brief Initializes the checkpointing functionality.
     *
     * This function needs to be called at least once before using the checkpointing functionality.
     *
     * @param[in] x Solution vector
     * @param[in] gv Grid variables
     * @param[in] tl Time loop
     */
    template<class SolutionVector, class GridVariables, class TimeLoop>
    void initializeCheckpoint(SolutionVector &x,
                              GridVariables &gv,
                              TimeLoop &tl);

    template<class SolutionVector, class GridVariables>
    void initializeCheckpoint(SolutionVector &x, GridVariables &gv);

    template<class SolutionVector>
    void initializeCheckpoint(SolutionVector &x);
    /*!
     * @brief Writes the solver state to a checkpoint if required by preCICE.
     *
     */
    bool writeCheckpointIfRequired();
    /*!
     * @brief Reads the solver state from a checkpoint if required by preCICE.
     */
    bool readCheckpointIfRequired();

    /*!
     * @brief Checks if the participant is required to provide initial data. If true, the participant needs to write initial data to defined vertices prior to calling initialize().
     *
     * @return true Initial coupling data has to be provided.
     * @return false No further action is needed.
     */
    bool requiresToWriteInitialData();

    /*!
     * @brief Adds mesh for coupling of solvers. With the mesh size, the data maps inside the adapter initialize the relevant data vector to size of meshSize*dataDimension.
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
    void setMesh(const std::string &meshName,
                 const std::vector<double> &positions);
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
    void createIndexMapping(const std::vector<FaceID> &dumuxFaceIDs);
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
     * @brief Checks whether the time window has completed.
     *
     * @return true Time window has completed.
     * @return false Time window is still ongoing.
     */
    bool isTimeWindowComplete();
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
                                   const FaceID faceID);
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
                                   const FaceID faceID,
                                   const double value);
    /*!
     * @brief Writes the quantity value vector into the data map.
     *
     * @param[in] meshName Name of the mesh.
     * @param[in] dataName Name of the data.
     * @param[in] values Value of the scalar or vector quantity.
     */
    void writeQuantityVector(const std::string &meshName,
                             const std::string &dataName,
                             const std::vector<double> &values);
    /*!
     * @brief Checks whether face with given identifier is part of coupling interface.
     *
     * @param[in] faceID Identifier of the face according to DuMuX' numbering.
     * @return true Face is part of coupling interface.
     * @return false Face is not part of coupling interface.
     */
    bool isCoupledEntity(const int faceID) const;
    /*!
     * @brief Prints status of coupling adapter to given output stream.
     *
     * @param os Output stream.
     */
    void print(std::ostream &os);
};

template<class SolutionVector, class GridVariables, class TimeLoop>
void CouplingAdapter::initializeCheckpoint(SolutionVector &x,
                                           GridVariables &gv,
                                           TimeLoop &tl)
{
    states_.emplace_back(
        std::make_unique<SolverStateGridVarTimeLoop<SolutionVector,
                                                    GridVariables, TimeLoop>>(
            x, gv, tl));
}

template<class SolutionVector, class GridVariables>
void CouplingAdapter::initializeCheckpoint(SolutionVector &x, GridVariables &gv)
{
    states_.emplace_back(
        std::make_unique<SolverStateGridVar<SolutionVector, GridVariables>>(
            x, gv));
}

template<class SolutionVector>
void CouplingAdapter::initializeCheckpoint(SolutionVector &x)
{
    states_.emplace_back(std::make_unique<SolverStateOnly<SolutionVector>>(x));
}
}  // namespace Dumux::Precice
#endif
