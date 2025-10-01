#ifndef DUMUXPRECICEINDEXWRAPPER_H
#define DUMUXPRECICEINDEXWRAPPER_H

#include <cassert>
#include <ostream>
#include <unordered_map>
#include <vector>

/*!
 * @brief Namespace of dumux-precice internals
 *
 * The namespace contains classes and implementations that should not
 * be needed by the users themselves.
 *
 */
namespace Dumux::Precice::Internal
{
/*!
 * @brief Mapping between preCICE vertex indices and DuMuX face indices.
 *
 * @tparam FaceID Data type of the Dumux face ID.
 * @tparam VertexID Data type of the preCICE vertex ID.
 */
template<typename FaceID, typename VertexID>
class DumuxPreciceIndexMapper
{
private:
    //!  Mapping from Dumux' face indices to preCICE's vertex indices.
    std::unordered_map<FaceID, VertexID> dumuxFaceIndexToPreciceIndex_;
    //!  Mapping from preCICE's vertex indices to Dumux' face indices.
    std::unordered_map<VertexID, FaceID> preciceVertexToDumuxFaceIndex_;

public:
    /*!
     * @brief Creates/computes the bijective mapping between preCICE's
     *        vertex and DuMuX' face indices.
     *
     * \note It is assumed that the order of vertex and face indices in the two
     *       input parameters is identical. That means the first entry in
     *       `dumuxIndices` should be connected to the first entry in
     *       `preciceIndices`.
     *
     * @param[in] dumuxIndices Vector of DuMuX' face indices.
     * @param[in] preciceIndices Vector of preCICE's vertex indices.
     */
    void createMapping(const std::vector<FaceID> &dumuxIndices,
                       const std::vector<VertexID> &preciceIndices)
    {
        assert(dumuxIndices.size() == preciceIndices.size());
        const size_t size_ = dumuxIndices.size();
        preciceVertexToDumuxFaceIndex_.reserve(size_);
        dumuxFaceIndexToPreciceIndex_.reserve(size_);

        for (size_t i = 0; i < size_; i++) {
            preciceVertexToDumuxFaceIndex_.emplace(preciceIndices[i],
                                                   dumuxIndices[i]);
            dumuxFaceIndexToPreciceIndex_.emplace(dumuxIndices[i],
                                                  preciceIndices[i]);
        }
    }
    /*!
     * @brief Gets preCICE's vertex index based on a DuMuX face index.
     *
     * @param[in] dumuxId DuMuX face index.
     * @return const T preCICE vertex index.
     */
    const VertexID getPreciceId(const FaceID dumuxId) const
    {
        assert(isDumuxIdMapped(dumuxId));
        return dumuxFaceIndexToPreciceIndex_.at(dumuxId);
    }
    /*!
     * @brief Gets DuMuX' face index basde on a preCICE vertex index.
     *
     * @param[in] preciceId preCICE vertex index.
     * @return const T DuMuX face index.
     */
    const FaceID getDumuxId(const VertexID preciceId) const
    {
        assert(isPreciceIdMapped(preciceId));
        return preciceVertexToDumuxFaceIndex_.at(preciceId);
    }
    /*!
     * @brief Checks if a DuMuX face index is mapped to a preCICE vertex index.
     *
     * @param[in] dumuxId DuMuX face index.
     * @return true DuMuX face index is mapped to a preCICE vertex index.
     * @return false No mapping for the given index available.
     */
    bool isDumuxIdMapped(const FaceID dumuxId) const
    {
        return dumuxFaceIndexToPreciceIndex_.count(dumuxId) == 1;
    }
    /*!
     * @brief Checkes if a preCICE vertex index is mapped to a DuMuX face index.
     *
     * @param[in] preciceId preCICE vertex index.
     * @return true preCICE vertex index is mapped to a DuMuX face index.
     * @return false  No mapping for the given index available.
     */
    bool isPreciceIdMapped(const VertexID preciceId) const
    {
        return preciceVertexToDumuxFaceIndex_.count(preciceId) == 1;
    }
    /*!
     * @brief Gets the size of the mapping table
     *
     * @return size_t Number of face/vertex indices mapped.
     */
    size_t getSize() const { return preciceVertexToDumuxFaceIndex_.size(); }
    /*!
     * @brief Destructor
     *
     */
    virtual ~DumuxPreciceIndexMapper() {}

    /*!
     * @brief Prints state of the DumuxPreciceIndexMapper object to the given outstream.
     *
     * @tparam D Data type of the Dumux indices.
     * @tparam P Data type of the preCICE indices.
     * @param os Outstream.
     * @param wrapper DumuxPreciceIndexMapper object that should be printed.
     * @return std::ostream& Outstream.
     */
    template<typename D, typename P>
    friend std::ostream &operator<<(
        std::ostream &os,
        const DumuxPreciceIndexMapper<D, P> &wrapper);
};

template<typename D, typename P>
std::ostream &operator<<(std::ostream &os,
                         const DumuxPreciceIndexMapper<D, P> &wrapper)
{
    os << "preCICE to DuMuX mapping "
       << "\n";
    for (const auto &v : wrapper.preciceVertexToDumuxFaceIndex_) {
        os << v.first << " -> " << wrapper.getDumuxId(v.first) << "\n";
    }

    os << "\n\n";
    os << "Dumux to preCICE mapping "
       << "\n";
    for (const auto &v : wrapper.dumuxFaceIndexToPreciceIndex_) {
        os << v.first << " -> " << wrapper.getPreciceId(v.first) << "\n";
    }

    return os;
}
}  // namespace Dumux::Precice::Internal

#endif  // DUMUXPRECICEINDEXWRAPPER_H
