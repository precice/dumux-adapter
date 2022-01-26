#ifndef DUMUXPRECICEINDEXWRAPPER_H
#define DUMUXPRECICEINDEXWRAPPER_H

#include <cassert>
#include <iostream>
#include <map>
#include <ostream>
#include <vector>

template<typename T>
class DumuxPreciceIndexMapper
{
private:
    // Use a boost bimap?
    std::map<T, T> dumuxFaceIndexToPreciceIndex_;
    std::map<T, T> preciceVertexToDumuxFaceIndex_;
    // AJ: I think we should not rely on the fact that the preCICE vertexIDs start at 0.
    // It might happen that we have more than one interface?! We should check the performance.
    //std::vector<T> preciceVertexToDumuxFaceIndex_;
    bool mappingWasCreated_;
    size_t size_;

public:
    DumuxPreciceIndexMapper() : mappingWasCreated_(false), size_(0) {}

    void createMapping(const std::vector<T> &dumuxIndices,
                       const std::vector<T> &preciceIndices)
    {
        assert(dumuxIndices.size() == preciceIndices.size());
        size_ = dumuxIndices.size();

        for (T i = 0; i < size_; i++) {
            preciceVertexToDumuxFaceIndex_.emplace(preciceIndices[i],
                                                   dumuxIndices[i]);
            dumuxFaceIndexToPreciceIndex_.emplace(dumuxIndices[i],
                                                  preciceIndices[i]);
        }
        mappingWasCreated_ = true;
    }

    bool isDumuxIdMapped(const T dumuxId) const
    {
        return dumuxFaceIndexToPreciceIndex_.count(dumuxId) == 1;
    }

    bool isPreciceIdMapped(const T preciceId) const
    {
        return preciceVertexToDumuxFaceIndex_.count(preciceId) == 1;
    }

    const T getPreciceId(const T dumuxId) const
    {
        assert(isDumuxIdMapped(dumuxId));
        return dumuxFaceIndexToPreciceIndex_.at(dumuxId);
    }

    const T getDumuxId(const T preciceId) const
    {
        assert(isPreciceIdMapped(preciceId));
        return preciceVertexToDumuxFaceIndex_.at(preciceId);
    }

    size_t getSize() const { return preciceVertexToDumuxFaceIndex_.size(); }

    virtual ~DumuxPreciceIndexMapper() {}

    template<typename U>
    friend std::ostream &operator<<(std::ostream &os,
                                    const DumuxPreciceIndexMapper<U> &wrapper);
};

template<typename T>
std::ostream &operator<<(std::ostream &os,
                         const DumuxPreciceIndexMapper<T> &wrapper)
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

#endif  // DUMUXPRECICEINDEXWRAPPER_H
