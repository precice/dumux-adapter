#include "thirdparty/catch2/catch.hpp"

#include <random>
#include <vector>
#include "../src/dumuxpreciceindexwrapper.hh"

template<typename T>
void fillVectors(std::vector<T> &preciceIndices,
                 std::vector<T> &solverIndices,
                 const T offset)
{
    REQUIRE(preciceIndices.size() == solverIndices.size());
    const size_t numIndices = preciceIndices.size();

    for (size_t i = 0; i < numIndices; ++i) {
        preciceIndices[i] = i;
        solverIndices[i] = i + offset;
    }
}

template<typename T>
void checkDumuxIdMapping(DumuxPreciceIndexMapper<T> &indexMapper,
                         std::vector<T> &preciceIndices,
                         std::vector<T> &solverIndices)
{
    const size_t numIndices = preciceIndices.size();

    for (size_t i = 0; i < numIndices; ++i) {
        REQUIRE(indexMapper.getDumuxId(solverIndices[i]) == preciceIndices[i]);
    }
}

template<typename T>
void checkPreciceIdMapping(DumuxPreciceIndexMapper<T> &indexMapper,
                           std::vector<T> &preciceIndices,
                           std::vector<T> &solverIndices)
{
    const size_t numIndices = preciceIndices.size();

    for (size_t i = 0; i < numIndices; ++i) {
        REQUIRE(indexMapper.getPreciceId(preciceIndices[i]) ==
                solverIndices[i]);
    }
}

template<typename T>
void checkMapper(DumuxPreciceIndexMapper<T> &indexMapper,
                 std::vector<T> &preciceIndices,
                 std::vector<T> &solverIndices,
                 const T numIndices)
{
    indexMapper.createMapping(preciceIndices, solverIndices);

    SECTION("Check size") { REQUIRE(indexMapper.getSize() == numIndices); }

    SECTION("Mapping with offset get solver ID")
    {
        checkDumuxIdMapping(indexMapper, preciceIndices, solverIndices);
    }

    SECTION("Mapping with offset get precice ID")
    {
        checkPreciceIdMapping(indexMapper, preciceIndices, solverIndices);
    }
}

TEST_CASE("DumuxPreciceIndexMapper mapping", "[DumuxPreciceIndexMapper]")
{
    DumuxPreciceIndexMapper<size_t> indexMapper;

    const size_t numIndices = 200;

    std::vector<size_t> preciceIndices(numIndices);
    std::vector<size_t> solverIndices(numIndices);

    REQUIRE(preciceIndices.size() == numIndices);
    REQUIRE(solverIndices.size() >= numIndices);

    SECTION("Mapping with offset 1")
    {
        const size_t offset = 1;

        fillVectors(preciceIndices, solverIndices, offset);
        checkMapper(indexMapper, preciceIndices, solverIndices, numIndices);
    }

    SECTION("Mapping with offset 5")
    {
        const size_t offset = 5;

        fillVectors(preciceIndices, solverIndices, offset);
        checkMapper(indexMapper, preciceIndices, solverIndices, numIndices);
    }

    SECTION("Mapping with offset 5")
    {
        const size_t offset = 3;

        fillVectors(preciceIndices, solverIndices, offset);
        std::mt19937 g(42);
        std::shuffle(preciceIndices.begin(), preciceIndices.end(), g);
        std::shuffle(solverIndices.begin(), solverIndices.end(), g);

        checkMapper(indexMapper, preciceIndices, solverIndices, numIndices);
    }
}
