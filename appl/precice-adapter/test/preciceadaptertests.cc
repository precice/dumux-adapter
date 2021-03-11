#include "../include/preciceadapter.hh"
#include "thirdparty/catch2/catch.hpp"

TEST_CASE("PreciceAdapter without data mapping", "[PreciceAdapter]")
{
    using namespace precice_adapter;
    PreciceAdapter &couplingInterface = PreciceAdapter::getInstance();

    couplingInterface.announceSolver("SolverA", "precice-basic-test-config.xml",
                                     0, 1);

    REQUIRE(couplingInterface.getDimensions() == 2);

    const size_t dim = couplingInterface.getDimensions();
    const size_t meshSize = 20;
    std::vector<double> coords(meshSize * dim);
    for (size_t i = 0; i < meshSize; ++i) {
        coords[dim * i] = double(i);
        coords[dim * i + 1] = 0.;
    }

    couplingInterface.setMesh("SolverAMesh", meshSize, coords);

    //  SECTION("Check: Get data ID from string")
    {
        REQUIRE(couplingInterface.getNumberOfVertices() == meshSize);
    }

    const auto fooId = couplingInterface.announceQuantity("Foo");
    const auto barId = couplingInterface.announceQuantity("Bar");

    //  SECTION("Check: Get data ID from string")
    {
        REQUIRE(couplingInterface.getIdFromName("Foo") == fooId);
        REQUIRE(couplingInterface.getIdFromName("Bar") == barId);
    }

    //  SECTION("Check: Get data name from ID")
    {
        REQUIRE(couplingInterface.getNameFromId(fooId) == "Foo");
        REQUIRE(couplingInterface.getNameFromId(barId) == "Bar");
    }

    //fooId data
    {
        //Should be zero since unitilialized
        const auto &readData = couplingInterface.getQuantityVector(fooId);

        REQUIRE(readData.size() == meshSize);
        for (size_t i = 0; i < meshSize; ++i) {
            REQUIRE(readData[i] == Approx(0.));
        }
    }

    {
        // Scalar data
        std::vector<double> data(meshSize);
        for (size_t i = 0; i < meshSize; ++i) {
            data[i] = double(i);
        }
        std::mt19937 g(42);
        std::shuffle(data.begin(), data.end(), g);

        couplingInterface.writeScalarQuantityVector(fooId, data);
        const auto &readData = couplingInterface.getQuantityVector(fooId);

        REQUIRE(readData.size() == data.size());
        for (size_t i = 0; i < meshSize; ++i) {
            REQUIRE(readData[i] == Approx(data[i]));
        }
    }

    //barId data
    {
        //Should be zero since unitilialized
        const auto &readData = couplingInterface.getQuantityVector(barId);

        REQUIRE(readData.size() == meshSize);
        for (size_t i = 0; i < meshSize; ++i) {
            REQUIRE(readData[i] == Approx(0.));
        }
    }

    {
        // Scalar data
        std::vector<double> data(meshSize);
        for (size_t i = 0; i < meshSize; ++i) {
            data[i] = double(i);
        }
        std::mt19937 g(42);
        std::shuffle(data.begin(), data.end(), g);

        couplingInterface.writeScalarQuantityVector(barId, data);
        const auto &readData = couplingInterface.getQuantityVector(barId);

        REQUIRE(readData.size() == data.size());
        for (size_t i = 0; i < meshSize; ++i) {
            REQUIRE(readData[i] == Approx(data[i]));
        }
    }

    couplingInterface.finalize();
}
