#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "thirdparty/catch2/catch.hpp"

unsigned int Factorial(unsigned int number)
{
    return number > 1 ? Factorial(number - 1) * number : 1;
}

TEST_CASE("Factorials are computed", "[factorial]")
{
    REQUIRE(Factorial(0) == 1);
    REQUIRE(Factorial(1) == 1);
    REQUIRE(Factorial(2) == 2);
    REQUIRE(Factorial(3) == 6);
    REQUIRE(Factorial(10) == 3628800);
}
