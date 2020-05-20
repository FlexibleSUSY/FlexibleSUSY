#include "string_conversion.hpp"
#include <limits>
#include <iostream>
#include <string>

int errors = 0;


#define CHECK(cond)                                                     \
   do {                                                                 \
      if (!(cond)) {                                                    \
         std::cout << "FAILED test in line " << __LINE__ << std::endl;  \
         errors++;                                                      \
      }                                                                 \
   } while (false)


#define CHECK_THROW(expr)                       \
   do {                                         \
      try {                                     \
         expr;                                  \
         std::cout << "FAILED test in line " << __LINE__ << std::endl;  \
         errors++;                              \
      } catch (...) {}                          \
   } while (false)


void test_to_int()
{
   using namespace flexiblesusy;

   CHECK(to_int("0")  ==  0);
   CHECK(to_int("1")  ==  1);
   CHECK(to_int("-1") == -1);
   CHECK(to_int(std::to_string(std::numeric_limits<int>::min()).c_str()) == std::numeric_limits<int>::min());
   CHECK(to_int(std::to_string(std::numeric_limits<int>::max()).c_str()) == std::numeric_limits<int>::max());

   CHECK_THROW(to_int(std::to_string(std::numeric_limits<long>::max()).c_str()));
   CHECK_THROW(to_int("a"));
   CHECK_THROW(to_int("1a"));
}


void test_to_double()
{
   using namespace flexiblesusy;

   CHECK(to_double("0")    ==  0.0);
   CHECK(to_double("0.0")  ==  0.0);
   CHECK(to_double("1")    ==  1.0);
   CHECK(to_double("1.0")  ==  1.0);
   CHECK(to_double("-1")   == -1.0);
   CHECK(to_double("-1.0") == -1.0);
   CHECK(to_double(std::to_string(std::numeric_limits<long>::min()).c_str()) == static_cast<double>(std::numeric_limits<long>::min()));
   CHECK(to_double(std::to_string(std::numeric_limits<long>::max()).c_str()) == static_cast<double>(std::numeric_limits<long>::max()));
   // CHECK(to_double(std::to_string(std::numeric_limits<double>::min()).c_str()) == std::numeric_limits<double>::min());
   CHECK(to_double(std::to_string(std::numeric_limits<double>::max()).c_str()) == std::numeric_limits<double>::max());

   CHECK(to_double("inf") == std::numeric_limits<double>::infinity());
   // CHECK(to_double("nan") == std::numeric_limits<double>::quiet_NaN());

   CHECK_THROW(to_double("a"));
   CHECK_THROW(to_double("1.0a"));
}


int main()
{
   test_to_int();
   test_to_double();

   std::cout << "==================================\n";
   if (errors) {
      std::cout << "FAILED tests: " << errors << std::endl;
   } else {
      std::cout << "All tests PASSED\n";
   }
   std::cout << "==================================\n";

   return errors;
}
