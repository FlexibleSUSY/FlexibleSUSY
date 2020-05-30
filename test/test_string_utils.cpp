#include "string_utils.hpp"
#include <iostream>
#include <vector>


int errors = 0;


#define CHECK_EQ(a,b)                                                   \
   do {                                                                 \
      if (!(a == b)) {                                                  \
         std::cout << "FAILED test in line " << __LINE__ << ": "        \
                   << a << " != " << b << std::endl;                    \
         errors++;                                                      \
      }                                                                 \
   } while (false)


void test_concat()
{
   using namespace flexiblesusy;

   CHECK_EQ(concat({"1", "2", "3"}),  "123");
}


void test_concat_with_separator()
{
   using namespace flexiblesusy;

   CHECK_EQ(concat({"1", "2", "3"}, ','), "1,2,3");
   CHECK_EQ(concat({"1", "2", "3"}, ","), "1,2,3");
}


void test_concat_ints_with_separator()
{
   using namespace flexiblesusy;

   std::vector<int> vec = {1, 2, 3};

   CHECK_EQ(concat(vec.cbegin(), vec.cend(), ","), "1,2,3");
   CHECK_EQ(concat(vec.cbegin(), vec.cend(), ','), "1,2,3");
}


int main()
{
   test_concat();
   test_concat_with_separator();
   test_concat_ints_with_separator();

   return errors;
}
