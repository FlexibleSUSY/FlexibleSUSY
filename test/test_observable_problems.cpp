#include "observable_problems.hpp"
#include <iostream>

int errors = 0;

#define CHECK(cond) do {                                                \
      if (!(cond)) {                                                    \
         std::cerr << "Error in line " << __LINE__ << std::endl;        \
         errors++;                                                      \
      }                                                                 \
   } while (false)


void print(flexiblesusy::Observable_problems& op)
{
   const auto strings = op.get_problem_strings();

   for (const auto& s: strings) {
      std::cout << s << std::endl;
   }
}


void test_empty()
{
   flexiblesusy::Observable_problems op;
   CHECK(!op.have_problem());
   CHECK(op.get_problem_strings().empty());
}


void test_error()
{
   flexiblesusy::Observable_problems op;
   op.a_muon.flag_non_perturbative_running(1.0);
   CHECK(op.have_problem());
   CHECK(op.get_problem_strings().size() == 1);
}


void test_print()
{
   flexiblesusy::Observable_problems op;
   op.a_muon.flag_non_perturbative_running(1.0);
   print(op);
}


int main()
{
   test_empty();
   test_error();
   test_print();

   return errors;
}
