#include "observable_problems.hpp"
#include "observable_problems_format.hpp"
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

int errors = 0;

#define CHECK(cond) do {                                                \
      if (!(cond)) {                                                    \
         std::cerr << "Error in line " << __LINE__ << std::endl;        \
         errors++;                                                      \
      }                                                                 \
   } while (false)


int count(const flexiblesusy::Observable_problems& op)
{
   std::vector<std::string> str;
   copy_problem_strings(op, std::back_insert_iterator<std::vector<std::string>>(str));
   return str.size();
}


void print(const flexiblesusy::Observable_problems& op)
{
   copy_problem_strings(op, std::ostream_iterator<std::string>(std::cout, "\n"));
}


void test_empty()
{
   flexiblesusy::Observable_problems op;
   CHECK(!op.have_problem());
   CHECK(count(op) == 0);
}


void test_error()
{
   flexiblesusy::Observable_problems op;
   op.a_muon.flag_non_perturbative_running(1.0);
   CHECK(op.have_problem());
   CHECK(count(op) == 1);
}


void test_print()
{
   flexiblesusy::Observable_problems op;
   op.a_muon.flag_non_perturbative_running(1.0);
   print(op);
}


void test_slha()
{
   flexiblesusy::Observable_problems op;
   op.a_muon.flag_non_perturbative_running(1.0);

   copy_problem_strings(op, flexiblesusy::observable_problems::SLHA_observable_problems_output_iterator<std::string>(std::cout, 1, 3));
}


int main()
{
   test_empty();
   test_error();
   test_print();
   test_slha();

   return errors;
}
