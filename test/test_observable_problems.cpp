#include "observable_problems.hpp"
#include "observable_problems_format_slha.hpp"
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


flexiblesusy::Observable_problems setup_errors()
{
   flexiblesusy::Observable_problems op;

   op.general.flag_thrown(nullptr);
   op.general.flag_non_perturbative_running(1.0);
   op.a_muon.flag_non_perturbative_running(1.0);

   return op;
}


int count_lines(const flexiblesusy::Observable_problems& op)
{
   std::vector<std::string> str;
   copy_problem_strings(op, std::back_insert_iterator<std::vector<std::string>>(str));
   return str.size();
}


void test_empty()
{
   flexiblesusy::Observable_problems op;
   CHECK(!op.have_problem());
   CHECK(op.number_of_problems() == 0);
   CHECK(count_lines(op) == 0);
}


void test_error()
{
   flexiblesusy::Observable_problems op = setup_errors();
   CHECK(op.have_problem());
   CHECK(op.number_of_problems() == 3);
   CHECK(count_lines(op) == 3);
}


void test_print()
{
   flexiblesusy::Observable_problems op = setup_errors();
   copy_problem_strings(op, std::ostream_iterator<std::string>(std::cout, "\n"));
}


void test_slha()
{
   flexiblesusy::Observable_problems op = setup_errors();
   slha_format_problems_and_warnings(op, std::ostream_iterator<std::string>(std::cout));
}


void test_ostream()
{
   flexiblesusy::Observable_problems op = setup_errors();
   std::cout << op;
}


int main()
{
   std::cout << "----- test empty struct ----------------------\n";
   test_empty();
   std::cout << "----- test problems --------------------------\n";
   test_error();
   std::cout << "----- test ostream_iterator ------------------\n";
   test_print();
   std::cout << "----- test SLHA output -----------------------\n";
   test_slha();
   std::cout << "----- test ostream output --------------------\n";
   test_ostream();

   std::cout << '\n';
   std::cout << "===== Summary ================================\n";
   if (errors == 0) {
      std::cout << "All tests PASSED\n";
   } else {
      std::cout << "Failed tests: " << errors << '\n';
   }
   std::cout << "==============================================\n";

   return errors;
}
