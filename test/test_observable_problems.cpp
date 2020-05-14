#include "observable_problems.hpp"
#include <iostream>
#include <iterator>
#include <string>
#include <vector>
#include <boost/format.hpp>

int errors = 0;

#define CHECK(cond) do {                                                \
      if (!(cond)) {                                                    \
         std::cerr << "Error in line " << __LINE__ << std::endl;        \
         errors++;                                                      \
      }                                                                 \
   } while (false)


std::vector<std::string> problem_strings(const flexiblesusy::Observable_problems& op)
{
   std::vector<std::string> str;
   op.copy_problem_strings(std::back_insert_iterator<std::vector<std::string>>(str));
   return str;
}

// should be moved to SLHA class
template <typename T>
class SLHA_output_iterator {
public:
   SLHA_output_iterator(std::ostream& ostr_, int idx1_, int idx2_)
      : ostr(ostr_), idx1(idx1_), idx2(idx2_) {}

   void operator=(const T& elem) {
      ostr << boost::format(" %5d %5d   %s\n") % idx1 % idx2 % elem;
   }
   void operator++(int) {}
private:
   std::ostream& ostr;
   int idx1, idx2;
};


void print(const flexiblesusy::Observable_problems& op)
{
   op.copy_problem_strings(std::ostream_iterator<const char*>(std::cout, "\n"));
}


void test_empty()
{
   flexiblesusy::Observable_problems op;
   CHECK(!op.have_problem());
   CHECK(problem_strings(op).empty());
}


void test_error()
{
   flexiblesusy::Observable_problems op;
   op.a_muon.flag_non_perturbative_running(1.0);
   CHECK(op.have_problem());
   CHECK(problem_strings(op).size() == 1);
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

   op.copy_problem_strings(SLHA_output_iterator<const char*>(std::cout, 1, 3));
}


int main()
{
   test_empty();
   test_error();
   test_print();
   test_slha();

   return errors;
}
