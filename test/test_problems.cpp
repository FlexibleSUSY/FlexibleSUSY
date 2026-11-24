
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_problems

#include <boost/test/unit_test.hpp>

#include "problems.hpp"
#include "names.hpp"

using namespace flexiblesusy;

class Dummy_names : public Names {
public:
   Dummy_names(int size) : s(size) {}
   virtual ~Dummy_names() = default;
   virtual const std::string& get(int) const { return name; }
   virtual int size() const { return s; }
private:
   int s{};
   std::string name{"P"};
};

BOOST_AUTO_TEST_CASE( test_initialization )
{
   const Dummy_names dummy_names(3);
   Problems problems("DummyModel", &dummy_names, &dummy_names);

   BOOST_CHECK(!problems.is_running_tachyon(0));
   BOOST_CHECK(!problems.is_running_tachyon(1));
   BOOST_CHECK(!problems.is_running_tachyon(2));
   BOOST_CHECK(!problems.have_tachyon());

   BOOST_CHECK(!problems.no_ewsb());
   BOOST_CHECK(!problems.no_perturbative());
}

BOOST_AUTO_TEST_CASE(test_failed_rg_integration)
{
   const Dummy_names dummy_names(3);
   Problems problems("DummyModel", &dummy_names, &dummy_names);

   BOOST_CHECK(!problems.have_failed_rg_integration());

   problems.flag_failed_rg_integration("too many steps", 1, 1e-14, 117.);

   BOOST_CHECK(problems.have_failed_rg_integration());
   BOOST_CHECK(problems.have_problem());
   BOOST_CHECK(!problems.no_perturbative());
   BOOST_CHECK(!problems.have_non_perturbative_parameter());
   BOOST_CHECK_EQUAL(problems.number_of_problems(), 1);

   const auto strings = problems.get_problem_strings();
   BOOST_REQUIRE_EQUAL(strings.size(), 1);
   BOOST_TEST_MESSAGE(strings.at(0));
   BOOST_CHECK(strings.at(0).find("RG integration failed (too many steps)") !=
               std::string::npos);
   BOOST_CHECK(strings.at(0).find("step size limited by P") !=
               std::string::npos);

   Problems other("DummyModel", &dummy_names, &dummy_names);
   other.add(problems);
   BOOST_CHECK(other.have_failed_rg_integration());

   problems.unflag_failed_rg_integration();
   BOOST_CHECK(!problems.have_failed_rg_integration());
   BOOST_CHECK(!problems.have_problem());

   other.clear();
   BOOST_CHECK(!other.have_failed_rg_integration());
}

BOOST_AUTO_TEST_CASE(test_out_of_range_parameter_index)
{
   const Dummy_names dummy_names(3);
   Problems problems("DummyModel", &dummy_names, &dummy_names);

   // an index out of range must not be used to look up a name
   problems.flag_non_perturbative_parameter(5, 1., 100.);
   problems.flag_failed_rg_integration("too many steps", 5, 1., 100.);

   const auto strings = problems.get_problem_strings();
   BOOST_REQUIRE_EQUAL(strings.size(), 2);
   BOOST_TEST_MESSAGE(strings.at(0));
   BOOST_TEST_MESSAGE(strings.at(1));
   BOOST_CHECK(strings.at(0).find("non-perturbative Q") != std::string::npos);
   BOOST_CHECK(strings.at(1).find("step size limited by") == std::string::npos);
}
