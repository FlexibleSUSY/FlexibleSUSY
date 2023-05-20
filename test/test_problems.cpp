
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
