#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_HiggsTools_CP

#include <boost/test/unit_test.hpp>

#include "Higgs/Predictions.hpp"

namespace HP = Higgs::predictions;

BOOST_AUTO_TEST_CASE( test_HiggsTools_CP )
{
   BOOST_CHECK (static_cast<HP::CP>(1) == HP::CP::even);
   BOOST_CHECK (static_cast<HP::CP>(-1) == HP::CP::odd);
   BOOST_CHECK (static_cast<HP::CP>(0) == HP::CP::undefined);
}
