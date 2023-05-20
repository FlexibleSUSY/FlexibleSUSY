// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

#ifndef TEST_RANDOM_SIGN_DATASET_H
#define TEST_RANDOM_SIGN_DATASET_H

#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <random>

namespace flexiblesusy
{
namespace boost_test_tools
{
class random_sign_dataset
{
public:
  using sample = int;
  enum { arity = 1 };

  class iterator
  {
    std::mt19937 generator;
    std::bernoulli_distribution distribution;
    
    int sign;
    int random_sign( void )
    { return distribution( generator ) ? -1 : 1; }
  public:
    iterator() :
    generator( std::random_device{}() ), distribution( 0.5 ),
    sign(random_sign()) {}

    int operator*() const { return sign; }
    void operator++()
    { sign = random_sign(); }
  };
  
  boost::unit_test::data::size_t size() const
  { return boost::unit_test::data::BOOST_TEST_DS_INFINITE_SIZE; }

  iterator begin() const { return iterator(); }
};
}
}

namespace boost
{
namespace unit_test
{
namespace data
{
namespace monomorphic
{
template<>
struct is_dataset<
	flexiblesusy::boost_test_tools::random_sign_dataset
> : boost::mpl::true_ {};
}
}
}
}

#endif
