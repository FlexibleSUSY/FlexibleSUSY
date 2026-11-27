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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_flexibledecay_problems

#include <boost/test/unit_test.hpp>

#include "decays/flexibledecay_problems.hpp"

#include <string>
#include <vector>

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE(test_no_problems)
{
   const FlexibleDecay_problems problems;

   BOOST_CHECK(!problems.have_problem());
   BOOST_CHECK(!problems.have_error());
   BOOST_CHECK(!problems.have_warning());
   BOOST_CHECK(problems.get_problem_strings().empty());
   BOOST_CHECK(problems.get_warning_strings().empty());
}

BOOST_AUTO_TEST_CASE(test_warning_is_not_a_problem)
{
   FlexibleDecay_problems problems;
   problems.add_warning("warning");

   BOOST_CHECK(!problems.have_problem());
   BOOST_CHECK(!problems.have_error());
   BOOST_CHECK(problems.have_warning());
   BOOST_CHECK(problems.get_warning_strings() ==
               std::vector<std::string>{"warning"});
}

// e.g. a channel set to zero because its width isn't finite
BOOST_AUTO_TEST_CASE(test_error_is_not_a_problem)
{
   FlexibleDecay_problems problems;
   problems.add_error("error");

   BOOST_CHECK(!problems.have_problem());
   BOOST_CHECK(problems.have_error());
   BOOST_CHECK(!problems.have_warning());
   BOOST_CHECK(problems.get_problem_strings() ==
               std::vector<std::string>{"error"});
}

// an error that interrupted the decay calculation
BOOST_AUTO_TEST_CASE(test_fatal_error_is_a_problem)
{
   FlexibleDecay_problems problems;
   problems.add_error("error");
   problems.add_fatal_error("fatal error");

   BOOST_CHECK(problems.have_problem());
   BOOST_CHECK(problems.have_error());
   BOOST_CHECK(!problems.have_warning());
   BOOST_CHECK((problems.get_problem_strings() ==
                std::vector<std::string>{"error", "fatal error"}));
}

BOOST_AUTO_TEST_CASE(test_clear)
{
   FlexibleDecay_problems problems;
   problems.add_warning("warning");
   problems.add_fatal_error("error");
   problems.clear();

   BOOST_CHECK(!problems.have_problem());
   BOOST_CHECK(!problems.have_error());
   BOOST_CHECK(!problems.have_warning());
   BOOST_CHECK(problems.get_problem_strings().empty());
   BOOST_CHECK(problems.get_warning_strings().empty());
}
