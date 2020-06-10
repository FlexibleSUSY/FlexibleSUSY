#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_gsl_vector

#include <boost/test/unit_test.hpp>

#include "gsl_vector.hpp"
#include "error.hpp"

#include <limits>

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_init_default )
{
   GSL_vector v;
   BOOST_CHECK_EQUAL(v.size(), 0);
}

BOOST_AUTO_TEST_CASE( test_init_0 )
{
   GSL_vector v(static_cast<std::size_t>(0));
   BOOST_CHECK_EQUAL(v.size(), 0);
}

BOOST_AUTO_TEST_CASE( test_init )
{
   GSL_vector v(3);
   BOOST_CHECK_EQUAL(v[0], 0.);
   BOOST_CHECK_EQUAL(v[1], 0.);
   BOOST_CHECK_EQUAL(v[2], 0.);
}

BOOST_AUTO_TEST_CASE( test_operator )
{
   GSL_vector v(3);
   v[0] = 1.;
   v[1] = 2.;
   v[2] = 3.;

   BOOST_CHECK_EQUAL(v[0], 1.);
   BOOST_CHECK_EQUAL(v[1], 2.);
   BOOST_CHECK_EQUAL(v[2], 3.);
   BOOST_CHECK_EQUAL(v(0), 1.);
   BOOST_CHECK_EQUAL(v(1), 2.);
   BOOST_CHECK_EQUAL(v(2), 3.);
}

BOOST_AUTO_TEST_CASE( test_copy )
{
   GSL_vector v(3);
   v[0] = 1.;
   v[1] = 2.;
   v[2] = 3.;

   GSL_vector v2(v);

   BOOST_CHECK_EQUAL(v2[0], 1.);
   BOOST_CHECK_EQUAL(v2[1], 2.);
   BOOST_CHECK_EQUAL(v2[2], 3.);
}

BOOST_AUTO_TEST_CASE( test_assign )
{
   GSL_vector v(3);
   v[0] = 1.;
   v[1] = 2.;
   v[2] = 3.;

   GSL_vector v2;
   v2 = v;

   BOOST_CHECK_EQUAL(v2[0], 1.);
   BOOST_CHECK_EQUAL(v2[1], 2.);
   BOOST_CHECK_EQUAL(v2[2], 3.);
}

BOOST_AUTO_TEST_CASE( test_gsl_assign )
{
   gsl_vector* v;
   v = gsl_vector_alloc(3);
   gsl_vector_set(v, 0, 1.);
   gsl_vector_set(v, 1, 2.);
   gsl_vector_set(v, 2, 3.);

   GSL_vector v2(3);
   v2 = v;

   gsl_vector_free(v);

   BOOST_CHECK_EQUAL(v2[0], 1.);
   BOOST_CHECK_EQUAL(v2[1], 2.);
   BOOST_CHECK_EQUAL(v2[2], 3.);
}

BOOST_AUTO_TEST_CASE( test_gsl_empty_assign )
{
   gsl_vector* v;
   v = gsl_vector_alloc(3);
   gsl_vector_set(v, 0, 1.);
   gsl_vector_set(v, 1, 2.);
   gsl_vector_set(v, 2, 3.);

   GSL_vector v2;
   v2 = v;

   gsl_vector_free(v);

   BOOST_CHECK_EQUAL(v2[0], 1.);
   BOOST_CHECK_EQUAL(v2[1], 2.);
   BOOST_CHECK_EQUAL(v2[2], 3.);
}

BOOST_AUTO_TEST_CASE( test_gsl_null_assign )
{
   gsl_vector* v = 0;
   GSL_vector v2(3);
   v2 = v;

   BOOST_CHECK_EQUAL(v2.size(), 0);
}

BOOST_AUTO_TEST_CASE( test_assign_chain )
{
   GSL_vector v(3);
   v[0] = 1.;
   v[1] = 2.;
   v[2] = 3.;

   GSL_vector v2, v3;
   v3 = v2 = v;

   BOOST_CHECK_EQUAL(v2[0], 1.);
   BOOST_CHECK_EQUAL(v2[1], 2.);
   BOOST_CHECK_EQUAL(v2[2], 3.);
   BOOST_CHECK_EQUAL(v3[0], 1.);
   BOOST_CHECK_EQUAL(v3[1], 2.);
   BOOST_CHECK_EQUAL(v3[2], 3.);
}

BOOST_AUTO_TEST_CASE( test_move_assign )
{
   GSL_vector v(3);
   v[0] = 1.;
   v[1] = 2.;
   v[2] = 3.;

   GSL_vector v2(2);
   v2 = std::move(v);

   BOOST_CHECK(!v.raw());
   BOOST_CHECK_EQUAL(v2[0], 1.);
   BOOST_CHECK_EQUAL(v2[1], 2.);
   BOOST_CHECK_EQUAL(v2[2], 3.);
}

BOOST_AUTO_TEST_CASE( test_move_assign_to_empty )
{
   GSL_vector v(3);
   v[0] = 1.;
   v[1] = 2.;
   v[2] = 3.;

   GSL_vector v2;
   v2 = std::move(v);

   BOOST_CHECK(!v.raw());
   BOOST_CHECK_EQUAL(v2[0], 1.);
   BOOST_CHECK_EQUAL(v2[1], 2.);
   BOOST_CHECK_EQUAL(v2[2], 3.);
}

BOOST_AUTO_TEST_CASE( test_move_assign_from_empty )
{
   GSL_vector v;
   GSL_vector v2;
   v2 = std::move(v);

   BOOST_CHECK(!v.raw());
   BOOST_CHECK(!v2.raw());
}

BOOST_AUTO_TEST_CASE( test_move_assign_chain )
{
   GSL_vector v(3);
   v[0] = 1.;
   v[1] = 2.;
   v[2] = 3.;

   GSL_vector v2, v3;
   v3 = v2 = std::move(v);

   BOOST_CHECK(!v.raw());
   BOOST_CHECK_EQUAL(v2[0], 1.);
   BOOST_CHECK_EQUAL(v2[1], 2.);
   BOOST_CHECK_EQUAL(v2[2], 3.);
   BOOST_CHECK_EQUAL(v3[0], 1.);
   BOOST_CHECK_EQUAL(v3[1], 2.);
   BOOST_CHECK_EQUAL(v3[2], 3.);
}

BOOST_AUTO_TEST_CASE( test_self_assign )
{
   GSL_vector v(3);
   v[0] = 1.;
   v[1] = 2.;
   v[2] = 3.;

   v = v;

   BOOST_CHECK_EQUAL(v[0], 1.);
   BOOST_CHECK_EQUAL(v[1], 2.);
   BOOST_CHECK_EQUAL(v[2], 3.);
}

BOOST_AUTO_TEST_CASE( test_initializer_list )
{
   GSL_vector v = {1., 2., 3.};

   BOOST_CHECK_EQUAL(v.size(), 3);
   BOOST_CHECK_EQUAL(v[0], 1.);
   BOOST_CHECK_EQUAL(v[1], 2.);
   BOOST_CHECK_EQUAL(v[2], 3.);
}

BOOST_AUTO_TEST_CASE( test_empty_initializer_list )
{
   GSL_vector v = {};

   BOOST_CHECK_EQUAL(v.size(), 0);
}

BOOST_AUTO_TEST_CASE( test_empty )
{
   GSL_vector v(3), v2;

   BOOST_CHECK(!v.empty());
   BOOST_CHECK(v2.empty());
}

BOOST_AUTO_TEST_CASE( test_set_all )
{
   GSL_vector v(3);
   v.set_all(1.);

   BOOST_CHECK_EQUAL(v[0], 1.);
   BOOST_CHECK_EQUAL(v[1], 1.);
   BOOST_CHECK_EQUAL(v[2], 1.);
}

BOOST_AUTO_TEST_CASE( test_set_all_empty )
{
   GSL_vector v;
   v.set_all(1.);

   BOOST_CHECK_EQUAL(v.size(), 0);
}

BOOST_AUTO_TEST_CASE( test_move )
{
   GSL_vector v(3);
   v[0] = 1.;
   v[1] = 2.;
   v[2] = 3.;

   GSL_vector v2(std::move(v));

   BOOST_CHECK(!v.raw());
   BOOST_CHECK_EQUAL(v2[0], 1.);
   BOOST_CHECK_EQUAL(v2[1], 2.);
   BOOST_CHECK_EQUAL(v2[2], 3.);
}

BOOST_AUTO_TEST_CASE( test_output )
{
   GSL_vector v(3);
   v[0] = 1.;
   v[1] = 2.;
   v[2] = 3.;

   BOOST_TEST_MESSAGE(v);
}

BOOST_AUTO_TEST_CASE( test_empty_bounds_check )
{
   GSL_vector v;

   BOOST_CHECK_THROW(v(0), OutOfBoundsError);
   BOOST_CHECK_THROW(v(1), OutOfBoundsError);
}

BOOST_AUTO_TEST_CASE( test_nonempty_bounds_check )
{
   GSL_vector v(3);

   BOOST_CHECK_NO_THROW(v(0));
   BOOST_CHECK_NO_THROW(v(1));
   BOOST_CHECK_NO_THROW(v(2));
   BOOST_CHECK_THROW(v(3), OutOfBoundsError);
}

BOOST_AUTO_TEST_CASE( test_alloc )
{
   gsl_set_error_handler_off();
   BOOST_CHECK_THROW(GSL_vector v(std::numeric_limits<std::size_t>::max()),
                     OutOfMemoryError);
}

BOOST_AUTO_TEST_CASE( test_begin )
{
   GSL_vector v(3);
   v[0] = 1.;
   v[1] = 2.;
   v[2] = 3.;

   BOOST_REQUIRE(begin(v));
   BOOST_CHECK_EQUAL(*begin(v), 1.);
}

BOOST_AUTO_TEST_CASE( test_end )
{
   GSL_vector v(3);
   v[0] = 1.;
   v[1] = 2.;
   v[2] = 3.;

   const double* e = begin(v);
   e++; e++; e++;

   BOOST_REQUIRE(end(v));
   BOOST_CHECK_EQUAL(e, end(v));
}

BOOST_AUTO_TEST_CASE( test_cbegin )
{
   GSL_vector v(3);
   v[0] = 1.;
   v[1] = 2.;
   v[2] = 3.;

   BOOST_REQUIRE(cbegin(v));
   BOOST_CHECK_EQUAL(*cbegin(v), 1.);
}

BOOST_AUTO_TEST_CASE( test_cend )
{
   GSL_vector v(3);
   v[0] = 1.;
   v[1] = 2.;
   v[2] = 3.;

   const double* e = begin(v);
   e++; e++; e++;

   BOOST_REQUIRE(cend(v));
   BOOST_CHECK_EQUAL(e, cend(v));
}

BOOST_AUTO_TEST_CASE( test_begin_end_empty )
{
   GSL_vector v;
   BOOST_CHECK(begin(v) == end(v));
}

BOOST_AUTO_TEST_CASE( test_range_for )
{
   GSL_vector v(3);
   v[0] = 1.;
   v[1] = 2.;
   v[2] = 3.;

   double sum = 0;

   for (auto& x: v)
      sum += x;

   BOOST_CHECK_EQUAL(sum, 6.);
}

BOOST_AUTO_TEST_CASE( test_release )
{
   GSL_vector v(3);
   v[0] = 1.;
   v[1] = 2.;
   v[2] = 3.;

   gsl::owner<gsl_vector>* vp = v.release();

   BOOST_CHECK(v.raw() == nullptr);
   BOOST_CHECK(vp != nullptr);

   gsl_vector_free(vp);
}
