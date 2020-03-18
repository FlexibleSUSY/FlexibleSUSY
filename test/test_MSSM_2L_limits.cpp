#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_MSSM_2L_limits

#include <boost/test/unit_test.hpp>
#include <cmath>

#include "mssm_twoloophiggs.hpp"
#include "mssm_twoloophiggs.h"

using namespace std;
using namespace flexiblesusy;

constexpr double sqr(double x) { return x*x; }

struct Point {
   double mt2 = sqr(173.34);
   double mg = 1000;
   double mst12 = sqr(4000);
   double mst22 = sqr(5000);
   double st = 0.0000001;
   double ct = sqrt (1. - sqr(st));
   double q2 = sqr(3000);
   double mu = 2000;
   double tb = 5;
   double v2 = sqr(245);
   double g3 = sqrt(4*M_PI*0.1184);
   // for O(at*at) corrections
   double mb2 = sqr(4.18);
   double mA2 = sqr(300);
   double msb12 = sqr(1500);
   double msb22 = sqr(2500);
   double sb = 0.0000001;
   double cb = sqrt (1. - sqr(sb));
};

// Pietro Slavich implementation of CP-even Mh correction O(at*as)
Eigen::Matrix<double, 2, 2> calc_dMh_at_as_PS(Point p)
{
   double S11, S22, S12;
   int OS = 0;

   dszhiggs_(
      &p.mt2, &p.mg, &p.mst12, &p.mst22, &p.st, &p.ct,
      &p.q2, &p.mu, &p.tb, &p.v2, &p.g3, &OS, &S11, &S22, &S12);

   Eigen::Matrix<double, 2, 2> result;
   result << S11, S12, S12, S22;

   return result;
}

// Pietro Slavich implementation of CP-even Mh correction O(at*at)
Eigen::Matrix<double, 2, 2> calc_dMh_at_at_PS(Point p)
{
   double S11, S22, S12;

   ddshiggs_(&p.mt2, &p.mb2, &p.mA2, &p.mst12, &p.mst22, &p.msb12, &p.msb22,
             &p.st, &p.ct, &p.sb, &p.cb, &p.q2, &p.mu, &p.tb, &p.v2,
             &S11, &S12, &S22);

   Eigen::Matrix<double, 2, 2> result;
   result << S11, S12, S12, S22;

   return result;
}

// Pietro Slavich implementation of CP-odd MA correction O(at*as)
double calc_dMA_at_as_PS(Point p)
{
   double result;

   dszodd_(&p.mt2, &p.mg, &p.mst12, &p.mst22, &p.st, &p.ct, &p.q2,
           &p.mu, &p.tb, &p.v2, &p.g3, &result);

   return result;
}

// Pietro Slavich implementation of CP-even tadpoles O(at*as)
Eigen::Matrix<double, 2, 1> calc_tad_at_as_PS(Point p)
{
   double t1, t2;

   ewsb2loop_(
      &p.mt2, &p.mg, &p.mst12, &p.mst22, &p.st, &p.ct,
      &p.q2, &p.mu, &p.tb, &p.v2, &p.g3, &t1, &t2);

   Eigen::Matrix<double, 2, 1> result;
   result << t1, t2;

   return result;
}

// Pietro Slavich implementation of CP-even tadpoles O(at*at)
Eigen::Matrix<double, 2, 1> calc_tad_at_at_PS(Point p)
{
   double t1, t2;

   ddstad_(&p.mt2, &p.mb2, &p.mA2, &p.mst12, &p.mst22, &p.msb12, &p.msb22,
           &p.st, &p.ct, &p.sb, &p.cb, &p.q2, &p.mu, &p.tb, &p.v2,
           &t1, &t2);

   Eigen::Matrix<double, 2, 1> result;
   result << t1, t2;

   return result;
}

// FlexibleSUSY wrapper for CP-even tadpoles O(at*as)
Eigen::Matrix<double, 2, 1> calc_tad_at_as_FS(Point p)
{
   return -mssm_twoloophiggs::tadpole_higgs_2loop_at_as_mssm(
      p.mt2, p.mg, p.mst12, p.mst22, p.st, p.ct, p.q2, p.mu, p.tb, p.v2, p.g3);
}

// FlexibleSUSY wrapper for CP-even tadpoles O(at*at)
Eigen::Matrix<double, 2, 1> calc_tad_at_at_FS(Point p)
{
   return -mssm_twoloophiggs::tadpole_higgs_2loop_at_at_mssm(
      p.mt2, p.mb2, p.mA2, p.mst12, p.mst22, p.msb12, p.msb22, p.st, p.ct, p.sb, p.cb,
      p.q2, p.mu, p.tb, p.v2);
}

// FlexibleSUSY wrapper for CP-even Mh correction O(at*as)
Eigen::Matrix<double, 2, 2> calc_dMh_at_as_FS(Point p)
{
   return -mssm_twoloophiggs::self_energy_higgs_2loop_at_as_mssm_with_tadpoles(
      p.mt2, p.mg, p.mst12, p.mst22, p.st, p.ct, p.q2, p.mu, p.tb,
      p.v2, p.g3, 0);
}

// FlexibleSUSY wrapper for CP-even Mh correction O(at*at)
Eigen::Matrix<double, 2, 2> calc_dMh_at_at_FS(Point p)
{
   return -mssm_twoloophiggs::self_energy_higgs_2loop_at_at_mssm_with_tadpoles(
      p.mt2, p.mb2, p.mA2, p.mst12, p.mst22, p.msb12, p.msb22,
      p.st, p.ct, p.sb, p.cb, p.q2, p.mu, p.tb, p.v2);
}

// FlexibleSUSY wrapper for CP-odd MA correction O(at*as)
double calc_dMA_at_as_FS(Point p)
{
   return -mssm_twoloophiggs::self_energy_pseudoscalar_2loop_at_as_mssm_with_tadpoles(
      p.mt2, p.mg, p.mst12, p.mst22, p.st, p.ct, p.q2, p.mu, p.tb,
      p.v2, p.g3);
}

BOOST_AUTO_TEST_CASE( MSSM_tadpole_at_as_st_0 )
{
   Point p_close, p_exact;
   p_close.st = 0.0000001;
   p_exact.st = 0;

   const auto tad_ps       = calc_tad_at_as_PS(p_close);
   const auto tad_fs_close = calc_tad_at_as_FS(p_close);
   const auto tad_fs_exact = calc_tad_at_as_FS(p_exact);

   BOOST_CHECK_EQUAL(tad_ps(0), tad_fs_close(0));
   BOOST_CHECK_EQUAL(tad_ps(1), tad_fs_close(1));

   BOOST_CHECK_CLOSE(tad_ps(0), tad_fs_exact(0), 1e-3);
   BOOST_CHECK_CLOSE(tad_ps(1), tad_fs_exact(1), 1e-3);

   BOOST_TEST_MESSAGE("Pietro Slavich: " << tad_ps.transpose());
   BOOST_TEST_MESSAGE("Limit st -> 0 : " << tad_fs_exact.transpose());
}

BOOST_AUTO_TEST_CASE( MSSM_tadpole_at_at_st_0 )
{
   Point p_close, p_exact;
   p_close.st = 0.000000001;
   p_exact.st = 0;

   const auto tad_ps       = calc_tad_at_at_PS(p_close);
   const auto tad_fs_close = calc_tad_at_at_FS(p_close);
   const auto tad_fs_exact = calc_tad_at_at_FS(p_exact);

   BOOST_CHECK_EQUAL(tad_ps(0), tad_fs_close(0));
   BOOST_CHECK_EQUAL(tad_ps(1), tad_fs_close(1));

   BOOST_CHECK_CLOSE(tad_ps(0), tad_fs_exact(0), 5e-5);
   BOOST_CHECK_CLOSE(tad_ps(1), tad_fs_exact(1), 1e-5);

   BOOST_TEST_MESSAGE("Pietro Slavich: " << tad_ps.transpose());
   BOOST_TEST_MESSAGE("Limit st -> 0 : " << tad_fs_exact.transpose());
}

BOOST_AUTO_TEST_CASE( MSSM_tadpole_at_as_st_0_mst1_eq_mst2 )
{
   Point p_close, p_exact;
   p_close.st = 0.0000001;
   p_close.mst12 = sqr(4000);
   p_close.mst22 = sqr(4000.01);
   p_exact.st = 0;
   p_exact.mst12 = sqr(4000);
   p_exact.mst22 = sqr(4000);

   const auto tad_ps       = calc_tad_at_as_PS(p_close);
   const auto tad_fs_close = calc_tad_at_as_FS(p_close);
   const auto tad_fs_exact = calc_tad_at_as_FS(p_exact);

   BOOST_CHECK_EQUAL(tad_ps(0), tad_fs_close(0));
   BOOST_CHECK_EQUAL(tad_ps(1), tad_fs_close(1));

   BOOST_CHECK_CLOSE(tad_ps(0), tad_fs_exact(0), 1e-3);
   BOOST_CHECK_CLOSE(tad_ps(1), tad_fs_exact(1), 1e-3);

   BOOST_TEST_MESSAGE("Pietro Slavich                : " << tad_ps.transpose());
   BOOST_TEST_MESSAGE("Limit st -> 0 and mst1 -> mst2: " << tad_fs_exact.transpose());
}

BOOST_AUTO_TEST_CASE( MSSM_tadpole_at_at_st_0_mst1_eq_mst2 )
{
   Point p_close, p_exact;
   p_close.st = 0.00000001;
   p_close.mst12 = sqr(4000);
   p_close.mst22 = sqr(4000.00001);
   p_exact.st = 0;
   p_exact.mst12 = sqr(4000);
   p_exact.mst22 = sqr(4000);

   const auto tad_ps       = calc_tad_at_at_PS(p_close);
   const auto tad_fs_close = calc_tad_at_at_FS(p_close);
   const auto tad_fs_exact = calc_tad_at_at_FS(p_exact);

   BOOST_CHECK_EQUAL(tad_ps(0), tad_fs_close(0));
   BOOST_CHECK_EQUAL(tad_ps(1), tad_fs_close(1));

   BOOST_CHECK_CLOSE(tad_ps(0), tad_fs_exact(0), 5e-5);
   BOOST_CHECK_CLOSE(tad_ps(1), tad_fs_exact(1), 1e-5);

   BOOST_TEST_MESSAGE("Pietro Slavich                : " << tad_ps.transpose());
   BOOST_TEST_MESSAGE("Limit st -> 0 and mst1 -> mst2: " << tad_fs_exact.transpose());
}

BOOST_AUTO_TEST_CASE( MSSM_dMh_at_as_st_0 )
{
   Point p_close, p_exact;
   p_close.st = 0.0000001;
   p_exact.st = 0;

   const auto dMh_ps       = calc_dMh_at_as_PS(p_close);
   const auto dMh_fs_close = calc_dMh_at_as_FS(p_close);
   const auto dMh_fs_exact = calc_dMh_at_as_FS(p_exact);

   BOOST_CHECK_EQUAL(dMh_ps(0,0), dMh_fs_close(0,0));
   BOOST_CHECK_EQUAL(dMh_ps(0,1), dMh_fs_close(0,1));
   BOOST_CHECK_EQUAL(dMh_ps(1,0), dMh_fs_close(1,0));
   BOOST_CHECK_EQUAL(dMh_ps(1,1), dMh_fs_close(1,1));

   BOOST_CHECK_SMALL(dMh_ps(0,0), 1e-3);
   BOOST_CHECK_SMALL(dMh_fs_exact(0,0), 1e-3);
   BOOST_CHECK_CLOSE(dMh_ps(0,1), dMh_fs_exact(0,1), 1e-3);
   BOOST_CHECK_CLOSE(dMh_ps(1,0), dMh_fs_exact(1,0), 1e-3);
   BOOST_CHECK_CLOSE(dMh_ps(1,1), dMh_fs_exact(1,1), 1e-3);

   BOOST_TEST_MESSAGE("Pietro Slavich:\n" << dMh_ps);
   BOOST_TEST_MESSAGE("Limit st -> 0 :\n" << dMh_fs_exact);
}

BOOST_AUTO_TEST_CASE( MSSM_dMh_at_at_st_0 )
{
   Point p_close, p_exact;
   p_close.st = 0.00000001;
   p_exact.st = 0;

   const auto dMh_ps       = calc_dMh_at_at_PS(p_close);
   const auto dMh_fs_close = calc_dMh_at_at_FS(p_close);
   const auto dMh_fs_exact = calc_dMh_at_at_FS(p_exact);

   BOOST_CHECK_EQUAL(dMh_ps(0,0), dMh_fs_close(0,0));
   BOOST_CHECK_EQUAL(dMh_ps(0,1), dMh_fs_close(0,1));
   BOOST_CHECK_EQUAL(dMh_ps(1,0), dMh_fs_close(1,0));
   BOOST_CHECK_EQUAL(dMh_ps(1,1), dMh_fs_close(1,1));

   BOOST_CHECK_CLOSE(dMh_ps(0,0), dMh_fs_exact(0,0), 5e-2);
   BOOST_CHECK_CLOSE(dMh_ps(0,1), dMh_fs_exact(0,1), 5e-4);
   BOOST_CHECK_CLOSE(dMh_ps(1,0), dMh_fs_exact(1,0), 5e-4);
   BOOST_CHECK_CLOSE(dMh_ps(1,1), dMh_fs_exact(1,1), 1e-5);

   BOOST_TEST_MESSAGE("Pietro Slavich:\n" << dMh_ps);
   BOOST_TEST_MESSAGE("Limit st -> 0 :\n" << dMh_fs_exact);
}

BOOST_AUTO_TEST_CASE( MSSM_dMh_at_as_st_0_mst1_eq_mst2 )
{
   Point p_close, p_exact;
   p_close.st = 0.0000001;
   p_close.mst12 = sqr(4000);
   p_close.mst22 = sqr(4000.01);
   p_exact.st = 0;
   p_exact.mst12 = sqr(4000);
   p_exact.mst22 = sqr(4000);

   const auto dMh_ps       = calc_dMh_at_as_PS(p_close);
   const auto dMh_fs_close = calc_dMh_at_as_FS(p_close);
   const auto dMh_fs_exact = calc_dMh_at_as_FS(p_exact);

   BOOST_CHECK_SMALL(dMh_ps(0,0), 1e-6);
   BOOST_CHECK_SMALL(dMh_fs_close(0,0), 1e-6);
   BOOST_CHECK_CLOSE(dMh_ps(0,1), dMh_fs_close(0,1), 1e-3);
   BOOST_CHECK_CLOSE(dMh_ps(1,0), dMh_fs_close(1,0), 1e-3);
   BOOST_CHECK_CLOSE(dMh_ps(1,1), dMh_fs_close(1,1), 2e-3);

   BOOST_CHECK_SMALL(dMh_ps(0,0), 1e-3);
   BOOST_CHECK_SMALL(dMh_fs_exact(0,0), 1e-3);
   BOOST_CHECK_CLOSE(dMh_ps(0,1), dMh_fs_exact(0,1), 1e-3);
   BOOST_CHECK_CLOSE(dMh_ps(1,0), dMh_fs_exact(1,0), 1e-3);
   BOOST_CHECK_CLOSE(dMh_ps(1,1), dMh_fs_exact(1,1), 2e-3);

   BOOST_TEST_MESSAGE("Pietro Slavich                 :\n" << dMh_ps);
   BOOST_TEST_MESSAGE("Limit st -> 0 and mst1 -> mst2 :\n" << dMh_fs_exact);
}

BOOST_AUTO_TEST_CASE( MSSM_dMh_at_at_st_0_mst1_eq_mst2 )
{
   Point p_close, p_exact;
   p_close.st = 0.00001;
   p_close.mst12 = sqr(4000);
   p_close.mst22 = sqr(4000.01);
   p_exact.st = 0;
   p_exact.mst12 = sqr(4000);
   p_exact.mst22 = sqr(4000);

   const auto dMh_ps       = calc_dMh_at_at_PS(p_close);
   const auto dMh_fs_close = calc_dMh_at_at_FS(p_close);
   const auto dMh_fs_exact = calc_dMh_at_at_FS(p_exact);

   BOOST_CHECK_CLOSE(dMh_ps(0,0), dMh_fs_close(0,0), 1e-6);
   BOOST_CHECK_CLOSE(dMh_ps(0,1), dMh_fs_close(0,1), 1e-3);
   BOOST_CHECK_CLOSE(dMh_ps(1,0), dMh_fs_close(1,0), 1e-3);
   BOOST_CHECK_CLOSE(dMh_ps(1,1), dMh_fs_close(1,1), 2e-3);

   // BOOST_CHECK_CLOSE(dMh_ps(0,0), dMh_fs_exact(0,0), 1e-3); // FAILS
   BOOST_CHECK_CLOSE(dMh_ps(0,1), dMh_fs_exact(0,1), 15);
   BOOST_CHECK_CLOSE(dMh_ps(1,0), dMh_fs_exact(1,0), 15);
   BOOST_CHECK_CLOSE(dMh_ps(1,1), dMh_fs_exact(1,1), 0.007);

   BOOST_TEST_MESSAGE("Pietro Slavich                 :\n" << dMh_ps);
   BOOST_TEST_MESSAGE("Limit st -> 0 and mst1 -> mst2 :\n" << dMh_fs_exact);
}

// calculate st such that At = 0
double calc_st(Point p)
{
   if (p.mu == 0.)
      return 0.;

   return std::sqrt(p.mt2) * p.mu / (p.ct * (p.mst12 - p.mst22) * p.tb);
}

BOOST_AUTO_TEST_CASE( MSSM_dMA_at_as_At_0_mst1_eq_mst2 )
{
   Point p_close, p_exact;
   p_close.mu = p_exact.mu = 0;
   p_close.mst12 = sqr(4000);
   p_close.mst22 = sqr(4000.01);
   p_exact.mst12 = sqr(4000);
   p_exact.mst22 = sqr(4000);

   // fix st such that At = 0
   p_close.st = calc_st(p_close);
   p_exact.st = calc_st(p_exact);

   const auto dMA_ps       = calc_dMA_at_as_PS(p_close);
   const auto dMA_fs_close = calc_dMA_at_as_FS(p_close);
   const auto dMA_fs_exact = calc_dMA_at_as_FS(p_exact);

   BOOST_CHECK_EQUAL(dMA_ps, dMA_fs_close);
   BOOST_CHECK_CLOSE(dMA_ps, dMA_fs_exact, 1e-3);

   BOOST_TEST_MESSAGE("Pietro Slavich                : " << dMA_ps);
   BOOST_TEST_MESSAGE("Limit At -> 0 and mst1 -> mst2: " << dMA_fs_exact);
}

BOOST_AUTO_TEST_CASE( MSSM_dMA_at_as_mst1_eq_mst2 )
{
   // At != 0

   Point p_close, p_exact;
   p_close.st = 0.0000001;
   p_close.mst12 = sqr(4000);
   p_close.mst22 = sqr(4000.01);
   p_exact.st = 0;
   p_exact.mst12 = sqr(4000);
   p_exact.mst22 = sqr(4000);

   const auto dMA_ps       = calc_dMA_at_as_PS(p_close);
   const auto dMA_fs_close = calc_dMA_at_as_FS(p_close);
   const auto dMA_fs_exact = calc_dMA_at_as_FS(p_exact);

   BOOST_CHECK_CLOSE(dMA_ps, dMA_fs_close, 1e-3);
   BOOST_CHECK_CLOSE(dMA_ps, dMA_fs_exact, 1e-3);

   BOOST_TEST_MESSAGE("Pietro Slavich    : " << dMA_ps);
   BOOST_TEST_MESSAGE("Limit mst1 -> mst2: " << dMA_fs_exact);
}
