
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_betafunction

#include <boost/test/unit_test.hpp>

#include "rge.h"
#include "betafunction.hpp"
#include "wrappers.hpp"
#include "conversion.hpp"
#include "error.hpp"
#include "logger.hpp"

using namespace flexiblesusy;
using namespace softsusy;

class RGE_model : public RGE {
public:
   RGE_model() : RGE(), pars(10) {
      setPars(10);
      setMu(100.);
   }
   virtual ~RGE_model() {}
   virtual const DoubleVector display() const { return pars; }
   virtual void set(const DoubleVector& s) { pars = s; }
   virtual DoubleVector beta() const {
      DoubleVector beta(pars.size());
      for (std::size_t i = 1; i <= pars.size(); i++)
         beta(i) = pars(i) * pars(i) * 0.1 + 0.05 * pars(pars.size() + 1 - i);
      return beta;
   }
private:
   DoubleVector pars;
};

class Eigen_model : public Beta_function {
public:
   Eigen_model() : Beta_function(), pars(Eigen::ArrayXd::Zero(10)) {
      set_scale(100.);
      set_number_of_parameters(10);
      set_loops(1);
   }
   virtual ~Eigen_model() {}
   virtual Eigen::ArrayXd get() const { return pars; }
   virtual void set(const Eigen::ArrayXd& s) { pars = s; }
   virtual Eigen::ArrayXd beta() const {
      Eigen::ArrayXd beta(pars.size());
      for (int i = 0; i < pars.rows(); i++)
         beta(i) = pars(i) * pars(i) * 0.1 + 0.05 * pars(pars.rows() - 1 - i);
      return beta;
   }
private:
   Eigen::ArrayXd pars;
};

namespace
{

/// compares element-wise; after running, the two implementations may
/// differ by rounding (e.g. with fused multiply-adds)
void check_close(const DoubleVector& a, const DoubleVector& b, double eps)
{
   BOOST_REQUIRE_EQUAL(a.size(), b.size());
   for (std::size_t i = 1; i <= a.size(); i++) {
      BOOST_CHECK_CLOSE_FRACTION(a(i), b(i), eps);
   }
}

} // anonymous namespace

BOOST_AUTO_TEST_CASE( test_running )
{
   DoubleVector pars(10);
   for (std::size_t i = 1; i <= pars.size(); i++)
      pars(i) = 0.01 * i;

   RGE_model rge;
   Eigen_model eig;

   rge.set(pars);
   eig.set(ToEigenArray(pars));

   BOOST_CHECK_EQUAL(rge.display(), ToDoubleVector(eig.get()));
   BOOST_CHECK_EQUAL(rge.displayMu(), eig.get_scale());
   BOOST_CHECK_EQUAL(rge.howMany(), eig.get_number_of_parameters());

   int status_rge = rge.runto(1.0e10);
   int status_beta = 0;

   try {
      eig.run_to(1.0e10);
   } catch (Error& error) {
      status_beta = 1;
   }

   BOOST_CHECK_EQUAL(status_rge, status_beta);
   check_close(rge.display(), ToDoubleVector(eig.get()), 1e-12);
   BOOST_CHECK_EQUAL(rge.displayMu(), eig.get_scale());
}

BOOST_AUTO_TEST_CASE( test_custom_precision )
{
   DoubleVector pars(10);
   for (std::size_t i = 1; i <= pars.size(); i++)
      pars(i) = 0.01 * i;

   RGE_model rge;
   Eigen_model eig;

   rge.set(pars);
   eig.set(ToEigenArray(pars));

   BOOST_CHECK_EQUAL(rge.display(), ToDoubleVector(eig.get()));
   BOOST_CHECK_EQUAL(rge.displayMu(), eig.get_scale());
   BOOST_CHECK_EQUAL(rge.howMany(), eig.get_number_of_parameters());

   const double custom_precision = 0.02;
   rge.runto(1.0e10, custom_precision);
   eig.run_to(1.0e10, custom_precision);

   check_close(rge.display(), ToDoubleVector(eig.get()), 1e-12);
   BOOST_CHECK_EQUAL(rge.displayMu(), eig.get_scale());
}
