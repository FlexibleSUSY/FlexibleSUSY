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

#include "weinberg_angle.hpp"
#include "config.h"
#include "ew_input.hpp"
#include "logger.hpp"
#include "numerics.h"
#include "numerics2.hpp"

#include <cmath>
#include <complex>
#include <limits>

#define WARN_IF_ZERO(p,fun)                     \
   if (is_zero(p))                              \
      WARNING(#fun ": " #p " is zero!");

#define QUIT_IF(condition,fun)                  \
   if (condition)                               \
      FATAL(#fun ": " #condition);

#define ROOT2 Electroweak_constants::root2

namespace flexiblesusy {

namespace {

constexpr double Pi = 3.141592653589793;

constexpr double oneOver16PiSqr = 1.0/(16.0*Pi*Pi);

double sqr(double x) noexcept { return x*x; }

} // anonymous namespace

namespace weinberg_angle {

using namespace softsusy;

Weinberg_angle::Data::Data()
   : scale(0.)
   , alpha_em_drbar(0.)
   , fermi_contant(0.)
   , self_energy_z_at_mz(0.)
   , self_energy_w_at_0(0.)
   , self_energy_w_at_mw(0.)
   , mw_pole(0.)
   , mz_pole(0.)
   , mt_pole(0.)
   , mh_drbar(0.)
   , hmix_12(0.)
   , msel_drbar(0.)
   , msmul_drbar(0.)
   , msve_drbar(0.)
   , msvm_drbar(0.)
   , mn_drbar()
   , mc_drbar()
   , zn()
   , um()
   , up()
   , gY(0.)
   , g2(0.)
   , g3(0.)
   , tan_beta(0.)
{
}

Weinberg_angle::Self_energy_data::Self_energy_data()
   : scale(0.)
   , mt_pole(0.)
   , mt_drbar(0.)
   , mb_drbar(0.)
   , gY(0.)
   , g2(0.)
{
}

/**
 * Sets the maximum number of iterations to 20 and
 * sets the precision goal to 1.0e-8 and
 * sets enables the SUSY contributions by default.
 */
Weinberg_angle::Weinberg_angle()
   : number_of_iterations(20)
   , number_of_loops(2)
   , precision_goal(1.0e-8)
   , rho_hat(0.)
   , sin_theta(0.)
   , data()
   , susy_contributions(true)
{
}

void Weinberg_angle::enable_susy_contributions()
{
   susy_contributions = true;
}

void Weinberg_angle::disable_susy_contributions()
{
   susy_contributions = false;
}

void Weinberg_angle::set_number_of_iterations(int n)
{
   number_of_iterations = n;
}

void Weinberg_angle::set_number_of_loops(int n)
{
   number_of_loops = n;
}

void Weinberg_angle::set_data(const Data& data_)
{
   data = data_;
}

void Weinberg_angle::set_precision_goal(double p)
{
   precision_goal = p;
}

double Weinberg_angle::get_rho_hat() const
{
   return rho_hat;
}

double Weinberg_angle::get_sin_theta() const
{
   return sin_theta;
}

/**
 * Calculates the DR-bar weak mixing angle \f$\sin\hat{\theta}_W\f$ as
 * defined in Eq. (C.3) from hep-ph/9606211 given the Fermi constant,
 * the Z-boson pole mass and the DR-bar electromagnetic coupling as
 * input.
 *
 * The function returns 1 if the iterative procedure to determine the
 * weak mixing angle does not converge.
 *
 * @param rho_start initial guess for the rho-hat-parameter
 * @param sin_start initial guess for the sinus of the weak mixing angle
 *
 * @return value != 0 if an error has occured
 */
int Weinberg_angle::calculate(double rho_start, double sin_start)
{
   const double alphaDRbar = data.alpha_em_drbar;
   const double mz_pole    = data.mz_pole;
   const double gfermi     = data.fermi_contant;

   int iteration = 0;
   bool not_converged = true;
   double rho_old = rho_start, sin_old = sin_start;
   double rho_new = rho_start, sin_new = sin_start;

   while (not_converged && iteration < number_of_iterations) {
      double deltaR
         = calculate_delta_r(rho_old, sin_old, data, susy_contributions,
                             number_of_loops);

      if (deltaR > 1.) {
#if defined(ENABLE_VERBOSE) || defined(ENABLE_DEBUG)
         WARNING("delta_r > 1");
#endif
         deltaR = 0.;
      }

      if (!std::isfinite(deltaR)) {
#if defined(ENABLE_VERBOSE) || defined(ENABLE_DEBUG)
         WARNING("delta_r non-finite");
#endif
         deltaR = 0.;
      }

      double sin2thetasqO4 = Pi * alphaDRbar /
         (ROOT2 * sqr(mz_pole) * gfermi * (1.0 - deltaR));

      if (sin2thetasqO4 >= 0.25)
         sin2thetasqO4 = 0.25;

      if (sin2thetasqO4 < 0.0)
         sin2thetasqO4 = 0.0;

      const double sin2theta = std::sqrt(4.0 * sin2thetasqO4);
      const double theta = 0.5 * std::asin(sin2theta);

      sin_new = std::sin(theta);

      double deltaRho
         = calculate_delta_rho(rho_old, sin_new, data, susy_contributions,
                               number_of_loops);

      if (!std::isfinite(deltaRho)) {
#if defined(ENABLE_VERBOSE) || defined(ENABLE_DEBUG)
         WARNING("delta_rho non-finite");
#endif
         deltaRho = 0.;
      }

      if (std::abs(deltaRho) < 1.0)
         rho_new = 1.0 / (1.0 - deltaRho);
      else
         rho_new = 1.0;

      const double precision
         = std::abs(rho_old / rho_new - 1.0) + std::abs(sin_old / sin_new - 1.0);

      VERBOSE_MSG("\t\tIteration step " << iteration
                  << ": prec=" << precision
                  << " dr=" << deltaR
                  << " drho=" << deltaRho
                  << " rho_new=" << rho_new
                  << " sin_new=" << sin_new);

      not_converged = precision >= precision_goal;

      rho_old = rho_new;
      sin_old = sin_new;
      iteration++;
   }

   rho_hat = rho_new;
   sin_theta = sin_new;

   const int no_convergence_error = iteration == number_of_iterations;

   return no_convergence_error;
}

/**
 * Calculates the \f$\Delta\hat{\rho}\f$ corrections as defined in
 * Eqs. (C.4), (C.6) from hep-ph/9606211 .
 *
 * @param rho rho-hat-parameter
 * @param sinThetaW sin(theta_W)
 * @param data data structure with model parameters
 * @param add_susy_contributions add SUSY particle contributions
 * @param number_of_loops loop order (0, 1, 2)
 *
 * @return \f$\Delta\hat{\rho}\f$ as defined in (C.5) and (C.5) from hep-ph/9606211
 */
double Weinberg_angle::calculate_delta_rho(
   double rho,
   double sinThetaW,
   const Data& data,
   bool add_susy_contributions,
   int number_of_loops
)
{
   const double mz = data.mz_pole;
   const double mw = data.mw_pole;
   const double mt = data.mt_pole;
   const double mh = data.mh_drbar;
   const double sinb = std::sin(std::atan(data.tan_beta));
   const double xt = 3.0 * data.fermi_contant * sqr(mt) * ROOT2 * oneOver16PiSqr;
   const double alphaDRbar = data.alpha_em_drbar;
   const double g3 = data.g3;
   const double pizztMZ = data.self_energy_z_at_mz;
   const double piwwtMW = data.self_energy_w_at_mw;
   const double hmix12 = data.hmix_12;

#if defined(ENABLE_VERBOSE) || defined(ENABLE_DEBUG)
   WARN_IF_ZERO(rho, calculate_delta_rho)
   WARN_IF_ZERO(sinThetaW, calculate_delta_rho)
   WARN_IF_ZERO(mz, calculate_delta_rho)
   WARN_IF_ZERO(mw, calculate_delta_rho)
   WARN_IF_ZERO(mt, calculate_delta_rho)
   WARN_IF_ZERO(mh, calculate_delta_rho)
   WARN_IF_ZERO(xt, calculate_delta_rho)
   WARN_IF_ZERO(alphaDRbar, calculate_delta_rho)
   WARN_IF_ZERO(pizztMZ, calculate_delta_rho)
   WARN_IF_ZERO(piwwtMW, calculate_delta_rho)
   if (add_susy_contributions) {
      WARN_IF_ZERO(sinb, calculate_delta_rho)
      WARN_IF_ZERO(hmix12, calculate_delta_rho)
   }
#endif

   double hmix_r = 1.0;
   if (add_susy_contributions)
      hmix_r = sqr(hmix12 / sinb);

   double deltaRhoOneLoop = 0.;

   if (number_of_loops > 0)
      deltaRhoOneLoop = pizztMZ / (rho * sqr(mz)) - piwwtMW / sqr(mw);

   double deltaRho2LoopSm = 0.;

   if (number_of_loops > 1) {
      deltaRho2LoopSm = alphaDRbar * sqr(g3) /
         (16.0 * Pi * sqr(Pi) * sqr(sinThetaW)) *
         (-2.145 * sqr(mt) / sqr(mw) + 1.262 * log(mt / mz) - 2.24
          - 0.85 * sqr(mz)
          / sqr(mt)) + sqr(xt) * hmix_r *
         rho_2(mh / mt) / 3.0;
   }

   const double deltaRho = deltaRhoOneLoop + deltaRho2LoopSm;

   return deltaRho;
}

/**
 * Calculates the \f$\Delta\hat{r}\f$ corrections as defined in Eqs. (C.3),
 * (C.5) from hep-ph/9606211 .
 *
 * @param rho rho-hat-parameter
 * @param sinThetaW sin(theta_W)
 * @param data data structure with model parameters
 * @param add_susy_contributions add SUSY particle contributions
 * @param number_of_loops loop order (0, 1, 2)
 *
 * @return \f$\Delta\hat{r}\f$ as defined in (C.5) and (C.5) from hep-ph/9606211
 */
double Weinberg_angle::calculate_delta_r(
   double rho,
   double sinThetaW,
   const Data& data,
   bool add_susy_contributions,
   int number_of_loops
)
{
   const double outcos = std::cos(std::asin(sinThetaW));
   const double sinb = std::sin(std::atan(data.tan_beta));
   const double mz = data.mz_pole;
   const double mw = data.mw_pole;
   const double mt = data.mt_pole;
   const double mh = data.mh_drbar;
   const double xt = 3.0 * data.fermi_contant * sqr(mt) * ROOT2 * oneOver16PiSqr;
   const double alphaDRbar = data.alpha_em_drbar;
   const double g3 = data.g3;
   const double pizztMZ = data.self_energy_z_at_mz;
   const double piwwt0 = data.self_energy_w_at_0;
   const double hmix12 = data.hmix_12;

#if defined(ENABLE_VERBOSE) || defined(ENABLE_DEBUG)
   WARN_IF_ZERO(rho, calculate_delta_r)
   WARN_IF_ZERO(sinThetaW, calculate_delta_r)
   WARN_IF_ZERO(mz, calculate_delta_r)
   WARN_IF_ZERO(mw, calculate_delta_r)
   WARN_IF_ZERO(mt, calculate_delta_r)
   WARN_IF_ZERO(mh, calculate_delta_r)
   WARN_IF_ZERO(xt, calculate_delta_r)
   WARN_IF_ZERO(alphaDRbar, calculate_delta_r)
   WARN_IF_ZERO(g3, calculate_delta_r)
   WARN_IF_ZERO(pizztMZ, calculate_delta_r)
   WARN_IF_ZERO(piwwt0, calculate_delta_r)
   if (add_susy_contributions) {
      WARN_IF_ZERO(sinb, calculate_delta_rho)
      WARN_IF_ZERO(hmix12, calculate_delta_rho)
   }
#endif

   double dvb = 0.;

   if (number_of_loops > 0)
      dvb = calculate_delta_vb(rho, sinThetaW, data, add_susy_contributions);

   double hmix_r = 1.0;
   if (add_susy_contributions)
      hmix_r = sqr(hmix12 / sinb);

   double deltaR = 0.;

   if (number_of_loops > 0)
      deltaR = rho * piwwt0 / sqr(mw) - pizztMZ / sqr(mz) + dvb;

   double deltaR2LoopSm = 0.;

   if (number_of_loops > 1) {
      deltaR2LoopSm = alphaDRbar * sqr(g3) /
         (16.0 * sqr(Pi) * Pi * sqr(sinThetaW) * sqr(outcos)) *
         (2.145 * sqr(mt) / sqr(mz) + 0.575 * log(mt / mz) - 0.224
          - 0.144 * sqr(mz) / sqr(mt)) -
         sqr(xt) * hmix_r *
         rho_2(mh / mt) * (1.0 - deltaR) * rho / 3.0;
   }

   const double deltaR_full = deltaR + deltaR2LoopSm;

   return deltaR_full;
}

/**
 * Calculates the vertex, box and external wave-function
 * renormalizations \f$\delta_{\text{VB}}\f$ as given in
 * Eqs. (C.11)-(C.16), (C.20) from hep-ph/9606211 .
 *
 * @param rho rho-hat-parameter
 * @param sinThetaW sin(theta_W)
 * @param data data structure with model parameters
 * @param add_susy_contributions add SUSY particle contributions
 *
 * @return \f$\delta_{\text{VB}}\f$ as defined in (C.11) from hep-ph/9606211
 */
double Weinberg_angle::calculate_delta_vb(
   double rho,
   double sinThetaW,
   const Data& data,
   bool add_susy_contributions)
{
   double deltaVbSusy = 0.;
   const double deltaVbSm = calculate_delta_vb_sm(rho, sinThetaW, data);

   if (add_susy_contributions)
      deltaVbSusy = calculate_delta_vb_susy(sinThetaW, data);

   const double deltaVb = deltaVbSm + deltaVbSusy;

   return deltaVb;
}

/**
 * Calculates the Standard Model vertex, box corrections
 * \f$\delta_{\text{VB}}^{\text{SM}}\f$ as given in Eqs. (C.12) from
 * hep-ph/9606211 .
 *
 * @param rho rho-hat-parameter
 * @param sinThetaW sin(theta_W)
 * @param data data structure with model parameters
 *
 * @return \f$\delta_{\text{VB}}^{\text{SM}}\f$ as defined in (C.11)
 * from hep-ph/9606211
 */
double Weinberg_angle::calculate_delta_vb_sm(
   double rho,
   double sinThetaW,
   const Data& data
)
{
  const double mz      = data.mz_pole;
  const double mw      = data.mw_pole;
  const double costh   = mw / mz;
  const double cw2     = sqr(costh);
  const double sw2     = 1.0 - cw2;
  const double sinThetaW2 = sqr(sinThetaW);
  const double outcos  = sqrt(1.0 - sinThetaW2);
  const double alphaDRbar = data.alpha_em_drbar;
  const double scale   = data.scale;

#if defined(ENABLE_VERBOSE) || defined(ENABLE_DEBUG)
   WARN_IF_ZERO(rho, calculate_delta_vb_sm)
   WARN_IF_ZERO(sinThetaW, calculate_delta_vb_sm)
   WARN_IF_ZERO(mz, calculate_delta_vb_sm)
   WARN_IF_ZERO(mw, calculate_delta_vb_sm)
   WARN_IF_ZERO(alphaDRbar, calculate_delta_vb_sm)
#endif

  const double deltaVbSm =
     rho * alphaDRbar / (4.0 * Pi * sinThetaW2) *
     (6.0 + log(cw2) / sw2 *
      (3.5 - 2.5 * sw2 - sinThetaW2 * (5.0 - 1.5 * cw2 / sqr(outcos)))
      - 4. * std::log(sqr(mz/scale)));

  return deltaVbSm;
}

/**
 * Calculates the SUSY vertex, box and external wave-function
 * renormalizations \f$\delta_{\text{VB}}^{\text{SUSY}}\f$ as given in
 * Eqs. (C.11)-(C.16), (C.20) from hep-ph/9606211 .
 *
 * @param sinThetaW sin(theta_W)
 * @param data data structure with model parameters
 *
 * @return \f$\delta_{\text{VB}}^{\text{SUSY}}\f$ as defined in (C.11)
 * from hep-ph/9606211
 */
double Weinberg_angle::calculate_delta_vb_susy(
   double sinThetaW,
   const Data& data
)
{
  const double g       = data.g2;
  const double gp      = data.gY;
  const double mz      = data.mz_pole;
  const double sinThetaW2 = sqr(sinThetaW);
  const double outcos  = sqrt(1.0 - sinThetaW2);
  const double q       = data.scale;
  const double alphaDRbar = data.alpha_em_drbar;
  const double mselL = data.msel_drbar;
  const double msmuL = data.msmul_drbar;
  const double msnue = data.msve_drbar;
  const double msnumu = data.msvm_drbar;
  const Eigen::ArrayXd& mneut(data.mn_drbar);
  const Eigen::ArrayXd& mch(data.mc_drbar);
  const Eigen::MatrixXcd& n(data.zn);
  const Eigen::MatrixXcd& u(data.um);
  const Eigen::MatrixXcd& v(data.up);

#if defined(ENABLE_VERBOSE) || defined(ENABLE_DEBUG)
   WARN_IF_ZERO(sinThetaW, calculate_delta_vb_susy)
   WARN_IF_ZERO(g, calculate_delta_vb_susy)
   WARN_IF_ZERO(gp, calculate_delta_vb_susy)
   WARN_IF_ZERO(mz, calculate_delta_vb_susy)
   WARN_IF_ZERO(q, calculate_delta_vb_susy)
   WARN_IF_ZERO(alphaDRbar, calculate_delta_vb_susy)
   WARN_IF_ZERO(mselL, calculate_delta_vb_susy)
   WARN_IF_ZERO(msmuL, calculate_delta_vb_susy)
   WARN_IF_ZERO(msnue, calculate_delta_vb_susy)
   WARN_IF_ZERO(msnumu, calculate_delta_vb_susy)
   QUIT_IF(mneut.rows() < 4, calculate_delta_vb_susy)
   QUIT_IF(mneut.cols() != 1, calculate_delta_vb_susy)
   QUIT_IF(mch.rows() < 2, calculate_delta_vb_susy)
   QUIT_IF(mch.cols() != 1, calculate_delta_vb_susy)
   QUIT_IF(mneut.rows() != n.rows(), calculate_delta_vb_susy)
   QUIT_IF(mneut.rows() != n.cols(), calculate_delta_vb_susy)
   QUIT_IF(mch.rows() != u.rows(), calculate_delta_vb_susy)
   QUIT_IF(mch.rows() != u.cols(), calculate_delta_vb_susy)
   QUIT_IF(mch.rows() != v.rows(), calculate_delta_vb_susy)
   QUIT_IF(mch.rows() != v.cols(), calculate_delta_vb_susy)
#endif

  const int dimN = mneut.rows();
  const int dimC = mch.rows();

  Eigen::VectorXd bPsi0NuNul(Eigen::VectorXd::Zero(dimN)),
     bPsicNuSell(Eigen::VectorXd::Zero(dimC));
  Eigen::VectorXd bPsi0ESell(Eigen::VectorXd::Zero(dimN)),
     aPsicESnul(Eigen::VectorXd::Zero(dimC));
  Eigen::VectorXcd bChi0NuNul, bChicNuSell;
  Eigen::VectorXcd bChi0ESell, aChicESnul;

  bPsicNuSell(0) = g;
  bPsi0NuNul(1) = ROOT2 * g * 0.5;
  bPsi0NuNul(0) = -gp / ROOT2;
  aPsicESnul(0) = g;
  bPsi0ESell(0) = -gp / ROOT2;
  bPsi0ESell(1) = -g * ROOT2 * 0.5;

  bChicNuSell = u * bPsicNuSell;
  bChi0ESell =  n * bPsi0ESell;
  bChi0NuNul = n * bPsi0NuNul;

  aChicESnul = v.conjugate() * aPsicESnul;

  double deltaZnue = 0.0, deltaZe = 0.0;
  for (int i = 0; i < dimN; i++) {
     deltaZnue += - sqr(std::abs(bChi0NuNul(i))) * b1(0.0, mneut(i), msnue, q);
     deltaZe   += - sqr(std::abs(bChi0ESell(i))) * b1(0.0, mneut(i), mselL, q);
  }
  for (int i = 0; i < dimC; i++) {
     deltaZnue += - sqr(std::abs(bChicNuSell(i))) * b1(0.0, mch(i), mselL, q);
     deltaZe   += - sqr(std::abs(aChicESnul(i))) * b1(0.0, mch(i), msnue, q);
  }

  Eigen::VectorXd bPsicNuSmul(Eigen::VectorXd::Zero(dimC));
  Eigen::VectorXd bPsi0MuSmul(Eigen::VectorXd::Zero(dimN)),
     aPsicMuSnul(Eigen::VectorXd::Zero(dimC));
  Eigen::VectorXcd bChicNuSmul;
  Eigen::VectorXcd bChi0MuSmul, aChicMuSnul;

  bPsicNuSmul(0) = g;
  aPsicMuSnul(0) = g;
  bPsi0MuSmul(0) = -gp / ROOT2;
  bPsi0MuSmul(1) = -g * ROOT2 * 0.5;

  bChicNuSmul = u * bPsicNuSmul;
  bChi0MuSmul =  n * bPsi0MuSmul;
  aChicMuSnul = v.conjugate() * aPsicMuSnul;

  double deltaZnumu = 0.0, deltaZmu = 0.0;
  for (int i = 0; i < dimN; i++) {
     deltaZnumu += - sqr(std::abs(bChi0NuNul(i))) * b1(0.0, mneut(i), msnumu, q);
     deltaZmu   += - sqr(std::abs(bChi0MuSmul(i))) * b1(0.0, mneut(i), msmuL, q);
  }
  for (int i = 0; i < dimC; i++) {
     deltaZnumu += - sqr(std::abs(bChicNuSmul(i))) * b1(0.0, mch(i), msmuL, q);
     deltaZmu   += - sqr(std::abs(aChicMuSnul(i))) * b1(0.0, mch(i), msnumu, q);
  }

  Eigen::MatrixXd aPsi0PsicW(Eigen::MatrixXd::Zero(dimN,dimC)),
     bPsi0PsicW(Eigen::MatrixXd::Zero(dimN,dimC)),
     fW(Eigen::MatrixXd::Zero(dimN,dimC)),
     gW(Eigen::MatrixXd::Zero(dimN,dimC));
  Eigen::MatrixXcd aChi0ChicW, bChi0ChicW;

  aPsi0PsicW(1, 0) = - g;
  bPsi0PsicW(1, 0) = - g;
  aPsi0PsicW(3, 1) = g / ROOT2;
  bPsi0PsicW(2, 1) = -g / ROOT2;

  aChi0ChicW = n.conjugate() * aPsi0PsicW * v.transpose();
  bChi0ChicW = n * bPsi0PsicW * u.adjoint();

  const double b0_0_mselL_msnue = b0(0.0, mselL, msnue, q);
  const double b0_0_mnmuL_msnumu = b0(0.0, msmuL, msnumu, q);

  std::complex<double> deltaVE;
  for (int i = 0; i < dimC; i++) {
     for (int j = 0; j < dimN; j++) {
        const double c0_mselL_mch_mneut = c0(mselL, mch(i), mneut(j));
        const double c0_msnue_mch_mneut = c0(msnue, mch(i), mneut(j));
        const double b0_0_mch_mneut = b0(0.0, mch(i), mneut(j), q);

        deltaVE += bChicNuSell(i) * std::conj(bChi0ESell(j)) *
           (- ROOT2 / g * aChi0ChicW(j, i) * mch(i) * mneut(j) *
            c0_mselL_mch_mneut + 1.0 / (ROOT2 * g) *
            bChi0ChicW(j, i) *
            (b0_0_mch_mneut + sqr(mselL) *
             c0_mselL_mch_mneut - 0.5));
        deltaVE += - aChicESnul(i) * bChi0NuNul(j) *
           (- ROOT2 / g * bChi0ChicW(j, i) * mch(i) * mneut(j) *
            c0_msnue_mch_mneut + 1.0 / (ROOT2 * g) *
            aChi0ChicW(j, i) *
            (b0_0_mch_mneut + sqr(msnue) *
             c0_msnue_mch_mneut - 0.5));
     }
  }

  for (int j = 0; j < dimN; j++) {
     deltaVE +=
        0.5 * std::conj(bChi0ESell(j)) * bChi0NuNul(j) *
        (b0_0_mselL_msnue + sqr(mneut(j)) *
         c0(mneut(j), mselL, msnue) + 0.5);
  }

  std::complex<double> deltaVMu;
  for (int i = 0; i < dimC; i++) {
     for (int j = 0; j < dimN; j++) {
        const double c0_msmuL_mch_mneut = c0(msmuL, mch(i), mneut(j));
        const double c0_msnumu_mch_mneut = c0(msnumu, mch(i), mneut(j));
        const double b0_0_mch_mneut = b0(0.0, mch(i), mneut(j), q);

        deltaVMu += bChicNuSmul(i) * std::conj(bChi0MuSmul(j)) *
           (- ROOT2 / g * aChi0ChicW(j, i) * mch(i) * mneut(j) *
            c0_msmuL_mch_mneut + 1.0 / (ROOT2 * g) *
            bChi0ChicW(j, i) *
            (b0_0_mch_mneut + sqr(msmuL) *
             c0_msmuL_mch_mneut - 0.5));
        deltaVMu += - aChicMuSnul(i) * bChi0NuNul(j) *
           (- ROOT2 / g * bChi0ChicW(j, i) * mch(i) * mneut(j) *
            c0_msnumu_mch_mneut + 1.0 / (ROOT2 * g) *
            aChi0ChicW(j, i) *
            (b0_0_mch_mneut + sqr(msnumu) *
             c0_msnumu_mch_mneut - 0.5));
     }
  }

  for (int j = 0; j < dimN; j++) {
     deltaVMu +=
        0.5 * std::conj(bChi0MuSmul(j)) * bChi0NuNul(j) *
        (b0_0_mnmuL_msnumu + sqr(mneut(j)) *
         c0(mneut(j), msmuL, msnumu) + 0.5);
  }

  std::complex<double> a1;
  for(int i = 0; i < dimC; i++) {
     for(int j = 0; j < dimN; j++) {
        a1 += 0.5 * aChicMuSnul(i) * std::conj(bChicNuSell(i)) *
           bChi0NuNul(j) * bChi0ESell(j) * mch(i) * mneut(j) *
           d0(mselL, msnumu, mch(i), mneut(j));
        a1 += 0.5 * std::conj(aChicESnul(i)) * bChicNuSmul(i) *
           std::conj(bChi0NuNul(j)) * std::conj(bChi0MuSmul(j)) * mch(i) * mneut(j) *
           d0(msmuL, msnue, mch(i), mneut(j));
        a1 += bChicNuSmul(i) * std::conj(bChicNuSell(i)) *
           std::conj(bChi0MuSmul(j)) * bChi0ESell(j) *
           d27(msmuL, mselL, mch(i), mneut(j));
        a1 += std::conj(aChicMuSnul(i)) * aChicESnul(i) *
           bChi0NuNul(j) * std::conj(bChi0NuNul(j)) *
           d27(msnumu, msnue, mch(i), mneut(j));
     }
  }

  const double deltaVbSusy =
     (-sinThetaW2 * sqr(outcos) / (2.0 * Pi * alphaDRbar) * sqr(mz)
      * a1.real() + deltaVE.real() + deltaVMu.real() +
      0.5 * (deltaZe + deltaZnue + deltaZmu + deltaZnumu) ) * oneOver16PiSqr;

  return deltaVbSusy;
}

/**
 * Calculates \f$\rho^{(2)}(r)\f$ as given in Eqs. (C.7)-(C.8) from
 * hep-ph/9606211 .
 *
 * @param r ratio of Higgs mass over top quark mass
 *
 * @return \f$\rho^{(2)}(r)\f$
 */
double Weinberg_angle::rho_2(double r)
{
   const double Pi2 = Pi * Pi;

   if (r <= std::numeric_limits<double>::epsilon()) {
#if defined(ENABLE_VERBOSE) || defined(ENABLE_DEBUG)
      WARNING("rho_2: value of r is invalid: r = " << r);
      WARNING("-> setting 2-loop corrections ~ xt^2 to 0");
#endif
      return 0.;
   }

   if (r <= 1.9) {
      const double r2 = sqr(r);
      return 19.0 - 16.5 * r + 43.0 * r2 / 12.0 + 7.0 / 120.0 * r2 * r -
         Pi * sqrt(r) * (4.0 - 1.5 * r + 3.0 / 32.0 * r2 + r2 * r /
                         256.0) - Pi2 * (2.0 - 2.0 * r + 0.5 * r2) -
         std::log(r) * (3.0 * r - 0.5 * r2);
   } else {
      const double rm1 = 1.0 / r, rm2 = sqr(rm1), rm3 = rm2 * rm1,
         rm4 = rm3 * rm1, rm5 = rm4 * rm1;
      return sqr(std::log(r)) * (1.5 - 9.0 * rm1 - 15.0 * rm2 - 48.0 * rm3 - 168.0
                            * rm4 - 612.0 * rm5) -
         std::log(r) * (13.5 + 4.0 * rm1 - 125.0 / 4.0 * rm2 - 558.0 / 5.0 * rm3 -
                   8307.0 / 20.0 * rm4 - 109321.0 / 70.0 * rm5)
         + Pi2 * (1.0 - 4.0 * rm1 - 5.0 * rm2 - 16.0 * rm3 -
                      56.0 * rm4 - 204.0 * rm5)
         + 49.0 / 4.0 + 2.0 / 3.0 * rm1 + 1613.0 / 48.0 * rm2 + 87.57 * rm3 +
         341959.0 / 1200.0 * rm4 + 9737663.0 / 9800.0 * rm5;
   }
}

/**
 * Calculates 1-loop top-quark contribution to W boson self-energy.
 *
 * @param p momentum
 * @param mt top-quark mass
 * @param data model parameters
 *
 * @return 1-loop top-quark contribution to W boson self-energy
 */
double Weinberg_angle::calculate_self_energy_w_top(
   double p, double mt, const Self_energy_data& data)
{
   const double q = data.scale;
   const double mb = data.mb_drbar;
   const double Nc = 3.0;
   const double g2 = data.g2;

   const double self_energy_w =
      0.5 * Nc * hfn(p, mt, mb, q) * sqr(g2) * oneOver16PiSqr;

   return self_energy_w;
}

/**
 * Calculates 1-loop top-quark contribution to Z boson self-energy.
 *
 * @param p momentum
 * @param mt top-quark mass
 * @param data model parameters
 *
 * @return 1-loop top-quark contribution to Z boson self-energy
 */
double Weinberg_angle::calculate_self_energy_z_top(
   double p, double mt, const Self_energy_data& data)
{
   const double q = data.scale;
   const double Nc = 3.0;
   const double gY = data.gY;
   const double g2 = data.g2;
   const double gY2 = sqr(gY);
   const double g22 = sqr(g2);
   const double sw2 = gY2 / (gY2 + g22);
   const double cw2 = 1.0 - sw2;
   const double guL = 0.5 - 2.0 * sw2 / 3.0;
   const double guR = 2.0 * sw2 / 3.0;

   const double self_energy_z =
      Nc * (hfn(p, mt, mt, q) * (sqr(guL) + sqr(guR))
             - 4.0 * guL * guR * sqr(mt) * b0(p, mt, mt, q))
      * sqr(g2) / cw2 * oneOver16PiSqr;

   return self_energy_z;
}

/**
 * Replaces the 1-loop top-quark DR-bar mass contribution in the given
 * W boson self-energy by the 1-loop top-quark pole mass contribution.
 *
 * @param self_energy_w full W boson self-energy
 * @param p momentum
 * @param data model parameters
 *
 * @return W self-energy with top DR-bar mass replaced by top pole mass
 */
double Weinberg_angle::replace_mtop_in_self_energy_w(
   double self_energy_w, double p, const Self_energy_data& data)
{
   const double mt_pole = data.mt_pole;
   const double mt_drbar = data.mt_drbar;

   const double self_energy_w_mt_drbar
      = calculate_self_energy_w_top(p, mt_drbar, data);
   const double self_energy_w_mt_pole
      = calculate_self_energy_w_top(p, mt_pole, data);

   const double self_energy_w_with_mt_pole
      = self_energy_w - self_energy_w_mt_drbar + self_energy_w_mt_pole;

   return self_energy_w_with_mt_pole;
}

/**
 * Replaces the 1-loop top-quark DR-bar mass contribution in the given
 * Z boson self-energy by the 1-loop top-quark pole mass contribution.
 *
 * @param self_energy_z full Z boson self-energy
 * @param p momentum
 * @param data model parameters
 *
 * @return Z self-energy with top DR-bar mass replaced by top pole mass
 */
double Weinberg_angle::replace_mtop_in_self_energy_z(
   double self_energy_z, double p, const Self_energy_data& data)
{
   const double mt_pole = data.mt_pole;
   const double mt_drbar = data.mt_drbar;

   const double self_energy_z_mt_drbar
      = calculate_self_energy_z_top(p, mt_drbar, data);
   const double self_energy_z_mt_pole
      = calculate_self_energy_z_top(p, mt_pole, data);

   const double self_energy_z_with_mt_pole
      = self_energy_z - self_energy_z_mt_drbar + self_energy_z_mt_pole;

   return self_energy_z_with_mt_pole;
}

} // namespace weinberg_angle

} // namespace flexiblesusy
