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

/** \file lowe.cpp
   - Project:     SOFTSUSY
   - Author:      Ben Allanach, Alexander Voigt
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/
*/

#include "lowe.h"
#include "eigen_utils.hpp"
#include "error.hpp"
#include "ew_input.hpp"
#include "string_format.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

namespace softsusy {

namespace {

constexpr double sqr(double a) noexcept { return a*a; }

// Given a value of mt, and alphas(MZ), find alphas(mt) to 1 loops in qcd:
// it's a very good approximation at these scales, better than 10^-3 accuracy
double getAsmt(double mtop, double alphasMz, double mz) {
  return alphasMz /
      (1.0 - 23.0 * alphasMz / (6.0 * M_PI) * std::log(mz / mtop));
}

// Input pole mass of top and alphaS(mt), outputs running mass mt(mt)
// including one-loop standard model correction only
double getRunMt(double poleMt, double asmt) {
  return poleMt / (1.0 + (4.0 / (3.0 * M_PI)) * asmt);
}

// Given pole mass and alphaS(MZ), returns running top mass -- one loop qcd
double getRunMtFromMz(double poleMt, double asMZ, double mz) {
  return getRunMt(poleMt, getAsmt(poleMt, asMZ, mz));
}

} // anonymous namespace

const std::array<std::string, NUMBER_OF_LOW_ENERGY_INPUT_PARAMETERS> QedQcd_input_parmeter_names = {
   "alpha_em_MSbar_at_MZ",
   "alpha_s_MSbar_at_MZ",
   "GFermi",
   "MZ_pole", "MW_pole",
   "Mv1_pole", "Mv2_pole", "Mv3_pole",
   "Me_pole", "Mm_pole", "Mtau_pole",
   "mu_2GeV", "ms_2GeV", "Mt_pole",
   "md_2GeV", "mc_mc", "mb_mb",
   "CKM_theta_12", "CKM_theta_13", "CKM_theta_23", "CKM_delta",
   "PMNS_theta_12", "PMNS_theta_13", "PMNS_theta_23", "PMNS_delta", "PMNS_alpha_1", "PMNS_alpha_2"
};

QedQcd::QedQcd()
   : mbPole(flexiblesusy::Electroweak_constants::PMBOTTOM)
{
  set_number_of_parameters(a.size() + mf.size());
  mf(0) = flexiblesusy::Electroweak_constants::MUP;
  mf(1) = flexiblesusy::Electroweak_constants::MCHARM;
  mf(2) = getRunMtFromMz(flexiblesusy::Electroweak_constants::PMTOP,
                         flexiblesusy::Electroweak_constants::alpha3,
                         flexiblesusy::Electroweak_constants::MZ);
  mf(3) = flexiblesusy::Electroweak_constants::MDOWN;
  mf(4) = flexiblesusy::Electroweak_constants::MSTRANGE;
  mf(5) = flexiblesusy::Electroweak_constants::MBOTTOM;
  mf(6) = flexiblesusy::Electroweak_constants::MELECTRON;
  mf(7) = flexiblesusy::Electroweak_constants::MMUON;
  mf(8) = flexiblesusy::Electroweak_constants::MTAU;
  a(0) = flexiblesusy::Electroweak_constants::aem;
  a(1) = flexiblesusy::Electroweak_constants::alpha3;
  input(alpha_em_MSbar_at_MZ) = flexiblesusy::Electroweak_constants::aem;
  input(alpha_s_MSbar_at_MZ) = flexiblesusy::Electroweak_constants::alpha3;
  input(Mt_pole) = flexiblesusy::Electroweak_constants::PMTOP;
  input(mb_mb) = flexiblesusy::Electroweak_constants::MBOTTOM;
  input(Mtau_pole) = flexiblesusy::Electroweak_constants::MTAU;
  input(Mm_pole) = flexiblesusy::Electroweak_constants::MMUON;
  input(Me_pole) = flexiblesusy::Electroweak_constants::MELECTRON;
  input(MW_pole) = flexiblesusy::Electroweak_constants::MW;
  input(MZ_pole) = flexiblesusy::Electroweak_constants::MZ;
  input(GFermi) = flexiblesusy::Electroweak_constants::gfermi;
  input(mc_mc) = flexiblesusy::Electroweak_constants::MCHARM;
  input(mu_2GeV) = flexiblesusy::Electroweak_constants::MUP;
  input(md_2GeV) = flexiblesusy::Electroweak_constants::MDOWN;
  input(ms_2GeV) = flexiblesusy::Electroweak_constants::MSTRANGE;
  set_scale(flexiblesusy::Electroweak_constants::MZ);
  set_loops(3);
  set_thresholds(1);
}

Eigen::ArrayXd QedQcd::get() const
{
   Eigen::ArrayXd y(a.size() + mf.size());
   y(0) = a(0);
   y(1) = a(1);
   for (int i = 0; i < mf.size(); i++) {
      y(i + 2) = mf(i);
   }
   return y;
}

void QedQcd::set(const Eigen::ArrayXd& y)
{
   a(0) = y(0);
   a(1) = y(1);
   for (int i = 0; i < mf.size(); i++) {
      mf(i) = y(i + 2);
   }
}

Eigen::ArrayXd QedQcd::beta() const
{
   Eigen::ArrayXd dydx(a.size() + mf.size());
   dydx(0) = qedBeta();
   dydx(1) = qcdBeta();
   const auto y = massBeta();
   for (int i = 0; i < y.size(); i++) {
      dydx(i + 2) = y(i);
   }
   return dydx;
}

void QedQcd::runto_safe(double scale, double eps)
{
   try {
      run_to(scale, eps);
   } catch (...) {
      throw flexiblesusy::NonPerturbativeRunningQedQcdError(
         "Non-perturbative running to Q = "
         + flexiblesusy::to_string(scale)
         + " during determination of the SM(5) parameters.");
   }
}

//  Active flavours at energy mu
int QedQcd::flavours(double mu) const {
  int k = 0;
  // if (mu > mf(mTop - 1)) k++;
  if (mu > mf(mCharm - 1)) { k++; }
  if (mu > mf(mUp - 1)) { k++; }
  if (mu > mf(mDown - 1)) { k++; }
  if (mu > mf(mBottom - 1)) { k++; }
  if (mu > mf(mStrange - 1)) { k++; }
  return k;
}

void QedQcd::setCKM(const flexiblesusy::CKM_parameters& ckm)
{
   input(CKM_theta_12) = ckm.theta_12;
   input(CKM_theta_13) = ckm.theta_13;
   input(CKM_theta_23) = ckm.theta_23;
   input(CKM_delta)    = ckm.delta;
}

void QedQcd::setPMNS(const flexiblesusy::PMNS_parameters& pmns)
{
   input(PMNS_theta_12) = pmns.theta_12;
   input(PMNS_theta_13) = pmns.theta_13;
   input(PMNS_theta_23) = pmns.theta_23;
   input(PMNS_delta)    = pmns.delta;
   input(PMNS_alpha_1)  = pmns.alpha_1;
   input(PMNS_alpha_2)  = pmns.alpha_2;
}

flexiblesusy::CKM_parameters QedQcd::displayCKM() const
{
   flexiblesusy::CKM_parameters ckm;
   ckm.theta_12 = input(CKM_theta_12);
   ckm.theta_13 = input(CKM_theta_13);
   ckm.theta_23 = input(CKM_theta_23);
   ckm.delta    = input(CKM_delta);
   return ckm;
}

flexiblesusy::PMNS_parameters QedQcd::displayPMNS() const
{
   flexiblesusy::PMNS_parameters pmns;
   pmns.theta_12 = input(PMNS_theta_12);
   pmns.theta_13 = input(PMNS_theta_13);
   pmns.theta_23 = input(PMNS_theta_23);
   pmns.delta    = input(PMNS_delta);
   pmns.alpha_1  = input(PMNS_alpha_1);
   pmns.alpha_2  = input(PMNS_alpha_2);
   return pmns;
}

std::ostream& operator<<(std::ostream &left, const QedQcd &m) {
  left << "mU: " << m.displayMass(mUp)
       << "  mC: " << m.displayMass(mCharm)
       << "  mt: " << m.displayMass(mTop)
       << "  mt^pole: " << m.displayPoleMt()
       << '\n';
  left << "mD: " << m.displayMass(mDown)
       << "  mS: " << m.displayMass(mStrange)
       << "  mB: " << m.displayMass(mBottom)
       << "  mb(mb):  " << m.displayMbMb()
       << '\n';
  left << "mE: " << m.displayMass(mElectron)
       << "  mM: " << m.displayMass(mMuon)
       <<  "  mT: " << m.displayMass(mTau)
       << "  mb^pole: " << m.displayPoleMb()
       << '\n';
  left << "aE: " << 1.0 / m.displayAlpha(ALPHA)
       << "  aS: " << m.displayAlpha(ALPHAS)
       << "   Q: " << m.get_scale()
       << "  mT^pole: " << m.displayPoleMtau()
       << '\n';
  left << "loops: " << m.get_loops()
       << "        thresholds: " << m.get_thresholds() << '\n';

  return left;
}

/// returns QED beta function in SM(5) (without the top quark)
double QedQcd::qedBeta() const {
  double x = 24.0 / 9.0;

  if (get_scale() > mf(mCharm - 1)) { x += 8.0 / 9.0; }
  // if (get_scale() > mf(mTop - 1)) { x += 8.0 / 9.0; }
  if (get_scale() > mf(mBottom - 1)) { x += 2.0 / 9.0; }
  if (get_scale() > mf(mTau - 1)) { x += 2.0 / 3.0; }
  if (get_scale() > displayPoleMW()) { x += -7.0 / 2.0; }

  return (x * sqr(a(ALPHA - 1)) / M_PI);
}

/// Returns QCD beta function to 3 loops in QCD for the SM(5). Note
/// that if quark masses are running, the number of active quarks will
/// be taken into account.
double QedQcd::qcdBeta() const {
  static const double INVPI = 1.0 / M_PI;
  const int quarkFlavours = flavours(get_scale());
  const double qb0 = (11.0e0 - (2.0e0 / 3.0e0 * quarkFlavours)) / 4.0;
  const double qb1 = (102.0e0 - (38.0e0 * quarkFlavours) / 3.0e0) / 16.0;
  const double qb2 = (2.857e3 * 0.5 - (5.033e3 * quarkFlavours) / 18.0  +
                      (3.25e2 * sqr(quarkFlavours) ) / 5.4e1) / 64;

  double qa0 = 0., qa1 = 0., qa2 = 0.;

  if (get_loops() > 0) {
     qa0 = qb0 * INVPI;
  }
  if (get_loops() > 1) {
     qa1 = qb1 * sqr(INVPI);
  }
  if (get_loops() > 2) {
     qa2 = qb2 * sqr(INVPI) * INVPI;
  }

  // add contributions of the one, two and three loop constributions resp.
  const double beta =
    -2.0 * sqr(displayAlpha(ALPHAS)) *
    (qa0 + qa1 * displayAlpha(ALPHAS) + qa2 *
     sqr(displayAlpha(ALPHAS)));

  return beta;
}

/// returns fermion mass beta functions
Eigen::Array<double,9,1> QedQcd::massBeta() const {
  static const double INVPI = 1.0 / M_PI, ZETA3 = 1.202056903159594;

  // qcd bits: 1,2,3 loop resp.
  double qg1 = 0., qg2 = 0., qg3 = 0.;
  const int quarkFlavours = flavours(get_scale());

  if (get_loops() > 0) {
     qg1 = INVPI;
  }
  if (get_loops() > 1) {
     qg2 = (202.0 / 3.0 - (20.0e0 * quarkFlavours) / 9.0) * sqr(INVPI) / 16.0;
  }
  if (get_loops() > 2) {
     qg3 = (1.249e3 - ((2.216e3 * quarkFlavours) / 27.0e0 +
                       1.6e2 * ZETA3 * quarkFlavours / 3.0e0) -
            140.0e0 * quarkFlavours * quarkFlavours / 81.0e0) * sqr(INVPI) *
        INVPI / 64.0;
  }

  const double qcd = -2.0 * a(ALPHAS - 1) * (
     qg1  + qg2 * a(ALPHAS - 1) + qg3 * sqr(a(ALPHAS - 1)));
  const double qed = -a(ALPHA - 1) * INVPI / 2;

  Eigen::Array<double,9,1> x(Eigen::Array<double,9,1>::Zero());

  for (int i = 0; i < 3; i++) {   // up quarks
     x(i) = (qcd + 4.0 * qed / 3.0) * mf(i);
  }
  for (int i = 3; i < 6; i++) {   // down quarks
     x(i) = (qcd + qed / 3.0) * mf(i);
  }
  for (int i = 6; i < 9; i++) {   // leptons
     x(i) = 3.0 * qed * mf(i);
  }

  // switch off relevant beta functions
  if (get_thresholds() > 0) {
     for(int i = 0; i < x.size(); i++) {
        if (get_scale() < mf(i)) {
           x(i) = 0.0;
        }
     }
  }
  // nowadays, u,d,s masses defined at 2 GeV: don't run them below that
  if (get_scale() < 2.0) {
     x(mUp - 1) = x(mDown - 1) = x(mStrange - 1) = 0.0;
  }

  return x;
}

/// Supposed to be done at mb(mb) -- MSbar, calculates pole mass
double QedQcd::extractPoleMb(double alphasMb)
{
  if (get_scale() != displayMass(mBottom)) {
    throw flexiblesusy::SetupError(
       "QedQcd::extractPoleMb called at scale "
       + flexiblesusy::to_string(get_scale()) + " instead of mb(mb)");
  }

  // Following is the MSbar correction from QCD, hep-ph/9912391
  double delta = 0.0;

  if (get_loops() > 0) {
     delta = delta + 4.0 / 3.0 * alphasMb / M_PI;
  }
  if (get_loops() > 1) {
     delta = delta + sqr(alphasMb / M_PI) *
        (9.2778 + (displayMass(mUp) + displayMass(mDown) + displayMass(mCharm) +
                   displayMass(mStrange)) / mbPole);
  }
  if (get_loops() > 2) {
     delta = delta + 94.4182 * alphasMb / M_PI * sqr(alphasMb / M_PI);
  }

  const double mbPole = displayMass(mBottom) * (1.0 + delta);

  return mbPole;
}

/// Takes QedQcd object created at MZ and spits it out at MZ
void QedQcd::toMz()
{
   to(displayPoleMZ());
}

/**
 * Calculates all running parameters in the SM w/o top quark at Q.
 * This function can be called multiple times, leading to the same
 * result.
 *
 * @param scale target renormalization scale
 * @param precision_goal precision goal
 * @param max_iterations maximum number of iterations
 */
void QedQcd::to(double scale, double precision_goal, int max_iterations) {
   int it = 0;
   bool converged = false;
   auto qedqcd_old(get()), qedqcd_new(get());
   const double running_precision = 0.1 * precision_goal;

   while (!converged && it < max_iterations) {
      // set alpha_i(MZ)
      runto_safe(displayPoleMZ(), running_precision);
      setAlpha(ALPHA, input(alpha_em_MSbar_at_MZ));
      setAlpha(ALPHAS, input(alpha_s_MSbar_at_MZ));

      // set mb(mb)
      runto_safe(displayMbMb(), running_precision);
      setMass(mBottom, displayMbMb());
      setPoleMb(extractPoleMb(displayAlpha(ALPHAS)));

      // set mc(mc)
      runto_safe(displayMcMc(), running_precision);
      setMass(mCharm, displayMcMc());

      // set mu, md, ms at 2 GeV
      runto_safe(2.0, running_precision);
      setMass(mUp, displayMu2GeV());
      setMass(mDown, displayMd2GeV());
      setMass(mStrange, displayMs2GeV());

      // set me, mm, ml at 2 GeV
      setMass(mElectron, displayPoleMel());
      setMass(mMuon, displayPoleMmuon());
      setMass(mTau, displayPoleMtau());

      // check convergence
      runto_safe(scale, running_precision);
      qedqcd_new = get();

      converged = flexiblesusy::is_equal_rel(qedqcd_old, qedqcd_new, precision_goal);

      qedqcd_old = qedqcd_new;

      it++;
   }

   // set alpha_i(MZ) on last time
   runto_safe(displayPoleMZ(), precision_goal);
   setAlpha(ALPHA, input(alpha_em_MSbar_at_MZ));
   setAlpha(ALPHAS, input(alpha_s_MSbar_at_MZ));

   runto_safe(scale, precision_goal);

   if (!converged && max_iterations > 0) {
      std::string msg =
         "Iteration to determine SM(5) parameters did not"
         " converge after " + std::to_string(max_iterations) +
         " iterations (precision goal: " + std::to_string(precision_goal)
         + ").";
      throw flexiblesusy::NoConvergenceError(max_iterations, msg);
   }
}

/**
 * Returns the three coupling constants of the Standard Model without
 * the top quark (SM(5)) at the given [scale].
 *
 * @note The returned alpha_1 is in GUT-normalized.
 *
 * @param scale output scale
 *
 * @return {alpha_1, alpha_2, alpha_3}
 */
Eigen::Array<double,3,1> QedQcd::guess_alpha_SM5(double scale) const
{
  auto oneset = *this;
  oneset.runto_safe(scale);

  const double aem = oneset.displayAlpha(ALPHA);
  const double MW = oneset.displayPoleMW();
  const double MZ = oneset.displayPoleMZ();
  const double sin2th = 1. - sqr(MW / MZ);

  Eigen::Array<double,3,1> alpha(Eigen::Array<double,3,1>::Zero());

  alpha(0) = 5.0 * aem / (3.0 * (1.0 - sin2th));
  alpha(1) = aem / sin2th;
  alpha(2) = oneset.displayAlpha(ALPHAS);

  return alpha;
}

std::array<std::string, NUMBER_OF_LOW_ENERGY_INPUT_PARAMETERS> QedQcd::display_input_parameter_names()
{
   return QedQcd_input_parmeter_names;
}

bool operator ==(const QedQcd& a, const QedQcd& b)
{
   const double eps = 1e-10;

   return
      std::fabs(a.get_scale() - b.get_scale()) < eps &&
      std::fabs(a.get_loops() - b.get_loops()) < eps &&
      std::fabs(a.get_thresholds() - b.get_thresholds()) < eps &&
      std::fabs(a.displayAlpha(ALPHA) - b.displayAlpha(ALPHA)) < eps &&
      std::fabs(a.displayAlpha(ALPHAS) - b.displayAlpha(ALPHAS)) < eps &&
      std::fabs(a.displayMass(mUp) - b.displayMass(mUp)) < eps &&
      std::fabs(a.displayMass(mCharm) - b.displayMass(mCharm)) < eps &&
      std::fabs(a.displayMass(mTop) - b.displayMass(mTop)) < eps &&
      std::fabs(a.displayMass(mDown) - b.displayMass(mDown)) < eps &&
      std::fabs(a.displayMass(mStrange) - b.displayMass(mStrange)) < eps &&
      std::fabs(a.displayMass(mBottom) - b.displayMass(mBottom)) < eps &&
      std::fabs(a.displayMass(mElectron) - b.displayMass(mElectron)) < eps &&
      std::fabs(a.displayMass(mMuon) - b.displayMass(mMuon)) < eps &&
      std::fabs(a.displayMass(mTau) - b.displayMass(mTau)) < eps &&
      std::fabs(a.displayNeutrinoPoleMass(1) - b.displayNeutrinoPoleMass(1)) < eps &&
      std::fabs(a.displayNeutrinoPoleMass(2) - b.displayNeutrinoPoleMass(2)) < eps &&
      std::fabs(a.displayNeutrinoPoleMass(3) - b.displayNeutrinoPoleMass(3)) < eps &&
      std::fabs(a.displayPoleMt() - b.displayPoleMt()) < eps &&
      std::fabs(a.displayPoleMb() - b.displayPoleMb()) < eps &&
      std::fabs(a.displayPoleMtau() - b.displayPoleMtau()) < eps &&
      std::fabs(a.displayPoleMW() - b.displayPoleMW()) < eps &&
      std::fabs(a.displayPoleMZ() - b.displayPoleMZ()) < eps &&
      std::fabs(a.displayFermiConstant() - b.displayFermiConstant()) < eps;
}

} // namespace softsusy
