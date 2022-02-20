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

#include "sm_mw.hpp"
#include "logger.hpp"
#include <cmath>

namespace flexiblesusy {
namespace sm_mw {

namespace {

constexpr double sqr(double x) noexcept { return x*x; }

} // anonymous namespace

/** 
 * Calculates the prediction for the W boson mass in the Standard
 * Model, using the fit formula Eq (45) from arXiv:1411.7040 (MS-bar
 * calculation).
 * 
 * @param mh SM Higgs boson pole mass
 * @param mt SM top quark pole mass
 * @param as MS-bar alpha_s(MZ) in the SM with 5 quark flavours
 * @param da5had hadronic contributison Delta alpha_{had}^{(5)}(MZ^2)
 *
 * The authors of arXiv:1411.7040 used the following input values (Table 1):
 *
 * mz = 91.1876 GeV (Z boson pole mass)
 * mh = 125.15 GeV (SM Higgs boson pole mass)
 * mt = 173.34 GeV
 * as = 0.1184
 * da5had = 0.02750
 *
 * @return W boson pole mass as predicted in the Standard Model (first
 * entry) and corresponding uncertainty (second entry)
 */
std::pair<double, double> calculate_mw_pole_SM_fit_MSbar(double mh, double mt, double as, double da5had) noexcept
{
   // Table 3, 2nd column, 124.42 <= mh <= 125.87
   const double p[8] = {
      80.35712, -0.06017, 0.0, 0.0, 0.52749, -0.00613, -0.08178, -0.50530
   };
   // Table 3, 3rd column, 50 <= mh <= 450
   const double q[8] = {
      80.35714, -0.06094, -0.00971, 0.00028, 0.52655, -0.00646, -0.08199, -0.50259
   };

   const double das = as/0.1184 - 1;      // Eq.(20)
   const double da5 = da5had/0.02750 - 1; // defined below Eq.(44)
   const double dh = sqr(mh/125.15) - 1;  // defined below Eq.(45)
   const double dH = std::log(mh/125.15); // Eq.(20)
   const double dt = sqr(mt/173.34) - 1;  // defined below Eq.(41)

   const double* w = (124.42 <= mh && mh <= 125.87) ? p : q;
   const double dmw = (124.42 <= mh && mh <= 125.87) ? 0.11e-3 : 0.5e-3; // below Eq.(45)

   // Eq.(45)
   const double mw = w[0] + w[1]*dH + w[2]*dH*dH + w[3]*dh + w[4]*dt
      + w[5]*dH*dt + w[6]*das + w[7]*da5;

   if (mw < 0) {
      WARNING("calculate_mw_pole_SM_fit_MSbar: Standard Model MW " << mw << " < 0");
   }

   return std::make_pair(std::abs(mw), dmw);
}

} // namespace sm_mw
} // namespace flexiblesusy
