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

#include "sm_fourloop_as.hpp"
#include <cmath>
#include <ostream>

namespace flexiblesusy {
namespace sm_fourloop_as {

namespace {
   constexpr double Pi = 3.1415926535897932384626433832795;
   template <typename T> constexpr T power2(T x) { return x*x; }
   template <typename T> constexpr T power3(T x) { return x*x*x; }
   template <typename T> constexpr T power4(T x) { return x*x*x*x; }
} // anonymous namespace

/**
 * 1-loop O(alpha_s) contributions to Delta alpha_s, Eq (23) of
 * [hep-ph/0004189]
 *
 * @param pars parameters
 *
 * @return 1-loop Delta alpha_s
 */
double delta_alpha_s_1loop_as(const Parameters& pars)
{
   const double as = pars.as;
   const double mt2 = power2(pars.mt);
   const double Q2 = power2(pars.Q);
   const double L = std::log(Q2/mt2);

   return as / Pi * (1./6. * L);
}

/**
 * 2-loop O(alpha_s^2) contributions to Delta alpha_s, Eq (23) of
 * [hep-ph/0004189]
 *
 * @param pars parameters
 *
 * @return 2-loop Delta alpha_s
 */
double delta_alpha_s_2loop_as_as(const Parameters& pars)
{
   const double as = pars.as;
   const double mt2 = power2(pars.mt);
   const double Q2 = power2(pars.Q);
   const double L = std::log(Q2/mt2);

   return power2(as / Pi) * (-11./72. + 11./24*L + 1./36. * power2(L));
}

/**
 * 3-loop O(alpha_s^3) contributions to Delta alpha_s, Eq (23) of
 * [hep-ph/0004189]
 *
 * @param pars parameters
 *
 * @return 3-loop Delta alpha_s
 */
double delta_alpha_s_3loop_as_as_as(const Parameters& pars)
{
   const double as = pars.as;
   const double mt2 = power2(pars.mt);
   const double Q2 = power2(pars.Q);
   const double L = std::log(Q2/mt2);
   const double L2 = power2(L);
   const double L3 = power3(L);
   const double nl = 5.;
   const double zeta3 = 1.202056903159594;

   return power3(as / Pi) * (
      - 564731./124416.
      + 82043./27648. * zeta3
      + 2645./1728. * L
      + 167./576. * L2
      + 1./216. * L3
      + nl * (2633./31104. - 67./576. * L + 1./36. * L2)
   );
}

/**
 * 4-loop O(alpha_s^4) contributions to Delta alpha_s, Eq (38) of
 * [hep-ph/0512060]
 *
 * @param pars parameters
 *
 * @return 4-loop Delta alpha_s
 */
double delta_alpha_s_4loop_as_as_as_as(const Parameters& pars)
{
   const double as = pars.as;
   const double mt2 = power2(pars.mt);
   const double Q2 = power2(pars.Q);
   const double L = std::log(Q2/mt2);
   const double L2 = power2(L);
   const double L3 = power3(L);
   const double L4 = power4(L);
   const double nl = 5.;
   const double nl2 = power2(nl);
   const double zeta3 = 1.202056903159594;
   const double delta_MS_4 =
      5.170346990805882 - 1.00993152453019 * nl - 0.0219783748689228 * nl2;

   return power4(as / Pi) * (
      + 121./1728.
      - delta_MS_4
      - 11093717./746496. * L
      + 3022001./165888. * zeta3 * L
      + 1837./1152. * L2
      + 2909./10368. * L3
      + 1./1296. * L4
      + nl * (
         + 141937./373248. * L
         - 110779./82944. * zeta3 * L
         + 277./10368. * L2
         + 271./5184. * L3
      )
      + nl2 * (
         - 6865./186624. * L
         + 77./20736. * L2
         - 1./324. * L3
      )
   );
}

/**
 * Calculates alpha_s(SM(6)) from alpha_s(SM(5)) in QCD using
 * (inverted) Eq (14) from [hep-ph/0512060].
 *
 * @param pars parameters
 *
 * @return alpha_s(SM(6))
 */
double calc_alpha_s(const Parameters& pars, int loops)
{
   const auto as = pars.as; // SM(5)
   const auto d1 = loops < 1 ? 0.0 : delta_alpha_s_1loop_as(pars);
   const auto d2 = loops < 2 ? 0.0 : delta_alpha_s_2loop_as_as(pars);
   const auto d3 = loops < 3 ? 0.0 : delta_alpha_s_3loop_as_as_as(pars);
   const auto d4 = loops < 4 ? 0.0 : delta_alpha_s_4loop_as_as_as_as(pars);

   return as * (1.0 + d1 + d2 + d3 + d4);
}


/**
 * Calculates alpha_s(SM(6)) from alpha_s(SM(5)) in QCD using
 * the expressions given in [hep-ph/0512060].
 *
 * @param pars parameters
 *
 * @return alpha_s(SM(6))
 */
double calc_alpha_s_alternative(const Parameters& pars, int loops)
{
   const auto as = pars.as; // SM(5)

   const auto d1 = loops < 1 ? 0.0 : delta_alpha_s_1loop_as(pars);
   const auto d2 = loops < 2 ? 0.0 : delta_alpha_s_2loop_as_as(pars);
   const auto d3 = loops < 3 ? 0.0 : delta_alpha_s_3loop_as_as_as(pars);
   const auto d4 = loops < 4 ? 0.0 : delta_alpha_s_4loop_as_as_as_as(pars);

   const auto del1 = d1;
   const auto del2 = d2 - power2(d1);
   const auto del3 = d3 + power3(d1) - 2.0*d1*d2;
   const auto del4 = d4 - 2.0*d1*d3 - power2(d2) + 3.0*power2(d1)*d2 - power4(d1);

   const auto delta = del1 + del2 + del3 + del4;

   return as / (1.0 - delta);
}

std::ostream& operator<<(std::ostream& out, const Parameters& pars)
{
   out <<
      "Delta alpha_s(SM) parameters:\n"
      "alpha_s(SM(5)) = " <<  pars.as   << '\n' <<
      "mt             = " <<  pars.mt   << '\n' <<
      "Q              = " <<  pars.Q    << '\n';

   return out;
}

} // namespace sm_fourloop_as
} // namespace flexiblesusy
