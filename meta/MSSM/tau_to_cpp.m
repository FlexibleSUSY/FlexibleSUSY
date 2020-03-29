(* :Copyright:

   ====================================================================
   This file is part of FlexibleSUSY.

   FlexibleSUSY is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published
   by the Free Software Foundation, either version 3 of the License,
   or (at your option) any later version.

   FlexibleSUSY is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with FlexibleSUSY.  If not, see
   <http://www.gnu.org/licenses/>.
   ====================================================================

*)

Get[FileNameJoin[{"meta", "TextFormatting.m"}]];

scale = Q;
At  = xt + MUE CB/SB;
Ab  = xb + MUE SB/CB;
Al  = xtau + MUE SB/CB;
YT  = yt; YB = yb; YL = ytau;
MW  = mw; MZ = mz;
mh0 = mh; mH0 = mH; mA0 = mA; mHp = mC;
MUE = mu; SB = sb; CB = cb; SA = sa; CA = ca;

mw /: mw^n_ := mw2^(n/2) /; EvenQ[n];
mz /: mz^n_ := mz2^(n/2) /; EvenQ[n];
mt /: mt^n_ := mt2^(n/2) /; EvenQ[n];
mtau /: mtau^n_ := mtau2^(n/2) /; EvenQ[n];
mh /: mh^n_ := mh2^(n/2) /; EvenQ[n];
mH /: mH^n_ := mH2^(n/2) /; EvenQ[n];
mA /: mA^n_ := mA2^(n/2) /; EvenQ[n];
mC /: mC^n_ := mC2^(n/2) /; EvenQ[n];
mu /: mu^n_ := mu2^(n/2) /; EvenQ[n];
Q  /: Q^n_  := Q2^(n/2)  /; EvenQ[n];

mst1 /: mst1^2 := mst12;
mst2 /: mst2^2 := mst22;
mst1 /: mst1^4 := mst14;
mst2 /: mst2^4 := mst24;
msb1 /: msb1^2 := msb12;
msb2 /: msb2^2 := msb22;
msb1 /: msb1^4 := msb14;
msb2 /: msb2^4 := msb24;
mstau1 /: mstau1^2 := mstau12;
mstau2 /: mstau2^2 := mstau22;
mstau1 /: mstau1^4 := mstau14;
mstau2 /: mstau2^4 := mstau24;
msntau /: msntau^2 := msntau2;
msntau /: msntau^(-2) := 1/msntau2;
msntau /: msntau^4 := msntau4;

tau2l = Get[FileNameJoin[{"meta", "MSSM", "dmtauas2.m"}]];

tau2lyl4    = Coefficient[tau2l, ytau^4];
tau2lyl2yt2 = Coefficient[tau2l, ytau^2 yt^2];
tau2lyl2yb2 = Coefficient[tau2l, ytau^2 yb^2];

fpart[0, m1_, m2_, scale_]    := Fin20[m1^2, m2^2, scale^2];
fpart[m1_, m2_, m3_, Q_]      := Fin3[m1^2, m2^2, m3^2, Q^2];
delta3[m1_, m2_, m3_]         := Delta[m1^2, m2^2, m3^2, -1]; 
Delta[m1_, m2_, m3_, -1]      := DeltaInv[m1,m2,m3];

Simp[expr_] := Collect[expr, {xt, xb, xtau, Fin3[__]}, Simplify[#,TimeConstraint->1200]&] //. {
        Power[x_,n_] /; n > 0 :> Symbol["power" <> ToString[n]][x],
        Power[x_,-2]          :> 1/Symbol["power" <> ToString[2]][x],
        Power[x_,-3]          :> 1/Symbol["power" <> ToString[3]][x],
        Power[x_,-4]          :> 1/Symbol["power" <> ToString[4]][x],
        Power[x_,-5]          :> 1/Symbol["power" <> ToString[5]][x],
        Power[x_,-6]          :> 1/Symbol["power" <> ToString[6]][x],
        Log[x_]               :> log[x]
    };

ToCPP[expr_] := ToString[Simp[expr], CForm];

headerName = "mssm_twoloop_mtau.hpp";
implName   = "mssm_twoloop_mtau.cpp";

header = "\
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

// This file has been generated at " <> DateString[] <> "
// with the script \"tau_to_cpp.m\".

#ifndef MSSM_TWO_LOOP_MTAU_H
#define MSSM_TWO_LOOP_MTAU_H

#include <iosfwd>

namespace flexiblesusy {
namespace mssm_twoloop_mtau {

using Real = long double;

struct Parameters {
    Real yt{};     ///< MSSM top Yukawa coupling DR-bar
    Real yb{};     ///< MSSM bottom Yukawa coupling DR-bar
    Real ytau{};   ///< MSSM tau Yukawa coupling DR-bar
    Real mt{};     ///< MSSM top mass DR-bar
    Real mb{};     ///< MSSM bottom mass DR-bar
    Real mtau{};   ///< MSSM tau mass DR-bar
    Real mst1{};   ///< MSSM light stop mass DR-bar
    Real mst2{};   ///< MSSM heavy stop mass DR-bar
    Real msb1{};   ///< MSSM light sbottom mass DR-bar
    Real msb2{};   ///< MSSM heavy sbottom mass DR-bar
    Real mstau1{}; ///< MSSM light stau mass DR-bar
    Real mstau2{}; ///< MSSM heavy stau mass DR-bar
    Real msntau{}; ///< MSSM tau sneutrino mass DR-bar
    Real xt{};     ///< MSSM stop mixing parameter DR-bar
    Real xb{};     ///< MSSM sbottom mixing parameter DR-bar
    Real xtau{};   ///< MSSM stau mixing parameter DR-bar
    Real mw{};     ///< MSSM W boson mass DR-bar
    Real mz{};     ///< MSSM Z boson mass DR-bar
    Real mh{};     ///< MSSM light CP-even Higgs mass DR-bar
    Real mH{};     ///< MSSM heavy CP-even Higgs mass DR-bar
    Real mC{};     ///< MSSM charged Higgs mass DR-bar
    Real mA{};     ///< MSSM CP-odd Higgs mass DR-bar
    Real mu{};     ///< MSSM mu superpotential parameter DR-bar
    Real tb{};     ///< MSSM tan(beta) DR-bar
    Real Q{};      ///< renormalization scale
};

/// 2-loop contribution to mtau O(alpha_tau^2) [0912.4652]
double delta_mtau_2loop_atau_atau(const Parameters&);

/// 2-loop contribution to mtau O(alpha_tau*alpha_t) [0912.4652]
double delta_mtau_2loop_atau_at(const Parameters&);

/// 2-loop contribution to mtau O(alpha_tau*alpha_b) [0912.4652]
double delta_mtau_2loop_atau_ab(const Parameters&);

std::ostream& operator<<(std::ostream&, const Parameters&);

} // namespace mssm_twoloop_mtau
} // namespace flexiblesusy

#endif
";

impl = "\
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

// This file has been generated at " <> DateString[] <> "
// with the script \"tau_to_cpp.m\".

#include \"" <> headerName <> "\"
#include \"dilog.hpp\"
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <ostream>

namespace flexiblesusy {
namespace mssm_twoloop_mtau {

namespace {
   const Real Pi = 3.1415926535897932384626433832795l;

   template <typename T> T power2(T x)  { return x*x; }
   template <typename T> T power3(T x)  { return x*x*x; }
   template <typename T> T power4(T x)  { return x*x*x*x; }
   template <typename T> T power5(T x)  { return x*x*x*x*x; }
   template <typename T> T power6(T x)  { return x*x*x*x*x*x; }
   template <typename T> T power7(T x)  { return x*x*x*x*x*x*x; }
   template <typename T> T power8(T x)  { return x*x*x*x*x*x*x*x; }
   template <typename T> T power10(T x) { return x*x*x*x*x * x*x*x*x*x; }
   template <typename T> T power12(T x) { return x*x*x*x*x * x*x*x*x*x * x*x; }
   template <typename T> T power14(T x) { return x*x*x*x*x * x*x*x*x*x * x*x*x*x; }

   const Real oneLoop = 1.l/power2(4*Pi);
   const Real twoLoop = power2(oneLoop);

   template <typename T>
   bool is_zero(T a, T prec = std::numeric_limits<T>::epsilon())
   {
      return std::fabs(a) < prec;
   }

   template <typename T>
   bool is_equal(T a, T b, T prec = std::numeric_limits<T>::epsilon())
   {
      return is_zero(a - b, prec);
   }

   template <typename T>
   bool is_equal_rel(T a, T b, T prec = std::numeric_limits<T>::epsilon())
   {
      if (is_equal(a, b, std::numeric_limits<T>::epsilon()))
         return true;

      if (std::abs(a) < std::numeric_limits<T>::epsilon() ||
          std::abs(b) < std::numeric_limits<T>::epsilon())
         return false;

      return std::abs((a - b)/a) < prec;
   }

   /**
    * Fin20[] function from twoloopbubble.m .
    *
    * @param mm1 squared mass \\f$m_1^2\\f$
    * @param mm2 squared mass \\f$m_2^2\\f$
    * @param mmu squared renormalization scale
    *
    * @return Fin20(m12, m22, mmu)
    */
   Real Fin20(Real mm1, Real mm2, Real mmu)
   {
      using std::log;

      return (6*(mm1*log(mm1/mmu) + mm2*log(mm2/mmu)) +
         (-mm1 - mm2)*(7 + power2(Pi)/6.) +
         (mm1 - mm2)*(2*dilog(1 - mm1/mm2) +
            power2(log(mm1/mm2))/2.) +
         ((mm1 + mm2)*power2(log(mm1/mm2)))/2. -
         2*(mm1*power2(log(mm1/mmu)) + mm2*power2(log(mm2/mmu))))/2.;
   }

   Real LambdaSquared(Real x, Real y)
   {
      return power2(1 - x - y) - 4*x*y;
   }

   /// ClausenCl[2,x]
   Real ClausenCl2(Real x)
   {
      using std::exp;
      const std::complex<Real> img(0.l,1.l);

      return std::imag(dilog(exp(img*x)));
   }

   /// x < 1 && y < 1, LambdaSquared(x,y) > 0
   Real PhiPos(Real x, Real y)
   {
      const Real lambda = std::sqrt(LambdaSquared(x,y));

      return (-(log(x)*log(y))
              + 2*log((1 - lambda + x - y)/2.)*log((1 - lambda - x + y)/2.)
              - 2*dilog((1 - lambda + x - y)/2.)
              - 2*dilog((1 - lambda - x + y)/2.)
              + power2(Pi)/3.)/lambda;
   }

   /// LambdaSquared(x,y) < 0
   Real PhiNeg(Real x, Real y)
   {
      using std::acos;
      using std::sqrt;
      const Real lambda = std::sqrt(-LambdaSquared(x,y));

      return 2*(+ ClausenCl2(2*acos((1 + x - y)/(2.*sqrt(x))))
                + ClausenCl2(2*acos((1 - x + y)/(2.*sqrt(y))))
                + ClausenCl2(2*acos((-1 + x + y)/(2.*sqrt(x*y)))))/lambda;
   }

   Real Phi(Real x, Real y)
   {
      const Real lambda = LambdaSquared(x,y);

      if (lambda > 0.)
         return PhiPos(x,y);

      return PhiNeg(x,y);
   }

   /**
    * Fin3[] function from twoloopbubble.m .
    *
    * @param mm1 squared mass \\f$m_1^2\\f$
    * @param mm2 squared mass \\f$m_2^2\\f$
    * @param mm3 squared mass \\f$m_3^2\\f$
    * @param mmu squared renormalization scale
    *
    * @return Fin3(m12, m22, m32, mmu)
    */
   Real Fin3(Real mm1, Real mm2, Real mm3, Real mmu)
   {
      using std::log;

      std::array<Real,3> masses = { mm1, mm2, mm3 };
      std::sort(masses.begin(), masses.end());

      const Real mm = masses[2];
      const Real x = masses[0]/mm;
      const Real y = masses[1]/mm;

      const Real lambda = LambdaSquared(x,y);

      if (is_zero(lambda, 1e-10l)) {
         return -(mm*(2*y*(-3 + 2*log(mm/mmu))*log(y)
                      + log(x)*(2*x*(-3 + 2*log(mm/mmu)) + (-1 + x + y)*log(y))
                      + (1 + x + y)*(7 - 6*log(mm/mmu) + power2(Pi)/6. + 2*power2(log(mm/mmu)))
                      + x*power2(log(x)) + y*power2(log(y))))/2.;
      }

      return mm*((-7 + 6*log(mm/mmu) + log(x)*log(y)
                  - lambda*Phi(x,y) - power2(Pi)/6. - 2*power2(log(mm/mmu)))/2.
                 - (x*(7 - 6*log(mm/mmu) + log(x)*(-6 + 4*log(mm/mmu) + log(y))
                       + power2(Pi)/6. + 2*power2(log(mm/mmu)) + power2(log(x))))/2.
                 - (y*(7 - 6*log(mm/mmu) + (
                     -6 + 4*log(mm/mmu) + log(x))*log(y) + power2(Pi)/6.
                       + 2*power2(log(mm/mmu)) + power2(log(y))))/2.);
   }

   /// Delta[m1,m2,m3,-1]
   Real DeltaInv(Real m1, Real m2, Real m3)
   {
      return 1./(power2(m1) + power2(m2) + power2(m3) - 2*(m1*m2 + m1*m3 + m2*m3));
   }

   /// calculates sin(theta)
   Real calc_sin_theta(Real mf, Real xf, Real msf12, Real msf22)
   {
      if (is_zero(mf, 1e-10l) || is_zero(xf, 1e-10l))
         return 0.;

      const Real sin_2theta = 2.0*mf*xf / (msf12 - msf22);
      const Real theta = 0.5*std::asin(sin_2theta);

      return std::sin(theta);
   }

   /// calculates Higgs mixing angle from squarde Higgs masses and tan(beta)
   Real calc_alpha(Real mh2, Real mH2, Real tb)
   {
      const Real beta = std::atan(tb);
      const Real sin_2alpha = -(mH2 + mh2)/(mH2 - mh2) * std::sin(2*beta);

      return 0.5*std::asin(sin_2alpha);
   }

} // anonymous namespace

double delta_mtau_2loop_atau_atau(const Parameters& pars)
{
   using std::log;
   const Real ytau  = pars.ytau;
   const Real xtau  = pars.xtau;
   const Real mstau1 = pars.mstau1;
   const Real mstau12 = power2(pars.mstau1);
   const Real mstau2  = pars.mstau2;
   const Real mstau22 = power2(pars.mstau2);
   const Real mstau24 = power4(pars.mstau2);
   const Real msntau2 = power2(pars.msntau);
   const Real mw2   = power2(pars.mw);
   const Real mz2   = power2(pars.mz);
   const Real mh2   = power2(pars.mh);
   const Real mH2   = power2(pars.mH);
   const Real mC2   = power2(pars.mC);
   const Real mA2   = power2(pars.mA);
   const Real mu    = pars.mu;
   const Real mu2   = power2(pars.mu);
   const Real tb    = pars.tb;
   const Real sb    = tb / std::sqrt(1. + power2(tb));
   const Real cb    = 1. / std::sqrt(1. + power2(tb));
   const Real Q2    = power2(pars.Q);
   const Real alpha = calc_alpha(mh2, mH2, tb);
   const Real sa    = std::sin(alpha);
   const Real ca    = std::cos(alpha);

   const double result =
" <> WrapText @ IndentText[ToCPP[tau2lyl4] <> ";"] <> "

   return result * power4(ytau) * twoLoop;
}

double delta_mtau_2loop_atau_at(const Parameters& pars)
{
   using std::log;
   const Real ytau  = pars.ytau;
   const Real yt    = pars.yt;
   const Real xt    = pars.xt;
   const Real xtau  = pars.xtau;
   const Real mt    = pars.mt;
   const Real mt2   = power2(pars.mt);
   const Real mst1  = pars.mst1;
   const Real mst12 = power2(pars.mst1);
   const Real mst14 = power4(pars.mst1);
   const Real mst2  = pars.mst2;
   const Real mst22 = power2(pars.mst2);
   const Real mst24 = power4(pars.mst2);
   const Real msb1  = pars.msb1;
   const Real msb12 = power2(pars.msb1);
   const Real msb14 = power4(pars.msb1);
   const Real mstau12 = power2(pars.mstau1);
   const Real mstau22 = power2(pars.mstau2);
   const Real msntau2 = power2(pars.msntau);
   const Real mh2   = power2(pars.mh);
   const Real mH2   = power2(pars.mH);
   const Real mC2   = power2(pars.mC);
   const Real mA2   = power2(pars.mA);
   const Real mu    = pars.mu;
   const Real mu2   = power2(pars.mu);
   const Real tb    = pars.tb;
   const Real sb    = tb / std::sqrt(1. + power2(tb));
   const Real cb    = 1. / std::sqrt(1. + power2(tb));
   const Real Q2    = power2(pars.Q);
   const Real snt   = calc_sin_theta(mt, xt, mst12, mst22);
   const Real alpha = calc_alpha(mh2, mH2, tb);
   const Real sa    = std::sin(alpha);
   const Real ca    = std::cos(alpha);

   const double result =
" <> WrapText @ IndentText[ToCPP[tau2lyl2yt2] <> ";"] <> "

   return result * power2(ytau) * power2(yt) * twoLoop;
}

double delta_mtau_2loop_atau_ab(const Parameters& pars)
{
   using std::log;
   const Real ytau  = pars.ytau;
   const Real yb    = pars.yb;
   const Real xt    = pars.xt;
   const Real xb    = pars.xb;
   const Real mt    = pars.mt;
   const Real mt2   = power2(pars.mt);
   const Real mst1  = pars.mst1;
   const Real mst12 = power2(pars.mst1);
   const Real mst14 = power4(pars.mst1);
   const Real mst2  = pars.mst2;
   const Real mst22 = power2(pars.mst2);
   const Real mst24 = power4(pars.mst2);
   const Real msb1  = pars.msb1;
   const Real msb12 = power2(pars.msb1);
   const Real msb14 = power4(pars.msb1);
   const Real msb2  = pars.msb2;
   const Real msb22 = power2(pars.msb2);
   const Real msb24 = power4(pars.msb2);
   const Real mstau12 = power2(pars.mstau1);
   const Real mstau22 = power2(pars.mstau2);
   const Real msntau2 = power2(pars.msntau);
   const Real msntau4 = power4(pars.msntau);
   const Real mh2   = power2(pars.mh);
   const Real mH2   = power2(pars.mH);
   const Real mC2   = power2(pars.mC);
   const Real mA2   = power2(pars.mA);
   const Real mu    = pars.mu;
   const Real mu2   = power2(pars.mu);
   const Real tb    = pars.tb;
   const Real sb    = tb / std::sqrt(1. + power2(tb));
   const Real cb    = 1. / std::sqrt(1. + power2(tb));
   const Real Q2    = power2(pars.Q);
   const Real snt   = calc_sin_theta(mt, xt, mst12, mst22);
   const Real alpha = calc_alpha(mh2, mH2, tb);
   const Real sa    = std::sin(alpha);
   const Real ca    = std::cos(alpha);

   const double result =
" <> WrapText @ IndentText[ToCPP[tau2lyl2yb2] <> ";"] <> "

   return result * power2(ytau) * power2(yb) * twoLoop;
}

std::ostream& operator<<(std::ostream& out, const Parameters& pars)
{
   out <<
      \"Delta m_tau 2L parameters:\\n\"
      \"yt     = \" << pars.yt     << '\\n' <<
      \"yb     = \" << pars.yb     << '\\n' <<
      \"ytau   = \" << pars.ytau   << '\\n' <<
      \"mt     = \" << pars.mt     << '\\n' <<
      \"mb     = \" << pars.mb     << '\\n' <<
      \"mtau   = \" << pars.mtau   << '\\n' <<
      \"mst1   = \" << pars.mst1   << '\\n' <<
      \"mst2   = \" << pars.mst2   << '\\n' <<
      \"msb1   = \" << pars.msb1   << '\\n' <<
      \"msb2   = \" << pars.msb2   << '\\n' <<
      \"mstau1 = \" << pars.mstau1 << '\\n' <<
      \"mstau2 = \" << pars.mstau2 << '\\n' <<
      \"msntau = \" << pars.msntau << '\\n' <<
      \"xt     = \" << pars.xt     << '\\n' <<
      \"xb     = \" << pars.xb     << '\\n' <<
      \"xtau   = \" << pars.xtau   << '\\n' <<
      \"mw     = \" << pars.mw     << '\\n' <<
      \"mz     = \" << pars.mz     << '\\n' <<
      \"mh     = \" << pars.mh     << '\\n' <<
      \"mH     = \" << pars.mH     << '\\n' <<
      \"mC     = \" << pars.mC     << '\\n' <<
      \"mA     = \" << pars.mA     << '\\n' <<
      \"mu     = \" << pars.mu     << '\\n' <<
      \"tb     = \" << pars.tb     << '\\n' <<
      \"Q      = \" << pars.Q      << '\\n';

   return out;
}

} // namespace mssm_twoloop_mtau
} // namespace flexiblesusy
";

Export[headerName, header, "String"];
Export[implName  , impl  , "String"];
