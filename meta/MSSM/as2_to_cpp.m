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

GS  = g3; mgl = mg; scale = Q;
YT  = yt; YB = yb;
MW  = mw; MZ = mz;
mh0 = mh; mH0 = mH; mA0 = mA; mHp = mC;
MUE = mu; SB = sb; CB = cb; SA = sa; CA = ca;

mw /: mw^n_ := mw2^(n/2) /; EvenQ[n];
mz /: mz^n_ := mz2^(n/2) /; EvenQ[n];
mt /: mt^n_ := mt2^(n/2) /; EvenQ[n];
mg /: mg^n_ := mg2^(n/2) /; EvenQ[n];
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
mst1 /: mst1^6 := mst16;
mst2 /: mst2^6 := mst26;
msb1 /: msb1^2 := msb12;
msb2 /: msb2^2 := msb22;
msb1 /: msb1^4 := msb14;
msb2 /: msb2^4 := msb24;
msb1 /: msb1^6 := msb16;
msb2 /: msb2^6 := msb26;
msd1 /: msd1^2 := msd12;
msd2 /: msd2^2 := msd22;
msd1 /: msd1^4 := msd14;
msd2 /: msd2^4 := msd24;
msd1 /: msd1^6 := msd16;
msd2 /: msd2^6 := msd26;

mg2  /: mg2^2  := mg4;
mg2  /: mg2^3  := mg6;

mst1 /: mst1^n_ := mst12^(n/2) /; EvenQ[n];
mst2 /: mst2^n_ := mst22^(n/2) /; EvenQ[n];
msb1 /: msb1^n_ := msb12^(n/2) /; EvenQ[n];
msb2 /: msb2^n_ := msb22^(n/2) /; EvenQ[n];

simpA = {
    At -> xt + MUE CB/SB,
    Ab -> xb + MUE SB/CB
};

a2l     = Get[FileNameJoin[{"meta", "MSSM", "das2.m"}]];
a2lsqcd = Coefficient[a2l, g3^4] /. simpA;
a2latas = Coefficient[a2l, g3^2 yt^2];
a2labas = Coefficient[a2l, g3^2 yb^2];

fpart[m1_, m2_, m3_, Q_]      := Fin3[m1^2, m2^2, m3^2, Q^2];
delta3[m1_, m2_, m3_]         := Delta[m1^2, m2^2, m3^2, -1]; 
Delta[m1_, m2_, m3_, -1]      := DeltaInv[m1,m2,m3];

Simp[expr_] :=
    Collect[expr //. {
                (mst12 - mst22)^n_:-1 /; n < 0 :> (invdmst)^(-n),

                (mst12 - mg2)^n_:-1 /; n < 0 :> ( invdmst1g)^(-n),
                (mg2 - mst12)^n_:-1 /; n < 0 :> (-invdmst1g)^(-n),
                (mst22 - mg2)^n_:-1 /; n < 0 :> ( invdmst2g)^(-n),
                (mg2 - mst22)^n_:-1 /; n < 0 :> (-invdmst2g)^(-n),

                (msb12 - mg2)^n_:-1 /; n < 0 :> ( invdmsb1g)^(-n),
                (mg2 - msb12)^n_:-1 /; n < 0 :> (-invdmsb1g)^(-n),
                (msb22 - mg2)^n_:-1 /; n < 0 :> ( invdmsb2g)^(-n),
                (mg2 - msb22)^n_:-1 /; n < 0 :> (-invdmsb2g)^(-n),

                (msd12 - mg2)^n_:-1 /; n < 0 :> ( invdmsd1g)^(-n),
                (mg2 - msd12)^n_:-1 /; n < 0 :> (-invdmsd1g)^(-n),
                (msd22 - mg2)^n_:-1 /; n < 0 :> ( invdmsd2g)^(-n),
                (mg2 - msd22)^n_:-1 /; n < 0 :> (-invdmsd2g)^(-n),

                (mC2  -  mt2)^n_:-1 /; n < 0 :> ( invdct)^(-n),
                (mt2  -  mC2)^n_:-1 /; n < 0 :> (-invdct)^(-n),

                (mt2  -  mw2)^n_:-1 /; n < 0 :> ( invdtw)^(-n),
                (mw2  -  mt2)^n_:-1 /; n < 0 :> (-invdtw)^(-n),

                Log[mst12/Q2] -> lmst12,
                Log[mst22/Q2] -> lmst22,
                Log[msb12/Q2] -> lmsb12,
                Log[msb22/Q2] -> lmsb22,
                Log[msd12/Q2] -> lmsd12,
                Log[msd22/Q2] -> lmsd22,
                Log[mg2/Q2]   -> lmg2,
                Log[mt2/Q2]   -> lmt2,
                Log[mw2/Q2]   -> lmw2,
                Log[mz2/Q2]   -> lmz2,
                Log[mu2/Q2]   -> lmu2,
                Log[mH2/Q2]   -> lmH2,
                Log[mA2/Q2]   -> lmA2,
                Log[mh2/Q2]   -> lmh2,
                Log[mC2/Q2]   -> lmC2
            },
            {
                Fin3[__],
                DeltaInv[__],
                xt,
                invdmst
            },
            Simplify[#, TimeConstraint -> 300]&
    ] //. {
        Power[x_,n_] /; n > 0 :> Symbol["power" <> ToString[n]][x],
        Power[x_,-2]          :> 1/Symbol["power2"][x],
        Power[x_,-3]          :> 1/Symbol["power3"][x],
        Power[x_,-4]          :> 1/Symbol["power4"][x],
        Power[x_,-5]          :> 1/Symbol["power5"][x],
        Power[x_,-6]          :> 1/Symbol["power6"][x],
        cb*mu - At*sb         -> -xt,
        -(cb*mu) + At*sb      -> xt
    };

ToCPP[expr_] := ToString[Simp[expr], CForm];

headerName = "mssm_twoloop_as.hpp";
implName   = "mssm_twoloop_as.cpp";

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
// with the script \"as2_to_cpp.m\".

#ifndef MSSM_TWO_LOOP_AS_H
#define MSSM_TWO_LOOP_AS_H

#include <iosfwd>

namespace flexiblesusy {
namespace mssm_twoloop_as {

using Real = long double;

struct Parameters {
    Real g3{};    ///< MSSM strong gauge coupling DR-bar
    Real yt{};    ///< MSSM top Yukawa coupling DR-bar
    Real yb{};    ///< MSSM bottom Yukawa coupling DR-bar
    Real mt{};    ///< MSSM top mass DR-bar
    Real mb{};    ///< MSSM bottom mass DR-bar
    Real mg{};    ///< MSSM gluino mass DR-bar
    Real mst1{};  ///< MSSM light stop mass DR-bar
    Real mst2{};  ///< MSSM heavy stop mass DR-bar
    Real msb1{};  ///< MSSM light sbottom mass DR-bar
    Real msb2{};  ///< MSSM heavy sbottom mass DR-bar
    Real msd1{};  ///< MSSM light sdown mass DR-bar
    Real msd2{};  ///< MSSM heavy sdown mass DR-bar
    Real xt{};    ///< MSSM stop mixing parameter DR-bar
    Real xb{};    ///< MSSM sbottom mixing parameter DR-bar
    Real mw{};    ///< MSSM W boson mass DR-bar
    Real mz{};    ///< MSSM Z boson mass DR-bar
    Real mh{};    ///< MSSM light CP-even Higgs mass DR-bar
    Real mH{};    ///< MSSM heavy CP-even Higgs mass DR-bar
    Real mC{};    ///< MSSM charged Higgs mass DR-bar
    Real mA{};    ///< MSSM CP-odd Higgs mass DR-bar
    Real mu{};    ///< MSSM mu superpotential parameter DR-bar
    Real tb{};    ///< MSSM tan(beta) DR-bar
    Real Q{};     ///< renormalization scale
};

/// 2-loop O(alpha_s^2) contributions to Delta g_3 [hep-ph/0509048,arXiv:0810.5101]
Real delta_alpha_s_2loop_as_as(const Parameters&);

/// 2-loop O(alpha_t*alpha_s) contributions to Delta g_3 [arXiv:1009.5455]
Real delta_alpha_s_2loop_at_as(const Parameters&);

/// 2-loop O(alpha_b*alpha_s) contributions to Delta g_3 [arXiv:1009.5455]
Real delta_alpha_s_2loop_ab_as(const Parameters&);

std::ostream& operator<<(std::ostream&, const Parameters&);

} // namespace mssm_twoloop_as
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
// with the script \"as2_to_cpp.m\".

#include \"" <> headerName <> "\"
#include \"dilog.hpp\"
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <ostream>

namespace flexiblesusy {
namespace mssm_twoloop_as {

namespace {
   const Real Pi = 3.1415926535897932384626433832795l;

   template <typename T> T power2(T x)  noexcept { return x*x; }
   template <typename T> T power3(T x)  noexcept { return x*x*x; }
   template <typename T> T power4(T x)  noexcept { return power2(power2(x)); }
   template <typename T> T power5(T x)  noexcept { return x*power4(x); }
   template <typename T> T power6(T x)  noexcept { return power2(power3(x)); }
   template <typename T> T power7(T x)  noexcept { return x*power6(x); }
   template <typename T> T power8(T x)  noexcept { return power2(power4(x)); }
   template <typename T> T power10(T x) noexcept { return power2(power5(x)); }
   template <typename T> T power12(T x) noexcept { return power3(power4(x)); }
   template <typename T> T power14(T x) noexcept { return power2(x)*power12(x); }

   const Real oneLoop = 1/power2(4*Pi);
   const Real twoLoop = power2(oneLoop);

   template <typename T>
   bool is_zero(T a, T prec = std::numeric_limits<T>::epsilon())
   {
      return std::abs(a) < prec;
   }

   template <typename T>
   bool is_equal(T a, T b, T prec = std::numeric_limits<T>::epsilon())
   {
      return is_zero(a - b, prec);
   }

   template <typename T>
   bool is_equal_rel(T a, T b, T prec = std::numeric_limits<T>::epsilon())
   {
      if (is_equal(a, b, std::numeric_limits<T>::epsilon())) {
         return true;
      }

      if (std::abs(a) < std::numeric_limits<T>::epsilon() ||
          std::abs(b) < std::numeric_limits<T>::epsilon()) {
         return false;
      }

      return std::abs((a - b)/a) < prec;
   }

   Real Log(Real x) noexcept
   {
      return std::log(x);
   }

   Real LambdaSquared(Real x, Real y) noexcept
   {
      return power2(1 - x - y) - 4*x*y;
   }

   /// ClausenCl[2,x]
   Real ClausenCl2(Real x) noexcept
   {
      const std::complex<Real> img(0.0l, 1.0l);

      return std::imag(dilog(std::exp(img*x)));
   }

   /// x < 1 && y < 1, LambdaSquared(x,y) > 0
   Real PhiPos(Real x, Real y) noexcept
   {
      const Real lambda = std::sqrt(LambdaSquared(x,y));

      return (-(std::log(x)*std::log(y))
              + 2*std::log((1 - lambda + x - y)/2)*std::log((1 - lambda - x + y)/2)
              - 2*dilog((1 - lambda + x - y)/2)
              - 2*dilog((1 - lambda - x + y)/2)
              + power2(Pi)/3)/lambda;
   }

   /// LambdaSquared(x,y) < 0
   Real PhiNeg(Real x, Real y) noexcept
   {
      const Real lambda = std::sqrt(-LambdaSquared(x,y));

      return 2*(+ ClausenCl2(2*std::acos((1 + x - y)/(2*std::sqrt(x))))
                + ClausenCl2(2*std::acos((1 - x + y)/(2*std::sqrt(y))))
                + ClausenCl2(2*std::acos((-1 + x + y)/(2*std::sqrt(x*y)))))/lambda;
   }

   Real Phi(Real x, Real y) noexcept
   {
      const Real lambda = LambdaSquared(x,y);

      if (lambda > 0) {
         return PhiPos(x,y);
      }

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
   Real Fin3(Real mm1, Real mm2, Real mm3, Real mmu) noexcept
   {
      std::array<Real,3> masses = { mm1, mm2, mm3 };
      std::sort(masses.begin(), masses.end());

      const Real mm = masses[2];
      const Real x = masses[0]/mm;
      const Real y = masses[1]/mm;
      const Real lambda = LambdaSquared(x,y);
      const Real logx = std::log(x);
      const Real logy = std::log(y);
      const Real logm = std::log(mm/mmu);

      if (is_zero(lambda, 1e-10l)) {
         return -(mm*(2*y*(-3 + 2*logm)*logy
                      + logx*(2*x*(-3 + 2*logm) + (-1 + x + y)*logy)
                      + (1 + x + y)*(7 - 6*logm + power2(Pi)/6 + 2*power2(logm))
                      + x*power2(logx) + y*power2(logy)))/2;
      }

      return mm*((-7 + 6*logm + logx*logy
                  - lambda*Phi(x,y) - power2(Pi)/6 - 2*power2(logm))/2
                 - (x*(7 - 6*logm + logx*(-6 + 4*logm + logy)
                       + power2(Pi)/6 + 2*power2(logm) + power2(logx)))/2
                 - (y*(7 - 6*logm + (
                     -6 + 4*logm + logx)*logy + power2(Pi)/6
                       + 2*power2(logm) + power2(logy)))/2);
   }

   /// Delta[m1,m2,m3,-1]
   Real DeltaInv(Real m1, Real m2, Real m3) noexcept
   {
      return 1/(power2(m1) + power2(m2) + power2(m3) - 2*(m1*m2 + m1*m3 + m2*m3));
   }

   /// calculates sin(theta)
   Real calc_sin_theta(Real mf, Real xf, Real msf12, Real msf22) noexcept
   {
      if (is_zero(mf, 1e-10l) || is_zero(xf, 1e-10l)) {
         return 0;
      }

      const Real sin_2theta = 2*mf*xf / (msf12 - msf22);
      const Real theta = 0.5l*std::asin(sin_2theta);

      return std::sin(theta);
   }

   /// calculates Higgs mixing angle from squarde Higgs masses and tan(beta)
   Real calc_alpha(Real mh2, Real mH2, Real tb) noexcept
   {
      const Real beta = std::atan(tb);
      const Real sin_2alpha = -(mH2 + mh2)/(mH2 - mh2) * std::sin(2*beta);

      return 0.5l*std::asin(sin_2alpha);
   }

   /// shift gluino mass away from mst1 and mst2 if too close
   Real shift_mg(Real mg, Real mst1, Real mst2) noexcept
   {
      if (is_equal_rel(std::min(mst1, mst2), mg, 0.01l)) {
         return mg * 0.98l;
      }

      if (is_equal_rel(std::max(mst1, mst2), mg, 0.01l)) {
         return mg * 1.02l;
      }

      return mg;
   }

} // anonymous namespace

/**
 * @brief 2-loop \\f$O(\\alpha_s^2)\\f$ contributions to \\f$\\Delta g_3^{(2L)}\\f$ [hep-ph/0509048, arXiv:0810.5101]
 *
 * The function returns \\f$\\Delta g_3^{(2L)}\\f$ which appears in
 * conversion of the strong \\f$\\overline{\\text{MS}}\\f$ gauge coupling
 * of the Standard Model with 5 active quark flavours,
 * \\f$g_3^{\\text{SM(5)},\\overline{\\text{MS}}}\\f$, to the strong
 * \\f$\\overline{\\text{DR}}\\f$ gauge coupling of the full MSSM
 * \\f$g_3^{\\text{MSSM},\\overline{\\text{DR}}}\\f$,
 * \\f{align*}{
  g_3^{\\text{SM(5)},\\overline{\\text{MS}}} =
  g_3^{\\text{MSSM},\\overline{\\text{DR}}} \\left[
     1 + \\Delta g_3^{(1L)} + \\Delta g_3^{(2L)}
  \\right]
 * \\f}
 */
Real delta_alpha_s_2loop_as_as(const Parameters& pars)
{
   const Real g3        = pars.g3;
   const Real xt        = pars.xt;
   const Real xb        = pars.xb;
   const Real mt        = pars.mt;
   const Real mt2       = power2(pars.mt);
   const Real mb        = pars.mb;
   const Real mg        = shift_mg(pars.mg, pars.mst1, pars.mst2);
   const Real mg2       = power2(mg);
   const Real mg4       = power2(mg2);
   const Real mg6       = mg2*mg4;
   const Real mst12     = power2(pars.mst1);
   const Real mst14     = power2(mst12);
   const Real mst16     = mst12*mst14;
   const Real mst22     = power2(pars.mst2);
   const Real mst24     = power2(mst22);
   const Real mst26     = mst22*mst24;
   const Real msb12     = power2(pars.msb1);
   const Real msb22     = power2(pars.msb2);
   const Real msd12     = power2(pars.msd1);
   const Real msd14     = power2(msd12);
   const Real msd16     = msd12*msd14;
   const Real msd22     = power2(pars.msd2);
   const Real msd24     = power2(msd22);
   const Real msd26     = msd22*msd24;
   const Real Q2        = power2(pars.Q);
   const Real snt       = calc_sin_theta(mt, xt, mst12, mst22);
   const Real snb       = calc_sin_theta(mb, xb, msb12, msb22);
   const Real invdmst   = 1/(mst12 - mst22);
   const Real invdmsb1g = 1/(msb12 - mg2);
   const Real invdmsb2g = 1/(msb22 - mg2);
   const Real invdmsd1g = 1/(msd12 - mg2);
   const Real invdmsd2g = 1/(msd22 - mg2);
   const Real lmst12    = std::log(mst12/Q2);
   const Real lmst22    = std::log(mst22/Q2);
   const Real lmsb12    = std::log(msb12/Q2);
   const Real lmsb22    = std::log(msb22/Q2);
   const Real lmsd12    = std::log(msd12/Q2);
   const Real lmsd22    = std::log(msd22/Q2);
   const Real lmg2      = std::log(mg2/Q2);
   const Real lmt2      = std::log(mt2/Q2);

   const Real result =
" <> WrapText @ IndentText[ToCPP[a2lsqcd] <> ";"] <> "

   return power4(g3) * result * twoLoop;
}

/**
 * @brief 2-loop \\f$O(\\alpha_t\\alpha_s)\\f$ contributions to \\f$\\Delta g_3^{(2L)}\\f$ [arXiv:1009.5455]
 *
 * The function returns \\f$\\Delta g_3^{(2L)}\\f$ which appears in
 * conversion of the strong \\f$\\overline{\\text{MS}}\\f$ gauge coupling
 * of the Standard Model with 5 active quark flavours,
 * \\f$g_3^{\\text{SM(5)},\\overline{\\text{MS}}}\\f$, to the strong
 * \\f$\\overline{\\text{DR}}\\f$ gauge coupling of the full MSSM
 * \\f$g_3^{\\text{MSSM},\\overline{\\text{DR}}}\\f$,
 * \\f{align*}{
  g_3^{\\text{SM(5)},\\overline{\\text{MS}}} =
  g_3^{\\text{MSSM},\\overline{\\text{DR}}} \\left[
     1 + \\Delta g_3^{(1L)} + \\Delta g_3^{(2L)}
  \\right]
 * \\f}
 */
Real delta_alpha_s_2loop_at_as(const Parameters& pars)
{
   const Real g3        = pars.g3;
   const Real yt        = pars.yt;
   const Real xt        = pars.xt;
   const Real xb        = pars.xb;
   const Real mt        = pars.mt;
   const Real mt2       = power2(pars.mt);
   const Real mb        = pars.mb;
   const Real mg        = shift_mg(pars.mg, pars.mst1, pars.mst2);
   const Real mg2       = power2(mg);
   const Real mg4       = power2(mg2);
   const Real mg6       = mg2*mg4;
   const Real mst12     = power2(pars.mst1);
   const Real mst14     = power2(mst12);
   const Real mst16     = mst12*mst14;
   const Real mst22     = power2(pars.mst2);
   const Real mst24     = power2(mst22);
   const Real mst26     = mst22*mst24;
   const Real msb12     = power2(pars.msb1);
   const Real msb14     = power2(msb12);
   const Real msb16     = msb12*msb14;
   const Real msb22     = power2(pars.msb2);
   const Real msb24     = power2(msb22);
   const Real msb26     = msb22*msb24;
   const Real msd12     = power2(pars.msd1);
   const Real msd14     = power2(msd12);
   const Real msd16     = msd12*msd14;
   const Real msd22     = power2(pars.msd2);
   const Real msd24     = power2(msd22);
   const Real msd26     = msd22*msd24;
   const Real mw2       = power2(pars.mw);
   const Real mz2       = power2(pars.mz);
   const Real mh2       = power2(pars.mh);
   const Real mH2       = power2(pars.mH);
   const Real mC2       = power2(pars.mC);
   const Real mA2       = power2(pars.mA);
   const Real mu        = pars.mu;
   const Real mu2       = power2(pars.mu);
   const Real tb        = pars.tb;
   const Real sb        = tb / std::sqrt(1. + power2(tb));
   const Real cb        = 1. / std::sqrt(1. + power2(tb));
   const Real Q2        = power2(pars.Q);
   const Real snt       = calc_sin_theta(mt, xt, mst12, mst22);
   const Real snb       = calc_sin_theta(mb, xb, msb12, msb22);
   const Real alpha     = calc_alpha(mh2, mH2, tb);
   const Real sa        = std::sin(alpha);
   const Real ca        = std::cos(alpha);
   const Real At        = xt + mu*cb/sb;
   const Real Ab        = xb + mu*sb/cb;
   const Real invdmst   = 1/(mst12 - mst22);
   const Real invdct    = 1/(mC2 - mt2);
   const Real invdtw    = 1/(mt2 - mw2);
   const Real lmst12    = std::log(mst12/Q2);
   const Real lmst22    = std::log(mst22/Q2);
   const Real lmsb12    = std::log(msb12/Q2);
   const Real lmsb22    = std::log(msb22/Q2);
   const Real lmt2      = std::log(mt2/Q2);
   const Real lmw2      = std::log(mw2/Q2);
   const Real lmz2      = std::log(mz2/Q2);
   const Real lmu2      = std::log(mu2/Q2);
   const Real lmH2      = std::log(mH2/Q2);
   const Real lmA2      = std::log(mA2/Q2);
   const Real lmh2      = std::log(mh2/Q2);
   const Real lmC2      = std::log(mC2/Q2);

   const Real result =
" <> WrapText @ IndentText[ToCPP[a2latas] <> ";"] <> "

   return power2(g3) * power2(yt) * result * twoLoop;
}

/**
 * @brief 2-loop \\f$O(\\alpha_b\\alpha_s)\\f$ contributions to \\f$\\Delta g_3^{(2L)}\\f$ [arXiv:1009.5455]
 *
 * The function returns \\f$\\Delta g_3^{(2L)}\\f$ which appears in
 * conversion of the strong \\f$\\overline{\\text{MS}}\\f$ gauge coupling
 * of the Standard Model with 5 active quark flavours,
 * \\f$g_3^{\\text{SM(5)},\\overline{\\text{MS}}}\\f$, to the strong
 * \\f$\\overline{\\text{DR}}\\f$ gauge coupling of the full MSSM
 * \\f$g_3^{\\text{MSSM},\\overline{\\text{DR}}}\\f$,
 * \\f{align*}{
  g_3^{\\text{SM(5)},\\overline{\\text{MS}}} =
  g_3^{\\text{MSSM},\\overline{\\text{DR}}} \\left[
     1 + \\Delta g_3^{(1L)} + \\Delta g_3^{(2L)}
  \\right]
 * \\f}
 */
Real delta_alpha_s_2loop_ab_as(const Parameters& pars)
{
   const Real g3        = pars.g3;
   const Real yb        = pars.yb;
   const Real xt        = pars.xt;
   const Real xb        = pars.xb;
   const Real mt        = pars.mt;
   const Real mt2       = power2(pars.mt);
   const Real mb        = pars.mb;
   const Real mg        = shift_mg(pars.mg, pars.mst1, pars.mst2);
   const Real mg2       = power2(mg);
   const Real mg4       = power2(mg2);
   const Real mg6       = mg2*mg4;
   const Real mst12     = power2(pars.mst1);
   const Real mst14     = power2(mst12);
   const Real mst16     = mst12*mst14;
   const Real mst22     = power2(pars.mst2);
   const Real mst24     = power2(mst22);
   const Real mst26     = mst22*mst24;
   const Real msb12     = power2(pars.msb1);
   const Real msb14     = power2(msb12);
   const Real msb16     = msb12*msb14;
   const Real msb22     = power2(pars.msb2);
   const Real msb24     = power2(msb22);
   const Real msb26     = msb22*msb24;
   const Real msd12     = power2(pars.msd1);
   const Real msd14     = power2(msd12);
   const Real msd16     = msd12*msd14;
   const Real msd22     = power2(pars.msd2);
   const Real msd24     = power2(msd22);
   const Real msd26     = msd22*msd24;
   const Real mw2       = power2(pars.mw);
   const Real mz2       = power2(pars.mz);
   const Real mh2       = power2(pars.mh);
   const Real mH2       = power2(pars.mH);
   const Real mC2       = power2(pars.mC);
   const Real mA2       = power2(pars.mA);
   const Real mu        = pars.mu;
   const Real mu2       = power2(pars.mu);
   const Real tb        = pars.tb;
   const Real sb        = tb / std::sqrt(1. + power2(tb));
   const Real cb        = 1. / std::sqrt(1. + power2(tb));
   const Real Q2        = power2(pars.Q);
   const Real snt       = calc_sin_theta(mt, xt, mst12, mst22);
   const Real snb       = calc_sin_theta(mb, xb, msb12, msb22);
   const Real alpha     = calc_alpha(mh2, mH2, tb);
   const Real sa        = std::sin(alpha);
   const Real ca        = std::cos(alpha);
   const Real At        = xt + mu*cb/sb;
   const Real Ab        = xb + mu*sb/cb;
   const Real invdmst   = 1/(mst12 - mst22);
   const Real invdct    = 1/(mC2 - mt2);
   const Real invdtw    = 1/(mt2 - mw2);
   const Real lmst12    = std::log(mst12/Q2);
   const Real lmst22    = std::log(mst22/Q2);
   const Real lmsb12    = std::log(msb12/Q2);
   const Real lmsb22    = std::log(msb22/Q2);
   const Real lmt2      = std::log(mt2/Q2);
   const Real lmw2      = std::log(mw2/Q2);
   const Real lmz2      = std::log(mz2/Q2);
   const Real lmu2      = std::log(mu2/Q2);
   const Real lmH2      = std::log(mH2/Q2);
   const Real lmA2      = std::log(mA2/Q2);
   const Real lmh2      = std::log(mh2/Q2);
   const Real lmC2      = std::log(mC2/Q2);

   const Real result =
" <> WrapText @ IndentText[ToCPP[a2labas] <> ";"] <> "

   return power2(g3) * power2(yb) * result * twoLoop;
}

std::ostream& operator<<(std::ostream& out, const Parameters& pars)
{
   out <<
      \"Delta alpha_s 2L parameters:\\n\"
      \"g3   = \" <<  pars.g3   << '\\n' <<
      \"yt   = \" <<  pars.yt   << '\\n' <<
      \"yb   = \" <<  pars.yb   << '\\n' <<
      \"mt   = \" <<  pars.mt   << '\\n' <<
      \"mb   = \" <<  pars.mb   << '\\n' <<
      \"mg   = \" <<  pars.mg   << '\\n' <<
      \"mst1 = \" <<  pars.mst1 << '\\n' <<
      \"mst2 = \" <<  pars.mst2 << '\\n' <<
      \"msb1 = \" <<  pars.msb1 << '\\n' <<
      \"msb2 = \" <<  pars.msb2 << '\\n' <<
      \"msd1 = \" <<  pars.msd1 << '\\n' <<
      \"msd2 = \" <<  pars.msd2 << '\\n' <<
      \"xt   = \" <<  pars.xt   << '\\n' <<
      \"xb   = \" <<  pars.xb   << '\\n' <<
      \"mw   = \" <<  pars.mw   << '\\n' <<
      \"mz   = \" <<  pars.mz   << '\\n' <<
      \"mh   = \" <<  pars.mh   << '\\n' <<
      \"mH   = \" <<  pars.mH   << '\\n' <<
      \"mC   = \" <<  pars.mC   << '\\n' <<
      \"mA   = \" <<  pars.mA   << '\\n' <<
      \"mu   = \" <<  pars.mu   << '\\n' <<
      \"tb   = \" <<  pars.tb   << '\\n' <<
      \"Q    = \" <<  pars.Q    << '\\n';

   return out;
}

} // namespace mssm_twoloop_as
} // namespace flexiblesusy
";

Export[headerName, header, "String"];
Export[implName  , impl  , "String"];
