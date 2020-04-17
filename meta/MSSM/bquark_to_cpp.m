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

b2l = Get[FileNameJoin[{"meta", "MSSM", "bquark_2loop_sqcd_decoupling.m"}]];

colorCF = 4/3; colorCA = 3;
GS = g3 = 1; scale = Q;
MT = mt; MB = mb; MGl = mgl;

mmgl  /: mmgl^2  := mmgl2;
mmgl  /: mmgl^3  := mmgl3;
mmsb1 /: mmsb1^2 := mmsb12;
mmsb2 /: mmsb2^2 := mmsb22;

Delta[m1_, m2_, m3_, -1] := DeltaInv[m1,m2,m3];
fin[0, m1_, m2_]         := Fin20[m1,m2,mmu];
fin[m1_, m2_, m3_]       := Fin3[m1,m2,m3,mmu];

Simp[expr_] := Collect[expr //. {
       mb^2                  -> mmb,
       1/mb^2                -> 1/mmb,
       (-mmgl + mmsb1)^n_:-1 /; n < 0 :> (invdgb1)^(-n),
       (-mmgl + mmsb2)^n_:-1 /; n < 0 :> (invdgb2)^(-n),
       (mmsb1 - mmsb2)^n_:-1 /; n < 0 :> (invdb12)^(-n),
       Power[x_,n_] /; n > 0 :> Symbol["pow" <> ToString[n]][x],
       Power[x_,-2]          :> 1/Symbol["pow2"][x],
       Power[x_,-3]          :> 1/Symbol["pow3"][x],
       Power[x_,-4]          :> 1/Symbol["pow4"][x],
       Power[x_,-5]          :> 1/Symbol["pow5"][x],
       Power[x_,-6]          :> 1/Symbol["pow6"][x],
       Log[mmgl/mmu]         -> lgu,
       Log[mmst1/mmu]        -> lt1u,
       Log[mmst2/mmu]        -> lt2u,
       Log[mmsb1/mmu]        -> lb1u,
       Log[mmsb2/mmu]        -> lb2u,
       Log[mmsusy/mmu]       -> lsu,
       Log[mmt/mmu]          -> ltu
    },
    { Fin20[__], Fin3[__], DeltaInv[__], s2t, s2b },
    Refine
];

ToCPP[expr_] := ToString[Simp[expr], CForm];

headerName = "mssm_twoloop_mb.hpp";
implName   = "mssm_twoloop_mb.cpp";

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
// with the script \"bquark_to_cpp.m\".

#ifndef MSSM_TWO_LOOP_SQCD_MB_H
#define MSSM_TWO_LOOP_SQCD_MB_H

namespace flexiblesusy {
namespace mssm_twoloop_mb {

struct Parameters {
    Parameters() = default;
    Parameters(double g3_, double mt_, double mb_, double mg_,
               double mst1_, double mst2_,
               double msb1_, double msb2_,
               double msusy_,
               double xt_, double xb_, double Q_)
       : g3(g3_), mt(mt_), mb(mb_), mg(mg_)
       , mst1(mst1_), mst2(mst2_)
       , msb1(msb1_), msb2(msb2_)
       , msusy(msusy_)
       , xt(xt_), xb(xb_), Q(Q_)
       {}

    double g3{};    ///< MSSM strong gauge coupling DR-bar
    double mt{};    ///< MSSM top mass DR-bar
    double mb{};    ///< SM   bottom mass DR-bar
    double mg{};    ///< MSSM gluino mass DR-bar
    double mst1{};  ///< MSSM light stop mass DR-bar
    double mst2{};  ///< MSSM heavy stop mass DR-bar
    double msb1{};  ///< MSSM light sbottom mass DR-bar
    double msb2{};  ///< MSSM heavy sbottom mass DR-bar
    double msusy{}; ///< MSSM light squark masses DR-bar
    double xt{};    ///< MSSM sbottom mixing parameter DR-bar
    double xb{};    ///< MSSM stop mixing parameter DR-bar
    double Q{};     ///< renormalization scale
};

/// 2-loop full SQCD contributions to mb [arXiv:0707.0650]
double delta_mb_2loop(const Parameters&);

} // namespace mssm_twoloop_mb
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
// with the script \"bquark_to_cpp.m\".

#include \"" <> headerName <> "\"
#include \"dilog.hpp\"
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <limits>

namespace flexiblesusy {
namespace mssm_twoloop_mb {

namespace {
   const double Pi  = 3.1415926535897932384626433832795;
   const double zt2 = 1.6449340668482264364724151666460;

   template <typename T> T pow2(T x) noexcept { return x*x; }
   template <typename T> T pow3(T x) noexcept { return x*x*x; }
   template <typename T> T pow4(T x) noexcept { return pow2(pow2(x)); }
   template <typename T> T pow5(T x) noexcept { return x*pow4(x); }

   const double oneLoop = 1/pow2(4*Pi);
   const double twoLoop = pow2(oneLoop);

   template <typename T>
   bool is_zero(T a, T prec = std::numeric_limits<T>::epsilon()) noexcept
   {
      return std::abs(a) < prec;
   }

   template <typename T>
   bool is_equal(T a, T b, T prec = std::numeric_limits<T>::epsilon()) noexcept
   {
      return is_zero(a - b, prec);
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
   double Fin20(double mm1, double mm2, double mmu) noexcept
   {
      const double log12 = std::log(mm1/mm2);
      const double log1u = std::log(mm1/mmu);
      const double log2u = std::log(mm2/mmu);

      return (6*(mm1*log1u + mm2*log2u) +
         (-mm1 - mm2)*(7 + pow2(Pi)/6) +
         (mm1 - mm2)*(2*dilog(1 - mm1/mm2) +
            pow2(log12)/2) +
         ((mm1 + mm2)*pow2(log12))/2 -
         2*(mm1*pow2(log1u) + mm2*pow2(log2u)))/2;
   }

   double LambdaSquared(double x, double y) noexcept
   {
      return pow2(1 - x - y) - 4*x*y;
   }

   /// ClausenCl[2,x]
   double ClausenCl2(double x) noexcept
   {
      const std::complex<double> img(0.0, 1.0);

      return std::imag(dilog(std::exp(img*x)));
   }

   /// x < 1 && y < 1, LambdaSquared(x,y) > 0
   double PhiPos(double x, double y) noexcept
   {
      const double lambda = std::sqrt(LambdaSquared(x,y));

      return (-(std::log(x)*std::log(y))
              + 2*std::log((1 - lambda + x - y)/2)*std::log((1 - lambda - x + y)/2)
              - 2*dilog((1 - lambda + x - y)/2)
              - 2*dilog((1 - lambda - x + y)/2)
              + pow2(Pi)/3)/lambda;
   }

   /// LambdaSquared(x,y) < 0
   double PhiNeg(double x, double y) noexcept
   {
      const double lambda = std::sqrt(-LambdaSquared(x,y));

      return 2*(+ ClausenCl2(2*std::acos((1 + x - y)/(2*std::sqrt(x))))
                + ClausenCl2(2*std::acos((1 - x + y)/(2*std::sqrt(y))))
                + ClausenCl2(2*std::acos((-1 + x + y)/(2*std::sqrt(x*y)))))/lambda;
   }

   double Phi(double x, double y) noexcept
   {
      const double lambda = LambdaSquared(x,y);

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
   double Fin3(double mm1, double mm2, double mm3, double mmu) noexcept
   {
      std::array<double,3> masses = { mm1, mm2, mm3 };
      std::sort(masses.begin(), masses.end());

      const double mm = masses[2];
      const double x = masses[0]/mm;
      const double y = masses[1]/mm;
      const double lambda = LambdaSquared(x,y);
      const double logx = std::log(x);
      const double logy = std::log(y);
      const double logm = std::log(mm/mmu);

      if (is_zero(lambda, 1e-10)) {
         return -(mm*(2*y*(-3 + 2*logm)*logy
                      + logx*(2*x*(-3 + 2*logm) + (-1 + x + y)*logy)
                      + (1 + x + y)*(7 - 6*logm + pow2(Pi)/6 + 2*pow2(logm))
                      + x*pow2(logx) + y*pow2(logy)))/2;
      }

      return mm*((-7 + 6*logm + logx*logy
                  - lambda*Phi(x,y) - pow2(Pi)/6 - 2*pow2(logm))/2
                 - (x*(7 - 6*logm + logx*(-6 + 4*logm + logy)
                       + pow2(Pi)/6 + 2*pow2(logm) + pow2(logx)))/2
                 - (y*(7 - 6*logm + (
                     -6 + 4*logm + logx)*logy + pow2(Pi)/6
                       + 2*pow2(logm) + pow2(logy)))/2);
   }

   /// Delta[m1,m2,m3,-1]
   double DeltaInv(double m1, double m2, double m3) noexcept
   {
      return 1/(pow2(m1) + pow2(m2) + pow2(m3) - 2*(m1*m2 + m1*m3 + m2*m3));
   }

} // anonymous namespace

/**
 * The function returns the 2-loop SQCD (QCD + SUSY) relation between
 * the Standard Model MS-bar bottom mass
 * \\f$m_b^{\\text{SM},\\overline{\\text{MS}}}\\f$ and the MSSM DR-bar
 * bottom mass \\f$m_b^{\\text{MSSM},\\overline{\\text{DR}}}\\f$,
 * Eq. (61) of [arxiv:0707.0650].  The relation has the form
 *
 * \\f[
    m_b^{\\text{SM},\\overline{\\text{MS}}} =
    m_b^{\\text{MSSM},\\overline{\\text{DR}}} \\left[
       1 + \\left(\\frac{\\Delta m_b}{m_b}\\right)_{1L}
         + \\left(\\frac{\\Delta m_b}{m_b}\\right)_{2L}
    \\right]
   \\f]
 *
 * The function returns \\f$(\\Delta m_b/m_b)_{2L}\\f$.
 */
double delta_mb_2loop(const Parameters& pars)
{
   const double g3     = pars.g3;
   const double Xt     = pars.xt;
   const double Xb     = pars.xb;
   const double mgl    = pars.mg;
   const double mt     = pars.mt;
   const double mb     = pars.mb;
   const double mmt    = pow2(pars.mt);
   const double mmb    = pow2(pars.mb);
   const double mmgl   = pow2(pars.mg);
   const double mmgl2  = pow2(mmgl);
   const double mmgl3  = mmgl*mmgl2;
   const double mmu    = pow2(pars.Q);
   const double mmst1  = pow2(pars.mst1);
   const double mmst2  = pow2(pars.mst2);
   const double mmsb1  = pow2(pars.msb1);
   const double mmsb2  = pow2(pars.msb2);
   const double mmsb12 = pow2(mmsb1);
   const double mmsb22 = pow2(mmsb2);
   const double mmsusy = pow2(pars.msusy);
   const double lgu    = std::log(mmgl/mmu);
   const double lt1u   = std::log(mmst1/mmu);
   const double lt2u   = std::log(mmst2/mmu);
   const double lb1u   = std::log(mmsb1/mmu);
   const double lb2u   = std::log(mmsb2/mmu);
   const double lsu    = std::log(mmsusy/mmu);
   const double ltu    = std::log(mmt/mmu);
   const double s2t    = 2*mt*Xt / (mmst1 - mmst2);
   const double s2b    = 2*mb*Xb / (mmsb1 - mmsb2);
   const double invdgb1 = 1/(-mmgl + mmsb1);
   const double invdgb2 = 1/(-mmgl + mmsb2);
   const double invdb12 = 1/(mmsb1 - mmsb2);

   const double result =
" <> WrapText @ IndentText[ToCPP[b2l] <> ";"] <> "

   return pow4(g3) * result * twoLoop;
}

} // namespace mssm_twoloop_mb
} // namespace flexiblesusy
";

Export[headerName, header, "String"];
Export[implName  , impl  , "String"];
