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
GS = g3; scale = Q;
MT = mt; MB = mb; MGl = mgl;
At = Xt + MUE CB / SB;
Ab = Xb + MUE SB / CB;
s2t = 2 mt Xt / (mmst1 - mmst2);
s2b = 2 mb Xb / (mmsb1 - mmsb2);

Delta[m1_, m2_, m3_, -1] := DeltaInv[m1,m2,m3];
fin[0, m1_, m2_]         := Fin20[m1,m2,mmu];
fin[m1_, m2_, m3_]       := Fin3[m1,m2,m3,mmu];

Simp[expr_] := Collect[expr, { g3, Xt, Xb }] //. {
        mb^2                  -> mmb,
        1/mb^2                -> 1/mmb,
        Power[x_,n_] /; n > 0 :> Symbol["pow" <> ToString[n]][x],
        Power[x_,-2]          :> 1/Symbol["pow" <> ToString[2]][x],
        Power[x_,-3]          :> 1/Symbol["pow" <> ToString[3]][x],
        Power[x_,-4]          :> 1/Symbol["pow" <> ToString[4]][x],
        Power[x_,-5]          :> 1/Symbol["pow" <> ToString[5]][x],
        Power[x_,-6]          :> 1/Symbol["pow" <> ToString[6]][x],
        Log[x_]               :> log[x],
        PolyLog[2,x_]         :> dilog[x]
    };

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

   template <typename T> T pow2(T x)  { return x*x; }
   template <typename T> T pow3(T x)  { return x*x*x; }
   template <typename T> T pow4(T x)  { return x*x*x*x; }
   template <typename T> T pow5(T x)  { return x*x*x*x*x; }

   const double oneLoop = 1./pow2(4*Pi);
   const double twoLoop = pow2(oneLoop);

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

   /**
    * Fin20[] function from twoloopbubble.m .
    *
    * @param mm1 squared mass \\f$m_1^2\\f$
    * @param mm2 squared mass \\f$m_2^2\\f$
    * @param mmu squared renormalization scale
    *
    * @return Fin20(m12, m22, mmu)
    */
   double Fin20(double mm1, double mm2, double mmu)
   {
      using std::log;
      const double PI = 3.14159265358979323846264338327950288;

      return (6*(mm1*log(mm1/mmu) + mm2*log(mm2/mmu)) +
         (-mm1 - mm2)*(7 + pow2(PI)/6.) +
         (mm1 - mm2)*(2*dilog(1 - mm1/mm2) +
            pow2(log(mm1/mm2))/2.) +
         ((mm1 + mm2)*pow2(log(mm1/mm2)))/2. -
         2*(mm1*pow2(log(mm1/mmu)) + mm2*pow2(log(mm2/mmu))))/2.;
   }

   double LambdaSquared(double x, double y)
   {
      return pow2(1 - x - y) - 4*x*y;
   }

   /// ClausenCl[2,x]
   double ClausenCl2(double x)
   {
      using std::exp;
      const std::complex<double> img(0.,1.);

      return std::imag(dilog(exp(img*x)));
   }

   /// x < 1 && y < 1, LambdaSquared(x,y) > 0
   double PhiPos(double x, double y)
   {
      const double lambda = std::sqrt(LambdaSquared(x,y));

      return (-(log(x)*log(y))
              + 2*log((1 - lambda + x - y)/2.)*log((1 - lambda - x + y)/2.)
              - 2*dilog((1 - lambda + x - y)/2.)
              - 2*dilog((1 - lambda - x + y)/2.)
              + pow2(Pi)/3.)/lambda;
   }

   /// LambdaSquared(x,y) < 0
   double PhiNeg(double x, double y)
   {
      using std::acos;
      using std::sqrt;
      const double lambda = std::sqrt(-LambdaSquared(x,y));

      return 2*(+ ClausenCl2(2*acos((1 + x - y)/(2.*sqrt(x))))
                + ClausenCl2(2*acos((1 - x + y)/(2.*sqrt(y))))
                + ClausenCl2(2*acos((-1 + x + y)/(2.*sqrt(x*y)))))/lambda;
   }

   double Phi(double x, double y)
   {
      const double lambda = LambdaSquared(x,y);

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
   double Fin3(double mm1, double mm2, double mm3, double mmu)
   {
      using std::log;

      std::array<double,3> masses = { mm1, mm2, mm3 };
      std::sort(masses.begin(), masses.end());

      const double mm = masses[2];
      const double x = masses[0]/mm;
      const double y = masses[1]/mm;

      const double lambda = LambdaSquared(x,y);

      if (is_zero(lambda, 1e-10)) {
         return -(mm*(2*y*(-3 + 2*log(mm/mmu))*log(y)
                      + log(x)*(2*x*(-3 + 2*log(mm/mmu)) + (-1 + x + y)*log(y))
                      + (1 + x + y)*(7 - 6*log(mm/mmu) + pow2(Pi)/6. + 2*pow2(log(mm/mmu)))
                      + x*pow2(log(x)) + y*pow2(log(y))))/2.;
      }

      return mm*((-7 + 6*log(mm/mmu) + log(x)*log(y)
                  - lambda*Phi(x,y) - pow2(Pi)/6. - 2*pow2(log(mm/mmu)))/2.
                 - (x*(7 - 6*log(mm/mmu) + log(x)*(-6 + 4*log(mm/mmu) + log(y))
                       + pow2(Pi)/6. + 2*pow2(log(mm/mmu)) + pow2(log(x))))/2.
                 - (y*(7 - 6*log(mm/mmu) + (
                     -6 + 4*log(mm/mmu) + log(x))*log(y) + pow2(Pi)/6.
                       + 2*pow2(log(mm/mmu)) + pow2(log(y))))/2.);
   }

   /// Delta[m1,m2,m3,-1]
   double DeltaInv(double m1, double m2, double m3)
   {
      return 1./(pow2(m1) + pow2(m2) + pow2(m3) - 2*(m1*m2 + m1*m3 + m2*m3));
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
   using std::log;
   const double g3     = pars.g3;
   const double Xt     = pars.xt;
   const double Xb     = pars.xb;
   const double mgl    = pars.mg;
   const double mmt    = pow2(pars.mt);
   const double mmb    = pow2(pars.mb);
   const double mmgl   = pow2(pars.mg);
   const double mmu    = pow2(pars.Q);
   const double mmst1  = pow2(pars.mst1);
   const double mmst2  = pow2(pars.mst2);
   const double mmsb1  = pow2(pars.msb1);
   const double mmsb2  = pow2(pars.msb2);
   const double mmsusy = pow2(pars.msusy);

   const double result =
" <> WrapText @ IndentText[ToCPP[b2l] <> ";"] <> "

   return result * twoLoop;
}

} // namespace mssm_twoloop_mb
} // namespace flexiblesusy
";

Export[headerName, header, "String"];
Export[implName  , impl  , "String"];
