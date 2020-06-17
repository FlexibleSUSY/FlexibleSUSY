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

t1l    = Get[FileNameJoin[{"meta", "MSSM", "tquark_1loop_strong.m"}]];
t1lqcd = Get[FileNameJoin[{"meta", "MSSM", "tquark_1loop_qcd.m"}]] /. GSY -> GS;
t2l    = Get[FileNameJoin[{"meta", "MSSM", "tquark_2loop_strong.m"}]];
t2lqcd = Get[FileNameJoin[{"meta", "MSSM", "tquark_2loop_qcd.m"}]] /. GSY -> GS;

colorCA = 3; colorCF = 4/3; Tf = 1/2; GS = g3; g3 = 1;
MGl = mgl; MT = mt;
(* SX = 2 mt Xt; s2t = SX / (mmst1 - mmst2); *)
fin[0, args__] := fin[args, mmu];

Simp[expr_] := Simplify[expr] //.
    {
        Power[x_,n_] /; n > 0 :> Symbol["pow" <> ToString[n]][x],
        Power[x_,-2]          :> 1/Symbol["pow" <> ToString[2]][x],
        Power[x_,-3]          :> 1/Symbol["pow" <> ToString[3]][x],
        Power[x_,-4]          :> 1/Symbol["pow" <> ToString[4]][x],
        Power[x_,-5]          :> 1/Symbol["pow" <> ToString[5]][x],
        Power[x_,-6]          :> 1/Symbol["pow" <> ToString[6]][x],
        Log[x_/y_]            :> Symbol["log" <> ToString[x] <> ToString[y]]
    };

ToCPP[expr_] := ToString[Simp[expr], CForm];

(* ******* calculate limits ******* *)

CollectTerms[expr_] := Collect[expr /. zt2 -> Zeta[2], {Log[__]}, Together];

loopFunctions = {
    Hmine[mm1_,mm2_] :> 2 * PolyLog[2, 1-mm1/mm2] + 1/2 * Log[mm1/mm2]^2,
    fin[mm1_,mm2_,mmu_] :>
       1/2 * ( - (mm1 + mm2) * ( 7 + Zeta[2] )
               + 6 * (mm1 * Log[mm1/mmu] + mm2 * Log[mm2/mmu])
               - 2 * (mm1 * Log[mm1/mmu]^2 + mm2 * Log[mm2/mmu]^2 )
               +1/2 * (mm1 + mm2) * Log[mm1/mm2]^2 + (mm1-mm2)*Hmine[mm1,mm2] )
};

t1lLimitS1S2MG = Normal[Series[(t1l - t1lqcd) //. {
    mmst1 -> mmgl + x * dst1,
    mmst2 -> mmgl + x * dst2
    }, {x,0,0}]];
t1lLimitS1S2 = Normal[Series[(t1l - t1lqcd) //. {
    mmst1 -> mmst2 + x * dst1
    }, {x,0,0}]];
t1lLimitS1MG = Normal[Series[(t1l - t1lqcd) /. mmst1 -> mmgl + x * dst1, {x,0,0}]];
t1lLimitS2MG = Normal[Series[(t1l - t1lqcd) /. mmst2 -> mmgl + x * dst1, {x,0,0}]];

t2lLimitS1S2MSMG = CollectTerms @ Normal[Series[(t2l - t2lqcd) //. loopFunctions //. {
    mmst1  -> mmsusy + x * dst1,
    mmst2  -> mmsusy + x * dst2,
    mmgl   -> mmsusy + x * dst3
    }, {x,0,0}]];

(* t2lLimitS1S2MG = CollectTerms @ Normal[Series[(t2l - t2lqcd) //. loopFunctions //. { *)
(*     mmst1 -> mmgl + x * dst1, *)
(*     mmst2 -> mmgl + x * dst2 *)
(*     }, {x,0,0}]]; *)
(* t2lLimitS1MSMG = CollectTerms @ Normal[Series[(t2l - t2lqcd) //. loopFunctions //. { *)
(*     mmst1  -> mmsusy + x * dst1, *)
(*     mmgl   -> mmsusy + x * dst2 *)
(*     }, {x,0,0}]]; *)
(* t2lLimitS2MSMG = CollectTerms @ Normal[Series[(t2l - t2lqcd) //. loopFunctions //. { *)
(*     mmst2  -> mmsusy + x * dst1, *)
(*     mmgl   -> mmsusy + x * dst2 *)
(*     }, {x,0,0}]]; *)
(* t2lLimitS1S2MS = CollectTerms @ Normal[Series[(t2l - t2lqcd) //. loopFunctions //. { *)
(*     mmst1 -> mmsusy + x * dst1, *)
(*     mmst2 -> mmsusy + x * dst2 *)
(*     }, {x,0,0}]]; *)
(* t2lLimitS1S2 = CollectTerms @ Normal[Series[(t2l - t2lqcd) //. loopFunctions //. { *)
(*     mmst1 -> mmst2 + x * dst1 *)
(*     }, {x,0,0}]]; *)
(* t2lLimitS1MG = CollectTerms @ Normal[Series[(t2l - t2lqcd) //. loopFunctions /. mmst1  -> mmgl   + x * dst1, {x,0,0}]]; *)
(* t2lLimitS1MS = CollectTerms @ Normal[Series[(t2l - t2lqcd) //. loopFunctions /. mmst1  -> mmsusy + x * dst1, {x,0,0}]]; *)
(* t2lLimitS2MG = CollectTerms @ Normal[Series[(t2l - t2lqcd) //. loopFunctions /. mmst2  -> mmgl   + x * dst1, {x,0,0}]]; *)
(* t2lLimitS2MS = CollectTerms @ Normal[Series[(t2l - t2lqcd) //. loopFunctions /. mmst2  -> mmsusy + x * dst1, {x,0,0}]]; *)
(* t2lLimitMSMG = CollectTerms @ Normal[Series[(t2l - t2lqcd) //. loopFunctions /. mmsusy -> mmgl   + x * dst1, {x,0,0}]]; *)

mt /: mt^2 = mmt;

headerName = "mssm_twoloop_mt.hpp";
implName   = "mssm_twoloop_mt.cpp";

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
// with the script \"tquark_to_cpp.m\".

#ifndef MSSM_TWO_LOOP_SQCD_MT_H
#define MSSM_TWO_LOOP_SQCD_MT_H

namespace flexiblesusy {
namespace mssm_twoloop_mt {

struct Parameters {
    Parameters() = default;
    Parameters(double g3_, double mt_, double mg_, double mst1_,
               double mst2_, double msusy_, double xt_, double Q_)
       : g3(g3_), mt(mt_), mg(mg_), mst1(mst1_)
       , mst2(mst2_), msusy(msusy_), xt(xt_), Q(Q_)
       {}

    double g3{};    ///< MSSM strong gauge coupling DR-bar
    double mt{};    ///< MSSM top mass DR-bar
    double mg{};    ///< MSSM gluino mass DR-bar
    double mst1{};  ///< MSSM light stop mass DR-bar
    double mst2{};  ///< MSSM heavy stop mass DR-bar
    double msusy{}; ///< MSSM SUSY particle mass scale
    double xt{};    ///< MSSM stop mixing parameter DR-bar
    double Q{};     ///< renormalization scale
};

/// 1-loop QCD contributions to Delta Mt over mt [hep-ph/0210258]
double dMt_over_mt_1loop_qcd(const Parameters&);
/// 1-loop SUSY contributions to Delta Mt over mt [hep-ph/0210258]
double dMt_over_mt_1loop_susy(const Parameters&);
/// 1-loop full SQCD contributions to Delta Mt over mt [hep-ph/0210258]
double dMt_over_mt_1loop(const Parameters&);

/// 2-loop QCD contributions to Delta Mt over mt [hep-ph/0507139]
double dMt_over_mt_2loop_qcd(const Parameters&);
/// 2-loop SUSY contributions to Delta Mt over mt [hep-ph/0507139]
double dMt_over_mt_2loop_susy(const Parameters&);
/// 2-loop full SQCD contributions to Delta Mt over mt [hep-ph/0507139]
double dMt_over_mt_2loop(const Parameters&);

} // namespace mssm_twoloop_mt
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
// with the script \"tquark_to_cpp.m\".

#include \"" <> headerName <> "\"
#include \"dilog.hpp\"
#include <cmath>
#include <limits>

namespace flexiblesusy {
namespace mssm_twoloop_mt {

namespace {
   constexpr double Pi      = 3.1415926535897932384626433832795;
   constexpr double zt2     = 1.6449340668482264364724151666460;    // Zeta[2]
   constexpr double zt3     = 1.2020569031595942853997381615114;    // Zeta[3]
   constexpr double log2    = 6.9314718055994530941723212145818e-1; // Log[2]
   constexpr double oneLoop = 6.3325739776461107152424664506080e-3; // 1/(4Pi)^2
   constexpr double twoLoop = 4.0101493182360684332628059637182e-5; // 1/(4Pi)^4

   constexpr double pow2(double x) noexcept { return x*x; }
   constexpr double pow3(double x) noexcept { return x*x*x; }
   constexpr double pow4(double x) noexcept { return pow2(pow2(x)); }
   constexpr double pow5(double x) noexcept { return x*pow4(x); }
   constexpr double pow6(double x) noexcept { return pow2(pow3(x)); }
   constexpr double pow7(double x) noexcept { return x*pow6(x); }
   constexpr double pow8(double x) noexcept { return pow2(pow4(x)); }

   bool is_zero(double a, double prec) noexcept
   {
      return std::abs(a) < prec;
   }

   bool is_equal(double a, double b, double prec) noexcept
   {
      return is_zero(a - b, prec);
   }

   bool is_equal_rel(double a, double b, double prec) noexcept
   {
      if (is_equal(a, b, std::numeric_limits<double>::epsilon()))
         return true;

      const double min = std::min(std::abs(a), std::abs(b));

      if (min < std::numeric_limits<double>::epsilon())
         return is_equal(a, b, prec);

      const double max = std::max(std::abs(a), std::abs(b));

      return is_equal(a, b, prec*max);
   }

   /**
    * fin[] function from arXiv:hep-ph/0507139 .
    *
    * @param mm1 squared mass \\f$m_1^2\\f$
    * @param mm2 squared mass \\f$m_2^2\\f$
    * @param mmu squared renormalization scale
    *
    * @return fin(m12, m22)
    */
   double fin(double mm1, double mm2, double mmu) noexcept
   {
      const double log1u = std::log(mm1/mmu);
      const double log2u = std::log(mm2/mmu);
      const double log12 = std::log(mm1/mm2);

      return (6*(mm1*log1u + mm2*log2u) +
         (-mm1 - mm2)*(7 + pow2(Pi)/6.) +
         (mm1 - mm2)*(2*dilog(1 - mm1/mm2) + pow2(log12)/2.) +
         ((mm1 + mm2)*pow2(log12))/2. -
         2*(mm1*pow2(log1u) + mm2*pow2(log2u)))/2.;
   }

   /// shift gluino mass away from mst1 and mst2 if too close
   double shift_mg(double mg, double mst1, double mst2) noexcept
   {
      if (is_equal_rel(std::min(mst1, mst2), mg, 0.0003))
         return mg * 0.9995;

      if (is_equal_rel(std::max(mst1, mst2), mg, 0.0003))
         return mg * 1.0005;

      return mg;
   }

} // anonymous namespace

/// 1-loop QCD contributions to Delta Mt over mt [hep-ph/0210258]
double dMt_over_mt_1loop_qcd(const Parameters& pars)
{
   const double g32 = pow2(pars.g3);
   const double mmt = pow2(pars.mt);
   const double mmu = pow2(pars.Q);
   const double logmmtmmu = std::log(mmt/mmu);

   const double result = " <> ToCPP[t1lqcd] <> ";

   return result * g32 * oneLoop;
}

/// 1-loop SUSY contributions to Delta Mt over mt [hep-ph/0210258]
double dMt_over_mt_1loop_susy(const Parameters& pars)
{
   const double g32    = pow2(pars.g3);
   const double mt     = pars.mt;
   const double mgl    = pars.mg;
   const double mmu    = pow2(pars.Q);
   const double mmgl   = pow2(pars.mg);
   const double mmst1  = pow2(pars.mst1);
   const double mmst2  = pow2(pars.mst2);
   const double SX     = 2*mt*pars.xt;
   const double s2t    = SX / (mmst1 - mmst2);

   if (is_equal(mmst1, mmst2, 1e-6) && is_equal(mmst1, mmgl, 1e-6)) {
      const double logmmglmmu = std::log(mmgl/mmu);
      const double result = " <> ToCPP[t1lLimitS1S2MG] <> ";

      return result * g32 * oneLoop;
   }

   if (is_equal(mmst1, mmst2, 1e-6)) {
      const double logmmglmmu  = std::log(mmgl/mmu);
      const double logmmst2mmu = std::log(mmst2/mmu);

      const double result =
" <> WrapLines @ IndentText @ IndentText[ToCPP[t1lLimitS1S2] <> ";"] <> "

      return result * g32 * oneLoop;
   }

   if (is_equal(mmgl, mmst1, 1e-6)) {
      const double logmmglmmu  = std::log(mmgl/mmu);
      const double logmmst2mmu = std::log(mmst2/mmu);

      const double result =
" <> WrapLines @ IndentText @ IndentText[ToCPP[t1lLimitS1MG] <> ";"] <> "

      return result * g32 * oneLoop;
   }

   if (is_equal(mmgl, mmst2, 1e-6)) {
      const double logmmglmmu = std::log(mmgl/mmu);
      const double logmmst1mmu = std::log(mmst1/mmu);

      const double result =
" <> WrapLines @ IndentText @ IndentText[ToCPP[t1lLimitS2MG] <> ";"] <> "
      return result * g32 * oneLoop;

   }

   const double logmmglmmu  = std::log(mmgl/mmu);
   const double logmmst1mmu = std::log(mmst1/mmu);
   const double logmmst2mmu = std::log(mmst2/mmu);

   const double result =
" <> WrapLines @ IndentText[ToCPP[t1l - t1lqcd] <> ";"] <> "

   return result * g32 * oneLoop;
}

/// 1-loop full SQCD contributions to Delta Mt over mt [hep-ph/0210258]
double dMt_over_mt_1loop(const Parameters& pars)
{
   return dMt_over_mt_1loop_qcd(pars) + dMt_over_mt_1loop_susy(pars);
}

/// 2-loop QCD contributions to Delta Mt over mt [hep-ph/0507139]
double dMt_over_mt_2loop_qcd(const Parameters& pars)
{
   const double g34 = pow4(pars.g3);
   const double mmt = pow2(pars.mt);
   const double mmu = pow2(pars.Q);
   const double logmmtmmu = std::log(mmt/mmu);

   const double result =
" <> WrapLines @ IndentText[ToCPP[t2lqcd] <> ";"] <> "

   return result * g34 * twoLoop;
}

/// 2-loop SUSY contributions to Delta Mt over mt [hep-ph/0507139]
double dMt_over_mt_2loop_susy(const Parameters& pars)
{
   const double g34    = pow4(pars.g3);
   const double Xt     = pars.xt;
   const double mt     = pars.mt;
   const double mgl    = shift_mg(pars.mg, pars.mst1, pars.mst2);
   const double mmu    = pow2(pars.Q);
   const double mmt    = pow2(pars.mt);
   const double mmgl   = pow2(mgl);
   const double mmst1  = pow2(pars.mst1);
   const double mmst2  = pow2(pars.mst2);
   const double mmsusy = pow2(pars.msusy);
   const double SX     = 2*mt*pars.xt;
   const double s2t    = SX / (mmst1 - mmst2);

   if (is_equal(mmst1, mmst2, mmt) && is_equal(mmst1, mmgl, mmt) &&
       is_equal(mmst1, mmsusy, mmt) && is_equal(std::abs(Xt), 0., 1e-1)) {
      const double logmmsusymmu = std::log(mmsusy/mmu);
      const double logmmtmmu    = std::log(mmt/mmu);

      const double result =
" <> WrapLines @ IndentText @ IndentText[ToCPP[t2lLimitS1S2MSMG] <> ";"] <> "

      return result * g34 * twoLoop;
   }

   const double logmmsusymmu = std::log(mmsusy/mmu);
   const double logmmtmmu    = std::log(mmt/mmu);
   const double logmmst1mmu  = std::log(mmst1/mmu);
   const double logmmst2mmu  = std::log(mmst2/mmu);
   const double logmmglmmu   = std::log(mmgl/mmu );

   const double result =
" <> WrapLines @ IndentText[ToCPP[t2l - t2lqcd] <> ";"] <> "

   return result * g34 * twoLoop;
}

/// 2-loop full SQCD contributions to Delta Mt over mt [hep-ph/0507139]
double dMt_over_mt_2loop(const Parameters& pars)
{
   return dMt_over_mt_2loop_qcd(pars) + dMt_over_mt_2loop_susy(pars);
}

} // namespace mssm_twoloop_mt
} // namespace flexiblesusy
";

Export[headerName, header, "String"];
Export[implName  , impl  , "String"];
