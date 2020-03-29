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

colorCA = 3; colorCF = 4/3; Tf = 1/2; GS = g3;
MGl = mgl; MT = mt; SX = 2 mt Xt; s2t = SX / (mmst1 - mmst2);
fin[0, args__] := fin[args, mmu];

Simp[expr_] := Collect[expr, { g3 }] //.
    {
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
   const double Pi  = 3.1415926535897932384626433832795;
   const double zt2 = 1.6449340668482264364724151666460;
   const double zt3 = 1.2020569031595942853997381615114;
   const double log2 = std::log(2.);

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
    * fin[] function from arXiv:hep-ph/0507139 .
    *
    * @param mm1 squared mass \\f$m_1^2\\f$
    * @param mm2 squared mass \\f$m_2^2\\f$
    * @param mmu squared renormalization scale
    *
    * @return fin(m12, m22)
    */
   double fin(double mm1, double mm2, double mmu)
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

   /// shift gluino mass away from mst1 and mst2 if too close
   double shift_mg(double mg, double mst1, double mst2)
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
   using std::log;
   const double g3  = pars.g3;
   const double mmt = pow2(pars.mt);
   const double mmu = pow2(pars.Q);

   const double result = " <> ToCPP[t1lqcd] <> ";

   return result * oneLoop;
}

/// 1-loop SUSY contributions to Delta Mt over mt [hep-ph/0210258]
double dMt_over_mt_1loop_susy(const Parameters& pars)
{
   using std::log;
   const double g3     = pars.g3;
   const double Xt     = pars.xt;
   const double mt     = pars.mt;
   const double mgl    = pars.mg;
   const double mmu    = pow2(pars.Q);
   const double mmgl   = pow2(pars.mg);
   const double mmst1  = pow2(pars.mst1);
   const double mmst2  = pow2(pars.mst2);

   if (is_equal(mmst1, mmst2, 1e-6) && is_equal(mmst1, mmgl, 1e-6)) {
      const double result = " <> ToCPP[t1lLimitS1S2MG] <> ";
      return result * oneLoop;
   }

   if (is_equal(mmst1, mmst2, 1e-6)) {
      const double result =
" <> WrapLines @ IndentText @ IndentText[ToCPP[t1lLimitS1S2] <> ";"] <> "

      return result * oneLoop;
   }

   if (is_equal(mmgl, mmst1, 1e-6)) {
      const double result =
" <> WrapLines @ IndentText @ IndentText[ToCPP[t1lLimitS1MG] <> ";"] <> "

      return result * oneLoop;
   }

   if (is_equal(mmgl, mmst2, 1e-6)) {
      const double result =
" <> WrapLines @ IndentText @ IndentText[ToCPP[t1lLimitS2MG] <> ";"] <> "
      return result * oneLoop;

   }

   const double result =
" <> WrapLines @ IndentText[ToCPP[t1l - t1lqcd] <> ";"] <> "

   return result * oneLoop;
}

/// 1-loop full SQCD contributions to Delta Mt over mt [hep-ph/0210258]
double dMt_over_mt_1loop(const Parameters& pars)
{
   return dMt_over_mt_1loop_qcd(pars) + dMt_over_mt_1loop_susy(pars);
}

/// 2-loop QCD contributions to Delta Mt over mt [hep-ph/0507139]
double dMt_over_mt_2loop_qcd(const Parameters& pars)
{
   using std::log;
   const double g3  = pars.g3;
   const double mmt = pow2(pars.mt);
   const double mmu = pow2(pars.Q);

   const double result =
" <> WrapLines @ IndentText[ToCPP[t2lqcd] <> ";"] <> "

   return result * twoLoop;
}

/// 2-loop SUSY contributions to Delta Mt over mt [hep-ph/0507139]
double dMt_over_mt_2loop_susy(const Parameters& pars)
{
   using std::log;
   const double g3     = pars.g3;
   const double Xt     = pars.xt;
   const double mt     = pars.mt;
   const double mgl    = shift_mg(pars.mg, pars.mst1, pars.mst2);
   const double mmu    = pow2(pars.Q);
   const double mmt    = pow2(pars.mt);
   const double mmgl   = pow2(mgl);
   const double mmst1  = pow2(pars.mst1);
   const double mmst2  = pow2(pars.mst2);
   const double mmsusy = pow2(pars.msusy);

   if (is_equal(mmst1, mmst2, mmt) && is_equal(mmst1, mmgl, mmt) &&
       is_equal(mmst1, mmsusy, mmt) && is_equal(std::abs(Xt), 0., 1e-1)) {
      const double result =
" <> WrapLines @ IndentText @ IndentText[ToCPP[t2lLimitS1S2MSMG] <> ";"] <> "

      return result * twoLoop;
   }

   const double result =
" <> WrapLines @ IndentText[ToCPP[t2l - t2lqcd] <> ";"] <> "

   return result * twoLoop;
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
