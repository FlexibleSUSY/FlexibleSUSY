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

(* TODO:
 * Implement limit mU3^2 = mU3^2 = m3^2 or numerically deviate.
 *)

Get[FileNameJoin[{"meta", "TextFormatting.m"}]];
Get[FileNameJoin[{"meta", "MSSM", "Mh2_3loop_DR_SQCD_program.m"}]];

M  = MR;
z2 = zt2;
z3 = zt3;
mt  /: mt^n_  := mt2^(n/2)  /; EvenQ[n];
m3  /: m3^n_  := m32^(n/2)  /; EvenQ[n];
mU3 /: mU3^n_ := mU32^(n/2) /; EvenQ[n];
mQ3 /: mQ3^n_ := mQ32^(n/2) /; EvenQ[n];
msq /: msq^n_ := msq2^(n/2) /; EvenQ[n];
MR  /: MR^n_  := MR2^(n/2)  /; EvenQ[n];

Simp[expr_] := expr //.
    {
        mt^2                  -> mt2,
        mt^-2                 -> 1/mt2,
        m3^2                  -> m32,
        m3^-2                 -> 1/m32,
        mU3^2                 -> mU32,
        mU3^-2                -> 1/mU32,
        mQ3^2                 -> mQ32,
        mQ3^-2                -> 1/mQ32,
        msq^2                 -> msq2,
        msq^-2                -> 1/msq2,
        MR^2                  -> MR2,
        MR^-2                 -> 1/MR2,
        Power[x_,n_] /; n > 0 :> Symbol["pow" <> ToString[n]][x],
        Power[x_,-2]          :> 1/Symbol["pow" <> ToString[2]][x],
        Power[x_,-3]          :> 1/Symbol["pow" <> ToString[3]][x],
        Power[x_,-4]          :> 1/Symbol["pow" <> ToString[4]][x],
        Power[x_,-5]          :> 1/Symbol["pow" <> ToString[5]][x],
        Power[x_,-6]          :> 1/Symbol["pow" <> ToString[6]][x],
        Log[x_]               :> log[x],
        PolyLog[2,x_]         :> dilog[x],
        PolyLog[4,1/2]        -> N[PolyLog[4,1/2]]
    };

ToCPP[expr_] := ToString[Simp[expr], CForm];

Corrections[alphaSOrder_, logOrder_] :=
    k^((alphaSOrder-1)/2) deltamh2[alphaSOrder, logOrder] /.
    SM3LoopConst /.
    Inactive -> Identity;

impl = "\
#include \"dilog.hpp\"
#include <cmath>
#include <limits>

namespace himalaya {
namespace mh2_eft {

namespace {
   const double Pi  = 3.1415926535897932384626433832795;
   const double zt2 = 1.6449340668482264364724151666460;
   const double zt3 = 1.2020569031595942853997381615114;

   template <typename T> T pow2(T x)  { return x*x; }
   template <typename T> T pow3(T x)  { return x*x*x; }
   template <typename T> T pow4(T x)  { return x*x*x*x; }
   template <typename T> T pow5(T x)  { return x*x*x*x*x; }
   template <typename T> T pow6(T x)  { return x*x*x*x*x*x; }
   template <typename T> T pow7(T x)  { return x*x*x*x*x*x*x; }
   template <typename T> T pow8(T x)  { return x*x*x*x*x*x*x*x; }
   template <typename T> T pow9(T x)  { return x*x*x*x*x*x*x*x*x; }
   template <typename T> T pow11(T x) { return pow2(x)*pow9(x); }
   template <typename T> T pow12(T x) { return pow2(pow6(x)); }
   template <typename T> T pow13(T x) { return pow4(x)*pow9(x); }
   template <typename T> T pow14(T x) { return pow2(pow7(x)); }
   template <typename T> T pow15(T x) { return pow6(x)*pow9(x); }
   template <typename T> T pow16(T x) { return pow2(pow8(x)); }
   template <typename T> T pow18(T x) { return pow2(pow9(x)); }
} // anonymous namespace

/// 1-loop coefficient O(at*log^0)
double coeff_as_0_log_0(double mQ32, double mU32, double Xt, double MR2)
{
   using std::log;

   const double result =
" <> WrapLines @ IndentText @ IndentText[ToCPP[Corrections[0,0]] <> ";"] <> "
   return result;
}

/// 1-loop coefficient O(at*log^1)
double coeff_as_0_log_1()
{
   return " <> ToCPP[Corrections[0,1]] <> ";
}

/// 2-loop coefficient O(at*as*log^0)
double coeff_as_1_log_0(double mQ32, double mU32, double Xt, double m3, double MR2)
{
   using std::log;
   const double m32 = pow2(m3);

   const double result =
" <> WrapLines @ IndentText @ IndentText[ToCPP[Corrections[1,0]] <> ";"] <> "
   return result;
}

/// 2-loop coefficient O(at*as*log^1)
double coeff_as_1_log_1(double mQ32, double mU32, double Xt, double m3, double MR2)
{
   using std::log;
   const double m32 = pow2(m3);

   const double result =
" <> WrapLines @ IndentText @ IndentText[ToCPP[Corrections[1,1]] <> ";"] <> "
   return result;
}

/// 2-loop coefficient O(at*as*log^2)
double coeff_as_1_log_2()
{
   return " <> ToCPP[Corrections[1,2]] <> ";
}

/// 3-loop coefficient O(at*as^2*log^0)
double coeff_as_2_log_0(double mQ32, double mU32, double Xt, double m3, double msq2, double MR2)
{
   using std::log;
   const double m32 = pow2(m3);
   const double dlatas2 = 0.; // 3L threshold correction Δλ^(3L) O(at*as^2)

   const double result =
" <> WrapLines @ IndentText @ IndentText[ToCPP[Corrections[2,0]] <> ";"] <> "
   return result;
}

/// 3-loop coefficient O(at*as^2*log^1)
double coeff_as_2_log_1(double mQ32, double mU32, double Xt, double m3, double msq2, double MR2)
{
   using std::log;
   const double m32 = pow2(m3);

   const double result =
" <> WrapLines @ IndentText @ IndentText[ToCPP[Corrections[2,1]] <> ";"] <> "
   return result;
}

/// 3-loop coefficient O(at*as^2*log^2)
double coeff_as_2_log_2(double mQ32, double mU32, double Xt, double m3, double msq2, double MR2)
{
   using std::log;
   const double m32 = pow2(m3);

   const double result =
" <> WrapLines @ IndentText @ IndentText[ToCPP[Corrections[2,2]] <> ";"] <> "
   return result;
}

/// 3-loop coefficient O(at*as^2*log^3)
double coeff_as_2_log_3()
{
   return " <> ToCPP[Corrections[2,3]] <> ";
}

/// 1-loop correction O(at), at = yt^2 Sin[beta]^2/(4 Pi)
double Mh2_EFT_1loop(
   double at, double mt, double mQ32, double mU32,
   double Xt, double MR2)
{
   const double mt2 = pow2(mt);
   const double L = std::log(MR2/mt2);
   const double oneLoop = 1./(4*Pi);

   return oneLoop * mt2 * at * (
      coeff_as_0_log_0(mQ32, mU32, Xt, MR2) +
      coeff_as_0_log_1() * L
   );
}

/// 2-loop correction O(at*as)
/// at = yt^2 Sin[beta]^2/(4 Pi), as = g3^2/(4 Pi), 
double Mh2_EFT_2loop(
   double at, double mt, double mQ32, double mU32,
   double Xt, double MR2, double g3, double m3)
{
   const double mt2 = pow2(mt);
   const double L = std::log(MR2/mt2);
   const double as = pow2(g3)/(4*Pi);
   const double twoLoop = 1./pow2(4*Pi);

   return twoLoop * mt2 * at * as * (
      coeff_as_1_log_0(mQ32, mU32, Xt, m3, MR2) +
      coeff_as_1_log_1(mQ32, mU32, Xt, m3, MR2) * L +
      coeff_as_1_log_2() * pow2(L)
   );
}

/// 3-loop correction O(at*as^2), without Δλ^(3L),
/// at = yt^2 Sin[beta]^2/(4 Pi), as = g3^2/(4 Pi), 
double Mh2_EFT_3loop(
   double at, double mt, double mQ32, double mU32,
   double Xt, double MR2, double g3, double m3, double msq2)
{
   const double mt2 = pow2(mt);
   const double L = std::log(MR2/mt2);
   const double as = pow2(g3)/(4*Pi);
   const double threeLoop = 1./pow3(4*Pi);

   return threeLoop * mt2 * at * pow2(as) * (
      coeff_as_2_log_0(mQ32, mU32, Xt, m3, msq2, MR2) +
      coeff_as_2_log_1(mQ32, mU32, Xt, m3, msq2, MR2) * L +
      coeff_as_2_log_2(mQ32, mU32, Xt, m3, msq2, MR2) * pow2(L) +
      coeff_as_2_log_3() * pow3(L)
   );
}

} // namespace mh2_eft
} // namespace himalaya
";

Export["Mh2_EFT.cpp",impl, "String"];
Print["FINISHED"];
