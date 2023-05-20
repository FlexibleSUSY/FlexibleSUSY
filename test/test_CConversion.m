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

Needs["TestSuite`", "TestSuite.m"];
Needs["CConversion`", "CConversion.m"];

Print["testing ConvertGreekLetters[] ..."];

TestEquality[CConversion`Private`ConvertGreekLetters[\[Alpha]], Alpha];
TestEquality[CConversion`Private`ConvertGreekLetters[\[Beta]], Betax]; (* Beta is already defined by Mathematica *)
TestEquality[CConversion`Private`ConvertGreekLetters[\[Mu]]   , Mu];
TestEquality[CConversion`Private`ConvertGreekLetters[\[Zeta]], Zetax]; (* Zeta is already defined by Mathematica *)
SARAH`Delta;
TestEquality[CConversion`Private`ConvertGreekLetters[\[Delta]], Deltax];

Print["testing ToValidCSymbol[] ..."];

TestEquality[ToValidCSymbol[1], 1];
TestEquality[ToValidCSymbol[1.1], 1.1];
TestEquality[ToValidCSymbol[a], a];
TestEquality[ToValidCSymbol[a[1]], a1];
TestEquality[ToValidCSymbol[a[b]], ab];
TestEquality[ToValidCSymbol[a[b,c]], abc];
TestEquality[ToValidCSymbol[a[b,c][d]], abcd];

(* Test that symbols, which consist of greek letters are correctly
   converted. *)
TestEquality[ToValidCSymbol[\[Lambda]3], Lambda3];
TestEquality[ToValidCSymbol[T[\[Lambda]3]], TLambda3];
TestEquality[ToValidCSymbol[\[Lambda]], Lambdax];
(* no unique name necessary *)
(* TestEquality[ToValidCSymbol[T[\[Lambda]]], TLambda]; *)

Print["testing ToValidCSymbolString[] ..."];

TestEquality[ToValidCSymbolString[1], "1"];
TestEquality[ToValidCSymbolString[1.1], "1.1"];
TestEquality[ToValidCSymbolString[a], "a"];
TestEquality[ToValidCSymbolString[a[b]], "ab"];
TestEquality[ToValidCSymbolString[a[b,c]], "abc"];
TestEquality[ToValidCSymbolString[a[b,c][d]], "abcd"];

Print["testing RValueToCFormString[] ..."];

TestEquality[RValueToCFormString[1], "1"];
TestEquality[RValueToCFormString[a], "a"];
TestEquality[RValueToCFormString[a+b], "a + b"];
TestEquality[RValueToCFormString[a^1], "a"];
TestEquality[RValueToCFormString[a^2], "Sqr(a)"];
TestEquality[RValueToCFormString[f[x]], "f(x)"];
TestEquality[RValueToCFormString[a MatMul[A]], "a*A"];
TestEquality[RValueToCFormString[a MatMul[A,B]], "a*(A*B)"];
TestEquality[RValueToCFormString[a MatMul[A,B,A]], "a*(A*B*A)"];
TestEquality[RValueToCFormString[a trace[A]], "a*(A).trace()"];
TestEquality[RValueToCFormString[a trace[A,B]], "a*(A*B).trace()"];
TestEquality[RValueToCFormString[a trace[A,B,A]], "a*(A*B*A).trace()"];

TestEquality[RValueToCFormString[MatMul[Adj[A]]], "A.adjoint()"];
TestEquality[RValueToCFormString[MatMul[A Adj[A]]], "A*A.adjoint()"];
TestEquality[RValueToCFormString[trace[Adj[A]]], "(A.adjoint()).trace()"];
TestEquality[RValueToCFormString[trace[A Adj[A]]], "(A*A.adjoint()).trace()"];

(* test greek symbol conversion *)
TestEquality[RValueToCFormString[\[Mu]], "Mu"];
TestEquality[RValueToCFormString[SARAH`B[\[Mu]]], "BMu"];
TestEquality[RValueToCFormString[A[\[Mu]]], "A(Mu)"];
TestEquality[RValueToCFormString[\[Mu][1]], "Mu(1)"];
TestEquality[RValueToCFormString[\[Mu][1,2]], "Mu(1,2)"];
TestEquality[RValueToCFormString[\[Mu][1,2,3]], "Mu(1,2,3)"];
TestEquality[RValueToCFormString[MACROSTRING[T[Yu]][0,0]/MACROSTRING[Yu][0,0]], "MACROSTRING(TYu)(0,0)/MACROSTRING(Yu)(0,0)"];
TestEquality[RValueToCFormString[MACROSTRING[T[\[Kappa]]][0,0]/MACROSTRING[\[Kappa]][0,0]], "MACROSTRING(TKappa)(0,0)/MACROSTRING(Kappa)(0,0)"];

TestEquality[CreateCBoolValue[True], "true"];
TestEquality[CreateCBoolValue[False], "false"];

Print["testing GetHead[] ..."];

TestEquality[GetHead[1], 1];
TestEquality[GetHead[a], a];
TestEquality[GetHead[f[x]], f];
TestEquality[GetHead[f[x][y]], f];
TestEquality[GetHead[f[x[y]]], f];

PrintTestSummary[];
