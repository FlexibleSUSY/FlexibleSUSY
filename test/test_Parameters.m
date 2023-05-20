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
Needs["Parameters`", "Parameters.m"];

Off[General::stop];

Print["testing GuessInputParameterType[] ..."];

inputPars = {m, n, Sign[o], FlexibleSUSY`Phase[p]};

TestEquality[Parameters`GuessInputParameterType[a],
             CConversion`ScalarType[CConversion`realScalarCType]];
TestEquality[Parameters`GuessInputParameterType[m],
             CConversion`ScalarType[CConversion`realScalarCType]];
TestEquality[Parameters`GuessInputParameterType[Sign[o]],
             CConversion`ScalarType[CConversion`integerScalarCType]];
TestEquality[Parameters`GuessInputParameterType[FlexibleSUSY`Phase[p]],
             CConversion`ScalarType[CConversion`complexScalarCType]];

Print["testing GuessExtraParameterType[] ..."];

extraPars = {t, u, FlexibleSUSY`Phase[v]};

TestEquality[Parameters`GuessExtraParameterType[t],
             CConversion`ScalarType[CConversion`realScalarCType]];
TestEquality[Parameters`GuessExtraParameterType[FlexibleSUSY`Phase[v]],
             CConversion`ScalarType[CConversion`complexScalarCType]];

Print["testing IsRealParameter[] ..."];

realPars = {a, b, c};
compPars = {x, y, z};
allPars  = Join[realPars, compPars];

Parameters`SetModelParameters[allPars];
Parameters`AddRealParameter[realPars];
SARAH`RealParameters = {};
Parameters`SetInputParameters[inputPars];
Parameters`AddInputParameters[{{q, CConversion`MatrixType[CConversion`realScalarCType, 3, 3]}}];
Parameters`AddRealParameter[{t}];
Parameters`SetExtraParameters[extraPars];
Parameters`AddExtraParameters[{{w, CConversion`TensorType[CConversion`realScalarCType, 2, 2, 3]}}];

TestEquality[Parameters`IsRealParameter[a], True];
TestEquality[Parameters`IsRealParameter[b], True];
TestEquality[Parameters`IsRealParameter[c], True];
TestEquality[Parameters`IsRealParameter[x], False];
TestEquality[Parameters`IsRealParameter[y], False];
TestEquality[Parameters`IsRealParameter[z], False];
TestEquality[Parameters`IsRealParameter[m], True];
TestEquality[Parameters`IsRealParameter[n], True];
TestEquality[Parameters`IsRealParameter[Sign[o]], True];
TestEquality[Parameters`IsRealParameter[FlexibleSUSY`Phase[p]], False];
TestEquality[Parameters`IsRealParameter[q], True];
TestEquality[Parameters`IsRealParameter[t], True];
TestEquality[Parameters`IsRealParameter[u], True];
TestEquality[Parameters`IsRealParameter[FlexibleSUSY`Phase[v]], False];
TestEquality[Parameters`IsRealParameter[w], True];

Print["testing IsComplexParameter[] ..."];

TestEquality[Parameters`IsComplexParameter[a], False];
TestEquality[Parameters`IsComplexParameter[x], True];
TestEquality[Parameters`IsComplexParameter[m], False];
TestEquality[Parameters`IsComplexParameter[FlexibleSUSY`Phase[p]], True];

Print["testing IsRealExpression[] ..."];

TestEquality[Parameters`IsRealExpression[a], True];
TestEquality[Parameters`IsRealExpression[a[1]], True];
TestEquality[Parameters`IsRealExpression[x], False];

TestEquality[Parameters`IsRealExpression[a^2], True];
TestEquality[Parameters`IsRealExpression[a b], True];
TestEquality[Parameters`IsRealExpression[a b c], True];

TestEquality[Parameters`IsRealExpression[Sin[a] b], True];

TestEquality[Parameters`IsRealExpression[a b c x], False];
TestEquality[Parameters`IsRealExpression[x^2], False];
TestEquality[Parameters`IsRealExpression[Conjugate[x] x], True];
TestEquality[Parameters`IsRealExpression[conj[x] x], True];
TestEquality[Parameters`IsRealExpression[Conj[x] x], True];

TestEquality[Parameters`IsRealExpression[trace[a]], True];
TestEquality[Parameters`IsRealExpression[trace[x]], False];

TestEquality[Parameters`IsRealExpression[trace[a b]], True];

TestEquality[Parameters`IsRealExpression[trace[Adj[a],a]], True];
TestEquality[Parameters`IsRealExpression[trace[Adj[x],x]], True];

TestEquality[Parameters`IsRealExpression[trace[Adj[x],x,Adj[y],y]], True];
TestEquality[Parameters`IsRealExpression[trace[Adj[x],x,Adj[y],y,Adj[z],z]], False];

TestEquality[Parameters`IsRealExpression[trace[conj[a],Tp[b]]], True];
TestEquality[Parameters`IsRealExpression[trace[a,Tp[b]]], True];
TestEquality[Parameters`IsRealExpression[trace[a,b]], True];

TestEquality[Parameters`IsRealExpression[trace[conj[x],x]], False];
TestEquality[Parameters`IsRealExpression[trace[conj[x],Tp[x]]], True];
TestEquality[Parameters`IsRealExpression[trace[x,Tp[y]]], False];
TestEquality[Parameters`IsRealExpression[trace[x,y]], False];

TestEquality[Parameters`IsRealExpression[trace[conj[x],Tp[x],Adj[y],y]], True];

(* make a hermitian*)
SARAH`getDimParameters[a] = {1};
TestEquality[Parameters`Private`IsHermitian[a], True];

TestEquality[Parameters`IsRealExpression[trace[Adj[x],x,a]], True];

(* make x hermitian*)
SARAH`getDimParameters[x] = {2,2};
SARAH`ListSoftBreakingScalarMasses = {x};

TestEquality[Parameters`Private`IsHermitian[x], True];
TestEquality[Parameters`IsRealExpression[trace[x]], True];

(* test sum and products *)
TestEquality[Parameters`IsRealExpression[c + b trace[x]], True];
TestEquality[Parameters`IsRealExpression[c - b trace[x] - 2 trace[a]], True];

Print["testing IsMatrix[] ..."];

TestEquality[Parameters`IsMatrix[x], True];
TestEquality[Parameters`IsMatrix[m], False];
TestEquality[Parameters`IsMatrix[q], True];
TestEquality[Parameters`IsMatrix[w], False];

Print["testing IsTensor[] ..."];

TestEquality[Parameters`IsTensor[w], True];

Print["testing FindAllParameters[] ..."];

modelParameters = { Mu, SARAH`B[Mu], WOp, SARAH`Q[WOp] };
Parameters`SetModelParameters[modelParameters];

expr = 2 * Mu SARAH`B[Mu] + WOp SARAH`Q[WOp] + a;

TestEquality[Sort[Parameters`FindAllParameters[expr]],
             Sort[modelParameters]
            ];

expr = 2 * Mu Re[SARAH`B[Mu]] + Im[WOp] SARAH`Q[WOp] + a;

TestEquality[Sort[Parameters`FindAllParameters[expr]],
             Sort[modelParameters]
            ];

expr = 2 * SARAH`B[Mu] + SARAH`Q[WOp];

TestEquality[Sort[Parameters`FindAllParameters[expr]],
             Sort[{SARAH`B[Mu], SARAH`Q[WOp]}]
            ];

expr = Which[WOp > 1, Sqrt[B[Mu]]];

TestEquality[Sort[Parameters`FindAllParameters[expr]],
             Sort[{B[Mu], WOp}]];

Print["testing GetType[] ..."];

TestEquality[Parameters`GetType[x],
             CConversion`MatrixType[CConversion`complexScalarCType, 2, 2]];
TestEquality[Parameters`GetType[Sign[o]],
             CConversion`ScalarType[CConversion`integerScalarCType]];
TestEquality[Parameters`GetType[FlexibleSUSY`Phase[p]],
             CConversion`ScalarType[CConversion`complexScalarCType]];
TestEquality[Parameters`GetType[m],
             CConversion`ScalarType[CConversion`realScalarCType]];
TestEquality[Parameters`GetType[q],
             CConversion`MatrixType[CConversion`realScalarCType, 3, 3]];
TestEquality[Parameters`GetType[t],
             CConversion`ScalarType[CConversion`realScalarCType]];
TestEquality[Parameters`GetType[FlexibleSUSY`Phase[v]],
             CConversion`ScalarType[CConversion`complexScalarCType]];
TestEquality[Parameters`GetType[w],
             CConversion`TensorType[CConversion`realScalarCType, 2, 2, 3]];

Print["testing GetParameterDimensions[] ..."];

TestEquality[Parameters`GetParameterDimensions[x], {2, 2}];
TestEquality[Parameters`GetParameterDimensions[m], {1}];
TestEquality[Parameters`GetParameterDimensions[Sign[o]], {1}];
TestEquality[Parameters`GetParameterDimensions[q], {3, 3}];
TestEquality[Parameters`GetParameterDimensions[t], {1}];
TestEquality[Parameters`GetParameterDimensions[w], {2, 2, 3}];

Print["testing GetRealTypeFromDimension[] ..."];

TestEquality[Parameters`GetRealTypeFromDimension[{1}],
             CConversion`ScalarType[CConversion`realScalarCType]];
TestEquality[Parameters`GetRealTypeFromDimension[{2}],
             CConversion`VectorType[CConversion`realScalarCType, 2]];
TestEquality[Parameters`GetRealTypeFromDimension[{0, 0}],
             CConversion`MatrixType[CConversion`realScalarCType, 0, 0]];
TestEquality[Parameters`GetRealTypeFromDimension[{1, 1}],
             CConversion`MatrixType[CConversion`realScalarCType, 1, 1]];
TestEquality[Parameters`GetRealTypeFromDimension[{1, 3}],
             CConversion`MatrixType[CConversion`realScalarCType, 1, 3]];
TestEquality[Parameters`GetRealTypeFromDimension[{3, 1}],
             CConversion`MatrixType[CConversion`realScalarCType, 3, 1]];
TestEquality[Parameters`GetRealTypeFromDimension[{2, 2}],
             CConversion`MatrixType[CConversion`realScalarCType, 2, 2]];
TestEquality[Parameters`GetRealTypeFromDimension[{1, 1, 1}],
             CConversion`TensorType[CConversion`realScalarCType, 1, 1, 1]];
TestEquality[Parameters`GetRealTypeFromDimension[{2, 3, 2}],
             CConversion`TensorType[CConversion`realScalarCType, 2, 3, 2]];
TestEquality[Parameters`GetRealTypeFromDimension[{5, 4, 3, 2}],
             CConversion`TensorType[CConversion`realScalarCType, 5, 4, 3, 2]];

PrintTestSummary[];
