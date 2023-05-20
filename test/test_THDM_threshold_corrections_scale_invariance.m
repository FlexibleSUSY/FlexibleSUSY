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
Needs["THDMThresholds1L`", FileNameJoin[{Directory[], "meta", "THDM", "Thresholds_1L_full.m"}]];

mssmGaugeRGEs = FileNameJoin[{Directory[], "Output", "MSSM", "RGEs", "BetaGauge.m"}];
thdmLambdaRGEs = FileNameJoin[{Directory[], "Output", "THDM-II", "RGEs", "BetaLijkl.m"}];
thdmGaugeRGEs = FileNameJoin[{Directory[], "Output", "THDM-II", "RGEs", "BetaGauge.m"}];

conj[0] = 0;
gRules = {g1 -> Sqrt[5/3] gY};

(* Note: Yukawa couplings are transposed, compared to SARAH *)
YuMat = Table[Yu[j, i], {i, 1, 3}, {j, 1, 3}];
YdMat = Table[Yd[j, i], {i, 1, 3}, {j, 1, 3}];
YeMat = Table[Ye[j, i], {i, 1, 3}, {j, 1, 3}];

expandTraces = {
    trace[Yu, Adj[Yu]] -> Tr[YuMat.ConjugateTranspose[YuMat]],
    trace[Yd, Adj[Yd]] -> Tr[YdMat.ConjugateTranspose[YdMat]],
    trace[Ye, Adj[Ye]] -> Tr[YeMat.ConjugateTranspose[YeMat]],
    trace[Yu, Adj[Yu], Yu, Adj[Yu]] -> Tr[YuMat.ConjugateTranspose[YuMat].YuMat.ConjugateTranspose[YuMat]],
    trace[Yd, Adj[Yd], Yd, Adj[Yd]] -> Tr[YdMat.ConjugateTranspose[YdMat].YdMat.ConjugateTranspose[YdMat]],
    trace[Ye, Adj[Ye], Ye, Adj[Ye]] -> Tr[YeMat.ConjugateTranspose[YeMat].YeMat.ConjugateTranspose[YeMat]],
    trace[Yd, Adj[Yu], Yu, Adj[Yd]] -> Tr[YdMat.ConjugateTranspose[YuMat].YuMat.ConjugateTranspose[YdMat]]
};

betag1MSSM = Cases[Get[mssmGaugeRGEs], {g1, b_, __} :> b][[1]];
betagYMSSM = (Sqrt[3/5] betag1MSSM /. gRules);
betag2MSSM = Cases[Get[mssmGaugeRGEs], {g2, b_, __} :> b][[1]] /. gRules;

betag1THDM = Cases[Get[thdmGaugeRGEs], {g1, b_, __} :> b][[1]];
betagYTHDM = (Sqrt[3/5] betag1THDM /. gRules);
betag2THDM = Cases[Get[thdmGaugeRGEs], {g2, b_, __} :> b][[1]] /. gRules;

(* tree-level in SARAH convention *)
lambdaTree = {
    1/8 (gY^2 + g2^2),
    1/8 (gY^2 + g2^2),
    1/4 (-gY^2 + g2^2),
    -1/2 g2^2,
    0, 0, 0
};

lambdaTreeRules = {
    Lambda1 -> lambdaTree[[1]],
    Lambda2 -> lambdaTree[[2]],
    Lambda3 -> lambdaTree[[3]],
    Lambda4 -> lambdaTree[[4]],
    Lambda5 -> lambdaTree[[5]],
    Lambda6 -> lambdaTree[[6]],
    Lambda7 -> lambdaTree[[7]]
};

(* beta functions of lambda_i in the MSSM *)
betaLambdaMSSM = (Dt[#] & /@ lambdaTree) /. {
    Dt[gY] -> betagYMSSM,
    Dt[g2] -> betag2MSSM
} // Simplify;

(* change in scale dependence since threshold corrections are
   expressed in terms of THDM gauge couplings *)
betaLambdaGaugeDiff = (Dt[#] & /@ lambdaTree) /. {
    Dt[gY] -> betagYMSSM - betagYTHDM,
    Dt[g2] -> betag2MSSM - betag2THDM
} // Simplify;

(* beta functions of lambda_i in the THDM *)
betaLambdaTHDM = {
    Cases[Get[thdmLambdaRGEs], {Lambda1, b_, __} :> b][[1]],
    Cases[Get[thdmLambdaRGEs], {Lambda2, b_, __} :> b][[1]],
    Cases[Get[thdmLambdaRGEs], {Lambda3, b_, __} :> b][[1]],
    Cases[Get[thdmLambdaRGEs], {Lambda4, b_, __} :> b][[1]],
    Cases[Get[thdmLambdaRGEs], {Lambda5, b_, __} :> b][[1]],
    Cases[Get[thdmLambdaRGEs], {Lambda6, b_, __} :> b][[1]],
    Cases[Get[thdmLambdaRGEs], {Lambda7, b_, __} :> b][[1]]
} /. gRules;

(* difference of the beta functions in the two models *)
betaDiff =
  Expand[(betaLambdaTHDM - betaLambdaMSSM + betaLambdaGaugeDiff) /. lambdaTreeRules //. expandTraces];

(* convert to SARAH convention *)
lamSARAH = GetTHDMThresholds1L[];
lamSARAH[[1]] = lamSARAH[[1]]/2;
lamSARAH[[2]] = lamSARAH[[2]]/2;

thresh = lamSARAH //. GetTHDMThresholds1LLoopFunctions[];

(* mu-dependence of threshold corrections *)
threshMuDep = Expand[Q D[thresh, Q] 16 Pi^2] /. Derivative[1][Re][_] -> 1;

TestEquality[Simplify[betaDiff - threshMuDep], Table[0, {i,1,7}]];

PrintTestSummary[];
