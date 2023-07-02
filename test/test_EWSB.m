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
Needs["EWSB`", "EWSB.m"];
Needs["Parameters`", "Parameters.m"];

Print["testing MSSM-like EWSB for Mu and BMu ..."];

FlexibleSUSY`FSSolveEWSBTimeConstraint = 120;

mssmEwsbEqs = {
    mu^2 + x^2 + x y + z + 5,
    Bmu  - x^2 + x y + z + 5
};

mssmEwsbOutputParameters = { mu, Bmu };

Parameters`SetRealParameters[mssmEwsbOutputParameters];

{mssmSolution, mssmFreePhases} = EWSB`FindSolutionAndFreePhases[mssmEwsbEqs, mssmEwsbOutputParameters];

TestEquality[mssmFreePhases, {Sign[mu]}];

mssmFullSolution = EWSB`Private`FindSolution[mssmEwsbEqs, mssmEwsbOutputParameters];

TestEquality[Sort /@ mssmFullSolution,
             Sort /@ { {mu -> -Sqrt[-5 - x^2 - x*y - z], B[mu] -> -5 + x^2 - x*y - z},
                       {mu ->  Sqrt[-5 - x^2 - x*y - z], B[mu] -> -5 + x^2 - x*y - z}
                     }];

TestEquality[Sort[mssmSolution],
             Sort[{ B[mu] -> -5 + x^2 - x*y - z,
                    mu -> Sign[mu] Sqrt[-5 - x^2 - x*y - z]
                  }]];

Print["testing NMSSM-like EWSB for Kappa, vS and mS2 ..."];

muEff = lambda s;

BmuEff = Alambda + kappa s;

nmssmEwsbEqs = {
    vu (mHu2 + muEff^2 + lambda^2 vd^2 + g^2 (vu^2 - vd^2)) - vd muEff BmuEff,
    vd (mHd2 + muEff^2 + lambda^2 vu^2 + g^2 (vd^2 - vu^2)) - vu muEff BmuEff,
    s (mS2 + X + kappa s Akappa + kappa^2 s^2) + Y
};

nmssmEwsbOutputParameters = { s, kappa, mS2 };

Parameters`SetRealParameters[nmssmEwsbOutputParameters];

{nmssmSolution, nmssmFreePhases} = EWSB`FindSolutionAndFreePhases[nmssmEwsbEqs, nmssmEwsbOutputParameters];

TestEquality[nmssmFreePhases, {Sign[s]}];

nmssmFullSolution = EWSB`Private`FindSolution[nmssmEwsbEqs, nmssmEwsbOutputParameters];

TestEquality[Sort /@ nmssmFullSolution,
             Sort /@ EWSB`Private`ToMathematicaSolutionFormat@
                     {{{s -> -(Sqrt[-(mHd2*vd^2) - g^2*vd^4 + mHu2*vu^2 + g^2*vu^4]/
                               Sqrt[lambda^2*vd^2 - lambda^2*vu^2])},
                       {s -> Sqrt[-(mHd2*vd^2) - g^2*vd^4 + mHu2*vu^2 + g^2*vu^4]/
                        Sqrt[lambda^2*vd^2 - lambda^2*vu^2]}},
                      {{kappa ->
                        (-(Alambda*lambda*s*vd) + mHu2*vu + lambda^2*s^2*vu - g^2*vd^2*vu +
                         lambda^2*vd^2*vu + g^2*vu^3)/(lambda*s^2*vd)}},
                      {{mS2 -> (-(Akappa*kappa*s^2) - kappa^2*s^3 - s*X - Y)/s}}
                     }];

TestEquality[Sort[nmssmSolution],
             Sort[{ kappa -> (-(Alambda*lambda*s*vd) + mHu2*vu + lambda^2*s^2*vu -
                              g^2*vd^2*vu + lambda^2*vd^2*vu + g^2*vu^3)/(lambda*s^2*vd),
                    mS2 -> (-(Akappa*kappa*s^2) - kappa^2*s^3 - s*X - Y)/s,
                    s -> (Sqrt[-(mHd2*vd^2) - g^2*vd^4 + mHu2*vu^2 + g^2*vu^4]*
                          Sign[s])/Sqrt[lambda^2*vd^2 - lambda^2*vu^2]
                  }]];

Print["testing UMSSM-like EWSB for mHu2, mHd2 and mS2 ..."];

umssmEwsbEqs = {
    vu mHu2 + g^2 (vu^2 - vd^2),
    vd mHd2 + g^2 (vd^2 - vu^2),
    s mS2 + X
};

umssmEwsbOutputParameters = { mHu2, mHd2, mS2 };

Parameters`SetRealParameters[umssmEwsbOutputParameters];

{umssmSolution, umssmFreePhases} = EWSB`FindSolutionAndFreePhases[umssmEwsbEqs, umssmEwsbOutputParameters];

TestEquality[umssmFreePhases, {}];

umssmFullSolution = EWSB`Private`FindSolution[umssmEwsbEqs, umssmEwsbOutputParameters];

TestEquality[Sort /@ umssmFullSolution,
             Sort /@ EWSB`Private`ToMathematicaSolutionFormat@
                     { {{mHu2 -> (g^2*(vd^2 - vu^2))/vu}},
                       {{mHd2 -> -(g^2*(vd^2 - vu^2))/vd}},
                       {{mS2 -> -(X/s)}}
                     }];

TestEquality[Sort /@ umssmSolution,
             Sort /@ { mHu2 -> (g^2*(vd^2 - vu^2))/vu,
                       mHd2 -> -(g^2*(vd^2 - vu^2))/vd,
                       mS2 -> -(X/s)
                     }];

Print["testing EWSB for vu, vd, s ..."];

(* This test ensures that the algorith works even in cases where only
   one independent sub-equation can be found. *)

oneIndependentSubeq = {
    vu mHu2 + X vd + s,
    vd mHd2 + Y vu + s,
    s mS2 + Z           (* this Eq. is independent of all the others*)
};

oneIndependentSubeqEwsbOutputParameters = { vu, vd, s };

Parameters`SetRealParameters[oneIndependentSubeqEwsbOutputParameters];

{oneIndependentSubeqSolution, oneIndependentSubeqFreePhases} =
    EWSB`FindSolutionAndFreePhases[oneIndependentSubeq, oneIndependentSubeqEwsbOutputParameters];

TestEquality[oneIndependentSubeqFreePhases, {}];

TestEquality[Sort /@ Simplify[oneIndependentSubeqSolution],
             Sort /@ { s -> -(Z/mS2),
                       vu -> -((s*(mHd2 - X))/(mHd2*mHu2 - X*Y)),
                       vd -> (-s - mHu2*vu)/X
                     }];

Print["testing NSM EWSB for mH2, mS2 ..."];

nsmEwsbOutputParameters = {mH2, mS2};

Parameters`SetRealParameters[nsmEwsbOutputParameters];

nsmEwsbEqs = {
    mH2*vH - vH^3*l1 - vH*vS^2*l3 - vH*vS*l4 - tadpole[1],
    2*mS2*vS - 4*vS^3*l2 - vH^2*vS*l3 - (vH^2*l4)/2 - 3*vS^2*l5 - tadpole[2]
};

nsmFullSolution = EWSB`Private`FindSolution[nsmEwsbEqs, nsmEwsbOutputParameters];

TestEquality[Sort /@ nsmFullSolution,
             Sort /@ EWSB`Private`ToMathematicaSolutionFormat@
                     { {{mH2 -> (vH^3*l1 + vH*vS^2*l3 + vH*vS*l4 + tadpole[1])/vH}},
                       {{mS2 -> (8*vS^3*l2 + 2*vH^2*vS*l3 + vH^2*l4 + 6*vS^2*l5 + 2*tadpole[2])/(4*vS)}}
                     }];

Print["testing cMSSM-like EWSB for |Mu| and BMu ..."];

cmssmEwsbEqs = {
    Susyno`LieGroups`conj[Mu] Mu + x^2 + x y + z + 5,
    B[Mu] - x^2 + x y + z + 5
};

cmssmEwsbOutputParameters = { Mu, B[Mu] };

TestEquality[Parameters`IsRealParameter[Mu], False];
TestEquality[Parameters`IsRealParameter[B[Mu]], False];

Print["\t calling FindSolution[] ..."];

cmssmFullSolution = EWSB`Private`FindSolution[cmssmEwsbEqs, cmssmEwsbOutputParameters];

TestEquality[Sort /@ cmssmFullSolution,
             Sort /@ EWSB`Private`ToMathematicaSolutionFormat@
                     { {{Mu -> -Sqrt[-5 - x^2 - x*y - z]},
                        {Mu -> Sqrt[-5 - x^2 - x*y - z]}},
                       {{B[Mu] -> -5 + x^2 - x*y - z}}
                     }];

Print["\t calling FindSolutionAndFreePhases[] ..."];

{cmssmSolution, cmssmFreePhases} = EWSB`FindSolutionAndFreePhases[cmssmEwsbEqs, cmssmEwsbOutputParameters];

TestEquality[cmssmFreePhases, {FlexibleSUSY`Phase[Mu]}];

TestEquality[Sort[cmssmSolution],
             Sort[{ B[Mu] -> -5 + x^2 - x*y - z,
                    Mu -> FlexibleSUSY`Phase[Mu] Sqrt[-5 - x^2 - x*y - z]
                  }]];

Print["testing CMSSMCPV-like EWSB for |Mu|, Re[BMu] and Im[BMu] ..."];

cmssmcpvEwsbEqs = {
    x - (E^(I*eta)*vu*B[Mu])/2 + vd*Mu*Susyno`LieGroups`conj[Mu] -
    (vu*Susyno`LieGroups`conj[B[Mu]])/(2*E^(I*eta)),
    y - (E^(I*eta)*vd*B[Mu])/2 + vu*Mu*Susyno`LieGroups`conj[Mu] -
    (vd*Susyno`LieGroups`conj[B[Mu]])/(2*E^(I*eta)),
    (-I/2)*E^(I*eta)*vd*B[Mu] + ((I/2)*vd*Susyno`LieGroups`conj[B[Mu]])/E^(I*eta)
};

cmssmcpvEwsbOutputParameters = { Mu, Re[B[Mu]], Im[B[Mu]] };

TestEquality[Parameters`IsRealParameter[Mu], False];
TestEquality[Parameters`IsRealParameter[B[Mu]], False];

Print["\t calling FindSolution[] ..."];

cmssmcpvFullSolution = EWSB`Private`FindSolution[cmssmcpvEwsbEqs, cmssmcpvEwsbOutputParameters];

solutionForMuMathematica7 = (Sqrt[-((vd*vu*x)/(vd^2 - vu^2)) - y + (vd^2*y)/(vd^2 - vu^2)]/Sqrt[vu]);
solutionForMuMathematica8 = (Sqrt[-(vd*x) + vu*y]/Sqrt[vd^2 - vu^2]);
solutionForMu = Which[$VersionNumber <= 7., solutionForMuMathematica7,
                      $VersionNumber  > 7., solutionForMuMathematica8
                     ];

TestEquality[Sort /@ cmssmcpvFullSolution,
             Sort /@ {{Mu        -> - solutionForMu,
                       Re[B[Mu]] -> ((1 + E^((2*I)*eta))*(-(vu*x) + vd*y))/(2*E^(I*eta)*(vd^2 - vu^2)),
                       Im[B[Mu]] -> (I/2*(-1 + E^((2*I)*eta))*(-(vu*x) + vd*y))/(E^(I*eta)*(vd^2 - vu^2))},
                      {Mu        -> + solutionForMu,
                       Re[B[Mu]] -> ((1 + E^((2*I)*eta))*(-(vu*x) + vd*y))/(2*E^(I*eta)*(vd^2 - vu^2)),
                       Im[B[Mu]] -> (I/2*(-1 + E^((2*I)*eta))*(-(vu*x) + vd*y))/(E^(I*eta)*(vd^2 - vu^2))}
                     }];

Print["\t calling FindSolutionAndFreePhases[] ..."];

{cmssmcpvSolution, cmssmcpvFreePhases} = EWSB`FindSolutionAndFreePhases[cmssmcpvEwsbEqs, cmssmcpvEwsbOutputParameters];

TestEquality[cmssmcpvFreePhases, {FlexibleSUSY`Phase[Mu]}];

TestEquality[Sort[cmssmcpvSolution],
             Sort[{Mu        -> solutionForMu * FlexibleSUSY`Phase[Mu],
                   Im[B[Mu]] -> (I/2*(-1 + E^((2*I)*eta))*(-(vu*x) + vd*y))/(E^(I*eta)*(vd^2 - vu^2)),
                   Re[B[Mu]] -> ((1 + E^((2*I)*eta))*(-(vu*x) + vd*y))/(2*E^(I*eta)*(vd^2 - vu^2))
                  }]];

Print["testing EWSB for mHu2, mHd2 ..."];

solution = EWSB`Private`TimeConstrainedSolve[{a + b - 2 == 0, a - b == 0}, {a,b}];

TestEquality[Sort[solution], Sort[{{a -> 1, b -> 1}}]];

solution = EWSB`Private`TimeConstrainedSolve[{a + b - 2 == 0, a - b == 0, c == 1}, {a,b}];

TestEquality[Sort[solution],
             Sort[{{a -> 1, b -> 1}}]
            ];

solution = EWSB`Private`TimeConstrainedSolve[{a - 2 == 0, b - 1 == 0, c == 1}, {a,b}];

TestEquality[Sort[solution],
             Sort[{{a -> 2, b -> 1}}]
            ];

(* test case for the MSSM/CPV:

   Here we have three linear independent equations, which we want to
   solve for two the parametes mHu2, mHd2
 *)

Print["testing MSSM/CPV EWSB for mHu2, mHd2"];

solution = EWSB`Private`TimeConstrainedSolve[
    EWSB`FilterOutIndependentEqs[
        {mHd2*vd + x - (E^(I*eta)*vu*B[Mu])/2 - (vu*Susyno`LieGroups`conj[B[Mu]])/(2*E^(I*eta)) == 0,
         mHu2*vu - y - (E^(I*eta)*vd*B[Mu])/2 - (vd*Susyno`LieGroups`conj[B[Mu]])/(2*E^(I*eta)) == 0,
         -I/2*E^(I*eta)*vd*B[Mu] + (I/2*vd*Susyno`LieGroups`conj[B[Mu]])/E^(I*eta) == 0},
        {mHd2, mHu2}
    ]
    ,
    {mHd2, mHu2}
];

TestEquality[Sort[solution],
             Sort[{{mHd2 -> (-2*E^(I*eta)*x + E^((2*I)*eta)*vu*B[Mu] + vu*conj[B[Mu]])/(2*E^(I*eta)*vd),
                    mHu2 -> ( 2*E^(I*eta)*y + E^((2*I)*eta)*vd*B[Mu] + vd*conj[B[Mu]])/(2*E^(I*eta)*vu)}}]
            ];

(* test case for the NMSSM/CPV *)

Print["testing NMSSM/CPV EWSB for mHd2, Im[T[\[Kappa]]], Re[T[\[Kappa]]], Im[T[\[Lambda]]], Re[T[\[Lambda]]]"];

nmssmcpvEWSBEqs =
{mHd2*vd + x - (E^(I*eta - (2*I)*etaS)*vS^2*vu*\[Lambda]*
    conj[\[Kappa]])/4 - (E^((-I)*eta + (2*I)*etaS)*vS^2*vu*\[Kappa]*
    conj[\[Lambda]])/4 + (vd*vS^2*\[Lambda]*conj[\[Lambda]])/2 +
  (vd*vu^2*\[Lambda]*conj[\[Lambda]])/2 -
  (E^((-I)*eta - I*etaS)*vS*vu*conj[T[\[Lambda]]])/(2*Sqrt[2]) -
  (E^(I*eta + I*etaS)*vS*vu*T[\[Lambda]])/(2*Sqrt[2]),
 mHu2*vu + y - (E^(I*eta - (2*I)*etaS)*vd*vS^2*\[Lambda]*conj[\[Kappa]])/
   4 - (E^((-I)*eta + (2*I)*etaS)*vd*vS^2*\[Kappa]*conj[\[Lambda]])/4 +
  (vd^2*vu*\[Lambda]*conj[\[Lambda]])/2 + (vS^2*vu*\[Lambda]*conj[\[Lambda]])/
   2 - (E^((-I)*eta - I*etaS)*vd*vS*conj[T[\[Lambda]]])/(2*Sqrt[2]) -
  (E^(I*eta + I*etaS)*vd*vS*T[\[Lambda]])/(2*Sqrt[2]),
 ms2*vS + vS^3*\[Kappa]*conj[\[Kappa]] -
  (E^(I*eta - (2*I)*etaS)*vd*vS*vu*\[Lambda]*conj[\[Kappa]])/2 -
  (E^((-I)*eta + (2*I)*etaS)*vd*vS*vu*\[Kappa]*conj[\[Lambda]])/2 +
  (vd^2*vS*\[Lambda]*conj[\[Lambda]])/2 + (vS*vu^2*\[Lambda]*conj[\[Lambda]])/
   2 + (vS^2*conj[T[\[Kappa]]])/(2*Sqrt[2]*E^((3*I)*etaS)) -
  (E^((-I)*eta - I*etaS)*vd*vu*conj[T[\[Lambda]]])/(2*Sqrt[2]) +
  (E^((3*I)*etaS)*vS^2*T[\[Kappa]])/(2*Sqrt[2]) -
  (E^(I*eta + I*etaS)*vd*vu*T[\[Lambda]])/(2*Sqrt[2]),
 (-I/4)*E^(I*eta - (2*I)*etaS)*vS^2*vu*\[Lambda]*conj[\[Kappa]] +
  (I/4)*E^((-I)*eta + (2*I)*etaS)*vS^2*vu*\[Kappa]*conj[\[Lambda]] +
  ((I/2)*E^((-I)*eta - I*etaS)*vS*vu*conj[T[\[Lambda]]])/Sqrt[2] -
  ((I/2)*E^(I*eta + I*etaS)*vS*vu*T[\[Lambda]])/Sqrt[2],
 (-I/4)*E^(I*eta - (2*I)*etaS)*vd*vS^2*\[Lambda]*conj[\[Kappa]] +
  (I/4)*E^((-I)*eta + (2*I)*etaS)*vd*vS^2*\[Kappa]*conj[\[Lambda]] +
  ((I/2)*E^((-I)*eta - I*etaS)*vd*vS*conj[T[\[Lambda]]])/Sqrt[2] -
  ((I/2)*E^(I*eta + I*etaS)*vd*vS*T[\[Lambda]])/Sqrt[2],
 (I/2)*E^(I*eta - (2*I)*etaS)*vd*vS*vu*\[Lambda]*conj[\[Kappa]] -
  (I/2)*E^((-I)*eta + (2*I)*etaS)*vd*vS*vu*\[Kappa]*conj[\[Lambda]] -
  ((I/2)*vS^2*conj[T[\[Kappa]]])/(Sqrt[2]*E^((3*I)*etaS)) +
  ((I/2)*E^((-I)*eta - I*etaS)*vd*vu*conj[T[\[Lambda]]])/Sqrt[2] +
  ((I/2)*E^((3*I)*etaS)*vS^2*T[\[Kappa]])/Sqrt[2] -
  ((I/2)*E^(I*eta + I*etaS)*vd*vu*T[\[Lambda]])/Sqrt[2]};

nmssmcpvEWSBOutputParameters = { mHd2, Im[T[\[Kappa]]], Re[T[\[Kappa]]], Im[T[\[Lambda]]], Re[T[\[Lambda]]]};

Parameters`SetRealParameters[{vS}];
TestEquality[Parameters`IsRealParameter[Re[\[Kappa]]], True];
TestEquality[Parameters`IsRealParameter[Im[\[Kappa]]], True];
TestEquality[Parameters`IsRealParameter[Re[T[\[Kappa]]]], True];
TestEquality[Parameters`IsRealParameter[Im[T[\[Kappa]]]], True];
TestEquality[Parameters`IsRealParameter[vS], True];

solution = EWSB`Private`FindSolution[nmssmcpvEWSBEqs, nmssmcpvEWSBOutputParameters];

TestNonEquality[solution, {}];
TestNonEquality[solution, {{}}];
TestEquality[Length[solution], 1];
TestEquality[Length[solution[[1]]], 5];

Print["testing NMSSM/CPV EWSB for mHd2, mHu2, ms2, Im[T[\[Kappa]]], Im[T[\[Lambda]]]"];

nmssmcpvEWSBOutputParameters = { mHd2, mHu2, ms2, Im[T[\[Kappa]]], Im[T[\[Lambda]]] };

Parameters`AddRealParameter[{mHd2, mHu2, ms2}];

nmssmcpvEWSBEqs = EWSB`FilterOutLinearDependentEqs[nmssmcpvEWSBEqs, nmssmcpvEWSBOutputParameters];
solution = EWSB`Private`FindSolution[nmssmcpvEWSBEqs, nmssmcpvEWSBOutputParameters];

TestNonEquality[solution, {}];
TestNonEquality[solution, {{}}];
TestEquality[Length[solution], 1];
TestEquality[Length[solution[[1]]], 5];

Print["testing EWSB substitutions ..."];

Parameters`SetModelParameters[{\[Mu], B[\[Mu]], mHd2, mHu2}];

subEwsbEqs = {
    \[Mu]^2 + x^2 + x y + z + 5,
    B[\[Mu]] - x^2 + x y + z + 5
};

ewsbSubs = {
   Rule[\[Mu], Sign[\[Mu]] Sqrt[MuSqr]]
};

subEwsbOutputParameters = { MuSqr, B[\[Mu]] };

Parameters`SetRealParameters[subEwsbOutputParameters];

{subSolution, subFreePhases} = EWSB`FindSolutionAndFreePhases[subEwsbEqs, subEwsbOutputParameters, "", ewsbSubs];

TestEquality[subFreePhases, {}];
TestEquality[Sort[Rule[#[[1]],Expand[#[[2]]]]& /@ subSolution],
             Sort[Rule[#[[1]],Expand[#[[2]]]]& /@ {MuSqr -> -(x^2 + x y + z + 5) / Sign[\[Mu]]^2,
                                                   B[\[Mu]] -> x^2 - x y - z - 5}]];

ewsbSubs = {
   Rule[\[Mu], Sign[\[Mu]] Sqrt[MuSqr]],
   Rule[B[\[Mu]], BMu0]
};

subEwsbOutputParameters = { MuSqr, BMu0 };

{subSolution, subFreePhases} = EWSB`FindSolutionAndFreePhases[subEwsbEqs, subEwsbOutputParameters, "", ewsbSubs];

TestEquality[subFreePhases, {}];
TestEquality[Sort[Rule[#[[1]],Expand[#[[2]]]]& /@ subSolution],
             Sort[Rule[#[[1]],Expand[#[[2]]]]& /@ {MuSqr -> -(x^2 + x y + z + 5) / Sign[\[Mu]]^2,
                                                   BMu0 -> x^2 - x y - z - 5}]];

Parameters`SetRealParameters[subEwsbOutputParameters];

subEwsbEqs = {
   mHd2 + x,
   mHu2 + y
};

ewsbSubs = {
   Rule[mHd2, m0^2 + m12^2],
   Rule[mHu2, m0^2 + AzeroSqr]
};

subEwsbOutputParameters = {m12, AzeroSqr};

Parameters`SetRealParameters[subEwsbOutputParameters];

{subSolution, subFreePhases} = EWSB`FindSolutionAndFreePhases[subEwsbEqs, subEwsbOutputParameters, "", ewsbSubs];

TestEquality[subFreePhases, {Sign[m12]}];
TestEquality[Sort[Rule[#[[1]],Expand[#[[2]]]]& /@ subSolution],
             Sort[Rule[#[[1]],Expand[#[[2]]]]& /@ {m12 -> Sign[m12] Sqrt[-(m0^2 + x)], AzeroSqr -> -(m0^2 + y)}]];

PrintTestSummary[];
