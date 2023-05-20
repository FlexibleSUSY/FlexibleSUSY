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
Needs["Vertices`", "Vertices.m"];

workingDirectory = Directory[];
SARAH`SARAH[OutputDirectory] = CreateDirectory[];
Print["Current working directory: ", workingDirectory];
Print["SARAH output directory: ", SARAH`SARAH[OutputDirectory]];

Start["MRSSM"];

On[Assert];
Off[EnforceCpColorStructures::cpext];

FlexibleSUSY`FSEigenstates = SARAH`EWSB;
Parameters`SetModelParameters[{
    Yd, Ye, LamTD, LamTU, LamSD, LamSU, Yu, \[Mu], MuD, MuU, g1, g2, g3,
    vd, vu, vT, vS, B[\[Mu]], B[MuD], B[MuU],
    mq2, ml2, mHd2, mHu2, md2, mu2, me2, mS2, mT2, moc2, mRd2, mRu2,
    MDBS, MDWBT, MDGoc}];

FH := Vertices`Private`FieldHead;
SFI := Vertices`StripFieldIndices;

CheckExternalField[extField_, cp:Cp[a_, b_, __]] /;
    FH[a] === FH[b] === extField &&
    SFI[a] === Susyno`LieGroups`conj@SFI[b] := (
    Print[extField, " correctly appears at the first two positions in ", cp];
    True
);

CheckExternalField[extField_, cp_Cp] := (
    Print[extField, " does not appear at the first two positions in ", cp];
    False
);

CheckExternalField[extField_, lst_List] :=
    And @@ (CheckExternalField[extField, #]& /@ Cases[lst, _Cp, Infinity]);

If[FreeQ[SARAH`Particles[EWSB], phiO | sigmaO] (* SARAH version < 4.5 *),

(* external fields are at wrong positions in following Cp[]'s *)

selfEnergySOc = {
    FSSelfEnergy[SOc,
     -C A0[Mass2[SOc]] Cp[conj[SOc], conj[SOc], SOc, SOc] -
      C sum[gI1, 1, 6,
        A0[Mass2[Sd[{gI1}]]] Cp[conj[Sd[{gI1}]], SOc, conj[SOc], Sd[{gI1}]]] -
      C sum[gI1, 1, 6,
        A0[Mass2[Su[{gI1}]]] Cp[conj[Su[{gI1}]], SOc, Su[{gI1}], conj[SOc]]]
    ]};

selfEnergySd = {
    FSSelfEnergy[Sd[gO1, gO2],
     -C A0[Mass2[SOc]]
	Cp[USd[{gO1}], conj[SOc], conj[USd[{gO2}]], SOc] -
      C sum[gI1, 1, 6,
        A0[Mass2[Sd[{gI1}]]]
	Cp[USd[{gO1}], conj[Sd[{gI1}]], Sd[{gI1}], conj[USd[{gO2}]]]] -
      C sum[gI1, 1, 6,
        A0[Mass2[Su[{gI1}]]]
	Cp[Su[{gI1}], conj[USd[{gO2}]], conj[Su[{gI1}]], USd[{gO1}]]]
    ]};

selfEnergySu = {
    FSSelfEnergy[Su[gO1, gO2],
     -C A0[Mass2[SOc]]
	Cp[conj[SOc], USu[{gO1}], conj[USu[{gO2}]], SOc] -
      C sum[gI1, 1, 6,
        A0[Mass2[Sd[{gI1}]]]
	Cp[Sd[{gI1}], conj[USu[{gO2}]], USu[{gO1}], conj[Sd[{gI1}]]]] -
      C sum[gI1, 1, 6,
        A0[Mass2[Su[{gI1}]]]
	Cp[Su[{gI1}], conj[Su[{gI1}]], conj[USu[{gO2}]], USu[{gO1}]]]
    ]};

massMatrices = {
    FSMassMatrix[dummy, SOc, Null],
    FSMassMatrix[dummy, Sd, ZD],
    FSMassMatrix[dummy, Su, ZU]
};

Print["testing EnforceCpColorStructures[] ..."];

selfEnergySOc = EnforceCpColorStructures[selfEnergySOc];
selfEnergySd  = EnforceCpColorStructures[selfEnergySd ];
selfEnergySu  = EnforceCpColorStructures[selfEnergySu ];

TestEquality[CheckExternalField[SOc, selfEnergySOc], True];
TestEquality[CheckExternalField[USd, selfEnergySd ], True];
TestEquality[CheckExternalField[USu, selfEnergySu ], True];

Print["testing color summation ..."];

vertexRulesSOc = VertexRules[selfEnergySOc, massMatrices];
vertexRulesSd  = VertexRules[selfEnergySd , massMatrices];
vertexRulesSu  = VertexRules[selfEnergySu , massMatrices];

TestEquality[Cp[SOc, conj[SOc], conj[SOc], SOc] /. vertexRulesSOc,
    -3*g3^2];

TestEquality[Cp[SOc, conj[SOc], conj[Sd[{gI1}]], Sd[{gI2}]] /. vertexRulesSOc,
    0];

TestEquality[Cp[SOc, conj[SOc], conj[Su[{gI1}]], Su[{gI2}]] /. vertexRulesSOc,
    0];

TestEquality[Cp[USd[{gO1}], conj[USd[{gO2}]], conj[SOc], SOc] /. vertexRulesSd,
    0];

TestEquality[Cp[USu[{gO1}], conj[USu[{gO2}]], conj[SOc], SOc] /. vertexRulesSu,
    0]
,
(*
 * SARAH version >= 4.5 eliminates from self-energy expressions
 * vertices which vanish after color summation
 *)

selfEnergyPhiO = {
    FSSelfEnergy[phiO,
     -(C*A0[Mass2[sigmaO]]*Cp[phiO, sigmaO, phiO, sigmaO])/2
    ]};

selfEnergySigmaO = {
    FSSelfEnergy[sigmaO,
     -(C*A0[Mass2[phiO]]*Cp[phiO, phiO, sigmaO, sigmaO])/2
    ]};

selfEnergySd = {
    FSSelfEnergy[Sd[gO1, gO2],
     -C sum[gI1, 1, 6,
        A0[Mass2[Sd[{gI1}]]]
	Cp[USd[{gO1}], conj[Sd[{gI1}]], Sd[{gI1}], conj[USd[{gO2}]]]] -
      C sum[gI1, 1, 6,
        A0[Mass2[Su[{gI1}]]]
	Cp[Su[{gI1}], conj[USd[{gO2}]], conj[Su[{gI1}]], USd[{gO1}]]]
    ]};

selfEnergySu = {
    FSSelfEnergy[Su[gO1, gO2],
     -C sum[gI1, 1, 6,
        A0[Mass2[Sd[{gI1}]]]
	Cp[Sd[{gI1}], conj[USu[{gO2}]], USu[{gO1}], conj[Sd[{gI1}]]]] -
      C sum[gI1, 1, 6,
        A0[Mass2[Su[{gI1}]]]
	Cp[Su[{gI1}], conj[Su[{gI1}]], conj[USu[{gO2}]], USu[{gO1}]]]
    ]};

massMatrices = {
    FSMassMatrix[dummy, SOc, Null],
    FSMassMatrix[dummy, Sd, ZD],
    FSMassMatrix[dummy, Su, ZU]
};

Print["testing EnforceCpColorStructures[] ..."];

selfEnergyPhiO   = EnforceCpColorStructures[selfEnergyPhiO  ];
selfEnergySigmaO = EnforceCpColorStructures[selfEnergySigmaO];
selfEnergySd  	 = EnforceCpColorStructures[selfEnergySd    ];
selfEnergySu  	 = EnforceCpColorStructures[selfEnergySu    ];

TestEquality[CheckExternalField[phiO  , selfEnergyPhiO  ], True];
TestEquality[CheckExternalField[sigmaO, selfEnergySigmaO], True];
TestEquality[CheckExternalField[USd   , selfEnergySd    ], True];
TestEquality[CheckExternalField[USu   , selfEnergySu    ], True];

Print["testing color summation ..."];

vertexRulesPhiO   = VertexRules[selfEnergyPhiO  , massMatrices];
vertexRulesSigmaO = VertexRules[selfEnergySigmaO, massMatrices];

TestEquality[Cp[phiO, phiO, sigmaO, sigmaO] /. vertexRulesPhiO,
    -6*g3^2];

TestEquality[Cp[sigmaO, sigmaO, phiO, phiO] /. vertexRulesSigmaO,
    -6*g3^2];
]

DeleteDirectory[SARAH`SARAH[OutputDirectory], DeleteContents -> True];

PrintTestSummary[];
