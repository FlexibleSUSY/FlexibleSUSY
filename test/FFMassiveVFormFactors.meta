(* ::Package:: *)

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

BeginPackage["FFMassiveVFormFactors`",
   {"SARAH`", "TextFormatting`", "TreeMasses`", "Vertices`", "CXXDiagrams`", "ColorMathInterface`"}
];

FFMassiveVFormFactorsCreateInterface::usage = "";
ff::usage = "";
MassiveVIndices::usage = "";

Begin["Private`"];

MyReIm[z_] := If[$VersionNumber >= 10.1,
      ReIm[z],
      {Re[z], Im[z]}
   ];

MyVertex[particles_List] := MyVertex[particles] = SARAH`Vertex[particles];

MassiveVIndices[V_] :=
   Module[{name = CXXNameOfField[V], dim = TreeMasses`GetDimension[V]},
      (* we expect that the gauge group from which massive vector bosons come from is broken,
         i.e. there is no index *)
      Assert[dim === 1];

      "template<>\n" <>
      "typename field_indices<" <> name <> ">::type\n" <>
      "default_indices_for_spectator<" <> name <> "> (void)\n" <>
      "{\n" <>
         "return {};\n" <>
      "}\n"
   ];

FFMassiveVFormFactorsCreateInterface[inFermion_ -> outFermion_, spectator_, loopParticles_List] :=
    Module[{prototype, definition,
            numberOfIndices1 = CXXDiagrams`NumberOfFieldIndices[inFermion],
            numberOfIndices2 = CXXDiagrams`NumberOfFieldIndices[outFermion],
            numberOfIndices3 = CXXDiagrams`NumberOfFieldIndices[spectator]},

      prototype =
         "std::valarray<std::complex<double>> calculate_" <> CXXNameOfField[inFermion] <>
            "_" <> CXXNameOfField[outFermion] <> "_" <> CXXNameOfField[spectator] <> "_form_factors" <>
            " (\n" <>
            If[TreeMasses`GetDimension[inFermion] =!= 1,
               "   int generationIndex1, ",
               " "
            ] <>
            If[TreeMasses`GetDimension[outFermion] =!= 1,
               " int generationIndex2, ",
               " "
            ] <>
            "const " <> FlexibleSUSY`FSModelName <> "_mass_eigenstates& model )";

      definition =
          prototype <> "{\n" <>
            IndentText[
               "context_base context {model};\n" <>
               "std::array<int, " <> ToString @ numberOfIndices1 <> "> indices1 = {" <>
                     (* TODO: Specify indices correctly *)
                       If[TreeMasses`GetDimension[inFermion] =!= 1,
                          " generationIndex1" <>
                          If[numberOfIndices1 =!= 1,
                             StringJoin @ Table[", 0", {numberOfIndices1-1}],
                             ""] <> " ",
                          If[numberOfIndices1 =!= 0,
                             StringJoin @ Riffle[Table[" 0", {numberOfIndices1}], ","] <> " ",
                             ""]
                         ] <> "};\n" <>
                   "std::array<int, " <> ToString @ numberOfIndices2 <>
                     "> indices2 = {" <>
                       If[TreeMasses`GetDimension[outFermion] =!= 1,
                          " generationIndex2" <>
                          If[numberOfIndices2 =!= 1,
                             StringJoin @ Table[", 0", {numberOfIndices2-1}],
                             ""] <> " ",
                          If[numberOfIndices2 =!= 0,
                             StringJoin @ Riffle[Table[" 0", {numberOfIndices2}], ","] <> " ",
                             ""]
                         ] <> "};\n\n" <>

               "std::valarray<std::complex<double>> val {0.0, 0.0};\n\n" <>

               StringJoin[
                  ("val += std::complex<double> " <> (ToString @ N[MyReIm@ColorN[#[[2,1]]], 16]) <> " * FFMassiveVVertexCorrectionFS<" <>
                   StringJoin @ Riffle[CXXDiagrams`CXXNameOfField /@ {inFermion, outFermion, spectator, #[[1,1]], #[[1,2]]}, ","]  <>
                   ">::value(indices1, indices2, context);\n") & /@ loopParticles
               ] <> "\n" <>

               "return val;"

            ] <> "\n}\n\n";

    {prototype <> ";", definition}
  ];

(* if t is TreeMasses`IsScalar then returns list of scalars and anti-scalar, etc. *)
getParticlesOfType[t_] :=
    DeleteDuplicates@Join[#, CXXDiagrams`LorentzConjugate /@ #]& @
      Select[TreeMasses`GetParticles[], t];

vertexNonZero[vertex_] :=
    Transpose[Drop[vertex, 1]][[1]] =!= {0,0};

vertexNonZeroS[vertex_] :=
    Transpose[Drop[vertex, 1]][[1]] =!= {0};

VertexIsNonZeroQ[vertex_] :=
   Module[{},
      True
   ];


singleMassiveDiagram[inFermion_, outFermion_, spectator_, F_?TreeMasses`IsFermion, S_?TreeMasses`IsScalar] :=
   Module[{FBarFjSBar, FiBarFS, SBarSVBar, FBarFVBar, FjBarFjVBar, FiBarFiVBar, v1, v2, v3, v4, v5, v6, colorIndexAssociation, p},

      On[Assert];
      (* calculation of color coefficients  for massive vector bosons is correct only if they are color singlets *)
      Assert[IsMassless[spectator] || !ColorChargedQ[spectator]];

      (* if the electric charge of an incomind particle doesn't equal to the sum of charges of outgoing ones,
         return an {} *)
      If[TreeMasses`GetElectricCharge[inFermion] =!= Plus @@ (TreeMasses`GetElectricCharge /@ {S,F}),
         Return[{}]
      ];

      colorIndexAssociation =
         If[TreeMasses`ColorChargedQ[#],
            If[TreeMasses`GetDimension[#] === 1,
               #[{Unique["ct"]}],
               #[{Unique["gt"], Unique["ct"]}]
            ], #
         ]& /@ {outFermion, F, S, F, S, inFermion, spectator};

      p = colorIndexAssociation;
      p = {p[[6]], p[[1]], p[[7]], p[[4]], p[[2]], p[[5]], p[[3]]};

      v1 = {CXXDiagrams`LorentzConjugate[F], inFermion, CXXDiagrams`LorentzConjugate[S]};
      FBarFjSBar = SARAHToColorMathSymbols@SARAH`Vertex[{CXXDiagrams`LorentzConjugate[p[[4]]], p[[1]], CXXDiagrams`LorentzConjugate[p[[6]]]}];
      v2 = {CXXDiagrams`LorentzConjugate[outFermion], F, S};
      FiBarFS = SARAHToColorMathSymbols@SARAH`Vertex[{CXXDiagrams`LorentzConjugate[p[[2]]], p[[5]], p[[7]]}];
      v3 = {CXXDiagrams`LorentzConjugate[S], S, CXXDiagrams`LorentzConjugate[spectator]};
      SBarSVBar = SARAHToColorMathSymbols@SARAH`Vertex[{CXXDiagrams`LorentzConjugate[p[[7]]], p[[6]], CXXDiagrams`LorentzConjugate[p[[3]]]}];
      v4 = {CXXDiagrams`LorentzConjugate[F], F, CXXDiagrams`LorentzConjugate[spectator]};
      FBarFVBar = SARAHToColorMathSymbols@SARAH`Vertex[{CXXDiagrams`LorentzConjugate[p[[5]]], p[[4]], CXXDiagrams`LorentzConjugate[p[[3]]]}];

      v5 = {CXXDiagrams`LorentzConjugate[inFermion], inFermion, CXXDiagrams`LorentzConjugate[spectator]};
      FjBarFjVBar = SARAHToColorMathSymbols@SARAH`Vertex[{CXXDiagrams`LorentzConjugate[inFermion], inFermion, CXXDiagrams`LorentzConjugate[spectator]}];

      v6 = {CXXDiagrams`LorentzConjugate[outFermion], outFermion, CXXDiagrams`LorentzConjugate[spectator]};
      FiBarFiVBar = SARAHToColorMathSymbols@SARAH`Vertex[{CXXDiagrams`LorentzConjugate[outFermion], outFermion, CXXDiagrams`LorentzConjugate[spectator]}];

      If[vertexNonZero[FBarFjSBar] && vertexNonZero[FiBarFS]
         && (vertexNonZeroS[SBarSVBar] || vertexNonZero[FBarFVBar] || vertexNonZero[FiBarFiVBar] || vertexNonZero[FjBarFjVBar]),
            Return[
               {StripSU3Generators[p[[1]], p[[2]], p[[3]],
                  ColorMath`CSimplify[CalculateColorFactor[{FBarFjSBar, FiBarFS, SBarSVBar}] ConnectColorLines[p[[5]], p[[4]]]]
                  ]
               , (*{v1, v2, v3}*){v1, v2, v3, v4}}
            ]
         ];

      Return[{}];
   ];


ff[inFermion_ -> outFermion_, spectator_] :=
   Module[{scalars, fermions, internalParticles = {}, temp},

      scalars = getParticlesOfType[TreeMasses`IsScalar];
      fermions = getParticlesOfType[TreeMasses`IsFermion];

      Map[
         (temp = singleMassiveDiagram[inFermion, outFermion, spectator, #[[1]], #[[2]]];
         If[temp =!= {},
            AppendTo[internalParticles, {#, temp}]
         ])&,
         Tuples[{fermions, scalars}]
      ];

      internalParticles
   ];

(* TODO: add other topologies? *)

End[];
EndPackage[];
