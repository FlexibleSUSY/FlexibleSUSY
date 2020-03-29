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

BeginPackage["SelfEnergies`", {"SARAH`", "TextFormatting`", "CConversion`", "TreeMasses`", "Parameters`", "Vertices`", "Utils`"}];

FSSelfEnergy::usage="self-energy head";
FSHeavySelfEnergy::usage="head for self-energy w/o BSM particles";
FSHeavyRotatedSelfEnergy::usage="head for self-energy w/o BSM particles in mass eigenstate basis";
Tadpole::usage="tadpole head";

GetField::usage="Returns field in self-energy or tadpole";

ConvertSarahSelfEnergies::usage="converts SARAH's self-energies to our
own format: SelfEnergies`FSSelfEnergy[particle, expression]";

ConvertSarahTadpoles::usage="converts SARAH's tadpoles to our own
format: SelfEnergies`Tadpole[particle, expression]";

CreateVertexExpressions::usage="creates C/C++ functions for the
given list of vertices";

CreateNPointFunctions::usage="creates C/C++ functions for the
given list of self-energies and tadpoles";

CreateSelfEnergyFunctionName::usage="creates self-energy function name
for a given field";

CreateHeavySelfEnergyFunctionName::usage="creates heavy self-energy
function name for a given field";

CreateHeavyRotatedSelfEnergyFunctionName::usage="creates heavy rotated
self-energy function name for a given field";

FillArrayWithLoopTadpoles::usage="add loop tadpoles to array"

FillArrayWithTwoLoopTadpoles::usage="add two-loop tadpoles to array"

DivideTadpoleByVEV::usage="Divides each tadpole by corresponding VEV";

CreateTwoLoopTadpolesMSSM::usage="Creates function prototypes and
definitions for two-loop tadpoles in the MSSM";

CreateTwoLoopTadpolesNMSSM::usage="Creates function prototypes and
definitions for two-loop tadpoles in the NMSSM";

CreateTwoLoopSelfEnergiesSM::usage="Creates function prototypes and
definitions for 2-loop Higgs self-energies in the SM";

CreateThreeLoopSelfEnergiesSM::usage="Creates function prototypes and
definitions for 3-loop Higgs self-energies in the SM";

CreateFourLoopSelfEnergiesSM::usage="Creates function prototypes and
definitions for 4-loop Higgs self-energies in the SM";

CreateTwoLoopSelfEnergiesMSSM::usage="Creates function prototypes and
definitions for two-loop Higgs self-energies in the MSSM";

CreateTwoLoopSelfEnergiesNMSSM::usage="Creates function prototypes and
definitions for two-loop Higgs self-energies in the NMSSM";

CreateThreeLoopSelfEnergiesMSSM::usage="Creates function prototypes
and definitions for three-loop Higgs contribution in the MSSM";

CreateThreeLoopSelfEnergiesSplit::usage="Creates function prototypes and
definitions for three-loop Higgs self-energies in split-SUSY";

SelfEnergyIsSymmetric::usage = "";
CreateCouplingSymbol::usage = "";

ReplaceGhosts::usage="";

Begin["`Private`"];

GetExpression[selfEnergy_SelfEnergies`FSSelfEnergy] :=
    selfEnergy[[2]];

GetExpression[selfEnergy_SelfEnergies`FSHeavySelfEnergy] :=
    selfEnergy[[2]];

GetExpression[selfEnergy_SelfEnergies`FSHeavyRotatedSelfEnergy] :=
    selfEnergy[[2]];

GetExpression[tadpole_SelfEnergies`Tadpole] :=
    tadpole[[2]];

GetField[selfEnergy_SelfEnergies`FSSelfEnergy] :=
    selfEnergy[[1]];

GetField[selfEnergy_SelfEnergies`FSHeavySelfEnergy] :=
    selfEnergy[[1]];

GetField[selfEnergy_SelfEnergies`FSHeavyRotatedSelfEnergy] :=
    selfEnergy[[1]];

GetField[tadpole_SelfEnergies`Tadpole] :=
    tadpole[[1]];

GetField[sym_] :=
    Module[{},
           Print["Error: GetField is not called with a SelfEnergies`FSSelfEnergy",
                 "or SelfEnergies`Tadpole head: ", sym];
           Quit[1];
          ];

SelfEnergyIsSymmetric[s_SelfEnergies`FSSelfEnergy] :=
    SelfEnergyIsSymmetric[GetField[s]];

SelfEnergyIsSymmetric[particle_] :=
    Length[Flatten[{FindMixingMatrixSymbolFor[particle]}]] === 1;

ExprContainsParticle[expr_, particle_List] :=
    Or @@ (ExprContainsParticle[expr,#]& /@ particle);

ExprContainsParticle[expr_, particle_] :=
    !FreeQ[expr,particle];

RemoveParticle[head_[p_,expr_], particle_] :=
    Module[{strippedExpr, a},
           strippedExpr = expr //. {
               SARAH`Cp[a__  /; ExprContainsParticle[{a},particle]][_] -> 0,
               SARAH`Cp[a__  /; ExprContainsParticle[{a},particle]] -> 0
                                   };
           head[p,strippedExpr]
          ];

RemoveSMParticles[SelfEnergies`FSSelfEnergy[p_,expr__], _] :=
    SelfEnergies`FSSelfEnergy[p,expr];

RemoveSMParticles[SelfEnergies`Tadpole[p_,expr__], _] :=
    SelfEnergies`Tadpole[p,expr];

ExprContainsNonOfTheseParticles[expr_, particles_List] :=
    And @@ (FreeQ[expr,#]& /@ particles);

RemoveSMParticles[head_[p_,expr_], removeGoldstones_:True, except_:{}] :=
    Module[{strippedExpr, keepParticles, a, goldstones,
            goldstonesWithoutIndex, g, i},
           keepParticles = Join[TreeMasses`GetSusyParticles[], except];
           goldstones = TreeMasses`GetSMGoldstoneBosons[];
           goldstonesWithoutIndex = goldstones /. a_[_] :> a;
           If[removeGoldstones,
              keepParticles = Complement[keepParticles, goldstonesWithoutIndex];,
              keepParticles = DeleteDuplicates[Join[keepParticles, goldstonesWithoutIndex]];
             ];
           strippedExpr = expr /. ReplaceGhosts[];
           strippedExpr = strippedExpr //. {
               SARAH`Cp[a__  /; ExprContainsNonOfTheseParticles[{a},keepParticles]][_] -> 0,
               SARAH`Cp[a__  /; ExprContainsNonOfTheseParticles[{a},keepParticles]] -> 0
                                           };
           (* remove goldstone bosons *)
           If[removeGoldstones,
              For[i = 1, i <= Length[goldstones], i++,
                  g = CConversion`GetHead[goldstones[[i]]];
                  strippedExpr = strippedExpr //.
                  SARAH`sum[idx_,_,endIdx_,expression_] /; !FreeQ[expression,g[{idx}]] :> SARAH`sum[idx,TreeMasses`GetDimensionStartSkippingSMGoldstones[g],endIdx,expression];
                 ];
             ];
           head[p,strippedExpr]
          ];

ReplaceUnrotatedFields[SelfEnergies`FSSelfEnergy[p_,expr__]] :=
    SelfEnergies`FSSelfEnergy[p,expr];

ReplaceUnrotatedFields[SelfEnergies`FSHeavySelfEnergy[p_,expr__]] :=
    SelfEnergies`FSHeavySelfEnergy[p,expr];

ReplaceUnrotatedFields[SelfEnergies`FSHeavyRotatedSelfEnergy[p_,expr__]] :=
    SelfEnergies`FSHeavyRotatedSelfEnergy[p, Sequence @@ (
        { expr } /. {
            SARAH`Cp[a__][l_] :> ReplaceUnrotatedFields[SARAH`Cp[a][l]],
            SARAH`Cp[a__]     :> ReplaceUnrotatedFields[SARAH`Cp[a]]
        }
    )];

CreateMassEigenstateReplacements[] :=
    Cases[Join[
             Flatten[SARAH`diracSubBack1 /@ SARAH`NameOfStates],
             Flatten[SARAH`diracSubBack2 /@ SARAH`NameOfStates]
          ],
          HoldPattern[Except[0] -> _]
    ];

AppendFieldIndices[lst_List, idx__] :=
    Module[{k, field, result = lst},
           For[k = 1, k <= Length[result], k++,
               field = GetField[result[[k]]];
               If[GetDimension[field] > 1,
                  result[[k,1]] = field[idx];
                 ];
              ];
           result
          ];

(* If the external field has dimension 1, remove it's indices.  For
   example in the Glu self-energy, terms appear of the form

      Cp[Glu[{gO2}], conj[Sd[{gI1}]], Fd[{gI2}]]

   Since Glu has dimension 1, the C variables Glu is a double and must
   therefore not be accessed in the form Glu(gO2).
 *)
Remove1DimensionalFieldIndices[lst_List] :=
    Module[{k, field, result = lst},
           For[k = 1, k <= Length[result], k++,
               field = GetHead[GetField[result[[k]]]];
               If[GetDimension[field] == 1,
                  result[[k,2]] = result[[k,2]] /. field[{__}] :> field;
                 ];
              ];
           result
          ];

(* decompose fermionic self-energies into L,R,S parts *)
SplitFermionSelfEnergies[lst_List] :=
    Module[{result = lst, k, field, expr, fermionSE},
           (* filter out all fermionic self-energies *)
           fermionSE = Cases[result, SelfEnergies`FSSelfEnergy[_, List[__]]];
           result = Select[result, (Head[GetExpression[#]] =!= List)&];
           (* decompose fermionic self-energies into L,R,S parts *)
           For[k = 1, k <= Length[fermionSE], k++,
               field = GetField[fermionSE[[k]]];
               expr = GetExpression[fermionSE[[k]]];
               If[Length[expr] != 3,
                  Print["Error: self-energy of ", field, " does not have 3 parts"];
                  Continue[];
                 ];
               AppendTo[result, SelfEnergies`FSSelfEnergy[field[1]       , expr[[1]]]];
               AppendTo[result, SelfEnergies`FSSelfEnergy[field[SARAH`PR], expr[[2]]]];
               AppendTo[result, SelfEnergies`FSSelfEnergy[field[SARAH`PL], expr[[3]]]];
              ];
           result
          ];

ConvertSarahTadpoles[DeleteLightFieldContrubtions[tadpoles_,_,_]] :=
    ConvertSarahTadpoles[tadpoles];

ConvertSarahTadpoles[tadpoles_List] :=
    Module[{result},
           result = (SelfEnergies`Tadpole @@@ tadpoles) /. CreateMassEigenstateReplacements[];
           result = AppendFieldIndices[result, SARAH`gO1];
           result /. SARAH`Mass -> FlexibleSUSY`M
          ];

ConvertSarahSelfEnergies[selfEnergies_List] :=
    Module[{result, heavySE,
            tQuark = TreeMasses`GetSMTopQuarkMultiplet[],
            bQuark = TreeMasses`GetSMBottomQuarkMultiplet[]
           },
           result = (SelfEnergies`FSSelfEnergy @@@ selfEnergies) /. CreateMassEigenstateReplacements[];
           result = AppendFieldIndices[result, SARAH`gO1, SARAH`gO2];
           result = SplitFermionSelfEnergies[result];
           result = Remove1DimensionalFieldIndices[result];
           (* Create Bottom, Tau self-energy with only SUSY
              particles and W and Z bosons in the loop *)
           heavySE = Cases[result, SelfEnergies`FSSelfEnergy[
               p:bQuark[__][_]|bQuark[_]|((particle_[__][_]|particle_[_]) /; TreeMasses`IsSMChargedLepton[particle]), expr__] :>
                           SelfEnergies`FSHeavyRotatedSelfEnergy[p, expr]];
           result = Join[result,
                         ReplaceUnrotatedFields /@ (RemoveSMParticles[#,False,{SARAH`VectorZ,SARAH`VectorW,SARAH`HiggsBoson}]& /@ heavySE)];
           (* Create rotated Top self-energy with only SUSY
              particles and W, Z and photon bosons in the loop *)
           heavySE = Cases[result, SelfEnergies`FSSelfEnergy[p:tQuark[__][_]|tQuark[_], expr__] :>
                           SelfEnergies`FSHeavyRotatedSelfEnergy[p, expr]];
           result = Join[result,
                         ReplaceUnrotatedFields /@ (RemoveParticle[#,
                                                                   If[FlexibleSUSY`UseMSSMYukawa2Loop === True,
                                                                      {SARAH`VectorG,SARAH`Gluino},
                                                                      SARAH`VectorG
                                                                     ]
                                                                  ]& /@ heavySE)];
           (* Create unrotated Top self-energy with only SUSY
              particles and W, Z and photon bosons in the loop *)
           heavySE = Cases[result, SelfEnergies`FSSelfEnergy[p:tQuark[__][_]|tQuark[_], expr__] :>
                           SelfEnergies`FSHeavySelfEnergy[p, expr]];
           result = Join[result, RemoveParticle[#,
                                                If[FlexibleSUSY`UseMSSMYukawa2Loop === True,
                                                   {SARAH`VectorG,SARAH`Gluino},
                                                   SARAH`VectorG
                                                  ]
                                               ]& /@ heavySE];
           result /. SARAH`Mass -> FlexibleSUSY`M
          ];

GetParticleIndicesInCoupling[Cp[a__]] := Flatten[Cases[{a}, List[__], Infinity]];

GetParticleIndicesInCoupling[Cp[a__][_]] := GetParticleIndicesInCoupling[Cp[a]];

CreateCouplingSymbol[coupling_] :=
    Module[{symbol, indices},
           indices = GetParticleIndicesInCoupling[coupling];
           symbol = ToValidCSymbol[coupling /. a_[List[__]] :> a];
           symbol[Sequence @@ indices]
          ];

indexCount = 0;
MakeUniqueIdx[] :=
    Symbol["id" <> ToString[indexCount++]];

(* creates a C++ function that calculates a coupling
 *
 * Return: {prototypes_String, definitions_String, rules_List}
 *
 * prototypes is a string that contains all coupling function
 * prototypes.  definitions is a string that contains all coupling
 * function definitions.
 *
 * rules is a list of replacement rules of the form
 * { Cp[bar[UCha[{gO2}]], VZ, Cha[{gI2}]][PR] :>
 *   CpbarUChaVZChaPR[gO2, gI2],
 *   ...
 * }
 *)
CreateCouplingFunction[coupling_, expr_, inModelClass_] :=
    Module[{symbol, prototype = "", definition = "",
            indices = {}, body = "", cFunctionName = "", i,
            type, typeStr},
           indices = GetParticleIndicesInCoupling[coupling];
           symbol = CreateCouplingSymbol[coupling];
           cFunctionName = ToValidCSymbolString[GetHead[symbol]];
           cFunctionName = cFunctionName <> "(";
           For[i = 1, i <= Length[indices], i++,
               If[i > 1, cFunctionName = cFunctionName <> ", ";];
               cFunctionName = cFunctionName <> "int ";
               (* variable names must not be integers *)
               If[!IntegerQ[indices[[i]]] && !FreeQ[expr, indices[[i]]],
                  cFunctionName = cFunctionName <> ToValidCSymbolString[indices[[i]]];
                 ];
              ];
           cFunctionName = cFunctionName <> ")";
           If[Parameters`IsRealExpression[expr],
              type = CConversion`ScalarType[CConversion`realScalarCType];,
              type = CConversion`ScalarType[CConversion`complexScalarCType];];
           typeStr = CConversion`CreateCType[type];
           prototype = typeStr <> " " <> cFunctionName <> " const;\n";
           definition = typeStr <> " CLASSNAME::" <> cFunctionName <> " const\n{\n";
           body = If[inModelClass,
                     Parameters`CreateLocalConstRefsForInputParameters[expr, "LOCALINPUT"],
                     Parameters`CreateLocalConstRefs[expr]
                    ] <> "\n" <>
                  "const " <> typeStr <> " result = " <>
                  Parameters`ExpressionToString[expr] <> ";\n\n" <>
                  "return result;\n";

           body = IndentText[WrapLines[body]];
           definition = definition <> body <> "}\n";
           {prototype, definition,
            RuleDelayed @@ {Vertices`ToCpPattern[coupling], symbol}}
          ];

GetParticleList[Cp[a__]] := {a};

GetParticleList[Cp[a__][_]] := {a};

IsUnrotated[SARAH`bar[field_]] := IsUnrotated[field];

IsUnrotated[Susyno`LieGroups`conj[field_]] := IsUnrotated[field];

IsUnrotated[field_[__]] := IsUnrotated[field];

IsUnrotated[field_Symbol] := StringTake[ToString[field],1] == "U";

ToRotatedField[field_Symbol] :=
    Symbol[StringReplace[ToString[field], StartOfString ~~ "U" ~~ rest_ :> rest]];

ToRotatedField[SARAH`bar[field_]] := SARAH`bar[ToRotatedField[field]];

ToRotatedField[Susyno`LieGroups`conj[field_]] := Susyno`LieGroups`conj[ToRotatedField[field]];

ToRotatedField[field_List] := ToRotatedField /@ field;

ToRotatedField[field_[indices__]] := ToRotatedField[field][indices];

GetUnrotatedFields[coupling_] := Select[GetParticleList[coupling], IsUnrotated];

ReplaceUnrotatedFields[SARAH`Cp[p__]] :=
    Cp[Sequence @@ ToRotatedField[{p}]];

ReplaceUnrotatedFields[SARAH`Cp[p__][lorentz_]] :=
    ReplaceUnrotatedFields[Cp[p]][lorentz];

CreateVertexExpressions[vertexRules_List, inModelClass_:True] :=
    Module[{k, prototypes = "", defs = "", rules, coupling, expr,
            p, d, r, MakeIndex},
           MakeIndex[i_Integer] := MakeUniqueIdx[];
           MakeIndex[i_] := i;
           rules = Table[0, {Length[vertexRules]}];
           Utils`StartProgressBar[Dynamic[k], Length[vertexRules]];
           For[k = 1, k <= Length[vertexRules], k++,
               coupling = Vertices`ToCp[vertexRules[[k,1]]] /. p_[{idx__}] :> p[MakeIndex /@ {idx}];
               expr = vertexRules[[k,2]];
               Utils`UpdateProgressBar[k, Length[vertexRules]];
               {p,d,r} = CreateCouplingFunction[coupling, expr, inModelClass];
               prototypes = prototypes <> p;
               defs = defs <> d <> "\n";
               rules[[k]] = r;
              ];
           Utils`StopProgressBar[Length[vertexRules]];
           {prototypes, defs, Flatten[rules]}
          ];

ReplaceGhosts[states_:FlexibleSUSY`FSEigenstates] :=
    Module[{vectorBosons = {}, ghostStr, ghostSym, ghostCSym, ghosts = {}, k},
           vectorBosons = GetVectorBosons[states];
           For[k = 1, k <= Length[vectorBosons], k++,
               ghostStr = ToString[vectorBosons[[k]]];
               ghostStr = StringReplace[ghostStr, StartOfString ~~ "V" ~~ x_ :> "g" <> x];
               ghostSym = Symbol[ghostStr];
               ghostCSym = Symbol[ghostStr <> "C"];
               AppendTo[ghosts, ghostSym -> vectorBosons[[k]]];
               AppendTo[ghosts, ghostCSym -> vectorBosons[[k]]];
              ];
           Return[ghosts];
          ];

DeclareFieldIndices[field_Symbol] := "";

DeclareFieldIndices[field_[ind1_, ind2_]] :=
    ", int " <> ToValidCSymbolString[ind1] <>
    ", int " <> ToValidCSymbolString[ind2];

DeclareFieldIndices[field_[PL]] := DeclareFieldIndices[field];
DeclareFieldIndices[field_[PR]] := DeclareFieldIndices[field];
DeclareFieldIndices[field_[1]]  := DeclareFieldIndices[field];
DeclareFieldIndices[field_[ind_]] :=
    "int " <> ToValidCSymbolString[ind];

ExtractChiraility[field_[idx1_,idx2_]] := ExtractChiraility[field];
ExtractChiraility[field_[PL]]          := "_PL";
ExtractChiraility[field_[PR]]          := "_PR";
ExtractChiraility[field_[1]]           := "_1";
ExtractChiraility[field_[idx_]]        := ExtractChiraility[field];
ExtractChiraility[field_]              := "";

ExtractFieldName[field_[idx1_,idx2_]] := ExtractFieldName[field];
ExtractFieldName[field_[PL]]          := ExtractFieldName[field];
ExtractFieldName[field_[PR]]          := ExtractFieldName[field];
ExtractFieldName[field_[1]]           := ExtractFieldName[field];
ExtractFieldName[field_[idx_]]        := ExtractFieldName[field];
ExtractFieldName[field_]              := ToValidCSymbolString[field];

CreateSelfEnergyFunctionName[field_, loops_] :=
    "self_energy_" <> ExtractFieldName[field] <> "_" <> ToString[loops] <> "loop" <> ExtractChiraility[field];

CreateHeavySelfEnergyFunctionName[field_, loops_] :=
    "self_energy_" <> ExtractFieldName[field] <> "_" <> ToString[loops] <> "loop" <> ExtractChiraility[field] <> "_heavy";

CreateHeavyRotatedSelfEnergyFunctionName[field_, loops_] :=
    "self_energy_" <> ExtractFieldName[field] <> "_" <> ToString[loops] <> "loop" <> ExtractChiraility[field] <> "_heavy_rotated";

CreateTadpoleFunctionName[field_, loops_] :=
    "tadpole_" <> ExtractFieldName[field] <> "_" <> ToString[loops] <> "loop" <> ExtractChiraility[field];

CreateFunctionName[selfEnergy_SelfEnergies`FSSelfEnergy, loops_] :=
    CreateSelfEnergyFunctionName[GetField[selfEnergy], loops];

CreateFunctionName[selfEnergy_SelfEnergies`FSHeavySelfEnergy, loops_] :=
    CreateHeavySelfEnergyFunctionName[GetField[selfEnergy], loops];

CreateFunctionName[selfEnergy_SelfEnergies`FSHeavyRotatedSelfEnergy, loops_] :=
    CreateHeavyRotatedSelfEnergyFunctionName[GetField[selfEnergy], loops];

CreateFunctionName[tadpole_SelfEnergies`Tadpole, loops_] :=
    CreateTadpoleFunctionName[GetField[tadpole], loops];

CreateFunctionPrototype[tadpole_SelfEnergies`Tadpole, loops_] :=
    CreateFunctionName[tadpole, loops] <>
    "(" <> DeclareFieldIndices[GetField[tadpole]] <> ") const";

CreateFunctionPrototype[selfEnergy_, loops_] :=
    CreateFunctionName[selfEnergy, loops] <>
    "(" <> CreateCType[CConversion`ScalarType[CConversion`realScalarCType]] <> " p " <> DeclareFieldIndices[GetField[selfEnergy]] <> ") const";

CreateFunctionPrototypeMatrix[s_, loops_] :=
    CreateFunctionName[s, loops] <> "(double p) const";

ExpressionToStringSequentially[expr_Plus, heads_, result_String] :=
    StringJoin[(result <> " += " <> ExpressionToString[#,heads] <> ";\n")& /@ (List @@ expr)];

ExpressionToStringSequentially[expr_, heads_, result_String] :=
    result <> " = " <> ExpressionToString[expr, heads] <> ";\n";

(* decreases literal indices in SARAH couplings *)
DecreaseLiteralCouplingIndices[expr_, num_:1] :=
    Module[{DecIdxLit},
           DecIdxLit[p_[idx_Integer]]   := p[idx - num];
           DecIdxLit[p_[{idx_Integer}]] := p[{idx - num}];
           DecIdxLit[p_]                := p;
           expr /. {
               SARAH`Cp[a__][b_] :> SARAH`Cp[Sequence @@ (DecIdxLit /@ {a})][b],
               SARAH`Cp[a__]     :> SARAH`Cp[Sequence @@ (DecIdxLit /@ {a})]
           }
          ];

CreateNPointFunction[nPointFunction_, vertexRules_List] :=
    Module[{decl, expr, prototype, body, functionName},
           expr = GetExpression[nPointFunction];
           functionName = CreateFunctionPrototype[nPointFunction, 1];
           type = CConversion`CreateCType[CConversion`ScalarType[CConversion`complexScalarCType]];
           prototype = type <> " " <> functionName <> ";\n";
           decl = "\n" <> type <> " CLASSNAME::" <> functionName <> "\n{\n";
           body = type <> " result;\n\n" <>
                  ExpressionToStringSequentially[
                                     DecreaseLiteralCouplingIndices[expr] /.
                                     vertexRules /.
                                     a_[List[i__]] :> a[i] /.
                                     ReplaceGhosts[FlexibleSUSY`FSEigenstates] /.
                                     C -> 1,
                                     TreeMasses`GetParticles[], "result"]  <>
                  "\nreturn result * oneOver16PiSqr;";
           body = IndentText[WrapLines[body]];
           decl = decl <> body <> "\n}\n";
           Return[{prototype, decl}];
          ];

CreateNPointFunctionMatrix[_SelfEnergies`Tadpole] := { "", "" };

FillHermitianSelfEnergyMatrix[nPointFunction_, sym_String] :=
    Module[{field = GetField[nPointFunction], dim, name},
           dim = GetDimension[field];
           name = CreateFunctionName[nPointFunction, 1];
           "\
for (int i = 0; i < " <> ToString[dim] <> "; i++)
   for (int k = i; k < " <> ToString[dim] <> "; k++)
      " <> sym <> "(i, k) = " <> name <> "(p, i, k);

Hermitianize(" <> sym <> ");
"
          ];

FillGeneralSelfEnergyFunction[nPointFunction_, sym_String] :=
    Module[{field = GetField[nPointFunction], dim, name},
           dim = GetDimension[field];
           name = CreateFunctionName[nPointFunction, 1];
           "\
for (int i = 0; i < " <> ToString[dim] <> "; i++)
   for (int k = 0; k < " <> ToString[dim] <> "; k++)
      " <> sym <> "(i, k) = " <> name <> "(p, i, k);
"
          ];

FillSelfEnergyMatrix[nPointFunction_, sym_String] :=
    Module[{particle = GetField[nPointFunction]},
           Which[(IsScalar[particle] || IsVector[particle]) && SelfEnergyIsSymmetric[particle],
                 FillHermitianSelfEnergyMatrix[nPointFunction, sym],
                 True,
                 FillGeneralSelfEnergyFunction[nPointFunction, sym]
                ]
          ];

CreateNPointFunctionMatrix[nPointFunction_] :=
    Module[{dim, functionName, type, prototype, def},
           dim = GetDimension[GetField[nPointFunction]];
           If[dim == 1, Return[{ "", "" }]];
           functionName = CreateFunctionPrototypeMatrix[nPointFunction, 1];
           type = CConversion`CreateCType[CConversion`MatrixType[CConversion`complexScalarCType, dim, dim]];
           prototype = type <> " " <> functionName <> ";\n";
           def = "
" <> type <> " CLASSNAME::" <> functionName <> "
{
   " <> type <> " self_energy;

" <> IndentText[FillSelfEnergyMatrix[nPointFunction, "self_energy"]] <> "
   return self_energy;
}
";
           { prototype, def }
          ];

CreateNPointFunctions[nPointFunctions_List, vertexRules_List] :=
    Module[{prototypes = "", defs = "", vertexFunctionNames = {}, p, d,
            relevantVertexRules},
           (* create coupling functions for all vertices in the list *)
           Print["Converting vertex functions ..."];
           (* extract vertex rules needed for the given nPointFunctions *)
           relevantVertexRules = Cases[vertexRules, r:(Rule[a_,b_] /; !FreeQ[nPointFunctions,a]) :> r];
           {prototypes, defs, vertexFunctionNames} = CreateVertexExpressions[relevantVertexRules];
           (* creating n-point functions *)
           Print["Converting self energies ..."];
           Utils`StartProgressBar[Dynamic[k], Length[nPointFunctions]];
           For[k = 1, k <= Length[nPointFunctions], k++,
               Utils`UpdateProgressBar[k, Length[nPointFunctions]];
               {p,d} = CreateNPointFunction[nPointFunctions[[k]], vertexFunctionNames];
               prototypes = prototypes <> p;
               defs = defs <> d;
               {p,d} = CreateNPointFunctionMatrix[nPointFunctions[[k]]];
               prototypes = prototypes <> p;
               defs = defs <> d;
              ];
           Utils`StopProgressBar[Length[nPointFunctions]];
           {prototypes, defs}
          ];

FillArrayWithLoopTadpoles[loopLevel_, higgsAndIdx_List, arrayName_String, sign_String:"-", struct_String:""] :=
    Module[{body = "", v, field, idx, head, functionName},
           For[v = 1, v <= Length[higgsAndIdx], v++,
               field = higgsAndIdx[[v,1]];
               idx = higgsAndIdx[[v,2]];
               head = CConversion`ToValidCSymbolString[higgsAndIdx[[v,3]]];
               functionName = CreateTadpoleFunctionName[field, loopLevel];
               If[TreeMasses`GetDimension[field] == 1,
                  body = body <> arrayName <> "[" <> ToString[v-1] <> "] " <> sign <> "= " <>
                         head <> "(" <> struct <> functionName <> "());\n";
                  ,
                  body = body <> arrayName <> "[" <> ToString[v-1] <> "] " <> sign <> "= " <>
                         head <> "(" <> struct <> functionName <>
                         "(" <> ToString[idx - 1] <> "));\n";
                 ];
              ];
           body
          ];

FillArrayWithTwoLoopTadpoles[higgsBoson_, arrayName_String, sign_String:"-", struct_String:""] :=
    Module[{body, v, field, functionName, dim, dimStr},
           functionName = CreateTadpoleFunctionName[higgsBoson,2];
           dim = GetDimension[higgsBoson];
           dimStr = ToString[dim];
           body = "const auto tadpole_2l(" <> struct <> functionName <> "());\n";
           For[v = 1, v <= dim, v++,
               body = body <> arrayName <> "[" <> ToString[v-1] <> "] " <> sign <> "= " <>
                      "tadpole_2l(" <> ToString[v-1] <> ");\n";
              ];
           body
          ];

DivideTadpoleByVEV[higgsAndVEV_List, arrayName_String] :=
    Module[{body = "", v, vev},
           For[v = 1, v <= Length[higgsAndVEV], v++,
               vev = higgsAndVEV[[v,3]] higgsAndVEV[[v,4]];
               If[vev === 0,
                  body = body <> arrayName <> "[" <> ToString[v-1] <> "] = 0.;\n";,
                  vev = CConversion`RValueToCFormString[vev];
                  body = body <> arrayName <> "[" <> ToString[v-1] <> "] /= " <> vev <> ";\n";
                 ];
              ];
           body
          ];

AssertFieldDimension[field_, dim_, model_] :=
    Block[{fieldDim},
          fieldDim = GetDimension[field];
          If[fieldDim =!= dim,
             Print["Error: Cannot use ", model, " two-loop routines for ",
                   field, " (multiplet size ", fieldDim, ").  Multiplet size ",
                   dim, " required!"];
             Quit[1];
            ];
         ];

GetTwoLoopTadpoleCorrections[model_String /; model === "MSSM"] :=
    Module[{g3Str, mtStr, mbStr, mtauStr,
            mTop, mBot, mTau,
            vev2Str, tanbStr, muStr, m3Str, mA0Str},
           AssertFieldDimension[SARAH`HiggsBoson, 2, model];
           mTop    = TreeMasses`GetMass[TreeMasses`GetUpQuark[3,True]];
           mBot    = TreeMasses`GetMass[TreeMasses`GetDownQuark[3,True]];
           mTau    = TreeMasses`GetMass[TreeMasses`GetDownLepton[3,True]];
           mtStr   = CConversion`RValueToCFormString[mTop];
           mbStr   = CConversion`RValueToCFormString[mBot];
           mtauStr = CConversion`RValueToCFormString[mTau];
           g3Str   = CConversion`RValueToCFormString[SARAH`strongCoupling];
           vev2Str = CConversion`RValueToCFormString[SARAH`VEVSM1^2 + SARAH`VEVSM2^2];
           tanbStr = CConversion`RValueToCFormString[SARAH`VEVSM2 / SARAH`VEVSM1];
           muStr   = CConversion`RValueToCFormString[-Parameters`GetEffectiveMu[]];
           m3Str   = CConversion`RValueToCFormString[FlexibleSUSY`M[SARAH`Gluino]];
           mA0Str  = TreeMasses`CallPseudoscalarHiggsMassGetterFunction[] <> "(0)";
"\
using namespace flexiblesusy::mssm_twoloophiggs;

// calculate 3rd generation sfermion masses and mixing angles
double mst_1, mst_2, theta_t;
double msb_1, msb_2, theta_b;
double mstau_1, mstau_2, theta_tau;
double msnu_1, msnu_2, theta_nu;

" <> TreeMasses`CallGenerationHelperFunctionName[3, SARAH`TopSquark, "mst_1", "mst_2", "theta_t"] <>
";
" <> TreeMasses`CallGenerationHelperFunctionName[3, SARAH`BottomSquark, "msb_1", "msb_2", "theta_b"] <>
";
" <> TreeMasses`CallGenerationHelperFunctionName[3, SARAH`Selectron, "mstau_1", "mstau_2", "theta_tau"] <>
";
" <> TreeMasses`CallGenerationHelperFunctionName[3, SARAH`Sneutrino, "msnu_1", "msnu_2", "theta_nu"] <>
";

const double mst1sq = Sqr(mst_1), mst2sq = Sqr(mst_2);
const double msb1sq = Sqr(msb_1), msb2sq = Sqr(msb_2);
const double mstau1sq = Sqr(mstau_1), mstau2sq = Sqr(mstau_2);
const double msnusq = Sqr(msnu_2);
const double sxt = Sin(theta_t), cxt = Cos(theta_t);
const double sxb = Sin(theta_b), cxb = Cos(theta_b);
const double sintau = Sin(theta_tau), costau = Cos(theta_tau);
const double gs = " <> g3Str <> ";
const double rmtsq = Sqr(" <> mtStr <> ");
const double scalesq = Sqr(get_scale());
const double vev2 = " <> vev2Str <> ";
const double tanb = " <> tanbStr <> ";
const double amu = Re(" <> muStr <> ");
const double mg = " <> m3Str <> ";
const double mAsq = Sqr(" <> mA0Str <> ");
const double cotbeta = 1.0 / tanb;
const double rmbsq = Sqr(" <> mbStr <> ");
const double rmtausq = Sqr(" <> mtauStr <> ");

" <> GetTadpoleVectorCType[2] <> " tadpole_2l(" <> GetTadpoleVectorCType[2] <> "::Zero());

if (HIGGS_2LOOP_CORRECTION_AT_AS) {
   tadpole_2l += tadpole_higgs_2loop_at_as_mssm(
      rmtsq, mg, mst1sq, mst2sq, sxt, cxt, scalesq,
      amu, tanb, vev2, gs);
}

if (HIGGS_2LOOP_CORRECTION_AT_AT) {
   tadpole_2l += tadpole_higgs_2loop_at_at_mssm(
      rmtsq, rmbsq, mAsq, mst1sq, mst2sq, msb1sq, msb2sq,
      sxt, cxt, sxb, cxb, scalesq, amu, tanb, vev2);
}

if (HIGGS_2LOOP_CORRECTION_AB_AS) {
   tadpole_2l += tadpole_higgs_2loop_ab_as_mssm(
      rmbsq, mg, msb1sq, msb2sq, sxb, cxb, scalesq,
      amu, cotbeta, vev2, gs);
}

if (HIGGS_2LOOP_CORRECTION_ATAU_ATAU) {
   tadpole_2l += tadpole_higgs_2loop_atau_atau_mssm(
      rmtausq, mAsq, msnusq, mstau1sq, mstau2sq, sintau,
      costau, scalesq, amu, tanb, vev2);
}

tadpole_2l(0) *= " <> CConversion`ToValidCSymbolString[SARAH`VEVSM1] <> ";
tadpole_2l(1) *= " <> CConversion`ToValidCSymbolString[SARAH`VEVSM2] <> ";

return tadpole_2l;"
          ];

GetTwoLoopTadpoleCorrections[model_String /; model === "NMSSM"] :=
    Module[{mTop, mBot, mTau,
            g3Str, mtStr, mbStr, mtauStr, lambdaStr,
            vev2Str, svevStr, tanbStr, muStr, m3Str, mA0Str},
           AssertFieldDimension[SARAH`HiggsBoson, 3, model];
           mTop    = TreeMasses`GetMass[TreeMasses`GetUpQuark[3,True]];
           mBot    = TreeMasses`GetMass[TreeMasses`GetDownQuark[3,True]];
           mTau    = TreeMasses`GetMass[TreeMasses`GetDownLepton[3,True]];
           mtStr   = CConversion`RValueToCFormString[mTop];
           mbStr   = CConversion`RValueToCFormString[mBot];
           mtauStr = CConversion`RValueToCFormString[mTau];
           g3Str   = CConversion`RValueToCFormString[SARAH`strongCoupling];
           vev2Str = CConversion`RValueToCFormString[SARAH`VEVSM1^2 + SARAH`VEVSM2^2];
           tanbStr = CConversion`RValueToCFormString[SARAH`VEVSM2 / SARAH`VEVSM1];
           muStr   = CConversion`RValueToCFormString[-Parameters`GetEffectiveMu[]];
           m3Str   = CConversion`RValueToCFormString[FlexibleSUSY`M[SARAH`Gluino]];
           mA0Str  = CConversion`RValueToCFormString[Parameters`GetEffectiveMASqr[]];
           svevStr = CConversion`RValueToCFormString[Parameters`GetParameterFromDescription["Singlet-VEV"]];
           lambdaStr = CConversion`RValueToCFormString[Parameters`GetParameterFromDescription["Singlet-Higgs-Interaction"]];
"\
using namespace flexiblesusy::mssm_twoloophiggs;
using namespace flexiblesusy::nmssm_twoloophiggs;

// calculate 3rd generation sfermion masses and mixing angles
double mst_1, mst_2, theta_t;
double msb_1, msb_2, theta_b;
double mstau_1, mstau_2, theta_tau;
double msnu_1, msnu_2, theta_nu;

" <> TreeMasses`CallGenerationHelperFunctionName[3, SARAH`TopSquark, "mst_1", "mst_2", "theta_t"] <>
";
" <> TreeMasses`CallGenerationHelperFunctionName[3, SARAH`BottomSquark, "msb_1", "msb_2", "theta_b"] <>
";
" <> TreeMasses`CallGenerationHelperFunctionName[3, SARAH`Selectron, "mstau_1", "mstau_2", "theta_tau"] <>
";
" <> TreeMasses`CallGenerationHelperFunctionName[3, SARAH`Sneutrino, "msnu_1", "msnu_2", "theta_nu"] <>
";

const double mst1sq = Sqr(mst_1), mst2sq = Sqr(mst_2);
const double msb1sq = Sqr(msb_1), msb2sq = Sqr(msb_2);
const double mstau1sq = Sqr(mstau_1), mstau2sq = Sqr(mstau_2);
const double msnusq = Sqr(msnu_2);
const double sxt = Sin(theta_t), cxt = Cos(theta_t);
const double sxb = Sin(theta_b), cxb = Cos(theta_b);
const double sintau = Sin(theta_tau), costau = Cos(theta_tau);
const double gs = " <> g3Str <> ";
const double rmtsq = Sqr(" <> mtStr <> ");
const double scalesq = Sqr(get_scale());
const double vev2 = " <> vev2Str <> ";
const double tanb = " <> tanbStr <> ";
const double amu = Re(" <> muStr <> ");
const double mg = " <> m3Str <> ";
const double mAsq = " <> mA0Str <> ";
const double cotbeta = 1.0 / tanb;
const double rmbsq = Sqr(" <> mbStr <> ");
const double rmtausq = Sqr(" <> mtauStr <> ");
const double lam = Re(" <> lambdaStr <> ");
const double svev = Abs(amu / lam);

" <> GetTadpoleVectorCType[3] <> " tadpole_2l(" <> GetTadpoleVectorCType[3] <> "::Zero());

if (HIGGS_2LOOP_CORRECTION_AT_AS) {
   tadpole_2l += tadpole_higgs_2loop_at_as_nmssm(
      rmtsq, mg, mst1sq, mst2sq, sxt, cxt, scalesq,
      amu, tanb, vev2, gs, svev);
}

if (HIGGS_2LOOP_CORRECTION_AT_AT) {
   tadpole_2l.head<2>() += tadpole_higgs_2loop_at_at_mssm(
      rmtsq, rmbsq, mAsq, mst1sq, mst2sq, msb1sq, msb2sq,
      sxt, cxt, sxb, cxb, scalesq, amu, tanb, vev2);
}

if (HIGGS_2LOOP_CORRECTION_AB_AS) {
   tadpole_2l += tadpole_higgs_2loop_ab_as_nmssm(
      rmbsq, mg, msb1sq, msb2sq, sxb, cxb, scalesq,
      amu, cotbeta, vev2, gs, svev);
}

if (HIGGS_2LOOP_CORRECTION_ATAU_ATAU) {
   tadpole_2l.head<2>() += tadpole_higgs_2loop_atau_atau_mssm(
      rmtausq, mAsq, msnusq, mstau1sq, mstau2sq, sintau,
      costau, scalesq, amu, tanb, vev2);
}

tadpole_2l(0) *= " <> CConversion`ToValidCSymbolString[SARAH`VEVSM1] <> ";
tadpole_2l(1) *= " <> CConversion`ToValidCSymbolString[SARAH`VEVSM2] <> ";
tadpole_2l(2) *= " <>  svevStr <> ";

return tadpole_2l;"
          ];

GetTwoLoopTadpoleCorrections[model_] :=
    Module[{},
           Print["Error: two-loop tadpoles for ", model, " not available"];
           ""
          ];

GetTadpoleVectorType[dim_] :=
    CConversion`VectorType[CConversion`realScalarCType, dim];

GetTadpoleVectorCType[dim_] :=
    CConversion`CreateCType[GetTadpoleVectorType[dim]];

CreateTwoLoopTadpoles[higgsBoson_, model_String] :=
    Module[{prototype, function, functionName, dim, dimStr, cType},
           dim = GetDimension[higgsBoson];
           dimStr = ToString[dim];
           functionName = CreateTadpoleFunctionName[higgsBoson,2];
           cType = GetTadpoleVectorCType[dim];
           prototype = cType <> " " <> functionName <> "() const;\n";
           body = GetTwoLoopTadpoleCorrections[model];
           function = cType <> " CLASSNAME::" <> functionName <> "() const\n{\n" <>
                      IndentText[body] <> "\n}\n";
           Return[{prototype, function}];
          ];

CreateTwoLoopTadpolesMSSM[higgsBoson_] :=
    CreateTwoLoopTadpoles[higgsBoson, "MSSM"];

CreateTwoLoopTadpolesNMSSM[higgsBoson_] :=
    CreateTwoLoopTadpoles[higgsBoson, "NMSSM"];

GetNLoopSelfEnergyCorrections[particle_ /; particle === SARAH`HiggsBoson,
                              model_String /; model === "SM", 2] :=
    Module[{mtStr, ytStr, mbStr, ybStr, mtauStr, ytauStr, g3Str},
           AssertFieldDimension[particle, 1, model];
           mtStr   = CConversion`RValueToCFormString[TreeMasses`GetMass[TreeMasses`GetUpQuark[3,True]]];
           ytStr   = CConversion`RValueToCFormString[Parameters`GetThirdGeneration[SARAH`UpYukawa]];
           mbStr   = CConversion`RValueToCFormString[TreeMasses`GetMass[TreeMasses`GetDownQuark[3,True]]];
           ybStr   = CConversion`RValueToCFormString[Parameters`GetThirdGeneration[SARAH`DownYukawa]];
           mtauStr = CConversion`RValueToCFormString[TreeMasses`GetMass[TreeMasses`GetDownLepton[3,True]]];
           ytauStr = CConversion`RValueToCFormString[Parameters`GetThirdGeneration[SARAH`ElectronYukawa]];
           g3Str   = CConversion`RValueToCFormString[SARAH`strongCoupling];
"\
using namespace flexiblesusy::sm_twoloophiggs;

const double p2 = Sqr(p);
const double mb = " <> mbStr <> ";
const double mt = " <> mtStr <> ";
const double mtau = " <> mtauStr <> ";
const double yb = " <> ybStr <> ";
const double yt = " <> ytStr <> ";
const double ytau = " <> ytauStr <> ";
const double gs = " <> g3Str <> ";
const double scale = get_scale();
double self_energy = 0.;

if (HIGGS_2LOOP_CORRECTION_AT_AT) {
   self_energy -= delta_mh_2loop_at_at_sm(p2, scale, mt, yt, mb);
}

if (HIGGS_2LOOP_CORRECTION_AT_AS) {
   self_energy -= delta_mh_2loop_at_as_sm(p2, scale, mt, yt, gs);
}

if (HIGGS_2LOOP_CORRECTION_AB_AS) {
   self_energy -= delta_mh_2loop_ab_as_sm(p2, scale, mb, yb, gs);
}

if (HIGGS_2LOOP_CORRECTION_ATAU_ATAU) {
   self_energy -= delta_mh_2loop_atau_atau_sm(p2, scale, mtau, ytau);
}

return self_energy;"
          ];

GetNLoopSelfEnergyCorrections[particle_ /; particle === SARAH`HiggsBoson,
                              model_String /; model === "SM", 3] :=
    Module[{mTop, mtStr, yt, ytStr, g3Str, mHiggs, mhStr},
           AssertFieldDimension[particle, 1, model];
           mTop    = TreeMasses`GetMass[TreeMasses`GetUpQuark[3,True]];
           mtStr   = CConversion`RValueToCFormString[mTop];
           yt      = Parameters`GetThirdGeneration[SARAH`UpYukawa];
           ytStr   = CConversion`RValueToCFormString[yt];
           g3Str   = CConversion`RValueToCFormString[SARAH`strongCoupling];
           mHiggs  = TreeMasses`GetMass[particle];
           mhStr   = CConversion`RValueToCFormString[mHiggs];
"\
using namespace flexiblesusy::sm_threeloophiggs;

const double mt = " <> mtStr <> ";
const double yt = " <> ytStr <> ";
const double gs = " <> g3Str <> ";
const double mh = " <> mhStr <> ";
const double scale = get_scale();
double self_energy = 0.;

if (HIGGS_3LOOP_CORRECTION_AT_AT_AT) {
   self_energy -= delta_mh_3loop_at_at_at_sm(scale, mt, yt, mh);
}

if (HIGGS_3LOOP_CORRECTION_AT_AT_AS) {
   self_energy -= delta_mh_3loop_at_at_as_sm(scale, mt, yt, gs);
}

if (HIGGS_3LOOP_CORRECTION_AT_AS_AS) {
   self_energy -= delta_mh_3loop_at_as_as_sm(scale, mt, yt, gs);
}

return self_energy;"
          ];

GetNLoopSelfEnergyCorrections[particle_ /; particle === SARAH`HiggsBoson,
                              model_String /; model === "SM", 4] :=
    Module[{mTop, mtStr, yt, ytStr, g3Str, mHiggs, mhStr},
           AssertFieldDimension[particle, 1, model];
           mTop    = TreeMasses`GetMass[TreeMasses`GetUpQuark[3,True]];
           mtStr   = CConversion`RValueToCFormString[mTop];
           yt      = Parameters`GetThirdGeneration[SARAH`UpYukawa];
           ytStr   = CConversion`RValueToCFormString[yt];
           g3Str   = CConversion`RValueToCFormString[SARAH`strongCoupling];
           mHiggs  = TreeMasses`GetMass[particle];
           mhStr   = CConversion`RValueToCFormString[mHiggs];
"\
using namespace flexiblesusy::sm_fourloophiggs;

const double mt = " <> mtStr <> ";
const double yt = " <> ytStr <> ";
const double gs = " <> g3Str <> ";
const double scale = get_scale();
double self_energy = 0.;

if (HIGGS_4LOOP_CORRECTION_AT_AS_AS_AS) {
   self_energy -= delta_mh_4loop_at_as_as_as_sm(scale, mt, yt, gs);
}

return self_energy;"
          ];

GetNLoopSelfEnergyCorrections[particle_ /; particle === SARAH`HiggsBoson,
                              model_String /; model === "Split", 3] :=
    Module[{mTop, mGluino, mtStr, mgStr, yt, ytStr, g3Str},
           AssertFieldDimension[particle, 1, model];
           mTop    = TreeMasses`GetMass[TreeMasses`GetUpQuark[3,True]];
           mtStr   = CConversion`RValueToCFormString[mTop];
           mGluino = TreeMasses`GetMass[Parameters`GetParticleFromDescription["Gluino"]];
           mgStr   = CConversion`RValueToCFormString[mGluino];
           yt      = Parameters`GetThirdGeneration[SARAH`UpYukawa];
           ytStr   = CConversion`RValueToCFormString[yt];
           g3Str   = CConversion`RValueToCFormString[SARAH`strongCoupling];
"\
using namespace flexiblesusy::splitmssm_threeloophiggs;

const double mt = " <> mtStr <> ";
const double mg = " <> mgStr <> ";
const double yt = " <> ytStr <> ";
const double gs = " <> g3Str <> ";
const double scale = get_scale();
double self_energy = 0.;

if (HIGGS_3LOOP_CORRECTION_AT_AS_AS) {
   self_energy -= delta_mh_3loop_gluino_split(scale, mt, yt, gs, mg);
}

return self_energy;"
          ];

GetNLoopSelfEnergyCorrections[particle_ /; particle === SARAH`HiggsBoson,
                              model_String /; model === "MSSM", 2] :=
    Module[{g3Str, mtStr, mbStr, mtauStr,
            mTop, mBot, mTau,
            vev2Str, vuStr, vdStr, tanbStr, muStr, m3Str, mA0Str},
           AssertFieldDimension[particle, 2, model];
           mTop    = TreeMasses`GetMass[TreeMasses`GetUpQuark[3,True]];
           mBot    = TreeMasses`GetMass[TreeMasses`GetDownQuark[3,True]];
           mTau    = TreeMasses`GetMass[TreeMasses`GetDownLepton[3,True]];
           mtStr   = CConversion`RValueToCFormString[mTop];
           mbStr   = CConversion`RValueToCFormString[mBot];
           mtauStr = CConversion`RValueToCFormString[mTau];
           g3Str   = CConversion`RValueToCFormString[SARAH`strongCoupling];
           vev2Str = CConversion`RValueToCFormString[SARAH`VEVSM1^2 + SARAH`VEVSM2^2];
           vdStr   = CConversion`RValueToCFormString[SARAH`VEVSM1];
           vuStr   = CConversion`RValueToCFormString[SARAH`VEVSM2];
           tanbStr = CConversion`RValueToCFormString[SARAH`VEVSM2 / SARAH`VEVSM1];
           muStr   = CConversion`RValueToCFormString[-Parameters`GetEffectiveMu[]];
           m3Str   = CConversion`RValueToCFormString[FlexibleSUSY`M[SARAH`Gluino]];
           mA0Str  = TreeMasses`CallPseudoscalarHiggsMassGetterFunction[] <> "(0)";
"\
using namespace flexiblesusy::mssm_twoloophiggs;

// calculate 3rd generation sfermion masses and mixing angles
double mst_1, mst_2, theta_t;
double msb_1, msb_2, theta_b;
double mstau_1, mstau_2, theta_tau;
double msnu_1, msnu_2, theta_nu;

" <> TreeMasses`CallGenerationHelperFunctionName[3, SARAH`TopSquark, "mst_1", "mst_2", "theta_t"] <>
";
" <> TreeMasses`CallGenerationHelperFunctionName[3, SARAH`BottomSquark, "msb_1", "msb_2", "theta_b"] <>
";
" <> TreeMasses`CallGenerationHelperFunctionName[3, SARAH`Selectron, "mstau_1", "mstau_2", "theta_tau"] <>
";
" <> TreeMasses`CallGenerationHelperFunctionName[3, SARAH`Sneutrino, "msnu_1", "msnu_2", "theta_nu"] <>
";

const double mst1sq = Sqr(mst_1), mst2sq = Sqr(mst_2);
const double msb1sq = Sqr(msb_1), msb2sq = Sqr(msb_2);
const double mstau1sq = Sqr(mstau_1), mstau2sq = Sqr(mstau_2);
const double msnusq = Sqr(msnu_2);
const double sxt = Sin(theta_t), cxt = Cos(theta_t);
const double sxb = Sin(theta_b), cxb = Cos(theta_b);
const double sintau = Sin(theta_tau), costau = Cos(theta_tau);
const double gs = " <> g3Str <> ";
const double rmtsq = Sqr(" <> mtStr <> ");
const double scalesq = Sqr(get_scale());
const double vev2 = " <> vev2Str <> ";
const double tanb = " <> tanbStr <> ";
const double amu = Re(" <> muStr <> ");
const double mg = " <> m3Str <> ";
const double mAsq = Sqr(" <> mA0Str <> ");
const double cotbeta = 1.0 / tanb;
const double rmbsq = Sqr(" <> mbStr <> ");
const double rmtausq = Sqr(" <> mtauStr <> ");

" <> CConversion`CreateCType[TreeMasses`GetMassMatrixType[SARAH`HiggsBoson]] <> " self_energy_2l(" <> CConversion`CreateCType[TreeMasses`GetMassMatrixType[SARAH`HiggsBoson]] <> "::Zero());

if (HIGGS_2LOOP_CORRECTION_AT_AS) {
   self_energy_2l += self_energy_higgs_2loop_at_as_mssm(
      rmtsq, mg, mst1sq, mst2sq, sxt, cxt, scalesq, amu,
      tanb, vev2, gs);
}

if (HIGGS_2LOOP_CORRECTION_AB_AS) {
   self_energy_2l += self_energy_higgs_2loop_ab_as_mssm(
      rmbsq, mg, msb1sq, msb2sq, sxb, cxb, scalesq, amu,
      cotbeta, vev2, gs);
}

if (HIGGS_2LOOP_CORRECTION_AT_AT) {
   self_energy_2l += self_energy_higgs_2loop_at_at_mssm(
      rmtsq, rmbsq, mAsq, mst1sq, mst2sq, msb1sq, msb2sq,
      sxt, cxt, sxb, cxb, scalesq, amu, tanb, vev2);
}

if (HIGGS_2LOOP_CORRECTION_ATAU_ATAU) {
   self_energy_2l += self_energy_higgs_2loop_atau_atau_mssm(
      rmtausq, mAsq, msnusq, mstau1sq, mstau2sq, sintau,
      costau, scalesq, amu, tanb, vev2);
}

return self_energy_2l;"
          ];

GetNLoopSelfEnergyCorrections[particle_ /; particle === SARAH`PseudoScalar,
                              model_String /; model === "MSSM", 2] :=
    Module[{g3Str, mtStr, mbStr, mtauStr,
            mTop, mBot, mTau,
            vev2Str, vuStr, vdStr, tanbStr, muStr, m3Str, mA0Str},
           AssertFieldDimension[particle, 2, model];
           mTop    = TreeMasses`GetMass[TreeMasses`GetUpQuark[3,True]];
           mBot    = TreeMasses`GetMass[TreeMasses`GetDownQuark[3,True]];
           mTau    = TreeMasses`GetMass[TreeMasses`GetDownLepton[3,True]];
           mtStr   = CConversion`RValueToCFormString[mTop];
           mbStr   = CConversion`RValueToCFormString[mBot];
           mtauStr = CConversion`RValueToCFormString[mTau];
           g3Str   = CConversion`RValueToCFormString[SARAH`strongCoupling];
           vev2Str = CConversion`RValueToCFormString[SARAH`VEVSM1^2 + SARAH`VEVSM2^2];
           vdStr   = CConversion`RValueToCFormString[SARAH`VEVSM1];
           vuStr   = CConversion`RValueToCFormString[SARAH`VEVSM2];
           tanbStr = CConversion`RValueToCFormString[SARAH`VEVSM2 / SARAH`VEVSM1];
           muStr   = CConversion`RValueToCFormString[-Parameters`GetEffectiveMu[]];
           m3Str   = CConversion`RValueToCFormString[FlexibleSUSY`M[SARAH`Gluino]];
           mA0Str  = TreeMasses`CallPseudoscalarHiggsMassGetterFunction[] <> "(0)";
"\
using namespace flexiblesusy::mssm_twoloophiggs;

// calculate 3rd generation sfermion masses and mixing angles
double mst_1, mst_2, theta_t;
double msb_1, msb_2, theta_b;
double mstau_1, mstau_2, theta_tau;
double msnu_1, msnu_2, theta_nu;

" <> TreeMasses`CallGenerationHelperFunctionName[3, SARAH`TopSquark, "mst_1", "mst_2", "theta_t"] <>
";
" <> TreeMasses`CallGenerationHelperFunctionName[3, SARAH`BottomSquark, "msb_1", "msb_2", "theta_b"] <>
";
" <> TreeMasses`CallGenerationHelperFunctionName[3, SARAH`Selectron, "mstau_1", "mstau_2", "theta_tau"] <>
";
" <> TreeMasses`CallGenerationHelperFunctionName[3, SARAH`Sneutrino, "msnu_1", "msnu_2", "theta_nu"] <>
";

const double mst1sq = Sqr(mst_1), mst2sq = Sqr(mst_2);
const double msb1sq = Sqr(msb_1), msb2sq = Sqr(msb_2);
const double mstau1sq = Sqr(mstau_1), mstau2sq = Sqr(mstau_2);
const double msnusq = Sqr(msnu_2);
const double sxt = Sin(theta_t), cxt = Cos(theta_t);
const double sxb = Sin(theta_b), cxb = Cos(theta_b);
const double sintau = Sin(theta_tau), costau = Cos(theta_tau);
const double gs = " <> g3Str <> ";
const double rmtsq = Sqr(" <> mtStr <> ");
const double scalesq = Sqr(get_scale());
const double vev2 = " <> vev2Str <> ";
const double tanb = " <> tanbStr <> ";
const double amu = Re(" <> muStr <> ");
const double mg = " <> m3Str <> ";
const double mAsq = Sqr(" <> mA0Str <> ");
const double cotbeta = 1.0 / tanb;
const double rmbsq = Sqr(" <> mbStr <> ");
const double rmtausq = Sqr(" <> mtauStr <> ");

" <> CConversion`CreateCType[TreeMasses`GetMassMatrixType[SARAH`PseudoScalar]] <> " self_energy_2l(" <> CConversion`CreateCType[TreeMasses`GetMassMatrixType[SARAH`PseudoScalar]] <> "::Zero());

if (HIGGS_2LOOP_CORRECTION_AT_AS) {
   self_energy_2l += self_energy_pseudoscalar_2loop_at_as_mssm(
      rmtsq, mg, mst1sq, mst2sq, sxt, cxt, scalesq, amu,
      tanb, vev2, gs);
}

if (HIGGS_2LOOP_CORRECTION_AB_AS) {
   self_energy_2l += self_energy_pseudoscalar_2loop_ab_as_mssm(
      rmbsq, mg, msb1sq, msb2sq, sxb, cxb, scalesq, amu,
      cotbeta, vev2, gs);
}

if (HIGGS_2LOOP_CORRECTION_AT_AT) {
   self_energy_2l += self_energy_pseudoscalar_2loop_at_at_mssm(
      rmtsq, rmbsq, mAsq, mst1sq, mst2sq, msb1sq, msb2sq,
      sxt, cxt, sxb, cxb, scalesq, amu, tanb, vev2);
}

if (HIGGS_2LOOP_CORRECTION_ATAU_ATAU) {
   self_energy_2l += self_energy_pseudoscalar_2loop_atau_atau_mssm(
      rmtausq, mAsq, msnusq, mstau1sq, mstau2sq, sintau,
      costau, scalesq, amu, tanb, vev2);
}

return self_energy_2l;"
          ];

GetNLoopSelfEnergyCorrections[particle_ /; particle === SARAH`HiggsBoson,
                              model_String /; model === "NMSSM", 2] :=
    Module[{g3Str, mtStr, mbStr, mtauStr,
            mTop, mBot, mTau,
            vev2Str, vuStr, vdStr, vsStr, tanbStr, muStr, m3Str, mA0Str,
            lambdaStr},
           AssertFieldDimension[particle, 3, model];
           mTop    = TreeMasses`GetMass[TreeMasses`GetUpQuark[3,True]];
           mBot    = TreeMasses`GetMass[TreeMasses`GetDownQuark[3,True]];
           mTau    = TreeMasses`GetMass[TreeMasses`GetDownLepton[3,True]];
           mtStr   = CConversion`RValueToCFormString[mTop];
           mbStr   = CConversion`RValueToCFormString[mBot];
           mtauStr = CConversion`RValueToCFormString[mTau];
           g3Str   = CConversion`RValueToCFormString[SARAH`strongCoupling];
           vev2Str = CConversion`RValueToCFormString[SARAH`VEVSM1^2 + SARAH`VEVSM2^2];
           vdStr   = CConversion`RValueToCFormString[SARAH`VEVSM1];
           vuStr   = CConversion`RValueToCFormString[SARAH`VEVSM2];
           tanbStr = CConversion`RValueToCFormString[SARAH`VEVSM2 / SARAH`VEVSM1];
           muStr   = CConversion`RValueToCFormString[-Parameters`GetEffectiveMu[]];
           m3Str   = CConversion`RValueToCFormString[FlexibleSUSY`M[SARAH`Gluino]];
           mA0Str  = CConversion`RValueToCFormString[Parameters`GetEffectiveMASqr[]];
           vsStr   = CConversion`RValueToCFormString[Parameters`GetParameterFromDescription["Singlet-VEV"]];
           lambdaStr = CConversion`RValueToCFormString[Parameters`GetParameterFromDescription["Singlet-Higgs-Interaction"]];
"\
using namespace flexiblesusy::mssm_twoloophiggs;
using namespace flexiblesusy::nmssm_twoloophiggs;

// calculate 3rd generation sfermion masses and mixing angles
double mst_1, mst_2, theta_t;
double msb_1, msb_2, theta_b;
double mstau_1, mstau_2, theta_tau;
double msnu_1, msnu_2, theta_nu;

" <> TreeMasses`CallGenerationHelperFunctionName[3, SARAH`TopSquark, "mst_1", "mst_2", "theta_t"] <>
";
" <> TreeMasses`CallGenerationHelperFunctionName[3, SARAH`BottomSquark, "msb_1", "msb_2", "theta_b"] <>
";
" <> TreeMasses`CallGenerationHelperFunctionName[3, SARAH`Selectron, "mstau_1", "mstau_2", "theta_tau"] <>
";
" <> TreeMasses`CallGenerationHelperFunctionName[3, SARAH`Sneutrino, "msnu_1", "msnu_2", "theta_nu"] <>
";

const double mst1sq = Sqr(mst_1), mst2sq = Sqr(mst_2);
const double msb1sq = Sqr(msb_1), msb2sq = Sqr(msb_2);
const double mstau1sq = Sqr(mstau_1), mstau2sq = Sqr(mstau_2);
const double msnusq = Sqr(msnu_2);
const double sxt = Sin(theta_t), cxt = Cos(theta_t);
const double sxb = Sin(theta_b), cxb = Cos(theta_b);
const double sintau = Sin(theta_tau), costau = Cos(theta_tau);
const double gs = " <> g3Str <> ";
const double as = Sqr(gs) / (4.0 * Pi);
const double rmt = " <> mtStr <> ";
const double rmtsq = Sqr(rmt);
const double scalesq = Sqr(get_scale());
const double vev2 = " <> vev2Str <> ";
const double tanb = " <> tanbStr <> ";
const double amu = Re(" <> muStr <> ");
const double mg = " <> m3Str <> ";
const double mAsq = " <> mA0Str <> ";
const double cotb = 1.0 / tanb;
const double rmb = " <> mbStr <> ";
const double rmbsq = Sqr(rmb);
const double rmtausq = Sqr(" <> mtauStr <> ");
const double lam = Re(" <> lambdaStr <> ");
const double svev = Abs(amu / lam);

" <> CConversion`CreateCType[TreeMasses`GetMassMatrixType[SARAH`HiggsBoson]] <> " self_energy_2l(" <> CConversion`CreateCType[TreeMasses`GetMassMatrixType[SARAH`HiggsBoson]] <> "::Zero());

if (HIGGS_2LOOP_CORRECTION_AT_AS) {
   self_energy_2l += self_energy_higgs_2loop_at_as_nmssm(
      rmt, mg, mst1sq, mst2sq, sxt, cxt,
      scalesq, tanb, vev2, lam, svev, as, amu);
}

if (HIGGS_2LOOP_CORRECTION_AB_AS) {
   self_energy_2l += self_energy_higgs_2loop_ab_as_nmssm(
      rmb, mg, msb1sq, msb2sq, sxb, cxb,
      scalesq, cotb, vev2, lam, svev, as, amu);
}

// Corrections as in MSSM, not corrected for NMSSM,
// should be OK for MSSM states when S state is close to decoupled

if (HIGGS_2LOOP_CORRECTION_AT_AT) {
   self_energy_2l.topLeftCorner<2,2>() += self_energy_higgs_2loop_at_at_mssm(
      rmtsq, rmbsq, mAsq, mst1sq, mst2sq, msb1sq,
      msb2sq, sxt, cxt, sxb, cxb, scalesq, amu, tanb, vev2);
}

if (HIGGS_2LOOP_CORRECTION_ATAU_ATAU) {
   self_energy_2l.topLeftCorner<2,2>() += self_energy_higgs_2loop_atau_atau_mssm(
      rmtausq, mAsq, msnusq, mstau1sq, mstau2sq, sintau,
      costau, scalesq, amu, tanb, vev2);
}

return self_energy_2l;"
          ];

GetNLoopSelfEnergyCorrections[particle_ /; particle === SARAH`PseudoScalar,
                              model_String /; model === "NMSSM", 2] :=
    Module[{g3Str, mtStr, mbStr, mtauStr,
            mTop, mBot, mTau,
            vev2Str, vuStr, vdStr, vsStr, tanbStr, muStr, m3Str, mA0Str,
            lambdaStr},
           AssertFieldDimension[particle, 3, model];
           mTop    = TreeMasses`GetMass[TreeMasses`GetUpQuark[3,True]];
           mBot    = TreeMasses`GetMass[TreeMasses`GetDownQuark[3,True]];
           mTau    = TreeMasses`GetMass[TreeMasses`GetDownLepton[3,True]];
           mtStr   = CConversion`RValueToCFormString[mTop];
           mbStr   = CConversion`RValueToCFormString[mBot];
           mtauStr = CConversion`RValueToCFormString[mTau];
           g3Str   = CConversion`RValueToCFormString[SARAH`strongCoupling];
           vev2Str = CConversion`RValueToCFormString[SARAH`VEVSM1^2 + SARAH`VEVSM2^2];
           vdStr   = CConversion`RValueToCFormString[SARAH`VEVSM1];
           vuStr   = CConversion`RValueToCFormString[SARAH`VEVSM2];
           tanbStr = CConversion`RValueToCFormString[SARAH`VEVSM2 / SARAH`VEVSM1];
           muStr   = CConversion`RValueToCFormString[-Parameters`GetEffectiveMu[]];
           m3Str   = CConversion`RValueToCFormString[FlexibleSUSY`M[SARAH`Gluino]];
           mA0Str  = CConversion`RValueToCFormString[Parameters`GetEffectiveMASqr[]];
           vsStr   = CConversion`RValueToCFormString[Parameters`GetParameterFromDescription["Singlet-VEV"]];
           lambdaStr = CConversion`RValueToCFormString[Parameters`GetParameterFromDescription["Singlet-Higgs-Interaction"]];
"\
using namespace flexiblesusy::mssm_twoloophiggs;
using namespace flexiblesusy::nmssm_twoloophiggs;

// calculate 3rd generation sfermion masses and mixing angles
double mst_1, mst_2, theta_t;
double msb_1, msb_2, theta_b;
double mstau_1, mstau_2, theta_tau;
double msnu_1, msnu_2, theta_nu;

" <> TreeMasses`CallGenerationHelperFunctionName[3, SARAH`TopSquark, "mst_1", "mst_2", "theta_t"] <>
";
" <> TreeMasses`CallGenerationHelperFunctionName[3, SARAH`BottomSquark, "msb_1", "msb_2", "theta_b"] <>
";
" <> TreeMasses`CallGenerationHelperFunctionName[3, SARAH`Selectron, "mstau_1", "mstau_2", "theta_tau"] <>
";
" <> TreeMasses`CallGenerationHelperFunctionName[3, SARAH`Sneutrino, "msnu_1", "msnu_2", "theta_nu"] <>
";

const double mst1sq = Sqr(mst_1), mst2sq = Sqr(mst_2);
const double msb1sq = Sqr(msb_1), msb2sq = Sqr(msb_2);
const double mstau1sq = Sqr(mstau_1), mstau2sq = Sqr(mstau_2);
const double msnusq = Sqr(msnu_2);
const double sxt = Sin(theta_t), cxt = Cos(theta_t);
const double sxb = Sin(theta_b), cxb = Cos(theta_b);
const double sintau = Sin(theta_tau), costau = Cos(theta_tau);
const double gs = " <> g3Str <> ";
const double as = Sqr(gs) / (4.0 * Pi);
const double rmt = " <> mtStr <> ";
const double rmtsq = Sqr(rmt);
const double scalesq = Sqr(get_scale());
const double vev2 = " <> vev2Str <> ";
const double tanb = " <> tanbStr <> ";
const double amu = Re(" <> muStr <> ");
const double mg = " <> m3Str <> ";
const double mAsq = " <> mA0Str <> ";
const double cotb = 1.0 / tanb;
const double rmb = " <> mbStr <> ";
const double rmbsq = Sqr(rmb);
const double rmtausq = Sqr(" <> mtauStr <> ");
const double lam = Re(" <> lambdaStr <> ");
const double svev = Abs(amu / lam);

" <> CConversion`CreateCType[TreeMasses`GetMassMatrixType[SARAH`PseudoScalar]] <> " self_energy_2l(" <> CConversion`CreateCType[TreeMasses`GetMassMatrixType[SARAH`PseudoScalar]] <> "::Zero());

if (HIGGS_2LOOP_CORRECTION_AT_AS) {
   self_energy_2l += self_energy_pseudoscalar_2loop_at_as_nmssm(
      rmt, mg, mst1sq, mst2sq, sxt, cxt,
      scalesq, tanb, vev2, lam, svev, as, amu);
}

if (HIGGS_2LOOP_CORRECTION_AB_AS) {
   self_energy_2l += self_energy_pseudoscalar_2loop_ab_as_nmssm(
      rmb, mg, msb1sq, msb2sq, sxb, cxb,
      scalesq, cotb, vev2, lam, svev, as, amu);
}

// Corrections as in MSSM, not corrected for NMSSM,
// should be OK for MSSM states when S state is close to decoupled

if (HIGGS_2LOOP_CORRECTION_AT_AT) {
   self_energy_2l.topLeftCorner<2,2>() += self_energy_pseudoscalar_2loop_at_at_mssm(
      rmtsq, rmbsq, mAsq, mst1sq, mst2sq, msb1sq, msb2sq,
      sxt, cxt, sxb, cxb, scalesq, amu, tanb, vev2);
}

if (HIGGS_2LOOP_CORRECTION_ATAU_ATAU) {
   self_energy_2l.topLeftCorner<2,2>() += self_energy_pseudoscalar_2loop_atau_atau_mssm(
      rmtausq, mAsq, msnusq, mstau1sq, mstau2sq, sintau,
      costau, scalesq, amu, tanb, vev2);
}

return self_energy_2l;"
          ];

GetNLoopSelfEnergyCorrections[particle_ /; particle === SARAH`HiggsBoson,
                              model_String /; model === "MSSM", 3] :=
    Module[{g3Str, mtStr, mbStr, meStr, mTop, mBot, mTau,
            vuStr, vdStr, muStr, m3Str, mA0Str,
            AtStr, AbStr, AeStr, mWStr, mZStr,
            mq2Str, md2Str, mu2Str, ml2Str, me2Str},
           AssertFieldDimension[particle, 2, model];
           mTop    = TreeMasses`GetMass[TreeMasses`GetUpQuark[3,True]];
           mBot    = TreeMasses`GetMass[TreeMasses`GetDownQuark[3,True]];
           mTau    = TreeMasses`GetMass[TreeMasses`GetDownLepton[3,True]];
           mtStr   = CConversion`RValueToCFormString[mTop];
           mbStr   = CConversion`RValueToCFormString[mBot];
           meStr   = CConversion`RValueToCFormString[mTau];
           g3Str   = CConversion`RValueToCFormString[SARAH`strongCoupling];
           vdStr   = CConversion`RValueToCFormString[SARAH`VEVSM1];
           vuStr   = CConversion`RValueToCFormString[SARAH`VEVSM2];
           muStr   = CConversion`RValueToCFormString[Parameters`GetEffectiveMu[]];
           m3Str   = CConversion`RValueToCFormString[FlexibleSUSY`M[SARAH`Gluino]];
           mA0Str  = TreeMasses`CallPseudoscalarHiggsMassGetterFunction[] <> "(0)";
           AtStr   = CConversion`RValueToCFormString[SARAH`TrilinearUp[2,2] / SARAH`UpYukawa[2,2]];
           AbStr   = CConversion`RValueToCFormString[SARAH`TrilinearDown[2,2] / SARAH`DownYukawa[2,2]];
           AeStr   = CConversion`RValueToCFormString[SARAH`TrilinearLepton[2,2] / SARAH`ElectronYukawa[2,2]];
           mWStr   = CConversion`RValueToCFormString[FlexibleSUSY`M[SARAH`VectorW]];
           mZStr   = CConversion`RValueToCFormString[FlexibleSUSY`M[SARAH`VectorZ]];
           mq2Str  = CConversion`RValueToCFormString[SARAH`SoftSquark];
           mu2Str  = CConversion`RValueToCFormString[SARAH`SoftUp];
           md2Str  = CConversion`RValueToCFormString[SARAH`SoftDown];
           ml2Str  = CConversion`RValueToCFormString[SARAH`SoftLeftLepton];
           me2Str  = CConversion`RValueToCFormString[SARAH`SoftRightLepton];
CConversion`CreateCType[TreeMasses`GetMassMatrixType[SARAH`HiggsBoson]] <> " self_energy_3l(" <> CConversion`CreateCType[TreeMasses`GetMassMatrixType[SARAH`HiggsBoson]] <> "::Zero());

#ifdef ENABLE_HIMALAYA
// calculate 3rd generation sfermion masses and mixing angles
double mst_1, mst_2, theta_t;
double msb_1, msb_2, theta_b;

" <> TreeMasses`CallGenerationHelperFunctionName[3, SARAH`TopSquark, "mst_1", "mst_2", "theta_t"] <>
";
" <> TreeMasses`CallGenerationHelperFunctionName[3, SARAH`BottomSquark, "msb_1", "msb_2", "theta_b"] <>
";

himalaya::Parameters pars;
pars.scale = get_scale();
pars.mu = Re(" <> muStr <> ");
pars.g3 = " <> g3Str <> ";
pars.vd = " <> vdStr <> ";
pars.vu = " <> vuStr <> ";
pars.mq2 = Re(" <> mq2Str <> ");
pars.md2 = Re(" <> md2Str <> ");
pars.mu2 = Re(" <> mu2Str <> ");
pars.MG = " <> m3Str <> ";
pars.MW = " <> mWStr <> ";
pars.MZ = " <> mZStr <> ";
pars.Mt = " <> mtStr <> ";
pars.Mb = " <> mbStr <> ";
pars.MA = " <> mA0Str <> ";
pars.MSt << mst_1, mst_2;
pars.MSb << msb_1, msb_2;
pars.s2t = Sin(2*theta_t);
pars.s2b = Sin(2*theta_b);

#if Himalaya_VERSION_MAJOR < 2
   pars.At = Re(" <> AtStr <> ");
   pars.Ab = Re(" <> AbStr <> ");
#else
   pars.Au(2,2) = Re(" <> AtStr <> ");
   pars.Ad(2,2) = Re(" <> AbStr <> ");
   pars.Ae(2,2) = Re(" <> AeStr <> ");
   pars.ml2 = Re(" <> ml2Str <> ");
   pars.me2 = Re(" <> me2Str <> ");
   pars.Mtau = " <> meStr <> ";
#endif

try {
   const auto ren_scheme = HIGGS_3LOOP_SCHEME;
   const bool verbose = false;
   himalaya::HierarchyCalculator hc(pars, verbose);

   if (HIGGS_3LOOP_CORRECTION_AT_AS_AS) {
#if Himalaya_VERSION_MAJOR < 2
      const auto hier = hc.calculateDMh3L(false, ren_scheme);
#else
      const auto hier = hc.calculateDMh3L(false);
#endif

      // calculate the 3-loop corrections
      self_energy_3l += - hier.getDMh(3);

#if Himalaya_VERSION_MAJOR < 2
      if (ren_scheme) {
         // calculate shift DR -> MDR
         self_energy_3l += - hier.getDRToMDRShift();
      }
#else
      if (ren_scheme == 1) {
         // calculate shift DR' -> MDR'
         self_energy_3l += - hier.getDMhDRbarPrimeToMDRbarPrimeShift();
      } else if (ren_scheme == 1) {
         // calculate shift DR' -> H3m
         self_energy_3l += - hier.getDMhDRbarPrimeToH3mShift();
      }
#endif
   }

   if (HIGGS_3LOOP_CORRECTION_AB_AS_AS) {
#if Himalaya_VERSION_MAJOR < 2
      const auto hier = hc.calculateDMh3L(true, ren_scheme);
#else
      const auto hier = hc.calculateDMh3L(true);
#endif

      // calculate the 3-loop corrections
      self_energy_3l += - hier.getDMh(3);

#if Himalaya_VERSION_MAJOR < 2
      if (ren_scheme) {
         // calculate the shift DR -> MDR
         self_energy_3l += - hier.getDRToMDRShift();
      }
#else
      if (ren_scheme == 1) {
         // calculate shift DR' -> MDR'
         self_energy_3l += - hier.getDMhDRbarPrimeToMDRbarPrimeShift();
      } else if (ren_scheme == 1) {
         // calculate shift DR' -> H3m
         self_energy_3l += - hier.getDMhDRbarPrimeToH3mShift();
      }
#endif
   }
} catch (const std::exception& e) {
   VERBOSE_MSG(e.what());
   VERBOSE_MSG(pars);
   throw HimalayaError(e.what());
}
#else // ENABLE_HIMALAYA
throw HimalayaError(\"The 3-loop corrections to Mh require Himalaya 1.0 \"
                    \"(or higher), but FlexibleSUSY has not been \"
                    \"configured with Himalaya!\");
#endif // ENABLE_HIMALAYA

return self_energy_3l;"
          ];

GetNLoopSelfEnergyCorrections[particle_, model_, loop_] :=
    Module[{},
           Print["Error: ", loop,"-loop self-energy corrections not available for ",
                 particle, " in the ", model];
           ""
          ];

CreateNLoopSelfEnergy[particle_, model_String, loop_, args_String] :=
    Module[{prototype, function, functionName, dim, dimStr, cType},
           dim = Parameters`NumberOfIndependentEntriesOfSymmetricMatrix[GetDimension[particle]];
           dimStr = ToString[dim];
           functionName = CreateSelfEnergyFunctionName[particle,loop];
           cType = CConversion`CreateCType[TreeMasses`GetMassMatrixType[particle]];
           prototype = cType <> " " <> functionName <> "(" <> args <> ") const;\n";
           body = GetNLoopSelfEnergyCorrections[particle, model, loop];
           function = cType <> " CLASSNAME::" <> functionName <> "(" <> args <> ") const\n{\n" <>
                      IndentText[body] <> "\n}\n";
           Return[{prototype, function}];
          ];

CreateNLoopSelfEnergies[particles_List, model_String, loop_, args_String:""] :=
    Module[{prototype = "", function = "", i, p, f},
           For[i = 1, i <= Length[particles], i++,
               {p, f} = CreateNLoopSelfEnergy[particles[[i]], model, loop, args];
               prototype = prototype <> p;
               function = function <> f <> "\n";
              ];
           Return[{prototype, function}];
          ];

CreateTwoLoopSelfEnergiesSM[particles_List] :=
    CreateNLoopSelfEnergies[particles, "SM", 2, "double p"];

CreateThreeLoopSelfEnergiesSM[particles_List] :=
    CreateNLoopSelfEnergies[particles, "SM", 3];

CreateFourLoopSelfEnergiesSM[particles_List] :=
    CreateNLoopSelfEnergies[particles, "SM", 4];

CreateTwoLoopSelfEnergiesMSSM[particles_List] :=
    CreateNLoopSelfEnergies[particles, "MSSM", 2];

CreateTwoLoopSelfEnergiesNMSSM[particles_List] :=
    CreateNLoopSelfEnergies[particles, "NMSSM", 2];

CreateThreeLoopSelfEnergiesMSSM[particles_List] :=
    CreateNLoopSelfEnergies[particles, "MSSM", 3];

CreateThreeLoopSelfEnergiesSplit[particles_List] :=
    CreateNLoopSelfEnergies[particles, "Split", 3];

End[];

EndPackage[];
