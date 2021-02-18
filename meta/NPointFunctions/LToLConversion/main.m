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

BeginPackage@"LToLConversion`";Quiet[

LToLConversion`create::usage =
"@brief Main entrance point for the calculation.";

];Begin@"`Private`";

setCxx[obs:`type`observable] := Module[{cxx = CConversion`ToValidCSymbolString},
   Unprotect@"LToLConversion`Private`cxx`*";
   `cxx`in = cxx@in;
   `cxx`out = cxx@out;

   `cxx`fields = StringJoin@@Riffle["fields::"<>#&/@ cxx/@
      {in, SARAH`UpQuark, SARAH`DownQuark, SARAH`Photon}, ", "];

   {`cxx`classU, `cxx`classD} = StringJoin["conversion_", cxx@in, #, "_to_",
      cxx@out, #, "_for_", SymbolName@con]&/@ cxx/@
         {SARAH`UpQuark, SARAH`DownQuark};

   `cxx`penguin = StringJoin["calculate_", cxx@in, "_", cxx@out, "_",
      cxx@SARAH`Photon, "_form_factors"];

   (*TODO this code is partially duplicated in CalculateObservable.*)
   `cxx`prototype = CConversion`CreateCType@Observables`GetObservableType@obs <>
      " calculate_"<>cxx@in<>"_to_"<>cxx@out<>"_for_"<>SymbolName@con<>"(\n"<>
      "   int in, int out,\n"<>
      "   const " <> FlexibleSUSY`FSModelName <>
         "_l_to_l_conversion::Nucleus nucleus,\n" <>
      "   const " <> FlexibleSUSY`FSModelName <>
      "_mass_eigenstates& model, const softsusy::QedQcd& qedqcd)";
   Protect@"LToLConversion`Private`cxx`*";
];
setCxx // Utils`MakeUnknownInputDefinition;
setCxx ~ SetAttributes ~ {Protected, Locked};

create[obs:`type`observable] :=
Module[{npfVertices, npfHeader, npfDefinition, calculateDefinition},
   setCxx@obs;
   {npfVertices, npfHeader, npfDefinition} = `npf`create@obs;

   calculateDefinition = `cxx`prototype <> " {
   return forge_conversion<
      "<>`cxx`fields<>",
      "<>`cxx`penguin<>",
      npointfunctions::"<>`cxx`classU<>",
      npointfunctions::"<>`cxx`classD<>"
   >(in, out, nucleus, model, qedqcd);\n}";

   {
      npfVertices,
      {npfHeader, npfDefinition},
      {`cxx`prototype <> ";", calculateDefinition}}];

create[list:{__}] := Module[{unique},
   unique = DeleteDuplicates[list, SameQ@@({##} /. _Integer :> Sequence)&];
   {  DeleteDuplicates[Join@@#[[All,1]]],
      {  #[[1,2,1]], StringJoin@Riffle[#[[All,2,2]], "\n\n"]},
      {  StringJoin@Riffle[#[[All,3,1]], "\n\n"],
         StringJoin@Riffle[#[[All,3,2]], "\n\n"]}}&[create/@unique]];

create // Utils`MakeUnknownInputDefinition;
create ~ SetAttributes ~ {Protected, Locked};

`npf`clean[npf:NPointFunctions`internal`type`npf] :=
npf /.
   {  SARAH`sum[__] -> 0,
      LoopTools`B0i[i_, _, mm__] :> LoopTools`B0i[i, 0, mm],
      LoopTools`C0i[i_, Repeated[_, {3}], mm__] :>
         LoopTools`C0i[i, Sequence@@Array[0&, 3], mm],
      LoopTools`D0i[i_, Repeated[_, {6}], mm__] :>
         LoopTools`D0i[i, Sequence@@Array[0&, 6], mm]};
`npf`clean // Utils`MakeUnknownInputDefinition;
`npf`clean // Protect;

`npf`parse[obs:`type`observable] :=
Module[{parsed},
   parsed = SymbolName/@If[Head@# === List, #, {#}]&@con;
   Switch[parsed,
      {"All"},
         {Vectors, Scalars, Boxes},
      {"NoScalars"},
         {Vectors, Boxes},
      {"Penguins"},
         {Vectors, Scalars},
      _,
         Symbol/@parsed]];
`npf`parse // Utils`MakeUnknownInputDefinition;
`npf`parse // Protect;

`npf`create[obs:`type`observable] :=
Module[{npfU, npfD, fields, keep, dim6, codeU, codeD},
   keep = `npf`parse@obs;
   Utils`FSFancyLine@"<";
   Print[      "Calculation for "<>Utils`StringJoinWithSeparator[
      keep, ",\n                ", SymbolName]<>" started"];
   {npfU, npfD} = NPointFunctions`NPointFunction[
      {in,#},{out,#},
      NPointFunctions`OnShellFlag -> True,
      NPointFunctions`UseCache -> False,
      NPointFunctions`ZeroExternalMomenta -> NPointFunctions`ExceptLoops,
      NPointFunctions`KeepProcesses -> keep,
      NPointFunctions`Observable -> obs] &/@ {SARAH`UpQuark, SARAH`DownQuark};
   {npfU, npfD} = `npf`clean/@{npfU, npfD};

   fields[SARAH`UpQuark] = Flatten@NPointFunctions`internal`getProcess@npfU;
   fields[SARAH`DownQuark] = Flatten@NPointFunctions`internal`getProcess@npfD;

   dim6[q_] := Module[{sp, dc, l = SARAH`Lorentz, R = 6, L = 7},
      sp[f_, n_] := SARAH`DiracSpinor[fields[f][[n]], 0, 0];
      dc[a_, b__, c_] := NPointFunctions`internal`dc[sp[q, a], b, sp[q, c]];
      {  "S_LL" -> dc[3,L,1] dc[4,L,2],
         "S_LR" -> dc[3,L,1] dc[4,R,2],
         "S_RL" -> dc[3,R,1] dc[4,L,2],
         "S_RR" -> dc[3,R,1] dc[4,R,2],
         "V_LL" -> dc[3,R,l@1,1] dc[4,R,l@1,2],
         "V_LR" -> dc[3,R,l@1,1] dc[4,L,l@1,2],
         "V_RL" -> dc[3,L,l@1,1] dc[4,R,l@1,2],
         "V_RR" -> dc[3,L,l@1,1] dc[4,L,l@1,2],
         "T_LL" -> dc[3,-L,l@1,l@2,1] dc[4,-L,l@1,l@2,2],
         "T_RR" -> dc[3,-R,l@1,l@2,1] dc[4,-R,l@1,l@2,2]}];

   npfU = WilsonCoeffs`InterfaceToMatching[npfU, dim6@SARAH`UpQuark];
   npfD = WilsonCoeffs`InterfaceToMatching[npfD, dim6@SARAH`DownQuark];

   codeU = NPointFunctions`CreateCXXFunctions[
      npfU, `cxx`classU, SARAH`Delta, dim6@SARAH`UpQuark][[2]];
   codeD = NPointFunctions`CreateCXXFunctions[
      npfD, `cxx`classD, SARAH`Delta, dim6@SARAH`DownQuark][[2]];
   Utils`FSFancyLine@">";

   {  DeleteDuplicates@Join[
         NPointFunctions`VerticesForNPointFunction@npfU,
         NPointFunctions`VerticesForNPointFunction@npfD],
      NPointFunctions`CreateCXXHeaders[],
      codeU<>"\n\n"<>codeD}];
`npf`create // Utils`MakeUnknownInputDefinition;
`npf`create ~ SetAttributes ~ {Locked,Protected};

End[];
EndPackage[];
$ContextPath = DeleteCases[$ContextPath, "LToLConversion`"];
Unprotect@$Packages;
$Packages = DeleteCases[$Packages, "LToLConversion`"];
Protect@$Packages;
