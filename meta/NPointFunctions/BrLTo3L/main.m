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

BeginPackage@"BrLTo3L`";Quiet[

BrLTo3L`create::usage = "
@brief An interface function to make C++ code for an observable.
@param obs A sigle observable to calculate.
@param list A list of observables to calculate.
@returns A set of required vertices, headers and definitions for npf functions,
         prototypes and definitions for calculate function.";

];Begin@"`Private`";

$boxes = "boxes";
$boxes // Protect;

setGlobals[obs:`type`observable] :=
Module[{cxx = CConversion`ToValidCSymbolString},
   Unprotect@"BrLTo3L`Private`$*";
   $fields = Utils`StringJoinWithSeparator["fields::"<>#&/@ cxx/@
      {lep, TreeMasses`GetPhoton[]}, ", "];
   $calculate = StringReplace[calculate[obs, Head], s_~~"("~~__ :> s];
   $penguin = StringJoin["calculate_", cxx@lep, "_", cxx@lep, "_",
      cxx@TreeMasses`GetPhoton[], "_form_factors"];
   $prototype = CConversion`CreateCType@Observables`GetObservableType@obs <>
      " " <> calculate[obs, Head];
   Protect@"BrLTo3L`Private`$*";];
setGlobals // Utils`MakeUnknownInputDefinition;
setGlobals // Protect;

create::usage = "
@brief Creates ingredients for the C++ code, including all required vertices.
       Boxes and penguins are treated differently: boxes generate all
       contributions and penguins will be taken for t-contribution only, because
       of the simple and predictable structure of the result. It means, that
       we can reduce the computation time by separating those contributions and
       using one general function to calculate boxes.
@param obs An observable definition.
@param list A list of all requested observables.
@returns A list of vertices, npf headers string, npf definition string,
         calculate function prototypes, calculate function definitions.
@note We assume that generations for leptons are given by integer numbers,
      rather then different field names.";
create[list:{__}] :=
Module[{unique, box, join = Utils`StringJoinWithSeparator},
   unique = DeleteDuplicates[list, SameQ@@({##} /. _Integer :> Sequence)&];
   box = `npf`box@unique;
   {  DeleteDuplicates[Join@@Append[#[[All, 1]], box[[1, 1]]]],
      box[[2]],
      join[Append[#[[All, 2]], box[[1, 2]]], "\n\n"],
      join[#[[All, 3]], "\n\n"],
      join[#[[All, 4]], "\n\n"]}&[create/@unique]];
create[obs:`type`observable] :=
Module[{pengVert = {}, pengDef = "", boxQ, calculateDefinition},
   setGlobals@obs;
   {{pengVert, pengDef}, boxQ} = `npf`assemble@obs;
   calculateDefinition = $prototype <> " {
   return forge<
      "<>$fields<>",
      "<>$penguin<>",
      npointfunctions::"<>If[pengDef =!= "", $calculate, "zero"]<>",
      npointfunctions::"<>If[boxQ, $boxes, "zero"]<>"
   >(nI, nO, nA, model, qedqcd);\n}";
   {  pengVert,
      pengDef,
      $prototype <> ";",
      calculateDefinition}];
create // Utils`MakeUnknownInputDefinition;
create // Protect;

`npf`box[list:{__}] :=
Module[{box = Boxes, res = {{}, ""}},
   If[!FreeQ[`npf`parse/@list, box],
      Utils`FSFancyLine@"<";
      Print["Calculation for "<>SymbolName@box<>" started"];
      Unprotect@$calculate;
      $calculate = $boxes;
      Protect@$calculate;
      res = `npf`code@`npf`match@`npf`clean@`npf`create[list[[1]], {box}];
      Utils`FSFancyLine@">";];
   {res, NPointFunctions`CreateCXXHeaders[]}];
`npf`box // Utils`MakeUnknownInputDefinition;
`npf`box // Protect;

`npf`parse@`type`observable :=
Module[{parsed},
   parsed = SymbolName/@If[Head@# === List, #, {#}]&@proc;
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

`npf`assemble[obs:`type`observable] :=
Module[{keep, peng, out, boxQ},
   boxQ = True;
   keep = `npf`parse@obs;
   peng = Complement[keep, {Boxes}];
   Utils`FSFancyLine@"<";
   Print[      "Calculation for "<>Utils`StringJoinWithSeparator[
      peng, ",\n                ", SymbolName]<>" started"];
   Switch[peng,
      (* no boxes *) keep,
         out = `npf`code@`npf`match@`npf`clean@`npf`create[obs, peng];
         boxQ = False;,
      (* only boxes *) {},
         Print["Boxes exist already! Skipping."];
         out = {{}, ""};,
      (* both *) _,
         out = `npf`code@`npf`match@`npf`clean@`npf`create[obs, peng];];
   Utils`FSFancyLine@">";
   {out, boxQ}];
`npf`assemble // Utils`MakeUnknownInputDefinition;
`npf`assemble // Protect;

`npf`create[obs:`type`observable, keep:{__}] :=
NPointFunctions`NPointFunction[{lep, lep}, {lep, lep},
   NPointFunctions`OnShellFlag -> True,
   NPointFunctions`UseCache -> False,
   NPointFunctions`ZeroExternalMomenta -> NPointFunctions`ExceptLoops,
   NPointFunctions`KeepProcesses -> keep,
   NPointFunctions`Observable -> obs];
`npf`create // Utils`MakeUnknownInputDefinition;
`npf`create // Protect;

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

`npf`match::usage = "
@note String names on the lhs. are representing the final names of coefficients
      for the C++ code, after applying some relations on the C++ level! Check
      appropriate file in templates directory to see them.";
`npf`match[npf:NPointFunctions`internal`type`npf] :=
Module[{fields, sp, dc, dim6},
   fields = Flatten@NPointFunctions`internal`getProcess@npf;
   sp[i_] := SARAH`DiracSpinor[fields[[i]], 0, 0];
   dc[a_, b__, c_] := NPointFunctions`internal`dc[sp@a, b, sp@c];
   dim6 = With[{l = SARAH`Lorentz, R = 6, L = 7},
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
   {WilsonCoeffs`InterfaceToMatching[npf, dim6], dim6}];
`npf`match // Utils`MakeUnknownInputDefinition;
`npf`match // Protect;

`npf`code[{npf:NPointFunctions`internal`type`npf, b:{Rule[_String, _]..}}] :=
{  DeleteDuplicates@NPointFunctions`VerticesForNPointFunction@npf,
   NPointFunctions`CreateCXXFunctions[npf, $calculate, Identity, b][[2]]};
`npf`code // Utils`MakeUnknownInputDefinition;
`npf`code // Protect;

End[];
EndPackage[];
$ContextPath = DeleteCases[$ContextPath, "BrLTo3L`"];
Unprotect@$Packages;
$Packages = DeleteCases[$Packages, "BrLTo3L`"];
Protect@$Packages;
