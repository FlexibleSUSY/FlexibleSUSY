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

setCxx[obs:`type`observable] := Module[{cxx = CConversion`ToValidCSymbolString},
   Unprotect@"BrLTo3L`Private`cxx`*";
   `cxx`lep = cxx@lep;
   `cxx`fields = StringJoin@@Riffle["fields::"<>#&/@ cxx/@
      {lep, SARAH`Photon}, ", "];
   `cxx`name = StringReplace["npf_"<>calculate[obs, Head], s_~~"("~~__ :> s];
   `cxx`penguin = StringJoin["calculate_", cxx@lep, "_", cxx@lep, "_",
      cxx@SARAH`Photon, "_form_factors"];
   `cxx`proto = CConversion`CreateCType@Observables`GetObservableType@obs <>
      " " <> calculate[obs, Head];
   Protect@"BrLTo3L`Private`cxx`*";];
setCxx // Utils`MakeUnknownInputDefinition;
setCxx ~ SetAttributes ~ {Protected, Locked};

create[obs:`type`observable] :=
Module[{npfVertices = {}, npfHeader = "", npfDefinition = "",
      calculateDefinition},
   setCxx@obs;
   {npfVertices, npfHeader, npfDefinition} = `npf`create@obs;

   calculateDefinition = `cxx`proto <> " {
   return forge(nI, nO, nA, model, qedqcd);\n}";

   {  npfVertices,
      {npfHeader, npfDefinition},
      {`cxx`proto <> ";", calculateDefinition}}];

create[list:{__}] :=
   {  DeleteDuplicates[Join@@#[[All,1]]],
      {  #[[1,2,1]], StringJoin@Riffle[#[[All,2,2]], "\n\n"]},
      {  StringJoin@Riffle[#[[All,3,1]], "\n\n"],
         StringJoin@Riffle[#[[All,3,2]], "\n\n"]}}&[create/@list];

create // Utils`MakeUnknownInputDefinition;
create ~ SetAttributes ~ {Protected, Locked};

`npf`parse@`type`observable := Switch[proc,
   All, {
      NPointFunctions`FourFermionMassiveVectorPenguins,
      NPointFunctions`FourFermionScalarPenguins,
      NPointFunctions`FourFermionFlavourChangingBoxes},
   NPointFunctions`noScalars, {
      NPointFunctions`FourFermionMassiveVectorPenguins,
      NPointFunctions`FourFermionFlavourChangingBoxes},
   NPointFunctions`Penguins, {
      NPointFunctions`FourFermionScalarPenguins,
      NPointFunctions`FourFermionMassiveVectorPenguins},
   _List, proc,
   _, {proc}];
`npf`parse // Utils`MakeUnknownInputDefinition;
`npf`parse ~ SetAttributes ~ {Protected, Locked};

`npf`create[obs:`type`observable] := Module[{npf, fields, sp, dc, dim6, code},
   Utils`FSFancyLine@"<";
   Print@StringReplace["Calculation of #-#+ to #-#+ started", "#"->`cxx`lep];
   npf = NPointFunctions`NPointFunction[
      {lep, SARAH`bar@lep}, {lep, SARAH`bar@lep},
      NPointFunctions`OnShellFlag -> True,
      NPointFunctions`UseCache -> False,
      NPointFunctions`ZeroExternalMomenta -> NPointFunctions`ExceptLoops,
      NPointFunctions`KeepProcesses -> `npf`parse@obs,
      NPointFunctions`Observable -> obs];
   npf = npf /. SARAH`sum[__] -> 0;
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
         "minusT_LL" -> dc[3,-L,l@1,l@2,1] dc[4,-L,l@1,l@2,2],
         "minusT_RR" -> dc[3,-R,l@1,l@2,1] dc[4,-R,l@1,l@2,2]}];
   npf = WilsonCoeffs`InterfaceToMatching[npf, dim6];
   code = NPointFunctions`CreateCXXFunctions[
      npf, `cxx`name, Identity, dim6][[2]];
   Print@StringReplace["Calculation of #-#+ to #-#+ finished", "#"->`cxx`lep];
   Utils`FSFancyLine@">";
   {  DeleteDuplicates@NPointFunctions`VerticesForNPointFunction@npf,
      NPointFunctions`CreateCXXHeaders[],
      code}];
`npf`create // Utils`MakeUnknownInputDefinition;
`npf`create ~ SetAttributes ~ {Locked,Protected};

End[];
EndPackage[];
$ContextPath = DeleteCases[$ContextPath, "BrLTo3L`"];
Unprotect@$Packages;
$Packages = DeleteCases[$Packages, "BrLTo3L`"];
Protect@$Packages;
