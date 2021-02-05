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
      cxx@out, #, "_for_", cxx@con, cxx@massless]&/@ cxx/@
         {SARAH`UpQuark, SARAH`DownQuark};

   `cxx`penguin = StringJoin["calculate_", cxx@in, "_", cxx@out, "_",
      cxx@SARAH`Photon, "_form_factors"];

   (*TODO this code is partially duplicated in CalculateObservable.*)
   `cxx`prototype = CConversion`CreateCType@Observables`GetObservableType@obs <>
      " calculate_"<>cxx@in<>"_to_"<>cxx@out<>"_for_"<>cxx@con<>cxx@massless<>"(\n"<>
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

create[list:{__}] := {
   DeleteDuplicates[ Join@@#[[All,1]] ],
   {
      #[[1,2,1]],
      StringJoin@Riffle[#[[All,2,2]], "\n\n"]},
   {
      StringJoin@Riffle[#[[All,3,1]], "\n\n"],
      StringJoin@Riffle[#[[All,3,2]], "\n\n"]}}& [create /@ list];
create // Utils`MakeUnknownInputDefinition;
create ~ SetAttributes ~ {Protected, Locked};

`npf`create[`type`observable] := Module[{
      parsedCon,
      npfU, npfD,
      l=SARAH`Lorentz, p=SARAH`Mom, m=SARAH`Mass,
      dc = NPointFunctions`internal`dc,
      fiG, foG, uiG, uoG, diG, doG, (*@note particle | inc/out | generation*)
      regulator, (*@note arbitrary sqr(3-momenta) of quarks*)
      inner, sp, dim6,
      codeU, codeD},

   parsedCon = Switch[con,
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
      _, con];

   Print["<<npf<< calculation for ",`cxx`in," to ",`cxx`out," conversion started ..."];

   {npfU, npfD} = NPointFunctions`NPointFunction[
      {in,#},{out,#},
      NPointFunctions`OnShellFlag -> True,
      NPointFunctions`UseCache -> False,
      NPointFunctions`ZeroExternalMomenta -> If[massless===True, NPointFunctions`ExceptLoops, NPointFunctions`OperatorsOnly],
      NPointFunctions`KeepProcesses -> parsedCon] &/@ {SARAH`UpQuark,SARAH`DownQuark};

   {fiG, uiG, foG, uoG} = Flatten@NPointFunctions`internal`getProcess@npfU;
   {fiG, diG, foG, doG} = Flatten@NPointFunctions`internal`getProcess@npfD;
   regulator = m@fiG^2;

   inner = SARAH`sum[i_,1,4,SARAH`g[i_,i_]*p[#1,i_]*p[#2,i_]]&;
   {npfU, npfD} = {npfU, npfD} //. {
         inner[fiG,foG] :> m@fiG^2,
         inner[uiG,fiG] :> m@fiG*Sqrt[m@uiG^2+regulator],
         inner[uoG,fiG] :> inner[uiG,fiG],
         inner[uoG,foG] :> m@fiG^2/2+inner[uiG,fiG],
         inner[diG,fiG] :> m@fiG*Sqrt[m@diG^2+regulator],
         inner[doG,fiG] :> inner[diG,fiG],
         inner[doG,foG] :> m@fiG^2/2+inner[diG,fiG]
      };

   (* @note If ZeroExternalMomenta is set to True, replace p and m to zeroes *)
   sp[particle:_,num:_Integer] := If[massless===True, SARAH`DiracSpinor[#, 0, 0], SARAH`DiracSpinor[#, p@num, m@#]] &@
      particle@{Symbol["SARAH`gt"<>ToString@num]};

   dim6[i_,o_,q_] := {
      (*@note 6 means PR, 7 means PL.*)
      "S_LL" -> dc[o~sp~3,7,i~sp~1] dc[q~sp~4,7,q~sp~2],
      "S_LR" -> dc[o~sp~3,7,i~sp~1] dc[q~sp~4,6,q~sp~2],
      "S_RL" -> dc[o~sp~3,6,i~sp~1] dc[q~sp~4,7,q~sp~2],
      "S_RR" -> dc[o~sp~3,6,i~sp~1] dc[q~sp~4,6,q~sp~2],
      (*@note names are correct, one just need to commute projectors with
       *Dirac matrices. It changes 6 to 7 or 7 to 6.*)
      "V_LL" -> dc[o~sp~3,6,l@1,i~sp~1] dc[q~sp~4,6,l@1,q~sp~2],
      "V_LR" -> dc[o~sp~3,6,l@1,i~sp~1] dc[q~sp~4,7,l@1,q~sp~2],
      "V_RL" -> dc[o~sp~3,7,l@1,i~sp~1] dc[q~sp~4,6,l@1,q~sp~2],
      "V_RR" -> dc[o~sp~3,7,l@1,i~sp~1] dc[q~sp~4,7,l@1,q~sp~2],
      (*@note Minus, because FormCalc`s -6,Lor[1],Lor[2] is ours
       *-I*sigma[1,2] (according to FC definition of antisymmetrization), when
       *taking this twice we get I*I=-1. FC cites [Ni05] for Fierz identities,
       *where our conventions are used, but in FC manual on the page 20
       *weird convention for sigma_munu is shown.*)
      "minus_T_LL" -> dc[o~sp~3,-7,l@1,l@2,i~sp~1] dc[q~sp~4,-7,l@1,l@2,q~sp~2],
      "minus_T_RR" -> dc[o~sp~3,-6,l@1,l@2,i~sp~1] dc[q~sp~4,-6,l@1,l@2,q~sp~2]
   };
   npfU = npfU~WilsonCoeffs`InterfaceToMatching~dim6[in,out,SARAH`UpQuark];
   npfD = npfD~WilsonCoeffs`InterfaceToMatching~dim6[in,out,SARAH`DownQuark];

   Print[">>npf>> calculation for ",`cxx`in," to ",`cxx`out," conversion done."];

   Print["<<npf<< c++ code calculation for ",`cxx`in," to ",`cxx`out," conversion started ..."];

   codeU = NPointFunctions`CreateCXXFunctions[
      npfU,
      `cxx`classU,
      SARAH`Delta,
      dim6[in,out,SARAH`UpQuark]
   ][[2]];
   codeD = NPointFunctions`CreateCXXFunctions[
      npfD,
      `cxx`classD,
      SARAH`Delta,
      dim6[in,out,SARAH`DownQuark]
   ][[2]];

   Print[">>npf>> c++ code calculation for ",`cxx`in," to ",`cxx`out," conversion done."];
   {
      DeleteDuplicates@Join[
         NPointFunctions`VerticesForNPointFunction@npfU,
         NPointFunctions`VerticesForNPointFunction@npfD
      ],
      NPointFunctions`CreateCXXHeaders[],
      codeU<>"\n\n"<>codeD
   }];
`npf`create // Utils`MakeUnknownInputDefinition;
`npf`create ~ SetAttributes ~ {Locked,Protected};

End[];
EndPackage[];
(* Allows define any name for symbols inside the package if $ContextPath
   contains only System during the package body.
 *)
$ContextPath = DeleteCases[$ContextPath, "LToLConversion`"];
(* Allows to separate the body of the package into different files for better
   readability.
 *)
Unprotect@$Packages;
$Packages = DeleteCases[$Packages, "LToLConversion`"];
Protect@$Packages;
