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

BeginPackage@"WilsonCoeffs`";

{InterfaceToMatching,neglectBasisElements};

Begin["`internal`"];
`type`npf = NPointFunctions`internal`type`npf;
getGenericSums = NPointFunctions`internal`getGenericSums;
getSubexpressions = NPointFunctions`internal`getSubexpressions;
getName = NPointFunctions`internal`getName;

neglectBasisElements::usage=
"@brief Deletes specified basis elements with not anymore used subexpressions.
@param <npf> npf object to modify.
@param <{Rule[_String,_]..}> operatorBasis list with
{string name,fermion chain multiplication} pairs.
@returns <npf> object";
neglectBasisElements[obj:`type`npf, operatorBasis:{Rule[_String,_]..}]:=
Module[
   {
      basis = Last /@ findFermionChains[getSubexpressions[obj],operatorBasis],
      objNew,positionsToDelete,
      subsOnly = DeleteDuplicates@Cases[#,Alternatives@@(getName@getSubexpressions@obj),Infinity]&
   },
   If[basis === {},Return@obj];
   objNew = ReplaceAll[obj,#->0 &/@ basis];
   positionsToDelete = If[Length@#===1,#[[1]]~Take~3,##&[]] &@ Position[objNew,#] &/@ Complement[getName@getSubexpressions@obj,subsOnly@objNew[[2,1,1]]];
   Delete[objNew,positionsToDelete]
];
SetAttributes[neglectBasisElements,{Locked,Protected}];

InterfaceToMatching::usage=
"@brief Transforms GenericSum accordig to a given basis.
@param NPF NPF object.
@param operatorBasis list with {string name,fermion chain multiplication} pairs.
@returns Corresponding GenericSum which matches to a given basis.
@note the name convention of the chiral basis follows the FormCalc convention.";
InterfaceToMatching::errUnknownInput =
"Correct input has the folliwing form:
InterfaceToMatching@@{ <npf object>, {<string -> basis expression>..} }
and not
InterfaceToMatching@@`1`.";
InterfaceToMatching[obj:`type`npf, operatorBasis:{Rule[_String,_]}] := obj;
InterfaceToMatching[obj:`type`npf, operatorBasis:{Rule[_String,_]..}] :=
Module[{basis, coefficientsWilson},
   basis = findFermionChains[getSubexpressions@obj, operatorBasis];
   coefficientsWilson = removeFermionChains[createNewNPF[obj, basis]];
   coefficientsWilson
];
InterfaceToMatching[x___] :=
Utils`AssertOrQuit[False,InterfaceToMatching::errUnknownInput,{x}];

findFermionChains::usage =
"@brief Searches the FermionChains in the abbreviations rules.
@param subs substitution rules of the form {Rule[_,_]..}.
@param chiralBasis  name basis element rules of the form {Rule[_String,_]..}.
@returns List of <string chain name>->NPointFunctions`internal`mat[F#] pairs.";
findFermionChains[subs:{Rule[_,_]..}, chiralBasis:{Rule[_String,_]..}] :=
Module[
   {
      warning = If[!$Notebooks,"\033[1;33mWarning\033[1;0m",Style["Warning",Yellow]],
      (*@note there is F# <-> chain correspondence*)
      basisPos = Position[subs, #]& /@ chiralBasis[[All, 2]],
      i
   },
   Table[
      If[basisPos[[i]] === {},
         Print[warning,": " <> chiralBasis[[i,1]] <> " is absent in GenericSums."];
         chiralBasis[[i,1]]->NPointFunctions`internal`mat[],
         (*else*)
         chiralBasis[[i,1]]->NPointFunctions`internal`mat[Extract[subs,{basisPos[[i,1,1]],basisPos[[i,1,2]]-1}]]
      ],
      {i,Length@basisPos}]
];

createNewNPF::usage =
"@brief Extracts the coefficients for a given basis and NPF object.";
createNewNPF[obj:`type`npf,
   chiralBasis:{Rule[_String,_NPointFunctions`internal`mat]..}
] :=
Module[
   {
      newSums,
      newNPF=obj
   },
   newSums = extractCoeffs[#,chiralBasis]& /@ getGenericSums[obj];
   newNPF[[2, 1, 1]] = newSums;
   newNPF
];

extractCoeffs::errRemainingExpression =
"Probably missing basis element for fermionic chain. During calculation expression
`1`
was not taken into account appropriately.";
extractCoeffs[
   NPointFunctions`GenericSum[{expr_},sumFields:{__}],
   operators:{Rule[_String,_NPointFunctions`internal`mat]..}
] :=
Module[
   {
      coefficients = Coefficient[expr, #]& /@ operators[[All,2]],
      check
   },
   check = Expand[expr-coefficients.(operators[[All,2]])];
   If[check=!=0,Utils`AssertOrQuit[False,extractCoeffs::errRemainingExpression,FullForm@check]];
   NPointFunctions`GenericSum[coefficients,sumFields]
];

removeFermionChains::usage =
"@brief Removes DiracChains from the abbreviations rules.";
removeFermionChains[npointExpression:`type`npf] :=
Module[{pos},
   pos = Take[#, 3]& /@ Position[npointExpression, NPointFunctions`internal`dc];
   Delete[npointExpression, pos]
];

End[];
EndPackage[];
