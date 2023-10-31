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

BeginPackage@"tools`";

zeroRules::usage = "
@brief Given a set of rules that map to zero and a set that does
       not map to zero, apply the zero-rules to the non-zero ones
       recursively until the non-zero rules do not change anymore.";

subWrite::usage = "";

Begin@"`Private`";

subWrite = NPointFunctions`Private`subWrite;
subWrite // Protect;

zeroRules[nonzeroRules:{Rule[_, _]...}, zeroRules:{Rule[_, 0]...}] :=
Module[{newNonzero, newZeroRules},
   newNonzero = Thread[
      Rule[nonzeroRules[[All,1]], nonzeroRules[[All,2]]/. zeroRules]];
   If[newNonzero === nonzeroRules, Return[{nonzeroRules, zeroRules}]];
   newZeroRules = Cases[newNonzero, HoldPattern[_-> 0]];
   newNonzero = Complement[newNonzero, newZeroRules];
   zeroRules[newNonzero, Join[zeroRules, newZeroRules]]
];
zeroRules // secure;

End[];
Block[{$ContextPath}, EndPackage[]];
