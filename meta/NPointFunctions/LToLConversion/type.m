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

Off[LToLConversion`arguments::shdw];
BeginPackage["LToLConversion`"];

namespace::usage = "
@brief Returns a namespace for C++ code of a given observable.
@note Convention: the name of a corresponding C++ file without ``.cpp``.";

arguments::usage = "
@param in A symbol for incoming lepton.
@param inN A symbol for incoming lepton generation.
@param out A symbol for outgoing lepton.
@param outN A symbol for outgoing lepton generation.
@param nucleus A symbol for a nucleus.
@param contribution A symbol for a desired contribution.
@return A sequence of arguments for the observable.";

Begin["`Private`"];

namespace[File] := "ltolconversion";
namespace[C] := FlexibleSUSY`FSModelName<>"_"<>namespace[File]<>"::";
namespace // Utils`MakeUnknownInputDefinition;
namespace ~ SetAttributes ~ {Protected, Locked};

Off@RuleDelayed::rhs;
arguments[in_Symbol[inN_Symbol],
          out_Symbol[outN_Symbol],
          nucleus_Symbol,
          contribution_Symbol,
          loopN_Symbol] :=
   Sequence[Rule[(in_Symbol?TreeMasses`IsLepton)[inN_Integer],
                 (out_Symbol?TreeMasses`IsLepton)[outN_Integer]],
            nucleus_Symbol,
            contribution:_Symbol|{__Symbol},
            loopN:0|1];
On@RuleDelayed::rhs;
arguments // Utils`MakeUnknownInputDefinition;
arguments ~ SetAttributes ~ {Protected, Locked};

`type`observable = FlexibleSUSYObservable`LToLConversion@
   arguments[in@iIn, out@iOut, nucleus, con, loopN];
`type`observable ~ SetAttributes ~ {Protected, Locked};

End[];
Block[{$ContextPath}, EndPackage[]];
