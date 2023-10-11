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

BeginPackage["LToLConversion`"];

namespace::usage = "Returns a namespace for C++ code of a given observable.";

Begin["`Private`"];

namespace[File] := "ltolconversion";
namespace[C] := FlexibleSUSY`FSModelName<>"_"<>namespace[File]<>"::";
namespace // Utils`MakeUnknownInputDefinition;
namespace ~ SetAttributes ~ {Protected, Locked};

`type`observable = FlexibleSUSYObservable`LToLConversion[
   in_@iIn_ -> out_@iOut_, nucleus_, con_, loopN_];
`type`observable ~ SetAttributes ~ {Protected, Locked};

End[];
Block[{$ContextPath}, EndPackage[]];
