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

BeginPackage["References`"];

AddReference::usage = "Add a reference by a bibtex \"key\"";
CreateCitation::usage = "Creates \\cite{} command with all added keys.";

Begin["`Private`"];

fsreferences = {};

AddReference[key_] := (
    fsreferences = DeleteDuplicates[Append[fsreferences, ToString[key]]];
);

AddReference[r__] := AddReference /@ { r };

CreateCitation[{}] := "";

CreateCitation[refs_List] :=
    "\\cite{" <> StringRiffle[refs, ","] <> "}";

CreateCitation[] := CreateCitation[fsreferences];

End[];

EndPackage[];
