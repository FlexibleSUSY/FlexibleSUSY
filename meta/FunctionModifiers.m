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

BeginPackage["FunctionModifiers`"];

MakeAbstract::usage = "Adds virtual keyword (function specifier) to
 functions and makes it abstract (function body not defined).";

MakeOverride::usage = "Adds override function specifier to
 functions.";

MakeVirtual::usage = "Adds virtual keyword (function specifier) to
 functions.";

Begin["`Private`"];

MakeAbstract[str_String] :=
    Module[{ws = Whitespace ...},
           StringReplace[MakeVirtual[str], {
               ws ~~ Shortest["{" ~~ ___ ~~ "}"] ~~ ws ~~ ";" ... ~~ ws ~~ EndOfLine -> " = 0;",
               ws ~~ ")" ~~ p : (ws ~~ (CharacterRange["a", "z"] | ws) ...) ~~ ws ~~ ";" ~~ ws ~~ EndOfLine :> ")" <> p <> " = 0;"
           }]
          ];

MakeOverride[str_String] :=
    Module[{ws = Whitespace ...},
           StringReplace[str,
                         {
                             ws ~~ p : Shortest["{" ~~ ___ ~~ "}"] :> " override " <> p,
                             ")" ~~ ws ~~ "const;" -> ") const override;",
                             ");" -> ") override;"
                         }]
          ];

MakeVirtual[str_String] :=
    Module[{ws = Whitespace ...},
           StringReplace[#, {
               Repeated[("virtual" ~~ ws ..), {2, Infinity}] -> "virtual "
           }]& @
           StringReplace[#, {
               StartOfLine ~~ w : ws ~~ p_ :> w ~~ "virtual " <> p
           }]& @
           StringReplace[#, {
               StartOfLine ~~ Whitespace.. ~~ EndOfLine -> ""
           }]& @ str
          ];

End[];

EndPackage[];
