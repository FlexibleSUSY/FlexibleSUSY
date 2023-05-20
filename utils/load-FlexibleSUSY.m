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

(*
   This script loads SARAH and FlexibleSUSY.
   Usage:

     Get["utils/load-FlexibleSUSY.m"];
*)

AppendTo[$Path, FileNameJoin[{Directory[], "meta"}]];

Needs["SARAH`"];
Needs["FlexibleSUSY`", FileNameJoin[{Directory[], "meta", "FlexibleSUSY.m"}]];

SARAH`SARAH[OutputDirectory] = FileNameJoin[{Directory[], "Output"}];

SARAH`SARAH[InputDirectories] = {
    FileNameJoin[{Directory[], "sarah"}],
    ToFileName[{$sarahDir, "Models"}]
};

Print["Current working directory: ", Directory[]];
Print["SARAH output directory   : ", SARAH`SARAH[OutputDirectory]];
Print["SARAH models directories : ", SARAH`SARAH[InputDirectories]];
