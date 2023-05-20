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

BeginPackage["FlexibleTower`", {"CConversion`"}];

GetModelTypes::usage = "Returns comma-separated list of model types";

Begin["`Private`"];

GetModelType[name_String, templatePar_String] :=
    name <> "<" <> templatePar <> ">";

GetSLHAModelType[name_String, templatePar_String] :=
    GetModelType[name, templatePar];

GetModelTypes[] :=
    Utils`StringJoinWithSeparator[
        Which[FlexibleSUSY`FlexibleEFTHiggs === True,
              { GetSLHAModelType[FlexibleSUSY`FSModelName, "Solver_type"],
                GetModelType["standard_model::StandardModel", "Two_scale"] },
              True,
              { GetSLHAModelType[FlexibleSUSY`FSModelName, "Solver_type"] }
             ],
        ", "
    ];

End[];

EndPackage[];
