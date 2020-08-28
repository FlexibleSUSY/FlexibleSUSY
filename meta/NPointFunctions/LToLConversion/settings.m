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

`settings`topologyReplacements = {
   FourFermionMassiveVectorPenguins -> (`topologyQ`pinguinT@#&),
   FourFermionScalarPenguins -> (`topologyQ`pinguinT@#&),
   FourFermionFlavourChangingBoxes -> (`topologyQ`box@#&)
};

`settings`diagrams[ds:`type`diagramSet] := {
   FourFermionMassiveVectorPenguins -> {
      {
         {
            `topologyQ`self1pinguinT,
            removeGenericInsertionsBy,
            FeynArts`Field[7|8] -> getField[ds,1] /. i:`type`indexGen:>Blank[],
            "t-penguins: remove leptons in initial SED loop"
         },
         {
            `topologyQ`self1pinguinT,
            FeynArts`DiagramSelect,
            FreeQ[#,FeynArts`Field[7|8] -> FeynArts`V]&,
            "t-penguins: remove vector bosons in initial SED loop"
         },
         {
            `topologyQ`self3pinguinT,
            removeGenericInsertionsBy,
            FeynArts`Field[7|8] -> getField[ds,3] /. i:`type`indexGen:>Blank[],
            "t-penguins: remove leptons in final SED loop"
         },
         {
            `topologyQ`self3pinguinT,
            FeynArts`DiagramSelect,
            FreeQ[#,FeynArts`Field[7|8] -> FeynArts`V]&,
            "t-penguins: remove vector bosons in final SED loop"
         },
         {
           `topologyQ`trianglepinguinT,
            removeGenericInsertionsBy,
            FeynArts`Field[6|7] -> getField[ds,1] /. i:`type`indexGen:>Blank[],
            "t-penguins: remove leptons in triangle loop"
         },
         {
           `topologyQ`trianglepinguinT,
            FeynArts`DiagramSelect,
            FreeQ[#, FeynArts`Field[6|7|8] -> FeynArts`V]&,
            "t-penguins: remove vector bosons in triangle loop"
         }
      },
      {
         {
            `topologyQ`pinguinT,
            FeynArts`DiagramSelect,
            FreeQ[#,FeynArts`Field@5 -> FeynArts`V]&,
            "t-penguins: remove tree-like vector bosons"
         }
      }
   },
   FourFermionScalarPenguins -> {
      {
         {
            `topologyQ`self1pinguinT,
            removeGenericInsertionsBy,
            FeynArts`Field[7|8] -> getField[ds,1] /. i:`type`indexGen:>Blank[],
            "t-penguins: remove leptons in initial SED loop"
         },
         {
            `topologyQ`self1pinguinT,
            FeynArts`DiagramSelect,
            FreeQ[#,FeynArts`Field[7|8] -> FeynArts`V]&,
            "t-penguins: remove vector bosons in initial SED loop"
         },
         {
            `topologyQ`self3pinguinT,
            removeGenericInsertionsBy,
            FeynArts`Field[7|8] -> getField[ds,3] /. i:`type`indexGen:>Blank[],
            "t-penguins: remove leptons in final SED loop"
         },
         {
            `topologyQ`self3pinguinT,
            FeynArts`DiagramSelect,
            FreeQ[#,FeynArts`Field[7|8] -> FeynArts`V]&,
            "t-penguins: remove vector bosons in final SED loop"
         },
         {
           `topologyQ`trianglepinguinT,
            removeGenericInsertionsBy,
            FeynArts`Field[6|7] -> getField[ds,1] /. i:`type`indexGen:>Blank[],
            "t-penguins: remove leptons in triangle loop"
         },
         {
           `topologyQ`trianglepinguinT,
            FeynArts`DiagramSelect,
            FreeQ[#, FeynArts`Field[6|7|8] -> FeynArts`V]&,
            "t-penguins: remove vector bosons in triangle loop"
         }
      },
      {
         {
            `topologyQ`pinguinT,
            FeynArts`DiagramSelect,
            FreeQ[#,FeynArts`Field@5 -> FeynArts`S]&,
            "t-penguins: remove tree-like scalar bosons"
         }
      }
   },
   FourFermionFlavourChangingBoxes -> {
      {
         {
            `topologyQ`boxS,
            removeGenericInsertionsBy,
            FeynArts`Field[6] -> getField[ds,1] /. i:`type`indexGen:>Blank[],
            "s-boxes: remove loops with initial lepton"
         },
         {
            `topologyQ`boxS,
            FeynArts`DiagramSelect,
            FreeQ[#, FeynArts`Field[5|6|7|8] -> FeynArts`V]&,
            "s-boxes: remove loops with vector bosons"
         },
         {
            `topologyQ`boxU,
            removeGenericInsertionsBy,
            FeynArts`Field[5] -> getField[ds,1] /. i:`type`indexGen:>Blank[],
            "u-boxes: remove loops with initial lepton"
         },
         {
            `topologyQ`boxU,
            FeynArts`DiagramSelect,
            FreeQ[#, FeynArts`Field[5|6|7|8] -> FeynArts`V]&,
            "u-boxes: remove loops with vector bosons"
         }
      },
      {}
   }
};

`settings`amplitudes = {
   FourFermionMassiveVectorPenguins -> {
      {
         {
            `topologyQ`pinguinT,
            {genericMass[FeynArts`V, 5], UnsameQ, 0},
            "t-pinguins: remove tree-like massless vector bosons"
         }
      },
      {}
   },
   FourFermionScalarPenguins -> {
      {},
      {}
   },
   FourFermionFlavourChangingBoxes -> {
      {},
      {}
   }
};

`settings`sum[ds:`type`diagramSet] := {
   {
      "t-pinguins: no external lepton in initial SED bridge",
      `topologyQ`self1pinguinT,
      {6, getField[ds,1]}
   },
   {
      "t-pinguins: no external lepton in final SED bridge",
      `topologyQ`self3pinguinT,
      {6, getField[ds,3]}
   }
};
