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
            FeynArts`Field[7|8] -> getField[ds,1] /. i:`type`indexGeneration:>Blank[],
            "initial SED: remove initial lepton in the loop"
         },
         {
            `topologyQ`self1pinguinT,
            FeynArts`DiagramSelect,
            FreeQ[#,FeynArts`Field[7|8] -> FeynArts`V]&,
            "initial SED: remove vector bosons in the loop"
         },
         {
            `topologyQ`self3pinguinT,
            removeGenericInsertionsBy,
            FeynArts`Field[7|8] -> -getField[ds,3] /. i:`type`indexGeneration:>Blank[],
            "final SED: remove final lepton in the loop"
         },
         {
            `topologyQ`self3pinguinT,
            FeynArts`DiagramSelect,
            FreeQ[#,FeynArts`Field[7|8] -> FeynArts`V]&,
            "final SED: remove vector bosons in the loop"
         },
         {
           `topologyQ`trianglepinguinT,
            removeGenericInsertionsBy,
            FeynArts`Field[6|7] -> getField[ds,1] /. i:`type`indexGeneration:>Blank[],
            "penguins: remove initial lepton in triangle loop"
         },
         {
           `topologyQ`trianglepinguinT,
            FeynArts`DiagramSelect,
            FreeQ[#, FeynArts`Field[6|7|8] -> FeynArts`V]&,
            "penguins: remove vector bosons in triangle loop"
         }
      },
      {
         {
            `topologyQ`pinguinT,
            FeynArts`DiagramSelect,
            FreeQ[#,FeynArts`Field@5 -> FeynArts`V]&,
            "penguins: remove tree-like vector bosons"
         }
      }
   },
   FourFermionScalarPenguins -> {
      {
         {
            `topologyQ`self1pinguinT,
            removeGenericInsertionsBy,
            FeynArts`Field[7|8] -> getField[ds,1] /. i:`type`indexGeneration:>Blank[],
            "initial SED: remove initial lepton in the loop"
         },
         {
            `topologyQ`self1pinguinT,
            FeynArts`DiagramSelect,
            FreeQ[#,FeynArts`Field[7|8] -> FeynArts`V]&,
            "initial SED: remove vector bosons in the loop"
         },
         {
            `topologyQ`self3pinguinT,
            removeGenericInsertionsBy,
            FeynArts`Field[7|8] -> -getField[ds,3] /. i:`type`indexGeneration:>Blank[],
            "final SED: remove final lepton in the loop"
         },
         {
            `topologyQ`self3pinguinT,
            FeynArts`DiagramSelect,
            FreeQ[#,FeynArts`Field[7|8] -> FeynArts`V]&,
            "final SED: remove vector bosons in the loop"
         },
         {
           `topologyQ`trianglepinguinT,
            removeGenericInsertionsBy,
            FeynArts`Field[6|7] -> getField[ds,1] /. i:`type`indexGeneration:>Blank[],
            "penguins: remove initial lepton in triangle loop"
         },
         {
           `topologyQ`trianglepinguinT,
            FeynArts`DiagramSelect,
            FreeQ[#, FeynArts`Field[6|7|8] -> FeynArts`V]&,
            "penguins: remove vector bosons in triangle loop"
         }
      },
      {
         {
            `topologyQ`pinguinT,
            FeynArts`DiagramSelect,
            FreeQ[#,FeynArts`Field@5 -> FeynArts`S]&,
            "penguins: remove tree-like scalar bosons"
         }
      }
   },
   FourFermionFlavourChangingBoxes -> {
      {
         {
            `topologyQ`boxS,
            removeGenericInsertionsBy,
            FeynArts`Field[6] -> getField[ds,1] /. i:`type`indexGeneration:>Blank[],
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
            FeynArts`Field[5] -> getField[ds,1] /. i:`type`indexGeneration:>Blank[],
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
            "pinguins: remove tree-like massless vector bosons",
            `topologyQ`pinguinT,
            {UnsameQ, genericMass[FeynArts`V, 5], 0}
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
      "initial SED: skip initial lepton in sum",
      `topologyQ`self1pinguinT,
      {6, getField[ds,1]}
   },
   {
      "final SED: skip final lepton in sum",
      `topologyQ`self3pinguinT,
      {6, getField[ds,3]}
   }
};

`settings`massless[ds:`type`diagramSet] := {
   {
      "initial SED: use explicit final lepton mass",
      `topologyQ`self1pinguinT,
      {Append, FeynArts`F[6] :> 3}
   },
   {
      "initial SED: keep initial lepton mass untouched",
      `topologyQ`self1pinguinT,
      {Hold, 1}
   }
};

`settings`momenta = {
   `topologyQ`pinguinT -> 2,
   `topologyQ`boxS -> 2,
   `topologyQ`boxU -> 2
};

`settings`regularization = {
   `topologyQ`boxS -> D,
   `topologyQ`boxU -> D
};

`settings`order = {3, 1, 4, 2};

`settings`chains = {
   {ExceptLoops, OperatorsOnly} ->
   {
      1[k@4]                        -> 0,
                     2[k@1]         -> 0,
      1[l@1]         2[k@1,l@1]     -> 0,
      1[k@4,l@1]     2[l@1]         -> 0,
      1[k@4,l@1]     2[k@1,l@1]     -> 0,
      1[k@4,l@1,l@2] 2[k@1,l@1,l@2] -> 0,
      1[k@1,k@4,l@1] 2[l@1]         -> 0,
      1[l@1]         2[k@1,k@4,l@1] -> 0,
      1[k@1,l@1]     2[l@1]         -> 0,
      1[l@1,l@2]     2[k@1,l@1,l@2] -> 0
   }
};

