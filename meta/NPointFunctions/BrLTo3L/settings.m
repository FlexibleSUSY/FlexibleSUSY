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

`settings`topologyReplacements =
{  MassiveVectorPenguins -> (`topologyQ`penguinT@#&),
   ScalarPenguins -> (`topologyQ`penguinT@#||`topologyQ`penguinU@#&),
   FlavourChangingBoxes -> (`topologyQ`box@#&)};

`settings`diagrams =
{  MassiveVectorPenguins ->
      {  {  "penguins: remove external leptons from loops"[
               `topologyQ`penguinT|`topologyQ`penguinU,
               FreeQ[FeynArts`LoopFields@##, fieldPattern[#3, 1|3]]&],
            "penguins: remove vector bosons from loops"[
               `topologyQ`penguinT|`topologyQ`penguinU,
               FreeQ[FeynArts`LoopFields@##, FeynArts`V]&]},
         {  "penguins: remove tree-like vector bosons"[
               `topologyQ`penguinT|`topologyQ`penguinU,
               FreeQ[FeynArts`TreeFields@##, FeynArts`V]&]}},
   ScalarPenguins ->
      {  {  "penguins: remove external leptons from loops"[
               `topologyQ`penguinT|`topologyQ`penguinU,
               FreeQ[FeynArts`LoopFields@##, fieldPattern[#3, 1|3]]&],
            "penguins: remove vector bosons from loops"[
               `topologyQ`penguinT|`topologyQ`penguinU,
               FreeQ[FeynArts`LoopFields@##, FeynArts`V]&]},
         {  "penguins: remove tree-like scalar bosons"[
               `topologyQ`penguinT|`topologyQ`penguinU,
               FreeQ[FeynArts`TreeFields@##, FeynArts`S]&]}},
   FlavourChangingBoxes ->
      {  {  "boxes: remove external leptons from loops"[
               `topologyQ`box,
               FreeQ[FeynArts`LoopFields@##, fieldPattern[#3, 1|2|3|4]]&],
            "boxes: remove vector bosons from loops"[
               `topologyQ`box,
               FreeQ[FeynArts`LoopFields@##, FeynArts`V]&]},
         {}}};

`settings`amplitudes =
{  MassiveVectorPenguins ->
      {  {  "penguins: remove tree-like massless vector bosons"[
               `topologyQ`penguinT,
               {UnsameQ, genericMass[FeynArts`V, 5], 0}],
            "penguins: remove tree-like massless vector bosons"[
               `topologyQ`penguinU,
               {UnsameQ, genericMass[FeynArts`V, 5], 0}]},
         {}}};

`settings`sum[ds:`type`diagramSet] :=
{  {  "SED T: skip initial lepton in sum",
      `topologyQ`self1penguinT,
      {6, getField[ds, 1]}},
   {  "SED T: skip final lepton in sum",
      `topologyQ`self3penguinT,
      {6, getField[ds, 3]}},
   {  "SED U: skip initial lepton in sum",
      `topologyQ`self1penguinU,
      {6, getField[ds, 1]}},
   {  "SED U: skip final lepton in sum",
      `topologyQ`self4penguinU,
      {6, getField[ds, 4]}}};

`settings`massless[ds:`type`diagramSet] :=
{  {  "SED T: use explicit final lepton mass",
      `topologyQ`self1penguinT,
      {Append, FeynArts`F[6] :> 3}},
   {  "SED T: keep initial lepton mass untouched",
      `topologyQ`self1penguinT,
      {Hold, 1}},
   {  "SED U: use explicit final lepton mass",
      `topologyQ`self1penguinU,
      {Append, FeynArts`F[6] :> 4}},
   {  "SED U: keep initial lepton mass untouched",
      `topologyQ`self1penguinU,
      {Hold, 1}}};

`settings`momenta =
{  `topologyQ`penguinT -> 2,
   `topologyQ`penguinU -> 4,
   `topologyQ`box -> 2};

`settings`regularization =
{  `topologyQ`penguinU -> D,
   `topologyQ`box -> D};

`settings`order = {3, 1, 4, 2};

`settings`chains =
{  {  ExceptLoops} ->
      {  1[k[4|2], ___] -> 0, 2[k[3|1], ___] -> 0}};
