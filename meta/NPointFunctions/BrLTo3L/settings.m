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

`settings`topology =
{  Vectors -> (`topologyQ`penguinT@#&),
   Scalars -> (`topologyQ`penguinT@#&),
   Boxes -> (`topologyQ`box@#&)};

`settings`diagrams =
{  Vectors ->
      {  {  "penguins: remove external leptons from loops"[
               `topologyQ`penguinT,
               FreeQ[FeynArts`LoopFields@##, fieldPattern[#3, 1|3]]&],
            "penguins: remove vector bosons from loops"[
               `topologyQ`penguinT,
               FreeQ[FeynArts`LoopFields@##, FeynArts`V]&]},
         {  "penguins: remove tree-like vector bosons"[
               `topologyQ`penguinT,
               FreeQ[FeynArts`TreeFields@##, FeynArts`V]&]}},
   Scalars ->
      {  {  "penguins: remove external leptons from loops"[
               `topologyQ`penguinT,
               FreeQ[FeynArts`LoopFields@##, fieldPattern[#3, 1|3]]&],
            "penguins: remove vector bosons from loops"[
               `topologyQ`penguinT,
               FreeQ[FeynArts`LoopFields@##, FeynArts`V]&]},
         {  "penguins: remove tree-like scalar bosons"[
               `topologyQ`penguinT,
               FreeQ[FeynArts`TreeFields@##, FeynArts`S]&]}},
   Boxes ->
      {  {  "boxes: remove external leptons from loops"[
               `topologyQ`box,
               FreeQ[FeynArts`LoopFields@##, fieldPattern[#3, 1|2|3|4]]&],
            "boxes: remove vector bosons from loops"[
               `topologyQ`box,
               FreeQ[FeynArts`LoopFields@##, FeynArts`V]&]},
         {}}};

`settings`amplitudes =
{  Vectors ->
      {  {  {  "penguins: remove tree-like massless vector bosons",
               `topologyQ`penguinT,
               FreeQ[#, genericMass[FeynArts`V, 5] -> 0]&}},
         {}}};

`settings`sum :=
{  "SED T: skip initial lepton in sum"[
      `topologyQ`self1penguinT,
      {6, getField[#, 1]&}],
   "SED T: skip final lepton in sum"[
      `topologyQ`self3penguinT,
      {6, getField[#, 3]&}]};

`settings`massless =
{  "SED T: use explicit final lepton mass"[
      `topologyQ`self1penguinT,
      {Append, FeynArts`F[6] :> 3}],
   "SED T: keep initial lepton mass untouched"[
      `topologyQ`self1penguinT,
      {Hold, 1}]};

`settings`momenta =
{  `topologyQ`penguinT -> 2,
   `topologyQ`box -> 2};

`settings`regularization =
{  `topologyQ`box -> D};

`settings`order = {3, 1, 4, 2};

`settings`chains =
{  {  ExceptLoops} ->
      {  1[k[4|2], ___] -> 0, 2[k[3|1], ___] -> 0}};
