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
      {  {  {  "penguins: remove external leptons from loops",
               `topologyQ`penguinT,
               FreeQ[LoopFields@##, fieldPattern[#3, 1|3]]&},
            {  "penguins: remove vector bosons from loops",
               `topologyQ`penguinT,
               FreeQ[LoopFields@##, FeynArts`V]&}},
         {  {  "penguins: remove tree-like vector bosons",
               `topologyQ`penguinT,
               FreeQ[TreeFields@##, FeynArts`V]&}}},
   Scalars ->
      {  {  {  "penguins: remove external leptons from loops",
               `topologyQ`penguinT,
               FreeQ[LoopFields@##, fieldPattern[#3, 1|3]]&},
            {  "penguins: remove vector bosons from loops",
               `topologyQ`penguinT,
               FreeQ[LoopFields@##, FeynArts`V]&}},
         {  {  "penguins: remove tree-like scalar bosons",
               `topologyQ`penguinT,
               FreeQ[TreeFields@##, FeynArts`S]&}}},
   Boxes ->
      {  {
            {  "boxes: remove external leptons from loops",
               `topologyQ`box,
               FreeQ[LoopFields@##, fieldPattern[#3, 1|3]]&},
            {  "boxes: remove vector bosons from loops",
               `topologyQ`box,
               FreeQ[LoopFields@##, FeynArts`V]&}},
         {}}};

`settings`amplitudes =
{  Vectors ->
      {  {  {  "penguins: remove tree-like massless vector bosons",
               `topologyQ`penguinT,
               FreeQ[#, genericMass[FeynArts`V, 5] -> 0]&}},
         {}}};

`settings`sum =
{  "initial SED: skip initial lepton in sum"[
      `topologyQ`self1penguinT,
      {6, getField[#, 1]&}],
   "final SED: skip final lepton in sum"[
      `topologyQ`self3penguinT,
      {6, getField[#, 3]&}]};

`settings`massless =
{  "initial SED: use explicit final lepton mass"[
      `topologyQ`self1penguinT,
      {Append, FeynArts`F[6] :> 3}],
   "initial SED: keep initial lepton mass untouched"[
      `topologyQ`self1penguinT,
      {Hold, 1}]};

`settings`momenta =
{  `topologyQ`penguinT -> 2,
   `topologyQ`boxS -> 2,
   `topologyQ`boxU -> 2};

`settings`regularization =
{  `topologyQ`boxS -> D,
   `topologyQ`boxU -> D};

`settings`order = {3, 1, 4, 2};

`settings`chains =
{  {ExceptLoops, OperatorsOnly} ->
      {  1[k@4, ___] -> 0, 2[k@1, ___] -> 0}};
