topologies[0] = {
   Vectors -> treeAll,
   Scalars -> treeAll
};

topologies[1] = {
   Vectors -> penguinT,
   Scalars -> penguinT,
   Boxes -> boxAll
};

diagrams[0, Minus] = {
   Vectors -> {
      {"remove vector bosons", treeAll, FreeQ[TreeFields@##, FeynArts`V]&}
   },
   Scalars -> {
      {"remove scalar bosons", treeAll, FreeQ[TreeFields@##, FeynArts`S]&}
   }
};

diagrams[1, Plus] = {
   Vectors -> {
      {"penguins: remove external leptons",
         penguinT,
         FreeQ[LoopFields@##, FieldPattern[#3, 1|3]]&
      },
      {"penguins: remove loop vectors",
         penguinT,
         FreeQ[LoopFields@##, FeynArts`V]&
      }
   },
   Scalars -> {
      {"penguins: remove external leptons",
         penguinT,
         FreeQ[LoopFields@##, FieldPattern[#3, 1|3]]&
      },
      {"penguins: remove loop vectors",
         penguinT,
         FreeQ[LoopFields@##, FeynArts`V]&
      }
   },
   Boxes -> {
      {"boxes: remove external leptons",
         boxAll,
         FreeQ[LoopFields@##, FieldPattern[#3, 1|2|3|4]]&}
      ,
      {"boxes: remove vector bosons from loops",
         boxAll,
         FreeQ[LoopFields@##, FeynArts`V]&
      }
   }
};

diagrams[1, Minus] = {
   Vectors -> {
      {"penguins: remove tree vectors", penguinT, FreeQ[TreeFields@##, FeynArts`V]&}
   },
   Scalars -> {
      {"penguins: remove tree scalars", penguinT, FreeQ[TreeFields@##, FeynArts`S]&}
   }
};

amplitudes[0, Plus] = {
   Vectors -> {
      {"remove photons", treeAll, FreeQ[#, genericMass[FeynArts`V, 5] -> 0]&}
   }
};

amplitudes[1, Plus] = {
   Vectors -> {
      {"penguins: remove tree photons",
         penguinT,
         FreeQ[#, genericMass[FeynArts`V, 5] -> 0]&
      }
   }
};

`settings`sum =
{  {  "initial SED: skip initial lepton in sum",
      inSelfT,
      {6, Field[#3, 1]&}},
   {  "final SED: skip final lepton in sum",
      outSelfT,
      {6, Field[#3, 3]&}}};

`settings`massless =
{  {  "initial SED: use explicit final lepton mass",
      inSelfT,
      {Append, FeynArts`F[6] :> 3}},
   {  "initial SED: keep initial lepton mass untouched",
      inSelfT,
      {Hold, 1}}};

`settings`momenta =
{  penguinT -> 2,
   boxS -> 2,
   boxU -> 2};

`settings`regularization =
{  boxS -> D,
   boxU -> D};

`settings`order = {3, 1, 4, 2};

`settings`chains =
{  {ExceptLoops, OperatorsOnly} ->
      {  1[k@4, ___] -> 0, 2[k@1, ___] -> 0}};
