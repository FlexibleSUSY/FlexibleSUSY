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
      {"boxes: remove loop vectors",
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
         FreeQ[#, genericMass[FeynArts`V, 5] -> 0]&}
   }
};

`settings`sum :=
{  {  "SED T: skip initial lepton in sum",
      inSelfT,
      {6, Field[#3, 1]&}},
   {  "SED T: skip final lepton in sum",
      outSelfT,
      {6, Field[#3, 3]&}}};

`settings`massless =
{  {  "SED T: use explicit final lepton mass",
      inSelfT,
      {Append, FeynArts`F[6] :> 3}},
   {  "SED T: keep initial lepton mass untouched",
      inSelfT,
      {Hold, 1}}};

`settings`momenta =
{  penguinT -> 2,
   boxAll -> 2};

`settings`regularization =
{  boxAll -> D};

`settings`order = {3, 1, 4, 2};

`settings`chains =
{  {  ExceptLoops} ->
      {  1[k[4|2], ___] -> 0, 2[k[3|1], ___] -> 0}};
