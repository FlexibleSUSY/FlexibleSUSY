topologies[0] = {
   {Vectors, Scalars} -> treeAll
};

topologies[1] = {
   {Vectors, Scalars} -> penguinT,
   Boxes -> boxAll
};

diagrams[0, Minus] = {
   Vectors -> {
      treeAll -> {"remove vector bosons",
         FreeQ[TreeFields@##, FeynArts`V]&}},
   Scalars -> {
      treeAll -> {"remove scalar bosons",
         FreeQ[TreeFields@##, FeynArts`S]&}}
};

diagrams[1, Plus] = {
   Vectors -> {
      penguinT -> {"penguins: remove external leptons",
         FreeQ[LoopFields@##, FieldPattern[#3, 1|3]]&},
      penguinT -> {"penguins: remove loop vectors",
         FreeQ[LoopFields@##, FeynArts`V]&}},
   Scalars -> {
      penguinT -> {"penguins: remove external leptons",
         FreeQ[LoopFields@##, FieldPattern[#3, 1|3]]&},
      penguinT -> {"penguins: remove loop vectors",
         FreeQ[LoopFields@##, FeynArts`V]&}},
   Boxes -> {
      boxAll -> {"boxes: remove external leptons",
         FreeQ[LoopFields@##, FieldPattern[#3, 1|2|3|4]]&},
      boxAll -> {"boxes: remove loop vectors",
         FreeQ[LoopFields@##, FeynArts`V]&}}
};

diagrams[1, Minus] = {
   Vectors -> {
      penguinT -> {"penguins: remove tree vectors",
         FreeQ[TreeFields@##, FeynArts`V]&}},
   Scalars -> {
      penguinT -> {"penguins: remove tree scalars",
         FreeQ[TreeFields@##, FeynArts`S]&}}
};

amplitudes[0, Plus] = {
   Vectors -> {
      treeAll -> {"remove photons",
         FreeQ[#, InternalMass[FeynArts`V, 5] -> 0]&}}
};

amplitudes[1, Plus] = {
   Vectors -> {
      penguinT -> {"penguins: remove tree photons",
         FreeQ[#, InternalMass[FeynArts`V, 5] -> 0]&}}
};

order[] = {3, 1, 4, 2};

regularization[1] = {
   boxAll -> D
};

momenta[1] = {
   penguinT -> 2,
   boxAll -> 2
};

sum[1] = {
   inSelfT -> {"in-sed: skip initial lepton",
      {6, Field[#3, 1]&}},
   outSelfT -> {"out-sed: skip final lepton",
      {6, Field[#3, 3]&}}
};

chains[1] = {
   ExceptLoops -> {1[k[4|2], ___] -> 0, 2[k[3|1], ___] -> 0}
};

mass[1] = {
   inSelfT -> {"in-sed: use explicit final lepton mass",
      {Append, InternalMass[FeynArts`F, 6] :> ExternalMass[3]}},
   inSelfT -> {"in-sed: keep initial lepton mass untouched",
      {Hold, ExternalMass[1]}}
};
