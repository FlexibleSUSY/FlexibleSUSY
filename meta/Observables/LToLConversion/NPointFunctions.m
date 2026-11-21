topologies[0] = {
   {Vectors, Scalars} -> treeAll
};

topologies[1] = {
   {Vectors, Scalars} -> penguinT,
   Boxes -> boxAll
};

diagrams[0, Absent] = {
   Vectors -> {
      treeAll -> {"remove vector bosons", FreeQ[TreeFields@##, FeynArts`V]&}
   },
   Scalars -> {
      treeAll -> {"remove scalar bosons", FreeQ[TreeFields@##, FeynArts`S]&}
   }
};

diagrams[1, Present] = {
   Vectors -> {
      penguinT -> {"remove external leptons", FreeQ[LoopFields@##, FieldPattern[#, 1|3]]&},
      penguinT -> {"remove loop vectors", FreeQ[LoopFields@##, FeynArts`V]&}
   },
   Scalars -> {
      penguinT -> {"remove external leptons", FreeQ[LoopFields@##, FieldPattern[#, 1|3]]&},
      penguinT -> {"remove loop vectors", FreeQ[LoopFields@##, FeynArts`V]&}
   },
   Boxes -> {
      boxAll -> {"remove external leptons", FreeQ[LoopFields@##, FieldPattern[#, 1|2|3|4]]&},
      boxAll -> {"remove vector bosons from loops", FreeQ[LoopFields@##, FeynArts`V]&}
   }
};

diagrams[1, Absent] = {
   Vectors -> {
      penguinT -> {"remove tree vectors", FreeQ[TreeFields@##, FeynArts`V]&}
   },
   Scalars -> {
      penguinT -> {"remove tree scalars", FreeQ[TreeFields@##, FeynArts`S]&}
   }
};

amplitudes[0, Present] = {
   Vectors -> {
      treeAll -> {"remove photons", FreeQ[{##}, InternalMass[FeynArts`V, 5] -> 0]&}
   }
};

amplitudes[1, Present] = {
   Vectors -> {
      penguinT -> {"remove tree photons", FreeQ[{##}, InternalMass[FeynArts`V, 5] -> 0]&}
   }
};

order[] = {3, 1, 4, 2};

regularization[1] = {
   boxS -> D,
   boxU -> D
};

momenta[1] = {
   penguinT -> 2,
   boxS -> 2,
   boxU -> 2
};

sum[1] = {
   inSelfT  -> {"skip initial lepton", {6, Field[#, 1]&}},
   outSelfT -> {"skip final lepton",   {6, Field[#, 3]&}}
};

chains[1] = {
   {ExceptLoops, OperatorsOnly} -> {1[k[4], ___] -> 0, 2[k[1], ___] -> 0}
};

mass[1] = {
   inSelfT -> {"explicit final lepton mass", InternalMass[FeynArts`F, 6] :> ExternalMass[3]},
   inSelfT -> {"keep initial lepton mass untouched", Hold :> ExternalMass[1]}
};
