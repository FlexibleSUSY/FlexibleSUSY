topologies[0] = {
   Vectors -> (treeAll@#&),
   Scalars -> (treeAll@#&)
};

topologies[1] = {
   Vectors -> (penguinT@#&),
   Scalars -> (penguinT@#&),
   Boxes -> (boxAll@#&)
};

`settings`diagrams = {
   0 -> {
      Vectors -> {
         List[],
         List[
            {
               "remove vector bosons",
               treeAll,
               FreeQ[TreeFields@##, FeynArts`V]&
            }
         ]
      },
      Scalars -> {
         List[],
         List[
            {
               "remove scalar bosons",
               treeAll,
               FreeQ[TreeFields@##, FeynArts`S]&
            }
         ]
      }
   },
   1 -> {
      Vectors -> {
         List[
            {
               "penguins: remove external leptons from loops",
               penguinT,
               FreeQ[LoopFields@##, FieldPattern[#3, 1|3]]&
            },
            {
               "penguins: remove vector bosons from loops",
               penguinT,
               FreeQ[LoopFields@##, FeynArts`V]&
            }
         ],
         List[
            {
               "penguins: remove tree-like vector bosons",
               penguinT,
               FreeQ[TreeFields@##, FeynArts`V]&
            }
         ]
      },
      Scalars -> {
         List[
            {
               "penguins: remove external leptons from loops",
               penguinT,
               FreeQ[LoopFields@##, FieldPattern[#3, 1|3]]&
            },
            {
               "penguins: remove vector bosons from loops",
               penguinT,
               FreeQ[LoopFields@##, FeynArts`V]&
            }
         ],
         List[
            {
               "penguins: remove tree-like scalar bosons",
               penguinT,
               FreeQ[TreeFields@##, FeynArts`S]&
            }
         ]
      },
      Boxes -> {
         List[
            {
               "boxes: remove external leptons from loops",
               boxAll,
               FreeQ[LoopFields@##, FieldPattern[#3, 1|2|3|4]]&
            },
            {
               "boxes: remove vector bosons from loops",
               boxAll,
               FreeQ[LoopFields@##, FeynArts`V]&
            }
         ],
         List[]
      }
   }
};

`settings`amplitudes = {
   0 -> {
      Vectors -> {
         List[
            {
               "remove photons",
               treeAll,
               FreeQ[#, genericMass[FeynArts`V, 5] -> 0]&
            }
         ],
         List[]
      }
   },
   1 -> {
      Vectors -> {
         List[
            {
               "penguins: remove tree photons",
               penguinT,
               FreeQ[#, genericMass[FeynArts`V, 5] -> 0]&
            }
         ],
         List[]
      }
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
