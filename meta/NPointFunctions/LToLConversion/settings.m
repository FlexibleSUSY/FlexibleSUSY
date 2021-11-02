`settings`topology = {
   0 -> {
      Vectors -> (`topologyQ`tree22@#&),
      Scalars -> (`topologyQ`tree22@#&)
   },
   1 -> {
      Vectors -> (`topologyQ`penguinT@#&),
      Scalars -> (`topologyQ`penguinT@#&),
      Boxes -> (`topologyQ`box@#&)
   }
};

`settings`diagrams = {
   0 -> {
      Vectors -> {
         List[],
         List[
            {
               "remove vector bosons",
               `topologyQ`tree22,
               FreeQ[TreeFields@##, FeynArts`V]&
            }
         ]
      },
      Scalars -> {
         List[],
         List[
            {
               "remove scalar bosons",
               `topologyQ`tree22,
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
               `topologyQ`penguinT,
               FreeQ[LoopFields@##, FieldPattern[#3, 1|3]]&
            },
            {
               "penguins: remove vector bosons from loops",
               `topologyQ`penguinT,
               FreeQ[LoopFields@##, FeynArts`V]&
            }
         ],
         List[
            {
               "penguins: remove tree-like vector bosons",
               `topologyQ`penguinT,
               FreeQ[TreeFields@##, FeynArts`V]&
            }
         ]
      },
      Scalars -> {
         List[
            {
               "penguins: remove external leptons from loops",
               `topologyQ`penguinT,
               FreeQ[LoopFields@##, FieldPattern[#3, 1|3]]&
            },
            {
               "penguins: remove vector bosons from loops",
               `topologyQ`penguinT,
               FreeQ[LoopFields@##, FeynArts`V]&
            }
         ],
         List[
            {
               "penguins: remove tree-like scalar bosons",
               `topologyQ`penguinT,
               FreeQ[TreeFields@##, FeynArts`S]&
            }
         ]
      },
      Boxes -> {
         List[
            {
               "boxes: remove external leptons from loops",
               `topologyQ`box,
               FreeQ[LoopFields@##, FieldPattern[#3, 1|3]]&
            },
            {
               "boxes: remove vector bosons from loops",
               `topologyQ`box,
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
               `topologyQ`tree22,
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
               `topologyQ`penguinT,
               FreeQ[#, genericMass[FeynArts`V, 5] -> 0]&
            }
         ],
         List[]
      }
   }
};

`settings`sum =
{  {  "initial SED: skip initial lepton in sum",
      `topologyQ`self1penguinT,
      {6, Field[#3, 1]&}},
   {  "final SED: skip final lepton in sum",
      `topologyQ`self3penguinT,
      {6, Field[#3, 3]&}}};

`settings`massless =
{  {  "initial SED: use explicit final lepton mass",
      `topologyQ`self1penguinT,
      {Append, FeynArts`F[6] :> 3}},
   {  "initial SED: keep initial lepton mass untouched",
      `topologyQ`self1penguinT,
      {Hold, 1}}};

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
