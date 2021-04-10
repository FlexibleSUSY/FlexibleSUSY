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
      {  {  {  "boxes: remove external leptons from loops",
               `topologyQ`box,
               FreeQ[LoopFields@##, fieldPattern[#3, 1|2|3|4]]&},
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

`settings`sum :=
{  {  "SED T: skip initial lepton in sum",
      `topologyQ`self1penguinT,
      {6, Field[#3, 1]&}},
   {  "SED T: skip final lepton in sum",
      `topologyQ`self3penguinT,
      {6, Field[#3, 3]&}}};

`settings`massless =
{  {  "SED T: use explicit final lepton mass",
      `topologyQ`self1penguinT,
      {Append, FeynArts`F[6] :> 3}},
   {  "SED T: keep initial lepton mass untouched",
      `topologyQ`self1penguinT,
      {Hold, 1}}};

`settings`momenta =
{  `topologyQ`penguinT -> 2,
   `topologyQ`box -> 2};

`settings`regularization =
{  `topologyQ`box -> D};

`settings`order = {3, 1, 4, 2};

`settings`chains =
{  {  ExceptLoops} ->
      {  1[k[4|2], ___] -> 0, 2[k[3|1], ___] -> 0}};
