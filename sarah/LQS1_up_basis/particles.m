ParticleDefinitions[EWSB] = {

{hh,   {Description-> "Higgs",
        PDG-> {25},
        PDG.IX-> {101000001}}},
{Ah,   {Description-> "Pseudo-Scalar Higgs",
        PDG-> {0},
        PDG.IX-> {0},
        Width-> {0},
        Mass -> {0}}},
{Hp,   {Description-> "Charged Higgs",
        PDG-> {0},
        PDG.IX-> {0},
        Width-> {0},
        Mass -> {0},
        ElectricCharge-> 1}},
{Fu,   {Description-> "Up-Quarks"}},
{Fd,   {Description-> "Down-Quarks"}},
{Fv,   {Description-> "Neutrinos" }},
{Fe,   {Description-> "Leptons" }},
{VP,   {Description-> "Photon"}},
{VZ,   {Description-> "Z-Boson",
        Goldstone-> Ah}},
{VG,   {Description-> "Gluon" }},
{gP,   {Description-> "Photon Ghost"}},
{gZ,   {Description-> "Z-Boson Ghost" }},
{VWp,  {Description-> "W+ - Boson",
        Goldstone-> Hp,
        ElectricCharge-> 1}},
{gWp,  {Description-> "Positive W+ - Boson Ghost"}},
{gWpC, {Description-> "Negative W+ - Boson Ghost" }},
{gG,   {Description-> "Gluon Ghost" }},
{S1c,  {Description-> "Scalar Leptoquark S1c",
        PDG-> {41},
        OutputName-> "S1c",
        ElectricCharge-> -1/3}}
};

ParticleDefinitions[GaugeES] = {

{gB,   {Description-> "B-Boson Ghost"}},
{gG,   {Description-> "Gluon Ghost" }},
{gWB,  {Description-> "W-Boson Ghost"}},
{VB,   {Description-> "B-Boson"}},
{VG,   {Description-> "Gluon"}},
{VWB,  {Description-> "W-Bosons"}},
{H0,   {PDG-> {0},
        FeynArtsNr-> 1,
        OutputName-> "H0" }},
{Hp,   {PDG-> {0},
        FeynArtsNr-> 2,
        OutputName-> "Hp" }}
};

(* We need to define at least one here. Otherwise FS complains. *)
WeylFermionAndIndermediate = {
{sc,    {LaTeX -> "S_1^*"}}
};



