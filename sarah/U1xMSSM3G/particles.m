(* SARAH generated template for particles.m file for U1xMSSM *)
(* File created at 13:02 on 31.1.2016  *) 

(* IMPORTANT: *)
(* check/adjust in particular the lines which contain "CHECK!" *)
(* the correct information is needed there to have correct results! *)


(* ###############            EWSB       ################# *) 

ParticleDefinitions[EWSB] = { 
 
(* ----------- Already defined particles in existing particles.m -------------- *) 

{Ah,{ 
     Description -> "Pseudo-Scalar Higgs", 
     FeynArtsNr -> 2, 
     LaTeX -> "A^0", 
     Mass -> {0, LesHouches}, 
     OutputName -> "Ah", 
     PDG -> {0, 0, 36, 9900036, 9900046}, 
     ElectricCharge -> 0, 
     Width -> {0, External}}}, 

{Cha,{ 
     Description -> "Charginos", 
     FeynArtsNr -> 12, 
     LaTeX -> {"\\tilde{\\chi}^-","\\tilde{\\chi}^+"}, 
     Mass -> LesHouches, 
     OutputName -> "C", 
     PDG -> {-1000024, -1000037}, 
     ElectricCharge -> -1, 
     Width -> Automatic}}, 

{Chi,{ 
     Description -> "Neutralinos", 
     FeynArtsNr -> 11, 
     LaTeX -> "\\tilde{\\chi}^0", 
     Mass -> LesHouches, 
     OutputName -> "N", 
     PDG -> {1000022, 1000023, 1000025, 1000035, 1000032, 1000036, 1000039, 1000049}, 
     ElectricCharge -> 0, 
     Width -> Automatic}}, 

{Fd,{ 
     Description -> "Down-Quarks", 
     FeynArtsNr -> 4, 
     LaTeX -> "d", 
     Mass -> {0.0035, 0.104, 4.2}, 
     OutputName -> "d", 
     PDG -> {1, 3, 5}, 
     ElectricCharge -> -1/3, 
     Width -> {0, 0, 0}}}, 

{Fe,{ 
     Description -> "Leptons", 
     FeynArtsNr -> 2, 
     LaTeX -> "e", 
     Mass -> {0.000511, 0.105, 1.776}, 
     OutputName -> "e", 
     PDG -> {11, 13, 15, 9900148,9900248,9900348}   (* adjusted number of PDGs to number of generations  *), 
     ElectricCharge -> -1, 
     Width -> {0, 0, 0}}}, 

{Fu,{ 
     Description -> "Up-Quarks", 
     FeynArtsNr -> 3, 
     LaTeX -> "u", 
     Mass -> {0.0015, 1.27, 171.2}, 
     OutputName -> "u", 
     PDG -> {2, 4, 6, 9900150, 9900250, 9900350}   (* adjusted number of PDGs to number of generations  *), 
     ElectricCharge -> 2/3, 
     Width -> {0, 0, 1.51}}}, 

{Fv,{ 
     Description -> "Neutrinos", 
     FeynArtsNr -> 1, 
     LaTeX -> "\\nu", 
     Mass -> Automatic, 
     OutputName -> "nu", 
     PDG -> {12, 14, 16, 112, 114, 116}, 
     ElectricCharge -> 0, 
     Width -> Automatic}}, 

{gG,{ 
     Description -> "Gluon Ghost", 
     FeynArtsNr -> 5, 
     LaTeX -> "\\eta^G", 
     Mass -> 0, 
     OutputName -> "gG", 
     PDG -> {0}, 
     ElectricCharge -> 0, 
     Width -> 0}}, 

{Glu,{ 
     Description -> "Gluino", 
     FeynArtsNr -> 15, 
     LaTeX -> "\\tilde{g}", 
     Mass -> LesHouches, 
     OutputName -> "go", 
     PDG -> {1000021}, 
     ElectricCharge -> 0, 
     Width -> Automatic}}, 

{gP,{ 
     Description -> "Photon Ghost", 
     FeynArtsNr -> 1, 
     LaTeX -> "\\eta^{\\gamma}", 
     Mass -> 0, 
     OutputName -> "gA", 
     PDG -> {0}   (* adjusted number of PDGs to number of generations  *), 
     ElectricCharge -> 0, 
     Width -> 0}}, 

{gWm,{ 
     Description -> "Negative W-Boson Ghost", 
     FeynArtsNr -> 3, 
     LaTeX -> "\\eta^-", 
     Mass -> Mass[VWm], 
     OutputName -> "gWm", 
     PDG -> {0}   (* adjusted number of PDGs to number of generations  *), 
     ElectricCharge -> -1, 
     Width -> Automatic}}, 

{gWmC,{ 
     Description -> "Positive W-Boson Ghost", 
     FeynArtsNr -> 4, 
     LaTeX -> "\\eta^+", 
     Mass -> Mass[VWm], 
     OutputName -> "gWpC", 
     PDG -> {0}   (* adjusted number of PDGs to number of generations  *), 
     ElectricCharge -> 1, 
     Width -> Automatic}}, 

{gZ,{ 
     Description -> "Z-Boson Ghost", 
     FeynArtsNr -> 2, 
     LaTeX -> "\\eta^Z", 
     Mass -> Mass[VZ], 
     OutputName -> "gZ", 
     PDG -> {0}   (* adjusted number of PDGs to number of generations  *), 
     ElectricCharge -> 0, 
     Width -> Automatic}}, 

{gZp,{ 
     Description -> "Z'-Ghost", 
     FeynArtsNr -> 10, 
     LaTeX -> "\\eta^{Z'}", 
     Mass -> Automatic, 
     OutputName -> "gZp", 
     PDG -> {0}   (* adjusted number of PDGs to number of generations  *), 
     ElectricCharge -> 0, 
     Width -> 0}}, 

{hh,{ 
     Description -> "Higgs", 
     FeynArtsNr -> 1, 
     LaTeX -> "h", 
     Mass -> LesHouches, 
     OutputName -> "h", 
     PDG -> {25, 35, 9900025, 9900035, 9900045}, 
     ElectricCharge -> 0, 
     Width -> Automatic}}, 

{Hpm,{ 
     Description -> "Charged Higgs", 
     FeynArtsNr -> 3, 
     LaTeX -> {"H^-","H^+"}, 
     Mass -> {0, LesHouches}, 
     OutputName -> {"Hm","Hp"}, 
     PDG -> {0, -37}, 
     ElectricCharge -> -1, 
     Width -> {0, External}}}, 

{Sd,{ 
     Description -> "Down-Squarks", 
     FeynArtsNr -> 14, 
     LaTeX -> "\\tilde{d}", 
     Mass -> LesHouches, 
     OutputName -> "sd", 
     PDG -> {1000001, 1000003, 1000005, 2000001, 2000003, 2000005}, 
     ElectricCharge -> -1/3, 
     Width -> Automatic}}, 

{Se,{ 
     Description -> "Sleptons", 
     FeynArtsNr -> 12, 
     LaTeX -> "\\tilde{e}", 
     Mass -> LesHouches, 
     OutputName -> "se", 
     PDG -> {1000011, 1000013, 1000015, 2000011, 2000013, 2000015, 9900153, 9900154, 9900253, 9900254, 9900353, 9900354}   (* adjusted number of PDGs to number of generations  *), 
     ElectricCharge -> -1, 
     Width -> Automatic}}, 

{Su,{ 
     Description -> "Up-Squarks", 
     FeynArtsNr -> 13, 
     LaTeX -> "\\tilde{u}", 
     Mass -> LesHouches, 
     OutputName -> "su", 
     PDG -> {1000002, 1000004, 1000006, 2000002, 2000004, 2000006, 9900157, 9900158, 9900257, 9900258, 9900357, 9900358}   (* adjusted number of PDGs to number of generations  *), 
     ElectricCharge -> 2/3, 
     Width -> Automatic}}, 

{SvIm,{ 
     Description -> "CP-odd Sneutrino", 
     FeynArtsNr -> 40, 
     LaTeX -> "\\nu^I", 
     Mass -> LesHouches, 
     OutputName -> "nI", 
     PDG -> {4000012, 4000014, 4000016, 5000012, 5000014, 5000016}, 
     ElectricCharge -> 0, 
     Width -> Automatic}}, 

{SvRe,{ 
     Description -> "CP-even Sneutrino", 
     FeynArtsNr -> 41, 
     LaTeX -> "\\nu^R", 
     Mass -> LesHouches, 
     OutputName -> "nR", 
     PDG -> {1000012, 1000014, 1000016, 2000012, 2000014, 2000016}, 
     ElectricCharge -> 0, 
     Width -> Automatic}}, 

{VG,{ 
     Description -> "Gluon", 
     FeynArtsNr -> 5, 
     LaTeX -> "g", 
     Mass -> 0, 
     OutputName -> "g", 
     PDG -> {21}, 
     ElectricCharge -> 0, 
     Width -> 0}}, 

{VP,{ 
     Description -> "Photon", 
     FeynArtsNr -> 1, 
     LaTeX -> "\\gamma", 
     Mass -> 0, 
     OutputName -> "A", 
     PDG -> {22}, 
     ElectricCharge -> 0, 
     Width -> 0}}, 

{VWm,{ 
     Description -> "W-Boson", 
     Goldstone -> Hpm[{1}], 
     FeynArtsNr -> 3, 
     LaTeX -> {"W^-","W^+"}, 
     Mass -> Dependence, 
     MassDependence -> Sqrt[Mass[VZ]^2/2 + Sqrt[-((Pi*Mass[VZ]^2)/(Sqrt[2]*aEWinv*Gf)) + Mass[VZ]^4/4]], 
     OutputName -> {"Wm","Wp"}, 
     PDG -> {-24}, 
     ElectricCharge -> -1, 
     Width -> 2.141}}, 

{VZ,{ 
     Description -> "Z-Boson", 
     Goldstone -> Ah[{1}], 
     FeynArtsNr -> 2, 
     LaTeX -> "Z", 
     Mass -> 91.1876, 
     OutputName -> "Z", 
     PDG -> {23}, 
     ElectricCharge -> 0, 
     Width -> 2.4952}}, 

{VZp,{ 
     Description -> "Z'-Boson", 
     Goldstone -> Ah[{2}], 
     FeynArtsNr -> 10, 
     LaTeX -> "{Z'}", 
     Mass -> LesHouches, 
     OutputName -> "Zp", 
     PDG -> {31}, 
     ElectricCharge -> 0, 
     Width -> Automatic}} 
}; 




(* ###############            GaugeES       ################# *) 

ParticleDefinitions[GaugeES] = { 
 
(* ----------- Already defined particles in existing particles.m -------------- *) 

{Bino,{ 
     Description -> "Bino", 
     FeynArtsNr -> 11, 
     LaTeX -> "\\tilde{B}", 
     Mass -> Automatic, 
     OutputName -> "b", 
     PDG -> {0}, 
     Width -> Automatic}}, 

{FB,{ 
     Description -> "Bino'", 
     FeynArtsNr -> 10, 
     LaTeX -> "{\\tilde{B}'}", 
     Mass -> LesHouches, 
     OutputName -> "FBp", 
     PDG -> {0}, 
     Width -> Automatic}}, 

{FC,{ 
     FeynArtsNr -> 20, 
     LaTeX -> "C", 
     Mass -> LesHouches, 
     OutputName -> "FC", 
     PDG -> {9910018}   (* adjusted number of PDGs to number of generations  *), 
     Width -> Automatic}}, 

{Fd1,{ 
     Description -> "Dirac Left Down-Quark", 
     FeynArtsNr -> 120, 
     LaTeX -> "d^1", 
     Mass -> LesHouches, 
     OutputName -> "d1", 
     PDG -> {9910022, 9910023, 9910024}   (* adjusted number of PDGs to number of generations  *), 
     Width -> Automatic}}, 

{Fd2,{ 
     Description -> "Dirac Right Down-Quark", 
     FeynArtsNr -> 121, 
     LaTeX -> "d^2", 
     Mass -> LesHouches, 
     OutputName -> "d2", 
     PDG -> {9910028, 9910029, 9910030}   (* adjusted number of PDGs to number of generations  *), 
     Width -> Automatic}}, 

{Fe1,{ 
     Description -> "Dirac Left Electron", 
     FeynArtsNr -> 130, 
     LaTeX -> "e^1", 
     Mass -> LesHouches, 
     OutputName -> "e1", 
     PDG -> {9910034, 9910035, 9910036}   (* adjusted number of PDGs to number of generations  *), 
     Width -> Automatic}}, 

{Fe2,{ 
     Description -> "Dirac Right Electron", 
     FeynArtsNr -> 131, 
     LaTeX -> "e^2", 
     Mass -> LesHouches, 
     OutputName -> "e2", 
     PDG -> {9910040, 9910041, 9910042}   (* adjusted number of PDGs to number of generations  *), 
     Width -> Automatic}}, 

{Fu1,{ 
     Description -> "Dirac Left Up-Quark", 
     FeynArtsNr -> 110, 
     LaTeX -> "u^1", 
     Mass -> LesHouches, 
     OutputName -> "u1", 
     PDG -> {9910046, 9910047, 9910048}   (* adjusted number of PDGs to number of generations  *), 
     Width -> Automatic}}, 

{Fu2,{ 
     Description -> "Dirac Right Up-Quark", 
     FeynArtsNr -> 111, 
     LaTeX -> "u^2", 
     Mass -> LesHouches, 
     OutputName -> "u2", 
     PDG -> {9910052, 9910053, 9910054}   (* adjusted number of PDGs to number of generations  *), 
     Width -> Automatic}}, 

{Fv1,{ 
     Description -> "Dirac Left Neutrino", 
     FeynArtsNr -> 140, 
     LaTeX -> "\\nu^1", 
     Mass -> LesHouches, 
     OutputName -> "v1", 
     PDG -> {9910058, 9910059, 9910060}   (* adjusted number of PDGs to number of generations  *), 
     Width -> Automatic}}, 

{Fv2,{ 
     Description -> "Dirac Right Neutrino", 
     FeynArtsNr -> 141, 
     LaTeX -> "\\nu^2", 
     Mass -> LesHouches, 
     OutputName -> "v2", 
     PDG -> {9910064, 9910065, 9910066}   (* adjusted number of PDGs to number of generations  *), 
     Width -> Automatic}}, 

{gB,{ 
     Description -> "B-Boson Ghost", 
     FeynArtsNr -> 1, 
     LaTeX -> "\\eta^B", 
     Mass -> 0, 
     OutputName -> "gB", 
     PDG -> {0}, 
     Width -> 0}}, 

{gBp,{ 
     Description -> "B'-Boson Ghost", 
     FeynArtsNr -> 10, 
     LaTeX -> "{g_{B'}}", 
     Mass -> LesHouches, 
     OutputName -> "gBp", 
     PDG -> {0}, 
     Width -> Automatic}}, 

{gG,{ 
     Description -> "Gluon Ghost", 
     FeynArtsNr -> 5, 
     LaTeX -> "\\eta^G", 
     Mass -> 0, 
     OutputName -> "gG", 
     PDG -> {0}, 
     ElectricCharge -> 0, 
     Width -> 0}}, 

{Glu,{ 
     Description -> "Gluino", 
     FeynArtsNr -> 15, 
     LaTeX -> "\\tilde{g}", 
     Mass -> LesHouches, 
     OutputName -> "go", 
     PDG -> {1000021}, 
     ElectricCharge -> 0, 
     Width -> Automatic}}, 

{gWB,{ 
     Description -> "W-Boson Ghost", 
     FeynArtsNr -> 2, 
     LaTeX -> "\\eta^W", 
     Mass -> 0, 
     OutputName -> "gW", 
     PDG -> {0}   (* adjusted number of PDGs to number of generations  *), 
     Width -> 0}}, 

{H0,{ 
     Description -> "Neutral Higgsinos", 
     FeynArtsNr -> 21, 
     LaTeX -> "\\tilde{H}^0", 
     Mass -> Automatic, 
     OutputName -> "H0", 
     PDG -> {0}, 
     Width -> Automatic}}, 

{HC,{ 
     Description -> "Charged Higgsinos", 
     FeynArtsNr -> 22, 
     LaTeX -> "\\tilde{H}^-", 
     Mass -> Automatic, 
     OutputName -> "HC", 
     PDG -> {0}, 
     Width -> Automatic}}, 

{SC10,{ 
     FeynArtsNr -> 102, 
     LaTeX -> "\\eta", 
     Mass -> Automatic, 
     OutputName -> "C1", 
     PDG -> {0}, 
     Width -> 0}}, 

{SC20,{ 
     FeynArtsNr -> 103, 
     LaTeX -> "\\bar{\\eta}", 
     Mass -> Automatic, 
     OutputName -> "C2", 
     PDG -> {0}, 
     Width -> 0}}, 

{SdL,{ 
     Description -> "Left Down-Squarks", 
     FeynArtsNr -> 14, 
     LaTeX -> "\\tilde{d}_L^0", 
     Mass -> Automatic, 
     OutputName -> "dL", 
     PDG -> {0, 0, 0}, 
     Width -> Automatic}}, 

{SdR,{ 
     Description -> "Right Down-Squarks", 
     FeynArtsNr -> 15, 
     LaTeX -> "\\tilde{d}_R^0", 
     Mass -> Automatic, 
     OutputName -> "dR", 
     PDG -> {0, 0, 0}, 
     Width -> Automatic}}, 

{SeL,{ 
     Description -> "Left Selectron", 
     FeynArtsNr -> 12, 
     LaTeX -> "\\tilde{e}_L", 
     Mass -> Automatic, 
     OutputName -> "eL", 
     PDG -> {0, 0, 0}, 
     Width -> Automatic}}, 

{SeR,{ 
     Description -> "Right Selectron", 
     FeynArtsNr -> 13, 
     LaTeX -> "\\tilde{e}_R", 
     Mass -> Automatic, 
     OutputName -> "eR", 
     PDG -> {0, 0, 0}, 
     Width -> Automatic}}, 

{SHd0,{ 
     Description -> "Neutral Down-Higgs", 
     FeynArtsNr -> 1, 
     LaTeX -> "H_d^0", 
     Mass -> Automatic, 
     OutputName -> "Hd0", 
     PDG -> {0}, 
     Width -> Automatic}}, 

{SHdm,{ 
     Description -> "Charged Down-Higgs", 
     FeynArtsNr -> 2, 
     LaTeX -> "H_d^-", 
     Mass -> Automatic, 
     OutputName -> "Hdm", 
     PDG -> {0}, 
     Width -> Automatic}}, 

{SHu0,{ 
     Description -> "Neutral Up-Higgs", 
     FeynArtsNr -> 3, 
     LaTeX -> "H_u^0", 
     Mass -> Automatic, 
     OutputName -> "Hu0", 
     PDG -> {0}, 
     Width -> Automatic}}, 

{SHup,{ 
     Description -> "Charged Up-Higgs", 
     FeynArtsNr -> 4, 
     LaTeX -> "H_u^+", 
     Mass -> Automatic, 
     OutputName -> "Hup", 
     PDG -> {0}, 
     Width -> Automatic}}, 

{SuL,{ 
     Description -> "Left Up-Squarks", 
     FeynArtsNr -> 16, 
     LaTeX -> "\\tilde{u}_L^0", 
     Mass -> Automatic, 
     OutputName -> "uL", 
     PDG -> {0, 0, 0}, 
     Width -> Automatic}}, 

{SuR,{ 
     Description -> "Right Up-Squarks", 
     FeynArtsNr -> 17, 
     LaTeX -> "\\tilde{u}_R^0", 
     Mass -> Automatic, 
     OutputName -> "uR", 
     PDG -> {0, 0, 0}, 
     Width -> Automatic}}, 

{SvL,{ 
     Description -> "Left Sneutrino", 
     FeynArtsNr -> 10, 
     LaTeX -> "\\tilde{\\nu}_L", 
     Mass -> Automatic, 
     OutputName -> "nL", 
     PDG -> {0, 0, 0}, 
     Width -> Automatic}}, 

{SvR,{ 
     Description -> "Right Sneutrino", 
     FeynArtsNr -> 18, 
     LaTeX -> "\\tilde{\\nu}_R", 
     Mass -> Automatic, 
     OutputName -> "vR", 
     PDG -> {0, 0, 0}, 
     Width -> 0}}, 

{VB,{ 
     Description -> "B-Boson", 
     FeynArtsNr -> 1, 
     LaTeX -> "B", 
     Mass -> 0, 
     OutputName -> "B", 
     PDG -> {0}, 
     Width -> 0}}, 

{VBp,{ 
     Description -> "B'-Boson", 
     FeynArtsNr -> 10, 
     LaTeX -> "{B'}", 
     Mass -> LesHouches, 
     OutputName -> "Bp", 
     PDG -> {0}, 
     Width -> Automatic}}, 

{VG,{ 
     Description -> "Gluon", 
     FeynArtsNr -> 5, 
     LaTeX -> "g", 
     Mass -> 0, 
     OutputName -> "g", 
     PDG -> {21}, 
     ElectricCharge -> 0, 
     Width -> 0}}, 

{VWB,{ 
     Description -> "W-Bosons", 
     FeynArtsNr -> 2, 
     LaTeX -> "W", 
     Mass -> 0, 
     OutputName -> "W", 
     PDG -> {0, 9910069, 9910070}   (* adjusted number of PDGs to number of generations  *), 
     Width -> 0}}, 

{Wino,{ 
     Description -> "Wino", 
     FeynArtsNr -> 12, 
     LaTeX -> "\\tilde{W}", 
     Mass -> Automatic, 
     OutputName -> "w", 
     PDG -> {0, 9910073, 9910074}   (* adjusted number of PDGs to number of generations  *), 
     Width -> Automatic}}, 



(*  ----------- Additional particles in the model ---------- *) 

{FeRx,{ 
     FeynArtsNr -> 142   (* auto generated FeynArts number *), 
     LaTeX -> "eRx"   (* auto generated LaTeX name *), 
     Mass -> LesHouches, 
     OutputName -> "ferx"   (* auto generated Output name *), 
     PDG -> {9910076}   (* auto generated PDGs  *), 
     Width -> Automatic}}, 

{FeRxp,{ 
     FeynArtsNr -> 144   (* auto generated FeynArts number *), 
     LaTeX -> "eRxp"   (* auto generated LaTeX name *), 
     Mass -> LesHouches, 
     OutputName -> "ferxp"   (* auto generated Output name *), 
     PDG -> {9910078}   (* auto generated PDGs  *), 
     Width -> Automatic}}, 

{Fsing,{ 
     FeynArtsNr -> 146   (* auto generated FeynArts number *), 
     LaTeX -> "sing"   (* auto generated LaTeX name *), 
     Mass -> LesHouches, 
     OutputName -> "fsing"   (* auto generated Output name *), 
     PDG -> {9910080}   (* auto generated PDGs  *), 
     Width -> Automatic}}, 

{FuRx,{ 
     FeynArtsNr -> 148   (* auto generated FeynArts number *), 
     LaTeX -> "uRx"   (* auto generated LaTeX name *), 
     Mass -> LesHouches, 
     OutputName -> "furx"   (* auto generated Output name *), 
     PDG -> {9910082}   (* auto generated PDGs  *), 
     Width -> Automatic}}, 

{FuRxp,{ 
     FeynArtsNr -> 150   (* auto generated FeynArts number *), 
     LaTeX -> "uRxp"   (* auto generated LaTeX name *), 
     Mass -> LesHouches, 
     OutputName -> "furxp"   (* auto generated Output name *), 
     PDG -> {9910084}   (* auto generated PDGs  *), 
     Width -> Automatic}}, 

{SeRx,{ 
     FeynArtsNr -> 152   (* auto generated FeynArts number *), 
     LaTeX -> "\\tilde{eRx}"   (* auto generated LaTeX name *), 
     Mass -> LesHouches, 
     OutputName -> "serx"   (* auto generated Output name *), 
     PDG -> {9910086}   (* auto generated PDGs  *), 
     Width -> Automatic}}, 

{SeRxp,{ 
     FeynArtsNr -> 154   (* auto generated FeynArts number *), 
     LaTeX -> "\\tilde{eRxp}"   (* auto generated LaTeX name *), 
     Mass -> LesHouches, 
     OutputName -> "serxp"   (* auto generated Output name *), 
     PDG -> {9910088}   (* auto generated PDGs  *), 
     Width -> Automatic}}, 

{Ssing,{ 
     FeynArtsNr -> 156   (* auto generated FeynArts number *), 
     LaTeX -> "\\tilde{sing}"   (* auto generated LaTeX name *), 
     Mass -> LesHouches, 
     OutputName -> "ssing"   (* auto generated Output name *), 
     PDG -> {9910090}   (* auto generated PDGs  *), 
     Width -> Automatic}}, 

{SuRx,{ 
     FeynArtsNr -> 158   (* auto generated FeynArts number *), 
     LaTeX -> "\\tilde{uRx}"   (* auto generated LaTeX name *), 
     Mass -> LesHouches, 
     OutputName -> "surx"   (* auto generated Output name *), 
     PDG -> {9910092}   (* auto generated PDGs  *), 
     Width -> Automatic}}, 

{SuRxp,{ 
     FeynArtsNr -> 160   (* auto generated FeynArts number *), 
     LaTeX -> "\\tilde{uRxp}"   (* auto generated LaTeX name *), 
     Mass -> LesHouches, 
     OutputName -> "surxp"   (* auto generated Output name *), 
     PDG -> {9910094}   (* auto generated PDGs  *), 
     Width -> Automatic}} 
}; 




(* ###############            Weyl Spinors and intermediate states       ################# *) 

WeylFermionAndIndermediate = { 
 
(* ----------- Already defined particles in existing particles.m -------------- *) 

{fB,{ 
     Description -> "Bino Weyl-Spinor", 
     LaTeX -> "\\lambda_{\\tilde{B}}"}}, 

{fWB,{ 
     Description -> "Wino Weyl-Spinor", 
     LaTeX -> "\\lambda_{\\tilde{W}}"}}, 

{fG,{ 
     Description -> "Gluino Weyl-Spinor", 
     LaTeX -> "\\lambda_{\\tilde{g}}"}}, 

{fBp,{ 
     LaTeX -> "\\lambda_{B'}"}}, 

{FdL,{ 
     Description -> "Left Down-Quark", 
     LaTeX -> "d_L"}}, 

{FuL,{ 
     Description -> "Left Up-Quark", 
     LaTeX -> "u_L"}}, 

{FeL,{ 
     Description -> "Left Electron", 
     LaTeX -> "e_L"}}, 

{FvL,{ 
     Description -> "Left Neutrino", 
     LaTeX -> "\\nu_L"}}, 

{FHd0,{ 
     Description -> "Neutral Down-Higgsino", 
     LaTeX -> "\\tilde{H}_d^0"}}, 

{FHdm,{ 
     Description -> "Charged Down-Higgsino", 
     LaTeX -> "\\tilde{H}_d^-"}}, 

{FHu0,{ 
     Description -> "Neutral Up-Higgsino", 
     LaTeX -> "\\tilde{H}_u^0"}}, 

{FHup,{ 
     Description -> "Charged Up-Higgsino", 
     LaTeX -> "\\tilde{H}_u^+"}}, 

{FC10,{ 
     LaTeX -> "\\tilde{\\eta}"}}, 

{FC20,{ 
     LaTeX -> "\\tilde{\\bar{\\eta}}"}}, 

{fWm,{ 
     Description -> "Negative Wino", 
     LaTeX -> "\\tilde{W}^-"}}, 

{fWp,{ 
     Description -> "Positive Wino", 
     LaTeX -> "\\tilde{W}^+"}}, 

{fW0,{ 
     Description -> "Neutral Wino", 
     LaTeX -> "\\tilde{W}^0"}}, 

{L0,{ 
     Description -> "Neutralino Weyl-Spinor", 
     LaTeX -> "\\lambda^0"}}, 

{Lm,{ 
     Description -> "Negative Chargino Weyl-Spinor", 
     LaTeX -> "\\lambda^-"}}, 

{Lp,{ 
     Description -> "Positive Chargino Weyl-Spinor", 
     LaTeX -> "\\lambda^+"}}, 

{Fvm,{ 
     LaTeX -> "\\lambda_\\nu"}}, 

{FEL,{ 
     Description -> "Rotated Left Electron", 
     LaTeX -> "E_L"}}, 

{FER,{ 
     Description -> "Rotated Right Electron", 
     LaTeX -> "E_R"}}, 

{FDL,{ 
     Description -> "Rotated Left Up-Quark", 
     LaTeX -> "D_L"}}, 

{FDR,{ 
     Description -> "Rotated Right Up-Quark", 
     LaTeX -> "D_R"}}, 

{FUL,{ 
     Description -> "Rotated Left Down-Quark", 
     LaTeX -> "U_L"}}, 

{FUR,{ 
     Description -> "Rotated Right Down-Quark", 
     LaTeX -> "U_R"}}, 

{sigmad,{ 
     Description -> "Pseudo Scalar Down", 
     LaTeX -> "\\sigma_{d}"}}, 

{phid,{ 
     Description -> "Scalar Down", 
     LaTeX -> "\\phi_{d}"}}, 

{sigmau,{ 
     Description -> "Pseudo Scalar Up", 
     LaTeX -> "\\sigma_{u}"}}, 

{phiu,{ 
     Description -> "Scalar Up", 
     LaTeX -> "\\phi_{u}"}}, 

{sigma1,{ 
     LaTeX -> "\\sigma_{\\eta}"}}, 

{phi1,{ 
     LaTeX -> "\\phi_{\\eta}"}}, 

{sigma2,{ 
     LaTeX -> "\\sigma_{\\bar{\\eta}}"}}, 

{phi2,{ 
     LaTeX -> "\\phi_{\\bar{\\eta}}"}}, 

{q,{ 
     Description -> "Left Quark Superfield", 
     LaTeX -> "\\hat{q}"}}, 

{l,{ 
     Description -> "left Lepton Superfield", 
     LaTeX -> "\\hat{l}"}}, 

{Hd,{ 
     Description -> "Down-Higgs Superfield", 
     LaTeX -> "\\hat{H}_d"}}, 

{Hu,{ 
     Description -> "Up-Higgs Superfield", 
     LaTeX -> "\\hat{H}_u"}}, 

{d,{ 
     Description -> "Right Down-Quark Superfield", 
     LaTeX -> "\\hat{d}"}}, 

{u,{ 
     Description -> "Right Up-Quark Superfield", 
     LaTeX -> "\\hat{u}"}}, 

{e,{ 
     Description -> "Right Electron Superfield", 
     LaTeX -> "\\hat{e}"}}, 

{vR,{ 
     Description -> "Right Neutrino Superfield", 
     LaTeX -> "\\hat{\\nu}"}}, 

{C1,{ 
     LaTeX -> "\\hat{\\eta}"}}, 

{C2,{ 
     LaTeX -> "\\hat{\\bar{\\eta}}"}}, 

{B,{ 
     Description -> "B Superfield", 
     LaTeX -> "\\hat{B}"}}, 

{WB,{ 
     Description -> "W Superfield", 
     LaTeX -> "\\hat{W}"}}, 

{G,{ 
     Description -> "Gluon Superfield", 
     LaTeX -> "\\hat{g}"}}, 

{Bp,{ 
     LaTeX -> "\\hat{B}'"}}, 

{Fq,{ 
     Description -> "Left Quarks", 
     LaTeX -> "q"}}, 

{Fl,{ 
     Description -> "Left Leptons", 
     LaTeX -> "l"}}, 

{FHd,{ 
     Description -> "Down-Higgsino", 
     LaTeX -> "\\tilde{H}_d"}}, 

{FHu,{ 
     Description -> "Up-Higgsino", 
     LaTeX -> "\\tilde{H}_u"}}, 

{FC10,{ 
     LaTeX -> "\\tilde{\\eta}"}}, 

{FC20,{ 
     LaTeX -> "\\tilde{\\bar{\\eta}}"}}, 

{fB,{ 
     Description -> "Bino Weyl-Spinor", 
     LaTeX -> "\\lambda_{\\tilde{B}}"}}, 

{fWB,{ 
     Description -> "Wino Weyl-Spinor", 
     LaTeX -> "\\lambda_{\\tilde{W}}"}}, 

{fG,{ 
     Description -> "Gluino Weyl-Spinor", 
     LaTeX -> "\\lambda_{\\tilde{g}}"}}, 

{fBp,{ 
     LaTeX -> "\\lambda_{B'}"}}, 

{Sq,{ 
     Description -> "Left Squark", 
     LaTeX -> "\\tilde{q}"}}, 

{Sl,{ 
     Description -> "Left Slepton", 
     LaTeX -> "\\tilde{l}"}}, 

{SHd,{ 
     Description -> "Down-Higgs", 
     LaTeX -> "H_d"}}, 

{SHu,{ 
     Description -> "Up-Higgs", 
     LaTeX -> "H_u"}}, 



(*  ----------- Additional particles in the model ---------- *) 

{FvR,{ 
     LaTeX -> "{v}_R"   (* auto generated LaTeX name *)}}, 

{Fsing,{ 
     LaTeX -> "sing"   (* auto generated LaTeX name *)}}, 

{FuRx,{ 
     LaTeX -> "uRx"   (* auto generated LaTeX name *)}}, 

{FuRxp,{ 
     LaTeX -> "uRxp"   (* auto generated LaTeX name *)}}, 

{FeRx,{ 
     LaTeX -> "eRx"   (* auto generated LaTeX name *)}}, 

{FeRxp,{ 
     LaTeX -> "eRxp"   (* auto generated LaTeX name *)}}, 

{sigmaS,{ 
     LaTeX -> "{\\sigma}_{s}"   (* auto generated LaTeX name *)}}, 

{phiS,{ 
     LaTeX -> "{\\phi}_{s}"   (* auto generated LaTeX name *)}}, 

{sigmaL,{ 
     LaTeX -> "{\\sigma}_{l}"   (* auto generated LaTeX name *)}}, 

{phiL,{ 
     LaTeX -> "{\\phi}_{l}"   (* auto generated LaTeX name *)}}, 

{sigmaR,{ 
     LaTeX -> "{\\sigma}_{r}"   (* auto generated LaTeX name *)}}, 

{phiR,{ 
     LaTeX -> "{\\phi}_{r}"   (* auto generated LaTeX name *)}}, 

{S,{ 
     LaTeX -> "\\hat{S}"   (* auto generated LaTeX name *)}}, 

{UX,{ 
     LaTeX -> "\\hat{UX}"   (* auto generated LaTeX name *)}}, 

{UXp,{ 
     LaTeX -> "\\hat{UXp}"   (* auto generated LaTeX name *)}}, 

{EX,{ 
     LaTeX -> "\\hat{EX}"   (* auto generated LaTeX name *)}}, 

{EXp,{ 
     LaTeX -> "\\hat{EXp}"   (* auto generated LaTeX name *)}}, 

{conj[FdR],{ 
     LaTeX -> "conjfdr"   (* auto generated LaTeX name *)}}, 

{conj[FuR],{ 
     LaTeX -> "conjfur"   (* auto generated LaTeX name *)}}, 

{conj[FeR],{ 
     LaTeX -> "conjfer"   (* auto generated LaTeX name *)}}, 

{conj[FvR],{ 
     LaTeX -> "conjfvr"   (* auto generated LaTeX name *)}}, 

{Fsing,{ 
     LaTeX -> "sing"   (* auto generated LaTeX name *)}}, 

{conj[FuRx],{ 
     LaTeX -> "conjfurx"   (* auto generated LaTeX name *)}}, 

{FuRxp,{ 
     LaTeX -> "uRxp"   (* auto generated LaTeX name *)}}, 

{conj[FeRx],{ 
     LaTeX -> "conjferx"   (* auto generated LaTeX name *)}}, 

{FeRxp,{ 
     LaTeX -> "eRxp"   (* auto generated LaTeX name *)}}, 

{VB,{ 
     LaTeX -> "B"   (* auto generated LaTeX name *)}}, 

{VWB,{ 
     LaTeX -> "WB"   (* auto generated LaTeX name *)}}, 

{VG,{ 
     LaTeX -> "G"   (* auto generated LaTeX name *)}}, 

{VBp,{ 
     LaTeX -> "Bp"   (* auto generated LaTeX name *)}}, 

{conj[SdR],{ 
     LaTeX -> "conjsdr"   (* auto generated LaTeX name *)}}, 

{conj[SuR],{ 
     LaTeX -> "conjsur"   (* auto generated LaTeX name *)}}, 

{conj[SeR],{ 
     LaTeX -> "conjser"   (* auto generated LaTeX name *)}}, 

{conj[SvR],{ 
     LaTeX -> "conjsvr"   (* auto generated LaTeX name *)}}, 

{SC10,{ 
     LaTeX -> "\\tilde{{C1}^0}"   (* auto generated LaTeX name *)}}, 

{SC20,{ 
     LaTeX -> "\\tilde{{C2}^0}"   (* auto generated LaTeX name *)}}, 

{Ssing,{ 
     LaTeX -> "\\tilde{sing}"   (* auto generated LaTeX name *)}}, 

{conj[SuRx],{ 
     LaTeX -> "conjsurx"   (* auto generated LaTeX name *)}}, 

{SuRxp,{ 
     LaTeX -> "\\tilde{uRxp}"   (* auto generated LaTeX name *)}}, 

{conj[SeRx],{ 
     LaTeX -> "conjserx"   (* auto generated LaTeX name *)}}, 

{SeRxp,{ 
     LaTeX -> "\\tilde{eRxp}"   (* auto generated LaTeX name *)}} 
}; 




