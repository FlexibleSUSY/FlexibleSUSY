(* SARAH generated template for parameters.m file for U1xMSSM*)
(* File created at 13:02 on 31.1.2016  *) 
(* IMPORTANT: *)
(* check/adjust in particular the lines which contain "CHECK!" *)
(* the correct information is needed there to have correct results! *)


ParameterDefinitions = { 
 
(* ----------- Already defined parameters in existing parameters.m ----------- *) 

{aEWinv,{ 
     Description -> "inverse weak coupling constant at mZ", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> True, 
     Value -> 137.035999679, 
     LesHouches -> {SMINPUTS, 1}, 
     LaTeX -> "\\alpha^{-1}", 
     OutputName -> aEWinv}}, 

{AlphaS,{ 
     Description -> "Alpha Strong", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> True, 
     Value -> 0.119, 
     LesHouches -> {SMINPUTS, 3}, 
     LaTeX -> "\\alpha_S", 
     OutputName -> aS}}, 


{e,{ 
     Description -> "electric charge", 
     Dependence -> None, 
     DependenceNum -> 2*Sqrt[aEWinv^(-1)]*Sqrt[Pi], 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> True, 
     LaTeX -> "e", 
     OutputName -> el}}, 

{g1,{ 
     Description -> "Hypercharge-Coupling", 
     Dependence -> None, 
     DependenceNum -> e*Sec[ThetaW], 
     DependenceOptional -> e*Sec[ThetaW], 
     DependenceSPheno -> None, 
     GUTnormalization -> Sqrt[3/5], 
     Real -> True, 
     LesHouches -> {gauge, 1}, 
     LaTeX -> "g_1", 
     OutputName -> g1}}, 

{g1X,{ 
     Description -> "Mixed Gauge Coupling 1", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> True, 
     LesHouches -> {gauge, 10}, 
     LaTeX -> "g_{Y X}", 
     OutputName -> gYX}}, 

{g2,{ 
     Description -> "Left-Coupling", 
     Dependence -> None, 
     DependenceNum -> e*Csc[ThetaW], 
     DependenceOptional -> e*Csc[ThetaW], 
     DependenceSPheno -> None, 
     Real -> True, 
     LesHouches -> {gauge, 2}, 
     LaTeX -> "g_2", 
     OutputName -> g2}}, 

{g3,{ 
     Description -> "Strong-Coupling", 
     Dependence -> None, 
     DependenceNum -> 2*Sqrt[AlphaS]*Sqrt[Pi], 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> True, 
     LesHouches -> {gauge, 3}, 
     LaTeX -> "g_3", 
     OutputName -> g3}}, 

{gX,{ 
     Description -> "B-L-Coupling", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     GUTnormalization -> 1, 
     Real -> True, 
     LesHouches -> {gauge, 4}, 
     LaTeX -> "g_{X}", 
     OutputName -> gX}}, 

{gX1,{ 
     Description -> "Mixed Gauge Coupling 2", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> True, 
     LesHouches -> {gauge, 11}, 
     LaTeX -> "g_{X Y}", 
     OutputName -> gXY}}, 

{Gf,{ 
     Description -> "Fermi's constant", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> True, 
     Value -> 0.0000116639, 
     LesHouches -> {SMINPUTS, 2}, 
     LaTeX -> "G_f", 
     OutputName -> Gf}}, 

{lw,{ 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> False, 
     LesHouches -> {SING, 13}, 
     LaTeX -> "lw"  (* auto generated LaTeX name *), 
     OutputName -> lw}}, 

{MassB,{ 
     Description -> "Bino Mass parameter", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Form -> Scalar, 
     Real -> False, 
     LesHouches -> {MSOFT, 1}, 
     LaTeX -> "M_1", 
     OutputName -> M1}}, 

{MassBBX,{ 
     Description -> "Mixed Gaugino Mass 1", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> False, 
     LesHouches -> {XSECTOR, 32}, 
     LaTeX -> "{M}_{B B'}", 
     OutputName -> MBBp}}, 

{MassBX,{ 
     Description -> "Bino' Mass", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> False, 
     LesHouches -> {XSECTOR, 31}, 
     LaTeX -> "{M}_{BL}", 
     OutputName -> MBp}}, 

{MassG,{ 
     Description -> "Gluino Mass parameter", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Form -> Scalar, 
     Real -> False, 
     LesHouches -> {MSOFT, 3}, 
     LaTeX -> "M_3", 
     OutputName -> M3}}, 

{MassWB,{ 
     Description -> "Wino Mass parameter", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Form -> Scalar, 
     Real -> False, 
     LesHouches -> {MSOFT, 2}, 
     LaTeX -> "M_2", 
     OutputName -> M2}}, 

{mC12,{ 
     Description -> "Bilepton 1 Soft-Breaking mass", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> True, 
     LesHouches -> {XSECTOR, 11}, 
     LaTeX -> "m_{\\eta}^2", 
     OutputName -> mC12}}, 

{mC22,{ 
     Description -> "Bilepton 2 Soft-Breaking mass", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> True, 
     LesHouches -> {XSECTOR, 12}, 
     LaTeX -> "m_{\\bar{\\eta}}^2", 
     OutputName -> mC22}}, 

{md2,{ 
     Description -> "Softbreaking right Down-Squark Mass", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> False, 
     LesHouches -> MSD2, 
     LaTeX -> "m_d^2", 
     OutputName -> md2}}, 

{me2,{ 
     Description -> "Softbreaking right Slepton Mass", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> False, 
     LesHouches -> MSE2, 
     LaTeX -> "m_e^2", 
     OutputName -> me2}}, 

{mHd2,{ 
     Description -> "Softbreaking Down-Higgs Mass", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> True, 
     LesHouches -> {MSOFT, 21}, 
     LaTeX -> "m_{H_d}^2", 
     OutputName -> mHd2}}, 

{mHu2,{ 
     Description -> "Softbreaking Up-Higgs Mass", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> True, 
     LesHouches -> {MSOFT, 22}, 
     LaTeX -> "m_{H_u}^2", 
     OutputName -> mHu2}}, 

{ml2,{ 
     Description -> "Softbreaking left Slepton Mass", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> False, 
     LesHouches -> MSL2, 
     LaTeX -> "m_l^2", 
     OutputName -> ml2}}, 

{mq2,{ 
     Description -> "Softbreaking left Squark Mass", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> False, 
     LesHouches -> MSQ2, 
     LaTeX -> "m_q^2", 
     OutputName -> mq2}}, 

{mS2,{ 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> True, 
     LesHouches -> {SING, 14}, 
     LaTeX -> "m^2_{S}"  (* auto generated LaTeX name *), 
     OutputName -> ms2}}, 

{mu2,{ 
     Description -> "Softbreaking right Up-Squark Mass", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> False, 
     LesHouches -> MSU2, 
     LaTeX -> "m_u^2", 
     OutputName -> mu2}}, 

{mvR2,{ 
     Description -> "Softbreaking right Sneutrino Mass", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> False, 
     LesHouches -> mv2, 
     LaTeX -> "m_{\\nu}^2", 
     OutputName -> mv2}}, 

{PhaseGlu,{ 
     Description -> "Gluino-Phase", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> False, 
     LesHouches -> {Phases, 1}, 
     LaTeX -> "\\phi_{\\tilde{g}}", 
     OutputName -> pG}}, 

{TBetaP,{ 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> True, 
     LesHouches -> {XSECTOR, 80}, 
     LaTeX -> "\\tan(\\beta')", 
     OutputName -> TBp}}, 

{ThetaW,{ 
     Description -> "Weinberg-Angle", 
     Dependence -> None, 
     DependenceNum -> ArcSin[Sqrt[1 - Mass[VWm]^2/Mass[VZ]^2]], 
     DependenceOptional -> None, 
     DependenceSPheno -> ArcCos[Abs[ZZ[1, 1]]], 
     Real -> True, 
     LaTeX -> "\\Theta_W", 
     OutputName -> TW}}, 

{ThetaWp,{ 
     Description -> "Theta'", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> ArcCos[Abs[ZZ[3, 3]]], 
     Real -> True, 
     LesHouches -> {ANGLES, 10}, 
     LaTeX -> "{\\Theta'}_W", 
     OutputName -> TWp}}, 

{UM,{ 
     Description -> "Chargino-minus Mixing-Matrix", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> False, 
     LesHouches -> UMIX, 
     LaTeX -> "U", 
     OutputName -> UM}}, 

{UP,{ 
     Description -> "Chargino-plus Mixing-Matrix", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> False, 
     LesHouches -> VMIX, 
     LaTeX -> "V", 
     OutputName -> UP}}, 

{UV,{ 
     Description -> "Neutrino-Mixing-Matrix", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> False, 
     LesHouches -> UVMIX, 
     LaTeX -> "U^V", 
     OutputName -> UV}}, 

{v,{ 
     Description -> "EW-VEV", 
     Dependence -> None, 
     DependenceNum -> 2*Sqrt[Mass[VWm]^2/g2^2], 
     DependenceOptional -> None, 
     DependenceSPheno -> Sqrt[vd^2 + vu^2], 
     Real -> True, 
     LesHouches -> {HMIX, 3}, 
     LaTeX -> "v", 
     OutputName -> v}}, 

{vd,{ 
     Description -> "Down-VEV", 
     Dependence -> None, 
     DependenceNum -> v*Cos[\[Beta]], 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> True, 
     LesHouches -> {HMIX, 102}, 
     LaTeX -> "v_d", 
     OutputName -> vd}}, 

{vu,{ 
     Description -> "Up-VEV", 
     Dependence -> None, 
     DependenceNum -> v*Sin[\[Beta]], 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> True, 
     LesHouches -> {HMIX, 103}, 
     LaTeX -> "v_u", 
     OutputName -> vu}}, 

{x1,{ 
     Description -> "Bilepton 1 VEV", 
     Dependence -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     DependenceNum -> None,      
     Real -> True, 
     LesHouches -> {XSECTOR, 41}, 
     LaTeX -> "v_{\\eta}", 
     OutputName -> x1}}, 

{x2,{ 
     Description -> "Bilepton 2 VEV", 
     Dependence -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     DependenceNum -> None,           
     Real -> True, 
     LesHouches -> {XSECTOR, 42}, 
     LaTeX -> "v_{\\bar{\\eta}}", 
     OutputName -> x2}}, 

{xS,{ 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> True, 
     LesHouches -> {SING, 15}, 
     LaTeX -> "xS"  (* auto generated LaTeX name *), 
     OutputName -> xS}}, 

{Yd,{ 
     Description -> "Down-Yukawa-Coupling", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> False, 
     LesHouches -> Yd, 
     LaTeX -> "Y_d", 
     OutputName -> Yd}}, 

{Ye,{ 
     Description -> "Lepton-Yukawa-Coupling", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> False, 
     LesHouches -> Ye, 
     LaTeX -> "Y_e", 
     OutputName -> Ye}}, 

{Yn,{ 
     Description -> "Neutrino-X-Yukawa-Coupling", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> False, 
     LesHouches -> Yx, 
     LaTeX -> "Y_x", 
     OutputName -> Yx}}, 

{Yu,{ 
     Description -> "Up-Yukawa-Coupling", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> False, 
     LesHouches -> Yu, 
     LaTeX -> "Y_u", 
     OutputName -> Yu}}, 

{Yv,{ 
     Description -> "Neutrino-Yukawa-Coupling", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> False, 
     LesHouches -> Yv, 
     LaTeX -> "Y_\\nu", 
     OutputName -> Yv}}, 

{ZA,{ 
     Description -> "Pseudo-Scalar-Mixing-Matrix", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> True, 
     LesHouches -> PSEUDOSCALARMIX, 
     LaTeX -> "Z^A", 
     OutputName -> ZA}}, 

{ZD,{ 
     Description -> "Down-Squark-Mixing-Matrix", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> False, 
     LesHouches -> DSQMIX, 
     LaTeX -> "Z^D", 
     OutputName -> ZD}}, 

{ZDL,{ 
     Description -> "Left-Down-Mixing-Matrix", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> False, 
     LesHouches -> UDLMIX, 
     LaTeX -> "U^d_L", 
     OutputName -> ZDL}}, 

{ZDR,{ 
     Description -> "Right-Down-Mixing-Matrix", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> False, 
     LesHouches -> UDRMIX, 
     LaTeX -> "U^d_R", 
     OutputName -> ZDR}}, 

{ZE,{ 
     Description -> "Slepton-Mixing-Matrix", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> False, 
     LesHouches -> SELMIX, 
     LaTeX -> "Z^E", 
     OutputName -> ZE}}, 

{ZEL,{ 
     Description -> "Left-Lepton-Mixing-Matrix", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> False, 
     LesHouches -> UELMIX, 
     LaTeX -> "U^e_L", 
     OutputName -> ZEL}}, 

{ZER,{ 
     Description -> "Right-Lepton-Mixing-Matrix", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> False, 
     LesHouches -> UERMIX, 
     LaTeX -> "U^e_R", 
     OutputName -> ZER}}, 

{ZfW,{ 
     Description -> "Wino Mixing Matrix", 
     Dependence -> {{1/Sqrt[2], 1/Sqrt[2], 0}, {(-I)/Sqrt[2], I/Sqrt[2], 0}, {0, 0, 1}}, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> False, 
     LaTeX -> "Z^{\\tilde{W}}", 
     OutputName -> ZfW}}, 

{ZH,{ 
     Description -> "Scalar-Mixing-Matrix", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> True, 
     LesHouches -> SCALARMIX, 
     LaTeX -> "Z^H", 
     OutputName -> ZH}}, 

{ZN,{ 
     Description -> "Neutralino Mixing-Matrix", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> False, 
     LesHouches -> NMIX, 
     LaTeX -> "N", 
     OutputName -> ZN}}, 

{ZP,{ 
     Description -> "Charged-Mixing-Matrix", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> {{-Cos[\[Beta]], Sin[\[Beta]]}, {Sin[\[Beta]], Cos[\[Beta]]}}, 
     DependenceSPheno -> None, 
     Real -> True, 
     LesHouches -> CHARGEMIX, 
     LaTeX -> "Z^+", 
     OutputName -> ZP}}, 

{ZU,{ 
     Description -> "Up-Squark-Mixing-Matrix", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> False, 
     LesHouches -> USQMIX, 
     LaTeX -> "Z^U", 
     OutputName -> ZU}}, 

{ZUL,{ 
     Description -> "Left-Up-Mixing-Matrix", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> False, 
     LesHouches -> UULMIX, 
     LaTeX -> "U^u_L", 
     OutputName -> ZUL}}, 

{ZUR,{ 
     Description -> "Right-Up-Mixing-Matrix", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> False, 
     LesHouches -> UURMIX, 
     LaTeX -> "U^u_R", 
     OutputName -> ZUR}}, 

{ZVI,{ 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> False, 
     LesHouches -> SNUMIXI, 
     LaTeX -> "Z^I", 
     OutputName -> ZVI}}, 

{ZVR,{ 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> False, 
     LesHouches -> SNUMIXR, 
     LaTeX -> "Z^R", 
     OutputName -> ZVR}}, 

{ZW,{ 
     Description -> "W Mixing Matrix", 
     Dependence -> {{1/Sqrt[2], 1/Sqrt[2]}, {(-I)/Sqrt[2], I/Sqrt[2]}}, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> False, 
     LaTeX -> "Z^{W}", 
     OutputName -> ZW}}, 

{ZZ,{ 
     Description -> "Photon-Z-Z' Mixing Matrix", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> True, 
     LaTeX -> "Z^{\\gamma Z Z'}", 
     OutputName -> ZZ}}, 

{\[Beta],{ 
     Description -> "Pseudo Scalar mixing angle", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> ArcSin[Abs[ZP[1, 2]]], 
     Real -> True, 
     LesHouches -> {HMIX, 10}, 
     LaTeX -> "\\beta", 
     OutputName -> betaH}}, 

{\[Mu],{ 
     Description -> "Mu-parameter", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> False, 
     LesHouches -> {HMIX, 1}, 
     LaTeX -> "\\mu", 
     OutputName -> Mu}}, 

{B[\[Mu]],{ 
     Description -> "Bmu-parameter", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> False, 
     LesHouches -> {HMIX, 101}, 
     LaTeX -> "B_{\\mu}", 
     OutputName -> Bmu}}, 

{L[lw],{ 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> False, 
     LesHouches -> {SING, 22}, 
     LaTeX -> "L_{lw}"  (* auto generated LaTeX name *), 
     OutputName -> Llw}}, 

{T[Yd],{ 
     Description -> "Trilinear-Down-Coupling", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> False, 
     LesHouches -> Td, 
     LaTeX -> "T_d", 
     OutputName -> Td}}, 

{T[Ye],{ 
     Description -> "Trilinear-Lepton-Coupling", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> False, 
     LesHouches -> Te, 
     LaTeX -> "T_e", 
     OutputName -> Te}}, 

{T[Yn],{ 
     Description -> "Trilinear-Neutrino-X-Coupling", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> False, 
     LesHouches -> TX, 
     LaTeX -> "T_x", 
     OutputName -> Tx}}, 

{T[Yu],{ 
     Description -> "Trilinear-Up-Coupling", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> False, 
     LesHouches -> Tu, 
     LaTeX -> "T_u", 
     OutputName -> Tu}}, 

{T[Yv],{ 
     Description -> "Trilinear-Neutrino-Coupling", 
     Dependence -> None, 
     DependenceNum -> None, 
     DependenceOptional -> None, 
     DependenceSPheno -> None, 
     Real -> False, 
     LesHouches -> Tv, 
     LaTeX -> "T_\\nu", 
     OutputName -> Tv}}, 



(*  ----------- Additional parameters in the model ---------- *) 

{MVE,{ 
     Real -> False, 
     LesHouches -> MVE  (* auto generated Les Houches entry *), 
     LaTeX -> "M_E"  (* auto generated LaTeX name *), 
     OutputName -> mve  (* auto generated Output name *)}}, 

{B[MVE],{ 
     Real -> False, 
     LesHouches -> BVE  (* auto generated Les Houches entry *), 
     LaTeX -> "B_{E}"  (* auto generated LaTeX name *), 
     OutputName -> Bmve  (* auto generated Output name *)}}, 

{MS,{ 
     Real -> False, 
     LesHouches -> {U1XMSSM,6}  (* auto generated Les Houches entry *), 
     LaTeX -> "M_S"  (* auto generated LaTeX name *), 
     OutputName -> ms  (* auto generated Output name *)}}, 

{B[MS],{ 
     Real -> False, 
     LesHouches -> {U1XMSSM,8}  (* auto generated Les Houches entry *), 
     LaTeX -> "B_{S}"  (* auto generated LaTeX name *), 
     OutputName -> Bms  (* auto generated Output name *)}}, 

{MVU,{ 
     Real -> False, 
     LesHouches ->MVU  (* auto generated Les Houches entry *), 
     LaTeX -> "M_U"  (* auto generated LaTeX name *), 
     OutputName -> mvu  (* auto generated Output name *)}}, 

{B[MVU],{ 
     Real -> False, 
     LesHouches -> BVU  (* auto generated Les Houches entry *), 
     LaTeX -> "B_{U}"  (* auto generated LaTeX name *), 
     OutputName -> Bmvu  (* auto generated Output name *)}}, 

{Yep,{ 
     Real -> False, 
     LesHouches -> YEP  (* auto generated Les Houches entry *), 
     LaTeX -> "{Y'_e}"  (* auto generated LaTeX name *), 
     OutputName -> yep  (* auto generated Output name *)}}, 

{T[Yep],{ 
     Real -> False, 
     LesHouches -> TYEP  (* auto generated Les Houches entry *), 
     LaTeX -> "{T'_{e}}"  (* auto generated LaTeX name *), 
     OutputName -> Tyep  (* auto generated Output name *)}}, 

{lambdaC,{ 
     Real -> False, 
     LesHouches -> {U1XMSSM,14}  (* auto generated Les Houches entry *), 
     LaTeX -> "{\\lambda}_{C}"  (* auto generated LaTeX name *), 
     OutputName -> lamc  (* auto generated Output name *)}}, 

{T[lambdaC],{ 
     Real -> False, 
     LesHouches -> {U1XMSSM,16}  (* auto generated Les Houches entry *), 
     LaTeX -> "T_{{\\lambda}_{C}}"  (* auto generated LaTeX name *), 
     OutputName -> Tlc  (* auto generated Output name *)}}, 

{lambdaE,{ 
     Real -> False, 
     LesHouches -> LAME  (* auto generated Les Houches entry *), 
     LaTeX -> "{\\lambda}_{E}"  (* auto generated LaTeX name *), 
     OutputName -> lame  (* auto generated Output name *)}}, 

{T[lambdaE],{ 
     Real -> False, 
     LesHouches -> TLAME  (* auto generated Les Houches entry *), 
     LaTeX -> "T_{{\\lambda}_{E}}"  (* auto generated LaTeX name *), 
     OutputName -> Tle  (* auto generated Output name *)}}, 

{lambdaH,{ 
     Real -> False, 
     LesHouches -> {U1XMSSM,22}  (* auto generated Les Houches entry *), 
     LaTeX -> "{\\lambda}_{H}"  (* auto generated LaTeX name *), 
     OutputName -> lamh  (* auto generated Output name *)}}, 

{T[lambdaH],{ 
     Real -> False, 
     LesHouches -> {U1XMSSM,24}  (* auto generated Les Houches entry *), 
     LaTeX -> "T_{{\\lambda}_{H}}"  (* auto generated LaTeX name *), 
     OutputName -> Tlh  (* auto generated Output name *)}}, 

{kappa,{ 
     Real -> False, 
     LesHouches -> {U1XMSSM,26}  (* auto generated Les Houches entry *), 
     LaTeX -> "\\kappa"  (* auto generated LaTeX name *), 
     OutputName -> kap  (* auto generated Output name *)}}, 

{T[kappa],{ 
     Real -> False, 
     LesHouches -> {U1XMSSM,28}  (* auto generated Les Houches entry *), 
     LaTeX -> "T_{\\kappa}"  (* auto generated LaTeX name *), 
     OutputName -> Tkap  (* auto generated Output name *)}}, 

{lambdaU,{ 
     Real -> False, 
     LesHouches -> LAMU  (* auto generated Les Houches entry *), 
     LaTeX -> "{\\lambda}_{U}"  (* auto generated LaTeX name *), 
     OutputName -> lamu  (* auto generated Output name *)}}, 

{T[lambdaU],{ 
     Real -> False, 
     LesHouches -> TLAMU (* auto generated Les Houches entry *), 
     LaTeX -> "T_{{\\lambda}_{U}}"  (* auto generated LaTeX name *), 
     OutputName -> Tlu  (* auto generated Output name *)}}, 

{Yup,{ 
     Real -> False, 
     LesHouches -> YUP  (* auto generated Les Houches entry *), 
     LaTeX -> "{Y'_u}"  (* auto generated LaTeX name *), 
     OutputName -> yup  (* auto generated Output name *)}}, 

{T[Yup],{ 
     Real -> False, 
     LesHouches -> TYUP  (* auto generated Les Houches entry *), 
     LaTeX -> "{T'_{u}}"  (* auto generated LaTeX name *), 
     OutputName -> Tyup  (* auto generated Output name *)}}, 

{muUX2,{ 
     Real -> False, 
     LesHouches -> MUUX2  (* auto generated Les Houches entry *), 
     LaTeX -> "m^2_{uUX}"  (* auto generated LaTeX name *), 
     OutputName -> muux2  (* auto generated Output name *)}}, 
     

{meEX2,{ 
     Real -> False, 
     LesHouches -> MEEX2  (* auto generated Les Houches entry *), 
     LaTeX -> "m^2_{eEX}"  (* auto generated LaTeX name *), 
     OutputName -> meex2  (* auto generated Output name *)}},      

{mUX2,{ 
     Real -> False, 
     LesHouches -> MUX2  (* auto generated Les Houches entry *), 
     LaTeX -> "m^2_{UX}"  (* auto generated LaTeX name *), 
     OutputName -> mux2  (* auto generated Output name *)}}, 

{mUXp2,{ 
     Real -> False, 
     LesHouches -> MUXP2  (* auto generated Les Houches entry *), 
     LaTeX -> "m^2_{UXp}"  (* auto generated LaTeX name *), 
     OutputName -> muxp2  (* auto generated Output name *)}}, 

{mEX2,{ 
     Real -> False, 
     LesHouches -> MEX2  (* auto generated Les Houches entry *), 
     LaTeX -> "m^2_{EX}"  (* auto generated LaTeX name *), 
     OutputName -> mex2  (* auto generated Output name *)}}, 

{mEXp2,{ 
     Real -> False, 
     LesHouches -> MEXP2  (* auto generated Les Houches entry *), 
     LaTeX -> "m^2_{EXp}"  (* auto generated LaTeX name *), 
     OutputName -> mexp2  (* auto generated Output name *)}}, 
     
{MtE,{ 
     Real -> False, 
     LesHouches -> MtE  (* auto generated Les Houches entry *), 
     LaTeX -> "\\tilde{M}_E"  (* auto generated LaTeX name *), 
     OutputName -> MtE  (* auto generated Output name *)}},      
     
{MtU,{ 
     Real -> False, 
     LesHouches -> MtU  (* auto generated Les Houches entry *), 
     LaTeX -> "\\tilde{M}_U"  (* auto generated LaTeX name *), 
     OutputName -> MtU  (* auto generated Output name *)}},  
     
{B[MtE],{ 
     Real -> False, 
     LesHouches -> BtE  (* auto generated Les Houches entry *), 
     LaTeX -> "\\tilde{B}_E"  (* auto generated LaTeX name *), 
     OutputName -> BtE  (* auto generated Output name *)}},      
     
{B[MtU],{ 
     Real -> False, 
     LesHouches -> BtU  (* auto generated Les Houches entry *), 
     LaTeX -> "\\tilde{B}_U"  (* auto generated LaTeX name *), 
     OutputName -> BtU  (* auto generated Output name *)}},        

{Mass[VWm],{ 
     Real -> False, 
     LesHouches -> {U1XMSSM,42}  (* auto generated Les Houches entry *), 
     LaTeX -> "M_{VWm}"  (* auto generated LaTeX name *), 
     OutputName -> Mvwm  (* auto generated Output name *)}} 
} 
