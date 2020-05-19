DIR          := meta
MODNAME      := meta

META_MSSM_SRC:= \
		$(DIR)/MSSM/BDSZHiggs.m \
		$(DIR)/MSSM/beta_BMu.m \
		$(DIR)/MSSM/beta_g1.m \
		$(DIR)/MSSM/beta_g2.m \
		$(DIR)/MSSM/beta_g3.m \
		$(DIR)/MSSM/beta_M1.m \
		$(DIR)/MSSM/beta_M2.m \
		$(DIR)/MSSM/beta_M3.m \
		$(DIR)/MSSM/beta_md2.m \
		$(DIR)/MSSM/beta_me2.m \
		$(DIR)/MSSM/beta_mHd2.m \
		$(DIR)/MSSM/beta_mHu2.m \
		$(DIR)/MSSM/beta_ml2.m \
		$(DIR)/MSSM/beta_mq2.m \
		$(DIR)/MSSM/beta_mu2.m \
		$(DIR)/MSSM/beta_Mu.m \
		$(DIR)/MSSM/beta_TYd.m \
		$(DIR)/MSSM/beta_TYe.m \
		$(DIR)/MSSM/beta_TYu.m \
		$(DIR)/MSSM/beta_Yd.m \
		$(DIR)/MSSM/beta_Ye.m \
		$(DIR)/MSSM/beta_Yu.m \
		$(DIR)/MSSM/gamma_SdR.m \
		$(DIR)/MSSM/gamma_SeR.m \
		$(DIR)/MSSM/gamma_SHd.m \
		$(DIR)/MSSM/gamma_SHu.m \
		$(DIR)/MSSM/gamma_SlL.m \
		$(DIR)/MSSM/gamma_SqL.m \
		$(DIR)/MSSM/gamma_SuR.m \
		$(DIR)/MSSM/Mh2_3loop_DR_SQCD.m \
		$(DIR)/MSSM/das2.m \
		$(DIR)/MSSM/dmtauas2.m \
		$(DIR)/MSSM/bquark_2loop_sqcd_decoupling.m \
		$(DIR)/MSSM/bquark_2loop_sqcd_pole.m \
		$(DIR)/MSSM/tquark_1loop_qcd.m \
		$(DIR)/MSSM/tquark_1loop_strong.m \
		$(DIR)/MSSM/tquark_2loop_qcd.m \
		$(DIR)/MSSM/tquark_2loop_strong.m \
		$(DIR)/MSSM/twoloopbubble.m

META_SM_SRC  := \
		$(DIR)/SM/beta_g1.m \
		$(DIR)/SM/beta_g2.m \
		$(DIR)/SM/beta_g3.m \
		$(DIR)/SM/beta_gb.m \
		$(DIR)/SM/beta_gtau.m \
		$(DIR)/SM/beta_gt.m \
		$(DIR)/SM/beta_lambda.m \
		$(DIR)/SM/beta_m2.m \
		$(DIR)/SM/HSSUSY_corrections.m \
		$(DIR)/SM/HSSUSY_scale_variation.m \
		$(DIR)/SM/HSSUSY_scale_variation_1L.m \
		$(DIR)/SM/HSSUSY_scale_variation_2L.m \
		$(DIR)/SM/HSSUSY_scale_variation_3L.m

META_THDM_SRC:= \
		$(DIR)/THDM/Thresholds_1L_full.m

META_SRC     := \
		$(DIR)/AMuon.m \
		$(DIR)/AnomalousDimension.m \
		$(DIR)/BetaFunction.m \
		$(DIR)/CConversion.m \
		$(DIR)/Constraint.m \
		$(DIR)/ConvergenceTester.m \
		$(DIR)/CXXDiagrams.m \
		$(DIR)/NPointFunctions.m \
		$(DIR)/NPointFunctions/internal.m \
		$(DIR)/NPointFunctions/createFAModelFile.m \
		$(DIR)/WilsonCoeffs.m \
		$(DIR)/EDM.m \
		$(DIR)/FFVFormFactors.m \
		$(DIR)/FToFConversionInNucleus.m \
		$(DIR)/BrLToLGamma.m \
		$(DIR)/BtoSGamma.m \
		$(DIR)/EffectiveCouplings.m \
		$(DIR)/EWSB.m \
		$(DIR)/FlexibleEFTHiggsMatching.m \
		$(DIR)/FlexibleSUSY.m \
		$(DIR)/FlexibleTower.m \
		$(DIR)/Format.m \
		$(DIR)/FSMathLink.m \
		$(DIR)/FunctionModifiers.m \
		$(DIR)/Himalaya.m \
		$(DIR)/LatticeUtils.m \
		$(DIR)/LoopFunctions.m \
		$(DIR)/LoopFunctionsZeroMomentum.m \
		$(DIR)/LoopMasses.m \
		$(DIR)/Observables.m \
		$(DIR)/Parameters.m \
		$(DIR)/Phases.m \
		$(DIR)/ReadSLHA.m \
		$(DIR)/References.m \
		$(DIR)/RGIntegrator.m \
		$(DIR)/SelfEnergies.m \
		$(DIR)/SemiAnalytic.m \
		$(DIR)/TestSuite.m \
		$(DIR)/TextFormatting.m \
		$(DIR)/TerminalFormatting.m \
		$(DIR)/ThreeLoopMSSM.m \
		$(DIR)/ThreeLoopQCD.m \
		$(DIR)/ThreeLoopSM.m \
		$(DIR)/ThresholdCorrections.m \
		$(DIR)/ThresholdCorrectionLoopFunctions.m \
		$(DIR)/Traces.m \
		$(DIR)/TreeMasses.m \
		$(DIR)/TwoLoopMSSM.m \
		$(DIR)/TwoLoopQCD.m \
		$(DIR)/TwoLoopSM.m \
		$(DIR)/Utils.m \
		$(DIR)/Vertices.m \
		$(DIR)/WeinbergAngle.m \
		$(DIR)/WriteOut.m \
		$(DIR)/writeRGE.m \
		$(DIR)/writeNRGE.m \
		$(META_MSSM_SRC) \
		$(META_SM_SRC) \
		$(META_THDM_SRC)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

all-$(MODNAME):
		@true

clean-$(MODNAME):
		@true

distclean-$(MODNAME): clean-$(MODNAME)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)
