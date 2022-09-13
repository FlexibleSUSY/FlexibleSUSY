DIR          := doc/manuals
MODNAME      := doc-manuals
MANUALS_MK   := $(DIR)/module.mk

MANUALS_DIR  := $(DIR)

PDF_OUTPUT_DIR  := $(MANUALS_DIR)

MANUALS_INSTALL_DIR := $(INSTALL_DIR)/$(MANUALS_DIR)

MANUALS_EXAMPLES_DIR := $(MANUALS_DIR)/examples-2.0

MANUALS_PLOTS_DIR := $(MANUALS_DIR)/plots

MANUALS_PLOTS_THDMIIMSSMBC := \
		$(MANUALS_PLOTS_DIR)/THDMIIMSSMBC/THDMIIMSSMBC_Mh_MS_TB_MA-200_Xt-sqrt6.pdf \
		$(MANUALS_PLOTS_DIR)/THDMIIMSSMBC/THDMIIMSSMBC_Mh_MS_TB_MA-200_Xt-0.pdf

MANUALS_PLOTS_HSSUSY := \
		$(MANUALS_PLOTS_DIR)/HSSUSY/Mh_MS_EFT_uncertainty.pdf \
		$(MANUALS_PLOTS_DIR)/HSSUSY/HSSUSY_thresholds.pdf \
		$(MANUALS_PLOTS_DIR)/HSSUSY/Mh_Xt.pdf \
		$(MANUALS_PLOTS_DIR)/HSSUSY/Mh_Xt_uncertainty.pdf

MANUALS_PLOTS_FlexibleMW := \
		$(MANUALS_PLOTS_DIR)/FlexibleMW/scan_MW_M12.pdf \
		$(MANUALS_PLOTS_DIR)/FlexibleMW/scan_MW_LamTU.pdf

MANUALS_PLOTS_CMSSM_multiple_solutions := \
		$(MANUALS_PLOTS_DIR)/CMSSM_multiple_solutions/CMSSM_multiple_solutions.pdf

MANUALS_PLOTS_CE6SSM := \
		$(MANUALS_PLOTS_DIR)/CE6SSM/ce6ssm_tb10-scan-l-k-vs_combined_m0_m12_ccA0.pdf \
		$(MANUALS_PLOTS_DIR)/CE6SSM/ce6ssm_tb10-scan-l-k-vs_combined_m0_m12_ccmsu1.pdf \
		$(MANUALS_PLOTS_DIR)/CE6SSM/ce6ssm_tb10-scan-l-k-vs_combined_m0_m12_cclam.pdf \
		$(MANUALS_PLOTS_DIR)/CE6SSM/ce6ssm_tb10-scan-l-k-vs_combined_m0_m12_cckap.pdf

MANUALS_PLOTS_CSE6SSM := \
		$(MANUALS_PLOTS_DIR)/CSE6SSM/CSE6SSM_m12_Azero_m0_handwritten.pdf \
		$(MANUALS_PLOTS_DIR)/CSE6SSM/CSE6SSM_m12_Azero_m0.pdf

MANUALS_PLOTS_CNMSSM := \
		$(MANUALS_PLOTS_DIR)/CNMSSM/cnmssm_4D-m0sqpos-full-m0_m12_sort-A0.pdf \
		$(MANUALS_PLOTS_DIR)/CNMSSM/cnmssm_4D-full-squark_mass_m0sqgtr0_m0_m12.pdf \
		$(MANUALS_PLOTS_DIR)/CNMSSM/cnmssm_4D-m0sqpos-full_mhZoom_lam_m12_sort-mh.pdf \
		$(MANUALS_PLOTS_DIR)/CNMSSM/cnmssm_4D-m0sqgtr0-full-msu6.pdf

MANUALS_PLOTS_semi_analytic_benchmark := \
		$(MANUALS_PLOTS_DIR)/semi_analytic_benchmark/i7-4702MQ.pdf

MANUALS_PLOTS_Himalaya := \
		$(MANUALS_PLOTS_DIR)/Himalaya/scan_Mh_MS_TB-5_Xt-0_delta_alpha.pdf \
		$(MANUALS_PLOTS_DIR)/Himalaya/scan_Mh_Xt_TB-5_MS-2000.pdf \
		$(MANUALS_PLOTS_DIR)/Himalaya/scan_Mh_MS_TB-5_Xt-0.pdf \
		$(MANUALS_PLOTS_DIR)/Himalaya/scan_Mh_Xt_TB-5_MS-2000_delta_alpha.pdf

MANUALS_PLOTS_CMSSMCPV := \
		$(MANUALS_PLOTS_DIR)/CMSSMCPV/scan_TYu33_Q_m0-500_M12-500_TB-10_A0-500_PhiMu-0_PhiA0-Pi4.pdf \
		$(MANUALS_PLOTS_DIR)/CMSSMCPV/scan_Mh_PhiMu_m0-500_M12-500_TB-10_A0-0.pdf

MANUALS_PLOTS_FlexibleEFTHiggs := \
		$(MANUALS_PLOTS_DIR)/FlexibleEFTHiggs/scan_Mh_Xt_TB-5_MS-2000.pdf \
		$(MANUALS_PLOTS_DIR)/FlexibleEFTHiggs/scan_Mh_MS_TB-5_Xt-0.pdf

MANUALS_PLOTS_CMSSM_solution_regions := \
		$(MANUALS_PLOTS_DIR)/CMSSM_solution_regions/CMSSM_higgs_mass.pdf \
		$(MANUALS_PLOTS_DIR)/CMSSM_solution_regions/CMSSM_solution_regions.pdf

MANUALS_PLOTS_semi_analytic_rel_errs := \
		$(MANUALS_PLOTS_DIR)/semi_analytic_rel_errs/single_step_errors.pdf

MANUALS_PLOTS_MSSMMuBMu := \
		$(MANUALS_PLOTS_DIR)/MSSMMuBMu/scan_MSSMMuBMu_Mh_Xt_TB-5_MS-2000.pdf \
		$(MANUALS_PLOTS_DIR)/MSSMMuBMu/scan_MSSMMuBMu_Mh_MS_TB-5_Xt-0.pdf

MANUALS_PLOTS_MSSMCPV := \
		$(MANUALS_PLOTS_DIR)/MSSMCPV/contour_de_PhiMuTe.pdf \
		$(MANUALS_PLOTS_DIR)/MSSMCPV/contour_de_PhiMuM12.pdf

MANUALS_FEYNMAN_DIR := $(DIR)/feynman

MANUALS_FEYNMAN := \
		$(MANUALS_FEYNMAN_DIR)/deltaVBVertex1.pdf \
		$(MANUALS_FEYNMAN_DIR)/amuFFS.pdf \
		$(MANUALS_FEYNMAN_DIR)/deltaVBWave.pdf \
		$(MANUALS_FEYNMAN_DIR)/deltaVBVertex2.pdf \
		$(MANUALS_FEYNMAN_DIR)/amuFSS.pdf \
		$(MANUALS_FEYNMAN_DIR)/deltaVBBox1.pdf

MANUALS_EXAMPLES := \
		$(MANUALS_EXAMPLES_DIR)/scan_FlexibleEFTHiggs_uncertainty.m \
		$(MANUALS_EXAMPLES_DIR)/scan_HSSUSY.m \
		$(MANUALS_EXAMPLES_DIR)/scan_HSSUSY_EFT_uncertainty.m \
		$(MANUALS_EXAMPLES_DIR)/scan_HSSUSY_SM_uncertainty.m

PAPER_PDF_1     := $(PDF_OUTPUT_DIR)/flexiblesusy-1.0.pdf
PAPER_PDF_2     := $(PDF_OUTPUT_DIR)/flexiblesusy-2.0.pdf
PAPER_PDF_3     := $(PDF_OUTPUT_DIR)/flexiblesusy-new_features.pdf
PAPER_PDF       := $(PAPER_PDF_1) $(PAPER_PDF_2) $(PAPER_PDF_3)
PAPER_SRC_1     := $(DIR)/flexiblesusy-1.0.tex
PAPER_SRC_2     := \
		$(DIR)/flexiblesusy-2.0.tex \
		$(DIR)/flexiblesusy-2.0.bib
PAPER_SRC_3     := \
		$(DIR)/flexiblesusy-new_features.tex \
		$(DIR)/flexiblesusy-new_features.bib
PAPER_SRC       := $(PAPER_SRC_1) $(PAPER_SRC_2) $(PAPER_SRC_3)
PAPER_STY       := $(DIR)/JHEP.bst $(DIR)/tikz-uml.sty

LATEX_TMP       := \
		$(patsubst %.pdf, %.aux, $(PAPER_PDF)) \
		$(patsubst %.pdf, %.bbl, $(PAPER_PDF)) \
		$(patsubst %.pdf, %.blg, $(PAPER_PDF)) \
		$(patsubst %.pdf, %.log, $(PAPER_PDF)) \
		$(patsubst %.pdf, %.toc, $(PAPER_PDF)) \
		$(patsubst %.pdf, %.out, $(PAPER_PDF)) \
		$(patsubst %.pdf, %.spl, $(PAPER_PDF))

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME) \
		$(MODNAME)-pdf

clean-$(MODNAME):
		$(Q)-rm -f $(LATEX_TMP)

distclean-$(MODNAME): clean-$(MODNAME)
		$(Q)-rm -f $(PAPER_PDF)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

ifneq ($(INSTALL_DIR),)
install-src::
		$(Q)install -d $(MANUALS_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(MANUALS_MK) $(MANUALS_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(PAPER_SRC) $(MANUALS_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(PAPER_STY) $(MANUALS_INSTALL_DIR)
		$(Q)install -d $(INSTALL_DIR)/$(MANUALS_EXAMPLES_DIR)
		$(Q)install -m u=rw,g=r,o=r $(MANUALS_EXAMPLES) $(INSTALL_DIR)/$(MANUALS_EXAMPLES_DIR)

		$(Q)install -d $(INSTALL_DIR)/$(MANUALS_FEYNMAN_DIR)
		$(Q)install -m u=rw,g=r,o=r $(MANUALS_FEYNMAN) $(INSTALL_DIR)/$(MANUALS_FEYNMAN_DIR)

		$(Q)install -d $(INSTALL_DIR)/$(MANUALS_PLOTS_DIR)/THDMIIMSSMBC
		$(Q)install -m u=rw,g=r,o=r $(MANUALS_PLOTS_THDMIIMSSMBC) $(INSTALL_DIR)/$(MANUALS_PLOTS_DIR)/THDMIIMSSMBC

		$(Q)install -d $(INSTALL_DIR)/$(MANUALS_PLOTS_DIR)/HSSUSY
		$(Q)install -m u=rw,g=r,o=r $(MANUALS_PLOTS_HSSUSY) $(INSTALL_DIR)/$(MANUALS_PLOTS_DIR)/HSSUSY

		$(Q)install -d $(INSTALL_DIR)/$(MANUALS_PLOTS_DIR)/FlexibleMW
		$(Q)install -m u=rw,g=r,o=r $(MANUALS_PLOTS_FlexibleMW) $(INSTALL_DIR)/$(MANUALS_PLOTS_DIR)/FlexibleMW

		$(Q)install -d $(INSTALL_DIR)/$(MANUALS_PLOTS_DIR)/CMSSM_multiple_solutions
		$(Q)install -m u=rw,g=r,o=r $(MANUALS_PLOTS_CMSSM_multiple_solutions) $(INSTALL_DIR)/$(MANUALS_PLOTS_DIR)/CMSSM_multiple_solutions

		$(Q)install -d $(INSTALL_DIR)/$(MANUALS_PLOTS_DIR)/CE6SSM
		$(Q)install -m u=rw,g=r,o=r $(MANUALS_PLOTS_CE6SSM) $(INSTALL_DIR)/$(MANUALS_PLOTS_DIR)/CE6SSM

		$(Q)install -d $(INSTALL_DIR)/$(MANUALS_PLOTS_DIR)/CSE6SSM
		$(Q)install -m u=rw,g=r,o=r $(MANUALS_PLOTS_CSE6SSM) $(INSTALL_DIR)/$(MANUALS_PLOTS_DIR)/CSE6SSM

		$(Q)install -d $(INSTALL_DIR)/$(MANUALS_PLOTS_DIR)/CNMSSM
		$(Q)install -m u=rw,g=r,o=r $(MANUALS_PLOTS_CNMSSM) $(INSTALL_DIR)/$(MANUALS_PLOTS_DIR)/CNMSSM

		$(Q)install -d $(INSTALL_DIR)/$(MANUALS_PLOTS_DIR)/semi_analytic_benchmark
		$(Q)install -m u=rw,g=r,o=r $(MANUALS_PLOTS_semi_analytic_benchmark) $(INSTALL_DIR)/$(MANUALS_PLOTS_DIR)/semi_analytic_benchmark

		$(Q)install -d $(INSTALL_DIR)/$(MANUALS_PLOTS_DIR)/Himalaya
		$(Q)install -m u=rw,g=r,o=r $(MANUALS_PLOTS_Himalaya) $(INSTALL_DIR)/$(MANUALS_PLOTS_DIR)/Himalaya

		$(Q)install -d $(INSTALL_DIR)/$(MANUALS_PLOTS_DIR)/CMSSMCPV
		$(Q)install -m u=rw,g=r,o=r $(MANUALS_PLOTS_CMSSMCPV) $(INSTALL_DIR)/$(MANUALS_PLOTS_DIR)/CMSSMCPV

		$(Q)install -d $(INSTALL_DIR)/$(MANUALS_PLOTS_DIR)/FlexibleEFTHiggs
		$(Q)install -m u=rw,g=r,o=r $(MANUALS_PLOTS_FlexibleEFTHiggs) $(INSTALL_DIR)/$(MANUALS_PLOTS_DIR)/FlexibleEFTHiggs

		$(Q)install -d $(INSTALL_DIR)/$(MANUALS_PLOTS_DIR)/CMSSM_solution_regions
		$(Q)install -m u=rw,g=r,o=r $(MANUALS_PLOTS_CMSSM_solution_regions) $(INSTALL_DIR)/$(MANUALS_PLOTS_DIR)/CMSSM_solution_regions

		$(Q)install -d $(INSTALL_DIR)/$(MANUALS_PLOTS_DIR)/semi_analytic_rel_errs
		$(Q)install -m u=rw,g=r,o=r $(MANUALS_PLOTS_semi_analytic_rel_errs) $(INSTALL_DIR)/$(MANUALS_PLOTS_DIR)/semi_analytic_rel_errs

		$(Q)install -d $(INSTALL_DIR)/$(MANUALS_PLOTS_DIR)/MSSMMuBMu
		$(Q)install -m u=rw,g=r,o=r $(MANUALS_PLOTS_MSSMMuBMu) $(INSTALL_DIR)/$(MANUALS_PLOTS_DIR)/MSSMMuBMu

		$(Q)install -d $(INSTALL_DIR)/$(MANUALS_PLOTS_DIR)/MSSMCPV
		$(Q)install -m u=rw,g=r,o=r $(MANUALS_PLOTS_MSSMCPV) $(INSTALL_DIR)/$(MANUALS_PLOTS_DIR)/MSSMCPV

endif

$(MODNAME)-pdf: $(PAPER_PDF)

$(PAPER_PDF_1): $(PAPER_SRC_1) $(PAPER_STY)
		$(Q)cd $(MANUALS_DIR) && pdflatex flexiblesusy-1.0.tex
		$(Q)cd $(MANUALS_DIR) && pdflatex flexiblesusy-1.0.tex
		$(Q)cd $(MANUALS_DIR) && pdflatex flexiblesusy-1.0.tex

$(PAPER_PDF_2): $(PAPER_SRC_2) $(PAPER_STY)
		$(Q)cd $(MANUALS_DIR) && pdflatex flexiblesusy-2.0.tex
		$(Q)cd $(MANUALS_DIR) && bibtex flexiblesusy-2.0
		$(Q)cd $(MANUALS_DIR) && pdflatex flexiblesusy-2.0.tex
		$(Q)cd $(MANUALS_DIR) && pdflatex flexiblesusy-2.0.tex

$(PAPER_PDF_3): $(PAPER_SRC_3) $(PAPER_STY)
		$(Q)cd $(MANUALS_DIR) && pdflatex flexiblesusy-new_features.tex
		$(Q)cd $(MANUALS_DIR) && bibtex flexiblesusy-new_features
		$(Q)cd $(MANUALS_DIR) && pdflatex flexiblesusy-new_features.tex
		$(Q)cd $(MANUALS_DIR) && pdflatex flexiblesusy-new_features.tex
