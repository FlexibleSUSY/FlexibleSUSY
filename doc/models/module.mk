DIR          := doc/models
MODNAME      := doc-models

DOC_MODELS_DIR := $(DIR)

DOC_MODELS_MK := $(DOC_MODELS_DIR)/module.mk

DOC_MODELS_INSTALL_DIR := $(INSTALL_DIR)/$(DOC_MODELS_DIR)

DOC_MODELS_RST := \
		$(DOC_MODELS_DIR)/HSSUSY.rst \
		$(DOC_MODELS_DIR)/MSSMEFTHiggs.rst \
		$(DOC_MODELS_DIR)/NUHMSSMNoFVHimalaya.rst

DOC_MODELS_IMAGE_DIR := $(DOC_MODELS_DIR)/images

DOC_MODELS_IMAGES := \
		$(DOC_MODELS_IMAGE_DIR)/HSSUSY_Mh_MS.png \
		$(DOC_MODELS_IMAGE_DIR)/HSSUSY_Mh_Xt.png \
		$(DOC_MODELS_IMAGE_DIR)/HSSUSY_tower.svg \
		$(DOC_MODELS_IMAGE_DIR)/MSSMEFTHiggs_Mh_MS.png \
		$(DOC_MODELS_IMAGE_DIR)/MSSMEFTHiggs_tower.svg \
		$(DOC_MODELS_IMAGE_DIR)/NUHMSSMNoFVHimalaya_tower.svg

DOC_MODELS_EXAMPLE_DIR := $(DOC_MODELS_DIR)/examples

DOC_MODELS_EXAMPLES := \
		$(DOC_MODELS_EXAMPLE_DIR)/HSSUSY_uncertainty_estimate.m \
		$(DOC_MODELS_EXAMPLE_DIR)/MSSMEFTHiggs_uncertainty_estimate.m \
		$(DOC_MODELS_EXAMPLE_DIR)/NUHMSSMNoFVHimalaya_uncertainty_estimate.m

ifneq ($(INSTALL_DIR),)
install-src::
		$(Q)install -d $(DOC_MODELS_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(DOC_MODELS_MK) $(DOC_MODELS_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(DOC_MODELS_RST) $(DOC_MODELS_INSTALL_DIR)
		$(Q)install -d $(INSTALL_DIR)/$(DOC_MODELS_IMAGE_DIR)
		$(Q)install -m u=rw,g=r,o=r $(DOC_MODELS_IMAGES) $(INSTALL_DIR)/$(DOC_MODELS_IMAGE_DIR)
		$(Q)install -d $(INSTALL_DIR)/$(DOC_MODELS_EXAMPLE_DIR)
		$(Q)install -m u=rw,g=r,o=r $(DOC_MODELS_EXAMPLES) $(INSTALL_DIR)/$(DOC_MODELS_EXAMPLE_DIR)
endif
