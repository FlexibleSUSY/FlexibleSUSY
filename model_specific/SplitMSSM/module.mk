DIR          := model_specific/SplitMSSM
MODNAME      := model_specific_SplitMSSM
WITH_$(MODNAME) := yes

LIB_model_specific_SplitMSSM_MK  := \
		$(DIR)/module.mk

LIB_model_specific_SplitMSSM_SRC := \
		$(DIR)/splitmssm_threeloophiggs.cpp \
		$(DIR)/splitmssm_thresholds.cpp

LIB_model_specific_SplitMSSM_HDR := \
		$(DIR)/splitmssm_threeloophiggs.hpp \
		$(DIR)/splitmssm_thresholds.hpp

LIB_model_specific_SplitMSSM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIB_model_specific_SplitMSSM_SRC)))

LIB_model_specific_SplitMSSM_DEP := $(LIB_model_specific_SplitMSSM_OBJ:.o=.d)

LIB_model_specific_SplitMSSM     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

LIB_model_specific_SplitMSSM_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-dep \
		clean-$(MODNAME)-lib clean-$(MODNAME)-obj distclean-$(MODNAME)

all-$(MODNAME): $(LIB_model_specific_SplitMSSM)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		$(Q)install -d $(LIB_model_specific_SplitMSSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIB_model_specific_SplitMSSM_SRC) $(LIB_model_specific_SplitMSSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIB_model_specific_SplitMSSM_HDR) $(LIB_model_specific_SplitMSSM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIB_model_specific_SplitMSSM_MK) $(LIB_model_specific_SplitMSSM_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		$(Q)-rm -f $(LIB_model_specific_SplitMSSM_DEP)

clean-$(MODNAME)-lib:
		$(Q)-rm -f $(LIB_model_specific_SplitMSSM)

clean-$(MODNAME)-obj:
		$(Q)-rm -f $(LIB_model_specific_SplitMSSM_OBJ)

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		@true

distclean-$(MODNAME): clean-$(MODNAME)

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(LIB_model_specific_SplitMSSM_DEP) $(LIB_model_specific_SplitMSSM_OBJ): CPPFLAGS += $(EIGENFLAGS) $(BOOSTFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIB_model_specific_SplitMSSM_DEP) $(LIB_model_specific_SplitMSSM_OBJ): CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LIB_model_specific_SplitMSSM): $(LIB_model_specific_SplitMSSM_OBJ)
		@$(MSG)
		$(Q)$(MODULE_MAKE_LIB_CMD) $@ $^

ALLDEP += $(LIB_model_specific_SplitMSSM_DEP)
ALLLIB += $(LIB_model_specific_SplitMSSM)
