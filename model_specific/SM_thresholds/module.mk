DIR          := model_specific/SM_thresholds
MODNAME      := model_specific_SM_thresholds
WITH_$(MODNAME) := yes

LIB_model_specific_SM_thresholds_MK  := \
		$(DIR)/module.mk

LIB_model_specific_SM_thresholds_SRC := \
		$(DIR)/sm_twoloop_mt.cpp

LIB_model_specific_SM_thresholds_HDR := \
		$(DIR)/sm_twoloop_mt.hpp

LIB_model_specific_SM_thresholds_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIB_model_specific_SM_thresholds_SRC)))

LIB_model_specific_SM_thresholds_DEP := $(LIB_model_specific_SM_thresholds_OBJ:.o=.d)

LIB_model_specific_SM_thresholds     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

LIB_model_specific_SM_thresholds_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-dep \
		clean-$(MODNAME)-lib clean-$(MODNAME)-obj distclean-$(MODNAME)

all-$(MODNAME): $(LIB_model_specific_SM_thresholds)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		$(Q)install -d $(LIB_model_specific_SM_thresholds_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIB_model_specific_SM_thresholds_SRC) $(LIB_model_specific_SM_thresholds_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIB_model_specific_SM_thresholds_HDR) $(LIB_model_specific_SM_thresholds_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIB_model_specific_SM_thresholds_MK) $(LIB_model_specific_SM_thresholds_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		$(Q)-rm -f $(LIB_model_specific_SM_thresholds_DEP)

clean-$(MODNAME)-lib:
		$(Q)-rm -f $(LIB_model_specific_SM_thresholds)

clean-$(MODNAME)-obj:
		$(Q)-rm -f $(LIB_model_specific_SM_thresholds_OBJ)

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		@true

distclean-$(MODNAME): clean-$(MODNAME)

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(LIB_model_specific_SM_thresholds_DEP) $(LIB_model_specific_SM_thresholds_OBJ): CPPFLAGS += $(EIGENFLAGS) $(BOOSTFLAGS) $(TSILFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIB_model_specific_SM_thresholds_DEP) $(LIB_model_specific_SM_thresholds_OBJ): CPPFLAGS += $(LOOPFUNCFLAGS)
endif

ifeq ($(ENABLE_SHARED_LIBS),yes)
$(LIB_model_specific_SM_thresholds): $(LIB_model_specific_SM_thresholds_OBJ)
		@$(MSG)
		$(Q)$(MODULE_MAKE_LIB_CMD) $@ $^ $(TSILLIBS)
else
$(LIB_model_specific_SM_thresholds): $(LIB_model_specific_SM_thresholds_OBJ)
		@$(MSG)
		$(Q)$(MODULE_MAKE_LIB_CMD) $@ $^
endif

ALLDEP += $(LIB_model_specific_SM_thresholds_DEP)
ALLLIB += $(LIB_model_specific_SM_thresholds)
