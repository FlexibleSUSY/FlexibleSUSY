DIR          := model_specific/MSSM_thresholds
MODNAME      := model_specific_MSSM_thresholds
WITH_$(MODNAME) := yes

LIB_model_specific_MSSM_thresholds_MK  := \
		$(DIR)/module.mk

LIB_model_specific_MSSM_thresholds_SRC := \
		$(DIR)/mssm_twoloop_as.cpp \
		$(DIR)/mssm_twoloop_mb.cpp \
		$(DIR)/mssm_twoloop_mt.cpp \
		$(DIR)/mssm_twoloop_mtau.cpp

LIB_model_specific_MSSM_thresholds_HDR := \
		$(DIR)/mssm_twoloop_as.hpp \
		$(DIR)/mssm_twoloop_mb.hpp \
		$(DIR)/mssm_twoloop_mt.hpp \
		$(DIR)/mssm_twoloop_mtau.hpp

LIB_model_specific_MSSM_thresholds_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIB_model_specific_MSSM_thresholds_SRC)))

LIB_model_specific_MSSM_thresholds_DEP := $(LIB_model_specific_MSSM_thresholds_OBJ:.o=.d)

LIB_model_specific_MSSM_thresholds     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

LIB_model_specific_MSSM_thresholds_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-dep \
		clean-$(MODNAME)-lib clean-$(MODNAME)-obj distclean-$(MODNAME)

all-$(MODNAME): $(LIB_model_specific_MSSM_thresholds)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		$(Q)install -d $(LIB_model_specific_MSSM_thresholds_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIB_model_specific_MSSM_thresholds_SRC) $(LIB_model_specific_MSSM_thresholds_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIB_model_specific_MSSM_thresholds_HDR) $(LIB_model_specific_MSSM_thresholds_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIB_model_specific_MSSM_thresholds_MK) $(LIB_model_specific_MSSM_thresholds_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		$(Q)-rm -f $(LIB_model_specific_MSSM_thresholds_DEP)

clean-$(MODNAME)-lib:
		$(Q)-rm -f $(LIB_model_specific_MSSM_thresholds)

clean-$(MODNAME)-obj:
		$(Q)-rm -f $(LIB_model_specific_MSSM_thresholds_OBJ)

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		@true

distclean-$(MODNAME): clean-$(MODNAME)

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(LIB_model_specific_MSSM_thresholds_DEP) $(LIB_model_specific_MSSM_thresholds_OBJ): CPPFLAGS += $(EIGENFLAGS) $(BOOSTFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIB_model_specific_MSSM_thresholds_DEP) $(LIB_model_specific_MSSM_thresholds_OBJ): CPPFLAGS += $(LOOPFUNCFLAGS)
endif

$(LIB_model_specific_MSSM_thresholds): $(LIB_model_specific_MSSM_thresholds_OBJ)
		@$(MSG)
		$(Q)$(MODULE_MAKE_LIB_CMD) $@ $^

ALLDEP += $(LIB_model_specific_MSSM_thresholds_DEP)
ALLLIB += $(LIB_model_specific_MSSM_thresholds)
