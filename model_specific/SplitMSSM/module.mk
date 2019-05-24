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
		install -d $(LIB_model_specific_SplitMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIB_model_specific_SplitMSSM_SRC) $(LIB_model_specific_SplitMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIB_model_specific_SplitMSSM_HDR) $(LIB_model_specific_SplitMSSM_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIB_model_specific_SplitMSSM_MK) $(LIB_model_specific_SplitMSSM_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIB_model_specific_SplitMSSM_DEP)

clean-$(MODNAME)-lib:
		-rm -f $(LIB_model_specific_SplitMSSM)

clean-$(MODNAME)-obj:
		-rm -f $(LIB_model_specific_SplitMSSM_OBJ)

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		@true

distclean-$(MODNAME): clean-$(MODNAME)

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(LIB_model_specific_SplitMSSM_DEP) $(LIB_model_specific_SplitMSSM_OBJ): CPPFLAGS += $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(SQLITEFLAGS) $(TSILFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIB_model_specific_SplitMSSM_DEP) $(LIB_model_specific_SplitMSSM_OBJ): CPPFLAGS += $(LOOPFUNCFLAGS)
endif

ifeq ($(ENABLE_SHARED_LIBS),yes)
$(LIB_model_specific_SplitMSSM): $(LIB_model_specific_SplitMSSM_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^ $(BOOSTTHREADLIBS) $(GSLLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS)
else
$(LIB_model_specific_SplitMSSM): $(LIB_model_specific_SplitMSSM_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^
endif

ALLDEP += $(LIB_model_specific_SplitMSSM_DEP)
ALLLIB += $(LIB_model_specific_SplitMSSM)
