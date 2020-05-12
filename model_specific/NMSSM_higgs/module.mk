DIR          := model_specific/NMSSM_higgs
MODNAME      := model_specific_NMSSM_higgs
WITH_$(MODNAME) := yes
MODmodel_specific_NMSSM_higgs_MOD := MSSM_higgs
MODmodel_specific_NMSSM_higgs_DEP := $(patsubst %,model_specific/%,$(MODmodel_specific_NMSSM_higgs_MOD))
MODmodel_specific_NMSSM_higgs_INC := $(patsubst %,-Imodel_specific/%,$(MODmodel_specific_NMSSM_higgs_MOD))
MODmodel_specific_NMSSM_higgs_LIB := $(foreach M,$(MODmodel_specific_NMSSM_higgs_MOD),model_specific/$M/libmodel_specific_$M$(MODULE_LIBEXT))

LIB_model_specific_NMSSM_higgs_MK  := \
		$(DIR)/module.mk

LIB_model_specific_NMSSM_higgs_SRC := \
		$(DIR)/nmssm_twoloophiggs.cpp \
		$(DIR)/nmssm2loop.f

LIB_model_specific_NMSSM_higgs_HDR := \
		$(DIR)/nmssm_twoloophiggs.hpp \
		$(DIR)/nmssm2loop.h

LIB_model_specific_NMSSM_higgs_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIB_model_specific_NMSSM_higgs_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIB_model_specific_NMSSM_higgs_SRC)))

LIB_model_specific_NMSSM_higgs_DEP := $(LIB_model_specific_NMSSM_higgs_OBJ:.o=.d)

LIB_model_specific_NMSSM_higgs     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

LIB_model_specific_NMSSM_higgs_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-dep \
		clean-$(MODNAME)-lib clean-$(MODNAME)-obj distclean-$(MODNAME)

all-$(MODNAME): $(LIB_model_specific_NMSSM_higgs)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		$(Q)install -d $(LIB_model_specific_NMSSM_higgs_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIB_model_specific_NMSSM_higgs_SRC) $(LIB_model_specific_NMSSM_higgs_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIB_model_specific_NMSSM_higgs_HDR) $(LIB_model_specific_NMSSM_higgs_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIB_model_specific_NMSSM_higgs_MK) $(LIB_model_specific_NMSSM_higgs_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		$(Q)-rm -f $(LIB_model_specific_NMSSM_higgs_DEP)

clean-$(MODNAME)-lib:
		$(Q)-rm -f $(LIB_model_specific_NMSSM_higgs)

clean-$(MODNAME)-obj:
		$(Q)-rm -f $(LIB_model_specific_NMSSM_higgs_OBJ)

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		@true

distclean-$(MODNAME): clean-$(MODNAME)

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(LIB_model_specific_NMSSM_higgs_DEP) $(LIB_model_specific_NMSSM_higgs_OBJ): \
	CPPFLAGS += $(MODmodel_specific_NMSSM_higgs_INC) $(EIGENFLAGS) $(BOOSTFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIB_model_specific_NMSSM_higgs_DEP) $(LIB_model_specific_NMSSM_higgs_OBJ): \
	CPPFLAGS += $(LOOPFUNCFLAGS)
endif

ifeq ($(ENABLE_SHARED_LIBS),yes)
$(LIB_model_specific_NMSSM_higgs): $(LIB_model_specific_NMSSM_higgs_OBJ)
		@$(MSG)
		$(Q)$(MODULE_MAKE_LIB_CMD) $@ $^ $(FLIBS)
else
$(LIB_model_specific_NMSSM_higgs): $(LIB_model_specific_NMSSM_higgs_OBJ)
		@$(MSG)
		$(Q)$(MODULE_MAKE_LIB_CMD) $@ $^
endif

ALLDEP += $(LIB_model_specific_NMSSM_higgs_DEP)
ALLLIB += $(LIB_model_specific_NMSSM_higgs)
ALLMODDEP += $(MODmodel_specific_NMSSM_higgs_DEP)
