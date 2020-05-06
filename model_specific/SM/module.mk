DIR          := model_specific/SM
MODNAME      := model_specific_SM
WITH_$(MODNAME) := yes

LIB_model_specific_SM_MK  := \
		$(DIR)/module.mk

LIB_model_specific_SM_SRC := \
		$(DIR)/sm_fourloophiggs.cpp \
		$(DIR)/sm_fourloop_as.cpp \
		$(DIR)/sm_threeloophiggs.cpp \
		$(DIR)/sm_twoloophiggs.cpp \
		$(DIR)/standard_model.cpp \
		$(DIR)/standard_model_effective_couplings.cpp \
		$(DIR)/standard_model_physical.cpp \
		$(DIR)/standard_model_two_scale_convergence_tester.cpp \
		$(DIR)/standard_model_two_scale_low_scale_constraint.cpp \
		$(DIR)/standard_model_two_scale_model.cpp \
		$(DIR)/weinberg_angle.cpp

LIB_model_specific_SM_HDR := \
		$(DIR)/sm_fourloophiggs.hpp \
		$(DIR)/sm_fourloop_as.hpp \
		$(DIR)/sm_threeloophiggs.hpp \
		$(DIR)/sm_twoloophiggs.hpp \
		$(DIR)/standard_model.hpp \
		$(DIR)/standard_model_convergence_tester.hpp \
		$(DIR)/standard_model_effective_couplings.hpp \
		$(DIR)/standard_model_low_scale_constraint.hpp \
		$(DIR)/standard_model_physical.hpp \
		$(DIR)/standard_model_two_scale_convergence_tester.hpp \
		$(DIR)/standard_model_two_scale_low_scale_constraint.hpp \
		$(DIR)/standard_model_two_scale_model.hpp \
		$(DIR)/weinberg_angle.hpp

LIB_model_specific_SM_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIB_model_specific_SM_SRC)))

LIB_model_specific_SM_DEP := $(LIB_model_specific_SM_OBJ:.o=.d)

LIB_model_specific_SM     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

LIB_model_specific_SM_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-dep \
		clean-$(MODNAME)-lib clean-$(MODNAME)-obj distclean-$(MODNAME)

all-$(MODNAME): $(LIB_model_specific_SM)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		$(Q)install -d $(LIB_model_specific_SM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIB_model_specific_SM_SRC) $(LIB_model_specific_SM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIB_model_specific_SM_HDR) $(LIB_model_specific_SM_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIB_model_specific_SM_MK) $(LIB_model_specific_SM_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		$(Q)-rm -f $(LIB_model_specific_SM_DEP)

clean-$(MODNAME)-lib:
		$(Q)-rm -f $(LIB_model_specific_SM)

clean-$(MODNAME)-obj:
		$(Q)-rm -f $(LIB_model_specific_SM_OBJ)

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		@true

distclean-$(MODNAME): clean-$(MODNAME)

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(LIB_model_specific_SM_DEP) $(LIB_model_specific_SM_OBJ): CPPFLAGS += $(EIGENFLAGS) $(BOOSTFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIB_model_specific_SM_DEP) $(LIB_model_specific_SM_OBJ): CPPFLAGS += $(LOOPFUNCFLAGS)
endif

ifeq ($(ENABLE_SHARED_LIBS),yes)
$(LIB_model_specific_SM): $(LIB_model_specific_SM_OBJ)
		@$(MSG)
		$(Q)$(MODULE_MAKE_LIB_CMD) $@ $^ $(THREADLIBS)
else
$(LIB_model_specific_SM): $(LIB_model_specific_SM_OBJ)
		@$(MSG)
		$(Q)$(MODULE_MAKE_LIB_CMD) $@ $^
endif

ALLDEP += $(LIB_model_specific_SM_DEP)
ALLLIB += $(LIB_model_specific_SM)
