DIR             := test/SOFTSUSY
MODNAME         := SOFTSUSY
WITH_$(MODNAME) := yes
MODSOFTSUSY_MOD := NMSSM_higgs MSSM_higgs
MODSOFTSUSY_DEP := $(patsubst %,model_specific/%,$(MODSOFTSUSY_MOD))
MODSOFTSUSY_INC := $(patsubst %,-Imodel_specific/%,$(MODSOFTSUSY_MOD))
MODSOFTSUSY_LIB := $(foreach M,$(MODSOFTSUSY_MOD),model_specific/$M/libmodel_specific_$M$(MODULE_LIBEXT))

LIBSOFTSUSY_HDR  := \
		$(DIR)/conversion.hpp \
		$(DIR)/def.h \
		$(DIR)/diagonalization.hpp \
		$(DIR)/gut_scale_calculator.hpp \
		$(DIR)/linalg.h \
		$(DIR)/lowe_legacy.h \
		$(DIR)/mycomplex.h \
		$(DIR)/numerics_legacy.h \
		$(DIR)/rge.h \
		$(DIR)/rk_legacy.hpp \
		$(DIR)/utils.h \
		$(DIR)/xpr-base.h \
		$(DIR)/xpr-matrix.h \
		$(DIR)/xpr-vector.h

LIBSOFTSUSY_SRC  := \
		$(DIR)/conversion.cpp \
		$(DIR)/def.cpp \
		$(DIR)/diagonalization.cpp \
		$(DIR)/flavoursoft.cpp \
		$(DIR)/linalg.cpp \
		$(DIR)/lowe_legacy.cpp \
		$(DIR)/mssmUtils.cpp \
		$(DIR)/nmssmUtils.cpp \
		$(DIR)/nmssmsoftpars.cpp \
		$(DIR)/nmssmsoftsusy.cpp \
		$(DIR)/nmssmsusy.cpp \
		$(DIR)/numerics_legacy.cpp \
		$(DIR)/physpars.cpp \
		$(DIR)/rge.cpp \
		$(DIR)/rk_legacy.cpp \
		$(DIR)/susy.cpp \
		$(DIR)/tensor.cpp \
		$(DIR)/utils.cpp

ifneq ($(findstring two_scale,$(SOLVERS)),)
LIBSOFTSUSY_SRC  += \
		$(DIR)/SoftsusyMSSM_two_scale.cpp \
		$(DIR)/SoftsusyMSSM_two_scale_convergence_tester.cpp \
		$(DIR)/SoftsusyMSSM_two_scale_initial_guesser.cpp \
		$(DIR)/SoftsusyMSSM_two_scale_susy_scale_constraint.cpp \
		$(DIR)/SoftsusyMSSM_two_scale_low_scale_constraint.cpp \
		$(DIR)/SoftsusyMSSM_two_scale_sugra_constraint.cpp \
		$(DIR)/SoftsusyNMSSM_two_scale.cpp \
		$(DIR)/SoftsusyNMSSM_two_scale_convergence_tester.cpp \
		$(DIR)/SoftsusyNMSSM_two_scale_initial_guesser.cpp \
		$(DIR)/SoftsusyNMSSM_two_scale_low_scale_constraint.cpp \
		$(DIR)/SoftsusyNMSSM_two_scale_sugra_constraint.cpp \
		$(DIR)/SoftsusyNMSSM_two_scale_susy_scale_constraint.cpp
endif

LIBSOFTSUSY_OBJ  := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBSOFTSUSY_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBSOFTSUSY_SRC)))

LIBSOFTSUSY_DEP  := $(LIBSOFTSUSY_OBJ:.o=.d)

LIBSOFTSUSY      := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

RUN_SOFTSUSY_SRC := $(DIR)/run_softsusy.cpp
RUN_SOFTSUSY_OBJ := $(patsubst %.cpp, %.o, $(filter %.cpp, $(RUN_SOFTSUSY_SRC)))
RUN_SOFTSUSY_DEP := $(RUN_SOFTSUSY_OBJ:.o=.d)
RUN_SOFTSUSY_EXE := $(RUN_SOFTSUSY_OBJ:.o=.x)

RUN_SOFTPOINT_SRC := $(DIR)/run_softpoint.cpp
RUN_SOFTPOINT_OBJ := $(patsubst %.cpp, %.o, $(filter %.cpp, $(RUN_SOFTPOINT_SRC)))
RUN_SOFTPOINT_DEP := $(RUN_SOFTPOINT_OBJ:.o=.d)
RUN_SOFTPOINT_EXE := $(RUN_SOFTPOINT_OBJ:.o=.x)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME) \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj

all-$(MODNAME): $(LIBSOFTSUSY)

clean-$(MODNAME)-dep:
		$(Q)-rm -f $(LIBSOFTSUSY_DEP)
		$(Q)-rm -f $(RUN_SOFTSUSY_DEP)
		$(Q)-rm -f $(RUN_SOFTPOINT_DEP)

clean-$(MODNAME)-lib:
		$(Q)-rm -f $(LIBSOFTSUSY)

clean-$(MODNAME)-obj:
		$(Q)-rm -f $(LIBSOFTSUSY_OBJ)
		$(Q)-rm -f $(RUN_SOFTSUSY_OBJ)
		$(Q)-rm -f $(RUN_SOFTPOINT_OBJ)

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-obj
		$(Q)-rm -f $(LIBSOFTSUSY)
		$(Q)-rm -f $(RUN_SOFTSUSY_EXE)
		$(Q)-rm -f $(RUN_SOFTPOINT_EXE)

distclean-$(MODNAME): clean-$(MODNAME)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(LIBSOFTSUSY_DEP) $(RUN_SOFTSUSY_DEP) $(RUN_SOFTPOINT_DEP) \
$(LIBSOFTSUSY_OBJ) $(RUN_SOFTSUSY_OBJ) $(RUN_SOFTPOINT_OBJ): \
	CPPFLAGS += $(MODSOFTSUSY_INC) $(BOOSTFLAGS) $(EIGENFLAGS)

$(LIBSOFTSUSY): $(LIBSOFTSUSY_OBJ)
		@$(MSG)
		$(Q)$(MODULE_MAKE_LIB_CMD) $@ $^

$(RUN_SOFTSUSY_EXE): $(RUN_SOFTSUSY_OBJ) $(LIBSOFTSUSY) $(MODSOFTSUSY_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		@$(MSG)
		$(Q)$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(FLIBS)

$(RUN_SOFTPOINT_EXE): $(RUN_SOFTPOINT_OBJ) $(LIBSOFTSUSY) $(MODSOFTSUSY_LIB) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		@$(MSG)
		$(Q)$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(FLIBS)

ALLDEP += $(LIBSOFTSUSY_DEP) $(RUN_SOFTSUSY_DEP) $(RUN_SOFTPOINT_DEP)
ALLLIB += $(LIBSOFTSUSY)
ALLEXE += $(RUN_SOFTSUSY_EXE) $(RUN_SOFTPOINT_EXE)
ALLMODDEP += $(MODSOFTSUSY_DEP)
