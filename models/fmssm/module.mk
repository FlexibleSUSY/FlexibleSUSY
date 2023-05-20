DIR          := models/fmssm
MODNAME      := fmssm

LIBFMSSM_SRC  :=
LIBFMSSM_GENERATED_SRC :=
LIBFMSSM_INC  :=

ifneq ($(findstring lattice,$(SOLVERS)),)
LIBFMSSM_GENERATED_SRC += \
		$(DIR)/fmssm_lattice_rge.f \
		$(DIR)/fmssm_lattice_constraints.f \
		$(DIR)/fmssm_lattice_numerical_constraints_functions.f \
		$(DIR)/fmssm_lattice_numerical_constraints_dependence.cpp

LIBFMSSM_INC  += \
		$(DIR)/fmssm_lattice_translator.inc

LIBFMSSM_SRC  += \
		$(DIR)/fmssm_lattice.cpp \
		$(DIR)/fmssm_lattice_mz_constraint.cpp \
		$(DIR)/fmssm_lattice_numerical_mz_constraint.cpp \
		$(DIR)/fmssm_lattice_msusy_constraint.cpp \
		$(DIR)/fmssm_lattice_numerical_msusy_constraint.cpp \
		$(DIR)/fmssm_lattice_mx_constraint.cpp \
		$(DIR)/fmssm_oneloop.cpp \
		$(DIR)/fmssm_lattice_numerical_constraints.cpp \
		$(LIBFMSSM_GENERATED_SRC)
endif

LIBFMSSM_OBJ  := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBFMSSM_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBFMSSM_SRC)))

LIBFMSSM_DEP  := \
		$(LIBFMSSM_OBJ:.o=.d)

LIBFMSSM      := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj \
		clean-$(MODNAME)-src distclean-$(MODNAME)

all-$(MODNAME): $(LIBFMSSM)

clean-$(MODNAME)-dep:
		-rm -f $(LIBFMSSM_DEP)

clean-$(MODNAME)-lib:
		-rm -f $(LIBFMSSM)

clean-$(MODNAME)-obj:
		-rm -f $(LIBFMSSM_OBJ)

clean-$(MODNAME)-src:
		-rm -f $(LIBFMSSM_GENERATED_SRC)
		-rm -f $(LIBFMSSM_INC)

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		@true

distclean-$(MODNAME): clean-$(MODNAME) clean-$(MODNAME)-src
		@true

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(DIR)/fmssm_lattice.o: $(DIR)/fmssm_lattice_translator.inc

$(DIR)/%.cpp : $(DIR)/%.cpp.m
	$(MATH) -run 'filename="$@"; << $<; Quit[]'

$(DIR)/%.f : $(DIR)/%.f.m
	$(MATH) -run 'filename="$@"; << $<; Quit[]'

$(DIR)/%.inc : $(DIR)/%.inc.m
	$(MATH) -run 'filename="$@"; << $<; Quit[]'

ifneq ($(findstring lattice,$(SOLVERS)),)
$(LIBFMSSM_DEP) $(LIBFMSSM_OBJ): CPPFLAGS += $(EIGENFLAGS) $(GSLFLAGS) $(BOOSTFLAGS)
endif

$(LIBFMSSM): $(LIBFMSSM_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^

ALLDEP += $(LIBFMSSM_DEP)
ALLLIB += $(LIBFMSSM)
