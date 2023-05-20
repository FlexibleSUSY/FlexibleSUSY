DIR          := models/fmssmn
MODNAME      := fmssmn

LIBFMSSMN_SRC  :=
LIBFMSSMN_GENERATED_SRC :=
LIBFMSSMN_INC  :=

ifneq ($(findstring lattice,$(SOLVERS)),)
LIBFMSSMN_GENERATED_SRC += \
		$(DIR)/fmssmn_lattice_rge.f \
		$(DIR)/fmssmn_lattice_constraints.f \
		$(DIR)/fmssmn_lattice_numerical_constraints_functions.f \
		$(DIR)/fmssmn_lattice_numerical_constraints_dependence.cpp \
		$(DIR)/fmssm_fmssmn_lattice_matchings.f \
		$(DIR)/fmssm_fmssmn_lattice_numerical_matchings_functions.f \
		$(DIR)/fmssm_fmssmn_lattice_numerical_matchings_dependence.cpp

LIBFMSSMN_INC  += \
		$(DIR)/fmssmn_lattice_translator.inc

LIBFMSSMN_SRC  += \
		$(DIR)/fmssmn_lattice.cpp \
		$(DIR)/fmssmn_lattice_mz_constraint.cpp \
		$(DIR)/fmssmn_lattice_msusy_constraint.cpp \
		$(DIR)/fmssmn_lattice_mx_constraint.cpp \
		$(DIR)/fmssmn_lattice_numerical_constraints.cpp \
		$(DIR)/fmssm_fmssmn_lattice_numerical_matchings.cpp \
		$(LIBFMSSMN_GENERATED_SRC)
endif

LIBFMSSMN_OBJ  := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBFMSSMN_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBFMSSMN_SRC)))

LIBFMSSMN_DEP  := \
		$(LIBFMSSMN_OBJ:.o=.d)

LIBFMSSMN      := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) \
		clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj \
		clean-$(MODNAME)-src distclean-$(MODNAME)

all-$(MODNAME): $(LIBFMSSMN)

clean-$(MODNAME)-dep:
		-rm -f $(LIBFMSSMN_DEP)

clean-$(MODNAME)-lib:
		-rm -f $(LIBFMSSMN)

clean-$(MODNAME)-obj:
		-rm -f $(LIBFMSSMN_OBJ)

clean-$(MODNAME)-src:
		-rm -f $(LIBFMSSMN_GENERATED_SRC)
		-rm -f $(LIBFMSSMN_INC)

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		@true

distclean-$(MODNAME): clean-$(MODNAME) clean-$(MODNAME)-src
		@true

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(DIR)/fmssmn_lattice.o: $(DIR)/fmssmn_lattice_translator.inc

$(DIR)/%.cpp : $(DIR)/%.cpp.m
	$(MATH) -run 'filename="$@"; << $<; Quit[]'

$(DIR)/%.f : $(DIR)/%.f.m
	$(MATH) -run 'filename="$@"; << $<; Quit[]'

$(DIR)/%.inc : $(DIR)/%.inc.m
	$(MATH) -run 'filename="$@"; << $<; Quit[]'

ifneq ($(findstring lattice,$(SOLVERS)),)
$(LIBFMSSMN_DEP) $(LIBFMSSMN_OBJ): CPPFLAGS += $(EIGENFLAGS) $(GSLFLAGS) $(BOOSTFLAGS)
endif

$(LIBFMSSMN): $(LIBFMSSMN_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^

ALLDEP += $(LIBFMSSMN_DEP)
ALLLIB += $(LIBFMSSMN)
