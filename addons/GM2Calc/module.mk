DIR          := addons/GM2Calc
MODNAME      := GM2Calc
WITH_$(MODNAME) := yes

LIBGM2Calc_MK  := $(DIR)/module.mk

# source files
LIBGM2Calc_SRC := \
		$(DIR)/ffunctions.cpp \
		$(DIR)/gm2_1loop.cpp \
		$(DIR)/gm2_2loop.cpp \
		$(DIR)/gm2_mb.cpp \
		$(DIR)/gm2_slha_io.cpp \
		$(DIR)/gm2_uncertainty.cpp \
		$(DIR)/gm2calc_interface.cpp \
		$(DIR)/MSSMNoFV_onshell.cpp \
		$(DIR)/MSSMNoFV_onshell_mass_eigenstates.cpp \
		$(DIR)/MSSMNoFV_onshell_physical.cpp \
		$(DIR)/MSSMNoFV_onshell_problems.cpp \
		$(DIR)/MSSMNoFV_onshell_soft_parameters.cpp \
		$(DIR)/MSSMNoFV_onshell_susy_parameters.cpp

# main()
EXEGM2Calc_SRC := \
		$(DIR)/gm2calc.cpp

# header files
LIBGM2Calc_HDR := \
		$(DIR)/ffunctions.hpp \
		$(DIR)/gm2_1loop.hpp \
		$(DIR)/gm2_2loop.hpp \
		$(DIR)/gm2_error.hpp \
		$(DIR)/gm2_mb.hpp \
		$(DIR)/gm2_slha_io.hpp \
		$(DIR)/gm2_uncertainty.hpp \
		$(DIR)/gm2calc_interface.hpp \
		$(DIR)/MSSMNoFV_onshell.hpp \
		$(DIR)/MSSMNoFV_onshell_mass_eigenstates.hpp \
		$(DIR)/MSSMNoFV_onshell_physical.hpp \
		$(DIR)/MSSMNoFV_onshell_problems.hpp \
		$(DIR)/MSSMNoFV_onshell_soft_parameters.hpp \
		$(DIR)/MSSMNoFV_onshell_susy_parameters.hpp

LIBGM2Calc_SLHA_INPUT := \
		$(DIR)/example.slha \
		$(DIR)/example.gm2

LIBGM2Calc_INFO := \
		$(DIR)/AUTHORS \
		$(DIR)/COPYING \
		$(DIR)/README

LIBGM2Calc_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBGM2Calc_SRC)))

EXEGM2Calc_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXEGM2Calc_SRC)))

LIBGM2Calc_DEP := \
		$(LIBGM2Calc_OBJ:.o=.d)

EXEGM2Calc_DEP := \
		$(EXEGM2Calc_OBJ:.o=.d)

EXEGM2Calc_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXEGM2Calc_SRC)))

LIBGM2Calc     := \
		$(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

LIBGM2Calc_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-dep \
		clean-$(MODNAME)-lib clean-$(MODNAME)-obj distclean-$(MODNAME)

all-$(MODNAME): $(LIBGM2Calc) $(EXEGM2Calc_EXE)
		@true

clean-$(MODNAME)-dep:
		$(Q)-rm -f $(LIBGM2Calc_DEP)
		$(Q)-rm -f $(EXEGM2Calc_DEP)

clean-$(MODNAME)-lib:
		$(Q)-rm -f $(LIBGM2Calc)

clean-$(MODNAME)-obj:
		$(Q)-rm -f $(LIBGM2Calc_OBJ)
		$(Q)-rm -f $(EXEGM2Calc_OBJ)

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		$(Q)-rm -f $(EXEGM2Calc_EXE)

distclean-$(MODNAME): clean-$(MODNAME)

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

ifneq ($(INSTALL_DIR),)
install-src::
		$(Q)install -d $(LIBGM2Calc_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBGM2Calc_SRC) $(LIBGM2Calc_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBGM2Calc_HDR) $(LIBGM2Calc_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(EXEGM2Calc_SRC) $(LIBGM2Calc_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBGM2Calc_MK) $(LIBGM2Calc_INSTALL_DIR)
endif

$(LIBGM2Calc_DEP) $(EXEGM2Calc_DEP) $(LIBGM2Calc_OBJ) $(EXEGM2Calc_OBJ): CPPFLAGS += $(EIGENFLAGS) $(BOOSTFLAGS)

$(LIBGM2Calc): $(LIBGM2Calc_OBJ)
		@$(MSG)
		$(Q)$(MODULE_MAKE_LIB_CMD) $@ $^

$(DIR)/%.x: $(DIR)/%.o $(LIBGM2Calc) $(LIBFLEXI)
		@$(MSG)
		$(Q)$(CXX) -o $@ $(call abspathx,$^)

ALLDEP += $(LIBGM2Calc_DEP) $(EXEGM2Calc_DEP)
ALLSRC += $(LIBGM2Calc_SRC) $(EXEGM2Calc_SRC)
ALLLIB += $(LIBGM2Calc)
ALLEXE += $(EXEGM2Calc_EXE)
