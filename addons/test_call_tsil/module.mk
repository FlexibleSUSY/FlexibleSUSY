DIR          := addons/test_call_tsil
MODNAME      := test_call_tsil
WITH_$(MODNAME) := yes

LIBtest_call_tsil_MK  := $(DIR)/module.mk

# source files
LIBtest_call_tsil_SRC := \
		$(DIR)/call_tsil.cpp

# main()
EXEtest_call_tsil_SRC := \
		$(DIR)/run.cpp

# header files
LIBtest_call_tsil_HDR := \
		$(DIR)/call_tsil.hpp

LIBtest_call_tsil_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBtest_call_tsil_SRC))) \
		$(patsubst %.c, %.o, $(filter %.c, $(LIBtest_call_tsil_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBtest_call_tsil_SRC)))

EXEtest_call_tsil_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXEtest_call_tsil_SRC))) \
		$(patsubst %.c, %.o, $(filter %.c, $(EXEtest_call_tsil_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(EXEtest_call_tsil_SRC)))

LIBtest_call_tsil_DEP := \
		$(LIBtest_call_tsil_OBJ:.o=.d)

EXEtest_call_tsil_DEP := \
		$(EXEtest_call_tsil_OBJ:.o=.d)

EXEtest_call_tsil_EXE := \
		$(patsubst %.cpp, %.x, $(filter %.cpp, $(EXEtest_call_tsil_SRC))) \
		$(patsubst %.f, %.x, $(filter %.f, $(EXEtest_call_tsil_SRC)))

LIBtest_call_tsil     := \
		$(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

LIBtest_call_tsil_INSTALL_DIR := \
		$(INSTALL_DIR)/$(DIR)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-dep \
		clean-$(MODNAME)-lib clean-$(MODNAME)-obj clean-$(MODNAME)-src \
		distclean-$(MODNAME)

all-$(MODNAME): $(LIBtest_call_tsil) $(EXEtest_call_tsil_EXE)
		@true

clean-$(MODNAME)-dep:
		-rm -f $(LIBtest_call_tsil_DEP)
		-rm -f $(EXEtest_call_tsil_DEP)

clean-$(MODNAME)-lib:
		-rm -f $(LIBtest_call_tsil)

clean-$(MODNAME)-obj:
		-rm -f $(LIBtest_call_tsil_OBJ)
		-rm -f $(EXEtest_call_tsil_OBJ)

# BEGIN: NOT EXPORTED ##########################################
# remove generated files
clean-$(MODNAME)-src:
		@true

clean-$(MODNAME): clean-$(MODNAME)-src
# END:   NOT EXPORTED ##########################################

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		-rm -f $(EXEtest_call_tsil_EXE)

distclean-$(MODNAME): clean-$(MODNAME)

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(LIBtest_call_tsil_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBtest_call_tsil_SRC) $(LIBtest_call_tsil_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBtest_call_tsil_HDR) $(LIBtest_call_tsil_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(EXEtest_call_tsil_SRC) $(LIBtest_call_tsil_INSTALL_DIR)
		$(INSTALL_STRIPPED) $(LIBtest_call_tsil_MK) $(LIBtest_call_tsil_INSTALL_DIR) -m u=rw,g=r,o=r
endif

$(LIBtest_call_tsil_DEP) $(EXEtest_call_tsil_DEP) $(LIBtest_call_tsil_OBJ) $(EXEtest_call_tsil_OBJ): CPPFLAGS += $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBtest_call_tsil_DEP) $(EXEtest_call_tsil_DEP) $(LIBtest_call_tsil_OBJ) $(EXEtest_call_tsil_OBJ): CPPFLAGS += $(LOOPFUNCFLAGS)
endif

ifeq ($(ENABLE_TSIL),yes)
$(LIBtest_call_tsil_DEP) $(EXEtest_call_tsil_DEP) $(LIBtest_call_tsil_OBJ) $(EXEtest_call_tsil_OBJ): CPPFLAGS += $(TSILFLAGS)
endif

$(LIBtest_call_tsil): $(LIBtest_call_tsil_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^

$(EXEtest_call_tsil_EXE): $(EXEtest_call_tsil_OBJ) $(LIBtest_call_tsil) $(LIBFLEXI) $(filter-out -%,$(LOOPFUNCLIBS))
		$(CXX) -o $@ $(call abspathx,$^) $(filter -%,$(LOOPFUNCLIBS)) $(GSLLIBS) $(BOOSTTHREADLIBS) $(THREADLIBS) $(TSILLIBS) $(FLIBS)

ALLDEP += $(LIBtest_call_tsil_DEP) $(EXEtest_call_tsil_DEP)
ALLSRC += $(LIBtest_call_tsil_SRC) $(EXEtest_call_tsil_SRC)
ALLLIB += $(LIBtest_call_tsil)
ALLEXE += $(EXEtest_call_tsil_EXE)
