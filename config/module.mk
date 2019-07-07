DIR          := config
MODNAME      := config

CONFIG_HDR := \
		$(DIR)/config.h

CONFIG_MK    := \
		$(DIR)/module.mk

CONFIG_TMPL  := \
		$(DIR)/abspathx.mk \
		$(DIR)/config.h.in \
		$(DIR)/flexiblesusy-config.in \
		$(DIR)/list_sarah_model_files.sh.in

DEPGEN_SRC   := $(DIR)/depgen.cpp

DEPGEN_OBJ   := $(DEPGEN_SRC:.cpp=.o)

DEPGEN_EXE   := $(DEPGEN_SRC:.cpp=.x)

CONFIG_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

MAKEFILE_IN  := \
		$(DIR)/Makefile.in

REQUIRED_SARAH_VERSION_FILE := \
		$(DIR)/required_sarah_version.m

FLEXIBLESUSY_VERSION_FILE := \
		$(DIR)/version

FLEXIBLESUSY_GIT_COMMIT_FILE := \
		$(DIR)/git_commit

REMOVE_EXPORT_MARKERS := \
		$(DIR)/remove_export_markers.sh

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-dep \
		clean-$(MODNAME)-lib clean-$(MODNAME)-obj distclean-$(MODNAME)

all-$(MODNAME):
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		$(Q)install -d $(CONFIG_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(CONFIG_TMPL) $(CONFIG_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(CONFIG_MK) $(CONFIG_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(DEPGEN_SRC) $(CONFIG_INSTALL_DIR)
		$(Q)install -m u=rwx,g=r,o=r $(MATH_INC_PATHS) $(CONFIG_INSTALL_DIR)
		$(Q)$(INSTALL_STRIPPED) $(MAKEFILE_IN) $(CONFIG_INSTALL_DIR) -m u=rw,g=r,o=r
		$(Q)$(INSTALL_STRIPPED) $(REMOVE_EXPORT_MARKERS) $(CONFIG_INSTALL_DIR) -m u=rwx,g=r,o=r
		$(Q)$(INSTALL_STRIPPED) $(INSTALL_STRIPPED) $(CONFIG_INSTALL_DIR) -m u=rwx,g=r,o=r
		$(Q)$(INSTALL_STRIPPED) $(CONVERT_DOS_PATHS) $(CONFIG_INSTALL_DIR) -m u=rwx,g=r,o=r
endif

clean-$(MODNAME)-dep:
		@true

clean-$(MODNAME)-lib:
		@true

clean-$(MODNAME)-obj:
		$(Q)-rm -f $(DEPGEN_OBJ)

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		$(Q)-rm -f $(DEPGEN_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		$(Q)-rm -f $(CONFIG_HDR)
		$(Q)-rm -f $(FLEXIBLESUSY_VERSION_FILE)
		$(Q)-rm -f $(FLEXIBLESUSY_GIT_COMMIT_FILE)
		$(Q)-rm -f $(REQUIRED_SARAH_VERSION_FILE)

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(DEPGEN_EXE): $(DEPGEN_OBJ)
		@$(MSG)
		$(Q)$(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS)

ALLEXE += $(DEPGEN_EXE)
