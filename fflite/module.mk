DIR          := fflite
MODNAME      := fflite
WITH_$(MODNAME) := yes

LIBFFLITE_HDR := \
		$(DIR)/defs.h \
		$(DIR)/externals.h \
		$(DIR)/ff.h \
		$(DIR)/fferr.h \
		$(DIR)/fflite.hpp \
		$(DIR)/ffwarn.h \
		$(DIR)/lt.h \
		$(DIR)/types.h

LIBFFLITE_MK := \
		$(DIR)/module.mk

LIBFFLITE_SRC := \
		$(DIR)/BcoeffAD.F \
		$(DIR)/ffca0.F \
		$(DIR)/ffcb0.F \
		$(DIR)/ffcb1.F \
		$(DIR)/ffcb2p.F \
		$(DIR)/ffcc0.F \
		$(DIR)/ffcli2.F \
		$(DIR)/ffinit.F \
		$(DIR)/ffxa0.F \
		$(DIR)/ffxb0.F \
		$(DIR)/ffxb1.F \
		$(DIR)/ffxb2p.F \
		$(DIR)/ffxli2.F \
		$(DIR)/ini.F

LIBFFLITE_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBFFLITE_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBFFLITE_SRC))) \
		$(patsubst %.F, %.o, $(filter %.F, $(LIBFFLITE_SRC)))

LIBFFLITE_DEP := \
		$(LIBFFLITE_OBJ:.o=.d)

LIBFFLITE     := $(DIR)/lib$(MODNAME)$(MODULE_LIBEXT)

LIBFFLITE_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-dep \
		clean-$(MODNAME)-lib clean-$(MODNAME)-obj distclean-$(MODNAME)

all-$(MODNAME): $(LIBFFLITE)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		$(Q)install -d $(LIBFFLITE_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBFFLITE_SRC) $(LIBFFLITE_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBFFLITE_HDR) $(LIBFFLITE_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBFFLITE_MK) $(LIBFFLITE_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		$(Q)-rm -f $(LIBFFLITE_DEP)

clean-$(MODNAME)-lib:
		$(Q)-rm -f $(LIBFFLITE)

clean-$(MODNAME)-obj:
		$(Q)-rm -f $(LIBFFLITE_OBJ)

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		@true

distclean-$(MODNAME): clean-$(MODNAME)

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(LIBFFLITE): $(LIBFFLITE_OBJ)
		@$(MSG)
		$(Q)$(MODULE_MAKE_LIB_CMD) $@ $^

ifeq ($(ENABLE_FFLITE),yes)
ALLDEP += $(LIBFFLITE_DEP)
ALLLIB += $(LIBFFLITE)
endif
