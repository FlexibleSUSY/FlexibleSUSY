DIR          := nlohmann
MODNAME      := nlohmann

NLOHMANN_HDR   := \
		$(DIR)/json.hpp

NLOHMANN_MK    := \
		$(DIR)/module.mk

NLOHMANN_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

all-$(MODNAME):
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		$(Q)install -d $(NLOHMANN_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(NLOHMANN_HDR) $(NLOHMANN_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(NLOHMANN_MK) $(NLOHMANN_INSTALL_DIR)
endif

clean-$(MODNAME):
		@true

distclean-$(MODNAME): clean-$(MODNAME)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)
