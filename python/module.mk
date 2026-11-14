DIR          := python
MODNAME      := python

PYTHON_PKG_SRC := \
		$(wildcard $(DIR)/flexiblesusy/*.py)

PYTHON_PKG_BACKENDS_SRC := \
		$(wildcard $(DIR)/flexiblesusy/backends/*.py)

PYTHON_PKG_TESTS := \
		$(wildcard $(DIR)/tests/*.py)

PYTHON_PKG_FILES := \
		$(DIR)/README.md \
		$(DIR)/pyproject.toml \
		$(DIR)/module.mk

PYTHON_PKG_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

all-$(MODNAME):
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		$(Q)install -d $(PYTHON_PKG_INSTALL_DIR)
		$(Q)install -d $(PYTHON_PKG_INSTALL_DIR)/flexiblesusy
		$(Q)install -d $(PYTHON_PKG_INSTALL_DIR)/flexiblesusy/backends
		$(Q)install -d $(PYTHON_PKG_INSTALL_DIR)/tests
		$(Q)install -m u=rw,g=r,o=r $(PYTHON_PKG_FILES) $(PYTHON_PKG_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(PYTHON_PKG_SRC) $(PYTHON_PKG_INSTALL_DIR)/flexiblesusy
		$(Q)install -m u=rw,g=r,o=r $(PYTHON_PKG_BACKENDS_SRC) $(PYTHON_PKG_INSTALL_DIR)/flexiblesusy/backends
		$(Q)install -m u=rw,g=r,o=r $(PYTHON_PKG_TESTS) $(PYTHON_PKG_INSTALL_DIR)/tests
endif

clean-$(MODNAME):
		$(Q)-rm -rf $(DIR)/flexiblesusy/__pycache__ $(DIR)/flexiblesusy/backends/__pycache__ \
			$(DIR)/tests/__pycache__

distclean-$(MODNAME): clean-$(MODNAME)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)
