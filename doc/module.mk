include doc/manuals/module.mk
include doc/models/module.mk

DIR          := doc
MODNAME      := doc

DOC_MK       := $(DIR)/module.mk

MODELS_DOC_DIR := $(DIR)/models

DOC_TMPL     := \
		$(DIR)/FlexibleEFTHiggs.rst \
		$(DIR)/librarylink.rst \
		$(DIR)/meta_code.rst \
		$(DIR)/mainpage.dox.in \
		$(DIR)/model_file.rst \
		$(DIR)/slha_input.rst

DOC_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

HTML_OUTPUT_DIR := $(DIR)/html
MAN_OUTPUT_DIR  := $(DIR)/man

IMAGE_DIR       := $(DIR)/images
IMAGES          := $(IMAGE_DIR)/FS-logo.png
INDEX_PAGE      := $(HTML_OUTPUT_DIR)/index.html
MAN_PAGE        := $(MAN_OUTPUT_DIR)/index.html
DOXYFILE        := $(DIR)/Doxyfile
DOXYGEN_MAINPAGE:= $(DIR)/mainpage.dox

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME) \
		$(INDEX_PAGE) $(MAN_PAGE) doc doc-html doc-man doc-pdf

doc: all-$(MODNAME) all-doc-manuals

doc-pdf: doc-manuals-pdf

doc-html: $(INDEX_PAGE)

doc-man: $(MAN_PAGE)

all-$(MODNAME): doc-html doc-man doc-pdf
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		$(Q)install -d $(DOC_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(DOC_TMPL) $(DOC_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(DOC_MK) $(DOC_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(DOXYFILE) $(DOC_INSTALL_DIR)
		$(Q)install -d $(INSTALL_DIR)/$(IMAGE_DIR)
		$(Q)install -m u=rw,g=r,o=r $(IMAGES) $(INSTALL_DIR)/$(IMAGE_DIR)
endif

clean-$(MODNAME): clean-doc-manuals

distclean-$(MODNAME): clean-$(MODNAME) distclean-doc-manuals
		$(Q)-rm -rf $(HTML_OUTPUT_DIR)
		$(Q)-rm -f $(DOXYGEN_MAINPAGE)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(INDEX_PAGE):
		( cat $(DOXYFILE) ; \
		  echo "INPUT = $(MODULES) $(README_FILE)" ; \
		  echo "OUTPUT_DIRECTORY = $(HTML_OUTPUT_DIR)" ; \
		  echo "EXCLUDE = $(ALLDEP) $(META_SRC) $(TEMPLATES) \
		        $(TEST_SRC) $(TEST_META)"; \
		  echo "EXCLUDE_PATTERNS = */meta/* */test/*"; \
		  echo "IMAGE_PATH = $(IMAGE_DIR)"; \
		  echo "INCLUDE_PATH = $(MODULES)"; \
		) | doxygen -

$(MAN_PAGE):
		( cat $(DOXYFILE) ; \
		  echo "INPUT = $(MODULES) $(README_FILE)" ; \
		  echo "OUTPUT_DIRECTORY = $(MAN_OUTPUT_DIR)" ; \
		  echo "EXCLUDE = $(ALLDEP) $(META_SRC) $(TEMPLATES) \
		        $(TEST_SRC) $(TEST_META)"; \
		  echo "EXCLUDE_PATTERNS = */meta/* */test/*"; \
		  echo "IMAGE_PATH = $(IMAGE_DIR)"; \
		  echo "INCLUDE_PATH = $(MODULES)"; \
		  echo "GENERATE_MAN = YES"; \
		  echo "GENERATE_HTML = NO"; \
		) | doxygen -
