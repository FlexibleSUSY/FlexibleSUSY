DIR      := examples
MODNAME  := examples

EXAMPLES_SRC :=

EXAMPLES_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(EXAMPLES_SRC)))

EXAMPLES_DEP := \
		$(EXAMPLES_OBJ:.o=.d)

EXAMPLES_EXE := \
		$(EXAMPLES_OBJ:.o=.x)

STANDALONE_DIR := \
		$(DIR)/standalone-model \
		$(DIR)/standalone-rge \
		$(DIR)/tower

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-dep \
		clean-$(MODNAME)-lib clean-$(MODNAME)-obj distclean-$(MODNAME)

all-$(MODNAME): $(EXAMPLES_EXE)
		@true

clean-$(MODNAME)-dep:
		$(Q)-rm -f $(EXAMPLES_DEP)

clean-$(MODNAME)-lib:
		@true

clean-$(MODNAME)-obj:
		$(Q)-rm -f $(EXAMPLES_OBJ)

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		$(Q)-rm -f $(EXAMPLES_EXE)

distclean-$(MODNAME): clean-$(MODNAME)
		-@for d in $(STANDALONE_DIR); do \
			(cd $$d && make distclean); \
		 done

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

ALLDEP += $(EXAMPLES_DEP)
ALLEXE += $(EXAMPLES_EXE)
