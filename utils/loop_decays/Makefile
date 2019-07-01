MATH	:= WolframKernel

META_DIR	:= meta
TEMPLATES_DIR	:= templates
OUTPUT_DIR	:= output

META_SRC	:= \
		$(META_DIR)/OneLoopDecaysUtils.m

CALC_AMPLITUDES_SRC := \
		$(META_DIR)/calc_amplitudes_SFF.m \
		$(META_DIR)/calc_amplitudes_SSS.m \
		$(META_DIR)/calc_amplitudes_SSV.m \
		$(META_DIR)/calc_amplitudes_SVV.m

CALC_AMPLITUDES_DIAGS := $(patsubst $(META_DIR)/calc_amplitudes_%.m,$(OUTPUT_DIR)/diagrams_%.m,$(CALC_AMPLITUDES_SRC))
CALC_AMPLITUDES_AMPS := $(patsubst $(META_DIR)/calc_amplitudes_%.m,$(OUTPUT_DIR)/amplitudes_%.m,$(CALC_AMPLITUDES_SRC))
CALC_AMPLITUDES_EXPRS := $(patsubst $(META_DIR)/calc_amplitudes_%.m,$(OUTPUT_DIR)/amplitudes_exprs_%.m,$(CALC_AMPLITUDES_SRC))
CALC_AMPLITUDES_FF := $(patsubst $(META_DIR)/calc_amplitudes_%.m,$(OUTPUT_DIR)/form_factors_%.m,$(CALC_AMPLITUDES_SRC))

CALC_AMPLITUDES_OUTPUT := \
		$(CALC_AMPLITUDES_DIAGS) \
		$(CALC_AMPLITUDES_AMPS) \
		$(CALC_AMPLITUDES_EXPRS) \
		$(CALC_AMPLITUDES_FF)

CALC_AMPLITUDES_LOG := $(patsubst $(META_DIR)/calc_amplitudes_%.m,$(OUTPUT_DIR)/calc_amplitudes_%.log,$(CALC_AMPLITUDES_SRC))

TEMPLATES_SRC := \
		$(TEMPLATES_DIR)/one_loop_decay_diagrams.hpp.in \
		$(TEMPLATES_DIR)/one_loop_decay_diagrams.cpp.in

GENERATE_META	:= $(META_DIR)/generate_loop_decays.m

GENERATE_META_LOG := $(patsubst $(META_DIR)/%.m,$(OUTPUT_DIR)/%.log,$(GENERATE_META))

GENERATED_SRC	:= $(patsubst $(TEMPLATES_DIR)/%.in,$(OUTPUT_DIR)/%,$(TEMPLATES_SRC))

.PHONY: all clean distclean

all: generate-src

clean:
	-rm -f $(CALC_AMPLITUDES_LOG)
	-rm -f $(GENERATE_META_LOG)

distclean: clean
	-rm -f $(GENERATED_SRC)
	-rm -f $(CALC_AMPLITUDES_OUTPUT)
	-rm -rf $(OUTPUT_DIR)

$(OUTPUT_DIR):
	mkdir $(OUTPUT_DIR)

$(CALC_AMPLITUDES_DIAGS): $(OUTPUT_DIR)/diagrams_%.m: $(META_DIR)/calc_amplitudes_%.m $(META_SRC) | $(OUTPUT_DIR)
	"$(MATH)" -run "resultsDir = \"$(OUTPUT_DIR)\"; Get[\"$<\"]; Quit[status];" > "$(OUTPUT_DIR)/calc_amplitudes_$*.log" 2>&1

$(CALC_AMPLITUDES_AMPS): $(OUTPUT_DIR)/amplitudes_%.m: $(META_DIR)/calc_amplitudes_%.m $(META_SRC) | $(OUTPUT_DIR)
	"$(MATH)" -run "resultsDir = \"$(OUTPUT_DIR)\"; Get[\"$<\"]; Quit[status];" > "$(OUTPUT_DIR)/calc_amplitudes_$*.log" 2>&1

$(CALC_AMPLITUDES_EXPRS): $(OUTPUT_DIR)/amplitudes_exprs_%.m: $(META_DIR)/calc_amplitudes_%.m $(META_SRC) | $(OUTPUT_DIR)
	"$(MATH)" -run "resultsDir = \"$(OUTPUT_DIR)\"; Get[\"$<\"]; Quit[status];" > "$(OUTPUT_DIR)/calc_amplitudes_$*.log" 2>&1

generate-src: $(TEMPLATES_SRC) $(GENERATE_META) $(CALC_AMPLITUDES_OUTPUT)
	"$(MATH)" -run "templatesDir = \"$(TEMPLATES_DIR)\"; resultsDir = \"$(OUTPUT_DIR)\"; Get[\"$(GENERATE_META)\"]; Quit[status];" > "$(GENERATE_META_LOG)" 2>&1
