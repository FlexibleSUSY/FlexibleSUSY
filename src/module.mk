DIR          := src
MODNAME      := src
WITH_$(MODNAME) := yes

LIBFLEXI_MK  := \
		$(DIR)/module.mk

LIBFLEXI_SRC := \
		$(DIR)/betafunction.cpp \
		$(DIR)/build_info.cpp \
		$(DIR)/bvp_solver_problems.cpp \
		$(DIR)/ckm.cpp \
		$(DIR)/command_line_options.cpp \
		$(DIR)/composite_convergence_tester.cpp \
		$(DIR)/database.cpp \
		$(DIR)/dilog.cpp \
		$(DIR)/dilogc.f \
		$(DIR)/effective_couplings.cpp \
		$(DIR)/global_thread_pool.cpp \
		$(DIR)/gsl_utils.cpp \
		$(DIR)/gsl_vector.cpp \
		$(DIR)/logger.cpp \
		$(DIR)/lowe.cpp \
		$(DIR)/sfermions.cpp \
		$(DIR)/numerics.cpp \
		$(DIR)/numerics2.cpp \
		$(DIR)/physical_input.cpp \
		$(DIR)/pmns.cpp \
		$(DIR)/problems.cpp \
		$(DIR)/pv.cpp \
		$(DIR)/rkf_integrator.cpp \
		$(DIR)/scan.cpp \
		$(DIR)/slha_io.cpp \
		$(DIR)/spectrum_generator_problems.cpp \
		$(DIR)/spectrum_generator_settings.cpp \
		$(DIR)/string_utils.cpp \
		$(DIR)/threshold_corrections.cpp \
		$(DIR)/threshold_loop_functions.cpp \
		$(DIR)/weinberg_angle.cpp \
		$(DIR)/wrappers.cpp

LIBFLEXI_HDR := \
		$(DIR)/array_view.hpp \
		$(DIR)/basic_rk_integrator.hpp \
		$(DIR)/betafunction.hpp \
		$(DIR)/build_info.hpp \
		$(DIR)/bvp_solver_problems.hpp \
		$(DIR)/cextensions.hpp \
		$(DIR)/ckm.hpp \
		$(DIR)/command_line_options.hpp \
		$(DIR)/composite_convergence_tester.hpp \
		$(DIR)/compound_constraint.hpp \
		$(DIR)/concatenate.hpp \
		$(DIR)/constraint.hpp \
		$(DIR)/convergence_tester.hpp \
		$(DIR)/convergence_tester_drbar.hpp \
		$(DIR)/coupling_monitor.hpp \
		$(DIR)/database.hpp \
		$(DIR)/derivative.hpp \
		$(DIR)/dilog.hpp \
		$(DIR)/effective_couplings.hpp \
		$(DIR)/eigen_utils.hpp \
		$(DIR)/eigen_tensor.hpp \
		$(DIR)/error.hpp \
		$(DIR)/ew_input.hpp \
		$(DIR)/ewsb_solver.hpp \
		$(DIR)/find_if.hpp \
		$(DIR)/fixed_point_iterator.hpp \
		$(DIR)/functors.hpp \
		$(DIR)/global_thread_pool.hpp \
		$(DIR)/gsl.hpp \
		$(DIR)/gsl_utils.hpp \
		$(DIR)/gsl_vector.hpp \
		$(DIR)/loop_corrections.hpp \
		$(DIR)/if.hpp \
		$(DIR)/initial_guesser.hpp \
		$(DIR)/linalg2.hpp \
		$(DIR)/logger.hpp \
		$(DIR)/lowe.h \
		$(DIR)/mathlink_utils.hpp \
		$(DIR)/minimizer.hpp \
		$(DIR)/model.hpp \
		$(DIR)/multiindex.hpp \
		$(DIR)/numerics.h \
		$(DIR)/numerics2.hpp \
		$(DIR)/physical_input.hpp \
		$(DIR)/pmns.hpp \
		$(DIR)/pp_map.hpp \
		$(DIR)/problems.hpp \
		$(DIR)/pv.hpp \
		$(DIR)/raii.hpp \
		$(DIR)/rg_flow.hpp \
		$(DIR)/rk.hpp \
		$(DIR)/rkf_integrator.hpp \
		$(DIR)/root_finder.hpp \
		$(DIR)/scan.hpp \
		$(DIR)/sfermions.hpp \
		$(DIR)/single_scale_constraint.hpp \
		$(DIR)/single_scale_matching.hpp \
		$(DIR)/slha_io.hpp \
		$(DIR)/spectrum_generator_problems.hpp \
		$(DIR)/spectrum_generator_settings.hpp \
		$(DIR)/string_utils.hpp \
		$(DIR)/sum.hpp \
		$(DIR)/thread_pool.hpp \
		$(DIR)/threshold_corrections.hpp \
		$(DIR)/threshold_loop_functions.hpp \
		$(DIR)/weinberg_angle.hpp \
		$(DIR)/which.hpp \
		$(DIR)/wrappers.hpp

ifneq ($(findstring two_scale,$(SOLVERS)),)
LIBFLEXI_SRC += \
		$(DIR)/two_scale_running_precision.cpp \
		$(DIR)/two_scale_solver.cpp

LIBFLEXI_HDR += \
		$(DIR)/two_scale_running_precision.hpp \
		$(DIR)/two_scale_solver.hpp
endif

ifneq ($(findstring semi_analytic,$(SOLVERS)),)
LIBFLEXI_SRC += \
		$(DIR)/semi_analytic_solver.cpp \
		$(DIR)/two_scale_running_precision.cpp \
		$(DIR)/two_scale_solver.cpp

LIBFLEXI_HDR += \
		$(DIR)/semi_analytic_solver.hpp \
		$(DIR)/two_scale_running_precision.hpp \
		$(DIR)/two_scale_solver.hpp
endif

# remove duplicates in case multiple solvers are used
LIBFLEXI_SRC := $(sort $(LIBFLEXI_SRC))
LIBFLEXI_HDR := $(sort $(LIBFLEXI_HDR))

LIBFLEXI_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBFLEXI_SRC))) \
		$(patsubst %.c, %.o, $(filter %.c, $(LIBFLEXI_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBFLEXI_SRC)))

LIBFLEXI_DEP := \
		$(LIBFLEXI_OBJ:.o=.d)

LIBFLEXI     := $(DIR)/libflexisusy$(MODULE_LIBEXT)

LIBFLEXI_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-dep \
		clean-$(MODNAME)-lib clean-$(MODNAME)-obj distclean-$(MODNAME)

all-$(MODNAME): $(LIBFLEXI)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		install -d $(LIBFLEXI_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBFLEXI_SRC) $(LIBFLEXI_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBFLEXI_HDR) $(LIBFLEXI_INSTALL_DIR)
		install -m u=rw,g=r,o=r $(LIBFLEXI_MK) $(LIBFLEXI_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		-rm -f $(LIBFLEXI_DEP)

clean-$(MODNAME)-lib:
		-rm -f $(LIBFLEXI)

clean-$(MODNAME)-obj:
		-rm -f $(LIBFLEXI_OBJ)

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj
		@true

distclean-$(MODNAME): clean-$(MODNAME)

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(LIBFLEXI_DEP) $(LIBFLEXI_OBJ): CPPFLAGS += $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(SQLITEFLAGS) $(TSILFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBFLEXI_DEP) $(LIBFLEXI_OBJ): CPPFLAGS += $(LOOPFUNCFLAGS)
endif

ifeq ($(ENABLE_SHARED_LIBS),yes)
$(LIBFLEXI): $(LIBFLEXI_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^ $(BOOSTTHREADLIBS) $(GSLLIBS) $(FLIBS) $(SQLITELIBS) $(TSILLIBS) $(THREADLIBS)
else
$(LIBFLEXI): $(LIBFLEXI_OBJ)
		$(MODULE_MAKE_LIB_CMD) $@ $^
endif

ALLDEP += $(LIBFLEXI_DEP)
ALLLIB += $(LIBFLEXI)
