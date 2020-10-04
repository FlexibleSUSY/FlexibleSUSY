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
		$(DIR)/coupling_monitor.cpp \
		$(DIR)/database.cpp \
		$(DIR)/decays/decay.cpp \
		$(DIR)/decays/decay_amplitudes.cpp \
		$(DIR)/decays/decay_functions.cpp \
		$(DIR)/dilog.cpp \
		$(DIR)/effective_couplings.cpp \
		$(DIR)/global_thread_pool.cpp \
		$(DIR)/gm2calc_interface.cpp \
		$(DIR)/gsl_multimin_fminimizer.cpp \
		$(DIR)/gsl_multiroot_fsolver.cpp \
		$(DIR)/gsl_utils.cpp \
		$(DIR)/gsl_vector.cpp \
		$(DIR)/Li4.cpp \
		$(DIR)/logger.cpp \
		$(DIR)/loop_libraries/library_softsusy.cpp \
		$(DIR)/lowe.cpp \
		$(DIR)/sfermions.cpp \
		$(DIR)/mixings.cpp \
		$(DIR)/model.cpp \
		$(DIR)/numerics.cpp \
		$(DIR)/numerics2.cpp \
		$(DIR)/decays/one_loop_decay_diagrams.cpp \
		$(DIR)/observables.cpp \
		$(DIR)/observable_problems.cpp \
		$(DIR)/observable_problems_format.cpp \
		$(DIR)/physical_input.cpp \
		$(DIR)/pmns.cpp \
		$(DIR)/problems.cpp \
		$(DIR)/rkf_integrator.cpp \
		$(DIR)/scan.cpp \
		$(DIR)/slha_format.cpp \
		$(DIR)/slha_io.cpp \
		$(DIR)/spectrum_generator_problems.cpp \
		$(DIR)/spectrum_generator_settings.cpp \
		$(DIR)/string_conversion.cpp \
		$(DIR)/string_format.cpp \
		$(DIR)/string_utils.cpp \
		$(DIR)/threshold_corrections.cpp \
		$(DIR)/threshold_loop_functions.cpp \
		$(DIR)/trilog.cpp \
		$(DIR)/wrappers.cpp

LIBFLEXI_HDR := \
		$(DIR)/array_view.hpp \
		$(DIR)/basic_rk_integrator.hpp \
		$(DIR)/betafunction.hpp \
		$(DIR)/build_info.hpp \
		$(DIR)/bvp_solver_problems.hpp \
		$(DIR)/bvp_solver_problems_format_mathlink.hpp \
		$(DIR)/cextensions.hpp \
		$(DIR)/ckm.hpp \
		$(DIR)/command_line_options.hpp \
		$(DIR)/complex.hpp \
		$(DIR)/composite_convergence_tester.hpp \
		$(DIR)/compound_constraint.hpp \
		$(DIR)/concatenate.hpp \
		$(DIR)/constraint.hpp \
		$(DIR)/convergence_tester.hpp \
		$(DIR)/convergence_tester_drbar.hpp \
		$(DIR)/coupling_monitor.hpp \
		$(DIR)/database.hpp \
		$(DIR)/decays/decay.hpp \
		$(DIR)/decays/decay_amplitudes.hpp \
		$(DIR)/decays/decays_corrections.hpp \
		$(DIR)/decays/decay_loop_functions.hpp \
		$(DIR)/decays/decays_problems.hpp \
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
		$(DIR)/for_each.hpp \
		$(DIR)/functors.hpp \
		$(DIR)/global_thread_pool.hpp \
		$(DIR)/gm2calc_interface.hpp \
		$(DIR)/gsl.hpp \
		$(DIR)/gsl_multimin_fminimizer.hpp \
		$(DIR)/gsl_multiroot_fsolver.hpp \
		$(DIR)/gsl_utils.hpp \
		$(DIR)/gsl_vector.hpp \
		$(DIR)/Li4.hpp \
		$(DIR)/loop_corrections.hpp \
		$(DIR)/if.hpp \
		$(DIR)/initial_guesser.hpp \
		$(DIR)/linalg2.hpp \
		$(DIR)/logger.hpp \
		$(DIR)/lowe.h \
		$(DIR)/loop_libraries/library_softsusy.hpp \
		$(DIR)/mathlink_utils.hpp \
		$(DIR)/minimizer.hpp \
		$(DIR)/mixings.hpp \
		$(DIR)/model.hpp \
		$(DIR)/multiindex.hpp \
		$(DIR)/names.hpp \
		$(DIR)/numerics.h \
		$(DIR)/numerics2.hpp \
		$(DIR)/decays/one_loop_decay_diagrams.hpp \
		$(DIR)/observables.hpp \
		$(DIR)/observable_problems.hpp \
		$(DIR)/observable_problems_format.hpp \
		$(DIR)/observable_problems_format_slha.hpp \
		$(DIR)/observable_problems_format_mathlink.hpp \
		$(DIR)/physical_input.hpp \
		$(DIR)/pmns.hpp \
		$(DIR)/pp_map.hpp \
		$(DIR)/problems.hpp \
		$(DIR)/problems_format_mathlink.hpp \
		$(DIR)/raii.hpp \
		$(DIR)/rg_flow.hpp \
		$(DIR)/rk.hpp \
		$(DIR)/rkf_integrator.hpp \
		$(DIR)/root_finder.hpp \
		$(DIR)/scan.hpp \
		$(DIR)/sfermions.hpp \
		$(DIR)/single_scale_constraint.hpp \
		$(DIR)/single_scale_matching.hpp \
		$(DIR)/slha_format.hpp \
		$(DIR)/slha_io.hpp \
		$(DIR)/spectrum_generator_problems.hpp \
		$(DIR)/spectrum_generator_settings.hpp \
		$(DIR)/string_conversion.hpp \
		$(DIR)/string_format.hpp \
		$(DIR)/string_utils.hpp \
		$(DIR)/sum.hpp \
		$(DIR)/thread_pool.hpp \
		$(DIR)/threshold_corrections.hpp \
		$(DIR)/threshold_loop_functions.hpp \
		$(DIR)/trilog.hpp \
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

LIBAUX := ''
LIBAUX_AUTOGEN := ''

# files which allow some useful things with fortran functions ##########
ifeq (yes, $(sort $(filter yes, $(ENABLE_FFLITE) $(ENABLE_COLLIER) $(ENABLE_LOOPTOOLS) )))

LIBFLEXI_SRC += \
	$(DIR)/fortran_utils.cpp

LIBFLEXI_HDR += \
	$(DIR)/fortran_utils.hpp

FUTI := $(DIR)/libfortran_utils

$(DIR)/fortran_utils.hpp $(DIR)/fortran_utils.cpp : $(FUTI).a

$(FUTI).a : $(FUTI).o
	$(Q)$(MSG)
	$(Q)$(MODULE_MAKE_LIB_CMD) $@ $(FUTI).o

$(FUTI).o : $(FUTI).f90
	$(Q)$(MSG)
	$(Q)$(FC) $(FFLAGS) $(FSTD) $(FMOD) src -c $< -o $(FUTI).o

$(FUTI).mod : $(FUTI).f90 $(FUTI).o
	@true

LIBAUX += \
	$(FUTI).a

endif

# loop library #########################################################
LOOP_DIR := $(DIR)/loop_libraries

LOOP_SRC := \
		$(LOOP_DIR)/loop_library.cpp

LOOP_HDR := \
		$(LOOP_DIR)/loop_library.hpp \
		$(LOOP_DIR)/loop_library_interface.hpp

LIBFLEXI_SRC += $(LOOP_SRC)
LIBFLEXI_HDR += $(LOOP_HDR)

ifeq ($(ENABLE_LOOPTOOLS),yes)
LIBFLEXI_SRC += \
		$(LOOP_DIR)/library_looptools.cpp
LIBFLEXI_HDR += \
		$(LOOP_DIR)/library_looptools.hpp
endif

ifeq ($(ENABLE_FFLITE),yes)
LIBFLEXI_SRC += \
		$(LOOP_DIR)/library_fflite.cpp
LIBFLEXI_HDR += \
		$(LOOP_DIR)/library_fflite.hpp
endif

ifeq ($(ENABLE_COLLIER),yes)
LIBFLEXI_SRC += \
		$(LOOP_DIR)/library_collier.cpp
LIBFLEXI_HDR += \
		$(LOOP_DIR)/library_collier.hpp

COLLWRAP := $(LOOP_DIR)/libcollier_wrapper

$(LOOP_HDR) $(LOOP_SRC) : $(COLLWRAP).a

$(COLLWRAP).a : $(COLLWRAP).o $(COLLWRAP).mod
	$(Q)$(MSG)
	$(Q)$(MODULE_MAKE_LIB_CMD) $@ $(COLLWRAP).o

$(COLLWRAP).o : $(COLLWRAP).f90
	$(Q)$(MSG)
	$(Q)$(FC) $(FFLAGS) $(FSTD) $(FMOD) $(LOOP_DIR) -c $< $(COLLIERFLAGS) -o $(COLLWRAP).o

$(COLLWRAP).mod :  $(COLLWRAP).f90 $(COLLWRAP).o
	@true

$(COLLWRAP).f90 : $(COLLWRAP).cpp
	$(Q)$(MSG)
	$(Q)$(CXX) -E $< | tr '@' '\n' > $@

LIBAUX += \
	$(COLLWRAP).a

LIBAUX_AUTOGEN += \
	$(COLLWRAP).f90

endif

LIBFLEXI_OBJ := \
		$(patsubst %.cpp, %.o, $(filter %.cpp, $(LIBFLEXI_SRC))) \
		$(patsubst %.c, %.o, $(filter %.c, $(LIBFLEXI_SRC))) \
		$(patsubst %.f, %.o, $(filter %.f, $(LIBFLEXI_SRC)))

LIBAUX_OBJ := \
		$(patsubst %.a, %.o, $(LIBAUX)) \
		$(patsubst %.a, %.mod, $(LIBAUX))

LIBFLEXI_DEP := \
		$(LIBFLEXI_OBJ:.o=.d)

LIBFLEXI     := $(DIR)/libflexisusy$(MODULE_LIBEXT)

LIBFLEXI_INSTALL_DIR := $(INSTALL_DIR)/$(DIR)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) clean-$(MODNAME)-dep \
		clean-$(MODNAME)-lib clean-$(MODNAME)-obj distclean-$(MODNAME) \
		clean-$(MODNAME)-autogen

all-$(MODNAME): $(LIBFLEXI)
		@true

ifneq ($(INSTALL_DIR),)
install-src::
		$(Q)install -d $(LIBFLEXI_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBFLEXI_SRC) $(LIBFLEXI_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBFLEXI_HDR) $(LIBFLEXI_INSTALL_DIR)
		$(Q)install -m u=rw,g=r,o=r $(LIBFLEXI_MK) $(LIBFLEXI_INSTALL_DIR)
endif

clean-$(MODNAME)-dep:
		$(Q)-rm -f $(LIBFLEXI_DEP)

clean-$(MODNAME)-lib:
		$(Q)-rm -f $(LIBFLEXI) $(LIBAUX)

clean-$(MODNAME)-obj:
		$(Q)-rm -f $(LIBFLEXI_OBJ) $(LIBAUX_OBJ)

clean-$(MODNAME)-autogen:
		$(Q)-rm -f $(LIBAUX_AUTOGEN)

clean-$(MODNAME): clean-$(MODNAME)-dep clean-$(MODNAME)-lib clean-$(MODNAME)-obj clean-$(MODNAME)-autogen
		@true

distclean-$(MODNAME): clean-$(MODNAME)

clean-obj::     clean-$(MODNAME)-obj

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)

$(LIBFLEXI_DEP) $(LIBFLEXI_OBJ): CPPFLAGS += $(GSLFLAGS) $(EIGENFLAGS) $(BOOSTFLAGS) $(SQLITEFLAGS) $(GM2CALCFLAGS) $(TSILFLAGS)

ifneq (,$(findstring yes,$(ENABLE_LOOPTOOLS)$(ENABLE_FFLITE)))
$(LIBFLEXI_DEP) $(LIBFLEXI_OBJ): CPPFLAGS += $(LOOPFUNCFLAGS)
endif

ifeq ($(ENABLE_SHARED_LIBS),yes)
$(LIBFLEXI): $(LIBFLEXI_OBJ)
		@$(MSG)
		$(Q)$(MODULE_MAKE_LIB_CMD) $@ $^ $(GSLLIBS) $(FLIBS) $(SQLITELIBS) $(GM2CALCLIBS) $(THREADLIBS)
else
$(LIBFLEXI): $(LIBFLEXI_OBJ)
		@$(MSG)
		$(Q)$(MODULE_MAKE_LIB_CMD) $@ $^
endif

ALLDEP += $(LIBFLEXI_DEP)
ALLLIB += $(LIBFLEXI)
