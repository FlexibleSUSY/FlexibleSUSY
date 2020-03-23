DIR          := templates
MODNAME      := templates

BASE_TEMPLATES := \
		$(DIR)/cxx_qft/qft.hpp.in \
		$(DIR)/cxx_qft/fields.hpp.in \
		$(DIR)/cxx_qft/vertices.hpp.in \
		$(DIR)/cxx_qft/vertices_.cpp.in \
		$(DIR)/cxx_qft/context_base.hpp.in \
		$(DIR)/cxx_qft/npointfunctions_wilsoncoeffs.hpp.in \
		$(DIR)/a_muon.hpp.in \
		$(DIR)/a_muon.cpp.in \
		$(DIR)/edm.hpp.in \
		$(DIR)/edm.cpp.in \
		$(DIR)/FFV_form_factors.hpp.in \
		$(DIR)/FFV_form_factors.cpp.in \
		$(DIR)/l_to_lgamma.hpp.in \
		$(DIR)/l_to_lgamma.cpp.in \
		$(DIR)/f_to_f_conversion.hpp.in \
		$(DIR)/f_to_f_conversion.cpp.in \
		$(DIR)/b_to_s_gamma.cpp.in \
		$(DIR)/b_to_s_gamma.hpp.in \
		$(DIR)/convergence_tester.hpp.in \
		$(DIR)/ewsb_solver.hpp.in \
		$(DIR)/ewsb_solver_interface.hpp.in \
		$(DIR)/effective_couplings.hpp.in \
		$(DIR)/effective_couplings.cpp.in \
		$(DIR)/high_scale_constraint.hpp.in \
		$(DIR)/info.hpp.in \
		$(DIR)/info.cpp.in \
		$(DIR)/initial_guesser.hpp.in \
		$(DIR)/input_parameters.hpp.in \
		$(DIR)/input_parameters.cpp.in \
		$(DIR)/librarylink.cpp.in \
		$(DIR)/librarylink.m.in \
		$(DIR)/low_scale_constraint.hpp.in \
		$(DIR)/mass_eigenstates.hpp.in \
		$(DIR)/mass_eigenstates.cpp.in \
		$(DIR)/mass_eigenstates_interface.hpp.in \
		$(DIR)/mass_eigenstates_decoupling_scheme.hpp.in \
		$(DIR)/mass_eigenstates_decoupling_scheme.cpp.in \
		$(DIR)/model.hpp.in \
		$(DIR)/model_slha.hpp.in \
		$(DIR)/observables.hpp.in \
		$(DIR)/observables.cpp.in \
		$(DIR)/physical.hpp.in \
		$(DIR)/physical.cpp.in \
		$(DIR)/plot_rgflow.gnuplot.in \
		$(DIR)/plot_spectrum.gnuplot.in \
		$(DIR)/run.cpp.in \
		$(DIR)/run.m.in \
		$(DIR)/run_cmd_line.cpp.in \
		$(DIR)/scan.cpp.in \
		$(DIR)/slha_io.hpp.in \
		$(DIR)/slha_io.cpp.in \
		$(DIR)/soft_beta_.cpp.in \
		$(DIR)/soft_parameters.hpp.in \
		$(DIR)/soft_parameters.cpp.in \
		$(DIR)/spectrum_generator.hpp.in \
		$(DIR)/spectrum_generator_interface.hpp.in \
		$(DIR)/standard_model_spectrum_generator_interface.hpp.in \
		$(DIR)/standard_model_matching.hpp.in \
		$(DIR)/standard_model_matching.cpp.in \
		$(DIR)/susy_beta_.cpp.in \
		$(DIR)/susy_parameters.hpp.in \
		$(DIR)/susy_parameters.cpp.in \
		$(DIR)/susy_scale_constraint.hpp.in \
		$(DIR)/utilities.hpp.in \
		$(DIR)/utilities.cpp.in

TWO_SCALE_TEMPLATES := \
		$(DIR)/standard_model_two_scale_high_scale_initial_guesser.cpp.in \
		$(DIR)/standard_model_two_scale_high_scale_initial_guesser.hpp.in \
		$(DIR)/standard_model_two_scale_high_scale_spectrum_generator.hpp.in \
		$(DIR)/standard_model_two_scale_high_scale_spectrum_generator.cpp.in \
		$(DIR)/standard_model_two_scale_low_scale_initial_guesser.cpp.in \
		$(DIR)/standard_model_two_scale_low_scale_initial_guesser.hpp.in \
		$(DIR)/standard_model_two_scale_low_scale_spectrum_generator.hpp.in \
		$(DIR)/standard_model_two_scale_low_scale_spectrum_generator.cpp.in \
		$(DIR)/standard_model_two_scale_matching.hpp.in \
		$(DIR)/standard_model_two_scale_matching.cpp.in \
		$(DIR)/two_scale_convergence_tester.hpp.in \
		$(DIR)/two_scale_convergence_tester.cpp.in \
		$(DIR)/two_scale_ewsb_solver.hpp.in \
		$(DIR)/two_scale_ewsb_solver.cpp.in \
		$(DIR)/two_scale_high_scale_constraint.hpp.in \
		$(DIR)/two_scale_high_scale_constraint.cpp.in \
		$(DIR)/two_scale_high_scale_initial_guesser.hpp.in \
		$(DIR)/two_scale_high_scale_initial_guesser.cpp.in \
		$(DIR)/two_scale_high_scale_spectrum_generator.hpp.in \
		$(DIR)/two_scale_high_scale_spectrum_generator.cpp.in \
		$(DIR)/two_scale_low_scale_constraint.hpp.in \
		$(DIR)/two_scale_low_scale_constraint.cpp.in \
		$(DIR)/two_scale_low_scale_initial_guesser.hpp.in \
		$(DIR)/two_scale_low_scale_initial_guesser.cpp.in \
		$(DIR)/two_scale_low_scale_spectrum_generator.hpp.in \
		$(DIR)/two_scale_low_scale_spectrum_generator.cpp.in \
		$(DIR)/two_scale_model.hpp.in \
		$(DIR)/two_scale_model.cpp.in \
		$(DIR)/two_scale_susy_scale_constraint.hpp.in \
		$(DIR)/two_scale_susy_scale_constraint.cpp.in \
		$(DIR)/weinberg_angle.hpp.in \
		$(DIR)/weinberg_angle.cpp.in

SEMI_ANALYTIC_TEMPLATES := \
		$(DIR)/semi_analytic_convergence_tester.hpp.in \
		$(DIR)/semi_analytic_convergence_tester.cpp.in \
		$(DIR)/semi_analytic_ewsb_solver.hpp.in \
		$(DIR)/semi_analytic_ewsb_solver.cpp.in \
		$(DIR)/semi_analytic_high_scale_constraint.hpp.in \
		$(DIR)/semi_analytic_high_scale_constraint.cpp.in \
		$(DIR)/semi_analytic_high_scale_initial_guesser.hpp.in \
		$(DIR)/semi_analytic_high_scale_initial_guesser.cpp.in \
		$(DIR)/semi_analytic_high_scale_spectrum_generator.hpp.in \
		$(DIR)/semi_analytic_high_scale_spectrum_generator.cpp.in \
		$(DIR)/semi_analytic_low_scale_constraint.hpp.in \
		$(DIR)/semi_analytic_low_scale_constraint.cpp.in \
		$(DIR)/semi_analytic_low_scale_initial_guesser.hpp.in \
		$(DIR)/semi_analytic_low_scale_initial_guesser.cpp.in \
		$(DIR)/semi_analytic_low_scale_spectrum_generator.hpp.in \
		$(DIR)/semi_analytic_low_scale_spectrum_generator.cpp.in \
		$(DIR)/semi_analytic_model.hpp.in \
		$(DIR)/semi_analytic_model.cpp.in \
		$(DIR)/semi_analytic_soft_parameters_constraint.hpp.in \
		$(DIR)/semi_analytic_soft_parameters_constraint.cpp.in \
		$(DIR)/semi_analytic_solutions.hpp.in \
		$(DIR)/semi_analytic_solutions.cpp.in \
		$(DIR)/semi_analytic_susy_convergence_tester.hpp.in \
		$(DIR)/semi_analytic_susy_convergence_tester.cpp.in \
		$(DIR)/semi_analytic_susy_scale_constraint.hpp.in \
		$(DIR)/semi_analytic_susy_scale_constraint.cpp.in \
		$(DIR)/soft_parameters_constraint.hpp.in \
		$(DIR)/susy_convergence_tester.hpp.in

TEMPLATES    := \
		$(BASE_TEMPLATES) \
		$(TWO_SCALE_TEMPLATES) \
		$(SEMI_ANALYTIC_TEMPLATES)

.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

all-$(MODNAME):
		@true

clean-$(MODNAME):
		@true

distclean-$(MODNAME): clean-$(MODNAME)

clean::         clean-$(MODNAME)

distclean::     distclean-$(MODNAME)
