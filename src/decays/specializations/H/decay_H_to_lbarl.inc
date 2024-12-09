template <>
double CLASSNAME::get_partial_width<Higgs, bar<ChargedLepton>::type, ChargedLepton>(
   const context_base& context,
   typename field_indices<Higgs>::type const& indexIn,
   typename field_indices<bar<ChargedLepton>::type>::type const& indexOut1,
   typename field_indices<ChargedLepton>::type const& indexOut2)
{

   const double mHOS = context.physical_mass<Higgs>(indexIn);
   const double mL1OS = context.physical_mass<bar<ChargedLepton>::type>(indexOut1);
   const double mL2OS = context.physical_mass<ChargedLepton>(indexOut2);

   // phase space without symmetry factor
   const auto xOS1 = Sqr(mL1OS/mHOS);
   const auto xOS2 = Sqr(mL2OS/mHOS);
   const auto betaOS2 = 1.0 - Sqr((mL1OS + mL2OS)/mHOS);
   const double ps = 1./(8.*Pi)*std::sqrt(KallenLambda(1., xOS1, xOS2));

   const auto amp_tree = calculate_amplitude_tree<Higgs, bar<ChargedLepton>::type, ChargedLepton>(context, indexIn, indexOut1, indexOut2);
   const auto amp_1l   = calculate_amplitude_1l<Higgs, bar<ChargedLepton>::type, ChargedLepton>(context, indexIn, indexOut1, indexOut2);
   standard_model_cxx_diagrams::context_base sm_context {sm};
   const auto se_e = self_energy_1loop_PR<Fe>(context, mL1OS).coeff(indexOut1.at(0), indexOut1.at(0)) - standard_model_cxx_diagrams::self_energy_1loop_PR<standard_model_cxx_diagrams::fields::Fe>(sm_context, mL1OS).coeff(indexOut1.at(0), indexOut1.at(0));
   const auto dse_e = self_energy_1loop_PR_deriv_p2<Fe>(context, mL1OS).coeff(indexOut1.at(0), indexOut1.at(0)) - standard_model_cxx_diagrams::self_energy_1loop_PR_deriv_p2<standard_model_cxx_diagrams::fields::Fe>(sm_context, mL1OS).coeff(indexOut1.at(0), indexOut1.at(0));

   // flux * phase space factor * symmetry factor * |matrix element|^2
   double res_tree2 = 0.5 * ps * amp_tree.square()/mHOS;

   double res = res_tree2;
   // higher order corrections
   if (flexibledecay_settings.get(FlexibleDecay_settings::include_higher_order_corrections)) {
      // 1-loop QED corrections
      res *= 1. + get_alpha(context)/Pi*17./4.;
      //res += 2.*amplitude_interference(amp_tree, amp_1l);
   }

   return res;
}
