template <>
double CLASSNAME::get_partial_width<AH, A, Z>(
      const context_base& context,
      const typename field_indices<AH>::type& in_idx,
      const typename field_indices<A>::type& out1_idx,
      const typename field_indices<Z>::type& out2_idx)
{
   if (in_idx.at(0) < info::number_of_neutral_goldstones) {
      throw OutOfBoundsError("Error in " + create_process_string<AH,A,Z>(in_idx, out1_idx, out2_idx) + " decay. Decaying particle is a Goldstone.");
   }

   const auto amp = calculate_amplitude<AH, A, Z>(context, in_idx, out1_idx, out2_idx);
   const double mAH = context.physical_mass<AH>(in_idx);
   const double mZ = qedqcd.displayPoleMZ();
   const double ps = 1./(8.*Pi) * std::sqrt(KallenLambda(1., 0., Sqr(mZ/mAH)));
   const double flux = 0.5/mAH;
   auto res = flux * ps * amp.square();

   // use alpha_em in the Thomson limit
   if (flexibledecay_settings.get(FlexibleDecay_settings::use_Thomson_alpha_in_Phigamgam_and_PhigamZ)) {
      const double alpha_em_0 = physical_input.get(Physical_input::alpha_em_0);
      const double alpha_em = get_alpha(context);
      res *= alpha_em_0/alpha_em;
   }

#ifdef ENABLE_HIGGSTOOLS
   if (static_cast<bool>(flexibledecay_settings.get(FlexibleDecay_settings::call_higgstools))) {
      higgstools_input.push_back(std::make_tuple(field_as_string<AH>(in_idx), 22, 23, mAH, std::sqrt(res)));
   }
#endif

   return res;
}
