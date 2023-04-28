template <>
double CLASSNAME::get_partial_width<PseudoscalarHiggs, Photon, ZBoson>(
      const context_base& context,
      const typename field_indices<PseudoscalarHiggs>::type& in_idx,
      const typename field_indices<Photon>::type& out1_idx,
      const typename field_indices<ZBoson>::type& out2_idx)
{
   if (in_idx.at(0) < info::number_of_neutral_goldstones) {
      throw OutOfBoundsError("Error in " + create_process_string<PseudoscalarHiggs,Photon,ZBoson>(in_idx, out1_idx, out2_idx) + " decay. Decaying particle is a Goldstone.");
   }

   const auto amp = calculate_amplitude<PseudoscalarHiggs, Photon, ZBoson>(context, in_idx, out1_idx, out2_idx);
   const double mAh = context.physical_mass<PseudoscalarHiggs>(in_idx);
   const double mZ = qedqcd.displayPoleMZ();
   const double ps = 1./(8.*Pi) * std::sqrt(KallenLambda(1., 0., Sqr(mZ/mAh)));
   const double flux = 0.5/mAh;
   auto res = flux * ps * amp.square();

   // use alpha_em in the Thomson limit
   if (flexibledecay_settings.get(FlexibleDecay_settings::use_Thomson_alpha_in_Phigamgam_and_PhigamZ)) {
      const double alpha_em_0 = physical_input.get(Physical_input::alpha_em_0);
      const double alpha_em = get_alpha(context);
      res *= alpha_em_0/alpha_em;
   }

#ifdef ENABLE_HIGGSTOOLS
   if (static_cast<bool>(flexibledecay_settings.get(FlexibleDecay_settings::call_higgstools))) {
      higgstools_input.add_coupling(field_as_string<PseudoscalarHiggs>(in_idx), {22, 23}, std::sqrt(res));
   }
#endif

   return res;
}
