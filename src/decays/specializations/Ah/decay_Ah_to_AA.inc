template <>
double CLASSNAME::get_partial_width<PseudoscalarHiggs, Photon, Photon>(
      const context_base& context,
      const typename field_indices<PseudoscalarHiggs>::type& in_idx,
      const typename field_indices<Photon>::type& out1_idx,
      const typename field_indices<Photon>::type& out2_idx) const
{
   if (in_idx.at(0) < info::number_of_neutral_goldstones) {
      throw OutOfBoundsError("Error in " + create_process_string<PseudoscalarHiggs,Photon,Photon>(in_idx, out1_idx, out2_idx) + " decay. Decaying particle is a Goldstone.");
   }

   const auto amp = calculate_amplitude<PseudoscalarHiggs, Photon, Photon>(context, in_idx, out1_idx, out2_idx);
   const double mAH = context.physical_mass<PseudoscalarHiggs>(in_idx);
   static constexpr double ps {1./(8.*Pi)};
   static constexpr double ps_symmetry {0.5};
   const double flux = 0.5/mAH;
   auto res = flux * ps * ps_symmetry * amp.square();

   // use alpha_em in the Thomson limit
   if (flexibledecay_settings.get(FlexibleDecay_settings::use_Thomson_alpha_in_Phigamgam_and_PhigamZ)) {
      const double alpha_em_0 = physical_input.get(Physical_input::alpha_em_0);
      const double alpha_em = get_alpha(context);
      res *= Sqr(alpha_em_0/alpha_em);
   }

   return res;
}
