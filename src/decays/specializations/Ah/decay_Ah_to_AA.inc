template <>
double CLASSNAME::get_partial_width<PseudoscalarHiggs, Photon, Photon>(
      const context_base& context,
      const typename field_indices<PseudoscalarHiggs>::type& in_idx,
      const typename field_indices<Photon>::type& out1_idx,
      const typename field_indices<Photon>::type& out2_idx)
{
   if (in_idx.at(0) < info::number_of_neutral_goldstones) {
      throw OutOfBoundsError("Error in " + create_process_string<PseudoscalarHiggs,Photon,Photon>(in_idx, out1_idx, out2_idx) + " decay. Decaying particle is a Goldstone.");
   }

   const auto amp = calculate_amplitude<PseudoscalarHiggs, Photon, Photon>(context, in_idx, out1_idx, out2_idx);
   const double mAh = context.physical_mass<PseudoscalarHiggs>(in_idx);
   static constexpr double ps {1./(8.*Pi)};
   static constexpr double ps_symmetry {0.5};
   const double flux = 0.5/mAh;
   auto res = flux * ps * ps_symmetry * amp.square();

   // use alpha_em in the Thomson limit
   if (flexibledecay_settings.get(FlexibleDecay_settings::use_Thomson_alpha_in_Phigamgam_and_PhigamZ)) {
      const double alpha_em_0 = physical_input.get(Physical_input::alpha_em_0);
      const double alpha_em = get_alpha(context);
      res *= Sqr(alpha_em_0/alpha_em);
   }

   if (flexibledecay_settings.get(FlexibleDecay_settings::print_effc_block)) {
      effhiggscouplings_block_input.push_back(
         {
            36 + 10*(in_idx.size()>0 ? in_idx.at(0)-1 : 0), 22, 22,
            std::sqrt(res/(flux * ps * ps_symmetry)/(0.5*Power4(mAh))),
            field_as_string<PseudoscalarHiggs>(in_idx) + "-" + field_as_string<Photon>(out1_idx) + "-" + field_as_string<Photon>(out2_idx)
         }
      );
   }

   if (static_cast<bool>(flexibledecay_settings.get(FlexibleDecay_settings::call_higgstools)) ||
       static_cast<bool>(flexibledecay_settings.get(FlexibleDecay_settings::call_lilith)) ||
       static_cast<bool>(flexibledecay_settings.get(FlexibleDecay_settings::calculate_normalized_effc))
   ) {
      neutral_higgs_effc.add_coupling(
         field_as_string<PseudoscalarHiggs>(in_idx),
         {22, 22},
         std::pair<std::string, double> {field_as_string<PseudoscalarHiggs>(in_idx) + "-" + field_as_string<Photon>(out1_idx) + "-" + field_as_string<Photon>(out2_idx), std::sqrt(res)}
      );
   }

   return res;
}
