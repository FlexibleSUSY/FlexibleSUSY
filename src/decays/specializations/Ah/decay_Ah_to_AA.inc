template <>
double CLASSNAME::get_partial_width<Ah, A, A>(
      const context_base& context,
      const typename field_indices<Ah>::type& in_idx,
      const typename field_indices<A>::type& out1_idx,
      const typename field_indices<A>::type& out2_idx)
{
   if (in_idx.at(0) < info::number_of_neutral_goldstones) {
      throw OutOfBoundsError("Error in " + create_process_string<Ah,A,A>(in_idx, out1_idx, out2_idx) + " decay. Decaying particle is a Goldstone.");
   }

   const auto amp = calculate_amplitude<Ah, A, A>(context, in_idx, out1_idx, out2_idx);
   const double mAH = context.physical_mass<Ah>(in_idx);
   static constexpr double ps {1./(8.*Pi)};
   static constexpr double ps_symmetry {1./2.};
   const double flux = 0.5/mAH;
   auto res = flux * ps * ps_symmetry * amp.square();

   // use alpha_em in the Thomson limit
   if (flexibledecay_settings.get(FlexibleDecay_settings::use_Thomson_alpha_in_Phigamgam_and_PhigamZ)) {
      const double alpha_em_0 = physical_input.get(Physical_input::alpha_em_0);
      const double alpha_em = get_alpha(context);
      res *= Sqr(alpha_em_0/alpha_em);
   }

#ifdef ENABLE_HIGGSTOOLS
   if (static_cast<bool>(flexibledecay_settings.get(FlexibleDecay_settings::call_higgstools))) {
      higgstools_input.add_coupling(field_as_string<Ah>(in_idx), {22, 22}, std::sqrt(res));
   }
#endif

   return res;
}
