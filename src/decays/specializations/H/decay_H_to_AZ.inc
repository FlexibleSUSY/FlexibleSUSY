template <>
double CLASSNAME::get_partial_width<Higgs, Photon, ZBoson>(
      const context_base& context,
      const typename field_indices<Higgs>::type& in_idx,
      const typename field_indices<Photon>::type& out1_idx,
      const typename field_indices<ZBoson>::type& out2_idx)
{
   const auto amp = calculate_amplitude_1l<Higgs, Photon, ZBoson>(context, in_idx, out1_idx, out2_idx);
   const double mH = context.physical_mass<Higgs>(in_idx);
   const double mZ = qedqcd.displayPoleMZ();
   const double ps = 1./(8.*Pi) * std::sqrt(KallenLambda(1., 0., Sqr(mZ/mH)));
   const double flux = 0.5/mH;
   auto res = flux * ps * amp.square();

   // use alpha_em in the Thomson limit
   if (flexibledecay_settings.get(FlexibleDecay_settings::use_Thomson_alpha_in_Phigamgam_and_PhigamZ)) {
      const double alpha_em_0 = physical_input.get(Physical_input::alpha_em_0);
      const double alpha_em = get_alpha(context);
      res *= alpha_em_0/alpha_em;
   }

   if (flexibledecay_settings.get(FlexibleDecay_settings::print_effc_block)) {
      effhiggscouplings_block_input.push_back(
         {
            25 + 10*(in_idx.size()>0 ? in_idx.at(0) : 0), 22, 23,
            std::sqrt(res/(flux * ps)/(0.5*Sqr(Sqr(mH)-Sqr(mZ)))),
            field_as_string<Higgs>(in_idx) + "-" + field_as_string<Photon>(out1_idx) + "-" + field_as_string<ZBoson>(out2_idx)
         }
      );
   }

   return res;
}
