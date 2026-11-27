template <>
double CLASSNAME::get_partial_width<Higgs, WpBoson, WmBoson>(
   const context_base& context, typename field_indices<Higgs>::type const& indexIn,
   typename field_indices<WpBoson>::type const& indexOut1,
   typename field_indices<WmBoson>::type const& indexOut2)
{

   const double mHOS = context.physical_mass<Higgs>(indexIn);
   const double mWOS = context.physical_mass<WpBoson>(indexOut1);
   const double x = Sqr(mWOS/mHOS);

   const auto indices = concatenate(indexIn, indexOut1, indexOut2);
   const auto ghWW =
      Vertex<Higgs, WpBoson, WmBoson>::evaluate(indices, context).value();

   // double off-shell decay H -> W* W* -> 4 fermions
   const auto four_body = [&] {
      static constexpr double GammaW = 2.085; //3*std::norm(ghWW)/(16.*Pi*mWOS);
      const auto integral = hVV_4body(mHOS, mWOS, GammaW);
      const double normalization = 1./(64.*Cube(Pi)) * std::norm(ghWW) * Cube(mHOS)/Power4(mWOS);
      return HVV_4body_result{integral.value * normalization,
                              integral.error * normalization};
   };

   // single off-shell decay H -> W W* -> W f f'
   bool bsm_final_states = false;
   const auto three_body = [&] {
      if (check_3body_Vff_decay<BSMForWdecay, WpBoson>(context, mHOS, indexOut1)) {
         bsm_final_states = true;
         return 0.;
      }

      // absolute value of baru d W+ vertex (no CKM and no PL projector)
      const double g2 = context.model.get_g2();
      const double gWud = g2/std::sqrt(2.0);

      // number of final states
      static constexpr double NLF = 3;  // number of lepton flavours
      static constexpr double Nc = 3;   // number of colors
      static constexpr double NQF = 2;  // number of quark flavours

      return 1. / (768. * Power3(Pi) * mHOS) * RT(x) / x *
             std::norm(ghWW * gWud) * (NLF + Nc * NQF);
   };

   // two-body decay H -> W W
   const auto two_body = [&] {
      const double flux = 1. / (2 * mHOS);
      // phase space without symmetry factor
      const double ps = 1./(8.*Pi)*std::sqrt(KallenLambda(1., x, x));

      // matrix element squared
      const auto mat_elem = calculate_amplitude<Higgs, WpBoson, WmBoson>(
         context, indexIn, indexOut1, indexOut2);
      const auto mat_elem_sq = mat_elem.square();

      // flux * phase space factor * matrix element squared
      return flux * ps * mat_elem_sq;
   };

   const int offshell_VV_decays =
      flexibledecay_settings.get(FlexibleDecay_settings::offshell_VV_decays);
   const auto width =
      hVV_width(x, offshell_VV_decays, four_body, three_body, two_body);

   if (bsm_final_states) {
      const std::string index_as_string =
         indexIn.size() == 0 ? "" : "(" + std::to_string(indexIn.at(0)) + ")";
      WARNING("Warning in H" + index_as_string +
              "->WW decays: Single off-shell decays H->Wff' assume no possible "
              "BSM particles in the final state. Turning off.");
      return 0.;
   }

   const double res = width.value;
   static constexpr double errorThresholdPerc = 0.1;
   if (width.error * 100 > errorThresholdPerc * res) {
      problems.add_warning(create_process_string<Higgs, WpBoson, WmBoson>(
                              indexIn, indexOut1, indexOut2) +
                           ": Relative integration error > " +
                           std::to_string(errorThresholdPerc) + "%");
   }

   if (flexibledecay_settings.get(FlexibleDecay_settings::call_higgstools) != 0 ||
       flexibledecay_settings.get(FlexibleDecay_settings::call_lilith) != 0 ||
       flexibledecay_settings.get(FlexibleDecay_settings::calculate_normalized_effc) != 0
   ) {
      neutral_higgs_effc.add_coupling(
         field_as_string<Higgs>(indexIn),
         {-24, 24},
         std::pair<std::string, double> {field_as_string<Higgs>(indexIn) + "-" + field_as_string<WpBoson>(indexOut1) + "-" + field_as_string<WmBoson>(indexOut2), std::sqrt(res)},
         res
      );
   }

   return res;
}
