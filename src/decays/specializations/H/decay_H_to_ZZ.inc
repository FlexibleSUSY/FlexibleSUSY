template <>
double CLASSNAME::get_partial_width<Higgs,ZBoson,ZBoson>(
   const context_base& context,
   typename field_indices<Higgs>::type const& indexIn,
   typename field_indices<ZBoson>::type const& indexOut1,
   typename field_indices<ZBoson>::type const& indexOut2)
{

   const double mHOS = context.physical_mass<Higgs>(indexIn);
   // There might be large differences between mZ from mass block
   // and one from slha input, especially in the decoupling limit
   // so we use the latter one. There might be a problem with
   // models where Z mixes with something else.
   // const double mZOS = context.physical_mass<Z>(indexOut1);
   const double mZOS = qedqcd.displayPoleMZ();
   const double x = Sqr(mZOS/mHOS);

   const auto indices = concatenate(indexIn, indexOut2, indexOut1);
   const auto ghZZ =
      Vertex<Higgs, ZBoson, ZBoson>::evaluate(indices, context).value();

   // double off-shell decay H -> Z* Z* -> 4 fermions
   const auto four_body = [&] {
      static constexpr double GammaZ = 2.4952;
      const auto integral = hVV_4body(mHOS, mZOS, GammaZ);
      const double normalization = 1./(64.*Cube(Pi)) * std::norm(ghZZ) * Cube(mHOS)/Power4(mZOS)/2.;
      return HVV_4body_result{integral.value * normalization,
                              integral.error * normalization};
   };

   // single off-shell decay H -> Z Z* -> Z f fbar
   bool bsm_final_states = false;
   const auto three_body = [&] {
      if (check_3body_Vff_decay<BSMForZdecay,ZBoson>(context, mHOS, indexOut1)) {
         bsm_final_states = true;
         return 0.;
      }

      const double sw2 = Sqr(std::sin(context.model.ThetaW()));
      const double deltaV = 7.0/12.0 - 10.0/9.0*sw2 + 40.0/27.0*Sqr(sw2);
      const double g2 = context.model.get_g2();

      return 3. / (512. * Power3(Pi)) * 1. / mHOS * deltaV * RT(x) / x *
             std::norm(ghZZ * g2) / (1 - sw2);
   };

   // two-body decay H -> Z Z
   const auto two_body = [&] {
      const double flux = 1. / (2 * mHOS);
      // phase space without symmetry factor
      const double ps = 1. / (8. * Pi) * std::sqrt(KallenLambda(1., x, x));

      // phase space symmetry factor
      const double ps_symmetry = 1./2.;

      // matrix element squared
      const auto mat_elem = calculate_amplitude<Higgs, ZBoson, ZBoson>(
         context, indexIn, indexOut1, indexOut2);
      const auto mat_elem_sq = mat_elem.square();

      // flux * phase space factor * matrix element squared
      return flux * ps * ps_symmetry * mat_elem_sq;
   };

   const int offshell_VV_decays =
      flexibledecay_settings.get(FlexibleDecay_settings::offshell_VV_decays);
   const auto width =
      hVV_width(x, offshell_VV_decays, four_body, three_body, two_body);

   if (bsm_final_states) {
      const std::string index_as_string =
         indexIn.size() == 0 ? "" : "(" + std::to_string(indexIn.at(0)) + ")";
      WARNING("Warning in H" + index_as_string +
              "->ZZ decays: Single off-shell decays H->Zff' assume no possible "
              "BSM particles in the final state. Turning off.");
      return 0.;
   }

   const double res = width.value;
   static constexpr double errorThresholdPerc = 0.1;
   if (width.error * 100 > errorThresholdPerc * res) {
      problems.add_warning(create_process_string<Higgs, ZBoson, ZBoson>(
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
         {23, 23},
         std::pair<std::string, double> {field_as_string<Higgs>(indexIn) + "-" + field_as_string<ZBoson>(indexOut1) + "-" + field_as_string<ZBoson>(indexOut2), std::sqrt(res)},
         res
      );
   }

   return res;
}
