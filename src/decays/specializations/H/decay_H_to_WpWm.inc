template <>
double CLASSNAME::get_partial_width<Higgs, WpBoson, WmBoson>(
   const context_base& context, typename field_indices<Higgs>::type const& indexIn,
   typename field_indices<WpBoson>::type const& indexOut1,
   typename field_indices<WmBoson>::type const& indexOut2)
{

   const double mHOS = context.physical_mass<Higgs>(indexIn);
   const double mWOS = context.physical_mass<WpBoson>(indexOut1);
   const double x = Sqr(mWOS/mHOS);
   double res = 0;

   const int offshell_VV_decays = flexibledecay_settings.get(FlexibleDecay_settings::offshell_VV_decays);
   if ((x > 1.0 && offshell_VV_decays != 0) ||
       (4.*x > 1.0 && offshell_VV_decays == 2)) {

      // integrand
      static constexpr double GammaW = 2.085; //3*std::norm(ghWW)/(16.*Pi*mWOS);
      struct hVV_4body_params params = {mHOS, mWOS, GammaW};
      gsl_monte_function G = { &hVV_4body, 2, &params };

      // setup integration
      gsl_rng_env_setup ();
      double xl[2] = {0, 0};
      double xu[2] = {Sqr(mHOS), Sqr(mHOS)};
      // this gives relative error < 0.05%
      static constexpr size_t calls = 1'000'000;
      double err;
      const gsl_rng_type *T = gsl_rng_default;
      gsl_rng *r = gsl_rng_alloc (T);
      gsl_monte_miser_state *s = gsl_monte_miser_alloc (2);
      gsl_monte_miser_integrate (&G, xl, xu, 2, calls, r, s,
                                 &res, &err);

      // clean-up
      gsl_monte_miser_free (s);
      gsl_rng_free (r);

      // prefactor
      const auto indices = concatenate(indexIn, indexOut1, indexOut2);
      const auto ghWW =
         Vertex<Higgs, WpBoson, WmBoson>::evaluate(indices, context).value();

      const double normalization = 1./(64.*Cube(Pi)) * std::norm(ghWW) * Cube(mHOS)/Power4(mWOS);
      res *= normalization;
      err *= normalization;
      const double rel_err = err/res;
      static constexpr double errorThresholdPerc = 0.1;
      if (rel_err*100 > errorThresholdPerc) {
         problems.add_warning(
            create_process_string<Higgs,WpBoson,WmBoson>(indexIn, indexOut1, indexOut2)
               + ": Relative integration error > " + std::to_string(errorThresholdPerc) + "%"
         );
      }
   }
   // three-body decays form mW < mH < 2*mW
   else if (4 * x > 1.0 && offshell_VV_decays != 0) {

      if (check_3body_Vff_decay<BSMForWdecay, WpBoson>(context, mHOS, indexOut1)) {
         const std::string index_as_string = indexIn.size() == 0 ? "" : "(" + std::to_string(indexIn.at(0)) + ")";
         WARNING("Warning in H" + index_as_string + "->WW decays: Single off-shell decays H->Wff' assume no possible BSM particles in the final state. Turning off.");
         return 0.;
      }
      res = 1./(768.*Power3(Pi)*mHOS) * RT(x)/x;

      const auto indices = concatenate(indexIn, indexOut1, indexOut2);
      const auto ghWW =
         Vertex<Higgs, WpBoson, WmBoson>::evaluate(indices, context).value();

      // absolute value of baru d W+ vertex (no CKM and no PL projector)
      const double g2 = context.model.get_g2();
      const double gWud = g2/std::sqrt(2.0);

      res *= std::norm(ghWW*gWud);

      // multiply by number of final states
      static constexpr double NLF = 3;  // number of lepton flavours
      static constexpr double Nc = 3;   // number of colors
      static constexpr double NQF = 2;  // number of quark flavours
      res *= NLF + Nc*NQF;

   // two-body decay for mH > 2 mW
   }
   else if (4.*x < 1.0) {
      const double flux = 1. / (2 * mHOS);
      // phase space without symmetry factor
      const double ps = 1./(8.*Pi)*std::sqrt(KallenLambda(1., x, x));

      // matrix element squared
      const auto mat_elem = calculate_amplitude<Higgs, WpBoson, WmBoson>(
         context, indexIn, indexOut1, indexOut2);
      const auto mat_elem_sq = mat_elem.square();

      // flux * phase space factor * matrix element squared
      res = flux * ps * mat_elem_sq;
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
