// special case for H -> W+ W-
// TODO: implement higher order corrections

template <>
double CLASSNAME::get_partial_width<H, conj<W>::type, W>(
   const context_base& context, typename field_indices<H>::type const& indexIn,
   typename field_indices<conj<W>::type>::type const& indexOut1,
   typename field_indices<W>::type const& indexOut2) const
{

   const double mHOS = context.physical_mass<H>(indexIn);
   const double mWOS = context.physical_mass<W>(indexOut1);
   const double x = Sqr(mWOS/mHOS);
   double res = 0;

   if ((x > 1.0 && flexibledecay_settings.get(FlexibleDecay_settings::offshell_VV_decays) != 0) ||
       (4.*x > 1.0 && flexibledecay_settings.get(FlexibleDecay_settings::offshell_VV_decays) == 2)) {

      // integrand
      constexpr double GammaW = 2.085; //3*std::norm(ghWW)/(16.*Pi*mWOS);
      struct my_f_params params = {mHOS, mWOS, GammaW};
      gsl_monte_function G = { &hVV_4body, 2, &params };

      // setup integration
      gsl_rng_env_setup ();
      double xl[2] = {0, 0};
      double xu[2] = {Sqr(mHOS), Sqr(mHOS)};
      // this gives relative error < 0.05%
      constexpr size_t calls = 1'000'000;
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
      const auto indices = concatenate(indexIn, indexOut2, indexOut1);
      const auto ghWW =
         Vertex<H, conj<W>::type, W>::evaluate(indices, context).value();

      res *= 1./(16.*Cube(Pi))* std::norm(ghWW) * 0.25*Cube(mHOS)/Power4(mWOS);
      err *= 1./(16.*Cube(Pi))* std::norm(ghWW) * 0.25*Cube(mHOS)/Power4(mWOS);
      const double rel_err = err/res;
      if (rel_err > 1e-3) {
         std::ostringstream oss;
         oss <<
            "Warning: result of H->WW is "
            << res << " +/- " << err << " GeV (" << 1e+2*rel_err
            << "%). Relative error > 0.1%.\n";
         WARNING (oss.str());
      }
   }
   // three-body decays form mW < mH < 2*mW
   else if (4 * x > 1.0 && flexibledecay_settings.get(FlexibleDecay_settings::offshell_VV_decays) != 0) {

      if (check_3body_Vff_decay<BSMForWdecay, W>(context, mHOS, indexOut1)) {
         const std::string index_as_string = indexIn.size() == 0 ? "" : "(" + std::to_string(indexIn.at(0)) + ")";
         WARNING("Warning in H" + index_as_string + "->WW decays: Single off-shell decays H->Wff' assume no possible BSM particles in the final state. Turning off.");
         return 0.;
      }
      res = 1./(768.*Power3(Pi)*mHOS) * RT(x)/x;

      const auto indices = concatenate(indexIn, indexOut2, indexOut1);
      const auto ghWW =
         Vertex<H, conj<W>::type, W>::evaluate(indices, context).value();

      // absolute value of baru d W+ vertex (no CKM and no PL projector)
      const double g2 = context.model.get_g2();
      // M_SQRT1_2 =	1/sqrt(2)
      const double gWud = g2*M_SQRT1_2;

      res *= std::norm(ghWW*gWud);

      // multiply by number of final states
      constexpr double NLF = 3;  // number of lepton flavours
      constexpr double Nc = 3;   // number of colors
      constexpr double NQF = 2;  // number of quark flavours
      res *= NLF + Nc*NQF;

   // two-body decay for mH > 2 mW
   }
   else if (4.*x < 1.0) {
      const double flux = 1. / (2 * mHOS);
      // phase space without symmetry factor
      const double ps = 1./(8.*Pi)*std::sqrt(KallenLambda(1., x, x));

      // matrix element squared
      const auto mat_elem = calculate_amplitude<H, conj<W>::type, W>(
         context, indexIn, indexOut1, indexOut2);
      const auto mat_elem_sq =mat_elem.square();

      // flux * phase space factor * matrix element squared
      res = flux * ps * mat_elem_sq;
   }

   return res;
}
