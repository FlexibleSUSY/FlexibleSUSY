template <>
double CLASSNAME::get_partial_width<PseudoscalarHiggs, Gluon, Gluon>(
      const context_base& context,
      const typename field_indices<PseudoscalarHiggs>::type& in_idx,
      const typename field_indices<Gluon>::type& out1_idx,
      const typename field_indices<Gluon>::type& out2_idx)
{
   if (in_idx.at(0) < info::number_of_neutral_goldstones) {
      throw OutOfBoundsError("Error in " + create_process_string<PseudoscalarHiggs,Gluon,Gluon>(in_idx, out1_idx, out2_idx) + " decay. Decaying particle is a Goldstone.");
   }

   const auto amp = calculate_amplitude<PseudoscalarHiggs, Gluon, Gluon>(context, in_idx, out1_idx, out2_idx);
   const double mAh = context.physical_mass<PseudoscalarHiggs>(in_idx);
   static constexpr double ps {1/(8.*Pi)};
   static constexpr double ps_symmetry {0.5};
   static constexpr double color_fact = squared_color_generator<PseudoscalarHiggs, Gluon, Gluon>();
   const double flux = 0.5/mAh;

   double result = flux * color_fact * ps * ps_symmetry * amp.square();

   const int Nf = number_of_active_flavours(qedqcd, mAh);
   if (Nf < 5) {
      problems.add_warning(
         create_process_string<PseudoscalarHiggs,Gluon,Gluon>(in_idx, out1_idx, out2_idx)
         + ": Cannot determine the number of active quark flavours. Disabling higher-order corrections."
      );
   }

   // higher order QCD corrections
   const double tau = Sqr(mAh/(2.*context.mass<UpTypeQuark>({2})));
   if (static_cast<int>(flexibledecay_settings.get(FlexibleDecay_settings::include_higher_order_corrections)) && Nf >= 5 && tau < 0.7) {
      // number of active light flavours
      static constexpr int Nf = 5;
      auto qedqcd_ = qedqcd;
      qedqcd_.to(mAh);
      // 5-flavour SM alpha_s
      const double alpha_s_5f = qedqcd_.displayAlpha(softsusy::ALPHAS);

      const auto indices = concatenate(std::array<int, 1> {2}, std::array<int, 1> {2}, in_idx);
      const auto AhGGVertex = Vertex<bar<UpTypeQuark>::type, UpTypeQuark, PseudoscalarHiggs>::evaluate(indices, context);
      std::complex<double> const AhGGVertexVal = 0.5*(-AhGGVertex.left() + AhGGVertex.right());

      const double tau = Sqr(mAh/(2.*context.mass<UpTypeQuark>({2})));

      const std::complex<double> A12_A = 2.*f(tau)/tau;
      // LO width comming only from the top-loop
      // agrees up to a full double precision with autmatically generated one
      const double Gamma_SM_LO_P = mAh/(18.*Power3(Pi))*std::norm(alpha_s_5f * AhGGVertexVal*sqrt(tau) * 3./4*A12_A);

      const double mu = mAh;
      const double LH = std::log(Sqr(mu/mAh));
      const double deltaNLO {
         97/4. - 7/6.*Nf + (33-2*Nf)/6.*LH
      };

      const double mt {context.mass<UpTypeQuark>({2})};
      const double Lt = std::log(Sqr(mu/mt));
      // eq. D10 of 2207.01032 and 23 of 9807241
      const double deltaNNLO {
         51959/96. - 363/8.*zeta2 - 495/8.*zeta3 + Nf*(-473/8. + 11/2.*zeta2 + 5/4.*zeta3 + Lt) + Sqr(Nf)*(251/216. - 1/6.*zeta2)
         + (3405/16. - 73/3*Nf + 7/12.*Sqr(Nf))*LH
         + (363/16. - 11/4.*Nf + 1/12.*Sqr(Nf))*Sqr(LH)
      };

      const double alpha_s_red = alpha_s_5f/Pi;

      double pseudoscalar_corr = 0.0;
      switch (static_cast<int>(flexibledecay_settings.get(FlexibleDecay_settings::include_higher_order_corrections))) {
         case 4:
            [[fallthrough]];
         case 3:
            [[fallthrough]];
         case 2:
            pseudoscalar_corr += deltaNNLO*alpha_s_red;
            [[fallthrough]];
         case 1:
            pseudoscalar_corr += deltaNLO;
            pseudoscalar_corr *= alpha_s_red/std::norm(0.5*A12_A);
            pseudoscalar_corr += 1. - Sqr(get_alphas(context)/alpha_s_5f);
            pseudoscalar_corr *= Gamma_SM_LO_P;
            break;
         default:
            WARNING("Unknow correcion in Phi->gg");
      }
      result += pseudoscalar_corr;
   }

   if (flexibledecay_settings.get(FlexibleDecay_settings::print_effc_block)) {
      effhiggscouplings_block_input.push_back(
         {
            36 + 10*(in_idx.size()>0 ? in_idx.at(0)-1 : 0), 21, 21,
            std::sqrt(result/(flux * color_fact * ps * ps_symmetry)/(0.5*Power4(mAh))),
            field_as_string<PseudoscalarHiggs>(in_idx) + "-" + field_as_string<Gluon>(out1_idx) + "-" + field_as_string<Gluon>(out2_idx)
         }
      );
   }

   if (static_cast<bool>(flexibledecay_settings.get(FlexibleDecay_settings::call_higgstools)) ||
       static_cast<bool>(flexibledecay_settings.get(FlexibleDecay_settings::call_lilith)) ||
       static_cast<bool>(flexibledecay_settings.get(FlexibleDecay_settings::calculate_normalized_effc))
   ) {
      neutral_higgs_effc.add_coupling(
         field_as_string<Ah>(in_idx),
         {21, 21},
         std::pair<std::string const&, double> {field_as_string<PseudoscalarHiggs>(in_idx) + "-" + field_as_string<Gluon>(out1_idx) + "-" + field_as_string<Gluon>(out2_idx), std::sqrt(result)}
      );
   }

   return result;
}
