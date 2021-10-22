template <>
double CLASSNAME::get_partial_width<AH, G, G>(
      const context_base& context,
      const typename field_indices<AH>::type& in_idx,
      const typename field_indices<G>::type& out1_idx,
      const typename field_indices<G>::type& out2_idx) const
{
   const auto amp = calculate_amplitude<AH, G, G>(context, in_idx, out1_idx, out2_idx);
   const double mAh = context.physical_mass<AH>(in_idx);
   constexpr double ps {1./(8.*Pi)};
   constexpr double ps_symmetry {1./2.};
   constexpr double color_fact = squared_color_generator<AH, G, G>();
   const double flux = 0.5/mAh;

   double result = flux * color_fact * ps * ps_symmetry * amp.square();

   // higher order QCD corrections
   const double tau = Sqr(mAh/(2.*context.mass<uq>({2})));
   if (static_cast<int>(flexibledecay_settings.get(FlexibleDecay_settings::include_higher_order_corrections)) && tau < 0.7) {
      // number of active light flavours
      constexpr int Nf = 5;
      auto qedqcd_ = qedqcd;
      qedqcd_.to(mAh);
      // 5-flavour SM alpha_s
      const double alpha_s_5f = qedqcd_.displayAlpha(softsusy::ALPHAS);

      const auto indices = concatenate(std::array<int, 1> {2}, std::array<int, 1> {2}, in_idx);
      const auto AHGGVertex = Vertex<bar<uq>::type, uq, AH>::evaluate(indices, context);
      std::complex<double> const AHGGVertexVal = 0.5*(-AHGGVertex.left() + AHGGVertex.right());

      const double tau = Sqr(mAh/(2.*context.mass<uq>({2})));

      const std::complex<double> A12_A = 2.*f(tau)/tau;
      // LO width comming only from the top-loop
      // agrees up to a full double precision with autmatically generated one
      const double Gamma_SM_LO_P = mAh/(18.*Power3(Pi))*std::norm(alpha_s_5f * AHGGVertexVal*sqrt(tau) * 3./4*A12_A);

      const double mu = mAh;
      const double LH = std::log(Sqr(mu/mAh));
      const double deltaNLO {
         97./4. - 7./6.*Nf + (33.-2*Nf)/6*LH
      };

      const double mtpole {qedqcd.displayPoleMt()};
      const double Lt = std::log(Sqr(mu/mtpole));
      const double deltaNNLO {
         237311./864. - 529./24.*zeta2 - 445./8.*zeta3 + 5.*Lt
      };

      const double alpha_s_red = alpha_s_5f/Pi;

      double pseudoscalar_corr = 0.0;
      switch (static_cast<int>(flexibledecay_settings.get(FlexibleDecay_settings::include_higher_order_corrections))) {
         case 4:
         case 3:
         case 2:
            pseudoscalar_corr += deltaNNLO*alpha_s_red;
         case 1:
            pseudoscalar_corr += deltaNLO;
            pseudoscalar_corr *= alpha_s_red/std::norm(0.5*A12_A);
            pseudoscalar_corr += 1. - Sqr(get_alphas(context)/alpha_s_5f);
            pseudoscalar_corr *= Gamma_SM_LO_P;
            break;
         default:
            WARNING("Unknow correcion in Phi->gg");
      }
   }

   return result;
}
