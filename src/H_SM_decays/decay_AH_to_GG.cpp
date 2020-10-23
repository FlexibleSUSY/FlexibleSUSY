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
   const double mtpole {qedqcd.displayPoleMt()};
   const double tau = Sqr(mAh/(2.*context.mass<uq>({2})));
   if (tau < 0.7) {

      // number of active light flavours
      unsigned int Nf;
      if (mAh > 5 && mAh < mtpole) {
         Nf = 5;
      } else if (mAh > mtpole) {
         Nf = 6;
      } else {
         throw std::runtime_error(
            "Error in "
            + create_process_string<AH, G, G>(in_idx, out1_idx, out2_idx)
            + ": Could not determine the number of active light quark flavours"
         );
      }

      const auto indices = concatenate(in_idx, std::array<int, 1> {2}, std::array<int, 1> {2});
      const auto AHGGVertex = Vertex<AH, bar<uq>::type, uq>::evaluate(indices, context);
      std::complex<double> const AHGGVertexVal = 0.5*(-AHGGVertex.left() + AHGGVertex.right());

      const double tau = Sqr(mAh/(2.*context.mass<uq>({2})));

      const std::complex<double> Ff = -2.*f(tau)/std::sqrt(tau);
      // LO width comming only from the top-loop
      // agrees up to a full double precision with autmatically generated one
      const double Gamma_SM_LO = mAh/(32.*Power3(Pi))*std::norm(get_alphas(context)*AHGGVertexVal*Ff);

      const double mu = mAh;
      const double LH = std::log(Sqr(mu/mAh));
      const double deltaNLO {
         97./4. - 7./6.*Nf + (33.-2*Nf)/6*LH
      };

      const double log_mAh2OverMT2 {std::log(Sqr(mAh/mtpole))};
      const double deltaNNLO {
         237311./864. - 529./24.*zeta2 - 445./8.*zeta3 + 5.*log_mAh2OverMT2
      };

      const double g3 = context.model.get_g3();
      const double alpha_s_red = Sqr(g3/(2*Pi));

      switch (include_higher_order_corrections) {
         case SM_higher_order_corrections::enable:
            result += Gamma_SM_LO*(deltaNLO*alpha_s_red + deltaNNLO*Sqr(alpha_s_red));
            break;
         case SM_higher_order_corrections::disable:
            break;
         default:
            break;
      }
   }

   return result;
}
