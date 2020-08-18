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
   if (mAh/(2.*mtpole) < 0.4) {

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

      const double deltaNLO {97./4. - 7./6.*Nf};

      const double log_mAh2OverMT2 {std::log(Sqr(mAh/mtpole))};

      const double deltaNNLO {
         237311./864. - 529./24.*zeta2 - 445./8.*zeta3 + 5.*log_mAh2OverMT2
      };

      const double g3 = context.model.get_g3();
      const double alpha_s_red = Sqr(g3)/(4*Sqr(Pi));

      switch (include_higher_order_corrections) {
         case SM_higher_order_corrections::enable:
            result *= 1. + deltaNLO*alpha_s_red + deltaNNLO*Sqr(alpha_s_red) /*+ deltaNNNLO*Cube(alpha_s_red)*/;
            break;
         case SM_higher_order_corrections::disable:
            break;
         default:
            break;
      }
   }

   return result;
}
