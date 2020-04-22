template <>
double CLASSNAME::get_partial_width<AH, G, G>(
      const context_base& context,
      const typename field_indices<H>::type& in_idx,
      const typename field_indices<G>::type& out1_idx,
      const typename field_indices<G>::type& out2_idx) const
{
   const auto amp = calculate_amplitude<AH, G, G>(context, in_idx, out1_idx, out2_idx);
   const double mH = context.physical_mass<AH>(in_idx);
   constexpr double ps {1./(8.*Pi)};
   constexpr double ps_symmetry {1./2.};
   const double color_fact = squared_color_generator<AH, G, G>();
   const double flux = 0.5/mH;

   double result = flux * color_fact * ps * ps_symmetry * amp.square();

   // higher order QCD corrections

   // number of active light flavours
   unsigned int Nf;
   const double mtpole {qedqcd.displayPoleMt()};
   if (mH > 5 && mH < mtpole) {
      Nf = 5;
   } else if (mH > mtpole) {
      Nf = 6;
   } else {
      throw;
   }

   // NLO, NNLO - eq. 11 of hep-ph/9705240
   const double deltaNLO {95./4. - 7./6.*Nf};

   const double mtpole2 {Sqr(mtpole)};
   using boost::math::zeta;
   const double deltaNNLO {
      149533./288. - 363./8.*zeta(2) - 495./8.*zeta(3) + 19./8.*log(mH*mH/mtpole2)
            + Nf*(-4157./72 + 11./2.*zeta(2) + 5./4.*zeta(3) + 2./3.*log(mH*mH/mtpole2))
            + Nf*Nf*(127./108. - zeta(2)/6.)
   };
   // eq. 4.20 from Adam's thesis
   const double deltaNNNLO {467.683620788 + 122.440972222*log(mH*mH/mtpole2) + 10.9409722222*Sqr(log(mH*mH/mtpole2))};


   const double g3 = context.model.get_g3();
   const double alpha_s_red = Sqr(g3)/(4*Sqr(Pi));

   switch (include_higher_order_corrections) {
      case HigherOrderSMCorrections::enable:
         result *= 1. + deltaNLO*alpha_s_red + deltaNNLO*Sqr(alpha_s_red) + deltaNNNLO*Cube(alpha_s_red);
      case HigherOrderSMCorrections::disable:
         result *= 1.;
      default:
         break;
   }

   if (result < 0) {
      throw std::runtime_error("Width < 0");
   } else {
      return result;
   }
}
