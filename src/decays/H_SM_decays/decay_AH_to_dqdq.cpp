
// specialization for the AH -> Fd Fd case

template<>
double CLASSNAME::get_partial_width<AH, bar<dq>::type, dq>(
   const context_base& context,
   typename field_indices<AH>::type const& indexIn,
   typename field_indices<bar<dq>::type>::type const& indexOut1,
   typename field_indices<dq>::type const& indexOut2
   ) const
{
   // TODO: should we take the off-diagonal case at all?
   //       or should this never happen and we should crash
   if(!boost::range::equal(indexOut1, indexOut2)) {
      WARNING("Flavour violating decays of Ah->ddbar are currently not implemented");
      return 0.;
   }

   const double mAH = context.physical_mass<AH>(indexIn);
   const double mdq = context.physical_mass<dq>(indexOut1);

   // TODO: add off-shell decays?
   if (mAH < 2.*mdq) {
      return 0.;
   }

   const double g3 = context.model.get_g3();
   const double Nf = number_of_active_flavours(mAH);
   const double alpha_s_red = Sqr(g3)/(4*Sqr(Pi));
   const double mtpole = qedqcd.displayPoleMt();

   const double deltaqq = calc_Deltaqq(alpha_s_red, Nf);

   const double lt = std::log(Sqr(mAH/mtpole));
   const double lq = std::log(Sqr(mdq/mAH));
   // eq. 2.4 of hep-ph/0503173
   const double deltaAH2 = Sqr(alpha_s_red) * (3.83 - lt + 1.0/6.0*Sqr(lq));

   const double flux = 1./(2.*mAH);
   const double phase_space = 1./(8.*Pi) * std::sqrt(KallenLambda(1., Sqr(mdq/mAH), Sqr(mdq/mAH)));
   constexpr double color_factor = 3;

   const double result = flux * phase_space * color_factor *
      amplitude_squared<AH, bar<dq>::type, dq>(context, indexIn, indexOut1, indexOut2) *
      (1. + deltaqq);

   return result;
}
