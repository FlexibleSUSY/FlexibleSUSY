// CP-even Higgs to charged leptons

template <>
double CLASSNAME::get_partial_width<H, bar<lep>::type, lep>(
   const context_base& context,
   typename field_indices<H>::type const& indexIn,
   typename field_indices<bar<lep>::type>::type const& indexOut1,
   typename field_indices<lep>::type const& indexOut2) const
{

   const double mHOS = context.physical_mass<H>(indexIn);
   const double mL1OS = context.physical_mass<bar<lep>::type>(indexOut1);
   const double mL2OS = context.physical_mass<lep>(indexOut2);

   // phase space without symmetry factor
   const auto xOS1 = Sqr(mL1OS/mHOS);
   const auto xOS2 = Sqr(mL2OS/mHOS);
   const double ps = 1./(8.*Pi)*std::sqrt(KallenLambda(1., xOS1, xOS2));

   // matrix element squared
   const auto amp =
      calculate_amplitude<hh, typename bar<lep>::type, lep>(context, indexIn, indexOut1, indexOut2);
   const auto amp2 = amp.square();

   // flux * phase space factor * symmetry factor * matrix element^2
   double res = 0.5 * ps * amp2/mHOS;

   // higher order corrections

   if (include_higher_order_corrections == SM_higher_order_corrections::enable) {
      // 1-loop QED corrections
      res *= 1. + get_alpha(context)/Pi*17./4.;
   }

   return res;
}
