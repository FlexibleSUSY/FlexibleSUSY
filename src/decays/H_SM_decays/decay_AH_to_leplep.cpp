// CP-odd Higgs to charged leptons

template <>
double CLASSNAME::get_partial_width<AH, bar<lep>::type, lep>(
   const context_base& context,
   typename field_indices<AH>::type const& indexIn,
   typename field_indices<bar<lep>::type>::type const& indexOut1,
   typename field_indices<lep>::type const& indexOut2) const
{

   const double mAHOS = context.physical_mass<AH>(indexIn);
   const double mL1OS = context.physical_mass<bar<lep>::type>(indexOut1);
   const double mL2OS = context.physical_mass<lep>(indexOut2);

   // phase space without symmetry factor
   const auto xOS1 = Sqr(mL1OS/mAHOS);
   const auto xOS2 = Sqr(mL2OS/mAHOS);
   const double ps = 1./(8.*Pi)*std::sqrt(KallenLambda(1., xOS1, xOS2));

   // matrix element squared
   const auto amp =
      calculate_amplitude<AH, typename bar<lep>::type, lep>(context, indexIn, indexOut1, indexOut2);
   const auto amp2 = amp.square();

   // flux * phase space factor * symmetry factor * |matrix element|^2
   double res = 0.5 * ps * amp2/mAHOS;

   // higher order corrections

   if (include_higher_order_corrections == SM_higher_order_corrections::enable) {
      // 1-loop QED corrections
      res *= 1. + get_alpha(context)/Pi*17./4.;
   }

   return res;
}
