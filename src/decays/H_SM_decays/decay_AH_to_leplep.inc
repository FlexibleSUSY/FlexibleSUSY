// CP-odd Higgs to charged leptons

template <>
double CLASSNAME::get_partial_width<PseudoscalarHiggs, bar<ChargedLepton>::type, ChargedLepton>(
   const context_base& context,
   typename field_indices<PseudoscalarHiggs>::type const& indexIn,
   typename field_indices<bar<ChargedLepton>::type>::type const& indexOut1,
   typename field_indices<ChargedLepton>::type const& indexOut2) const
{
   if (indexIn.at(0) < info::number_of_neutral_goldstones) {
      throw OutOfBoundsError("Error in " + create_process_string<PseudoscalarHiggs,bar<ChargedLepton>::type, ChargedLepton>(indexIn, indexOut1, indexOut2) + " decay. Decaying particle is a Goldstone.");
   }

   const double mAHOS = context.physical_mass<PseudoscalarHiggs>(indexIn);
   const double mL1OS = context.physical_mass<bar<ChargedLepton>::type>(indexOut1);
   const double mL2OS = context.physical_mass<ChargedLepton>(indexOut2);

   // phase space without symmetry factor
   const auto xOS1 = Sqr(mL1OS/mAHOS);
   const auto xOS2 = Sqr(mL2OS/mAHOS);
   const double ps = 1./(8.*Pi)*std::sqrt(KallenLambda(1., xOS1, xOS2));

   // matrix element squared
   const auto amp =
      calculate_amplitude<PseudoscalarHiggs, typename bar<ChargedLepton>::type, ChargedLepton>(context, indexIn, indexOut1, indexOut2);
   const auto amp2 = amp.square();

   // flux * phase space factor * symmetry factor * |matrix element|^2
   double res = 0.5 * ps * amp2/mAHOS;

   // higher order corrections

   if (FlexibleDecay_settings::include_higher_order_corrections) {
      // 1-loop QED corrections
      res *= 1. + get_alpha(context)/Pi*17./4.;
   }

   return res;
}
