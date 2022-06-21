// CP-odd Higgs to charged leptons

template <>
double CLASSNAME::get_partial_width<AH, bar<lep>::type, lep>(
   const context_base& context,
   typename field_indices<AH>::type const& indexIn,
   typename field_indices<bar<lep>::type>::type const& indexOut1,
   typename field_indices<lep>::type const& indexOut2)
{
   if (indexIn.at(0) < info::number_of_neutral_goldstones) {
      throw OutOfBoundsError("Error in " + create_process_string<AH,bar<lep>::type, lep>(indexIn, indexOut1, indexOut2) + " decay. Decaying particle is a Goldstone.");
   }

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

   if (FlexibleDecay_settings::include_higher_order_corrections) {
      // 1-loop QED corrections
      res *= 1. + get_alpha(context)/Pi*17./4.;
   }

#ifdef ENABLE_HIGGSTOOLS
   if (static_cast<bool>(flexibledecay_settings.get(FlexibleDecay_settings::call_higgstools))) {
      const std::array<int, 3> clep_pdgs = {11, 13, 15};
      higgstools_input.push_back(
         std::make_tuple(
            field_as_string<AH>(indexIn),
            clep_pdgs.at(indexOut1.at(0)), -clep_pdgs.at(indexOut2.at(0)),
            mAHOS,
            1i*std::sqrt(res)
         )
      );
   }
#endif

   return res;
}
