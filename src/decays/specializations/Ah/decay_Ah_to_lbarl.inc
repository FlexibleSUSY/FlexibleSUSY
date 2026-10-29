template <>
double CLASSNAME::get_partial_width<PseudoscalarHiggs, bar<ChargedLepton>::type, ChargedLepton>(
   const context_base& context,
   typename field_indices<PseudoscalarHiggs>::type const& indexIn,
   typename field_indices<bar<ChargedLepton>::type>::type const& indexOut1,
   typename field_indices<ChargedLepton>::type const& indexOut2)
{
   if (indexIn.at(0) < info::number_of_neutral_goldstones) {
      throw OutOfBoundsError("Error in " + create_process_string<PseudoscalarHiggs,bar<ChargedLepton>::type, ChargedLepton>(indexIn, indexOut1, indexOut2) + " decay. Decaying particle is a Goldstone.");
   }

   const double mAhOS = context.physical_mass<PseudoscalarHiggs>(indexIn);
   const double mL1OS = context.physical_mass<bar<ChargedLepton>::type>(indexOut1);
   const double mL2OS = context.physical_mass<ChargedLepton>(indexOut2);

   // phase space without symmetry factor
   const auto xOS1 = Sqr(mL1OS/mAhOS);
   const auto xOS2 = Sqr(mL2OS/mAhOS);
   const double ps = 1./(8.*Pi)*std::sqrt(KallenLambda(1., xOS1, xOS2));
   const double flux = 0.5/mAhOS;
   const auto betaOS2_P = 1.0 - Sqr((mL1OS - mL2OS)/mAhOS);

   // matrix element squared
   const auto indices = concatenate(indexOut1, indexOut2, indexIn);
   const auto HLLbarVertexDR = Vertex<bar<ChargedLepton>::type, ChargedLepton, PseudoscalarHiggs>::evaluate(indices, context);
   const auto HLLbarVertexDR_P = 0.5*(HLLbarVertexDR.right() - HLLbarVertexDR.left());
   const double amp2 = Sqr(mAhOS) * betaOS2_P * 2.*std::norm(HLLbarVertexDR_P);

   // flux * phase space factor * symmetry factor * |matrix element|^2
   double res = flux * ps * amp2;

   // higher order corrections

   if (flexibledecay_settings.get_decay_corrections().sm > 0) {
      // 1-loop QED corrections
      res *= 1. + get_alpha(context)/Pi*17./4.;
   }

   if (flexibledecay_settings.get(FlexibleDecay_settings::call_higgstools) != 0 ||
       flexibledecay_settings.get(FlexibleDecay_settings::call_lilith) != 0 ||
       flexibledecay_settings.get(FlexibleDecay_settings::calculate_normalized_effc) != 0
   ) {
      const double remove_normalization = flux*ps*2*Sqr(mAhOS);
      neutral_higgs_effc.add_coupling(
            field_as_string<PseudoscalarHiggs>(indexIn),
            {
               boost::hana::unpack(bar<ChargedLepton>::pdgids, _to_array<bar<ChargedLepton>::numberOfGenerations>).at(indexOut1.at(0)),
               boost::hana::unpack(ChargedLepton::pdgids, _to_array<ChargedLepton::numberOfGenerations>).at(indexOut2.at(0))
            },
         std::pair<std::string, std::complex<double>> {
            field_as_string<PseudoscalarHiggs>(indexIn) + "-" + field_as_string<bar<ChargedLepton>::type>(indexOut1) + "-" + field_as_string<ChargedLepton>(indexOut2),
            1i*std::sqrt(res/(remove_normalization*betaOS2_P))
         },
         res
      );
   }

   return res;
}
