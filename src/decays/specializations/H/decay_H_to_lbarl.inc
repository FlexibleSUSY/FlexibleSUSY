template <>
double CLASSNAME::get_partial_width<Higgs, bar<ChargedLepton>::type, ChargedLepton>(
   const context_base& context,
   typename field_indices<Higgs>::type const& indexIn,
   typename field_indices<bar<ChargedLepton>::type>::type const& indexOut1,
   typename field_indices<ChargedLepton>::type const& indexOut2)
{

   const auto indices = concatenate(indexOut1, indexOut2, indexIn);
   const auto HLLbarVertexDR = Vertex<bar<ChargedLepton>::type, ChargedLepton, Higgs>::evaluate(indices, context);
   const auto HLLbarVertexDR_S = 0.5*(HLLbarVertexDR.left() + HLLbarVertexDR.right());
   const auto HLLbarVertexDR_P = 0.5*(HLLbarVertexDR.right() - HLLbarVertexDR.left());

   if (is_zero(HLLbarVertexDR_S) && is_zero(HLLbarVertexDR_S)) {
      return 0.;
   }

   const double mHOS = context.physical_mass<Higgs>(indexIn);
   const double mL1OS = context.physical_mass<bar<ChargedLepton>::type>(indexOut1);
   const double mL2OS = context.physical_mass<ChargedLepton>(indexOut2);

   // phase space without symmetry factor
   const auto xOS1 = Sqr(mL1OS/mHOS);
   const auto xOS2 = Sqr(mL2OS/mHOS);
   const auto betaOS2_S = 1.0 - Sqr((mL1OS + mL2OS)/mHOS);
   const auto betaOS2_P = 1.0 - Sqr((mL1OS - mL2OS)/mHOS);
   const double ps = 1./(8.*Pi)*std::sqrt(KallenLambda(1., xOS1, xOS2));
   const double flux = 1./(2.*mHOS);

   const double amp2_S = Sqr(mHOS) * betaOS2_S * 2*std::norm(HLLbarVertexDR_S);

   const double amp2_P =
      info::is_CP_violating_Higgs_sector ? Sqr(mHOS) * betaOS2_P * 2.*std::norm(HLLbarVertexDR_P) : 0.;

   // flux * phase space factor * symmetry factor * |matrix element|^2
   double res_S = flux * ps * amp2_S;
   double res_P = flux * ps * amp2_P;

   // higher order corrections
   if (flexibledecay_settings.get(FlexibleDecay_settings::include_higher_order_corrections)) {
      // 1-loop QED corrections
      res_S *= 1. + get_alpha(context)/Pi*17./4.;
      res_P *= 1. + get_alpha(context)/Pi*17./4.;
   }

   if (flexibledecay_settings.get(FlexibleDecay_settings::call_higgstools) != 0 ||
       flexibledecay_settings.get(FlexibleDecay_settings::call_lilith) != 0 ||
       flexibledecay_settings.get(FlexibleDecay_settings::calculate_normalized_effc) != 0
   ) {
      const double remove_normalization = flux*ps*2*Sqr(mHOS);
      neutral_higgs_effc.add_coupling(
         field_as_string<Higgs>(indexIn),
            {
               boost::hana::unpack(bar<ChargedLepton>::pdgids, _to_array<bar<ChargedLepton>::numberOfGenerations>).at(indexOut1.at(0)),
               boost::hana::unpack(ChargedLepton::pdgids, _to_array<ChargedLepton::numberOfGenerations>).at(indexOut2.at(0))
            },
         std::pair<std::string, std::complex<double>> {
            field_as_string<Higgs>(indexIn) + "-" + field_as_string<bar<ChargedLepton>::type>(indexOut1) + "-" + field_as_string<ChargedLepton>(indexOut2),
            std::sqrt(res_S/(remove_normalization*betaOS2_S)) + 1i*std::sqrt(res_P/(remove_normalization*betaOS2_P))
         },
         res_S + res_P
      );
   }

   return res_S + res_P;
}
