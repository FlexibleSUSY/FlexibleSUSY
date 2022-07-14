template <>
double CLASSNAME::get_partial_width<H, bar<l>::type, l>(
   const context_base& context,
   typename field_indices<H>::type const& indexIn,
   typename field_indices<bar<l>::type>::type const& indexOut1,
   typename field_indices<l>::type const& indexOut2)
{

   const auto indices = concatenate(indexOut1, indexOut2, indexIn);
   const auto HLLbarVertexDR = Vertex<bar<l>::type, l, H>::evaluate(indices, context);
   const auto HLLbarVertexDR_S = 0.5*(HLLbarVertexDR.left() + HLLbarVertexDR.right());
   const auto HLLbarVertexDR_P = 0.5*(HLLbarVertexDR.right() - HLLbarVertexDR.left());

   if (is_zero(HLLbarVertexDR_S) && is_zero(HLLbarVertexDR_S)) {
      return 0.;
   }

   const double mHOS = context.physical_mass<H>(indexIn);
   const double mL1OS = context.physical_mass<bar<l>::type>(indexOut1);
   const double mL2OS = context.physical_mass<l>(indexOut2);

   // phase space without symmetry factor
   const auto xOS1 = Sqr(mL1OS/mHOS);
   const auto xOS2 = Sqr(mL2OS/mHOS);
   const auto betaOS2 = 1.0 - Sqr((mL1OS + mL2OS)/mHOS);
   const double ps = 1./(8.*Pi)*std::sqrt(KallenLambda(1., xOS1, xOS2));

   const double amp2_S = Sqr(mHOS) * betaOS2 * 2*std::norm(HLLbarVertexDR_S);

   const double amp2_P =
      info::is_CP_violating_Higgs_sector ? 2. * (Sqr(mHOS) - Sqr(mL1OS -mL2OS)) * std::norm(HLLbarVertexDR_P) : 0.;

   // flux * phase space factor * symmetry factor * |matrix element|^2
   double res_S = 0.5 * ps * amp2_S/mHOS;
   double res_P = 0.5 * ps * amp2_P/mHOS;

   // higher order corrections
   if (flexibledecay_settings.get(FlexibleDecay_settings::include_higher_order_corrections)) {
      // 1-loop QED corrections
      res_S *= 1. + get_alpha(context)/Pi*17./4.;
      res_P *= 1. + get_alpha(context)/Pi*17./4.;
   }

#ifdef ENABLE_HIGGSTOOLS
   if (static_cast<bool>(flexibledecay_settings.get(FlexibleDecay_settings::call_higgstools))) {
      static constexpr std::array<int, 3> clep_pdgs = {11, 13, 15};
      higgstools_input.add_coupling(
         field_as_string<H>(indexIn), {-clep_pdgs.at(indexOut1.at(0)), clep_pdgs.at(indexOut2.at(0))},
         std::sqrt(res_S) + 1i*std::sqrt(res_P)
      );
   }
#endif

   return res_S + res_P;
}
