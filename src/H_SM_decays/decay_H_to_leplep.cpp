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
   const double mL1DR = context.mass<bar<lep>::type>(indexOut1);
   const double mL2OS = context.physical_mass<lep>(indexOut2);

   // phase space without symmetry factor
   const double ps = 1. / (8. * Pi) * std::sqrt(KallenLambda(1., Sqr(mL1OS/mHOS), Sqr(mL2OS/mHOS)));

   // matrix element squared
   const auto xOS = std::pow(mL1OS / mHOS, 2);
   const auto betaOS = sqrt(1. - 4. * xOS);
   const auto indices = concatenate(indexIn, indexOut1, indexOut2);
   const auto HLLbarVertexDR =
      Vertex<H, bar<lep>::type, lep>::evaluate(indices, context);
   const auto amp2OS = std::pow(mHOS, 2) * std::pow(betaOS, 2) * 2. *
                       std::norm(HLLbarVertexDR.left()) *
                       std::pow(mL1OS / mL1DR, 2);

   // flux * phase space factor * symmetry factor * matrix element^2
   const double res = 0.5 * ps * amp2OS / mHOS;

   if (res < 0) {
      throw std::runtime_error("Error in H->ll. Partial width < 0.");
   } else {
      return res;
   }
}
