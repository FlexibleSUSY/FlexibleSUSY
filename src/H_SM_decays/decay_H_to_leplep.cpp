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
   const auto xOS = Sqr(mL1OS / mHOS);
   const auto betaOS = sqrt(1. - 4. * xOS);
   const auto indices = concatenate(indexIn, indexOut1, indexOut2);
   const auto HLLbarVertexDR =
      Vertex<H, bar<lep>::type, lep>::evaluate(indices, context);
   const auto amp2OS = Sqr(mHOS) * Sqr(betaOS) * 2. *
                       std::norm(HLLbarVertexDR.left()) *
                       Sqr(mL1OS / mL1DR);

   // flux * phase space factor * symmetry factor * matrix element^2
   double res = 0.5 * ps * amp2OS / mHOS;

   // higher order corrections

   if (include_higher_order_corrections == SM_higher_order_corrections::enable &&
       !info::is_CP_violating_Higgs_sector) {
      // 1-loop QED corrections
      res *= 1. + get_alpha(context)/Pi*17./4.;
   }

   return res;
}
