// template specialization for the H -> Fu Fu case

template<>
double CLASSNAME::get_partial_width<H, bar<uq>::type, uq>(
   const context_base& context,
   typename field_indices<H>::type const& indexIn,
   typename field_indices<uq>::type const& indexOut1,
   typename field_indices<uq>::type const& indexOut2
   ) const
{
   // TODO: should we take the off-diagonal case at all?
   //       or should this never happen and we should crash
   if(!boost::range::equal(indexOut1, indexOut2)) {
      WARNING("Flavour violating decays of H->uubar currently not implemented!");
      return 0.;
   }

   const double mHOS = context.physical_mass<H>(indexIn);
   const double muqDR = context.mass<uq>(indexOut1);
   const double muqOS = context.physical_mass<uq>(indexOut1);
   if(is_zero(muqDR) || is_zero(muqOS)) {
      throw std::runtime_error("Error in H->uubar: down quarks cannot be massless");
   }
   const auto xOS = Sqr(muqOS/mHOS);
   const auto xDR = Sqr(muqDR/mHOS);

   // TODO: add off-shell decays?
   if (4.*std::max(xDR, xOS) > 1.) {
      return 0.;
   }

   const auto betaOS = std::sqrt(1.-4.*xOS);
   const auto betaDR = std::sqrt(1.-4.*xDR);

   const double flux = 1./(2.*mHOS);
   constexpr double color_factor = 3;
   const double phase_spaceDR = 1./(8.*Pi) * std::sqrt(KallenLambda(1., xDR, xDR));
   const double phase_spaceOS = 1./(8.*Pi) * std::sqrt(KallenLambda(1., xOS, xOS));

   // get HBBbar vertex
   // we don't use amplitude_squared here because we need both this vertex
   // both with running and pole masses
   const auto indices = concatenate(indexIn, indexOut1, indexOut2);
   const auto HBBbarVertexDR = Vertex<H, bar<uq>::type, uq>::evaluate(indices, context);

   const auto amp2DR = Sqr(mHOS) * Sqr(betaDR) *
               2.*std::norm(HBBbarVertexDR.left());
   const auto amp2OS = Sqr(mHOS) * Sqr(betaOS) *
                2.*std::norm(HBBbarVertexDR.left()) * Sqr(muqOS / muqDR);

   // low x limit
   double result_DR =
      flux * color_factor * phase_spaceDR * amp2DR;
   // high x limit
   double result_OS =
      flux * color_factor * phase_spaceOS * amp2OS;

   if (include_higher_order_corrections == SM_higher_order_corrections::enable) {

      const double alpha_s_red = get_alphas(context)/Pi;
      const double Nf = number_of_active_flavours(mHOS);
      const double alpha_red = get_alpha(context)/Pi;
      const double mtpole = qedqcd.displayPoleMt();

      double deltaqqOS = 0.;
      double deltaqqDR = 0.;
      double deltaqqDRQED = 0.;
      double deltaqqOSQED = 0.;
      // chirality breaking corrections
      double deltaH2 = 0.;

      if(!info::is_CP_violating_Higgs_sector) {
         deltaqqOS =
            4./3. * alpha_s_red * calc_DeltaH(betaOS);
         deltaqqDR =
            2.*(1. - 10.*xDR)/(1-4.*xDR)*(4./3. - std::log(xDR))*alpha_s_red +
            4./3. * alpha_s_red * calc_DeltaH(betaDR) +
            calc_deltaqq(alpha_s_red, Nf);

         deltaqqDRQED =
            17./4. * alpha_red * Sqr(uq::electric_charge);
         deltaqqOSQED =
            alpha_red * Sqr(uq::electric_charge) * calc_DeltaH(betaOS);

         const double lt = std::log(Sqr(mHOS/mtpole));
         const double lq = std::log(xDR);
         deltaH2 = Sqr(alpha_s_red) * (1.57 - 2.0/3.0*lt + 1.0/9.0*Sqr(lq));
      }

      result_DR *= 1. + deltaqqDR + deltaqqDRQED + deltaH2;
      result_OS *= 1. + deltaqqOS + deltaqqOSQED;
   }

   return (1-4.*xOS)*result_DR + 4*xOS*result_OS;
}
