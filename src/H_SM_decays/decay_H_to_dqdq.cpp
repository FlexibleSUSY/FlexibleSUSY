// template specialization for the H -> Fd Fd case

template<>
double CLASSNAME::get_partial_width<H,bar<dq>::type,dq>(
   const context_base& context,
   typename field_indices<H>::type const& indexIn,
   typename field_indices<dq>::type const& indexOut1,
   typename field_indices<dq>::type const& indexOut2
   ) const
{
   // TODO: should we take the off-diagonal case at all?
   //       or should this never happen and we should crash
   if(!boost::range::equal(indexOut1, indexOut2))
      return 0.;
//    BOOST_ASSERT_MSG(boost::range::equal(indexOut1, indexOut2),
      // "Template specialization for H -> Fd1 bar[Fd2] is only valid for Fd1 = Fd2"
//    );

   const double mHOS = context.physical_mass<H>(indexIn);
   const double mdqDR = context.mass<dq>(indexOut1);
   const double mdqOS = context.physical_mass<dq>(indexOut1);
   if(is_zero(mdqDR) || is_zero(mdqOS)) {
      throw std::runtime_error("Error in H->ddbar: down quarks cannot be massless");
   }
   const auto xOS = std::pow(mdqOS/mHOS, 2);
   const auto xDR = std::pow(mdqDR/mHOS, 2);
   const auto betaOS = sqrt(1. - 4.*xOS);
   const auto betaDR = sqrt(1. - 4.*xDR);

   // TODO: add off-shell decays?
   if (mHOS < 2.*mdqDR) {
      return 0.;
   }

   const double alpha_s = get_alphas(context);
   const double alpha_s_red = alpha_s/Pi;
   const double Nf = number_of_active_flavours(mHOS);
   const double alpha = get_alpha(context);
   const double alpha_red = alpha/Pi;
   const double mtpole = qedqcd.displayPoleMt();

   const double deltaqqOS =
      4./3. * alpha_s_red * calc_DeltaH(betaOS);
   const double deltaqqDR =
      2.*(1. - 10.*xDR)/(1-4.*xDR)*(4./3. - std::log(xDR))*alpha_s_red +
      4./3. * alpha_s_red * calc_DeltaH(betaDR) +
      calc_deltaqq(alpha_s_red, Nf);

   const double deltaqqDRQED =
      17./4. * alpha_red * pow(dq::electric_charge, 2);
   const double deltaqqOSQED =
      alpha_red * pow(dq::electric_charge, 2) * calc_DeltaH(betaOS);

   // chiral breaking correctios
   // TODO: probably shouldn't be applied in case of CP-breaking
   double deltaH2 = 0.;
   if(!info::CPViolationInHiggsSector) {
      const double lt = std::log(Sqr(mHOS/mtpole));
      const double lq = std::log(Sqr(mdqDR/mHOS));
      deltaH2 = Sqr(alpha_s_red) * (1.57 - 2.0/3.0*lt + 1.0/9.0*Sqr(lq));
   }

   const double flux = 1./(2.*mHOS);
   const double color_factor = 3;
   const double phase_spaceDR = 1./(8.*Pi) * std::sqrt(KallenLambda(1., Sqr(mdqDR/mHOS), Sqr(mdqDR/mHOS)));
   const double phase_spaceOS = 1./(8.*Pi) * std::sqrt(KallenLambda(1., Sqr(mdqOS/mHOS), Sqr(mdqOS/mHOS)));

   // get HBBbar vertex
   // we don't use amplitude_squared here because we need both this vertex
   // both with running and pole masses
   const auto indices = concatenate(indexIn, indexOut1, indexOut2);
   const auto HBBbarVertexDR = Vertex<H, bar<dq>::type, dq>::evaluate(indices, context);

   const auto amp2DR = std::pow(mHOS, 2) * std::pow(betaDR, 2) *
               2.*std::norm(HBBbarVertexDR.left());
   const auto amp2OS = std::pow(mHOS, 2) * std::pow(betaOS, 2) *
                2.*std::norm(HBBbarVertexDR.left()) * std::pow(mdqOS / mdqDR, 2);

   return flux * color_factor *
          (
             // low x limit
             (1 - 4. * xOS) * phase_spaceDR * amp2DR *
                (1. + deltaqqDR + deltaqqDRQED +  deltaH2) +
             // high x limit
             4 * xOS * phase_spaceOS * amp2OS * (1. + deltaqqOS + deltaqqOSQED));
}
