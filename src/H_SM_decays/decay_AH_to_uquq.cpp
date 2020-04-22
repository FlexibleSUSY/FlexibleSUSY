// special case for H -> Fu Fu
template<>
double CLASSNAME::get_partial_width<AH,bar<uq>::type,uq>(
   const context_base& context,
   typename field_indices<AH>::type const& indexIn,
   typename field_indices<uq>::type const& indexOut1,
   typename field_indices<uq>::type const& indexOut2
   ) const
{
   // TODO: should we take the off-diagonal case at all?
   //       or should this never happen and we should crash
   if(!boost::range::equal(indexOut1, indexOut2))
      return 0.;
//    BOOST_ASSERT_MSG(boost::range::equal(indexOut1, indexOut2), 
      // "Template specialization for H -> Fu1 bar[Fu2] is only valid for Fu1 = Fu2"
//    );

   const double mAH = context.mass<H>(indexIn);
   const double muq = context.mass<uq>(indexOut1);

   // @todo: add off-shell decays?
   if (mAH < 2.*muq) {
      return 0.;
   }

   // SM expression + pure BSM 1L corrections

   const double g3 = context.model.get_g3();
   const double alpha_s_red = Sqr(g3)/(4*Sqr(Pi));
   const double Nf = number_of_active_flavours(mAH);
   const double mtpole = qedqcd.displayPoleMt();

   double result = 0.;

   const double color_factor = squared_color_generator<H, bar<uq>::type, uq>();
   const double flux = 1./(2.*mAH);
   const double phase_space = 1./(8.*Pi) * std::sqrt(KallenLambda(1., Sqr(muq/mAH), Sqr(muq/mAH)));

   // top-quark needs special treatment
   const double x = 4.*Sqr(muq/mAH);
   if(indexOut1[1] == 2) {
     const double betaT = Sqrt(1 - x); // sqrt(1 - 4*Sqr(mtpole/mass));
     const double Abeta = (1 + Sqr(betaT))
                        * (4*PolyLog(2, (1-betaT)/(1+betaT))
                          + 2*PolyLog(2, (betaT-1)/(1+betaT))
                          - 3*Log((1+betaT)/(1-betaT))*Log(2.0/(1+betaT))
                          - 2*Log((1+betaT)/(1-betaT))*Log(betaT))
                        - 3*betaT*Log(4.0/(1-Sqr(betaT)))
                        - 4*betaT*Log(betaT);

     const double deltaHt = 4.0/3.0 * alpha_s_red * (Abeta/betaT
                          + (3 + 34*Sqr(betaT) - 13*Power(betaT,4))
                            * Log((1+betaT)/(1-betaT)) / (16*Power(betaT,3))
                          + 3.0/(8*Sqr(betaT)) * (7*Sqr(betaT) - 1));

     // @todo: check numerical prefactors
     result =
        flux * phase_space * color_factor
        * Power3(betaT)
        * amplitude_squared<H, bar<uq>::type, uq>(context, indexIn, indexOut1, indexOut2)
        // multiplicative 1-loop correction
        * (1 + deltaHt);

   } else {

   const double deltaqq = calc_deltaqq(alpha_s_red, Nf);
   const double lt = Log(Sqr(mAH/mtpole));
   const double lq = Log(Sqr(muq/mAH));
   const double deltaAH2 = Sqr(alpha_s_red) * (3.83 - lt + 1.0/6.0*Sqr(lq));

   result =
      flux * phase_space * color_factor
      * amplitude_squared<AH, bar<uq>::type, uq>(context, indexIn, indexOut1, indexOut2)
      * (1 + deltaqq + deltaAH2);
   }

   if (result < 0) {
      throw std::runtime_error("Error in Ah->uquq. Partial width < 0.");
   } else {
      return result;
   }
}
