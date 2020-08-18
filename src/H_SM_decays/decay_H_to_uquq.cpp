// template specialization for the H -> Fu Fu case

template<>
double CLASSNAME::get_partial_width<H,bar<uq>::type,uq>(
   const context_base& context,
   typename field_indices<H>::type const& indexIn,
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

   const double mHOS = context.physical_mass<H>(indexIn);
   const double muq = context.mass<uq>(indexOut1);

   // TODO: add off-shell decays?
   if (mHOS < 2.*muq) {
      return 0.;
   }

   // SM expression + pure BSM 1L corrections
//    return amplitude_squared<H, bar<uq>::type, uq>(context, indexIn, indexOut1, indexOut2);

   const double g3 = context.model.get_g3();
   const double alpha_s_red = Sqr(g3)/(4*Sqr(Pi));
   const double Nf = number_of_active_flavours(mHOS);
   const double mtpole = qedqcd.displayPoleMt();

   double result = 0.;

   const double flux = 1./(2.*mHOS);
   const double phase_space = 1./(8.*Pi) * std::sqrt(KallenLambda(mHOS*mHOS, muq*muq, muq*muq))/(mHOS*mHOS);
   constexpr double color_factor = squared_color_generator<H, bar<uq>::type, uq>();

   // top-quark needs special treatment
   const double x = 4.*Sqr(muq/mHOS);
   if(indexOut1[1] == 2) {
     const double betaT = Sqrt(1 - x);//Sqrt(1 - 4*Sqr(mtpole/mass));
     const double log_ratio {std::log((1 + betaT) / (1 - betaT))};
     const double Abeta = (1 + Sqr(betaT))
                        * (4*dilog((1-betaT)/(1+betaT))
                          + 2*dilog((betaT-1)/(1+betaT))
                          - 3*log_ratio*std::log(2.0/(1+betaT))
                          - 2*log_ratio*std::log(betaT))
                        - 3*betaT*std::log(4.0/(1-Sqr(betaT)))
                        - 4*betaT*std::log(betaT);

     const double deltaHt = 4.0/3.0 * alpha_s_red * (Abeta/betaT
                          + (3 + 34*Sqr(betaT) - 13*Power(betaT,4))
                            * log_ratio / (16*Power(betaT,3))
                          + 3.0/(8*Sqr(betaT)) * (7*Sqr(betaT) - 1));

     // @todo: check numerical prefactors
     result =
        flux * phase_space * color_factor
        * Power(betaT, 3)
        * calculate_amplitude<H, bar<uq>::type, uq>(context, indexIn, indexOut1, indexOut2).square()
        // multiplicative 1-loop correction
        * (1 + deltaHt);
   } else {

   const double deltaqq = calc_deltaqq(alpha_s_red, Nf);
   const double lt = std::log(Sqr(mHOS/mtpole));
   const double lq = std::log(Sqr(muq/mHOS));
   const double deltaH2 = Sqr(alpha_s_red) * (1.57 - 2.0/3.0*lt + 1.0/9.0*Sqr(lq));

   result = flux * phase_space * color_factor *
      amplitude_squared<H, bar<uq>::type, uq>(context, indexIn, indexOut1, indexOut2)
      * (1 + deltaqq + deltaH2);
   }

   return result;
}
