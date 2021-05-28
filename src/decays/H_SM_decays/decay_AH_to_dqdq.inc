
// specialization for the AH -> Fd Fd case

template<>
double CLASSNAME::get_partial_width<AH, bar<dq>::type, dq>(
   const context_base& context,
   typename field_indices<AH>::type const& indexIn,
   typename field_indices<bar<dq>::type>::type const& indexOut1,
   typename field_indices<dq>::type const& indexOut2
   ) const
{
   // get HBBbar vertex
   // we don't use amplitude_squared here because we need both this vertex
   // both with running and pole masses
   const auto indices = concatenate(indexOut1, indexOut2, indexIn);
   const auto HBBbarVertexDR = Vertex<bar<dq>::type, dq, AH>::evaluate(indices, context);

   const double mAOS = context.physical_mass<AH>(indexIn);
   const double flux = 1./(2.*mAOS);

   constexpr double color_factor = squared_color_generator<AH, bar<dq>::type, dq>();

   if(!boost::range::equal(indexOut1, indexOut2)) {
      if (!is_zero(HBBbarVertexDR.left()) || !is_zero(HBBbarVertexDR.right())) {
         const double mdqOS1 = context.physical_mass<dq>(indexOut1);
         const double mdqOS2 = context.physical_mass<dq>(indexOut2);
         const auto xOS1 = Sqr(mdqOS1/mAOS);
         const auto xOS2 = Sqr(mdqOS2/mAOS);
         const double phase_space = 1./(8.*Pi) * std::sqrt(KallenLambda(1., xOS1, xOS2));
         return flux * phase_space * color_factor * amplitude_squared<AH, bar<dq>::type, dq>(context, indexIn, indexOut1, indexOut2);
      }
      return 0.;
   }

   const double mdqDR = context.mass<dq>(indexOut1);
   const double mdqOS = context.physical_mass<dq>(indexOut1);
   if(is_zero(mdqDR) || is_zero(mdqOS)) {
      throw std::runtime_error("Error in A->ddbar: down quarks cannot be massless");
   }
   const auto xOS = Sqr(mdqOS/mAOS);
   const auto xDR = Sqr(mdqDR/mAOS);

   // TODO: add off-shell decays?
   if (4.*std::max(xDR, xOS) > 1.) {
      return 0.;
   }

   const auto betaOS = std::sqrt(1.-4.*xOS);
   const auto betaDR = std::sqrt(1.-4.*xDR);

   const double phase_spaceDR = 1./(8.*Pi) * std::sqrt(KallenLambda(1., xDR, xDR));
   const double phase_spaceOS = 1./(8.*Pi) * std::sqrt(KallenLambda(1., xOS, xOS));

   const auto HBBbarVertexDR_P = 0.5*(HBBbarVertexDR.right() - HBBbarVertexDR.left());

   double amp2DR_P = 0;
   double amp2OS_P = 0;
   amp2DR_P = Sqr(mAOS) *
              2*std::norm(HBBbarVertexDR_P);
   amp2OS_P = Sqr(mAOS) *
              2*std::norm(HBBbarVertexDR_P) * Sqr(mdqOS / mdqDR);

   if (static_cast<int>(flexibledecay_settings.get(FlexibleDecay_settings::include_higher_order_corrections)) > 0) {
         const int Nf = number_of_active_flavours(qedqcd, mAOS);
         double alpha_s_red;
         double Y_conversion = 1.;
         switch (Nf) {
            case 5: {
               auto qedqcd_ = qedqcd;
               qedqcd_.to(mAOS);
               alpha_s_red = qedqcd_.displayAlpha(softsusy::ALPHAS)/Pi;
               Y_conversion = Sqr(sm_down_quark_masses(qedqcd_, indexOut1.at(0))/mdqDR);
               break;
            }
            case 6:
               alpha_s_red = get_alphas(context)/Pi;
               break;
            default:
               throw std::runtime_error("Error in A->ddbar: Cannot determine the number of active flavours");
         }

         double deltaqq_QCD_DR_P = calc_Deltaqq(alpha_s_red, Nf, flexibledecay_settings);

         // 1L QED correction - eq. 17 in FD paper
         const double alpha_red = get_alpha(context)/Pi;
         const double deltaqq_QED_DR = 17./4.*Sqr(dq::electric_charge)*alpha_red;

         deltaqq_QCD_DR_P +=
            2.*(1. - 6.*xDR)/(1-4.*xDR)*(4./3. - std::log(xDR))*alpha_s_red +
            4./3.*alpha_s_red*calc_DeltaAH(betaDR);

         const double deltaqq_QCD_OS_P =
               4./3. * alpha_s_red * calc_DeltaAH(betaOS);

         const double deltaqq_QED_OS_P =
               alpha_red * Sqr(dq::electric_charge) * calc_DeltaAH(betaOS);

         const auto gbHoVEV_P = HBBbarVertexDR_P/context.mass<dq>(indexOut1);
         double deltaPhi2_P = 0.;
         double deltaqq_QCDxQED_DR = 0.;

         if (static_cast<int>(flexibledecay_settings.get(FlexibleDecay_settings::include_higher_order_corrections)) > 1) {
            deltaqq_QCDxQED_DR =
               (691/24. - 6*zeta3 - Sqr(Pi))*Sqr(dq::electric_charge)*alpha_red*alpha_s_red;
            if (!is_zero(gbHoVEV_P)) {
               const double mtpole = qedqcd.displayPoleMt();
               const double lt = std::log(Sqr(mAOS/mtpole));
               const double lq = std::log(xDR);
               const auto Httindices = concatenate(std::array<int, 1> {2}, std::array<int, 1> {2}, indexIn);
               const auto Httbar = Vertex<bar<uq>::type, uq, H>::evaluate(Httindices, context);
               const auto Httbar_P = 0.5*(Httbar.right() - Httbar.left());
               const auto gtHoVEV_P = Httbar_P/context.mass<uq>({2});
               deltaPhi2_P = Sqr(alpha_s_red) * std::real(gtHoVEV_P/gbHoVEV_P) * (23/6. - lt + 1.0/6.0*Sqr(lq));
            }
         }

         amp2DR_P *= Y_conversion*(1. + deltaqq_QCD_DR_P + deltaqq_QED_DR + deltaqq_QCDxQED_DR + deltaPhi2_P);
         amp2OS_P *= 1. + deltaqq_QCD_OS_P + deltaqq_QED_OS_P;
   }

   // low x limit
   double result_DR =
      flux * color_factor * phase_spaceDR * amp2DR_P;
   // high x limit
   double result_OS =
      flux * color_factor * phase_spaceOS * amp2OS_P;

   return (1-4.*xOS)*result_DR + 4*xOS*result_OS;
}
