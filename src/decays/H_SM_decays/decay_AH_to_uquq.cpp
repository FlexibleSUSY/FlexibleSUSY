// template specialization for the AH -> Fu Fu case

template<>
double CLASSNAME::get_partial_width<AH, bar<uq>::type, uq>(
   const context_base& context,
   typename field_indices<AH>::type const& indexIn,
   typename field_indices<uq>::type const& indexOut1,
   typename field_indices<uq>::type const& indexOut2
   ) const
{
   // get AHBBbar vertex
   // we don't use amplitude_squared here because we need both this vertex
   // both with running and pole masses
   const auto indices = concatenate(indexOut1, indexOut2, indexIn);
   const auto AHBBbarVertexDR = Vertex<bar<uq>::type, uq, AH>::evaluate(indices, context);

   const double mAHOS = context.physical_mass<AH>(indexIn);
   const double flux = 1./(2.*mAHOS);

   constexpr double color_factor = squared_color_generator<AH, bar<uq>::type, uq>();

   if(!boost::range::equal(indexOut1, indexOut2)) {
      if (!is_zero(AHBBbarVertexDR.left()) || !is_zero(AHBBbarVertexDR.right())) {
         const double muqOS1 = context.physical_mass<uq>(indexOut1);
         const double muqOS2 = context.physical_mass<uq>(indexOut2);
         const auto xOS1 = Sqr(muqOS1/mAHOS);
         const auto xOS2 = Sqr(muqOS2/mAHOS);
         const double phase_space = 1./(8.*Pi) * std::sqrt(KallenLambda(1., xOS1, xOS2));
         return flux * phase_space * color_factor * amplitude_squared<AH, bar<uq>::type, uq>(context, indexIn, indexOut1, indexOut2);
      }
      return 0.;
   }

   const double muqDR = context.mass<uq>(indexOut1);
   const double muqOS = context.physical_mass<uq>(indexOut1);
   if(is_zero(muqDR) || is_zero(muqOS)) {
      throw std::runtime_error("Error in AH->uubar: down quarks cannot be massless");
   }
   const auto xOS = Sqr(muqOS/mAHOS);
   const auto xDR = Sqr(muqDR/mAHOS);

   // TODO: add off-shell decays?
   if (4.*std::max(xDR, xOS) > 1.) {
      return 0.;
   }

   const auto betaOS = std::sqrt(1.-4.*xOS);
   const auto betaDR = std::sqrt(1.-4.*xDR);

   const double phase_spaceDR = 1./(8.*Pi) * std::sqrt(KallenLambda(1., xDR, xDR));
   const double phase_spaceOS = 1./(8.*Pi) * std::sqrt(KallenLambda(1., xOS, xOS));

   const auto AHBBbarVertexDR_P = 0.5*(AHBBbarVertexDR.right() - AHBBbarVertexDR.left());

   double amp2DR_P = 0;
   double amp2OS_P = 0;
   amp2DR_P = Sqr(mAHOS) *
              2*std::norm(AHBBbarVertexDR_P);
   amp2OS_P = Sqr(mAHOS) *
              2*std::norm(AHBBbarVertexDR_P) * Sqr(muqOS / muqDR);

   if (static_cast<int>(flexibledecay_settings.get(FlexibleDecay_settings::include_higher_order_corrections)) > 0) {
         const int Nf = number_of_active_flavours(qedqcd, mAHOS);
         double alpha_s_red;
         double Y_conversion = 1.;
         switch (Nf) {
            case 5: {
               auto qedqcd_ = qedqcd;
               qedqcd_.to(mAHOS);
               alpha_s_red = qedqcd_.displayAlpha(softsusy::ALPHAS)/Pi;
               Y_conversion = Sqr(sm_up_quark_masses(qedqcd_, indexOut1.at(0))/muqDR);
               break;
            }
            case 6:
               alpha_s_red = get_alphas(context)/Pi;
               break;
            default:
               throw std::runtime_error ("Error in AH->uubar: Cannot determine the number of active flavours");
         }
         double deltaqq_QCD_DR_P =
               calc_Deltaqq(alpha_s_red, Nf, flexibledecay_settings)
               + 2.*(1. - 6.*xDR)/(1-4.*xDR)*(4./3. - std::log(xDR))*alpha_s_red
               + 4./3.*alpha_s_red*calc_DeltaAH(betaDR);

         // 1L QED correction - eq. 21 in FD manual
         const double alpha_red = get_alpha(context)/Pi;
         const double deltaqq_QED_DR = 17./4.*Sqr(uq::electric_charge)*alpha_red;

         const double deltaqq_QCD_OS_P =
               4./3. * alpha_s_red * calc_DeltaAH(betaOS);

         const double deltaqq_QED_OS_P =
               alpha_red * Sqr(dq::electric_charge) * calc_DeltaAH(betaOS);

         double deltaqq_QCDxQED_DR = 0.;
         double deltaPhi2_P = 0.;
         if (static_cast<int>(flexibledecay_settings.get(FlexibleDecay_settings::include_higher_order_corrections)) > 1) {
            deltaqq_QCDxQED_DR =
               (691/24. - 6*zeta3 - Sqr(Pi))*Sqr(dq::electric_charge)*alpha_red*alpha_s_red;
            if ((indexOut1.at(0) < 2 || indexOut2.at(0) < 2)) {
               const double mtpole = qedqcd.displayPoleMt();
               const double lt = std::log(Sqr(mAHOS/mtpole));
               const double lq = std::log(xDR);
               // eq. 28 of hep-ph/9505358
               const auto AHttindices = concatenate(std::array<int, 1> {2}, std::array<int, 1> {2}, indexIn);
               const auto AHttbar = Vertex<bar<uq>::type, uq, AH>::evaluate(AHttindices, context);
               const auto CSuu = AHBBbarVertexDR_P/context.mass<uq>(indexOut1);
               if (!is_zero(CSuu)) {
                  const auto AHttbar_P = 0.5*(AHttbar.right() - AHttbar.left());
                  const auto CStu = AHttbar_P/context.mass<Fu>({2});
                  deltaPhi2_P = Sqr(alpha_s_red) * std::real(CStu/CSuu) * (23/6. - lt + 1.0/6.0*Sqr(lq));
               }
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
