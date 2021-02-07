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
   if(!boost::range::equal(indexOut1, indexOut2)) {
      WARNING("Flavour violating decays of H->ddbar currently not implemented!");
      return 0.;
   }

   const double mHOS = context.physical_mass<H>(indexIn);
   const double mdqDR = context.mass<dq>(indexOut1);
   const double mdqOS = context.physical_mass<dq>(indexOut1);
   if(is_zero(mdqDR) || is_zero(mdqOS)) {
      throw std::runtime_error("Error in H->ddbar: down quarks cannot be massless");
   }
   const auto xOS = Sqr(mdqOS/mHOS);
   const auto xDR = Sqr(mdqDR/mHOS);

   // TODO: add off-shell decays?
   if (4.*std::max(xDR, xOS) > 1.) {
      return 0.;
   }

   const auto betaOS = std::sqrt(1.-4.*xOS);
   const auto betaDR = std::sqrt(1.-4.*xDR);

   const double flux = 1./(2.*mHOS);
   constexpr double color_factor = squared_color_generator<H,bar<dq>::type,dq>();
   const double phase_spaceDR = 1./(8.*Pi) * std::sqrt(KallenLambda(1., xDR, xDR));
   const double phase_spaceOS = 1./(8.*Pi) * std::sqrt(KallenLambda(1., xOS, xOS));

   // get HBBbar vertex
   // we don't use amplitude_squared here because we need both this vertex
   // both with running and pole masses
   const auto indices = concatenate(indexOut1, indexOut2, indexIn);
   const auto HBBbarVertexDR = Vertex<bar<dq>::type, dq, H>::evaluate(indices, context);
   const auto HBBbarVertexDRV = HBBbarVertexDR.left() + HBBbarVertexDR.right();

   const auto amp2DR = Sqr(mHOS) * Sqr(betaDR) *
               2.*std::norm(HBBbarVertexDR.left());
   const auto amp2OS = Sqr(mHOS) * Sqr(betaOS) *
                2.*std::norm(HBBbarVertexDR.left()) * Sqr(mdqOS / mdqDR);

   // low x limit
   double result_DR =
      flux * color_factor * phase_spaceDR * amp2DR;
   // high x limit
   double result_OS =
      flux * color_factor * phase_spaceOS * amp2OS;

   switch (include_higher_order_corrections) {
      case SM_higher_order_corrections::enable: {
         double deltaqq_QCD_OS = 0.;
         const int Nf = number_of_active_flavours(mHOS);
         double alpha_s_red;
         switch (Nf) {
            case 5: {
               auto qedqcd_ = qedqcd;
               qedqcd_.to(mHOS);
               alpha_s_red = qedqcd_.displayAlpha(softsusy::ALPHAS)/Pi;
               break;
            }
            case 6:
               alpha_s_red = get_alphas(context)/Pi;
               break;
            default:
               ERROR("Error in H->ddbar: Cannot determine the number of active flavours");
               exit(1);
         }
         double deltaqq_QCD_DR = calc_Deltaqq(alpha_s_red, Nf);

         // eq. 21 in FD manual
         const double alpha_red = get_alpha(context)/Pi;
         const double deltaqq_QED_DR = 17./4.*Sqr(dq::electric_charge)*alpha_red;

         double deltaqq_QED_OS = 0.;
         // chirality breaking corrections
         double deltaH2 = 0.;

         if(!info::is_CP_violating_Higgs_sector) {

            deltaqq_QCD_OS =
               4./3. * alpha_s_red * calc_DeltaH(betaOS);
            deltaqq_QCD_DR +=
               2.*(1. - 10.*xDR)/(1-4.*xDR)*(4./3. - std::log(xDR))*alpha_s_red +
               4./3.*alpha_s_red*calc_DeltaH(betaDR);

            deltaqq_QED_OS =
               alpha_red * Sqr(dq::electric_charge) * calc_DeltaH(betaOS);

            const double mtpole = qedqcd.displayPoleMt();
            const double lt = std::log(Sqr(mHOS/mtpole));
            const double lq = std::log(xDR);
            // eq. 28 of hep-ph/9505358
            const auto Httindices = concatenate(std::array<int, 1> {2}, std::array<int, 1> {2}, indexIn);
            const auto Httbar = Vertex<bar<uq>::type, uq, H>::evaluate(Httindices, context);
            const auto HttbarV = Httbar.left() + Httbar.right();
            const auto gtHoVEV = HttbarV/context.mass<uq>({2});
            const auto gbHoVEV = HBBbarVertexDRV/context.mass<dq>(indexOut1);
            deltaH2 = Sqr(alpha_s_red) * std::real(gtHoVEV/gbHoVEV) * (1.57 - 2.0/3.0*lt + 1.0/9.0*Sqr(lq));
         }

         result_DR *= 1. + deltaqq_QCD_DR + deltaqq_QED_DR + deltaH2;
         result_OS *= 1. + deltaqq_QCD_OS + deltaqq_QED_OS;
         break;
      }
      case SM_higher_order_corrections::disable:
         break;
      default:
         WARNING("Unhandled option in H->ddbar decay");
   }

   return (1-4.*xOS)*result_DR + 4*xOS*result_OS;
}
