
// template specialization for the H -> Fd Fd case
template<>
double CLASSNAME::get_partial_width<H,bar<dq>::type,dq>(
   const context_base& context,
   typename field_indices<H>::type const& indexIn,
   typename field_indices<dq>::type const& indexOut1,
   typename field_indices<dq>::type const& indexOut2
   ) const
{
   // get HBBbar vertex
   // we don't use amplitude_squared here because we need this vertex
   // both with running and pole masses
   const auto indices = concatenate(indexOut1, indexOut2, indexIn);
   const auto HBBbarVertexDR = Vertex<bar<dq>::type, dq, H>::evaluate(indices, context);

   const double mHOS = context.physical_mass<H>(indexIn);
   const double flux = 1./(2.*mHOS);

   constexpr double color_factor = squared_color_generator<H,bar<dq>::type,dq>();

   if(!boost::range::equal(indexOut1, indexOut2)) {
      if (!is_zero(HBBbarVertexDR.left()) || !is_zero(HBBbarVertexDR.right())) {
         const double mdqOS1 = context.physical_mass<dq>(indexOut1);
         const double mdqOS2 = context.physical_mass<dq>(indexOut2);
         const auto xOS1 = Sqr(mdqOS1/mHOS);
         const auto xOS2 = Sqr(mdqOS2/mHOS);
         const double phase_space = 1./(8.*Pi) * std::sqrt(KallenLambda(1., xOS1, xOS2));
         return flux * phase_space * color_factor * amplitude_squared<H, bar<dq>::type, dq>(context, indexIn, indexOut1, indexOut2);
      }
      return 0.;
   }

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

   const double phase_spaceDR = 1./(8.*Pi) * std::sqrt(KallenLambda(1., xDR, xDR));
   const double phase_spaceOS = 1./(8.*Pi) * std::sqrt(KallenLambda(1., xOS, xOS));

   const auto HBBbarVertexDR_S = 0.5*(HBBbarVertexDR.left() + HBBbarVertexDR.right());
   const auto HBBbarVertexDR_P = 0.5*(HBBbarVertexDR.right() - HBBbarVertexDR.left());

   double amp2DR_S = Sqr(mHOS) * Sqr(betaDR) *
                     2*std::norm(HBBbarVertexDR_S);
   double amp2OS_S = Sqr(mHOS) * Sqr(betaOS) *
                     2*std::norm(HBBbarVertexDR_S) * Sqr(mdqOS / mdqDR);

   double amp2DR_P = 0;
   double amp2OS_P = 0;
   if (info::is_CP_violating_Higgs_sector) {
      amp2DR_P = Sqr(mHOS) *
                 2*std::norm(HBBbarVertexDR_P);
      amp2OS_P = Sqr(mHOS) *
                 2*std::norm(HBBbarVertexDR_P) * Sqr(mdqOS / mdqDR);
   }

   switch (include_higher_order_corrections) {
      case SM_higher_order_corrections::enable: {
         double deltaqq_QCD_OS = 0.;
         const int Nf = number_of_active_flavours(qedqcd, mHOS);
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
               throw std::runtime_error ("Error in H->ddbar: Cannot determine the number of active flavours");
         }
         double deltaqq_QCD_DR_S = calc_Deltaqq(alpha_s_red, Nf);
         double deltaqq_QCD_DR_P = deltaqq_QCD_DR_S;

         // 1L QED correction - eq. 21 in FD manual
         const double alpha_red = get_alpha(context)/Pi;
         const double deltaqq_QED_DR = 17./4.*Sqr(dq::electric_charge)*alpha_red;

         deltaqq_QCD_DR_S +=
            2.*(1. - 10.*xDR)/(1-4.*xDR)*(4./3. - std::log(xDR))*alpha_s_red +
            4./3.*alpha_s_red*calc_DeltaH(betaDR);

         const double deltaqq_QCD_OS_S =
            4./3. * alpha_s_red * calc_DeltaH(betaOS);

         const double deltaqq_QED_OS_S =
            alpha_red * Sqr(dq::electric_charge) * calc_DeltaH(betaOS);

         const double mtpole = qedqcd.displayPoleMt();
         const double lt = std::log(Sqr(mHOS/mtpole));
         const double lq = std::log(xDR);
         // eq. 28 of hep-ph/9505358
         const auto Httindices = concatenate(std::array<int, 1> {2}, std::array<int, 1> {2}, indexIn);
         const auto Httbar = Vertex<bar<uq>::type, uq, H>::evaluate(Httindices, context);
         const auto Httbar_S = 0.5*(Httbar.left() + Httbar.right());
         const auto gtHoVEV = Httbar_S/context.mass<uq>({2});
         const auto gbHoVEV = HBBbarVertexDR_S/context.mass<dq>(indexOut1);
         const double deltaPhi2_S = Sqr(alpha_s_red) * std::real(gtHoVEV/gbHoVEV) * (1.57 - 2.0/3.0*lt + 1.0/9.0*Sqr(lq));

         double deltaqq_QCD_OS_P = 0.;
         double deltaqq_QED_OS_P = 0.;
         double deltaPhi2_P = 0.;
         // don't waste time computing it in models without CPV
         if (info::is_CP_violating_Higgs_sector) {

            deltaqq_QCD_DR_P +=
               2.*(1. - 6.*xDR)/(1-4.*xDR)*(4./3. - std::log(xDR))*alpha_s_red +
               4./3.*alpha_s_red*calc_DeltaAH(betaDR);

            deltaqq_QCD_OS_P =
               4./3. * alpha_s_red * calc_DeltaAH(betaOS);

            deltaqq_QED_OS_P =
               alpha_red * Sqr(dq::electric_charge) * calc_DeltaAH(betaOS);

            const auto gbHoVEV_P = HBBbarVertexDR_P/context.mass<dq>(indexOut1);
            if (!is_zero(gbHoVEV_P)) {
               const auto Httbar_P = 0.5*(Httbar.right() - Httbar.left());
               const auto gtHoVEV_P = Httbar_P/context.mass<uq>({2});
               deltaPhi2_P = Sqr(alpha_s_red) * std::real(gtHoVEV_P/gbHoVEV_P) * (3.83 - lt + 1.0/6.0*Sqr(lq));
            }
         }

         amp2DR_S *= 1. + deltaqq_QCD_DR_S + deltaqq_QED_DR + deltaPhi2_S;
         amp2DR_P *= 1. + deltaqq_QCD_DR_P + deltaqq_QED_DR + deltaPhi2_P;
         amp2OS_S *= 1. + deltaqq_QCD_OS_S + deltaqq_QED_OS_S;
         amp2OS_P *= 1. + deltaqq_QCD_OS_P + deltaqq_QED_OS_P;
         break;
      }
      case SM_higher_order_corrections::disable:
         break;
      default:
         WARNING("Unhandled option in H->ddbar decay");
   }

   // low x limit
   double result_DR =
      flux * color_factor * phase_spaceDR * (amp2DR_S + amp2DR_P);
   // high x limit
   double result_OS =
      flux * color_factor * phase_spaceOS * (amp2OS_S + amp2OS_P);

   return (1-4.*xOS)*result_DR + 4*xOS*result_OS;
}
