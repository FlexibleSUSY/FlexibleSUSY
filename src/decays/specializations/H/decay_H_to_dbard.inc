template <>
double CLASSNAME::get_partial_width<Higgs,bar<DownTypeQuark>::type,DownTypeQuark>(
   const context_base& context,
   typename field_indices<Higgs>::type const& indexIn,
   typename field_indices<DownTypeQuark>::type const& indexOut1,
   typename field_indices<DownTypeQuark>::type const& indexOut2)
{
   // get HQQbar vertex
   // we don't use amplitude_squared here because we need this vertex
   // both with running and pole masses
   const auto indices = concatenate(indexOut1, indexOut2, indexIn);
   const auto HQQbarVertexDR = Vertex<bar<DownTypeQuark>::type, DownTypeQuark, Higgs>::evaluate(indices, context);

   const double mHOS = context.physical_mass<Higgs>(indexIn);
   const double flux = 1./(2.*mHOS);

   static constexpr double color_factor = squared_color_generator<Higgs,bar<DownTypeQuark>::type,DownTypeQuark>();

   if(!boost::range::equal(indexOut1, indexOut2)) {
      if (!is_zero(HQQbarVertexDR.left()) || !is_zero(HQQbarVertexDR.right())) {
         const double mdqOS1 = context.physical_mass<DownTypeQuark>(indexOut1);
         const double mdqOS2 = context.physical_mass<DownTypeQuark>(indexOut2);
         const auto xOS1 = Sqr(mdqOS1/mHOS);
         const auto xOS2 = Sqr(mdqOS2/mHOS);
         const double phase_space = 1./(8.*Pi) * std::sqrt(KallenLambda(1., xOS1, xOS2));
         auto amp = calculate_amplitude_tree<Higgs, bar<DownTypeQuark>::type, DownTypeQuark>(context, indexIn, indexOut1, indexOut2);
         return flux * phase_space * color_factor * amp.square();
      }
      return 0.;
   }

   const double mdqDR = context.mass<DownTypeQuark>(indexOut1);
   const double mdqOS = context.physical_mass<DownTypeQuark>(indexOut1);
   if(is_zero(mdqDR) || is_zero(mdqOS)) {
      throw std::runtime_error(
         create_process_string<Higgs,bar<DownTypeQuark>::type, DownTypeQuark>(indexIn, indexOut1, indexOut2)
            + ": Down-type quark cannot be massless. Aborting."
      );
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

   const auto HQQbarVertexDR_S = 0.5*(HQQbarVertexDR.left() + HQQbarVertexDR.right());
   const auto HQQbarVertexDR_P = 0.5*(HQQbarVertexDR.right() - HQQbarVertexDR.left());

   double amp2DR_S = Sqr(mHOS) * Sqr(betaDR) *
                     2*std::norm(HQQbarVertexDR_S);
   double amp2OS_S = Sqr(mHOS) * Sqr(betaOS) *
                     2*std::norm(HQQbarVertexDR_S) * Sqr(mdqOS / mdqDR);

   double amp2DR_P = 0;
   double amp2OS_P = 0;
   if (info::is_CP_violating_Higgs_sector) {
      amp2DR_P = Sqr(mHOS) *
                 2*std::norm(HQQbarVertexDR_P);
      amp2OS_P = Sqr(mHOS) *
                 2*std::norm(HQQbarVertexDR_P) * Sqr(mdqOS / mdqDR);
   }

   const int Nf = number_of_active_flavours(qedqcd, mHOS);
   if (Nf < 5) {
      problems.add_warning(
         create_process_string<Higgs,bar<DownTypeQuark>::type, DownTypeQuark>(indexIn, indexOut1, indexOut2)
            + ": Cannot determine the number of active quark flavours. Disabling higher-order corrections."
      );
   }

   if (static_cast<int>(flexibledecay_settings.get(FlexibleDecay_settings::include_higher_order_corrections)) > 0 && Nf >= 5) {
         double alpha_s_red;
         double Y_conversion = 1.;
         switch (Nf) {
            case 5: {
               auto qedqcd_ = qedqcd;
               qedqcd_.to(mHOS);
               alpha_s_red = qedqcd_.displayAlpha(softsusy::ALPHAS)/Pi;
               Y_conversion = Sqr(qedqcd_.displayDownQuarkRunningMass(indexOut1.at(0))/mdqDR);
               break;
            }
            case 6:
               alpha_s_red = get_alphas(context)/Pi;
               break;
            default:
               throw std::runtime_error(
                  create_process_string<Higgs,bar<DownTypeQuark>::type, DownTypeQuark>(indexIn, indexOut1, indexOut2)
                     + ": Cannot determine the number of active quark flavours"
               );
         }
         double deltaqq_QCD_DR_S = calc_Deltaqq(alpha_s_red, Nf, flexibledecay_settings);
         double deltaqq_QCD_DR_P = deltaqq_QCD_DR_S;

         // 1L QED correction - eq. 17 in FD paper
         const double alpha_red = get_alpha(context)/Pi;
         const double deltaqq_QED_DR = 17./4.*Sqr(DownTypeQuark::electricCharge)*alpha_red;

         deltaqq_QCD_DR_S +=
            2.*(1. - 10.*xDR)/(1-4.*xDR)*(4./3. - std::log(xDR))*alpha_s_red +
            4./3.*alpha_s_red*calc_DeltaH(betaDR);

         const double deltaqq_QCD_OS_S =
            4./3. * alpha_s_red * calc_DeltaH(betaOS);

         const double deltaqq_QED_OS_S =
            alpha_red * Sqr(DownTypeQuark::electricCharge) * calc_DeltaH(betaOS);

         double deltaPhi2_S = 0.;
         double deltaqq_QCD_OS_P = 0.;
         double deltaqq_QED_OS_P = 0.;
         double deltaPhi2_P = 0.;
         double deltaqq_QCDxQED_DR = 0.;
         if (static_cast<int>(flexibledecay_settings.get(FlexibleDecay_settings::include_higher_order_corrections)) > 1) {

            deltaqq_QCDxQED_DR =
               deltaqq_QCDxQED*Sqr(DownTypeQuark::electricCharge)*alpha_red*alpha_s_red;

            const double mtpole = qedqcd.displayPoleMt();
            const double lt = std::log(Sqr(mHOS/mtpole));
            const double lq = std::log(xDR);
            const auto Httindices = concatenate(std::array<int, 1> {2}, std::array<int, 1> {2}, indexIn);
            const auto Httbar = Vertex<bar<UpTypeQuark>::type, UpTypeQuark, Higgs>::evaluate(Httindices, context);
            const auto gbHoVEV = HQQbarVertexDR_S/mdqDR;
            if (!is_zero(gbHoVEV)) {
               // eq. 28 of hep-ph/9505358
               const auto Httbar_S = 0.5*(Httbar.left() + Httbar.right());
               const auto gtHoVEV = Httbar_S/context.mass<UpTypeQuark>({2});
               deltaPhi2_S = Sqr(alpha_s_red) * std::real(gtHoVEV/gbHoVEV) * (8/3. - Sqr(Pi/3.) - 2.0/3.0*lt + 1.0/9.0*Sqr(lq));
            }

            // don't waste time computing it in models without CPV
            if (info::is_CP_violating_Higgs_sector) {

               deltaqq_QCD_DR_P +=
                  2.*(1. - 6.*xDR)/(1-4.*xDR)*(4./3. - std::log(xDR))*alpha_s_red +
                  4./3.*alpha_s_red*calc_DeltaAh(betaDR);

               deltaqq_QCD_OS_P =
                  4./3. * alpha_s_red * calc_DeltaAh(betaOS);

               deltaqq_QED_OS_P =
                  alpha_red * Sqr(DownTypeQuark::electricCharge) * calc_DeltaAh(betaOS);

               const auto gbHoVEV_P = HQQbarVertexDR_P/mdqDR;
               if (!is_zero(gbHoVEV_P)) {
                  const auto Httbar_P = 0.5*(Httbar.right() - Httbar.left());
                  const auto gtHoVEV_P = Httbar_P/context.mass<UpTypeQuark>({2});
                  deltaPhi2_P = Sqr(alpha_s_red) * std::real(gtHoVEV_P/gbHoVEV_P) * (23/6. - lt + 1.0/6.0*Sqr(lq));
               }
            }
         }

         amp2DR_S *= Y_conversion*(1. + deltaqq_QCD_DR_S + deltaqq_QED_DR + deltaqq_QCDxQED_DR + deltaPhi2_S);
         amp2DR_P *= Y_conversion*(1. + deltaqq_QCD_DR_P + deltaqq_QED_DR + deltaqq_QCDxQED_DR + deltaPhi2_P);
         amp2OS_S *= 1. + deltaqq_QCD_OS_S + deltaqq_QED_OS_S;
         amp2OS_P *= 1. + deltaqq_QCD_OS_P + deltaqq_QED_OS_P;
   }

   // low x limit
   const double result_DR_S =
      flux * color_factor * phase_spaceDR * amp2DR_S;
   const double result_DR_P =
      flux * color_factor * phase_spaceDR * amp2DR_P;
   // high x limit
   const double result_OS_S =
      flux * color_factor * phase_spaceOS * amp2OS_S;
   const double result_OS_P =
      flux * color_factor * phase_spaceOS * amp2OS_P;

   const double result_S = (1-4.*xOS)*result_DR_S + 4*xOS*result_OS_S;
   const double result_P = (1-4.*xOS)*result_DR_P + 4*xOS*result_OS_P;

   if (flexibledecay_settings.get(FlexibleDecay_settings::call_higgstools) != 0 ||
       flexibledecay_settings.get(FlexibleDecay_settings::call_lilith) != 0 ||
       flexibledecay_settings.get(FlexibleDecay_settings::calculate_normalized_effc) != 0
   ) {
      const double remove_normalization = flux*phase_spaceOS*2*Sqr(mHOS);
      neutral_higgs_effc.add_coupling(field_as_string<Higgs>(indexIn),
         {
            boost::hana::unpack(bar<DownTypeQuark>::pdgids, _to_array<bar<DownTypeQuark>::numberOfGenerations>).at(indexOut1.at(0)),
            boost::hana::unpack(DownTypeQuark::pdgids, _to_array<DownTypeQuark::numberOfGenerations>).at(indexOut2.at(0))
         },
         std::pair<std::string, std::complex<double>> {
            field_as_string<Higgs>(indexIn) + "-" + field_as_string<bar<DownTypeQuark>::type>(indexOut1) + "-" + field_as_string<DownTypeQuark>(indexOut2),
            std::sqrt(result_S/(remove_normalization*Sqr(betaOS))) + std::sqrt(result_P/(remove_normalization))*1i
         },
         result_S + result_P
      );
   }

   return result_S + result_P;
}
