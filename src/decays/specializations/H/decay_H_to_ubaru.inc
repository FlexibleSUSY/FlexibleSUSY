template <>
double CLASSNAME::get_partial_width<Higgs, bar<UpTypeQuark>::type, UpTypeQuark>(
   const context_base& context,
   typename field_indices<Higgs>::type const& indexIn,
   typename field_indices<UpTypeQuark>::type const& indexOut1,
   typename field_indices<UpTypeQuark>::type const& indexOut2)
{
   // get HQQbar vertex
   // we don't use amplitude_squared here because we need both this vertex
   // both with running and pole masses
   const auto indices = concatenate(indexOut1, indexOut2, indexIn);
   const auto HQQbarVertexDR = Vertex<bar<UpTypeQuark>::type, UpTypeQuark, Higgs>::evaluate(indices, context);

   const double mHOS = context.physical_mass<Higgs>(indexIn);
   const double flux = 1./(2.*mHOS);

   static constexpr double color_factor = squared_color_generator<Higgs, bar<UpTypeQuark>::type, UpTypeQuark>();

   if(!boost::range::equal(indexOut1, indexOut2)) {
      const double muqOS1 = context.physical_mass<UpTypeQuark>(indexOut1);
      const double muqOS2 = context.physical_mass<UpTypeQuark>(indexOut2);
      if (mHOS > muqOS1 + muqOS2 && (!is_zero(HQQbarVertexDR.left()) || !is_zero(HQQbarVertexDR.right()))) {
         const auto xOS1 = Sqr(muqOS1/mHOS);
         const auto xOS2 = Sqr(muqOS2/mHOS);
         const double phase_space = 1./(8.*Pi) * std::sqrt(KallenLambda(1., xOS1, xOS2));
         return flux * phase_space * color_factor * amplitude_squared<Higgs, bar<UpTypeQuark>::type, UpTypeQuark>(context, indexIn, indexOut1, indexOut2);
      }
      return 0.;
   }

   const double muqDR = context.mass<UpTypeQuark>(indexOut1);
   const double muqOS = context.physical_mass<UpTypeQuark>(indexOut1);
   if(is_zero(muqDR) || is_zero(muqOS)) {
      throw std::runtime_error(
         create_process_string<Higgs,bar<UpTypeQuark>::type, UpTypeQuark>(indexIn, indexOut1, indexOut2)
            + ": Up-type quark cannot be massless. Aborting."
      );
   }
   const auto xOS = Sqr(muqOS/mHOS);
   const auto xDR = Sqr(muqDR/mHOS);

   const auto HQQbarVertexDR_S = 0.5*(HQQbarVertexDR.left() + HQQbarVertexDR.right());
   const auto HQQbarVertexDR_P = 0.5*(HQQbarVertexDR.right() - HQQbarVertexDR.left());

   double result_S = 0;
   double result_P = 0;

   if (4.*std::max(xDR, xOS) > 1.) {
   }
   else {
      const auto betaOS = std::sqrt(1.-4.*xOS);
      const auto betaDR = std::sqrt(1.-4.*xDR);

      const double phase_spaceDR = 1./(8.*Pi) * std::sqrt(KallenLambda(1., xDR, xDR));
      const double phase_spaceOS = 1./(8.*Pi) * std::sqrt(KallenLambda(1., xOS, xOS));

      double amp2DR_S = Sqr(mHOS) * Sqr(betaDR) *
                        2*std::norm(HQQbarVertexDR_S);
      double amp2OS_S = Sqr(mHOS) * Sqr(betaOS) *
                        2*std::norm(HQQbarVertexDR_S) * Sqr(muqOS / muqDR);

      double amp2DR_P = 0;
      double amp2OS_P = 0;
      if (info::is_CP_violating_Higgs_sector) {
         amp2DR_P = Sqr(mHOS) *
                    2*std::norm(HQQbarVertexDR_P);
         amp2OS_P = Sqr(mHOS) *
                    2*std::norm(HQQbarVertexDR_P) * Sqr(muqOS / muqDR);
      }

      const int Nf = number_of_active_flavours(qedqcd, mHOS);
      if (Nf < 5) {
         problems.add_warning(
            create_process_string<Higgs,bar<UpTypeQuark>::type, UpTypeQuark>(indexIn, indexOut1, indexOut2)
               + ": Cannot determine the number of active quark flavours. Disabling higher-order corrections."
         );
      }

      if (static_cast<int>(flexibledecay_settings.get(FlexibleDecay_settings::include_higher_order_corrections)) > 0 && Nf >= 5)  {
            double alpha_s_red;
            double Y_conversion = 1.;
            switch (Nf) {
               case 5: {
                  auto qedqcd_ = qedqcd;
                  qedqcd_.to(mHOS);
                  alpha_s_red = qedqcd_.displayAlpha(softsusy::ALPHAS)/Pi;
                  Y_conversion = Sqr(qedqcd_.displayUpQuarkRunningMass(indexOut1.at(0))/muqDR);
                  break;
               }
               case 6:
                  alpha_s_red = get_alphas(context)/Pi;
                  break;
               default:
                  throw std::runtime_error(
                     create_process_string<Higgs,bar<UpTypeQuark>::type, UpTypeQuark>(indexIn, indexOut1, indexOut2)
                        + ": Cannot determine the number of active quark flavours"
                  );
            }
            double deltaqq_QCD_DR_S = calc_Deltaqq(alpha_s_red, Nf, flexibledecay_settings);
            double deltaqq_QCD_DR_P = deltaqq_QCD_DR_S;

            // 1L QED correction - eq. 17 in FD paper
            const double alpha_red = get_alpha(context)/Pi;
            const double deltaqq_QED_DR = 17./4.*Sqr(UpTypeQuark::electricCharge)*alpha_red;

            deltaqq_QCD_DR_S +=
               2.*(1. - 10.*xDR)/(1-4.*xDR)*(4./3. - std::log(xDR))*alpha_s_red +
               4./3.*alpha_s_red*calc_DeltaH(betaDR);

            const double deltaqq_QCD_OS_S =
               4./3. * alpha_s_red * calc_DeltaH(betaOS);

            const double deltaqq_QED_OS_S =
               alpha_red * Sqr(UpTypeQuark::electricCharge) * calc_DeltaH(betaOS);

            double deltaqq_QCD_OS_P = 0.;
            double deltaqq_QED_OS_P = 0.;
            // don't waste time computing it in models without CPV
            if (info::is_CP_violating_Higgs_sector) {
               deltaqq_QCD_DR_P +=
                  2.*(1. - 6.*xDR)/(1-4.*xDR)*(4./3. - std::log(xDR))*alpha_s_red +
                  4./3.*alpha_s_red*calc_DeltaAh(betaDR);

               deltaqq_QCD_OS_P =
                  4./3. * alpha_s_red * calc_DeltaAh(betaOS);

               deltaqq_QED_OS_P =
                  alpha_red * Sqr(UpTypeQuark::electricCharge) * calc_DeltaAh(betaOS);
            }

            double deltaPhi2_S = 0.;
            double deltaPhi2_P = 0.;
            double deltaqq_QCDxQED_DR = 0.;
            if (static_cast<int>(flexibledecay_settings.get(FlexibleDecay_settings::include_higher_order_corrections)) > 1) {
               deltaqq_QCDxQED_DR =
                  deltaqq_QCDxQED*Sqr(UpTypeQuark::electricCharge)*alpha_red*alpha_s_red;
               if ((indexOut1.at(0) < 2 || indexOut2.at(0) < 2)) {
                  const double mtpole = qedqcd.displayPoleMt();
                  const double lt = std::log(Sqr(mHOS/mtpole));
                  const double lq = std::log(xDR);
                  // eq. 28 of hep-ph/9505358
                  const auto Httindices = concatenate(std::array<int, 1> {2}, std::array<int, 1> {2}, indexIn);
                  const auto Httbar = Vertex<bar<UpTypeQuark>::type, UpTypeQuark, Higgs>::evaluate(Httindices, context);
                  const auto gbHoVEV = HQQbarVertexDR_S/muqDR;
                  if (!is_zero(gbHoVEV)) {
                     const auto Httbar_S = 0.5*(Httbar.left() + Httbar.right());
                     const auto gtHoVEV = Httbar_S/context.mass<UpTypeQuark>({2});
                     deltaPhi2_S = Sqr(alpha_s_red) * std::real(gtHoVEV/gbHoVEV) * (8/3. - Sqr(Pi/3.) - 2.0/3.0*lt + 1.0/9.0*Sqr(lq));
                  }
                  if (info::is_CP_violating_Higgs_sector) {
                     const auto CSuu = HQQbarVertexDR_S/muqDR;
                     if (!is_zero(CSuu)) {
                        const auto Httbar_P = 0.5*(Httbar.right() - Httbar.left());
                        const auto CStu = Httbar_P/context.mass<Fu>({2});
                        deltaPhi2_P = Sqr(alpha_s_red) * std::real(CStu/CSuu) * (23/6. - lt + 1.0/6.0*Sqr(lq));
                     }
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

      result_S = (1-4.*xOS)*result_DR_S + 4*xOS*result_OS_S;
      result_P = (1-4.*xOS)*result_DR_P + 4*xOS*result_OS_P;
   }

#ifdef ENABLE_HIGGSTOOLS
   if (static_cast<bool>(flexibledecay_settings.get(FlexibleDecay_settings::call_higgstools))) {
      higgstools_input.add_coupling(field_as_string<Higgs>(indexIn),
         {
            boost::hana::unpack(bar<UpTypeQuark>::pdgids, _to_array<bar<UpTypeQuark>::numberOfGenerations>).at(indexOut1.at(0)),
            boost::hana::unpack(UpTypeQuark::pdgids, _to_array<UpTypeQuark::numberOfGenerations>).at(indexOut2.at(0))
         },
         4.*std::max(xDR, xOS) > 1. ? std::real(HQQbarVertexDR_S) + std::imag(HQQbarVertexDR_P)*1i : std::sqrt(result_S) + std::sqrt(result_P)*1i);
   }
#endif

   return result_S + result_P;
}
