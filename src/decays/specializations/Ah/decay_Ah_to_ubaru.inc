template <>
double CLASSNAME::get_partial_width<PseudoscalarHiggs, bar<UpTypeQuark>::type, UpTypeQuark>(
   const context_base& context,
   typename field_indices<PseudoscalarHiggs>::type const& indexIn,
   typename field_indices<UpTypeQuark>::type const& indexOut1,
   typename field_indices<UpTypeQuark>::type const& indexOut2)
{
   if (indexIn.at(0) < info::number_of_neutral_goldstones) {
      throw OutOfBoundsError("Error in " + create_process_string<PseudoscalarHiggs,bar<UpTypeQuark>::type, UpTypeQuark>(indexIn, indexOut1, indexOut2) + " decay. Decaying particle is a Goldstone.");
   }

   // get AhQQbar vertex
   // we don't use amplitude_squared here because we need both this vertex
   // both with running and pole masses
   const auto indices = concatenate(indexOut1, indexOut2, indexIn);
   const auto AhQQbarVertexDR = Vertex<bar<UpTypeQuark>::type, UpTypeQuark, PseudoscalarHiggs>::evaluate(indices, context);

   const double mAhOS = context.physical_mass<PseudoscalarHiggs>(indexIn);
   const double flux = 1./(2.*mAhOS);

   static constexpr double color_factor = squared_color_generator<PseudoscalarHiggs, bar<UpTypeQuark>::type, UpTypeQuark>();

   // flavour-violating decay
   if(!boost::range::equal(indexOut1, indexOut2)) {
      const double muqOS1 = context.physical_mass<UpTypeQuark>(indexOut1);
      const double muqOS2 = context.physical_mass<UpTypeQuark>(indexOut2);
      if (mAhOS > muqOS1 + muqOS2 && (!is_zero(AhQQbarVertexDR.left()) || !is_zero(AhQQbarVertexDR.right()))) {
         const auto xOS1 = Sqr(muqOS1/mAhOS);
         const auto xOS2 = Sqr(muqOS2/mAhOS);
         const double phase_space = 1./(8.*Pi) * std::sqrt(KallenLambda(1., xOS1, xOS2));
         return flux * phase_space * color_factor * amplitude_squared<PseudoscalarHiggs, bar<UpTypeQuark>::type, UpTypeQuark>(context, indexIn, indexOut1, indexOut2);
      }
      return 0.;
   }

   const double muqDR = context.mass<UpTypeQuark>(indexOut1);
   const double muqOS = context.physical_mass<UpTypeQuark>(indexOut1);
   if(is_zero(muqDR) || is_zero(muqOS)) {
      throw std::runtime_error(
         create_process_string<PseudoscalarHiggs,bar<UpTypeQuark>::type, UpTypeQuark>(indexIn, indexOut1, indexOut2)
            + ": Up-type quark cannot be massless. Aborting."
      );
   }
   const auto xOS = Sqr(muqOS/mAhOS);
   const auto xDR = Sqr(muqDR/mAhOS);

   double result = 0.;

   const auto AhQQbarVertexDR_P = 0.5*(AhQQbarVertexDR.right() - AhQQbarVertexDR.left());

   const double phase_spaceOS = 1./(8.*Pi) * std::sqrt(KallenLambda(1., xOS, xOS));

   // off-shell decay
   if (4.*std::max(xDR, xOS) > 1.) {
   }
   else {
      const auto betaOS = std::sqrt(1.-4.*xOS);
      const auto betaDR = std::sqrt(1.-4.*xDR);

      const double phase_spaceDR = 1./(8.*Pi) * std::sqrt(KallenLambda(1., xDR, xDR));

      double amp2DR_P = 0;
      double amp2OS_P = 0;
      amp2DR_P = Sqr(mAhOS) *
                 2*std::norm(AhQQbarVertexDR_P);
      amp2OS_P = Sqr(mAhOS) *
                 2*std::norm(AhQQbarVertexDR_P) * Sqr(muqOS / muqDR);

      const int Nf = number_of_active_flavours(qedqcd, mAhOS);
      if (Nf < 5) {
         problems.add_warning(
            create_process_string<PseudoscalarHiggs,bar<UpTypeQuark>::type, UpTypeQuark>(indexIn, indexOut1, indexOut2)
               + ": Cannot determine the number of active quark flavours. Disabling higher-order corrections."
         );
      }

      if (flexibledecay_settings.get_decay_corrections().sm > 0 && Nf >= 5) {
            double alpha_s_red;
            double Y_conversion = 1.;
            switch (Nf) {
               case 5: {
                  auto qedqcd_ = qedqcd;
                  qedqcd_.to(mAhOS);
                  alpha_s_red = qedqcd_.displayAlpha(softsusy::ALPHAS)/Pi;
                  Y_conversion = Sqr(qedqcd_.displayUpQuarkRunningMass(indexOut1.at(0))/muqDR);
                  break;
               }
               case 6:
                  alpha_s_red = get_alphas(context)/Pi;
                  break;
               default:
                  throw std::runtime_error(
                     create_process_string<PseudoscalarHiggs,bar<UpTypeQuark>::type, UpTypeQuark>(indexIn, indexOut1, indexOut2)
                        + ": Cannot determine the number of active quark flavours"
                  );
            }
            double deltaqq_QCD_DR_P =
                  calc_Deltaqq(alpha_s_red, Nf, flexibledecay_settings)
                  + 2.*(1. - 6.*xDR)/(1-4.*xDR)*(4./3. - std::log(xDR))*alpha_s_red
                  + 4./3.*alpha_s_red*calc_DeltaAh(betaDR);

            // 1L QED correction - eq. 17 in FD paper
            const double alpha_red = get_alpha(context)/Pi;
            const double deltaqq_QED_DR = 17./4.*Sqr(UpTypeQuark::electricCharge)*alpha_red;

            const double deltaqq_QCD_OS_P =
                  4./3. * alpha_s_red * calc_DeltaAh(betaOS);

            const double deltaqq_QED_OS_P =
                  alpha_red * Sqr(UpTypeQuark::electricCharge) * calc_DeltaAh(betaOS);

            double deltaqq_QCDxQED_DR = 0.;
            double deltaPhi2_P = 0.;
            if (static_cast<int>(flexibledecay_settings.get(FlexibleDecay_settings::include_higher_order_corrections)) > 1) {
               deltaqq_QCDxQED_DR =
                  deltaqq_QCDxQED*Sqr(UpTypeQuark::electricCharge)*alpha_red*alpha_s_red;
               if ((indexOut1.at(0) < 2 || indexOut2.at(0) < 2)) {
                  const double mtpole = qedqcd.displayPoleMt();
                  const double lt = std::log(Sqr(mAhOS/mtpole));
                  const double lq = std::log(xDR);
                  // eq. 28 of hep-ph/9505358
                  const auto Ahttindices = concatenate(std::array<int, 1> {2}, std::array<int, 1> {2}, indexIn);
                  const auto Ahttbar = Vertex<bar<UpTypeQuark>::type, UpTypeQuark, PseudoscalarHiggs>::evaluate(Ahttindices, context);
                  const auto CSuu = AhQQbarVertexDR_P/context.mass<UpTypeQuark>(indexOut1);
                  if (!is_zero(CSuu)) {
                     const auto Ahttbar_P = 0.5*(Ahttbar.right() - Ahttbar.left());
                     const auto CStu = Ahttbar_P/context.mass<Fu>({2});
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

      result = (1-4.*xOS)*result_DR + 4*xOS*result_OS;
   }

   if (flexibledecay_settings.get(FlexibleDecay_settings::call_higgstools) != 0 ||
       flexibledecay_settings.get(FlexibleDecay_settings::call_lilith) != 0 ||
       flexibledecay_settings.get(FlexibleDecay_settings::calculate_normalized_effc) != 0
   ) {
      const double remove_normalization = flux*phase_spaceOS*2*Sqr(mAhOS);
      neutral_higgs_effc.add_coupling(field_as_string<PseudoscalarHiggs>(indexIn),
         {
            boost::hana::unpack(bar<UpTypeQuark>::pdgids, _to_array<bar<UpTypeQuark>::numberOfGenerations>).at(indexOut1.at(0)),
            boost::hana::unpack(UpTypeQuark::pdgids, _to_array<UpTypeQuark::numberOfGenerations>).at(indexOut2.at(0))
         },
         std::pair<std::string, std::complex<double>> {
            field_as_string<PseudoscalarHiggs>(indexIn) + "-" + field_as_string<bar<UpTypeQuark>::type>(indexOut1) + "-" + field_as_string<UpTypeQuark>(indexOut2),
            // if the decay is not kinematically allowed return the tree-level coupling,
            // otherwise a sqrt of a partial width
            4.*std::max(xDR, xOS) > 1. ? std::imag(AhQQbarVertexDR_P)*1i : std::sqrt(result/remove_normalization)*1i
         },
         result
         );
   }

   return result;
}
