template <>
double CLASSNAME::get_partial_width<PseudoscalarHiggs, bar<DownTypeQuark>::type, DownTypeQuark>(
   const context_base& context,
   typename field_indices<PseudoscalarHiggs>::type const& indexIn,
   typename field_indices<bar<DownTypeQuark>::type>::type const& indexOut1,
   typename field_indices<DownTypeQuark>::type const& indexOut2
   )
{
   if (indexIn.at(0) < info::number_of_neutral_goldstones) {
      throw OutOfBoundsError("Error in " + create_process_string<PseudoscalarHiggs,bar<DownTypeQuark>::type, DownTypeQuark>(indexIn, indexOut1, indexOut2) + " decay. Decaying particle is a Goldstone.");
   }

   // get HBBbar vertex
   // we don't use amplitude_squared here because we need both this vertex
   // both with running and pole masses
   const auto indices = concatenate(indexOut1, indexOut2, indexIn);
   const auto HBBbarVertexDR = Vertex<bar<DownTypeQuark>::type, DownTypeQuark, PseudoscalarHiggs>::evaluate(indices, context);

   const double mAOS = context.physical_mass<PseudoscalarHiggs>(indexIn);
   const double flux = 1./(2.*mAOS);

   static constexpr double color_factor = noexcept(squared_color_generator<PseudoscalarHiggs, bar<DownTypeQuark>::type, DownTypeQuark>());

   if(!boost::range::equal(indexOut1, indexOut2)) {
      if (!is_zero(HBBbarVertexDR.left()) || !is_zero(HBBbarVertexDR.right())) {
         const double mdqOS1 = context.physical_mass<DownTypeQuark>(indexOut1);
         const double mdqOS2 = context.physical_mass<DownTypeQuark>(indexOut2);
         const auto xOS1 = Sqr(mdqOS1/mAOS);
         const auto xOS2 = Sqr(mdqOS2/mAOS);
         const double phase_space = 1./(8.*Pi) * std::sqrt(KallenLambda(1., xOS1, xOS2));
         return flux * phase_space * color_factor * amplitude_squared<PseudoscalarHiggs, bar<DownTypeQuark>::type, DownTypeQuark>(context, indexIn, indexOut1, indexOut2);
      }
      return 0.;
   }

   const double mdqDR = context.mass<DownTypeQuark>(indexOut1);
   const double mdqOS = context.physical_mass<DownTypeQuark>(indexOut1);
   if(is_zero(mdqDR) || is_zero(mdqOS)) {
      throw std::runtime_error(create_process_string<PseudoscalarHiggs,bar<DownTypeQuark>::type, DownTypeQuark>(indexIn, indexOut1, indexOut2) + ": Down-type quark cannot be massless. Aborting.");
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

   const int Nf = number_of_active_flavours(qedqcd, mAOS);
   if (Nf < 5) {
      problems.add_warning(create_process_string<PseudoscalarHiggs,bar<DownTypeQuark>::type, DownTypeQuark>(indexIn, indexOut1, indexOut2) + ": Cannot determine the number of active quark flavours. Disabling higher-order corrections.");
   }

   if (static_cast<int>(flexibledecay_settings.get(FlexibleDecay_settings::include_higher_order_corrections)) > 0 && Nf >= 5) {
         double alpha_s_red;
         double Y_conversion = 1.;
         switch (Nf) {
            case 5: {
               auto qedqcd_ = qedqcd;
               qedqcd_.to(mAOS);
               alpha_s_red = qedqcd_.displayAlpha(softsusy::ALPHAS)/Pi;
               Y_conversion = Sqr(qedqcd_.displayDownQuarkRunningMass(indexOut1.at(0))/mdqDR);
               break;
            }
            case 6:
               alpha_s_red = get_alphas(context)/Pi;
               break;
         }

         double deltaqq_QCD_DR_P = calc_Deltaqq(alpha_s_red, Nf, flexibledecay_settings);

         // 1L QED correction - eq. 17 in FD paper
         const double alpha_red = get_alpha(context)/Pi;
         const double deltaqq_QED_DR = 17./4.*Sqr(DownTypeQuark::electricCharge)*alpha_red;

         deltaqq_QCD_DR_P +=
            2.*(1. - 6.*xDR)/(1-4.*xDR)*(4./3. - std::log(xDR))*alpha_s_red +
            4./3.*alpha_s_red*calc_DeltaAh(betaDR);

         const double deltaqq_QCD_OS_P =
               4./3. * alpha_s_red * calc_DeltaAh(betaOS);

         const double deltaqq_QED_OS_P =
               alpha_red * Sqr(DownTypeQuark::electricCharge) * calc_DeltaAh(betaOS);

         const auto gbHoVEV_P = HBBbarVertexDR_P/context.mass<DownTypeQuark>(indexOut1);
         double deltaPhi2_P = 0.;
         double deltaqq_QCDxQED_DR = 0.;

         if (static_cast<int>(flexibledecay_settings.get(FlexibleDecay_settings::include_higher_order_corrections)) > 1) {
            deltaqq_QCDxQED_DR =
               deltaqq_QCDxQED*Sqr(DownTypeQuark::electricCharge)*alpha_red*alpha_s_red;
            if (!is_zero(gbHoVEV_P)) {
               const double mtpole = qedqcd.displayPoleMt();
               const double lt = std::log(Sqr(mAOS/mtpole));
               const double lq = std::log(xDR);
               const auto Ahttindices = concatenate(std::array<int, 1> {2}, std::array<int, 1> {2}, indexIn);
               const auto Ahttbar = Vertex<bar<UpTypeQuark>::type, UpTypeQuark, PseudoscalarHiggs>::evaluate(Ahttindices, context);
               const auto Ahttbar_P = 0.5*(Ahttbar.right() - Ahttbar.left());
               const auto gtHoVEV_P = Ahttbar_P/context.mass<UpTypeQuark>({2});
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

   const double result = (1-4.*xOS)*result_DR + 4*xOS*result_OS;

#ifdef ENABLE_HIGGSTOOLS
   if (static_cast<bool>(flexibledecay_settings.get(FlexibleDecay_settings::call_higgstools))) {
      higgstools_input.add_coupling(field_as_string<PseudoscalarHiggs>(indexIn),
         {
            boost::hana::unpack(bar<DownTypeQuark>::pdgids, _to_array<bar<DownTypeQuark>::numberOfGenerations>).at(indexOut1.at(0)),
            boost::hana::unpack(DownTypeQuark::pdgids, _to_array<DownTypeQuark::numberOfGenerations>).at(indexOut2.at(0))
         },
         1i*std::sqrt(result)
      );
   }
#endif

   return result;
}
