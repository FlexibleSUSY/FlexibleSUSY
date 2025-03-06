template <>
double CLASSNAME::get_partial_width<Higgs, bar<ChargedLepton>::type, ChargedLepton>(
   const context_base& context,
   typename field_indices<Higgs>::type const& indexIn,
   typename field_indices<bar<ChargedLepton>::type>::type const& indexOut1,
   typename field_indices<ChargedLepton>::type const& indexOut2)
{

   const auto indices = concatenate(indexOut1, indexOut2, indexIn);
   const auto HLLbarVertexDR = Vertex<bar<ChargedLepton>::type, ChargedLepton, Higgs>::evaluate(indices, context);
   const auto HLLbarVertexDR_S = 0.5*(HLLbarVertexDR.left() + HLLbarVertexDR.right());
   const auto HLLbarVertexDR_P = 0.5*(HLLbarVertexDR.right() - HLLbarVertexDR.left());

   if (is_zero(HLLbarVertexDR_S) && is_zero(HLLbarVertexDR_S)) {
      return 0.;
   }

   const double mHOS = context.physical_mass<Higgs>(indexIn);
   const double mL1OS = context.physical_mass<bar<ChargedLepton>::type>(indexOut1);
   const double mL2OS = context.physical_mass<ChargedLepton>(indexOut2);

   // phase space without symmetry factor
   const auto xOS1 = Sqr(mL1OS/mHOS);
   const auto xOS2 = Sqr(mL2OS/mHOS);
   const auto betaOS2 = 1.0 - Sqr((mL1OS + mL2OS)/mHOS);
   const double ps = 1./(8.*Pi)*std::sqrt(KallenLambda(1., xOS1, xOS2));

   const double amp2_S = Sqr(mHOS) * betaOS2 * 2*std::norm(HLLbarVertexDR_S);

   const double amp2_P =
      info::is_CP_violating_Higgs_sector ? 2. * (Sqr(mHOS) - Sqr(mL1OS -mL2OS)) * std::norm(HLLbarVertexDR_P) : 0.;

   // flux * phase space factor * symmetry factor * |matrix element|^2
   double res_S = 0.5 * ps * amp2_S/mHOS;
   double res_P = 0.5 * ps * amp2_P/mHOS;

   const auto amp_tree =
      calculate_amplitude_tree<Higgs, bar<ChargedLepton>::type, ChargedLepton>(context, indexIn, indexOut1, indexOut2);

   double res = res_S + res_P;

   // higher order SM corrections
   if (flexibledecay_settings.get(FlexibleDecay_settings::include_higher_order_corrections)) {
      // 1-loop QED corrections
      res_S *= 1. + get_alpha(context)/Pi*17./4.;
      res_P *= 1. + get_alpha(context)/Pi*17./4.;
   }

   if (true /* higher order BSM corrections */) {

      namespace smns = standard_model_cxx_diagrams;

      smns::context_base sm_context {sm};

   //**************************************************************************
   // Self-Energies for the different particles in the decoupling definition
   // that is the difference between BSM-and SM self-energy
   //**************************************************************************

   //**************************************************************************
   // Lepton Self-Energies

   const auto se_PL_lepton_BSM = self_energy_1loop_PL<Fe>(context, mL1OS).coeff(indexOut1.at(0), indexOut2.at(0)) - smns::self_energy_1loop_PL<smns::fields::Fe>(sm_context, mL1OS).coeff(indexOut1.at(0), indexOut2.at(0));
   const auto se_PR_lepton_BSM = self_energy_1loop_PR<Fe>(context, mL1OS).coeff(indexOut1.at(0), indexOut2.at(0)) - smns::self_energy_1loop_PR<smns::fields::Fe>(sm_context, mL1OS).coeff(indexOut1.at(0), indexOut2.at(0));
   const auto se_1_lepton_BSM = self_energy_1loop_1<Fe>(context, mL1OS).coeff(indexOut1.at(0), indexOut2.at(0)) - smns::self_energy_1loop_1<smns::fields::Fe>(sm_context, mL1OS).coeff(indexOut1.at(0), indexOut2.at(0));

   // End Lepton Self-Energies
   //**************************************************************************

   //**************************************************************************
   // Higgs Self-Energy

   const auto se_Higgs_BSM = self_energy_1loop<hh>(context, mHOS) - smns::self_energy_1loop<smns::fields::hh>(sm_context, mHOS);

   // End Higgs Self-Energy

   //**************************************************************************
   // W and Z Self-Energy

   const double mWOS = context.physical_mass<WpBoson>({});
   const double mZOS = context.physical_mass<ZBoson>({});

   const auto se_W_BSM = self_energy_1loop<VWp>(context, mWOS) - smns::self_energy_1loop<smns::fields::VWp>(sm_context, mWOS);
   const auto se_Z_BSM = self_energy_1loop<VZ>(context, mZOS) - smns::self_energy_1loop<smns::fields::VZ>(sm_context, mZOS);

   // End Z and W Self-Energy
   //**************************************************************************


   //**************************************************************************
   // renormalization constants for the Higgs-Lepton-Lepton counter term
   // the renromalization constants are evaluated with then on-shell masses
   //**************************************************************************  

   //**************************************************************************    
   // Lepton Renormalization constants in the Decoupling Renormalization Scheme
   // if there is a problem the first thing to check would be the definition of the scalar part of the fermion self-energy

   // for real couplings we can introduce the sum of left-and right handed vector part in the BSM model
   const auto sum_LR_BSM = self_energy_1loop_PL<Fe>(context, mL1OS).coeff(indexOut1.at(0), indexOut2.at(0)) + self_energy_1loop_PR<Fe>(context, mL1OS).coeff(indexOut1.at(0), indexOut1.at(0));

   // derivative
   const auto sum_LR_BSM_deriv = self_energy_1loop_PL_deriv_p2<Fe>(context, mL1OS).coeff(indexOut1.at(0), indexOut2.at(0)) + self_energy_1loop_PR_deriv_p2<Fe>(context, mL1OS).coeff(indexOut1.at(0), indexOut2.at(0));

   // the same helper functions for the pure Standard model
   const auto sum_LR_SM = smns::self_energy_1loop_PL<smns::fields::Fe>(sm_context, mL1OS).coeff(indexOut1.at(0), indexOut1.at(0)) + smns::self_energy_1loop_PR<smns::fields::Fe>(sm_context, mL1OS).coeff(indexOut1.at(0), indexOut1.at(0));

   // derivative
   const auto sum_LR_SM_deriv = smns::self_energy_1loop_PL_deriv_p2<smns::fields::Fe>(sm_context, mL1OS).coeff(indexOut1.at(0), indexOut1.at(0)) + smns::self_energy_1loop_PR_deriv_p2<smns::fields::Fe>(sm_context, mL1OS).coeff(indexOut1.at(0), indexOut1.at(0));

   // left field renormalization constant in the decoupling scheme
   const auto delta_ZL_lepton = -(self_energy_1loop_PL<Fe>(context, mL1OS).coeff(indexOut1.at(0), indexOut1.at(0)) - smns::self_energy_1loop_PL<smns::fields::Fe>(sm_context, mL1OS).coeff(indexOut1.at(0), indexOut1.at(0))) - mL1OS * mL1OS * (sum_LR_BSM_deriv - sum_LR_SM_deriv) - mL1OS * mL1OS * (self_energy_1loop_1_deriv_p2<Fe>(context, mL1OS).coeff(indexOut1.at(0), indexOut1.at(0)) - smns::self_energy_1loop_1_deriv_p2<smns::fields::Fe>(sm_context, mL1OS).coeff(indexOut1.at(0), indexOut1.at(0)));

   // right field renormalization constant in the decoupling scheme
   const auto delta_ZR_lepton = -(self_energy_1loop_PR<Fe>(context, mL1OS).coeff(indexOut1.at(0), indexOut1.at(0)) - smns::self_energy_1loop_PR<smns::fields::Fe>(sm_context, mL1OS).coeff(indexOut1.at(0), indexOut1.at(0))) - mL1OS * mL1OS * (sum_LR_BSM_deriv - sum_LR_SM_deriv) - mL1OS * mL1OS * (self_energy_1loop_1_deriv_p2<Fe>(context, mL1OS).coeff(indexOut1.at(0), indexOut1.at(0)) - smns::self_energy_1loop_1_deriv_p2<smns::fields::Fe>(sm_context, mL1OS).coeff(indexOut1.at(0), indexOut1.at(0)));

   // mass renormalization constant in the decoupling scheme
   const auto delta_Zm_lepton = sum_LR_BSM - sum_LR_SM + (self_energy_1loop_1<Fe>(context, mL1OS).coeff(indexOut1.at(0), indexOut1.at(0)) - smns::self_energy_1loop_1<smns::fields::Fe>(sm_context, mL1OS).coeff(indexOut1.at(0), indexOut1.at(0)));

   // End Fermion Renormalization constants
   //**************************************************************************


   //**************************************************************************
   // Higgs Renormalization constants in the Decoupling Renormalization Scheme
   // Here we only list the SM-like Higgs renormalization constants

   // Higgs Field renormalization constant in the decoupling scheme
   const auto delta_Z_higgs = self_energy_1loop_deriv_p2<hh>(context, mHOS) - smns::self_energy_1loop_deriv_p2<smns::fields::hh>(sm_context, mHOS);

   // Higgs Mass renormalization constant in the decoupling scheme
   const auto delta_Zm_higgs = self_energy_1loop<hh>(context, mHOS)/Sqr(mHOS) - smns::self_energy_1loop<smns::fields::hh>(sm_context, mHOS)/Sqr(mHOS);

   // End Higgs Renormalization constants
   //**************************************************************************

   //**************************************************************************
   // Gauge Field Renormalization constants in the Decoupling Scheme


   // W Field renormalization constant in the decoupling scheme
   const auto delta_Z_W = -(self_energy_1loop_deriv_p2<VWp>(context, mWOS) - smns::self_energy_1loop_deriv_p2<smns::fields::VWp>(sm_context, mWOS));

   // W mass renormalization constant in the decoupling scheme
   const auto delta_Zm_W = self_energy_1loop<VWp>(context, mWOS)/Sqr(mWOS) - smns::self_energy_1loop<smns::fields::VWp>(sm_context, mWOS)/Sqr(mWOS);

   // Z Field renormalization constant in the decoupling scheme
   const auto delta_Z_Z = -(self_energy_1loop_deriv_p2<VZ>(context, mZOS) - smns::self_energy_1loop_deriv_p2<smns::fields::VZ>(sm_context, mZOS));

   // Z mass renormalization constant in the decoupling scheme
   const auto delta_Zm_Z = self_energy_1loop<VZ>(context, mZOS)/Sqr(mZOS) - smns::self_energy_1loop<smns::fields::VZ>(sm_context, mZOS)/Sqr(mZOS);

   // photon field renormalization constant in the decoupling scheme
   const auto delta_Z_A = -(self_energy_1loop_deriv_p2<VP>(context, 0) - smns::self_energy_1loop<smns::fields::hh>(sm_context, 0));


   // there is still the A-Z mixing fieldrenormalization constants missing


   // End Gauge Field Renormalization constants
   //**************************************************************************

   //**************************************************************************
   // Additional Renormalization constants that can be expressed through the Renormalization constants above

   // charge renormalization constant in the decoupling scheme
   const auto delta_e = -0.5*delta_Z_A; //- tan(theta_W) / 2 * delta_Z_ZA

   // it is still necessary to decine the W-and Z-mass as well as the vev
   const double vev = 246;
   //vev renormalization constant in the decoupling scheme
   const auto delta_Z_vev = 0; //(1 - 1 / (2 * sin(theta_W)) * delta_Zm_W) + cot(theta_W) / 2 * delta_Zm_Z - delta_e

   // End Additional Renormalization constants
   //**************************************************************************

   //**************************************************************************
   // Counter Terms for the Left-and Right Form-Factor

   // counter term for the left form factor in the decoupling scheme
   const auto counter_term_FL = 0.5/vev * (mL1OS * delta_ZL_lepton + mL2OS * delta_ZR_lepton) - mL1OS / vev * delta_Zm_lepton - mL1OS * (delta_Z_higgs - delta_Z_vev);

   // counter term for the right form factor in the decoupling scheme
   const auto counter_term_FR = 0.5/vev * (mL1OS * delta_ZR_lepton + mL2OS * delta_ZL_lepton) - mL1OS / vev * delta_Zm_lepton - mL1OS * (delta_Z_higgs - delta_Z_vev);

   //End Left-and Right Form_Factors

      const auto amp_1l = calculate_amplitude_1l<Higgs, bar<ChargedLepton>::type, ChargedLepton>(context, indexIn, indexOut1, indexOut2);

      if (indexOut1.at(0) == 1 && indexOut2.at(0) == 1) {

         std::cout << "\n BSM Lepton self-energy contributions left: " << se_PL_lepton_BSM << '\n';
         std::cout << " BSM Lepton self-energy contributions right: " << se_PR_lepton_BSM << '\n';
         std::cout << " BSM Lepton self-energy contributions scalar: " << se_1_lepton_BSM << '\n';

         std::cout << "\n BSM self-energy contributions Higgs: " << se_Higgs_BSM << '\n';
         std::cout << "\n BSM self-energy contributions W: " << se_W_BSM << '\n';
         std::cout << " BSM self-energy contributions Z: " << se_Z_BSM << '\n';

         std::cout << "\n Lepton Field Renormalization Left: " << delta_ZL_lepton << '\n';
         std::cout << " Lepton Field Renormalization Right: " << delta_ZR_lepton << '\n';
         std::cout << " Lepton Mass Renormalization: " << delta_Zm_lepton << '\n';

         std::cout << "\n Counter Term left form-factor: " << counter_term_FL << '\n';
         std::cout << " Counter Term right form-factor: " << counter_term_FR << '\n';

         std::cout << "\n form factor left: " << amp_1l.form_factor_left << '\n';
         std::cout << " form factor right: " << amp_1l.form_factor_right << '\n';
         std::cout << "\n 1l* conj(tree): " << amplitude_interference(amp_1l, amp_tree) << '\n';
         std::cout << " |tree|^2: " << amp_tree.square() << '\n';
         std::cout << " |1l|^2: " << amp_1l.square() << '\n';
      }

      res += 2.*std::real(amplitude_interference(amp_tree, amp_1l)) /* + ren const */;
   }

   if (flexibledecay_settings.get(FlexibleDecay_settings::call_higgstools) != 0 ||
       flexibledecay_settings.get(FlexibleDecay_settings::call_lilith) != 0 ||
       flexibledecay_settings.get(FlexibleDecay_settings::calculate_normalized_effc) != 0
   ) {
      neutral_higgs_effc.add_coupling(
         field_as_string<Higgs>(indexIn),
            {
               boost::hana::unpack(bar<ChargedLepton>::pdgids, _to_array<bar<ChargedLepton>::numberOfGenerations>).at(indexOut1.at(0)),
               boost::hana::unpack(ChargedLepton::pdgids, _to_array<ChargedLepton::numberOfGenerations>).at(indexOut2.at(0))
            },
         std::pair<std::string, std::complex<double>> {
            field_as_string<Higgs>(indexIn) + "-" + field_as_string<bar<ChargedLepton>::type>(indexOut1) + "-" + field_as_string<ChargedLepton>(indexOut2),
            std::sqrt(res_S) + 1i*std::sqrt(res_P)
         }
      );
   }

   return res;
}
