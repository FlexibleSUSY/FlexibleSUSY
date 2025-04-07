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

   const auto lepton_PL_all_contribution = self_energy_1loop_PL<ChargedLepton>(context, mL1OS).coeff(indexOut1.at(0), indexOut2.at(0));
   const auto lepton_PL_sm_contribution = smns::self_energy_1loop_PL<smns::fields::Fe>(sm_context, mL1OS).coeff(indexOut1.at(0), indexOut2.at(0));

   const auto lepton_PR_all_contribution = self_energy_1loop_PR<ChargedLepton>(context, mL1OS).coeff(indexOut1.at(0), indexOut2.at(0));
   const auto lepton_PR_sm_contribution = smns::self_energy_1loop_PR<smns::fields::Fe>(sm_context, mL1OS).coeff(indexOut1.at(0), indexOut2.at(0));

   const auto lepton_1_all_contribution = self_energy_1loop_1<ChargedLepton>(context, mL1OS).coeff(indexOut1.at(0), indexOut2.at(0));
   const auto lepton_1_sm_contribution = smns::self_energy_1loop_1<smns::fields::Fe>(sm_context, mL1OS).coeff(indexOut1.at(0), indexOut2.at(0));

   const auto lepton_PL_se_dec = lepton_PL_all_contribution - lepton_PL_sm_contribution;
   const auto lepton_PR_se_dec = lepton_PR_all_contribution - lepton_PR_sm_contribution;
   const auto lepton_1_se_dec = lepton_1_all_contribution - lepton_1_sm_contribution;

   //derivative

   const auto lepton_PL_all_contribution_deriv_p2 = self_energy_1loop_PL_deriv_p2<ChargedLepton>(context, mL1OS).coeff(indexOut1.at(0), indexOut2.at(0));
   const auto lepton_PL_sm_contribution_deriv_p2 = smns::self_energy_1loop_PL_deriv_p2<smns::fields::Fe>(sm_context, mL1OS).coeff(indexOut1.at(0), indexOut2.at(0));

   const auto lepton_PR_all_contribution_deriv_p2 = self_energy_1loop_PR_deriv_p2<ChargedLepton>(context, mL1OS).coeff(indexOut1.at(0), indexOut2.at(0));
   const auto lepton_PR_sm_contribution_deriv_p2 = smns::self_energy_1loop_PR_deriv_p2<smns::fields::Fe>(sm_context, mL1OS).coeff(indexOut1.at(0), indexOut2.at(0));

   const auto lepton_1_all_contribution_deriv_p2 = self_energy_1loop_1_deriv_p2<ChargedLepton>(context, mL1OS).coeff(indexOut1.at(0), indexOut2.at(0));
   const auto lepton_1_sm_contribution_deriv_p2 = smns::self_energy_1loop_1_deriv_p2<smns::fields::Fe>(sm_context, mL1OS).coeff(indexOut1.at(0), indexOut2.at(0));

   const auto lepton_PL_se_dec_deriv_p2 = lepton_PL_all_contribution_deriv_p2 - lepton_PL_sm_contribution_deriv_p2;
   const auto lepton_PR_se_dec_deriv_p2 = lepton_PR_all_contribution_deriv_p2 - lepton_PR_sm_contribution_deriv_p2;
   const auto lepton_1_se_dec_deriv_p2 = lepton_1_all_contribution_deriv_p2 - lepton_1_sm_contribution_deriv_p2;

   // End Lepton Self-Energies
   //**************************************************************************

   //**************************************************************************
   // Higgs Self-Energy

   const auto higgs_all_contribution = self_energy_1loop<Higgs>(context, mHOS);

   const auto higgs_sm_contribution = smns::self_energy_1loop<smns::fields::hh>(sm_context, mHOS);

   const auto higgs_se_dec = higgs_all_contribution - higgs_sm_contribution;

   //derivative

   const auto higgs_all_contribution_deriv_p2 = self_energy_1loop_deriv_p2<Higgs>(context, mHOS);

   const auto higgs_sm_contribution_deriv_p2 = smns::self_energy_1loop_deriv_p2<smns::fields::hh>(sm_context, mHOS);

   const auto higgs_se_dec_deriv_p2 = higgs_all_contribution_deriv_p2 - higgs_sm_contribution_deriv_p2;

   // End Higgs Self-Energy
   //**************************************************************************

   //**************************************************************************
   // Photon Self-Energy

   const auto P_all_contribution = self_energy_1loop<Photon>(context, 0);
   const auto P_sm_contribution = smns::self_energy_1loop<smns::fields::VP>(sm_context, 0);

   const auto P_se_dec = P_all_contribution - P_sm_contribution;

   //derivative
   const auto P_all_contribution_deriv_p2 = self_energy_1loop_deriv_p2<Photon>(context, 0);
   const auto P_sm_contribution_deriv_p2 = smns::self_energy_1loop_deriv_p2<smns::fields::VP>(sm_context, 0);

   const auto P_se_dec_deriv_p2 = P_all_contribution_deriv_p2 - P_sm_contribution_deriv_p2;

   // End Z and W Self-Energy
   //**************************************************************************

   //**************************************************************************
   // Z Self-Energy

   const double mZOS = context.physical_mass<ZBoson>({});

   const auto Z_all_contribution = self_energy_1loop<ZBoson>(context, mZOS);
   const auto Z_sm_contribution = smns::self_energy_1loop<smns::fields::VZ>(sm_context, mZOS);

   const auto Z_se_dec = Z_all_contribution - Z_sm_contribution;

   //derivative
   const auto Z_all_contribution_deriv_p2 = self_energy_1loop_deriv_p2<ZBoson>(context, mZOS);
   const auto Z_sm_contribution_deriv_p2 = smns::self_energy_1loop_deriv_p2<smns::fields::VZ>(sm_context, mZOS);

   const auto Z_se_dec_deriv_p2 = Z_all_contribution_deriv_p2 - Z_sm_contribution_deriv_p2;

   // End Z Self-Energy
   //**************************************************************************

   //**************************************************************************
   // W Self-Energy

   const double mWOS = context.physical_mass<WpBoson>({});

   const auto W_all_contribution = self_energy_1loop<WpBoson>(context, mWOS);
   const auto W_sm_contribution = smns::self_energy_1loop<smns::fields::VWp>(sm_context, mWOS);

   const auto W_se_dec = W_all_contribution - W_sm_contribution;


   //derivative
   const auto W_all_contribution_deriv_p2 = self_energy_1loop_deriv_p2<WpBoson>(context, mWOS);
   const auto W_sm_contribution_deriv_p2 = smns::self_energy_1loop_deriv_p2<smns::fields::VWp>(sm_context, mWOS);


   const auto W_se_dec_deriv_p2 = W_all_contribution_deriv_p2 - W_sm_contribution_deriv_p2;

   // End W Self-Energy
   //**************************************************************************


   //**************************************************************************
   // renormalization constants for the Higgs-Lepton-Lepton counter term
   // the renromalization constants are evaluated with then on-shell masses
   //**************************************************************************

   //**************************************************************************
   // Lepton Renormalization constants in the Decoupling Renormalization Scheme
   // if there is a problem the first thing to check would be the definition of the scalar part of the fermion self-energy

   // left field renormalization constant in the decoupling scheme
   const auto delta_ZL_lepton_dec = -2 * lepton_PL_se_dec - mL1OS * mL1OS * (2 * lepton_PL_se_dec_deriv_p2 + 2 * lepton_PR_se_dec_deriv_p2) - mL1OS * (lepton_1_se_dec_deriv_p2 + lepton_1_se_dec_deriv_p2);

   // right field renormalization constant in the decoupling scheme
   const auto delta_ZR_lepton_dec = - 2 * lepton_PR_se_dec - mL1OS * mL1OS * (2 * lepton_PL_se_dec_deriv_p2 + 2 * lepton_PR_se_dec_deriv_p2) - mL1OS * (lepton_1_se_dec_deriv_p2 + lepton_1_se_dec_deriv_p2);

   // mass renormalization constant in the decoupling scheme
   const auto delta_Zm_lepton_dec = 0.5 * (2 * lepton_PL_se_dec + 2 * lepton_PR_se_dec) + 0.5 / mL1OS * (lepton_1_se_dec + lepton_1_se_dec);

   // End Fermion Renormalization constants
   //**************************************************************************


   //**************************************************************************
   // Higgs Renormalization constants in the Decoupling Renormalization Scheme
   // Here we only list the SM-like Higgs renormalization constants

   // Higgs Field renormalization constant in the decoupling scheme
   const auto delta_Z_higgs_dec = -higgs_se_dec_deriv_p2;

   // Higgs Mass renormalization constant in the decoupling scheme
   const auto delta_Zm_higgs_dec = higgs_se_dec/Sqr(mHOS);

   // End Higgs Renormalization constants
   //**************************************************************************

   //**************************************************************************
   // Gauge Field Renormalization constants in the Decoupling Scheme


   // W Field renormalization constant in the decoupling scheme
   const auto delta_Z_W_dec = -W_se_dec_deriv_p2;

   // W mass renormalization constant in the decoupling scheme
   const auto delta_Zm_W_dec = W_se_dec/Sqr(mWOS);

   // Z Field renormalization constant in the decoupling scheme
   const auto delta_Z_Z_dec = -Z_se_dec_deriv_p2;

   // Z mass renormalization constant in the decoupling scheme
   const auto delta_Zm_Z_dec = Z_se_dec/Sqr(mZOS);

   // photon field renormalization constant in the decoupling scheme
   const auto delta_Z_A_dec = -P_se_dec_deriv_p2;


   // there is still the A-Z mixing fieldrenormalization constants missing


   // End Gauge Field Renormalization constants
   //**************************************************************************

   //**************************************************************************
   // Additional Renormalization constants that can be expressed through the Renormalization constants above

   const double vev = 248.9397360685197;
   const double theta_W = context.model.ThetaW();

   // charge renormalization constant in the decoupling scheme
   const auto delta_e_dec = -0.5*delta_Z_A_dec; //- tan(theta_W) / 2 * delta_Z_ZA

   // it is still necessary to decine the W-and Z-mass as well as the vev
   //vev renormalization constant in the decoupling scheme
   const auto delta_Z_vev_dec = (1 - 1 / (2 * Sqr(sin(theta_W)))) * delta_Zm_W_dec + 0.5 * Sqr(cos(theta_W) / sin(theta_W)) * delta_Zm_Z_dec - delta_e_dec;

   // End Additional Renormalization constants
   //**************************************************************************

   //**************************************************************************
   // Counter Terms for the Left-and Right Form-Factor

   // counter term for the left form factor in the decoupling scheme
   const auto counter_term_FL1 = 0.5/vev * (mL1OS * delta_ZL_lepton_dec + mL2OS * delta_ZR_lepton_dec) + mL1OS / vev * delta_Zm_lepton_dec ;
   const auto counter_term_FL2 = mL1OS/vev * (delta_Z_higgs_dec - delta_Z_vev_dec);

   const auto counter_term_FL = counter_term_FL1 + counter_term_FL2;

   // counter term for the right form factor in the decoupling scheme
   const auto counter_term_FR1 = 0.5/vev * (mL1OS * delta_ZR_lepton_dec + mL2OS * delta_ZL_lepton_dec) + mL1OS / vev * delta_Zm_lepton_dec ;
   const auto counter_term_FR2 = mL1OS/vev * (delta_Z_higgs_dec - delta_Z_vev_dec);

   const auto counter_term_FR = counter_term_FR1 + counter_term_FR2;

   //End Left-and Right Form_Factors

      const auto amp_1l = calculate_amplitude_1l<Higgs, bar<ChargedLepton>::type, ChargedLepton>(context, indexIn, indexOut1, indexOut2);

      // 1loop contributions squared
      double res_1l2 = 0.5 * ps * amp_1l.square() / mHOS;

      // mix term between tree and 1loop amplitude
      double res_1l_tree = 0.5 * ps * 2.*std::real(amplitude_interference(amp_tree, amp_1l)) / mHOS;

      // add the 1loop contributions to the decay width --> this is now the MSbar renormalized decay width!!!
      res += res_1l2;
      res += res_1l_tree;

      // copy the structure of a generic Amplitude from the tree-level amplitude --> A ~ FL * PL + FR * PR
      Decay_amplitude_SFF amp_counterterm = amp_tree;


      // adjust the form factors to hold the counter terms for the left and right form-factors
      amp_counterterm.form_factor_left = counter_term_FL;
      amp_counterterm.form_factor_right = counter_term_FR;

      if (indexOut1.at(0) == 1 && indexOut2.at(0) == 1) {

         //std::cout << "\n model lam: " << context.model.get_Lambdax() << '\n';
         //std::cout << " sm lam: " << sm_context.model.get_Lambdax() << '\n';

         //std::cout << "\n model mH: " << context.physical_mass<Higgs>(indexIn) << '\n';
         //std::cout << "\n sm mH: " << sm_context.physical_mass<smns::fields::hh>(indexIn) << '\n';

         //std::cout << "\n Autogenerated Self-Energies: " << '\n';
         //std::cout << "\n Higgs Self-Energies: \n";

         //std::cout << "\n Higgs SM: " << higgs_all_contribution << '\n';
         //std::cout << " Higgs SM deriv: " << higgs_all_contribution_deriv_p2 << '\n';

         //std::cout << "\n Fermion Self-Energies: \n";

         //std::cout << "\n Fermion SM PL: " << lepton_PL_all_contribution << '\n';
         //std::cout << " Fermion SM PR: " << lepton_PR_all_contribution << '\n';
         //std::cout << " Fermion SM 1: " << lepton_1_all_contribution << '\n';

         //std::cout << "\n Fermion SM PL deriv: " << fermion_PL_all_contribution_deriv_p2 << '\n';
         //std::cout << " Fermion SM PR deriv: " << fermion_PR_all_contribution_deriv_p2 << '\n';
         //std::cout << " Fermion SM 1 deriv: " << fermion_1_all_contribution_deriv_p2 << '\n';

         //std::cout << "\n Z Self-Energies: \n";
         //std::cout << "\n Z SM: " << Z_all_contribution << '\n';
         //std::cout << " Z SM deriv: " << Z_all_contribution_deriv_p2 << '\n';

         //std::cout << "\n Photon Self-Energies: \n";
         //std::cout << "\n P SM: " << P_all_contribution << '\n';
         //std::cout << " P SM deriv: " << P_all_contribution_deriv_p2 << '\n';

         // Build-In Self-Energies
         //std::cout << "\n Build-in SM Self-Energies: " << '\n';
         //std::cout << "\n Higgs Self-Energies: \n";

         //std::cout << "\n Higgs SM: " << higgs_sm_contribution << '\n';
         //std::cout << " Higgs SM deriv: " << higgs_sm_contribution_deriv_p2 << '\n';

         //std::cout << "\n Fermion Self-Energies: \n";

         //std::cout << "\n Fermion SM PL: " << fermion_PL_sm_contribution << '\n';
         //std::cout << " Fermion SM PR: " << fermion_PR_sm_contribution << '\n';
         //std::cout << " Fermion SM 1: " << fermion_1_sm_contribution << '\n';

         //std::cout << "\n Fermion SM PL deriv: " << fermion_PL_sm_contribution_deriv_p2 << '\n';
         //std::cout << " Fermion SM PR deriv: " << fermion_PR_sm_contribution_deriv_p2 << '\n';
         //std::cout << " Fermion SM 1 deriv: " << fermion_1_sm_contribution_deriv_p2 << '\n';

         //std::cout << "\n Z Self-Energies: \n";
         //std::cout << "\n Z SM: " << Z_sm_contribution << '\n';
         //std::cout << " Z SM deriv: " << Z_sm_contribution_deriv_p2 << '\n';

         //std::cout << "\n Photon Self-Energies: \n";
         //std::cout << "\n P SM: " << P_sm_contribution << '\n';
         //std::cout << " P SM deriv: " << P_sm_contribution_deriv_p2 << '\n';

         // Renormalization Constants

         std::cout << "\n Decoupling Higgs Field Renormalization: " << delta_Z_higgs_dec << '\n';
         std::cout << " Decoupling Higgs Mass Renormalization: " << delta_Zm_higgs_dec << '\n';

         std::cout << "\n Decoupling Fermion PL Field Renormalization: " << delta_ZL_lepton_dec << '\n';
         std::cout << " Decoupling Fermion PR Field Renormalization: " << delta_ZR_lepton_dec << '\n';
         std::cout << " Decoupling Fermion Mass Renormalization: " << delta_Zm_lepton_dec << '\n';

         std::cout << "\n Decoupling Z Field Renormalization: " << delta_Z_Z_dec << '\n';
         std::cout << " Decoupling Z Mass Renormalization: " << delta_Zm_Z_dec << '\n';

         std::cout << "\n Decoupling W Field Renormalization: " << delta_Z_W_dec << '\n';
         std::cout << " Decoupling W Mass Renormalization: " << delta_Zm_W_dec << '\n';

         std::cout << "\n Decoupling Charge Renormalization: " << delta_e_dec << '\n';
         std::cout << "\n Decoupling vev Renormaliztion: " << delta_Z_vev_dec << '\n';


         //std::cout << "\n Counter Term 1 left form-factor: " << counter_term_FL1 << '\n';
         //std::cout << " Counter Term 1 right form-factor: " << counter_term_FR1 << '\n';

         //std::cout << "\n Counter Term 2 left form-factor: " << counter_term_FL2 << '\n';
         //std::cout << " Counter Term 2 right form-factor: " << counter_term_FR2 << '\n';

         //std::cout << "\n counter term left: " << counter_term_FL << '\n';
         //std::cout << " counter term right: " << counter_term_FR << '\n';

         //std::cout << "\n Counter Term left form-factor struct: " << amp_counterterm.form_factor_left << '\n';
         //std::cout << " Counter Term right form-factor struct: " << amp_counterterm.form_factor_right << '\n';

         std::cout << "\n counter term left struct: " << amp_counterterm.form_factor_left << '\n';
         std::cout << "\n counter term right struct: " << amp_counterterm.form_factor_right << '\n';

         std::cout << "\n form factor left tree: " << amp_tree.form_factor_left << '\n';
         std::cout << " form factor right tree: " << amp_tree.form_factor_right << '\n';

         std::cout << "\n form factor left 1loop: " << amp_1l.form_factor_left << '\n';
         std::cout << " form factor right 1loop: " << amp_1l.form_factor_right << '\n';

         std::cout << "\n form factor left dec: " << amp_1l.form_factor_left + counter_term_FL << '\n';
         std::cout << " form factor right dec: " << amp_1l.form_factor_right + counter_term_FR << '\n';

         std::cout << "\n 1l* conj(tree): " << amplitude_interference(amp_1l, amp_tree) << '\n';
         std::cout << " |tree|^2: " << amp_tree.square() << '\n';
         std::cout << " |1l|^2: " << amp_1l.square() << '\n';
         std::cout << " |A|^2: " << amp_tree.square() + amp_1l.square() + 2.0 * std::real(amplitude_interference(amp_1l, amp_tree)) << '\n';

         std::cout << "\n Decay Width tree: " << 0.5 * ps * amp_tree.square() / mHOS << '\n';
         std::cout << " Decay Width MSbar: " << res << '\n';
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
