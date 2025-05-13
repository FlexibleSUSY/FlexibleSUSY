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
   const auto betaOS2_S = 1.0 - Sqr((mL1OS + mL2OS)/mHOS);
   const auto betaOS2_P = 1.0 - Sqr((mL1OS - mL2OS)/mHOS);
   const double ps = 1./(8.*Pi)*std::sqrt(KallenLambda(1., xOS1, xOS2));
   const double flux = 1./(2.*mHOS);

   const double amp2_S = Sqr(mHOS) * betaOS2_S * 2*std::norm(HLLbarVertexDR_S);

   const double amp2_P =
      info::is_CP_violating_Higgs_sector ? Sqr(mHOS) * betaOS2_P * 2.*std::norm(HLLbarVertexDR_P) : 0.;

   // flux * phase space factor * symmetry factor * |matrix element|^2
   double res_S = flux * ps * amp2_S;
   double res_P = flux * ps * amp2_P;

   const auto amp_tree =
      calculate_amplitude_tree<Higgs, bar<ChargedLepton>::type, ChargedLepton>(context, indexIn, indexOut1, indexOut2);

   // higher order SM corrections
   if (flexibledecay_settings.get_decay_corrections().sm > 0) {
      // 1-loop QED corrections
      res_S *= 1. + get_alpha(context)/Pi*17./4.;
      res_P *= 1. + get_alpha(context)/Pi*17./4.;
   }

   double res = res_S + res_P;

#if DECAY_BSM_ORDER > 0
   if (flexibledecay_settings.get_decay_corrections().bsm > 0) {

      namespace smns = standard_model_cxx_diagrams;
      smns::context_base sm_context {sm};

      // self-energies, to be moved outside

      // lepton self-energy [PL]
      const auto lepton_se_all_PL =
         self_energy_1loop_PL<ChargedLepton>(context, mL1OS, indexOut1, indexOut2);
      const auto lepton_se_sm_PL =
         smns::self_energy_1loop_PL<smns::fields::Fe>(sm_context, mL1OS).coeff(indexOut1.at(0), indexOut2.at(0));
      const auto lepton_se_dec_PL = lepton_se_all_PL - lepton_se_sm_PL;

      // lepton self-energy [PR]
      const auto lepton_se_all_PR =
         self_energy_1loop_PR<ChargedLepton>(context, mL1OS, indexOut1, indexOut2);
      const auto lepton_se_sm_PR =
         smns::self_energy_1loop_PR<smns::fields::Fe>(sm_context, mL1OS).coeff(indexOut1.at(0), indexOut2.at(0));
      const auto lepton_se_dec_PR =
         lepton_se_all_PR - lepton_se_sm_PR;

      // lepton self-energy [1]
      const auto lepton_1_all =
         self_energy_1loop_1<ChargedLepton>(context, mL1OS, indexOut1, indexOut2);
      const auto lepton_1_sm =
         smns::self_energy_1loop_1<smns::fields::Fe>(sm_context, mL1OS).coeff(indexOut1.at(0), indexOut2.at(0));
      const auto lepton_1_se_dec =
         lepton_1_all - lepton_1_sm;

      // lepton self-energy derivative [PL]
      const auto lepton_se_deriv_all_PL =
         self_energy_1loop_PL_deriv_p2<ChargedLepton>(context, mL1OS, indexOut1, indexOut2);
      const auto lepton_se_deriv_sm_PL =
         smns::self_energy_1loop_PL_deriv_p2<smns::fields::Fe>(sm_context, mL1OS).coeff(indexOut1.at(0), indexOut2.at(0));
      const auto lepton_se_deriv_dec_PL =
         lepton_se_deriv_all_PL - lepton_se_deriv_sm_PL;

      // lepton self-energy derivative [PR]
      const auto lepton_se_deriv_all_PR =
         self_energy_1loop_PR_deriv_p2<ChargedLepton>(context, mL1OS, indexOut1, indexOut2);
      const auto lepton_se_deriv_sm_PR =
         smns::self_energy_1loop_PR_deriv_p2<smns::fields::Fe>(sm_context, mL1OS).coeff(indexOut1.at(0), indexOut2.at(0));
      const auto lepton_se_deriv_dec_PR =
         lepton_se_deriv_all_PR - lepton_se_deriv_sm_PR;

      // lepton self-energy derivative [1]
      const auto lepton_se_deriv_all_1 =
         self_energy_1loop_1_deriv_p2<ChargedLepton>(context, mL1OS, indexOut1, indexOut2);
      const auto lepton_se_deriv_sm_1 =
         smns::self_energy_1loop_1_deriv_p2<smns::fields::Fe>(sm_context, mL1OS).coeff(indexOut1.at(0), indexOut2.at(0));
      const auto lepton_se_deriv_dec_1 =
         lepton_se_deriv_all_1 - lepton_se_deriv_sm_1;

      // higgs se
      const auto higgs_se_all = self_energy_1loop<Higgs>(context, mHOS, indexIn, indexIn);
      const auto higgs_se_sm = smns::self_energy_1loop<smns::fields::hh>(sm_context, mHOS);
      const auto higgs_se_dec = higgs_se_all - higgs_se_sm;

      // higgs se derivative
      const auto higgs_se_deriv_all = self_energy_1loop_deriv_p2<Higgs>(context, mHOS, indexIn, indexIn);
      const auto higgs_se_deriv_sm = smns::self_energy_1loop_deriv_p2<smns::fields::hh>(sm_context, mHOS);
      const auto higgs_se_deriv_dec = higgs_se_deriv_all - higgs_se_deriv_sm;

      // photon se
      //const auto photon_se_all = self_energy_1loop<Photon>(context, 0, {}, {});
      //const auto photon_se_sm = smns::self_energy_1loop<smns::fields::VP>(sm_context, 0);
      //const auto photon_se_dec = photon_se_all - photon_se_sm;

      // photon se derivative
      const auto photon_se_deriv_all = self_energy_1loop_deriv_p2<Photon>(context, 0, {}, {});
      const auto photon_se_deriv_sm = smns::self_energy_1loop_deriv_p2<smns::fields::VP>(sm_context, 0);
      const auto photon_se_deriv_dec = photon_se_deriv_all - photon_se_deriv_sm;

      // Z se
      const double mZOS = context.physical_mass<ZBoson>({});
      const auto Z_se_all = self_energy_1loop<ZBoson>(context, mZOS, {}, {});
      const auto Z_se_sm = smns::self_energy_1loop<smns::fields::VZ>(sm_context, mZOS);
      const auto Z_se_dec = Z_se_all - Z_se_sm;

      // Z se derivative
      const auto Z_se_deriv_all = self_energy_1loop_deriv_p2<ZBoson>(context, mZOS, {}, {});
      const auto Z_se_deriv_sm = smns::self_energy_1loop_deriv_p2<smns::fields::VZ>(sm_context, mZOS);
      const auto Z_se_deriv_dec = Z_se_deriv_all - Z_se_deriv_sm;

      // W se
      const double mWOS = context.physical_mass<WpBoson>({});
      const auto W_se_all = self_energy_1loop<WpBoson>(context, mWOS, {}, {});
      const auto W_se_sm = smns::self_energy_1loop<smns::fields::VWp>(sm_context, mWOS);
      const auto W_se_dec = W_se_all - W_se_sm;

      // W se derivative
      const auto W_se_deriv_all = self_energy_1loop_deriv_p2<WpBoson>(context, mWOS, {}, {});
      const auto W_se_deriv_sm = smns::self_energy_1loop_deriv_p2<smns::fields::VWp>(sm_context, mWOS);
      const auto W_se_dec_deriv_p2 = W_se_deriv_all - W_se_deriv_sm;

      // --------------------------------------------------------------------------------------

      const auto delta_ZL_lepton_dec =
         -2*lepton_se_dec_PL - 2*mL1OS*mL1OS*(lepton_se_deriv_dec_PL + lepton_se_deriv_dec_PR) - mL1OS*(lepton_se_deriv_dec_1 + lepton_se_deriv_dec_1);
      const auto delta_ZR_lepton_dec =
         -2*lepton_se_dec_PR - 2*mL1OS*mL1OS*(lepton_se_deriv_dec_PL + lepton_se_deriv_dec_PR) - mL1OS*(lepton_se_deriv_dec_1 + lepton_se_deriv_dec_1);
      const auto delta_Zm_lepton_dec =
         (lepton_se_dec_PR + lepton_se_dec_PR) + 0.5/mL1OS*(lepton_1_se_dec + lepton_1_se_dec);

      const auto delta_Z_higgs_dec = -higgs_se_deriv_dec;
      const auto delta_Zm_higgs_dec = higgs_se_dec/Sqr(mHOS);

      const auto delta_Z_W_dec = -W_se_dec_deriv_p2;
      const auto delta_Zm_W_dec = W_se_dec/Sqr(mWOS);

      const auto delta_Z_Z_dec = -Z_se_deriv_dec;
      const auto delta_Zm_Z_dec = Z_se_dec/Sqr(mZOS);

      const auto delta_Z_A_dec = -photon_se_deriv_dec;

      const double theta_W = context.model.ThetaW();

      const auto delta_e_dec = -0.5*delta_Z_A_dec; //- tan(theta_W) / 2 * delta_Z_ZA

      // it is still necessary to decine the W-and Z-mass as well as the vev
      // vev renormalization constant in the decoupling scheme
      const auto delta_Z_vev_dec =
         (1 - 0.5/Sqr(std::sin(theta_W)))*delta_Zm_W_dec + 0.5/Sqr(std::tan(theta_W))*delta_Zm_Z_dec - delta_e_dec;


      // --------------------------------------------------------------------------------------

      const double vev = 248.9397360685197;
      const auto ct_FL1 =
         amp_tree.form_factor_left * (0.5*(delta_ZL_lepton_dec + delta_ZR_lepton_dec) + delta_Zm_lepton_dec);
      const auto ct_FL2 =
         amp_tree.form_factor_left * (delta_Z_higgs_dec - delta_Z_vev_dec);
      const auto ct_FL =
         ct_FL1 + ct_FL2;

      const auto ct_FR1 =
         amp_tree.form_factor_right * (0.5*(delta_ZR_lepton_dec + delta_ZL_lepton_dec) + delta_Zm_lepton_dec);
      const auto ct_FR2 =
        amp_tree.form_factor_right * (delta_Z_higgs_dec - delta_Z_vev_dec);
      const auto ct_FR = ct_FR1 + ct_FR2;

      // --------------------------------------------------------------------------------------

      const auto amp_1l = calculate_amplitude_1l<Higgs, bar<ChargedLepton>::type, ChargedLepton>(context, indexIn, indexOut1, indexOut2);
      const auto amp_1l_sm = sm_decays.calculate_amplitude_1l<smns::fields::hh, bar<smns::fields::Fe>::type, smns::fields::Fe>(sm_context, indexIn, indexOut1, indexOut2);
      const auto amp_1l_dec = amp_1l-amp_1l_sm;

      // mix term between tree and 1loop amplitude
      double amp_1l_tree = 2.*std::real(amplitude_interference(amp_tree, amp_1l_dec));

      // if the interference needs to be addded or subtracted dependes on the convention and still needs to be figured out
      res -= flux * ps * amp_1l_tree;

      // copy the structure of a generic Amplitude from the tree-level amplitude --> A ~ FL * PL + FR * PR
      Decay_amplitude_SFF amp_counterterm = amp_tree;

      // adjust the form factors to hold the counter terms for the left and right form-factors
      amp_counterterm.form_factor_left = ct_FL;
      amp_counterterm.form_factor_right = ct_FR;

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


         //std::cout << "\n Counter Term 1 left form-factor: " << ct_FL1 << '\n';
         //std::cout << " Counter Term 1 right form-factor: " << ct_FR1 << '\n';

         //std::cout << "\n Counter Term 2 left form-factor: " << ct_FL2 << '\n';
         //std::cout << " Counter Term 2 right form-factor: " << ct_FR2 << '\n';

         //std::cout << "\n counter term left: " << ct_FL << '\n';
         //std::cout << " counter term right: " << ct_FR << '\n';

         //std::cout << "\n Counter Term left form-factor struct: " << amp_counterterm.form_factor_left << '\n';
         //std::cout << " Counter Term right form-factor struct: " << amp_counterterm.form_factor_right << '\n';

         std::cout << "\n counter term left struct: " << amp_counterterm.form_factor_left << '\n';
         std::cout << "\n counter term right struct: " << amp_counterterm.form_factor_right << '\n';

         std::cout << "\n form factor left tree: " << amp_tree.form_factor_left << '\n';
         std::cout << " form factor right tree: " << amp_tree.form_factor_right << '\n';

         std::cout << "\n form factor left 1loop: " << amp_1l.form_factor_left << '\n';
         std::cout << " form factor right 1loop: " << amp_1l.form_factor_right << '\n';

         std::cout << "\n form factor left dec: " << amp_1l.form_factor_left + ct_FL << '\n';
         std::cout << " form factor right dec: " << amp_1l.form_factor_right + ct_FR << '\n';

         std::cout << "\n 1l* conj(tree): " << amplitude_interference(amp_1l, amp_tree) << '\n';
         std::cout << " |tree|^2: " << amp_tree.square() << '\n';
         std::cout << " |1l|^2: " << amp_1l.square() << '\n';
         std::cout << " |A|^2: " << amp_tree.square() + amp_1l.square() + 2.0 * std::real(amplitude_interference(amp_1l, amp_tree)) << '\n';

         std::cout << "\n Decay Width tree: " << 0.5 * ps * amp_tree.square() / mHOS << '\n';
         std::cout << " Decay Width MSbar: " << res << '\n';
      }

      res += flux * ps * 2.*std::real(amplitude_interference(amp_tree, amp_counterterm));
   }
#endif

   if (flexibledecay_settings.get(FlexibleDecay_settings::call_higgstools) != 0 ||
       flexibledecay_settings.get(FlexibleDecay_settings::call_lilith) != 0 ||
       flexibledecay_settings.get(FlexibleDecay_settings::calculate_normalized_effc) != 0
   ) {
      const double remove_normalization = flux*ps*2*Sqr(mHOS);
      neutral_higgs_effc.add_coupling(
         field_as_string<Higgs>(indexIn),
            {
               boost::hana::unpack(bar<ChargedLepton>::pdgids, _to_array<bar<ChargedLepton>::numberOfGenerations>).at(indexOut1.at(0)),
               boost::hana::unpack(ChargedLepton::pdgids, _to_array<ChargedLepton::numberOfGenerations>).at(indexOut2.at(0))
            },
         std::pair<std::string, std::complex<double>> {
            field_as_string<Higgs>(indexIn) + "-" + field_as_string<bar<ChargedLepton>::type>(indexOut1) + "-" + field_as_string<ChargedLepton>(indexOut2),
            std::sqrt(res_S/(remove_normalization*betaOS2_S)) + 1i*std::sqrt(res_P/(remove_normalization*betaOS2_P))
         },
         res_S + res_P
      );
   }

   return res;
}
