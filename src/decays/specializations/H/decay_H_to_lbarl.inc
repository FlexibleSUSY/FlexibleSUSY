#include <iostream>

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

   std::cout << std::setprecision(20) << '\n';
   //std::cout << " Tree Level Amp: " << amp_tree << '\n';

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

      const auto [gamma_se_1, gamma_se_PL, gamma_se_PR] = subtract_1l_QED_se(context, indexOut1.at(0), mHOS);
      const auto [gamma_se_deriv_1, gamma_se_deriv_PL, gamma_se_deriv_PR] = subtract_1l_QED_se_deriv(context, indexOut1.at(0), mHOS);
      const auto gamma_triangle = subtract_1l_QED_triangle<Higgs>(context, indexIn, indexOut1.at(0), indexOut2.at(0)) * amp_tree;

      //-----------------------------------------------------------------------------------------------------------------
      // Explicit photon diagrams (Build-in SM)
      //-----------------------------------------------------------------------------------------------------------------

      auto& lib = Loop_library::get();
      auto& libSM = Loop_library::get();
      looplibrary::Bcoeff_t B_coeffsm;
      looplibrary::Ccoeff_t C_coeffsm;

      const double scaleSM = sm.get_scale();
      const double mHOSSM = sm.get_physical().Mhh;
      const double mL1OSSM = sm.get_physical().MFe(indexOut1.at(0));
      const double mL2OSSM = sm.get_physical().MFe(indexOut2.at(0));
      const double mL1TreeSM = sm.get_MFe(indexOut1.at(0));
      const double mL2TreeSM = sm.get_MFe(indexOut2.at(0));

      // B,C functions for triangle diagram, evaluated at physical and tree masses
      libSM.B(B_coeffsm, Sqr(mL1OSSM), 0, Sqr(mL1TreeSM), Sqr(scaleSM));
      libSM.C(C_coeffsm, Sqr(mL1OSSM), Sqr(mL2OSSM), Sqr(mHOSSM), Sqr(mL1TreeSM), 0, Sqr(mL2TreeSM), Sqr(scaleSM));

      const auto b0SM = B_coeffsm.at(0);
      const auto b1SM = B_coeffsm.at(1);
      const auto db0SM = lib.DB0(Sqr(mL1OSSM), 0, Sqr(mL2TreeSM), Sqr(scaleSM));
      const auto db1SM = lib.DB1(Sqr(mL1OSSM), 0, Sqr(mL2TreeSM), Sqr(scaleSM));

      const auto c0SM = C_coeffsm.at(0);
      const auto c1SM = C_coeffsm.at(1);
      const auto c2SM = C_coeffsm.at(2);

      const auto chargeSM = unit_charge(sm_context);

      const auto gamma_se_PL_sm = oneOver16PiSqr * Sqr(chargeSM) * (2 * b0SM + 2 * b1SM - 1);
      const auto gamma_se_PR_sm = oneOver16PiSqr * Sqr(chargeSM) * (2 * b0SM + 2 * b1SM - 1);
      const auto gamma_se_1_sm = oneOver16PiSqr * Sqr(chargeSM) * (-2) * (2 * b0SM - 1) * mL1TreeSM;

      const auto gamma_se_deriv_PL_sm = oneOver16PiSqr * Sqr(chargeSM) * (2 * db0SM + 2 * b1SM);
      const auto gamma_se_deriv_PR_sm = oneOver16PiSqr * Sqr(chargeSM) * (2 * db0SM + 2 * db1SM);
      const auto gamma_se__deriv_1_sm = oneOver16PiSqr * Sqr(chargeSM) * (-2) * (2 * db0SM) * mL1TreeSM;

      // amp_tree must still be changed
      const auto gamma_triangle_sm = oneOver16PiSqr * Sqr(chargeSM) * (-2) * amp_tree
                              * (-1 + 2 * b0SM
                                 + c2SM * (2 * Sqr(mHOSSM) + (mL1TreeSM + mL2TreeSM) * (mL1OSSM - mL2OSSM))
                                 + c0SM * mL1TreeSM * (2 * mL1TreeSM + 2 * mL2TreeSM + mL1OSSM - mL2OSSM)
                                 + c1SM * (Sqr(mHOSSM) + mL1OSSM * (mL1TreeSM - mL2TreeSM) + Sqr(mL1OSSM) - Sqr(mL2OSSM)));


      //---------------------------------------------------------------------------------------------------


      //std::cout << " SM scale: " << sm.get_scale() << '\n';

      //std::cout << "\n SM physical Higgs Mass: " << sm.get_physical().Mhh << '\n';
      //std::cout << " SM physical Z Mass: " << sm.get_physical().MVZ << '\n';
      //std::cout << " SM physical W Mass: " << sm.get_physical().MVWp << '\n';
      //std::cout << " SM physical e Mass: " << sm.get_physical().MFe(0) << '\n';
      //std::cout << " SM physical mu Mass: " << sm.get_physical().MFe(1) << '\n';
      //std::cout << " SM physical tau Mass: " << sm.get_physical().MFe(2) << '\n';

      //std::cout << "\n SM Higgs Mass: " << sm.get_Mhh() << '\n';
      //std::cout << " SM Z Mass: " << sm.get_MVZ() << '\n';
      //std::cout << " SM W Mass: " << sm.get_MVWp() << '\n';
      //std::cout << " SM lepton Mass: " << sm.get_MFe() << '\n';
      //std::cout << " SM up Mass: " << sm.get_MFu() << '\n';
      //std::cout << " SM down Mass: " << sm.get_MFd() << '\n';

      //std::cout  << " SM g1: " << sm.get_g1() << '\n';
      //std::cout  << " SM g2: " << sm.get_g2() << '\n';
      //std::cout  << " SM v: " << sm.get_v() << '\n';


      const double mHSMOS = sm.get_physical().Mhh;
      const double mL1SMOS = sm.get_physical().MFe(indexOut1.at(0));
      const double mL2SMOS = sm.get_physical().MFe(indexOut2.at(0));

      // self-energies, to be moved outside

      // lepton self-energy [PL]
      auto lepton_se_all_PL =
         self_energy_1loop_PL<ChargedLepton,ChargedLepton>(context, mL1OS, indexOut1, indexOut2);

      lepton_se_all_PL -= 0.5 * gamma_se_PL; // Removal of photon contributions

      auto lepton_se_sm_PL =
         sm.self_energy_Fe_1loop_PL(mL1SMOS, indexOut1.at(0), indexOut2.at(0));

      lepton_se_sm_PL -= 0.5 * gamma_se_PL_sm; // removal of photon contribution

      const auto lepton_se_dec_PL = lepton_se_all_PL - lepton_se_sm_PL;

      // lepton self-energy [PR]
      auto lepton_se_all_PR =
         self_energy_1loop_PR<ChargedLepton,ChargedLepton>(context, mL1OS, indexOut1, indexOut2);

      lepton_se_all_PR -= 0.5 * gamma_se_PR; // removal of photon contribution

      auto lepton_se_sm_PR =
         sm.self_energy_Fe_1loop_PR(mL1SMOS, indexOut1.at(0), indexOut2.at(0));

      lepton_se_sm_PR -= 0.5 * gamma_se_PR_sm; // removal of photon contribution

      const auto lepton_se_dec_PR =
         lepton_se_all_PR - lepton_se_sm_PR;

      // lepton self-energy [1]
      auto lepton_1_all =
         self_energy_1loop_1<ChargedLepton,ChargedLepton>(context, mL1OS, indexOut1, indexOut2);

      lepton_1_all -= gamma_se_1; // removal of photon contibution

      auto lepton_1_sm =
         sm.self_energy_Fe_1loop_1(mL1SMOS, indexOut1.at(0), indexOut2.at(0));

      lepton_1_sm -= gamma_se_1_sm; // removal of photon contribution

      const auto lepton_1_se_dec =
         lepton_1_all - lepton_1_sm;

      // lepton self-energy derivative [PL]
      auto lepton_se_deriv_all_PL =
         self_energy_1loop_PL_deriv_p2<ChargedLepton>(context, mL1OS, indexOut1, indexOut2);

      lepton_se_deriv_all_PL -= 0.5 * gamma_se_deriv_PL; // removal of photon contribution

      auto lepton_se_deriv_sm_PL =
         sm.self_energy_Fe_1loop_PL_deriv_p2(mL1SMOS, indexOut1.at(0), indexOut2.at(0));

      lepton_se_deriv_sm_PL -= 0.5 * gamma_se_deriv_PL_sm; // removal of photon contribution

      const auto lepton_se_deriv_dec_PL =
         lepton_se_deriv_all_PL - lepton_se_deriv_sm_PL;

      // lepton self-energy derivative [PR]
      auto lepton_se_deriv_all_PR =
         self_energy_1loop_PR_deriv_p2<ChargedLepton>(context, mL1OS, indexOut1, indexOut2);

      lepton_se_deriv_all_PR -= 0.5 * gamma_se_deriv_PR; // removal of photon contribution

      auto lepton_se_deriv_sm_PR =
         sm.self_energy_Fe_1loop_PR_deriv_p2(mL1SMOS, indexOut1.at(0), indexOut2.at(0));

      lepton_se_deriv_sm_PR -= 0.5 * gamma_se_deriv_PR_sm; // removal of photon contribution

      const auto lepton_se_deriv_dec_PR =
         lepton_se_deriv_all_PR - lepton_se_deriv_sm_PR;

      // lepton self-energy derivative [1]
      auto lepton_se_deriv_all_1 =
         self_energy_1loop_1_deriv_p2<ChargedLepton>(context, mL1OS, indexOut1, indexOut2);

      lepton_se_deriv_all_1 -= gamma_se_deriv_1; // removal of photon contribution

      auto lepton_se_deriv_sm_1 =
         sm.self_energy_Fe_1loop_1_deriv_p2(mL1SMOS, indexOut1.at(0), indexOut2.at(0));

      lepton_se_deriv_sm_1 -= gamma_se__deriv_1_sm; // removal of photon contribution

      const auto lepton_se_deriv_dec_1 =
         lepton_se_deriv_all_1 - lepton_se_deriv_sm_1;

      // higgs se
      const auto higgs_se_all = self_energy_1loop<Higgs,Higgs>(context, mHOS, indexIn, indexIn);
      const auto higgs_se_sm = sm.self_energy_hh_1loop(mHSMOS);
      const auto higgs_se_dec = higgs_se_all - higgs_se_sm;

      // higgs se derivative
      const auto higgs_se_deriv_all = self_energy_1loop_deriv_p2<Higgs>(context, mHOS, indexIn, indexIn);
      const auto higgs_se_deriv_sm = sm.self_energy_hh_1loop_deriv_p2(mHSMOS);
      const auto higgs_se_deriv_dec = higgs_se_deriv_all - higgs_se_deriv_sm;

      // photon se
      const auto photon_se_all = self_energy_1loop<Photon,Photon>(context, 0, {}, {});
      const auto photon_se_sm = sm.self_energy_VWp_1loop(0);
      const auto photon_se_dec = photon_se_all - photon_se_sm;

      // photon se derivative
      const auto photon_se_deriv_all = self_energy_1loop_deriv_p2<Photon>(context, 0, {}, {});
      const auto photon_se_deriv_sm = sm.self_energy_VP_1loop_deriv_p2(0);
      const auto photon_se_deriv_dec = photon_se_deriv_all - photon_se_deriv_sm;

      // Z se
      const double mZOS = context.physical_mass<ZBoson>({});
      const double mZSMOS = sm.get_physical().MVZ;
      const auto Z_se_all = self_energy_1loop<ZBoson,ZBoson>(context, mZOS, {}, {});
      const auto Z_se_sm = sm.self_energy_VZ_1loop(mZSMOS);
      const auto Z_se_dec = Z_se_all - Z_se_sm;

      // Z se derivative
      const auto Z_se_deriv_all = self_energy_1loop_deriv_p2<ZBoson>(context, mZOS, {}, {});
      const auto Z_se_deriv_sm =  sm.self_energy_VZ_1loop_deriv_p2(mZSMOS);
      const auto Z_se_deriv_dec = Z_se_deriv_all - Z_se_deriv_sm;

      // W se
      const double mWOS = context.physical_mass<WpBoson>({});
      const double mWSMOS = sm.get_physical().MVWp;
      const auto W_se_all = self_energy_1loop<WpBoson,WpBoson>(context, mWOS, {}, {});
      const auto W_se_sm = sm.self_energy_VWp_1loop(mWSMOS);
      const auto W_se_dec = W_se_all - W_se_sm;

      // W se derivative
      const auto W_se_deriv_all = self_energy_1loop_deriv_p2<WpBoson>(context, mWOS, {}, {});
      const auto W_se_deriv_sm = sm.self_energy_VWp_1loop_deriv_p2(mWSMOS);
      const auto W_se_dec_deriv_p2 = W_se_deriv_all - W_se_deriv_sm;

      // photon-Z se
      
      const auto PZ_se_all = self_energy_1loop<Photon,ZBoson>(context, 0, {}, {});
      const auto PZ_se_sm = sm.self_energy_VPVZ_1loop(0);

      // --------------------------------------------------------------------------------------
      // SM On-Shell Renormalization Constants
      // --------------------------------------------------------------------------------------

      const double mLSMTree = sm.get_MFe(indexOut2.at(0));

      const auto delta_ZL_lepton_sm =
         -2*lepton_se_sm_PL - 2*mLSMTree*mLSMTree*(lepton_se_deriv_sm_PL + lepton_se_deriv_sm_PR) - mLSMTree*(lepton_se_deriv_sm_1 + lepton_se_deriv_sm_1);
      const auto delta_ZR_lepton_sm =
         -2*lepton_se_sm_PR - 2*mLSMTree*mLSMTree*(lepton_se_deriv_sm_PL + lepton_se_deriv_sm_PR) - mLSMTree*(lepton_se_deriv_sm_1 + lepton_se_deriv_sm_1);
      const auto delta_Zm_lepton_sm =
         (lepton_se_sm_PR + lepton_se_sm_PR) + 0.5/mLSMTree*(lepton_1_sm + lepton_1_sm);

      const auto delta_Z_higgs_sm = -higgs_se_deriv_sm;
      const auto delta_Zm_higgs_sm = higgs_se_sm/Sqr(mHSMOS);

      const auto delta_Z_W_sm = -W_se_deriv_sm;
      const auto delta_Zm_W_sm = W_se_sm/Sqr(mWSMOS);

      const auto delta_Z_Z_sm = -Z_se_deriv_sm;
      const auto delta_Zm_Z_sm = Z_se_sm/Sqr(mZSMOS);

      const auto delta_Z_A_sm = -photon_se_deriv_sm;
      
      // still needs the correct build-in ZA self energy
      const auto delta_Z_ZA_sm = 2 * PZ_se_sm / Sqr(mZSMOS);

      const double theta_W_sm = sm.ThetaW();

      const auto delta_e_sm = -0.5*delta_Z_A_sm - tan(theta_W_sm) / 2 * delta_Z_ZA_sm;

      // it is still necessary to decine the W-and Z-mass as well as the vev
      // vev renormalization constant in the decoupling scheme
      const auto delta_Z_vev_sm =
         (1 - 0.5/Sqr(std::sin(theta_W_sm)))*delta_Zm_W_sm + 0.5/Sqr(std::tan(theta_W_sm))*delta_Zm_Z_sm - delta_e_sm;
      
      // --------------------------------------------------------------------------------------
      // BSM On-Shell Renormalization Constants
      // --------------------------------------------------------------------------------------

      const double mLTree = context.model.get_MFe(indexOut2.at(0));

      const auto delta_ZL_lepton =
         -2*lepton_se_all_PL - 2*mLTree*mLTree*(lepton_se_deriv_all_PL + lepton_se_deriv_all_PR) - mLTree*(lepton_se_deriv_all_1 + lepton_se_deriv_all_1);
      const auto delta_ZR_lepton =
         -2*lepton_se_all_PR - 2*mLTree*mLTree*(lepton_se_deriv_all_PL + lepton_se_deriv_all_PR) - mLTree*(lepton_se_deriv_all_1 + lepton_se_deriv_all_1);
      const auto delta_Zm_lepton =
         (lepton_se_all_PR + lepton_se_all_PR) + 0.5/mLTree*(lepton_1_all + lepton_1_all);

      const auto delta_Z_higgs = -higgs_se_deriv_all;
      const auto delta_Zm_higgs = higgs_se_all/Sqr(mHOS);

      const auto delta_Z_W = -W_se_deriv_all;
      const auto delta_Zm_W = W_se_all/Sqr(mWOS);

      const auto delta_Z_Z = -Z_se_deriv_all;
      const auto delta_Zm_Z = Z_se_all/Sqr(mZOS);

      const auto delta_Z_A = -photon_se_deriv_all;
      const auto delta_Z_ZA = 2 * PZ_se_all / Sqr(mZOS);

      const double theta_W = context.model.ThetaW();

      const auto delta_e = -0.5*delta_Z_A - tan(theta_W) / 2 * delta_Z_ZA;

      // it is still necessary to decine the W-and Z-mass as well as the vev
      // vev renormalization constant in the decoupling scheme
      const auto delta_Z_vev =
         (1 - 0.5/Sqr(std::sin(theta_W)))*delta_Zm_W + 0.5/Sqr(std::tan(theta_W))*delta_Zm_Z - delta_e;


      // --------------------------------------------------------------------------------------
      // Decoupling Renormalization Constants
      // --------------------------------------------------------------------------------------

      const auto delta_ZL_lepton_dec = delta_ZL_lepton - delta_ZL_lepton_sm;
      const auto delta_ZR_lepton_dec = delta_ZR_lepton - delta_ZR_lepton_sm;
      const auto delta_Zm_lepton_dec = delta_Zm_lepton - delta_Zm_lepton_sm;

      const auto delta_Z_higgs_dec = delta_Z_higgs - delta_Z_higgs_sm;
      const auto delta_Zm_higgs_dec = delta_Zm_higgs - delta_Zm_higgs_sm;

      const auto delta_Z_W_dec = delta_Z_W - delta_Z_W_sm;
      const auto delta_Zm_W_dec = delta_Zm_W - delta_Zm_W_sm;

      const auto delta_Z_Z_dec = delta_Z_Z - delta_Z_Z_sm;
      const auto delta_Zm_Z_dec = delta_Zm_Z - delta_Zm_Z_sm;

      const auto delta_Z_A_dec = delta_Z_A - delta_Z_A_sm;

      const auto delta_e_dec = delta_e - delta_e_sm; //- tan(theta_W) / 2 * delta_Z_ZA

      // it is still necessary to decine the W-and Z-mass as well as the vev
      // vev renormalization constant in the decoupling scheme
      const auto delta_Z_vev_dec = (1 - 0.5/Sqr(std::sin(theta_W)))*delta_Zm_W_dec + 0.5/Sqr(std::tan(theta_W))*delta_Zm_Z_dec - delta_e_dec;


      // --------------------------------------------------------------------------------------

      const auto ct_FL1 =
         amp_tree.form_factor_left * (0.5*(delta_ZL_lepton + delta_ZR_lepton) + delta_Zm_lepton_dec);
      const auto ct_FL2 =
         amp_tree.form_factor_left * (delta_Z_higgs - delta_Z_vev_dec);
      const auto ct_FL =
         ct_FL1 + ct_FL2;

      const auto ct_FR1 =
         amp_tree.form_factor_right * (0.5*(delta_ZR_lepton + delta_ZL_lepton) + delta_Zm_lepton_dec);
      const auto ct_FR2 =
        amp_tree.form_factor_right * (delta_Z_higgs - delta_Z_vev_dec);
      const auto ct_FR = ct_FR1 + ct_FR2;

      // --------------------------------------------------------------------------------------

      auto amp_1l = calculate_amplitude_1l<Higgs, bar<ChargedLepton>::type, ChargedLepton>(context, indexIn, indexOut1, indexOut2);

      amp_1l = amp_1l - gamma_triangle; // removal of photon triangle

      auto amp_1l_sm = sm_decays.calculate_amplitude_1l<smns::fields::hh, bar<smns::fields::Fe>::type, smns::fields::Fe>(sm_context, {}, indexOut1, indexOut2);

      amp_1l_sm = amp_1l_sm - gamma_triangle_sm;

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

         std::cout << "\n model lam: " << context.model.get_Lambdax() << '\n';
         std::cout << " sm lam: " << sm_context.model.get_Lambdax() << '\n';

         //std::cout << "\n model lam: " << context.model.get_Lambdax() << '\n';
         //std::cout << " sm lam: " << sm.get_Lambdax() << '\n';

         //std::cout << "\n model mH: " << context.physical_mass<Higgs>(indexIn) << '\n';
         //std::cout << "\n sm mH: " << sm_context.physical_mass<smns::fields::hh>(indexIn) << '\n';

         std::cout << "\n Autogenerated Self-Energies: " << '\n';
         std::cout << "\n Higgs Self-Energies: \n";

         std::cout << "\n Higgs: " << higgs_se_all << '\n';
         std::cout << " Higgs deriv: " << higgs_se_deriv_all << '\n';

         std::cout << "\n Fermion Self-Energies: \n";

         std::cout << "\n Fermion PL: " << lepton_se_all_PL << '\n';
         std::cout << " Fermion PR: " << lepton_se_all_PR << '\n';
         std::cout << " Fermion 1: " << lepton_1_all << '\n';

         std::cout << "\n Fermion PL deriv: " << lepton_se_deriv_all_PL << '\n';
         std::cout << " Fermion PR deriv: " << lepton_se_deriv_all_PR << '\n';
         std::cout << " Fermion 1 deriv: " << lepton_se_deriv_all_1 << '\n';

         std::cout << "\n Z Self-Energies: \n";
         std::cout << "\n Z: " << Z_se_all << '\n';
         std::cout << " Z deriv: " << Z_se_deriv_all << '\n';

         std::cout << "\n Photon Self-Energies: \n";
         std::cout << "\n P: " << photon_se_all << '\n';
         std::cout << " P deriv: " << photon_se_deriv_all << '\n';

         std::cout << "\n W Self-Energies: \n";
         std::cout << "\n W: " << W_se_all << '\n';
         std::cout << " W deriv: " << W_se_deriv_all << '\n';

         std::cout << "\n Photon-Z Self-Energies: \n";
         std::cout << "\n BSM: " << PZ_se_all << '\n';
         std::cout << "\n SM:  " << PZ_se_sm << '\n';

         // Build-In Self-Energies
         std::cout << "\n Build-in SM Self-Energies: " << '\n';
         std::cout << "\n Higgs Self-Energies: \n";

         std::cout << "\n Higgs SM: " << higgs_se_sm << '\n';
         std::cout << " Higgs SM deriv: " << higgs_se_deriv_sm << '\n';

         std::cout << "\n Fermion Self-Energies: \n";

         std::cout << "\n Fermion SM PL: " << lepton_se_sm_PL << '\n';
         std::cout << " Fermion SM PR: " << lepton_se_sm_PR << '\n';
         std::cout << " Fermion SM 1: " << lepton_1_sm << '\n';

         std::cout << "\n Fermion SM PL deriv: " << lepton_se_deriv_sm_PL << '\n';
         std::cout << " Fermion SM PR deriv: " << lepton_se_deriv_sm_PR << '\n';
         std::cout << " Fermion SM 1 deriv: " << lepton_se_deriv_sm_1 << '\n';

         std::cout << "\n Z Self-Energies: \n";
         std::cout << "\n Z SM: " << Z_se_sm << '\n';
         std::cout << " Z SM deriv: " << Z_se_deriv_sm << '\n';

         std::cout << "\n Photon Self-Energies: \n";
         std::cout << "\n P SM: " << photon_se_sm << '\n';
         std::cout << " P SM deriv: " << photon_se_deriv_sm << '\n';
         
         std::cout << "\n W Self-Energies: \n";
         std::cout << "\n W SM: " << W_se_sm << '\n';
         std::cout << " W SM deriv: " << W_se_deriv_sm << '\n';

         // Renormalization Constants
         std::cout << "\n SM Renormalization Constants" << '\n';
         std::cout << "\n Higgs Field Renormalization: " << delta_Z_higgs_sm << '\n';
         std::cout << " Higgs Mass Renormalization: " << delta_Zm_higgs_sm << '\n';

         std::cout << "\n Fermion PL Field Renormalization: " << delta_ZL_lepton_sm << '\n';
         std::cout << " Fermion PR Field Renormalization: " << delta_ZR_lepton_sm << '\n';
         std::cout << " Fermion Mass Renormalization: " << delta_Zm_lepton_sm << '\n';

         std::cout << "\n Z Field Renormalization: " << delta_Z_Z_sm << '\n';
         std::cout << " Z Mass Renormalization: " << delta_Zm_Z_sm << '\n';

         std::cout << "\n W Field Renormalization: " << delta_Z_W_sm << '\n';
         std::cout << " W Mass Renormalization: " << delta_Zm_W_sm << '\n';

         std::cout << "\n Charge Renormalization: " << delta_e_sm << '\n';
         std::cout << "\n ev Renormaliztion: " << delta_Z_vev_sm << '\n';
         
         std::cout << "\n Decoupling Renormalization Constants" << '\n';
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


         std::cout << "\n Counter Term 1 left form-factor: " << ct_FL1 << '\n';
         std::cout << " Counter Term 1 right form-factor: " << ct_FR1 << '\n';

         std::cout << "\n Counter Term 2 left form-factor: " << ct_FL2 << '\n';
         std::cout << " Counter Term 2 right form-factor: " << ct_FR2 << '\n';

         std::cout << "\n counter term left: " << ct_FL << '\n';
         std::cout << " counter term right: " << ct_FR << '\n';

         std::cout << "\n form factor left tree: " << amp_tree.form_factor_left << '\n';
         std::cout << " form factor right tree: " << amp_tree.form_factor_right << '\n';

         std::cout << "\n form factor left 1loop SM: " << amp_1l_sm.form_factor_left << '\n';
         std::cout << " form factor right 1loop SM: " << amp_1l_sm.form_factor_right << '\n';

         std::cout << "\n form factor left 1loop: " << amp_1l.form_factor_left << '\n';
         std::cout << " form factor right 1loop: " << amp_1l.form_factor_right << '\n';

         std::cout << "\n form factor left dec: " << amp_1l_dec.form_factor_left - ct_FL << '\n';
         std::cout << " form factor right dec: " << amp_1l_dec.form_factor_right - ct_FR << '\n';

         std::cout << "\n 1l* conj(tree): " << amplitude_interference(amp_1l, amp_tree) << '\n';
         std::cout << " |tree|^2: " << amp_tree.square() << '\n';
         std::cout << " |1l|^2: " << amp_1l.square() << '\n';
         std::cout << " |A|^2: " << amp_tree.square() + amp_1l.square() + 2.0 * std::real(amplitude_interference(amp_1l, amp_tree)) << '\n';

         std::cout << "\n Decay Width tree: " << flux * ps * amp_tree.square() << '\n';
         std::cout << " Decay Width MSbar: " << res << '\n';
      }

      res += flux * ps * 2.*std::real(amplitude_interference(amp_tree, amp_counterterm));

      if (indexOut1.at(0) == 1 && indexOut2.at(0) == 1) {
         std::cout << " Decay Width Dec: " << res << '\n';
      }
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
