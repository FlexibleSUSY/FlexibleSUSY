
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_CMSSM_effective_couplings

#include <boost/test/unit_test.hpp>

#include "effective_couplings.hpp"
#include "CMSSM_mass_eigenstates.hpp"
#include "CMSSM_effective_couplings.hpp"

#include "effective_couplings.hpp"
#include "ew_input.hpp"
#include "linalg2.hpp"
#include "wrappers.hpp"

using namespace flexiblesusy;
using namespace effective_couplings;

class MSSM_loop_decays {
public:
   MSSM_loop_decays(const CMSSM_mass_eigenstates& model_,
                    const softsusy::QedQcd& qedqcd_,
                    const Physical_input& physical_inputs_)
      : model(model_)
      , qedqcd(qedqcd_)
      , physical_inputs(physical_inputs_)
      {}

   void do_rg_improve(bool flag) { rg_improve = flag; }
   void set_model(const CMSSM_mass_eigenstates& model_) { model = model_; }

   std::complex<double> get_eff_CphhVPVP(int);
   std::complex<double> get_eff_CphhVGVG(int);
   std::complex<double> get_eff_CpAhVPVP();
   std::complex<double> get_eff_CpAhVGVG();

private:
   CMSSM_mass_eigenstates model{};
   softsusy::QedQcd qedqcd{};
   Physical_input physical_inputs{};
   bool rg_improve{true};

   double get_higgs_mixing_angle() const;
   void run_SM_strong_coupling_to(double);

   Eigen::Matrix<std::complex<double>,2,2> get_up_type_squark_couplings(
      int, int) const;
   Eigen::Matrix<std::complex<double>,2,2> get_down_type_squark_couplings(
      int, int) const;
   Eigen::Matrix<std::complex<double>,2,2> get_slepton_couplings(
      int, int) const;

   std::complex<double> get_hhVPVP_VWm_contribution(int) const;
   std::complex<double> get_hhVPVP_Hpm_contribution(int) const;
   std::complex<double> get_hhVPVP_Cha_contribution(int) const;
   std::complex<double> get_hhVPVP_Su_contribution(int, int) const;
   std::complex<double> get_hhVPVP_Sd_contribution(int, int) const;
   std::complex<double> get_hhVPVP_Se_contribution(int, int) const;
   std::complex<double> get_hhVPVP_Fu_contribution(int, int) const;
   std::complex<double> get_hhVPVP_Fd_contribution(int, int) const;
   std::complex<double> get_hhVPVP_Fe_contribution(int, int) const;
   std::complex<double> get_hhVGVG_Su_contribution(int, int) const;
   std::complex<double> get_hhVGVG_Sd_contribution(int, int) const;
   std::complex<double> get_hhVGVG_Fu_contribution(int, int) const;
   std::complex<double> get_hhVGVG_Fd_contribution(int, int) const;
   std::complex<double> get_AhVPVP_Cha_contribution() const;
   std::complex<double> get_AhVPVP_Fu_contribution(int) const;
   std::complex<double> get_AhVPVP_Fd_contribution(int) const;
   std::complex<double> get_AhVPVP_Fe_contribution(int) const;
   std::complex<double> get_AhVGVG_Fu_contribution(int) const;
   std::complex<double> get_AhVGVG_Fd_contribution(int) const;
};

double MSSM_loop_decays::get_higgs_mixing_angle() const
{
   const auto ZH = model.get_physical().ZH;
   return ArcTan(ZH(1,1) / ZH(1,0));
}

void MSSM_loop_decays::run_SM_strong_coupling_to(double scale)
{
   standard_model::Standard_model sm;

   sm.set_loops(2);
   sm.set_thresholds(2);
   sm.set_physical_input(physical_inputs);

   sm.initialise_from_input(qedqcd);

   sm.run_to(scale);

   model.set_g3(sm.get_g3());
}

Eigen::Matrix<std::complex<double>,2,2> MSSM_loop_decays::get_up_type_squark_couplings(
   int idx, int gen) const
{
   const double Qf = 2. / 3.;
   const double I3Lf = 0.5;
   const double gY = Sqrt(3. / 5.) * model.get_g1();
   const double g2 = model.get_g2();
   const double gbar = Sqrt(Sqr(g2) + Sqr(gY));
   const double sw = gY / gbar;
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double v = Sqrt(Sqr(vd) + Sqr(vu));
   const double beta = ArcTan(vu / vd);

   const double vf = vu;
   const double yf = model.get_Yu(gen, gen);
   const double Tf = model.get_TYu(gen, gen);
   const double mu = model.get_Mu();
   const double mf = yf * vf / Sqrt(2.);
   const double mz = 0.5 * gbar * v;
   const double alpha = get_higgs_mixing_angle();

   Eigen::Matrix<std::complex<double>,2,2> Cff;
   if (idx == 0) {
      const double s1q = Cos(alpha) / Sin(beta);
      const double s2q = Sin(alpha) / Sin(beta);
      Cff(0,0) = -(I3Lf - Qf* Sqr(sw)) * Sqr(mz) * Sin(beta + alpha)
         + Sqr(mf) * s1q;
      Cff(0,1) = 0.5 * vf * (Tf * s1q + yf * mu * s2q) / Sqrt(2.);
      Cff(1,0) = Cff(0,1);
      Cff(1,1) = -Qf * Sqr(sw) * Sqr(mz) * Sin(beta + alpha)
         + Sqr(mf) * s1q;
   } else {
      const double r1q = Sin(alpha) / Sin(beta);
      const double r2q = -Cos(alpha) / Sin(beta);
      Cff(0,0) = (I3Lf - Qf * Sqr(sw)) * Sqr(mz) * Cos(beta + alpha)
         + Sqr(mf) * r1q;
      Cff(0,1) = 0.5 * vf * (Tf * r1q + yf * mu * r2q) / Sqrt(2.);
      Cff(1,0) = Cff(0,1);
      Cff(1,1) = Qf * Sqr(sw) * Sqr(mz) * Cos(beta + alpha)
         + Sqr(mf) * r1q;
   }

   return Cff;
}

Eigen::Matrix<std::complex<double>,2,2>  MSSM_loop_decays::get_down_type_squark_couplings(
   int idx, int gen) const
{
   const double Qf = -1. / 3.;
   const double I3Lf = -0.5;
   const double gY = Sqrt(3. / 5.) * model.get_g1();
   const double g2 = model.get_g2();
   const double gbar = Sqrt(Sqr(g2) + Sqr(gY));
   const double sw = gY / gbar;
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double v = Sqrt(Sqr(vd) + Sqr(vu));
   const double beta = ArcTan(vu / vd);

   const double vf = vd;
   const double yf = model.get_Yd(gen, gen);
   const double Tf = model.get_TYd(gen, gen);
   const double mu = model.get_Mu();
   const double mf = yf * vf / Sqrt(2.);
   const double mz = 0.5 * gbar * v;
   const double alpha = get_higgs_mixing_angle();

   Eigen::Matrix<std::complex<double>,2,2> Cff;
   if (idx == 0) {
      const double s1q = -Sin(alpha) / Cos(beta);
      const double s2q = -Cos(alpha) / Cos(beta);
      Cff(0,0) = -(I3Lf - Qf * Sqr(sw)) * Sqr(mz) * Sin(beta + alpha)
         + Sqr(mf) * s1q;
      Cff(0,1) = 0.5 * vf * (Tf * s1q + yf * mu * s2q) / Sqrt(2.);
      Cff(1,0) = Cff(0,1);
      Cff(1,1) = -Qf * Sqr(sw) * Sqr(mz) * Sin(beta + alpha)
         + Sqr(mf) * s1q;
   } else {
      const double r1q = Cos(alpha) / Cos(beta);
      const double r2q = -Sin(alpha) / Cos(beta);
      Cff(0,0) = (I3Lf - Qf * Sqr(sw)) * Sqr(mz) * Cos(beta + alpha)
         + Sqr(mf) * r1q;
      Cff(0,1) = 0.5 * vf * (Tf * r1q + yf * mu * r2q) / Sqrt(2.);
      Cff(1,0) = Cff(0,1);
      Cff(1,1) = Qf * Sqr(sw) * Sqr(mz) * Cos(beta + alpha)
         + Sqr(mf) * r1q;
   }

   return Cff;
}

Eigen::Matrix<std::complex<double>,2,2>  MSSM_loop_decays::get_slepton_couplings(
   int idx, int gen) const
{
   const double Qf = -1.;
   const double I3Lf = -0.5;
   const double gY = Sqrt(3. / 5.) * model.get_g1();
   const double g2 = model.get_g2();
   const double gbar = Sqrt(Sqr(g2) + Sqr(gY));
   const double sw = gY / gbar;
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double v = Sqrt(Sqr(vd) + Sqr(vu));
   const double beta = ArcTan(vu / vd);

   const double vf = vd;
   const double yf = model.get_Ye(gen, gen);
   const double Tf = model.get_TYe(gen, gen);
   const double mu = model.get_Mu();
   const double mf = yf * vf / Sqrt(2.);
   const double mz = 0.5 * gbar * v;
   const double alpha = get_higgs_mixing_angle();

   Eigen::Matrix<std::complex<double>,2,2> Cff;
   if (idx == 0) {
      const double s1q = -Sin(alpha) / Cos(beta);
      const double s2q = -Cos(alpha) / Cos(beta);
      Cff(0,0) = -(I3Lf - Qf* Sqr(sw)) * Sqr(mz) * Sin(beta + alpha)
         + Sqr(mf) * s1q;
      Cff(0,1) = 0.5 * vf * (Tf * s1q + yf * mu * s2q) / Sqrt(2.);
      Cff(1,0) = Cff(0,1);
      Cff(1,1) = -Qf * Sqr(sw) * Sqr(mz) * Sin(beta + alpha)
         + Sqr(mf) * s1q;
   } else {
      const double r1q = Cos(alpha) / Cos(beta);
      const double r2q = -Sin(alpha) / Cos(beta);
      Cff(0,0) = (I3Lf - Qf * Sqr(sw)) * Sqr(mz) * Cos(beta + alpha)
         + Sqr(mf) * r1q;
      Cff(0,1) = 0.5 * vf * (Tf * r1q + yf * mu * r2q) / Sqrt(2.);
      Cff(1,0) = Cff(0,1);
      Cff(1,1) = Qf * Sqr(sw) * Sqr(mz) * Cos(beta + alpha)
         + Sqr(mf) * r1q;
   }

   return Cff;
}

std::complex<double> MSSM_loop_decays::get_hhVPVP_VWm_contribution(
   int idx) const
{
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double beta = ArcTan(vu / vd);
   const double alpha = get_higgs_mixing_angle();

   std::complex<double> coupling;
   if (idx == 0) {
      coupling = Sin(beta - alpha);
   } else {
      coupling = Cos(beta - alpha);
   }

   const double MVWm = model.get_MVWm();
   const auto Mhh = model.get_physical().Mhh;
   const double tau = 0.25 * Sqr(Mhh(idx)) / Sqr(MVWm);

   return coupling * AS1(tau);
}

std::complex<double> MSSM_loop_decays::get_hhVPVP_Hpm_contribution(
   int idx) const
{
   const double g1 = model.get_g1();
   const double g2 = model.get_g2();
   const double cw = g2 / Sqrt(Sqr(g2) + 0.6 * Sqr(g1));
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double beta = ArcTan(vu / vd);
   const double alpha = get_higgs_mixing_angle();

   std::complex<double> coupling;
   if (idx == 0) {
      coupling = Cos(2. * beta) * Sin(beta + alpha)
         + 2. * Sqr(cw) * Sin(beta - alpha);
   } else {
      coupling = -Cos(2. * beta) * Cos(beta + alpha)
         + 2. * Sqr(cw) * Cos(beta - alpha);
   }

   const double MHpm = model.get_MHpm(1);
   const auto Mhh = model.get_physical().Mhh;
   const double tau = 0.25 * Sqr(Mhh(idx)) / Sqr(MHpm);

   return 0.125 * (Sqr(g2) + 0.6 * Sqr(g1)) * (Sqr(vd) + Sqr(vu))
      * coupling * AS0(tau) / Sqr(MHpm);
}

std::complex<double> MSSM_loop_decays::get_hhVPVP_Cha_contribution(
   int idx) const
{
   const double g2 = model.get_g2();
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double v = Sqrt(Sqr(vd) + Sqr(vu));

   const double alpha = get_higgs_mixing_angle();

   const auto U = (model.get_UM()).conjugate();
   const auto V = model.get_UP();
   std::complex<double> coupling_light;
   std::complex<double> coupling_heavy;
   if (idx == 0) {
      coupling_light = (-Sin(alpha) * V(0,0) * U(0,1)
                        + Cos(alpha) * V(0,1) * U(0,0)) / Sqrt(2.);
      coupling_heavy = (-Sin(alpha) * V(1,0) * U(1,1)
                        + Cos(alpha) * V(1,1) * U(1,0)) / Sqrt(2.);
   } else {
      coupling_light = (Cos(alpha) * V(0,0) * U(0,1)
                        + Sin(alpha) * V(0,1) * U(0,0)) / Sqrt(2.);
      coupling_heavy = (Cos(alpha) * V(1,0) * U(1,1)
                        + Sin(alpha) * V(1,1) * U(1,0)) / Sqrt(2.);
   }

   const auto MCha = model.get_MCha();
   const auto Mhh = model.get_physical().Mhh;
   const double tau_light = 0.25 * Sqr(Mhh(idx)) / Sqr(MCha(0));
   const double tau_heavy = 0.25 * Sqr(Mhh(idx)) / Sqr(MCha(1));

   return g2 * v * coupling_light * AS12(tau_light) / MCha(0)
      + g2 * v * coupling_heavy * AS12(tau_heavy) / MCha(1);
}

std::complex<double> MSSM_loop_decays::get_hhVPVP_Su_contribution(
   int idx, int gen) const
{
   const int Nc = 3;
   const double Qf = 2. / 3.;
   const double I3Lf = 0.5;
   const double gY = Sqrt(3. / 5.) * model.get_g1();
   const double g2 = model.get_g2();
   const double gbar = Sqrt(Sqr(g2) + Sqr(gY));
   const double sw = gY / gbar;
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double v = Sqrt(Sqr(vd) + Sqr(vu));
   const double beta = ArcTan(vu / vd);

   const double vf = vu;
   const double vo = vd;
   const double mLf2 = model.get_mq2(gen, gen);
   const double mRf2 = model.get_mu2(gen, gen);
   const double yf = model.get_Yu(gen, gen);
   const double Tf = model.get_TYu(gen, gen);
   const double mu = model.get_Mu();
   const double mf = yf * vf / Sqrt(2.);
   const double mix = (Tf * vf - yf * mu * vo) / Sqrt(2.);
   const double mz = 0.5 * gbar * v;

   Eigen::Matrix<std::complex<double>,2,2> Cff
      = get_up_type_squark_couplings(idx, gen);

   const double mLL2 = mLf2 + (I3Lf - Qf * Sqr(sw)) * Sqr(mz)
      * Cos(2. * beta);
   const double mRR2 = mRf2 + Qf * Sqr(sw) * Sqr(mz) * Cos(2. * beta);

   Eigen::Matrix<double,2,2> mass_matrix;
   mass_matrix << mLL2 + Sqr(mf), mix, mix, mRR2 + Sqr(mf);

   Eigen::Array<double,2,1> masses;
   Eigen::Matrix<std::complex<double>,2,2> mixings;
   fs_diagonalize_symmetric(mass_matrix, masses, mixings);

   const double msf12 = masses(0);
   const double msf22 = masses(1);

   const Eigen::Matrix<std::complex<double>,2,2> gff =
      mixings * Cff * mixings.transpose();

   const auto Mhh = model.get_physical().Mhh;
   const double tau_light = 0.25 * Sqr(Mhh(idx)) / msf12;
   const double tau_heavy = 0.25 * Sqr(Mhh(idx)) / msf22;

   return Nc * Sqr(Qf) * gff(0,0) * AS0(tau_light) / msf12
      + Nc * Sqr(Qf) * gff(1,1) * AS0(tau_heavy) / msf22;
}

std::complex<double> MSSM_loop_decays::get_hhVPVP_Sd_contribution(
   int idx, int gen) const
{
   const int Nc = 3;
   const double Qf = -1. / 3.;
   const double I3Lf = -0.5;
   const double gY = Sqrt(3. / 5.) * model.get_g1();
   const double g2 = model.get_g2();
   const double gbar = Sqrt(Sqr(g2) + Sqr(gY));
   const double sw = gY / gbar;
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double v = Sqrt(Sqr(vd) + Sqr(vu));
   const double beta = ArcTan(vu / vd);

   const double vf = vd;
   const double vo = vu;
   const double mLf2 = model.get_mq2(gen, gen);
   const double mRf2 = model.get_md2(gen, gen);
   const double yf = model.get_Yd(gen, gen);
   const double Tf = model.get_TYd(gen, gen);
   const double mu = model.get_Mu();
   const double mf = yf * vf / Sqrt(2.);
   const double mix = (Tf * vf - yf * mu * vo) / Sqrt(2.);
   const double mz = 0.5 * gbar * v;

   Eigen::Matrix<std::complex<double>,2,2> Cff
      = get_down_type_squark_couplings(idx, gen);

   const double mLL2 = mLf2 + (I3Lf - Qf * Sqr(sw)) * Sqr(mz)
      * Cos(2. * beta);
   const double mRR2 = mRf2 + Qf * Sqr(sw) * Sqr(mz) * Cos(2. * beta);

   Eigen::Matrix<double,2,2> mass_matrix;
   mass_matrix << mLL2 + Sqr(mf), mix, mix, mRR2 + Sqr(mf);

   Eigen::Array<double,2,1> masses;
   Eigen::Matrix<std::complex<double>,2,2> mixings;
   fs_diagonalize_symmetric(mass_matrix, masses, mixings);

   const double msf12 = masses(0);
   const double msf22 = masses(1);

   const Eigen::Matrix<std::complex<double>,2,2> gff =
      mixings * Cff * mixings.transpose();

   const auto Mhh = model.get_physical().Mhh;
   const double tau_light = 0.25 * Sqr(Mhh(idx)) / msf12;
   const double tau_heavy = 0.25 * Sqr(Mhh(idx)) / msf22;

   return Nc * Sqr(Qf) * gff(0,0) * AS0(tau_light) / msf12
      + Nc * Sqr(Qf) * gff(1,1) * AS0(tau_heavy) / msf22;
}

std::complex<double> MSSM_loop_decays::get_hhVPVP_Se_contribution(
   int idx, int gen) const
{
   const int Nc = 1;
   const double Qf = -1.;
   const double I3Lf = -0.5;
   const double gY = Sqrt(3. / 5.) * model.get_g1();
   const double g2 = model.get_g2();
   const double gbar = Sqrt(Sqr(g2) + Sqr(gY));
   const double sw = gY / gbar;
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double v = Sqrt(Sqr(vd) + Sqr(vu));
   const double beta = ArcTan(vu / vd);

   const double vf = vd;
   const double vo = vu;
   const double mLf2 = model.get_ml2(gen, gen);
   const double mRf2 = model.get_me2(gen, gen);
   const double yf = model.get_Ye(gen, gen);
   const double Tf = model.get_TYe(gen, gen);
   const double mu = model.get_Mu();
   const double mf = yf * vf / Sqrt(2.);
   const double mix = (Tf * vf - yf * mu * vo) / Sqrt(2.);
   const double mz = 0.5 * gbar * v;

   Eigen::Matrix<std::complex<double>,2,2> Cff
      = get_slepton_couplings(idx, gen);

   const double mLL2 = mLf2 + (I3Lf - Qf * Sqr(sw)) * Sqr(mz)
      * Cos(2. * beta);
   const double mRR2 = mRf2 + Qf * Sqr(sw) * Sqr(mz) * Cos(2. * beta);

   Eigen::Matrix<double,2,2> mass_matrix;
   mass_matrix << mLL2 + Sqr(mf), mix, mix, mRR2 + Sqr(mf);

   Eigen::Array<double,2,1> masses;
   Eigen::Matrix<std::complex<double>,2,2> mixings;
   fs_diagonalize_symmetric(mass_matrix, masses, mixings);

   const double msf12 = masses(0);
   const double msf22 = masses(1);

   const Eigen::Matrix<std::complex<double>,2,2> gff =
      mixings * Cff * mixings.transpose();

   const auto Mhh = model.get_physical().Mhh;
   const double tau_light = 0.25 * Sqr(Mhh(idx)) / msf12;
   const double tau_heavy = 0.25 * Sqr(Mhh(idx)) / msf22;

   return Nc * Sqr(Qf) * gff(0,0) * AS0(tau_light) / msf12
      + Nc * Sqr(Qf) * gff(1,1) * AS0(tau_heavy) / msf22;
}

std::complex<double> MSSM_loop_decays::get_hhVPVP_Fu_contribution(
   int idx, int gen) const
{
   const int Nc = 3;
   const double Qf = 2. / 3.;
   const double alpha = get_higgs_mixing_angle();
   const double beta = ArcTan(model.get_vu() / model.get_vd());

   std::complex<double> coupling;
   if (idx == 0) {
      coupling = Cos(alpha) / Sin(beta);
   } else {
      coupling = Sin(alpha) / Sin(beta);
   }

   const auto Mhh = model.get_physical().Mhh;
   const auto MFu = model.get_MFu();
   const double tau = 0.25 * Sqr(Mhh(idx)) / Sqr(MFu(gen));

   return Nc * Sqr(Qf) * coupling * AS12(tau);
}

std::complex<double> MSSM_loop_decays::get_hhVPVP_Fd_contribution(
   int idx, int gen) const
{
   const int Nc = 3;
   const double Qf = -1. / 3.;
   const double alpha = get_higgs_mixing_angle();
   const double beta = ArcTan(model.get_vu() / model.get_vd());

   std::complex<double> coupling;
   if (idx == 0) {
      coupling = -Sin(alpha) / Cos(beta);
   } else {
      coupling = Cos(alpha) / Cos(beta);
   }

   const auto Mhh = model.get_physical().Mhh;
   const auto MFd = model.get_MFd();
   const double tau = 0.25 * Sqr(Mhh(idx)) / Sqr(MFd(gen));

   return Nc * Sqr(Qf) * coupling * AS12(tau);
}

std::complex<double> MSSM_loop_decays::get_hhVPVP_Fe_contribution(
   int idx, int gen) const
{
   const double Qf = -1.;
   const double alpha = get_higgs_mixing_angle();
   const double beta = ArcTan(model.get_vu() / model.get_vd());

   std::complex<double> coupling;
   if (idx == 0) {
      coupling = -Sin(alpha) / Cos(beta);
   } else {
      coupling = Cos(alpha) / Cos(beta);
   }

   const auto Mhh = model.get_physical().Mhh;
   const auto MFe = model.get_MFe();
   const double tau = 0.25 * Sqr(Mhh(idx)) / Sqr(MFe(gen));

   return Sqr(Qf) * coupling * AS12(tau);
}

std::complex<double> MSSM_loop_decays::get_hhVGVG_Su_contribution(
   int idx, int gen) const
{
   const double Qf = 2. / 3.;
   const double I3Lf = 0.5;
   const double gY = Sqrt(3. / 5.) * model.get_g1();
   const double g2 = model.get_g2();
   const double gbar = Sqrt(Sqr(g2) + Sqr(gY));
   const double sw = gY / gbar;
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double v = Sqrt(Sqr(vd) + Sqr(vu));
   const double beta = ArcTan(vu / vd);

   const double vf = vu;
   const double vo = vd;
   const double mLf2 = model.get_mq2(gen, gen);
   const double mRf2 = model.get_mu2(gen, gen);
   const double yf = model.get_Yu(gen, gen);
   const double Tf = model.get_TYu(gen, gen);
   const double mu = model.get_Mu();
   const double mf = yf * vf / Sqrt(2.);
   const double mix = (Tf * vf - yf * mu * vo) / Sqrt(2.);
   const double mz = 0.5 * gbar * v;

   Eigen::Matrix<std::complex<double>,2,2> Cff
      = get_up_type_squark_couplings(idx, gen);

   const double mLL2 = mLf2 + (I3Lf - Qf * Sqr(sw)) * Sqr(mz)
      * Cos(2. * beta);
   const double mRR2 = mRf2 + Qf * Sqr(sw) * Sqr(mz) * Cos(2. * beta);

   Eigen::Matrix<double,2,2> mass_matrix;
   mass_matrix << mLL2 + Sqr(mf), mix, mix, mRR2 + Sqr(mf);

   Eigen::Array<double,2,1> masses;
   Eigen::Matrix<std::complex<double>,2,2> mixings;
   fs_diagonalize_symmetric(mass_matrix, masses, mixings);

   const double msf12 = masses(0);
   const double msf22 = masses(1);

   const Eigen::Matrix<std::complex<double>,2,2> gff =
      mixings * Cff * mixings.transpose();

   const auto Mhh = model.get_physical().Mhh;
   const double tau_light = 0.25 * Sqr(Mhh(idx)) / msf12;
   const double tau_heavy = 0.25 * Sqr(Mhh(idx)) / msf22;

   return gff(0,0) * AS0(tau_light) / msf12
      + gff(1,1) * AS0(tau_heavy) / msf22;
}

std::complex<double> MSSM_loop_decays::get_hhVGVG_Sd_contribution(
   int idx, int gen) const
{
   const double Qf = -1. / 3.;
   const double I3Lf = -0.5;
   const double gY = Sqrt(3. / 5.) * model.get_g1();
   const double g2 = model.get_g2();
   const double gbar = Sqrt(Sqr(g2) + Sqr(gY));
   const double sw = gY / gbar;
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double v = Sqrt(Sqr(vd) + Sqr(vu));
   const double beta = ArcTan(vu / vd);

   const double vf = vd;
   const double vo = vu;
   const double mLf2 = model.get_mq2(gen, gen);
   const double mRf2 = model.get_md2(gen, gen);
   const double yf = model.get_Yd(gen, gen);
   const double Tf = model.get_TYd(gen, gen);
   const double mu = model.get_Mu();
   const double mf = yf * vf / Sqrt(2.);
   const double mix = (Tf * vf - yf * mu * vo) / Sqrt(2.);
   const double mz = 0.5 * gbar * v;

   Eigen::Matrix<std::complex<double>,2,2> Cff
      = get_down_type_squark_couplings(idx, gen);

   const double mLL2 = mLf2 + (I3Lf - Qf * Sqr(sw)) * Sqr(mz)
      * Cos(2. * beta);
   const double mRR2 = mRf2 + Qf * Sqr(sw) * Sqr(mz) * Cos(2. * beta);

   Eigen::Matrix<double,2,2> mass_matrix;
   mass_matrix << mLL2 + Sqr(mf), mix, mix, mRR2 + Sqr(mf);

   Eigen::Array<double,2,1> masses;
   Eigen::Matrix<std::complex<double>,2,2> mixings;
   fs_diagonalize_symmetric(mass_matrix, masses, mixings);

   const double msf12 = masses(0);
   const double msf22 = masses(1);

   const Eigen::Matrix<std::complex<double>,2,2> gff =
      mixings * Cff * mixings.transpose();

   const auto Mhh = model.get_physical().Mhh;
   const double tau_light = 0.25 * Sqr(Mhh(idx)) / msf12;
   const double tau_heavy = 0.25 * Sqr(Mhh(idx)) / msf22;

   return gff(0,0) * AS0(tau_light) / msf12
      + gff(1,1) * AS0(tau_heavy) / msf22;
}

std::complex<double> MSSM_loop_decays::get_hhVGVG_Fu_contribution(
   int idx, int gen) const
{
   const double alpha = get_higgs_mixing_angle();
   const double beta = ArcTan(model.get_vu() / model.get_vd());

   std::complex<double> coupling;
   if (idx == 0) {
      coupling = Cos(alpha) / Sin(beta);
   } else {
      coupling = Sin(alpha) / Sin(beta);
   }

   const auto Mhh = model.get_physical().Mhh;
   const auto MFu = model.get_MFu();
   const double tau = 0.25 * Sqr(Mhh(idx)) / Sqr(MFu(gen));

   return coupling * AS12(tau);
}

std::complex<double> MSSM_loop_decays::get_hhVGVG_Fd_contribution(
   int idx, int gen) const
{
   const double alpha = get_higgs_mixing_angle();
   const double beta = ArcTan(model.get_vu() / model.get_vd());

   std::complex<double> coupling;
   if (idx == 0) {
      coupling = -Sin(alpha) / Cos(beta);
   } else {
      coupling = Cos(alpha) / Cos(beta);
   }

   const auto Mhh = model.get_physical().Mhh;
   const auto MFd = model.get_MFd();
   const double tau = 0.25 * Sqr(Mhh(idx)) / Sqr(MFd(gen));

   return coupling * AS12(tau);
}

std::complex<double> MSSM_loop_decays::get_AhVPVP_Cha_contribution() const
{
   const double g2 = model.get_g2();
   const double vd = model.get_vd();
   const double vu = model.get_vu();
   const double v = Sqrt(Sqr(vd) + Sqr(vu));

   const auto ZA = model.get_physical().ZA;
   const auto U = (model.get_UM()).conjugate();
   const auto V = model.get_UP();
   const std::complex<double> coupling_light
      = (-ZA(1,0) * V(0,0) * U(0,1)
         - ZA(1,1) * V(0,1) * U(0,0)) / Sqrt(2.)
      * std::complex<double>(0.,1.);
   std::complex<double> coupling_heavy
      = (-ZA(1,0) * V(1,0) * U(1,1)
         - ZA(1,1) * V(1,1) * U(1,0)) / Sqrt(2.)
      * std::complex<double>(0.,1.);

   const auto MCha = model.get_MCha();
   const auto MAh = model.get_physical().MAh;
   const double tau_light = 0.25 * Sqr(MAh(1)) / Sqr(MCha(0));
   const double tau_heavy = 0.25 * Sqr(MAh(1)) / Sqr(MCha(1));

   return g2 * v * coupling_light * 2. * AP12(tau_light) / MCha(0)
      + g2 * v * coupling_heavy * 2. * AP12(tau_heavy) / MCha(1);
}

std::complex<double> MSSM_loop_decays::get_AhVPVP_Fu_contribution(
   int gen) const
{
   const int Nc = 3;
   const double Qf = 2. / 3.;
   const auto ZA = model.get_physical().ZA;
   const double beta = ArcTan(model.get_vu() / model.get_vd());
   const std::complex<double> coupling = ZA(1,1) / Sin(beta)
      * std::complex<double>(0., 1.);

   const auto MAh = model.get_physical().MAh;
   const auto MFu = model.get_MFu();
   const double tau = 0.25 * Sqr(MAh(1)) / Sqr(MFu(gen));

   return Nc * Sqr(Qf) * coupling * 2. * AP12(tau);
}

std::complex<double> MSSM_loop_decays::get_AhVPVP_Fd_contribution(
   int gen) const
{
   const int Nc = 3;
   const double Qf = -1. / 3.;
   const auto ZA = model.get_physical().ZA;
   const double beta = ArcTan(model.get_vu() / model.get_vd());
   const std::complex<double> coupling = ZA(1,0) / Cos(beta)
      * std::complex<double>(0., 1.);

   const auto MAh = model.get_physical().MAh;
   const auto MFd = model.get_MFd();
   const double tau = 0.25 * Sqr(MAh(1)) / Sqr(MFd(gen));

   return Nc * Sqr(Qf) * coupling * 2. * AP12(tau);
}

std::complex<double> MSSM_loop_decays::get_AhVPVP_Fe_contribution(
   int gen) const
{
   const double Qf = -1.;
   const auto ZA = model.get_physical().ZA;
   const double beta = ArcTan(model.get_vu() / model.get_vd());
   const std::complex<double> coupling = ZA(1,0) / Cos(beta)
      * std::complex<double>(0., 1.);

   const auto MAh = model.get_physical().MAh;
   const auto MFe = model.get_MFe();
   const double tau = 0.25 * Sqr(MAh(1)) / Sqr(MFe(gen));

   return Sqr(Qf) * coupling * 2. * AP12(tau);
}

std::complex<double> MSSM_loop_decays::get_AhVGVG_Fu_contribution(
   int gen) const
{
   const auto ZA = model.get_physical().ZA;
   const double beta = ArcTan(model.get_vu() / model.get_vd());
   const std::complex<double> coupling = ZA(1,1) / Sin(beta)
      * std::complex<double>(0., 1.);

   const auto MAh = model.get_physical().MAh;
   const auto MFu = model.get_MFu();
   const double tau = 0.25 * Sqr(MAh(1)) / Sqr(MFu(gen));

   return coupling * 2. * AP12(tau);
}

std::complex<double> MSSM_loop_decays::get_AhVGVG_Fd_contribution(
   int gen) const
{
   const auto ZA = model.get_physical().ZA;
   const double beta = ArcTan(model.get_vu() / model.get_vd());
   const std::complex<double> coupling = ZA(1,0) / Cos(beta)
      * std::complex<double>(0., 1.);

   const auto MAh = model.get_physical().MAh;
   const auto MFd = model.get_MFd();
   const double tau = 0.25 * Sqr(MAh(1)) / Sqr(MFd(gen));

   return coupling * 2. * AP12(tau);
}

std::complex<double> MSSM_loop_decays::get_eff_CphhVPVP(int idx)
{
   const auto Mhh = model.get_physical().Mhh;
   if (rg_improve && model.get_scale() > Mhh(idx)) {
      model.run_to(Mhh(idx));
   }
   model.calculate_DRbar_masses();

   std::complex<double> result = get_hhVPVP_VWm_contribution(idx);
   result += get_hhVPVP_Hpm_contribution(idx);
   result += get_hhVPVP_Cha_contribution(idx);
   for (int gen = 0; gen < 3; ++gen) {
      result += get_hhVPVP_Su_contribution(idx, gen);
      result += get_hhVPVP_Sd_contribution(idx, gen);
      result += get_hhVPVP_Se_contribution(idx, gen);
      result += get_hhVPVP_Fu_contribution(idx, gen);
      result += get_hhVPVP_Fd_contribution(idx, gen);
      result += get_hhVPVP_Fe_contribution(idx, gen);
   }
   const double prefactor = 0.1892681907127351
      * physical_inputs.get(Physical_input::alpha_em_0)
      * Sqrt(qedqcd.displayFermiConstant());

   return prefactor * result;
}

std::complex<double> MSSM_loop_decays::get_eff_CphhVGVG(int idx)
{
   const double old_g3 = model.get_g3();
   const auto Mhh = model.get_physical().Mhh;
   if (rg_improve && model.get_scale() > Mhh(idx)) {
      model.run_to(Mhh(idx));
   }
   run_SM_strong_coupling_to(Mhh(idx));
   model.calculate_DRbar_masses();

   std::complex<double> result = std::complex<double>(0.,0.);
   for (int gen = 0; gen < 3; ++gen) {
      result += get_hhVGVG_Su_contribution(idx, gen);
      result += get_hhVGVG_Sd_contribution(idx, gen);
      result += get_hhVGVG_Fu_contribution(idx, gen);
      result += get_hhVGVG_Fd_contribution(idx, gen);
   }
   result *= 0.75;

   const double alpha_s = 0.07957747154594767*Sqr(model.get_g3());
   const double prefactor = 0.12617879380849006 * alpha_s
      * Sqrt(qedqcd.displayFermiConstant());

   model.set_g3(old_g3);

   return prefactor * result;
}

std::complex<double> MSSM_loop_decays::get_eff_CpAhVPVP()
{
   const auto MAh = model.get_physical().MAh;
   if (rg_improve && model.get_scale() > MAh(1)) {
      model.run_to(MAh(1));
   }
   model.calculate_DRbar_masses();

   std::complex<double> result = get_AhVPVP_Cha_contribution();
   for (int gen = 0; gen < 3; ++gen) {
      result += get_AhVPVP_Fu_contribution(gen);
      result += get_AhVPVP_Fd_contribution(gen);
      result += get_AhVPVP_Fe_contribution(gen);
   }
   const double prefactor = 0.1892681907127351
      * physical_inputs.get(Physical_input::alpha_em_0)
      * Sqrt(qedqcd.displayFermiConstant());

   return prefactor * result;
}

std::complex<double> MSSM_loop_decays::get_eff_CpAhVGVG()
{
   const double old_g3 = model.get_g3();
   const auto MAh = model.get_physical().MAh;
   if (rg_improve && model.get_scale() > MAh(1)) {
      model.run_to(MAh(1));
   }
   run_SM_strong_coupling_to(MAh(1));
   model.calculate_DRbar_masses();

   std::complex<double> result = std::complex<double>(0.,0.);
   for (int gen = 0; gen < 3; ++gen) {
      result += get_AhVGVG_Fu_contribution(gen);
      result += get_AhVGVG_Fd_contribution(gen);
   }
   result *= 0.75;

   const double alpha_s = 0.07957747154594767*Sqr(model.get_g3());
   const double prefactor = 0.12617879380849006 * alpha_s
      * Sqrt(qedqcd.displayFermiConstant());

   model.set_g3(old_g3);

   return prefactor * result;
}
void set_test_model_parameters(CMSSM_mass_eigenstates& model,
                               CMSSM_input_parameters& input,
                               const softsusy::QedQcd& qedqcd)
{
   input.m0 = 125.;
   input.m12 = 500.;
   input.Azero = 0.;
   input.TanBeta = 10.0;

   const double scale = 864.196350;

   const double g1 = Sqrt(5. / 3.) * 0.362388869;
   const double g2 = 0.643064079;
   const double g3 = 1.06292111;
   const double Mu = 623.398163;
   const double BMu = 53774.7440;
   const double vd = 25.0988057;
   const double vu = 242.826920;

   Eigen::Matrix<double,3,3> Yu(Eigen::Matrix<double,3,3>::Zero());
   Eigen::Matrix<double,3,3> Yd(Eigen::Matrix<double,3,3>::Zero());
   Eigen::Matrix<double,3,3> Ye(Eigen::Matrix<double,3,3>::Zero());

   Yu(0,0) = 7.32901347e-6;
   Yu(1,1) = 3.35149154e-3;
   Yu(2,2) = qedqcd.displayPoleMt()* Sqrt(2.) / vu;

   Yd(0,0) = 1.41183630e-4;
   Yd(1,1) = 3.09118035e-3;
   Yd(2,2) = 1.33272717e-1;

   Ye(0,0) = 2.89560942e-5;
   Ye(1,1) = 5.98720208e-3;
   Ye(2,2) = 1.00698298e-1;

   Eigen::Matrix<double,3,3> AYu(Eigen::Matrix<double,3,3>::Zero());
   Eigen::Matrix<double,3,3> AYd(Eigen::Matrix<double,3,3>::Zero());
   Eigen::Matrix<double,3,3> AYe(Eigen::Matrix<double,3,3>::Zero());

   AYu(0,0) = -1130.62590;
   AYu(1,1) = -1130.62089;
   AYu(2,2) = -872.799086;

   AYd(0,0) = -1384.23875;
   AYd(1,1) = -1384.23410;
   AYd(2,2) = -1293.85753;

   AYe(0,0) = -299.168792;
   AYe(1,1) = -299.163003;
   AYe(2,2) = -297.531417;

   Eigen::Matrix<double,3,3> TYu = (Yu.array() * AYu.array()).matrix();
   Eigen::Matrix<double,3,3> TYd = (Yd.array() * AYd.array()).matrix();
   Eigen::Matrix<double,3,3> TYe = (Ye.array() * AYe.array()).matrix();

   Eigen::Matrix<double,3,3> mq2(Eigen::Matrix<double,3,3>::Zero());
   Eigen::Matrix<double,3,3> ml2(Eigen::Matrix<double,3,3>::Zero());
   Eigen::Matrix<double,3,3> md2(Eigen::Matrix<double,3,3>::Zero());
   Eigen::Matrix<double,3,3> mu2(Eigen::Matrix<double,3,3>::Zero());
   Eigen::Matrix<double,3,3> me2(Eigen::Matrix<double,3,3>::Zero());

   mq2(0,0) = 1.01310643e6;
   mq2(1,1) = 1.01310128e6;
   mq2(2,2) = 8.61071428e5;

   ml2(0,0) = 1.24941646e5;
   ml2(1,1) = 1.24939130e5;
   ml2(2,2) = 1.24230841e5;

   md2(0,0) = 9.27220763e5;
   md2(1,1) = 9.27215601e5;
   md2(2,2) = 9.18044968e5;

   mu2(0,0) = 9.36437027e5;
   mu2(1,1) = 9.36431757e5;
   mu2(2,2) = 6.34856659e5;

   me2(0,0) = 4.93863949e4;
   me2(1,1) = 4.93812597e4;
   me2(2,2) = 4.79356909e4;

   const double mHd2 = 1.09303956e5;
   const double mHu2 = -3.76666907e5;
   const double MassB = 209.081139;
   const double MassWB = 387.961620;
   const double MassG = 1114.07969;

   model.set_input_parameters(input);
   model.set_scale(scale);
   model.set_loops(2);
   model.set_thresholds(3);
   model.set_g1(g1);
   model.set_g2(g2);
   model.set_g3(g3);
   model.set_Yu(Yu);
   model.set_Yd(Yd);
   model.set_Ye(Ye);
   model.set_Mu(Mu);
   model.set_vd(vd);
   model.set_vu(vu);
   model.set_BMu(BMu);
   model.set_MassB(MassB);
   model.set_MassG(MassWB);
   model.set_MassWB(MassG);
   model.set_mq2(mq2);
   model.set_ml2(ml2);
   model.set_md2(md2);
   model.set_mu2(mu2);
   model.set_me2(me2);
   model.set_mHd2(mHd2);
   model.set_mHu2(mHu2);
   model.set_TYu(TYu);
   model.set_TYd(TYd);
   model.set_TYe(TYe);

   model.calculate_DRbar_masses();
   model.solve_ewsb();
   model.calculate_spectrum();
}

BOOST_AUTO_TEST_CASE( test_LO_scalar_diphoton_couplings )
{
   softsusy::QedQcd qedqcd;
   Physical_input physical_inputs;

   CMSSM_input_parameters input;
   CMSSM_mass_eigenstates model;
   set_test_model_parameters(model, input, qedqcd);

   CMSSM_effective_couplings eff_cp(model, qedqcd, physical_inputs);
   eff_cp.do_run_couplings(false);
   eff_cp.do_include_qcd_corrections(false);
   eff_cp.calculate_effective_couplings();

   const std::complex<double> obtained_cphhVPVP_0 = eff_cp.get_eff_CphhVPVP(0);
   const std::complex<double> obtained_cphhVPVP_1 = eff_cp.get_eff_CphhVPVP(1);

   MSSM_loop_decays mssm_tester(model, qedqcd, physical_inputs);
   mssm_tester.do_rg_improve(false);

   const std::complex<double> expected_cphhVPVP_0
      = mssm_tester.get_eff_CphhVPVP(0);
   const std::complex<double> expected_cphhVPVP_1
      = mssm_tester.get_eff_CphhVPVP(1);

   BOOST_CHECK_CLOSE_FRACTION(-Re(expected_cphhVPVP_0),
                              Re(obtained_cphhVPVP_0), 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(-Im(expected_cphhVPVP_0),
                              Im(obtained_cphhVPVP_0), 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(-Re(expected_cphhVPVP_1),
                              Re(obtained_cphhVPVP_1), 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(-Im(expected_cphhVPVP_1),
                              Im(obtained_cphhVPVP_1), 1.0e-10);
}

BOOST_AUTO_TEST_CASE( test_LO_scalar_digluon_couplings )
{
   softsusy::QedQcd qedqcd;
   Physical_input physical_inputs;

   CMSSM_input_parameters input;
   CMSSM_mass_eigenstates model;
   set_test_model_parameters(model, input, qedqcd);

   CMSSM_effective_couplings eff_cp(model, qedqcd, physical_inputs);
   eff_cp.do_run_couplings(false);
   eff_cp.do_include_qcd_corrections(false);
   eff_cp.calculate_effective_couplings();

   const std::complex<double> obtained_cphhVGVG_0 = eff_cp.get_eff_CphhVGVG(0);
   const std::complex<double> obtained_cphhVGVG_1 = eff_cp.get_eff_CphhVGVG(1);

   MSSM_loop_decays mssm_tester(model, qedqcd, physical_inputs);
   mssm_tester.do_rg_improve(false);

   const std::complex<double> expected_cphhVGVG_0
      = mssm_tester.get_eff_CphhVGVG(0);
   const std::complex<double> expected_cphhVGVG_1
      = mssm_tester.get_eff_CphhVGVG(1);

   BOOST_CHECK_CLOSE_FRACTION(-Re(expected_cphhVGVG_0),
                              Re(obtained_cphhVGVG_0), 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(-Im(expected_cphhVGVG_0),
                              Im(obtained_cphhVGVG_0), 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(-Re(expected_cphhVGVG_1),
                              Re(obtained_cphhVGVG_1), 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(-Im(expected_cphhVGVG_1),
                              Im(obtained_cphhVGVG_1), 1.0e-10);
}

BOOST_AUTO_TEST_CASE( test_LO_pseudoscalar_diphoton_couplings )
{
   softsusy::QedQcd qedqcd;
   Physical_input physical_inputs;

   CMSSM_input_parameters input;
   CMSSM_mass_eigenstates model;
   set_test_model_parameters(model, input, qedqcd);

   CMSSM_effective_couplings eff_cp(model, qedqcd, physical_inputs);
   eff_cp.do_run_couplings(false);
   eff_cp.do_include_qcd_corrections(false);
   eff_cp.calculate_effective_couplings();

   const std::complex<double> obtained_cpAhVPVP_1 = eff_cp.get_eff_CpAhVPVP(1);

   MSSM_loop_decays mssm_tester(model, qedqcd, physical_inputs);
   mssm_tester.do_rg_improve(false);

   const std::complex<double> expected_cpAhVPVP_1
      = mssm_tester.get_eff_CpAhVPVP();

   BOOST_CHECK_CLOSE_FRACTION(-Re(expected_cpAhVPVP_1),
                              Re(obtained_cpAhVPVP_1), 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(-Im(expected_cpAhVPVP_1),
                              Im(obtained_cpAhVPVP_1), 1.0e-10);
}

BOOST_AUTO_TEST_CASE( test_LO_pseudoscalar_digluon_couplings )
{
   softsusy::QedQcd qedqcd;
   Physical_input physical_inputs;

   CMSSM_input_parameters input;
   CMSSM_mass_eigenstates model;
   set_test_model_parameters(model, input, qedqcd);

   CMSSM_effective_couplings eff_cp(model, qedqcd, physical_inputs);
   eff_cp.do_run_couplings(false);
   eff_cp.do_include_qcd_corrections(false);
   eff_cp.calculate_effective_couplings();

   const std::complex<double> obtained_cpAhVGVG_1 = eff_cp.get_eff_CpAhVGVG(1);

   MSSM_loop_decays mssm_tester(model, qedqcd, physical_inputs);
   mssm_tester.do_rg_improve(false);

   const std::complex<double> expected_cpAhVGVG_1
      = mssm_tester.get_eff_CpAhVGVG();

   BOOST_CHECK_CLOSE_FRACTION(-Re(expected_cpAhVGVG_1),
                              Re(obtained_cpAhVGVG_1), 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(-Im(expected_cpAhVGVG_1),
                              Im(obtained_cpAhVGVG_1), 1.0e-10);
}
