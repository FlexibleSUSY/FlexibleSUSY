// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

/** \file lowe.h
   - Project:     SOFTSUSY
   - Author:      Ben Allanach, Alexander Voigt
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/
   \brief QedQcd object contains Standard Model quark and lepton
   masses. It integrates them using 3 loop qcd x 1 loop qed effective theory.
*/

#ifndef LOWE_H
#define LOWE_H

#include "betafunction.hpp"
#include "ckm.hpp"
#include "pmns.hpp"
#include <array>
#include <iosfwd>
#include <Eigen/Core>

namespace softsusy {

/// used to give order of quark masses stored
enum EMass {mUp=1, mCharm, mTop, mDown, mStrange, mBottom, mElectron,
           mMuon, mTau};
/// order of gauge couplings stored in QedQcd
enum EGauge {ALPHA=1, ALPHAS};

enum QedQcd_input_parmeters : int {
   alpha_em_MSbar_at_MZ,
   alpha_s_MSbar_at_MZ,
   GFermi,
   MZ_pole, MW_pole,
   Mv1_pole, Mv2_pole, Mv3_pole,
   Me_pole, Mm_pole, Mtau_pole,
   mu_2GeV, ms_2GeV, Mt_pole,
   md_2GeV, mc_mc, mb_mb,
   CKM_theta_12, CKM_theta_13, CKM_theta_23, CKM_delta,
   PMNS_theta_12, PMNS_theta_13, PMNS_theta_23, PMNS_delta, PMNS_alpha_1, PMNS_alpha_2,
   NUMBER_OF_LOW_ENERGY_INPUT_PARAMETERS
};

extern const std::array<std::string, NUMBER_OF_LOW_ENERGY_INPUT_PARAMETERS> QedQcd_input_parmeter_names;

/// Quark and lepton masses and gauge couplings in QEDxQCD effective theory
class QedQcd: public flexiblesusy::Beta_function
{
public:
  using Input_t = Eigen::Array<double,NUMBER_OF_LOW_ENERGY_INPUT_PARAMETERS,1>;

private:
  Eigen::Array<double,2,1> a{Eigen::Array<double,2,1>::Zero()};  ///< gauge couplings
  Eigen::Array<double,9,1> mf{Eigen::Array<double,9,1>::Zero()}; ///< fermion running masses
  Input_t input{Input_t::Zero()}; ///< SLHA input parmeters
  double mbPole;    ///< pole masses of third family quarks

  double qedBeta() const;   ///< QED beta function
  double qcdBeta() const;   ///< QCD beta function
  Eigen::Array<double,9,1> massBeta() const; ///< beta functions of masses
  void runto_safe(double, double eps = -1.0); ///< throws if non-perturbative error occurs

  int flavours(double) const;  /// returns number of active flavours

  /// calculates pole bottom mass given alpha_s(Mb)^{MSbar} from running b mass
  double extractPoleMb(double asMb);

public:
  QedQcd();
  QedQcd(const QedQcd&) = default;
  QedQcd(QedQcd&&) = default;
  QedQcd& operator=(const QedQcd&) = default;
  QedQcd& operator=(QedQcd&&) = default;
  virtual ~QedQcd() = default;

  // Beta_function interface
  virtual Eigen::ArrayXd get() const override;
  virtual void set(const Eigen::ArrayXd&) override;
  virtual Eigen::ArrayXd beta() const override;

  void setPoleMt(double mt) { input(Mt_pole) = mt; } ///< set pole top mass
  void setPoleMb(double mb) { mbPole = mb; } ///< set pole bottom mass
  void setPoleMtau(double mtau) { input(Mtau_pole) = mtau; } ///< set pole tau mass
  void setPoleMmuon(double m) { input(Mm_pole) = m; } ///< set pole muon mass
  void setPoleMel(double m) { input(Me_pole) = m;  } ///< set pole electron mass
  void setMbMb(double mb)   { input(mb_mb) = mb;   } ///< set mb(mb)
  void setMcMc(double mc)   { input(mc_mc) = mc;   } ///< set mc(mc)
  void setMu2GeV(double mu) { input(mu_2GeV) = mu; } ///< set mu(2 GeV)
  void setMd2GeV(double md) { input(md_2GeV) = md; } ///< set md(2 GeV)
  void setMs2GeV(double ms) { input(ms_2GeV) = ms; } ///< set ms(2 GeV)
  void setPoleMW(double mw) { input(MW_pole) = mw; } ///< set W boson pole mass
  void setPoleMZ(double mz) { input(MZ_pole) = mz; } ///< set Z boson pole mass
  /// sets running quark masses
  void setMasses(const Eigen::Array<double,9,1>& m) { mf = m; }
  /// sets a running quark mass
  void setMass(EMass mno, double m) { mf(mno - 1) = m; }
  /// sets a neutrino pole mass
  void setNeutrinoPoleMass(int i, double m) { input(Mv1_pole + i - 1) = m; }
  /// sets QED and QCD structure constants
  void setAlphas(const Eigen::Array<double,2,1>& o) { a = o; }
  /// sets QED or QCD structure constant
  void setAlpha(EGauge ai, double ap) { a(ai - 1) = ap; }
  /// set input value of alpha_em(MZ)
  void setAlphaEmInput(double a) { input(alpha_em_MSbar_at_MZ) = a; }
  /// set input value of alpha_s(MZ)
  void setAlphaSInput(double a) { input(alpha_s_MSbar_at_MZ) = a; }
  /// sets CKM parameters (in the MS-bar scheme at MZ)
  void setCKM(const flexiblesusy::CKM_parameters&);
  /// sets PMNS parameters (in the MS-bar scheme at MZ)
  void setPMNS(const flexiblesusy::PMNS_parameters&);
  /// sets Fermi constant
  void setFermiConstant(double gf) { input(GFermi) = gf; }
  /// sets all input parameters
  void set_input(const Input_t& i) { input = i; }

  /// Displays input parameters
  Input_t displayInput() const { return input; }
  /// Display pole top mass
  double displayPoleMt() const { return input(Mt_pole); }
  /// Display pole tau mass
  double displayPoleMtau() const { return input(Mtau_pole); }
  /// Display pole muon mass
  double displayPoleMmuon() const { return input(Mm_pole); }
  /// Display pole electron mass
  double displayPoleMel() const { return input(Me_pole); }
  /// Returns bottom "pole" mass
  double displayPoleMb() const { return mbPole; }
  /// Returns W boson pole mass
  double displayPoleMW() const { return input(MW_pole); }
  /// Returns Z boson pole mass
  double displayPoleMZ() const { return input(MZ_pole); }
  /// Returns a vector of running fermion masses
  auto displayMass() const -> decltype(mf) { return mf; }
  /// Returns a single running mass
  double displayMass(EMass mno) const { return mf(mno - 1); }
  /// Returns a single neutrino pole mass
  double displayNeutrinoPoleMass(int i) const { return input(Mv1_pole + i - 1); }
  /// Returns a single gauge structure constant
  double displayAlpha(EGauge ai) const { return a(ai - 1); }
  /// Returns gauge structure constants
  auto displayAlphas() const -> decltype(a) { return a; }
  /// Returns input value alpha_em(MZ)
  double displayAlphaEmInput() const { return input(alpha_em_MSbar_at_MZ); }
  /// Returns input value alpha_s(MZ)
  double displayAlphaSInput() const { return input(alpha_s_MSbar_at_MZ); }
  /// Returns Fermi constant
  double displayFermiConstant() const { return input(GFermi); }
  /// returns vector of all input parameters
  Input_t display_input() const { return input; }
  /// returns vector of all parameter names
  static std::array<std::string, NUMBER_OF_LOW_ENERGY_INPUT_PARAMETERS> display_input_parameter_names();
  /// Returns mb(mb) MSbar
  double displayMbMb() const { return input(mb_mb); }
  /// Returns mc(mc) MSbar
  double displayMcMc() const { return input(mc_mc); }
  /// Returns mu(2 GeV)
  double displayMu2GeV() const { return input(mu_2GeV); }
  /// Returns md(2 GeV)
  double displayMd2GeV() const { return input(md_2GeV); }
  /// Returns ms(2 GeV)
  double displayMs2GeV() const { return input(ms_2GeV); }
  /// returns CKM parameters
  flexiblesusy::CKM_parameters displayCKM() const;
  /// Returns real CKM matrix
  Eigen::Matrix<double,3,3> get_real_ckm() const { return displayCKM().get_real_ckm(); }
  /// Returns complex CKM matrix
  Eigen::Matrix<std::complex<double>,3,3> get_complex_ckm() const { return displayCKM().get_complex_ckm(); }
  /// returns PMNS parameters
  flexiblesusy::PMNS_parameters displayPMNS() const;
  /// Returns real PMNS matrix
  Eigen::Matrix<double,3,3> get_real_pmns() const { return displayPMNS().get_real_pmns(); }
  /// Returns complex PMNS matrix
  Eigen::Matrix<std::complex<double>,3,3> get_complex_pmns() const { return displayPMNS().get_complex_pmns(); }

  /// Evolves object to MZ
  void toMz();
  /// Evolves object to given scale.
  void to(double scale, double precision_goal = 1e-5, int max_iterations = 20);
  /// guess coupling constants {alpha_1, alpha_2, alpha_3} in SM(5)
  Eigen::Array<double,3,1> guess_alpha_SM5(double scale) const;
};

/// Formatted output from QedQcd object
std::ostream & operator<<(std::ostream &, const QedQcd &);

bool operator ==(const QedQcd&, const QedQcd&);

} // namespace softsusy

#endif
