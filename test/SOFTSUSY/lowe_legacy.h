
/** \file lowe_legacy.h
   - Project:     SOFTSUSY
   - Author:      Ben Allanach
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/
   \brief QedQcd_legacy object contains Standard Model quark and lepton
   masses. It integrates them using 3 loop qcd x 1 loop qed effective theory.

*/

#ifndef LOWE_LEGACY_H
#define LOWE_LEGACY_H

#include "ew_input.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include "def.h"
#include "utils.h"
#include "linalg.h"
#include "rge.h"
#include "ckm.hpp"
#include "pmns.hpp"
#include <array>
#include <Eigen/Core>

namespace softsusy {

namespace legacy {

/// used to give order of quark masses stored
enum mass {mUp=1, mCharm, mTop, mDown, mStrange, mBottom, mElectron,
           mMuon, mTau};
/// order of gauge couplings stored in QedQcd
enum leGauge {ALPHA=1, ALPHAS};

enum qedqcd_input_parmeters : int {
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

} // namespace legacy

const double MUP = flexiblesusy::Electroweak_constants::MUP; ///< default running quark mass from PDG
const double MDOWN = flexiblesusy::Electroweak_constants::MDOWN; ///< default running quark mass from PDG
const double MSTRANGE = flexiblesusy::Electroweak_constants::MSTRANGE; ///< default running quark mass from PDG
const double MCHARM = flexiblesusy::Electroweak_constants::MCHARM; ///< default running quark mass from PDG
const double MBOTTOM = flexiblesusy::Electroweak_constants::MBOTTOM; ///< default running quark mass from PDG
const double MTOP = flexiblesusy::Electroweak_constants::MTOP; ///< default running quark mass from PDG
/// default pole lepton mass from PDG
const double MELECTRON = flexiblesusy::Electroweak_constants::MELECTRON;
const double MMUON = flexiblesusy::Electroweak_constants::MMUON; ///< default pole lepton mass from PDG
const double MTAU = flexiblesusy::Electroweak_constants::MTAU; ///< default pole lepton mass from PDG
const double ALPHASMZ = flexiblesusy::Electroweak_constants::alpha3; ///< default running mass from PDG
const double ALPHAMZ = flexiblesusy::Electroweak_constants::aem; ///< default running alpha(MZ) from PDG

const double PMTOP = flexiblesusy::Electroweak_constants::PMTOP; ///< default pole mass from CDF/D0 Run II 1207.1069
const double PMBOTTOM = flexiblesusy::Electroweak_constants::PMBOTTOM; ///< default pole mass from PDG
/// default central values of CKM matrix elements from PDG 2006 in radians
const double THETA12CKM = flexiblesusy::Electroweak_constants::CKM_THETA12; ///< From Vus/Vud in global CKM fit, PDG
const double THETA13CKM = flexiblesusy::Electroweak_constants::CKM_THETA13; ///< From Vub in global CKM fit, PDG
const double THETA23CKM = flexiblesusy::Electroweak_constants::CKM_THETA23; ///< From Vcb/Vtb in global CKM fit, PDG

/// Quark and lepton masses and gauge couplings in QEDxQCD effective theory
class QedQcd_legacy: public RGE
{
private:
  DoubleVector a;   ///< gauge couplings
  DoubleVector mf;  ///< fermion running masses
  Eigen::ArrayXd input; ///< SLHA input parmeters
  double mbPole;    ///< pole masses of third family quarks
  flexiblesusy::CKM_parameters ckm; ///< CKM parameters (in the MS-bar scheme at MZ)
  flexiblesusy::PMNS_parameters pmns; ///< PMNS parameters (in the MS-bar scheme at MZ)

  DoubleVector runSMGauge(double, const DoubleVector&);
  void runto_safe(double, double); ///< throws if non-perturbative error occurs

public:
  QedQcd_legacy(); ///< Initialises with default values defined in lowe_legacy.h
  QedQcd_legacy(const QedQcd_legacy &); ///< Initialises object with another
  const QedQcd_legacy& operator=(const QedQcd_legacy & m); ///< Sets two objects equal
  virtual ~QedQcd_legacy() {};
  
  void setPoleMt(double mt) { input(legacy::Mt_pole) = mt; }; ///< set pole top mass
  void setPoleMb(double mb) { mbPole = mb; }; ///< set pole bottom mass
  void setPoleMtau(double mtau) { input(legacy::Mtau_pole) = mtau; }; ///< set pole tau mass
  void setPoleMmuon(double m) { input(legacy::Mm_pole) = m; } ///< set pole muon mass
  void setPoleMel(double m) { input(legacy::Me_pole) = m; } ///< set pole electron mass
  void setMbMb(double mb)   { input(legacy::mb_mb) = mb;   }; ///< set mb(mb)
  void setMcMc(double mc)   { input(legacy::mc_mc) = mc;   }  ///< set mc(mc)
  void setMu2GeV(double mu) { input(legacy::mu_2GeV) = mu; } ///< set mu(2 GeV)
  void setMd2GeV(double md) { input(legacy::md_2GeV) = md; } ///< set md(2 GeV)
  void setMs2GeV(double ms) { input(legacy::ms_2GeV) = ms; } ///< set ms(2 GeV)
  void setPoleMW(double mw) { input(legacy::MW_pole) = mw; } ///< set W boson pole mass
  void setPoleMZ(double mz) { input(legacy::MZ_pole) = mz; } ///< set Z boson pole mass
  /// sets running quark masses
  void setMasses(const DoubleVector& m) { mf = m; };
  /// sets a running quark mass
  void setMass(legacy::mass mno, double m) { mf(mno) = m; };
  /// sets a neutrino pole mass
  void setNeutrinoPoleMass(int i, double m) { input(legacy::Mv1_pole + i - 1) = m; }
  /// sets QED and QCD structure constants
  void setAlphas(const DoubleVector& o) { a = o; }
  /// sets QED or QCD structure constant
  void setAlpha(legacy::leGauge ai, double ap) { a(ai) = ap; }
  /// set input value of alpha_em(MZ)
  void setAlphaEmInput(double a) { input(legacy::alpha_em_MSbar_at_MZ) = a; }
  /// set input value of alpha_s(MZ)
  void setAlphaSInput(double a) { input(legacy::alpha_s_MSbar_at_MZ) = a; }
  /// sets CKM parameters (in the MS-bar scheme at MZ)
  void setCKM(const flexiblesusy::CKM_parameters& ckm_) { ckm = ckm_; }
  /// sets PMNS parameters (in the MS-bar scheme at MZ)
  void setPMNS(const flexiblesusy::PMNS_parameters& pmns_) { pmns = pmns_; }
  /// sets Fermi constant
  void setFermiConstant(double gf) { input(legacy::GFermi) = gf; }
  /// For exporting beta functions to Runge-Kutta
  void set(const DoubleVector &);
  /// sets all input parameters
  void set_input(const Eigen::ArrayXd&);

  /// Display pole top mass
  double displayPoleMt() const { return input(legacy::Mt_pole); };
  /// Display pole tau mass
  double displayPoleMtau() const { return input(legacy::Mtau_pole); };
  /// Display pole muon mass
  double displayPoleMmuon() const { return input(legacy::Mm_pole); };
  /// Display pole electron mass
  double displayPoleMel() const { return input(legacy::Me_pole); };
  /// Returns bottom "pole" mass
  double displayPoleMb() const { return mbPole; };
  /// Returns W boson pole mass
  double displayPoleMW() const { return input(legacy::MW_pole); }
  /// Returns Z boson pole mass
  double displayPoleMZ() const { return input(legacy::MZ_pole); }
  /// Returns a vector of running fermion masses
  const DoubleVector & displayMass() const { return mf; };
  /// Returns a single running mass
  double displayMass(legacy::mass mno) const { return mf.display(mno); };
  /// Returns a single neutrino pole mass
  double displayNeutrinoPoleMass(int i) const { return input(legacy::Mv1_pole + i - 1); }
  /// Returns a single gauge structure constant
  double displayAlpha(legacy::leGauge ai) const { return a.display(ai); };
  /// Returns all single gauge structure constants
  const DoubleVector& displayAlphas() const { return a; }
  /// Returns input value alpha_em(MZ)
  double displayAlphaEmInput() const { return input(legacy::alpha_em_MSbar_at_MZ); }
  /// Returns input value alpha_s(MZ)
  double displayAlphaSInput() const { return input(legacy::alpha_s_MSbar_at_MZ); }
  /// Returns Fermi constant
  double displayFermiConstant() const { return input(legacy::GFermi); }
  /// Obgligatory: returns vector of all running parameters
  const DoubleVector display() const;
  /// returns vector of all input parameters
  Eigen::ArrayXd display_input() const;
  /// returns vector of all parameter names
  static std::array<std::string, legacy::NUMBER_OF_LOW_ENERGY_INPUT_PARAMETERS> display_input_parameter_names();
  /// Returns mb(mb) MSbar
  double displayMbMb() const { return input(legacy::mb_mb); }
  /// Returns mc(mc) MSbar
  double displayMcMc() const { return input(legacy::mc_mc); }
  /// Returns mu(2 GeV)
  double displayMu2GeV() const { return input(legacy::mu_2GeV); }
  /// Returns md(2 GeV)
  double displayMd2GeV() const { return input(legacy::md_2GeV); }
  /// Returns ms(2 GeV)
  double displayMs2GeV() const { return input(legacy::ms_2GeV); }
  /// returns CKM parameters
  flexiblesusy::CKM_parameters displayCKM() const { return ckm; }
  /// Returns real CKM matrix
  Eigen::Matrix<double,3,3> get_real_ckm() const { return ckm.get_real_ckm(); }
  /// Returns complex CKM matrix
  Eigen::Matrix<std::complex<double>,3,3> get_complex_ckm() const { return ckm.get_complex_ckm(); }
  /// returns PMNS parameters
  flexiblesusy::PMNS_parameters displayPMNS() const { return pmns; }
  /// Returns real PMNS matrix
  Eigen::Matrix<double,3,3> get_real_pmns() const { return pmns.get_real_pmns(); }
  /// Returns complex PMNS matrix
  Eigen::Matrix<std::complex<double>,3,3> get_complex_pmns() const { return pmns.get_complex_pmns(); }

  int flavours(double) const;  /// returns number of active flavours

  double qedBeta() const;   ///< QED beta function
  double qcdBeta() const;   ///< QCD beta function
  void massBeta(DoubleVector &) const; ///< beta functions of masses
  /// Beta functions of both beta-functions and all MSbar masses
  DoubleVector beta() const;

  /// Does not run the masses, just gauge couplings from start to end
  void runGauge(double start, double end);
  /// calculates pole bottom mass given alpha_s(Mb)^{MSbar} from running b mass
  double extractPoleMb(double asMb);
  /// Done at pole mb: extracts running mb(polemb)
  double extractRunningMb(double asMb);
  /// calculates running bottom mass given alpha_s(Mb)^{MSbar} from pole m_b
  void calcRunningMb();
  /// Calculates the pole mass from the running mass, which should be defined
  /// at mb
  void calcPoleMb();

  /// Evolves object to running top mass
  void toMt();
  /// Evolves object to MZ
  void toMz();
  /// Evolves object to given scale.  This implementation can be called multiple times
  void to(double, double tol = 1e-5, unsigned max_iterations = 20);
  /// This will calculate the three gauge couplings of the Standard Model at
  /// the scale m2.
  /// It's a simple one-loop calculation only and no
  /// thresholds are assumed. Range of validity is electroweak to top scale.
  // alpha1 is in the GUT normalisation. sinth = sin^2 thetaW(Q) in MSbar
  // scheme
  DoubleVector  getGaugeMu(const double m2, const
		     double sinth) const;
};

/// Input numbers into the object: by file stream
std::ostream & operator <<(std::ostream &, const QedQcd_legacy &);
/// Formatted output from QedQcd_legacy object
std::istream & operator >>(std::istream &left, QedQcd_legacy &m);

/// Reads in a QedQed-type object and returns it in oneset.
/// Call with fname "" if you want it to come from standard input
/// "massIn" is an example of a data initialisation file:
void readIn(QedQcd_legacy & oneset, const char fname[80]);

inline QedQcd_legacy::QedQcd_legacy(const QedQcd_legacy &m)
   : RGE()
   , a(m.a)
   , mf(m.mf)
   , input(m.input)
   , mbPole(m.mbPole)
   , ckm(m.ckm)
   , pmns(m.pmns)
{
  setPars(11);
  setMu(m.displayMu());
  setLoops(m.displayLoops());
  setThresholds(m.displayThresholds());
}

/// Returns diagonal fermion mass matrices given input object r
// void massFermions(const QedQcd & r, DoubleMatrix & mDon,
//		  DoubleMatrix & mUpq, DoubleMatrix & mEle);

/// Returns diagonal fermion mass matrices given input object r
void massFermions(const QedQcd_legacy & r, DoubleMatrix & mDon,
		  DoubleMatrix & mUpq, DoubleMatrix & mEle);

bool operator ==(const QedQcd_legacy&, const QedQcd_legacy&);

} // namespace softsusy

#endif
