
/** \file susy.h
   - Project:     SOFTSUSY 
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/

   \brief MssmSusy contains all SUSY couplings and tan beta, as well as their
   beta functions 
*/

#ifndef SUSY_H
#define SUSY_H

#include <iostream>
#include <cmath>
#include <fstream>
#include "linalg.h"
#include "lowe_legacy.h"
#include "rge.h"
#include "utils.h"

namespace softsusy {

class Tensor;

const static int numSusyPars = 33;

/// For accessing up, down and charged lepton Yukawa matrices respectively
typedef enum {YU=1, YD, YE} yukawa;

/// Contains data needed in beta function calculation to make it faster
struct sBrevity {
  /// dt=\f$Y_D^T\f$, ut=\f$Y_U^T\f$, et=\f$Y_E^T\f$, u2=\f$Y_U Y_U^T\f$,
  /// d2=\f$Y_D Y_D^T\f$, e2=\f$Y_E Y_E^T\f$, 
  /// u2t=\f$Y_U^T Y_U\f$, d2t=\f$Y_D^T Y_D\f$, e2t=\f$Y_E^T Y_E\f$
  DoubleMatrix dt, ut, et, u2, d2, e2, u2t, e2t, d2t;
  /// gsq=\f$g_i^2\f$, g3=\f$g_i^3\f$, g4=\f$g_i^4\f$
  DoubleVector  gsq, g3, g4;
  /// uuT=\f$Tr(Y_U^T Y_U)\f$, ddt=\f$Tr(Y_D^T Y_D)\f$, eet=\f$Tr(Y_E^T Y_E)\f$
  double uuT, ddT, eeT;
  /// d1=\f$Y_D\f$, u1=\f$Y_U\f$, e1=\f$Y_E\f$
  DoubleMatrix u1, d1, e1;

  /// Constructor fills sruct with zeroes by default
  sBrevity();
  /// Constructor sets struct to be equal to another
  sBrevity(const sBrevity &);
  /// Sets struct to be equal to another
  const sBrevity & operator=(const sBrevity &);
  
  /// Calculates all the data given Yukawa matrices yu, yd and ye and gauge
  /// couplings g
  void calculate(const DoubleMatrix & yu, const DoubleMatrix & yd, const
		 DoubleMatrix & ye, const DoubleVector  & g);
};

inline sBrevity::sBrevity()
  : dt(3, 3), ut(3, 3), et(3, 3), u2(3, 3), 
    d2(3, 3), e2(3, 3), u2t(3, 3), 
    e2t(3, 3), d2t(3, 3), gsq(1, 3), g3(1, 3), g4(1, 3), uuT(0.0),
    ddT(0.0), eeT(0.0), u1(3, 3), d1(3, 3), e1(3, 3)
{}

inline sBrevity::sBrevity(const sBrevity &s)
  : dt(s.dt), ut(s.ut), et(s.ut), u2(s.u2), d2(s.d2),
    e2(s.e2), u2t(s.u2t),
    e2t(s.e2t), d2t(s.d2t), gsq(s.gsq), g3(s.g3), g4(s.g4),
    uuT(s.uuT), ddT(s.ddT), eeT(s.eeT), u1(s.u1), d1(s.d1), e1(s.e1) 
{}

/// Contains all supersymmetric RPC-MSSM parameters and RGEs 
class MssmSusy: public RGE {
private:
  DoubleMatrix u, d, e; ///< Yukawa matrices for ups, downs and leptons
  DoubleVector g; ///< Gauge couplings in GUT normalisation (for g1)
  /// Bilinear Higgs superpotential parameter and ratio of Higgs VEVs,
  /// \f$ v_1/v_2 \f$ 
  double smu, tanb, hVev; 
  
public:
  MssmSusy(); ///< Constructor fills object with zeroes by default
  /// Constructor sets object to be equal to another
  MssmSusy(const MssmSusy &); 
  /// Constructor given Yukawa matrices u,d,e, gauge couplings v, mu
  /// parameter=m, tan beta=tb, renormalisation scale MU, number of loops in
  /// RG evolution l and thresholds parameter t
  MssmSusy(const DoubleMatrix & u, const DoubleMatrix & d, const
	   DoubleMatrix & e, const DoubleVector & v, double m,
	   double tb, double MU, int l, int t, double h);
  virtual ~MssmSusy() {}; ///< Default destructor

  /// sets object to be equal to another
  const MssmSusy & operator=(const MssmSusy & s);
  /// sets object to be equal to another
  void setSusy(const MssmSusy &s);

  /// Sets DRbar running Higgs vev
  void setHvev(double h);
  /// Copies Yukawa matrices and gauge couplings from s only
  void setSomePars(const MssmSusy & s);
  /// Sets one element of a Yukawa matrix
  void setYukawaElement(yukawa, int, int, double);
  /// Sets whole Yukawa matrix 
  void setYukawaMatrix(yukawa, const DoubleMatrix &);
  /// Set a single gauge coupling
  void setGaugeCoupling(int, double);
  /// Set all gauge couplings
  void setAllGauge(const DoubleVector  &);
  /// Sets superpotential mu parameter
  void setSusyMu(double);
  /// Sets all RGE parameters to elements of vector
  void set(const DoubleVector &);
  /// Sets tan beta
  void setTanb(double);
  
  /// Returns DRbar running Higgs vev
  double displayHvev() const;
  /// Returns whole object as a const
  inline const MssmSusy & displaySusy() const;
  /// Returns a single Yukawa matrix element
  double displayYukawaElement(yukawa, int, int) const;
  /// Returns a whole Yukawa matrix
  const DoubleMatrix & displayYukawaMatrix(yukawa) const;
  /// Returns a single gauge coupling
  double displayGaugeCoupling(int) const;
  /// Returns all gauge couplings
  DoubleVector  displayGauge() const;
  /// Returns superpotential mu parameter
  double displaySusyMu() const;
  /// Returns all parameters as elements of a vector
  const DoubleVector display() const;
  /// Returns tan beta
  double displayTanb() const;
  
  /// outputs object QedQcd_legacy & r valid at 1 GeV from SUSY data at mt, from
  /// diagonal elements of Yukawa couplings and Higgs VEV vev. 
  void getMasses(QedQcd_legacy & r, double vev) const;
  /// This turns diagonal Yukawa couplings at MZ into CKM mixed ones
  /// Takes diagonal quark Yukawa matrices and mixes them up
  /// according to the CKM matrix assuming:
  /// mix=2, all mixing is in down sector
  /// mix=1, all mixing is in up sector
  void quarkMixing(const DoubleMatrix & CKM, int mix);
  /// Sets diagonal Yukawa couplings according to data in QedQcd_legacy input and
  /// Higgs VEV parameter vev=\f$v_1^2+v_2^2\f$
  void setDiagYukawas(const QedQcd_legacy &, double vev);
  /// Defines mixed Yukawa matrices from data input in form of CKM matrix and
  /// r, vev. If mix=2, all mixing is in down sector
  /// mix=1, all mixing is in up sector
  void getQuarkMixedYukawas(const QedQcd_legacy & r, const DoubleMatrix &
				    CKM, int mix, double vev);
  /// Calculate beta functions of SUSY preserving parameters of RPC MSSM
  DoubleVector beta() const;
  /// Calculate beta functions of SUSY preserving parameters of RPC MSSM
  MssmSusy beta(sBrevity &) const;
  /// Outputs one-loop anomlous dimensions gii given matrix inputs.
  /// for RH leptons, LH leptons, LH quarks, RH downs, RH ups, H1 and H2
  /// respectively. Note that we use the convention (for matrices in terms of
  /// gamma's): gamma^Li_Lj = M_ij for LH fields and
  /// gamma^Rj_Ri = M_ij for RH fields (since they are really the complex
  /// conjugates of the RH fields). a should already be defined.
  void getOneLpAnom(DoubleMatrix & gEE, DoubleMatrix & gLL,
		    DoubleMatrix & gQQ, DoubleMatrix & gDD,
		    DoubleMatrix & gUU, double & gH1H1, double &
		    gH2H2, sBrevity & a) const; 
  /// Outputs two-loop anomlous dimensions gii given matrix inputs.
  /// for RH leptons, LH leptons, LH quarks, RH downs, RH ups, H1 and H2
  /// respectively. Note that we use the convention (for matrices in terms of
  /// gamma's): gamma^Li_Lj = M_ij for LH fields and
  /// gamma^Rj_Ri = M_ij for RH fields (since they are really the complex
  /// conjugates of the RH fields). a should already be defined.
  void getTwoLpAnom(DoubleMatrix & gEE, DoubleMatrix & gLL,
		    DoubleMatrix & gQQ, DoubleMatrix & gDD,
		    DoubleMatrix & gUU, double & gH1H1, double &
		    gH2H2, sBrevity & a) const; 
  /// Calculates three-loop anomalous dimensions in the dominant third family
  /// approximation and adds them
  void getThreeLpAnom(DoubleMatrix & gEE, DoubleMatrix & gLL,
                      DoubleMatrix & gQQ, DoubleMatrix & gDD,
                      DoubleMatrix & gUU, double & gH1H1, double &
                      gH2H2, sBrevity & a) const;
  /// Outputs wave function renormalisation for SUSY parameters and gauge beta
  /// functions up to 2 loops. Also calculates and outputs a.
  /// IO parameters: RH leptons, LH leptons, LH quarks, RH downs, RH ups, H1
  /// and H2 respectively. 
  /// g^Li_Lj = m_{ij} for LH fields
  /// g^Ei_Ej = m_{ji} for RH fields
  void anomalousDimension(DoubleMatrix & gEE, DoubleMatrix & gLL,
			  DoubleMatrix & gQQ, DoubleMatrix & gUU,
			  DoubleMatrix & gDD, DoubleVector & dg, 
			  double & gH1H1, double & gH2H2, 
			  sBrevity & a) const; 
  /// Rotates to quark mass basis, returning the mixing matrices defined as 
  /// \f$(Y_U)_{diag}\f$ = vul \f$Y_U\f$ vur^+,  
  /// \f$(Y_D)_{diag}\f$ = vdl \f$Y_D\f$ vdr^+  
  /// All matrices should be 3 by 3
  void diagQuarkBasis(DoubleMatrix & vdl, DoubleMatrix & vdr, 
			  DoubleMatrix & vul, DoubleMatrix & vur) const;
};

/// Formatted output
std::ostream & operator <<(std::ostream &, const MssmSusy &);

/// Outputs beta function coefficients for MSSM gauge coupling evolution in
/// arguments. 
void setBetas(DoubleMatrix &, DoubleVector  &, DoubleVector  &, DoubleVector
	       &, DoubleVector  &);

void setBetasThreeLoop(Tensor &, DoubleMatrix &, DoubleMatrix &, 
                       DoubleMatrix &, DoubleVector &, DoubleVector &, DoubleVector &,
                       DoubleVector &, DoubleVector &, DoubleVector &, DoubleVector &,  
                       DoubleVector &); 

inline const MssmSusy & MssmSusy::displaySusy() const { return *this; }

inline void MssmSusy::setGaugeCoupling(int i, double f) { g(i) = f; }

inline void MssmSusy::setAllGauge(const DoubleVector & v) { 
  if (v.displayStart() != 1 || v.displayEnd() !=3) {
    std::ostringstream ii;
    ii << 
      "Initialising SUSY params gauge function with vector NOT 1..3\n" <<
      v;
    throw ii.str();
  }
  g = v; 
}

inline double MssmSusy::displayHvev() const { return hVev; } 

inline void MssmSusy::setHvev(double h) { hVev = h; }
inline void MssmSusy::setSusyMu(double f) { smu = f; }
inline void MssmSusy::setTanb(double f) { tanb = f; }
inline DoubleVector MssmSusy::displayGauge() const { return g; }
inline double MssmSusy::displayGaugeCoupling(int i) const { 
  return g.display(i); 
}
inline double MssmSusy::displaySusyMu() const { return smu; }

} // namespace softsusy

#endif
