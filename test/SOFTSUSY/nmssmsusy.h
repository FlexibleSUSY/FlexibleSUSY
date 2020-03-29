
/** \file nmssmsusy.h
   - Project:     SOFTSUSY
   - Author:      Ben Allanach, Peter Athron, Alexander Voigt
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/

   \brief nMssmSusy contains all NMSSM SUSY couplings and tan beta, as
   well as their beta functions
*/

#ifndef NMSSMSUSY_H
#define NMSSMSUSY_H

#include "susy.h"

namespace softsusy {

// with the 5 new nmssm parameters we now have 38 in total
const static int numNMssmPars = numSusyPars + 5;

/// Contains data needed in beta function calculation to make it faster
struct nmsBrevity : public sBrevity {
  double lsq, ksq, l4, k4;

  /// Constructor fills struct with zeroes by default
  nmsBrevity();
  /// Constructor sets struct to be equal to another
  nmsBrevity(const nmsBrevity&);
  /// Sets struct to be equal to another
  const nmsBrevity & operator=(const nmsBrevity &);
  /// Fill MSSM part of the struct
  const nmsBrevity & operator=(const sBrevity &);

  /// Calculates all the data given Yukawa matrices yu, yd and ye and gauge
  /// couplings g, lambda and kappa
  void calculate(const DoubleMatrix & yu, const DoubleMatrix & yd, const
		 DoubleMatrix & ye, const DoubleVector & g,
                 const double & lam, const double & kap);
};

inline nmsBrevity::nmsBrevity()
   : sBrevity()
   , lsq(0.0), ksq(0.0), l4(0.0), k4(0.0)
{}

inline nmsBrevity::nmsBrevity(const nmsBrevity &s)
   : sBrevity(s)
   , lsq(s.lsq), ksq(s.ksq), l4(s.l4), k4(s.k4)
{}


/// Contains all supersymmetric RPC-MSSM parameters and RGEs
class NmssmSusy: public MssmSusy
{
private:
  /// new nmssm parameters, lambda, kappa appearing as superpotential
  /// terms, lambda S H_u H_d and \frac{1}{3} kappa S^3 and sVev is
  /// the vev of the singlet superfield defined by sVev = \sqrt{2}<S>.
  double lambda, kappa, sVev;
  /// new parameters appearing the general MSSM, but not the ususal
  /// form of the NMSSM which solves the mu-problem nor the MSSM.
  /// They appear in superpotential terms xiF S and \frac{1}{2} mupr
  /// S^2
  double xiF, mupr;

public:
  NmssmSusy(); ///< Constructor fills object with zeroes by default
  /// Constructor sets object to be equal to another
  NmssmSusy(const NmssmSusy &);
  /// PA: Constructor given Yukawa matrices u,d,e, gauge couplings v, mu
  /// parameter=m, tan beta=tb, lambda, kappa, mupr, xiF,
  //renormalisation scale MU, number of loops in
  /// RG evolution l and thresholds parameter t
  NmssmSusy(const DoubleMatrix & u, const DoubleMatrix & d, const
            DoubleMatrix & e, const DoubleVector & v, double m,
            double tb, double MU, int l, int t, double h, double s,
            double lambda, double kappa, double xiF, double mupr);
  
  NmssmSusy(const MssmSusy &m);
  
  virtual ~NmssmSusy(); ///< Default destructor

  /// sets object to be equal to another
  const NmssmSusy & operator=(const NmssmSusy & s);
  /// set the MssmSusy part
  const NmssmSusy & operator=(const MssmSusy & s);
  /// sets object to be equal to another
  void setSusy(const NmssmSusy &s);

  /// sets DRbar running singlet vev.
  void setSvev(double s);
  /// Copies MSSM Yukawa matrices and gauge couplings from s only
  void setSomePars(const NmssmSusy & s);
  /// sets the \lambda S H_u H_d coupling
  void setLambda(double);
  /// sets the \frac{1}{3} \kappa S^3 coupling
  void setKappa(double);
  /// sets the \frac{1}{2} mupr S^2 coupling
  void setMupr(double);
   /// sets the xiF S coupling
  void setXiF(double);
  /// Sets all RGE parameters to elements of vector
  void set(const DoubleVector &);

  /// returns DRbar running Singlet Higgs vev
  double displaySvev() const;
  /// Returns whole object as a const
  inline const NmssmSusy & displaySusy() const;
  /// returns superpotential parameter lambda
  double displayLambda() const;
  /// returns superpotential parameter lambda
  double displayKappa() const;
  /// returns mupr superpotential parameter
  double displayMupr() const;
  /// returns xiF superpotential parameter
  double displayXiF() const;
  /// Returns all parameters as elements of a vector
  const DoubleVector display() const;

  /// Calculate beta functions of SUSY preserving parameters of RPC MSSM
  DoubleVector beta() const;
  /// Calculate beta functions of SUSY preserving parameters of RPC MSSM
  NmssmSusy beta(nmsBrevity &) const;
  /// Outputs one-loop anomlous dimensions gii given matrix inputs.
  /// for RH leptons, LH leptons, LH quarks, RH downs, RH ups, H1 and H2
  /// respectively. Note that we use the convention (for matrices in terms of
  /// gamma's): gamma^Li_Lj = M_ij for LH fields and
  /// gamma^Rj_Ri = M_ij for RH fields (since they are really the complex
  /// conjugates of the RH fields). a should already be defined.
  void getOneLpAnom(DoubleMatrix & gEE, DoubleMatrix & gLL,
		    DoubleMatrix & gQQ, DoubleMatrix & gDD,
		    DoubleMatrix & gUU, double & gH1H1, double &
		    gH2H2, double & gSS, nmsBrevity & a) const;
  /// Outputs two-loop anomlous dimensions gii given matrix inputs.
  /// for RH leptons, LH leptons, LH quarks, RH downs, RH ups, H1 and H2
  /// respectively. Note that we use the convention (for matrices in terms of
  /// gamma's): gamma^Li_Lj = M_ij for LH fields and
  /// gamma^Rj_Ri = M_ij for RH fields (since they are really the complex
  /// conjugates of the RH fields). a should already be defined.
  void getTwoLpAnom(DoubleMatrix & gEE, DoubleMatrix & gLL,
		    DoubleMatrix & gQQ, DoubleMatrix & gDD,
		    DoubleMatrix & gUU, double & gH1H1, double &
		    gH2H2, double & gSS, nmsBrevity & a) const;
  /// Outputs wave function renormalisation for SUSY parameters and gauge beta
  /// functions up to 2 loops. Also calculates and outputs a.
  /// IO parameters: RH leptons, LH leptons, LH quarks, RH downs, RH ups, H1
  /// and H2 respectively.
  /// g^Li_Lj = m_{ij} for LH fields
  /// g^Ei_Ej = m_{ji} for RH fields
  void anomalousDimension(DoubleMatrix & gEE, DoubleMatrix & gLL,
			  DoubleMatrix & gQQ, DoubleMatrix & gUU,
			  DoubleMatrix & gDD, DoubleVector & dg,
			  double & gH1H1, double & gH2H2, double & gSS,
			  nmsBrevity & a) const;
};

/// Formatted output
std::ostream & operator <<(std::ostream &, const NmssmSusy &);

/// Outputs beta function coefficients for MSSM gauge coupling evolution in
/// arguments.
void nmsetBetas(DoubleMatrix &, DoubleVector  &, DoubleVector  &, DoubleVector
                &, DoubleVector  &, DoubleVector  & );

/// Outputs beta function coefficients for lambda.
void setBetaLambda(DoubleVector&);

inline const NmssmSusy & NmssmSusy::displaySusy() const { return *this; }

inline double NmssmSusy::displaySvev() const { return sVev; }

inline void NmssmSusy::setSvev(double h) { sVev = h; }
inline void NmssmSusy::setMupr(double f) { mupr = f; }
inline void NmssmSusy::setXiF(double z) { xiF = z; }

inline void NmssmSusy::setLambda(double l) { lambda = l; }
inline void NmssmSusy::setKappa(double k) { kappa = k; }

inline double NmssmSusy::displayMupr() const { return mupr; }
inline double NmssmSusy::displayXiF() const { return xiF; }
inline double NmssmSusy::displayLambda() const { return lambda; }
inline double NmssmSusy::displayKappa() const { return kappa; }

} // namespace softsusy

#endif
