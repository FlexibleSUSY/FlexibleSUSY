
/** \file softpars.h
   - Project:     SOFTSUSY
   - Author:      Ben Allanach
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/
   \brief Soft SUSY breaking parameters

*/

#ifndef NMSSMSOFTPARS_H
#define NMSSMSOFTPARS_H

#include <cmath>
#include "def.h"
#include "linalg.h"
#include "utils.h"
#include "numerics.h"
#include "nmssmsusy.h"
#include "softpars.h"

namespace softsusy {

/// PA: Number of parameters contained in the NMSSM RGEs (5 more soft paraemeters than the MSSM).
const static int numSoftParsNmssm = 83 + numNMssmPars;

/// Soft SUSY breaking parameters and beta functions.
class SoftParsNmssm: public SoftPars<NmssmSusy, nmsBrevity>
{
private:
  // NMSSM Trilinear soft terms not present in MSSM
  double alambda, akappa;
  // NMSSM soft breaking masses \f$ m_S^2 |S|^2 + \frac{1}{2}m{'2}_S S^2 + \xiS S + h.c. \f$ repectively.
  double mSsq, mSpsq, xiS;

public:
  using SoftPars<NmssmSusy, nmsBrevity>::displayGaugino;
  using SoftPars<NmssmSusy, nmsBrevity>::setM32;
  using SoftPars<NmssmSusy, nmsBrevity>::setM3Squared;
  using SoftPars<NmssmSusy, nmsBrevity>::displayM3Squared;
  using SoftPars<NmssmSusy, nmsBrevity>::setMh1Squared;
  using SoftPars<NmssmSusy, nmsBrevity>::displayMh1Squared;
  using SoftPars<NmssmSusy, nmsBrevity>::setMh2Squared;
  using SoftPars<NmssmSusy, nmsBrevity>::displayMh2Squared;
  using SoftPars<NmssmSusy, nmsBrevity>::yTildes;
  using SoftPars<NmssmSusy, nmsBrevity>::u1R_PQflip;
  using SoftPars<NmssmSusy, nmsBrevity>::universalGauginos;
  using SoftPars<NmssmSusy, nmsBrevity>::displayTrilinear;
  using SoftPars<NmssmSusy, nmsBrevity>::setTrilinearElement;
  using SoftPars<NmssmSusy, nmsBrevity>::setTrilinearMatrix;
  using SoftPars<NmssmSusy, nmsBrevity>::displaySoftMassSquared;
  using SoftPars<NmssmSusy, nmsBrevity>::displaySoftA;

  /// Default constructor fills object with zeroes
  SoftParsNmssm();
  /// Constructor fills SUSY conserving parts with another object, all
  /// SUSY breaking parameters set to zero
  SoftParsNmssm(const NmssmSusy &);
  /// Constructor sets all parameters equal to those in another object
  SoftParsNmssm(const SoftParsNmssm &);
  /// Constructor sets all parameters equal to those in the base class
  SoftParsNmssm(const SoftPars<NmssmSusy, nmsBrevity> &);
  /// Sets all parameters equal to those in another object
  const SoftParsNmssm & operator=(const SoftParsNmssm & s);
  /// Constructor sets RPC SUSY parameters to s, gaugino masses to mG,
  /// trilinears to aU, aD, aE for au, ad, ae
  /// trilnears respectively,  \f$m_Q^2\f$=mQl, \f$m_U^2\f$=mUr,
  /// \f$m_D^2\f$=mDr, \f$m_L^2\f$=mLl, \f$m_E^2\f$=mEr, \f$ m_3^2\f$=m3sq,
  /// \f$m_{H_1}^2\f$=mH1sq, \f$m_{H_2}^2\f$=mH2sq, mu parameter, number of
  /// loops=l, and threshold parameter=t
  SoftParsNmssm(const NmssmSusy & s, const DoubleVector & mG, const
	       DoubleMatrix & aU, const DoubleMatrix & aD, const DoubleMatrix
	       & aE, const double & alambda, const double & akappa,
               const DoubleMatrix & mQl, const DoubleMatrix & mUr, const
	       DoubleMatrix & mDr, const DoubleMatrix & mLl, const
	       DoubleMatrix & mEr, double m3sq, double mH1sq, double mH2sq,
               double mSsq, double mSpsq, double xiS,
	       double mGravitino, double mu, int l, int t);

  /// Returns whole object as a const
  const SoftParsNmssm & displaySoftPars() const;

  //PA: Return trilinear soft mass $a_\lambda S H_u H_d$
  double displayTrialambda() const;
  //PA: Return trillinear soft mass $a_\kappa S S S$
  double displayTriakappa() const;
  //PA: Return trilinear soft mass in "SUGRA style" $ \lambda A_\lambda S H_u H_d$
  double displaySoftAlambda() const;
  //PA: Return trillinear soft mass in "SUGRA style" $\kappa A_\kappa S S S$
  double displaySoftAkappa() const;

  double displayMsSquared() const; ///< Return \f$m_{S}^2\f$=mSsq
  double displayMspSquared() const; ///< Return \f$m_{S'}^2\f$=mSpsq
  double displayXiS() const; ///< Return xiS i.e \f$xiS S\f$
  /// Return contents of object in a vector: for RG evolution
  virtual const DoubleVector display() const;

  /// Sets whole thing equal to another object
  void setSoftPars(SoftParsNmssm const &);
  //PA: Set trilinear SUSY breaking parameter alambda
  void setTrialambda(double al);
  //PA: Set trilinear SUSY breaking parameter akappa
  void setTriakappa(double ak);
  void setMsSquared(double); //PA: <sets soft scalar mass \f$m_S^2 |S|^2\f$
  void setMspSquared(double); //PA: <sets soft scalar mass \f$m_S^2 S^2\f$
  void setXiS(double); //PA: <sets soft breaking \f$ xiS S\f$
  /// Sets total set of RGE parameters equal to elements of a vector
  void set(const DoubleVector &);
  //  void setSusy(const NmssmSusy &);

  /// Returns double vector containing numerical beta functions of parameters
  DoubleVector beta() const;
  /// Returns numerical beta functions of parameters
  SoftParsNmssm beta2() const;
  /// Returns derivatives of anomalous dimensions of fields with respect to
  /// renormalisation scale in MSSM for: RH leptons, LH leptons, LH quarks, RH
  /// up quarks, RH down quarks, H1 and H2 respectively
  void anomalousDeriv(DoubleMatrix & gEE, DoubleMatrix & gLL,
		      DoubleMatrix & gQQ, DoubleMatrix & gUU,
		      DoubleMatrix & gDD,
		      double & gH1H1, double & gH2H2, double & gSS) const;

  //adds amsb contribution to soft masses for NMSSM.
  void addAmsb(double maux);
  /// Reads in universal boundary conditions at the current scale: m0, M1/2, A
  void standardSugra(double m0, double m12, double a0);
  /// Reads in universal boundary conditions at the current scale: m0, M1/2, A, mS, Al and Ak
  void standardsemiSugra(double m0,  double m12, double a0, double Al, 
			 double Ak, double mS=0.0);
  /// Sets all flavour-diagonal SUSY breaking scalar masses to m0
  void universalScalars(double m0);
  /// Sets singet soft mass to mS and all other flavour-diagonal SUSY breaking scalar masses to m0
  void semiuniversalScalars(double m0);
  /// Sets all SUSY breaking trilinear couplings to a0
  void universalTrilinears(double a0);
  /// Sets singlet trilinears to ak  and al and all SUSY breaking trilinear couplings to a0
  void semiuniversalTrilinears(double a0, double al, double ak);
};

/// Formatted ouput of whole object
std::ostream & operator <<(std::ostream &left, const SoftParsNmssm &s);

inline SoftParsNmssm::SoftParsNmssm()
  : SoftPars<NmssmSusy, nmsBrevity>(), alambda(0.0), akappa(0.0),
  mSsq(0.0), mSpsq(0.0), xiS(0.0) {
  setPars(numSoftParsNmssm);
  setMu(0.0);
  setLoops(0);
  setThresholds(0);
}

inline SoftParsNmssm::SoftParsNmssm(const SoftParsNmssm & s)
  : SoftPars<NmssmSusy, nmsBrevity>(s),
  alambda(s.displayTrialambda()), akappa(s.displayTriakappa()),
  mSsq(s.displayMsSquared()), mSpsq(s.displayMspSquared()),
  xiS(s.displayXiS()) {
  setPars(numSoftParsNmssm);
  setMu(s.displayMu());
  setLoops(s.displayLoops());
  setThresholds(s.displayThresholds());
}

inline SoftParsNmssm::SoftParsNmssm(const NmssmSusy &s)
    : SoftPars<NmssmSusy, nmsBrevity>(s), alambda(0.0), akappa(0.0),
    mSsq(0.0), mSpsq(0.0), xiS(0.0) {
      setPars(numSoftParsNmssm);
      setMu(s.displayMu());
      setLoops(s.displayLoops());
      setThresholds(s.displayThresholds());
}

inline SoftParsNmssm::SoftParsNmssm(const SoftPars<NmssmSusy, nmsBrevity> & s)
    : SoftPars<NmssmSusy, nmsBrevity>(s), alambda(0.0), akappa(0.0),
    mSsq(0.0), mSpsq(0.0), xiS(0.0) {
      setPars(numSoftParsNmssm);
      setMu(s.displayMu());
      setLoops(s.displayLoops());
      setThresholds(s.displayThresholds());
}

inline SoftParsNmssm::SoftParsNmssm
(const NmssmSusy & s, const DoubleVector & mG, const
 DoubleMatrix & aU, const DoubleMatrix & aD, const DoubleMatrix & aE,
 const double & aL, const double & aK, const DoubleMatrix & mQl,
 const DoubleMatrix & mUr, const DoubleMatrix & mDr, const
 DoubleMatrix & mLl, const DoubleMatrix & mEr, double m3sqn, double mH1sq,
 double mH2sq, double mSsq, double mSpsq, double xiS,
 double mg, double mu, int l, int t)
   : SoftPars<NmssmSusy, nmsBrevity>(s, mG, aU, aD, aE, mQl, mUr, mDr, mLl, mEr, m3sqn, mH1sq, mH2sq, mg, mu, l, t)
   , alambda(aL), akappa(aK), mSsq(mSsq), mSpsq(mSpsq), xiS(xiS) {
      setPars(numSoftParsNmssm);
      setMu(mu);
      setLoops(l);
      setThresholds(t);
}

inline const SoftParsNmssm & SoftParsNmssm::displaySoftPars() const { return *this; }

inline double SoftParsNmssm::displayMsSquared() const { return mSsq; }

inline double SoftParsNmssm::displayMspSquared() const { return mSpsq; }

inline double SoftParsNmssm::displayXiS() const { return xiS; }

inline double SoftParsNmssm::displayTrialambda() const { return alambda; }

inline double SoftParsNmssm::displayTriakappa() const { return akappa; }

inline void SoftParsNmssm::setMsSquared(double f) { mSsq = f; }
inline void SoftParsNmssm::setMspSquared(double f) { mSpsq = f; }
inline void SoftParsNmssm::setXiS(double f) { xiS = f; }
inline void SoftParsNmssm::setSoftPars(SoftParsNmssm const & s) { *this = s; }

} // namespace softsusy

#endif
