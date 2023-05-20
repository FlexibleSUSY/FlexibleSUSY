
/** \file softpars.h
   - Project:     SOFTSUSY
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/
   \brief Soft SUSY breaking parameters

*/

#ifndef SOFTPARS_H
#define SOFTPARS_H

#include <cmath>
#include "susy.h"
#include "def.h"
#include "linalg.h"
#include "utils.h"
#include "numerics.h"
#include "numerics_legacy.h"

namespace softsusy {

/// SUSY breaking soft mass squared parameters
typedef enum {mQl=1, mUr, mDr, mLl, mEr} softMasses;
/// SUSY breaking trilinear parameters
typedef enum {UA=1, DA, EA} trilinears;

/// Number of parameters contained in RGEs
const static int numSoftParsMssm = 78 + numSusyPars;

/// Soft SUSY breaking parameters and beta functions.
template<class Susy, class Brevity>
class SoftPars: public Susy
{
private:
  DoubleVector mGaugino; ///< Gaugino masses, see ::beta for definitions
  DoubleMatrix ua, da, ea; ///< Trilinear soft terms..
  /// soft mass squared matrices of \f$ m_Q^2, m_U^2, m_D^2, m_L^2, m_E^2 \f$
  /// respectively.
  DoubleMatrix mQLsq, mURsq, mDRsq, mLLsq, mSEsq; 
  /// Bilinear Higgs soft parameters: \f$ m_3^2, m_{H_1}^2, m_{H_2}^2 \f$
  /// respectively. 
  double m3sq, mH1sq, mH2sq; 
  double m32;         ///< Gravitino mass
public:
  using Susy::setPars;
  using Susy::setMu;
  using Susy::setLoops;
  using Susy::setThresholds;
  using Susy::setSusy;
  using Susy::setSusyMu;
  using Susy::displayGaugeCoupling;
  using Susy::displayYukawaElement;
  using Susy::displayYukawaMatrix;
  using Susy::displayLoops;
  using Susy::displaySusyMu;
  using Susy::displaySusy;
  using Susy::displayMu;
  using Susy::displayThresholds;
  using Susy::displayGauge;

  typedef Susy susy_type;

  /// Default constructor fills object with zeroes
  SoftPars();
  /// Constructor fills SUSY conserving parts with another object, all
  /// SUSY breaking parameters set to zero
  SoftPars(const Susy &);
  /// Constructor sets all parameters equal to those in another object
  SoftPars(const SoftPars &);
  /// Sets all parameters equal to those in another object
  const SoftPars & operator=(const SoftPars & s);
  /// Constructor sets RPC SUSY parameters to s, gaugino masses to mG,
  /// trilinears to aU, aD, aE for au, ad, ae
  /// trilnears respectively,  \f$m_Q^2\f$=mQl, \f$m_U^2\f$=mUr,
  /// \f$m_D^2\f$=mDr, \f$m_L^2\f$=mLl, \f$m_E^2\f$=mEr, \f$ m_3^2\f$=m3sq,
  /// \f$m_{H_1}^2\f$=mH1sq, \f$m_{H_2}^2\f$=mH2sq, mu parameter, number of
  /// loops=l, and threshold parameter=t
  SoftPars(const Susy & s, const DoubleVector & mG, const
	       DoubleMatrix & aU, const DoubleMatrix & aD, const DoubleMatrix
	       & aE, const DoubleMatrix & mQl, const DoubleMatrix & mUr, const
	       DoubleMatrix & mDr, const DoubleMatrix & mLl, const
	       DoubleMatrix & mEr, double m3sq, double mH1sq, double mH2sq,
	       double mGravitino, double mu, int l, int t);

  /// Returns whole object as a const
  const SoftPars & displaySoftPars() const;

  /// Return a trilinear coupling matrix
  const DoubleMatrix & displayTrilinear(trilinears) const;
  /// Return a trilinear element
  double displayTrilinear(trilinears, int i, int j) const;
  /// Return a trilinear element in "SUGRA style"
  double displaySoftA(trilinears, int, int) const;
  /// Return a soft mass squared matrix
  const DoubleMatrix & displaySoftMassSquared(softMasses) const;
  /// Return a soft mass squared element
  double displaySoftMassSquared(softMasses, int i, int j) const;

  double displayGravitino() const; ///< Returns the gravitino mass
  double displayM3Squared() const;     ///< Return \f$ m_3^2\f$
  double displayMh1Squared() const;    ///< Return \f$m_{H_1}^2\f$
  double displayMh2Squared() const;    ///< Return \f$m_{H_2}^2\f$=mH2sq
  DoubleVector displayGaugino() const; ///< Return \f$M_{G_i}\f$
  double displayGaugino(int i) const;  ///< Return \f$M_{G_i}\f$
  /// Return contents of object in a vector: for RG evolution
  virtual const DoubleVector display() const; 
  
  /// Sets gravitino mass
  void setM32(double);
  /// Sets whole thing equal to another object
  void setSoftPars(SoftPars const &);
  /// Set one element of a soft mass squared matrix 
  void setSoftMassElement(softMasses, int, int, double);
  /// Set whole of a soft mass squared matrix
  void setSoftMassMatrix(softMasses, const DoubleMatrix &);
  /// Set whole of a trilinear SUSY breaking parameter matrix
  void setTrilinearMatrix(trilinears, const DoubleMatrix &);
  /// Set one element of a trilinear SUSY breaking parameter matrix
  void setTrilinearElement(trilinears k, int i, int j, double a);
  /// Set one gaugino mass
  void setGauginoMass(int, double);
  /// Set all gaugino masses
  void setAllGauginos(const DoubleVector &);
  void setM3Squared(double);  ///< Sets \f$ m_3^2\f$
  void setMh1Squared(double); ///< Sets \f$ m_{H_1}^2\f$
  void setMh2Squared(double); ///< Sets \f$ m_{H_2}^2\f$
  /// Sets total set of RGE parameters equal to elements of a vector
  void set(const DoubleVector &);
  //  void setSusy(const Susy &);

  /// Returns double vector containing numerical beta functions of parameters
  DoubleVector beta() const; 
  /// Returns numerical beta functions of parameters  
  SoftPars beta2() const;
  /// Returns numerical beta functions of parameters and Brevity
  SoftPars beta2(Brevity&) const;
  /// Returns derivatives of anomalous dimensions of fields with respect to
  /// renormalisation scale in MSSM for: RH leptons, LH leptons, LH quarks, RH
  /// up quarks, RH down quarks, H1 and H2 respectively
  void anomalousDeriv(DoubleMatrix & gEE, DoubleMatrix & gLL,
		      DoubleMatrix & gQQ, DoubleMatrix & gUU,
		      DoubleMatrix & gDD, 
		      double & gH1H1, double & gH2H2) const;
  /// Ytilde quantities are for calculational brevity in beta functions.
  void yTildes(DoubleMatrix & yu, DoubleMatrix & yd, DoubleMatrix &ye) const;
  
  /// Reads in universal boundary conditions at the current scale:
  /// m0, M1/2, A0, B-parameter and mu
  void universal(double m0,  double m12,  double a0,  double mu,
		 double m3sq);
  /// Give it a SUSY object and a value of M3/2, and it will return a soft
  /// object with AMSB soft breaking terms. Note that the sleptons will be
  /// tachyonic, ie nothing has been done to fix that problem.
  /// Note that in the following, we are neglecting all Yukawa couplings except
  /// that of the third family.
  void addAmsb(double m32);
  /// This flips the signs of various parameters (including mu) with no
  /// physical effect - apparently, through a U(1)_{R-PQ} transformation. See
  /// hep-ph/0312378 
  void u1R_PQflip();
  /// Reads in universal boundary conditions at the current scale: m0, M1/2, A0
  void standardSugra(double m0,  double m12, double a0);
  /// Sets all flavour-diagonal SUSY breaking scalar masses to m0
  void universalScalars(double m0);
  /// Sets all flavour-diagonal SUSY breaking gaugino masses to m12
  void universalGauginos(double m12);
  /// Sets all SUSY breaking trilinear couplings to a0
  void universalTrilinears(double a0);
  /// Boundary conditions to be applied at messenger scale for Gauge mediated
  /// SUSY breaking (see hep-ph/9703211 for example), n5 is the number of
  /// 5-plets, mMess is the messenger scale and lambda is the GMSB scale
  void minimalGmsb(int n5, double lambda, double mMess, double cgrav);  
};

typedef SoftPars<MssmSusy, sBrevity> SoftParsMssm;

/// Formatted ouput of whole object
template<class Susy, class Brevity>
std::ostream & operator <<(std::ostream &left, const SoftPars<Susy, Brevity> &s);

template<class Susy, class Brevity>
SoftPars<Susy, Brevity>::SoftPars()
  : Susy(), mGaugino(3), ua(3, 3), da(3, 3), ea(3, 3),
  mQLsq(3, 3), mURsq(3, 3), mDRsq(3, 3), mLLsq(3, 3), mSEsq(3, 3), m3sq(0.0),
  mH1sq(0.0), mH2sq(0.0), m32(1.e19) {      
  
  setPars(numSoftParsMssm);
  setMu(0.0);
  setLoops(0);
  setThresholds(0);
}

template<class Susy, class Brevity>
SoftPars<Susy, Brevity>::SoftPars(const SoftPars<Susy, Brevity> & s)
  : Susy(s.displaySusy()),
  mGaugino(s.displayGaugino()), ua(s.displayTrilinear(UA)),
  da(s.displayTrilinear(DA)), ea(s.displayTrilinear(EA)),
  mQLsq(s.displaySoftMassSquared(mQl)), 
  mURsq(s.displaySoftMassSquared(mUr)), 
  mDRsq(s.displaySoftMassSquared(mDr)),
  mLLsq(s.displaySoftMassSquared(mLl)),
  mSEsq(s.displaySoftMassSquared(mEr)), 
  m3sq(s.displayM3Squared()), mH1sq(s.displayMh1Squared()),
  mH2sq(s.displayMh2Squared()), m32(s.displayGravitino()) {
  
  setPars(numSoftParsMssm);
  setMu(s.displayMu()); 
  setLoops(s.displayLoops());
  setThresholds(s.displayThresholds());
}

template<class Susy, class Brevity>
SoftPars<Susy, Brevity>::SoftPars(const Susy &s)
  : Susy(s), mGaugino(3), ua(3, 3), da(3, 3), ea(3, 3),
    mQLsq(3, 3), mURsq(3, 3), mDRsq(3, 3), mLLsq(3, 3), mSEsq(3, 3), m3sq(0.0),
    mH1sq(0.0),  mH2sq(0.0), m32(1.e19) { 
      setPars(numSoftParsMssm);
      setMu(s.displayMu()); 
      setLoops(s.displayLoops());
      setThresholds(s.displayThresholds());
}

template<class Susy, class Brevity>
SoftPars<Susy, Brevity>::SoftPars
(const Susy & s, const DoubleVector & mG, const
 DoubleMatrix & aU, const DoubleMatrix & aD, const DoubleMatrix & aE, const
 DoubleMatrix & mQl, const DoubleMatrix & mUr, const DoubleMatrix & mDr, const
 DoubleMatrix & mLl, const DoubleMatrix & mEr, double m3sqn, double mH1sq,
 double mH2sq, double mg, double mu, int l, int t)
  : Susy(s), mGaugino(mG), ua(aU), da(aD), ea(aE),
    mQLsq(mQl), mURsq(mUr), mDRsq(mDr), mLLsq(mLl), mSEsq(mEr), m3sq(m3sqn),
    mH1sq(mH1sq), mH2sq(mH2sq), m32(mg) {
      setPars(numSoftParsMssm);
      setMu(mu);
      setLoops(l);
      setThresholds(t);
}

template<class Susy, class Brevity>
const SoftPars<Susy, Brevity> & SoftPars<Susy, Brevity>::displaySoftPars() const {
   return *this;
}

template<class Susy, class Brevity>
double SoftPars<Susy, Brevity>::displayM3Squared() const {
   return m3sq;
}

template<class Susy, class Brevity>
double SoftPars<Susy, Brevity>::displayMh1Squared() const {
   return mH1sq;
}

template<class Susy, class Brevity>
double SoftPars<Susy, Brevity>::displayMh2Squared() const {
   return mH2sq;
}

template<class Susy, class Brevity>
DoubleVector SoftPars<Susy, Brevity>::displayGaugino() const {
   return mGaugino;
}

template<class Susy, class Brevity>
double SoftPars<Susy, Brevity>::displayGaugino(int i) const {
  return mGaugino.display(i);
}

template<class Susy, class Brevity>
double SoftPars<Susy, Brevity>::displayGravitino() const { return m32; }

template<class Susy, class Brevity>
void SoftPars<Susy, Brevity>::setGauginoMass(int i, double f) {
  mGaugino(i) = f;
}

template<class Susy, class Brevity>
void SoftPars<Susy, Brevity>::setM3Squared(double f) { m3sq = f; }

template<class Susy, class Brevity>
void SoftPars<Susy, Brevity>::setMh1Squared(double f) { mH1sq = f; }

template<class Susy, class Brevity>
void SoftPars<Susy, Brevity>::setMh2Squared(double f) { mH2sq = f; }

template<class Susy, class Brevity>
void SoftPars<Susy, Brevity>::setSoftPars(SoftPars<Susy, Brevity> const & s) { *this = s; }

template<class Susy, class Brevity>
void SoftPars<Susy, Brevity>::setM32(double a) { m32 = a; }

#include "softpars.cpp"

} // namespace softsusy

#endif
