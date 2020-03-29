
/** \file physpars.h
   - Project:     SOFTSUSY 
   - File:        physpars.h
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/

   \brief Flags for potential problems in sProblem structure, and structure
   for containing physical MSSM parameters such as pole masses in sPhysical. 
*/

#ifndef PHYSPARSH
#define PHYSPARSH

#include <iostream>
#include "linalg.h"

namespace softsusy {
  typedef enum {none=0, selectron, smuon, stau, sup, scharm, stop, sdown, 
             sstrange, sbottom, h0, A0, hpm, snue, snumu, snutau, W, Z} 
  tachyonType;

  const static char* tachyonNames[18] = {
   "none", "selectron", "smuon", "stau", "sup", "scharm", "stop", "sdown", 
   "sstrange", "sbottom", "h0", "A0", "hpm", "snue", "snumu", "snutau", 
   "W", "Z"};

/// Masses of the physical particles. 
struct sPhysical {
  DoubleVector mh0, mA0;
  double mHpm;
  /// sneutrino masses
  DoubleVector msnu;
  /// chargino/neutralino masses: ordered
  DoubleVector mch, mneut;
  /// Gluino mass
  double mGluino;
  /// neutralino mixing
  DoubleMatrix mixNeut;
  /// chargino and third family mixing angles
  double thetaL, thetaR, thetat, thetab, thetatau;
  /// sparticle masses in order (i=L/R, family)
  DoubleMatrix mu, md, me;
  /// Higgs mixing angle (alpha)
   double thetaH, thetaA0;
  /// Higgs mixing matrices
  DoubleMatrix mixh0;
  
  sPhysical(); ///< Constructor: initialises with zeroes
  sPhysical(const sPhysical &); ///< Constructor copies another object
  /// Displays the object in situations where a const is required
  const sPhysical & displaysPhysical() const { return *this; };
  /// Sets whole contents to those of another object s
  void setsPhysical(const sPhysical &s) { *this = s; };
  /// Sets whole contents to those of another object s
  const sPhysical & operator = (const sPhysical & s);
  
  /// Displays contents in a C-style convention *a (starts with index zero)
  void display(double *a) const;

  bool operator==(const sPhysical& o) const {
     return (mh0(1) == o.mh0(1))
        && (mA0(1) == o.mA0(1))
        && (mh0(2) == o.mh0(2))
        && (mHpm == o.mHpm)
        && (msnu == o.msnu)
        && (mch == o.mch)
        && (mneut == o.mneut)
        && (mGluino == o.mGluino)
        && (mixNeut == o.mixNeut)
        && (thetaL == o.thetaL)
        && (thetaR == o.thetaR)
        && (thetat == o.thetat)
        && (thetab == o.thetab)
        && (thetatau == o.thetatau)
        && (mu == o.mu)
        && (md == o.md)
        && (me == o.me)
        && (thetaH == o.thetaH);
  }
}; 

/// Formatted printout
std::ostream & operator <<(std::ostream &, const sPhysical &); 

/// Various boolean values to flag any problems with the parameter point
struct sProblem {
  bool mgutOutOfBounds;
  bool badConvergence; ///< Nowhere near a decent solution
  bool irqfp; ///< Infra-red quasi fixed point breached
  bool noRhoConvergence; ///< Couldn't calculate electroweak rho parameter
  bool noConvergence; ///< Iteration did not converge: not always serious
  tachyonType tachyon; ///< Tachyonic point
  bool muSqWrongSign; ///< mu^2 came out with wrong sign; no REWSB
  bool m3sq; ///< m3sq came out with wrong sign; no REWSB
  bool higgsUfb; ///< Higgs potential inconsistent with a good minimum
  bool nonperturbative; ///< Running went non-perturbative
  bool noMuConvergence; ///< mu couldn't be calculated
  /// Higgs mass is potentially inaccurate and cant be trusted
  bool inaccurateHiggsMass; 
  bool problemThrown;
  /// returns true if there's any problem 
  bool test() const 
  {return (mgutOutOfBounds || irqfp || noConvergence || tachyon || 
	   muSqWrongSign || higgsUfb || nonperturbative || noRhoConvergence || 
	   noMuConvergence || m3sq || badConvergence || inaccurateHiggsMass ||
	   problemThrown);}; 
  /// Only returns true if there's a serious problem
  bool testSeriousProblem() const 
  {return (irqfp || tachyon || muSqWrongSign || higgsUfb || nonperturbative 
	   || noRhoConvergence || noMuConvergence || m3sq || badConvergence ||
	   mgutOutOfBounds || problemThrown);}; 

  inline sProblem(); ///< constructor full of false values
  /// Constructor that sets flags equal to those of s
  inline sProblem(const sProblem & s);
  /// Sets flags equal to those of s
  const sProblem & operator = (const sProblem &);
  
};
/// Formatted output, but won't print unflagged problems
std::ostream & operator <<(std::ostream &st, const sProblem & p);

/// DRbar values of masses and mixings in MSSM
struct drBarPars: public sPhysical { 
  double mz, mw;       /// Running electroweak gauge boson masses
  double mt, mb, mtau; /// Running top, bottom and tau mass
  double ht, hb, htau; /// Running diagonal top, bottom, tau Yukawa couplings
  double ut, ub, utau; /// Running diagonal top, bottom, tau trlinear couplings
  /// BPMZ convention mixing matrices for neutralinos and charginos
  ComplexMatrix nBpmz, uBpmz, vBpmz; 
  /// positive definite masses for neutralinos and charginos
  softsusy::DoubleVector mnBpmz, mchBpmz; 

  inline drBarPars(); ///< Initialises with zero values
  inline drBarPars(const drBarPars &); ///< Initialises with another object
  /// Sets contents equal to those of another object
  const drBarPars & operator = (const drBarPars &s);

  /// Returns mixing matrix o and neutralino masses mn in the MPZ convention
  /// (hep-ph/9606211), n is 4 by 4 and mneut is 1->4, ie
  /// Calculates mnBpmz, nBpmz
  void mpzNeutralinos(); 
  /// Calculates uBpmz, vBpmz, mchBpmz, ie
  /// Returns mixing matrices u,v and neutralino masses mneut in the MPZ
  /// convention (hep-ph/9606211),  u+v are (2,2) and mch is 1->2.
  void mpzCharginos();
};

std::ostream & operator <<(std::ostream &, const drBarPars &); 

// ---------------------- inline class members ------------------------
inline drBarPars::drBarPars()
  : sPhysical(), mz(0.0), mw(0.0), 
    mt(0.0), mb(0.0), mtau(0.0), ht(0.0), hb(0.0), htau(0.0), 
    ut(0.0), ub(0.0), utau(0.0),
    nBpmz(4, 4), uBpmz(2, 2), 
    vBpmz(2, 2), mnBpmz(4), mchBpmz(2)
{}

inline drBarPars::drBarPars(const drBarPars &s)
  : sPhysical(s.displaysPhysical()), mz(s.mz), mw(s.mw),
    mt(s.mt), mb(s.mb), mtau(s.mtau), ht(s.ht), hb(s.hb), htau(s.htau),
    ut(s.ut), ub(s.ub), utau(s.utau),
    nBpmz(s.nBpmz), uBpmz(s.uBpmz), vBpmz(s.vBpmz), mnBpmz(s.mnBpmz), 
    mchBpmz(s.mchBpmz)
{}

inline sPhysical::sPhysical()
  : mh0(2), mA0(1), mHpm(0.), msnu(3), mch(2), mneut(4), 
    mGluino(0.0),
    mixNeut(4, 4), thetaL(0.0), thetaR(0.0), thetat(0.0), thetab(0.0),
                 thetatau(0.0), mu(2, 3), md(2, 3), me(2, 3), thetaH(0.0), 
                 thetaA0(0.0), mixh0(2,2)
{}

inline sPhysical::sPhysical(const sPhysical & s)
  : mh0(s.mh0), mA0(s.mA0), mHpm(s.mHpm), msnu(s.msnu), mch(s.mch), 
    mneut(s.mneut), mGluino(s.mGluino), mixNeut(s.mixNeut), thetaL(s.thetaL),
    thetaR(s.thetaR), thetat(s.thetat), thetab(s.thetab),
    thetatau(s.thetatau), mu(s.mu), md(s.md), me(s.me), thetaH(s.thetaH),
    thetaA0(s.thetaA0), mixh0(s.mixh0)
{}

inline sProblem::sProblem()
  : mgutOutOfBounds(false), badConvergence(false), 
    irqfp(false), noRhoConvergence(false), noConvergence(false),
    tachyon(none), muSqWrongSign(false), m3sq(false), higgsUfb(false), 
    nonperturbative(false), noMuConvergence(false),     
    inaccurateHiggsMass(false), problemThrown(false)
{}

inline sProblem::sProblem(const sProblem & s)
  : mgutOutOfBounds(s.mgutOutOfBounds), badConvergence(s.badConvergence), 
    irqfp(s.irqfp), noRhoConvergence(s.noRhoConvergence), 
    noConvergence(s.noConvergence),
    tachyon(s.tachyon), muSqWrongSign(s.muSqWrongSign), m3sq(s.m3sq),
    higgsUfb(s.higgsUfb), nonperturbative(s.nonperturbative), 
    noMuConvergence(s.noMuConvergence), 
    inaccurateHiggsMass(s.inaccurateHiggsMass), problemThrown(s.problemThrown)
{}

} // namespace softsusy

/// Returns the relative difference between the orders of magnitude of two
/// numbers unless they are less than one, in which case it returns the value
/// of the absolute difference
double sTfn(double sTins, double sTouts);

/// LCT: Calculates fractional difference in Drbar masses between a and b
namespace softsusy {
   class drBarPars;
   class DoubleVector;
}
void sumTol(const softsusy::drBarPars & a, const softsusy::drBarPars & b, softsusy::DoubleVector & sT);

#endif
