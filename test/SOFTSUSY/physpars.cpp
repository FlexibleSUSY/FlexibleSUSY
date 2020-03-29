
/** \file physpars.cpp
   - Project:     SOFTSUSY 
   - File:        physpars.cpp
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/

*/

#include "physpars.h"
#include "utils.h"

namespace softsusy {

const drBarPars & drBarPars::operator=(const drBarPars &s) {
  if (this == &s) return *this;
  mz = s.mz; mw = s.mw;
  mt = s.mt; mb = s.mb; mtau = s.mtau;
  ht = s.ht; hb = s.hb; htau = s.htau;
  ut = s.ut; ub = s.ub; utau = s.utau;
  nBpmz = s.nBpmz; uBpmz = s.uBpmz; vBpmz = s.vBpmz; 
  mnBpmz = s.mnBpmz; mchBpmz = s.mchBpmz;
  setsPhysical(s.displaysPhysical());
  
  return *this;
}

const sPhysical & sPhysical::operator=(const sPhysical &s) {
  if (this == &s) return *this;
  mh0 = s.mh0; mA0 = s.mA0; mHpm = s.mHpm;
  msnu = s.msnu; 
  mch = s.mch; mneut = s.mneut; mixNeut = s.mixNeut;
  thetaL = s.thetaL; thetaR = s.thetaR; mGluino = s.mGluino;
  thetat = s.thetat; thetab = s.thetab; thetatau = s.thetatau;
  mu = s.mu; md = s.md; me = s.me; thetaH = s.thetaH;  
  thetaA0 = s.thetaA0; mixh0 = s.mixh0; 
  return *this;
}

// a should be in C convention ie start from index zero
void sPhysical::display(double *a) const {
  std::size_t k = 0;
  mh0.fillArray(a,k);   k += mh0.size();
  mA0.fillArray(a,k);   k += mA0.size();
  a[k++] = mHpm;
  msnu.fillArray(a,k);  k += msnu.size();
  mch.fillArray(a,k);   k += mch.size();
  mneut.fillArray(a,k); k += mneut.size();
  a[k++] = mGluino;
  mixNeut.fillArray(a,k); k += mixNeut.size();
  a[k++] = thetaL;
  a[k++] = thetaR;
  a[k++] = thetat;
  a[k++] = thetab;
  a[k++] = thetatau;
  mu.fillArray(a,k); k += mu.size();
  md.fillArray(a,k); k += md.size();
  me.fillArray(a,k); k += me.size();
  a[k++] = thetaH;
  mixh0.fillArray(a,k); k += mixh0.size();
  a[k++] = thetaA0;
}

std::ostream & operator <<(std::ostream & left, const drBarPars &s) {
  left << s.displaysPhysical();
  left << "BPMZ conventions, N" << s.nBpmz << "U" << s.uBpmz << "V" 
       << s.vBpmz;
  left << "mt: "  << s.mt << " mb: " << s.mb << " mtau: " << s.mtau << std::endl;
  left << "ht: "  << s.ht << " hb: " << s.hb << " htau: " << s.htau << std::endl;
  left << "Ut: "  << s.ut << " Ub: " << s.ub << " Utau: " << s.utau << std::endl;
  left << "mz: "  << s.mz << " mw: " << s.mw << std::endl;

  return left;
}

std::ostream & operator <<(std::ostream & left, const sPhysical &s) {
  left << "mh^0: " << s.mh0 << "mA^0: " << s.mA0
       << "mH^+-: " << s.mHpm << "\n";
  left << "alpha: " << s.thetaH << "\n";
  left << "sneutrinos" << s.msnu; 
  left << "mU~" << s.mu << "mD~" << s.md << "mE~" << s.me;
  left << "thetat: " << s.thetat << " thetab: " << s.thetab << 
    " thetatau: " << s.thetatau << "\n";
  left << "mGluino:  " << s.mGluino << "\n";
  left << "charginos" << s.mch;
  left << "thetaL: " << s.thetaL << " thetaR: " << s.thetaR << "\n";
  left << "neutralinos" << s.mneut;
  left << "neutralino mixing matrix " << s.mixNeut;
  left << "CP even Higgs mixing matrix: " << s.mixh0;
  left << "CP odd mixing angle: " << s.thetaA0;
  return left;
}

std::ostream & operator <<(std::ostream &st, const sProblem & p) {
  if (!p.test()) return st;
  st << "[ ";
  if (p.mgutOutOfBounds) st << "GUT scale too high or too low ";
  if (p.badConvergence) st << "No acceptable solution found ";
  if (p.irqfp) st << "Quasi-fixed point breached ";
  if (p.noMuConvergence) st << "No mu convergence ";
  if (p.noRhoConvergence) st << "No rho convergence ";
  if (p.nonperturbative) st << "Non-perturbative ";
  if (p.noConvergence) st << "No convergence ";
  if (p.tachyon) {
    /*    switch (p.tachyonType) {
          case selectron:
    case smuon:
      if (p.tachyonType == )
      }*/
    st << tachyonNames[p.tachyon] << " tachyon ";
  }
  if (p.muSqWrongSign) st << "MuSqWrongsign ";
  if (p.m3sq) st << "m3sq-problem ";
  if (p.higgsUfb) st << "Higgs potential ufb ";
  if (p.inaccurateHiggsMass) st << "Inaccurate Higgs mass ";
  if (p.problemThrown) st << "Numerical problemThrown ";
  st << "]";
  return st;
}

const sProblem & sProblem::operator=(const sProblem &s) {
  if (this == &s) return *this;
  mgutOutOfBounds = s.mgutOutOfBounds;
  irqfp = s.irqfp;
  badConvergence = s.badConvergence;
  noMuConvergence = s.noMuConvergence;
  noRhoConvergence = s.noRhoConvergence;
  nonperturbative = s.nonperturbative;
  noConvergence = s.noConvergence;
  tachyon = s.tachyon;
  muSqWrongSign = s.muSqWrongSign;
  higgsUfb = s.higgsUfb;
  m3sq = s.m3sq;
  problemThrown = s.problemThrown;
  return *this;
}

// Returns mixing matrix o and neutralino masses mn in the MPZ convention
// (hep-ph/9606211), n is 4 by 4 and mneut is 1->4.
void drBarPars::mpzNeutralinos() { 
  // We want to change the PHASES of the neutralino mixing matrix in order to
  // produce POSITIVE neutralino masses, a la Matchev, Pierce and Zhang
  const int rank = mneut.displayEnd();
  DoubleVector temp(mneut);
  
  ComplexMatrix K(rank, rank);
  int i; for (i=1; i<=rank; i++) 
    if (mneut.display(i) < 0.0) K(i, i) = Complex(0.0, -1.0);
    else
      K(i, i) = Complex(1.0, 0.0);
  
  mnBpmz = temp.apply(fabs);
  nBpmz = K * mixNeut.transpose();
}

// Returns mixing matrices u,v and neutralino masses mneut in the MPZ
// convention (hep-ph/9606211),  u+v are (2,2) and mch is 1->2.
void drBarPars::mpzCharginos() {
  // We want to change the PHASES of the neutralino mixing matrix in order to
  // produce POSITIVE neutralino masses, a la Matchev, Pierce and Zhang
  ComplexMatrix u(2, 2), v(2, 2);
  positivise(thetaL, thetaR, mch, u, v);
  uBpmz = u; vBpmz = v;
  mchBpmz = mch.apply(fabs); 
}

} // namespace softsusy

double sTfn(double sTins, double sTouts) {
  double sTin  = fabs(sTins);
  double sTout = fabs(sTouts);
  if (sTin < 1. && sTout < 1.) return fabs(sTin - sTout);
  else return fabs(1.0 - minimum(sTin, sTout) / maximum(sTin, sTout));
}

/// LCT: Difference between two drBarPars objects
void sumTol(const softsusy::drBarPars & a, const softsusy::drBarPars & b, softsusy::DoubleVector & sT) {
  int k = 1;
  
  sT(k) = sTfn(a.mGluino, b.mGluino); k++;
  int i; for (i=1; i<=a.mh0.displayEnd(); i++) {
    sT(k) = sTfn(a.mh0(i), b.mh0(i)); k++;
  }
  for (i=1; i<=a.mA0.displayEnd(); i++) {
    sT(k) = sTfn(a.mA0(i), b.mA0(i)); k++;
  }
  sT(k) = sTfn(a.mHpm, b.mHpm); k++;
  for (i=1; i<=3; i++) {
    sT(k) = sTfn(a.msnu(i), b.msnu(i)); k++;
  }
  for (i=1; i<=2; i++) {
    sT(k) = sTfn(a.mch(i), b.mch(i)); k++;
  }
  for (i=1; i<=a.mneut.displayEnd(); i++) {
    sT(k) = sTfn(a.mneut(i), b.mneut(i)); k++;
  }
  int j; for (j=1; j<=3; j++)
    for(i=1; i<=2; i++) {
      sT(k) = sTfn(a.mu(i, j), b.mu(i, j)); k++;
      sT(k) = sTfn(a.md(i, j), b.md(i, j)); k++;
      sT(k) = sTfn(a.me(i, j), b.me(i, j)); k++;
    }
}
