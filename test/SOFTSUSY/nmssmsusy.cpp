
/** \file nmssmsusy.cpp
   - Project:     SOFTSUSY
   - Author:      Ben Allanach, Peter Athron, Alexander Voigt
   - Manual:      hep-ph/0104145
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/
   - Description: All NMSSM SUSY respecting parameters
*/

#include "nmssmsusy.h"
#include <cassert>

namespace softsusy {

const nmsBrevity & nmsBrevity::operator=(const sBrevity &s) {
  if (this == &s) return *this;
  sBrevity::operator=(s);
  return *this;
}

const nmsBrevity & nmsBrevity::operator=(const nmsBrevity &s) {
  if (this == &s) return *this;
  sBrevity::operator=(s);
  lsq = s.lsq; ksq = s.ksq; l4 = s.l4; k4 = s.k4;
  return *this;
}

void nmsBrevity::calculate(const DoubleMatrix & yu, const DoubleMatrix & yd,
                         const DoubleMatrix & ye, const DoubleVector & g,
                         const double & lam, const double & kap) {
  sBrevity::calculate(yu, yd, ye, g);
  lsq = lam * lam; ksq = kap * kap;
  l4 = lsq * lsq; k4 = ksq * ksq;
}


NmssmSusy::NmssmSusy()
   : MssmSusy()
   , lambda(0.0), kappa(0.0), sVev(0.0), xiF(0.0), mupr(0.0)
{
    setPars(numNMssmPars);
}

NmssmSusy::NmssmSusy(const NmssmSusy &s)
   : MssmSusy(s)
   , lambda(s.lambda), kappa(s.kappa), sVev(s.sVev), xiF(s.xiF), mupr(s.mupr)
{
    setPars(numNMssmPars);
}

NmssmSusy::NmssmSusy(const MssmSusy &m)
   : MssmSusy(m)
   , lambda(0), kappa(0), sVev(0), xiF(0), mupr(0)
{
    setPars(numNMssmPars);
}


NmssmSusy::NmssmSusy(const DoubleMatrix & u, const DoubleMatrix & d, const
		     DoubleMatrix & e, const DoubleVector & v, double m,
		     double tb, double MU, int l, int t, double hv, double sv,
                     double lam, double kap, double mpr, double z)
   : MssmSusy(u, d, e, v, m, tb, MU, l, t, hv)
   , lambda(lam), kappa(kap), sVev(sv), xiF(z), mupr(mpr)
{
    setPars(numNMssmPars);
}


NmssmSusy::~NmssmSusy() {
}

const NmssmSusy & NmssmSusy::operator=(const NmssmSusy & s) {
  if (this == &s) return *this;
  MssmSusy::operator=(s);
  sVev = s.sVev;
  lambda = s.lambda;
  kappa = s.kappa;
  mupr = s.mupr;
  xiF = s.xiF;
  return *this;
}

const NmssmSusy & NmssmSusy::operator=(const MssmSusy & s) {
  if (this == &s) return *this;
  MssmSusy::operator=(s);
  return *this;
}

void NmssmSusy::setSomePars(const NmssmSusy & s) {
  MssmSusy::setSomePars(s);
}

const DoubleVector NmssmSusy::display() const {
  DoubleVector y(MssmSusy::display());
  assert(y.displayStart() == 1 && y.displayEnd() == numSusyPars);
  y.setEnd(numNMssmPars);
  y(34) = sVev;
  y(35) = lambda;
  y(36) = kappa;
  y(37) = mupr;
  y(38) = xiF;
  return y;
}

void NmssmSusy::set(const DoubleVector & y) {
  assert(y.displayEnd() - y.displayStart() + 1 >= numNMssmPars);
  MssmSusy::set(y);
  sVev = y.display(34);
  lambda = y.display(35);
  kappa = y.display(36);
  mupr = y.display(37);
  xiF = y.display(38);
}

std::ostream & operator <<(std::ostream &left, const NmssmSusy &s) {
  left << static_cast<MssmSusy>(s)
       << "singlet VEV: " << s.displaySvev()
       << " lambda: " << s.displayLambda()
       << " kappa: " << s.displayKappa()
       << " smu: " << s.displaySusyMu()
       << " mupr: " << s.displayMupr()
       << " xiF: " << s.displayXiF()
       << '\n';
  return left;
}

void NmssmSusy::setSusy(const NmssmSusy & s) {
  MssmSusy::setSusy(s);
  setLambda(s.displayLambda());
  setKappa(s.displayKappa());
  setSvev(s.displaySvev());
  setMupr(s.displayMupr());
  setXiF(s.displayXiF());
}

// Outputs derivatives (DRbar scheme) in the form of ds. a contains the
// matrices calculated that are handy for computation.
// W=  LL Y^E H1 ER + QL Y^D H1 DR + QL Y^U H2 UR + smu H2 H1
// is the superpotential. Consistent with Allanach, Dedes, Dreiner
// hep-ph/9902251 and Barger, Berger and Ohmann hep-ph/9209232, 9311269
// EXCEPT for the sign of smu, which is opposite. These equations are also
// valid for W=  - LL Y^E H1 ER - QL Y^D H1 DR + QL Y^U H2 UR + smu H2 H1, the
// New SOFTSUSY convention
NmssmSusy NmssmSusy::beta(nmsBrevity & a) const {
  // Wave function renormalisations: convention for g**(i, j) is that i is the
  // LOWER index and j the upper in our paper hep-ph/9902251
  static DoubleMatrix gEE(3, 3), gLL(3, 3), gQQ(3, 3), gDD(3, 3),
    gUU(3, 3);

  double gH1H1=0.0, gH2H2=0.0, gSS = 0.0;
  static DoubleVector dg(1,3);

  // keep this option in order to interface with RPVSUSY
  anomalousDimension(gEE, gLL, gQQ, gUU, gDD, dg, gH1H1, gH2H2, gSS, a);

  // To keep this a const function
  const DoubleMatrix &u1 = displayYukawaMatrix(YU).display(),
     &d1 = displayYukawaMatrix(YD).display(),
     &e1 = displayYukawaMatrix(YE).display();

  // contain derivatives of up, down quarks and leptons
  static DoubleMatrix du(3, 3), dd(3, 3), de(3, 3);
  static double dl, dk, dz, dmupr;
  // mu parameter derivatives
  double dmu;

  // RGEs of SUSY parameters
  du = u1 * (gUU + gH2H2) + gQQ * u1;
  dd = d1 * (gDD + gH1H1) + gQQ * d1;
  de = e1 * (gEE + gH1H1) + gLL * e1;
  dl = lambda * (gH1H1 + gH2H2 + gSS);
  dk = kappa * (3.0 * gSS);
  dmu = displaySusyMu() * (gH1H1 + gH2H2);
  dmupr = 2.0 * mupr * gSS;
  dz = xiF * gSS;

  // Following is from hep-ph/9308335: scalar H anomalous dimensions (as
  // opposed to the chiral superfield one - see hep-ph/0111209).
  // Additional contribution from Feynman gauge running at two-loops of tan
  // beta: we need this to link up with BPMZ: hep-ph/0112251
  const double &uuT = a.uuT, &ddT = a.ddT, &eeT = a.eeT;
  const DoubleVector &gsq=a.gsq, &g4 = a.g4;
  const DoubleMatrix &u2=a.u2, &d2=a.d2, &e2=a.e2, &d2t=a.d2t;
  const double t = (d2 * u2).trace();
  const double &lsq = a.lsq, &ksq = a.ksq, &l4 = a.l4, &k4 = a.k4;
  static const double oneLoop = 1.0 / (16.0 * sqr(PI));
  double sH1H1 = oneLoop * (3.0 * ddT + eeT + lsq);
  double sH2H2 = oneLoop * (3.0 * uuT + lsq);

  const static double twolp = 4.010149318236068e-5; // 1/(16 pi^2)^2
  if (displayLoops() > 1) {
    const double g4terms = 1.035 * g4(1) + 0.45 * gsq(1) * gsq(2) + 5.875 * g4(2);
    sH1H1 = sH1H1 + twolp *
      (-(3.0 * (e2 * e2).trace() + 9.0 * (d2t * d2t).trace() + 3.0 * t) +
       (16 * gsq(3) - 0.4 * gsq(1)) * ddT + 1.2 * gsq(1) * eeT
       + g4terms
       - 2 * ksq * lsq - 3 * l4 - 3 * lsq * uuT);
    sH2H2 = sH2H2 + twolp *
      (- (9.0 * (u2 * u2).trace() + 3.0 * t) +
       (16 * gsq(3) + 0.8 * gsq(1)) * uuT
       + g4terms
       - 2 * ksq * lsq - 3 * l4 - lsq * (3 * ddT + eeT));
  }

  double cosb2 = sqr(cos(atan(displayTanb()))), sinb2 = 1.0 - cosb2;
  double feynman = 1.5 * gsq(2) + 0.3 * gsq(1);
  /// One-loop RGEs in Feynman gauge
  double dt = displayTanb() * (sH1H1 - sH2H2);
  double dHvev = displayHvev() *
    (cosb2 * (-sH1H1 + feynman * oneLoop) +
     sinb2 * (-sH2H2 + feynman * oneLoop));
  double dSvev = - sVev * oneLoop * 2 * (lsq + ksq);

  if (displayLoops() > 1) {
    /// Two-loop pieces
    dt = dt + displayTanb() * twolp * (3.0 * ddT + eeT - 3.0 * uuT) * feynman;
    dHvev = dHvev - displayHvev() * twolp * (cosb2 * (3.0 * ddT + eeT) +
				    sinb2 * 3.0 * uuT + lsq) * feynman
                  + displayHvev() * twolp * 4.5 * g4(2);
    dSvev = dSvev + sVev * twolp * (8 * k4 + 8 * ksq * lsq
            + 2 * lsq * (-2 * feynman + 2 * lsq + 3 * ddT + eeT + 3 * uuT));
  }

  // Contains all susy derivatives:
  NmssmSusy ds(du, dd, de, dg, dmu, dt, displayMu(), displayLoops(),
	       displayThresholds(), dHvev, dSvev, dl, dk, dmupr, dz);

  return ds;
}

void nmsetBetas(DoubleMatrix & babBeta, DoubleVector &cuBeta, DoubleVector
                & cdBeta, DoubleVector & ceBeta, DoubleVector & clBeta,
                DoubleVector & bBeta) {
   setBetas(babBeta, cuBeta, cdBeta, ceBeta, bBeta);
   setBetaLambda(clBeta);
}

void setBetaLambda(DoubleVector& clBeta) {
   clBeta(1) = 6.0 / 5.0;  clBeta(2) = 2.0; clBeta(3) = 0.0;
}

// outputs one-loop anomlous dimensions gii given matrix inputs
// Note that we use the convention (for matrices in terms of gamma's)
// gamma^Li_Lj = M_ij for LH fields and
// gamma^Rj_Ri = M_ij for RH fields (since they are really the complex
// conjugates of the RH fields): CHECKED 23/5/02
void NmssmSusy::getOneLpAnom(DoubleMatrix & gEE, DoubleMatrix & gLL,
				DoubleMatrix & gQQ, DoubleMatrix & gDD,
				DoubleMatrix & gUU, double & gH1H1, double &
                                gH2H2, double & gSS, nmsBrevity & a) const {
  const double ksq = a.ksq, lsq = a.lsq;
  static const double oneO16Pisq = 1.0 / (16.0 * sqr(PI));

  MssmSusy::getOneLpAnom(gEE, gLL, gQQ, gDD, gUU, gH1H1, gH2H2, a);

  gH1H1 += oneO16Pisq * lsq;
  gH2H2 += oneO16Pisq * lsq;
  gSS    = oneO16Pisq * (2.0 * lsq + 2.0 * ksq);
}


// adds two-loop anomalous dimension contribution to gii given matrix inputs
// g^Li_Lj = m_{ij} for LH fields
// g^Ei_Ej = m_{ji} for RH fields CHECKED: 23/5/02
void NmssmSusy::getTwoLpAnom(DoubleMatrix & gEE, DoubleMatrix & gLL,
				DoubleMatrix & gQQ, DoubleMatrix & gDD,
				DoubleMatrix & gUU, double & gH1H1, double &
                                gH2H2, double & gSS, nmsBrevity & a) const {
  // For calculational brevity
  DoubleMatrix &u2=a.u2, &d2=a.d2, &e2=a.e2,
    &u2t=a.u2t, &d2t=a.d2t, &e2t=a.e2t;
  double &uuT = a.uuT, &ddT = a.ddT, &eeT = a.eeT;
  double &lsq = a.lsq, &ksq = a.ksq, &l4 = a.l4, &k4 = a.k4;
  DoubleVector &gsq = a.gsq;

  // Everything gets the (1/16pi^2)^2 factor at the bottom
  double h1h1, h2h2;

  // Two-loop pure gauge anom dimensions
  double ss = 1.2 * lsq * gsq(1) + 6.0 * lsq * gsq(2);

  // Two-loop pure Yukawa contributions
  const double s = (eeT + 3.0 * ddT);

  const DoubleMatrix ll(- lsq * e2);
  const DoubleMatrix ee(- 2.0 * lsq * e2t);
  const DoubleMatrix qq(- lsq * (u2 + d2));
  const DoubleMatrix dd(- 2.0 * lsq * d2t);
  const DoubleMatrix uu(- 2.0 * lsq * u2t);
  h1h1 = - (3.0 * l4 + 2.0 * lsq * ksq + 3.0 * lsq * uuT);
  h2h2 = - (3.0 * l4 + 2.0 * lsq * ksq + lsq * s);
  ss = ss - (4.0 * l4 + 8.0 * k4 + 8.0 * lsq * ksq + 2.0 * lsq * s
             + 6.0 * lsq * uuT);

  const static double twolp = 4.010149318236068e-5; // 1/(16 pi^2)^2

  MssmSusy::getTwoLpAnom(gEE, gLL, gQQ, gDD, gUU, gH1H1, gH2H2, a);

  gLL = gLL + twolp * ll;
  gEE = gEE + twolp * ee;
  gQQ = gQQ + twolp * qq;
  gDD = gDD + twolp * dd;
  gUU = gUU + twolp * uu;
  gH1H1 = gH1H1 + twolp * h1h1;
  gH2H2 = gH2H2 + twolp * h2h2;
  gSS = gSS + twolp * ss;
}

// Outputs wave function renormalisation for SUSY parameters and gauge beta
// functions up to 2 loops.
void NmssmSusy::anomalousDimension(DoubleMatrix & gEE, DoubleMatrix & gLL,
				    DoubleMatrix & gQQ, DoubleMatrix & gUU,
				    DoubleMatrix & gDD, DoubleVector & dg,
				    double & gH1H1, double & gH2H2,
                                    double & gSS,
				    nmsBrevity & a)  const {
  // Constants for gauge running
  static DoubleVector bBeta(3), cuBeta(3), cdBeta(3), ceBeta(3), clBeta(3);
  static DoubleMatrix babBeta(3, 3);
  if (bBeta(1) < 1.0e-5) // Constants not set yet
     nmsetBetas(babBeta, cuBeta, cdBeta, ceBeta, clBeta, bBeta);

  // nmsBrevity a contains all of the shortcutted matrices etc;
  a.calculate(displayYukawaMatrix(YU).display(),
              displayYukawaMatrix(YD).display(),
              displayYukawaMatrix(YE).display(),
              displayGauge().display(), lambda, kappa);

  // For calculational brevity
  double &uuT = a.uuT, &ddT = a.ddT, &eeT = a.eeT, &lsq = a.lsq;
  DoubleVector &gsq=a.gsq, &g3=a.g3;

  // 1 loop contributions:
  if (displayLoops() > 0) {
    static const double oneO16Pisq = 1.0 / (16.0 * sqr(PI));
    getOneLpAnom(gEE, gLL, gQQ, gDD, gUU, gH1H1, gH2H2, gSS, a);
    dg = oneO16Pisq * g3 * bBeta;
  }

  if (displayLoops() > 1) {
    getTwoLpAnom(gEE, gLL, gQQ, gDD, gUU, gH1H1, gH2H2, gSS, a);
    const static double twolp = 4.010149318236068e-5;
    dg = dg + g3 * (babBeta * gsq - cuBeta * uuT - cdBeta *
		    ddT - ceBeta * eeT - clBeta * lsq) * twolp;
    //PA: checked numerically 14/9/2012
  }
}

// Outputs derivatives vector y[n] for SUSY parameters: interfaces to
// integration routines
DoubleVector NmssmSusy::beta() const {
  static nmsBrevity a;

  // calculate the derivatives
  static NmssmSusy ds;

  ds = beta(a);

  return ds.display(); // convert to a long vector
}

} // namespace softsusy
