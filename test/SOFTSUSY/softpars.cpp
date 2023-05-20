
/** \file softpars.cpp
   - Project:     SOFTSUSY 
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/
   - Description: includes all MSSM parameters in the Lagrangian

*/

#ifdef SOFTPARS_H

using std::fabs;
using std::log;
using std::cout;
using std::endl;

template<class Susy, class Brevity>
const SoftPars<Susy, Brevity> & SoftPars<Susy, Brevity>::operator=(const SoftPars<Susy, Brevity> & s) {
  if (this == &s) return *this;
  mGaugino = s.mGaugino;
  ua = s.ua;
  da = s.da;
  ea = s.ea;
  m32 = s.m32;  
  mQLsq = s.mQLsq;
  mURsq = s.mURsq;
  mDRsq = s.mDRsq;
  mLLsq = s.mLLsq;
  mSEsq = s.mSEsq;
  m3sq = s.m3sq;
  mH1sq = s.mH1sq;
  mH2sq = s.mH2sq;
  setSusy(s.displaySusy());
  setMu(s.displayMu());
  setLoops(s.displayLoops());
  setThresholds(s.displayThresholds());
  return *this;
}

template<class Susy, class Brevity>
const DoubleMatrix & SoftPars<Susy, Brevity>::displayTrilinear(trilinears k) const {
  switch(k) {
  case UA: return ua; break;
  case DA: return da; break;
  case EA: return ea; break;
  default: 
    std::ostringstream ii;
    ii << "In SoftPars::displayTrilinear, called with illegal argument";
    ii << " " << k << endl;
    throw ii.str();
    break;
  }
}

template<class Susy, class Brevity>
double SoftPars<Susy, Brevity>::displayTrilinear(trilinears k, int i, int j) 
  const {
  switch(k) {
  case UA: return ua.display(i, j); break;
  case DA: return da.display(i, j); break;
  case EA: return ea.display(i, j); break;
  default: 
    std::ostringstream ii;
    ii << "In SoftPars::displayTrilinear, called with illegal argument";
    ii << " " << k << endl;
    throw ii.str();
    break;
  }
}

// Will display zero if that's what the A term is regardless of what the
// Yukawa coupling is (even if it's zero).
template<class Susy, class Brevity>
double SoftPars<Susy, Brevity>::displaySoftA(trilinears k, int i, int j) const {
  double am = displayTrilinear(k, i, j);
  if (fabs(am) < EPSTOL) return 0.0;
  
  if (fabs(displayYukawaElement(yukawa(k), i, j)) < 1.0e-100) {
    std::ostringstream ii;
    ii << "WARNING: asking for SoftPars::displaySoftA(" << int(k) << "," 
	 << i << "," << j << "), where Yukawa coupling is " <<
      fabs(displayYukawaElement(yukawa(k), i, j)) << endl;
    throw ii.str();
  }
  double temp = displayTrilinear(k, i, j) 
    / displayYukawaElement(yukawa(k), i, j); 
  return temp;
}

template<class Susy, class Brevity>
const DoubleMatrix & SoftPars<Susy, Brevity>::displaySoftMassSquared(softMasses k) const {
  switch(k) {
  case mQl: return mQLsq; break;
  case mUr: return mURsq; break;
  case mDr: return mDRsq; break;
  case mLl: return mLLsq; break;
  case mEr: return mSEsq; break;
  default: 
    std::ostringstream ii;
    ii << "SoftPars::displaySoftMassSquared with illegal argument ";
    ii << k << endl;
    throw ii.str();
    break;
  }
}

template<class Susy, class Brevity>
double SoftPars<Susy, Brevity>::displaySoftMassSquared(softMasses k, int i, int j) 
  const {
  switch(k) {
  case mQl: return mQLsq.display(i, j); break;
  case mUr: return mURsq.display(i, j); break;
  case mDr: return mDRsq.display(i, j); break;
  case mLl: return mLLsq.display(i, j); break;
  case mEr: return mSEsq.display(i, j); break;
  default: 
    std::ostringstream ii;
    ii << "SoftPars::displaySoftMassSquared with illegal argument ";
    ii << k << endl;
    throw ii.str();
    break;
  }
}

template<class Susy, class Brevity>
const DoubleVector SoftPars<Susy, Brevity>::display() const {
  DoubleVector y(Susy::display());
  y.setEnd(numSoftParsMssm);
  int i, j, k=numSusyPars;
  for (i=1; i<=3; i++) {
    k++;
    y(k) = displayGaugino(i);
  }
  
  for (i=1; i<=3; i++)    
    for (j=1; j<=3; j++) {
      k++;
      y(k) = displayTrilinear(UA, i, j);
      y(k+9) = displayTrilinear(DA, i, j);
      y(k+18) = displayTrilinear(EA, i, j);
      y(k+27) = displaySoftMassSquared(mQl, i, j);
      y(k+36) = displaySoftMassSquared(mUr, i, j);
      y(k+45) = displaySoftMassSquared(mDr, i, j);
      y(k+54) = displaySoftMassSquared(mLl, i, j);
      y(k+63) = displaySoftMassSquared(mEr, i, j);
    }
  y(k+64) = m3sq;
  y(k+65) = mH1sq;
  y(k+66) = mH2sq;
  return y;
}

template<class Susy, class Brevity>
void SoftPars<Susy, Brevity>::setSoftMassElement(softMasses k, int i, int j, 
					  double f) { 
  switch(k) {
  case mQl: mQLsq(i, j) = f; break;
  case mUr: mURsq(i, j) = f; break;
  case mDr: mDRsq(i, j) = f; break;
  case mLl: mLLsq(i, j) = f; break;
  case mEr: mSEsq(i, j) = f; break;
  }
}

template<class Susy, class Brevity>
void SoftPars<Susy, Brevity>::setSoftMassMatrix(softMasses k, const DoubleMatrix & m) { 
  switch(k) {
  case mQl: mQLsq = m; break;
  case mUr: mURsq = m; break;
  case mDr: mDRsq = m; break;
  case mLl: mLLsq = m; break;
  case mEr: mSEsq = m; break;
  }
}

template<class Susy, class Brevity>
void SoftPars<Susy, Brevity>::setTrilinearMatrix(trilinears k, const DoubleMatrix & m) { 
  switch(k) {
  case UA: ua = m; break;
  case DA: da = m; break;
  case EA: ea = m; break;
  }
}

template<class Susy, class Brevity>
void SoftPars<Susy, Brevity>::setTrilinearElement(trilinears k, int i, int j, 
					  double m) { 
  switch(k) {
  case UA: ua(i, j) = m; break;
  case DA: da(i, j) = m; break;
  case EA: ea(i, j) = m; break;
  }
}

template<class Susy, class Brevity>
void SoftPars<Susy, Brevity>::setAllGauginos(const DoubleVector & v) { 
  if (v.displayStart() != 1 || v.displayEnd() !=3) {
    std::ostringstream ii;
    ii << "Initialising SoftPars::setAllGauginos with vector"
	 << v;
    throw ii.str();
  }
  mGaugino = v; 
}


template<class Susy, class Brevity>
void SoftPars<Susy, Brevity>::set(const DoubleVector & y) {
  Susy::set(y);
  int i, j, k=numSusyPars;
  for (i=1; i<=3; i++) {
    k++;
    mGaugino(i) = y.display(k);
  }
  
  for (i=1; i<=3; i++)
    for (j=1; j<=3; j++) {
      k++;
      ua(i, j) = y.display(k);
      da(i, j) = y.display(k+9);
      ea(i, j) = y.display(k+18);
      mQLsq(i, j) = y.display(k+27);
      mURsq(i, j) = y.display(k+36);
      mDRsq(i, j) = y.display(k+45);
      mLLsq(i, j) = y.display(k+54);
      mSEsq(i, j) = y.display(k+63);
    }
  m3sq = y.display(k+64);
  mH1sq = y.display(k+65);
  mH2sq = y.display(k+66);
}

template<class Susy, class Brevity>
SoftPars<Susy, Brevity> SoftPars<Susy, Brevity>::beta2() const {
  static Brevity a;
  return beta2(a);
}

// Outputs derivatives (DRbar scheme) in the form of dsoft
// thresholds = 0 and NOTHING is decoupled.
// thresholds = 2 and SUSY params/gauge couplings are decoupled at sparticle
// thresholds.
// CHECKED: 24/05/02
template<class Susy, class Brevity>
SoftPars<Susy, Brevity> SoftPars<Susy, Brevity>::beta2(Brevity& a) const {

  // Constants for gauge running
  static DoubleVector bBeta(3), cuBeta(3), cdBeta(3), ceBeta(3);
  static DoubleMatrix babBeta(3, 3);
  
  if (bBeta(1) < 1.0e-5) // Constants not set yet
    setBetas(babBeta, cuBeta, cdBeta, ceBeta, bBeta);
  
  // convert to beta functions
  static Susy dsb;
  
  // calculate derivatives for full SUSY spectrum. Brevity calculations come
  // out encoded in a
  dsb = Susy::beta(a);
  
  // To keep this a const function: TIME SAVINGS
  DoubleMatrix &u1=a.u1, &d1=a.d1, &e1=a.e1;
  DoubleMatrix &u2=a.u2, &d2=a.d2, &e2=a.e2,
    &u2t=a.u2t, &d2t=a.d2t, &e2t=a.e2t, &dt=a.dt, &ut=a.ut, &et=a.et;      
  double &uuT = a.uuT, &ddT = a.ddT, &eeT = a.eeT;
  DoubleVector &gsq=a.gsq;
  
  // Best to make these const DoubleMatrix & type
  const DoubleMatrix &hu = displayTrilinear(UA), &hd = displayTrilinear(DA),
    &he = displayTrilinear(EA);

  static DoubleMatrix hut(3, 3), hdt(3, 3), het(3, 3);
  static DoubleMatrix hu2(3, 3), hd2(3, 3), he2(3, 3);
  static DoubleMatrix hu2t(3, 3), hd2t(3, 3), he2t(3, 3);
  const DoubleMatrix &mq = displaySoftMassSquared(mQl),
    &mu = displaySoftMassSquared(mUr), &md = displaySoftMassSquared(mDr),
    &ml = displaySoftMassSquared(mLl), &me = displaySoftMassSquared(mEr); 
  static DoubleVector mG(displayGaugino()); 
  static DoubleVector msq(1, 3), gsqM(1, 3), gMsq(1, 3);
  
  hut = hu.transpose(); hdt = hd.transpose(); het = he.transpose();
  hu2 = hu * hut; hd2 = hd * hdt; he2 = he * het;
  double huT = hu2.trace(), hdT = hd2.trace(), heT = he2.trace();
  hu2t = hut * hu; hd2t = hdt * hd; he2t = het * he;
  double mqT = mq.trace(), muT = mu.trace(), mdT = md.trace(), meT =
    me.trace(), mlT = ml.trace();
  
  static DoubleMatrix huut(3,3), hddt(3,3), heet(3,3);
  huut = hu * ut; hddt = hd * dt; heet = he * et; // AVB
  double huuT = (hu * ut).trace(), hddT = (hd * dt).trace(), 
    heeT = (he * et).trace(); 
  mG = mGaugino; msq = mG * mG; gsqM = gsq * mG; gMsq = gsq * msq;
  
  double m3sq = displayM3Squared();

  // derivatives of soft parameters
  static DoubleVector dmG(1, 3);
  static double dmH1sq, dmH2sq, dm3sq;
  static DoubleMatrix dmq(3, 3), dmu(3, 3), dmd(3, 3), dme(3, 3), dml(3, 3);
  static DoubleMatrix dhu(3, 3), dhd(3, 3), dhe(3, 3);
  
  const double ONEO16Pisq = 1.0 / (16.0 * sqr(PI));
  if (displayLoops() > 0) {
    const static double sixteenO3 = 16.0 / 3.0, oneO15 = 1.0 /
      15.0;
    
    dmG = gsqM * 2.0 * bBeta;
    
    double curlyS = mH2sq - mH1sq + mqT - mlT - 2.0 * muT + mdT + meT;
    
    dm3sq = 2.0 * displaySusyMu() * (0.6 * gsqM(1) + 3.0 * gsqM(2) + 
			3.0 * huuT  +  3.0 * hddT + heeT) 
      + m3sq * (3.0 * uuT + 3.0 * ddT + eeT - 0.6 * gsq(1) - 3.0 * gsq(2));

    dmH1sq = 2.0 * 
      (-0.3 * gsq(1) * curlyS - 0.6 * gMsq(1) - 3.0 * gMsq(2) +
       3.0 * (mH1sq * ddT + (d2 * mq).trace() + (d2t * md).trace() + hdT) +
       (mH1sq * eeT + (e2 * ml).trace() + (e2t * me).trace() + heT));
    
    dmH2sq = 2.0 * 
      (0.3 * gsq(1) * curlyS - 0.6 * gMsq(1) - 3.0 * gMsq(2) +
       3.0 * (mH2sq * uuT + (u2 * mq).trace() + (u2t * mu).trace() + huT));

    dmq = 2.0 *
      (0.1 * gsq(1) * curlyS - oneO15 * gMsq(1) - 3.0 * gMsq(2) -
       sixteenO3 * gMsq(3) +
       0.5 * (u2 * mq + mq * u2 + 2.0 * (u1 * mu * ut + mH2sq * u2 + hu2))
       + 0.5 * (d2 * mq + mq * d2 + 2.0 * (d1 * md * dt + mH1sq * d2 +
					   hd2)));
    
    dml = 2.0 *
      ( -0.3 * gsq(1) * curlyS - 0.6 * gMsq(1) - 3.0 * gMsq(2) + 
	0.5 * (e2 * ml + ml * e2 + 2.0 * (e1 * me * et + mH1sq * e2 + he2)));
    
    dmu = 2.0 *
      (-0.4 * gsq(1) * curlyS - 16.0 * oneO15 * gMsq(1) - sixteenO3 *
       gMsq(3) + 
       (u2t * mu + mu * u2t + 2.0 * (ut * mq * u1 + mH2sq * u2t + hu2t)));
    
    dmd = 2.0 * (0.2 * gsq(1) * curlyS 
		 - 4.0 * oneO15 * gMsq(1) - sixteenO3 * gMsq(3) + 
		 (d2t * md + md * d2t + 2.0 * (dt * mq * d1 + mH1sq * d2t +
					       hd2t)));
    
    dme = 2.0 *
      (0.6 * gsq(1) * curlyS - 36.0 * oneO15 * gMsq(1) +
       e2t * me + me * e2t + 2.0 * (et * ml * e1 + mH1sq * e2t + he2t));
    
    dhu = (3.0 * uuT 
	   - sixteenO3 * gsq(3) - 3.0 * gsq(2) - 13.0 * oneO15 * gsq(1))
      * hu
      + 2.0 * (3.0 * huuT + 13.0 * oneO15 * gsqM(1) + 3 * gsqM(2) + sixteenO3
	     * gsqM(3)) * u1 +
      4.0 * hu * u2t + 5.0 * u2 * hu + 2.0 * hd * dt * u1 + d2 * hu;
    
    dhd = (eeT + 3.0 * ddT -sixteenO3 * gsq(3) - 3.0 * gsq(2) - 
	   7.0 * oneO15 * gsq(1)) * hd 
      + 2.0 * (7.0 * oneO15 * gsqM(1) + 3.0 * gsqM(2) + sixteenO3
	     * gsqM(3) + heeT + 3.0 * hddT) * d1 +
      4.0 * hd * d2t + 5.0 * d2 * hd + 
      2.0 * hu * ut * d1 + u2 * hd;
    
    dhe = (eeT + 3.0 * ddT - 3.0 * gsq(2) - 1.8 * gsq(1)) * he
      + (3.6 * gsqM(1) + 6.0 * gsqM(2) + 2.0 * heeT + 6.0 * hddT) * e1 +
      4.0 * he * e2t + 5.0 * e2 * he;

    // convert to proper derivatives: 
    dmG    = ONEO16Pisq * dmG;
    dm3sq  *= ONEO16Pisq;
    dmH1sq *= ONEO16Pisq;
    dmH2sq *= ONEO16Pisq;
    dmq    *= ONEO16Pisq;
    dml    *= ONEO16Pisq;
    dmu    *= ONEO16Pisq;
    dmd    *= ONEO16Pisq;
    dme    *= ONEO16Pisq;
    dhu    *= ONEO16Pisq;
    dhd    *= ONEO16Pisq;
    dhe    *= ONEO16Pisq;
  }
  
  // two-loop contributions. I got these from hep-ph/9311340. WIth respect to
  // their notation, the Yukawas and h's are TRANSPOSED. Gaugino masses are
  // identical, as are soft masses. B(mV) = mu(BBO) B(BBO)
  if (displayLoops() > 1) {
    static DoubleVector dmG2(1, 3);
    static DoubleMatrix dmq2(3, 3), dmu2(3, 3), dmd2(3, 3), dme2(3, 3), 
      dml2(3, 3);
    static DoubleMatrix dhu2(3, 3), dhd2(3, 3), dhe2(3, 3);
    static DoubleVector sigma(3);

    const double oneO16Pif = sqr(ONEO16Pisq);
    
    DoubleVector &g4=a.g4;

    double mq3 = displaySoftMassSquared(mQl, 3, 3);
    double mu3 = displaySoftMassSquared(mUr, 3, 3);
    double md3 = displaySoftMassSquared(mDr, 3, 3);
    double ml3 = displaySoftMassSquared(mLl, 3, 3);
    double me3 = displaySoftMassSquared(mEr, 3, 3);
    double ht = displayYukawaElement(YU, 3, 3), ht2 = sqr(ht), ht4 = sqr(ht2);
    double htau = displayYukawaElement(YE, 3, 3), htau2 = sqr(htau), 
      htau4 = sqr(htau2);
    double hb = displayYukawaElement(YD, 3, 3), hb2 = sqr(hb), hb4 = sqr(hb2);
    double Ut = displayTrilinear(UA, 3, 3), Ut2 = sqr(Ut);
    double Ub = displayTrilinear(DA, 3, 3), Ub2 = sqr(Ub);
    double Utau = displayTrilinear(EA, 3, 3), Utau2 = sqr(Utau);

    double sP;
    if (MIXING < 1) {
      sP = ///< dominant 3rd family approximation
	-(3.0 * mH2sq + mq3 - 4.0 * mu3) * ht2 +  
	(3.0 * mH1sq - mq3 - 2.0 * md3) * hb2 + 
	(mH1sq + ml3 - 2.0 * me3) * htau2 +
	(1.5 * gsq(2) + 0.3 * gsq(1)) * (mH2sq - mH1sq - mlT) +
	(8.0 / 3.0 * gsq(3) + 1.5 * gsq(2) + gsq(1) / 30.0) * mqT -
	(16.0 / 3.0 * gsq(3) + 16.0 / 15.0 * gsq(1)) * muT +
	(8.0 / 3.0 * gsq(3) + 2.0 / 15.0 * gsq(1)) * mdT +       
	1.2 * gsq(1) * meT;  //checked
    } else {
      sP = ///< fully mixed case -- slower
	-(3.0 * mH2sq * uuT + (mq * u2).trace()) 
	+ 4.0 * (u1 * mu * ut).trace() + 
	3.0 * mH1sq * ddT - (mq * d2).trace() - 2.0 * (d1 * md * dt).trace() + 
	mH1sq * eeT + (ml * e2).trace() -
	(2.0 * e1 * me * et).trace() +
	(1.5 * gsq(2) + 0.3 * gsq(1)) * (mH2sq - mH1sq - mlT) +
	(8.0 / 3.0 * gsq(3) + 1.5 * gsq(2) + gsq(1) / 30.0) * mqT -
	(16.0 / 3.0 * gsq(3) + 16.0 / 15.0 * gsq(1)) * muT +
	(8.0 / 3.0 * gsq(3) + 2.0 / 15.0 * gsq(1)) * mdT +       
	1.2 * gsq(1) * meT;  //checked
    }
    
    sigma(1) = 0.2 * gsq(1) *
      (3.0 * (mH1sq + mH2sq) + mqT + 3.0 * mlT + 8.0 * muT + 2.0 * mdT + 6.0 *
       meT);
    sigma(2) = gsq(2) * (mH1sq + mH2sq + (3.0 * mqT + mlT));
    sigma(3) = gsq(3) * (2.0 * mqT + muT + mdT); // these 3 checked

    DoubleMatrix u4(u2 * u2), d4(d2 * d2), u4t(u2t * u2t), d4t(d2t * d2t);

    if (INCLUDE_2_LOOP_SCALAR_CORRECTIONS) {
    /* new dominant 3-family version */
      if (MIXING < 1) {
	dmq2 = -(2.0 * mq + 8.0 * mH2sq) * u4 - 
	  4.0 * u1 * mu * u2t * ut - 
	  4.0 * u2 * mq * u2 - 4.0 * u2 * u1 * mu * ut - 2.0 * u4 * mq -
	  (2.0 * mq + 8.0 * mH1sq) * d4 - 4.0 * d1 * md * d2t * dt - 
	  4.0 * d2 * mq * d2 - 4.0 * d2 * d1 * md * dt - 2.0 * d4 * mq -
	  ((mq + 4.0 * mH2sq) * u2 + 2.0 * u1 * mu * ut + u2 * mq) * uuT 
	  * 3.0 -
	  ((mq + 4.0 * mH1sq) * d2 + 2.0 * d1 * md * dt + d2 * mq) * 
	  (3.0 * ddT + eeT) -
	  6.0 * u2 * ((mq * u2).trace() + (u1 * mu * ut).trace()) -
	  d2 * ((6.0 * mq * d2).trace() + (6.0 * d1 * md * dt).trace() + 
		(2.0 * ml * e2).trace() + (2.0 * e1 * me * et).trace()) -
	  4.0 * (u2 * hu2 + hu2 * u2 + u1 * hu2t * ut + hu * u2t * hut) -
	  4.0 * (d2 * hd2 + hd2 * d2 + d1 * hd2t * dt + hd * d2t * hdt) -
	  hu2 * uuT * 6.0 - u2 * 6.0 * Ut2 - hu * ut * 6.0 * huuT -
	  u1 * hut * 6.0 * huuT - hd2 * (6.0 * ddT + 2.0 * eeT) - 
	  d2 * (6.0 * Ub2 + 2.0 * Utau2) -
	  hd * dt * (6.0 * hddT + 2.0 * heeT) - 
	  d1 * hdt * (6.0 * hddT + 2.0 * heeT) + 0.4 * gsq(1) * 
	  ((2.0 * mq +  4.0 * mH2sq) * u2 + 4.0 * u1 * mu * ut +
	   2.0 * u2 * mq + 4.0 * hu2 - 4.0 * mG(1) * hu * ut - 4.0
	   * mG(1) * u1 * hut + 8.0 * msq(1) * u2 +
	   (mq + 2.0 * mH1sq) * d2 + 2.0 * d1 * md * dt +
	   d2 * mq + 2.0 * hd2 - 2.0 * mG(1) * hd * dt - 
	   2.0 * mG(1) * d1 * hdt + 4.0 * msq(1) * d2) +
	  0.4 * gsq(1) * sP - 128.0 / 3.0 * g4(3) * msq(3) + 32.0 * gsq(3) *
	  gsq(2) * (msq(3) + msq(2) + mG(2) * mG(3)) +
	  32.0 / 45.0 * gsq(3) * gsq(1) * 
	  (msq(3) + msq(1) + mG(1) * mG(3)) + 33.0
	  * g4(2) * msq(2) + 0.4 * gsq(2) * gsq(1) * (msq(2) + msq(1) + mG(1) *
						      mG(2)) +
	  199.0 / 75.0 * g4(1) * msq(1) + 16.0 / 3.0 * gsq(3) * sigma(3) +
	  3.0 * gsq(2) * sigma(2) + gsq(1) * sigma(1) / 15.0; // checked

	dml2 = -(2.0 * ml + 8.0 * mH1sq) * e2 * e2 - 4.0 * e1 * me * e2t * et -
	  4.0 * e2 * ml * e2 - 4.0 * e2 * e1 * me * et - 2.0 * e2 * e2 * ml -
	  ((ml + 4.0 * mH1sq) * e2 + 2.0 * e1 * me * et + e2 * ml) * 
	  (3.0 * ddT + eeT) -
	  e2 * (6.0 * mq * d2 + 6.0 * d1 * md * dt + 2.0 * ml * e2 + 
		2.0 * e1 * me * et).trace() -
	  4.0 * (e2 * he2 + he2 * e2 + e1 * he2t * et + he * e2t * het) -
	  he2 * (6.0 * ddT + 2.0 * eeT) - e2 * (6.0 * Ub2 + 2.0 * Utau2) -
	  he * et * (6.0 * hddT + 2.0 * heeT) - 
	  e1 * het * (6.0 * hddT + 2.0 * heeT) +
	  1.2 * gsq(2) * 
	  ((ml + 2.0 * mH1sq) * e2 + 2.0 * e1 * me * et + e2 * ml +
	   2.0 * he2 - 2.0 * mG(1) * he * et - 2.0 * mG(1) * e1 *
			  het + 4.0 * msq(1) * e2) -
	  1.2 * gsq(1) * sP + 33.0 * g4(2) * msq(2) + 3.6 * gsq(2) * gsq(1) *
	  (msq(2) + msq(1) + mG(1) * mG(2)) + 621.0 / 25.0 * g4(1) * msq(1) +
	  3.0 * gsq(2) * sigma(2) + 0.6 * gsq(1) * sigma(1); //checked

	dmu2 = -(2.0 * mu + 8.0 * mH2sq) * u4t - 4.0 * ut * mq * u2 * u1 
	  - 4.0 * u2t * mu * u2t - 4.0 * u2t * ut * mq * u1 - 
	  2.0 * u4t * mu - (2.0 * mu + 4.0 * mH2sq + 4.0 * mH1sq) 
	  * ut * d2 * u1 - 4.0 * ut * mq * d2 * u1 - 
	  4.0 * ut * d1 * md * dt * u1 -
	  4.0 * ut * d2 * mq * u1 - 2.0 * ut * d2 * u1 * mu -
	  ((mu + 4.0 * mH2sq) * u2t + 2.0 * ut * mq * u1 + u2t * mu) * 
	  6.0 * uuT -
	  12.0 * u2t * (mq * u2 + u1 * mu * ut).trace() -
	  4.0 * (hu2t * u2t + u2t * hu2t + hut * u2 * hu + ut * hu2 * u1) -
	  4.0 * (hut * hd * dt * u1 + ut * d1 * hdt * hu + hut * d2 * hu + 
	     ut * hd2 * u1) -
	  12.0 * (hu2t * uuT + u2t * Ut2 + hut * u1 * huuT +
		  ut * hu * huuT) +
	  (6.0 * gsq(2) - 0.4 * gsq(1)) * 
	  ((mu + 2.0 * mH2sq) * u2t + 2.0 * ut * mq * u1 + u2t * mu + 
	   2.0 * hu2t) +
	  12.0 * gsq(2) * (2.0 * msq(2) * u2t - mG(2) * hut * u1 - 
			   mG(2) * ut * hu) -
	  0.8 * gsq(1) * (2.0 * msq(1) * u2t - mG(1) * hut * u1 - 
			  mG(1) * ut * hu)-
	  1.6 * gsq(1) * sP - 128.0 / 3.0 * g4(3) * msq(3) + 
	  512.0 / 45.0 * gsq(1) * gsq(3) * (msq(3) + msq(1) + mG(1) * mG(3)) +
	  3424.0 / 75.0 * g4(1) * msq(1) + 16.0 / 3.0 * gsq(3) * sigma(3) +
	  16.0 / 15.0 * gsq(1) * sigma(1); // checked

	dmd2 = -(2.0 * md + 8.0 * mH1sq) * d4t - 
	  4.0 * dt * mq * d1 * d2t -
	  4.0 * d2t * md * d2t - 4.0 * d2t * dt * mq * d1 - 
	  2.0 * d4t * md -
	  (2.0 * md + 4.0 * mH2sq + 4.0 * mH1sq) * dt * u2 * d1 -
	  4.0 * dt * mq * u2 * d1 - 4.0 * dt * u1 * mu * ut * d1 - 
	  4.0 * dt * u2 * mq * d1 - 2.0 * dt * u2 * d1 * md - 
	  ((md + 4.0 * mH1sq) * d2t + 2.0 * dt * mq * d1 + d2t * md) * 
	  (6.0 * ddT + 2.0 * eeT) -
	  4.0 * d2t * (3.0 * mq * d2 + 3.0 * d1 * md * dt + ml * e2 
		       + e1 * me * et).trace() -
	  4.0 * (hd2t * d2t + d2t * hd2t + hdt * d2 * hd + dt * hd2 * d1) -
	  4.0 * (hdt * hu * ut * d1 + dt * u1 * hut * hd + hdt * u2 * hd + 
		 dt * hu2 * d1) -
	  4.0 * hd2t * (3.0 * ddT + eeT) - 4.0 * d2t * (3.0 * hdT + heT) -
	  4.0 * hdt * d1 * (3.0 * hddT + heeT) - 
	  4.0 * dt * hd * (3.0 * hddT + heeT) +
	  (6.0 * gsq(2) + 0.4 * gsq(1)) * 
	  ((md + 2.0 * mH1sq) * d2t + 2.0 * dt * mq * d1 + d2t * md + 
	   2.0 * hd2t) +
	  12.0 * gsq(2) * (2.0 * msq(2) * d2t - mG(2) * hdt * d1 
			   - mG(2) * dt * hd) +
	  0.8 * gsq(1) * (2.0 * msq(1) * d2t - mG(1) * hdt * d1 - 
			  mG(1) * dt * hd)+
	  0.8 * gsq(1) * sP - 128.0 / 3.0 * g4(3) * msq(3) + 
	  128.0 / 45.0 * gsq(3) * gsq(1) * (msq(3) + msq(1) + mG(1) * mG(3)) +
	  808.0 / 75.0 * g4(1) * msq(1) + 16.0 / 3.0 * gsq(3) * sigma(3) +
	  4.0 / 15.0 * gsq(1) * sigma(1); // checked

	dme2 = -(2.0 * me + 8.0 * mH1sq) * e2t * e2t - 
	  4.0 * et * ml * e2 * e1 -
	  4.0 * e2t * me * e2t - 4.0 * e2t * et * ml * e1 - 
	  2.0 * e2t * e2t * me -
	  ((me + 4.0 * mH1sq) * e2t + 2.0 * et * ml * e1 + e2t * me) *
	  (6.0 * ddT + 2.0 * eeT) - 4.0 * e2t * 
	  (3.0 * mq * d2 + 3.0 * d1 * md * dt + ml * e2 + 
	   e1 * me * et).trace() -
	  4.0 * (he2t * e2t + e2t * he2t + het * e2 * he + et * he2 * e1) -
	  4.0 * he2t * (3.0 * ddT + eeT) - 4.0 * e2t * (3.0 * Ub2 + Utau2) -
	  4.0 * het * e1 * (3.0 * hddT + heeT) - 4.0 * et * he * 
	  (3.0 * hddT + heeT) + (6.0 * gsq(2) - 1.2 * gsq(1)) * 
	  ((me + 2.0 * mH1sq) * e2t + 2.0 * et * ml * e1 + e2t * me + 
	   2.0 * he2t) +
	  12.0 * gsq(2) * (2.0 * msq(2) * e2t - mG(2) * het * e1 
			   - mG(2) * et * he) -
	  2.4 * gsq(1) * (2.0 * msq(1) * e2t - mG(1) * het * e1 - 
			  mG(1) * et * he)
	  + 2.4 * gsq(1) * sP + 2808. / 25.0 * g4(1) * msq(1) + 
	  2.4 * gsq(1) * sigma(1); // checked

	dhu2 = 
	  (-3.0 * (3.0 * ht4 + ht2 * hb2) - d2 * (3.0 * ddT + eeT)
	   -15.0 * u2 * uuT - 6.0 * u4 - 2.0 * d4 -
	   4.0 * u2 * d2 + (16.0 * gsq(3) + 0.8 * gsq(1)) * uuT +
	   12.0 * gsq(2) * u2 + 0.4 * gsq(1) * d2 - 16.0 / 9.0 * g4(3) +
	   8.0 * gsq(3) * gsq(2) + 136.0 / 45.0 * gsq(3) * gsq(1) + 
	   7.5 * g4(2) + 
	   gsq(2) * gsq(1) + 2743.0 / 450.0 * g4(1)) * hu +
	  (-6.0 * (6.0 * Ut * ht2 * ht + Ut * hb2 * ht + Ub * ht2 * hb) -
	   18.0 * u2 * huuT - d2 * (6.0 * hddT + 2.0 * heeT) - 
	   12.0 * hu * ut * uuT - hd * dt * (6.0 * ddT + 2.0 * eeT) - 
	   6.0 * hu * u2t * ut - 8.0 * u2 * hu * ut - 4.0 * hd * d2t * dt -
	   4.0 * d2 * hd * dt - 2.0 * hu * ut * d2 - 4.0 * u2 * hd * dt +
	   (32.0 * gsq(3) + 1.6 * gsq(1)) * huuT + 
	   (6.0 * gsq(2) + 1.2 * gsq(1)) * hu * ut + 0.8 * gsq(1) * hd * dt -
	   (32.0 * gsq(3) * mG(3) + 1.6 * gsq(1) * mG(1)) * uuT -
	   (12.0 * gsq(2) * mG(2) + 0.8 * gsq(1) * mG(1)) * u2 -
	   0.8 * gsq(1) * mG(1) * d2 + 64.0 / 9.0 * g4(3) * mG(3) - 
	   16.0 * gsq(3) * gsq(2) * (mG(3) + mG(2)) -
	   272.0 / 45.0 * gsq(3) * gsq(1) * (mG(1) + mG(3)) -
	   30.0 * g4(2) * mG(2) - 2.0 * gsq(2) * gsq(1) * (mG(2) + mG(1)) -
	   5486.0 / 225.0 * g4(1) * mG(1)) * u1; // checked

	dhd2 =  // checked + debugged 19.05.04 by M R Ramage (thanks Mike!)
	  (-3.0 * (3.0 * hb4 + ht2 * hb2 + htau4) -
	   3.0 * u2 * uuT - 5.0 * d2 * (3.0 * ddT + eeT) - 6.0 * d4 -
	   2.0 * u4 - 4.0 * d2 * u2 + (16.0 * gsq(3) - 
					    0.4 * gsq(1)) * ddT +
	   1.2 * gsq(1) * eeT + 0.8 * gsq(1) * u2 + 
	   (12.0 * gsq(2) + 1.2 * gsq(1)) * d2 - 16.0 / 9.0 * g4(3) + 
	   8.0 * gsq(3) * gsq(2) + 8.0 / 9.0 * gsq(3) * gsq(1) + 7.5 * g4(2) +
	   gsq(2) * gsq(1) + 287.0 / 90.0 * g4(1)) * hd +
	  (-6.0 * (6.0 * Ub * hb2 * hb + Ut * hb2 * ht + Ub * ht2 * hb +
		   2.0 * Utau * htau2 * htau) - 6.0 * u2 * huuT -
	   6.0 * d2 * (3.0 * hddT + heeT) - 6.0 * hu * ut * uuT -
	   4.0 * hd * dt * (3.0 * ddT + eeT) - 6.0 * hd * d2t * dt -
	   8.0 * d2 * hd * dt - 4.0 * u2 * hu * ut - 4.0 * hu * u2t * ut -
	   4.0 * d2 * hu * ut - 2.0 * hd * dt * u2 + 
	   (32.0 * gsq(3) - 0.8 * gsq(1)) * hddT + 2.4 * gsq(1) * heeT +
	   1.6 * gsq(1) * hu * ut + (6.0 * gsq(2) + 1.2 * gsq(1)) * hd * dt -
	   (32.0 * gsq(3) * mG(3) - 0.8 * gsq(1) * mG(1)) * ddT - 
	   2.4 * gsq(1) * mG(1) * eeT - 
	   (12.0 * gsq(2) * mG(2) + 1.6 * gsq(1) * mG(1)) * d2 -
	   1.6 * gsq(1) * mG(1) * u2 + 64.0 / 9.0 * g4(3) * mG(3) - 
	   16.0 * gsq(3) * gsq(2) * (mG(3) + mG(2)) - 
	   16.0 / 9.0 * gsq(3) * gsq(1) * (mG(3) + mG(1)) - 
	   30.0 * g4(2) * mG(2) - 
	   2.0 * gsq(2) * gsq(1) * (mG(2) + mG(1)) - 
	   574.0 / 45.0 * g4(1) * mG(1)) * d1;

	dhe2 = 
	  (-3.0 * (3.0 * hb4 + ht2 * hb2 + htau4) -
	   5.0 * e2 * (3.0 * ddT + eeT) - 6.0 * e2 * e2 + 
	   (16.0 * gsq(3) - 0.4 * gsq(1)) * ddT + 1.2 * gsq(1) * eeT +
	   (12.0 * gsq(2) - 1.2 * gsq(1)) * eeT + 7.5 * g4(2) + 
	   1.8 * gsq(2) * gsq(1) + 13.5 * g4(1)) * he +
	  (-6.0 * (6.0 * Ub * hb2 * hb + Ut * hb2 * ht + Ub * ht2 * hb +
		   2.0 * Utau * htau2 * htau) - 
	   4.0 * he * et * (3.0 * ddT + eeT) - 6.0 * e2 * (3.0 * hddT + heeT) -
	   6.0 * he * e2t * et - 8.0 * e2 * he * et + 
	   (32.0 * gsq(3) - 0.8 * gsq(1)) * hddT + 2.4 * gsq(1) * heeT + 
	   (6.0 * gsq(2) + 1.2 * gsq(1)) * he * et -
	   (32.0 * gsq(3) * mG(3) - 0.8 * gsq(1) * mG(1)) * ddT - 
	   2.4 * gsq(1) * mG(1) * eeT - 12.0 * gsq(2) * mG(2) * e2 -
	   30.0 * g4(2) * mG(2) - 3.6 * gsq(2) * gsq(1) * (mG(1) + mG(2)) -
	   54.0 * g4(1) * mG(1)) * e1; // checked
      } else {
	// NB you should introduce speedy computation here. For example sum
	// traces, not add then trace. Also surely some of the products are
	// repeated, eg u2 * u2 = u4, d2 * d2 = d4, u1 * mu * ut, d1 * md *
	// dt, d2 * d1, d2 * mq
	dmq2 = -(2.0 * mq + 8.0 * mH2sq) * u4 - 4.0 * u1 * mu * u2t * ut 
	  - 4.0 * u2 * mq * u2 - 4.0 * u2 * u1 * mu * ut - 2.0 * u4 * mq -
	  (2.0 * mq + 8.0 * mH1sq) * d4 - 4.0 * d1 * md * d2t * dt - 
	  4.0 * d2 * mq * d2 - 4.0 * d2 * d1 * md * dt - 2.0 * d4 * mq -
	  ((mq + 4.0 * mH2sq) * u2 + 2.0 * u1 * mu * ut + u2 * mq) * uuT * 3.0 
	  -((mq + 4.0 * mH1sq) * d2 + 2.0 * d1 * md * dt + d2 * mq) * 
	  (3.0 * ddT + eeT) -
	  6.0 * u2 * ((mq * u2).trace() + (u1 * mu * ut).trace()) -
	  d2 * (6.0 * (mq * d2).trace() + 6.0 * (d1 * md * dt).trace() + 
		2.0 * (ml * e2).trace() + 2.0 * (e1 * me * et).trace()) -
	  4.0 * (u2 * hu2 + hu2 * u2 + u1 * hu2t * ut + hu * u2t * hut) -
	  4.0 * (d2 * hd2 + hd2 * d2 + d1 * hd2t * dt + hd * d2t * hdt) -
	  hu2 * uuT * 6.0 - u2 * 6.0 * huT - 
	  hu * ut * 6.0 * huuT - u1 * hut * 6.0 * huuT - hd2 * 
	  (6.0 * ddT + 2.0 * eeT) - 
	  d2 * (6.0 * hdT + 2.0 * heT) -
	  hd * dt * (6.0 * hddT + 2.0 * heeT) - 
	  d1 * hdt * (6.0 * hddT + 2.0 * heeT) +
	  0.4 * gsq(1) * 
	  ((2.0 * mq +  4.0 * mH2sq) * u2 + 4.0 * u1 * mu * ut +
	   2.0 * u2 * mq + 4.0 * hu2 - 4.0 * mG(1) * hu * ut - 4.0
	   * mG(1) * u1 * hut + 8.0 * msq(1) * u2 +
	   (mq + 2.0 * mH1sq) * d2 + 2.0 * d1 * md * dt +
	   d2 * mq + 2.0 * hd2 - 2.0 * mG(1) * hd * dt - 
	   2.0 * mG(1) * d1 * hdt + 4.0 * msq(1) * d2) +
	  0.4 * gsq(1) * sP - 128.0 / 3.0 * g4(3) * msq(3) + 32.0 * gsq(3) *
	  gsq(2) * (msq(3) + msq(2) + mG(2) * mG(3)) +
	  32.0 / 45.0 * gsq(3) * gsq(1) * (msq(3) + msq(1) + 
					   mG(1) * mG(3)) + 33.0
	  * g4(2) * msq(2) + 0.4 * gsq(2) * gsq(1) * (msq(2) + msq(1) + mG(1) *
						      mG(2)) +
	  199.0 / 75.0 * g4(1) * msq(1) + 16.0 / 3.0 * gsq(3) * sigma(3) +
	  3.0 * gsq(2) * sigma(2) + gsq(1) * sigma(1) / 15.0; // checked

	dml2 = -(2.0 * ml + 8.0 * mH1sq) * e2 * e2 - 4.0 * e1 * me * e2t * et -
	  4.0 * e2 * ml * e2 - 4.0 * e2 * e1 * me * et - 2.0 * e2 * e2 * ml -
	  ((ml + 4.0 * mH1sq) * e2 + 2.0 * e1 * me * et + e2 * ml) * 
	  (3.0 * ddT + eeT) -
	  e2 * (6.0 * mq * d2 + 6.0 * d1 * md * dt + 2.0 * ml * e2 + 
		2.0 * e1 * me * et).trace() -
	  4.0 * (e2 * he2 + he2 * e2 + e1 * he2t * et + he * e2t * het) -
	  he2 * (6.0 * ddT + 2.0 * eeT) - 
	  e2 * (6.0 * hdT + 2.0 * heT) -
	  he * et * (6.0 * hddT + 2.0 * heeT) - 
	  e1 * het * (6.0 * hddT + 2.0 * heeT) +
	  1.2 * gsq(1) * ((ml + 2.0 * mH1sq) * e2 + 2.0 * e1 * me * et + 
			  e2 * ml +
			  2.0 * he2 - 2.0 * mG(1) * he * et - 
			  2.0 * mG(1) * e1 *
			  het + 4.0 * msq(1) * e2) -
	  1.2 * gsq(1) * sP + 33.0 * g4(2) * msq(2) + 3.6 * gsq(2) * gsq(1) *
	  (msq(2) + msq(1) + mG(1) * mG(2)) + 621.0 / 25.0 * g4(1) * msq(1) +
	  3.0 * gsq(2) * sigma(2) + 0.6 * gsq(1) * sigma(1); //checked
	
	dmu2 = -(2.0 * mu + 8.0 * mH2sq) * u4t - 
	  4.0 * ut * mq * u2 * u1 -
	  4.0 * u2t * mu * u2t - 4.0 * u2t * ut * mq * u1 - 
	  2.0 * u4t * mu - (2.0 * mu + 4.0 * mH2sq + 
				  4.0 * mH1sq) * ut * d2
	  * u1 - 4.0 * ut * mq * d2 * u1 - 4.0 * ut * d1 * md * dt * u1 -
	  4.0 * ut * d2 * mq * u1 - 2.0 * ut * d2 * u1 * mu -
	  ((mu + 4.0 * mH2sq) * u2t + 2.0 * ut * mq * u1 + 
	   u2t * mu) * 6.0 * uuT -
	  12.0 * u2t * (mq * u2 + u1 * mu * ut).trace() -
	  4.0 * (hu2t * u2t + u2t * hu2t + hut * u2 * hu + ut * hu2 * u1) -
	  4.0 * (hut * hd * dt * u1 + ut * d1 * hdt * hu + hut * d2 * hu + 
		 ut * hd2 * u1) -
	  12.0 * (hu2t * uuT + u2t * huT + hut * u1 * huuT +
		  ut * hu * huuT) +
	  (6.0 * gsq(2) - 0.4 * gsq(1)) * 
	  ((mu + 2.0 * mH2sq) * u2t + 2.0 * ut * mq * u1 + u2t * mu + 
	   2.0 * hu2t) +
	  12.0 * gsq(2) * (2.0 * msq(2) * u2t - mG(2) * hut * u1 - 
			   mG(2) * ut * hu) -
	  0.8 * gsq(1) * (2.0 * msq(1) * u2t - mG(1) * hut * u1 - 
			  mG(1) * ut * hu)-
	  1.6 * gsq(1) * sP - 128.0 / 3.0 * g4(3) * msq(3) + 
	  512.0 / 45.0 * gsq(1) * gsq(3) * (msq(3) + msq(1) + mG(1) * mG(3)) +
	  3424.0 / 75.0 * g4(1) * msq(1) + 16.0 / 3.0 * gsq(3) * sigma(3) +
	  16.0 / 15.0 * gsq(1) * sigma(1); // checked

	dmd2 = -(2.0 * md + 8.0 * mH1sq) * d4t - 
	  4.0 * dt * mq * d1 * d2t -
	  4.0 * d2t * md * d2t - 4.0 * d2t * dt * mq * d1 - 
	  2.0 * d4t * md -
	  (2.0 * md + 4.0 * mH2sq + 4.0 * mH1sq) * dt * u2 * d1 -
	  4.0 * dt * mq * u2 * d1 - 4.0 * dt * u1 * mu * ut * d1 - 
	  4.0 * dt * u2 * mq * d1 - 2.0 * dt * u2 * d1 * md - 
	  ((md + 4.0 * mH1sq) * d2t + 2.0 * dt * mq * d1 + d2t * md) * 
	  (6.0 * ddT + 2.0 * eeT) -
	  4.0 * d2t * (3.0 * mq * d2 + 3.0 * d1 * md * dt + ml * e2 
		       + e1 * me * et).trace() -
	  4.0 * (hd2t * d2t + d2t * hd2t + hdt * d2 * hd + dt * hd2 * d1) -
	  4.0 * (hdt * hu * ut * d1 + dt * u1 * hut * hd + hdt * u2 * hd + 
		 dt * hu2 * d1) -
	  4.0 * hd2t * (3.0 * ddT + eeT) - 4.0 * d2t * (3.0 * hdT + heT) -
	  4.0 * hdt * d1 * (3.0 * hddT + heeT) - 
	  4.0 * dt * hd * (3.0 * hddT + heeT) +
	  (6.0 * gsq(2) + 0.4 * gsq(1)) * 
	  ((md + 2.0 * mH1sq) * d2t + 2.0 * dt * mq * d1 + d2t * md + 
	   2.0 * hd2t) +
	  12.0 * gsq(2) * (2.0 * msq(2) * d2t - mG(2) * hdt * d1 
			   - mG(2) * dt * hd) +
	  0.8 * gsq(1) * (2.0 * msq(1) * d2t - mG(1) * hdt * d1 - 
			  mG(1) * dt * hd)+
	  0.8 * gsq(1) * sP - 128.0 / 3.0 * g4(3) * msq(3) + 
	  128.0 / 45.0 * gsq(3) * gsq(1) * (msq(3) + msq(1) + mG(1) * mG(3)) +
	  808.0 / 75.0 * g4(1) * msq(1) + 16.0 / 3.0 * gsq(3) * sigma(3) +
	  4.0 / 15.0 * gsq(1) * sigma(1); // checked

	dme2 = -(2.0 * me + 8.0 * mH1sq) * e2t * e2t - 
	  4.0 * et * ml * e2 * e1 -
	  4.0 * e2t * me * e2t - 4.0 * e2t * et * ml * e1 - 
	  2.0 * e2t * e2t * me -
	  ((me + 4.0 * mH1sq) * e2t + 2.0 * et * ml * e1 + e2t * me) *
	  (6.0 * ddT + 2.0 * eeT) - 4.0 * e2t * 
	  (3.0 * mq * d2 + 3.0 * d1 * md * dt + ml * e2 + 
	   e1 * me * et).trace() -
	  4.0 * (he2t * e2t + e2t * he2t + het * e2 * he + et * he2 * e1) -
	  4.0 * he2t * (3.0 * ddT + eeT) - 4.0 * e2t * (3.0 * hdT + heT) -
	  4.0 * het * e1 * (3.0 * hddT + heeT) - 4.0 * et * he * 
	  (3.0 * hddT + heeT) + (6.0 * gsq(2) - 1.2 * gsq(1)) * 
	  ((me + 2.0 * mH1sq) * e2t + 2.0 * et * ml * e1 + e2t * me + 
	   2.0 * he2t) +
	  12.0 * gsq(2) * (2.0 * msq(2) * e2t - mG(2) * het * e1 
			   - mG(2) * et * he) -
	  2.4 * gsq(1) * (2.0 * msq(1) * e2t - mG(1) * het * e1 - 
			  mG(1) * et * he)
	  + 2.4 * gsq(1) * sP + 2808. / 25.0 * g4(1) * msq(1) + 
	  2.4 * gsq(1) * sigma(1); // checked

    /* Old full 3-family version 
    dhd2 =  // checked + debugged 19.05.04 by M R Ramage (thanks Mike!)
      (-3.0 * (3.0 * d2t * d2t + ut * d2 * u1 + e2t * e2t).trace() -
       3.0 * u2 * uuT - 5.0 * d2 * (3.0 * ddT + eeT) - 6.0 * d2 * d2 -
       2.0 * u2 * u2 - 4.0 * d2 * u2 + (16.0 * gsq(3) - 0.4 * gsq(1)) * ddT +
       1.2 * gsq(1) * eeT + 0.8 * gsq(1) * u2 + 
       (12.0 * gsq(2) + 1.2 * gsq(1)) * d2 - 16.0 / 9.0 * g4(3) + 
       8.0 * gsq(3) * gsq(2) + 8.0 / 9.0 * gsq(3) * gsq(1) + 7.5 * g4(2) +
       gsq(2) * gsq(1) + 287.0 / 90.0 * g4(1)) * hd +
      (-6.0 * (6.0 * hdt * d2 * d1 + hut * d2 * u1 + hdt * u2 * d1 +
	       2.0 * het * e2 * e1).trace() - 6.0 * u2 * huuT -
       6.0 * d2 * (3.0 * hddT + heeT) - 6.0 * hu * ut * uuT -
       4.0 * hd * dt * (3.0 * ddT + eeT) - 6.0 * hd * d2t * dt -
       8.0 * d2 * hd * dt - 4.0 * u2 * hu * ut - 4.0 * hu * u2t * ut -
       4.0 * d2 * hu * ut - 2.0 * hd * dt * u2 + 
       (32.0 * gsq(3) - 0.8 * gsq(1)) * hddT + 2.4 * gsq(1) * heeT +
       1.6 * gsq(1) * hu * ut + (6.0 * gsq(2) + 1.2 * gsq(1)) * hd * dt -
       (32.0 * gsq(3) * mG(3) - 0.8 * gsq(1) * mG(1)) * ddT - 
       2.4 * gsq(1) * mG(1) * eeT - 
       (12.0 * gsq(2) * mG(2) + 1.6 * gsq(1) * mG(1)) * d2 -
       1.6 * gsq(1) * mG(1) * u2 + 64.0 / 9.0 * g4(3) * mG(3) - 
       16.0 * gsq(3) * gsq(2) * (mG(3) + mG(2)) - 
       16.0 / 9.0 * gsq(3) * gsq(1) * (mG(3) + mG(1)) - 30.0 * g4(2) * mG(2) - 
       2.0 * gsq(2) * gsq(1) * (mG(2) + mG(1)) - 574.0 / 45.0 * g4(1) * mG(1)) 
       * d1; 
    */      

	dhu2 = 
	  (-3.0 * (3.0 * u4 + ut * d2 * u1).trace() - 
	   d2 * (3.0 * ddT + eeT)
	   -15.0 * u2 * uuT - 6.0 * u4 - 2.0 * d4 -
	   4.0 * u2 * d2 + (16.0 * gsq(3) + 0.8 * gsq(1)) * uuT +
	   12.0 * gsq(2) * u2 + 0.4 * gsq(1) * d2 - 16.0 / 9.0 * g4(3) +
	   8.0 * gsq(3) * gsq(2) + 136.0 / 45.0 * gsq(3) * gsq(1) + 
	   7.5 * g4(2) + 
	   gsq(2) * gsq(1) + 2743.0 / 450.0 * g4(1)) * hu +
	  (-6.0 * (6.0 * hut * u2 * u1 + hut * d2 * u1 + 
		   hdt * u2 * d1).trace() -
	   18.0 * u2 * huuT - d2 * (6.0 * hddT + 2.0 * heeT) - 
	   12.0 * hu * ut * uuT - hd * dt * (6.0 * ddT + 2.0 * eeT) - 
	   6.0 * hu * u2t * ut - 8.0 * u2 * hu * ut - 4.0 * hd * d2t * dt -
	   4.0 * d2 * hd * dt - 2.0 * hu * ut * d2 - 4.0 * u2 * hd * dt +
	   (32.0 * gsq(3) + 1.6 * gsq(1)) * huuT + 
	   (6.0 * gsq(2) + 1.2 * gsq(1)) * hu * ut + 0.8 * gsq(1) * hd * dt -
	   (32.0 * gsq(3) * mG(3) + 1.6 * gsq(1) * mG(1)) * uuT -
	   (12.0 * gsq(2) * mG(2) + 0.8 * gsq(1) * mG(1)) * u2 -
	   0.8 * gsq(1) * mG(1) * d2 + 64.0 / 9.0 * g4(3) * mG(3) - 
	   16.0 * gsq(3) * gsq(2) * (mG(3) + mG(2)) - 
	   272.0 / 45.0 * gsq(3) * gsq(1) * (mG(1) + mG(3)) - 
	   30.0 * g4(2) * mG(2) - 2.0 * gsq(2) * gsq(1) * (mG(2) + mG(1)) - 
	   5486.0 / 225.0 * g4(1) * mG(1)) * u1; // checked

	dhd2 =  // checked + debugged 19.05.04 by M R Ramage (thanks Mike!)
	  (-3.0 * (3.0 * d4t + ut * d2 * u1 + e2t * e2t).trace() -
	   3.0 * u2 * uuT - 5.0 * d2 * (3.0 * ddT + eeT) - 6.0 * d4 -
	   2.0 * u4 - 4.0 * d2 * u2 + (16.0 * gsq(3) 
					    - 0.4 * gsq(1)) * ddT +
	   1.2 * gsq(1) * eeT + 0.8 * gsq(1) * u2 + 
	   (12.0 * gsq(2) + 1.2 * gsq(1)) * d2 - 16.0 / 9.0 * g4(3) + 
	   8.0 * gsq(3) * gsq(2) + 8.0 / 9.0 * gsq(3) * gsq(1) + 7.5 * g4(2) +
	   gsq(2) * gsq(1) + 287.0 / 90.0 * g4(1)) * hd +
	  (-6.0 * (6.0 * hdt * d2 * d1 + hut * d2 * u1 + hdt * u2 * d1 +
		   2.0 * het * e2 * e1).trace() - 6.0 * u2 * huuT -
	   6.0 * d2 * (3.0 * hddT + heeT) - 6.0 * hu * ut * uuT -
	   4.0 * hd * dt * (3.0 * ddT + eeT) - 6.0 * hd * d2t * dt -
	   8.0 * d2 * hd * dt - 4.0 * u2 * hu * ut - 4.0 * hu * u2t * ut -
	   4.0 * d2 * hu * ut - 2.0 * hd * dt * u2 + 
	   (32.0 * gsq(3) - 0.8 * gsq(1)) * hddT + 2.4 * gsq(1) * heeT +
	   1.6 * gsq(1) * hu * ut + (6.0 * gsq(2) + 1.2 * gsq(1)) * hd * dt -
	   (32.0 * gsq(3) * mG(3) - 0.8 * gsq(1) * mG(1)) * ddT -
	   2.4 * gsq(1) * mG(1) * eeT -
	   (12.0 * gsq(2) * mG(2) + 1.6 * gsq(1) * mG(1)) * d2 -
	   1.6 * gsq(1) * mG(1) * u2 + 64.0 / 9.0 * g4(3) * mG(3) -
	   16.0 * gsq(3) * gsq(2) * (mG(3) + mG(2)) - 
	   16.0 / 9.0 * gsq(3) * gsq(1) * (mG(3) + mG(1)) -
	   30.0 * g4(2) * mG(2) -
	   2.0 * gsq(2) * gsq(1) * (mG(2) + mG(1)) -
	   574.0 / 45.0 * g4(1) * mG(1)) 
	  * d1; 

	dhe2 = 
	  (-3.0 * (3.0 * d4t + ut * d2 * u1 + e2t * e2t).trace() -
	   5.0 * e2 * (3.0 * ddT + eeT) - 6.0 * e2 * e2 + 
	   (16.0 * gsq(3) - 0.4 * gsq(1)) * ddT + 1.2 * gsq(1) * eeT +
	   (12.0 * gsq(2) - 1.2 * gsq(1)) * e2 + 7.5 * g4(2) + 
	   1.8 * gsq(2) * gsq(1) + 13.5 * g4(1)) * he +
	  (-6.0 * (6.0 * hdt * d2 * d1 + hut * d2 * u1 + hdt * u2 * d1 +
		   2.0 * het * e2 * e1).trace() - 
	   4.0 * he * et * (3.0 * ddT + eeT) - 6.0 * e2 * (3.0 * hddT + heeT) -
	   6.0 * he * e2t * et - 8.0 * e2 * he * et +
	   (32.0 * gsq(3) - 0.8 * gsq(1)) * hddT + 2.4 * gsq(1) * heeT +
	   (6.0 * gsq(2) + 1.2 * gsq(1)) * he * et -
	   (32.0 * gsq(3) * mG(3) - 0.8 * gsq(1) * mG(1)) * ddT -
	   2.4 * gsq(1) * mG(1) * eeT - 12.0 * gsq(2) * mG(2) * e2 -
	   30.0 * g4(2) * mG(2) - 3.6 * gsq(2) * gsq(1) * (mG(1) + mG(2)) -
	   54.0 * g4(1) * mG(1)) * e1; // checked
      }	
    
    // add onto one-loop beta functions
    dmq2 *= oneO16Pif; dmq += dmq2;
    dml2 *= oneO16Pif; dml += dml2;
    dmu2 *= oneO16Pif; dmu += dmu2;
    dmd2 *= oneO16Pif; dmd += dmd2;
    dme2 *= oneO16Pif; dme += dme2;
    dhu2 *= oneO16Pif; dhu += dhu2;
    dhd2 *= oneO16Pif; dhd += dhd2;
    dhe2 *= oneO16Pif; dhe += dhe2;
    }

    // Default is to include these 2-loop corrections anyhow because they can 
    // be so important: gauginos + higgs 
    dmG2 = 2.0 * gsq * 
      (babBeta * gsqM + mG * (babBeta * gsq) +
       cuBeta * huuT - cuBeta * mG * uuT + 
       cdBeta * hddT - cdBeta * mG * ddT + 
       ceBeta * heeT - ceBeta * mG * eeT); // checked

    double dmH1sq2, dm3sq2, dmH2sq2;
    if (MIXING < 1) {
      /// The following are valid in the third-family approximation -- they are
      /// much faster, and should be a good approximation to the no-mixed case
      dm3sq2 = m3sq * 
	(-3.0 * (3.0 * ht4 + 3.0 * hb4 + 2.0 * hb2 * ht2 + htau4) +
	 (16.0 * gsq(3) + 0.8 * gsq(1)) * uuT + 
	 (16.0 * gsq(3) - 0.4 * gsq(1)) * ddT + 
	 1.2 * gsq(1) * eeT + 7.5 * g4(2) + 1.8 * gsq(1) * gsq(2) + 
	 207. / 50. * g4(1)) +
	displaySusyMu() * 
	(-12.0 * (3.0 * Ut * ht * ht2 + 3.0 * Ub * hb * hb2 + Ut * hb2 * ht +
		  Ub * ht2 * hb + Utau * htau2 * htau) +
	 (32.0 * gsq(3) + 1.6 * gsq(1)) * huuT + 
	 (32.0 * gsq(3) - 0.8 * gsq(1)) * hddT + 2.4 * gsq(1) * heeT -
	 (32.0 * gsq(3) * mG(3) + 1.6 * gsq(1) * mG(1)) * uuT -
	 (32.0 * gsq(3) * mG(3) - 0.8 * gsq(1) * mG(1)) * ddT -
	 2.4 * gsq(1) * mG(1) * eeT - 30.0 * g4(2) * mG(2) -
	 18.0 / 5.0 * gsq(1)
	 * gsq(2) * (mG(1) + mG(2)) - 414. / 25.0 * g4(1) * mG(1)); // checked
      
      dmH2sq2 = -6.0 * 
	(6.0 * (mH2sq + mq3 + mu3) * ht4 + 
	 hb2 * ht2 * (mH1sq + mH2sq + 2.0 * mq3 + mu3 + md3) + 
	 ht2 * 12.0 * Ut2 + ht2 * Ub2 + hb2 * Ut2 + 2.0 * hb * ht * Ut * Ub) +
	(32.0 * gsq(3) + 1.6 * gsq(1)) * 
	((mH2sq + mq3 + mu3) * ht2 + Ut2) +
	32.0 * gsq(3) * (2.0 * msq(3) * ht2 - 2.0 * mG(3) * ht * Ut) +
	1.6 * gsq(1) * (2.0 * msq(1) * ht2 - 2.0 * mG(1) * ht * Ut) +
	1.2 * gsq(1) * sP + 33.0 * g4(2) * msq(2) + 
	18.0 / 5.0 * gsq(2) * gsq(1) * (msq(2) + msq(1) + mG(1) * mG(2)) +
	621.0 / 25.0 * g4(1) * msq(1) + 3.0 * gsq(2) * sigma(2) +
	0.6 * gsq(1) * sigma(1); // checked
    
      dmH1sq2 = -6.0 *
	(6.0 * (mH1sq + mq3 + md3) * hb4 + 
	 (mH1sq + mH2sq + 2.0 * mq3 + mu3 + md3) * ht2 * hb2 + 
	 2.0 * (mH1sq + ml3 + me3) * htau4 +
	 12.0 * Ub2 * hb2 + hb2 * Ut2 + ht2 * Ub2 + 2.0 * Ut * ht * Ub * hb +
	 4.0 * htau2 * Utau2) +
	(32.0 * gsq(3) - 0.8 * gsq(1)) * ((mH1sq + mq3 + md3) * hb2 + Ub2) +
	32.0 * gsq(3) * (2.0 * msq(3) * ddT - 2.0 * mG(3) * hddT) -
	0.8 * gsq(1) * (2.0 * msq(1) * ddT - 2.0 * mG(1) * hddT) +
	2.4 * gsq(1) * (((mH1sq + ml3 + me3) * htau2 + Utau2) +
			2.0 * msq(1) * eeT - 2.0 * mG(1) * heeT) 
	- 1.2 * gsq(1) * sP +
	33.0 * g4(2) * msq(2) + 
	3.6 * gsq(2) * gsq(1) * (msq(2) + msq(1) + mG(1) * mG(2)) +
	621.0 / 25.0 * g4(1) * msq(1) + 3.0 * gsq(2) * sigma(2) + 
	0.6 * gsq(1) * sigma(1); // checked
    } else {
      /// In the mixed case, we need to use the slower full 3-family
      /// expressions
      dmH2sq2 = -6.0 * 
	(6.0 * (mH2sq + mq) * u4 + 6.0 * u1 * mu * u2t * ut +
	 (mH1sq + mH2sq + mq) * u2 * d2 + u1 * mu * ut * d2 +
	 u2 * mq * d2 + u2 * d1 * md * dt + 6.0 * hu2 * u2 +
	 6.0 * hu * u2t * hut + hd2 * u2 + d2 * hu2 + 
	 hd * dt * u1 * hut + d1 * hdt * hu * ut).trace() +
	(32.0 * gsq(3) + 1.6 * gsq(1)) * 
	((mH2sq + mq) * u2 + u1 * mu * ut + hu2).trace() +

	32.0 * gsq(3) * (2.0 * msq(3) * uuT - 2.0 * mG(3) * huuT) +
	1.6 * gsq(1) * (2.0 * msq(1) * uuT - 2.0 * mG(1) * huuT) +
	1.2 * gsq(1) * sP + 33.0 * g4(2) * msq(2) + 
	18.0 / 5.0 * gsq(2) * gsq(1) * (msq(2) + msq(1) + mG(1) * mG(2)) +
	621.0 / 25.0 * g4(1) * msq(1) + 3.0 * gsq(2) * sigma(2) +
	0.6 * gsq(1) * sigma(1); // checked
      
      dmH1sq2 = -6.0 *
	(6.0 * (mH1sq + mq) * d4 + 6.0 * dt * md * d2 * d1 +
	 (mH1sq + mH2sq + mq) * u2 * d2 + u1 * mu * ut * d2 +
	 u2 * mq * d2 + u2 * d1 * md * dt + 2.0 * (mH1sq + ml) * e2 * e2 +
	 2.0 * e1 * me * e2t * et + 6.0 * hd2 * d2 + 6.0 * hd * d2t * hdt + 
	 hu2 * d2 + u2 * hd2 + hu * ut * d1 * hdt + u1 * hut * hd * dt +
	 2.0 * he2 * e2 + 2.0 * he * e2t * het).trace() +
	(32.0 * gsq(3) - 0.8 * gsq(1)) * ((mH1sq + mq) * d2 + d1 * md * dt +
					  hd2).trace() +
	32.0 * gsq(3) * (2.0 * msq(3) * ddT - 2.0 * mG(3) * hddT) -
	0.8 * gsq(1) * (2.0 * msq(1) * ddT - 2.0 * mG(1) * hddT) +
	2.4 * gsq(1) * (((mH1sq + ml) * e2 + e1 * me * et + he2).trace() +
			2.0 * msq(1) * eeT - 2.0 * mG(1) * heeT) - 
	1.2 * gsq(1) * sP +
	33.0 * g4(2) * msq(2) + 
	3.6 * gsq(2) * gsq(1) * (msq(2) + msq(1) + mG(1) * mG(2)) +
	621.0 / 25.0 * g4(1) * msq(1) + 3.0 * gsq(2) * sigma(2) + 
	0.6 * gsq(1) *
	sigma(1); // checked+corrected 27/9/05

      dm3sq2 = m3sq * 
	(-3.0 * (3.0 * u4t + 3.0 * d4t + 2.0 * ut * d2
		 * u1 + e2t * e2t).trace() +
	 (16.0 * gsq(3) + 0.8 * gsq(1)) * uuT + 
	 (16.0 * gsq(3) - 0.4 * gsq(1)) * ddT + 
	 1.2 * gsq(1) * eeT + 7.5 * g4(2) + 1.8 * gsq(1) * gsq(2) + 
	 207. / 50. * g4(1)) +
	displaySusyMu() * 
	(-12.0 * (3.0 * hut * u2 * u1 + 3.0 * hdt * d2 * d1 + hut * d2 * u1 +
		  hdt * u2 * d1 + het * e2 * e1).trace() +
	 (32.0 * gsq(3) + 1.6 * gsq(1)) * huuT + 
	 (32.0 * gsq(3) - 0.8 * gsq(1)) * hddT + 2.4 * gsq(1) * heeT -
	 (32.0 * gsq(3) * mG(3) + 1.6 * gsq(1) * mG(1)) * uuT -
	 (32.0 * gsq(3) * mG(3) - 0.8 * gsq(1) * mG(1)) * ddT -
	 2.4 * gsq(1) * mG(1) * eeT - 30.0 * g4(2) * mG(2) 
	 - 18.0 / 5.0 * gsq(1)
	 * gsq(2) * (mG(1) + mG(2)) - 414. / 25.0 * g4(1) * mG(1)); // checked
    }    

    dmG = dmG + dmG2 * oneO16Pif;
    dm3sq = dm3sq + dm3sq2 * oneO16Pif;
    dmH1sq = dmH1sq + dmH1sq2 * oneO16Pif;
    dmH2sq = dmH2sq + dmH2sq2 * oneO16Pif;

    if (displayLoops() > 2) {
      
      const static double threelp = 2.53945672191370e-7; // 1/(16 pi^2)^3
      const static double kz = 7.21234141895757; // 6 Zeta(3)
      
      // (powers of) gaugino masses 
      double M1 = mG(1), M12 = sqr(M1);
      double M2 = mG(2), M22 = sqr(M2);
      double M3 = mG(3), M32 = sqr(M3);
      
      //powers of gauge couplings
      double a1 = gsq(1),      a2  = gsq(2),       a3  = gsq(3);   
      double a12 = a1*gsq(1),  a22 = a2*gsq(2),    a32 = a3*gsq(3);   
      double a13 = a12*gsq(1), a23 = a22*gsq(2),   a33 = a32*gsq(3);   
      
      double dmH2sq3, dmH1sq3, dm3sq3, dM13, dM23, dM33;
      double smu_ = displaySusyMu();
      
      // three-loop contribution 3rd family approximation
      if (MIXING < 0) {
	
	static const double O750= .00133333333333333333 ;
	static const double O27= .03703703703703703703 ;
	static const double O90= .01111111111111111111 ;
	static const double O45= .02222222222222222222 ;
	static const double O6= .16666666666666666666 ;
	static const double O6750= .00014814814814814814 ;
	static const double O375= .00266666666666666666 ;
	static const double O300= .00333333333333333333 ;
	static const double O75= .01333333333333333333 ;
	static const double O9= .11111111111111111111 ;
	static const double O1125= .00088888888888888888 ;
	static const double O150= .00666666666666666666 ;
	static const double O12= .08333333333333333333 ;
	static const double O30= .03333333333333333333 ;
	static const double O900= .00111111111111111111 ;
	static const double O180= .00555555555555555555 ;
	static const double O225= .00444444444444444444 ;
	static const double O450= .00222222222222222222 ;
	static const double O15= .06666666666666666666 ;
	static const double O3= .33333333333333333333 ;
	static const double O60= .01666666666666666666 ;
	
	double dmq3, dmt3, dmb3, dmtau3, dml3;
	double dht3, dhb3, dhl3; 
	
	double mq3 = displaySoftMassSquared(mQl, 3, 3);
	double mu3 = displaySoftMassSquared(mUr, 3, 3);
	double md3 = displaySoftMassSquared(mDr, 3, 3);
	double ml3 = displaySoftMassSquared(mLl, 3, 3);
	double me3 = displaySoftMassSquared(mEr, 3, 3);
	double ht = displayYukawaElement(YU, 3, 3), ht2 = sqr(ht), ht3 = ht2*ht, ht4 = sqr(ht2), ht5 = ht4*ht, ht6 = ht2*ht4;
	double htau = displayYukawaElement(YE, 3, 3), htau2 = sqr(htau), htau3 = htau2*htau,
	  htau4 = sqr(htau2), htau5=htau4*htau, htau6 = htau2*htau4;
	double hb = displayYukawaElement(YD, 3, 3), hb2 = sqr(hb), hb3 = hb2*hb, hb4 = sqr(hb2), hb5 = hb4*hb, hb6 = hb2*hb4;
	double Ut = displayTrilinear(UA, 3, 3), Ut2 = sqr(Ut);
	double Ub = displayTrilinear(DA, 3, 3), Ub2 = sqr(Ub);
	double Utau = displayTrilinear(EA, 3, 3), Utau2 = sqr(Utau);
	dmq3 =  
	  + htau4*Ub2 * (+ 10)
	  + ht4*Ut2 * (+ 270+ 18*kz)
	  + ht4*Ub2 * (+ 20)
	  + ht6 * (+ 90*mu3+ 6*mu3*kz+ 90*mq3+ 6*mq3*kz+ 90*mH2sq+ 6*mH2sq*kz)
	  + hb*htau3*Ub*Utau * (+ 40)
	  + hb*ht3*Ub*Ut * (+ 80)
	  + hb2*htau2*Utau2 * (+ 40)
	  + hb2*htau2*Ub2 * (- 32)
	  + hb2*htau4 * (+ 20*ml3+ 20*me3+ 10*md3+ 10*mq3+ 30*mH1sq)
	  + hb2*ht2*Ut2 * (+ 80)
	  + hb2*ht2*Ub2 * (+ 80)
	  + hb2*ht4 * (+ 40*mu3+ 20*md3+ 60*mq3+ 40*mH2sq+ 20*mH1sq)
	  + hb3*htau*Ub*Utau * (- 32)
	  + hb3*ht*Ub*Ut * (+ 80)
	  + hb4*Utau2 * (- 8)
	  + hb4*Ut2 * (+ 20)
	  + hb4*Ub2 * (+ 270+ 18*kz)
	  + hb4*htau2 * (- 8*ml3- 8*me3- 16*md3- 16*mq3- 24*mH1sq)
	  + hb4*ht2 * (+ 20*mu3+ 40*md3+ 60*mq3+ 20*mH2sq+ 40*mH1sq)
	  + hb6 * (+ 90*md3+ 6*md3*kz+ 90*mq3+ 6*mq3*kz+ 90*mH1sq+ 6*mH1sq*kz)
	  + a3*htau2*Ub2 * (+ 16)
	  + a3*ht2*Ut2 * (+ 320*O3+ 64*kz)
	  + a3*ht4 * (+ 160*O3*mu3+ 32*mu3*kz+ 160*O3*mq3+ 32*mq3*kz+ 160*O3*mH2sq+ 32*mH2sq*kz)
	  + a3*hb*htau*Ub*Utau * (+ 32)
	  + a3*hb2*Utau2 * (+ 16)
	  + a3*hb2*Ub2 * (+ 320*O3+ 64*kz)
	  + a3*hb2*htau2 * (+ 16*ml3+ 16*me3+ 16*md3+ 16*mq3+ 32*mH1sq)
	  + a3*hb4 * (+ 160*O3*md3+ 32*md3*kz+ 160*O3*mq3+ 32*mq3*kz+ 160*O3*mH1sq+ 32*mH1sq*kz)
	  + a3*M3*ht3*Ut * (- 320*O3- 64*kz)
	  + a3*M3*hb*htau2*Ub * (- 32)
	  + a3*M3*hb2*htau*Utau * (- 32)
	  + a3*M3*hb3*Ub * (- 320*O3- 64*kz)
	  + a3*M32*ht4 * (+ 160*O3+ 32*kz)
	  + a3*M32*hb2*htau2 * (+ 32)
	  + a3*M32*hb4 * (+ 160*O3+ 32*kz)
	  + a32*Ut2 * (- 176*O3- 272*O9*kz)
	  + a32*Ub2 * (- 176*O3- 272*O9*kz)
	  + a32*ht2 * (- 176*O3*mu3- 272*O9*mu3*kz- 176*O3*mq3- 272*O9*mq3*kz- 176*O3*mH2sq- 272*O9*mH2sq*kz)
	  + a32*hb2 * (- 176*O3*md3- 272*O9*md3*kz- 176*O3*mq3- 272*O9*mq3*kz- 176*O3*mH1sq- 272*O9*mH1sq*kz)
	  + a32*M3*ht*Ut * (+ 192+ 1088*O9*kz)
	  + a32*M3*hb*Ub * (+ 192+ 1088*O9*kz)
	  + a32*M32*ht2 * (- 288- 544*O3*kz)
	  + a32*M32*hb2 * (- 288- 544*O3*kz)
	  + a33 * (+ 320*O9*mu3+ 320*O9*md3+ 640*O9*mq3)
	  + a33*M32 * (+ 20512*O9+ 1280*kz)
	  + a2*htau2*Ub2 * (+ 12)
	  + a2*ht2*Ut2 * (+ 120+ 24*kz)
	  + a2*ht4 * (+ 60*mu3+ 12*mu3*kz+ 60*mq3+ 12*mq3*kz+ 60*mH2sq+ 12*mH2sq*kz)
	  + a2*hb*htau*Ub*Utau * (+ 24)
	  + a2*hb2*Utau2 * (+ 12)
	  + a2*hb2*Ub2 * (+ 120+ 24*kz)
	  + a2*hb2*htau2 * (+ 12*ml3+ 12*me3+ 12*md3+ 12*mq3+ 24*mH1sq)
	  + a2*hb4 * (+ 60*md3+ 12*md3*kz+ 60*mq3+ 12*mq3*kz+ 60*mH1sq+ 12*mH1sq*kz)
	  + a2*M2*ht3*Ut * (- 120- 24*kz)
	  + a2*M2*hb*htau2*Ub * (- 24)
	  + a2*M2*hb2*htau*Utau * (- 24)
	  + a2*M2*hb3*Ub * (- 120- 24*kz)
	  + a2*M22*ht4 * (+ 60+ 12*kz)
	  + a2*M22*hb2*htau2 * (+ 24)
	  + a2*M22*hb4 * (+ 60+ 12*kz)
	  + a2*a3*Ut2 * (- 8)
	  + a2*a3*Ub2 * (- 8)
	  + a2*a3*ht2 * (- 8*mu3- 8*mq3- 8*mH2sq)
	  + a2*a3*hb2 * (- 8*md3- 8*mq3- 8*mH1sq)
	  + a2*a3*M3*ht*Ut * (+ 16)
	  + a2*a3*M3*hb*Ub * (+ 16)
	  + a2*a3*M32*ht2 * (- 16)
	  + a2*a3*M32*hb2 * (- 16)
	  + a2*a3*M2*ht*Ut * (+ 16)
	  + a2*a3*M2*hb*Ub * (+ 16)
	  + a2*a3*M2*M3*ht2 * (- 16)
	  + a2*a3*M2*M3*hb2 * (- 16)
	  + a2*a3*M22*ht2 * (- 16)
	  + a2*a3*M22*hb2 * (- 16)
	  + a2*a32 * (- 16*mu3- 16*md3- 32*mq3)
	  + a2*a32*M32 * (+ 192- 144*kz)
	  + a2*a32*M2*M3 * (+ 64- 96*kz)
	  + a2*a32*M22 * (+ 80- 48*kz)
	  + a22*Utau2 * (- 18)
	  + a22*Ut2 * (- 76.5- 10.5*kz)
	  + a22*Ub2 * (- 76.5- 10.5*kz)
	  + a22*htau2 * (- 18*ml3- 18*me3- 18*mH1sq)
	  + a22*ht2 * (- 76.5*mu3- 10.5*mu3*kz- 76.5*mq3- 10.5*mq3*kz- 76.5*mH2sq- 10.5*mH2sq*kz)
	  + a22*hb2 * (- 76.5*md3- 10.5*md3*kz- 76.5*mq3- 10.5*mq3*kz- 76.5*mH1sq- 10.5*mH1sq*kz)
	  + a22*M2*htau*Utau * (+ 60)
	  + a22*M2*ht*Ut * (+ 270+ 42*kz)
	  + a22*M2*hb*Ub * (+ 270+ 42*kz)
	  + a22*M22*htau2 * (- 90)
	  + a22*M22*ht2 * (- 405- 63*kz)
	  + a22*M22*hb2 * (- 405- 63*kz)
	  + a22*a3 * (- 16*ml3- 48*mq3- 16*mH2sq- 16*mH1sq)
	  + a22*a3*M32 * (+ 272- 72*kz)
	  + a22*a3*M2*M3 * (+ 400- 144*kz)
	  + a22*a3*M22 * (+ 664- 216*kz)
	  + a23 * (- 3*ml3- 9*mq3- 3*mH2sq- 3*mH1sq)
	  + a23*M22 * (+ 2157+ 630*kz)
	  + a1*htau2*Ub2 * (- 3.2+ 0.8*kz)
	  + a1*ht2*Ut2 * (+ 136*O3- 14.4*kz)
	  + a1*ht4 * (+ 68*O3*mu3- 7.2*mu3*kz+ 68*O3*mq3- 7.2*mq3*kz+ 68*O3*mH2sq- 7.2*mH2sq*kz)
	  + a1*hb*htau*Ub*Utau * (- 6.4+ 1.6*kz)
	  + a1*hb2*Utau2 * (- 3.2+ 0.8*kz)
	  + a1*hb2*Ub2 * (+ 88*O3- 8*kz)
	  + a1*hb2*htau2 * (- 3.2*ml3+ 0.8*ml3*kz- 3.2*me3+ 0.8*me3*kz- 3.2*md3+ 0.8*md3*kz- 3.2*mq3+ 0.8*mq3*kz- 6.4*mH1sq+ 1.6*mH1sq*kz)
	  + a1*hb4 * (+ 44*O3*md3- 4*md3*kz+ 44*O3*mq3- 4*mq3*kz+ 44*O3*mH1sq- 4*mH1sq*kz)
	  + a1*M1*ht3*Ut * (- 136*O3+ 14.4*kz)
	  + a1*M1*hb*htau2*Ub * (+ 6.4- 1.6*kz)
	  + a1*M1*hb2*htau*Utau * (+ 6.4- 1.6*kz)
	  + a1*M1*hb3*Ub * (- 88*O3+ 8*kz)
	  + a1*M12*ht4 * (+ 68*O3- 7.2*kz)
	  + a1*M12*hb2*htau2 * (- 6.4+ 1.6*kz)
	  + a1*M12*hb4 * (+ 44*O3- 4*kz)
	  + a1*a3*Ut2 * (- 27.2+ 128*O45*kz)
	  + a1*a3*Ub2 * (- 152*O15+ 128*O45*kz)
	  + a1*a3*ht2 * (- 27.2*mu3+ 128*O45*mu3*kz- 27.2*mq3+ 128*O45*mq3*kz- 27.2*mH2sq+ 128*O45*mH2sq*kz)
	  + a1*a3*hb2 * (- 152*O15*md3+ 128*O45*md3*kz- 152*O15*mq3+ 128*O45*mq3*kz- 152*O15*mH1sq+ 128*O45*mH1sq*kz)
	  + a1*a3*M3*ht*Ut * (+ 54.4- 256*O45*kz)
	  + a1*a3*M3*hb*Ub * (+ 304*O15- 256*O45*kz)
	  + a1*a3*M32*ht2 * (- 54.4+ 256*O45*kz)
	  + a1*a3*M32*hb2 * (- 304*O15+ 256*O45*kz)
	  + a1*a3*M1*ht*Ut * (+ 54.4- 256*O45*kz)
	  + a1*a3*M1*hb*Ub * (+ 304*O15- 256*O45*kz)
	  + a1*a3*M1*M3*ht2 * (- 54.4+ 256*O45*kz)
	  + a1*a3*M1*M3*hb2 * (- 304*O15+ 256*O45*kz)
	  + a1*a3*M12*ht2 * (- 54.4+ 256*O45*kz)
	  + a1*a3*M12*hb2 * (- 304*O15+ 256*O45*kz)
	  + a1*a32 * (- 16*O45*mu3- 16*O45*md3- 32*O45*mq3)
	  + a1*a32*M32 * (+ 2464*O15- 35.2*kz)
	  + a1*a32*M1*M3 * (+ 4864*O45- 352*O15*kz)
	  + a1*a32*M12 * (+ 592*O9- 176*O15*kz)
	  + a1*a2*Ut2 * (- 11.8+ 3*kz)
	  + a1*a2*Ub2 * (- 8.2+ 0.6*kz)
	  + a1*a2*ht2 * (- 11.8*mu3+ 3*mu3*kz- 11.8*mq3+ 3*mq3*kz- 11.8*mH2sq+ 3*mH2sq*kz)
	  + a1*a2*hb2 * (- 8.2*md3+ 0.6*md3*kz- 8.2*mq3+ 0.6*mq3*kz- 8.2*mH1sq+ 0.6*mH1sq*kz)
	  + a1*a2*M2*ht*Ut * (+ 23.6- 6*kz)
	  + a1*a2*M2*hb*Ub * (+ 16.4- 1.2*kz)
	  + a1*a2*M22*ht2 * (- 23.6+ 6*kz)
	  + a1*a2*M22*hb2 * (- 16.4+ 1.2*kz)
	  + a1*a2*M1*ht*Ut * (+ 23.6- 6*kz)
	  + a1*a2*M1*hb*Ub * (+ 16.4- 1.2*kz)
	  + a1*a2*M1*M2*ht2 * (- 23.6+ 6*kz)
	  + a1*a2*M1*M2*hb2 * (- 16.4+ 1.2*kz)
	  + a1*a2*M12*ht2 * (- 23.6+ 6*kz)
	  + a1*a2*M12*hb2 * (- 16.4+ 1.2*kz)
	  + a1*a2*a3*M32 * (- 6.4)
	  + a1*a2*a3*M2*M3 * (- 6.4)
	  + a1*a2*a3*M22 * (- 6.4)
	  + a1*a2*a3*M1*M3 * (- 6.4)
	  + a1*a2*a3*M1*M2 * (- 6.4)
	  + a1*a2*a3*M12 * (- 6.4)
	  + a1*a22 * (- 0.2*ml3- 0.6*mq3- 0.2*mH2sq- 0.2*mH1sq)
	  + a1*a22*M22 * (+ 75.8- 16.2*kz)
	  + a1*a22*M1*M2 * (+ 50- 10.8*kz)
	  + a1*a22*M12 * (+ 30.4- 5.4*kz)
	  + a12*Utau2 * (- 0.72)
	  + a12*Ut2 * (- 3923*O150+ 143*O450*kz)
	  + a12*Ub2 * (- 13.22+ 7*O90*kz)
	  + a12*htau2 * (- 0.72*ml3- 0.72*me3- 0.72*mH1sq)
	  + a12*ht2 * (- 0.96*ml3- 1.92*me3- 4307*O150*mu3+ 143*O450*mu3*kz- 0.64*md3- 3971*O150*mq3+ 143*O450*mq3*kz- 4067*O150*mH2sq+ 143*O450*mH2sq*kz- 0.96*mH1sq)
	  + a12*hb2 * (- 0.48*ml3- 0.96*me3- 1.28*mu3- 13.54*md3+ 7*O90*md3*kz- 13.38*mq3+ 7*O90*mq3*kz- 0.48*mH2sq- 13.7*mH1sq+ 7*O90*mH1sq*kz)
	  + a12*M1*htau*Utau * (+ 2.4)
	  + a12*M1*ht*Ut * (+ 103.92- 286*O225*kz)
	  + a12*M1*hb*Ub * (+ 3938*O75- 14*O45*kz)
	  + a12*M12*htau2 * (- 3.6)
	  + a12*M12*ht2 * (- 155.88+ 143*O75*kz)
	  + a12*M12*hb2 * (- 78.76+ 7*O15*kz)
	  + a12*a3 * (- 16*O75*ml3- 32*O75*me3- 128*O225*mu3- 32*O225*md3- 16*O225*mq3- 16*O75*mH2sq- 16*O75*mH1sq)
	  + a12*a3*M32 * (+ 208*O45- 88*O75*kz)
	  + a12*a3*M1*M3 * (+ 1552*O225- 176*O75*kz)
	  + a12*a3*M12 * (+ 776*O75- 3.52*kz)
	  + a12*a2 * (- 0.12*ml3- 0.24*me3- 0.32*mu3- 0.08*md3- 0.04*mq3- 0.12*mH2sq- 0.12*mH1sq)
	  + a12*a2*M22 * (+ 0.8- 0.36*kz)
	  + a12*a2*M1*M2 * (+ 0.88- 0.72*kz)
	  + a12*a2*M12 * (+ 1.32- 1.08*kz)
	  + a13 * (- 199*O375*ml3- 398*O375*me3- 1592*O1125*mu3- 398*O1125*md3- 199*O1125*mq3- 199*O375*mH2sq- 199*O375*mH1sq)
	  + a13*M12 * (+ 57511*O1125- 3.184*kz);
	
	dmt3 = 
	  + ht2*htau2*Ub2 * (+ 8)
	  + ht4*Ut2 * (+ 540+ 36*kz)
	  + ht4*Ub2 * (+ 4)
	  + ht6 * (+ 180*mu3+ 12*mu3*kz+ 180*mq3+ 12*mq3*kz+ 180*mH2sq+ 12*mH2sq*kz)
	  + hb*ht*htau2*Ub*Ut * (+ 16)
	  + hb*ht2*htau*Ub*Utau * (+ 16)
	  + hb*ht3*Ub*Ut * (+ 16)
	  + hb2*htau2*Ut2 * (+ 8)
	  + hb2*ht*htau*Ut*Utau * (+ 16)
	  + hb2*ht2*Utau2 * (+ 8)
	  + hb2*ht2*Ut2 * (+ 16)
	  + hb2*ht2*Ub2 * (+ 144)
	  + hb2*ht2*htau2 * (+ 8*ml3+ 8*me3+ 8*mu3+ 8*md3+ 16*mq3+ 8*mH2sq+ 16*mH1sq)
	  + hb2*ht4 * (+ 8*mu3+ 4*md3+ 12*mq3+ 8*mH2sq+ 4*mH1sq)
	  + hb3*ht*Ub*Ut * (+ 144)
	  + hb4*Ut2 * (+ 36)
	  + hb4*ht2 * (+ 36*mu3+ 72*md3+ 108*mq3+ 36*mH2sq+ 72*mH1sq)
	  + a3*ht2*Ut2 * (+ 128*O3+ 128*kz)
	  + a3*ht2*Ub2 * (+ 128*O3)
	  + a3*ht4 * (+ 64*O3*mu3+ 64*mu3*kz+ 64*O3*mq3+ 64*mq3*kz+ 64*O3*mH2sq+ 64*mH2sq*kz)
	  + a3*hb*ht*Ub*Ut * (+ 256*O3)
	  + a3*hb2*Ut2 * (+ 128*O3)
	  + a3*hb2*ht2 * (+ 128*O3*mu3+ 128*O3*md3+ 256*O3*mq3+ 128*O3*mH2sq+ 128*O3*mH1sq)
	  + a3*M3*ht3*Ut * (- 128*O3- 128*kz)
	  + a3*M3*hb*ht2*Ub * (- 256*O3)
	  + a3*M3*hb2*ht*Ut * (- 256*O3)
	  + a3*M32*ht4 * (+ 64*O3+ 64*kz)
	  + a3*M32*hb2*ht2 * (+ 256*O3)
	  + a32*Ut2 * (- 160*O3- 544*O9*kz)
	  + a32*Ub2 * (- 64)
	  + a32*ht2 * (- 160*O3*mu3- 544*O9*mu3*kz- 160*O3*mq3- 544*O9*mq3*kz- 160*O3*mH2sq- 544*O9*mH2sq*kz)
	  + a32*hb2 * (- 64*md3- 64*mq3- 64*mH1sq)
	  + a32*M3*ht*Ut * (+ 512*O3+ 2176*O9*kz)
	  + a32*M3*hb*Ub * (+ 640*O3)
	  + a32*M32*ht2 * (- 256- 1088*O3*kz)
	  + a32*M32*hb2 * (- 320)
	  + a33 * (+ 320*O9*mu3+ 320*O9*md3+ 640*O9*mq3)
	  + a33*M32 * (+ 20512*O9+ 1280*kz)
	  + a2*ht2*Ut2 * (+ 288- 96*kz)
	  + a2*ht2*Ub2 * (+ 18- 6*kz)
	  + a2*ht4 * (+ 144*mu3- 48*mu3*kz+ 144*mq3- 48*mq3*kz+ 144*mH2sq- 48*mH2sq*kz)
	  + a2*hb*ht*Ub*Ut * (+ 36- 12*kz)
	  + a2*hb2*Ut2 * (+ 18- 6*kz)
	  + a2*hb2*ht2 * (+ 18*mu3- 6*mu3*kz+ 18*md3- 6*md3*kz+ 36*mq3- 12*mq3*kz+ 18*mH2sq- 6*mH2sq*kz+ 18*mH1sq- 6*mH1sq*kz)
	  + a2*M2*ht3*Ut * (- 288+ 96*kz)
	  + a2*M2*hb*ht2*Ub * (- 36+ 12*kz)
	  + a2*M2*hb2*ht*Ut * (- 36+ 12*kz)
	  + a2*M22*ht4 * (+ 144- 48*kz)
	  + a2*M22*hb2*ht2 * (+ 36- 12*kz)
	  + a2*a3*Ut2 * (- 176+ 32*kz)
	  + a2*a3*ht2 * (- 176*mu3+ 32*mu3*kz- 176*mq3+ 32*mq3*kz- 176*mH2sq+ 32*mH2sq*kz)
	  + a2*a3*M3*ht*Ut * (+ 352- 64*kz)
	  + a2*a3*M32*ht2 * (- 352+ 64*kz)
	  + a2*a3*M2*ht*Ut * (+ 352- 64*kz)
	  + a2*a3*M2*M3*ht2 * (- 352+ 64*kz)
	  + a2*a3*M22*ht2 * (- 352+ 64*kz)
	  + a2*a32*M32 * (+ 720- 144*kz)
	  + a2*a32*M2*M3 * (+ 480- 96*kz)
	  + a2*a32*M22 * (+ 288- 48*kz)
	  + a22*Ut2 * (- 87- 3*kz)
	  + a22*ht2 * (- 12*ml3- 87*mu3- 3*mu3*kz- 123*mq3- 3*mq3*kz- 99*mH2sq- 3*mH2sq*kz- 12*mH1sq)
	  + a22*M2*ht*Ut * (+ 348+ 12*kz)
	  + a22*M22*ht2 * (- 474- 18*kz)
	  + a1*ht2*Ut2 * (+ 160*O3+ 19.2*kz)
	  + a1*ht2*Ub2 * (+ 38*O15+ 1.2*kz)
	  + a1*ht4 * (+ 80*O3*mu3+ 9.6*mu3*kz+ 80*O3*mq3+ 9.6*mq3*kz+ 80*O3*mH2sq+ 9.6*mH2sq*kz)
	  + a1*hb*ht*Ub*Ut * (+ 76*O15+ 2.4*kz)
	  + a1*hb2*Ut2 * (+ 38*O15+ 1.2*kz)
	  + a1*hb2*ht2 * (+ 38*O15*mu3+ 1.2*mu3*kz+ 38*O15*md3+ 1.2*md3*kz+ 76*O15*mq3+ 2.4*mq3*kz+ 38*O15*mH2sq+ 1.2*mH2sq*kz+ 38*O15*mH1sq+ 1.2*mH1sq*kz)
	  + a1*M1*ht3*Ut * (- 160*O3- 19.2*kz)
	  + a1*M1*hb*ht2*Ub * (- 76*O15- 2.4*kz)
	  + a1*M1*hb2*ht*Ut * (- 76*O15- 2.4*kz)
	  + a1*M12*ht4 * (+ 80*O3+ 9.6*kz)
	  + a1*M12*hb2*ht2 * (+ 76*O15+ 2.4*kz)
	  + a1*a3*Ut2 * (- 16*O15- 224*O45*kz)
	  + a1*a3*ht2 * (- 16*O15*mu3- 224*O45*mu3*kz- 16*O15*mq3- 224*O45*mq3*kz- 16*O15*mH2sq- 224*O45*mH2sq*kz)
	  + a1*a3*M3*ht*Ut * (+ 32*O15+ 448*O45*kz)
	  + a1*a3*M32*ht2 * (- 32*O15- 448*O45*kz)
	  + a1*a3*M1*ht*Ut * (+ 32*O15+ 448*O45*kz)
	  + a1*a3*M1*M3*ht2 * (- 32*O15- 448*O45*kz)
	  + a1*a3*M12*ht2 * (- 32*O15- 448*O45*kz)
	  + a1*a32 * (- 256*O45*mu3- 256*O45*md3- 512*O45*mq3)
	  + a1*a32*M32 * (- 176*O15- 35.2*kz)
	  + a1*a32*M1*M3 * (- 1376*O45- 352*O15*kz)
	  + a1*a32*M12 * (- 32*O9- 176*O15*kz)
	  + a1*a2*Ut2 * (- 26.8+ 5.2*kz)
	  + a1*a2*ht2 * (- 26.8*mu3+ 5.2*mu3*kz- 26.8*mq3+ 5.2*mq3*kz- 26.8*mH2sq+ 5.2*mH2sq*kz)
	  + a1*a2*M2*ht*Ut * (+ 53.6- 10.4*kz)
	  + a1*a2*M22*ht2 * (- 53.6+ 10.4*kz)
	  + a1*a2*M1*ht*Ut * (+ 53.6- 10.4*kz)
	  + a1*a2*M1*M2*ht2 * (- 53.6+ 10.4*kz)
	  + a1*a2*M12*ht2 * (- 53.6+ 10.4*kz)
	  + a12*Utau2 * (- 11.52)
	  + a12*Ut2 * (- 48.6- 247*O225*kz)
	  + a12*Ub2 * (- 8.96)
	  + a12*htau2 * (- 11.52*ml3- 11.52*me3- 11.52*mH1sq)
	  + a12*ht2 * (+ 0.48*ml3+ 0.96*me3- 47.32*mu3- 247*O225*mu3*kz+ 0.32*md3- 48.44*mq3- 247*O225*mq3*kz- 48.12*mH2sq- 247*O225*mH2sq*kz+ 0.48*mH1sq)
	  + a12*hb2 * (- 8.96*md3- 8.96*mq3- 8.96*mH1sq)
	  + a12*M1*htau*Utau * (+ 38.4)
	  + a12*M1*ht*Ut * (+ 13748*O75+ 988*O225*kz)
	  + a12*M1*hb*Ub * (+ 448*O15)
	  + a12*M12*htau2 * (- 57.6)
	  + a12*M12*ht2 * (- 274.96- 494*O75*kz)
	  + a12*M12*hb2 * (- 44.8)
	  + a12*a3 * (- 256*O75*ml3- 512*O75*me3- 2048*O225*mu3- 512*O225*md3- 256*O225*mq3- 256*O75*mH2sq- 256*O75*mH1sq)
	  + a12*a3*M32 * (+ 512*O9- 1408*O75*kz)
	  + a12*a3*M1*M3 * (+ 17152*O225- 2816*O75*kz)
	  + a12*a3*M12 * (+ 8576*O75- 56.32*kz)
	  + a12*a2*M22 * (+ 34.56- 5.76*kz)
	  + a12*a2*M1*M2 * (+ 57.6- 11.52*kz)
	  + a12*a2*M12 * (+ 86.4- 17.28*kz)
	  + a13 * (- 3424*O375*ml3- 6848*O375*me3- 27392*O1125*mu3- 6848*O1125*md3- 3424*O1125*mq3- 3424*O375*mH2sq- 3424*O375*mH1sq)
	  + a13*M12 * (+ 864496*O1125- 50.944*kz);
	
	dmb3 = 
	  + htau4*Ub2 * (+ 20)
	  + ht2*htau2*Ub2 * (- 4)
	  + ht4*Ub2 * (+ 36)
	  + hb*htau3*Ub*Utau * (+ 80)
	  + hb*ht*htau2*Ub*Ut * (- 8)
	  + hb*ht2*htau*Ub*Utau * (- 8)
	  + hb*ht3*Ub*Ut * (+ 144)
	  + hb2*htau2*Utau2 * (+ 80)
	  + hb2*htau2*Ut2 * (- 4)
       + hb2*htau2*Ub2 * (- 80)
	  + hb2*htau4 * (+ 40*ml3+ 40*me3+ 20*md3+ 20*mq3+ 60*mH1sq)
	  + hb2*ht*htau*Ut*Utau * (- 8)
	  + hb2*ht2*Utau2 * (- 4)
	  + hb2*ht2*Ut2 * (+ 144)
	  + hb2*ht2*Ub2 * (+ 16)
	  + hb2*ht2*htau2 * (- 4*ml3- 4*me3- 4*mu3- 4*md3- 8*mq3- 4*mH2sq- 8*mH1sq)
	  + hb2*ht4 * (+ 72*mu3+ 36*md3+ 108*mq3+ 72*mH2sq+ 36*mH1sq)
	  + hb3*htau*Ub*Utau * (- 80)
	  + hb3*ht*Ub*Ut * (+ 16)
	  + hb4*Utau2 * (- 20)
	  + hb4*Ut2 * (+ 4)
	  + hb4*Ub2 * (+ 540+ 36*kz)
	  + hb4*htau2 * (- 20*ml3- 20*me3- 40*md3- 40*mq3- 60*mH1sq)
	  + hb4*ht2 * (+ 4*mu3+ 8*md3+ 12*mq3+ 4*mH2sq+ 8*mH1sq)
	  + hb6 * (+ 180*md3+ 12*md3*kz+ 180*mq3+ 12*mq3*kz+ 180*mH1sq+ 12*mH1sq*kz)
	  + a3*htau2*Ub2 * (+ 32)
	  + a3*ht2*Ub2 * (+ 128*O3)
	  + a3*hb*htau*Ub*Utau * (+ 64)
	  + a3*hb*ht*Ub*Ut * (+ 256*O3)
	  + a3*hb2*Utau2 * (+ 32)
	  + a3*hb2*Ut2 * (+ 128*O3)
	  + a3*hb2*Ub2 * (+ 128*O3+ 128*kz)
	  + a3*hb2*htau2 * (+ 32*ml3+ 32*me3+ 32*md3+ 32*mq3+ 64*mH1sq)
	  + a3*hb2*ht2 * (+ 128*O3*mu3+ 128*O3*md3+ 256*O3*mq3+ 128*O3*mH2sq+ 128*O3*mH1sq)
	  + a3*hb4 * (+ 64*O3*md3+ 64*md3*kz+ 64*O3*mq3+ 64*mq3*kz+ 64*O3*mH1sq+ 64*mH1sq*kz)
	  + a3*M3*hb*htau2*Ub * (- 64)
	  + a3*M3*hb*ht2*Ub * (- 256*O3)
	  + a3*M3*hb2*htau*Utau * (- 64)
	  + a3*M3*hb2*ht*Ut * (- 256*O3)
	  + a3*M3*hb3*Ub * (- 128*O3- 128*kz)
	  + a3*M32*hb2*htau2 * (+ 64)
	  + a3*M32*hb2*ht2 * (+ 256*O3)
	  + a3*M32*hb4 * (+ 64*O3+ 64*kz)
	  + a32*Ut2 * (- 64)
	  + a32*Ub2 * (- 160*O3- 544*O9*kz)
	  + a32*ht2 * (- 64*mu3- 64*mq3- 64*mH2sq)
	  + a32*hb2 * (- 160*O3*md3- 544*O9*md3*kz- 160*O3*mq3- 544*O9*mq3*kz- 160*O3*mH1sq- 544*O9*mH1sq*kz)
	  + a32*M3*ht*Ut * (+ 640*O3)
	  + a32*M3*hb*Ub * (+ 512*O3+ 2176*O9*kz)
	  + a32*M32*ht2 * (- 320)
	  + a32*M32*hb2 * (- 256- 1088*O3*kz)
	  + a33 * (+ 320*O9*mu3+ 320*O9*md3+ 640*O9*mq3)
	  + a33*M32 * (+ 20512*O9+ 1280*kz)
	  + a2*htau2*Ub2 * (+ 18- 6*kz)
	  + a2*ht2*Ub2 * (+ 18- 6*kz)
	  + a2*hb*htau*Ub*Utau * (+ 36- 12*kz)
	  + a2*hb*ht*Ub*Ut * (+ 36- 12*kz)
	  + a2*hb2*Utau2 * (+ 18- 6*kz)
	  + a2*hb2*Ut2 * (+ 18- 6*kz)
	  + a2*hb2*Ub2 * (+ 288- 96*kz)
	  + a2*hb2*htau2 * (+ 18*ml3- 6*ml3*kz+ 18*me3- 6*me3*kz+ 18*md3- 6*md3*kz+ 18*mq3- 6*mq3*kz+ 36*mH1sq- 12*mH1sq*kz)
	  + a2*hb2*ht2 * (+ 18*mu3- 6*mu3*kz+ 18*md3- 6*md3*kz+ 36*mq3- 12*mq3*kz+ 18*mH2sq- 6*mH2sq*kz+ 18*mH1sq- 6*mH1sq*kz)
	  + a2*hb4 * (+ 144*md3- 48*md3*kz+ 144*mq3- 48*mq3*kz+ 144*mH1sq- 48*mH1sq*kz)
	  + a2*M2*hb*htau2*Ub * (- 36+ 12*kz)
	  + a2*M2*hb*ht2*Ub * (- 36+ 12*kz)
	  + a2*M2*hb2*htau*Utau * (- 36+ 12*kz)
	  + a2*M2*hb2*ht*Ut * (- 36+ 12*kz)
	  + a2*M2*hb3*Ub * (- 288+ 96*kz)
	  + a2*M22*hb2*htau2 * (+ 36- 12*kz)
	  + a2*M22*hb2*ht2 * (+ 36- 12*kz)
	  + a2*M22*hb4 * (+ 144- 48*kz)
	  + a2*a3*Ub2 * (- 176+ 32*kz)
	  + a2*a3*hb2 * (- 176*md3+ 32*md3*kz- 176*mq3+ 32*mq3*kz- 176*mH1sq+ 32*mH1sq*kz)
	  + a2*a3*M3*hb*Ub * (+ 352- 64*kz)
	  + a2*a3*M32*hb2 * (- 352+ 64*kz)
	  + a2*a3*M2*hb*Ub * (+ 352- 64*kz)
	  + a2*a3*M2*M3*hb2 * (- 352+ 64*kz)
	  + a2*a3*M22*hb2 * (- 352+ 64*kz)
	  + a2*a32*M32 * (+ 720- 144*kz)
	  + a2*a32*M2*M3 * (+ 480- 96*kz)
	  + a2*a32*M22 * (+ 288- 48*kz)
	  + a22*Ub2 * (- 87- 3*kz)
	  + a22*hb2 * (- 12*ml3- 87*md3- 3*md3*kz- 123*mq3- 3*mq3*kz- 12*mH2sq- 99*mH1sq- 3*mH1sq*kz)
	  + a22*M2*hb*Ub * (+ 348+ 12*kz)
	  + a22*M22*hb2 * (- 474- 18*kz)
	  + a1*htau2*Ub2 * (- 6+ 2*kz)
	  + a1*ht2*Ub2 * (- 58*O15+ 1.2*kz)
	  + a1*hb*htau*Ub*Utau * (- 12+ 4*kz)
	  + a1*hb*ht*Ub*Ut * (- 116*O15+ 2.4*kz)
	  + a1*hb2*Utau2 * (- 6+ 2*kz)
	  + a1*hb2*Ut2 * (- 58*O15+ 1.2*kz)
	  + a1*hb2*Ub2 * (+ 160*O3- 6.4*kz)
	  + a1*hb2*htau2 * (- 6*ml3+ 2*ml3*kz- 6*me3+ 2*me3*kz- 6*md3+ 2*md3*kz- 6*mq3+ 2*mq3*kz- 12*mH1sq+ 4*mH1sq*kz)
	  + a1*hb2*ht2 * (- 58*O15*mu3+ 1.2*mu3*kz- 58*O15*md3+ 1.2*md3*kz- 116*O15*mq3+ 2.4*mq3*kz- 58*O15*mH2sq+ 1.2*mH2sq*kz- 58*O15*mH1sq+ 1.2*mH1sq*kz)
	  + a1*hb4 * (+ 80*O3*md3- 3.2*md3*kz+ 80*O3*mq3- 3.2*mq3*kz+ 80*O3*mH1sq- 3.2*mH1sq*kz)
	  + a1*M1*hb*htau2*Ub * (+ 12- 4*kz)
	  + a1*M1*hb*ht2*Ub * (+ 116*O15- 2.4*kz)
	  + a1*M1*hb2*htau*Utau * (+ 12- 4*kz)
	  + a1*M1*hb2*ht*Ut * (+ 116*O15- 2.4*kz)
	  + a1*M1*hb3*Ub * (- 160*O3+ 6.4*kz)
	  + a1*M12*hb2*htau2 * (- 12+ 4*kz)
	  + a1*M12*hb2*ht2 * (- 116*O15+ 2.4*kz)
	  + a1*M12*hb4 * (+ 80*O3- 3.2*kz)
	  + a1*a3*Ub2 * (- 9.6+ 32*O9*kz)
	  + a1*a3*hb2 * (- 9.6*md3+ 32*O9*md3*kz- 9.6*mq3+ 32*O9*mq3*kz- 9.6*mH1sq+ 32*O9*mH1sq*kz)
	  + a1*a3*M3*hb*Ub * (+ 19.2- 64*O9*kz)
	  + a1*a3*M32*hb2 * (- 19.2+ 64*O9*kz)
	  + a1*a3*M1*hb*Ub * (+ 19.2- 64*O9*kz)
	  + a1*a3*M1*M3*hb2 * (- 19.2+ 64*O9*kz)
	  + a1*a3*M12*hb2 * (- 19.2+ 64*O9*kz)
	  + a1*a32 * (- 64*O45*mu3- 64*O45*md3- 128*O45*mq3)
	  + a1*a32*M32 * (+ 1936*O15- 35.2*kz)
	  + a1*a32*M1*M3 * (+ 3616*O45- 352*O15*kz)
	  + a1*a32*M12 * (+ 2336*O45- 176*O15*kz)
	  + a1*a2*Ub2 * (- 17.2+ 2.8*kz)
	  + a1*a2*hb2 * (- 17.2*md3+ 2.8*md3*kz- 17.2*mq3+ 2.8*mq3*kz- 17.2*mH1sq+ 2.8*mH1sq*kz)
	  + a1*a2*M2*hb*Ub * (+ 34.4- 5.6*kz)
	  + a1*a2*M22*hb2 * (- 34.4+ 5.6*kz)
	  + a1*a2*M1*hb*Ub * (+ 34.4- 5.6*kz)
	  + a1*a2*M1*M2*hb2 * (- 34.4+ 5.6*kz)
	  + a1*a2*M12*hb2 * (- 34.4+ 5.6*kz)
	  + a12*Utau2 * (- 2.88)
	  + a12*Ut2 * (- 4.16)
	  + a12*Ub2 * (- 1853*O75- 7*O225*kz)
	  + a12*htau2 * (- 2.88*ml3- 2.88*me3- 2.88*mH1sq)
	  + a12*ht2 * (- 4.16*mu3- 4.16*mq3- 4.16*mH2sq)
	  + a12*hb2 * (- 0.48*ml3- 0.96*me3- 1.28*mu3- 1877*O75*md3- 7*O225*md3*kz- 373*O15*mq3- 7*O225*mq3*kz- 0.48*mH2sq- 1889*O75*mH1sq- 7*O225*mH1sq*kz)
	  + a12*M1*htau*Utau * (+ 9.6)
	  + a12*M1*ht*Ut * (+ 208*O15)
	  + a12*M1*hb*Ub * (+ 292*O3+ 28*O225*kz)
	  + a12*M12*htau2 * (- 14.4)
	  + a12*M12*ht2 * (- 20.8)
	  + a12*M12*hb2 * (- 146- 14*O75*kz)
	  + a12*a3 * (- 64*O75*ml3- 128*O75*me3- 512*O225*mu3- 128*O225*md3- 64*O225*mq3- 64*O75*mH2sq- 64*O75*mH1sq)
	  + a12*a3*M32 * (+ 3968*O225- 352*O75*kz)
	  + a12*a3*M1*M3 * (+ 5824*O225- 704*O75*kz)
	  + a12*a3*M12 * (+ 2912*O75- 14.08*kz)
	  + a12*a2*M22 * (+ 8.64- 1.44*kz)
	  + a12*a2*M1*M2 * (+ 14.4- 2.88*kz)
	  + a12*a2*M12 * (+ 21.6- 4.32*kz)
	  + a13 * (- 808*O375*ml3- 1616*O375*me3- 6464*O1125*mu3- 1616*O1125*md3- 808*O1125*mq3- 808*O375*mH2sq- 808*O375*mH1sq)
	  + a13*M12 * (+ 227548*O1125- 12.736*kz);
	
	dml3 = 
	  + htau4*Utau2 * (+ 126+ 18*kz)
	  + htau6 * (+ 42*ml3+ 6*ml3*kz+ 42*me3+ 6*me3*kz+ 42*mH1sq+ 6*mH1sq*kz)
	  + ht2*htau2*Ub2 * (+ 12)
	  + hb*ht*htau2*Ub*Ut * (+ 24)
	  + hb*ht2*htau*Ub*Utau * (+ 24)
	  + hb2*htau2*Ut2 * (+ 12)
	  + hb2*htau2*Ub2 * (+ 72)
	  + hb2*ht*htau*Ut*Utau * (+ 24)
	  + hb2*ht2*Utau2 * (+ 12)
	  + hb2*ht2*htau2 * (+ 12*ml3+ 12*me3+ 12*mu3+ 12*md3+ 24*mq3+ 12*mH2sq+ 24*mH1sq)
	  + hb3*htau*Ub*Utau * (+ 72)
	  + hb4*Utau2 * (+ 18)
	  + hb4*htau2 * (+ 18*ml3+ 18*me3+ 36*md3+ 36*mq3+ 54*mH1sq)
	  + a3*htau2*Ub2 * (- 64+ 16*kz)
	  + a3*hb*htau*Ub*Utau * (- 128+ 32*kz)
	  + a3*hb2*Utau2 * (- 64+ 16*kz)
	  + a3*hb2*htau2 * (- 64*ml3+ 16*ml3*kz- 64*me3+ 16*me3*kz- 64*md3+ 16*md3*kz- 64*mq3+ 16*mq3*kz- 128*mH1sq+ 32*mH1sq*kz)
	  + a3*M3*hb*htau2*Ub * (+ 128- 32*kz)
	  + a3*M3*hb2*htau*Utau * (+ 128- 32*kz)
	  + a3*M32*hb2*htau2 * (- 128+ 32*kz)
	  + a2*htau2*Utau2 * (+ 24+ 24*kz)
	  + a2*htau2*Ub2 * (+ 36)
	  + a2*htau4 * (+ 12*ml3+ 12*ml3*kz+ 12*me3+ 12*me3*kz+ 12*mH1sq+ 12*mH1sq*kz)
	  + a2*hb*htau*Ub*Utau * (+ 72)
	  + a2*hb2*Utau2 * (+ 36)
	  + a2*hb2*htau2 * (+ 36*ml3+ 36*me3+ 36*md3+ 36*mq3+ 72*mH1sq)
	  + a2*M2*htau3*Utau * (- 24- 24*kz)
	  + a2*M2*hb*htau2*Ub * (- 72)
	  + a2*M2*hb2*htau*Utau * (- 72)
	  + a2*M22*htau4 * (+ 12+ 12*kz)
	  + a2*M22*hb2*htau2 * (+ 72)
	  + a22*Utau2 * (- 40.5- 10.5*kz)
	  + a22*Ut2 * (- 54)
	  + a22*Ub2 * (- 54)
	  + a22*htau2 * (- 40.5*ml3- 10.5*ml3*kz- 40.5*me3- 10.5*me3*kz- 40.5*mH1sq- 10.5*mH1sq*kz)
	  + a22*ht2 * (- 54*mu3- 54*mq3- 54*mH2sq)
	  + a22*hb2 * (- 54*md3- 54*mq3- 54*mH1sq)
	  + a22*M2*htau*Utau * (+ 150+ 42*kz)
	  + a22*M2*ht*Ut * (+ 180)
	  + a22*M2*hb*Ub * (+ 180)
	  + a22*M22*htau2 * (- 225- 63*kz)
	  + a22*M22*ht2 * (- 270)
	  + a22*M22*hb2 * (- 270)
	  + a22*a3*M32 * (+ 432- 72*kz)
	  + a22*a3*M2*M3 * (+ 720- 144*kz)
	  + a22*a3*M22 * (+ 1080- 216*kz)
	  + a23 * (- 3*ml3- 9*mq3- 3*mH2sq- 3*mH1sq)
	  + a23*M22 * (+ 2157+ 630*kz)
	  + a1*htau2*Utau2 * (+ 72- 14.4*kz)
	  + a1*htau2*Ub2 * (+ 16- 4*kz)
	  + a1*htau4 * (+ 36*ml3- 7.2*ml3*kz+ 36*me3- 7.2*me3*kz+ 36*mH1sq- 7.2*mH1sq*kz)
	  + a1*hb*htau*Ub*Utau * (+ 32- 8*kz)
	  + a1*hb2*Utau2 * (+ 16- 4*kz)
	  + a1*hb2*htau2 * (+ 16*ml3- 4*ml3*kz+ 16*me3- 4*me3*kz+ 16*md3- 4*md3*kz+ 16*mq3- 4*mq3*kz+ 32*mH1sq- 8*mH1sq*kz)
	  + a1*M1*htau3*Utau * (- 72+ 14.4*kz)
	  + a1*M1*hb*htau2*Ub * (- 32+ 8*kz)
	  + a1*M1*hb2*htau*Utau * (- 32+ 8*kz)
	  + a1*M12*htau4 * (+ 36- 7.2*kz)
	  + a1*M12*hb2*htau2 * (+ 32- 8*kz)
	  + a1*a2*Utau2 * (- 16.2+ 5.4*kz)
	  + a1*a2*htau2 * (- 16.2*ml3+ 5.4*ml3*kz- 16.2*me3+ 5.4*me3*kz- 16.2*mH1sq+ 5.4*mH1sq*kz)
	  + a1*a2*M2*htau*Utau * (+ 32.4- 10.8*kz)
	  + a1*a2*M22*htau2 * (- 32.4+ 10.8*kz)
	  + a1*a2*M1*htau*Utau * (+ 32.4- 10.8*kz)
	  + a1*a2*M1*M2*htau2 * (- 32.4+ 10.8*kz)
	  + a1*a2*M12*htau2 * (- 32.4+ 10.8*kz)
	  + a1*a22 * (- 1.8*ml3- 5.4*mq3- 1.8*mH2sq- 1.8*mH1sq)
	  + a1*a22*M22 * (+ 34.2- 16.2*kz)
	  + a1*a22*M1*M2 * (+ 18- 10.8*kz)
	  + a1*a22*M12 * (+ 14.4- 5.4*kz)
	  + a12*Utau2 * (- 61.38+ 0.54*kz)
	  + a12*Ut2 * (- 9.36)
	  + a12*Ub2 * (- 5.04)
	  + a12*htau2 * (- 62.82*ml3+ 0.54*ml3*kz- 64.26*me3+ 0.54*me3*kz- 3.84*mu3- 0.96*md3- 0.48*mq3- 1.44*mH2sq- 62.82*mH1sq+ 0.54*mH1sq*kz)
	  + a12*ht2 * (- 9.36*mu3- 9.36*mq3- 9.36*mH2sq)
	  + a12*hb2 * (- 5.04*md3- 5.04*mq3- 5.04*mH1sq)
	  + a12*M1*htau*Utau * (+ 241.2- 2.16*kz)
	  + a12*M1*ht*Ut * (+ 31.2)
	  + a12*M1*hb*Ub * (+ 16.8)
	  + a12*M12*htau2 * (- 361.8+ 3.24*kz)
	  + a12*M12*ht2 * (- 46.8)
	  + a12*M12*hb2 * (- 25.2)
	  + a12*a3*M32 * (+ 63.36- 10.56*kz)
	  + a12*a3*M1*M3 * (+ 105.6- 21.12*kz)
	  + a12*a3*M12 * (+ 158.4- 31.68*kz)
	  + a12*a2 * (- 1.08*ml3- 2.16*me3- 2.88*mu3- 0.72*md3- 0.36*mq3- 1.08*mH2sq- 1.08*mH1sq)
	  + a12*a2*M22 * (+ 4.32- 3.24*kz)
	  + a12*a2*M1*M2 * (+ 2.16- 6.48*kz)
	  + a12*a2*M12 * (+ 3.24- 9.72*kz)
	  + a13 * (- 4.968*ml3- 9.936*me3- 13.248*mu3- 3.312*md3- 1.656*mq3- 4.968*mH2sq- 4.968*mH1sq)
	  + a13*M12 * (+ 446.136- 28.656*kz);
   
    dmtau3 = 
      + htau4*Utau2 * (+ 324+ 36*kz)
      + htau4*Ub2 * (- 12)
      + htau6 * (+ 108*ml3+ 12*ml3*kz+ 108*me3+ 12*me3*kz+ 108*mH1sq+ 12*mH1sq*kz)
      + ht2*htau2*Ub2 * (+ 24)
      + hb*htau3*Ub*Utau * (- 48)
      + hb*ht*htau2*Ub*Ut * (+ 48)
      + hb*ht2*htau*Ub*Utau * (+ 48)
      + hb2*htau2*Utau2 * (- 48)
      + hb2*htau2*Ut2 * (+ 24)
      + hb2*htau2*Ub2 * (+ 144)
      + hb2*htau4 * (- 24*ml3- 24*me3- 12*md3- 12*mq3- 36*mH1sq)
      + hb2*ht*htau*Ut*Utau * (+ 48)
      + hb2*ht2*Utau2 * (+ 24)
      + hb2*ht2*htau2 * (+ 24*ml3+ 24*me3+ 24*mu3+ 24*md3+ 48*mq3+ 24*mH2sq+ 48*mH1sq)
      + hb3*htau*Ub*Utau * (+ 144)
      + hb4*Utau2 * (+ 36)
      + hb4*htau2 * (+ 36*ml3+ 36*me3+ 72*md3+ 72*mq3+ 108*mH1sq)
      + a3*htau2*Ub2 * (- 128+ 32*kz)
      + a3*hb*htau*Ub*Utau * (- 256+ 64*kz)
      + a3*hb2*Utau2 * (- 128+ 32*kz)
      + a3*hb2*htau2 * (- 128*ml3+ 32*ml3*kz- 128*me3+ 32*me3*kz- 128*md3+ 32*md3*kz- 128*mq3+ 32*mq3*kz- 256*mH1sq+ 64*mH1sq*kz)
      + a3*M3*hb*htau2*Ub * (+ 256- 64*kz)
      + a3*M3*hb2*htau*Utau * (+ 256- 64*kz)
      + a3*M32*hb2*htau2 * (- 256+ 64*kz)
      + a2*htau2*Utau2 * (+ 144- 48*kz)
      + a2*htau2*Ub2 * (+ 54- 18*kz)
      + a2*htau4 * (+ 72*ml3- 24*ml3*kz+ 72*me3- 24*me3*kz+ 72*mH1sq- 24*mH1sq*kz)
      + a2*hb*htau*Ub*Utau * (+ 108- 36*kz)
      + a2*hb2*Utau2 * (+ 54- 18*kz)
      + a2*hb2*htau2 * (+ 54*ml3- 18*ml3*kz+ 54*me3- 18*me3*kz+ 54*md3- 18*md3*kz+ 54*mq3- 18*mq3*kz+ 108*mH1sq- 36*mH1sq*kz)
      + a2*M2*htau3*Utau * (- 144+ 48*kz)
      + a2*M2*hb*htau2*Ub * (- 108+ 36*kz)
      + a2*M2*hb2*htau*Utau * (- 108+ 36*kz)
      + a2*M22*htau4 * (+ 72- 24*kz)
      + a2*M22*hb2*htau2 * (+ 108- 36*kz)
      + a22*Utau2 * (- 87- 3*kz)
      + a22*htau2 * (- 99*ml3- 3*ml3*kz- 87*me3- 3*me3*kz- 36*mq3- 12*mH2sq- 99*mH1sq- 3*mH1sq*kz)
      + a22*M2*htau*Utau * (+ 348+ 12*kz)
      + a22*M22*htau2 * (- 474- 18*kz)
      + a1*htau2*Utau2 * (+ 28.8+ 28.8*kz)
      + a1*htau2*Ub2 * (+ 42.8+ 2.8*kz)
      + a1*htau4 * (+ 14.4*ml3+ 14.4*ml3*kz+ 14.4*me3+ 14.4*me3*kz+ 14.4*mH1sq+ 14.4*mH1sq*kz)
      + a1*hb*htau*Ub*Utau * (+ 85.6+ 5.6*kz)
      + a1*hb2*Utau2 * (+ 42.8+ 2.8*kz)
      + a1*hb2*htau2 * (+ 42.8*ml3+ 2.8*ml3*kz+ 42.8*me3+ 2.8*me3*kz+ 42.8*md3+ 2.8*md3*kz+ 42.8*mq3+ 2.8*mq3*kz+ 85.6*mH1sq+ 5.6*mH1sq*kz)
      + a1*M1*htau3*Utau * (- 28.8- 28.8*kz)
      + a1*M1*hb*htau2*Ub * (- 85.6- 5.6*kz)
      + a1*M1*hb2*htau*Utau * (- 85.6- 5.6*kz)
      + a1*M12*htau4 * (+ 14.4+ 14.4*kz)
      + a1*M12*hb2*htau2 * (+ 85.6+ 5.6*kz)
      + a1*a2*Utau2 * (- 54+ 10.8*kz)
      + a1*a2*htau2 * (- 54*ml3+ 10.8*ml3*kz- 54*me3+ 10.8*me3*kz- 54*mH1sq+ 10.8*mH1sq*kz)
      + a1*a2*M2*htau*Utau * (+ 108- 21.6*kz)
      + a1*a2*M22*htau2 * (- 108+ 21.6*kz)
      + a1*a2*M1*htau*Utau * (+ 108- 21.6*kz)
      + a1*a2*M1*M2*htau2 * (- 108+ 21.6*kz)
      + a1*a2*M12*htau2 * (- 108+ 21.6*kz)
      + a12*Utau2 * (- 86.04- 5.4*kz)
      + a12*Ut2 * (- 37.44)
      + a12*Ub2 * (- 20.16)
      + a12*htau2 * (- 84.6*ml3- 5.4*ml3*kz- 83.16*me3- 5.4*me3*kz+ 3.84*mu3+ 0.96*md3+ 0.48*mq3+ 1.44*mH2sq- 84.6*mH1sq- 5.4*mH1sq*kz)
      + a12*ht2 * (- 37.44*mu3- 37.44*mq3- 37.44*mH2sq)
      + a12*hb2 * (- 20.16*md3- 20.16*mq3- 20.16*mH1sq)
      + a12*M1*htau*Utau * (+ 326.88+ 21.6*kz)
      + a12*M1*ht*Ut * (+ 124.8)
      + a12*M1*hb*Ub * (+ 67.2)
      + a12*M12*htau2 * (- 490.32- 32.4*kz)
      + a12*M12*ht2 * (- 187.2)
      + a12*M12*hb2 * (- 100.8)
      + a12*a3*M32 * (+ 253.44- 42.24*kz)
      + a12*a3*M1*M3 * (+ 422.4- 84.48*kz)
      + a12*a3*M12 * (+ 633.6- 126.72*kz)
      + a12*a2*M22 * (+ 77.76- 12.96*kz)
      + a12*a2*M1*M2 * (+ 129.6- 25.92*kz)
      + a12*a2*M12 * (+ 194.4- 38.88*kz)
      + a13 * (- 22.464*ml3- 44.928*me3- 59.904*mu3- 14.976*md3- 7.488*mq3- 22.464*mH2sq- 22.464*mH1sq)
      + a13*M12 * (+ 1535.71- 114.624*kz);
    
    dmH1sq3 =
      + htau4*Utau2 * (+ 126+ 18*kz)
      + htau4*Ub2 * (+ 36)
      + htau6 * (+ 42*ml3+ 6*ml3*kz+ 42*me3+ 6*me3*kz+ 42*mH1sq+ 6*mH1sq*kz)
      + ht4*Ub2 * (+ 54)
      + hb*htau3*Ub*Utau * (+ 144)
      + hb*ht3*Ub*Ut * (+ 216)
      + hb2*htau2*Utau2 * (+ 144)
      + hb2*htau2*Ub2 * (+ 144)
      + hb2*htau4 * (+ 72*ml3+ 72*me3+ 36*md3+ 36*mq3+ 108*mH1sq)
      + hb2*ht2*Ut2 * (+ 216)
      + hb2*ht4 * (+ 108*mu3+ 54*md3+ 162*mq3+ 108*mH2sq+ 54*mH1sq)
      + hb3*htau*Ub*Utau * (+ 144)
      + hb4*Utau2 * (+ 36)
      + hb4*Ub2 * (+ 1026+ 54*kz)
      + hb4*htau2 * (+ 36*ml3+ 36*me3+ 72*md3+ 72*mq3+ 108*mH1sq)
      + hb6 * (+ 342*md3+ 18*md3*kz+ 342*mq3+ 18*mq3*kz+ 342*mH1sq+ 18*mH1sq*kz)
      + a3*ht2*Ub2 * (+ 48- 16*kz)
      + a3*hb*ht*Ub*Ut * (+ 96- 32*kz)
      + a3*hb2*Ut2 * (+ 48- 16*kz)
      + a3*hb2*Ub2 * (+ 576- 192*kz)
      + a3*hb2*ht2 * (+ 48*mu3- 16*mu3*kz+ 48*md3- 16*md3*kz+ 96*mq3- 32*mq3*kz+ 48*mH2sq- 16*mH2sq*kz+ 48*mH1sq- 16*mH1sq*kz)
      + a3*hb4 * (+ 288*md3- 96*md3*kz+ 288*mq3- 96*mq3*kz+ 288*mH1sq- 96*mH1sq*kz)
      + a3*M3*hb*ht2*Ub * (- 96+ 32*kz)
      + a3*M3*hb2*ht*Ut * (- 96+ 32*kz)
      + a3*M3*hb3*Ub * (- 576+ 192*kz)
      + a3*M32*hb2*ht2 * (+ 96- 32*kz)
      + a3*M32*hb4 * (+ 288- 96*kz)
      + a32*Ub2 * (- 320*O3- 16*O3*kz)
      + a32*hb2 * (- 32*mu3- 416*O3*md3- 16*O3*md3*kz- 512*O3*mq3- 16*O3*mq3*kz- 320*O3*mH1sq- 16*O3*mH1sq*kz)
      + a32*M3*hb*Ub * (+ 1280*O3+ 64*O3*kz)
      + a32*M32*hb2 * (- 448- 32*kz)
      + a2*htau2*Utau2 * (+ 24+ 24*kz)
      + a2*htau4 * (+ 12*ml3+ 12*ml3*kz+ 12*me3+ 12*me3*kz+ 12*mH1sq+ 12*mH1sq*kz)
      + a2*ht2*Ub2 * (+ 36)
      + a2*hb*ht*Ub*Ut * (+ 72)
      + a2*hb2*Ut2 * (+ 36)
      + a2*hb2*Ub2 * (+ 72+ 72*kz)
      + a2*hb2*ht2 * (+ 36*mu3+ 36*md3+ 72*mq3+ 36*mH2sq+ 36*mH1sq)
      + a2*hb4 * (+ 36*md3+ 36*md3*kz+ 36*mq3+ 36*mq3*kz+ 36*mH1sq+ 36*mH1sq*kz)
      + a2*M2*htau3*Utau * (- 24- 24*kz)
      + a2*M2*hb*ht2*Ub * (- 72)
      + a2*M2*hb2*ht*Ut * (- 72)
      + a2*M2*hb3*Ub * (- 72- 72*kz)
      + a2*M22*htau4 * (+ 12+ 12*kz)
      + a2*M22*hb2*ht2 * (+ 72)
      + a2*M22*hb4 * (+ 36+ 36*kz)
      + a2*a3*Ub2 * (- 264+ 48*kz)
      + a2*a3*hb2 * (- 264*md3+ 48*md3*kz- 264*mq3+ 48*mq3*kz- 264*mH1sq+ 48*mH1sq*kz)
      + a2*a3*M3*hb*Ub * (+ 528- 96*kz)
      + a2*a3*M32*hb2 * (- 528+ 96*kz)
      + a2*a3*M2*hb*Ub * (+ 528- 96*kz)
      + a2*a3*M2*M3*hb2 * (- 528+ 96*kz)
      + a2*a3*M22*hb2 * (- 528+ 96*kz)
      + a22*Utau2 * (- 40.5- 10.5*kz)
      + a22*Ut2 * (- 54)
      + a22*Ub2 * (- 121.5- 31.5*kz)
      + a22*htau2 * (- 40.5*ml3- 10.5*ml3*kz- 40.5*me3- 10.5*me3*kz- 40.5*mH1sq- 10.5*mH1sq*kz)
      + a22*ht2 * (- 54*mu3- 54*mq3- 54*mH2sq)
      + a22*hb2 * (- 121.5*md3- 31.5*md3*kz- 121.5*mq3- 31.5*mq3*kz- 121.5*mH1sq- 31.5*mH1sq*kz)
      + a22*M2*htau*Utau * (+ 150+ 42*kz)
      + a22*M2*ht*Ut * (+ 180)
      + a22*M2*hb*Ub * (+ 450+ 126*kz)
      + a22*M22*htau2 * (- 225- 63*kz)
      + a22*M22*ht2 * (- 270)
      + a22*M22*hb2 * (- 675- 189*kz)
      + a22*a3*M32 * (+ 432- 72*kz)
      + a22*a3*M2*M3 * (+ 720- 144*kz)
      + a22*a3*M22 * (+ 1080- 216*kz)
      + a23 * (- 3*ml3- 9*mq3- 3*mH2sq- 3*mH1sq)
      + a23*M22 * (+ 2157+ 630*kz)
      + a1*htau2*Utau2 * (+ 72- 14.4*kz)
      + a1*htau4 * (+ 36*ml3- 7.2*ml3*kz+ 36*me3- 7.2*me3*kz+ 36*mH1sq- 7.2*mH1sq*kz)
      + a1*ht2*Ub2 * (- 4.8+ 2.8*kz)
      + a1*hb*ht*Ub*Ut * (- 9.6+ 5.6*kz)
      + a1*hb2*Ut2 * (- 4.8+ 2.8*kz)
      + a1*hb2*Ub2 * (+ 24+ 14.4*kz)
      + a1*hb2*ht2 * (- 4.8*mu3+ 2.8*mu3*kz- 4.8*md3+ 2.8*md3*kz- 9.6*mq3+ 5.6*mq3*kz- 4.8*mH2sq+ 2.8*mH2sq*kz- 4.8*mH1sq+ 2.8*mH1sq*kz)
      + a1*hb4 * (+ 12*md3+ 7.2*md3*kz+ 12*mq3+ 7.2*mq3*kz+ 12*mH1sq+ 7.2*mH1sq*kz)
      + a1*M1*htau3*Utau * (- 72+ 14.4*kz)
      + a1*M1*hb*ht2*Ub * (+ 9.6- 5.6*kz)
      + a1*M1*hb2*ht*Ut * (+ 9.6- 5.6*kz)
      + a1*M1*hb3*Ub * (- 24- 14.4*kz)
      + a1*M12*htau4 * (+ 36- 7.2*kz)
      + a1*M12*hb2*ht2 * (- 9.6+ 5.6*kz)
      + a1*M12*hb4 * (+ 12+ 7.2*kz)
      + a1*a3*Ub2 * (- 568*O15+ 112*O15*kz)
      + a1*a3*hb2 * (- 568*O15*md3+ 112*O15*md3*kz- 568*O15*mq3+ 112*O15*mq3*kz- 568*O15*mH1sq+ 112*O15*mH1sq*kz)
      + a1*a3*M3*hb*Ub * (+ 1136*O15- 224*O15*kz)
      + a1*a3*M32*hb2 * (- 1136*O15+ 224*O15*kz)
      + a1*a3*M1*hb*Ub * (+ 1136*O15- 224*O15*kz)
      + a1*a3*M1*M3*hb2 * (- 1136*O15+ 224*O15*kz)
      + a1*a3*M12*hb2 * (- 1136*O15+ 224*O15*kz)
      + a1*a2*Utau2 * (- 16.2+ 5.4*kz)
      + a1*a2*Ub2 * (- 0.6- 3*kz)
      + a1*a2*htau2 * (- 16.2*ml3+ 5.4*ml3*kz- 16.2*me3+ 5.4*me3*kz- 16.2*mH1sq+ 5.4*mH1sq*kz)
      + a1*a2*hb2 * (- 0.6*md3- 3*md3*kz- 0.6*mq3- 3*mq3*kz- 0.6*mH1sq- 3*mH1sq*kz)
      + a1*a2*M2*htau*Utau * (+ 32.4- 10.8*kz)
      + a1*a2*M2*hb*Ub * (+ 1.2+ 6*kz)
      + a1*a2*M22*htau2 * (- 32.4+ 10.8*kz)
      + a1*a2*M22*hb2 * (- 1.2- 6*kz)
      + a1*a2*M1*htau*Utau * (+ 32.4- 10.8*kz)
      + a1*a2*M1*hb*Ub * (+ 1.2+ 6*kz)
      + a1*a2*M1*M2*htau2 * (- 32.4+ 10.8*kz)
      + a1*a2*M1*M2*hb2 * (- 1.2- 6*kz)
      + a1*a2*M12*htau2 * (- 32.4+ 10.8*kz)
      + a1*a2*M12*hb2 * (- 1.2- 6*kz)
      + a1*a22 * (- 1.8*ml3- 5.4*mq3- 1.8*mH2sq- 1.8*mH1sq)
      + a1*a22*M22 * (+ 34.2- 16.2*kz)
      + a1*a22*M1*M2 * (+ 18- 10.8*kz)
      + a1*a22*M12 * (+ 14.4- 5.4*kz)
      + a12*Utau2 * (- 61.38+ 0.54*kz)
      + a12*Ut2 * (- 9.36)
      + a12*Ub2 * (- 4501*O150- 77*O150*kz)
      + a12*htau2 * (- 62.82*ml3+ 0.54*ml3*kz- 64.26*me3+ 0.54*me3*kz- 3.84*mu3- 0.96*md3- 0.48*mq3- 1.44*mH2sq- 62.82*mH1sq+ 0.54*mH1sq*kz)
      + a12*ht2 * (- 9.36*mu3- 9.36*mq3- 9.36*mH2sq)
      + a12*hb2 * (+ 0.48*ml3+ 0.96*me3+ 1.28*mu3- 4453*O150*md3- 77*O150*md3*kz- 4477*O150*mq3- 77*O150*mq3*kz+ 0.48*mH2sq- 4429*O150*mH1sq- 77*O150*mH1sq*kz)
      + a12*M1*htau*Utau * (+ 241.2- 2.16*kz)
      + a12*M1*ht*Ut * (+ 31.2)
      + a12*M1*hb*Ub * (+ 350*O3+ 154*O75*kz)
      + a12*M12*htau2 * (- 361.8+ 3.24*kz)
      + a12*M12*ht2 * (- 46.8)
      + a12*M12*hb2 * (- 175- 3.08*kz)
      + a12*a3*M32 * (+ 63.36- 10.56*kz)
      + a12*a3*M1*M3 * (+ 105.6- 21.12*kz)
      + a12*a3*M12 * (+ 158.4- 31.68*kz)
      + a12*a2 * (- 1.08*ml3- 2.16*me3- 2.88*mu3- 0.72*md3- 0.36*mq3- 1.08*mH2sq- 1.08*mH1sq)
      + a12*a2*M22 * (+ 4.32- 3.24*kz)
      + a12*a2*M1*M2 * (+ 2.16- 6.48*kz)
      + a12*a2*M12 * (+ 3.24- 9.72*kz)
      + a13 * (- 4.968*ml3- 9.936*me3- 13.248*mu3- 3.312*md3- 1.656*mq3- 4.968*mH2sq- 4.968*mH1sq)
      + a13*M12 * (+ 446.136- 28.656*kz);
    
    
    dmH2sq3 = 
      + ht2*htau2*Ub2 * (+ 12)
      + ht4*Ut2 * (+ 1026+ 54*kz)
      + ht6 * (+ 342*mu3+ 18*mu3*kz+ 342*mq3+ 18*mq3*kz+ 342*mH2sq+ 18*mH2sq*kz)
      + hb*ht*htau2*Ub*Ut * (+ 24)
      + hb*ht2*htau*Ub*Utau * (+ 24)
      + hb2*htau2*Ut2 * (+ 12)
      + hb2*ht*htau*Ut*Utau * (+ 24)
      + hb2*ht2*Utau2 * (+ 12)
      + hb2*ht2*Ub2 * (+ 216)
      + hb2*ht2*htau2 * (+ 12*ml3+ 12*me3+ 12*mu3+ 12*md3+ 24*mq3+ 12*mH2sq+ 24*mH1sq)
      + hb3*ht*Ub*Ut * (+ 216)
      + hb4*Ut2 * (+ 54)
      + hb4*ht2 * (+ 54*mu3+ 108*md3+ 162*mq3+ 54*mH2sq+ 108*mH1sq)
      + a3*ht2*Ut2 * (+ 576- 192*kz)
      + a3*ht2*Ub2 * (+ 48- 16*kz)
      + a3*ht4 * (+ 288*mu3- 96*mu3*kz+ 288*mq3- 96*mq3*kz+ 288*mH2sq- 96*mH2sq*kz)
      + a3*hb*ht*Ub*Ut * (+ 96- 32*kz)
      + a3*hb2*Ut2 * (+ 48- 16*kz)
      + a3*hb2*ht2 * (+ 48*mu3- 16*mu3*kz+ 48*md3- 16*md3*kz+ 96*mq3- 32*mq3*kz+ 48*mH2sq- 16*mH2sq*kz+ 48*mH1sq- 16*mH1sq*kz)
      + a3*M3*ht3*Ut * (- 576+ 192*kz)
      + a3*M3*hb*ht2*Ub * (- 96+ 32*kz)
      + a3*M3*hb2*ht*Ut * (- 96+ 32*kz)
      + a3*M32*ht4 * (+ 288- 96*kz)
      + a3*M32*hb2*ht2 * (+ 96- 32*kz)
      + a32*Ut2 * (- 320*O3- 16*O3*kz)
      + a32*ht2 * (- 416*O3*mu3- 16*O3*mu3*kz- 32*md3- 512*O3*mq3- 16*O3*mq3*kz- 320*O3*mH2sq- 16*O3*mH2sq*kz)
      + a32*M3*ht*Ut * (+ 1280*O3+ 64*O3*kz)
      + a32*M32*ht2 * (- 448- 32*kz)
      + a2*ht2*Ut2 * (+ 72+ 72*kz)
      + a2*ht2*Ub2 * (+ 36)
      + a2*ht4 * (+ 36*mu3+ 36*mu3*kz+ 36*mq3+ 36*mq3*kz+ 36*mH2sq+ 36*mH2sq*kz)
      + a2*hb*ht*Ub*Ut * (+ 72)
      + a2*hb2*Ut2 * (+ 36)
      + a2*hb2*ht2 * (+ 36*mu3+ 36*md3+ 72*mq3+ 36*mH2sq+ 36*mH1sq)
      + a2*M2*ht3*Ut * (- 72- 72*kz)
      + a2*M2*hb*ht2*Ub * (- 72)
      + a2*M2*hb2*ht*Ut * (- 72)
      + a2*M22*ht4 * (+ 36+ 36*kz)
      + a2*M22*hb2*ht2 * (+ 72)
      + a2*a3*Ut2 * (- 264+ 48*kz)
      + a2*a3*ht2 * (- 264*mu3+ 48*mu3*kz- 264*mq3+ 48*mq3*kz- 264*mH2sq+ 48*mH2sq*kz)
      + a2*a3*M3*ht*Ut * (+ 528- 96*kz)
      + a2*a3*M32*ht2 * (- 528+ 96*kz)
      + a2*a3*M2*ht*Ut * (+ 528- 96*kz)
      + a2*a3*M2*M3*ht2 * (- 528+ 96*kz)
      + a2*a3*M22*ht2 * (- 528+ 96*kz)
      + a22*Utau2 * (- 18)
      + a22*Ut2 * (- 121.5- 31.5*kz)
      + a22*Ub2 * (- 54)
      + a22*htau2 * (- 18*ml3- 18*me3- 18*mH1sq)
      + a22*ht2 * (- 121.5*mu3- 31.5*mu3*kz- 121.5*mq3- 31.5*mq3*kz- 121.5*mH2sq- 31.5*mH2sq*kz)
      + a22*hb2 * (- 54*md3- 54*mq3- 54*mH1sq)
      + a22*M2*htau*Utau * (+ 60)
      + a22*M2*ht*Ut * (+ 450+ 126*kz)
      + a22*M2*hb*Ub * (+ 180)
      + a22*M22*htau2 * (- 90)
      + a22*M22*ht2 * (- 675- 189*kz)
      + a22*M22*hb2 * (- 270)
      + a22*a3*M32 * (+ 432- 72*kz)
      + a22*a3*M2*M3 * (+ 720- 144*kz)
      + a22*a3*M22 * (+ 1080- 216*kz)
      + a23 * (- 3*ml3- 9*mq3- 3*mH2sq- 3*mH1sq)
      + a23*M22 * (+ 2157+ 630*kz)
      + a1*ht2*Ut2 * (+ 91.2- 4.8*kz)
      + a1*ht2*Ub2 * (+ 2.4+ 0.4*kz)
      + a1*ht4 * (+ 45.6*mu3- 2.4*mu3*kz+ 45.6*mq3- 2.4*mq3*kz+ 45.6*mH2sq- 2.4*mH2sq*kz)
      + a1*hb*ht*Ub*Ut * (+ 4.8+ 0.8*kz)
      + a1*hb2*Ut2 * (+ 2.4+ 0.4*kz)
      + a1*hb2*ht2 * (+ 2.4*mu3+ 0.4*mu3*kz+ 2.4*md3+ 0.4*md3*kz+ 4.8*mq3+ 0.8*mq3*kz+ 2.4*mH2sq+ 0.4*mH2sq*kz+ 2.4*mH1sq+ 0.4*mH1sq*kz)
      + a1*M1*ht3*Ut * (- 91.2+ 4.8*kz)
      + a1*M1*hb*ht2*Ub * (- 4.8- 0.8*kz)
      + a1*M1*hb2*ht*Ut * (- 4.8- 0.8*kz)
      + a1*M12*ht4 * (+ 45.6- 2.4*kz)
      + a1*M12*hb2*ht2 * (+ 4.8+ 0.8*kz)
      + a1*a3*Ut2 * (- 248*O3+ 208*O15*kz)
      + a1*a3*ht2 * (- 248*O3*mu3+ 208*O15*mu3*kz- 248*O3*mq3+ 208*O15*mq3*kz- 248*O3*mH2sq+ 208*O15*mH2sq*kz)
      + a1*a3*M3*ht*Ut * (+ 496*O3- 416*O15*kz)
      + a1*a3*M32*ht2 * (- 496*O3+ 416*O15*kz)
      + a1*a3*M1*ht*Ut * (+ 496*O3- 416*O15*kz)
      + a1*a3*M1*M3*ht2 * (- 496*O3+ 416*O15*kz)
      + a1*a3*M12*ht2 * (- 496*O3+ 416*O15*kz)
      + a1*a2*Ut2 * (- 11.4+ 4.2*kz)
      + a1*a2*ht2 * (- 11.4*mu3+ 4.2*mu3*kz- 11.4*mq3+ 4.2*mq3*kz- 11.4*mH2sq+ 4.2*mH2sq*kz)
      + a1*a2*M2*ht*Ut * (+ 22.8- 8.4*kz)
      + a1*a2*M22*ht2 * (- 22.8+ 8.4*kz)
      + a1*a2*M1*ht*Ut * (+ 22.8- 8.4*kz)
      + a1*a2*M1*M2*ht2 * (- 22.8+ 8.4*kz)
      + a1*a2*M12*ht2 * (- 22.8+ 8.4*kz)
      + a1*a22 * (- 1.8*ml3- 5.4*mq3- 1.8*mH2sq- 1.8*mH1sq)
      + a1*a22*M22 * (+ 34.2- 16.2*kz)
      + a1*a22*M1*M2 * (+ 18- 10.8*kz)
      + a1*a22*M12 * (+ 14.4- 5.4*kz)
      + a12*Utau2 * (- 6.48)
      + a12*Ut2 * (- 10849*O150- 13*O30*kz)
      + a12*Ub2 * (- 5.04)
      + a12*htau2 * (- 6.48*ml3- 6.48*me3- 6.48*mH1sq)
      + a12*ht2 * (- 0.96*ml3- 1.92*me3- 11233*O150*mu3- 13*O30*mu3*kz- 0.64*md3- 10897*O150*mq3- 13*O30*mq3*kz- 10993*O150*mH2sq- 13*O30*mH2sq*kz- 0.96*mH1sq)
      + a12*hb2 * (- 5.04*md3- 5.04*mq3- 5.04*mH1sq)
      + a12*M1*htau*Utau * (+ 21.6)
      + a12*M1*ht*Ut * (+ 4246*O15+ 26*O15*kz)
      + a12*M1*hb*Ub * (+ 16.8)
      + a12*M12*htau2 * (- 32.4)
      + a12*M12*ht2 * (- 424.6- 2.6*kz)
      + a12*M12*hb2 * (- 25.2)
      + a12*a3*M32 * (+ 63.36- 10.56*kz)
      + a12*a3*M1*M3 * (+ 105.6- 21.12*kz)
      + a12*a3*M12 * (+ 158.4- 31.68*kz)
      + a12*a2 * (- 1.08*ml3- 2.16*me3- 2.88*mu3- 0.72*md3- 0.36*mq3- 1.08*mH2sq- 1.08*mH1sq)
      + a12*a2*M22 * (+ 4.32- 3.24*kz)
      + a12*a2*M1*M2 * (+ 2.16- 6.48*kz)
      + a12*a2*M12 * (+ 3.24- 9.72*kz)
      + a13 * (- 4.968*ml3- 9.936*me3- 13.248*mu3- 3.312*md3- 1.656*mq3- 4.968*mH2sq- 4.968*mH1sq)
      + a13*M12 * (+ 446.136- 28.656*kz);
    
    dm3sq3 = 
      + htau5*Utau*smu_ * (+ 42+ 6*kz)
      + htau6*m3sq * (+ 7+ kz)
      + ht5*Ut*smu_ * (+ 342+ 18*kz)
      + ht6*m3sq * (+ 57+ 3*kz)
      + hb*htau4*Ub*smu_ * (+ 36)
      + hb*ht2*htau2*Ub*smu_ * (+ 12)
      + hb*ht4*Ub*smu_ * (+ 54)
      + hb2*htau3*Utau*smu_ * (+ 72)
      + hb2*htau4*m3sq * (+ 18)
      + hb2*ht*htau2*Ut*smu_ * (+ 12)
      + hb2*ht2*htau*Utau*smu_ * (+ 12)
      + hb2*ht2*htau2*m3sq * (+ 6)
      + hb2*ht3*Ut*smu_ * (+ 108)
      + hb2*ht4*m3sq * (+ 27)
      + hb3*htau2*Ub*smu_ * (+ 72)
      + hb3*ht2*Ub*smu_ * (+ 108)
      + hb4*htau*Utau*smu_ * (+ 36)
      + hb4*htau2*m3sq * (+ 18)
      + hb4*ht*Ut*smu_ * (+ 54)
      + hb4*ht2*m3sq * (+ 27)
      + hb5*Ub*smu_ * (+ 342+ 18*kz)
      + hb6*m3sq * (+ 57+ 3*kz)
      + a3*ht3*Ut*smu_ * (+ 288- 96*kz)
      + a3*ht4*m3sq * (+ 72- 24*kz)
      + a3*hb*ht2*Ub*smu_ * (+ 96- 32*kz)
      + a3*hb2*ht*Ut*smu_ * (+ 96- 32*kz)
      + a3*hb2*ht2*m3sq * (+ 48- 16*kz)
      + a3*hb3*Ub*smu_ * (+ 288- 96*kz)
      + a3*hb4*m3sq * (+ 72- 24*kz)
      + a3*M3*ht4*smu_ * (- 144+ 48*kz)
      + a3*M3*hb2*ht2*smu_ * (- 96+ 32*kz)
      + a3*M3*hb4*smu_ * (- 144+ 48*kz)
      + a32*ht*Ut*smu_ * (- 320*O3- 16*O3*kz)
      + a32*ht2*m3sq * (- 160*O3- 8*O3*kz)
      + a32*hb*Ub*smu_ * (- 320*O3- 16*O3*kz)
      + a32*hb2*m3sq * (- 160*O3- 8*O3*kz)
      + a32*M3*ht2*smu_ * (+ 640*O3+ 32*O3*kz)
      + a32*M3*hb2*smu_ * (+ 640*O3+ 32*O3*kz)
      + a2*htau3*Utau*smu_ * (+ 12+ 12*kz)
      + a2*htau4*m3sq * (+ 3+ 3*kz)
      + a2*ht3*Ut*smu_ * (+ 36+ 36*kz)
      + a2*ht4*m3sq * (+ 9+ 9*kz)
      + a2*hb*ht2*Ub*smu_ * (+ 72)
      + a2*hb2*ht*Ut*smu_ * (+ 72)
      + a2*hb2*ht2*m3sq * (+ 36)
      + a2*hb3*Ub*smu_ * (+ 36+ 36*kz)
      + a2*hb4*m3sq * (+ 9+ 9*kz)
      + a2*M2*htau4*smu_ * (- 6- 6*kz)
      + a2*M2*ht4*smu_ * (- 18- 18*kz)
      + a2*M2*hb2*ht2*smu_ * (- 72)
      + a2*M2*hb4*smu_ * (- 18- 18*kz)
      + a2*a3*ht*Ut*smu_ * (- 264+ 48*kz)
      + a2*a3*ht2*m3sq * (- 132+ 24*kz)
      + a2*a3*hb*Ub*smu_ * (- 264+ 48*kz)
      + a2*a3*hb2*m3sq * (- 132+ 24*kz)
      + a2*a3*M3*ht2*smu_ * (+ 264- 48*kz)
      + a2*a3*M3*hb2*smu_ * (+ 264- 48*kz)
      + a2*a3*M2*ht2*smu_ * (+ 264- 48*kz)
      + a2*a3*M2*hb2*smu_ * (+ 264- 48*kz)
      + a22*htau*Utau*smu_ * (- 52.5- 10.5*kz)
      + a22*htau2*m3sq * (- 26.25- 5.25*kz)
      + a22*ht*Ut*smu_ * (- 157.5- 31.5*kz)
      + a22*ht2*m3sq * (- 78.75- 15.75*kz)
      + a22*hb*Ub*smu_ * (- 157.5- 31.5*kz)
      + a22*hb2*m3sq * (- 78.75- 15.75*kz)
      + a22*M2*htau2*smu_ * (+ 105+ 21*kz)
      + a22*M2*ht2*smu_ * (+ 315+ 63*kz)
      + a22*M2*hb2*smu_ * (+ 315+ 63*kz)
      + a22*a3*m3sq * (+ 180- 36*kz)
      + a22*a3*M3*smu_ * (- 360+ 72*kz)
      + a22*a3*M2*smu_ * (- 720+ 144*kz)
      + a23*m3sq * (+ 172.5+ 52.5*kz)
      + a23*M2*smu_ * (- 1035- 315*kz)
      + a1*htau3*Utau*smu_ * (+ 36- 7.2*kz)
      + a1*htau4*m3sq * (+ 9- 1.8*kz)
      + a1*ht3*Ut*smu_ * (+ 45.6- 2.4*kz)
      + a1*ht4*m3sq * (+ 11.4- 0.6*kz)
      + a1*hb*ht2*Ub*smu_ * (- 2.4+ 3.2*kz)
      + a1*hb2*ht*Ut*smu_ * (- 2.4+ 3.2*kz)
      + a1*hb2*ht2*m3sq * (- 1.2+ 1.6*kz)
      + a1*hb3*Ub*smu_ * (+ 12+ 7.2*kz)
      + a1*hb4*m3sq * (+ 3+ 1.8*kz)
      + a1*M1*htau4*smu_ * (- 18+ 3.6*kz)
      + a1*M1*ht4*smu_ * (- 22.8+ 1.2*kz)
      + a1*M1*hb2*ht2*smu_ * (+ 2.4- 3.2*kz)
      + a1*M1*hb4*smu_ * (- 6- 3.6*kz)
      + a1*a3*ht*Ut*smu_ * (- 248*O3+ 208*O15*kz)
      + a1*a3*ht2*m3sq * (- 124*O3+ 104*O15*kz)
      + a1*a3*hb*Ub*smu_ * (- 568*O15+ 112*O15*kz)
      + a1*a3*hb2*m3sq * (- 284*O15+ 56*O15*kz)
      + a1*a3*M3*ht2*smu_ * (+ 248*O3- 208*O15*kz)
      + a1*a3*M3*hb2*smu_ * (+ 568*O15- 112*O15*kz)
      + a1*a3*M1*ht2*smu_ * (+ 248*O3- 208*O15*kz)
      + a1*a3*M1*hb2*smu_ * (+ 568*O15- 112*O15*kz)
      + a1*a2*htau*Utau*smu_ * (- 16.2+ 5.4*kz)
      + a1*a2*htau2*m3sq * (- 8.1+ 2.7*kz)
      + a1*a2*ht*Ut*smu_ * (- 11.4+ 4.2*kz)
      + a1*a2*ht2*m3sq * (- 5.7+ 2.1*kz)
      + a1*a2*hb*Ub*smu_ * (- 0.6- 3*kz)
      + a1*a2*hb2*m3sq * (- 0.3- 1.5*kz)
      + a1*a2*M2*htau2*smu_ * (+ 16.2- 5.4*kz)
      + a1*a2*M2*ht2*smu_ * (+ 11.4- 4.2*kz)
      + a1*a2*M2*hb2*smu_ * (+ 0.6+ 3*kz)
       + a1*a2*M1*htau2*smu_ * (+ 16.2- 5.4*kz)
      + a1*a2*M1*ht2*smu_ * (+ 11.4- 4.2*kz)
      + a1*a2*M1*hb2*smu_ * (+ 0.6+ 3*kz)
      + a1*a22*m3sq * (+ 4.5- 2.7*kz)
      + a1*a22*M2*smu_ * (- 18+ 10.8*kz)
      + a1*a22*M1*smu_ * (- 9+ 5.4*kz)
      + a12*htau*Utau*smu_ * (- 65.7+ 0.54*kz)
      + a12*htau2*m3sq * (- 32.85+ 0.27*kz)
      + a12*ht*Ut*smu_ * (- 2357*O30- 13*O30*kz)
      + a12*ht2*m3sq * (- 2357*O60- 13*O60*kz)
      + a12*hb*Ub*smu_ * (- 1001*O30- 77*O150*kz)
      + a12*hb2*m3sq * (- 1001*O60- 77*O300*kz)
      + a12*M1*htau2*smu_ * (+ 131.4- 1.08*kz)
      + a12*M1*ht2*smu_ * (+ 2357*O15+ 13*O15*kz)
      + a12*M1*hb2*smu_ * (+ 1001*O15+ 77*O75*kz)
      + a12*a3*m3sq * (+ 26.4- 5.28*kz)
      + a12*a3*M3*smu_ * (- 52.8+ 10.56*kz)
      + a12*a3*M1*smu_ * (- 105.6+ 21.12*kz)
      + a12*a2*m3sq * (+ 0.54- 1.62*kz)
      + a12*a2*M2*smu_ * (- 1.08+ 3.24*kz)
      + a12*a2*M1*smu_ * (- 2.16+ 6.48*kz)
      + a13*m3sq * (+ 36.78- 2.388*kz)
      + a13*M1*smu_ * (- 220.68+ 14.328*kz);
    
    dht3 = 
      + ht6*Ut * (+ 714+ 42*kz)
      + hb*ht*htau4*Ub * (+ 10)
      + hb*ht3*htau2*Ub * (+ 20)
      + hb*ht5*Ub * (+ 24)
      + hb2*htau4*Ut * (+ 5)
      + hb2*ht*htau3*Utau * (+ 20)
      + hb2*ht2*htau2*Ut * (+ 30)
      + hb2*ht3*htau*Utau * (+ 20)
      + hb2*ht4*Ut * (+ 60)
      + hb3*ht*htau2*Ub * (- 16)
      + hb3*ht3*Ub * (+ 220)
      + hb4*htau2*Ut * (- 4)
      + hb4*ht*htau*Utau * (- 8)
      + hb4*ht2*Ut * (+ 165)
      + hb5*ht*Ub * (+ 90+ 6*kz)
      + hb6*Ut * (+ 15+ kz)
      + a3*ht4*Ut * (+ 1360*O3)
      + a3*hb*ht*htau2*Ub * (+ 16)
      + a3*hb*ht3*Ub * (+ 272*O3- 16*kz)
      + a3*hb2*htau2*Ut * (+ 8)
      + a3*hb2*ht*htau*Utau * (+ 16)
      + a3*hb2*ht2*Ut * (+ 136- 24*kz)
      + a3*hb3*ht*Ub * (+ 160*O3+ 32*kz)
      + a3*hb4*Ut * (+ 40*O3+ 8*kz)
      + a3*M3*ht5 * (- 544*O3)
      + a3*M3*hb2*ht*htau2 * (- 16)
      + a3*M3*hb2*ht3 * (- 272*O3+ 16*kz)
      + a3*M3*hb4*ht * (- 80*O3- 16*kz)
      + a32*ht2*Ut * (- 296- 144*kz)
      + a32*hb*ht*Ub * (- 304*O3- 272*O9*kz)
      + a32*hb2*Ut * (- 152*O3- 136*O9*kz)
      + a32*M3*ht3 * (+ 1184*O3+ 192*kz)
      + a32*M3*hb2*ht * (+ 608*O3+ 544*O9*kz)
      + a33*Ut * (+ 5440*O27+ 320*O3*kz)
      + a33*M3*ht * (- 10880*O9- 640*kz)
      + a2*ht4*Ut * (+ 300)
      + a2*hb*ht*htau2*Ub * (+ 12)
      + a2*hb*ht3*Ub * (+ 54- 6*kz)
      + a2*hb2*htau2*Ut * (+ 6)
      + a2*hb2*ht*htau*Utau * (+ 12)
      + a2*hb2*ht2*Ut * (+ 81- 9*kz)
      + a2*hb3*ht*Ub * (+ 60+ 12*kz)
      + a2*hb4*Ut * (+ 15+ 3*kz)
      + a2*M2*ht5 * (- 120)
      + a2*M2*hb2*ht*htau2 * (- 12)
      + a2*M2*hb2*ht3 * (- 54+ 6*kz)
      + a2*M2*hb4*ht * (- 30- 6*kz)
      + a2*a3*ht2*Ut * (- 672+ 120*kz)
      + a2*a3*hb*ht*Ub * (- 8)
      + a2*a3*hb2*Ut * (- 4)
      + a2*a3*M3*ht3 * (+ 448- 80*kz)
      + a2*a3*M3*hb2*ht * (+ 8)
      + a2*a3*M2*ht3 * (+ 448- 80*kz)
      + a2*a3*M2*hb2*ht * (+ 8)
      + a2*a32*Ut * (+ 68- 24*kz)
      + a2*a32*M3*ht * (- 272+ 96*kz)
      + a2*a32*M2*ht * (- 136+ 48*kz)
      + a22*htau2*Ut * (- 15)
      + a22*ht*htau*Utau * (- 30)
      + a22*ht2*Ut * (- 400.5- 67.5*kz)
      + a22*hb*ht*Ub * (- 112.5- 10.5*kz)
      + a22*hb2*Ut * (- 56.25- 5.25*kz)
      + a22*M2*ht*htau2 * (+ 60)
      + a22*M2*ht3 * (+ 534+ 90*kz)
      + a22*M2*hb2*ht * (+ 225+ 21*kz)
      + a22*a3*Ut * (+ 140- 36*kz)
      + a22*a3*M3*ht * (- 280+ 72*kz)
      + a22*a3*M2*ht * (- 560+ 144*kz)
      + a23*Ut * (+ 172.5+ 52.5*kz)
      + a23*M2*ht * (- 1035- 315*kz)
      + a1*ht4*Ut * (+ 356*O3)
      + a1*hb*ht*htau2*Ub * (- 3.2+ 0.8*kz)
      + a1*hb*ht3*Ub * (+ 74*O15+ 1.6*kz)
      + a1*hb2*htau2*Ut * (- 1.6+ 0.4*kz)
      + a1*hb2*ht*htau*Utau * (- 3.2+ 0.8*kz)
      + a1*hb2*ht2*Ut * (+ 7.4+ 2.4*kz)
      + a1*hb3*ht*Ub * (+ 44*O3- 4*kz)
      + a1*hb4*Ut * (+ 11*O3- kz)
      + a1*M1*ht5 * (- 712*O15)
      + a1*M1*hb2*ht*htau2 * (+ 3.2- 0.8*kz)
      + a1*M1*hb2*ht3 * (- 74*O15- 1.6*kz)
      + a1*M1*hb4*ht * (- 22*O3+ 2*kz)
      + a1*a3*ht2*Ut * (- 166.4+ 17.6*kz)
      + a1*a3*hb*ht*Ub * (- 152*O15+ 128*O45*kz)
      + a1*a3*hb2*Ut * (- 76*O15+ 64*O45*kz)
      + a1*a3*M3*ht3 * (+ 1664*O15- 176*O15*kz)
      + a1*a3*M3*hb2*ht * (+ 152*O15- 128*O45*kz)
      + a1*a3*M1*ht3 * (+ 1664*O15- 176*O15*kz)
      + a1*a3*M1*hb2*ht * (+ 152*O15- 128*O45*kz)
      + a1*a32*Ut * (+ 436*O45- 88*O15*kz)
      + a1*a32*M3*ht * (- 1744*O45+ 352*O15*kz)
      + a1*a32*M1*ht * (- 872*O45+ 176*O15*kz)
      + a1*a2*ht2*Ut * (- 75+ 18.6*kz)
      + a1*a2*hb*ht*Ub * (- 8.2+ 0.6*kz)
      + a1*a2*hb2*Ut * (- 4.1+ 0.3*kz)
      + a1*a2*M2*ht3 * (+ 50- 12.4*kz)
      + a1*a2*M2*hb2*ht * (+ 8.2- 0.6*kz)
      + a1*a2*M1*ht3 * (+ 50- 12.4*kz)
      + a1*a2*M1*hb2*ht * (+ 8.2- 0.6*kz)
      + a1*a2*a3*Ut * (- 1.6)
      + a1*a2*a3*M3*ht * (+ 3.2)
      + a1*a2*a3*M2*ht * (+ 3.2)
      + a1*a2*a3*M1*ht * (+ 3.2)
      + a1*a22*Ut * (+ 8.5- 2.7*kz)
      + a1*a22*M2*ht * (- 34+ 10.8*kz)
      + a1*a22*M1*ht * (- 17+ 5.4*kz)
      + a12*htau2*Ut * (- 7.8)
      + a12*ht*htau*Utau * (- 15.6)
      + a12*ht2*Ut * (- 213.86- 1.82*kz)
      + a12*hb*ht*Ub * (- 3719*O150+ 7*O90*kz)
      + a12*hb2*Ut * (- 3719*O300+ 7*O180*kz)
      + a12*M1*ht*htau2 * (+ 31.2)
      + a12*M1*ht3 * (+ 21386*O75+ 182*O75*kz)
      + a12*M1*hb2*ht * (+ 3719*O75- 7*O45*kz)
      + a12*a3*Ut * (+ 5308*O225- 572*O75*kz)
      + a12*a3*M3*ht * (- 10616*O225+ 1144*O75*kz)
      + a12*a3*M1*ht * (- 21232*O225+ 2288*O75*kz)
      + a12*a2*Ut * (+ 7.58- 2.34*kz)
      + a12*a2*M2*ht * (- 15.16+ 4.68*kz)
      + a12*a2*M1*ht * (- 30.32+ 9.36*kz)
      + a13*Ut * (+ 352097*O6750- 2587*O750*kz)
      + a13*M1*ht * (- 352097*O1125+ 20.696*kz);
    
    dhb3 = 
      + htau6*Ub * (+ 7+ kz)
      + ht6*Ub * (+ 15+ kz)
      + hb*htau5*Utau * (+ 42+ 6*kz)
      + hb*ht5*Ut * (+ 90+ 6*kz)
      + hb2*htau4*Ub * (+ 99)
      + hb2*ht2*htau2*Ub * (- 6)
      + hb2*ht4*Ub * (+ 165)
      + hb3*htau3*Utau * (+ 132)
      + hb3*ht*htau2*Ut * (- 4)
      + hb3*ht2*htau*Utau * (- 4)
      + hb3*ht3*Ut * (+ 220)
      + hb4*htau2*Ub * (+ 20)
      + hb4*ht2*Ub * (+ 60)
      + hb5*htau*Utau * (+ 8)
      + hb5*ht*Ut * (+ 24)
      + hb6*Ub * (+ 714+ 42*kz)
      + a3*ht4*Ub * (+ 40*O3+ 8*kz)
      + a3*hb*ht3*Ut * (+ 160*O3+ 32*kz)
      + a3*hb2*htau2*Ub * (+ 72)
      + a3*hb2*ht2*Ub * (+ 136- 24*kz)
      + a3*hb3*htau*Utau * (+ 48)
      + a3*hb3*ht*Ut * (+ 272*O3- 16*kz)
      + a3*hb4*Ub * (+ 1360*O3)
      + a3*M3*hb*ht4 * (- 80*O3- 16*kz)
      + a3*M3*hb3*htau2 * (- 48)
      + a3*M3*hb3*ht2 * (- 272*O3+ 16*kz)
      + a3*M3*hb5 * (- 544*O3)
      + a32*ht2*Ub * (- 152*O3- 136*O9*kz)
      + a32*hb*ht*Ut * (- 304*O3- 272*O9*kz)
      + a32*hb2*Ub * (- 296- 144*kz)
      + a32*M3*hb*ht2 * (+ 608*O3+ 544*O9*kz)
      + a32*M3*hb3 * (+ 1184*O3+ 192*kz)
      + a33*Ub * (+ 5440*O27+ 320*O3*kz)
      + a33*M3*hb * (- 10880*O9- 640*kz)
      + a2*htau4*Ub * (+ 3+ 3*kz)
      + a2*ht4*Ub * (+ 15+ 3*kz)
      + a2*hb*htau3*Utau * (+ 12+ 12*kz)
      + a2*hb*ht3*Ut * (+ 60+ 12*kz)
      + a2*hb2*htau2*Ub * (+ 45- 9*kz)
      + a2*hb2*ht2*Ub * (+ 81- 9*kz)
      + a2*hb3*htau*Utau * (+ 30- 6*kz)
      + a2*hb3*ht*Ut * (+ 54- 6*kz)
      + a2*hb4*Ub * (+ 300)
      + a2*M2*hb*htau4 * (- 6- 6*kz)
      + a2*M2*hb*ht4 * (- 30- 6*kz)
      + a2*M2*hb3*htau2 * (- 30+ 6*kz)
      + a2*M2*hb3*ht2 * (- 54+ 6*kz)
      + a2*M2*hb5 * (- 120)
      + a2*a3*ht2*Ub * (- 4)
      + a2*a3*hb*ht*Ut * (- 8)
      + a2*a3*hb2*Ub * (- 672+ 120*kz)
      + a2*a3*M3*hb*ht2 * (+ 8)
      + a2*a3*M3*hb3 * (+ 448- 80*kz)
      + a2*a3*M2*hb*ht2 * (+ 8)
      + a2*a3*M2*hb3 * (+ 448- 80*kz)
      + a2*a32*Ub * (+ 68- 24*kz)
      + a2*a32*M3*hb * (- 272+ 96*kz)
      + a2*a32*M2*hb * (- 136+ 48*kz)
      + a22*htau2*Ub * (- 26.25- 5.25*kz)
      + a22*ht2*Ub * (- 56.25- 5.25*kz)
      + a22*hb*htau*Utau * (- 52.5- 10.5*kz)
      + a22*hb*ht*Ut * (- 112.5- 10.5*kz)
      + a22*hb2*Ub * (- 400.5- 67.5*kz)
      + a22*M2*hb*htau2 * (+ 105+ 21*kz)
      + a22*M2*hb*ht2 * (+ 225+ 21*kz)
      + a22*M2*hb3 * (+ 534+ 90*kz)
      + a22*a3*Ub * (+ 140- 36*kz)
      + a22*a3*M3*hb * (- 280+ 72*kz)
      + a22*a3*M2*hb * (- 560+ 144*kz)
      + a23*Ub * (+ 172.5+ 52.5*kz)
      + a23*M2*hb * (- 1035- 315*kz)
      + a1*htau4*Ub * (+ 9- 1.8*kz)
      + a1*ht4*Ub * (+ 17*O3- 1.8*kz)
      + a1*hb*htau3*Utau * (+ 36- 7.2*kz)
      + a1*hb*ht3*Ut * (+ 68*O3- 7.2*kz)
      + a1*hb2*htau2*Ub * (- 13.8+ 4.2*kz)
      + a1*hb2*ht2*Ub * (- 13+ 6*kz)
      + a1*hb3*htau*Utau * (- 9.2+ 2.8*kz)
      + a1*hb3*ht*Ut * (- 26*O3+ 4*kz)
      + a1*hb4*Ub * (+ 200*O3)
      + a1*M1*hb*htau4 * (- 18+ 3.6*kz)
      + a1*M1*hb*ht4 * (- 34*O3+ 3.6*kz)
      + a1*M1*hb3*htau2 * (+ 9.2- 2.8*kz)
      + a1*M1*hb3*ht2 * (+ 26*O3- 4*kz)
      + a1*M1*hb5 * (- 80*O3)
      + a1*a3*ht2*Ub * (- 13.6+ 64*O45*kz)
      + a1*a3*hb*ht*Ut * (- 27.2+ 128*O45*kz)
      + a1*a3*hb2*Ub * (- 86.4+ 20.8*kz)
      + a1*a3*M3*hb*ht2 * (+ 27.2- 128*O45*kz)
      + a1*a3*M3*hb3 * (+ 57.6- 208*O15*kz)
      + a1*a3*M1*hb*ht2 * (+ 27.2- 128*O45*kz)
      + a1*a3*M1*hb3 * (+ 57.6- 208*O15*kz)
      + a1*a32*Ub * (+ 212*O9- 88*O15*kz)
      + a1*a32*M3*hb * (- 848*O9+ 352*O15*kz)
      + a1*a32*M1*hb * (- 424*O9+ 176*O15*kz)
      + a1*a2*htau2*Ub * (- 8.1+ 2.7*kz)
      + a1*a2*ht2*Ub * (- 5.9+ 1.5*kz)
      + a1*a2*hb*htau*Utau * (- 16.2+ 5.4*kz)
      + a1*a2*hb*ht*Ut * (- 11.8+ 3*kz)
      + a1*a2*hb2*Ub * (- 39+ 0.6*kz)
      + a1*a2*M2*hb*htau2 * (+ 16.2- 5.4*kz)
      + a1*a2*M2*hb*ht2 * (+ 11.8- 3*kz)
      + a1*a2*M2*hb3 * (+ 26- 0.4*kz)
      + a1*a2*M1*hb*htau2 * (+ 16.2- 5.4*kz)
      + a1*a2*M1*hb*ht2 * (+ 11.8- 3*kz)
      + a1*a2*M1*hb3 * (+ 26- 0.4*kz)
      + a1*a2*a3*Ub * (- 1.6)
      + a1*a2*a3*M3*hb * (+ 3.2)
      + a1*a2*a3*M2*hb * (+ 3.2)
      + a1*a2*a3*M1*hb * (+ 3.2)
      + a1*a22*Ub * (+ 8.5- 2.7*kz)
      + a1*a22*M2*hb * (- 34+ 10.8*kz)
      + a1*a22*M1*hb * (- 17+ 5.4*kz)
      + a12*htau2*Ub * (- 31.65+ 0.27*kz)
      + a12*ht2*Ub * (- 5587*O300+ 143*O900*kz)
      + a12*hb*htau*Utau * (- 63.3+ 0.54*kz)
      + a12*hb*ht*Ut * (- 5587*O150+ 143*O450*kz)
      + a12*hb2*Ub * (- 99.94- 0.7*kz)
      + a12*M1*hb*htau2 * (+ 126.6- 1.08*kz)
      + a12*M1*hb*ht2 * (+ 5587*O75- 143*O225*kz)
      + a12*M1*hb3 * (+ 9994*O75+ 14*O15*kz)
      + a12*a3*Ub * (+ 3892*O225- 308*O75*kz)
      + a12*a3*M3*hb * (- 7784*O225+ 616*O75*kz)
      + a12*a3*M1*hb * (- 15568*O225+ 1232*O75*kz)
      + a12*a2*Ub * (+ 2.18- 1.26*kz)
      + a12*a2*M2*hb * (- 4.36+ 2.52*kz)
      + a12*a2*M1*hb * (- 8.72+ 5.04*kz)
      + a13*Ub * (+ 194651*O6750- 1393*O750*kz)
      + a13*M1*hb * (- 194651*O1125+ 11.144*kz);
    
    dhl3 = 
      + htau6*Utau * (+ 224+ 28*kz)
      + hb*htau5*Ub * (+ 24)
      + hb*ht2*htau3*Ub * (+ 36)
      + hb*ht4*htau*Ub * (+ 54)
      + hb2*htau4*Utau * (+ 60)
      + hb2*ht*htau3*Ut * (+ 36)
      + hb2*ht2*htau2*Utau * (+ 54)
      + hb2*ht3*htau*Ut * (+ 108)
      + hb2*ht4*Utau * (+ 27)
      + hb3*htau3*Ub * (+ 180)
      + hb4*htau2*Utau * (+ 135)
      + hb5*htau*Ub * (+ 342+ 18*kz)
      + hb6*Utau * (+ 57+ 3*kz)
      + a3*hb*htau3*Ub * (- 192+ 48*kz)
      + a3*hb*ht2*htau*Ub * (+ 48- 16*kz)
      + a3*hb2*htau2*Utau * (- 288+ 72*kz)
      + a3*hb2*ht*htau*Ut * (+ 48- 16*kz)
      + a3*hb2*ht2*Utau * (+ 24- 8*kz)
      + a3*hb3*htau*Ub * (+ 288- 96*kz)
      + a3*hb4*Utau * (+ 72- 24*kz)
      + a3*M3*hb2*htau3 * (+ 192- 48*kz)
      + a3*M3*hb2*ht2*htau * (- 48+ 16*kz)
      + a3*M3*hb4*htau * (- 144+ 48*kz)
      + a32*hb*htau*Ub * (- 320*O3- 16*O3*kz)
      + a32*hb2*Utau * (- 160*O3- 8*O3*kz)
      + a32*M3*hb2*htau * (+ 640*O3+ 32*O3*kz)
      + a2*htau4*Utau * (+ 120)
      + a2*hb*htau3*Ub * (+ 90- 18*kz)
      + a2*hb*ht2*htau*Ub * (+ 36)
      + a2*hb2*htau2*Utau * (+ 135- 27*kz)
      + a2*hb2*ht*htau*Ut * (+ 36)
      + a2*hb2*ht2*Utau * (+ 18)
      + a2*hb3*htau*Ub * (+ 36+ 36*kz)
      + a2*hb4*Utau * (+ 9+ 9*kz)
      + a2*M2*htau5 * (- 48)
      + a2*M2*hb2*htau3 * (- 90+ 18*kz)
      + a2*M2*hb2*ht2*htau * (- 36)
      + a2*M2*hb4*htau * (- 18- 18*kz)
      + a2*a3*hb*htau*Ub * (- 264+ 48*kz)
      + a2*a3*hb2*Utau * (- 132+ 24*kz)
      + a2*a3*M3*hb2*htau * (+ 264- 48*kz)
      + a2*a3*M2*hb2*htau * (+ 264- 48*kz)
      + a22*htau2*Utau * (- 243- 36*kz)
      + a22*ht*htau*Ut * (- 90)
      + a22*ht2*Utau * (- 45)
      + a22*hb*htau*Ub * (- 157.5- 31.5*kz)
      + a22*hb2*Utau * (- 78.75- 15.75*kz)
      + a22*M2*htau3 * (+ 324+ 48*kz)
      + a22*M2*ht2*htau * (+ 180)
      + a22*M2*hb2*htau * (+ 315+ 63*kz)
      + a22*a3*Utau * (+ 180- 36*kz)
      + a22*a3*M3*htau * (- 360+ 72*kz)
      + a22*a3*M2*htau * (- 720+ 144*kz)
      + a23*Utau * (+ 172.5+ 52.5*kz)
      + a23*M2*htau * (- 1035- 315*kz)
      + a1*htau4*Utau * (+ 108)
      + a1*hb*htau3*Ub * (+ 58.8- 1.2*kz)
      + a1*hb*ht2*htau*Ub * (- 4.8+ 2.8*kz)
      + a1*hb2*htau2*Utau * (+ 88.2- 1.8*kz)
      + a1*hb2*ht*htau*Ut * (- 4.8+ 2.8*kz)
      + a1*hb2*ht2*Utau * (- 2.4+ 1.4*kz)
      + a1*hb3*htau*Ub * (+ 12+ 7.2*kz)
      + a1*hb4*Utau * (+ 3+ 1.8*kz)
      + a1*M1*htau5 * (- 43.2)
      + a1*M1*hb2*htau3 * (- 58.8+ 1.2*kz)
      + a1*M1*hb2*ht2*htau * (+ 4.8- 2.8*kz)
      + a1*M1*hb4*htau * (- 6- 3.6*kz)
      + a1*a3*hb*htau*Ub * (- 568*O15+ 112*O15*kz)
      + a1*a3*hb2*Utau * (- 284*O15+ 56*O15*kz)
      + a1*a3*M3*hb2*htau * (+ 568*O15- 112*O15*kz)
      + a1*a3*M1*hb2*htau * (+ 568*O15- 112*O15*kz)
      + a1*a2*htau2*Utau * (- 129.6+ 32.4*kz)
      + a1*a2*hb*htau*Ub * (- 0.6- 3*kz)
      + a1*a2*hb2*Utau * (- 0.3- 1.5*kz)
      + a1*a2*M2*htau3 * (+ 86.4- 21.6*kz)
      + a1*a2*M2*hb2*htau * (+ 0.6+ 3*kz)
      + a1*a2*M1*htau3 * (+ 86.4- 21.6*kz)
      + a1*a2*M1*hb2*htau * (+ 0.6+ 3*kz)
      + a1*a22*Utau * (+ 4.5- 2.7*kz)
      + a1*a22*M2*htau * (- 18+ 10.8*kz)
      + a1*a22*M1*htau * (- 9+ 5.4*kz)
      + a12*htau2*Utau * (- 303.48- 6.48*kz)
      + a12*ht*htau*Ut * (- 46.8)
      + a12*ht2*Utau * (- 23.4)
      + a12*hb*htau*Ub * (- 301*O6- 77*O150*kz)
      + a12*hb2*Utau * (- 301*O12- 77*O300*kz)
      + a12*M1*htau3 * (+ 404.64+ 8.64*kz)
      + a12*M1*ht2*htau * (+ 93.6)
      + a12*M1*hb2*htau * (+ 301*O3+ 77*O75*kz)
      + a12*a3*Utau * (+ 79.2- 15.84*kz)
      + a12*a3*M3*htau * (- 158.4+ 31.68*kz)
      + a12*a3*M1*htau * (- 316.8+ 63.36*kz)
      + a12*a2*Utau * (+ 16.74- 4.86*kz)
      + a12*a2*M2*htau * (- 33.48+ 9.72*kz)
      + a12*a2*M1*htau * (- 66.96+ 19.44*kz)
      + a13*Utau * (+ 99.972- 7.164*kz)
      + a13*M1*htau * (- 599.832+ 42.984*kz);
    
    dM13 =
      a1*htau3*Utau * (- 62.4)
      + a1*ht3*Ut * (- 139.2)
      + a1*hb*htau2*Ub * (- 33.6)
      + a1*hb*ht2*Ub * (- 23.2)
      + a1*hb2*htau*Utau * (- 33.6)
      + a1*hb2*ht*Ut * (- 23.2)
      + a1*hb3*Ub * (- 72)
      + a1*M1*htau4 * (+ 31.2)
      + a1*M1*ht4 * (+ 69.6)
      + a1*M1*hb2*htau2 * (+ 33.6)
      + a1*M1*hb2*ht2 * (+ 23.2)
      + a1*M1*hb4 * (+ 36)
      + a1*a3*ht*Ut * (+ 704*O15)
      + a1*a3*hb*Ub * (+ 512*O15)
      + a1*a3*M3*ht2 * (- 704*O15)
      + a1*a3*M3*hb2 * (- 512*O15)
      + a1*a3*M1*ht2 * (- 704*O15)
      + a1*a3*M1*hb2 * (- 512*O15)
      + a1*a32*M3 * (+ 1936*O15)
      + a1*a32*M1 * (+ 968*O15)
      + a1*a2*htau*Utau * (+ 25.2)
      + a1*a2*ht*Ut * (+ 34.8)
      + a1*a2*hb*Ub * (+ 13.2)
      + a1*a2*M2*htau2 * (- 25.2)
      + a1*a2*M2*ht2 * (- 34.8)
      + a1*a2*M2*hb2 * (- 13.2)
      + a1*a2*M1*htau2 * (- 25.2)
      + a1*a2*M1*ht2 * (- 34.8)
      + a1*a2*M1*hb2 * (- 13.2)
      + a1*a2*a3*M3 * (- 9.6)
      + a1*a2*a3*M2 * (- 9.6)
      + a1*a2*a3*M1 * (- 9.6)
      + a1*a22*M2 * (- 64.8)
      + a1*a22*M1 * (- 32.4)
      + a12*htau*Utau * (+ 6.48)
      + a12*ht*Ut * (+ 338*O75)
      + a12*hb*Ub * (+ 98*O75)
      + a12*M1*htau2 * (- 12.96)
      + a12*M1*ht2 * (- 676*O75)
      + a12*M1*hb2 * (- 196*O75)
      + a12*a3*M3 * (- 2192*O75)
      + a12*a3*M1 * (- 4384*O75)
      + a12*a2*M2 * (- 5.52)
      + a12*a2*M1 * (- 11.04)
      + a13*M1 * (- 513.872);
    
    dM23 =
      + a2*htau3*Utau * (- 40)
      + a2*ht3*Ut * (- 168)
      + a2*hb*htau2*Ub * (- 24)
      + a2*hb*ht2*Ub * (- 24)
      + a2*hb2*htau*Utau * (- 24)
      + a2*hb2*ht*Ut * (- 24)
      + a2*hb3*Ub * (- 168)
      + a2*M2*htau4 * (+ 20)
      + a2*M2*ht4 * (+ 84)
      + a2*M2*hb2*htau2 * (+ 24)
      + a2*M2*hb2*ht2 * (+ 24)
      + a2*M2*hb4 * (+ 84)
      + a2*a3*ht*Ut * (+ 64)
      + a2*a3*hb*Ub * (+ 64)
      + a2*a3*M3*ht2 * (- 64)
      + a2*a3*M3*hb2 * (- 64)
      + a2*a3*M2*ht2 * (- 64)
      + a2*a3*M2*hb2 * (- 64)
      + a2*a32*M3 * (+ 176)
      + a2*a32*M2 * (+ 88)
      + a22*htau*Utau * (+ 22)
      + a22*ht*Ut * (+ 66)
      + a22*hb*Ub * (+ 66)
      + a22*M2*htau2 * (- 44)
      + a22*M2*ht2 * (- 132)
      + a22*M2*hb2 * (- 132)
      + a22*a3*M3 * (+ 48)
      + a22*a3*M2 * (+ 96)
      + a23*M2 * (+ 210)
      + a1*a2*htau*Utau * (+ 8.4)
      + a1*a2*ht*Ut * (+ 11.6)
      + a1*a2*hb*Ub * (+ 4.4)
      + a1*a2*M2*htau2 * (- 8.4)
      + a1*a2*M2*ht2 * (- 11.6)
      + a1*a2*M2*hb2 * (- 4.4)
      + a1*a2*M1*htau2 * (- 8.4)
      + a1*a2*M1*ht2 * (- 11.6)
      + a1*a2*M1*hb2 * (- 4.4)
      + a1*a2*a3*M3 * (- 3.2)
      + a1*a2*a3*M2 * (- 3.2)
      + a1*a2*a3*M1 * (- 3.2)
      + a1*a22*M2 * (+ 7.2)
      + a1*a22*M1 * (+ 3.6)
      + a12*a2*M2 * (- 36.56)
      + a12*a2*M1 * (- 73.12);
    
    dM33 =
      + a3*ht3*Ut * (- 120)
      + a3*hb*htau2*Ub * (- 12)
      + a3*hb*ht2*Ub * (- 16)
      + a3*hb2*htau*Utau * (- 12)
      + a3*hb2*ht*Ut * (- 16)
      + a3*hb3*Ub * (- 120)
      + a3*M3*ht4 * (+ 60)
      + a3*M3*hb2*htau2 * (+ 12)
      + a3*M3*hb2*ht2 * (+ 16)
      + a3*M3*hb4 * (+ 60)
      + a32*ht*Ut * (+ 208*O3)
      + a32*hb*Ub * (+ 208*O3)
      + a32*M3*ht2 * (- 416*O3)
      + a32*M3*hb2 * (- 416*O3)
      + a33*M3 * (+ 694)
      + a2*a3*ht*Ut * (+ 24)
      + a2*a3*hb*Ub * (+ 24)
      + a2*a3*M3*ht2 * (- 24)
      + a2*a3*M3*hb2 * (- 24)
      + a2*a3*M2*ht2 * (- 24)
      + a2*a3*M2*hb2 * (- 24)
      + a2*a32*M3 * (+ 24)
      + a2*a32*M2 * (+ 12)
      + a22*a3*M3 * (- 54)
      + a22*a3*M2 * (- 108)
      + a1*a3*ht*Ut * (+ 88*O15)
      + a1*a3*hb*Ub * (+ 64*O15)
      + a1*a3*M3*ht2 * (- 88*O15)
      + a1*a3*M3*hb2 * (- 64*O15)
      + a1*a3*M1*ht2 * (- 88*O15)
      + a1*a3*M1*hb2 * (- 64*O15)
      + a1*a32*M3 * (+ 88*O15)
      + a1*a32*M1 * (+ 44*O15)
      + a1*a2*a3*M3 * (- 1.2)
      + a1*a2*a3*M2 * (- 1.2)
      + a1*a2*a3*M1 * (- 1.2)
      + a12*a3*M3 * (- 3404*O75)
      + a12*a3*M1 * (- 6808*O75);
    
    // add onto one-loop beta functions
    dmq3 *= threelp; dmq(3,3) += dmq3;
    dmt3 *= threelp; dmu(3,3) += dmt3;
    dmb3 *= threelp; dmd(3,3) += dmb3;
    dml3 *= threelp; dml(3,3) += dml3;
    dmtau3 *= threelp; dme(3,3) += dmtau3;
    dmH1sq3 *= threelp; dmH1sq += dmH1sq3;
    dmH2sq3 *= threelp; dmH2sq += dmH2sq3;
    dm3sq3 *= threelp; dm3sq += dm3sq3;
    dM13 *=threelp; dmG(1) += dM13;
    dM23 *=threelp; dmG(2) += dM23;
    dM33 *=threelp; dmG(3) += dM33;
    dht3 *=threelp; dhu(3,3) += dht3;
    dhb3 *=threelp; dhd(3,3) += dhb3;
    dhl3 *=threelp; dhe(3,3) += dhl3;
    
    // if (getenv("SOFTSUSY_THREELOOP_RGE_DEBUG")!=NULL) {
    //   dout << "dmq3 = " << dmq3 << endl
    //        << "dmt3 = " << dmt3 << endl	 
    //        << "dmb3 = " << dmb3 << endl	 
    //        << "dml3 = " << dml3 << endl	 
    //        << "dmtau3 = " << dmtau3 << endl	 
    //        << "dmH1sq3 = " << dmH1sq3 << endl	 
    //        << "dmH2sq3 = " << dmH2sq3 << endl	 
    //        << "dm3sq3 = " << dm3sq3 << endl	 
    //        << "dM13 = " << dM13 << endl	 
    //        << "dM23 = " << dM23 << endl	 
    //        << "dM33 = " << dM33 << endl	 
    //        << "dht3 = " << dht3 << endl	 
    //        << "dhb3 = " << dhb3 << endl	 
    //        << "dhl3 = " << dhl3 << endl;	 
    // }
    
      } else { // full matrix rge
	
	static DoubleMatrix dmq3(3, 3), 
	  dmu3(3, 3), 
	  dmd3(3, 3), 
	  dme3(3, 3), 
	  dml3(3, 3);
	static DoubleMatrix dhu3(3, 3), dhd3(3, 3), dhe3(3, 3);
	//static DoubleMatrix dmq3(3,3), dmt3(3,3), dmb3(3,3), dmtau3(3,3), dml3(3,3), dht3(3,3), dhb3(3,3), dhl3(3,3); 
	
	// aux matrices
	const DoubleMatrix u6(u4*u2), d6(d4*d2), e4(e2*e2), e6(e4*e2),
	  u6t(u4t*u2t), d6t(d4t*d2t), e4t(e2t*e2t), e6t(e4t*e2t);
	const DoubleMatrix utd(ut*d1), // Yu^+ Yd
	  dtu(dt*u1), // Yd^+ Yu
	  utd2u(utd*dtu), // Yu^+ Yd Yd^+ Yu
	  dtu2d(dtu*utd),  // Yd^+ Yu Yu^+ Yd
	  umuut(u1*mu*ut),  //Yu mu Yu^+
	  dmddt(d1*md*dt),  //Yd md Yd^+
	  emeet(e1*me*et),  //Ye me Ye^+
	  utmqu(ut*mq*u1),  //Yu^+ mq Yu * checked * 
	  dtmqd(dt*mq*d1),  //Yd^+ mq Yd * found BUG <- was d1*mq*d1
	  utmqd(ut*mq*d1),  //Yu^+ mq Yd * checked * 
	  dtmqu(dt*mq*u1),  //Yd^+ mq Yu * found BUG <- was d1*mq*u1
	  etmle(et*me*e1),  //Ye^+ ml Ye 
	  uhut(huut.transpose()), // Yu hu^+
	  uthu(ut*hu),  // Yu^+ hu
	  hutu(uthu.transpose()), // hu^+ Yu
	  dhdt(hddt.transpose()), // Yd hd^+
	  dthd(dt*hd),  // Yd^+ hd
	  hdtd(dthd.transpose()), // hd^+ Yd
	  ehet(heet.transpose()), // Ye he^+
	  ethe(et*he),  // Ye^+ he
	  hete(ethe.transpose()), // he^+ Ye
	  u2td2t(u2t*d2t), // Yu^+ Yu Yd^+ Yd
	  d2tu2t(d2t*u2t), // Yd^+ Yd Yu^+ Yu
	  u2d2(u2*d2), // Yu Yu^+ Yd Yd^+
	  d2u2(d2*u2), // Yd Yd^+ Yu Yu^+
	  ud2tut(u1*d2t*ut), // Yu Yd^+ Yd Yu^+
	  du2tdt(d1*u2t*dt), // Yd Yu^+ Yu Yd^+
	  //?
	  hutu2u(hutu*u2t), // hu^+ Yu Yu^+ Yu
	  uhutd2(uhut*d2), // Yu hu^+ Yd Yd^+	
	  dhdtd2(dhdt*d2), // Yd hd^+ Yd Yd^+	
	  dhdtu2(dhdt*u2), // Yd hd^+ Yu Yu^+	
	  ehete2(ehet*e2), // Ye he^+ Ye Ye^+
	  uhu2tut(uhut*huut),   // Yu hu^+ hu Yu^+
	  hd2u2(hd2*u2), // hd hd^+ Yu Yu^+
	  hu2d2(hu2*d2), // hu hu^+ Yd Yu^+
	  uhuthddt(uhut*hddt), // Yu hu^+ hd Yd^+
	  dhd2tdt(dhdt*hddt), // Yd hd^+ hd Yd^+
	  ehe2tet(ehet*heet), // Ye he^+ he Ye^+
	  u4mq(u4*mq), // Yu Yu^+ Yu Yu^+ mq
	  d4mq(d4*mq), // Yd Yd^+ Yd Yd^+ mq
	  u2d2mq(u2d2*mq), // Yu Yu^+ Yd Yd^+ mq
	  d2u2mq(d2u2*mq), // Yd Yd^+ Yu Yu^+ mq
	  e4ml(e4*ml), // Ye Ye^+ Ye Ye^+ ml 
	  u4tmu(u4t*mu), // Yu^+ Yu Yu^+ Yu mt 
	  d2umuut(d2*umuut), // Yd Yd ^+ Yu mu Yu^+
	  u2dmddt(u2*dmddt), // Yu Yu^+ Yd md Yd^+ 
	  d2dmddt(d2*dmddt), //  Yd Yd^+ Yd md Yd^+ 
	  e2emeet(e2*emeet), // Ye Ye^+ Ye me Ye^+  
	  //u4d2(u4*d2),  // Yu Yu^+ Yu^+ Yu Yd Yd^+
	  //d4u2(d4*u2),  // Yd Yd^+ Yd^+ Yd Yu Yu^+
	  hu2tu4t(hu2t*u4t), //  hu^+ hu Yu^+ Yu Yu^+ Yu 
	  u2dhdthddt(u2*dhdt*hddt), //  Yu Yu^+ Yd hd^+ hd  Yd^+ 
	  hd2d2u2(hd2*d2u2), // hd hd^+ Yd Yd Yu Yu^+ 
	  hu2d2u2(hu2*d2u2), // hu hu^+ Yd Yd Yu Yu^+ 
	  uhuthuutd2(uhut*huut*d2), // Yu hu^+ hu Yu^+ Yd Yd^+	
	  hu2d4(hu2*d4), //  hu hu^+ Yd Yd^+ Yd Yd^+ 
	  hu2u4(hu2*u4), //  hu hu^+ Yu Yu^+ Yu Yu^+ 
	  hd2u4(hd2*u4), //  hd hd^+ Yu Yu^+ Yu Yu^+ 
	  huutuhutd2(huut*uhut*d2), //  hu Yu^+ Yu hu^+ Yd Yd^+
	  hd2d4(hd2*d4), //  hd hd^+ Yd Yd^+ Yd Yd^+
	  hd2td4t(hd2t*d4t), //  hd^+ hd Yd^+ Yd Yd^+ Yd 	
	  uhutu2huut(uhut*u2*huut), //  Yu hu^+ Yu Yu^+ hu Yu^+
	  uhutd2hddt(uhutd2*hddt), //  Yu hu^+ Yd Yd^+ hd Yd^+ 
	  uhutu2hddt(uhut*u2*hddt), //  Yu hu^+ Yu Yu^+ hd Yd^+ 
	  dhdtd2hddt(dhdtd2*hddt), //  Yd hd^+ Yd Yd^+ hd Yd^+ 	
	  he2te4t(he2t*e4t), //  he^+ he Ye^+ Ye Ye^+ Ye  
	  he2e4(he2*e4),  //  he he^+ Ye Ye^+ Ye Ye^+  
	  ehete2heet(ehete2*heet), //  Ye he^+ Ye Ye^+ 
	  d2uhuthddt(d2*uhuthddt), // Yd Yd^+ Yu hu^+ hd Yd^+ 
	  u2uhuthddt(u2*uhuthddt), // Yu Yu^+ Yu hu^+ hd Yd^+  
	dhdtu2hddt(dhdtu2*hddt), // Yd hd^+ Yu Yu^+ hd Yd^+ 
	  uhutd2u2(uhutd2*u2), //  Yu hu^+ Yd Yd^+ Yu Yu^+
	  dhdtd2u2(dhdtd2*u2), // Yd hd^+ Yd Yd^+ Yu Yu^+  
	  dhdtu2d2(dhdtu2*d2), //  Yd hd^+ Yt Yt^+ Yd Yd^+  
	  dhdtu4(dhdt*u4), // Yd hd^+ Yt Yt^+ Yt Yt^+  
	  uhutu2d2(uhut*u2d2), //  Yu hu^+ Yt Yt^+ Yd Yd^+   
	  dhdtd4(dhdt*d4), //  Yd hd^+ Yd Yd^+ Yd Yd^+   
	  uhutu4(uhut*u4), // Yu hu^+ Yu Yu^+ Yu Yu^+    
	  uhutd4(uhut*d4), // Yu hu^+ Yd Yd^+ Yd Yd^+    
	  ehete4(ehet*e4), // Ye he^+ Ye Ye^+ Ye Ye^+     
	  u6mq(u6*mq),  
	  d6mq(d6*mq),  
	  u2d4(u2d2*d2), // Yu Yu^+ Yd^+ Yd Yd^+ Yd 	 
	  u2d4mq(u2d4*mq),   
	  d2u4(d2u2*u2), // Yd Yd^+ Yu^+ Yu Yu^+ Yu 	 
	  d2u4mq(d2u4*mq),   
	  d2u2d2(d2u2*d2), // Yd Yd^+ Yu^+ Yu Yd^+ Yd 	 
	  d2u2d2mq(d2u2d2*mq),
	  d4u2(d4*u2), // Yd Yd^+ Yd^+ Yd Yu^+ Yu 	 
	  d4u2mq(d4u2*mq), // Yd Yd^+ Yd^+ Yd Yu^+ Yu mq 
	  u4d2(u4*d2), // Yu Yu^+ Yu^+ Yu Yd^+ Yd 	 
	  u4d2mq(u4d2*mq), // Yu Yu^+ Yu^+ Yu Yd^+ Yd mq 
	  u2d2u2(u2d2*u2), // Yu Yu^+ Yd^+ Yd Yu^+ Yu 	 
	  u2d2u2mq(u2d2u2*mq), // Yu Yu^+ Yd^+ Yd Yu^+ Yu mq 
	  e6ml(e6*ml), 
	  u6tmu(u6t*mu),
	  u2d2dmddt(u2d2*dmddt), // Yu Yu^+ Yd Yd^+ Yd md Yd^+
	  d2u2dmddt(d2u2*dmddt), // Yd Yd^+ Yu Yu^+ Yd md Yd^+
	  u4dmddt(u4*dmddt), // Yu Yu^+ Yu Yu^+ Yd md Yd^+
	  d6tmd(d6t*md), //  Yd^+ Yd Yd^+ Yd Yd^+ Yd md  
	  u2d2umuut(u2d2*umuut), // Yu Yu^+ Yd Yd^+ Yu md Yu^+
	  d2u2umuut(d2u2*umuut), // Yu Yu^+ Yd Yd^+ Yu md Yu^+
	  d4umuut(d4*umuut), // Yd Yd^+ Yd Yd^+ Yu md Yu^+ 
	  e6tme(e6t*me);  // Ye^+ Ye Ye^+ Ye Ye^+ Ye me  
	
	
	// Traces
	// Old ones
	const double & u2T  = uuT, // Tr ( Yu Yu^+ ) 
	  & d2T  = ddT, // Tr ( Yd Yd^+ ) 
	  & e2T  = eeT, // Tr ( Ye Ye^+ )  
	  & hutuT = huuT, // Tr ( hu^+ Yu ) = Tr( Yu^+ hu ),  Im h = Im Y = 0 assumed
	  & hdtdT = hddT, // Tr ( hd^+ Yd )
	  & heteT = heeT, // Tr ( he^+ Ye )
	  & hu2T = huT, // Tr ( hu^+ hu )
	  & hd2T = hdT, // Tr ( hd^+ hd )
	  & he2T = heT; // Tr ( he^+ he )
	// New ones 
	const double u2mqT = utmqu.trace(), // Tr ( Yu^+ mq Yu ) 
	  d2mqT = dtmqd.trace(), // Tr ( Yd^+ mq Yd )  
	  e2mlT = etmle.trace(), // Tr ( Ye^+ ml Ye )  
	  u2tmuT = umuut.trace(), // Tr ( Yu mu Yu^+ )  
	  d2tmdT = dmddt.trace(), // Tr ( Yd md Yd^+ )  
	  e2tmeT = emeet.trace(), // Tr ( Ye me Ye^+ )  
	  u4T = u4.trace(),
	  d4T = d4.trace(),
	  u2d2T = u2d2.trace(), // Tr( Yu Yu^+ Yd Yd^+ )
	  e4T = e4.trace(),
	  hutu2uT = (hutu2u).trace(),  // Tr ( hu^+ Yu Yu^+ Yu)
	  hutd2uT = (uhutd2).trace(),  // Tr ( Yu hu^+ Yd Yd^+) 
	  hdtd2dT = (dhdtd2).trace(),  // Tr ( hd^+ Yd Yd^+ Yd)  = Tr( Yd hd^+ Yd Yd^+) 
	  hdtu2dT = (dhdtu2).trace(), // Tr ( Yd hd^+ Yu Yu^+ ) 
	  hete2eT = (ehete2).trace(), // Tr ( he^+ Ye Ye^+ Ye) = Tr (Ye he^+ Ye Ye^+)
	  hu2tu2tT = (uhu2tut).trace(),	    // Tr( Yu hu^+ hu Yu^+)
	  hd2u2T = (hd2u2).trace(),   // tr( hd hd^+ Yu Yu^+ )
	  hu2d2T = (hu2d2).trace(),   // tr( hu hu^+ Yd Yd^+ )
	  uhuthddtT = (uhuthddt).trace(), // Tr (Yt hu^+ hd Yd^+)
	  hd2td2tT = (dhd2tdt).trace(), // Tr (Yd hd^+ hd Yd^+ )
	  he2te2tT = (ehe2tet).trace(), // Tr ( Ye he^+ he Ye^+)
	  u4mqT = (u4mq).trace(), // Tr ( Yu Yu^+ Yu Yu^+ mq)
	  d4mqT = (d4mq).trace(), // Tr ( Yd Yd^+ Yd Yd^+ mq)
	  u2d2mqT = (u2d2mq).trace(), // Tr ( Yu Yu^+ Yd Yd^+ mq)
	  d2u2mqT = (d2u2mq).trace(), // Tr ( Yd Yd^+ Yu Yu^+ mq)
	  e2e2mlT = (e4ml).trace(), // Tr (Ye Ye^+ Ye Ye^+ ml )
	  u2umuutT = (u4tmu).trace(), // Tr ( Yu Yu^+ Yu mu Yu^+ ) = Tr (Yu^+ Yu Yu^+ Yu mu) 
	  d2umuutT = (d2umuut).trace(), // Tr ( Yd Yd ^+ Yu mu Yu^+)
	  u2dmddtT = (u2dmddt).trace(), // Tr (Yu Yu^+ Yd md Yd^+ )
	  d2dmddtT = (d2dmddt).trace(), //  Tr (Yd Yd^+ Yd md Yd^+ )
	  e2emeetT = (e2emeet).trace(), // Tr (Ye Ye^+ Ye me Ye^+ ) 
	  ///
	  u6T = u6.trace(), 
	  d6T = d6.trace(), 
	  e6T = e6.trace(), 
	  u4d2T = (u4d2).trace(),  // Tr ( Yu Yu^+ Yu^+ Yu Yd Yd^+)
	  d4u2T = (d4u2).trace(),  // Tr ( Yd Yd^+ Yd^+ Yd Yu Yu^+)
	  hu2tu4tT = (hu2tu4t).trace(), // Tr ( hu^+ hu Yu^+ Yu Yu^+ Yu )
	  hd2dtu2dT = (u2dhdthddt).trace(), // Tr ( Yu Yu^+ Yd hd^+ hd  Yd^+ ) = Tr ( hd^+ hd Yd^+ Yu Yu^+ Yd)
	  hd2d2u2T = (hd2d2u2).trace(), // Tr (Yu Yu^+ hd hd^+ Yd Yd) 
	  hu2d2u2T = (hu2d2u2).trace(), // Tr (Yu Yu^+ hu hu^+ Yd Yd) 
	  hu2tutd2uT = (uhuthuutd2).trace(), // Tr (hu^+ hu Yu^+ Yd Yd^+ Yu) = Tr (Yu hu^+ hu Yu^+ Yd Yd^+)	
	  hu2u4T = (hu2u4).trace(), // Tr ( hu hu^+ Yu Yu^+ Yu Yu^+) 
	  hu2d4T = (hu2d4).trace(), // Tr ( hu hu^+ Yd Yd^+ Yd Yd^+) 
	  hd2u4T = (hd2u4).trace(), // Tr ( hd hd^+ Yu Yu^+ Yu Yu^+) 
	  hutd2huu2tT = (huutuhutd2).trace(), // Tr ( hu Yu^+ Yu hu^+ Yd Yd^+)
	  hd2d4T = (hd2d4).trace(), // Tr ( hd hd^+ Yd Yd^+ Yd Yd^+)
	  hd2td4tT = (hd2td4t).trace(), // Tr ( hd^+ hd Yd^+ Yd Yd^+ Yd) 	
	  hutu2huu2tT = (uhutu2huut).trace(), // Tr ( Yu hu^+ Yu Yu^+ hu Yu^+)
	  hutd2hddtuT = (uhutd2hddt).trace(), // Tr ( Yu hu^+ Yd Yd^+ hd Yd^+) 
	  hutu2hddtuT = (uhutu2hddt).trace(), // Tr ( Yu hu^+ Yu Yu^+ hd Yd^+) 
	  hdtd2hdd2tT = (dhdtd2hddt).trace(), // Tr ( Yd hd^+ Yd Yd^+ hd Yd^+) 	
	  he2te4tT = (he2te4t).trace(), // Tr ( he^+ he Ye^+ Ye Ye^+ Ye ) 
	  he2e4T = (he2e4).trace(), // Tr ( he he^+ Ye Ye^+ Ye Ye^+ ) 
	  hete2hee2tT = (ehete2heet).trace(), // Tr( Ye he^+ Ye Ye^+ he Ye^+ ) 
	  huthdd2tdtuT = (d2uhuthddt).trace(), // Tr (Yd Yd^+ Yu hu^+ hd Yd^+) 
	  huthddtuu2tT = (u2uhuthddt).trace(), // Tr (Yu Yu^+ Yu hu^+ hd Yd^+)  
	  hdtu2hdd2tT  = (dhdtu2hddt).trace(), // Tr (Yd hd^+ Yu Yu^+ hd Yd^+ )
	  hutd2u2uT = (uhutd2u2).trace(), // Tr ( Yu hu^+ Yd Yd^+ Yu Yu^+ )
	  hdtd2u2dT = (dhdtd2u2).trace(), // Tr ( Yd hd^+ Yd Yd^+ Yu Yu^+ ) 
	  hdtu2d2dT = (dhdtu2d2).trace(), // Tr ( Yd hd^+ Yt Yt^+ Yd Yd^+ ) 
	  hdtu4dT = (dhdtu4).trace(), // Tr ( Yd hd^+ Yt Yt^+ Yt Yt^+ ) 
	  hutu2d2uT = (uhutu2d2).trace(), // Tr ( Yu hu^+ Yt Yt^+ Yd Yd^+ )  
	  hdtd4dT = (dhdtd4).trace(), // Tr ( Yd hd^+ Yd Yd^+ Yd Yd^+ )  
	  hutu4uT = (uhutu4).trace(), // Tr (Yu hu^+ Yu Yu^+ Yu Yu^+ )   
	  hutd4uT = (uhutd4).trace(), // Tr (Yu hu^+ Yd Yd^+ Yd Yd^+ )   
	  hete4eT = (ehete4).trace(), // Tr (Ye he^+ Ye Ye^+ Ye Ye^+ )    
	  u6mqT = (u6mq).trace(),  
	  d6mqT = (d6mq).trace(),  
	  u2d4mqT = (u2d4mq).trace(),  
	  d2u4mqT = (d2u4mq).trace(),  
	  d2u2d2mqT = (d2u2d2mq).trace(), 
	  d4u2mqT = (d4u2mq).trace(), 
	  u4d2mqT = (u4d2mq).trace(), 
	  u2d2u2mqT = (u2d2u2mq).trace(), 
	  e6mlT = (e6ml).trace(),
	  u6tmuT = (u6tmu).trace(),
	  dtu2d2dmdT  = (u2d2dmddt).trace(),
	  dtd2u2dmdT  = (d2u2dmddt).trace(), // Tr ( Yd^+ Yd Yd^+ Yu Yu^+ Yd md ) = Tr (Yd Yd^+ Yu Yu^+ Yd md Yd^+)
	  dtu4dmdT  = (u4dmddt).trace(), // Tr ( Yd^+ Yu Yu^+ Yu Yu^+ Yd md ) = Tr (Yu Yu^+ Yu Yu^+ Yd md Yd^+)
	  d6tmdT  = (d6tmd).trace(), // Tr ( Yd^+ Yd Yd^+ Yd Yd^+ Yd md ) 
	  utu2d2umuT = (u2d2umuut).trace(), // Tr ( Yu^+ Yu Yu^+ Yd Yd^+ Yu md ) = Tr (Yu Yu^+ Yd Yd^+ Yu md Yu^+)
	  utd2u2umuT = (d2u2umuut).trace(), // Tr ( Yu^+ Yu Yu^+ Yd Yd^+ Yu md ) = Tr (Yu Yu^+ Yd Yd^+ Yu md Yu^+)
	  utd4umuT = (d4umuut).trace(), // Tr ( Yu^+ Yd Yd^+ Yd Yd^+ Yu md ) = Tr (Yd Yd^+ Yd Yd^+ Yu md Yu^+) 
	  e6tmeT = (e6tme).trace();  //Tr ( Ye^+ Ye Ye^+ Ye Ye^+ Ye me ) 
	
	dM13=-513.872*a13*M1+a12*(1.3066666666666666*hdtdT+6.48*heteT+4.506666666666667*hutuT+a2*(-11.04*M1-5.52*M2)+a3*(-58.45333333333333*M1-29.226666666666667*M3)+M1*(-2.6133333333333333*d2T-12.96*e2T-9.013333333333334*u2T))+a1*(-28.8*d2T*hdtdT-43.2*(hdtd2dT+hete2eT)-19.2*e2T*heteT-33.6*(e2T*hdtdT+d2T*heteT)-23.2*(hdtu2dT+hutd2uT)-67.2*hutu2uT+a22*(-32.4*M1-64.8*M2)+a32*(64.53333333333333*M1+129.06666666666666*M3)+a3*(34.13333333333333*hdtdT+46.93333333333333*hutuT+M1*(-34.13333333333333*d2T-46.93333333333333*u2T)+M3*(-34.13333333333333*d2T-46.93333333333333*u2T))+a2*(13.2*hdtdT+25.2*heteT+34.8*hutuT+a3*(-9.6*M1-9.6*M2-9.6*M3)+M1*(-13.2*d2T-25.2*e2T-34.8*u2T)+M2*(-13.2*d2T-25.2*e2T-34.8*u2T))-72.*hutuT*u2T+M1*(14.4*d2T*d2T+9.6*e2T*e2T+21.6*(d4T+e4T)+23.2*u2d2T+36.*u2T*u2T+33.6*(d2T*e2T+u4T))); 
	
	dM23=a12*a2*(-73.12*M1-36.56*M2)+210.*a23*M2+a1*(a22*(3.6*M1+7.2*M2)+a2*(4.4*hdtdT+8.4*heteT+11.6*hutuT+a3*(-3.2*M1-3.2*M2-3.2*M3)+M1*(-4.4*d2T-8.4*e2T-11.6*u2T)+M2*(-4.4*d2T-8.4*e2T-11.6*u2T)))+a22*(22.*heteT+66.*(hdtdT+hutuT)+a3*(96.*M2+48.*M3)+M2*(-44.*e2T-132.*(d2T+u2T)))+a2*(-32.*hete2eT-e2T*(24.*hdtdT+8.*heteT)-24.*(hdtu2dT+d2T*heteT+hutd2uT)-96.*(hdtd2dT+hutu2uT)+a32*(88.*M2+176.*M3)-72.*(d2T*hdtdT+hutuT*u2T)+a3*(64.*(hdtdT+hutuT)-64.*M2*(d2T+u2T)-64.*M3*(d2T+u2T))+M2*(24.*d2T*e2T+4.*e2T*e2T+16.*e4T+24.*u2d2T+36.*(d2T*d2T+u2T*u2T)+48.*(d4T+u4T)));
	
	dM33=a22*a3*(-108.*M2-54.*M3)+a12*a3*(-90.77333333333333*M1-45.38666666666666*M3)+694.*a33*M3+a1*(a2*a3*(-1.2*M1-1.2*M2-1.2*M3)+a32*(2.933333333333333*M1+5.866666666666666*M3)+a3*(4.266666666666667*hdtdT+5.866666666666666*hutuT+M1*(-4.266666666666667*d2T-5.866666666666666*u2T)+M3*(-4.266666666666667*d2T-5.866666666666666*u2T)))+a32*(69.33333333333333*(hdtdT+hutuT)-138.66666666666666*M3*(d2T+u2T))+a2*(a32*(12.*M2+24.*M3)+a3*(24.*(hdtdT+hutuT)-24.*M2*(d2T+u2T)-24.*M3*(d2T+u2T)))+a3*(-12.*(e2T*hdtdT+d2T*heteT)-16.*(hdtu2dT+hutd2uT)-48.*(hdtd2dT+hutu2uT)-72.*(d2T*hdtdT+hutuT*u2T)+M3*(12.*d2T*e2T+16.*u2d2T+36.*(d2T*d2T+u2T*u2T)+24.*(d4T+u4T)));
	
	dmH1sq3=216.*d2dmddtT*d2T+216.*d2T*d4mqT+72.*d2T*e2e2mlT+72.*d2T*e2emeetT+36.*d4T*e2mlT+72.*d2dmddtT*e2T+72.*d4mqT*e2T+24.*e2e2mlT*e2T+24.*e2emeetT*e2T+36.*d4T*e2tmeT+36.*d2mqT*e4T+36.*d2tmdT*e4T+12.*e2mlT*e4T+12.*e2tmeT*e4T+36.*e4T*hd2T+108.*d4T*(d2mqT+d2tmdT+hd2T)+144.*e2T*hd2td2tT+432.*(d2T*hd2td2tT+hdtd2dT*hdtdT)+36.*d4T*he2T+12.*e4T*he2T+144.*d2T*he2te2tT+48.*e2T*he2te2tT+144.*hdtdT*hete2eT+144.*hdtd2dT*heteT+48.*hete2eT*heteT+36.*hu2d2u2T+18.*(d2u4mqT+d6mqT+d6tmdT+dtu4dmdT+hd2d4T+hd2td4tT+hd2u4T+hdtd2hdd2tT+hu2tutd2uT+hutd2huu2tT)+36.*huthddtuu2tT+36.*hutu2hddtuT+72.*hdtu2dT*hutuT+72.*hutd2uT*hutuT+6.*(e6mlT+e6tmeT+he2e4T+he2te4tT+hete2hee2tT)*(1+kz)+a23*((2157+630.*kz)*M22-3.*mH1sq-3.*mH2sq-3.*mlT-9.*mqT)+a13*((55767./125-28.656*kz)*M12-3.312*mdT-9.936*meT-4.968*mH1sq-4.968*mH2sq-4.968*mlT-1.656*mqT-13.248*muT)+a32*(-((d2mqT+d2tmdT+hd2T)*(106.66666666666667+5.333333333333333*kz))+hdtdT*(426.6666666666667+21.333333333333332*kz)*M3-d2T*(448.+32.*kz)*M32-d2T*(106.66666666666667+5.333333333333333*kz)*mH1sq-64.*d2T*mqT-32.*d2T*(mdT+muT))+36.*hu2T*u2d2T+36.*u2d2T*u2mqT+36.*d2u2mqT*u2T+36.*d2umuutT*u2T+36.*hd2u2T*u2T+36.*hu2d2T*u2T+36.*u2d2mqT*u2T+36.*u2dmddtT*u2T+36.*u2d2T*u2tmuT+a22*(-121.5*(d2mqT+d2tmdT+hd2T)-40.5*he2T-31.5*d2mqT*kz-31.5*d2tmdT*kz-31.5*hd2T*kz-10.5*he2T*kz-e2mlT*(40.5+10.5*kz)-e2tmeT*(40.5+10.5*kz)+(180.*hutuT+heteT*(150.+42.*kz)+hdtdT*(450.+126.*kz))*M2+a3*((1080-216.*kz)*M22+(720-144.*kz)*M2*M3+(432-72.*kz)*M32)+(-(e2T*(40.5+10.5*kz))-d2T*(121.5+31.5*kz))*mH1sq+M22*(-(e2T*(225.+63.*kz))-d2T*(675.+189.*kz)-270.*u2T)-54.*mH2sq*u2T-54.*(hu2T+u2mqT+u2tmuT))+a12*(-30.006666666666668*hd2T-61.38*(e2mlT+e2tmeT+he2T)-d2mqT*(30.006666666666668+0.5133333333333333*kz)-d2tmdT*(30.006666666666668+0.5133333333333333*kz)+0.54*e2mlT*kz+0.54*e2tmeT*kz-0.5133333333333333*hd2T*kz+0.54*he2T*kz+(31.2*hutuT+heteT*(241.2-2.16*kz)+hdtdT*(116.66666666666667+2.0533333333333332*kz))*M1+a3*((792./5-31.68*kz)*M12+(528./5-21.12*kz)*M1*M3+(1584./25-10.56*kz)*M32)+(-(d2T*(29.526666666666667+0.5133333333333333*kz))+e2T*(-62.82+0.54*kz))*mH1sq+a2*((81./25-9.72*kz)*M12+(54./25-6.48*kz)*M1*M2+(108./25-3.24*kz)*M22-0.72*mdT-2.16*meT-1.08*mH1sq-1.08*mH2sq-1.08*mlT-0.36*mqT-2.88*muT)+d2T*(0.32*mdT+0.96*meT+0.48*mlT+0.16*mqT+1.28*muT)-e2T*(0.96*mdT+2.88*meT+1.44*mlT+0.48*mqT+3.84*muT)+M12*(-361.8*e2T+3.24*e2T*kz-d2T*(175.+3.08*kz)-46.8*u2T)+mH2sq*(0.48*d2T-1.44*e2T-9.36*u2T)-9.36*(hu2T+u2mqT+u2tmuT))+mH2sq*(72.*u2d2T*u2T+36.*u4d2T)+mH1sq*(324.*d2T*d4T+108.*d4T*e2T+108.*d2T*e4T+36.*e2T*e4T+6.*e6T*(1+kz)+36.*u2d2T*u2T+18.*(d6T+d6T*kz+u4d2T))+72.*u2T*uhuthddtT+a3*(288.*d2dmddtT+288.*d4mqT+576.*hd2td2tT+48.*hd2u2T+48.*hu2d2T+d2u2mqT*(48.-16.*kz)+d2umuutT*(48.-16.*kz)-96.*d2dmddtT*kz-96.*d4mqT*kz-192.*hd2td2tT*kz-16.*hd2u2T*kz-16.*hu2d2T*kz+(-576.*hdtd2dT-96.*hdtu2dT-96.*hutd2uT+192.*hdtd2dT*kz+32.*hdtu2dT*kz+32.*hutd2uT*kz)*M3+48.*u2d2mqT-16.*kz*u2d2mqT+(48.-16.*kz)*mH2sq*u2d2T+mH1sq*(288.*d4T+48.*u2d2T-kz*(96.*d4T+16.*u2d2T))+M32*(288.*d4T+96.*u2d2T-kz*(96.*d4T+32.*u2d2T))+48.*u2dmddtT-16.*kz*u2dmddtT+96.*uhuthddtT-32.*kz*uhuthddtT)+a1*(12.*d2dmddtT-4.8*d2u2mqT-4.8*d2umuutT+12.*d4mqT+36.*e2e2mlT+36.*e2emeetT+24.*hd2td2tT-4.8*hd2u2T+72.*he2te2tT-4.8*hu2d2T+7.2*d2dmddtT*kz+2.8*d2u2mqT*kz+2.8*d2umuutT*kz+7.2*d4mqT*kz-7.2*e2e2mlT*kz-7.2*e2emeetT*kz+14.4*hd2td2tT*kz+2.8*hd2u2T*kz-14.4*he2te2tT*kz+2.8*hu2d2T*kz+(-24.*hdtd2dT-72.*hete2eT+hdtu2dT*(9.6-5.6*kz)+hutd2uT*(9.6-5.6*kz)-14.4*hdtd2dT*kz+14.4*hete2eT*kz)*M1+a3*((d2mqT+d2tmdT+hd2T)*(-37.86666666666667+7.466666666666667*kz)+d2T*(-75.73333333333333+14.933333333333334*kz)*M12+hdtdT*(75.73333333333333-14.933333333333334*kz)*M3+d2T*(-75.73333333333333+14.933333333333334*kz)*M32+M1*(hdtdT*(75.73333333333333-14.933333333333334*kz)+d2T*(-75.73333333333333+14.933333333333334*kz)*M3)+d2T*(-37.86666666666667+7.466666666666667*kz)*mH1sq)+a2*(-0.6*hd2T-16.2*(e2mlT+e2tmeT+he2T)+5.4*e2mlT*kz+5.4*e2tmeT*kz-3.*hd2T*kz+5.4*he2T*kz-d2mqT*(0.6+3.*kz)-d2tmdT*(0.6+3.*kz)+(-(d2T*(1.2+6.*kz))+e2T*(-32.4+10.8*kz))*M12+(heteT*(32.4-10.8*kz)+hdtdT*(1.2+6.*kz))*M2+(-(d2T*(1.2+6.*kz))+e2T*(-32.4+10.8*kz))*M22+M1*(heteT*(32.4-10.8*kz)+hdtdT*(1.2+6.*kz)+(-(d2T*(1.2+6.*kz))+e2T*(-32.4+10.8*kz))*M2)+(-(d2T*(0.6+3.*kz))+e2T*(-16.2+5.4*kz))*mH1sq)+a22*((72./5-5.4*kz)*M12+(18-10.8*kz)*M1*M2+(171./5-16.2*kz)*M22-1.8*mH1sq-1.8*mH2sq-1.8*mlT-5.4*mqT)-4.8*u2d2mqT+2.8*kz*u2d2mqT+(-4.8+2.8*kz)*mH2sq*u2d2T+mH1sq*(12.*d4T+36.*e4T+7.2*(d4T-e4T)*kz-4.8*u2d2T+2.8*kz*u2d2T)+M12*(12.*d4T+36.*e4T+7.2*(d4T-e4T)*kz-9.6*u2d2T+5.6*kz*u2d2T)-4.8*u2dmddtT+2.8*kz*u2dmddtT-9.6*uhuthddtT+5.6*kz*uhuthddtT)+a2*(72.*hd2td2tT+24.*he2te2tT+36.*(d2dmddtT+d2u2mqT+d2umuutT+d4mqT+hd2u2T+hu2d2T)+36.*d2dmddtT*kz+36.*d4mqT*kz+24.*he2te2tT*kz+12.*e2e2mlT*(1+kz)+12.*e2emeetT*(1+kz)+(-24.*hete2eT*(1+kz)-72.*(hdtd2dT+hdtu2dT+hutd2uT+hdtd2dT*kz))*M2+a3*((d2mqT+d2tmdT+hd2T)*(-264.+48.*kz)+d2T*(-528.+96.*kz)*M22+hdtdT*(528.-96.*kz)*M3+d2T*(-528.+96.*kz)*M32+M2*(hdtdT*(528.-96.*kz)+d2T*(-528.+96.*kz)*M3)+d2T*(-264.+48.*kz)*mH1sq)+36.*u2d2mqT+36.*mH2sq*u2d2T+M22*((36.*d4T+12.*e4T)*(1+kz)+72.*u2d2T)+mH1sq*(12.*e4T*(1+kz)+36.*(d4T+d4T*kz+u2d2T))+36.*u2dmddtT+72.*(hd2td2tT*kz+uhuthddtT))+18.*((d6mqT+d6tmdT+hd2d4T+hd2td4tT+hdtd2hdd2tT)*kz+u2d2u2mqT+u4d2mqT+utd2u2umuT+utu2d2umuT);
	
	dmH2sq3=36.*d2T*d2u2mqT+36.*d2T*d2umuutT+12.*d2u2mqT*e2T+12.*d2umuutT*e2T+36.*hd2d2u2T+36.*d2T*hd2u2T+12.*e2T*hd2u2T+36.*d2T*hu2d2T+12.*e2T*hu2d2T+36.*hutd2hddtuT+(72.*hdtdT+24.*heteT)*(hdtu2dT+hutd2uT)+36.*huthdd2tdtuT+a23*((2157+630.*kz)*M22-3.*mH1sq-3.*mH2sq-3.*mlT-9.*mqT)+a13*((55767./125-28.656*kz)*M12-3.312*mdT-9.936*meT-4.968*mH1sq-4.968*mH2sq-4.968*mlT-1.656*mqT-13.248*muT)+36.*d2T*u2d2mqT+12.*e2T*u2d2mqT+36.*d2mqT*u2d2T+36.*d2tmdT*u2d2T+12.*e2mlT*u2d2T+12.*e2tmeT*u2d2T+36.*hd2T*u2d2T+12.*he2T*u2d2T+mH1sq*(36.*d4u2T+72.*d2T*u2d2T+24.*e2T*u2d2T)+36.*d2T*u2dmddtT+12.*e2T*u2dmddtT+432.*(hutu2uT*hutuT+hu2tu2tT*u2T)+a12*(-5.04*(d2mqT+d2tmdT+hd2T)-6.48*(e2mlT+e2tmeT+he2T)+(16.8*hdtdT+21.6*heteT+hutuT*(283.06666666666666+1.7333333333333334*kz))*M1+a3*((792./5-31.68*kz)*M12+(528./5-21.12*kz)*M1*M3+(1584./25-10.56*kz)*M32)+a2*((81./25-9.72*kz)*M12+(54./25-6.48*kz)*M1*M2+(108./25-3.24*kz)*M22-0.72*mdT-2.16*meT-1.08*mH1sq-1.08*mH2sq-1.08*mlT-0.36*mqT-2.88*muT)+mH1sq*(-5.04*d2T-6.48*e2T-0.96*u2T)-(73.28666666666666+0.43333333333333335*kz)*mH2sq*u2T-(0.64*mdT+1.92*meT+0.96*mlT+0.32*mqT+2.56*muT)*u2T+M12*(-25.2*d2T-32.4*e2T-(424.6+2.6*kz)*u2T)-(72.32666666666667+0.43333333333333335*kz)*(hu2T+u2mqT+u2tmuT))+a32*(hutuT*(426.6666666666667+21.333333333333332*kz)*M3-(448.+32.*kz)*M32*u2T-(106.66666666666667+5.333333333333333*kz)*mH2sq*u2T-64.*mqT*u2T-32.*(mdT+muT)*u2T-(106.66666666666667+5.333333333333333*kz)*(hu2T+u2mqT+u2tmuT))+a22*(-54.*(d2mqT+d2tmdT+hd2T)-18.*(e2mlT+e2tmeT+he2T)+(180.*hdtdT+60.*heteT+hutuT*(450.+126.*kz))*M2+a3*((1080-216.*kz)*M22+(720-144.*kz)*M2*M3+(432-72.*kz)*M32)+(-54.*d2T-18.*e2T)*mH1sq-(121.5+31.5*kz)*mH2sq*u2T+M22*(-270.*d2T-90.*e2T-(675.+189.*kz)*u2T)-(121.5+31.5*kz)*(hu2T+u2mqT+u2tmuT))+216.*u2T*(u2umuutT+u4mqT)+108.*(hu2T+u2mqT+u2tmuT)*u4T+mH2sq*(36.*d2T*u2d2T+12.*e2T*u2d2T+324.*u2T*u4T+18.*(d4u2T+u6T+kz*u6T))+72.*d2T*uhuthddtT+24.*e2T*uhuthddtT+a3*(48.*hd2u2T+48.*hu2d2T+576.*hu2tu2tT+d2u2mqT*(48.-16.*kz)+d2umuutT*(48.-16.*kz)-16.*hd2u2T*kz-16.*hu2d2T*kz-192.*hu2tu2tT*kz+(-96.*hdtu2dT-96.*hutd2uT-576.*hutu2uT+32.*hdtu2dT*kz+32.*hutd2uT*kz+192.*hutu2uT*kz)*M3+48.*u2d2mqT-16.*kz*u2d2mqT+(48.-16.*kz)*mH1sq*u2d2T+48.*u2dmddtT-16.*kz*u2dmddtT+288.*u2umuutT-96.*kz*u2umuutT+288.*u4mqT-96.*kz*u4mqT+mH2sq*(48.*u2d2T+288.*u4T-kz*(16.*u2d2T+96.*u4T))+M32*(96.*u2d2T+288.*u4T-kz*(32.*u2d2T+96.*u4T))+96.*uhuthddtT-32.*kz*uhuthddtT)+a1*(2.4*hd2u2T+2.4*hu2d2T+91.2*hu2tu2tT+d2u2mqT*(2.4+0.4*kz)+d2umuutT*(2.4+0.4*kz)+0.4*hd2u2T*kz+0.4*hu2d2T*kz-4.8*hu2tu2tT*kz+(-(hdtu2dT*(4.8+0.8*kz))-hutd2uT*(4.8+0.8*kz)+hutu2uT*(-91.2+4.8*kz))*M1+a22*((72./5-5.4*kz)*M12+(18-10.8*kz)*M1*M2+(171./5-16.2*kz)*M22-1.8*mH1sq-1.8*mH2sq-1.8*mlT-5.4*mqT)+2.4*u2d2mqT+0.4*kz*u2d2mqT+(2.4+0.4*kz)*mH1sq*u2d2T+2.4*u2dmddtT+0.4*kz*u2dmddtT+a2*(hutuT*(22.8-8.4*kz)*M2+(-22.8+8.4*kz)*M12*u2T+(-22.8+8.4*kz)*M22*u2T+(-11.4+4.2*kz)*mH2sq*u2T+M1*(hutuT*(22.8-8.4*kz)+(-22.8+8.4*kz)*M2*u2T)+(-11.4+4.2*kz)*(hu2T+u2mqT+u2tmuT))+a3*(hutuT*(165.33333333333334-27.733333333333334*kz)*M3+(-165.33333333333334+27.733333333333334*kz)*M12*u2T+(-165.33333333333334+27.733333333333334*kz)*M32*u2T+(-82.66666666666667+13.866666666666667*kz)*mH2sq*u2T+M1*(hutuT*(165.33333333333334-27.733333333333334*kz)+(-165.33333333333334+27.733333333333334*kz)*M3*u2T)+(-82.66666666666667+13.866666666666667*kz)*(hu2T+u2mqT+u2tmuT))+45.6*u2umuutT-2.4*kz*u2umuutT+45.6*u4mqT-2.4*kz*u4mqT+mH2sq*(2.4*u2d2T+0.4*kz*u2d2T+45.6*u4T-2.4*kz*u4T)+M12*(4.8*u2d2T+0.8*kz*u2d2T+45.6*u4T-2.4*kz*u4T)+4.8*uhuthddtT+0.8*kz*uhuthddtT)+a2*(-72.*(hdtu2dT+hutd2uT+hutu2uT+hutu2uT*kz)*M2+36.*mH1sq*u2d2T+a3*(hutuT*(528.-96.*kz)*M3+(-528.+96.*kz)*M22*u2T+(-528.+96.*kz)*M32*u2T+(-264.+48.*kz)*mH2sq*u2T+M2*(hutuT*(528.-96.*kz)+(-528.+96.*kz)*M3*u2T)+(-264.+48.*kz)*(hu2T+u2mqT+u2tmuT))+36.*(d2u2mqT+d2umuutT+hd2u2T+hu2d2T+u2d2mqT+u2dmddtT+(1+kz)*(u2umuutT+u4mqT))+36.*mH2sq*(u2d2T+u4T+kz*u4T)+M22*(72.*u2d2T+36.*(1+kz)*u4T)+72.*(hu2tu2tT+hu2tu2tT*kz+uhuthddtT))+18.*(d2u2d2mqT+d4u2mqT+dtd2u2dmdT+dtu2d2dmdT+hd2dtu2dT+hdtu2hdd2tT+hu2d4T+hu2u4T+hutu2huu2tT+hu2tu4tT*(1+kz)+u2d4mqT+u6mqT+u6tmuT+kz*(hu2u4T+hutu2huu2tT+u6mqT+u6tmuT)+utd4umuT);
	
	dhu3=(9.*d4u2T+a13*(52.16251851851852-3.449333333333333*kz)+a23*(172.5+52.5*kz)+a33*(201.4814814814815+106.66666666666667*kz)+18.*d2T*u2d2T+6.*e2T*u2d2T+a22*(-45.*d2T-15.*e2T+a3*(140.-36.*kz)-78.75*u2T-15.75*kz*u2T)+a32*(-53.333333333333336*d2T-106.66666666666667*u2T-2.6666666666666665*kz*u2T)+a12*(-6.066666666666666*d2T-7.8*e2T+a3*(23.59111111111111-7.626666666666667*kz)+a2*(7.58-2.34*kz)-42.75*u2T-0.21666666666666667*kz*u2T)+54.*u2T*u4T+a1*(a32*(9.688888888888888-5.866666666666666*kz)+a22*(8.5-2.7*kz)+1.2*u2d2T+0.2*kz*u2d2T+a3*(-41.333333333333336+6.933333333333334*kz)*u2T+a2*(-1.6*a3+(-5.7+2.1*kz)*u2T)+11.4*u4T-0.6*kz*u4T)+a2*(a32*(68.-24.*kz)+18.*u2d2T+a3*(-132.+24.*kz)*u2T+9.*(1+kz)*u4T)+a3*(24.*u2d2T+72.*u4T-kz*(8.*u2d2T+24.*u4T))+3.*(1+kz)*u6T)*hu+(18.*(hdtd2u2dT+hdtu2d2dT+hutd4uT+hutu4uT+hutu4uT*kz)+a13*(-312.9751111111111+20.696*kz)*M1+a23*(-1035.-315.*kz)*M2+a33*(-1208.888888888889-640.*kz)*M3+36.*(d2T*(hdtu2dT+hutd2uT)+hdtdT*u2d2T)+12.*(e2T*(hdtu2dT+hutd2uT)+heteT*u2d2T)+216.*hutu2uT*u2T+a12*(-12.133333333333333*hdtdT-15.6*heteT-85.5*hutuT-0.43333333333333335*hutuT*kz+a2*((-30.32+9.36*kz)*M1+(-15.16+4.68*kz)*M2)+a3*((-94.36444444444444+30.506666666666668*kz)*M1+(-47.18222222222222+15.253333333333334*kz)*M3)+M1*(24.266666666666666*d2T+31.2*e2T+171.*u2T+0.8666666666666667*kz*u2T))+a32*(-106.66666666666667*hdtdT-213.33333333333334*hutuT-5.333333333333333*hutuT*kz+M3*(213.33333333333334*d2T+426.6666666666667*u2T+10.666666666666666*kz*u2T))+a22*(-90.*hdtdT-30.*heteT-157.5*hutuT-31.5*hutuT*kz+a3*((-560.+144.*kz)*M2+(-280.+72.*kz)*M3)+M2*(180.*d2T+60.*e2T+315.*u2T+63.*kz*u2T))+108.*hutuT*u4T+a1*(hutu2uT*(45.6-2.4*kz)+hdtu2dT*(2.4+0.4*kz)+hutd2uT*(2.4+0.4*kz)+a22*((-17.+5.4*kz)*M1+(-34.+10.8*kz)*M2)+a32*((-19.377777777777776+11.733333333333333*kz)*M1+(-38.75555555555555+23.466666666666665*kz)*M3)+a2*(hutuT*(-11.4+4.2*kz)+a3*(3.2*M1+3.2*M2+3.2*M3)+(11.4-4.2*kz)*M1*u2T+(11.4-4.2*kz)*M2*u2T)+a3*(hutuT*(-82.66666666666667+13.866666666666667*kz)+(82.66666666666667-13.866666666666667*kz)*M1*u2T+(82.66666666666667-13.866666666666667*kz)*M3*u2T)+M1*(-2.4*u2d2T-0.4*kz*u2d2T-22.8*u4T+1.2*kz*u4T))+a3*(hutu2uT*(288.-96.*kz)+hdtu2dT*(48.-16.*kz)+hutd2uT*(48.-16.*kz)+M3*(-48.*u2d2T+16.*kz*u2d2T-144.*u4T+48.*kz*u4T))+a2*(36.*(hdtu2dT+hutd2uT+hutu2uT+hutu2uT*kz)+a32*((-136.+48.*kz)*M2+(-272.+96.*kz)*M3)+a3*(hutuT*(-264.+48.*kz)+(264.-48.*kz)*M2*u2T+(264.-48.*kz)*M3*u2T)+M2*(-36.*u2d2T-18.*(1+kz)*u4T)))*u1+(-18.*d2T*d2T+36.*d4T-12.*d2T*e2T-2.*e2T*e2T+a2*(-8.*a3+36.*d2T+12.*e2T)+16.*a3*(e2T+d2T*(-1+kz))+a32*(5.333333333333333-30.22222222222222*kz)+a22*(-22.5-10.5*kz)+a12*(-12.66+0.07777777777777778*kz)+a1*(6.4*d2T-3.2*e2T+a2*(-8.2+0.6*kz)-1.6*d2T*kz+0.8*e2T*kz+a3*(-10.133333333333333+2.8444444444444446*kz))+12.*(e4T+u2d2T))*hddt*u1+(a32*(10.666666666666666-60.44444444444444*kz)+a22*(-66.-12.*kz)+a12*(-41.093333333333334-0.2311111111111111*kz)+24.*u2d2T+32.*a3*(-1+kz)*u2T-36.*u2T*u2T+a2*(a3*(-96.+16.*kz)+(63.-9.*kz)*u2T)+a1*(a3*(-27.733333333333334+0.35555555555555557*kz)+a2*(-25.2+5.6*kz)+(11.-0.2*kz)*u2T)+72.*u4T)*hu*u2t+(-9.*d2T*d2T+18.*d4T-6.*d2T*e2T-1.*e2T*e2T+a2*(-4.*a3+18.*d2T+6.*e2T)+8.*a3*(e2T+d2T*(-1+kz))+a32*(2.6666666666666665-15.11111111111111*kz)+a22*(-11.25-5.25*kz)+a12*(-6.33+0.03888888888888889*kz)+a1*(3.2*d2T-1.6*e2T+a2*(-4.1+0.3*kz)-0.8*d2T*kz+0.4*e2T*kz+a3*(-5.066666666666666+1.4222222222222223*kz))+6.*(e4T+u2d2T))*d2*hu+(72.*hdtd2dT-36.*d2T*hdtdT+24.*hete2eT-e2T*(12.*hdtdT+4.*heteT)+12.*(hdtu2dT-d2T*heteT+hutd2uT)+a12*(25.32-0.15555555555555556*kz)*M1+a22*(45.+21.*kz)*M2+a32*(-10.666666666666666+60.44444444444444*kz)*M3+a3*(16.*(heteT+hdtdT*(-1+kz))-16.*(e2T+d2T*(-1+kz))*M3)+a2*(36.*hdtdT+12.*heteT+(-36.*d2T-12.*e2T)*M2+a3*(8.*M2+8.*M3))+a1*(6.4*hdtdT-3.2*heteT-1.6*hdtdT*kz+0.8*heteT*kz+(-6.4*d2T+3.2*e2T+1.6*d2T*kz-0.8*e2T*kz)*M1+a2*((8.2-0.6*kz)*M1+(8.2-0.6*kz)*M2)+a3*((10.133333333333333-2.8444444444444446*kz)*M1+(10.133333333333333-2.8444444444444446*kz)*M3)))*d2*u1+(a32*(13.333333333333334-75.55555555555556*kz)+a22*(-98.25-8.25*kz)+a12*(-44.516666666666666-0.9388888888888889*kz)+30.*u2d2T+40.*a3*(-1+kz)*u2T-45.*u2T*u2T+a2*(a3*(-180.+32.*kz)+(72.-18.*kz)*u2T)+a1*(a3*(-14.666666666666666-3.5555555555555554*kz)+a2*(-32.7+6.7*kz)+(16.+2.*kz)*u2T)+90.*u4T)*u2*hu+(36.*(hdtu2dT+hutd2uT)+216.*hutu2uT+a12*(114.14666666666666+1.56*kz)*M1+a22*(219.+27.*kz)*M2+a32*(-32.+181.33333333333334*kz)*M3-108.*hutuT*u2T+a1*(hutuT*(18.+1.2*kz)+a2*((38.6-8.2*kz)*M1+(38.6-8.2*kz)*M2)+a3*((28.266666666666666+2.1333333333333333*kz)*M1+(28.266666666666666+2.1333333333333333*kz)*M3)-(18.+1.2*kz)*M1*u2T)+a2*(hutuT*(90.-18.*kz)+a3*((184.-32.*kz)*M2+(184.-32.*kz)*M3)+(-90.+18.*kz)*M2*u2T)+a3*(48.*hutuT*(-1+kz)-48.*(-1+kz)*M3*u2T))*u1*u2t+(42.666666666666664*a3+12.*d2T+4.*e2T+a1*(0.9333333333333333-0.4*kz)+6.*a2*(-1+kz))*hddt*d1*dtu+(21.333333333333332*a3+12.*d2T+4.*e2T+a2*(9.-3.*kz)+a1*(1.2666666666666666+0.6*kz)-6.*u2T)*hu*utd2u+(64.*a3+a1*(7.-1.*kz)+3.*a2*(1+kz)+18.*u2T)*hu*u4t+(42.666666666666664*a3+12.*d2T+4.*e2T+a1*(0.9333333333333333-0.4*kz)+6.*a2*(-1+kz))*d2*hddt*u1+(21.333333333333332*a3+6.*d2T+2.*e2T+a1*(0.4666666666666667-0.2*kz)+3.*a2*(-1+kz))*d4*hu+(12.*hdtdT+4.*heteT+a1*(-0.9333333333333333+0.4*kz)*M1-6.*a2*(-1+kz)*M2-42.666666666666664*a3*M3)*d4*u1+(42.666666666666664*a3+24.*d2T+8.*e2T+a2*(18.-6.*kz)+a1*(2.533333333333333+1.2*kz)-12.*u2T)*u2*hddt*u1+(6.666666666666667*a1+12.*a2+85.33333333333333*a3+24.*u2T)*u2*hu*u2t+(42.666666666666664*a3+24.*d2T+8.*e2T+a2*(18.-6.*kz)+a1*(2.533333333333333+1.2*kz)-12.*u2T)*u2*d2*hu+(24.*hdtdT+8.*heteT-12.*hutuT+a1*(-2.533333333333333-1.2*kz)*M1+a2*(-18.+6.*kz)*M2-42.666666666666664*a3*M3)*u1*utd2u+(64.*a3+a2*(15.-3.*kz)+a1*(3.+kz)+18.*u2T)*u4*hu+(24.*hutuT-6.666666666666667*a1*M1-12.*a2*M2-85.33333333333333*a3*M3)*u1*u4t+2.*kz*hddt*d4*u1+8.*hddt*u1*utd2u+6.*huut*d4*u1+6.*huut*d2*u1*u2t-2.*huut*u1*utd2u+(6+4.*kz)*hu*u6t+2.*kz*d2*hddt*d1*dtu+8.*d2*hu*utd2u+2.*kz*d4*hddt*u1+kz*d6*hu+8.*d2*u2*hddt*u1+4.*d2*u2d2*hu+12.*u2*hddt*d1*dtu+4.*u2*hddt*u1*u2t-4.*u2*hu*utd2u+(12+6.*kz)*u2*hu*u4t+12.*u2*d1*dthd*dtu+4.*u2*d2*hu*u2t+12.*u2*d4*hu-4.*u4*hddt*u1+(12+6.*kz)*u4*hu*u2t-4.*u4d2*hu+(12+5.*kz)*u6*hu;
	
	dhd3=(54.*d2T*d4T+18.*d4T*e2T+18.*d2T*e4T+6.*e2T*e4T+a13*(28.837185185185184-1.8573333333333333*kz)+3.*d6T*(1+kz)+e6T*(1+kz)+a23*(172.5+52.5*kz)+a33*(201.4814814814815+106.66666666666667*kz)+a2*(a32*(68.-24.*kz)+(9.*d4T+3.*e4T)*(1+kz)+a3*d2T*(-132.+24.*kz)+18.*u2d2T)+a1*(3.*d4T+9.*e4T+a32*(23.555555555555557-5.866666666666666*kz)+a22*(8.5-2.7*kz)+1.8*(d4T-e4T)*kz+a3*d2T*(-18.933333333333334+3.7333333333333334*kz)+a2*(-1.6*a3-0.3*d2T-8.1*e2T-1.5*d2T*kz+2.7*e2T*kz)-2.4*u2d2T+1.4*kz*u2d2T)+a3*(72.*d4T+24.*u2d2T-kz*(24.*d4T+8.*u2d2T))+a32*(-106.66666666666667*d2T-2.6666666666666665*d2T*kz-53.333333333333336*u2T)+a22*(-78.75*d2T-26.25*e2T+a3*(140.-36.*kz)-15.75*d2T*kz-5.25*e2T*kz-45.*u2T)+a12*(-15.75*d2T-31.65*e2T+a3*(17.297777777777778-4.1066666666666665*kz)+a2*(2.18-1.26*kz)-0.25666666666666665*d2T*kz+0.27*e2T*kz-6.066666666666666*u2T)+18.*u2d2T*u2T+9.*u4d2T)*hd+(216.*d2T*hdtd2dT+72.*e2T*hdtd2dT+108.*d4T*hdtdT+72.*d2T*hete2eT+24.*e2T*hete2eT+e4T*(36.*hdtdT+12.*heteT)+6.*hete4eT*(1+kz)+18.*(hdtd4dT+hdtu4dT+hutd2u2uT+hutu2d2uT+hdtd4dT*kz)+a13*(-173.02311111111112+11.144*kz)*M1+a23*(-1035.-315.*kz)*M2+a33*(-1208.888888888889-640.*kz)*M3+a2*(12.*hete2eT*(1+kz)+36.*(hdtd2dT+hdtu2dT+hutd2uT+hdtd2dT*kz)+a3*(hdtdT*(-264.+48.*kz)+d2T*(264.-48.*kz)*M2+d2T*(264.-48.*kz)*M3)+a32*((-136.+48.*kz)*M2+(-272.+96.*kz)*M3)+M2*(-((18.*d4T+6.*e4T)*(1+kz))-36.*u2d2T))+a1*(12.*hdtd2dT-4.8*hdtu2dT+36.*hete2eT-4.8*hutd2uT+7.2*hdtd2dT*kz+2.8*hdtu2dT*kz-7.2*hete2eT*kz+2.8*hutd2uT*kz+a22*((-17.+5.4*kz)*M1+(-34.+10.8*kz)*M2)+a3*(hdtdT*(-37.86666666666667+7.466666666666667*kz)+d2T*(37.86666666666667-7.466666666666667*kz)*M1+d2T*(37.86666666666667-7.466666666666667*kz)*M3)+a32*((-47.111111111111114+11.733333333333333*kz)*M1+(-94.22222222222223+23.466666666666665*kz)*M3)+a2*(-0.6*hdtdT-16.2*heteT-3.*hdtdT*kz+5.4*heteT*kz+(e2T*(16.2-5.4*kz)+d2T*(0.6+3.*kz))*M1+(e2T*(16.2-5.4*kz)+d2T*(0.6+3.*kz))*M2+a3*(3.2*M1+3.2*M2+3.2*M3))+M1*(-6.*d4T-18.*e4T+3.6*(-d4T+e4T)*kz+4.8*u2d2T-2.8*kz*u2d2T))+a3*(hdtd2dT*(288.-96.*kz)+hdtu2dT*(48.-16.*kz)+hutd2uT*(48.-16.*kz)+M3*(-144.*d4T+48.*d4T*kz-48.*u2d2T+16.*kz*u2d2T))+36.*(d4T*heteT+hutuT*u2d2T+(hdtu2dT+hutd2uT)*u2T)+a12*(-31.5*hdtdT-63.3*heteT-12.133333333333333*hutuT-0.5133333333333333*hdtdT*kz+0.54*heteT*kz+a2*((-8.72+5.04*kz)*M1+(-4.36+2.52*kz)*M2)+a3*((-69.19111111111111+16.426666666666666*kz)*M1+(-34.595555555555556+8.213333333333333*kz)*M3)+M1*(e2T*(126.6-1.08*kz)+d2T*(63.+1.0266666666666666*kz)+24.266666666666666*u2T))+a22*(-157.5*hdtdT-52.5*heteT-90.*hutuT-31.5*hdtdT*kz-10.5*heteT*kz+a3*((-560.+144.*kz)*M2+(-280.+72.*kz)*M3)+M2*(e2T*(105.+21.*kz)+d2T*(315.+63.*kz)+180.*u2T))+a32*(-213.33333333333334*hdtdT-106.66666666666667*hutuT-5.333333333333333*hdtdT*kz+M3*(d2T*(426.6666666666667+10.666666666666666*kz)+213.33333333333334*u2T)))*d1+(-36.*d2T*d2T+72.*d4T-24.*d2T*e2T-4.*e2T*e2T+32.*a3*(e2T+d2T*(-1+kz))+a32*(10.666666666666666-60.44444444444444*kz)+a22*(-66.-12.*kz)+a12*(-23.893333333333334+0.06222222222222222*kz)+a1*(13.4*d2T-6.2*e2T-2.6*d2T*kz+1.8*e2T*kz+a2*(-16.8+2.*kz)+a3*(-14.933333333333334+4.622222222222222*kz))+a2*(d2T*(63.-9.*kz)+e2T*(21.-3.*kz)+a3*(-96.+16.*kz))+24.*(e4T+u2d2T))*hd*d2t+(a32*(5.333333333333333-30.22222222222222*kz)+a22*(-22.5-10.5*kz)+a12*(-25.113333333333333+0.31777777777777777*kz)+12.*u2d2T+16.*a3*(-1+kz)*u2T-18.*u2T*u2T+a2*(-8.*a3+36.*u2T)+a1*(a3*(-27.2+2.8444444444444446*kz)+a2*(-11.8+3.*kz)+(4.-1.6*kz)*u2T)+36.*u4T)*huut*d1+(-45.*d2T*d2T+90.*d4T-30.*d2T*e2T-5.*e2T*e2T+40.*a3*(e2T+d2T*(-1+kz))+a32*(13.333333333333334-75.55555555555556*kz)+a22*(-98.25-8.25*kz)+a12*(-28.796666666666667+0.0077777777777777776*kz)+a1*(17.2*d2T-7.6*e2T-2.8*d2T*kz+2.4*e2T*kz+a2*(-21.3+3.1*kz)+a3*(-14.666666666666666+4.977777777777778*kz))+a2*(d2T*(72.-18.*kz)+e2T*(24.-6.*kz)+a3*(-180.+32.*kz))+30.*(e4T+u2d2T))*d2*hd+(216.*hdtd2dT-108.*d2T*hdtdT+72.*hete2eT-e2T*(36.*hdtdT+12.*heteT)+36.*(hdtu2dT-d2T*heteT+hutd2uT)+a12*(70.25333333333333-0.09333333333333334*kz)*M1+a22*(219.+27.*kz)*M2+a32*(-32.+181.33333333333334*kz)*M3+a3*(48.*(heteT+hdtdT*(-1+kz))-48.*(e2T+d2T*(-1+kz))*M3)+a2*(hdtdT*(90.-18.*kz)+heteT*(30.-6.*kz)+(-90.*d2T-30.*e2T+18.*d2T*kz+6.*e2T*kz)*M2+a3*((184.-32.*kz)*M2+(184.-32.*kz)*M3))+a1*(20.4*hdtdT-9.2*heteT-3.6*hdtdT*kz+2.8*heteT*kz+(-20.4*d2T+9.2*e2T+3.6*d2T*kz-2.8*e2T*kz)*M1+a2*((25.4-3.4*kz)*M1+(25.4-3.4*kz)*M2)+a3*((19.733333333333334-6.4*kz)*M1+(19.733333333333334-6.4*kz)*M3)))*d1*d2t+(a32*(2.6666666666666665-15.11111111111111*kz)+a22*(-11.25-5.25*kz)+a12*(-12.556666666666667+0.15888888888888889*kz)+6.*u2d2T+8.*a3*(-1+kz)*u2T-9.*u2T*u2T+a2*(-4.*a3+18.*u2T)+a1*(a3*(-13.6+1.4222222222222223*kz)+a2*(-5.9+1.5*kz)+(2.-0.8*kz)*u2T)+18.*u4T)*u2*hd+(12.*(hdtu2dT+hutd2uT)+72.*hutu2uT+a12*(50.22666666666667-0.6355555555555555*kz)*M1+a22*(45.+21.*kz)*M2+a32*(-10.666666666666666+60.44444444444444*kz)*M3-36.*hutuT*u2T+a1*(hutuT*(4.-1.6*kz)+a2*((11.8-3.*kz)*M1+(11.8-3.*kz)*M2)+a3*((27.2-2.8444444444444446*kz)*M1+(27.2-2.8444444444444446*kz)*M3)+(-4.+1.6*kz)*M1*u2T)+a2*(36.*hutuT+a3*(8.*M2+8.*M3)-36.*M2*u2T)+a3*(16.*hutuT*(-1+kz)-16.*(-1+kz)*M3*u2T))*u2*d1+(64.*a3+18.*d2T+6.*e2T+a1*(0.6-0.2*kz)+3.*a2*(1+kz))*hd*d4t+(21.333333333333332*a3-6.*d2T-2.*e2T+a2*(9.-3.*kz)+a1*(-1.9333333333333333+0.6*kz)+12.*u2T)*hd*dtu2d+(42.666666666666664*a3+a1*(7.333333333333333-2.*kz)+6.*a2*(-1+kz)+12.*u2T)*huut*u1*utd+(0.26666666666666666*a1+12.*a2+85.33333333333333*a3+24.*d2T+8.*e2T)*d2*hd*d2t+(42.666666666666664*a3-12.*d2T-4.*e2T+a2*(18.-6.*kz)+a1*(-3.8666666666666667+1.2*kz)+24.*u2T)*d2*huut*d1+(64.*a3+18.*d2T+6.*e2T+a2*(15.-3.*kz)+0.2*a1*(-1+kz))*d4*hd+(24.*hdtdT+8.*heteT-0.26666666666666666*a1*M1-12.*a2*M2-85.33333333333333*a3*M3)*d1*d4t+(42.666666666666664*a3-12.*d2T-4.*e2T+a2*(18.-6.*kz)+a1*(-3.8666666666666667+1.2*kz)+24.*u2T)*d2*u2*hd+(-12.*hdtdT-4.*heteT+24.*hutuT+a1*(3.8666666666666667-1.2*kz)*M1+a2*(-18.+6.*kz)*M2-42.666666666666664*a3*M3)*d1*dtu2d+(42.666666666666664*a3+a1*(7.333333333333333-2.*kz)+6.*a2*(-1+kz)+12.*u2T)*u2*huut*d1+(21.333333333333332*a3+a1*(3.6666666666666665-1.*kz)+3.*a2*(-1+kz)+6.*u2T)*u4*hd+(12.*hutuT+a1*(-7.333333333333333+2.*kz)*M1-6.*a2*(-1+kz)*M2-42.666666666666664*a3*M3)*u4*d1+(6+4.*kz)*hd*d6t-2.*hddt*d1*dtu2d+6.*hddt*u2*d1*d2t+6.*hddt*u4*d1+8.*huut*d1*dtu2d+2.*kz*huut*u4*d1+(12+6.*kz)*d2*hd*d4t-4.*d2*hd*dtu2d+4.*d2*huut*d1*d2t+12.*d2*huut*u1*utd+(12+6.*kz)*d4*hd*d2t-4.*d4*huut*d1+(12+5.*kz)*d6*hd-4.*d4u2*hd+4.*d2*u2*hd*d2t+12.*d2*u1*uthu*utd+12.*d2*u4*hd+8.*u2*hd*dtu2d+2.*kz*u2*huut*u1*utd+8.*u2*d2*huut*d1+4.*u2*d2u2*hd+2.*kz*u4*huut*d1+kz*u6*hd;
	
	dhe3=(54.*d2T*d4T+18.*d4T*e2T+18.*d2T*e4T+6.*e2T*e4T+a13*(99.972-7.164*kz)+3.*d6T*(1+kz)+e6T*(1+kz)-a32*d2T*(53.333333333333336+2.6666666666666665*kz)+a23*(172.5+52.5*kz)+a2*((9.*d4T+3.*e4T)*(1+kz)+a3*d2T*(-132.+24.*kz)+18.*u2d2T)+a1*(3.*d4T+9.*e4T+a22*(4.5-2.7*kz)+1.8*(d4T-e4T)*kz+a3*d2T*(-18.933333333333334+3.7333333333333334*kz)+a2*(-0.3*d2T-8.1*e2T-1.5*d2T*kz+2.7*e2T*kz)-2.4*u2d2T+1.4*kz*u2d2T)+a3*(72.*d4T+24.*u2d2T-kz*(24.*d4T+8.*u2d2T))+a22*(-78.75*d2T-26.25*e2T+a3*(180.-36.*kz)-15.75*d2T*kz-5.25*e2T*kz-45.*u2T)+a12*(-25.083333333333332*d2T-43.65*e2T+a3*(79.2-15.84*kz)+a2*(16.74-4.86*kz)-0.25666666666666665*d2T*kz+0.27*e2T*kz-23.4*u2T)+18.*u2d2T*u2T+9.*u4d2T)*he+(216.*d2T*hdtd2dT+72.*e2T*hdtd2dT+108.*d4T*hdtdT+72.*d2T*hete2eT+24.*e2T*hete2eT+e4T*(36.*hdtdT+12.*heteT)+6.*hete4eT*(1+kz)+18.*(hdtd4dT+hdtu4dT+hutd2u2uT+hutu2d2uT+hdtd4dT*kz)+a13*(-599.832+42.984*kz)*M1+a23*(-1035.-315.*kz)*M2+a32*(-(hdtdT*(106.66666666666667+5.333333333333333*kz))+d2T*(213.33333333333334+10.666666666666666*kz)*M3)+a2*(12.*hete2eT*(1+kz)+36.*(hdtd2dT+hdtu2dT+hutd2uT+hdtd2dT*kz)+a3*(hdtdT*(-264.+48.*kz)+d2T*(264.-48.*kz)*M2+d2T*(264.-48.*kz)*M3)+M2*(-((18.*d4T+6.*e4T)*(1+kz))-36.*u2d2T))+a1*(12.*hdtd2dT-4.8*hdtu2dT+36.*hete2eT-4.8*hutd2uT+7.2*hdtd2dT*kz+2.8*hdtu2dT*kz-7.2*hete2eT*kz+2.8*hutd2uT*kz+a22*((-9.+5.4*kz)*M1+(-18.+10.8*kz)*M2)+a2*(-0.6*hdtdT-16.2*heteT-3.*hdtdT*kz+5.4*heteT*kz+(e2T*(16.2-5.4*kz)+d2T*(0.6+3.*kz))*M1+(e2T*(16.2-5.4*kz)+d2T*(0.6+3.*kz))*M2)+a3*(hdtdT*(-37.86666666666667+7.466666666666667*kz)+d2T*(37.86666666666667-7.466666666666667*kz)*M1+d2T*(37.86666666666667-7.466666666666667*kz)*M3)+M1*(-6.*d4T-18.*e4T+3.6*(-d4T+e4T)*kz+4.8*u2d2T-2.8*kz*u2d2T))+a3*(hdtd2dT*(288.-96.*kz)+hdtu2dT*(48.-16.*kz)+hutd2uT*(48.-16.*kz)+M3*(-144.*d4T+48.*d4T*kz-48.*u2d2T+16.*kz*u2d2T))+36.*(d4T*heteT+hutuT*u2d2T+(hdtu2dT+hutd2uT)*u2T)+a12*(-50.166666666666664*hdtdT-87.3*heteT-46.8*hutuT-0.5133333333333333*hdtdT*kz+0.54*heteT*kz+a2*((-66.96+19.44*kz)*M1+(-33.48+9.72*kz)*M2)+a3*((-316.8+63.36*kz)*M1+(-158.4+31.68*kz)*M3)+M1*(e2T*(174.6-1.08*kz)+d2T*(100.33333333333333+1.0266666666666666*kz)+93.6*u2T))+a22*(-157.5*hdtdT-52.5*heteT-90.*hutuT-31.5*hdtdT*kz-10.5*heteT*kz+a3*((-720.+144.*kz)*M2+(-360.+72.*kz)*M3)+M2*(e2T*(105.+21.*kz)+d2T*(315.+63.*kz)+180.*u2T)))*e1+(-36.*d2T*d2T+72.*d4T-24.*d2T*e2T-4.*e2T*e2T+a2*(d2T*(63.-9.*kz)+e2T*(21.-3.*kz))+a22*(-66.-12.*kz)+a12*(-84.96-2.16*kz)+a3*d2T*(-128.+32.*kz)+a1*(d2T*(37.4-2.6*kz)+1.8*e2T*(1+kz)+a2*(-43.2+10.8*kz))+24.*(e4T+u2d2T))*he*e2t+(-45.*d2T*d2T+90.*d4T-30.*d2T*e2T-5.*e2T*e2T+a2*(d2T*(72.-18.*kz)+e2T*(24.-6.*kz))+a22*(-98.25-8.25*kz)+a12*(-87.57-5.13*kz)+a3*d2T*(-160.+40.*kz)+a1*(d2T*(50.8+0.8*kz)+3.6*e2T*(1+kz)+a2*(-62.1+13.5*kz))+30.*(e4T+u2d2T))*e2*he+(216.*hdtd2dT-108.*d2T*hdtdT+72.*hete2eT-e2T*(36.*hdtdT+12.*heteT)+36.*(hdtu2dT-d2T*heteT+hutd2uT)+a12*(230.04+9.72*kz)*M1+a22*(219.+27.*kz)*M2+a2*(hdtdT*(90.-18.*kz)+heteT*(30.-6.*kz)+(-90.*d2T-30.*e2T+18.*d2T*kz+6.*e2T*kz)*M2)+a1*(hdtdT*(58.8-1.2*kz)+3.6*heteT*(1+kz)+(-3.6*e2T*(1+kz)+d2T*(-58.8+1.2*kz))*M1+a2*((70.2-16.2*kz)*M1+(70.2-16.2*kz)*M2))+a3*(hdtdT*(-192.+48.*kz)+d2T*(192.-48.*kz)*M3))*e1*e2t+(18.*d2T+6.*e2T+a1*(19.8-1.8*kz)+3.*a2*(1+kz))*he*e4t+(21.6*a1+12.*a2+24.*d2T+8.*e2T)*e2*he*e2t+(18.*d2T+6.*e2T+a2*(15.-3.*kz)+a1*(12.6+1.8*kz))*e4*he+(24.*hdtdT+8.*heteT-21.6*a1*M1-12.*a2*M2)*e1*e4t+(6+4.*kz)*he*e6t+(12+6.*kz)*e2*he*e4t+(12+6.*kz)*e4*he*e2t+(12+5.*kz)*e6*he;
	
	dmq3=a23*((2157.+630.*kz)*M22-3.*mH1sq-3.*mH2sq-3.*mlT-9.*mqT)+a13*((51.12088888888889-3.184*kz)*M12-0.3537777777777778*mdT-1.0613333333333332*meT-0.5306666666666666*mH1sq-0.5306666666666666*mH2sq-0.5306666666666666*mlT-0.1768888888888889*mqT-1.4151111111111112*muT)+a2*a32*((80.-48.*kz)*M22+(64.-96.*kz)*M2*M3+(192.-144.*kz)*M32-32.*mqT-16.*(mdT+muT))+a33*((2279.1111111111113+1280.*kz)*M32+71.11111111111111*mqT+35.55555555555556*(mdT+muT))+a1*(a2*a3*(-6.4*M12-6.4*M22+M1*(-6.4*M2-6.4*M3)-6.4*M2*M3-6.4*M32)+a22*((30.4-5.4*kz)*M12+(50.-10.8*kz)*M1*M2+(75.8-16.2*kz)*M22-0.2*mH1sq-0.2*mH2sq-0.2*mlT-0.6*mqT)+a32*((65.77777777777777-11.733333333333333*kz)*M12+(108.08888888888889-23.466666666666665*kz)*M1*M3+(164.26666666666668-35.2*kz)*M32-0.7111111111111111*mqT-0.35555555555555557*(mdT+muT)))+a12*(-0.56*(d2mqT+d2tmdT+hd2T)-0.72*(e2mlT+e2tmeT+he2T)+(1.8666666666666667*hdtdT+2.4*heteT+3.466666666666667*hutuT)*M1+(-0.56*d2T-0.72*e2T)*mH1sq+a3*((10.346666666666666-3.52*kz)*M12+(6.897777777777778-2.3466666666666667*kz)*M1*M3+(4.622222222222222-1.1733333333333333*kz)*M32-0.14222222222222222*mdT-0.4266666666666667*meT-0.21333333333333335*mH1sq-0.21333333333333335*mH2sq-0.21333333333333335*mlT-0.07111111111111111*mqT-0.5688888888888889*muT)+a2*((1.32-1.08*kz)*M12+(0.88-0.72*kz)*M1*M2+(0.8-0.36*kz)*M22-0.08*mdT-0.24*meT-0.12*mH1sq-0.12*mH2sq-0.12*mlT-0.04*mqT-0.32*muT)+M12*(-2.8*d2T-3.6*e2T-5.2*u2T)-1.04*mH2sq*u2T-1.04*(hu2T+u2mqT+u2tmuT))+a32*(213.33333333333334*(hdtdT+hutuT)*M3-64.*d2T*mH1sq-64.*mH2sq*u2T-320.*M32*(d2T+u2T)-64.*(d2mqT+d2tmdT+hd2T+hu2T+u2mqT+u2tmuT))+a22*(-18.*(e2mlT+e2tmeT+he2T)+(60.*heteT+180.*(hdtdT+hutuT))*M2+(-54.*d2T-18.*e2T)*mH1sq+a3*((664.-216.*kz)*M22+(400.-144.*kz)*M2*M3+(272.-72.*kz)*M32-16.*mH1sq-16.*mH2sq-16.*mlT-48.*mqT)-54.*mH2sq*u2T+M22*(-90.*e2T-270.*(d2T+u2T))-54.*(d2mqT+d2tmdT+hd2T+hu2T+u2mqT+u2tmuT))+(-18.*d2T*d2T+36.*d4T-12.*d2T*e2T-2.*e2T*e2T+a2*(-8.*a3+36.*d2T+12.*e2T)+16.*a3*(e2T+d2T*(-1+kz))+a32*(5.333333333333333-30.22222222222222*kz)+a22*(-22.5-10.5*kz)+a12*(-12.66+0.07777777777777778*kz)+a1*(6.4*d2T-3.2*e2T+a2*(-8.2+0.6*kz)-1.6*d2T*kz+0.8*e2T*kz+a3*(-10.133333333333333+2.8444444444444446*kz))+12.*(e4T+u2d2T))*hd2+(72.*hdtd2dT-36.*d2T*hdtdT+24.*hete2eT-e2T*(12.*hdtdT+4.*heteT)+12.*(hdtu2dT-d2T*heteT+hutd2uT)+a12*(25.32-0.15555555555555556*kz)*M1+a22*(45.+21.*kz)*M2+a32*(-10.666666666666666+60.44444444444444*kz)*M3+a3*(16.*(heteT+hdtdT*(-1+kz))-16.*(e2T+d2T*(-1+kz))*M3)+a2*(36.*hdtdT+12.*heteT+(-36.*d2T-12.*e2T)*M2+a3*(8.*M2+8.*M3))+a1*(6.4*hdtdT-3.2*heteT-1.6*hdtdT*kz+0.8*heteT*kz+(-6.4*d2T+3.2*e2T+1.6*d2T*kz-0.8*e2T*kz)*M1+a2*((8.2-0.6*kz)*M1+(8.2-0.6*kz)*M2)+a3*((10.133333333333333-2.8444444444444446*kz)*M1+(10.133333333333333-2.8444444444444446*kz)*M3)))*hddt+(a32*(5.333333333333333-30.22222222222222*kz)+a22*(-22.5-10.5*kz)+a12*(-25.113333333333333+0.31777777777777777*kz)+12.*u2d2T+16.*a3*(-1+kz)*u2T-18.*u2T*u2T+a2*(-8.*a3+36.*u2T)+a1*(a3*(-27.2+2.8444444444444446*kz)+a2*(-11.8+3.*kz)+(4.-1.6*kz)*u2T)+36.*u4T)*hu2+(12.*(hdtu2dT+hutd2uT)+72.*hutu2uT+a12*(50.22666666666667-0.6355555555555555*kz)*M1+a22*(45.+21.*kz)*M2+a32*(-10.666666666666666+60.44444444444444*kz)*M3-36.*hutuT*u2T+a1*(hutuT*(4.-1.6*kz)+a2*((11.8-3.*kz)*M1+(11.8-3.*kz)*M2)+a3*((27.2-2.8444444444444446*kz)*M1+(27.2-2.8444444444444446*kz)*M3)+(-4.+1.6*kz)*M1*u2T)+a2*(36.*hutuT+a3*(8.*M2+8.*M3)-36.*M2*u2T)+a3*(16.*hutuT*(-1+kz)-16.*(-1+kz)*M3*u2T))*huut+(72.*hdtd2dT-36.*d2T*hdtdT+24.*hete2eT-e2T*(12.*hdtdT+4.*heteT)+12.*(hdtu2dT-d2T*heteT+hutd2uT)+a12*(25.32-0.15555555555555556*kz)*M1+a22*(45.+21.*kz)*M2+a32*(-10.666666666666666+60.44444444444444*kz)*M3+a3*(16.*(heteT+hdtdT*(-1+kz))-16.*(e2T+d2T*(-1+kz))*M3)+a2*(36.*hdtdT+12.*heteT+(-36.*d2T-12.*e2T)*M2+a3*(8.*M2+8.*M3))+a1*(6.4*hdtdT-3.2*heteT-1.6*hdtdT*kz+0.8*heteT*kz+(-6.4*d2T+3.2*e2T+1.6*d2T*kz-0.8*e2T*kz)*M1+a2*((8.2-0.6*kz)*M1+(8.2-0.6*kz)*M2)+a3*((10.133333333333333-2.8444444444444446*kz)*M1+(10.133333333333333-2.8444444444444446*kz)*M3)))*dhdt+(72.*(d2dmddtT+d4mqT)+24.*(e2e2mlT+e2emeetT)-e2mlT*(12.*d2T+4.*e2T)-36.*d2T*(d2mqT+d2tmdT+hd2T)-12.*e2T*(d2mqT+d2tmdT+hd2T)+144.*hd2td2tT-36.*hdtdT*hdtdT-4.*e2T*(e2tmeT+he2T)+48.*he2te2tT-24.*hdtdT*heteT-4.*heteT*heteT+a2*(36.*(d2mqT+d2tmdT+hd2T)+12.*(e2mlT+e2tmeT+he2T)+(-72.*hdtdT-24.*heteT)*M2+(72.*d2T+24.*e2T)*M22+a3*(-16.*M22-16.*M2*M3-16.*M32-8.*mH1sq)+(72.*d2T+24.*e2T)*mH1sq)+a3*(16.*(e2mlT+e2tmeT+he2T+d2mqT*(-1+kz)+d2tmdT*(-1+kz)+hd2T*(-1+kz))-32.*(heteT+hdtdT*(-1+kz))*M3+32.*(e2T+d2T*(-1+kz))*M32+32.*(e2T+d2T*(-1+kz))*mH1sq)+a32*((32.-181.33333333333334*kz)*M32+(5.333333333333333-30.22222222222222*kz)*mH1sq)+a22*((-135.-63.*kz)*M22+(-22.5-10.5*kz)*mH1sq)+a1*(-3.2*e2mlT-3.2*e2tmeT+6.4*(d2mqT+d2tmdT+hd2T)-3.2*he2T-1.6*d2mqT*kz-1.6*d2tmdT*kz+0.8*e2mlT*kz+0.8*e2tmeT*kz-1.6*hd2T*kz+0.8*he2T*kz+(-12.8*hdtdT+6.4*heteT+3.2*hdtdT*kz-1.6*heteT*kz)*M1+(12.8*d2T-6.4*e2T-3.2*d2T*kz+1.6*e2T*kz)*M12+(12.8*d2T-6.4*e2T-3.2*d2T*kz+1.6*e2T*kz)*mH1sq+a2*((-16.4+1.2*kz)*M12+(-16.4+1.2*kz)*M1*M2+(-16.4+1.2*kz)*M22+(-8.2+0.6*kz)*mH1sq)+a3*((-20.266666666666666+5.688888888888889*kz)*M12+(-20.266666666666666+5.688888888888889*kz)*M1*M3+(-20.266666666666666+5.688888888888889*kz)*M32+(-10.133333333333333+2.8444444444444446*kz)*mH1sq))+a12*((-75.96+0.4666666666666667*kz)*M12-0.32*mdT-0.96*meT+(-13.14+0.07777777777777778*kz)*mH1sq-0.48*mH2sq-0.48*mlT-0.16*mqT-1.28*muT)+12.*mH2sq*u2d2T+mH1sq*(-54.*d2T*d2T+108.*d4T-36.*d2T*e2T-6.*e2T*e2T+36.*e4T+24.*u2d2T)+12.*(d2u2mqT+d2umuutT+hd2u2T-d2T*(e2tmeT+he2T)+hu2d2T+u2d2mqT+u2dmddtT)+24.*uhuthddtT)*d2+(12.*(hdtu2dT+hutd2uT)+72.*hutu2uT+a12*(50.22666666666667-0.6355555555555555*kz)*M1+a22*(45.+21.*kz)*M2+a32*(-10.666666666666666+60.44444444444444*kz)*M3-36.*hutuT*u2T+a1*(hutuT*(4.-1.6*kz)+a2*((11.8-3.*kz)*M1+(11.8-3.*kz)*M2)+a3*((27.2-2.8444444444444446*kz)*M1+(27.2-2.8444444444444446*kz)*M3)+(-4.+1.6*kz)*M1*u2T)+a2*(36.*hutuT+a3*(8.*M2+8.*M3)-36.*M2*u2T)+a3*(16.*hutuT*(-1+kz)-16.*(-1+kz)*M3*u2T))*uhut+(144.*hu2tu2tT+a32*((32.-181.33333333333334*kz)*M32+(5.333333333333333-30.22222222222222*kz)*mH2sq)+a22*((-135.-63.*kz)*M22+(-22.5-10.5*kz)*mH2sq)+a12*((-150.68+1.9066666666666667*kz)*M12-0.64*mdT-1.92*meT-0.96*mH1sq+(-26.073333333333334+0.31777777777777777*kz)*mH2sq-0.96*mlT-0.32*mqT-2.56*muT)+12.*mH1sq*u2d2T+12.*(d2u2mqT+d2umuutT+hd2u2T+hu2d2T+u2d2mqT+u2dmddtT)+a2*(-72.*hutuT*M2+a3*(-16.*M22-16.*M2*M3-16.*M32-8.*mH2sq)+72.*M22*u2T+72.*mH2sq*u2T+36.*(hu2T+u2mqT+u2tmuT))+a3*(-32.*hutuT*(-1+kz)*M3+32.*(-1+kz)*M32*u2T+32.*(-1+kz)*mH2sq*u2T+16.*(-1+kz)*(hu2T+u2mqT+u2tmuT))+a1*(hutuT*(-8.+3.2*kz)*M1+a3*((-54.4+5.688888888888889*kz)*M12+(-54.4+5.688888888888889*kz)*M1*M3+(-54.4+5.688888888888889*kz)*M32+(-27.2+2.8444444444444446*kz)*mH2sq)+a2*((-23.6+6.*kz)*M12+(-23.6+6.*kz)*M1*M2+(-23.6+6.*kz)*M22+(-11.8+3.*kz)*mH2sq)+(8.-3.2*kz)*M12*u2T+(8.-3.2*kz)*mH2sq*u2T-(-4.+1.6*kz)*(hu2T+u2mqT+u2tmuT))-36.*(hutuT*hutuT+u2T*(hu2T+u2mqT+u2tmuT))+72.*(u2umuutT+u4mqT)+mH2sq*(24.*u2d2T-54.*u2T*u2T+108.*u4T)+24.*uhuthddtT)*u2+(-9.*d2T*d2T+18.*d4T-6.*d2T*e2T-1.*e2T*e2T+a2*(-4.*a3+18.*d2T+6.*e2T)+8.*a3*(e2T+d2T*(-1+kz))+a32*(2.6666666666666665-15.11111111111111*kz)+a22*(-11.25-5.25*kz)+a12*(-6.33+0.03888888888888889*kz)+a1*(3.2*d2T-1.6*e2T+a2*(-4.1+0.3*kz)-0.8*d2T*kz+0.4*e2T*kz+a3*(-5.066666666666666+1.4222222222222223*kz))+6.*(e4T+u2d2T))*mq*d2+(a32*(2.6666666666666665-15.11111111111111*kz)+a22*(-11.25-5.25*kz)+a12*(-12.556666666666667+0.15888888888888889*kz)+6.*u2d2T+8.*a3*(-1+kz)*u2T-9.*u2T*u2T+a2*(-4.*a3+18.*u2T)+a1*(a3*(-13.6+1.4222222222222223*kz)+a2*(-5.9+1.5*kz)+(2.-0.8*kz)*u2T)+18.*u4T)*mq*u2+(-18.*d2T*d2T+36.*d4T-12.*d2T*e2T-2.*e2T*e2T+a2*(-8.*a3+36.*d2T+12.*e2T)+16.*a3*(e2T+d2T*(-1+kz))+a32*(5.333333333333333-30.22222222222222*kz)+a22*(-22.5-10.5*kz)+a12*(-12.66+0.07777777777777778*kz)+a1*(6.4*d2T-3.2*e2T+a2*(-8.2+0.6*kz)-1.6*d2T*kz+0.8*e2T*kz+a3*(-10.133333333333333+2.8444444444444446*kz))+12.*(e4T+u2d2T))*dmddt+(-9.*d2T*d2T+18.*d4T-6.*d2T*e2T-1.*e2T*e2T+a2*(-4.*a3+18.*d2T+6.*e2T)+8.*a3*(e2T+d2T*(-1+kz))+a32*(2.6666666666666665-15.11111111111111*kz)+a22*(-11.25-5.25*kz)+a12*(-6.33+0.03888888888888889*kz)+a1*(3.2*d2T-1.6*e2T+a2*(-4.1+0.3*kz)-0.8*d2T*kz+0.4*e2T*kz+a3*(-5.066666666666666+1.4222222222222223*kz))+6.*(e4T+u2d2T))*d2*mq+(a32*(5.333333333333333-30.22222222222222*kz)+a22*(-22.5-10.5*kz)+a12*(-25.113333333333333+0.31777777777777777*kz)+12.*u2d2T+16.*a3*(-1+kz)*u2T-18.*u2T*u2T+a2*(-8.*a3+36.*u2T)+a1*(a3*(-27.2+2.8444444444444446*kz)+a2*(-11.8+3.*kz)+(4.-1.6*kz)*u2T)+36.*u4T)*umuut+(a32*(2.6666666666666665-15.11111111111111*kz)+a22*(-11.25-5.25*kz)+a12*(-12.556666666666667+0.15888888888888889*kz)+6.*u2d2T+8.*a3*(-1+kz)*u2T-9.*u2T*u2T+a2*(-4.*a3+18.*u2T)+a1*(a3*(-13.6+1.4222222222222223*kz)+a2*(-5.9+1.5*kz)+(2.-0.8*kz)*u2T)+18.*u4T)*u2*mq+(42.666666666666664*a3+12.*d2T+4.*e2T+a1*(0.9333333333333333-0.4*kz)+6.*a2*(-1+kz))*hd2*d2+(42.666666666666664*a3+12.*d2T+4.*e2T+a1*(0.9333333333333333-0.4*kz)+6.*a2*(-1+kz))*hddt*dhdt+(12.*hdtdT+4.*heteT+a1*(-0.9333333333333333+0.4*kz)*M1-6.*a2*(-1+kz)*M2-42.666666666666664*a3*M3)*hddt*d2+(42.666666666666664*a3+a1*(7.333333333333333-2.*kz)+6.*a2*(-1+kz)+12.*u2T)*hu2*u2+(42.666666666666664*a3+a1*(7.333333333333333-2.*kz)+6.*a2*(-1+kz)+12.*u2T)*huut*uhut+(12.*hutuT+a1*(-7.333333333333333+2.*kz)*M1-6.*a2*(-1+kz)*M2-42.666666666666664*a3*M3)*huut*u2+(42.666666666666664*a3+12.*d2T+4.*e2T+a1*(0.9333333333333333-0.4*kz)+6.*a2*(-1+kz))*dhd2tdt+(12.*hdtdT+4.*heteT+a1*(-0.9333333333333333+0.4*kz)*M1-6.*a2*(-1+kz)*M2-42.666666666666664*a3*M3)*dhdtd2+(42.666666666666664*a3+12.*d2T+4.*e2T+a1*(0.9333333333333333-0.4*kz)+6.*a2*(-1+kz))*d2*hd2+(12.*hdtdT+4.*heteT+a1*(-0.9333333333333333+0.4*kz)*M1-6.*a2*(-1+kz)*M2-42.666666666666664*a3*M3)*d2*hddt+(12.*hdtdT+4.*heteT+a1*(-0.9333333333333333+0.4*kz)*M1-6.*a2*(-1+kz)*M2-42.666666666666664*a3*M3)*d2*dhdt+(12.*(d2mqT+d2tmdT+hd2T)+4.*(e2mlT+e2tmeT+he2T)+(36.*d2T+12.*e2T)*mH1sq+a3*(85.33333333333333*M32+85.33333333333333*mH1sq)+a1*((1.8666666666666667-0.8*kz)*M12+(1.8666666666666667-0.8*kz)*mH1sq)+a2*(12.*(-1+kz)*M22+12.*(-1+kz)*mH1sq))*d4+(42.666666666666664*a3+a1*(7.333333333333333-2.*kz)+6.*a2*(-1+kz)+12.*u2T)*uhu2tut+(12.*hutuT+a1*(-7.333333333333333+2.*kz)*M1-6.*a2*(-1+kz)*M2-42.666666666666664*a3*M3)*uhut*u2+(42.666666666666664*a3+a1*(7.333333333333333-2.*kz)+6.*a2*(-1+kz)+12.*u2T)*u2*hu2+(12.*hutuT+a1*(-7.333333333333333+2.*kz)*M1-6.*a2*(-1+kz)*M2-42.666666666666664*a3*M3)*u2*huut+(12.*hutuT+a1*(-7.333333333333333+2.*kz)*M1-6.*a2*(-1+kz)*M2-42.666666666666664*a3*M3)*u2*uhut+(a3*(85.33333333333333*M32+85.33333333333333*mH2sq)+a1*((14.666666666666666-4.*kz)*M12+(14.666666666666666-4.*kz)*mH2sq)+a2*(12.*(-1+kz)*M22+12.*(-1+kz)*mH2sq)+36.*mH2sq*u2T+12.*(hu2T+u2mqT+u2tmuT))*u4+(21.333333333333332*a3+6.*d2T+2.*e2T+a1*(0.4666666666666667-0.2*kz)+3.*a2*(-1+kz))*mq*d4+(21.333333333333332*a3+a1*(3.6666666666666665-1.*kz)+3.*a2*(-1+kz)+6.*u2T)*mq*u4+(42.666666666666664*a3+12.*d2T+4.*e2T+a1*(0.9333333333333333-0.4*kz)+6.*a2*(-1+kz))*dmddt*d2+(42.666666666666664*a3+12.*d2T+4.*e2T+a1*(0.9333333333333333-0.4*kz)+6.*a2*(-1+kz))*d2*mq*d2+(42.666666666666664*a3+12.*d2T+4.*e2T+a1*(0.9333333333333333-0.4*kz)+6.*a2*(-1+kz))*d2dmddt+(21.333333333333332*a3+6.*d2T+2.*e2T+a1*(0.4666666666666667-0.2*kz)+3.*a2*(-1+kz))*d4mq+(42.666666666666664*a3+a1*(7.333333333333333-2.*kz)+6.*a2*(-1+kz)+12.*u2T)*umuut*u2+(42.666666666666664*a3+a1*(7.333333333333333-2.*kz)+6.*a2*(-1+kz)+12.*u2T)*u2*mq*u2+(42.666666666666664*a3+a1*(7.333333333333333-2.*kz)+6.*a2*(-1+kz)+12.*u2T)*u2*umuut+(21.333333333333332*a3+a1*(3.6666666666666665-1.*kz)+3.*a2*(-1+kz)+6.*u2T)*u4mq+2.*kz*hd2d4+8.*hd2u2*d2+2.*kz*hddt*dhdtd2+2.*kz*hddt*d2*dhdt+8.*hddt*uhutd2+8.*hddt*u2*dhdt+8.*hu2d2*u2+2.*kz*hu2u4+8.*huut*dhdtu2+8.*huut*d2*uhut+2.*kz*huut*uhut*u2+2.*kz*huut*u2*uhut+2.*kz*dhd2tdt*d2+8.*dhdt*huut*d2+2.*kz*dhdt*d1*dthd*dt+8.*dhdtu2hddt+2.*kz*d2*hd2*d2+2.*kz*d2*hddt*dhdt+8.*d2*hu2d2+8.*d2*huut*dhdt+2.*kz*d2*dhd2tdt+2.*kz*d4*hd2+6.*kz*mH1sq*d6+8.*d2uhuthddt+8.*d2*u2*hd2+(16.*mH1sq+8.*mH2sq)*d2*u2*d2+8.*uhuthddt*u2+2.*kz*uhu2tut*u2+8.*uhutd2*huut+2.*kz*uhutu2huut+8.*u2*hd2u2+8.*u2*hddt*uhut+2.*kz*u2*hu2*u2+2.*kz*u2*huut*uhut+8.*u2*dhdt*huut+8.*u2*d2*hu2+(8.*mH1sq+16.*mH2sq)*u2*d2*u2+2.*kz*u2*uhu2tut+2.*kz*u4*hu2+6.*kz*mH2sq*u6+kz*mq*d6+4.*mq*d2u2d2+4.*mq*u2d2u2+kz*mq*u6+2.*kz*dmddt*d4+8.*dmddt*u2*d2+2.*kz*d2*mq*d4+8.*d2*mq*u2*d2+2.*kz*d2dmddt*d2+2.*kz*d4mq*d2+2.*kz*d4*dmddt+kz*d6mq+8.*d2umuut*d2+8.*d2*u2*mq*d2+8.*d2*u2dmddt+4.*d2*u2d2mq+8.*umuut*d2*u2+2.*kz*umuut*u4+8.*u2*mq*d2*u2+2.*kz*u2*mq*u4+8.*u2dmddt*u2+8.*u2*d2*mq*u2+8.*u2*d2umuut+4.*u2*d2u2mq+2.*kz*u2*umuut*u2+2.*kz*u4mq*u2+2.*kz*u4*umuut+kz*u6mq;
	
	dml3=a23*((2157.+630.*kz)*M22-3.*mH1sq-3.*mH2sq-3.*mlT-9.*mqT)+a1*a22*((14.4-5.4*kz)*M12+(18.-10.8*kz)*M1*M2+(34.2-16.2*kz)*M22-1.8*mH1sq-1.8*mH2sq-1.8*mlT-5.4*mqT)+a13*((446.136-28.656*kz)*M12-3.312*mdT-9.936*meT-4.968*mH1sq-4.968*mH2sq-4.968*mlT-1.656*mqT-13.248*muT)+a12*(-5.04*(d2mqT+d2tmdT+hd2T)-6.48*(e2mlT+e2tmeT+he2T)+(16.8*hdtdT+21.6*heteT+31.2*hutuT)*M1+a3*((158.4-31.68*kz)*M12+(105.6-21.12*kz)*M1*M3+(63.36-10.56*kz)*M32)+(-5.04*d2T-6.48*e2T)*mH1sq+a2*((3.24-9.72*kz)*M12+(2.16-6.48*kz)*M1*M2+(4.32-3.24*kz)*M22-0.72*mdT-2.16*meT-1.08*mH1sq-1.08*mH2sq-1.08*mlT-0.36*mqT-2.88*muT)+M12*(-25.2*d2T-32.4*e2T-46.8*u2T)-9.36*mH2sq*u2T-9.36*(hu2T+u2mqT+u2tmuT))+a22*(-18.*(e2mlT+e2tmeT+he2T)+(60.*heteT+180.*(hdtdT+hutuT))*M2+a3*((1080.-216.*kz)*M22+(720.-144.*kz)*M2*M3+(432.-72.*kz)*M32)+(-54.*d2T-18.*e2T)*mH1sq-54.*mH2sq*u2T+M22*(-90.*e2T-270.*(d2T+u2T))-54.*(d2mqT+d2tmdT+hd2T+hu2T+u2mqT+u2tmuT))+(-18.*d2T*d2T+36.*d4T-12.*d2T*e2T-2.*e2T*e2T+a2*(36.*d2T+12.*e2T)+a22*(-22.5-10.5*kz)+a12*(-54.9+0.54*kz)+a3*d2T*(-64.+16.*kz)+a1*(d2T*(16.-4.*kz)+a2*(-16.2+5.4*kz))+12.*(e4T+u2d2T))*he2+(72.*hdtd2dT-36.*d2T*hdtdT+24.*hete2eT-e2T*(12.*hdtdT+4.*heteT)+12.*(hdtu2dT-d2T*heteT+hutd2uT)+a12*(109.8-1.08*kz)*M1+a22*(45.+21.*kz)*M2+a2*(36.*hdtdT+12.*heteT+(-36.*d2T-12.*e2T)*M2)+a1*(hdtdT*(16.-4.*kz)+d2T*(-16.+4.*kz)*M1+a2*((16.2-5.4*kz)*M1+(16.2-5.4*kz)*M2))+a3*(hdtdT*(-64.+16.*kz)+d2T*(64.-16.*kz)*M3))*heet+(72.*hdtd2dT-36.*d2T*hdtdT+24.*hete2eT-e2T*(12.*hdtdT+4.*heteT)+12.*(hdtu2dT-d2T*heteT+hutd2uT)+a12*(109.8-1.08*kz)*M1+a22*(45.+21.*kz)*M2+a2*(36.*hdtdT+12.*heteT+(-36.*d2T-12.*e2T)*M2)+a1*(hdtdT*(16.-4.*kz)+d2T*(-16.+4.*kz)*M1+a2*((16.2-5.4*kz)*M1+(16.2-5.4*kz)*M2))+a3*(hdtdT*(-64.+16.*kz)+d2T*(64.-16.*kz)*M3))*ehet+(72.*(d2dmddtT+d4mqT)+24.*(e2e2mlT+e2emeetT)-e2mlT*(12.*d2T+4.*e2T)-36.*d2T*(d2mqT+d2tmdT+hd2T)-12.*e2T*(d2mqT+d2tmdT+hd2T)+144.*hd2td2tT-36.*hdtdT*hdtdT-4.*e2T*(e2tmeT+he2T)+48.*he2te2tT-24.*hdtdT*heteT-4.*heteT*heteT+a2*(36.*(d2mqT+d2tmdT+hd2T)+12.*(e2mlT+e2tmeT+he2T)+(-72.*hdtdT-24.*heteT)*M2+(72.*d2T+24.*e2T)*M22+(72.*d2T+24.*e2T)*mH1sq)+a22*((-135.-63.*kz)*M22+(-22.5-10.5*kz)*mH1sq)+a3*((d2mqT+d2tmdT+hd2T)*(-64.+16.*kz)+hdtdT*(128.-32.*kz)*M3+d2T*(-128.+32.*kz)*M32+d2T*(-128.+32.*kz)*mH1sq)+a1*(-((d2mqT+d2tmdT+hd2T)*(-16.+4.*kz))+hdtdT*(-32.+8.*kz)*M1+d2T*(32.-8.*kz)*M12+d2T*(32.-8.*kz)*mH1sq+a2*((-32.4+10.8*kz)*M12+(-32.4+10.8*kz)*M1*M2+(-32.4+10.8*kz)*M22+(-16.2+5.4*kz)*mH1sq))+a12*((-329.4+3.24*kz)*M12-0.96*mdT-2.88*meT+(-56.34+0.54*kz)*mH1sq-1.44*mH2sq-1.44*mlT-0.48*mqT-3.84*muT)+12.*mH2sq*u2d2T+mH1sq*(-54.*d2T*d2T+108.*d4T-36.*d2T*e2T-6.*e2T*e2T+36.*e4T+24.*u2d2T)+12.*(d2u2mqT+d2umuutT+hd2u2T-d2T*(e2tmeT+he2T)+hu2d2T+u2d2mqT+u2dmddtT)+24.*uhuthddtT)*e2+(-9.*d2T*d2T+18.*d4T-6.*d2T*e2T-1.*e2T*e2T+a2*(18.*d2T+6.*e2T)+a22*(-11.25-5.25*kz)+a12*(-27.45+0.27*kz)+a3*d2T*(-32.+8.*kz)+a1*(d2T*(8.-2.*kz)+a2*(-8.1+2.7*kz))+6.*(e4T+u2d2T))*ml*e2+(-18.*d2T*d2T+36.*d4T-12.*d2T*e2T-2.*e2T*e2T+a2*(36.*d2T+12.*e2T)+a22*(-22.5-10.5*kz)+a12*(-54.9+0.54*kz)+a3*d2T*(-64.+16.*kz)+a1*(d2T*(16.-4.*kz)+a2*(-16.2+5.4*kz))+12.*(e4T+u2d2T))*emeet+(-9.*d2T*d2T+18.*d4T-6.*d2T*e2T-1.*e2T*e2T+a2*(18.*d2T+6.*e2T)+a22*(-11.25-5.25*kz)+a12*(-27.45+0.27*kz)+a3*d2T*(-32.+8.*kz)+a1*(d2T*(8.-2.*kz)+a2*(-8.1+2.7*kz))+6.*(e4T+u2d2T))*e2*ml+(12.*d2T+4.*e2T+a1*(18.-3.6*kz)+6.*a2*(-1+kz))*he2*e2+(12.*d2T+4.*e2T+a1*(18.-3.6*kz)+6.*a2*(-1+kz))*heet*ehet+(12.*hdtdT+4.*heteT+a1*(-18.+3.6*kz)*M1-6.*a2*(-1+kz)*M2)*heet*e2+(12.*d2T+4.*e2T+a1*(18.-3.6*kz)+6.*a2*(-1+kz))*ehe2tet+(12.*hdtdT+4.*heteT+a1*(-18.+3.6*kz)*M1-6.*a2*(-1+kz)*M2)*ehete2+(12.*d2T+4.*e2T+a1*(18.-3.6*kz)+6.*a2*(-1+kz))*e2*he2+(12.*hdtdT+4.*heteT+a1*(-18.+3.6*kz)*M1-6.*a2*(-1+kz)*M2)*e2*heet+(12.*hdtdT+4.*heteT+a1*(-18.+3.6*kz)*M1-6.*a2*(-1+kz)*M2)*e2*ehet+(12.*(d2mqT+d2tmdT+hd2T)+4.*(e2mlT+e2tmeT+he2T)+(36.*d2T+12.*e2T)*mH1sq+a1*((36.-7.2*kz)*M12+(36.-7.2*kz)*mH1sq)+a2*(12.*(-1+kz)*M22+12.*(-1+kz)*mH1sq))*e4+(6.*d2T+2.*e2T+a1*(9.-1.8*kz)+3.*a2*(-1+kz))*ml*e4+(12.*d2T+4.*e2T+a1*(18.-3.6*kz)+6.*a2*(-1+kz))*emeet*e2+(12.*d2T+4.*e2T+a1*(18.-3.6*kz)+6.*a2*(-1+kz))*e2*ml*e2+(12.*d2T+4.*e2T+a1*(18.-3.6*kz)+6.*a2*(-1+kz))*e2emeet+(6.*d2T+2.*e2T+a1*(9.-1.8*kz)+3.*a2*(-1+kz))*e4ml+2.*kz*he2e4+2.*kz*heet*ehete2+2.*kz*heet*e2*ehet+2.*kz*ehe2tet*e2+2.*kz*ehete2heet+2.*kz*e2*he2*e2+2.*kz*e2*heet*ehet+2.*kz*e2*ehe2tet+2.*kz*e4*he2+6.*kz*mH1sq*e6+kz*ml*e6+2.*kz*emeet*e4+2.*kz*e2*ml*e4+2.*kz*e2emeet*e2+2.*kz*e4ml*e2+2.*kz*e4*emeet+kz*e6ml;
	
	dmd3=a2*a32*((288.-48.*kz)*M22+(480.-96.*kz)*M2*M3+(720.-144.*kz)*M32)+a13*((202.2648888888889-12.736*kz)*M12-1.4364444444444444*mdT-4.309333333333333*meT-2.1546666666666665*mH1sq-2.1546666666666665*mH2sq-2.1546666666666665*mlT-0.7182222222222222*mqT-5.745777777777778*muT)+a1*a32*((51.91111111111111-11.733333333333333*kz)*M12+(80.35555555555555-23.466666666666665*kz)*M1*M3+(129.06666666666666-35.2*kz)*M32-2.8444444444444446*mqT-1.4222222222222223*(mdT+muT))+a33*((2279.1111111111113+1280.*kz)*M32+71.11111111111111*mqT+35.55555555555556*(mdT+muT))+a12*(-2.24*(d2mqT+d2tmdT+hd2T)-2.88*(e2mlT+e2tmeT+he2T)+(7.466666666666667*hdtdT+9.6*heteT+13.866666666666667*hutuT)*M1+a2*((21.6-4.32*kz)*M12+(14.4-2.88*kz)*M1*M2+(8.64-1.44*kz)*M22)+(-2.24*d2T-2.88*e2T)*mH1sq+a3*((38.82666666666667-14.08*kz)*M12+(25.884444444444444-9.386666666666667*kz)*M1*M3+(17.635555555555555-4.693333333333333*kz)*M32-0.5688888888888889*mdT-1.7066666666666668*meT-0.8533333333333334*mH1sq-0.8533333333333334*mH2sq-0.8533333333333334*mlT-0.28444444444444444*mqT-2.2755555555555556*muT)+M12*(-11.2*d2T-14.4*e2T-20.8*u2T)-4.16*mH2sq*u2T-4.16*(hu2T+u2mqT+u2tmuT))+a32*(213.33333333333334*(hdtdT+hutuT)*M3-64.*d2T*mH1sq-64.*mH2sq*u2T-320.*M32*(d2T+u2T)-64.*(d2mqT+d2tmdT+hd2T+hu2T+u2mqT+u2tmuT))+(-36.*d2T*d2T+72.*d4T-24.*d2T*e2T-4.*e2T*e2T+32.*a3*(e2T+d2T*(-1+kz))+a32*(10.666666666666666-60.44444444444444*kz)+a22*(-87.-3.*kz)+a12*(-22.466666666666665-0.03111111111111111*kz)+a1*(14.*d2T-6.*e2T-2.*d2T*kz+2.*e2T*kz+a2*(-17.2+2.8*kz)+a3*(-9.6+3.5555555555555554*kz))+a2*(d2T*(54.-18.*kz)+e2T*(18.-6.*kz)+a3*(-176.+32.*kz))+24.*(e4T+u2d2T))*hd2t+(144.*hdtd2dT-72.*d2T*hdtdT+48.*hete2eT-e2T*(24.*hdtdT+8.*heteT)+24.*(hdtu2dT-d2T*heteT+hutd2uT)+a12*(44.93333333333333+0.06222222222222222*kz)*M1+a22*(174.+6.*kz)*M2+a32*(-21.333333333333332+120.88888888888889*kz)*M3+a3*(32.*(heteT+hdtdT*(-1+kz))-32.*(e2T+d2T*(-1+kz))*M3)+a2*(hdtdT*(54.-18.*kz)+heteT*(18.-6.*kz)+(-54.*d2T-18.*e2T+18.*d2T*kz+6.*e2T*kz)*M2+a3*((176.-32.*kz)*M2+(176.-32.*kz)*M3))+a1*(14.*hdtdT-6.*heteT-2.*hdtdT*kz+2.*heteT*kz+(-14.*d2T+6.*e2T+2.*d2T*kz-2.*e2T*kz)*M1+a2*((17.2-2.8*kz)*M1+(17.2-2.8*kz)*M2)+a3*((9.6-3.5555555555555554*kz)*M1+(9.6-3.5555555555555554*kz)*M3)))*hdtd+(144.*hdtd2dT-72.*d2T*hdtdT+48.*hete2eT-e2T*(24.*hdtdT+8.*heteT)+24.*(hdtu2dT-d2T*heteT+hutd2uT)+a12*(44.93333333333333+0.06222222222222222*kz)*M1+a22*(174.+6.*kz)*M2+a32*(-21.333333333333332+120.88888888888889*kz)*M3+a3*(32.*(heteT+hdtdT*(-1+kz))-32.*(e2T+d2T*(-1+kz))*M3)+a2*(hdtdT*(54.-18.*kz)+heteT*(18.-6.*kz)+(-54.*d2T-18.*e2T+18.*d2T*kz+6.*e2T*kz)*M2+a3*((176.-32.*kz)*M2+(176.-32.*kz)*M3))+a1*(14.*hdtdT-6.*heteT-2.*hdtdT*kz+2.*heteT*kz+(-14.*d2T+6.*e2T+2.*d2T*kz-2.*e2T*kz)*M1+a2*((17.2-2.8*kz)*M1+(17.2-2.8*kz)*M2)+a3*((9.6-3.5555555555555554*kz)*M1+(9.6-3.5555555555555554*kz)*M3)))*dthd+(144.*(d2dmddtT+d4mqT)+48.*(e2e2mlT+e2emeetT)-e2mlT*(24.*d2T+8.*e2T)-72.*d2T*(d2mqT+d2tmdT+hd2T)-24.*e2T*(d2mqT+d2tmdT+hd2T)+288.*hd2td2tT-72.*hdtdT*hdtdT-8.*e2T*(e2tmeT+he2T)+96.*he2te2tT-48.*hdtdT*heteT-8.*heteT*heteT+a3*(32.*(e2mlT+e2tmeT+he2T+d2mqT*(-1+kz)+d2tmdT*(-1+kz)+hd2T*(-1+kz))-64.*(heteT+hdtdT*(-1+kz))*M3+64.*(e2T+d2T*(-1+kz))*M32+64.*(e2T+d2T*(-1+kz))*mH1sq)+a32*((64.-362.6666666666667*kz)*M32+(10.666666666666666-60.44444444444444*kz)*mH1sq)+a1*(-6.*e2mlT-6.*e2tmeT+14.*(d2mqT+d2tmdT+hd2T)-6.*he2T-2.*d2mqT*kz-2.*d2tmdT*kz+2.*e2mlT*kz+2.*e2tmeT*kz-2.*hd2T*kz+2.*he2T*kz+(-28.*hdtdT+12.*heteT+4.*hdtdT*kz-4.*heteT*kz)*M1+(28.*d2T-12.*e2T-4.*d2T*kz+4.*e2T*kz)*M12+(28.*d2T-12.*e2T-4.*d2T*kz+4.*e2T*kz)*mH1sq+a2*((-34.4+5.6*kz)*M12+(-34.4+5.6*kz)*M1*M2+(-34.4+5.6*kz)*M22+(-17.2+2.8*kz)*mH1sq)+a3*((-19.2+7.111111111111111*kz)*M12+(-19.2+7.111111111111111*kz)*M1*M3+(-19.2+7.111111111111111*kz)*M32+(-9.6+3.5555555555555554*kz)*mH1sq))+a2*(54.*(d2mqT+d2tmdT+hd2T)+18.*he2T+e2mlT*(18.-6.*kz)+e2tmeT*(18.-6.*kz)-18.*d2mqT*kz-18.*d2tmdT*kz-18.*hd2T*kz-6.*he2T*kz+(-108.*hdtdT-36.*heteT+36.*hdtdT*kz+12.*heteT*kz)*M2+(d2T*(108.-36.*kz)+e2T*(36.-12.*kz))*M22+(d2T*(108.-36.*kz)+e2T*(36.-12.*kz))*mH1sq+a3*((-352.+64.*kz)*M22+(-352.+64.*kz)*M2*M3+(-352.+64.*kz)*M32+(-176.+32.*kz)*mH1sq))+a22*((-474.-18.*kz)*M22+(-99.-3.*kz)*mH1sq-12.*mH2sq-12.*mlT-36.*mqT)+a12*((-134.8-0.18666666666666668*kz)*M12-0.32*mdT-0.96*meT+(-22.946666666666665-0.03111111111111111*kz)*mH1sq-0.48*mH2sq-0.48*mlT-0.16*mqT-1.28*muT)+24.*mH2sq*u2d2T+mH1sq*(-108.*d2T*d2T+216.*d4T-72.*d2T*e2T-12.*e2T*e2T+72.*e4T+48.*u2d2T)+24.*(d2u2mqT+d2umuutT+hd2u2T-d2T*(e2tmeT+he2T)+hu2d2T+u2d2mqT+u2dmddtT)+48.*uhuthddtT)*d2t+(-18.*d2T*d2T+36.*d4T-12.*d2T*e2T-2.*e2T*e2T+16.*a3*(e2T+d2T*(-1+kz))+a32*(5.333333333333333-30.22222222222222*kz)+a22*(-43.5-1.5*kz)+a12*(-11.233333333333333-0.015555555555555555*kz)+a1*(d2T*(7.-1.*kz)+e2T*(-3.+kz)+a2*(-8.6+1.4*kz)+a3*(-4.8+1.7777777777777777*kz))+a2*(d2T*(27.-9.*kz)+e2T*(9.-3.*kz)+a3*(-88.+16.*kz))+12.*(e4T+u2d2T))*md*d2t+(-36.*d2T*d2T+72.*d4T-24.*d2T*e2T-4.*e2T*e2T+32.*a3*(e2T+d2T*(-1+kz))+a32*(10.666666666666666-60.44444444444444*kz)+a22*(-87.-3.*kz)+a12*(-22.466666666666665-0.03111111111111111*kz)+a1*(14.*d2T-6.*e2T-2.*d2T*kz+2.*e2T*kz+a2*(-17.2+2.8*kz)+a3*(-9.6+3.5555555555555554*kz))+a2*(d2T*(54.-18.*kz)+e2T*(18.-6.*kz)+a3*(-176.+32.*kz))+24.*(e4T+u2d2T))*dtmqd+(-18.*d2T*d2T+36.*d4T-12.*d2T*e2T-2.*e2T*e2T+16.*a3*(e2T+d2T*(-1+kz))+a32*(5.333333333333333-30.22222222222222*kz)+a22*(-43.5-1.5*kz)+a12*(-11.233333333333333-0.015555555555555555*kz)+a1*(d2T*(7.-1.*kz)+e2T*(-3.+kz)+a2*(-8.6+1.4*kz)+a3*(-4.8+1.7777777777777777*kz))+a2*(d2T*(27.-9.*kz)+e2T*(9.-3.*kz)+a3*(-88.+16.*kz))+12.*(e4T+u2d2T))*d2t*md+(42.666666666666664*a3+12.*d2T+4.*e2T+a2*(18.-6.*kz)+a1*(-0.6666666666666666+0.4*kz))*hd2t*d2t+(42.666666666666664*a3-12.*d2T-4.*e2T+a2*(18.-6.*kz)+a1*(-3.8666666666666667+1.2*kz)+24.*u2T)*hdt*huut*d1+(42.666666666666664*a3+12.*d2T+4.*e2T+a2*(18.-6.*kz)+a1*(-0.6666666666666666+0.4*kz))*hdtd*dthd+(12.*hdtdT+4.*heteT+a1*(0.6666666666666666-0.4*kz)*M1+a2*(-18.+6.*kz)*M2-42.666666666666664*a3*M3)*hdtd*d2t+(42.666666666666664*a3-12.*d2T-4.*e2T+a2*(18.-6.*kz)+a1*(-3.8666666666666667+1.2*kz)+24.*u2T)*hdt*u2*hd+(-12.*hdtdT-4.*heteT+24.*hutuT+a1*(3.8666666666666667-1.2*kz)*M1+a2*(-18.+6.*kz)*M2-42.666666666666664*a3*M3)*hdt*u1*utd+(42.666666666666664*a3+12.*d2T+4.*e2T+a2*(18.-6.*kz)+a1*(-0.6666666666666666+0.4*kz))*dthd*hdtd+(12.*hdtdT+4.*heteT+a1*(0.6666666666666666-0.4*kz)*M1+a2*(-18.+6.*kz)*M2-42.666666666666664*a3*M3)*dthd*d2t+(42.666666666666664*a3-12.*d2T-4.*e2T+a2*(18.-6.*kz)+a1*(-3.8666666666666667+1.2*kz)+24.*u2T)*dt*hu2*d1+(-12.*hdtdT-4.*heteT+24.*hutuT+a1*(3.8666666666666667-1.2*kz)*M1+a2*(-18.+6.*kz)*M2-42.666666666666664*a3*M3)*dt*huut*d1+(42.666666666666664*a3+12.*d2T+4.*e2T+a2*(18.-6.*kz)+a1*(-0.6666666666666666+0.4*kz))*d2t*hd2t+(12.*hdtdT+4.*heteT+a1*(0.6666666666666666-0.4*kz)*M1+a2*(-18.+6.*kz)*M2-42.666666666666664*a3*M3)*d2t*hdtd+(12.*hdtdT+4.*heteT+a1*(0.6666666666666666-0.4*kz)*M1+a2*(-18.+6.*kz)*M2-42.666666666666664*a3*M3)*d2t*dthd+(12.*(d2mqT+d2tmdT+hd2T)+4.*(e2mlT+e2tmeT+he2T)+(36.*d2T+12.*e2T)*mH1sq+a3*(85.33333333333333*M32+85.33333333333333*mH1sq)+a2*((36.-12.*kz)*M22+(36.-12.*kz)*mH1sq)+a1*((-1.3333333333333333+0.8*kz)*M12+(-1.3333333333333333+0.8*kz)*mH1sq))*d4t+(42.666666666666664*a3-12.*d2T-4.*e2T+a2*(18.-6.*kz)+a1*(-3.8666666666666667+1.2*kz)+24.*u2T)*dt*uhut*hd+(-12.*hdtdT-4.*heteT+24.*hutuT+a1*(3.8666666666666667-1.2*kz)*M1+a2*(-18.+6.*kz)*M2-42.666666666666664*a3*M3)*dt*uhut*d1+(-12.*hdtdT-4.*heteT+24.*hutuT+a1*(3.8666666666666667-1.2*kz)*M1+a2*(-18.+6.*kz)*M2-42.666666666666664*a3*M3)*dtu*ut*hd+(-12.*(d2mqT+d2tmdT+hd2T)-4.*(e2mlT+e2tmeT+he2T)+a3*(85.33333333333333*M32+42.666666666666664*mH1sq+42.666666666666664*mH2sq)+a2*((36.-12.*kz)*M22+(18.-6.*kz)*mH1sq+(18.-6.*kz)*mH2sq)+a1*((-7.733333333333333+2.4*kz)*M12+(-3.8666666666666667+1.2*kz)*mH1sq+(-3.8666666666666667+1.2*kz)*mH2sq)+mH2sq*(-12.*d2T-4.*e2T+48.*u2T)+mH1sq*(-8.*e2T+24.*(-d2T+u2T))+24.*(hu2T+u2mqT+u2tmuT))*dtu2d+(21.333333333333332*a3+6.*d2T+2.*e2T+a2*(9.-3.*kz)+a1*(-0.3333333333333333+0.2*kz))*md*d4t+(21.333333333333332*a3-6.*d2T-2.*e2T+a2*(9.-3.*kz)+a1*(-1.9333333333333333+0.6*kz)+12.*u2T)*md*dtu2d+(42.666666666666664*a3+12.*d2T+4.*e2T+a2*(18.-6.*kz)+a1*(-0.6666666666666666+0.4*kz))*dtmqd*d2t+(42.666666666666664*a3-12.*d2T-4.*e2T+a2*(18.-6.*kz)+a1*(-3.8666666666666667+1.2*kz)+24.*u2T)*dtmqu*utd+(42.666666666666664*a3+12.*d2T+4.*e2T+a2*(18.-6.*kz)+a1*(-0.6666666666666666+0.4*kz))*d2t*md*d2t+(42.666666666666664*a3+12.*d2T+4.*e2T+a2*(18.-6.*kz)+a1*(-0.6666666666666666+0.4*kz))*d2t*dtmqd+(21.333333333333332*a3+6.*d2T+2.*e2T+a2*(9.-3.*kz)+a1*(-0.3333333333333333+0.2*kz))*d4t*md+(42.666666666666664*a3-12.*d2T-4.*e2T+a2*(18.-6.*kz)+a1*(-3.8666666666666667+1.2*kz)+24.*u2T)*dt*umuut*d1+(42.666666666666664*a3-12.*d2T-4.*e2T+a2*(18.-6.*kz)+a1*(-3.8666666666666667+1.2*kz)+24.*u2T)*dtu*utmqd+(21.333333333333332*a3-6.*d2T-2.*e2T+a2*(9.-3.*kz)+a1*(-1.9333333333333333+0.6*kz)+12.*u2T)*dtu2d*md+(12+4.*kz)*hd2td4t-4.*hd2t*dtu2d-4.*hdt*huut*d1*d2t+12.*hdt*huut*u1*utd+(12+4.*kz)*hdtd*dthd*d2t-4.*hdtd*dt*huut*d1+(12+4.*kz)*hdtd*d2t*dthd-4.*hdtd*dt*u2*hd-4.*hdt*u2*hd*d2t+12.*hdt*u1*uthu*utd-4.*hdt*u2d2*hd+12.*hdt*u4*hd+(12+4.*kz)*dthd*hdtd*d2t-4.*dthd*hdt*u1*utd+(12+4.*kz)*dthd*d2t*hdtd-4.*dthd*dt*uhut*d1-4.*dt*hu2*d1*d2t+12.*dt*hu*hutu*utd-4.*dt*huut*d1*hdtd+12.*dt*huut*uhut*d1+(12+4.*kz)*d2t*hd2t*d2t-4.*d2t*hdt*huut*d1+(12+4.*kz)*d2t*hdtd*dthd-4.*d2t*hdt*u2*hd+(12+4.*kz)*d2t*dthd*hdtd-4.*d2t*dt*hu2*d1+(12+4.*kz)*d4t*hd2t+(36.+12.*kz)*mH1sq*d6t-4.*d2t*dt*uhut*hd+(-8.*mH1sq-4.*mH2sq)*d2t*dtu2d-4.*dt*uhut*hd*d2t+12.*dt*uhu2tut*d1-4.*dt*uhut*d1*dthd+12.*dtu*hutu*ut*hd-4.*dtu*ut*hd*hdtd+12.*dtu*uthu*hut*d1-4.*dtu2d*hd2t+(-8.*mH1sq-4.*mH2sq)*dtu2d*d2t+12.*dtu*ut*uhut*hd+(12.*mH1sq+24.*mH2sq)*dt*u4*d1+(6+2.*kz)*md*d6t-2.*md*d2t*dtu2d-2.*md*dtu2d*d2t+6.*md*dt*u4*d1+(12+4.*kz)*dtmqd*d4t-4.*dtmqd*dtu2d-4.*dtmqu*ut*d1*d2t+12.*dtmqu*u2t*utd+(12+4.*kz)*d2t*md*d4t-4.*d2t*md*dtu2d+(12+4.*kz)*d2t*dtmqd*d2t-4.*d2t*dtmqu*utd+(12+4.*kz)*d4t*md*d2t+(12+4.*kz)*d4t*dtmqd+(6+2.*kz)*d6t*md-4.*d2t*dt*umuut*d1-4.*d2t*dt*u1*utmqd-2.*d2t*dtu2d*md-4.*dt*umuut*d1*d2t+12.*dt*umuut*u1*utd-4.*dtu*utmqd*d2t+12.*dtu*utmqu*utd-4.*dtu2d*md*d2t-4.*dtu2d*dtmqd-2.*dtu2d*d2t*md+12.*dtu*ut*umuut*d1+12.*dtu*u2t*utmqd+6.*dt*u4*d1*md;
	
	dmu3=a2*a32*((288.-48.*kz)*M22+(480.-96.*kz)*M2*M3+(720.-144.*kz)*M32)+a13*((768.4408888888889-50.944*kz)*M12-6.087111111111111*mdT-18.261333333333333*meT-9.130666666666666*mH1sq-9.130666666666666*mH2sq-9.130666666666666*mlT-3.0435555555555553*mqT-24.348444444444443*muT)+a1*a32*((-3.5555555555555554-11.733333333333333*kz)*M12+(-30.57777777777778-23.466666666666665*kz)*M1*M3+(-11.733333333333333-35.2*kz)*M32-11.377777777777778*mqT-5.688888888888889*(mdT+muT))+a33*((2279.1111111111113+1280.*kz)*M32+71.11111111111111*mqT+35.55555555555556*(mdT+muT))+a12*(-8.96*(d2mqT+d2tmdT+hd2T)-11.52*(e2mlT+e2tmeT+he2T)+(29.866666666666667*hdtdT+38.4*heteT+55.46666666666667*hutuT)*M1+a2*((86.4-17.28*kz)*M12+(57.6-11.52*kz)*M1*M2+(34.56-5.76*kz)*M22)+(-8.96*d2T-11.52*e2T)*mH1sq+a3*((114.34666666666666-56.32*kz)*M12+(76.2311111111111-37.54666666666667*kz)*M1*M3+(56.888888888888886-18.773333333333333*kz)*M32-2.2755555555555556*mdT-6.826666666666667*meT-3.4133333333333336*mH1sq-3.4133333333333336*mH2sq-3.4133333333333336*mlT-1.1377777777777778*mqT-9.102222222222222*muT)+M12*(-44.8*d2T-57.6*e2T-83.2*u2T)-16.64*mH2sq*u2T-16.64*(hu2T+u2mqT+u2tmuT))+a32*(213.33333333333334*(hdtdT+hutuT)*M3-64.*d2T*mH1sq-64.*mH2sq*u2T-320.*M32*(d2T+u2T)-64.*(d2mqT+d2tmdT+hd2T+hu2T+u2mqT+u2tmuT))+(a32*(10.666666666666666-60.44444444444444*kz)+a22*(-87.-3.*kz)+a12*(-31.96-1.0977777777777777*kz)+24.*u2d2T+32.*a3*(-1+kz)*u2T-36.*u2T*u2T+a2*(a3*(-176.+32.*kz)+(54.-18.*kz)*u2T)+a1*(a3*(-1.0666666666666667-4.977777777777778*kz)+a2*(-26.8+5.2*kz)+(14.+2.8*kz)*u2T)+72.*u4T)*hu2t+(24.*(hdtu2dT+hutd2uT)+144.*hutu2uT+a12*(63.92+2.1955555555555555*kz)*M1+a22*(174.+6.*kz)*M2+a32*(-21.333333333333332+120.88888888888889*kz)*M3-72.*hutuT*u2T+a1*(hutuT*(14.+2.8*kz)+a2*((26.8-5.2*kz)*M1+(26.8-5.2*kz)*M2)+a3*((1.0666666666666667+4.977777777777778*kz)*M1+(1.0666666666666667+4.977777777777778*kz)*M3)-(14.+2.8*kz)*M1*u2T)+a2*(hutuT*(54.-18.*kz)+a3*((176.-32.*kz)*M2+(176.-32.*kz)*M3)+(-54.+18.*kz)*M2*u2T)+a3*(32.*hutuT*(-1+kz)-32.*(-1+kz)*M3*u2T))*hutu+(24.*(hdtu2dT+hutd2uT)+144.*hutu2uT+a12*(63.92+2.1955555555555555*kz)*M1+a22*(174.+6.*kz)*M2+a32*(-21.333333333333332+120.88888888888889*kz)*M3-72.*hutuT*u2T+a1*(hutuT*(14.+2.8*kz)+a2*((26.8-5.2*kz)*M1+(26.8-5.2*kz)*M2)+a3*((1.0666666666666667+4.977777777777778*kz)*M1+(1.0666666666666667+4.977777777777778*kz)*M3)-(14.+2.8*kz)*M1*u2T)+a2*(hutuT*(54.-18.*kz)+a3*((176.-32.*kz)*M2+(176.-32.*kz)*M3)+(-54.+18.*kz)*M2*u2T)+a3*(32.*hutuT*(-1+kz)-32.*(-1+kz)*M3*u2T))*uthu+(288.*hu2tu2tT+a32*((64.-362.6666666666667*kz)*M32+(10.666666666666666-60.44444444444444*kz)*mH2sq)+a22*((-474.-18.*kz)*M22-12.*mH1sq+(-99.-3.*kz)*mH2sq-12.*mlT-36.*mqT)+a12*((-191.76-6.586666666666667*kz)*M12+0.32*mdT+0.96*meT+0.48*mH1sq+(-31.48-1.0977777777777777*kz)*mH2sq+0.48*mlT+0.16*mqT+1.28*muT)+24.*mH1sq*u2d2T+24.*(d2u2mqT+d2umuutT+hd2u2T+hu2d2T+u2d2mqT+u2dmddtT)+a3*(-64.*hutuT*(-1+kz)*M3+64.*(-1+kz)*M32*u2T+64.*(-1+kz)*mH2sq*u2T+32.*(-1+kz)*(hu2T+u2mqT+u2tmuT))+a1*(-(hutuT*(28.+5.6*kz)*M1)+a3*((-2.1333333333333333-9.955555555555556*kz)*M12+(-2.1333333333333333-9.955555555555556*kz)*M1*M3+(-2.1333333333333333-9.955555555555556*kz)*M32+(-1.0666666666666667-4.977777777777778*kz)*mH2sq)+a2*((-53.6+10.4*kz)*M12+(-53.6+10.4*kz)*M1*M2+(-53.6+10.4*kz)*M22+(-26.8+5.2*kz)*mH2sq)+(28.+5.6*kz)*M12*u2T+(28.+5.6*kz)*mH2sq*u2T+(14.+2.8*kz)*(hu2T+u2mqT+u2tmuT))+a2*(hutuT*(-108.+36.*kz)*M2+a3*((-352.+64.*kz)*M22+(-352.+64.*kz)*M2*M3+(-352.+64.*kz)*M32+(-176.+32.*kz)*mH2sq)+(108.-36.*kz)*M22*u2T+(108.-36.*kz)*mH2sq*u2T-(-54.+18.*kz)*(hu2T+u2mqT+u2tmuT))-72.*(hutuT*hutuT+u2T*(hu2T+u2mqT+u2tmuT))+144.*(u2umuutT+u4mqT)+mH2sq*(48.*u2d2T-108.*u2T*u2T+216.*u4T)+48.*uhuthddtT)*u2t+(a32*(5.333333333333333-30.22222222222222*kz)+a22*(-43.5-1.5*kz)+a12*(-15.98-0.5488888888888889*kz)+12.*u2d2T+16.*a3*(-1+kz)*u2T-18.*u2T*u2T+a2*(a3*(-88.+16.*kz)+(27.-9.*kz)*u2T)+a1*(a3*(-0.5333333333333333-2.488888888888889*kz)+a2*(-13.4+2.6*kz)+(7.+1.4*kz)*u2T)+36.*u4T)*mu*u2t+(a32*(10.666666666666666-60.44444444444444*kz)+a22*(-87.-3.*kz)+a12*(-31.96-1.0977777777777777*kz)+24.*u2d2T+32.*a3*(-1+kz)*u2T-36.*u2T*u2T+a2*(a3*(-176.+32.*kz)+(54.-18.*kz)*u2T)+a1*(a3*(-1.0666666666666667-4.977777777777778*kz)+a2*(-26.8+5.2*kz)+(14.+2.8*kz)*u2T)+72.*u4T)*utmqu+(a32*(5.333333333333333-30.22222222222222*kz)+a22*(-43.5-1.5*kz)+a12*(-15.98-0.5488888888888889*kz)+12.*u2d2T+16.*a3*(-1+kz)*u2T-18.*u2T*u2T+a2*(a3*(-88.+16.*kz)+(27.-9.*kz)*u2T)+a1*(a3*(-0.5333333333333333-2.488888888888889*kz)+a2*(-13.4+2.6*kz)+(7.+1.4*kz)*u2T)+36.*u4T)*u2t*mu+(42.666666666666664*a3+24.*d2T+8.*e2T+a2*(18.-6.*kz)+a1*(2.533333333333333+1.2*kz)-12.*u2T)*hut*hddt*u1+(42.666666666666664*a3+a2*(18.-6.*kz)+a1*(-0.6666666666666666+2.*kz)+12.*u2T)*hu2t*u2t+(42.666666666666664*a3+24.*d2T+8.*e2T+a2*(18.-6.*kz)+a1*(2.533333333333333+1.2*kz)-12.*u2T)*hut*d2*hu+(24.*hdtdT+8.*heteT-12.*hutuT+a1*(-2.533333333333333-1.2*kz)*M1+a2*(-18.+6.*kz)*M2-42.666666666666664*a3*M3)*hut*d1*dtu+(42.666666666666664*a3+a2*(18.-6.*kz)+a1*(-0.6666666666666666+2.*kz)+12.*u2T)*hutu*uthu+(12.*hutuT+a1*(0.6666666666666666-2.*kz)*M1+a2*(-18.+6.*kz)*M2-42.666666666666664*a3*M3)*hutu2u+(42.666666666666664*a3+24.*d2T+8.*e2T+a2*(18.-6.*kz)+a1*(2.533333333333333+1.2*kz)-12.*u2T)*ut*hd2*u1+(24.*hdtdT+8.*heteT-12.*hutuT+a1*(-2.533333333333333-1.2*kz)*M1+a2*(-18.+6.*kz)*M2-42.666666666666664*a3*M3)*ut*hddt*u1+(42.666666666666664*a3+a2*(18.-6.*kz)+a1*(-0.6666666666666666+2.*kz)+12.*u2T)*uthu*hutu+(12.*hutuT+a1*(0.6666666666666666-2.*kz)*M1+a2*(-18.+6.*kz)*M2-42.666666666666664*a3*M3)*uthu*u2t+(42.666666666666664*a3+24.*d2T+8.*e2T+a2*(18.-6.*kz)+a1*(2.533333333333333+1.2*kz)-12.*u2T)*ut*dhdt*hu+(24.*hdtdT+8.*heteT-12.*hutuT+a1*(-2.533333333333333-1.2*kz)*M1+a2*(-18.+6.*kz)*M2-42.666666666666664*a3*M3)*ut*dhdt*u1+(24.*hdtdT+8.*heteT-12.*hutuT+a1*(-2.533333333333333-1.2*kz)*M1+a2*(-18.+6.*kz)*M2-42.666666666666664*a3*M3)*utd*dt*hu+(24.*(d2mqT+d2tmdT+hd2T)+8.*(e2mlT+e2tmeT+he2T)+a3*(85.33333333333333*M32+42.666666666666664*mH1sq+42.666666666666664*mH2sq)+a2*((36.-12.*kz)*M22+(18.-6.*kz)*mH1sq+(18.-6.*kz)*mH2sq)+a1*((5.066666666666666+2.4*kz)*M12+(2.533333333333333+1.2*kz)*mH1sq+(2.533333333333333+1.2*kz)*mH2sq)+mH2sq*(8.*e2T+24.*(d2T-u2T))+mH1sq*(48.*d2T+16.*e2T-12.*u2T)-12.*(hu2T+u2mqT+u2tmuT))*utd2u+(42.666666666666664*a3+a2*(18.-6.*kz)+a1*(-0.6666666666666666+2.*kz)+12.*u2T)*u2t*hu2t+(12.*hutuT+a1*(0.6666666666666666-2.*kz)*M1+a2*(-18.+6.*kz)*M2-42.666666666666664*a3*M3)*u2t*hutu+(12.*hutuT+a1*(0.6666666666666666-2.*kz)*M1+a2*(-18.+6.*kz)*M2-42.666666666666664*a3*M3)*u2t*uthu+(a3*(85.33333333333333*M32+85.33333333333333*mH2sq)+a2*((36.-12.*kz)*M22+(36.-12.*kz)*mH2sq)+a1*((-1.3333333333333333+4.*kz)*M12+(-1.3333333333333333+4.*kz)*mH2sq)+36.*mH2sq*u2T+12.*(hu2T+u2mqT+u2tmuT))*u4t+(21.333333333333332*a3+12.*d2T+4.*e2T+a2*(9.-3.*kz)+a1*(1.2666666666666666+0.6*kz)-6.*u2T)*mu*utd2u+(21.333333333333332*a3+a2*(9.-3.*kz)+a1*(-0.3333333333333333+kz)+6.*u2T)*mu*u4t+(42.666666666666664*a3+24.*d2T+8.*e2T+a2*(18.-6.*kz)+a1*(2.533333333333333+1.2*kz)-12.*u2T)*utmqd*dtu+(42.666666666666664*a3+a2*(18.-6.*kz)+a1*(-0.6666666666666666+2.*kz)+12.*u2T)*utmqu*u2t+(42.666666666666664*a3+24.*d2T+8.*e2T+a2*(18.-6.*kz)+a1*(2.533333333333333+1.2*kz)-12.*u2T)*ut*dmddt*u1+(42.666666666666664*a3+24.*d2T+8.*e2T+a2*(18.-6.*kz)+a1*(2.533333333333333+1.2*kz)-12.*u2T)*utd*dtmqu+(21.333333333333332*a3+12.*d2T+4.*e2T+a2*(9.-3.*kz)+a1*(1.2666666666666666+0.6*kz)-6.*u2T)*utd2u*mu+(42.666666666666664*a3+a2*(18.-6.*kz)+a1*(-0.6666666666666666+2.*kz)+12.*u2T)*u2t*mu*u2t+(42.666666666666664*a3+a2*(18.-6.*kz)+a1*(-0.6666666666666666+2.*kz)+12.*u2T)*u2t*utmqu+(21.333333333333332*a3+a2*(9.-3.*kz)+a1*(-0.3333333333333333+kz)+6.*u2T)*u4t*mu+12.*hut*hddt*d1*dtu-4.*hut*hddt*u1*u2t-4.*hu2t*utd2u+(12+4.*kz)*hu2tu4t+12.*hut*d1*dthd*dtu-4.*hut*d2*hu*u2t+12.*hut*d4*hu-4.*hut*d2*u1*uthu-4.*hutu*ut*hddt*u1+(12+4.*kz)*hutu*uthu*u2t-4.*hutu*ut*d2*hu+(12+4.*kz)*hutu2u*uthu+12.*ut*hd*hdtd*dtu-4.*ut*hd2*u1*u2t+12.*ut*hddt*dhdt*u1-4.*ut*hddt*u1*hutu-4.*uthu*hut*d1*dtu+(12+4.*kz)*uthu*hutu2u-4.*uthu*ut*dhdt*u1+(12+4.*kz)*uthu*u2t*hutu+12.*ut*dhd2tdt*u1-4.*ut*dhdt*hu*u2t+12.*utd*hdtd*dt*hu-4.*ut*dhdt*u1*uthu+12.*utd*dthd*hdt*u1-4.*utd*dt*hu*hutu+12.*utd*dt*dhdt*hu+(24.*mH1sq+12.*mH2sq)*ut*d4*u1-4.*utd2u*hu2t+(-4.*mH1sq-8.*mH2sq)*utd2u*u2t-4.*u2t*hut*hddt*u1+(12+4.*kz)*u2t*hu2t*u2t-4.*u2t*hut*d2*hu+(12+4.*kz)*u2t*hutu*uthu-4.*u2t*ut*hd2*u1+(12+4.*kz)*u2t*uthu*hutu-4.*u2t*ut*dhdt*hu+(-4.*mH1sq-8.*mH2sq)*u2t*utd2u+(12+4.*kz)*u4t*hu2t+(36.+12.*kz)*mH2sq*u6t+6.*mu*ut*d4*u1-2.*mu*utd2u*u2t-2.*mu*u2t*utd2u+(6+2.*kz)*mu*u6t+12.*utmqd*d2t*dtu-4.*utmqd*dt*u1*u2t-4.*utmqu*utd2u+(12+4.*kz)*utmqu*u4t+12.*ut*dmddt*d1*dtu-4.*ut*dmddt*u1*u2t+12.*utd*dtmqd*dtu-4.*utd*dtmqu*u2t+12.*utd*dt*dmddt*u1+12.*utd*d2t*dtmqu+6.*ut*d4*u1*mu-4.*utd2u*mu*u2t-4.*utd2u*utmqu-2.*utd2u*u2t*mu-4.*u2t*mu*utd2u+(12+4.*kz)*u2t*mu*u4t-4.*u2t*utmqd*dtu+(12+4.*kz)*u2t*utmqu*u2t-4.*u2t*ut*dmddt*u1-4.*u2t*ut*d1*dtmqu-2.*u2t*utd2u*mu+(12+4.*kz)*u4t*mu*u2t+(12+4.*kz)*u4t*utmqu+(6+2.*kz)*u6t*mu;
	
	dme3=a13*((1535.712-114.624*kz)*M12-14.976*mdT-44.928*meT-22.464*mH1sq-22.464*mH2sq-22.464*mlT-7.488*mqT-59.904*muT)+a12*(-20.16*(d2mqT+d2tmdT+hd2T)-25.92*(e2mlT+e2tmeT+he2T)+(67.2*hdtdT+86.4*heteT+124.8*hutuT)*M1+a2*((194.4-38.88*kz)*M12+(129.6-25.92*kz)*M1*M2+(77.76-12.96*kz)*M22)+a3*((633.6-126.72*kz)*M12+(422.4-84.48*kz)*M1*M3+(253.44-42.24*kz)*M32)+(-20.16*d2T-25.92*e2T)*mH1sq+M12*(-100.8*d2T-129.6*e2T-187.2*u2T)-37.44*mH2sq*u2T-37.44*(hu2T+u2mqT+u2tmuT))+(-36.*d2T*d2T+72.*d4T-24.*d2T*e2T-4.*e2T*e2T+a2*(d2T*(54.-18.*kz)+e2T*(18.-6.*kz))+a12*(-60.12-5.4*kz)+a22*(-87.-3.*kz)+a3*d2T*(-128.+32.*kz)+a1*(3.6*e2T*(1+kz)+d2T*(42.8+2.8*kz)+a2*(-54.+10.8*kz))+24.*(e4T+u2d2T))*he2t+(144.*hdtd2dT-72.*d2T*hdtdT+48.*hete2eT-e2T*(24.*hdtdT+8.*heteT)+24.*(hdtu2dT-d2T*heteT+hutd2uT)+a12*(120.24+10.8*kz)*M1+a22*(174.+6.*kz)*M2+a2*(hdtdT*(54.-18.*kz)+heteT*(18.-6.*kz)+(-54.*d2T-18.*e2T+18.*d2T*kz+6.*e2T*kz)*M2)+a1*(3.6*heteT*(1+kz)+hdtdT*(42.8+2.8*kz)+(-3.6*e2T*(1+kz)-d2T*(42.8+2.8*kz))*M1+a2*((54.-10.8*kz)*M1+(54.-10.8*kz)*M2))+a3*(hdtdT*(-128.+32.*kz)+d2T*(128.-32.*kz)*M3))*hete+(144.*hdtd2dT-72.*d2T*hdtdT+48.*hete2eT-e2T*(24.*hdtdT+8.*heteT)+24.*(hdtu2dT-d2T*heteT+hutd2uT)+a12*(120.24+10.8*kz)*M1+a22*(174.+6.*kz)*M2+a2*(hdtdT*(54.-18.*kz)+heteT*(18.-6.*kz)+(-54.*d2T-18.*e2T+18.*d2T*kz+6.*e2T*kz)*M2)+a1*(3.6*heteT*(1+kz)+hdtdT*(42.8+2.8*kz)+(-3.6*e2T*(1+kz)-d2T*(42.8+2.8*kz))*M1+a2*((54.-10.8*kz)*M1+(54.-10.8*kz)*M2))+a3*(hdtdT*(-128.+32.*kz)+d2T*(128.-32.*kz)*M3))*ethe+(144.*(d2dmddtT+d4mqT)+48.*(e2e2mlT+e2emeetT)-e2mlT*(24.*d2T+8.*e2T)-72.*d2T*(d2mqT+d2tmdT+hd2T)-24.*e2T*(d2mqT+d2tmdT+hd2T)+288.*hd2td2tT-72.*hdtdT*hdtdT-8.*e2T*(e2tmeT+he2T)+96.*he2te2tT-48.*hdtdT*heteT-8.*heteT*heteT+a2*(54.*(d2mqT+d2tmdT+hd2T)+18.*he2T+e2mlT*(18.-6.*kz)+e2tmeT*(18.-6.*kz)-18.*d2mqT*kz-18.*d2tmdT*kz-18.*hd2T*kz-6.*he2T*kz+(-108.*hdtdT-36.*heteT+36.*hdtdT*kz+12.*heteT*kz)*M2+(d2T*(108.-36.*kz)+e2T*(36.-12.*kz))*M22+(d2T*(108.-36.*kz)+e2T*(36.-12.*kz))*mH1sq)+a3*((d2mqT+d2tmdT+hd2T)*(-128.+32.*kz)+hdtdT*(256.-64.*kz)*M3+d2T*(-256.+64.*kz)*M32+d2T*(-256.+64.*kz)*mH1sq)+a1*(42.8*(d2mqT+d2tmdT+hd2T)+3.6*(e2mlT+e2tmeT+he2T)+2.8*(d2mqT+d2tmdT+hd2T)*kz+3.6*(e2mlT+e2tmeT+he2T)*kz+(-7.2*heteT*(1+kz)-hdtdT*(85.6+5.6*kz))*M1+(7.2*e2T*(1+kz)+d2T*(85.6+5.6*kz))*M12+(7.2*e2T*(1+kz)+d2T*(85.6+5.6*kz))*mH1sq+a2*((-108.+21.6*kz)*M12+(-108.+21.6*kz)*M1*M2+(-108.+21.6*kz)*M22+(-54.+10.8*kz)*mH1sq))+a22*((-474.-18.*kz)*M22+(-99.-3.*kz)*mH1sq-12.*mH2sq-12.*mlT-36.*mqT)+a12*((-360.72-32.4*kz)*M12+0.96*mdT+2.88*meT+(-58.68-5.4*kz)*mH1sq+1.44*mH2sq+1.44*mlT+0.48*mqT+3.84*muT)+24.*mH2sq*u2d2T+mH1sq*(-108.*d2T*d2T+216.*d4T-72.*d2T*e2T-12.*e2T*e2T+72.*e4T+48.*u2d2T)+24.*(d2u2mqT+d2umuutT+hd2u2T-d2T*(e2tmeT+he2T)+hu2d2T+u2d2mqT+u2dmddtT)+48.*uhuthddtT)*e2t+(-18.*d2T*d2T+36.*d4T-12.*d2T*e2T-2.*e2T*e2T+a2*(d2T*(27.-9.*kz)+e2T*(9.-3.*kz))+a12*(-30.06-2.7*kz)+a22*(-43.5-1.5*kz)+a3*d2T*(-64.+16.*kz)+a1*(1.8*e2T*(1+kz)+d2T*(21.4+1.4*kz)+a2*(-27.+5.4*kz))+12.*(e4T+u2d2T))*me*e2t+(-36.*d2T*d2T+72.*d4T-24.*d2T*e2T-4.*e2T*e2T+a2*(d2T*(54.-18.*kz)+e2T*(18.-6.*kz))+a12*(-60.12-5.4*kz)+a22*(-87.-3.*kz)+a3*d2T*(-128.+32.*kz)+a1*(3.6*e2T*(1+kz)+d2T*(42.8+2.8*kz)+a2*(-54.+10.8*kz))+24.*(e4T+u2d2T))*etmle+(-18.*d2T*d2T+36.*d4T-12.*d2T*e2T-2.*e2T*e2T+a2*(d2T*(27.-9.*kz)+e2T*(9.-3.*kz))+a12*(-30.06-2.7*kz)+a22*(-43.5-1.5*kz)+a3*d2T*(-64.+16.*kz)+a1*(1.8*e2T*(1+kz)+d2T*(21.4+1.4*kz)+a2*(-27.+5.4*kz))+12.*(e4T+u2d2T))*e2t*me+(12.*d2T+4.*e2T+a2*(18.-6.*kz)+3.6*a1*(1+kz))*he2t*e2t+(12.*d2T+4.*e2T+a2*(18.-6.*kz)+3.6*a1*(1+kz))*hete*ethe+(12.*hdtdT+4.*heteT-3.6*a1*(1+kz)*M1+a2*(-18.+6.*kz)*M2)*hete*e2t+(12.*d2T+4.*e2T+a2*(18.-6.*kz)+3.6*a1*(1+kz))*ethe*hete+(12.*hdtdT+4.*heteT-3.6*a1*(1+kz)*M1+a2*(-18.+6.*kz)*M2)*ethe*e2t+(12.*d2T+4.*e2T+a2*(18.-6.*kz)+3.6*a1*(1+kz))*e2t*he2t+(12.*hdtdT+4.*heteT-3.6*a1*(1+kz)*M1+a2*(-18.+6.*kz)*M2)*e2t*hete+(12.*hdtdT+4.*heteT-3.6*a1*(1+kz)*M1+a2*(-18.+6.*kz)*M2)*e2t*ethe+(12.*(d2mqT+d2tmdT+hd2T)+4.*(e2mlT+e2tmeT+he2T)+(36.*d2T+12.*e2T)*mH1sq+a2*((36.-12.*kz)*M22+(36.-12.*kz)*mH1sq)+a1*(7.2*(1+kz)*M12+7.2*(1+kz)*mH1sq))*e4t+(6.*d2T+2.*e2T+a2*(9.-3.*kz)+1.8*a1*(1+kz))*me*e4t+(12.*d2T+4.*e2T+a2*(18.-6.*kz)+3.6*a1*(1+kz))*etmle*e2t+(12.*d2T+4.*e2T+a2*(18.-6.*kz)+3.6*a1*(1+kz))*e2t*me*e2t+(12.*d2T+4.*e2T+a2*(18.-6.*kz)+3.6*a1*(1+kz))*e2t*etmle+(6.*d2T+2.*e2T+a2*(9.-3.*kz)+1.8*a1*(1+kz))*e4t*me+(12+4.*kz)*he2te4t+(12+4.*kz)*hete*ethe*e2t+(12+4.*kz)*hete*e2t*ethe+(12+4.*kz)*ethe*hete*e2t+(12+4.*kz)*ethe*et*ehet*e1+(12+4.*kz)*e2t*he2t*e2t+(12+4.*kz)*e2t*hete*ethe+(12+4.*kz)*e2t*ethe*hete+(12+4.*kz)*e4t*he2t+(36.+12.*kz)*mH1sq*e6t+(6+2.*kz)*me*e6t+(12+4.*kz)*etmle*e4t+(12+4.*kz)*e2t*me*e4t+(12+4.*kz)*e2t*etmle*e2t+(12+4.*kz)*e4t*me*e2t+(12+4.*kz)*e4t*etmle+(6+2.*kz)*e6tme;

dm3sq3=m3sq*(18.*d4T*e2T+e6T*(1+kz)+6.*e2T*(e4T+u2d2T)+18.*(d2T*(e4T+u2d2T)+u2d2T*u2T)+9.*(d4u2T+u4d2T)+54.*(d2T*d4T+u2T*u4T)+3.*(1+kz)*(d6T+u6T))+(108.*d4T*hdtdT+36.*e4T*hdtdT+36.*d2T*hdtu2dT+72.*d2T*hete2eT+36.*d4T*heteT+12.*e4T*heteT+36.*d2T*hutd2uT+e2T*(72.*hdtd2dT+24.*hete2eT+12.*(hdtu2dT+hutd2uT))+6.*hete4eT*(1+kz)+18.*(hdtd2u2dT+hdtd4dT+hdtu2d2dT+hdtu4dT+hutd2u2uT+hutd4uT+hutu2d2uT+hutu4uT+(hdtd4dT+hutu4uT)*kz)+36.*hdtdT*u2d2T+12.*heteT*u2d2T+36.*hutuT*u2d2T+36.*hdtu2dT*u2T+36.*hutd2uT*u2T+216.*(d2T*hdtd2dT+hutu2uT*u2T)+108.*hutuT*u4T)*(smu_)+a13*((1839./50-2.388*kz)*m3sq+(-220.68+14.328*kz)*M1*(smu_))+a23*((345./2+52.5*kz)*m3sq-(1035.+315.*kz)*M2*(smu_))+a32*(-((53.333333333333336+2.6666666666666665*kz)*m3sq*(d2T+u2T))-(hdtdT+hutuT)*(106.66666666666667+5.333333333333333*kz)*(smu_)+(213.33333333333334+10.666666666666666*kz)*M3*(d2T+u2T)*(smu_))+a3*(m3sq*(48.*u2d2T+72.*(d4T+u4T)-kz*(16.*u2d2T+24.*(d4T+u4T)))+(hdtu2dT*(96.-32.*kz)+hutd2uT*(96.-32.*kz)-(hdtd2dT+hutu2uT)*(-288.+96.*kz))*(smu_)+M3*(-96.*u2d2T+32.*kz*u2d2T-144.*(d4T+u4T)+48.*kz*(d4T+u4T))*(smu_))+a12*(-(m3sq*(16.683333333333334*d2T+32.85*e2T+kz*(0.25666666666666665*d2T-0.27*e2T+0.21666666666666667*u2T)+39.28333333333333*u2T))-(33.36666666666667*hdtdT+65.7*heteT+hutuT*(78.56666666666666+0.43333333333333335*kz)+0.5133333333333333*hdtdT*kz-0.54*heteT*kz)*(smu_)+M1*(66.73333333333333*d2T+131.4*e2T+kz*(1.0266666666666666*d2T-1.08*e2T+0.8666666666666667*u2T)+157.13333333333333*u2T)*(smu_)+a2*((27./50-1.62*kz)*m3sq+(-2.16+6.48*kz)*M1*(smu_)+(-1.08+3.24*kz)*M2*(smu_))+a3*((132./5-5.28*kz)*m3sq+(-105.6+21.12*kz)*M1*(smu_)+(-52.8+10.56*kz)*M3*(smu_)))+a22*(-(m3sq*(e2T*(26.25+5.25*kz)+(78.75+15.75*kz)*(d2T+u2T)))-(heteT*(52.5+10.5*kz)+(hdtdT+hutuT)*(157.5+31.5*kz))*(smu_)+M2*(e2T*(105.+21.*kz)+(315.+63.*kz)*(d2T+u2T))*(smu_)+a3*((180-36.*kz)*m3sq+(-720.+144.*kz)*M2*(smu_)+(-360.+72.*kz)*M3*(smu_)))+a2*(36.*m3sq*u2d2T+(1+kz)*m3sq*(3.*e4T+9.*(d4T+u4T))+(72.*(hdtu2dT+hutd2uT)+(12.*hete2eT+36.*(hdtd2dT+hutu2uT))*(1+kz))*(smu_)-M2*(72.*u2d2T+(1+kz)*(6.*e4T+18.*(d4T+u4T)))*(smu_)+a3*((-132.+24.*kz)*m3sq*(d2T+u2T)+(hdtdT+hutuT)*(-264.+48.*kz)*(smu_)-(-264.+48.*kz)*M2*(d2T+u2T)*(smu_)-(-264.+48.*kz)*M3*(d2T+u2T)*(smu_)))+a1*(m3sq*(3.*d4T+9.*e4T-1.2*u2d2T+kz*(1.8*(d4T-e4T)+1.6*u2d2T-0.6*u4T)+11.4*u4T)+(12.*hdtd2dT-2.4*hdtu2dT+36.*hete2eT-2.4*hutd2uT+45.6*hutu2uT+7.2*hdtd2dT*kz+3.2*hdtu2dT*kz-7.2*hete2eT*kz+3.2*hutd2uT*kz-2.4*hutu2uT*kz)*(smu_)+M1*(-6.*d4T-18.*e4T+2.4*u2d2T-22.8*u4T+kz*(3.6*(-d4T+e4T)-3.2*u2d2T+1.2*u4T))*(smu_)+a22*((9./2-2.7*kz)*m3sq+(-9.+5.4*kz)*M1*(smu_)+(-18.+10.8*kz)*M2*(smu_))+a3*(m3sq*(-18.933333333333334*d2T+3.7333333333333334*d2T*kz-41.333333333333336*u2T+6.933333333333334*kz*u2T)+(-37.86666666666667*hdtdT-82.66666666666667*hutuT+7.466666666666667*hdtdT*kz+13.866666666666667*hutuT*kz)*(smu_)+M1*(d2T*(37.86666666666667-7.466666666666667*kz)+(82.66666666666667-13.866666666666667*kz)*u2T)*(smu_)+M3*(d2T*(37.86666666666667-7.466666666666667*kz)+(82.66666666666667-13.866666666666667*kz)*u2T)*(smu_))+a2*(-(m3sq*(8.1*e2T+d2T*(0.3+1.5*kz)+5.7*u2T-kz*(2.7*e2T+2.1*u2T)))-(heteT*(16.2-5.4*kz)+hutuT*(11.4-4.2*kz)+hdtdT*(0.6+3.*kz))*(smu_)+M1*(16.2*e2T+d2T*(0.6+3.*kz)+11.4*u2T-kz*(5.4*e2T+4.2*u2T))*(smu_)+M2*(16.2*e2T+d2T*(0.6+3.*kz)+11.4*u2T-kz*(5.4*e2T+4.2*u2T))*(smu_)));
 
 dmG(1) = dmG(1) + dM13 * threelp;
 dmG(2) = dmG(2) + dM23 * threelp;
 dmG(3) = dmG(3) + dM33 * threelp;
 dm3sq = dm3sq + dm3sq3 * threelp;
 dmH1sq = dmH1sq + dmH1sq3 * threelp;
 dmH2sq = dmH2sq + dmH2sq3 * threelp;
 
 // add onto one-loop beta functions
 dmq3 *= threelp; dmq += dmq3;
 dml3 *= threelp; dml += dml3;
 dmu3 *= threelp; dmu += dmu3;
 dmd3 *= threelp; dmd += dmd3;
 dme3 *= threelp; dme += dme3;
 dhu3 *= threelp; dhu += dhu3;
 dhd3 *= threelp; dhd += dhd3;
 dhe3 *= threelp; dhe += dhe3;
      }
    } ///< end of 3-loop 
  }

  SoftPars<Susy, Brevity> dsoft(dsb, dmG, dhu, dhd, dhe, dmq, dmu,
		     dmd, dml, dme, dm3sq, dmH1sq, dmH2sq, 
		     displayGravitino(), displayMu(),
		     displayLoops(), displayThresholds());

  return dsoft;
}

// Outputs derivatives vector y[109] for SUSY parameters: interfaces to
// integration routines
template<class Susy, class Brevity>
DoubleVector SoftPars<Susy, Brevity>::beta() const {
  // calculate the derivatives
  static SoftPars<Susy, Brevity> dsoft; dsoft = beta2();

  return dsoft.display(); // convert to a long vector
}

// Outputs derivatives of anomalous dimensions, from which the running can be
// derived. 
template<class Susy, class Brevity>
void SoftPars<Susy, Brevity>::anomalousDeriv(DoubleMatrix & gEE, DoubleMatrix & gLL,
				  DoubleMatrix & gQQ, DoubleMatrix & gUU,
				  DoubleMatrix & gDD, 
				  double & gH1H1, double & gH2H2)  const {
  // Constants for gauge running
  static DoubleVector bBeta(3), cuBeta(3), cdBeta(3), ceBeta(3);
  static DoubleMatrix babBeta(3, 3);
  
  if (bBeta(1) < 1.0e-5) // Constants not set yet
    setBetas(babBeta, cuBeta, cdBeta, ceBeta, bBeta);
  
  // For calculational brevity: 
  static Brevity a;
  // convert to beta functions
  static Susy dsb;
  
  // calculate derivatives for full SUSY spectrum. Brevity calculations come
  // out encoded in a
  dsb = Susy::beta(a);
  
  static DoubleVector g1(3);
  g1 = displayGauge();
  
  // To keep this a const function: TIME SAVINGS
  DoubleMatrix &dt=a.dt, &ut=a.ut, &et=a.et;      
  DoubleVector &gsq=a.gsq; // &g3=a.g3, &g4=a.g4;
  
  static DoubleMatrix hu(3, 3), hd(3, 3), he(3, 3);
  static DoubleMatrix hut(3, 3), hdt(3, 3), het(3, 3);
  static DoubleVector gsqM(1, 3);
  
  hu = displayTrilinear(UA); hd = displayTrilinear(DA); 
  he = displayTrilinear(EA); 
  hut = hu.transpose(); hdt = hd.transpose(); het = he.transpose();
  
  double huuT = (hu * ut).trace(), hddT = (hd * dt).trace(), 
    heeT = (he * et).trace(); 
  gsqM = gsq * mGaugino; 
  
  if (displayLoops() > 0) { // CHECKED: agrees with our conventions
    const static double eightO3 = 8.0 / 3.0, oneO30 = 1.0 / 30.0;
    gLL = - (he * et + 0.3 * gsqM(1) + 1.5 * gsqM(2));
    gEE = - (2.0 * (et * he) + 1.2 * gsqM(1));
    gQQ = - (hd * dt + hu * ut + oneO30 * gsqM(1) + 1.5 * gsqM(2) +
	     eightO3 * gsqM(3));
    gDD = - (2.0 * dt * hd + 4.0 * oneO30 * gsqM(1) +
	     eightO3 * gsqM(3));
    gUU = - (2.0 * ut * hu + 16.0 * oneO30 * gsqM(1) +
	     eightO3 * gsqM(3));
    gH1H1 = - (3.0 * hddT + heeT + 0.3 * gsqM(1) + 1.5 * gsqM(2));
    gH2H2 = - (3.0 * huuT + 0.3 * gsqM(1) + 1.5 * gsqM(2));

    const static double oneO16Pisq = 1.0 / (16.0 * sqr(PI));
    gH1H1 = gH1H1 * oneO16Pisq;
    gH2H2 = gH2H2 * oneO16Pisq;
    gEE  *= oneO16Pisq;    
    gLL  *= oneO16Pisq;
    gQQ  *= oneO16Pisq;
    gUU  *= oneO16Pisq;
    gDD  *= oneO16Pisq;
  }
}

// Gives the ytilde terms relevant for the soft mass running: CHECKED 23/5/02
template<class Susy, class Brevity>
void SoftPars<Susy, Brevity>::yTildes(DoubleMatrix & yu, DoubleMatrix & yd, DoubleMatrix
			   &ye) const {
  ye = displaySoftMassSquared(mLl) * displayYukawaMatrix(YE) + 
    displayYukawaMatrix(YE) * displayMh1Squared() + 
    displayYukawaMatrix(YE) * displaySoftMassSquared(mEr);

  yd = displaySoftMassSquared(mQl) * displayYukawaMatrix(YD) + 
    displayYukawaMatrix(YD) * displayMh1Squared() + 
    displayYukawaMatrix(YD) * displaySoftMassSquared(mDr);

  yu = displaySoftMassSquared(mQl) * displayYukawaMatrix(YU) + 
    displayYukawaMatrix(YU) * displayMh2Squared() + 
    displayYukawaMatrix(YU) * displaySoftMassSquared(mUr);
}

/* 
   Give it a SUSY object and a value of M3/2, and it will return a soft
   object with AMSB soft breaking terms. Note that the sleptons will be
   tachyonic, ie nothing has been done to fix that problem.
   Note that in the following, we are neglecting all Yukawa couplings except
   that of the third family.
   
   THE CURRENT STATE OF PLAY:
   Two loop additions are possible, but a pain.
   */
template<class Susy, class Brevity>
void SoftPars<Susy, Brevity>::addAmsb(double maux) {
  Susy run(displaySusy());
  const double ONEO16pisq = 1.0 / (16. * sqr(PI));
  const double ONEO16pif = sqr(ONEO16pisq);
  double     g1   = run.displayGaugeCoupling(1), 
    g2   = run.displayGaugeCoupling(2), 
    g3   = run.displayGaugeCoupling(3); 

  // For calculational brevity: 
  static Brevity a;
  static Susy dsb;
  
  // calculate derivatives for full SUSY spectrum. Brevity calculations come
  // out encoded in a
  dsb = run.beta(a);

  /** three family version -- will be out by some two-loop terms. See
      hep-ph/9904378 */
  mQLsq = mQLsq + sqr(maux) * 
    (ONEO16pif * (-11. / 50. * sqr(sqr(g1)) - 1.5 * sqr(sqr(g2)) + 
		  8.0 * sqr(sqr(g3))) +
     ONEO16pisq * 0.5 *  
     (dsb.displayYukawaMatrix(YD) * displayYukawaMatrix(YD).transpose() +
      displayYukawaMatrix(YD) * dsb.displayYukawaMatrix(YD).transpose() +
      displayYukawaMatrix(YU) * dsb.displayYukawaMatrix(YU).transpose() +
      dsb.displayYukawaMatrix(YU) * displayYukawaMatrix(YU).transpose()));
  mLLsq = mLLsq + sqr(maux) * 
    (ONEO16pif * (-99. / 50. * sqr(sqr(g1)) - 1.5 * sqr(sqr(g2))) + 
     ONEO16pisq * 0.5 * 
     (dsb.displayYukawaMatrix(YE) * displayYukawaMatrix(YE).transpose() +
      displayYukawaMatrix(YE) * dsb.displayYukawaMatrix(YE).transpose())); 
  mURsq = mURsq + sqr(maux) * 
    (ONEO16pif * (-88. / 25. * sqr(sqr(g1)) + 8 * sqr(sqr(g3))) + 
     ONEO16pisq * 
     (dsb.displayYukawaMatrix(YU).transpose() * displayYukawaMatrix(YU) +
      displayYukawaMatrix(YU).transpose() * dsb.displayYukawaMatrix(YU)));
  mDRsq = mDRsq + sqr(maux) * 
    (ONEO16pif * (-22. / 25. * sqr(sqr(g1)) + 8 * sqr(sqr(g3))) + 
     ONEO16pisq * 
     (dsb.displayYukawaMatrix(YD).transpose() * displayYukawaMatrix(YD) +
      displayYukawaMatrix(YD).transpose() * dsb.displayYukawaMatrix(YD)));
  mSEsq = mSEsq + sqr(maux) * 
    (ONEO16pif * (-198. / 25. * sqr(sqr(g1))) + 
     ONEO16pisq * 
     (dsb.displayYukawaMatrix(YE).transpose() * displayYukawaMatrix(YE) +
      displayYukawaMatrix(YE).transpose() * dsb.displayYukawaMatrix(YE)));
  mH1sq = mH1sq + sqr(maux) * 
    (ONEO16pif * (-99. / 50. * sqr(sqr(g1)) - 1.5 * sqr(sqr(g2))) + 
     ONEO16pisq * (3.0 * (dsb.displayYukawaMatrix(YD) *
			  displayYukawaMatrix(YD).transpose()).trace() 
		   + (dsb.displayYukawaMatrix(YE) * 
		      displayYukawaMatrix(YE).transpose()).trace()));
  mH2sq = mH2sq + sqr(maux) * 
    (ONEO16pif * (-99. / 50. * sqr(sqr(g1)) - 1.5 * sqr(sqr(g2))) + 
     ONEO16pisq * (3.0 * (dsb.displayYukawaMatrix(YU) * 
			  displayYukawaMatrix(YU).transpose()).trace()));

  // AMSB spectrum
  DoubleVector amsbGaugino(3);
  DoubleMatrix temp(3, 3);

  int i; for (i=1; i<=3; i++)
    amsbGaugino(i) = dsb.displayGaugeCoupling(i) * maux / 
      run.displayGaugeCoupling(i);
  mGaugino = mGaugino + amsbGaugino;
  
  ua = ua - dsb.displayYukawaMatrix(YU) * maux;
  da = da - dsb.displayYukawaMatrix(YD) * maux;
  ea = ea - dsb.displayYukawaMatrix(YE) * maux;
  
  m32 = maux;
  
  //  u1R_PQflip();
  return;
}


template<class Susy, class Brevity>
void SoftPars<Susy, Brevity>::u1R_PQflip() {
  setSusyMu(-displaySusyMu());
  mGaugino = -1. * mGaugino;
  ua = -1. * ua;
  da = -1. * da;
  ea = -1. * ea;
}

// Reads in universal boundary conditions at the current scale:
// m0, M1/2, A0, B and sign of mu
template<class Susy, class Brevity>
void SoftPars<Susy, Brevity>::universal(double m0,  double m12,  double a0,  double mu,
			      double m3sq) {
  standardSugra(m0, m12, a0);  
  setSusyMu(mu);
  setM3Squared(m3sq);
}

template<class Susy, class Brevity>
void SoftPars<Susy, Brevity>::universalScalars(double m0) {
  // scalar masses
  DoubleMatrix ID(3, 3), mm0(3, 3);
  int i; for (i=1; i<=3; i++) ID(i, i) = 1.0;
  mm0 = ID * sqr(m0);
  setSoftMassMatrix(mQl, mm0); setSoftMassMatrix(mUr, mm0);
  setSoftMassMatrix(mDr, mm0); setSoftMassMatrix(mLl, mm0);
  setSoftMassMatrix(mEr, mm0);
  setMh1Squared(sqr(m0)); setMh2Squared(sqr(m0));
}

template<class Susy, class Brevity>
void SoftPars<Susy, Brevity>::universalGauginos(double m12) {  
  // gaugino masses
  int i; for (i=1; i<=3; i++) setGauginoMass(i, m12);
}

template<class Susy, class Brevity>
void SoftPars<Susy, Brevity>::universalTrilinears(double a0)  {  
  // trilinears
  setTrilinearMatrix(UA, a0 * displayYukawaMatrix(YU)); 
  setTrilinearMatrix(DA, a0 * displayYukawaMatrix(YD));
  setTrilinearMatrix(EA, a0 * displayYukawaMatrix(YE));
}

// Input m0, NOT m0 squared.
template<class Susy, class Brevity>
void SoftPars<Susy, Brevity>::standardSugra(double m0,  double m12, double a0) {
  /*  if (m0 < 0.0) {
    std::ostringstream ii;
    ii << "m0=" << m0 << " passed to universal boundary" <<
      "conditions illegally negative.";
    throw ii.str();
    } Deleted on request from A Pukhov */
  universalScalars(m0);
  universalGauginos(m12);
  universalTrilinears(a0);
}

template<class Susy, class Brevity>
std::ostream & operator <<(std::ostream &left, const SoftPars<Susy, Brevity> &s) {
  left << "SUSY breaking MSSM parameters at Q: " << s.displayMu() << endl;
  left << " UA" << s.displayTrilinear(UA) 
       << " UD" << s.displayTrilinear(DA) 
       << " UE" << s.displayTrilinear(EA); 
  left << " mQLsq" << s.displaySoftMassSquared(mQl) 
       << " mURsq" << s.displaySoftMassSquared(mUr) 
       << " mDRsq" << s.displaySoftMassSquared(mDr) 
       << " mLLsq" << s.displaySoftMassSquared(mLl) 
       << " mSEsq" << s.displaySoftMassSquared(mEr);
  left << "m3sq: " << s.displayM3Squared() << " mH1sq: " <<
    s.displayMh1Squared() << " mH2sq: " << s.displayMh2Squared() << '\n';
  left << "Gaugino masses" << s.displayGaugino();
  left << s.displaySusy();
  return left;
}

// Boundary conditions to be applied at messenger scale for Gauge mediated
// SUSY breaking (see hep-ph/9703211 for example)
template<class Susy, class Brevity>
void SoftPars<Susy, Brevity>::minimalGmsb(int n5, double lambda, double mMess, 
			       double cgrav) {

// Modified thresholds by JEL 1-26-04 to accomodate numerical infinities

  const double epstol = 1.0e-4;
  double x = lambda / mMess;

  double f, g;

  if(fabs(x) < epstol) {
    g = 1.0 + x*x/6.0 + sqr(x*x)/15.0;
    f = 1.0 + x*x/36.0 + 11.0*sqr(x*x)/450.0;
  }
  else if(fabs(x-1.0) < 0.0001) {
    g  =  log(4.0);
    f  = -sqr(PI)/6.0 + log(4.0) + 0.5*sqr(log(4.0));
    g -=  0.0008132638905771205626;
    f -= -0.0049563838821509165200;
  }
  else {
    g = 1.0 / sqr(x) * 
      ((1.0 + x) * log(1.0 + x) + (1.0 - x) * log(1.0 - x));
    f = (1.0 + x) / sqr(x) * 
    (log(1.0 + x) - 2.0 * dilog(x / (1.0 + x)) + 0.5 * 
     dilog(2.0 * x / (1.0 + x))) + 
     (1.0 - x) / sqr(x) * (log(1.0 - x) - 2.0 * dilog(-x / (1.0 - x)) +
			 0.5 * dilog(-2.0 * x / (1.0 - x)));
  }

  double n5d = double(n5);

  /// There is a relative minus in the mGMSB conditions for gaugino masses,
  /// since these equations are for L=-M/2 lambda lambda. See hep-ph/9801271:
  /// BCA 27/7/12
  double m1, m2, m3;
  m1 = n5d * sqr(displayGaugeCoupling(1)) / (16.0 * sqr(PI)) * lambda * g; 
  m2 = n5d * sqr(displayGaugeCoupling(2)) / (16.0 * sqr(PI)) * lambda * g; 
  m3 = n5d * sqr(displayGaugeCoupling(3)) / (16.0 * sqr(PI)) * lambda * g; 
  setGauginoMass(1, m1);   setGauginoMass(2, m2);   setGauginoMass(3, m3);

  setM32(2.37e-19 * lambda * mMess * cgrav);

  double g1f = sqr(sqr(displayGaugeCoupling(1)));
  double g2f = sqr(sqr(displayGaugeCoupling(2)));
  double g3f = sqr(sqr(displayGaugeCoupling(3)));

  double mursq, mdrsq, mersq, mqlsq, mllsq;
  mursq = 2.0 * f * sqr(lambda) * n5d * 
    (4.0 / 3.0 * g3f + 0.6 * 4.0 / 9.0 * g1f) 
    / sqr(16.0 * sqr(PI));
  mdrsq = 2.0 * f * sqr(lambda) * n5d * 
    (4.0 / 3.0 * g3f + 0.6 * 1.0 / 9.0 * g1f) 
    / sqr(16.0 * sqr(PI));
  mersq = 2.0 * f * sqr(lambda) * n5d * 
    (0.6 * g1f) 
    / sqr(16.0 * sqr(PI));
  mqlsq = 2.0 * f * sqr(lambda) * n5d * 
    (4.0 / 3.0 * g3f + 0.75 * g2f + 0.6 * g1f / 36.0) 
    / sqr(16.0 * sqr(PI));
  mllsq = 2.0 * f * sqr(lambda) * n5d * 
    (                  0.75 * g2f + 0.6 * 0.25 * g1f) 
    / sqr(16.0 * sqr(PI));

  // You need Higgs masses too!

  DoubleMatrix id(3, 3);
  id(1, 1) = 1.0; id(2, 2) = 1.0; id(3, 3) = 1.0;

  setSoftMassMatrix(mQl, mqlsq * id);
  setSoftMassMatrix(mUr, mursq * id);
  setSoftMassMatrix(mDr, mdrsq * id);
  setSoftMassMatrix(mLl, mllsq * id);  
  setMh1Squared(mllsq);
  setMh2Squared(mllsq);
  setSoftMassMatrix(mEr, mersq * id);

  universalTrilinears(0.0);
}

#endif
