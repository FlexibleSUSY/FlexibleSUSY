/** \file softsusy.cpp 
    Project: SOFTSUSY 
    Author: Ben Allanach 
    Manual: hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
    Webpage: http://hepforge.cedar.ac.uk/softsusy/
*/

#ifdef SOFTSUSY_H

using std::fabs;
using std::sqrt;
using std::exp;
using std::log;
using std::cout;
using std::endl;

extern double sw2, gnuL, guL, gdL, geL, guR, gdR, geR, yuL, yuR, ydL,
  ydR, yeL, yeR, ynuL;


template<class SoftPars>
const Softsusy<SoftPars>& Softsusy<SoftPars>::operator=(const Softsusy<SoftPars>& s) {
  if (this == &s) return *this;
  setSoftPars(s.displaySoftPars());
  setAltEwsbMssm(s.displayAltEwsbMssm());
  physpars = s.physpars;
  problem = s.problem;
  mw = s.mw;
  minV = s.minV;
  forLoops = s.forLoops;
  msusy = s.msusy;
  dataSet = s.dataSet;
  fracDiff = s.fracDiff;
  setMu(s.displayMu());
  setLoops(s.displayLoops());
  setThresholds(s.displayThresholds());
  setSetTbAtMX(s.displaySetTbAtMX());
  altEwsb = s.altEwsb;
  predMzSq = s.displayPredMzSq();
  t1OV1Ms = s.displayTadpole1Ms(); 
  t2OV2Ms = s.displayTadpole2Ms(); 
  t1OV1Ms1loop = s.displayTadpole1Ms1loop(); 
  t2OV2Ms1loop = s.displayTadpole2Ms1loop(); 
  alternativeMs = s.alternativeMs;
  mxBC = s.displayMxBC();

  return *this;
}

/// Returns mu from rewsb requirement. 
/// returns 1 if there's a problem. Call at MSusy
template<class SoftPars>
int Softsusy<SoftPars>::rewsbMu(int sgnMu, double & mu) const {
  int flag = 0;
   if (abs(sgnMu) != 1) {
    std::ostringstream ii;     
    ii << "Error: sign mu = " << sgnMu << "\n";
    throw ii.str();
  }
  double mH1sq = displayMh1Squared(), mH2sq = displayMh2Squared(), tanb =
    displayTanb(), tanb2 =  sqr(tanb);
  
  /// Tree-level relation
  double musq = (mH1sq - mH2sq * tanb2) 
    / (tanb2 - 1.0) - 0.5 * sqr(displayMz());
  
  if (musq < 0.0) {
    mu = sgnMu * sqrt(fabs(musq)); flag = 1; /// mu has incorrect sign
  }
  else
    mu = sgnMu * sqrt(musq); 

  return flag;
}

/// returns 1 if mu < 1.0e-9
template<class SoftPars>
int Softsusy<SoftPars>::rewsbM3sq(double mu, double & m3sq) const {
  int flag = 0;

  if (fabs(mu) < 1.0e-9) 
    { flag = 1; m3sq = 0.0; }
  else
    m3sq =  0.5 * 
	 (displayMh1Squared() + displayMh2Squared() - displayTadpole2Ms() -
	  displayTadpole1Ms() + 2.0 * sqr(mu)) * 
      sin(2 * atan(displayTanb()));
  
  /// Following means no good rewsb
  if (m3sq < 0.0) flag = 1;
  else if (testNan(m3sq)) {
    flag = 1; m3sq = EPSTOL;
  }

  return flag;
}


/// Predicts tan beta once mu and soft terms are predicted at low energy
/// Useful for fine-tuning calculation. Call at MSusy only.
/*double MssmSoftsusy::predTanb() const  {
  double sin2t = 2.0 * displayM3Squared() / 
    (displayMh1Squared() - displayTadpole1Ms() + 
     displayMh2Squared() - displayTadpole2Ms() + 2.0 *
     sqr(susyMu)); 
  
  /// Note: we want to take inverse sine so that fundamental domain is greater
  /// than pi/4. sin(pi - 2 beta)=sin 2 beta should achieve this.
  /// we also use tan (pi/2 - theta) = 1/tan(theta)
  double theta;
  if (fabs(sin2t) < 1.0) theta = asin(sin2t) * 0.5;
  else return 0.0;
  
  return 1.0 / tan(theta);
  }*/
/// Predicts tan beta once mu and soft terms are predicted at low energy
/// Useful for fine-tuning calculation. Call at MSusy only.
template<class SoftPars>
double Softsusy<SoftPars>::predTanb(double susyMu) const  {
  if (susyMu < -6.e66) susyMu = displaySusyMu();

  double sin2t = 2.0 * displayM3Squared() / 
    (displayMh1Squared() - displayTadpole1Ms() + 
     displayMh2Squared() - displayTadpole2Ms() + 2.0 *
     sqr(susyMu)); 
  
  /// Note: we want to take inverse sine so that fundamental domain is greater
  /// than pi/4. sin(pi - 2 beta)=sin 2 beta should achieve this.
  /// we also use tan (pi/2 - theta) = 1/tan(theta)
  double theta;
  if (fabs(sin2t) < 1.0) theta = asin(sin2t) * 0.5;
  else return 0.0;
  
  return 1.0 / tan(theta);
}


template<class SoftPars>
void Softsusy<SoftPars>::doTadpoles(double mt, double sinthDRbar) {

    calcTadpole1Ms1loop(mt, sinthDRbar);
    calcTadpole2Ms1loop(mt, sinthDRbar);
    
    t1OV1Ms = t1OV1Ms1loop;
    t2OV2Ms = t2OV2Ms1loop;

    /// tachyons tend to screw up this, so only calculate if we don't have them
    if (numRewsbLoops > 1) {
      /// add the two-loop terms, prepare inputs
      double s1s = 0., s2s = 0., s1t = 0., s2t = 0.,
	gs = displayGaugeCoupling(3), 
	rmtsq = sqr(forLoops.mt), scalesq = sqr(displayMu()), 
	vev2 = sqr(displayHvev()), tanb = displayTanb(), 
	amu = -displaySusyMu(), mg = displayGaugino(3), 
	mAsq = sqr(forLoops.mA0(1)); 
      
      double sxt = sin(forLoops.thetat), cxt = cos(forLoops.thetat);
      double mst1sq = sqr(forLoops.mu(1, 3)), 
	mst2sq = sqr(forLoops.mu(2, 3));
      /// two-loop Higgs corrections: alpha_s alpha_b
      double sxb = sin(forLoops.thetab), 
	cxb = cos(forLoops.thetab);
      double sintau = sin(forLoops.thetatau), 
	costau = cos(forLoops.thetatau);
      double msnusq = sqr(forLoops.msnu(3));
      double msb1sq = sqr(forLoops.md(1, 3)), 
	msb2sq = sqr(forLoops.md(2, 3));
      double mstau1sq = sqr(forLoops.me(1, 3)), 
	mstau2sq = sqr(forLoops.me(2, 3));
      double cotbeta = 1.0 / tanb;
      double rmbsq = sqr(forLoops.mb);
      double rmtausq = sqr(forLoops.mtau);
      double s1b = 0.0, s2b = 0.0, s1tau = 0.0, s2tau = 0.0;

      ewsb2loop_(&rmtsq, &mg, &mst1sq, &mst2sq, &sxt, &cxt, &scalesq, 
        	 &amu, &tanb, &vev2, &gs, &s1s, &s2s);
      ddstad_(&rmtsq, &rmbsq, &mAsq, &mst1sq, &mst2sq, &msb1sq, &msb2sq, 
              &sxt, &cxt, &sxb, &cxb, &scalesq, &amu, &tanb, &vev2, &s1t, 
      	      &s2t);
      ewsb2loop_(&rmbsq, &mg, &msb1sq, &msb2sq, &sxb, &cxb, &scalesq,
        	 &amu, &cotbeta, &vev2, &gs, &s2b, &s1b);
      tausqtad_(&rmtausq, &mAsq, &msnusq, &mstau1sq, &mstau2sq, &sintau, 
        	&costau, &scalesq, &amu, &tanb, &vev2, &s1tau, &s2tau);

      if (!testNan(s1s * s1t * s1b * s1tau * s2s * s2t * s2b * s2tau)) {
	t1OV1Ms += - s1s - s1t - s1b - s1tau;
	t2OV2Ms += - s2s - s2t - s2b - s2tau; 
	/// end of 2-loop bit
      }
      else  {
	flagNoMuConvergence(true);
	if (PRINTOUT > 1) cout << "2-loop tadpoles are nans\n";
      }
    }
}

//PA: fixes trilnear H1-sfermion-sfermion couplings 
//for use in doCalcTadpole1oneLoop
template<class SoftPars>
void Softsusy<SoftPars>::H1SfSfCouplings(DoubleMatrix & lTS1Lr, DoubleMatrix & lBS1Lr, DoubleMatrix & lTauS1Lr, double gmzOcthW, double mu, double cosb, double v1) const{
  double ht = displayDrBarPars().ht;  
  double hbsq = sqr(displayDrBarPars().hb); 
  double htausq = sqr(displayDrBarPars().htau);
    
    lTS1Lr(1, 1) = gmzOcthW * guL * cosb;
    lTS1Lr(1, 2) = ht * mu / root2;
    lTS1Lr(2, 1) = lTS1Lr(1, 2);
    lTS1Lr(2, 2) = gmzOcthW * guR * cosb;
    
    lBS1Lr(1, 1) = gmzOcthW * gdL * cosb + hbsq * v1;
    lBS1Lr(1, 2) = forLoops.ub / root2;
    lBS1Lr(2, 1) = lBS1Lr(1, 2);
    lBS1Lr(2, 2) = gmzOcthW * gdR * cosb + hbsq * v1;
  
    lTauS1Lr(1, 1) = gmzOcthW * geL * cosb + htausq * v1;
    lTauS1Lr(1, 2) = forLoops.utau / root2;
    lTauS1Lr(2, 1) = lTauS1Lr(1, 2);
    lTauS1Lr(2, 2) = gmzOcthW * geR * cosb + htausq * v1;

 }


//PA: fixes trilnear H1-sfermion-sfermion couplings 
//for use in doCalcTadpole1oneLoop
template<class SoftPars>
void Softsusy<SoftPars>::H2SfSfCouplings(DoubleMatrix & lTS2Lr, DoubleMatrix & lBS2Lr, DoubleMatrix & lTauS2Lr, double gmzOcthW, double mu, double sinb) const{
  double hb = displayDrBarPars().hb;  
  double htsq = sqr(displayDrBarPars().ht); 
  double htau = displayDrBarPars().htau;
  double v2 = displayHvev() * sinb;
  lTS2Lr(1, 1) = - gmzOcthW * guL * sinb + htsq * v2;
  lTS2Lr(1, 2) = forLoops.ut / root2;
  lTS2Lr(2, 1) = lTS2Lr(1, 2);
  lTS2Lr(2, 2) = - gmzOcthW * guR * sinb + htsq * v2;
  
  lBS2Lr(1, 1) = - gmzOcthW * gdL * sinb;
  lBS2Lr(1, 2) = hb * mu / root2;
  lBS2Lr(2, 1) = lBS2Lr(1, 2);
  lBS2Lr(2, 2) = - gmzOcthW * gdR * sinb;
  
  lTauS2Lr(1, 1) = - gmzOcthW * geL * sinb;
  lTauS2Lr(1, 2) = htau * mu / root2;
  lTauS2Lr(2, 1) = lTauS2Lr(1, 2);
  lTauS2Lr(2, 2) = - gmzOcthW * geR * sinb;  
    
 }


//PA: routine to calculate sfermiom contributions to (16 \pi^2) t1 / v1
template<class SoftPars>
double Softsusy<SoftPars>::doCalcTad1Sfermions(DoubleMatrix lTS1Lr, DoubleMatrix lBS1Lr, DoubleMatrix lTauS1Lr, double costhDRbar) const {
  double g = displayGaugeCoupling(2), mz = displayMzRun();
  double tanb = displayTanb(), cosb = cos(atan(tanb));
  double q = displayMu(); 
  double gO2mwcosb = g / (2.0 * displayMwRun() * cosb);
  /// sneutrino coupling
  double lSnu = g * mz / costhDRbar * 0.5 * cosb;

   /// stop, sbottom, stau couplings
  DoubleMatrix lTS112(2, 2),lBS112(2, 2),lTauS112(2, 2), rotate(2, 2);

  rotate = rot2d(forLoops.thetat);
  lTS112 = rotate * lTS1Lr * rotate.transpose();
  
  rotate = rot2d(forLoops.thetab);
  lBS112 = rotate * lBS1Lr * rotate.transpose();
  
  rotate = rot2d(forLoops.thetatau);
  lTauS112 = rotate * lTauS1Lr * rotate.transpose();
   
  double stops = 0., sbots = 0.;
  /// third generation squarks
  stops = stops + 3.0 * gO2mwcosb * lTS112(1, 1) *
    a0(forLoops.mu(1, 3), q);
  stops = stops + 3.0 * gO2mwcosb * lTS112(2, 2) *
    a0(forLoops.mu(2, 3), q);
  sbots = sbots + 3.0 * gO2mwcosb * lBS112(1, 1) *
    a0(forLoops.md(1, 3), q);
  sbots = sbots + 3.0 * gO2mwcosb * lBS112(2, 2) *
    a0(forLoops.md(2, 3), q);

  /// third generation sleptons
  double staus = gO2mwcosb * lTauS112(1, 1) *
    a0(forLoops.me(1, 3), q);
  staus = staus + gO2mwcosb * lTauS112(2, 2) *
    a0(forLoops.me(2, 3), q);
  
  double sneuts = 0.0;
  sneuts = sneuts + gO2mwcosb * lSnu * a0(forLoops.msnu(3), q);
  
  int family; 
  double gsqMZO2mwcw = sqr(g) * mz * 0.5 / (displayMwRun() * costhDRbar);
  double sups = 0.0, sdowns = 0.0;
  /// first two families of squarks
  for (family = 1; family <=2; family++) {
    sups = sups + 3.0 * gsqMZO2mwcw * 
      ( guL * a0(forLoops.mu(1, family), q) + 
	guR * a0(forLoops.mu(2, family), q)); 
    sdowns = sdowns + 3.0 * gsqMZO2mwcw * 
      ( gdL * a0(forLoops.md(1, family), q) +
	gdR * a0(forLoops.md(2, family), q));
  }
  
  double sleps = 0.;
  /// sleptons
  for (family = 1; family <=2; family++) {
    sleps = sleps + gsqMZO2mwcw * 
      (geL * a0(forLoops.me(1, family), q) +
       geR * a0(forLoops.me(2, family), q));
    sneuts = sneuts + gO2mwcosb * lSnu * a0(forLoops.msnu(family), q);
  }
  
   double sfermions = stops + sbots + staus + sneuts + sleps + sups + sdowns;
   return sfermions;

 }

//PA: routine to calculate sfermiom contributions to (16 \pi^2) t1 / v1
template<class SoftPars>
double Softsusy<SoftPars>::doCalcTad2Sfermions(DoubleMatrix lTS2Lr, DoubleMatrix lBS2Lr, DoubleMatrix lTauS2Lr, double costhDRbar) const {
  double g = displayGaugeCoupling(2), mz = displayMzRun();
  double tanb = displayTanb(), sinb = sin(atan(tanb));
  double q = displayMu(); 
  DoubleMatrix lTS212(2, 2),  lBS212(2, 2), lTauS212(2, 2), rotate(2, 2);
  rotate = rot2d(forLoops.thetat);
  lTS212 = rotate * lTS2Lr * rotate.transpose();
  rotate = rot2d(forLoops.thetab);
  lBS212 = rotate * lBS2Lr * rotate.transpose(); 
  rotate = rot2d(forLoops.thetatau);
  lTauS212 = rotate * lTauS2Lr * rotate.transpose();
  /// sneutrino coupling
  double lSnu = - g * mz / costhDRbar * 0.5 * sinb;
  /// third generation squarks
  double sfermions = 3.0 * g * lTS212(1, 1) / (2.0 * displayMwRun() * sinb) *
    a0(forLoops.mu(1, 3), q);
  sfermions = sfermions + 3.0 * g * lTS212(2, 2) / 
    (2.0 * displayMwRun() * sinb) *
    a0(forLoops.mu(2, 3), q);
  sfermions = sfermions + 3.0 * g * lBS212(1, 1) / 
    (2.0 * displayMwRun() * sinb) *
    a0(forLoops.md(1, 3), q);
  sfermions = sfermions + 3.0 * g * lBS212(2, 2) / (2.0 * displayMwRun() * sinb) *
    a0(forLoops.md(2, 3), q);
  /// third generation sleptons
  sfermions = sfermions + g * lTauS212(1, 1) / (2.0 * displayMwRun() * sinb) *
    a0(forLoops.me(1, 3), q);
  sfermions = sfermions + g * lTauS212(2, 2) / (2.0 * displayMwRun() * sinb) *
    a0(forLoops.me(2, 3), q);
  sfermions = sfermions + 
    g * lSnu / (2.0 * displayMwRun() * sinb) * a0(forLoops.msnu(3), q);
  
  for (int family = 1; family <=2; family++)
    {
      sfermions = sfermions + 3.0 * g * 
	(- g * mz / costhDRbar * guL * sinb * a0(forLoops.mu(1, family), q) 
	 - g * mz / costhDRbar * guR * sinb * a0(forLoops.mu(2, family), q)) 
	/ (2.0 * displayMwRun() * sinb);
      sfermions = sfermions + 3.0 * g * 
	(- g * mz / costhDRbar * gdL * sinb * a0(forLoops.md(1, family), q) 
	 - g * mz / costhDRbar * gdR * sinb * a0(forLoops.md(2, family), q)) 
	/ (2.0 * displayMwRun() * sinb);	
    }
  
  /// sleptons
  for (int family = 1; family <=2; family++)
    {
      sfermions = sfermions + g * 
	(- g * mz / costhDRbar * geL * sinb * a0(forLoops.me(1, family), q) 
	 - g * mz / costhDRbar * geR * sinb * a0(forLoops.me(2, family), q)) 
	/ (2.0 * displayMwRun() * sinb);
      sfermions = sfermions + g * lSnu * a0(forLoops.msnu(family), q)
	/ (2.0 * displayMwRun() * sinb);
    }
  
  return sfermions;

}
//PA: fixes trilnear H1-fermion-fermion couplings 
//for use in doCalcTadpole1oneLoop  
template<class SoftPars>
double Softsusy<SoftPars>::doCalcTad1fermions(double q, double v1) const {
  double mb = forLoops.mb;
  double mtau = forLoops.mtau;
  double hb = displayDrBarPars().hb, htau = forLoops.htau;
  /// bottom quark and tau, ignore others - factor (10^-2)^3 down
  /// I have included the bottom pole mass in the propagators and the Yukawa
  /// for the coupling, hence BPMZ's hb is written mb * root2 / v1
  double  fermions = - 6.0 * hb * mb * root2 / v1 * a0(mb, q) 
      - 2.0 * htau * mtau * root2 / v1 * a0(mtau, q);
  return fermions;
}

//PA: fixes trilnear H2-fermion-fermion couplings 
//for use in doCalcTadpole1oneLoop  
template<class SoftPars> 
double Softsusy<SoftPars>::doCalcTad2fermions(double q) const {
  double mtop = forLoops.mt, ht = displayDrBarPars().ht;
  /// top quark, ignore others - factor (10^-2)^3 down
  double fermions = - 6.0 * sqr(ht) * a0(mtop, q);
  return fermions;
}

//one loop H1 tadpole contributions from Higgs bosons in the loops
// Follwing BPMZ Goldstone bosons are not included in this.
template<class SoftPars>
double Softsusy<SoftPars>::doCalcTad1Higgs(double q, double costhDRbar2, 
                                           double g, double tanb) const {
  double mA = forLoops.mA0(1),  mh = forLoops.mh0(1), 
     mH0 = forLoops.mh0(2), mHp = forLoops.mHpm;
  double alpha = forLoops.thetaH, sina2 = sqr(sin(alpha)), 
    cosa2 = 1.0 - sina2, cos2b = cos(2.0 * atan(tanb)), 
    sin2a = sin(2.0 * alpha);
  double higgs = 0.0;
  higgs = higgs - sqr(g) * cos2b / (8.0 * costhDRbar2) *
    (a0(mA, q) + 2.0 * a0(mHp, q)) +
    sqr(g) * a0(mHp, q) * 0.5 +
    sqr(g) / (8.0 * costhDRbar2) * a0(mh, q) *
    (3.0 * sina2 - cosa2 + sin2a * tanb) +
    sqr(g) / (8.0 * costhDRbar2) * a0(mH0, q) *
    (3.0 * cosa2 - sina2 - sin2a * tanb); 
  return higgs;
}

//one loop H2 tadpole contributions from Higgs bosons in the loops
// Follwing BPMZ Goldstone bosons are not included in this.
template<class SoftPars>
double Softsusy<SoftPars>::doCalcTad2Higgs(double q, double costhDRbar2, 
                                           double g, double tanb) const {
    /// Higgs
  double alpha = forLoops.thetaH, sina2 = sqr(sin(alpha)), cosa2 = 1.0 -
    sina2, cos2b = cos(2.0 * atan(tanb)), 
    mh = forLoops.mh0(1), sin2a = sin(2.0 * alpha);
  double mA = forLoops.mA0(1), mH0 = forLoops.mh0(2), mHp =
    forLoops.mHpm;
    
  double higgs = 0.0;
  higgs = sqr(g) * cos2b / (8.0 * costhDRbar2) *
    (a0(mA, q) + 2.0 * a0(mHp, q)) +
    sqr(g) * a0(mHp, q) * 0.5 +
    sqr(g) / (8.0 * costhDRbar2) * a0(mh, q) *
    (3.0 * cosa2 - sina2 + sin2a / tanb) +
    sqr(g) / (8.0 * costhDRbar2) * a0(mH0, q) *
    (3.0 * sina2 - cosa2 - sin2a / tanb);
  
  return higgs;
}
template<class SoftPars>
double Softsusy<SoftPars>::doCalcTad1Neutralinos(double q, double costhDRbar, 
                                                 double g, double cosb) const {
  ComplexMatrix n(forLoops.nBpmz);
  DoubleVector mneut(forLoops.mnBpmz);
  
  double neutralinos = 0.0;
  double tanthDRbar = tan(acos(costhDRbar));
  for (int family = 1; family <= 4; family++)
    neutralinos = neutralinos - 
      sqr(g) * mneut(family) / (displayMwRun() * cosb) *
      (n(family, 3) * (n(family, 2) - n(family, 1) * tanthDRbar)).real() * 
      a0(mneut(family), q);

  return neutralinos;
}

template<class SoftPars>
double Softsusy<SoftPars>::doCalcTad2Neutralinos(double q, double costhDRbar, 
                                                 double g, double sinb) const {
  ComplexMatrix n(forLoops.nBpmz);
  DoubleVector mneut(forLoops.mnBpmz);
  double tanthDRbar = tan(acos(costhDRbar));
  double neutralinos = 0.0;
  for (int family = 1; family <= 4; family++)
    neutralinos = neutralinos + sqr(g) * mneut(family) / 
      (displayMwRun() * sinb) *
      (n(family, 4) * (n(family, 2) - n(family, 1) * tanthDRbar)).real() * 
      a0(mneut(family), q); 

  return neutralinos;
}

template<class SoftPars>
double Softsusy<SoftPars>::doCalcTad1Charginos(double q, double g, 
                                                double cosb) const {
  ComplexMatrix u(forLoops.uBpmz), v(forLoops.vBpmz); 
  DoubleVector mch(forLoops.mchBpmz); 
  double charginos = 0.0;
  for (int family=1; family<=2; family++)
    charginos = charginos - root2 * sqr(g) / (displayMwRun() * cosb)
      * mch(family) * (v(family, 1) * u(family, 2)).real()
	 * a0(mch(family), q);
  return charginos;
}

template<class SoftPars>
double Softsusy<SoftPars>::doCalcTad2Charginos(double q, double g, 
                                                double sinb) const {
ComplexMatrix u(forLoops.uBpmz), v(forLoops.vBpmz); 
  DoubleVector mch(forLoops.mchBpmz); 
  double charginos = 0.0;
  for (int family = 1; family <= 2; family++)
    charginos = charginos - root2 * sqr(g) * mch(family) /
      (displayMwRun() * sinb)
      * (v(family, 2) * u(family, 1)).real() * a0(mch(family), q);

  return charginos;
 }

template<class SoftPars>
double Softsusy<SoftPars>::doCalcTad1GaugeBosons(double q, double costhDRbar2, 
                                                 double g, double tanb) const {
  double cos2b = cos(2.0 * atan(tanb));
  double mz = displayMzRun();
  double mw = displayMwRun();

  double gaugeBosons = 0.0;
  gaugeBosons = gaugeBosons + 3.0 * sqr(g) / 4.0 * 
    (2.0 * a0(mw, q) + a0(mz, q) / costhDRbar2)
    + sqr(g) * cos2b / (8.0 * costhDRbar2) * (2.0 * a0(mw, q) + 
					      a0(mz, q));
  return gaugeBosons;
}

  template<class SoftPars>
double Softsusy<SoftPars>::doCalcTad2GaugeBosons(double q, double costhDRbar2, 
                                                 double g, double tanb) const {
  /// Weak bosons
  double cos2b = cos(2.0 * atan(tanb));
  double mz = displayMzRun();
  double mw = displayMwRun();
  double gaugeBosons = 3.0 * sqr(g) / 4.0 * 
    (2.0 * a0(mw, q) + a0(mz, q) / costhDRbar2)
    - sqr(g) * cos2b / (8.0 * costhDRbar2) * 
    (2.0 * a0(mw, q) + a0(mz, q));
  
  return gaugeBosons;
}

/// From hep-ph/9606211's appendix. It should be done at MSusy to minimize the
/// 1-loop contributions. Only call if you've calculated drbarpars.
/// inputs are running top/bottom masses: call at MSusy only
template<class SoftPars>
double Softsusy<SoftPars>::doCalcTadpole1oneLoop(double /* mt */, double sinthDRbar) {

 if (forLoops.mu(1, 3) == 0.0 || forLoops.mu(2, 3) == 0.0) {
   if (PRINTOUT > 1)
    cout << "Trying to calculate tadpole without having first calculated"
	 << " the DRbar masses.\n";
   return 0.0; 
  }

  double g = displayGaugeCoupling(2), mz = displayMzRun(),  
     costhDRbar = cos(asin(sinthDRbar)), costhDRbar2 = sqr(costhDRbar), 
     tanb = displayTanb(), cosb = cos(atan(tanb)), 
    mu = -displaySusyMu(), q = displayMu();
  double beta = atan(displayTanb());
  double v1 = displayHvev() * cos(beta);
  double gmzOcthW =  g * mz / costhDRbar;

  double fermions =  doCalcTad1fermions(q, v1);
  /// PA: stop, sbottom, stau, couplings in the left right basis 
  // will be stored in these matrices 
  DoubleMatrix lTS1Lr(2, 2), lBS1Lr(2, 2), lTauS1Lr(2, 2);
  H1SfSfCouplings(lTS1Lr, lBS1Lr, lTauS1Lr, gmzOcthW, mu, cosb, v1);
  //PA: Now we take these couplings and obtain sfermion contributions
  double sfermions =  doCalcTad1Sfermions(lTS1Lr, lBS1Lr, lTauS1Lr, 
					  costhDRbar);
 /// Higgs
  double higgs = doCalcTad1Higgs(q, costhDRbar2, g, tanb);
  /// Neutralinos
  double neutralinos = doCalcTad1Neutralinos(q, costhDRbar, g, cosb);
  /// Charginos
  double charginos = doCalcTad1Charginos(q, g, cosb);
  /// Weak bosons
  double gaugeBosons = doCalcTad1GaugeBosons(q, costhDRbar2, g, tanb);
  double delta = fermions + sfermions + higgs + charginos + neutralinos + 
    gaugeBosons;

  return delta / (16.0 * sqr(PI));
}

/// checked
template<class SoftPars>
void Softsusy<SoftPars>::calcTadpole1Ms1loop(double mt, double sinthDRbar) { 
  t1OV1Ms1loop = doCalcTadpole1oneLoop(mt, sinthDRbar);
  if (testNan(t1OV1Ms1loop)) {
    t1OV1Ms1loop = 0.0;
    flagNoMuConvergence(true);
  }
}

template<class SoftPars>
double Softsusy<SoftPars>::displayMwRun() const {
  double costhDRbar = cos(asin(calcSinthdrbar()));
  return displayMzRun() * costhDRbar;
}

/// From hep-ph/9311269's appendix. It should be done at MSusy to minimize the
/// 1-loop contributions. Only call if you've calculated physpars
/// inputs are running top/bottom masses. Call at MSusy
template<class SoftPars>
double Softsusy<SoftPars>::doCalcTadpole2oneLoop(double /* mt */, double sinthDRbar) {
/// CHECKED
 if (forLoops.mu(1, 3) == 0.0 || forLoops.mu(2, 3) == 0.0) {
   if (PRINTOUT > 1)
    cout << "Trying to calculate tadpole without having first calculated"
	 << " the DRbar masses.\n";
   return 0.0; 
  }
  
  double g = displayGaugeCoupling(2), 
    costhDRbar = cos(asin(sinthDRbar)), costhDRbar2 = sqr(costhDRbar),
    tanb = displayTanb(), sinb = sin(atan(tanb)), 
    mu = -displaySusyMu(), q = displayMu(),
    mz = displayMzRun();
  const double gmzOcthW =  g * mz / costhDRbar;
  /// Sfermion couplings
  DoubleMatrix lTS2Lr(2, 2),  lBS2Lr(2, 2),  lTauS2Lr(2, 2);
  H2SfSfCouplings(lTS2Lr, lBS2Lr, lTauS2Lr, gmzOcthW, mu, sinb);
  double fermions = doCalcTad2Sfermions(lTS2Lr, lBS2Lr, lTauS2Lr, costhDRbar);
  double sfermions = doCalcTad2fermions(q);
  double higgs = doCalcTad2Higgs(q, costhDRbar2, g, tanb);
  /// Neutralinos
  double neutralinos = doCalcTad2Neutralinos(q, costhDRbar, g, sinb);
  double charginos = doCalcTad2Charginos(q, g, sinb);
  double gaugeBosons = doCalcTad2GaugeBosons(q, costhDRbar2, g, tanb);
  
  double delta = fermions + sfermions + higgs + charginos + neutralinos + 
    gaugeBosons;

  return delta / (16.0 * sqr(PI));
}

template<class SoftPars>
void Softsusy<SoftPars>::calcTadpole2Ms1loop(double mt, double sinthDRbar) {/// CHECKED
  t2OV2Ms1loop = doCalcTadpole2oneLoop(mt, sinthDRbar); 
  if (testNan(t2OV2Ms1loop)) {
    flagNoMuConvergence(true);
    t2OV2Ms1loop = 0.0;
  }
}

/// Apply at scale MSusy: checked 19.12.2000
/// Displays PHYSICAL MZ, ie MZ(q) - piZz^T(q)
/// Fixed pizztMS to resummed version 6/1/13
template<class SoftPars>
double Softsusy<SoftPars>::predMzsq(double & tanb, double muOld, double eps) {

  if (fabs(displayTadpole1Ms()) < EPSTOL && 
      fabs(displayTadpole2Ms()) < EPSTOL) {
    double sinthDRbar = calcSinthdrbar();
    calcDrBarPars(); 
    double mt = forLoops.mt;
    doTadpoles(mt, sinthDRbar);
  }

  double susyMu = displaySusyMu();
  tanb = predTanb(susyMu);
  if (muOld > -6.e66) susyMu = susyMu / eps - muOld * (1. / eps - 1.);

  double pizztMS = sqr(displayMzRun()) - sqr(displayMz()); ///< resums logs
  double MZsq = 2.0 *
    ((displayMh1Squared() - displayTadpole1Ms() - 
      (displayMh2Squared() - displayTadpole2Ms()) *
      sqr(tanb)) / (sqr(tanb) - 1.0) - sqr(susyMu)) - 
    pizztMS;

  return MZsq;
}


/// Used to get useful information into ftCalc
static MssmSoftsusy *tempSoft1;
static int ftFunctionality;
static DoubleVector ftPars(3); 
static void (*ftBoundaryCondition)(MssmSoftsusy &, const DoubleVector &);

inline double ftCalc(double x) {
  /// Stores running parameters in a vector
  DoubleVector storeObject(tempSoft1->display());
  double initialMu = tempSoft1->displayMu();
  drBarPars saveDrBar(tempSoft1->displayDrBarPars());

  if (PRINTOUT > 1) cout << "#";
  
  if (ftFunctionality <= ftPars.displayEnd()) {
    ftPars(ftFunctionality) = x;
    ftBoundaryCondition(*tempSoft1, ftPars);
    if (PRINTOUT > 1) cout << "p" << ftFunctionality << "=";
  }
  else if (ftFunctionality == ftPars.displayEnd() + 1) { 
    tempSoft1->setSusyMu(x); if (PRINTOUT > 1) cout << "mu= "; 
  }
  else if (ftFunctionality == ftPars.displayEnd() + 2) { 
    tempSoft1->setM3Squared(x); if (PRINTOUT > 1) cout << "m3sq= "; 
  } 
  else if (ftFunctionality == ftPars.displayEnd() + 3) {
    tempSoft1->setYukawaElement(YU, 3, 3, x); if (PRINTOUT > 1) cout << "ht= ";
  }
  else {
    std::ostringstream ii;
    ii << "MssmSoftsusy:ftCalc called with incorrect functionality=" <<
      ftFunctionality << endl;
    throw ii.str();
  }
  
  double referenceMzsq, predTanb;
  
  tempSoft1->runto(tempSoft1->calcMs());
  tempSoft1->calcDrBarPars();
  tempSoft1->runto(tempSoft1->calcMs());
  
  tempSoft1->setHvev(tempSoft1->getVev());
  double mt = tempSoft1->displayDrBarPars().mt;
  double sinthDRbar = tempSoft1->calcSinthdrbar();
  /// We miss two-loop terms in our calculation of fine-tuning...
  tempSoft1->calcTadpole2Ms1loop(mt, sinthDRbar); 
  tempSoft1->calcTadpole1Ms1loop(mt, sinthDRbar);
  
  referenceMzsq = tempSoft1->predMzsq(predTanb);  
  
  if (PRINTOUT > 1) cout << x << " MZ=" << sqrt(fabs(referenceMzsq)) 
			 << " tanb=" << predTanb << "\n";

  /// Restore initial parameters at correct scale
  tempSoft1->setMu(initialMu);
  tempSoft1->set(storeObject);
  tempSoft1->setDrBarPars(saveDrBar);

  return referenceMzsq;
}

template<class SoftPars>
double Softsusy<SoftPars>::it1par(int numPar, const DoubleVector & bcPars) {
   throw "it1par not implemented for non-MSSM models";
}

/// Give it a GUT scale object consistent with rewsb
/// and it'll return the fine tuning by varying m32, mu and m3sq at the high
/// scale
template<>
inline double Softsusy<SoftParsMssm>::it1par(int numPar, const DoubleVector & bcPars) {
  
  double ftParameter = 0.0, err, h = 0.01;
  
  tempSoft1 = this;

  /// Stores running parameters in a vector
  DoubleVector storeObject(display());
  double initialMu = displayMu();
  sPhysical savePhys(displayPhys());
  
  ftFunctionality = numPar;
  double x;
  
  /// Defines starting value to calculate derivative from
  if (numPar > 0 && numPar <= bcPars.displayEnd()) {
    x = bcPars.display(numPar); h = 0.01 * x; 
  }
  else if (numPar == bcPars.displayEnd() + 1) {//mu
    x = displaySusyMu(); h = 0.01 * x; 
  }
  else if (numPar == bcPars.displayEnd() + 2) {//m3sq
    x = displayM3Squared(); h = 0.01 * x; 
  }
  else if (numPar == bcPars.displayEnd() + 3) {//ht  
    x = displayYukawaElement(YU, 3, 3); h = x * 0.0005; 
  }
  else {
    std::ostringstream ii;
    ii << "it1par called with functionality " << ftFunctionality << 
      " out of range.\n";
    throw ii.str();
  }
  
  ftPars.setEnd(bcPars.displayEnd()); ftPars = bcPars;
  
  /// Fine tuning is a / MZ^2 d MZ^2 / d a for any parameter a
  double derivative;
  if (fabs(x) < 1.0e-10) {
    ftParameter = 0.0; derivative = 0.0; err = 0.0;
  }
  else {
    derivative = calcDerivative(ftCalc, x, h, &err);
    ftParameter = x * derivative / sqr(displayMz());
  }
  
  if (PRINTOUT > 1)
    cout << "derivative=" << derivative << " error=" << err << "\n";
  
  /// High error: if can't find a derivative, error comes back with 1.0e30
  if (ftParameter > TOLERANCE && fabs(err / derivative) > 1.0) 
    return numberOfTheBeast;

  /// Restore initial parameters at correct scale
  setMu(initialMu);
  set(storeObject);
  setPhys(savePhys);

  return fabs(ftParameter);
}

/// Pass it an object and it'll return the fine tuning parameters
template<>
inline DoubleVector Softsusy<SoftParsMssm>::fineTune
(void (*boundaryCondition)(Softsusy<SoftParsMssm> &, const DoubleVector &),
  const DoubleVector  & bcPars, double mx, bool doTop) {

  /// Stores running parameters in a vector
  DoubleVector storeObject(display());
  double initialMu = displayMu();
  sPhysical savePhys(displayPhys());
  
  runto(mx);

  ftBoundaryCondition = boundaryCondition;

  int numPars = bcPars.displayEnd() + 2;  
  if (doTop) numPars = numPars + 1;
  
  DoubleVector tempFineTuning(numPars);

  int firstParam = 1;
  if (boundaryCondition == &gmsbBcs) firstParam = 2;
  
  int i; for (i = firstParam; i<= numPars; i++) {
    tempFineTuning(i) = it1par(i, bcPars);
    /// flag problem FT calculation with NaN...too inaccurate
    if (tempFineTuning(i) > 1.0e66) tempFineTuning(i) = asin(2.0);
  }
  
  /// Restore initial parameters at correct scale
  setMu(initialMu);
  set(storeObject);
  setPhys(savePhys);

  return tempFineTuning;
}

/// Pass it an object and it'll return the fine tuning parameters
template<class SoftPars>
DoubleVector Softsusy<SoftPars>::fineTune
(void (*boundaryCondition)(Softsusy<SoftPars> &, const DoubleVector &),
 const DoubleVector  & bcPars, double mx, bool doTop) {
   throw "fineTune not implemented for non-MSSM models";
}

/// Obtains solution of one-loop effective potential minimisation via iteration
/// technique
/// err is 1 if no iteration reached
/// 2 if incorrect rewsb

template<class SoftPars>
void Softsusy<SoftPars>::iterateMu(double & muold, int sgnMu,
			     double mt, int maxTries, double pizzMS,
			     double sinthDRbar, double tol, int & err) {
  static int numTries = 0;
  static double munew = 0.;
  
  if (numTries - 1 > maxTries) { 
    if (PRINTOUT) cout << "iterateMu reached maxtries\n"; 
    numTries = 0; munew = 0.0;
    err = 1; return;
  }
  /// How close to convergence are we?
  double c = 1.0 - minimum(fabs(muold), fabs(munew)) / 
    maximum(fabs(muold), fabs(munew));
  
  if (PRINTOUT > 2) cout << " diff=" << c;

  if (c < tol) { 
    muold = munew; //err = 0;
    numTries = 0; munew = 0.0;
    if (PRINTOUT > 2) cout << " mu converged\n";
    return; 
  }

  numTries = numTries + 1;
  
  muold = munew;
  
  double mH1sq = displayMh1Squared(), mH2sq = displayMh2Squared(), 
    tanb = displayTanb(), tanb2 = sqr(tanb);
  
  double treeLevelMusq = (mH1sq - mH2sq * tanb2) / (tanb2 - 1.0) - 0.5 *
    sqr(displayMz());
  
  try {
    calcDrBarPars();
    double oneLoopMusq = 0.;
    /// calculate the new one-loop tadpoles with old value of mu
    if (numRewsbLoops > 0) {  
      doTadpoles(mt, sinthDRbar);
      oneLoopMusq = treeLevelMusq - 0.5 * pizzMS + 
	(displayTadpole2Ms() * sqr(tanb) - displayTadpole1Ms()) /
	(sqr(tanb) - 1.0);
    }

    /// Error in rewsb
    if (oneLoopMusq < 0.0) {
      err = 2; 
      if (PRINTOUT > 1) cout << " mu^2<0 ";
    }  
    
    munew = sgnMu * sqrt(fabs(oneLoopMusq));
    setSusyMu(munew); double m3sqnew;
    if (rewsbM3sq(munew, m3sqnew) == 0) flagM3sq(false);
    else flagM3sq(true);
    setM3Squared(m3sqnew);
  }
  catch(const char *a) {
    numTries = 0;
    throw a;
  }
  catch(const std::string &a) {
    numTries = 0;
    throw a;
  }
  
  if (PRINTOUT > 2) cout << " mu=" << munew;

  iterateMu(muold, sgnMu, mt, maxTries, pizzMS, sinthDRbar, tol, err);
}


template<class SoftPars>
void Softsusy<SoftPars>::alternativeEwsb(double mt) {
  calcDrBarPars();
  double sinthDRbarMS = calcSinthdrbar();
  double tanb = displayTanb(), beta = atan(tanb);
  double mzRun = displayMzRun();
  doTadpoles(mt, sinthDRbarMS);
  
  double piaa = piAA(displayMaCond(), displayMu());
  
  double mAsq =  sqr(displayDrBarPars().mA0(1));
  
  double gstrong = displayGaugeCoupling(3), 
    rmtsq = sqr(displayDrBarPars().mt), scalesq = sqr(displayMu()), 
    vev2 = sqr(displayHvev()), tbeta = displayTanb(), 
    amu = -displaySusyMu(), mg = displayGaugino()(3);
  double newMh1sq = 0., newMh2sq = 0.;

  double p2s = 0., p2w = 0., p2b = 0., p2tau = 0.;
  if (numHiggsMassLoops > 1) {
    /// two-loop Higgs corrections
    double sintau = sin(displayDrBarPars().thetatau), 
      costau = cos(displayDrBarPars().thetatau);
    double msnusq = sqr(displayDrBarPars().msnu.display(3));
    double sxb = sin(displayDrBarPars().thetab), 
      cxb = cos(displayDrBarPars().thetab);
    double msb1sq = sqr(displayDrBarPars().md.display(1, 3)), 
      msb2sq = sqr(displayDrBarPars().md.display(2, 3));
    double mstau1sq = sqr(displayDrBarPars().me.display(1, 3)), 
      mstau2sq = sqr(displayDrBarPars().me.display(2, 3));
    double cotbeta = 1.0 / tbeta;
    double rmbsq = sqr(displayDrBarPars().mb);
    double rmtausq = sqr(displayDrBarPars().mtau);      
    double sxt = sin(displayDrBarPars().thetat), 
      cxt = cos(displayDrBarPars().thetat);
    double mst1sq = sqr(displayDrBarPars().mu.display(1, 3)), 
      mst2sq = sqr(displayDrBarPars().mu.display(2, 3));
      
    dszodd_(&rmtsq, &mg, &mst1sq, &mst2sq, &sxt, &cxt, &scalesq, &amu,
            &tbeta, &vev2, &gstrong, &p2s); 
    ddsodd_(&rmtsq, &rmbsq, &mAsq, &mst1sq, &mst2sq, &msb1sq, &msb2sq, 
	    &sxt, &cxt, &sxb, &cxb, &scalesq, &amu, &tanb, &vev2, 
	    &p2w);

    dszodd_(&rmbsq, &mg, &msb1sq, &msb2sq, &sxb, &cxb, &scalesq, &amu,
            &cotbeta, &vev2, &gstrong, &p2b); 
    tausqodd_(&rmtausq, &mAsq, &msnusq, &mstau1sq, &mstau2sq, &sintau,
	      &costau, &scalesq, &amu, &tanb, &vev2, &p2tau);
    
    double dMA = p2s + p2b + p2w + p2tau;
    
    newMh1sq = sqr(sin(beta)) * (sqr(displayMaCond()) + piaa + 
				 sqr(mzRun) - dMA) 
      - (sqr(displaySusyMu()) + 0.5 * sqr(mzRun)) +
      displayTadpole1Ms() - sqr(sqr(sin(beta))) * 
      t1OV1Ms1loop -
      t2OV2Ms1loop * sqr(sin(beta)) * sqr(cos(beta));
    
    newMh2sq = sqr(cos(beta)) * (sqr(displayMaCond()) + piaa + 
				 sqr(mzRun) - dMA) 
      - (sqr(displaySusyMu()) + 0.5 * sqr(mzRun)) -
      t1OV1Ms1loop * sqr(sin(beta)) * sqr(cos(beta)) +
      displayTadpole2Ms() - sqr(sqr(cos(beta))) * 
      t2OV2Ms1loop;
  }
  setMh1Squared(newMh1sq);
  setMh2Squared(newMh2sq);
  
  double m3sqnew = 0.;
  if (rewsbM3sq(displaySusyMu(), m3sqnew) == 0) flagM3sq(false);
  else flagM3sq(true);
  setM3Squared(m3sqnew);
  
  if ((displayMh1Squared() + 2.0 * sqr(displaySusyMu()) +
       displayMh2Squared() - 2.0 * fabs(displayM3Squared())) < 0.0 )
    flagHiggsufb(true);
  else 
    flagHiggsufb(false);
}

template<class SoftPars>
double Softsusy<SoftPars>::treeLevelMuSq() {
  return (displayMh1Squared() - displayMh2Squared() * sqr(displayTanb())) / 
    (sqr(displayTanb()) - 1.0) - 0.5 * sqr(displayMz());
}

template<class SoftPars>
void Softsusy<SoftPars>::rewsbTreeLevel(int sgnMu) {
  if (altEwsb) {
    setSusyMu(displayMuCond());
    double newMh1sq, newMh2sq;
    double beta = atan(displayTanb());
    newMh1sq = sqr(sin(beta)) * (sqr(displayMaCond()) + sqr(MZ)) 
      - (sqr(displaySusyMu()) + 0.5 * sqr(MZ));
    
    newMh2sq = sqr(cos(beta)) * (sqr(displayMaCond()) + sqr(MZ)) 
      - (sqr(displaySusyMu()) + 0.5 * sqr(MZ));
    
    setMh1Squared(newMh1sq);
    setMh2Squared(newMh2sq);
    
    double m3sq;
    if (rewsbM3sq(displaySusyMu(), m3sq)) flagM3sq(true);
    else flagM3sq(false); 
    
    setM3Squared(m3sq);
    return;
  } else {
    double mu, m3sq;
    if (rewsbMu(sgnMu, mu)) flagMusqwrongsign(true);
    else flagMusqwrongsign(false);
    setSusyMu(mu);
    if (rewsbM3sq(mu, m3sq)) flagM3sq(true);
    else flagM3sq(false); 
    
    setM3Squared(m3sq);
    
    if ((displayMh1Squared() + 2.0 * sqr(displaySusyMu()) +
	 displayMh2Squared() - 2.0 * fabs(displayM3Squared())) < 0.0 )
      flagHiggsufb(true);
    else 
      flagHiggsufb(false);
    
    return;
  }
}

/// Organises rewsb: call it at the low scale MS^2 = sqrt(0.5 * (mT1^2 +
/// mT2^2)) is best, or below if it's decoupled from there. 
/// Call with zero, or no mt if you want tree level
template<class SoftPars>
void Softsusy<SoftPars>::rewsb(int sgnMu, double mt, const DoubleVector & /* pars */,
			 double muOld, double eps) {
  if (altEwsb) {
    alternativeEwsb(mt);
    return;
  } else {
    double munew, m3sqnew;
    
    double sinthDRbarMS = calcSinthdrbar();
    
    calcTadpole1Ms1loop(mt, sinthDRbarMS);  
    calcTadpole2Ms1loop(mt, sinthDRbarMS); 
    
    munew = displaySusyMu();
    
    /// Iterate to get a self-consistent mu solution
    int maxTries = 20, err = 0;
    double tol = TOLERANCE * 1.0e-6;
    
    double pizztMS = sqr(displayMzRun()) - sqr(displayMz()); ///< resums logs
    
    iterateMu(munew, sgnMu, mt, maxTries, pizztMS, sinthDRbarMS,
	      tol, err); 
    
    if (err == 2) {
      flagMusqwrongsign(true);
      if (PRINTOUT > 2) cout << " mu^2<0 ";
    }
    else flagMusqwrongsign(false); 
    if (err == 1) {
      flagNoMuConvergence(true);
      if (PRINTOUT > 2) cout << " no mu convergence ";
    }
    else setSusyMu(munew);
    
    /// average mu with the input value of muOld, if it isn't the number of the
    /// beast   
    if (muOld > -6.e66) {
      munew = (munew * eps + muOld * (1. - eps));
      setSusyMu(munew);
    }
    
    if (rewsbM3sq(munew, m3sqnew) == 0) flagM3sq(false);
    else flagM3sq(true);

    setM3Squared(m3sqnew);

    if ((displayMh1Squared() + 2.0 * sqr(displaySusyMu()) +
	 displayMh2Squared() - 2.0 * fabs(displayM3Squared())) < 0.0 )
      flagHiggsufb(true);
    else 
      flagHiggsufb(false);
  }
}


/// Gives a summary of important properties of a SUSY object: mu, m3sq, mH0,
/// mChi0, lightest stau, lightest mGluino, lightest stop, lightest chargino,
/// lightest tau sneutrino, lightest sbottom, down squark, up squark and
/// selectron masses, minimum of potential (if calculated) and fine-tuning
/// parameter (if passed)
inline void printShortInitialise() {
  cout <<
    "     mu     " << "   m3sq     " << 
    "   mstau1   " << "   msbott   " << "  mstop1    " <<
    " mSquark   " << "   mSlep   " << "  msnu1     " <<
    "    mhl0    " << "   mchi10   " << " mchargino1 " <<
    "  mGluino  " << std::flush;
  }

/// Prints mu, B and important spectral information
template<class SoftPars>
std::string Softsusy<SoftPars>::printShort() const {
  
  std::ostringstream a;
  const double problemFlag = -1.0;

  double mu, m3sq;
  mu = displaySusyMu(); m3sq = displayM3Squared(); 
  int pos;
  sPhysical s(this->displayPhys());
  
  if (displayProblem().muSqWrongSign || displayProblem().m3sq ||
      displayProblem().higgsUfb || displayProblem().noRhoConvergence) 
    mu = -1.0 *  mu;
  if (displayProblem().tachyon) m3sq = -1.0 * m3sq;
  
  if (displayProblem().nonperturbative || displayProblem().noMuConvergence) {
    a << " ";
    int i; for (i=1; i<=12; i++) a << problemFlag << " ";
  }
  else {
    a << " " << mu << " " << m3sq << " "
	   << minimum(s.me(1, 3), s.me(2, 3)) << " " 
	   << minimum(s.md(1, 3), s.md(2, 3)) << " " 
	   << minimum(s.mu(1, 3), s.mu(2, 3)) << " " 
	 << minimum(minimum(s.md(1, 1), s.md(2, 1)), 
		    minimum(s.mu(1, 1), s.mu(2, 1)))
	   << " " << minimum (s.me(1, 1), s.me(2, 1)) << " "
	   << s.msnu.min(pos) << " " 
	   << s.mh0(1) << " " 
	   << (s.mneut.apply(fabs)).min(pos) << " "
	   << (s.mch.apply(fabs)).min(pos) << " "
	   << s.mGluino << " ";
    }
  
  a << std::flush;
  return a.str();
}

template<class SoftPars>
std::string Softsusy<SoftPars>::printLong() {
  /// output:
  ///  1  2     3      4   5   6   7    8  9  10 11   12    13    
  /// mu  m3sq mH1sq mH2sq g1 g2 mt(mt) mh mA mH mH+ alphaH msnu3
  ///  14     15     16     17    18     19    20   
  /// msnu1 mstopL mstopR msupL msupR msbotL msbotR
  ///  21   22   23     24     25     26     27    28     29     
  /// msdL msdR mstauL mstauR mselL mselR thetat thetab thetatau
  /// 30   31   32   33   34      35      36     37     38    
  /// mgl mch1 mch2 thetaL thetaR mneut1 mneut2 mneut3 mneut4
  ///   39    40     41 
  /// sinthW t1ov1 t2ov2 
  std::ostringstream a;
  double mu = displaySusyMu();
  if (displayProblem().muSqWrongSign || displayProblem().m3sq ||
      displayProblem().higgsUfb) 
    mu = -1.0 *  mu;
  
  a << " " << mu << " " 
       << displayM3Squared() << " " 
       << displayMh1Squared() << " " <<
    displayMh2Squared() << " " << 
    displayGaugeCoupling(1) << " " <<
    displayGaugeCoupling(2) << " " <<
    calcRunningMt() << " " <<
    displayPhys().mh0(1) << " " <<
    displayPhys().mA0(1) << " " <<
    displayPhys().mh0(2) << " " <<
    displayPhys().mHpm << " " <<
    displayPhys().thetaH << " " <<
    displayPhys().msnu.display(3) << " " <<
    displayPhys().msnu.display(1) << " " <<
    displayPhys().mu.display(1, 3) << " " <<
    displayPhys().mu.display(2, 3) << " " <<
    displayPhys().mu.display(1, 1) << " " <<
    displayPhys().mu.display(2, 1) << " " <<
    displayPhys().md.display(1, 3) << " " <<
    displayPhys().md.display(2, 3) << " " <<
    displayPhys().md.display(1, 1) << " " <<
    displayPhys().md.display(2, 1) << " " <<
    displayPhys().me.display(1, 3) << " " <<
    displayPhys().me.display(2, 3) << " " <<
    displayPhys().me.display(1, 1) << " " <<
    displayPhys().me.display(2, 1) << " " <<
    displayPhys().thetat << " " <<
    displayPhys().thetab << " " <<
    displayPhys().thetatau << " " <<
    displayPhys().mGluino << " " <<
    displayPhys().mch.display(1) << " " <<
    displayPhys().mch.display(2) << " " <<
    displayPhys().thetaL << " " <<
    displayPhys().thetaR << " " <<
    displayPhys().mneut.display(1) << " " <<
    displayPhys().mneut.display(2) << " " <<
    displayPhys().mneut.display(3) << " " <<
    displayPhys().mneut.display(4) << " " <<
    calcSinthdrbar() << " " <<
    displayTadpole1Ms() << " " <<
    displayTadpole2Ms() << " ";

  if (displayProblem().test()) a << "%" << displayProblem();
  a << std::flush;
  return a.str();
}



template<class SoftPars>
bool Softsusy<SoftPars>::higgs(int accuracy, double piwwtMS,
                               double /* pizztMS */) {

  double tanb = displayTanb();
  double beta = atan(tanb);
  double sinb = sin(beta), cosb = cos(beta);
  double sinb2 = sqr(sinb), cosb2 = sqr(cosb), mzPole = displayMz(), 
    mzRun2 = sqr(displayMzRun());
  double mApole = physpars.mA0(1); /// physical value
  ///  double mApole2 = sqr(mApole);

  /// There'll be trouble if B has the opp sign to mu. This isn't really
  /// tree-level since it includes some one loops correrctions, namely
  /// tadpoles and Z loops 
  /*  double mAsq =  (displayMh2Squared() - displayTadpole2Ms() - 
		  displayMh1Squared() + displayTadpole1Ms()) / 
		  cos(2.0 * beta) - mzRun2; */
  double mAsq = displayM3Squared() / (sinb * cosb);

  DoubleMatrix mHtree(2, 2);
  mHtree(1, 1) = mAsq * sinb2 + mzRun2 * cosb2;
  mHtree(1, 2) = - sinb * cosb * (mAsq + mzRun2); 
  mHtree(2, 2) = mAsq * cosb2 + mzRun2 * sinb2; 
  mHtree(2, 1) = mHtree(1 ,2); 

  DoubleMatrix mhAtmh(mHtree), mhAtmH(mHtree);
  DoubleMatrix sigmaMh(2, 2), sigmaMH(2, 2);

  double q = displayMu(), p;     
  double p2s = 0.0, p2w = 0.0, p2b = 0.0, p2tau = 0.0, dMA = 0.;
  /// radiative corrections:
  if (accuracy > 0) {
    /// one-loop radiative corrections included in sigma
    p = physpars.mh0(1);
    sigmaMh(1, 1) = pis1s1(p, q); 
    sigmaMh(1, 2) = pis1s2(p, q); 
    sigmaMh(2, 2) = pis2s2(p, q); 

    p = physpars.mh0(2);
    sigmaMH(1, 1) = pis1s1(p, q); 
    sigmaMH(1, 2) = pis1s2(p, q); 
    sigmaMH(2, 2) = pis2s2(p, q);

    if (numHiggsMassLoops > 1) {
      double s11s = 0., s22s = 0., s12s = 0., 
	gstrong = displayGaugeCoupling(3), 
	rmtsq = sqr(forLoops.mt), scalesq = sqr(displayMu()), 
	vev2 = sqr(displayHvev()), tbeta = displayTanb(), 
	amu = -displaySusyMu(), mg = displayGaugino(3);
      
      double sxt = sin(forLoops.thetat), 
	cxt = cos(forLoops.thetat);
      double mst1sq = sqr(forLoops.mu(1, 3)), 
	mst2sq = sqr(forLoops.mu(2, 3));
      /// two-loop Higgs corrections: alpha_s alpha_b
      double sxb = sin(forLoops.thetab), 
      cxb = cos(forLoops.thetab);
      double msb1sq = sqr(forLoops.md(1, 3)), 
	msb2sq = sqr(forLoops.md(2, 3));
      double cotbeta = 1.0 / tbeta;
      double rmbsq = sqr(forLoops.mb);

      double s11b = 0.0, s12b = 0.0, s22b = 0.0;
      double s11tau = 0.0, s12tau = 0.0, s22tau = 0.0;
      double s11w = 0.0, s22w = 0.0, s12w = 0.0;
      int kkk = 0; /// chooses DR-bar scheme from slavich et al

      double fmasq = fabs(mAsq);
      /// two-loop Higgs corrections: alpha_s alpha_t, alpha_s alpha_b and
      /// alpha_b^2, alpha_t*2, alpha_b alpha_t
      dszhiggs_(&rmtsq, &mg, &mst1sq, &mst2sq, &sxt, &cxt, &scalesq, &amu, 
        	&tbeta, &vev2, &gstrong, &kkk, &s11s, &s22s, &s12s);
      dszodd_(&rmtsq, &mg, &mst1sq, &mst2sq, &sxt, &cxt, &scalesq, &amu,
              &tbeta, &vev2, &gstrong, &p2s); 
      dszhiggs_(&rmbsq, &mg, &msb1sq, &msb2sq, &sxb, &cxb, &scalesq, &amu, 
        	&cotbeta, &vev2, &gstrong, &kkk, &s22b, &s11b, &s12b);
      dszodd_(&rmbsq, &mg, &msb1sq, &msb2sq, &sxb, &cxb, &scalesq, &amu,
              &cotbeta, &vev2, &gstrong, &p2b); 
      ddshiggs_(&rmtsq, &rmbsq, &fmasq, &mst1sq, &mst2sq, &msb1sq, &msb2sq, 
              &sxt, &cxt, &sxb, &cxb, &scalesq, &amu, &tanb, &vev2, &s11w, 
      	      &s12w, &s22w);
      ddsodd_(&rmtsq, &rmbsq, &fmasq, &mst1sq, &mst2sq, &msb1sq, &msb2sq, 
              &sxt, &cxt, &sxb, &cxb, &scalesq, &amu, &tanb, &vev2, 
      	      &p2w);
       
      /// In hep-ph/0406277 found the lambda_tau^2 and lambda_tau lambda_b
      /// corrections to be negligible. Have commented them here for
      /// calculational efficiency. Uncomment to see their effect.
      double sintau = sin(forLoops.thetatau),
	costau = cos(forLoops.thetatau);
      double rmtausq = sqr(forLoops.mtau);
      int OS = 0;
      double mstau1sq = sqr(forLoops.me(1, 3)), 
	mstau2sq = sqr(forLoops.me(2, 3));
      double msnusq = sqr(forLoops.msnu(3));
      tausqhiggs_(&rmtausq, &fmasq, &msnusq, &mstau1sq, &mstau2sq, &sintau,
        	  &costau, &scalesq, &amu, &tanb, &vev2, &OS, &s11tau, 
        	  &s22tau, &s12tau);
      tausqodd_(&rmtausq, &fmasq, &msnusq, &mstau1sq, &mstau2sq, &sintau,
        	&costau, &scalesq, &amu, &tanb, &vev2, &p2tau);
      
      sigmaMh(1, 1) = sigmaMh(1, 1) - s11s - s11w - s11b - s11tau;
      sigmaMH(1, 1) = sigmaMH(1, 1) - s11s - s11w - s11b - s11tau;
      sigmaMh(1, 2) = sigmaMh(1, 2) - s12s - s12w - s12b - s12tau;
      sigmaMH(1, 2) = sigmaMH(1, 2) - s12s - s12w - s12b - s12tau;
      sigmaMh(2, 2) = sigmaMh(2, 2) - s22s - s22w - s22b - s22tau;
      sigmaMH(2, 2) = sigmaMH(2, 2) - s22s - s22w - s22b - s22tau;
 
    }

    sigmaMh(2, 1) = sigmaMh(1, 2);
    sigmaMH(2, 1) = sigmaMH(1, 2);
    /*
      As the calculation stands without the two-loop terms, BPMZ have
      obviously organised PI_Sij (CP-even loop corrections) so that their pole
      masses come out correctly to one-loop order, hence there is no need to
      add any one-loop self energies of the PI^AA terms to DMA. Secondly, the
      one loop terms involving mA in the PI_Sij are not proportional to either
      alpha_s, alpha_t or alpha_b. So, within PI_Sij(1 loop), if I consider
      the 2-loop order difference due to not using mA(pole), it will be (some
      other coupling) x (either alpha_s OR alpha_t OR alpha_b) ie it will not
      be of order alpha_s alpha_t, alpha_b alpha_s or alpha_t^2, thus
      remaining consistent with your two-loop terms. They have also performed
      the calculation where their stuff "adds on" to the 1-loop BPMZ stuff. 
      The tadpoles therefore appear only as 1-loop pieces, except where
      minimisation conditions are explicitly used (in the calculation of mA)
     */
    dMA = p2s + p2w + p2b + p2tau;
    mhAtmh(1, 1) = mHtree(1, 1) + t1OV1Ms1loop + dMA * sqr(sin(beta));
    mhAtmh(1, 2) = mHtree(1, 2) - dMA * sin(beta) * cos(beta);
    mhAtmh(2, 2) = mHtree(2, 2) + t2OV2Ms1loop + dMA * sqr(cos(beta));
    mhAtmh(2, 1) = mhAtmh(1 ,2);
    mhAtmh = mhAtmh - sigmaMh;

    mhAtmH(1, 1) = mHtree(1, 1) + t1OV1Ms1loop + dMA * sqr(sin(beta));
    mhAtmH(1, 2) = mHtree(1, 2) - dMA * sin(beta) * cos(beta);
    mhAtmH(2, 2) = mHtree(2, 2) + t2OV2Ms1loop + dMA * sqr(cos(beta));
    mhAtmH(2, 1) = mhAtmH(1 ,2);
    mhAtmH = mhAtmH - sigmaMH;
    
  }
  

  DoubleVector temp(2);  
  double theta;

  temp = mhAtmh.sym2by2(theta);

  bool h0Htachyon = false;
  if (temp(1) < 0.0 || temp(2) < 0.0) {
    h0Htachyon = true;
    if (PRINTOUT > 2) cout << " h0/H tachyon: m^2=" << temp;
  }
  temp = temp.apply(ccbSqrt);

  /// If certain DRbar ratios are large, they can cause massive higher order
  /// corrections in the higgs mass, making it have O(1) uncertainties. 
  /// In these cases, you should switch to an OS calculation (eg by using
  /// FEYNHIGGS) for the Higgs mass (but there are other points at high
  /// tan beta where the DRbar scheme does better).
  double mstop1 = minimum(displayDrBarPars().mu(1, 3), 
			  displayDrBarPars().mu(2, 3));
  double mgluino = displayGaugino(3);
  if (sqr(mgluino / mstop1) > (16.0 * sqr(PI)) ||
      sqr(displayMu() / mstop1) > (16.0 * sqr(PI))) 
    flagInaccurateHiggsMass(true);
    
  /// Definitions are such that theta should diagonalise the matrix like
  /// O = [ cos  sin ]
  ///     [-sin  cos ] 
  /// and
  /// O m O^T = [ m2^2      ]
  ///           [      m1^2 ]
  /// where m1 < m2, therefore if they come out in the wrong order, 
  /// add pi/2 onto theta.
  if (temp(2) > temp(1)) theta = theta + PI * 0.5; 

  physpars.thetaH = theta; /// theta defined for p=mh  
  int i; double littleMh = temp.apply(fabs).min(i);

  temp = mhAtmH.sym2by2(theta);

  if (temp(1) < 0.0 && temp(2) < 0.0) {
    h0Htachyon = true;
    if (PRINTOUT > 2) cout << " h0/H tachyon: m^2=" << temp;
  }
  temp = temp.apply(ccbSqrt);
  double bigMh = temp.max();


double piaa = piAA(mApole, displayMu()); 
  //  double piaa = piAA(displayDrBarPars().mA0, displayMu());
  double poleMasq = (displayMh2Squared() - displayMh1Squared() )
    / cos(2.0 * beta) - sqr(mzPole);
  
  if (accuracy > 0) {
      poleMasq = 
	(displayMh2Squared() - displayTadpole2Ms() - 
	 displayMh1Squared() + displayTadpole1Ms()) / 
	cos(2.0 * beta) - mzRun2 - piaa +
	sqr(sin(beta)) * t1OV1Ms1loop + sqr(cos(beta)) *
	t2OV2Ms1loop + dMA;
    }

  double pihphm = piHpHm(physpars.mHpm, displayMu());

  double poleMhcSq = poleMasq + sqr(displayMw()) + piaa + piwwtMS - pihphm;


  //PA:  below is the rearranged approach to calculating charged and cpodd Higgs used for comaprisons vs nmssmsoftsusy.cpp.   This leaves the initial accuracy == 0 calculation changed and while the actual calculation for accuarcy > 0 is unchanged this alters the final result for charged Higgs due to the numerical precion of the Higgs iteration. 
  /* ------------------------------------------------------------------------*/
  // double piaa = 0.0; 
  // if(accuracy > 0) piaa = piAA(mApole, displayMu()); 
  // //  double piaa = piAA(displayDrBarPars().mA0(1), displayMu());
  // double poleMasq = (displayMh2Squared()-  displayTadpole2Ms() 
  // 		     - displayMh1Squared()  + displayTadpole1Ms())
  //   / cos(2.0 * beta) -  sqr(displayMz());

  // if (accuracy > 0) {
  //   poleMasq = 
  //     (displayMh2Squared() - displayTadpole2Ms() - 
  //      displayMh1Squared() + displayTadpole1Ms()) / 
  //     cos(2.0 * beta) - mzRun2 - piaa +
  //       sqr(sin(beta)) * t1OV1Ms1loop + sqr(cos(beta)) *
  //     t2OV2Ms1loop + dMA;
  // }
  
  // double pihphm = 0.0;
  // if(accuracy > 0) pihphm = piHpHm(physpars.mHpm, displayMu()); 
  // // double poleMhcSq = poleMasq + sqr(displayMw()) + piaa + piwwtMS - pihphm;
  // double poleMasq2 = poleMasq;
  // if (accuracy > 0) {
  //   poleMasq2 = 
  //     (displayMh2Squared() - displayTadpole2Ms() - 
  //      displayMh1Squared() + displayTadpole1Ms()) / 
  //     cos(2.0 * beta) - mzRun2 - piaa +
  //       sqr(sin(beta)) * t1OV1Ms1loop + sqr(cos(beta)) *
  //     t2OV2Ms1loop;
  // }

  // double poleMhcSq = poleMasq2 + sqr(displayMw()) + piaa + piwwtMS - pihphm;
  /*-------------------------------------------------------------------------*/

  physpars.mh0(1) = littleMh;
  physpars.mA0(1) = ccbSqrt(poleMasq);
  physpars.mh0(2) = bigMh;
  physpars.mHpm = ccbSqrt(poleMhcSq);

  if (poleMhcSq > 0. && poleMasq > 0. && !h0Htachyon) return false;
  else {
    if (PRINTOUT) cout << " mA(phys)^2=" << poleMasq  
		       << " mHc(phys)^2=" << poleMhcSq 
		       << " but may be first iteration" << endl;
    return true;
  }
}
//PA: adds sfermion contribitions to the left right and scalar parts 
//of the self energy
template<class SoftPars>
void Softsusy<SoftPars>::addChaLoopSfermion(double p, DoubleMatrix & sigmaL, DoubleMatrix & sigmaR, DoubleMatrix & sigmaS) const {
  double q = displayMu();
  double g = displayGaugeCoupling(2), 
    ht = displayDrBarPars().ht,
    htau = displayDrBarPars().htau,
    hb = displayDrBarPars().hb;
  double tanb = displayTanb(), beta = atan(tanb);
   double mt = ht * displayHvev() * sin(beta) / root2,
     mb = hb * displayHvev() * cos(beta) / root2,
     mtau = htau * displayHvev() * cos(beta) / root2;
 /// basis is (psi1+ psi2+, L R)
  DoubleMatrix aPsicDu(2, 2), aPsicUd(2, 2), aPsicENu(2, 2), aPsicNue(2, 2);
  DoubleMatrix aPsicBt(2, 2), aPsicTb(2, 2), aPsicTauNu(2, 2), 
    aPsicNuTau(2, 2);
  DoubleMatrix aPsicBtm(2, 2), aPsicTbm(2, 2), aPsicTauNum(2, 2), 
    aPsicNuTaum(2, 2);

  DoubleMatrix bPsicDu(2, 2), bPsicUd(2, 2), bPsicENu(2, 2), bPsicNue(2, 2);
  DoubleMatrix bPsicBt(2, 2), bPsicTb(2, 2), bPsicTauNu(2, 2), 
    bPsicNuTau(2, 2);
  DoubleMatrix bPsicBtm(2, 2), bPsicTbm(2, 2), bPsicTauNum(2, 2), 
    bPsicNuTaum(2, 2);

  ComplexMatrix aChicDu(2, 2), aChicUd(2, 2), aChicENu(2, 2), aChicNue(2, 2);
  ComplexMatrix aChicBt(2, 2), aChicTb(2, 2), aChicTauNu(2, 2), 
    aChicNuTau(2, 2);

  ComplexMatrix bChicDu(2, 2), bChicUd(2, 2), bChicENu(2, 2), bChicNue(2, 2);
  ComplexMatrix bChicBt(2, 2), bChicTb(2, 2), bChicTauNu(2, 2), 
    bChicNuTau(2, 2);

  aPsicDu(1, 1) = g;       aPsicENu(1, 1) = g;
  bPsicUd(1, 1) = g;       bPsicNue(1, 1) = g;
  aPsicBt(1, 1) = g;       aPsicTauNu(1, 1) = g;
  bPsicTb(1, 1) = g;       bPsicNuTau(1, 1) = g;
  bPsicBt(2, 1) = - hb;    bPsicTauNu(2, 1) = - htau;
  bPsicTb(2, 2) = - hb;    bPsicNuTau(2, 2) = - htau;
  aPsicTb(2, 1) = - ht;
  aPsicBt(2, 2) = - ht;	     
 
  /// mix 3rd generation sfermions
  DoubleMatrix O(2, 2); DoubleVector t1(2), t2(2), tt(2);
  O = rot2d(forLoops.thetat); int i;
  for (i=1; i<=2; i++) {
    tt(1) = aPsicBt(i, 1); tt(2) = aPsicBt(i, 2);
    t1 = O * tt;
    tt(1) = bPsicBt(i, 1); tt(2) = bPsicBt(i, 2);
    t2 = O * tt;
    aPsicBtm(i, 1) = t1(1);     aPsicBtm(i, 2) = t1(2); 
    bPsicBtm(i, 1) = t2(1);     bPsicBtm(i, 2) = t2(2); 
  }

  O = rot2d(forLoops.thetab);
  for (i=1; i<=2; i++) {
    tt(1) = aPsicTb(i, 1); tt(2) = aPsicTb(i, 2);
    t1 = O * tt;
    tt(1) = bPsicTb(i, 1); tt(2) = bPsicTb(i, 2);
    t2 = O * tt;
    aPsicTbm(i, 1) = t1(1);     aPsicTbm(i, 2) = t1(2); 
    bPsicTbm(i, 1) = t2(1);     bPsicTbm(i, 2) = t2(2); 
  }

  O = rot2d(forLoops.thetatau);
  for (i=1; i<=2; i++) {
    tt(1) = aPsicNuTau(i, 1); tt(2) = aPsicNuTau(i, 2);
    t1 = O * tt;
    tt(1) = bPsicNuTau(i, 1); tt(2) = bPsicNuTau(i, 2);
    t2 = O * tt;
    aPsicNuTaum(i, 1) = t1(1);     aPsicNuTaum(i, 2) = t1(2); 
    bPsicNuTaum(i, 1) = t2(1);     bPsicNuTaum(i, 2) = t2(2); 
  }

 DoubleVector msup(2), msdown(2), msel(2), mscharm(2), msstrange(2), 
    msmuon(2), msnumu(2), mstop(2), msbot(2), mstau(2), msnutau(2), msnue(2);
  msup(1) = forLoops.mu(1, 1);      msup(2) = forLoops.mu(2, 1); 
  mscharm(1) = forLoops.mu(1, 2);   mscharm(2) = forLoops.mu(2, 2); 
  mstop(1) = forLoops.mu(1, 3);     mstop(2) = forLoops.mu(2, 3); 
  msdown(1) = forLoops.md(1, 1);    msdown(2) = forLoops.md(2, 1); 
  msstrange(1) = forLoops.md(1, 2); msstrange(2) = forLoops.md(2, 2); 
  msbot(1) = forLoops.md(1, 3);     msbot(2) = forLoops.md(2, 3); 
  msel(1) = forLoops.me(1, 1);      msel(2) = forLoops.me(2, 1); 
  msmuon(1) = forLoops.me(1, 2);    msmuon(2) = forLoops.me(2, 2); 
  mstau(1) = forLoops.me(1, 3);     mstau(2) = forLoops.me(2, 3); 
  msnue(1) = forLoops.msnu(1); 
  msnumu(1) = forLoops.msnu(2);
  msnutau(1) = forLoops.msnu(3);
 for (int i=1; i<=2; i++) 
    for (int j=1; j<=2; j++) 
      for (int k=1; k<=2; k++) {
	sigmaL(i, j) = sigmaL(i, j) + 0.5 * 
	  (3.0 * aPsicDu(i, k) * aPsicDu(j, k) * b1(p, 0., msup(k), q) +
	   3.0 * aPsicUd(i, k) * aPsicUd(j, k) * b1(p, 0., msdown(k), q) +
	   aPsicNue(i, k) * aPsicNue(j, k) * b1(p, 0., msel(k), q) +
	   aPsicENu(i, k) * aPsicENu(j, k) * b1(p, 0., msnue(k), q));
	sigmaL(i, j) = sigmaL(i, j) + 0.5 * 
	  (3.0 * aPsicDu(i, k) * aPsicDu(j, k) * b1(p, 0., mscharm(k), q) +
	   3.0 * aPsicUd(i, k) * aPsicUd(j, k) * b1(p, 0., msstrange(k), q) +
	   aPsicNue(i, k) * aPsicNue(j, k) * b1(p, 0., msmuon(k), q) +
	   aPsicENu(i, k) * aPsicENu(j, k) * b1(p, 0., msnumu(k), q));
	sigmaL(i, j) = sigmaL(i, j) + 0.5 * 
	  (3.0 * aPsicBtm(i, k) * aPsicBtm(j, k) * b1(p, mb, mstop(k), q) +
	   3.0 * aPsicTbm(i, k) * aPsicTbm(j, k) * b1(p, mt, msbot(k), q) +
	   aPsicNuTaum(i, k) * aPsicNuTaum(j, k) * b1(p, 0., mstau(k), q) +
	   aPsicTauNu(i, k) * aPsicTauNu(j, k) * b1(p, mtau, msnutau(k), q));
	sigmaR(i, j) = sigmaR(i, j) + 0.5 * 
	  (3.0 * bPsicDu(i, k) * bPsicDu(j, k) * b1(p, 0., msup(k), q) +
	   3.0 * bPsicUd(i, k) * bPsicUd(j, k) * b1(p, 0., msdown(k), q) +
	   bPsicNue(i, k) * bPsicNue(j, k) * b1(p, 0., msel(k), q) +
	   bPsicENu(i, k) * bPsicENu(j, k) * b1(p, 0., msnue(k), q));
	sigmaR(i, j) = sigmaR(i, j) + 0.5 * 
	  (3.0 * bPsicDu(i, k) * bPsicDu(j, k) * b1(p, 0., mscharm(k), q) +
	   3.0 * bPsicUd(i, k) * bPsicUd(j, k) * b1(p, 0., msstrange(k), q) +
	   bPsicNue(i, k) * bPsicNue(j, k) * b1(p, 0., msmuon(k), q) +
	   bPsicENu(i, k) * bPsicENu(j, k) * b1(p, 0., msnumu(k), q));
	sigmaR(i, j) = sigmaR(i, j) + 0.5 * 
	  (3.0 * bPsicBtm(i, k) * bPsicBtm(j, k) * b1(p, mb, mstop(k), q) +
	   3.0 * bPsicTbm(i, k) * bPsicTbm(j, k) * b1(p, mt, msbot(k), q) +
	   bPsicNuTaum(i, k) * bPsicNuTaum(j, k) * b1(p, 0., mstau(k), q) +
	   bPsicTauNu(i, k) * bPsicTauNu(j, k) * b1(p, mtau, msnutau(k), q));
	sigmaS(i, j) = sigmaS(i, j) + 
	  (3.0 * bPsicBtm(i, k) * aPsicBtm(j, k) * mb * b0(p, mb, mstop(k), q) +
	   3.0 * bPsicTbm(i, k) * aPsicTbm(j, k) * mt * b0(p, mt, msbot(k), q) +
	   bPsicTauNu(i, k) * aPsicTauNu(j, k) * mtau *
	   b0(p, mtau, msnutau(k), q));
      }
}

//PA: adds gauge boson contribitions to the left right and scalar parts 
//of the chargino self energy
template<class SoftPars>
void Softsusy<SoftPars>::addChaLoopGauge(double /* p */, DoubleMatrix & sigmaL, DoubleMatrix & sigmaR, DoubleMatrix & sigmaS, DoubleMatrix b1pCha, DoubleMatrix b0pCha, DoubleMatrix b1pNeut, DoubleMatrix b0pNeut) const {
  double g = displayGaugeCoupling(2);
  double e = g * calcSinthdrbar();
  
  ComplexMatrix n(forLoops.nBpmz);
  DoubleVector mneut(forLoops.mnBpmz);
  ComplexMatrix u(forLoops.uBpmz), v(forLoops.vBpmz); 
  DoubleVector mch(forLoops.mchBpmz); 
  const int dimN = mneut.displayEnd();

  /// checked and corrected
  ComplexMatrix aPsi0PsicW(dimN, 2), bPsi0PsicW(dimN, 2), aChi0PsicW(dimN, 2),
    bChi0PsicW(dimN, 2);
  aPsi0PsicW(2, 1) = - g;
  bPsi0PsicW(2, 1) = - g;
  aPsi0PsicW(4, 2) = g / root2;		     
  bPsi0PsicW(3, 2) = -g / root2;		     
  aChi0PsicW = n * aPsi0PsicW;
  bChi0PsicW = n.complexConjugate() * bPsi0PsicW;
  
  /// checked 
  ComplexMatrix aPsiPsiZ(2, 2), bPsiPsiZ(2, 2), aPsiChiZ(2, 2), bPsiChiZ(2, 2);
  double sinthW = calcSinthdrbar();
  double thetaWDRbar = asin(sinthW);
  aPsiPsiZ(1, 1) = g * cos(thetaWDRbar);
  aPsiPsiZ(2, 2) = g * cos(2.0 * thetaWDRbar) / (2.0 * cos(thetaWDRbar));
  bPsiPsiZ = aPsiPsiZ;
  aPsiChiZ = aPsiPsiZ * v.transpose();
  bPsiChiZ = bPsiPsiZ * u.hermitianConjugate();

  /// checked and corrected 4/10/12
  ComplexMatrix aPsiChiGam(2, 2), bPsiChiGam(2, 2);
  aPsiChiGam = e * v.transpose(); 
  bPsiChiGam = e * u.hermitianConjugate();

  for (int i=1; i<=2; i++) 
    for (int j=1; j<=2; j++) 
      for (int k=1; k<=dimN; k++) {
	//W
	sigmaL(i, j) = sigmaL(i, j) +
	  (aChi0PsicW(k, i).conj() * aChi0PsicW(k, j) * b1pNeut(k,1)).real();
	sigmaR(i, j) = sigmaR(i, j) +
	  (bChi0PsicW(k, i).conj() * bChi0PsicW(k, j) * b1pNeut(k,1)).real();
	sigmaS(i, j) = sigmaS(i, j) - 4.0 * mneut(k) * 
	  (bChi0PsicW(k, i).conj() * aChi0PsicW(k, j) * b0pNeut(k,1)).real();
	
		if (k <= 2) {
	  /// Z0
	  sigmaL(i, j) = sigmaL(i, j) +
	    (aPsiChiZ(i, k).conj() * aPsiChiZ(j, k) * b1pCha(k,2)).real();
	  sigmaR(i, j) = sigmaR(i, j) +
	    (bPsiChiZ(i, k).conj() * bPsiChiZ(j, k) * b1pCha(k,2)).real();
	  sigmaS(i, j) = sigmaS(i, j) - 4.0 * mch(k) * 
	    (bPsiChiZ(i, k).conj() * aPsiChiZ(j, k) * b0pCha(k,2)).real();
	  //photon
	  sigmaL(i, j) = sigmaL(i, j) +
	    (aPsiChiGam(i, k).conj() * aPsiChiGam(j, k) * b1pCha(k,1)).real();
	  sigmaR(i, j) = sigmaR(i, j) +
	    (bPsiChiGam(i, k).conj() * bPsiChiGam(j, k) * b1pCha(k,1)).real();
	  sigmaS(i, j) = sigmaS(i, j) - 4.0 * mch(k) *
	    (bPsiChiGam(i, k).conj() * aPsiChiGam(j, k) * b0pCha(k,1)).real();

		}

      }
}

//PA: adds gauge boson contribitions to the left right and scalar parts 
//of the chargino self energy
template<class SoftPars>
void Softsusy<SoftPars>::addChaLoopHiggs(double /* p */, DoubleMatrix & sigmaL, DoubleMatrix & sigmaR, DoubleMatrix & sigmaS, DoubleMatrix b1pCha, DoubleMatrix b0pCha, DoubleMatrix b1pNeut, DoubleMatrix b0pNeut) const{
  double g = displayGaugeCoupling(2), 
    gp = displayGaugeCoupling(1) * sqrt(0.6), 
    tanb = displayTanb();
  double beta = atan(tanb);
  double sinb = sin(beta), cosb = cos(beta);

  ComplexMatrix n(forLoops.nBpmz);
  DoubleVector mneut(forLoops.mnBpmz);
  ComplexMatrix u(forLoops.uBpmz), v(forLoops.vBpmz); 
  DoubleVector mch(forLoops.mchBpmz); 
  const int dimN =  mneut.displayEnd();
  /// checked and corrected
  DoubleMatrix aPsiPsiHc1(dimN, 2), bPsiPsiHc1(dimN, 2);
  DoubleMatrix aPsiPsiHc2(dimN, 2), bPsiPsiHc2(dimN, 2);
  aPsiPsiHc1(1, 2) = gp / root2;
  bPsiPsiHc2(1, 2) = aPsiPsiHc1(1, 2);
  aPsiPsiHc1(2, 2) = g / root2;
  bPsiPsiHc2(2, 2) = g / root2;
  aPsiPsiHc1(3, 1) = -g;
  bPsiPsiHc2(4, 1) = g;
  ComplexMatrix aChiPsiHc1(dimN, 2), bChiPsiHc1(dimN, 2);
  ComplexMatrix aChiPsiHc2(dimN, 2), bChiPsiHc2(dimN, 2);
  aChiPsiHc1 = n.complexConjugate() * bPsiPsiHc1;
  bChiPsiHc1 = n * aPsiPsiHc1;
  aChiPsiHc2 = n.complexConjugate() * bPsiPsiHc2;
  bChiPsiHc2 = n * aPsiPsiHc2;
  ComplexMatrix aChiPsiHHp(dimN, 2), bChiPsiHHp(dimN, 2);
  ComplexMatrix aChiPsiHGp(dimN, 2), bChiPsiHGp(dimN, 2);
  for (int i=1; i<=dimN; i++)
    for (int j=1; j<=2; j++) {
      aChiPsiHGp(i, j) = cosb * aChiPsiHc1(i, j) + sinb * aChiPsiHc2(i, j);
      bChiPsiHGp(i, j) = cosb * bChiPsiHc1(i, j) + sinb * bChiPsiHc2(i, j);
      aChiPsiHHp(i, j) =-sinb * aChiPsiHc1(i, j) + cosb * aChiPsiHc2(i, j);
      bChiPsiHHp(i, j) =-sinb * bChiPsiHc1(i, j) + cosb * bChiPsiHc2(i, j);
    }

  /// checked this block
  ComplexMatrix aPsiPsis1(2, 2), aPsiPsis2(2, 2), 
    aPsiPsip1(2, 2), aPsiPsip2(2, 2);
  ComplexMatrix bPsiPsis1(2, 2), bPsiPsis2(2, 2), 
    bPsiPsip1(2, 2), bPsiPsip2(2, 2);
  aPsiPsis1(1, 2) = g / root2;
  aPsiPsis2(2, 1) = g / root2;
  aPsiPsip1(1, 2) = g / root2;
  aPsiPsip2(2, 1) =-g / root2;
  /// checked and corrected 2/12/08
  bPsiPsis1 = aPsiPsis1.transpose(); 
  bPsiPsis2 = aPsiPsis2.transpose();
  /// end of correction 2/12/08
  bPsiPsip1 = -1.0 * aPsiPsip1.transpose();
  bPsiPsip2 = -1.0 * aPsiPsip2.transpose();
  ComplexMatrix aPsiChis1(2, 2), aPsiChis2(2, 2), 
    aPsiChip1(2, 2), aPsiChip2(2, 2);
  ComplexMatrix bPsiChis1(2, 2), bPsiChis2(2, 2), 
    bPsiChip1(2, 2), bPsiChip2(2, 2);
  aPsiChis1 = aPsiPsis1 * u.hermitianConjugate();
  aPsiChis2 = aPsiPsis2 * u.hermitianConjugate();
  aPsiChip1 = aPsiPsip1 * u.hermitianConjugate();
  aPsiChip2 = aPsiPsip2 * u.hermitianConjugate();
  bPsiChis1 = bPsiPsis1 * v.transpose();
  bPsiChis2 = bPsiPsis2 * v.transpose();
  bPsiChip1 = bPsiPsip1 * v.transpose();
  bPsiChip2 = bPsiPsip2 * v.transpose();
  ComplexMatrix aPsiChiH(2, 2), aPsiChih(2, 2), aPsiChiG(2, 2), aPsiChiA(2, 2);
  ComplexMatrix bPsiChiH(2, 2), bPsiChih(2, 2), bPsiChiG(2, 2), 
    bPsiChiA(2, 2);
  double cosa = cos(forLoops.thetaH), sina = sin(forLoops.thetaH);
  aPsiChiH = cosa * aPsiChis1 + sina * aPsiChis2;
  aPsiChih =-sina * aPsiChis1 + cosa * aPsiChis2;
  aPsiChiG = cosb * aPsiChip1 + sinb * aPsiChip2;
  aPsiChiA =-sinb * aPsiChip1 + cosb * aPsiChip2;
  bPsiChiH = cosa * bPsiChis1 + sina * bPsiChis2;
  bPsiChih =-sina * bPsiChis1 + cosa * bPsiChis2;
  bPsiChiG = cosb * bPsiChip1 + sinb * bPsiChip2;
  bPsiChiA =-sinb * bPsiChip1 + cosb * bPsiChip2;

  /// checked and corrected  
  for (int i=1; i<=2; i++) 
    for (int j=1; j<=2; j++) 
      for (int k=1; k<=dimN; k++) {
	/// G+ 
	sigmaL(i, j) = sigmaL(i, j) + 0.5 *
	  (aChiPsiHGp(k, i).conj() * aChiPsiHGp(k, j) * b1pNeut(k,1)).real();
	sigmaR(i, j) = sigmaR(i, j) + 0.5 *
	  (bChiPsiHGp(k, i).conj() * bChiPsiHGp(k, j) * b1pNeut(k,1)).real();
	sigmaS(i, j) = sigmaS(i, j) + mneut(k) * 
	  (bChiPsiHGp(k, i).conj() * aChiPsiHGp(k, j) * b0pNeut(k,1)).real();

	/// H+	
	sigmaL(i, j) = sigmaL(i, j) + 0.5 *
	  (aChiPsiHHp(k, i).conj() * aChiPsiHHp(k, j) * b1pNeut(k,2)).real();
	sigmaR(i, j) = sigmaR(i, j) + 0.5 *
	  (bChiPsiHHp(k, i).conj() * bChiPsiHHp(k, j) * b1pNeut(k,2)).real();
	sigmaS(i, j) = sigmaS(i, j) + mneut(k) * 
	  (bChiPsiHHp(k, i).conj() * aChiPsiHHp(k, j) * b0pNeut(k,2)).real();

	if (k <= 2) {
	  /// H
	  sigmaL(i, j) = sigmaL(i, j) + 0.5 *
	    (aPsiChiH(i, k).conj() * aPsiChiH(j, k) * b1pCha(k,5)).real();
	  sigmaR(i, j) = sigmaR(i, j) + 0.5 *
	    (bPsiChiH(i, k).conj() * bPsiChiH(j, k) * b1pCha(k,5)).real();
	  sigmaS(i, j) = sigmaS(i, j) + mch(k) * 
	    (bPsiChiH(i, k).conj() * aPsiChiH(j, k) * b0pCha(k,5)).real();
	  
	  /// h
	  sigmaL(i, j) = sigmaL(i, j) + 0.5 *
	    (aPsiChih(i, k).conj() * aPsiChih(j, k) * b1pCha(k,4)).real();
	  sigmaR(i, j) = sigmaR(i, j) + 0.5 *
	    (bPsiChih(i, k).conj() * bPsiChih(j, k) * b1pCha(k,4)).real();
	  sigmaS(i, j) = sigmaS(i, j) + mch(k) * 
	    (bPsiChih(i, k).conj() * aPsiChih(j, k) * b0pCha(k,4)).real();

	  /// G0
	  sigmaL(i, j) = sigmaL(i, j) + 0.5 *
	    (aPsiChiG(i, k).conj() * aPsiChiG(j, k) * b1pCha(k,2)).real();
	  sigmaR(i, j) = sigmaR(i, j) + 0.5 *
	    (bPsiChiG(i, k).conj() * bPsiChiG(j, k) * b1pCha(k,2)).real();
	  sigmaS(i, j) = sigmaS(i, j) + mch(k) * 
	    (bPsiChiG(i, k).conj() * aPsiChiG(j, k) * b0pCha(k,2)).real();

	  /// A0
	  sigmaL(i, j) = sigmaL(i, j) + 0.5 *
	    (aPsiChiA(i, k).conj() * aPsiChiA(j, k) * b1pCha(k,3)).real();
	  sigmaR(i, j) = sigmaR(i, j) + 0.5 *
	    (bPsiChiA(i, k).conj() * bPsiChiA(j, k) * b1pCha(k,3)).real();
	  sigmaS(i, j) = sigmaS(i, j) + mch(k) * 
	    (bPsiChiA(i, k).conj() * aPsiChiA(j, k) * b0pCha(k,3)).real();
	}
      }
}

template<class SoftPars>
void Softsusy<SoftPars>::addCharginoLoop(double p, DoubleMatrix & mass) {
  ///  double p = sqrt(forLoops.mchBpmz(1) * forLoops.mchBpmz(2)); 
  double mz = displayMzRun(), mw = displayMwRun(), q = displayMu();
  ComplexMatrix n(forLoops.nBpmz);
  DoubleVector mneut(forLoops.mnBpmz);
  ComplexMatrix u(forLoops.uBpmz), v(forLoops.vBpmz); 
  DoubleVector mch(forLoops.mchBpmz); 
  const int dimN =  mneut.displayEnd();
  DoubleMatrix sigmaL(2, 2), sigmaR(2, 2), sigmaS(2, 2);
  //PA: calculate P-V's fpr gauge and Higgs at outset 
  //avoids repeated calls for same function
  //1st indice runs over neutralino or charginos
  //2nd runs over higgs and gaige bosons in bases:
  //for Cha - (gamma,Z,A,h,H)
  //for Neut - (W, Hpm) 
  DoubleMatrix b1pCha(2,5), b0pCha(2,5), b1pNeut(dimN,2), b0pNeut(dimN,2); 
  for (int k=1; k<=2; k++) {
    b1pCha(k,1) = b1(p, mch(k), 0., q);	
    b0pCha(k,1) = b0(p, mch(k), 0., q);	
    b1pCha(k,2) = b1(p, mch(k), mz, q);
    b0pCha(k,2) = b0(p, mch(k), mz, q);
    b1pCha(k,3) = b1(p, mch(k), forLoops.mA0(1), q);
    b0pCha(k,3) = b0(p, mch(k), forLoops.mA0(1), q);
    b1pCha(k,4) = b1(p, mch(k), forLoops.mh0(1), q);
    b0pCha(k,4) = b0(p, mch(k), forLoops.mh0(1), q);
    b1pCha(k,5) = b1(p, mch(k), forLoops.mh0(2), q);
    b0pCha(k,5) = b0(p, mch(k), forLoops.mh0(2), q);
  }
  
  for (int k=1; k<=dimN; k++) {
    b1pNeut(k,1) = b1(p, mneut(k), mw, q);
    b0pNeut(k,1) = b0(p, mneut(k), mw, q);
    b1pNeut(k,2) = b1(p, mneut(k), forLoops.mHpm, q);
    b0pNeut(k,2) = b0(p, mneut(k), forLoops.mHpm, q);
  }
  
  //PA: sfermion contributions to sigmaL, sigmaR and sigmaS   
  addChaLoopSfermion(p, sigmaL, sigmaR, sigmaS);
  //Gauge bosons
  addChaLoopGauge(p, sigmaL, sigmaR, sigmaS, b1pCha, b0pCha, b1pNeut, b0pNeut);
  //Higgs
  addChaLoopHiggs(p, sigmaL, sigmaR, sigmaS, b1pCha, b0pCha, b1pNeut, b0pNeut);
 
  mass = mass - 1.0 / (16.0 * sqr(PI)) * 
    (sigmaR * mass + mass * sigmaL + sigmaS);
}

/// checked
template<class SoftPars>
void Softsusy<SoftPars>::charginos(int accuracy, double piwwtMS) {
  double tanb = displayTanb(), smu = displaySusyMu();
  DoubleMatrix mCh(2, 2);
  double m2 = displayGaugino(2); 

  double mwOneLarg = sqr(displayMw()) + piwwtMS;

  /// tree level mass matrix
  mCh(1, 1) = m2;
  mCh(2, 1) = root2 * sqrt(fabs(mwOneLarg)) * cos(atan(tanb)); 
  mCh(1, 2) = mCh(2, 1) * tanb;
  mCh(2, 2) = smu;
  if (accuracy == 0) {
    physpars.mch = mCh.asy2by2(physpars.thetaL, physpars.thetaR);
    return;
  }
  DoubleMatrix mCh2(mCh);

  double p1 = fabs(forLoops.mch(1)), p2 = fabs(forLoops.mch(2));
  addCharginoLoop(p1, mCh);
  addCharginoLoop(p2, mCh2);
  
  double x = 0., y = 0.;
  DoubleVector mch1(mCh.asy2by2(physpars.thetaL, physpars.thetaR));
  DoubleVector mch2(mCh2.asy2by2(x, y));
  physpars.mch(1) = mch1(1);
  /// You should take the sign of the chargino mass to be the same as
  /// got from the chargino_1 determination. Otherwise, if there's a
  /// difference, this will screw things up...
  double sgn_mass = mch1(2) / abs(mch1(2));
  physpars.mch(2) = sgn_mass * abs(mch2(2));
}

template<class SoftPars>
double Softsusy<SoftPars>::thet(double a, double b, double c) {
  double yy = maximum(sqr(a), sqr(b));
  yy = maximum(yy, sqr(c));
  return log(yy / sqr(displayMu()));
}

template<class SoftPars>
void Softsusy<SoftPars>::addNeutLoopSfermion(double p, DoubleMatrix & sigmaL, DoubleMatrix & sigmaR, DoubleMatrix & sigmaS) {
  double    g        = displayGaugeCoupling(2);
  double    gp       = displayGaugeCoupling(1) * sqrt(0.6);
  double    ht       = displayDrBarPars().ht;
  double    htau     = displayDrBarPars().htau;
  double    hb       = displayDrBarPars().hb;
  double    q        = displayMu();
  double    beta     = atan(displayTanb());
  double    mt       = ht * displayHvev() * sin(beta) / root2;
  double    mb       = hb * displayHvev() * cos(beta) / root2;
  double    mtau     = htau * displayHvev() * cos(beta) / root2;
  double    thetat   = forLoops.thetat;
  double    thetab   = forLoops.thetab;
  double    thetatau = forLoops.thetatau;

  DoubleVector msup(2), msdown(2), msel(2), mscharm(2), msstrange(2), 
    msmuon(2), msnumu(2), mstop(2), msbot(2), mstau(2), msnutau(2), msnue(2);
  msup(1) = forLoops.mu(1, 1);      msup(2) = forLoops.mu(2, 1); 
  mscharm(1) = forLoops.mu(1, 2);   mscharm(2) = forLoops.mu(2, 2); 
  mstop(1) = forLoops.mu(1, 3);     mstop(2) = forLoops.mu(2, 3); 
  msdown(1) = forLoops.md(1, 1);    msdown(2) = forLoops.md(2, 1); 
  msstrange(1) = forLoops.md(1, 2); msstrange(2) = forLoops.md(2, 2); 
  msbot(1) = forLoops.md(1, 3);     msbot(2) = forLoops.md(2, 3); 
  msel(1) = forLoops.me(1, 1);      msel(2) = forLoops.me(2, 1); 
  msmuon(1) = forLoops.me(1, 2);    msmuon(2) = forLoops.me(2, 2); 
  mstau(1) = forLoops.me(1, 3);     mstau(2) = forLoops.me(2, 3); 
  msnue(1) = forLoops.msnu(1); 
  msnumu(1) = forLoops.msnu(2);
  msnutau(1) = forLoops.msnu(3);

  ComplexMatrix n(forLoops.nBpmz);
  DoubleVector mneut(forLoops.mnBpmz);
  ComplexMatrix u(forLoops.uBpmz), v(forLoops.vBpmz); 
  DoubleVector mch(forLoops.mchBpmz); 

  const int rank = mneut.displayEnd();

  DoubleMatrix aPsiUu(rank, 2), aPsiDd(rank, 2), 
    aPsiEe(rank, 2), aPsiNuNu(rank, 2);
  DoubleMatrix bPsiUu(rank, 2), bPsiDd(rank, 2), 
    bPsiEe(rank, 2), bPsiNuNu(rank, 2);
  DoubleMatrix aPsiTt(rank, 2), aPsiBb(rank, 2), 
    aPsiTauTau(rank, 2), aPsiNutNut(rank, 2);
  DoubleMatrix bPsiTt(rank, 2), bPsiBb(rank, 2), 
    bPsiTauTau(rank, 2), bPsiNutNut(rank, 2);

  aPsiUu(1, 2) = gp / root2 * yuR;
  bPsiUu(1, 1) = gp / root2 * yuL;
  bPsiUu(2, 1) = g * root2 * 0.5;
  aPsiTt(1, 2) = gp / root2 * yuR;
  bPsiTt(1, 1) = gp / root2 * yuL;
  bPsiTt(2, 1) = g * root2 * 0.5;
  aPsiTt(4, 1) = ht;
  bPsiTt(4, 2) = ht;

  aPsiDd(1, 2) = gp / root2 * ydR;
  bPsiDd(1, 1) = gp / root2 * ydL;
  bPsiDd(2, 1) = -g * root2 * 0.5;
  aPsiBb(1, 2) = gp / root2 * ydR;
  bPsiBb(1, 1) = gp / root2 * ydL;
  bPsiBb(2, 1) = -g * root2 * 0.5;
  aPsiBb(3, 1) = hb;
  bPsiBb(3, 2) = hb;

  bPsiNuNu(1, 1) = gp / root2 * ynuL;
  bPsiNuNu(2, 1) = g * root2 * 0.5;
  bPsiNutNut(1, 1) = gp / root2 * ynuL;
  bPsiNutNut(2, 1) = g * root2 * 0.5;

  aPsiEe(1, 2) = gp / root2 * yeR;
  bPsiEe(1, 1) = gp / root2 * yeL;
  bPsiEe(2, 1) = -g * root2 * 0.5;
  aPsiTauTau(1, 2) = gp / root2 * yeR;
  bPsiTauTau(1, 1) = gp / root2 * yeL;
  bPsiTauTau(2, 1) = -g * root2 * 0.5;
  aPsiTauTau(3, 1) = htau;
  bPsiTauTau(3, 2) = htau;

  /// mix up the third family sfermions
  DoubleMatrix O(2, 2);
  O = rot2d(thetat);
  DoubleVector t1(2), t2(2), tt(2);
  int i; for (i=1; i<=rank; i++) {
    tt(1) = aPsiTt(i, 1); tt(2) = aPsiTt(i, 2);      
    t1 = O * tt;
    tt(1) = bPsiTt(i, 1); tt(2) = bPsiTt(i, 2);      
    t2 = O * tt;
    aPsiTt(i, 1) = t1(1);     aPsiTt(i, 2) = t1(2); 
    bPsiTt(i, 1) = t2(1);     bPsiTt(i, 2) = t2(2); 
  }

  O = rot2d(thetab);
  for (i=1; i<=rank; i++) {
    tt(1) = aPsiBb(i, 1); tt(2) = aPsiBb(i, 2);      
    t1 = O * tt;
    tt(1) = bPsiBb(i, 1); tt(2) = bPsiBb(i, 2);      
    t2 = O * tt;
    aPsiBb(i, 1) = t1(1);     aPsiBb(i, 2) = t1(2); 
    bPsiBb(i, 1) = t2(1);     bPsiBb(i, 2) = t2(2); 
  }

  O = rot2d(thetatau);
  for (i=1; i<=rank; i++) {
    tt(1) = aPsiTauTau(i, 1); tt(2) = aPsiTauTau(i, 2);      
    t1 = O * tt;
    tt(1) = bPsiTauTau(i, 1); tt(2) = bPsiTauTau(i, 2);      
    t2 = O * tt;
    aPsiTauTau(i, 1) = t1(1);     aPsiTauTau(i, 2) = t1(2); 
    bPsiTauTau(i, 1) = t2(1);     bPsiTauTau(i, 2) = t2(2); 
  }

  int j, k;
 for (i=1; i<=rank; i++) 
    for (j=1; j<=rank; j++) 
      for (k=1; k<=2; k++) {
	sigmaL(i, j) = sigmaL(i, j) + 
	  (3.0 * aPsiUu(i, k) * aPsiUu(j, k) * b1(p, 0., msup(k), q) +
	   3.0 * aPsiDd(i, k) * aPsiDd(j, k) * b1(p, 0., msdown(k), q) +
	   aPsiEe(i, k) * aPsiEe(j, k) * b1(p, 0., msel(k), q) +
	   aPsiNuNu(i, k) * aPsiNuNu(j, k) * b1(p, 0., msnue(k), q));
	sigmaL(i, j) = sigmaL(i, j) + 
	  (3.0 * aPsiUu(i, k) * aPsiUu(j, k) * b1(p, 0., mscharm(k), q) +
	   3.0 * aPsiDd(i, k) * aPsiDd(j, k) * b1(p, 0., msstrange(k), q) +
	   aPsiEe(i, k) * aPsiEe(j, k) * b1(p, 0., msmuon(k), q) +
	   aPsiNuNu(i, k) * aPsiNuNu(j, k) * b1(p, 0., msnumu(k), q));
	sigmaL(i, j) = sigmaL(i, j) + 
	  (3.0 * aPsiTt(i, k) * aPsiTt(j, k) * b1(p, mt, mstop(k), q) +
	   3.0 * aPsiBb(i, k) * aPsiBb(j, k) * b1(p, mb, msbot(k), q) +
	   aPsiTauTau(i, k) * aPsiTauTau(j, k) * b1(p, mtau, mstau(k), q) +
	   aPsiNutNut(i, k) * aPsiNutNut(j, k) * b1(p, 0., msnutau(k), q));
	sigmaR(i, j) = sigmaR(i, j) + 
	  (3.0 * bPsiUu(i, k) * bPsiUu(j, k) * b1(p, 0., msup(k), q) +
	   3.0 * bPsiDd(i, k) * bPsiDd(j, k) * b1(p, 0., msdown(k), q) +
	   bPsiEe(i, k) * bPsiEe(j, k) * b1(p, 0., msel(k), q) +
	   bPsiNuNu(i, k) * bPsiNuNu(j, k) * b1(p, 0., msnue(k), q));
	sigmaR(i, j) = sigmaR(i, j) + 
	  (3.0 * bPsiUu(i, k) * bPsiUu(j, k) * b1(p, 0., mscharm(k), q) +
	   3.0 * bPsiDd(i, k) * bPsiDd(j, k) * b1(p, 0., msstrange(k), q) +
	   bPsiEe(i, k) * bPsiEe(j, k) * b1(p, 0., msmuon(k), q) +
	   bPsiNuNu(i, k) * bPsiNuNu(j, k) * b1(p, 0., msnumu(k), q));
	sigmaR(i, j) = sigmaR(i, j) + 
	  (3.0 * bPsiTt(i, k) * bPsiTt(j, k) * b1(p, mt, mstop(k), q) +
	   3.0 * bPsiBb(i, k) * bPsiBb(j, k) * b1(p, mb, msbot(k), q) +
	   bPsiTauTau(i, k) * bPsiTauTau(j, k) * b1(p, mtau, mstau(k), q) +
	   bPsiNutNut(i, k) * bPsiNutNut(j, k) * b1(p, 0., msnutau(k), q));
	sigmaS(i, j) = sigmaS(i, j) + 2.0 * 
	  (3.0 * bPsiTt(i, k) * aPsiTt(j, k) * mt * b0(p, mt, mstop(k), q) +
	   3.0 * bPsiBb(i, k) * aPsiBb(j, k) * mb * b0(p, mb, msbot(k), q) +
	   bPsiTauTau(i, k) * aPsiTauTau(j, k) * mtau *
	   b0(p, mtau, mstau(k), q));
      }
}

template<class SoftPars>
void Softsusy<SoftPars>::getNeutPassarinoVeltman(double p, double q, DoubleMatrix & b0fn, DoubleMatrix & b1fn) {
  drBarPars forLoops(displayDrBarPars());
  double    mz = displayMzRun(), mw = displayMwRun();
  DoubleVector mneut(forLoops.mnBpmz);
  DoubleVector mch(forLoops.mchBpmz); 

  const int rank = mneut.displayEnd();

  int k;
     for (k=1; k<=rank; k++) {
       if (k<=2){
	 /// LCT: W+
	 b0fn(k, 1) = b0(p, mch(k), mw, q);
	 b1fn(k, 1) = b1(p, mch(k), mw, q);
	 /// LCT: H+
	 b0fn(k, 2) = b0(p, mch(k), forLoops.mHpm, q); 
	 b1fn(k, 2) = b1(p, mch(k), forLoops.mHpm, q);
       }
 	
       /// LCT: Z0
       b0fn(k, 3) = b0(p, mneut(k), mz, q);
       b1fn(k, 3) = b1(p, mneut(k), mz, q);

       /// LCT: A0
       b0fn(k, 4) = b0(p, mneut(k), forLoops.mA0(1), q);
       b1fn(k, 4) = b1(p, mneut(k), forLoops.mA0(1), q);

       /// LCT: h0
       b0fn(k, 5) = b0(p, mneut(k), forLoops.mh0(1), q);
       b1fn(k, 5) = b1(p, mneut(k), forLoops.mh0(1), q);
       b0fn(k, 6) = b0(p, mneut(k), forLoops.mh0(2), q);
       b1fn(k, 6) = b1(p, mneut(k), forLoops.mh0(2), q);
     }
}

template<class SoftPars>
void Softsusy<SoftPars>::addNeutLoopGauge(double p, DoubleMatrix & sigmaL, 
DoubleMatrix & sigmaR, DoubleMatrix & sigmaS) {
  double g = displayGaugeCoupling(2); 
  double q = displayMu();

  ComplexMatrix n(forLoops.nBpmz);
  DoubleVector mneut(forLoops.mnBpmz);
  ComplexMatrix u(forLoops.uBpmz), v(forLoops.vBpmz); 
  DoubleVector mch(forLoops.mchBpmz); 

  const int rank = mneut.displayEnd();

 /// checked
  ComplexMatrix aPsi0PsicW(rank, 2), bPsi0PsicW(rank, 2), aPsi0ChicW(rank, 2),
    bPsi0ChicW(rank, 2);
  aPsi0PsicW(2, 1) = - g;
  bPsi0PsicW(2, 1) = - g;
  aPsi0PsicW(4, 2) = g / root2;		     
  bPsi0PsicW(3, 2) = -g / root2;		     
  aPsi0ChicW = aPsi0PsicW * v.transpose();
  bPsi0ChicW = bPsi0PsicW * u.hermitianConjugate();

  /// checked
  ComplexMatrix aPsiPsiZ(rank, rank), bPsiPsiZ(rank, rank), 
    aPsiChiZ(rank, rank), bPsiChiZ(rank, rank);
  double sinthW = calcSinthdrbar();
  double thetaWDRbar = asin(sinthW);
  double costh = cos(thetaWDRbar);
  aPsiPsiZ(3, 3) = g / (2.0 * costh);
  aPsiPsiZ(4, 4) =-g / (2.0 * costh);
  bPsiPsiZ = -1.0 * aPsiPsiZ;
  aPsiChiZ = aPsiPsiZ * n.transpose();
  bPsiChiZ = bPsiPsiZ * n.hermitianConjugate();

  DoubleMatrix b0fn(rank, 6), b1fn(rank, 6);
  getNeutPassarinoVeltman(p, q, b0fn, b1fn);

  int i,j,k;
  for (i=1; i<=rank; i++) 
    for (j=1; j<=rank; j++) 
      for (k=1; k<=rank; k++) {
	if (k<=2) {
	  sigmaL(i, j) = sigmaL(i, j) + 2.0 * 
	    (aPsi0ChicW(i, k).conj() * aPsi0ChicW(j, k) * b1fn(k, 1)).real();
	  sigmaR(i, j) = sigmaR(i, j) + 2.0 * 
	    (bPsi0ChicW(i, k).conj() * bPsi0ChicW(j, k) * b1fn(k, 1)).real();
	  sigmaS(i, j) = sigmaS(i, j) - 8.0 * mch(k) * 
	    (bPsi0ChicW(i, k).conj() * aPsi0ChicW(j, k) * b0fn(k, 1)).real();
	}
	
	sigmaL(i, j) = sigmaL(i, j) +
	  (aPsiChiZ(i, k).conj() * aPsiChiZ(j, k) * b1fn(k, 3)).real();
	sigmaR(i, j) = sigmaR(i, j) +
	  (bPsiChiZ(i, k).conj() * bPsiChiZ(j, k) * b1fn(k, 3)).real();
	sigmaS(i, j) = sigmaS(i, j) - 4.0 * mneut(k) * 
	  (bPsiChiZ(i, k).conj() * aPsiChiZ(j, k) * b0fn(k, 3)).real();
      }
}

template<class SoftPars>
void Softsusy<SoftPars>::addNeutLoopHiggs(double p, DoubleMatrix & sigmaL, 
DoubleMatrix & sigmaR, DoubleMatrix & sigmaS) {
  double g    = displayGaugeCoupling(2);
  double gp   = displayGaugeCoupling(1) * sqrt(0.6);
  double q    = displayMu();
  double beta = atan(displayTanb());
  double sinb = sin(beta), cosb = cos(beta);

  ComplexMatrix n(forLoops.nBpmz);
  DoubleVector mneut(forLoops.mnBpmz);
  ComplexMatrix u(forLoops.uBpmz), v(forLoops.vBpmz); 
  DoubleVector mch(forLoops.mchBpmz);

 /// checked
  DoubleMatrix aPsiPsiHc1(4, 2), bPsiPsiHc1(4, 2);
  DoubleMatrix aPsiPsiHc2(4, 2), bPsiPsiHc2(4, 2);
  aPsiPsiHc1(1, 2) = gp / root2;
  bPsiPsiHc2(1, 2) = aPsiPsiHc1(1, 2);
  aPsiPsiHc1(2, 2) = g / root2;
  bPsiPsiHc2(2, 2) = g / root2;
  aPsiPsiHc1(3, 1) = -g;
  bPsiPsiHc2(4, 1) = g;
  ComplexMatrix aPsiChiHc1(4, 2), bPsiChiHc1(4, 2);
  ComplexMatrix aPsiChiHc2(4, 2), bPsiChiHc2(4, 2);
  aPsiChiHc1 = aPsiPsiHc1 * u.hermitianConjugate();
  aPsiChiHc2 = aPsiPsiHc2 * u.hermitianConjugate();
  bPsiChiHc1 = bPsiPsiHc1 * v.transpose();
  bPsiChiHc2 = bPsiPsiHc2 * v.transpose();
  ComplexMatrix aPsiChiHHp(4, 2), bPsiChiHHp(4, 2);
  ComplexMatrix aPsiChiHGp(4, 2), bPsiChiHGp(4, 2);
  int i,j,k; for (i=1; i<=4; i++)
    for (j=1; j<=2; j++) {
      aPsiChiHGp(i, j) = cosb * aPsiChiHc1(i, j) + sinb * aPsiChiHc2(i, j);
      bPsiChiHGp(i, j) = cosb * bPsiChiHc1(i, j) + sinb * bPsiChiHc2(i, j);
      aPsiChiHHp(i, j) =-sinb * aPsiChiHc1(i, j) + cosb * aPsiChiHc2(i, j);
      bPsiChiHHp(i, j) =-sinb * bPsiChiHc1(i, j) + cosb * bPsiChiHc2(i, j);
    }
  
  /// checked this block
  ComplexMatrix aPsiPsis1(4, 4), aPsiPsis2(4, 4), 
    aPsiPsip1(4, 4), aPsiPsip2(4, 4);
  ComplexMatrix bPsiPsis1(4, 4), bPsiPsis2(4, 4), 
    bPsiPsip1(4, 4), bPsiPsip2(4, 4);
  ComplexMatrix aPsiChis1(4, 4), aPsiChis2(4, 4), 
    aPsiChip1(4, 4), aPsiChip2(4, 4);
  ComplexMatrix bPsiChis1(4, 4), bPsiChis2(4, 4), 
    bPsiChip1(4, 4), bPsiChip2(4, 4);
  aPsiPsis1(1, 3) = -gp * 0.5;
  aPsiPsis1(2, 3) = g * 0.5;
  aPsiPsis2(2, 4) = -g * 0.5;
  aPsiPsis2(1, 4) = gp * 0.5;
  aPsiPsip1(1, 3) = -gp * 0.5;
  aPsiPsip1(2, 3) = g * 0.5;
  aPsiPsip2(2, 4) = g * 0.5;
  aPsiPsip2(1, 4) = -gp * 0.5;
  aPsiPsis1.symmetrise();
  aPsiPsis2.symmetrise();
  aPsiPsip1.symmetrise();
  aPsiPsip2.symmetrise();
  bPsiPsis1 = aPsiPsis1;
  bPsiPsis2 = aPsiPsis2;
  bPsiPsip1 =-1.0 * aPsiPsip1;
  bPsiPsip2 =-1.0 * aPsiPsip2;
  aPsiChis1 = aPsiPsis1 * n.hermitianConjugate();
  aPsiChis2 = aPsiPsis2 * n.hermitianConjugate();
  aPsiChip1 = aPsiPsip1 * n.hermitianConjugate();
  aPsiChip2 = aPsiPsip2 * n.hermitianConjugate();
  bPsiChis1 = bPsiPsis1 * n.transpose();
  bPsiChis2 = bPsiPsis2 * n.transpose();
  bPsiChip1 = bPsiPsip1 * n.transpose();
  bPsiChip2 = bPsiPsip2 * n.transpose();
  ComplexMatrix aPsiChiH(4, 4), aPsiChih(4, 4), aPsiChiG(4, 4), aPsiChiA(4, 4);
  ComplexMatrix bPsiChiH(4, 4), bPsiChih(4, 4), bPsiChiG(4, 4), 
    bPsiChiA(4, 4);
  double cosa = cos(forLoops.thetaH), sina = sin(forLoops.thetaH);
  aPsiChiH = cosa * aPsiChis1 + sina * aPsiChis2;
  aPsiChih =-sina * aPsiChis1 + cosa * aPsiChis2;
  aPsiChiG = cosb * aPsiChip1 + sinb * aPsiChip2;
  aPsiChiA =-sinb * aPsiChip1 + cosb * aPsiChip2;
  bPsiChiH = cosa * bPsiChis1 + sina * bPsiChis2;
  bPsiChih =-sina * bPsiChis1 + cosa * bPsiChis2;
  bPsiChiG = cosb * bPsiChip1 + sinb * bPsiChip2;
  bPsiChiA =-sinb * bPsiChip1 + cosb * bPsiChip2;

  DoubleMatrix b0fn(4, 6), b1fn(4, 6);
  getNeutPassarinoVeltman(p, q, b0fn, b1fn);

  for (i=1; i<=4; i++) 
    for (j=1; j<=4; j++) 
      for (k=1; k<=4; k++) {
	if (k<=2) {
	  /// G+
	  sigmaL(i, j) = sigmaL(i, j) + 
	    (aPsiChiHGp(i, k).conj() * aPsiChiHGp(j, k) * b1fn(k, 1)).real();
	  sigmaR(i, j) = sigmaR(i, j) + 
	    (bPsiChiHGp(i, k).conj() * bPsiChiHGp(j, k) * b1fn(k, 1)).real();
	  sigmaS(i, j) = sigmaS(i, j) + 2.0 * mch(k) * 
	    (bPsiChiHGp(i, k).conj() * aPsiChiHGp(j, k) * b0fn(k, 1)).real();

	  /// H+
	  sigmaL(i, j) = sigmaL(i, j) + 
	    (aPsiChiHHp(i, k).conj() * aPsiChiHHp(j, k) * b1fn(k, 2)).real();
	  sigmaR(i, j) = sigmaR(i, j) + 
	    (bPsiChiHHp(i, k).conj() * bPsiChiHHp(j, k) * b1fn(k, 2)).real();
	  sigmaS(i, j) = sigmaS(i, j) + 2.0 * mch(k) * 
	    (bPsiChiHHp(i, k).conj() * aPsiChiHHp(j, k) * b0fn(k, 2)).real();
	}
	
	/// H
	sigmaL(i, j) = sigmaL(i, j) + 0.5 *
	  (aPsiChiH(i, k).conj() * aPsiChiH(j, k) * b1fn(k, 6)).real();
	sigmaR(i, j) = sigmaR(i, j) + 0.5 *  
	  (bPsiChiH(i, k).conj() * bPsiChiH(j, k) * b1fn(k, 6)).real();
	sigmaS(i, j) = sigmaS(i, j) + mneut(k) * 
	  (bPsiChiH(i, k).conj() * aPsiChiH(j, k) * b0fn(k, 6)).real();

	/// h
	sigmaL(i, j) = sigmaL(i, j) + 0.5 *
	  (aPsiChih(i, k).conj() * aPsiChih(j, k) * b1fn(k, 5)).real();
	sigmaR(i, j) = sigmaR(i, j) + 0.5 *
	  (bPsiChih(i, k).conj() * bPsiChih(j, k) * b1fn(k, 5)).real();
	sigmaS(i, j) = sigmaS(i, j) + mneut(k) * 
	  (bPsiChih(i, k).conj() * aPsiChih(j, k) * b0fn(k, 5)).real();

	/// G0
	sigmaL(i, j) = sigmaL(i, j) + 0.5 *
	  (aPsiChiG(i, k).conj() * aPsiChiG(j, k) * b1fn(k, 3)).real();
	sigmaR(i, j) = sigmaR(i, j) + 0.5 * 
	  (bPsiChiG(i, k).conj() * bPsiChiG(j, k) * b1fn(k, 3)).real();
	sigmaS(i, j) = sigmaS(i, j) + mneut(k) * 
	  (bPsiChiG(i, k).conj() * aPsiChiG(j, k) * b0fn(k, 3)).real();

	/// A0
	sigmaL(i, j) = sigmaL(i, j) + 0.5 *
	  (aPsiChiA(i, k).conj() * aPsiChiA(j, k) * b1fn(k, 4)).real();
	sigmaR(i, j) = sigmaR(i, j) + 0.5 *
	  (bPsiChiA(i, k).conj() * bPsiChiA(j, k) * b1fn(k, 4)).real();
	sigmaS(i, j) = sigmaS(i, j) + mneut(k) * 
	  (bPsiChiA(i, k).conj() * aPsiChiA(j, k) * b0fn(k, 4)).real();
	}
}

template<class SoftPars>
void Softsusy<SoftPars>::addNeutralinoLoop(double p, DoubleMatrix & mass) {

  DoubleMatrix sigmaL(4, 4), sigmaR(4, 4), sigmaS(4, 4);

  /// LCT: Sfermion contribution
  addNeutLoopSfermion(p, sigmaL, sigmaR, sigmaS);
  /// LCT: Gauge contributions
  addNeutLoopGauge(p, sigmaL, sigmaR, sigmaS);
  /// LCT: Higgs contribution
  addNeutLoopHiggs(p, sigmaL, sigmaR, sigmaS);

  DoubleMatrix deltaM(4, 4);
  deltaM = -sigmaR * mass - mass * sigmaL - sigmaS;

  deltaM = (deltaM + deltaM.transpose()) / (32.0 * sqr(PI));

  mass = mass + deltaM;
}
  


/// mixNeut set to diagonal = mixNeut^T mNeutralino mixNeut: checked
template<class SoftPars>
void Softsusy<SoftPars>::neutralinos(int accuracy, double /* piwwtMS */,
                                     double /* pizztMS */) {
  double tanb = displayTanb();
  double cosb = cos(atan(tanb));
  DoubleMatrix mNeut(4, 4);
  
  double m1 = displayGaugino(1);
  double m2 = displayGaugino(2); 
  double smu = displaySusyMu();

  /// tree level
  mNeut(1, 1) = m1;
  mNeut(2, 2) = m2;
  mNeut(1, 3) = - displayMzRun() * cosb * calcSinthdrbar();
  mNeut(1, 4) = - mNeut(1, 3) * tanb;
  mNeut(2, 3) = displayMwRun() * cosb;

  mNeut(2, 4) = - mNeut(2, 3) * tanb;
  mNeut(3, 4) = - smu;

  /// symmetrise tree-level
  mNeut.symmetrise();

  if (accuracy == 0) {
    mNeut.diagonaliseSym(physpars.mixNeut, physpars.mneut);
    return;
  }

  DoubleMatrix mNeut4(mNeut), mNeut2(mNeut), mNeut3(mNeut);

  addNeutralinoLoop(fabs(forLoops.mneut(1)), mNeut);
  addNeutralinoLoop(fabs(forLoops.mneut(2)), mNeut2);
  addNeutralinoLoop(fabs(forLoops.mneut(3)), mNeut3);
  addNeutralinoLoop(fabs(forLoops.mneut(4)), mNeut4);
  
  DoubleVector mneut(4), mneut2(4), mneut3(4), mneut4(4);

  DoubleMatrix dummyMix(4, 4);
  double acceptableTol = TOLERANCE * 1.0e-3;
  
  if (mNeut.diagonaliseSym(physpars.mixNeut, mneut) > acceptableTol ||
      mNeut2.diagonaliseSym(dummyMix, mneut2) > acceptableTol ||
      mNeut3.diagonaliseSym(dummyMix, mneut3) > acceptableTol ||
      mNeut4.diagonaliseSym(dummyMix, mneut4) > acceptableTol) { 
    std::ostringstream ii;
    ii << "accuracy bad in neutralino diagonalisation"<< std::flush;
    ii << "diagonalising " << physpars.mneut << " with "   
       << physpars.mixNeut;
    throw ii.str(); 
  }

  /// We should choose sign conventions from the case where the mixing is
  /// defined, in case there is a difference 
  physpars.mneut(1) = mneut(1); 
  physpars.mneut(2) = mneut(2) / abs(mneut(2)) * abs(mneut2(2));
  physpars.mneut(3) = mneut(3) / abs(mneut(3)) * abs(mneut3(3)); 
  physpars.mneut(4) = mneut(4) / abs(mneut(4)) * abs(mneut4(4));
}

/// One loop corrections to gluino pole mass: hep-ph/9606211
/// BUG fixed to use g3 at current scale 8.1.2001
/// Changed to resummed version 11.05.2001
template<class SoftPars>
void Softsusy<SoftPars>::gluino(int accuracy) {

  if (accuracy == 0) {
    physpars.mGluino = displayGaugino(3);
    return;
  }

  if (forLoops.mGluino == 0.0) forLoops.mGluino = displayGaugino(3);

  double p = fabs(displayGaugino(3)), Q = displayMu();

  /// SUSY QCD radiative correction (gluon/gluino)
  double delta = 15.0 + 9.0 * log(sqr(Q / p));
  
  /// Quark/squark correction
  delta = delta -
    (b1(p, dataSet.displayMass(legacy::mUp), forLoops.mu(1, 1), Q) +
     b1(p, dataSet.displayMass(legacy::mCharm), forLoops.mu(1, 2), Q) +
     b1(p, forLoops.mt, forLoops.mu(1, 3), Q) + 
     b1(p, dataSet.displayMass(legacy::mUp), forLoops.mu(2, 1), Q) + 
     b1(p, dataSet.displayMass(legacy::mCharm), forLoops.mu(2, 2), Q) + 
     b1(p, forLoops.mt, forLoops.mu(2, 3), Q) + 
     b1(p, dataSet.displayMass(legacy::mDown), forLoops.md(1, 1), Q) +
     b1(p, dataSet.displayMass(legacy::mStrange), forLoops.md(1, 2), Q) +
     b1(p, forLoops.mb, forLoops.md(1, 3), Q) +
     b1(p, dataSet.displayMass(legacy::mDown), forLoops.md(2, 1), Q) +
     b1(p, dataSet.displayMass(legacy::mStrange), forLoops.md(2, 2), Q) + 
     b1(p, forLoops.mb, forLoops.md(2, 3), Q) );

  /// Third family mixing contribution: NB changed sign of these 4/6/10
  /// Matchev says BPMZ may be wrong! Fixed again by dividing by M3 on 26/8/11
  /// Corrected this incorrect "fix" 2/10/12
  delta = delta + forLoops.mt * sin(2.0 * forLoops.thetat) / displayGaugino(3) *
    (b0(p, forLoops.mt, forLoops.mu(1, 3), Q) - 
     b0(p, forLoops.mt, forLoops.mu(2, 3), Q));
  delta = delta + forLoops.mb * sin(2.0 * forLoops.thetab) / displayGaugino(3) *
    (b0(p, forLoops.mb, forLoops.md(1, 3), Q) - 
     b0(p, forLoops.mb, forLoops.md(2, 3), Q)); 
  
  delta = delta * sqr(displayGaugeCoupling(3)) / (16.0 * sqr(PI));

  if (testNan(delta)) {
    if (PRINTOUT > 2) cout << "Nan in gluino loop\n";
    flagProblemThrown(true);
    physpars.mGluino = displayGaugino(3);
    return;
  } 

  physpars.mGluino = displayGaugino(3) * (1.0 + delta); 
}

template<class SoftPars>
double Softsusy<SoftPars>::calcRunMtQCD() const {
  double mt = forLoops.mt;
/// 1 loop QCD only -- DRbar 10% correction
  double qcd = - (5.0 + 6.0 * log(displayMu() / mt)) * 4.0 *
    sqr(displayGaugeCoupling(3)) / 3.0;
  
  /// 2 loop QCD: hep-ph/0210258 -- debugged 15-6-03
  double l = 2.0 * log(mt / displayMu());

double  twoLoopQcd = sqr(sqr(displayGaugeCoupling(3))) * 
    (-0.5383144424082562 + 0.1815337873591885 * l - 
     0.03799544386587666 * sqr(l));

  return qcd + twoLoopQcd;
}

template<class SoftPars>
double Softsusy<SoftPars>::calcRunMtStopGluino() const {
  double    mstop1  = forLoops.mu(1,3);
  double    mstop2  = forLoops.mu(2,3);
  double    mg      = forLoops.mGluino; 
  double    thetat  = forLoops.thetat ;
  double    mtpole  = dataSet.displayPoleMt();
  double p = mtpole;
  double q = displayMu();
  double mt = forLoops.mt;
  /// stop/gluino correction 6% correction
  double  stopGluino = 4.0 * sqr(displayGaugeCoupling(3)) / 3.0 *
    (b1(p, mg, mstop1, q) + 
     b1(p, mg, mstop2, q) -
     sin(2.0 * thetat) * mg / mtpole * //PA: should be running mass?
     (b0(p, mg, mstop1, q) - 
      b0(p, mg, mstop2, q)));



/// 2 loop QCD involving MSSM sparticles -- hep-ph/0210258, in the
  /// approximation that all squarks and the gluino 
  /// have mass mSUSY: a few per mille error induced at SPS1a.
  double twoLoopMssm = 0.0;
  if (includeTwoLoopMssmCorrectionsToMt) {
  const static double cf = 4.0 / 3.0, ca = 3.0;
  /// colour weighted average mass scale of squarks and gluino
  double m = (3.0 * (forLoops.mu(1, 3) + forLoops.mu(2, 3) + 
		     forLoops.mu(1, 2) + forLoops.mu(2, 2) + 
		     forLoops.mu(1, 1) + forLoops.mu(2, 1) + 
		     forLoops.md(1, 3) + forLoops.md(2, 3) + 
		     forLoops.md(1, 2) + forLoops.md(2, 2) + 
		     forLoops.md(1, 1) + forLoops.md(2, 1)) +
	      8.0 * mg) / 44.0;
  double aq = displaySoftA(UA, 3, 3) - displaySusyMu() / displayTanb();
  double logMoQsq = 2.0 * log(m / q);
 twoLoopMssm = -cf * sqr(sqr(displayGaugeCoupling(3))) / 
    (16.0 * sqr(PI)) *
    (47.0 / 3.0 + 20.0 * logMoQsq + 12.0 * logMoQsq * log(m / mt) +
     cf * (23.0 / 24.0 - 13.0 / 6.0 * logMoQsq + 0.5 * sqr(logMoQsq) -
	   6.0 * logMoQsq * log(mt / q)) + 
     ca * (175.0 / 72.0 + 41.0 / 6.0 * logMoQsq - 0.5 * sqr(logMoQsq) -
	   4.0 * logMoQsq * log(mt / q)) +
     aq / m * (-4.0 - 8.0 * logMoQsq) + 
     cf * aq / m * (7.0 / 3.0 - 11.0 / 3.0 * logMoQsq + 6.0 * log(mt / q)) +
     ca * aq / m * (-8.0 / 3.0 + 4.0 * logMoQsq));
  }


  return stopGluino + twoLoopMssm;

}

template<class SoftPars>
double Softsusy<SoftPars>::calcRunMtHiggs() const {

  double    mH      = forLoops.mh0(2); 
  double    alpha   = forLoops.thetaH ;
  double    mh0     = forLoops.mh0(1);
  double    mA      = forLoops.mA0(1);
  double    mHc     = forLoops.mHpm;
  double    beta    = atan(displayTanb());
  double    g       = displayGaugeCoupling(2);
  double    e       = g * calcSinthdrbar();
  double    ht      = forLoops.ht;
  double    hb      = forLoops.hb;
  double    mtpole  = dataSet.displayPoleMt();
  double    mt      = forLoops.mt;
  double    mb      = forLoops.mb;
  double    q       = displayMu();
  double    mz      = displayMzRun();
  double    thetaWDRbar = asin(calcSinthdrbar());
  double    cw2DRbar    = sqr(cos(thetaWDRbar));
  double p = mtpole;
  
  const double  costh   = (displayMw() / displayMz());
  const double    cw2   = sqr(costh) ;
  const double    sw2   = (1.0 - cw2);
  double gtL = 0.5 - 2.0 * sw2 / 3.0, gtR = 2.0 * sw2 / 3.0;

  double higgs = sqr(ht) / 2.0 * 
    (sqr(sin(alpha)) * (b1(p, mt, mH, q) + b0(p, mt, mH, q))
     + sqr(cos(alpha)) * (b1(p, mt, mh0, q) + 
			  b0(p, mt, mh0, q))
     + sqr(cos(beta)) * (b1(p, mt, mA, q) - 
			 b0(p, mt, mA, q)) 
     + sqr(sin(beta)) * (b1(p, mt, mz, q) - 
			 b0(p, mt, mz, q))) + 
    0.5 * ( (sqr(hb) * sqr(sin(beta)) + sqr(ht) * sqr(cos(beta))) * 
	    b1(p, mb, mHc, q) +
	    (sqr(g) + sqr(hb) * sqr(cos(beta)) + sqr(ht) * sqr(sin(beta))) * 
	    b1(p, mb, displayMwRun(), q)) +
    sqr(hb) * sqr(cos(beta)) * (b0(p, mb, mHc, q) 
				- b0(p, mb, displayMwRun(), q)) -
    sqr(e) * 4.0 / 9.0 * (5.0 + 6.0 * log(q / mt)) +
    sqr(g) / cw2DRbar * ( (sqr(gtL) + sqr(gtR)) * b1(p, mt, mz, q) +
			  4.0 * gtL * gtR * b0(p, mt, mz, q) );

  return higgs;

}

template<class SoftPars>
double Softsusy<SoftPars>::calcRunMtNeutralinos() const {

  double    q       = displayMu();
  double    mtpole  = dataSet.displayPoleMt();
  double    ht      = forLoops.ht;
  double    g       = displayGaugeCoupling(2);
  double    gp      = displayGaugeCoupling(1) * sqrt(0.6);
  double    thetat  = forLoops.thetat ;
  double p = mtpole;
  DoubleMatrix neutralinoContribution(4, 2);
  
  DoubleVector aPsi0TStopr(4), bPsi0TStopr(4), aPsi0TStopl(4),
    bPsi0TStopl(4); 
  aPsi0TStopr(1) = - 4 * gp / (3.0 * root2);
  bPsi0TStopl(1) = gp / (3.0 * root2);
  bPsi0TStopl(2) = g / root2;
  aPsi0TStopl(4) = ht;
  bPsi0TStopr(4) = ht;

  ComplexVector aChi0TStopl(4), bChi0TStopl(4), aChi0TStopr(4),
    bChi0TStopr(4);

  /// Neutralinos
  ComplexMatrix n(forLoops.nBpmz);
  DoubleVector mneut(forLoops.mnBpmz);
  
  aChi0TStopl = n.complexConjugate() * aPsi0TStopl;
  bChi0TStopl = n * bPsi0TStopl;
  aChi0TStopr = n.complexConjugate() * aPsi0TStopr;
  bChi0TStopr = n * bPsi0TStopr;

  ComplexMatrix aNeutTStop(4, 2), bNeutTStop(4, 2);
  DoubleMatrix fNeutTStop(4, 2), gNeutTStop(4, 2);
  double neutralinos = 0.0;
  int i, j;
  DoubleMatrix O(2, 2);
  ComplexVector tt(2), t1(2), t2(2);
  O = rot2d(thetat);
  for (i=1; i<=4; i++) {
    tt(1) = aChi0TStopl(i); tt(2) = aChi0TStopr(i);      
    t1 = O * tt;

    tt(1) = bChi0TStopl(i); tt(2) = bChi0TStopr(i);      
    t2 = O * tt;    
    for (j=1; j<=2; j++) {
      aNeutTStop(i, j) = t1(j);
      bNeutTStop(i, j) = t2(j);
      /// functions of couplings needed for loops
      fNeutTStop(i, j) = sqr(aNeutTStop(i, j).mod()) + 
	sqr(bNeutTStop(i, j).mod());

      gNeutTStop(i, j) = 2.0 * 
	(aNeutTStop(i, j) * bNeutTStop(i, j).conj()).real(); 
      
      neutralinoContribution(i, j) = (fNeutTStop(i, j) * 
	 b1(p, mneut(i), forLoops.mu(j, 3), q) + 
	 gNeutTStop(i, j) * mneut(i) /  mtpole *  
	 b0(p, mneut(i), forLoops.mu(j, 3), q)) * 0.5;

      neutralinos = neutralinos + neutralinoContribution(i, j);
    }
  }

  return neutralinos;

}

template<class SoftPars>
double Softsusy<SoftPars>::calcRunMtCharginos() const {
  double    mtpole  = dataSet.displayPoleMt();
  double    q       = displayMu();
  double    ht      = forLoops.ht;
  double    hb      = forLoops.hb;
  double    thetab  = forLoops.thetab ;
  double    g       = displayGaugeCoupling(2);
  double p = mtpole;
  
  ComplexMatrix u(forLoops.uBpmz), v(forLoops.vBpmz); 
  DoubleVector mch(forLoops.mchBpmz); 
  DoubleVector bPsicTSbotl(2), bPsicTSbotr(2), aPsicTSbotl(2); 

  bPsicTSbotl(1) = g;
  bPsicTSbotr(2) = -hb;
  aPsicTSbotl(2) = -ht;
  
  DoubleVector aPsicTSbotr(2), aPsicCSbotl(2);
  ComplexVector aChicTSbotr(2), aChicTSbotl(2), bChicTSbotl(2),
      bChicTSbotr(2);
  ComplexMatrix aChTSbot(2, 2), bChTSbot(2, 2);
  DoubleMatrix fChTSbot(2, 2), gChTSbot(2, 2); 

  aChicTSbotl = v.complexConjugate() * aPsicTSbotl;
  bChicTSbotl = u * bPsicTSbotl;
  aChicTSbotr = v.complexConjugate() * aPsicTSbotr;
  bChicTSbotr = u * bPsicTSbotr;
       
  double charginoContribution = 0.0;
  DoubleMatrix O(2, 2);
  ComplexVector tt(2), t1(2), t2(2);
  for(int i=1; i<=2; i++) {  
    O = rot2d(thetab);      
    tt(1) = aChicTSbotl(i); tt(2) = aChicTSbotr(i);
    t1 = O * tt;
    tt(1) = bChicTSbotl(i); tt(2) = bChicTSbotr(i);      
    t2 = O * tt;
    for (int j=1; j<=2; j++) {
	aChTSbot(i, j) = t1(j);
	bChTSbot(i, j) = t2(j);

      	fChTSbot(i, j) = sqr(aChTSbot(i, j).mod()) + 
	  sqr(bChTSbot(i, j).mod());
	gChTSbot(i, j) = 2.0 * (aChTSbot(i, j) * 
				bChTSbot(i, j).conj()).real(); 
      
      	charginoContribution = charginoContribution + 
	  (fChTSbot(i, j) * 
	   b1(p, mch(i), fabs(forLoops.md(j, 3)),
	      q) +
	   gChTSbot(i, j) * mch(i) / mtpole * 
	   b0(p, mch(i), forLoops.md(j, 3), q)) * 0.5;
    }            
  }

  return charginoContribution;

}
///  Formulae from hep-ph/9801365: checked but should be checked again!
/// Implicitly calculates at the current scale.
template<class SoftPars>
double Softsusy<SoftPars>::calcRunningMt() {
  double    mtpole  = dataSet.displayPoleMt();
  double resigmat = 0.0; 
  double qcd = 0.0, stopGluino = 0.0, higgs = 0.0; 
  //one and two loop qcd
  qcd =  calcRunMtQCD();
  resigmat = resigmat + qcd;
  /// stop/gluino correction 6% correction
  stopGluino = calcRunMtStopGluino();
  resigmat = resigmat + stopGluino;
  /// rest are extra bits from Matchev et al: 2% corrections  
  //Higgs contribution
  higgs = calcRunMtHiggs();
  resigmat = resigmat + higgs;
  /// Neutralino contribution
  double neutralinos = calcRunMtNeutralinos();
  resigmat = resigmat + neutralinos;
  // Chargino contribution
  double charginoContribution = calcRunMtCharginos();
  resigmat = resigmat + charginoContribution; 
    
  resigmat = resigmat * mtpole / (16.0 * sqr(PI));  

  return mtpole + resigmat;
}

template<class SoftPars>
double Softsusy<SoftPars>::calcRunMbDrBarConv() const {
   double    g       = displayGaugeCoupling(2);
   double    gp      = displayGaugeCoupling(1) * sqrt(0.6);
   double alphasMZ = sqr(displayGaugeCoupling(3)) / (4.0 * PI);
   double conversion = (1.0 - alphasMZ / (3.0 * PI) 
                        + 3.0 * sqr(g) / (128.0 * sqr(PI)) +
                        13.0 * sqr(gp) / (1152. * sqr(PI)));
   return conversion;
}

template<class SoftPars>
double Softsusy<SoftPars>::calcRunMbSquarkGluino() const {
  double    msbot1  = displayDrBarPars().md(1,3);
  double    msbot2  = displayDrBarPars().md(2,3);
  double    mg      = displayDrBarPars().mGluino;
  double    thetab  = displayDrBarPars().thetab;
  double mbMZ = dataSet.displayMass(legacy::mBottom),  
     alphasMZ = sqr(displayGaugeCoupling(3)) / (4.0 * PI);
  double p = mbMZ;
  double mbMSSM  = displayDrBarPars().mb;
  double deltaSquarkGluino = - alphasMZ / (3.0 * PI) *
     (b1(p, mg, msbot1, displayMu()) + 
      b1(p, mg, msbot2, displayMu()) - 
      sin(2.0 * thetab) * mg / mbMSSM *  
      (b0(p, mg, msbot1, displayMu()) - 
       b0(p, mg, msbot2, displayMu())));
  
  return deltaSquarkGluino;
}


template<class SoftPars>
double Softsusy<SoftPars>::calcRunMbChargino() const {
  double mbMZ = dataSet.displayMass(legacy::mBottom);
  double p = mbMZ;
  double q = displayMu();
  double   mbMSSM  = displayDrBarPars().mb;
  double g  = displayGaugeCoupling(2);
  double thetat = displayDrBarPars().thetat;
  DoubleVector bPsicBstopl(2), bPsicBstopr(2), 
    aPsicBstopl(2), aPsicBstopr(2); 

  aPsicBstopl(1) = g;
  aPsicBstopr(2) = -displayDrBarPars().ht;
  bPsicBstopl(2) = -displayDrBarPars().hb;
  
  DoubleVector aPsicCStopl(2);
  ComplexVector aChicBstopr(2), aChicBstopl(2), bChicBstopl(2),
      bChicBstopr(2);
  ComplexMatrix aChBstop(2, 2), bChBstop(2, 2);
  DoubleMatrix fChBstop(2, 2), gChBstop(2, 2); 

  ComplexMatrix u(displayDrBarPars().uBpmz), v(displayDrBarPars().vBpmz); 
  DoubleVector mch(displayDrBarPars().mchBpmz); 

  aChicBstopl = v.complexConjugate() * aPsicBstopl;
  bChicBstopl = u * bPsicBstopl;
  aChicBstopr = v.complexConjugate() * aPsicBstopr;
  bChicBstopr = u * bPsicBstopr;
       
  double charginoContribution = 0.0;

  int i, j; DoubleMatrix O(2, 2); ComplexVector tt(2), t1(2), t2(2);
  for(i=1; i<=2; i++) {  
    O = rot2d(thetat);      
    tt(1) = aChicBstopl(i); tt(2) = aChicBstopr(i);
    t1 = O * tt;
    tt(1) = bChicBstopl(i); tt(2) = bChicBstopr(i);      
    t2 = O * tt;
    for (j=1; j<=2; j++) {
	aChBstop(i, j) = t1(j);
	bChBstop(i, j) = t2(j);

      	fChBstop(i, j) = sqr(aChBstop(i, j).mod()) + 
	  sqr(bChBstop(i, j).mod());
	gChBstop(i, j) = 2.0 * (aChBstop(i, j) * 
				bChBstop(i, j).conj()).real(); 
      
      	charginoContribution = charginoContribution + 
	  (fChBstop(i, j) * 
	   b1(p, mch(i), fabs(displayDrBarPars().mu(j, 3)),
	      q) +
	   gChBstop(i, j) * mch(i) / mbMSSM * 
	   b0(p, mch(i), displayDrBarPars().mu(j, 3), q)) * 0.5;
    }            
  }
  double deltaSquarkChargino = -charginoContribution / (16.0 * sqr(PI));
  
  return deltaSquarkChargino;

}

template<class SoftPars>
double Softsusy<SoftPars>::calcRunMbHiggs() const {
  double mbMZ = dataSet.displayMass(legacy::mBottom);
  double p = mbMZ;
  double q = displayMu();
  double deltaHiggs = 0.;
  double hb = displayDrBarPars().hb;
  double ht = displayDrBarPars().ht;
  double mb = displayDrBarPars().mb;
  double mt = displayDrBarPars().mt;
  double  mh      = displayDrBarPars().mh0(1);
  double  mA      = displayDrBarPars().mA0(1);
  double  mH      = displayDrBarPars().mh0(2);
  double  mHp     = displayDrBarPars().mHpm;
  double  ca      = cos(displayDrBarPars().thetaH);
  double  sa      = sin(displayDrBarPars().thetaH);
  double  cosb    = cos(atan(displayTanb()));
  double  sinb    = sin(atan(displayTanb()));
  double  mz = displayMzRun();
  double  mw = displayMwRun();
  double  thetaWDRbar = asin(calcSinthdrbar());
  double  cw2DRbar    = sqr(cos(thetaWDRbar));
  double g  = displayGaugeCoupling(2);
  deltaHiggs = 0.5 * sqr(hb) * 
    (sqr(ca) * (b1(p, mb, mH, q) + b0(p, mb, mH, q)) + 
     sqr(sa) * (b1(p, mb, mh, q) + b0(p, mb, mh, q)) + 
     sqr(sinb) * (b1(p, mb, mA, q) - b0(p, mb, mA, q)) + 
     sqr(cosb) * (b1(p, mb, mz, q) - b0(p, mb, mz, q))) +
    0.5 * ((sqr(ht) * sqr(cosb) + sqr(hb) * sqr(sinb)) * b1(p, mt, mHp, q) +
	   (sqr(g) + sqr(ht) * sqr(sinb) + sqr(hb) * sqr(cosb)) * 
	   b1(p, mt, mw, q)) +
    sqr(ht) * sqr(sinb) * (b0(p, mt, mHp, q) - b0(p, mt, mw, q)) +
    sqr(g) / cw2DRbar * ((sqr(gdL) + sqr(gdR)) * b1(p, mb, mz, q) +
				    4.0 * gdL * gdR * b0(p, mb, mz, q));
  
  deltaHiggs = - deltaHiggs / (16.0 * sqr(PI));

  return deltaHiggs;
}

template<class SoftPars>
double Softsusy<SoftPars>::calcRunMbNeutralinos() const {
  double mbMZ = dataSet.displayMass(legacy::mBottom);
  double p = mbMZ;
  double q = displayMu();
  double thetab  = displayDrBarPars().thetab;
  double g       = displayGaugeCoupling(2);
  double gp      = displayGaugeCoupling(1) * sqrt(0.6);
  double hb      = displayDrBarPars().hb;
  double mbMSSM  = displayDrBarPars().mb;
  ComplexMatrix n(displayDrBarPars().nBpmz);
  DoubleVector mneut(displayDrBarPars().mnBpmz);
  DoubleVector aPsi0Bsbotr(4), bPsi0Bsbotr(4), aPsi0Bsbotl(4),
    bPsi0Bsbotl(4); 
  aPsi0Bsbotr(1) = gp / (root2 * 3.0) * 2.0;
  bPsi0Bsbotl(1) = gp / (root2 * 3.0);
  bPsi0Bsbotl(2) = -root2 * g * 0.5;
  aPsi0Bsbotl(3) = hb;
  bPsi0Bsbotr(3) = hb;

  ComplexVector aChi0Bsbotl(4), bChi0Bsbotl(4), aChi0Bsbotr(4),
    bChi0Bsbotr(4);

  double deltaNeutralino = 0.;
  aChi0Bsbotl = n.complexConjugate() * aPsi0Bsbotl;
  bChi0Bsbotl = n * bPsi0Bsbotl;
  aChi0Bsbotr = n.complexConjugate() * aPsi0Bsbotr;
  bChi0Bsbotr = n * bPsi0Bsbotr;

  ComplexMatrix aNeutBsbot(4, 2), bNeutBsbot(4, 2);
  DoubleMatrix fNeutBsbot(4, 2), gNeutBsbot(4, 2), 
    neutralinoContribution(4, 2);
  int i, j; DoubleMatrix O(2, 2); ComplexVector tt(2), t1(2), t2(2);
  O = rot2d(thetab);
  for (i=1; i<=4; i++) {
    tt(1) = aChi0Bsbotl(i); tt(2) = aChi0Bsbotr(i);      
    t1 = O * tt;

    tt(1) = bChi0Bsbotl(i); tt(2) = bChi0Bsbotr(i);      
    t2 = O * tt;    
    for (j=1; j<=2; j++) {
      aNeutBsbot(i, j) = t1(j);
      bNeutBsbot(i, j) = t2(j);
      /// functions of couplings needed for loops
      fNeutBsbot(i, j) = sqr(aNeutBsbot(i, j).mod()) + 
	sqr(bNeutBsbot(i, j).mod());

      gNeutBsbot(i, j) = 2.0 * 
	(aNeutBsbot(i, j) * bNeutBsbot(i, j).conj()).real(); 
      
      neutralinoContribution(i, j) = (fNeutBsbot(i, j) * 
	 b1(p, mneut(i), displayDrBarPars().md(j, 3), q) + 
	 gNeutBsbot(i, j) * mneut(i) /  mbMSSM *  
	 b0(p, mneut(i), displayDrBarPars().md(j, 3), q)) * 0.5;

      deltaNeutralino = deltaNeutralino + neutralinoContribution(i, j);
    }
  }

  deltaNeutralino = -deltaNeutralino / (16.0 * sqr(PI));

return deltaNeutralino; 
}


template<class SoftPars>
double Softsusy<SoftPars>::calcRunningMb() const {

  if (displayMu() != displayMz()) {
    std::ostringstream ii;
    ii << "Softsusy<SoftPars>::calcRunningMb called with mu=" <<
      displayMu() << endl; 
    throw ii.str();
  }
  
  double mbMZ = dataSet.displayMass(legacy::mBottom);
  /// First convert mbMZ into DRbar value from hep-ph/9703293,0207126,9701308
  /// (SM gauge boson contributions)
  mbMZ = mbMZ * calcRunMbDrBarConv(); 

  double deltaSquarkGluino = calcRunMbSquarkGluino();
  //Chargino-squark loops
  double deltaSquarkChargino = calcRunMbChargino();
  /// Higgs
  double deltaHiggs = calcRunMbHiggs();
  /// Neutralinos
  double deltaNeutralino = calcRunMbNeutralinos();

  /// it's NOT clear if this resummation is reliable in the full 1-loop scheme
  /// but it's at least valid to 1 loop. Warning though: if you add higher
  /// loops, you'll have to re-arrange.
  return mbMZ / (1.0 + deltaSquarkGluino + deltaSquarkChargino + deltaHiggs
		 + deltaNeutralino);
}
template<class SoftPars>
double Softsusy<SoftPars>::calcRunMtauDrBarConv() const {
 double gp = displayGaugeCoupling(1) * sqrt(0.6);
 double conv = (1.0 - 3.0 * (sqr(gp) - sqr(displayGaugeCoupling(2))) / (128.0 * sqr(PI)));
 return conv;
}

template<class SoftPars>
double Softsusy<SoftPars>::calcRunMtauCharginos(double mTauSMMZ) const {
   
  double g        = displayGaugeCoupling(2);
  double htau     = displayDrBarPars().htau;
  double msnutau  = fabs(displayDrBarPars().msnu(3));
  double mTauPole = MTAU;
  double p = mTauPole;
  double q = displayMu();
  DoubleVector aPsicTauSnul(2), bPsicTauSnul(2); 

  aPsicTauSnul(1) = g;
  bPsicTauSnul(2) = -htau;

  ComplexVector aChicTauSnul(2), bChicTauSnul(2);
  DoubleVector fChiTauSnu(2), gChiTauSnu(2); 

  ComplexMatrix u(displayDrBarPars().uBpmz), v(displayDrBarPars().vBpmz); 
  DoubleVector mch(displayDrBarPars().mchBpmz); 

  /// Mass eignebasis of charginos
  aChicTauSnul = v.complexConjugate() * aPsicTauSnul;
  bChicTauSnul = u * bPsicTauSnul;
       
  DoubleVector charg(2);
  int i; 
  for(i=1; i<=2; i++) {  
    fChiTauSnu(i) = sqr(aChicTauSnul(i).mod()) + 
      sqr(bChicTauSnul(i).mod());
    gChiTauSnu(i) = 2.0 * (aChicTauSnul(i) * bChicTauSnul(i).conj()).real(); 
      
    charg(i) = 
      (fChiTauSnu(i) * b1(p, mch(i), msnutau, q) +
       gChiTauSnu(i) * mch(i) / mTauSMMZ * b0(p, mch(i), msnutau, q)) * 0.5;
    }            
  
  double sigmaChargino = (charg(1) + charg(2)) / (16.0 * sqr(PI));
  /// checked charginos

  return sigmaChargino; 

}

template<class SoftPars>
double Softsusy<SoftPars>::calcRunMtauHiggs() const {
  double mTauPole = MTAU;
  double p = mTauPole;
  double q = displayMu();
  double mh  = displayDrBarPars().mh0(1);
  double mA  = displayDrBarPars().mA0(1);
  double mH  = displayDrBarPars().mh0(2);
  double mHp = displayDrBarPars().mHpm;
  double mz  = displayMzRun();
  double mw  = displayMwRun();
  double ca  = cos(displayDrBarPars().thetaH);
  double sa  = sin(displayDrBarPars().thetaH);
  double cosb = cos(atan(displayTanb()));
  double sinb = sin(atan(displayTanb()));
  double htau = displayDrBarPars().htau;
  double mtau = displayDrBarPars().mtau;
  double g    = displayGaugeCoupling(2);
  double thetaWDRbar = asin(calcSinthdrbar());
  double cw2DRbar    = sqr(cos(thetaWDRbar));
  double mnu = 0.;
  double sigmaHiggs = 0.5 * sqr(htau) * 
    (sqr(ca) * (b1(p, mtau, mH, q) + b0(p, mtau, mH, q)) + 
     sqr(sa) * (b1(p, mtau, mh, q) + b0(p, mtau, mh, q)) + 
     sqr(sinb) * (b1(p, mtau, mA, q) - b0(p, mtau, mA, q)) + 
     sqr(cosb) * (b1(p, mtau, mz, q) - b0(p, mtau, mz, q))) +
    0.5 * (sqr(htau) * sqr(sinb) * b1(p, mnu, mHp, q) + 
	   (sqr(g) + sqr(htau) * sqr(cosb)) * b1(p, mnu, mw, q)) +
    sqr(g) / cw2DRbar * ((sqr(geL) + sqr(geR)) * b1(p, mtau, mz, q) + 
    4.0 * geL * geR * b0(p, mtau, mz, q));
  
  sigmaHiggs = sigmaHiggs / (16.0 * sqr(PI));
  
  return sigmaHiggs;
}

template<class SoftPars>
double Softsusy<SoftPars>::calcRunMtauNeutralinos(double mTauSMMZ) const {
  double thetatau = displayDrBarPars().thetatau;
  double g        = displayGaugeCoupling(2);
  double gp       = displayGaugeCoupling(1) * sqrt(0.6);
  double htau = displayDrBarPars().htau;
  double mTauPole = MTAU;
  double p = mTauPole;
  double q = displayMu();
  
  ComplexMatrix n(displayDrBarPars().nBpmz);
  DoubleVector mneut(displayDrBarPars().mnBpmz);
  DoubleVector aPsi0TauStaur(4), bPsi0TauStaur(4), aPsi0TauStaul(4),
    bPsi0TauStaul(4); 
  aPsi0TauStaur(1) = gp / root2 * 2.0;
  bPsi0TauStaul(1) = -gp / root2;
  bPsi0TauStaul(2) = -root2 * g * 0.5;
  aPsi0TauStaul(3) = htau;
  bPsi0TauStaur(3) = htau;

  ComplexVector aChi0TauStaul(4), bChi0TauStaul(4), aChi0TauStaur(4),
    bChi0TauStaur(4);

  double sigmaNeutralino = 0.;
  aChi0TauStaul = n.complexConjugate() * aPsi0TauStaul;
  bChi0TauStaul = n * bPsi0TauStaul;
  aChi0TauStaur = n.complexConjugate() * aPsi0TauStaur;
  bChi0TauStaur = n * bPsi0TauStaur;

  ComplexMatrix aNeutTauStau(4, 2), bNeutTauStau(4, 2);
  DoubleMatrix fNeutTauStau(4, 2), gNeutTauStau(4, 2), 
    neutralinoContribution(4, 2);

  DoubleMatrix O(2, 2); O = rot2d(thetatau);
  ComplexVector t1(2), t2(2), tt(2);
  for (int i=1; i<=4; i++) {
    tt(1) = aChi0TauStaul(i); tt(2) = aChi0TauStaur(i);      
    t1 = O * tt;

    tt(1) = bChi0TauStaul(i); tt(2) = bChi0TauStaur(i);      
    t2 = O * tt;    
    int j;
    for (j=1; j<=2; j++) {
      aNeutTauStau(i, j) = t1(j);
      bNeutTauStau(i, j) = t2(j);
      /// functions of couplings needed for loops
      fNeutTauStau(i, j) = sqr(aNeutTauStau(i, j).mod()) + 
	sqr(bNeutTauStau(i, j).mod());

      gNeutTauStau(i, j) = 2.0 * 
	(aNeutTauStau(i, j) * bNeutTauStau(i, j).conj()).real(); 
      
      neutralinoContribution(i, j) = (fNeutTauStau(i, j) * 
	 b1(p, mneut(i), displayDrBarPars().me(j, 3), q) + 
	 gNeutTauStau(i, j) * mneut(i) /  mTauSMMZ *  
	 b0(p, mneut(i), displayDrBarPars().me(j, 3), q)) * 0.5;

      sigmaNeutralino = sigmaNeutralino + neutralinoContribution(i, j);
    }
  }

  sigmaNeutralino = sigmaNeutralino / (16.0 * sqr(PI));
  return sigmaNeutralino;
}
/// Full BPMZ expression
template<class SoftPars>
double Softsusy<SoftPars>::calcRunningMtau() const {
  /// MSbar value
  double mTauSMMZ = displayDataSet().displayMass(legacy::mTau);
  // double mTauPole = MTAU;
  /// conversion to DRbar
  mTauSMMZ = mTauSMMZ * calcRunMtauDrBarConv();
  /// Chargino contribution  
  double sigmaChargino = calcRunMtauCharginos(mTauSMMZ);
  /// Higgs
  double sigmaHiggs = calcRunMtauHiggs();
  /// Neutralinos
  double  sigmaNeutralino = calcRunMtauNeutralinos(mTauSMMZ);

  /// old calculation of tau mass
  /**  double delta = sqr(displayGaugeCoupling(2)) / (16 * sqr(PI)) *
    (-displaySusyMu()) * displayGaugino(2) * displayTanb() /
    (sqr(displaySusyMu()) - sqr(displayGaugino(2))) *
    (b0(mTauPole, displayGaugino(2), 
	displayDrBarPars().msnu(3), displayMu()) -
	b0(mTauPole, -displaySusyMu(), displayDrBarPars().msnu(3), displayMu()));*/

  /// From hep-ph/9912516
  return mTauSMMZ * (1.0 + sigmaNeutralino + sigmaChargino + sigmaHiggs);
}

template<class SoftPars>
void Softsusy<SoftPars>::treeUpSquark(DoubleMatrix & mass, double mtrun, 
                                      double /* pizztMS */, double sinthDRbarMS,
                                      int family) { 
  const double cu = 2.0 / 3.0;
  double mz2 = sqr(displayMzRun()), mt2 = sqr(mtrun);
  double beta = atan(displayTanb()), mu = displaySusyMu(),
    c2b = cos(2 * beta), tanb = displayTanb(), sinth2 = sqr(sinthDRbarMS);

  mass(1, 1) = displaySoftMassSquared(mQl, family, family) +
    mz2 * (0.5 - cu * sinth2) * c2b;
  mass(2, 2) = displaySoftMassSquared(mUr, family, family) +
    mz2 * cu * sinth2 * c2b;

  if (family != 3) mass(1, 2) = 0.0;
  else {
    if (fabs(displayDrBarPars().ht) < EPSTOL) 
      mass(1, 2) = mtrun * (displaySoftA(UA, 3, 3) - mu / tanb);
    else 
      mass(1, 2) = mtrun * (displayDrBarPars().ut / displayDrBarPars().ht - mu / tanb);

    mass(1, 1) = mass(1, 1) + mt2;
    mass(2, 2) = mass(2, 2) + mt2;
  }

  mass(2, 1) = mass(1, 2);
}

/// Now these are calculated at the squark scale
template<class SoftPars>
void Softsusy<SoftPars>::addSquarkCorrection(DoubleMatrix & mass) {

/// No point adding radiative corrections to tachyonic particles
  if (mass(1, 1) < 0.0 || mass(2, 2) < 0.0) { 
    flagTachyon(sup); flagTachyon(sdown); flagTachyon(scharm); 
    flagTachyon(sstrange);
    return;
  }

  DoubleVector x(2), delta(2);
  int i; for (i=1; i<=2; i++)  
    { 
      x(i)= sqr(forLoops.mGluino) / mass(i, i);

      if (close(x(i), 1.0, EPSTOL))
	delta(i) = sqr(displayGaugeCoupling(3)) / (6.0 * sqr(PI)) * 
	  (1.0 + 3.0 * x(i) 
	   - sqr(x(i)) * log(fabs(x(i))) + 2.0 * x(i) * 
	   log(sqr(displayMu()) / mass(i, i)));
      else
	delta(i) = sqr(displayGaugeCoupling(3)) / (6.0 * sqr(PI)) * 
	  (1.0 + 3.0 * x(i) + sqr(x(i) - 1.0) * log(fabs(x(i) - 1.0))
	   - sqr(x(i)) * log(fabs(x(i))) + 2.0 * x(i) * 
	   log(sqr(displayMu()) / mass(i, i)));

      mass(i, i) = mass(i, i) * (1.0 + delta(i));
    }  
}

template<class SoftPars>
void Softsusy<SoftPars>::addSnuTauSfermion(double /* p */, double & /* stop */, double & sbottom) {
  double  thetatau  = forLoops.thetatau;
  double  ctau      = cos(thetatau);
  double  stau      = sin(thetatau);;
  double  q         = displayMu();
  double  htau      = forLoops.htau, htausq = sqr(htau);
  DoubleVector mstau(2);
  mstau(1)          = forLoops.me(1, 3);
  mstau(2)          = forLoops.me(2, 3);

  sbottom = 
    htausq * (sqr(stau) * a0(mstau(1), q) + sqr(ctau) * a0(mstau(2), q));
}

template<class SoftPars>
double Softsusy<SoftPars>::addSnuTauHiggs(double p, double & higgs) {
  double    sinthDrbar  = calcSinthdrbar();
  double    costhDrbar  = sqrt(1.0 - sqr(sinthDrbar));
  double    costhDrbar2 = 1.0 - sqr(sinthDrbar);
  double    alpha       = forLoops.thetaH;
  double    g	        = displayGaugeCoupling(2);
  double    beta        = atan(displayTanb());
  double    sinb        = sin(beta), cosb = cos(beta);
  double    htau        = forLoops.htau, htausq = sqr(htau);
  double    thetatau    = forLoops.thetatau;
  double    smu         = -displaySusyMu();
  double    q           = displayMu();
  double    v1          = displayHvev() * cos(beta);
  double    mz          = displayMzRun();
  double    mtau        = htau * displayHvev() / root2 * cos(beta);
  DoubleVector msnu(3);
  msnu(1)               = forLoops.msnu(1);
  msnu(2)               = forLoops.msnu(2);
  msnu(3)               = forLoops.msnu(3);
  DoubleVector mstau(2);
  mstau(1)              = forLoops.me(1, 3);
  mstau(2)              = forLoops.me(2, 3);

  /// Define Higgs vector of masses in 't-Hooft Feynman gauge, and couplings:
  DoubleVector higgsm(4), higgsc(2), dnu(4), dnd(4), cn(4);
  assignHiggsSfermions(higgsm, higgsc, dnu, dnd, cn, beta);

  DoubleVector lsSnuLSnuLR(4);
  /// Order (s1 s2 G A, L R)
  lsSnuLSnuLR(1) = g * mz * gnuL * cosb / costhDrbar;
  lsSnuLSnuLR(2) = -g * mz * gnuL * sinb / costhDrbar;

  DoubleVector lHSnuLSnu12(lsSnuLSnuLR);
  DoubleVector temp(2), temp2(2);
  /// Mix CP-even Higgses up
  temp(1) = lsSnuLSnuLR(1);
  temp(2) = lsSnuLSnuLR(2);
  temp2 = rot2d(alpha) * temp;
  lHSnuLSnu12(1) = temp2(1);
  lHSnuLSnu12(2) = temp2(2);
  
  /// Charged Higgs Feynman rules
  DoubleMatrix lChHsnuLstauLR(2, 2); /// (H+ G+, L R) basis
  lChHsnuLstauLR(1, 1) = (g * displayMwRun() * sin(2.0 * beta) / root2
    - htau * mtau * sinb);
  lChHsnuLstauLR(1, 2) = (smu * htau * cosb - 
			   forLoops.utau * sinb);
  lChHsnuLstauLR(2, 1) = (-g * displayMwRun() * cos(2.0 * beta) 
    + htausq * v1 * cosb) / root2;
  lChHsnuLstauLR(2, 2) = htau * smu * sinb + forLoops.utau * cosb;

  DoubleMatrix lChHsnuLstau12(2, 2);
  temp(1) = lChHsnuLstauLR(1, 1);
  temp(2) = lChHsnuLstauLR(1, 2);
  temp2 = rot2d(thetatau) * temp;
  lChHsnuLstau12(1, 1) = temp2(1);
  lChHsnuLstau12(1, 2) = temp2(2);
  temp(1) = lChHsnuLstauLR(2, 1);
  temp(2) = lChHsnuLstauLR(2, 2);
  temp2 = rot2d(thetatau) * temp;
  lChHsnuLstau12(2, 1) = temp2(1);
  lChHsnuLstau12(2, 2) = temp2(2);

  int i;
 for (i=1; i<=4; i++) {
    higgs = higgs + 
      0.5 * (- sqr(g) * gnuL * 0.5 / sqr(costhDrbar) * cn(i)) 
      * a0(higgsm(i), q);
  }

  for (i=3; i<=4; i++) { 
    double a0p = a0(higgsc(i - 2), q);
    higgs = higgs + 
      (htausq * dnu(i) + sqr(g) * 
       (gnuL * 0.5 / costhDrbar2 - 0.5) * cn(i))* a0p;
  }

  int j; for(i=1; i<=4; i++) {
    double b0p = b0(p, higgsm(i), msnu(3), q);
   higgs = higgs + sqr(lHSnuLSnu12(i)) * b0p;
  }

 for(i=1; i<=2; i++)
    for (j=1; j<=2; j++) {
      double b0p = b0(p, higgsc(i), mstau(j), q); 
      higgs = higgs + sqr(lChHsnuLstau12(i, j)) * b0p;
    }

  return higgs;
}

template<class SoftPars>
double Softsusy<SoftPars>::addSnuTauEweak(double p, double & electroweak) {
  double    sinthDrbar  = calcSinthdrbar();
  double    costhDrbar  = sqrt(1.0 - sqr(sinthDrbar));
  double    costhDrbar2 = 1.0 - sqr(sinthDrbar);
  double    g  	        = displayGaugeCoupling(2);
  double    gp          = displayGaugeCoupling(1) * sqrt(0.6);
  double    thetat      = forLoops.thetat;
  double    thetab      = forLoops.thetab;
  double    thetatau    = forLoops.thetatau;
  double    ct          = cos(thetat);
  double    st          = sin(thetat);
  double    cb          = cos(thetab);
  double    sb          = sin(thetab);
  double    ctau        = cos(thetatau);
  double    stau        = sin(thetatau);
  double    q           = displayMu();
  double    mz          = displayMzRun();
  DoubleVector msnu(3);
  msnu(1)               = forLoops.msnu(1);
  msnu(2)               = forLoops.msnu(2);
  msnu(3)               = forLoops.msnu(3);
  DoubleVector mstau(2);
  mstau(1)              = forLoops.me(1, 3);
  mstau(2)              = forLoops.me(2, 3);

  /// EW bosons
  electroweak = electroweak + 
    4.0 * sqr(g) / costhDrbar2 * sqr(gnuL) * a0(mz, q) + 
    2.0 * sqr(g) * a0(displayMwRun(), q) + 
    sqr(g * gnuL / costhDrbar) * ffn(p, msnu(3), mz, q) +
    sqr(g) * 0.5 * (sqr(ctau) * ffn(p, mstau(1), displayMwRun(), q) + 
		    sqr(stau) * ffn(p, mstau(2), displayMwRun(), q)) +
    sqr(g) * 0.25 * 
    (a0(msnu(3), q) + 2.0 *
     (sqr(ctau) * a0(mstau(1), q) + sqr(stau) * a0(mstau(2), q))) +
    sqr(g) * 0.5 * 
    (1.5 * a0(forLoops.mu(1, 1), q) + 
     1.5 * a0(forLoops.mu(1, 2), q) +
     1.5 * (sqr(ct) * a0(forLoops.mu(1, 3), q) + 
	    sqr(st) * a0(forLoops.mu(2, 3), q)) -
     1.5 * a0(forLoops.md(1, 1), q) -
     1.5 * a0(forLoops.md(1, 2), q) -
     1.5 * (sqr(cb) * a0(forLoops.md(1, 3), q) +
	    sqr(sb) * a0(forLoops.md(2, 3), q)) +
     0.5 * (a0(forLoops.msnu(1), q) + 
	    a0(forLoops.msnu(2), q) +
	    a0(forLoops.msnu(3), q)) -
     0.5 * (a0(forLoops.me(1, 1), q) + a0(forLoops.me(1, 2), q) +
	    sqr(ctau) * a0(forLoops.me(1, 3), q) +
	    sqr(stau) * a0(forLoops.me(2, 3), q))) +
    sqr(gp) * 0.25 * sqr(ynuL) * a0(msnu(3), q) +
    sqr(gp) * 0.25 * ynuL * 
    (3.0 * yuL * (a0(forLoops.mu(1, 1), q) + 
		  a0(forLoops.mu(1, 2), q) + 
		  sqr(ct) * a0(forLoops.mu(1, 3), q) + 
		  sqr(st) * a0(forLoops.mu(2, 3), q)) +
     3.0 * yuR * (a0(forLoops.mu(2, 1), q) + 
		  a0(forLoops.mu(2, 2), q) + 
		  sqr(st) * a0(forLoops.mu(1, 3), q) + 
		  sqr(ct) * a0(forLoops.mu(2, 3), q)) +
     3.0 * ydL * (a0(forLoops.md(1, 1), q) + 
		  a0(forLoops.md(1, 2), q) + 
		  sqr(cb) * a0(forLoops.md(1, 3), q) + 
		  sqr(sb) * a0(forLoops.md(2, 3), q)) +
     3.0 * ydR * (a0(forLoops.md(2, 1), q) + 
		  a0(forLoops.md(2, 2), q) + 
		  sqr(sb) * a0(forLoops.md(1, 3), q) + 
		  sqr(cb) * a0(forLoops.md(2, 3), q)) +
     yeL * (a0(forLoops.me(1, 1), q) + 
	    a0(forLoops.me(1, 2), q) + 
	    sqr(ctau) * a0(forLoops.me(1, 3), q) + 
	    sqr(stau) * a0(forLoops.me(2, 3), q)) +
     yeR * (a0(forLoops.me(2, 1), q) + 
	    a0(forLoops.me(2, 2), q) + 
	    sqr(stau) * a0(forLoops.me(1, 3), q) + 
	    sqr(ctau) * a0(forLoops.me(2, 3), q)) +
     ynuL * (a0(forLoops.msnu(1), q) + 
	     a0(forLoops.msnu(2), q) +
	     a0(forLoops.msnu(3), q)));

  return electroweak;
}

template<class SoftPars>
void Softsusy<SoftPars>::addSnuTauGaugino(double p, double & chargino, double & neutralino) {
  double    g	    = displayGaugeCoupling(2);
  double    gp      = displayGaugeCoupling(1) * sqrt(0.6);
  double    beta    = atan(displayTanb());
  double    htau    = forLoops.htau;
  double    mtau    = htau * displayHvev() / root2 * cos(beta);
  double    q       = displayMu();

  /// Neutralino Feynman rules
  ComplexMatrix n(forLoops.nBpmz);
  DoubleVector mneut(forLoops.mnBpmz);
  ComplexMatrix u(forLoops.uBpmz), v(forLoops.vBpmz); 
  DoubleVector mch(forLoops.mchBpmz); 

  const int rank = mneut.displayEnd();

  DoubleVector aPsi0TSnul(rank), bPsi0TSnul(rank); 
  bPsi0TSnul(1) = gp * ynuL / root2;
  bPsi0TSnul(2) = g / root2;

  ComplexVector aChi0TSnul(rank), bChi0TSnul(rank);

  aChi0TSnul = n.complexConjugate() * aPsi0TSnul;
  bChi0TSnul = n * bPsi0TSnul;

  DoubleVector gChi0NuSnuLL(rank), fChi0NuSnuLL(rank);
  int i; for (i=1; i<=rank; i++) {
    fChi0NuSnuLL(i) = (aChi0TSnul(i) * aChi0TSnul(i).conj() + 
      bChi0TSnul(i) * bChi0TSnul(i).conj()).real();
    gChi0NuSnuLL(i) = (bChi0TSnul(i).conj() * aChi0TSnul(i) + 
      bChi0TSnul(i) * aChi0TSnul(i).conj()).real();
  }

  /// Chargino Feynman Rules
  DoubleVector bPsicBSnul(2), aPsicBSnul(2);
  aPsicBSnul(1) = g;
  bPsicBSnul(2) = -htau;
  
  DoubleVector aPsicCStaul(2);
  ComplexVector aChicBSnur(2), aChicBSnul(2), bChicBSnul(2),
      bChicBSnur(2);
  ComplexMatrix aChBSnu(2, 2), bChBSnu(2, 2);

  aChicBSnul = v.complexConjugate() * aPsicBSnul;
  bChicBSnul = u * bPsicBSnul;

  DoubleVector fChBSnuLL(2), gChBSnuLL(2) ;
  for (i=1; i<=2; i++) {
    fChBSnuLL(i) = (aChicBSnul(i).conj() * aChicBSnul(i) +
		      bChicBSnul(i).conj() * bChicBSnul(i)).real();
    gChBSnuLL(i) = (bChicBSnul(i).conj() * aChicBSnul(i) +
		      aChicBSnul(i).conj() * bChicBSnul(i)).real();
  }

  for (i=1; i<=2; i++) {
    double one = gfn(p, mch(i), mtau, q);
    double two = mch(i) * mtau * b0(p, mch(i), mtau, q) * 2.0;
    chargino = chargino + fChBSnuLL(i) * one - gChBSnuLL(i) * two;
  }

  for (i=1; i<=rank; i++) {
    double one = gfn(p, mneut(i), 0., q);
    neutralino = neutralino + fChi0NuSnuLL(i) * one;
  }
}

template<class SoftPars>
void Softsusy<SoftPars>::addSnuTauCorrection(double & mass) {

  /// No point adding radiative corrections to tachyonic particles
  if (mass < 0.0) { 
    flagTachyon(snutau);
    mass = EPSTOL;
    return;
  }

  double p = sqrt(mass);

  /// one-loop correction matrix
  double piSq; 
	
  double stop = 0., sbottom = 0., higgs = 0., electroweak = 0.,
    chargino = 0., neutralino = 0.;

  /// LCT: Sfermion contribution
  addSnuTauSfermion(p, stop, sbottom);
  /// LCT: Higgs contribution
  addSnuTauHiggs(p, higgs);
  /// LCT: Electroweak contribution
  addSnuTauEweak(p, electroweak);
  /// LCT: Gaugino contributions    
  addSnuTauGaugino(p, chargino, neutralino);

  piSq = 1.0 / (16.0 * sqr(PI)) * 
    (stop + sbottom + higgs + electroweak + chargino + neutralino);

  mass = mass - piSq;	  
}

template<class SoftPars>
double Softsusy<SoftPars>::addSnuHiggs(double p, int family, double & higgs) {
  double    sinthDrbar  = calcSinthdrbar();
  double    costhDrbar  = sqrt(1.0 - sqr(sinthDrbar));
  double    costhDrbar2 = 1.0 - sqr(sinthDrbar);
  double    alpha       = forLoops.thetaH;
  double    g	        = displayGaugeCoupling(2);
  double    beta        = atan(displayTanb());
  double    sinb        = sin(beta), cosb = cos(beta);
  double    q           = displayMu();
  double    mz          = displayMzRun();
  DoubleVector msnu(3);
  msnu(1)               = forLoops.msnu(1);
  msnu(2)               = forLoops.msnu(2);
  msnu(3)               = forLoops.msnu(3);

  /// Define Higgs vector of masses in 't-Hooft Feynman gauge, and couplings:
  DoubleVector higgsm(4), higgsc(2), dnu(4), dnd(4), cn(4);
  assignHiggsSfermions(higgsm, higgsc, dnu, dnd, cn, beta);

  DoubleVector lsSnuLSnuLR(4);
  /// Order (s1 s2 G A, L R)
  lsSnuLSnuLR(1) = g * mz * gnuL * cosb / costhDrbar;
  lsSnuLSnuLR(2) = -g * mz * gnuL * sinb / costhDrbar;

  DoubleVector lHSnuLSnu12(lsSnuLSnuLR);
  DoubleVector temp(2), temp2(2);
  /// Mix CP-even Higgses up
  temp(1) = lsSnuLSnuLR(1);
  temp(2) = lsSnuLSnuLR(2);
  temp2 = rot2d(alpha) * temp;
  lHSnuLSnu12(1) = temp2(1);
  lHSnuLSnu12(2) = temp2(2);
  
  /// Charged Higgs Feynman rules
  DoubleMatrix lChHsnuLstauLR(2, 2); /// (H+ G+, L R) basis
  lChHsnuLstauLR(2, 1) = -g * displayMwRun() * cos(2.0 * beta) / root2;
  lChHsnuLstauLR(1, 1) = (g * displayMwRun() * sin(2.0 * beta) / root2);

  DoubleMatrix lChHsnuLstau12(2, 2);
  temp(1) = lChHsnuLstauLR(1, 1);
  temp(2) = lChHsnuLstauLR(1, 2);
  temp2 = temp;
  lChHsnuLstau12(1, 1) = temp2(1);
  lChHsnuLstau12(1, 2) = temp2(2);
  temp(1) = lChHsnuLstauLR(2, 1);
  temp(2) = lChHsnuLstauLR(2, 2);
  temp2 = temp;
  lChHsnuLstau12(2, 1) = temp2(1);
  lChHsnuLstau12(2, 2) = temp2(2);

  int i;
 for (i=1; i<=4; i++) {
    higgs = higgs + 
      0.5 * (- sqr(g) * gnuL * 0.5 / sqr(costhDrbar) * cn(i)) 
      * a0(higgsm(i), q);
  }

  for (i=3; i<=4; i++) { 
    double a0p = a0(higgsc(i - 2), q);
    higgs = higgs +
      (gnuL * 0.5 / costhDrbar2 - 0.5) * cn(i)* a0p * sqr(g);
  }

  int j; for(i=1; i<=4; i++) {
    double b0p = b0(p, higgsm(i), msnu(family), q);
   higgs = higgs + sqr(lHSnuLSnu12(i)) * b0p;
  }

  DoubleVector msel(2); 
  msel(1) = forLoops.me(1, family);
  msel(2) = forLoops.me(2, family);
  for(i=1; i<=2; i++) 
    for (j=1; j<=2; j++) {
      double b0p = b0(p, higgsc(i), msel(j), q); 
      higgs = higgs + sqr(lChHsnuLstau12(i, j)) * b0p;
  }

  return higgs;
}

template<class SoftPars>
double Softsusy<SoftPars>::addSnuEweak(double p, int family, double & electroweak) {
  double    sinthDrbar  = calcSinthdrbar();
  double    costhDrbar  = sqrt(1.0 - sqr(sinthDrbar));
  double    costhDrbar2 = 1.0 - sqr(sinthDrbar);
  double    g	        = displayGaugeCoupling(2);
  double    gp          = displayGaugeCoupling(1) * sqrt(0.6);
  double    thetat      = forLoops.thetat;
  double    thetab      = forLoops.thetab;
  double    thetatau    = forLoops.thetatau;
  double    ct          = cos(thetat) ;
  double    st          = sin(thetat) ;
  double    cb          = cos(thetab) ;
  double    sb          = sin(thetab) ;
  double    ctau        = cos(thetatau);
  double    stau        = sin(thetatau);
  double    q           = displayMu();
  double    mz          = displayMzRun();
  double    meL         = forLoops.me(1, family);
  DoubleVector msnu(3);
  msnu(1)               = forLoops.msnu(1);
  msnu(2)               = forLoops.msnu(2);
  msnu(3)               = forLoops.msnu(3);

  /// LCT: EW contributions
  electroweak = electroweak + 
    4.0 * sqr(g) / costhDrbar2 * sqr(gnuL) * a0(mz, q) + 
    2.0 * sqr(g) * a0(displayMwRun(), q) + 
    sqr(g * gnuL / costhDrbar) * ffn(p, msnu(family), mz, q) +
    sqr(g) * 0.5 * ffn(p, meL, displayMwRun(), q);

  electroweak = electroweak + 
    sqr(g) * 0.25 * 
    (a0(msnu(family), q) + 2.0 * a0(meL, q)) + sqr(g) * 0.5 * 
    (1.5 * a0(forLoops.mu(1, 1), q) + 
     1.5 * a0(forLoops.mu(1, 2), q) +
     1.5 * (sqr(ct) * a0(forLoops.mu(1, 3), q) + 
	    sqr(st) * a0(forLoops.mu(2, 3), q)) -
     1.5 * a0(forLoops.md(1, 1), q) -
     1.5 * a0(forLoops.md(1, 2), q) -
     1.5 * (sqr(cb) * a0(forLoops.md(1, 3), q) +
	    sqr(sb) * a0(forLoops.md(2, 3), q)) +
     0.5 * (a0(forLoops.msnu(1), q) + 
	    a0(forLoops.msnu(2), q) +
	    a0(forLoops.msnu(3), q)) -
     0.5 * (a0(forLoops.me(1, 1), q) + a0(forLoops.me(1, 2), q) +
	    sqr(ctau) * a0(forLoops.me(1, 3), q) +
	    sqr(stau) * a0(forLoops.me(2, 3), q))) +
    sqr(gp) * 0.25 * sqr(ynuL) * a0(msnu(family), q) +
    sqr(gp) * 0.25 * ynuL * 
    (3.0 * yuL * (a0(forLoops.mu(1, 1), q) + 
		  a0(forLoops.mu(1, 2), q) + 
		  sqr(ct) * a0(forLoops.mu(1, 3), q) + 
		  sqr(st) * a0(forLoops.mu(2, 3), q)) +
     3.0 * yuR * (a0(forLoops.mu(2, 1), q) + 
		  a0(forLoops.mu(2, 2), q) + 
		  sqr(st) * a0(forLoops.mu(1, 3), q) + 
		  sqr(ct) * a0(forLoops.mu(2, 3), q)) +
     3.0 * ydL * (a0(forLoops.md(1, 1), q) + 
		  a0(forLoops.md(1, 2), q) + 
		  sqr(cb) * a0(forLoops.md(1, 3), q) + 
		  sqr(sb) * a0(forLoops.md(2, 3), q)) +
     3.0 * ydR * (a0(forLoops.md(2, 1), q) + 
		  a0(forLoops.md(2, 2), q) + 
		  sqr(sb) * a0(forLoops.md(1, 3), q) + 
		  sqr(cb) * a0(forLoops.md(2, 3), q)) +
     yeL * (a0(forLoops.me(1, 1), q) + 
	    a0(forLoops.me(1, 2), q) + 
	    sqr(ctau) * a0(forLoops.me(1, 3), q) + 
	    sqr(stau) * a0(forLoops.me(2, 3), q)) +
     yeR * (a0(forLoops.me(2, 1), q) + 
	    a0(forLoops.me(2, 2), q) + 
	    sqr(stau) * a0(forLoops.me(1, 3), q) + 
	    sqr(ctau) * a0(forLoops.me(2, 3), q)) +
     ynuL * (a0(forLoops.msnu(1), q) + 
	     a0(forLoops.msnu(2), q) +
	     a0(forLoops.msnu(3), q)));

  return electroweak;
}

template<class SoftPars>
void Softsusy<SoftPars>::addSnuGaugino(double p, int /* family */, double & chargino, double & neutralino) {
  double    g	    = displayGaugeCoupling(2);
  double    gp      = displayGaugeCoupling(1) * sqrt(0.6);
  double    q       = displayMu();

  DoubleVector mneut(forLoops.mnBpmz);
  ComplexMatrix u(forLoops.uBpmz), v(forLoops.vBpmz); 
  DoubleVector mch(forLoops.mchBpmz); 

  const int rank = mneut.displayEnd();

 /// Neutralino Feynman rules
  DoubleVector aPsi0TSnul(rank), bPsi0TSnul(rank); 
  bPsi0TSnul(1) = gp * ynuL / root2;
  bPsi0TSnul(2) = g / root2;

  ComplexVector aChi0TSnul(rank), bChi0TSnul(rank);
  ComplexMatrix n(forLoops.nBpmz);

  aChi0TSnul = n.complexConjugate() * aPsi0TSnul;
  bChi0TSnul = n * bPsi0TSnul;

  DoubleVector gChi0NuSnuLL(rank), fChi0NuSnuLL(rank);
  int i; for (i=1; i<=rank; i++) {
    fChi0NuSnuLL(i) = (aChi0TSnul(i) * aChi0TSnul(i).conj() + 
      bChi0TSnul(i) * bChi0TSnul(i).conj()).real();
    gChi0NuSnuLL(i) = (bChi0TSnul(i).conj() * aChi0TSnul(i) + 
      bChi0TSnul(i) * aChi0TSnul(i).conj()).real();
  }

  /// Chargino Feynman Rules
  DoubleVector bPsicBSnul(2), aPsicBSnul(2);
  aPsicBSnul(1) = g;
  
  DoubleVector aPsicCStaul(2);
  ComplexVector aChicBSnur(2), aChicBSnul(2), bChicBSnul(2),
      bChicBSnur(2);
  ComplexMatrix aChBSnu(2, 2), bChBSnu(2, 2);

  aChicBSnul = v.complexConjugate() * aPsicBSnul;
  bChicBSnul = u * bPsicBSnul;

  DoubleVector fChBSnuLL(2), gChBSnuLL(2) ;
  for (i=1; i<=2; i++) {
    fChBSnuLL(i) = (aChicBSnul(i).conj() * aChicBSnul(i) +
		      bChicBSnul(i).conj() * bChicBSnul(i)).real();
    gChBSnuLL(i) = (bChicBSnul(i).conj() * aChicBSnul(i) +
		      aChicBSnul(i).conj() * bChicBSnul(i)).real();
  }

  for (i=1; i<=2; i++) {
    double one = gfn(p, mch(i), 0., q);
    chargino = chargino + fChBSnuLL(i) * one;
  }

  for (i=1; i<=rank; i++) {
    double one = gfn(p, mneut(i), 0., q);
    neutralino = neutralino + fChi0NuSnuLL(i) * one;
  }

}

/// Found+fixed bug 7/09/06. Thanks to J Kersten.
template<class SoftPars>
void Softsusy<SoftPars>::addSnuCorrection(double & mass, int family) {

  /// No point adding radiative corrections to tachyonic particles
  if (mass < 0.0) { 
    if (family == 1) flagTachyon(snue);
    else if (family == 2) flagTachyon(snumu);
    return;
  }

  double p = sqrt(mass);

  /// one-loop correction matrix
  double piSq; 
	
  double higgs = 0., electroweak = 0., chargino = 0., neutralino = 0.;

  /// LCT: Higg contribution
  addSnuHiggs(p, family, higgs);
  /// LCT: Electroweak contributions
  addSnuEweak(p, family, electroweak);
  /// LCT: Gaugino contributions
  addSnuGaugino(p, family, chargino, neutralino);
 
  piSq = 1.0 / (16.0 * sqr(PI)) * 
    (higgs + electroweak + chargino + neutralino);
  
  mass = mass - piSq;	  
}

template<class SoftPars>
DoubleMatrix Softsusy<SoftPars>::addStopQCD(double p, double mt, DoubleMatrix & strong) {
  double    mg      = forLoops.mGluino;
  double    thetat  = forLoops.thetat;
  double    ct      = cos(thetat);
  double    st      = sin(thetat);
  double    q       = displayMu();
  double    g3sq    = sqr(displayGaugeCoupling(3));
  DoubleVector mstop(2);
  mstop(1)          = forLoops.mu(1, 3);
  mstop(2)          = forLoops.mu(2, 3);
  double    a0t1    = a0(mstop(1), q), a0t2 = a0(mstop(2), q);
  double    ft1     = ffn(p, mstop(1), 0.0, q), ft2 = ffn(p, mstop(2), 0.0, q);
  double    ggt     = gfn(p, mg, mt, q);

  strong(1, 1) = 4.0 * g3sq / 3.0 *
    (2.0 * ggt + sqr(ct) * (ft1 + a0t1) + sqr(st) * (ft2 + a0t2));
  strong(2, 2) = 4.0 * g3sq / 3.0 *
    (2.0 * ggt + sqr(st) * (ft1 + a0t1) + sqr(ct) * (ft2 + a0t2));
  strong(1, 2) = 4.0 * g3sq / 3.0 *
    (4.0 * mg * mt * b0(p, mg, mt, q) + 
     st * ct * (ft1 - a0t1 - ft2 + a0t2));

  return strong;
}

template<class SoftPars>
DoubleMatrix Softsusy<SoftPars>::addStopStop(double /* p */, double /* mt */, DoubleMatrix & stop) {
  double thetat     = forLoops.thetat, ct = cos(thetat), st = sin(thetat);
  double q          = displayMu();
  DoubleVector mstop(2);
  mstop(1)          = forLoops.mu(1, 3);
  mstop(2)          = forLoops.mu(2, 3);
  double a0t1       = a0(mstop(1), q), a0t2 = a0(mstop(2), q);
  double ht         = forLoops.ht, htsq = sqr(ht);

  stop(1, 1) = htsq * (sqr(st) * a0t1 + sqr(ct) * a0t2);
  stop(2, 2) = htsq * (sqr(ct) * a0t1 + sqr(st) * a0t2);
  stop(1, 2) = htsq * ct * st * 3.0 * (a0t1 - a0t2);

  return stop;
}

template<class SoftPars>
DoubleMatrix Softsusy<SoftPars>::addStopSbottom(double /* p */, double /* mt */,
                                                DoubleMatrix & sbottom) {
  double thetab     = forLoops.thetab, cb = cos(thetab), sb = sin(thetab) ;
  double q          = displayMu();
  double ht         = forLoops.ht, htsq = sqr(ht);
  double hb         = forLoops.hb, hbsq = sqr(hb);
  DoubleVector msbot(2);
  msbot(1)          = forLoops.md(1, 3);
  msbot(2)          = forLoops.md(2, 3);

  sbottom(1, 1) = 
    hbsq * (sqr(sb) * a0(msbot(1), q) + sqr(cb) * a0(msbot(2), q));
  sbottom(2, 2) = 
    htsq * (sqr(cb) * a0(msbot(1), q) + sqr(sb) * a0(msbot(2), q));

  return sbottom;
}

template<class SoftPars>
DoubleMatrix Softsusy<SoftPars>::addStopHiggs(double p, double /* mt */,
                                              DoubleMatrix & higgs) {
  double    sinthDrbar  = calcSinthdrbar();
  double    costhDrbar  = sqrt(1.0 - sqr(sinthDrbar));
  double    costhDrbar2 = 1.0 - sqr(sinthDrbar);
  double    alpha       = forLoops.thetaH;
  double    g	        = displayGaugeCoupling(2);
  double    beta        = atan(displayTanb());
  double    sinb        = sin(beta), cosb = cos(beta);
  double    smu         = -displaySusyMu(); 
  double    thetat      = forLoops.thetat;
  double    thetab      = forLoops.thetab;
  double    q           = displayMu(); 
  double    ht          = forLoops.ht, htsq = sqr(ht);
  double    hb          = forLoops.hb, hbsq = sqr(hb);
  double    v1          = displayHvev() * cos(beta);
  double    v2          = displayHvev() * sin(beta);
  double    mz          = displayMzRun();
  DoubleVector msbot(2);
  msbot(1)              = forLoops.md(1, 3);
  msbot(2)              = forLoops.md(2, 3);
  DoubleVector mstop(2);
  mstop(1)              = forLoops.mu(1, 3);
  mstop(2)              = forLoops.mu(2, 3);

  /// define Higgs vector in 't-Hooft Feynman gauge, and couplings:
  DoubleVector higgsm(4), higgsc(2);
  DoubleVector dnu(4), dnd(4), cn(4);
  assignHiggsSfermions(higgsm, higgsc, dnu, dnd, cn, beta);
  setNeutCurrCouplings(sinthDrbar, sw2, guL, gdL, geL, guR, gdR, geR);

  /// Higgs-sfermion-sfermion couplings
  DoubleMatrix lsStopLStopLR(4, 2), lsStopLStop12(4, 2);
  DoubleMatrix lsStopRStopLR(4, 2), lsStopRStop12(4, 2);

  /// LCT: CP-{even,odd}-sfermion-sfermion Feynman rules
  /// Order (s1 s2 G A, L R)
  lsStopLStopLR(1, 1) = g * mz * guL * cosb / costhDrbar;
  lsStopLStopLR(1, 2) = ht * smu / root2;
  lsStopLStopLR(3, 2) = 1.0 / root2 * 
    (smu * cosb * ht + forLoops.ut * sinb);
  lsStopLStopLR(4, 2) = -1.0 / root2 * 
    (smu * sinb * ht - forLoops.ut * cosb);
  lsStopLStopLR(2, 1) = -g * mz * guL * sinb / costhDrbar
    + htsq * v2;
  lsStopLStopLR(2, 2) = forLoops.ut / root2;

  lsStopRStopLR(1, 1) = lsStopLStopLR(1, 2);
  lsStopRStopLR(1, 2) = g * mz * guR * cosb / costhDrbar;
  lsStopRStopLR(2, 1) = lsStopLStopLR(2, 2);
  lsStopRStopLR(2, 2) = -g * mz * guR * sinb / costhDrbar
    + htsq * v2;
  lsStopRStopLR(3, 1) = -lsStopLStopLR(3, 2);
  lsStopRStopLR(4, 1) = -lsStopLStopLR(4, 2);

  /// Mix stops up
  int i;
  DoubleVector temp(2), temp2(2);
  for (i=1; i<=4; i++) {
    temp(1) = lsStopLStopLR(i, 1);
    temp(2) = lsStopLStopLR(i, 2);
    temp2 = rot2d(thetat) * temp;
    lsStopLStop12(i, 1) = temp2(1);
    lsStopLStop12(i, 2) = temp2(2);
    temp(1) = lsStopRStopLR(i, 1);
    temp(2) = lsStopRStopLR(i, 2);
    temp2 = rot2d(thetat) * temp;
    lsStopRStop12(i, 1) = temp2(1);
    lsStopRStop12(i, 2) = temp2(2);
 }

  DoubleMatrix lHStopLStop12(lsStopLStop12), lHStopRStop12(lsStopRStop12);
  /// Mix CP-even Higgses up
  for (i=1; i<=2; i++) { /// i is the s1/s2 label
    temp(1) = lsStopLStop12(1, i);
    temp(2) = lsStopLStop12(2, i);
    temp2 = rot2d(alpha) * temp;
    lHStopLStop12(1, i) = temp2(1);
    lHStopLStop12(2, i) = temp2(2);
    temp(1) = lsStopRStop12(1, i);
    temp(2) = lsStopRStop12(2, i);
    temp2 = rot2d(alpha) * temp;
    lHStopRStop12(1, i) = temp2(1);
    lHStopRStop12(2, i) = temp2(2);
  }

  /// Charged Higgs Feynman rules  (H+ G+, L R) basis
  DoubleMatrix lChHstopLsbotLR(2, 2), lChHstopRsbotLR(2, 2);
  lChHstopLsbotLR(1, 1) = (g * displayMwRun() * sin(2.0 * beta) 
    - htsq * v2 * cosb - hbsq * v1 * sinb) / root2;
  lChHstopLsbotLR(1, 2) = (smu * hb * cosb - 
			   forLoops.ub * sinb);
  lChHstopLsbotLR(2, 1) = (-g * displayMwRun() * cos(2.0 * beta) 
    - htsq * v2 * sinb + hbsq * v1 * cosb) / root2;
  lChHstopLsbotLR(2, 2) = hb * smu * sinb + forLoops.ub * cosb;

  lChHstopRsbotLR(1, 1) = ht * smu * sinb - forLoops.ut * cosb;
  lChHstopRsbotLR(1, 2) = ht * hb * (- v1 * cosb - v2 * sinb) / root2;
  lChHstopRsbotLR(2, 1) = -ht * smu * cosb - forLoops.ut * sinb;

  /// LCT: Rotate sfermions to (1, 2) mass basis
  DoubleMatrix lChHstopLsbot12(2, 2);
  temp(1) = lChHstopLsbotLR(1, 1);
  temp(2) = lChHstopLsbotLR(1, 2);
  temp2 = rot2d(thetab) * temp;
  lChHstopLsbot12(1, 1) = temp2(1);
  lChHstopLsbot12(1, 2) = temp2(2);
  temp(1) = lChHstopLsbotLR(2, 1);
  temp(2) = lChHstopLsbotLR(2, 2);
  temp2 = rot2d(thetab) * temp;
  lChHstopLsbot12(2, 1) = temp2(1);
  lChHstopLsbot12(2, 2) = temp2(2);

  /// LCT: Rotate sfermions to (1, 2) mass basis
  DoubleMatrix lChHstopRsbot12(2, 2);
  temp(1) = lChHstopRsbotLR(1, 1);
  temp(2) = lChHstopRsbotLR(1, 2);
  temp2 = rot2d(thetab) * temp;
  lChHstopRsbot12(1, 1) = temp2(1);
  lChHstopRsbot12(1, 2) = temp2(2);
  temp(1) = lChHstopRsbotLR(2, 1);
  temp(2) = lChHstopRsbotLR(2, 2);
  temp2 = rot2d(thetab) * temp;
  lChHstopRsbot12(2, 1) = temp2(1);
  lChHstopRsbot12(2, 2) = temp2(2);

  /// LCT: Contributions start here
  for (i=3; i<=4; i++) { 
    double a0p = a0(higgsc(i - 2), q);
    higgs(1, 1) += 
      (hbsq * dnu(i) + sqr(g) * (guL * 0.5 / costhDrbar2 - 0.5) * cn(i))* a0p;
    higgs(2, 2) += 
      (htsq * dnd(i) + sqr(g) * guR * 0.5 / costhDrbar2 * cn(i))* a0p;    
  }

  for (i=1; i<=4; i++) {
     higgs(1, 1) += 
        0.5 * (htsq * dnu(i) - sqr(g) * guL * 0.5 / sqr(costhDrbar) * cn(i)) 
        * a0(higgsm(i), q);
     higgs(2, 2) += 
        0.5 * (htsq * dnu(i) - sqr(g) * guR * 0.5 / sqr(costhDrbar) * cn(i)) 
        * a0(higgsm(i), q);
  }

  int j; for(i=1; i<=4; i++)
    for (j=1; j<=2; j++) {
      double b0p = b0(p, higgsm(i), mstop(j), q);
      higgs(1, 1) += sqr(lHStopLStop12(i, j)) * b0p;
      higgs(1, 2) += 
  	lHStopLStop12(i, j) * lHStopRStop12(i, j) * b0p;
      higgs(2, 2) += sqr(lHStopRStop12(i, j)) * b0p;
    }

  for(i=1; i<=2; i++)
    for (j=1; j<=2; j++) {
      double b0p = b0(p, higgsc(i), msbot(j), q); 
      higgs(1, 1) += sqr(lChHstopLsbot12(i, j)) * b0p;
      higgs(1, 2) += 
        lChHstopLsbot12(i, j) * lChHstopRsbot12(i, j) * b0p;
      higgs(2, 2) = higgs(2, 2) + sqr(lChHstopRsbot12(i, j)) * b0p;
    }

  return higgs;
}

template<class SoftPars>
DoubleMatrix Softsusy<SoftPars>::addStopEweak(double p, DoubleMatrix & electroweak) {
  double    sinthDrbar  = calcSinthdrbar();
  double    costhDrbar  = sqrt(1.0 - sqr(sinthDrbar));
  double    costhDrbar2 = 1.0 - sqr(sinthDrbar);
  double    g           = displayGaugeCoupling(2);
  double    gp          = displayGaugeCoupling(1) * sqrt(0.6);
  double    e           = g * sinthDrbar; /// DRbar value of e
  double    thetat      = forLoops.thetat;
  double    thetab      = forLoops.thetab;
  double    thetatau    = forLoops.thetatau;
  double    ct          = cos(thetat);
  double    st          = sin(thetat);
  double    cb          = cos(thetab);
  double    sb          = sin(thetab);
  double    ctau        = cos(thetatau);
  double    stau        = sin(thetatau);
  double    q           = displayMu();
  double    mz          = displayMzRun();
  DoubleVector mstop(2);
  mstop(1)              = forLoops.mu(1, 3);
  mstop(2)              = forLoops.mu(2, 3);
  double    msbot1      = forLoops.md(1, 3), msbot2 = forLoops.md(2, 3);

  electroweak(1, 1) += 
    4.0 * sqr(g) / costhDrbar2 * sqr(guL) * a0(mz, q) + 
    2.0 * sqr(g) * a0(displayMwRun(), q) + sqr(2.0 / 3.0 * e) * 
    (sqr(ct) * ffn(p, mstop(1), 0.0, q) + sqr(st) * ffn(p, mstop(2), 0.0, q))
    + sqr(g * guL / costhDrbar) * 
    (sqr(ct) * ffn(p, mstop(1), mz, q) + sqr(st) * ffn(p, mstop(2), mz, q)) +
    sqr(g) * 0.5 * (sqr(cb) * ffn(p, msbot1, displayMwRun(), q) + sqr(sb) *
  		    ffn(p, msbot2, displayMwRun(), q)) +
    sqr(g) * 0.25 * 
    (sqr(ct) * a0(mstop(1), q) + sqr(st) * a0(mstop(2), q) + 2.0 *
     (sqr(cb) * a0(msbot1, q) + sqr(sb) * a0(msbot2, q))) +
    sqr(g) * 0.5 * 
    (1.5 * a0(forLoops.mu(1, 1), q) + 
     1.5 * a0(forLoops.mu(1, 2), q) +
     1.5 * (sqr(ct) * a0(forLoops.mu(1, 3), q) + 
  	    sqr(st) * a0(forLoops.mu(2, 3), q)) -
     1.5 * a0(forLoops.md(1, 1), q) -
     1.5 * a0(forLoops.md(1, 2), q) -
     1.5 * (sqr(cb) * a0(forLoops.md(1, 3), q) +
  	    sqr(sb) * a0(forLoops.md(2, 3), q)) +
     0.5 * (a0(forLoops.msnu(1), q) + 
  	    a0(forLoops.msnu(2), q) +
  	    a0(forLoops.msnu(3), q)) -
     0.5 * (a0(forLoops.me(1, 1), q) + a0(forLoops.me(1, 2), q) +
  	    sqr(ctau) * a0(forLoops.me(1, 3), q) +
  	    sqr(stau) * a0(forLoops.me(2, 3), q))) +
    sqr(gp) * 0.25 * sqr(yuL) * 
    (sqr(ct) * a0(mstop(1), q) + sqr(st) * a0(mstop(2), q)) +
    sqr(gp) * 0.25 * yuL * 
    (3.0 * yuL * (a0(forLoops.mu(1, 1), q) + 
  		  a0(forLoops.mu(1, 2), q) + 
  		  sqr(ct) * a0(forLoops.mu(1, 3), q) + 
  		  sqr(st) * a0(forLoops.mu(2, 3), q)) +
     3.0 * yuR * (a0(forLoops.mu(2, 1), q) + 
  		  a0(forLoops.mu(2, 2), q) + 
  		  sqr(st) * a0(forLoops.mu(1, 3), q) + 
  		  sqr(ct) * a0(forLoops.mu(2, 3), q)) +
     3.0 * ydL * (a0(forLoops.md(1, 1), q) + 
  		  a0(forLoops.md(1, 2), q) + 
  		  sqr(cb) * a0(forLoops.md(1, 3), q) + 
  		  sqr(sb) * a0(forLoops.md(2, 3), q)) +
     3.0 * ydR * (a0(forLoops.md(2, 1), q) + 
  		  a0(forLoops.md(2, 2), q) + 
  		  sqr(sb) * a0(forLoops.md(1, 3), q) + 
  		  sqr(cb) * a0(forLoops.md(2, 3), q)) +
     yeL * (a0(forLoops.me(1, 1), q) + 
  	    a0(forLoops.me(1, 2), q) + 
  	    sqr(ctau) * a0(forLoops.me(1, 3), q) + 
  	    sqr(stau) * a0(forLoops.me(2, 3), q)) +
     yeR * (a0(forLoops.me(2, 1), q) + 
  	    a0(forLoops.me(2, 2), q) + 
  	    sqr(stau) * a0(forLoops.me(1, 3), q) + 
  	    sqr(ctau) * a0(forLoops.me(2, 3), q)) +
     ynuL * (a0(forLoops.msnu(1), q) + 
  	     a0(forLoops.msnu(2), q) +
  	     a0(forLoops.msnu(3), q)));
     
  electroweak(2, 2) += 
    4.0 * sqr(g) / costhDrbar2 * sqr(guR) * a0(mz, q) + 
    sqr(2.0 / 3.0 * e) * 
    (sqr(st) * ffn(p, mstop(1), 0.0, q) + sqr(ct) * ffn(p, mstop(2), 0.0, q))
    + sqr(g * guR / costhDrbar) * 
    (sqr(st) * ffn(p, mstop(1), mz, q) + sqr(ct) * ffn(p, mstop(2), mz, q)) +
    sqr(gp) * 0.25 * sqr(yuR) * 
    (sqr(st) * a0(mstop(1), q) + sqr(ct) * a0(mstop(2), q)) +
    sqr(gp) * 0.25 * yuR * 
    (3.0 * yuL * (a0(forLoops.mu(1, 1), q) + 
  		  a0(forLoops.mu(1, 2), q) + 
  		  sqr(ct) * a0(forLoops.mu(1, 3), q) + 
  		  sqr(st) * a0(forLoops.mu(2, 3), q)) +
     3.0 * yuR * (a0(forLoops.mu(2, 1), q) + 
  		  a0(forLoops.mu(2, 2), q) + 
  		  sqr(st) * a0(forLoops.mu(1, 3), q) + 
  		  sqr(ct) * a0(forLoops.mu(2, 3), q)) +
     3.0 * ydL * (a0(forLoops.md(1, 1), q) + 
  		  a0(forLoops.md(1, 2), q) + 
  		  sqr(cb) * a0(forLoops.md(1, 3), q) + 
  		  sqr(sb) * a0(forLoops.md(2, 3), q)) +
     3.0 * ydR * (a0(forLoops.md(2, 1), q) + 
  		  a0(forLoops.md(2, 2), q) + 
  		  sqr(sb) * a0(forLoops.md(1, 3), q) + 
  		  sqr(cb) * a0(forLoops.md(2, 3), q)) +
     yeL * (a0(forLoops.me(1, 1), q) + 
  	    a0(forLoops.me(1, 2), q) + 
  	    sqr(ctau) * a0(forLoops.me(1, 3), q) + 
  	    sqr(stau) * a0(forLoops.me(2, 3), q)) +
     yeR * (a0(forLoops.me(2, 1), q) + 
  	    a0(forLoops.me(2, 2), q) + 
  	    sqr(stau) * a0(forLoops.me(1, 3), q) + 
  	    sqr(ctau) * a0(forLoops.me(2, 3), q)) +
     ynuL * (a0(forLoops.msnu(1), q) + 
  	     a0(forLoops.msnu(2), q) +
  	     a0(forLoops.msnu(3), q)));

  electroweak(1, 2) += 
    sqr(gp) * 0.25 * yuL * yuR * st * ct *
    (a0(mstop(1), q) - a0(mstop(2), q)) +
    sqr(2.0 / 3.0 * e) * st * ct * 
    (ffn(p, mstop(1), 0.0, q) - ffn(p, mstop(2), 0.0, q)) -
    sqr(g) / costhDrbar2 * guL * guR * st * ct *
    (ffn(p, mstop(1), mz, q) - ffn(p, mstop(2), mz, q));

  return electroweak;
}

template<class SoftPars>
DoubleMatrix Softsusy<SoftPars>::addStopNeutralino(double p, double mt, DoubleMatrix & neutralino) {
  double q  = displayMu();
  double g  = displayGaugeCoupling(2);
  double gp = displayGaugeCoupling(1) * sqrt(0.6);
  double ht = forLoops.ht;

  ComplexMatrix n(forLoops.nBpmz);
  DoubleVector mneut(forLoops.mnBpmz);

  int rank = mneut.displayEnd();

  /// Neutralino Feynman rules
  DoubleVector aPsi0TStopr(rank), bPsi0TStopr(rank), 
    aPsi0TStopl(rank), bPsi0TStopl(rank);
  aPsi0TStopr(1) = - 4.0 * gp / (3.0 * root2);
  bPsi0TStopl(1) = gp / (3.0 * root2);
  bPsi0TStopl(2) = g / root2;
  aPsi0TStopl(4) = ht;
  bPsi0TStopr(4) = ht;

  ComplexVector aChi0TStopl(rank), bChi0TStopl(rank), 
    aChi0TStopr(rank), bChi0TStopr(rank);
  aChi0TStopl = n.complexConjugate() * aPsi0TStopl;
  bChi0TStopl = n * bPsi0TStopl;
  aChi0TStopr = n.complexConjugate() * aPsi0TStopr;
  bChi0TStopr = n * bPsi0TStopr;

  DoubleVector gChi0TopStopLL(rank), fChi0TopStopLL(rank);
  DoubleVector gChi0TopStopLR(rank), fChi0TopStopLR(rank);
  DoubleVector gChi0TopStopRR(rank), fChi0TopStopRR(rank);

  int i;
  for (i=1; i<=rank; i++) {
    fChi0TopStopLL(i) = (aChi0TStopl(i) * aChi0TStopl(i).conj() + 
      bChi0TStopl(i) * bChi0TStopl(i).conj()).real();
    gChi0TopStopLL(i) = (bChi0TStopl(i).conj() * aChi0TStopl(i) + 
      bChi0TStopl(i) * aChi0TStopl(i).conj()).real();
    fChi0TopStopRR(i) = (aChi0TStopr(i) * aChi0TStopr(i).conj() + 
      bChi0TStopr(i) * bChi0TStopr(i).conj()).real();
    gChi0TopStopRR(i) = (bChi0TStopr(i).conj() * aChi0TStopr(i) + 
      bChi0TStopr(i) * aChi0TStopr(i).conj()).real();
    fChi0TopStopLR(i) = (aChi0TStopr(i) * aChi0TStopl(i).conj() + 
      bChi0TStopr(i) * bChi0TStopl(i).conj()).real();
    gChi0TopStopLR(i) = (bChi0TStopl(i).conj() * aChi0TStopr(i) + 
      bChi0TStopr(i) * aChi0TStopl(i).conj()).real();
  }

  for (i=1; i<=rank; i++) {
    double one = gfn(p, mneut(i), mt, q);
    double two = 2.0 * mneut(i) * mt * b0(p, mneut(i), mt, q);
    neutralino(1, 1) = neutralino(1, 1) +
      fChi0TopStopLL(i) * one - gChi0TopStopLL(i) * two;
    neutralino(2, 2) = neutralino(2, 2) +
      fChi0TopStopRR(i) * one - gChi0TopStopRR(i) * two;
    neutralino(1, 2) = neutralino(1, 2) +
      fChi0TopStopLR(i) * one - gChi0TopStopLR(i) * two;
  }

  return neutralino;
}

template<class SoftPars>
DoubleMatrix Softsusy<SoftPars>::addStopChargino(double p, DoubleMatrix & chargino) {
  double g  = displayGaugeCoupling(2);
  double ht = forLoops.ht;
  double hb = forLoops.hb;
  double mb = forLoops.mb;
  double q  = displayMu();

  /// Chargino Feynman Rules
  DoubleVector bPsicBStopl(2), bPsicBStopr(2), aPsicBStopl(2), aPsicBStopr(2);
  aPsicBStopl(1) = g;
  aPsicBStopr(2) = -ht;
  bPsicBStopl(2) = -hb;
  
  /// LCT: Define mixing matrices and mass vector
  ComplexMatrix u(forLoops.uBpmz), v(forLoops.vBpmz); 
  DoubleVector mch(forLoops.mchBpmz); 

  DoubleVector aPsicCSbotl(2);
  ComplexVector aChicBStopr(2), aChicBStopl(2), bChicBStopl(2), bChicBStopr(2);
  ComplexMatrix aChBStop(2, 2), bChBStop(2, 2);

  aChicBStopl = v.complexConjugate() * aPsicBStopl;
  bChicBStopl = u * bPsicBStopl;
  aChicBStopr = v.complexConjugate() * aPsicBStopr;
  bChicBStopr = u * bPsicBStopr;

  DoubleVector fChBStopLL(2), gChBStopLL(2) ;
  DoubleVector fChBStopLR(2), gChBStopLR(2); 
  DoubleVector fChBStopRR(2), gChBStopRR(2); 

  int i;
  for (i=1; i<=2; i++) {
    fChBStopLL(i) = (aChicBStopl(i).conj() * aChicBStopl(i) +
		      bChicBStopl(i).conj() * bChicBStopl(i)).real();
    gChBStopLL(i) = (bChicBStopl(i).conj() * aChicBStopl(i) +
		      aChicBStopl(i).conj() * bChicBStopl(i)).real();
    fChBStopLR(i) = (aChicBStopl(i).conj() * aChicBStopr(i) +
		      bChicBStopl(i).conj() * bChicBStopr(i)).real();
    gChBStopLR(i) = (bChicBStopl(i).conj() * aChicBStopr(i) +
		      aChicBStopl(i).conj() * bChicBStopr(i)).real();
    fChBStopRR(i) = (aChicBStopr(i).conj() * aChicBStopr(i) +
		      bChicBStopr(i).conj() * bChicBStopr(i)).real();
    gChBStopRR(i) = (bChicBStopr(i).conj() * aChicBStopr(i) +
		      aChicBStopr(i).conj() * bChicBStopr(i)).real();
  }

  for (i=1; i<=2; i++) {
    double one = gfn(p, mch(i), mb, q);
    double two = mch(i) * mb * b0(p, mch(i), mb, q) * 2.0;
    chargino(1, 1) = chargino(1, 1) + fChBStopLL(i) * one -
      gChBStopLL(i) * two;
    chargino(1, 2) = chargino(1, 2) + fChBStopLR(i) * one -
      gChBStopLR(i) * two;
    chargino(2, 2) = chargino(2, 2) + fChBStopRR(i) * one - 
      gChBStopRR(i) * two;
  }
  return chargino;
}

/// As in BPMZ appendix, INCLUDING weak boson loops.
template<class SoftPars>
void Softsusy<SoftPars>::addStopCorrection(double p, DoubleMatrix & mass, 
				     double mt) {

/// No point adding radiative corrections to tachyonic particles
  if (mass(1, 1) < 0.0 || mass(2, 2) < 0.0) { 
    flagTachyon(stop);
    if (mass(1, 1) < 0.) mass(1, 1) = EPSTOL;
    else mass(2, 2) = EPSTOL;
    return;
  }

  /// one-loop correction matrix
  DoubleMatrix piSq(2, 2); /// Self-energy matrix

  /// Corrections themselves start here
  DoubleMatrix strong(2, 2), stop(2, 2), sbottom(2, 2), 
    higgs(2, 2), electroweak(2, 2), chargino(2, 2), neutralino(2, 2);
  /// LCT: Corrections from strong interactions
  addStopQCD(p, mt, strong);
  /// LCT: Corrections from stops
  addStopStop(p, mt, stop);
  /// LCT: Corrections from sbottoms
  addStopSbottom(p, mt, sbottom);
  /// LCT: Corrections from Higgses
  addStopHiggs(p, mt, higgs);
  /// LCT: Electroweak corrections
  addStopEweak(p, electroweak);
  /// LCT: Chargino contribution
  addStopChargino(p, chargino);
  /// LCT: Neutralino contribution
  addStopNeutralino(p, mt, neutralino);
  
  piSq = 1.0 / (16.0 * sqr(PI)) * 
    (strong + stop + sbottom + higgs + electroweak + chargino + neutralino);

  piSq(2, 1) = piSq(1, 2);
  mass = mass - piSq;	  
}

template<class SoftPars>
void Softsusy<SoftPars>::assignHiggs(DoubleVector & higgsm, DoubleVector & higgsc)
  const {
  drBarPars f(displayDrBarPars());

  higgsm(1) = f.mh0(2);
  higgsm(2) = f.mh0(1);
  higgsm(3) = displayMzRun();
  higgsm(4) = f.mA0(1);
  higgsc(1) = displayMwRun();
  higgsc(2) = f.mHpm;
}

template<class SoftPars>
void Softsusy<SoftPars>::assignHiggs(DoubleVector & higgsm, DoubleVector & higgsc, 
			       DoubleVector & dnu, DoubleVector & dnd, 
			       DoubleVector & cn, double beta) const {
  double sinb = sin(beta), cosb = cos(beta);

  assignHiggs(higgsm, higgsc);

  cn(1)  = - cos(2.0 * forLoops.thetaH);
  cn(2)  = - cn(1);
  cn(3)  = - cos(2.0 * beta);
  cn(4)  = - cn(3);
  dnu(1) = sqr(sin(displayDrBarPars().thetaH));
  dnu(2) = sqr(cos(displayDrBarPars().thetaH));
  dnu(3) = sqr(sinb);
  dnu(4) = sqr(cosb);
  dnd(1) = dnu(2);
  dnd(2) = dnu(1);
  dnd(3) = dnu(4);
  dnd(4) = dnu(3);
}

/// some switches due to BPMZ's different conventions
template<class SoftPars>
void Softsusy<SoftPars>::assignHiggsSfermions(DoubleVector & higgsm, 
					DoubleVector & higgsc, 
					DoubleVector & dnu, 
					DoubleVector & dnd, 
					DoubleVector & cn, double beta) const {

  assignHiggs(higgsm, higgsc, dnu, dnd, cn, beta);

  /// swap ordering of charged Higgs' for sfermions due to BPMZ conventions
  double higgsc1 = higgsc(1);
  higgsc(1) = higgsc(2);
  higgsc(2) = higgsc1;
}

template<class SoftPars>
DoubleMatrix Softsusy<SoftPars>::addSlepHiggs(double p1, double p2, int family, DoubleMatrix & higgs) {
  double    mz         = displayMzRun();
  double    sinthDrbar = calcSinthdrbar();
  double    costhDrbar = sqrt(1.0 - sqr(sinthDrbar));
  double    alpha      = forLoops.thetaH;
  double    g	       = displayGaugeCoupling(2);
  double    beta       = atan(displayTanb());
  double    sinb       = sin(beta), cosb = cos(beta);
  double    q          = displayMu();
  DoubleVector msel(2);
  msel(1)              = forLoops.me(1, family);
  msel(2)              = forLoops.me(2, family);
  DoubleVector msnu(3);
  msnu(1)              = forLoops.msnu(1);
  msnu(2)              = forLoops.msnu(2);
  msnu(3)              = forLoops.msnu(3);

  /// Define Higgs vector of masses in 't-Hooft Feynman gauge, and couplings:
  DoubleVector higgsm(4), higgsc(2), dnu(4), dnd(4), cn(4);
  assignHiggsSfermions(higgsm, higgsc, dnu, dnd, cn, beta);

  /// stau now stands for generation 1 or 2
  DoubleMatrix lsStauLStauLR(4, 2), lsStauLStau12(4, 2);
  DoubleMatrix lsStauRStauLR(4, 2), lsStauRStau12(4, 2);
  /// Order (s1 s2 G A, L R)
  lsStauLStauLR(1, 1) = g * mz * geL * cosb / costhDrbar;
  lsStauLStauLR(2, 1) = -g * mz * geL * sinb / costhDrbar;

  lsStauRStauLR(1, 2) = g * mz * geR * cosb / costhDrbar;
  lsStauRStauLR(2, 2) = -g * mz * geR * sinb / costhDrbar;
  lsStauRStauLR(3, 1) = -lsStauLStauLR(3, 2);
  lsStauRStauLR(4, 1) = -lsStauLStauLR(4, 2);

  DoubleMatrix lHStauLStau12(lsStauLStauLR), lHStauRStau12(lsStauRStauLR);
  DoubleVector temp(2), temp2(2);
  /// Mix CP-even Higgses up
  int i; for (i=1; i<=2; i++) { /// i is the s1/s2 label
    temp(1) = lsStauLStauLR(1, i);
    temp(2) = lsStauLStauLR(2, i);
    temp2 = rot2d(alpha) * temp;
    lHStauLStau12(1, i) = temp2(1);
    lHStauLStau12(2, i) = temp2(2);
    temp(1) = lsStauRStauLR(1, i);
    temp(2) = lsStauRStauLR(2, i);
    temp2 = rot2d(alpha) * temp;
    lHStauRStau12(1, i) = temp2(1);
    lHStauRStau12(2, i) = temp2(2);
  }

  /// Charged Higgs Feynman rules
  DoubleMatrix lChHstauLsnu12(2, 2); /// (H+ G+, L R) basis
  lChHstauLsnu12(2, 1) = -g * displayMwRun() * cos(2.0 * beta) / root2;
  lChHstauLsnu12(1, 1) = g * displayMwRun() * sin(2.0 * beta) / root2;

 for (i=1; i<=4; i++) {
    higgs(1, 1) = higgs(1, 1) + 
      0.5 * ( - sqr(g) * geL * cn(i) 
	     / (2.0 * sqr(costhDrbar))) * a0(higgsm(i), q);
    higgs(2, 2) = higgs(2, 2) + 
      0.5 * (- sqr(g) * geR * cn(i) 
	     / (2.0 * sqr(costhDrbar))) * a0(higgsm(i), q);
  }

  for (i=3; i<=4; i++) {
    double a0p = a0(higgsc(i - 2), q);
    higgs(1, 1) = higgs(1, 1) + a0p * 
	  (sqr(g) * (geL / (2.0 * sqr(costhDrbar)) + 0.5) * cn(i));
    higgs(2, 2) = higgs(2, 2) + a0p *
      (sqr(g) * geR / (2.0 * sqr(costhDrbar)) * cn(i));
  }

  int j; for(i=1; i<=4; i++)
    for (j=1; j<=2; j++) {
      double b0p1 = b0(p1, higgsm(i), msel(j), q);
      double b0p2 = b0(p2, higgsm(i), msel(j), q);
      higgs(1, 1) = higgs(1, 1) + sqr(lHStauLStau12(i, j)) * b0p1;
      higgs(2, 2) = higgs(2, 2) + sqr(lHStauRStau12(i, j)) * b0p2;
    }

  for(i=1; i<=2; i++) {
      double b0p = b0(p1, msnu(family), higgsc(i), q); 
      higgs(1, 1) = higgs(1, 1) + sqr(lChHstauLsnu12(i, 1)) * b0p;
    }

  return higgs;
}

template<class SoftPars>
DoubleMatrix Softsusy<SoftPars>::addSlepEweak(double p1, double p2, int family, DoubleMatrix & electroweak) {
  double    mw         = displayMwRun();
  double    mz         = displayMzRun();
  double    sinthDrbar = calcSinthdrbar();
  double    costhDrbar = sqrt(1.0 - sqr(sinthDrbar));
  double    g	       = displayGaugeCoupling(2);
  double    gp         = displayGaugeCoupling(1) * sqrt(0.6);
  double    e          = g * sinthDrbar;
  double    thetat     = forLoops.thetat;
  double    thetab     = forLoops.thetab;
  double    thetatau   = forLoops.thetatau;
  double    ct         = cos(thetat);
  double    st         = sin(thetat);
  double    cb         = cos(thetab);
  double    sb         = sin(thetab);
  double    ctau       = cos(thetatau);
  double    stau       = sin(thetatau);
  double    q          = displayMu();
  DoubleVector msbot(2);
  msbot(1)             = forLoops.md(1, 3);
  msbot(2)             = forLoops.md(2, 3);
  DoubleVector mstop(2);
  mstop(1)             = forLoops.mu(1, 3);
  mstop(2)             = forLoops.mu(2, 3);
  DoubleVector mstau(2);
  mstau(1)             = forLoops.me(1, 3);
  mstau(2)             = forLoops.me(2, 3);
  DoubleVector msel(2);
  msel(1)              = forLoops.me(1, family);
  msel(2)              = forLoops.me(2, family);
  DoubleVector msnu(3);
  msnu(1)              = forLoops.msnu(1);
  msnu(2)              = forLoops.msnu(2);
  msnu(3)              = forLoops.msnu(3);

  /// LCT: EW contributions
  electroweak(1, 1) = electroweak(1, 1) + 
    4.0 * sqr(g) / sqr(costhDrbar) * sqr(geL) * a0(mz, q) +
    2.0 * sqr(g) * a0(mw, q) + sqr(e) * 
    ffn(p1, msel(1), 0., q);
  electroweak(1, 1) = electroweak(1, 1) +   
    sqr(g) / sqr(costhDrbar) * sqr(geL) * 
    ffn(p1, msel(1), mz, q) +
    sqr(g) * 0.5 * ffn(p1, msnu(family), mw, q);

  electroweak(1, 1) = electroweak(1, 1) +   
    sqr(g) * 0.25 * (a0(msel(1), q) + 2.0 * a0(msnu(family), q));
  electroweak(1, 1) = electroweak(1, 1) + sqr(g) * 
    (-3.0 * 0.25 * (sqr(ct) * a0(mstop(1), q) + sqr(st) * a0(mstop(2), q))
     +3.0 * 0.25 * (sqr(cb) * a0(msbot(1), q) + sqr(sb) * a0(msbot(2), q))
     +      0.25 * (sqr(ctau) * a0(mstau(1), q) + sqr(stau) * a0(mstau(2), q))
     -      0.25 * a0(msnu(3), q)
     -3.0 * 0.25 * (a0(forLoops.mu(1, 2), q) + a0(forLoops.mu(1, 1), q))
     +3.0 * 0.25 * (a0(forLoops.md(1, 2), q) + a0(forLoops.md(1, 1), q))
     +      0.25 * (a0(forLoops.me(1, 2), q) + a0(forLoops.me(1, 1), q))
     -      0.25 * (a0(msnu(2), q) + a0(msnu(1), q)))
    + sqr(gp) * 0.25 * sqr(yeL) * a0(msel(1), q);
  electroweak(1, 1) = electroweak(1, 1) + sqr(gp) * 0.25 * yeL *
    (+3.0 * yuL * (sqr(ct) * a0(mstop(1), q) + sqr(st) * a0(mstop(2), q))
     +3.0 * ydL * (sqr(cb) * a0(msbot(1), q) + sqr(sb) * a0(msbot(2), q))
     +      yeL * (sqr(ctau) * a0(mstau(1), q) + sqr(stau) * a0(mstau(2), q))
     +      ynuL * a0(msnu(3), q)
     +3.0 * yuL * (a0(forLoops.mu(1, 2), q) + a0(forLoops.mu(1, 1), q))
     +3.0 * ydL * (a0(forLoops.md(1, 2), q) + a0(forLoops.md(1, 1), q))
     +      yeL * (a0(forLoops.me(1, 2), q) + a0(forLoops.me(1, 1), q))
     +      ynuL * (a0(msnu(1), q) + a0(msnu(2), q))
     //
     +3.0 * yuR * (sqr(st) * a0(mstop(1), q) + sqr(ct) * a0(mstop(2), q))
     +3.0 * ydR * (sqr(sb) * a0(msbot(1), q) + sqr(cb) * a0(msbot(2), q))
     +yeR * (sqr(stau) * a0(mstau(1), q) + sqr(ctau) * a0(mstau(2), q))
     +3.0 * yuR * (a0(forLoops.mu(2, 2), q) + a0(forLoops.mu(2, 1), q))
     +3.0 * ydR * (a0(forLoops.md(2, 2), q) + a0(forLoops.md(2, 1), q))
     +      yeR * (a0(forLoops.me(2, 2), q) + a0(forLoops.me(2, 1), q)));


  electroweak(2, 2) = electroweak(2, 2) + 
    4.0 * sqr(g) / sqr(costhDrbar) * sqr(geR) * a0(mz, q) +
    sqr(e) * ffn(p2, msel(2), 0., q);
  electroweak(2, 2) = electroweak(2, 2) + 
    sqr(g) / sqr(costhDrbar) * sqr(geR) * ffn(p2, msel(2), mz, q);
  electroweak(2, 2) = electroweak(2, 2) + 
    sqr(gp) * 0.25 * sqr(yeR) * a0(msel(2), q);
  electroweak(2, 2) = electroweak(2, 2) + sqr(gp) * 0.25 * yeR *
    (+3.0 * yuL * (sqr(ct) * a0(mstop(1), q) + sqr(st) * a0(mstop(2), q))
     +3.0 * ydL * (sqr(cb) * a0(msbot(1), q) + sqr(sb) * a0(msbot(2), q))
     +      yeL * (sqr(ctau) * a0(mstau(1), q) + sqr(stau) * a0(mstau(2), q))
     +      ynuL * a0(msnu(3), q) 
     +3.0 * yuL * (a0(forLoops.mu(1, 2), q) + a0(forLoops.mu(1, 1), q))
     +3.0 * ydL * (a0(forLoops.md(1, 2), q) + a0(forLoops.md(1, 1), q))
     +      yeL * (a0(forLoops.me(1, 2), q) + a0(forLoops.me(1, 1), q))
     +      ynuL * (a0(msnu(1), q) + a0(msnu(2), q))
     //
     +3.0 * yuR * (sqr(st) * a0(mstop(1), q) + sqr(ct) * a0(mstop(2), q))
     +3.0 * ydR * (sqr(sb) * a0(msbot(1), q) + sqr(cb) * a0(msbot(2), q))
     +      yeR * (sqr(stau) * a0(mstau(1), q) + sqr(ctau) * a0(mstau(2), q))
     +3.0 * yuR * (a0(forLoops.mu(2, 2), q) + a0(forLoops.mu(2, 1), q))
     +3.0 * ydR * (a0(forLoops.md(2, 2), q) + a0(forLoops.md(2, 1), q))
     +      yeR * (a0(forLoops.me(2, 2), q) + a0(forLoops.me(2, 1), q)));

  return electroweak;
}

template<class SoftPars>
void Softsusy<SoftPars>::addSlepGaugino(double p1, double p2, int /* family */, DoubleMatrix & chargino, DoubleMatrix & neutralino) {
  double    g	    = displayGaugeCoupling(2);
  double    gp      = displayGaugeCoupling(1) * sqrt(0.6);
  double    q       = displayMu();

  /// Neutralino Feynman rules
  ComplexMatrix n(forLoops.nBpmz);
  DoubleVector mneut(forLoops.mnBpmz);
  ComplexMatrix u(forLoops.uBpmz), v(forLoops.vBpmz); 
  DoubleVector mch(forLoops.mchBpmz);

  const int rank = mneut.displayEnd();

  DoubleVector aPsi0TauStaur(rank), bPsi0TauStaur(rank), aPsi0TauStaul(rank),
    bPsi0TauStaul(rank); 
  aPsi0TauStaur(1) = gp * root2;
  bPsi0TauStaul(1) = -gp / root2;
  bPsi0TauStaul(2) = -g / root2;

  ComplexVector aChi0TauStaul(rank), bChi0TauStaul(rank), aChi0TauStaur(rank),
    bChi0TauStaur(rank); 

  aChi0TauStaul = n.complexConjugate() * aPsi0TauStaul;
  bChi0TauStaul = n * bPsi0TauStaul;
  aChi0TauStaur = n.complexConjugate() * aPsi0TauStaur;
  bChi0TauStaur = n * bPsi0TauStaur;

  DoubleVector gChi0TauStauLL(rank), fChi0TauStauLL(rank);
  DoubleVector gChi0TauStauRR(rank), fChi0TauStauRR(rank);
  int i;
  for (i=1; i<=rank; i++) {
    fChi0TauStauLL(i) = (aChi0TauStaul(i) * aChi0TauStaul(i).conj() + 
      bChi0TauStaul(i) * bChi0TauStaul(i).conj()).real();
    gChi0TauStauLL(i) = (bChi0TauStaul(i).conj() * aChi0TauStaul(i) + 
      bChi0TauStaul(i) * aChi0TauStaul(i).conj()).real();
    fChi0TauStauRR(i) = (aChi0TauStaur(i) * aChi0TauStaur(i).conj() + 
      bChi0TauStaur(i) * bChi0TauStaur(i).conj()).real();
    gChi0TauStauRR(i) = (bChi0TauStaur(i).conj() * aChi0TauStaur(i) + 
      bChi0TauStaur(i) * aChi0TauStaur(i).conj()).real();
  }

  /// Chargino Feynman Rules
  DoubleVector bPsicNuStaul(2);//, bPsicNuStaur(2);

  bPsicNuStaul(1) = g;
  
  ComplexVector bChicNuStaul(2);
  ComplexMatrix bChNuStau(2, 2);

  bChicNuStaul = u * bPsicNuStaul;

  DoubleVector fChNuStauLL(2) ;
  DoubleVector fChNuStauRR(2); 
  for (i=1; i<=2; i++) {
    fChNuStauLL(i) = (bChicNuStaul(i).conj() * bChicNuStaul(i)).real();
  }

  for (i=1; i<=2; i++) {
    double one = gfn(p1, mch(i), 0., q);
    chargino(1, 1) = chargino(1, 1) + fChNuStauLL(i) * one;
  }

  for (i=1; i<=rank; i++) {
    double one = gfn(p1, mneut(i), 0., q);
    double one1 = gfn(p2, mneut(i), 0., q);
    neutralino(1, 1) = neutralino(1, 1) + fChi0TauStauLL(i) * one;
    neutralino(2, 2) = neutralino(2, 2) + fChi0TauStauRR(i) * one1;
  }

}

/// 16.09.05 checked. 
template<class SoftPars>
void Softsusy<SoftPars>::addSlepCorrection(DoubleMatrix & mass, int family) {

/// No point adding radiative corrections to tachyonic particles
  if (mass(1, 1) < 0.0 || mass(2, 2) < 0.0) { 
    if (mass(1, 1) < 0.) mass(1, 1) = EPSTOL;
    else mass(2, 2) = EPSTOL;
    if (family == 1) flagTachyon(selectron);
    if (family == 2) flagTachyon(smuon);
    return;
  }

  /// one-loop correction matrix
  DoubleMatrix piSq(2, 2); /// Self-energy matrix

  DoubleVector msel(2);
  msel(1)           = forLoops.me(1, family);
  msel(2)           = forLoops.me(2, family);	
  double p1         = msel(1), p2 = msel(2);
  ///  double    p       = sqrt(msel(1) * msel(2));

  DoubleMatrix higgs(2, 2), chargino(2, 2), neutralino(2, 2), electroweak(2, 2);

  /// LCT: Higgs contribution
  addSlepHiggs(p1, p2, family, higgs);
  /// LCT: Electroweak contribution
  addSlepEweak(p1, p2, family, electroweak);
  /// LCT: Chargino and neutralino contribution
  addSlepGaugino(p1, p2, family, chargino, neutralino);

  piSq = 1.0 / (16.0 * sqr(PI)) * 
    (higgs + chargino + neutralino + electroweak);

  mass = mass - piSq;	
}

template<class SoftPars>
void Softsusy<SoftPars>::addStauSfermion(double /* p */, double /* mtau */, DoubleMatrix & stop, DoubleMatrix & sbottom) {
  double    thetab   = forLoops.thetab;
  double    thetatau = forLoops.thetatau;
  double    cb       = cos(thetab);
  double    sb       = sin(thetab);
  double    ctau     = cos(thetatau);
  double    stau     = sin(thetatau);
  double    q        = displayMu(); 
  double    hb       = forLoops.hb;
  double    htau     = forLoops.htau;
  double    htausq   = sqr(htau);
  DoubleVector msbot(2);
  msbot(1)           = forLoops.md(1, 3);
  msbot(2)           = forLoops.md(2, 3);
  DoubleVector mstau(2);
  mstau(1)           = forLoops.me(1, 3);
  mstau(2)           = forLoops.me(2, 3);
  DoubleVector msnu(3);
  msnu(1)            = forLoops.msnu(1);
  msnu(2)            = forLoops.msnu(2);
  msnu(3)            = forLoops.msnu(3);

  double     a0t1    = a0(mstau(1), q), a0t2 = a0(mstau(2), q);

  sbottom(2, 2) = htausq * a0(msnu(3), q);

  stop(1, 1) = htausq * (sqr(stau) * a0t1 + sqr(ctau) * a0t2);
  stop(2, 2) = htausq * (sqr(ctau) * a0t1 + sqr(stau) * a0t2);
  stop(1, 2) = htausq * ctau * stau * (a0t1 - a0t2)
     + 3.0 * htau * hb * cb * sb * (a0(msbot(1), q) - a0(msbot(2), q));
}

template<class SoftPars>
DoubleMatrix Softsusy<SoftPars>::addStauHiggs(double p, double mtau, DoubleMatrix & higgs) {
  double    mz         = displayMzRun();
  double    sinthDrbar = calcSinthdrbar();
  double    costhDrbar = sqrt(1.0 - sqr(sinthDrbar));
  double    alpha      = forLoops.thetaH;
  double    g	       = displayGaugeCoupling(2);
  double    beta       = atan(displayTanb());
  double    sinb       = sin(beta), cosb = cos(beta);
  double    smu        = -displaySusyMu();
  double    q          = displayMu();
  double    htau       = forLoops.htau, htausq = sqr(htau);
  double    thetatau   = forLoops.thetatau;
  DoubleVector mstau(2);
  mstau(1)             = forLoops.me(1, 3);
  mstau(2)             = forLoops.me(2, 3);
  DoubleVector msnu(3);
  msnu(1)              = forLoops.msnu(1);
  msnu(2)              = forLoops.msnu(2);
  msnu(3)              = forLoops.msnu(3);

  /// Define Higgs vector of masses in 't-Hooft Feynman gauge, and couplings:
  DoubleVector higgsm(4), higgsc(2);
  DoubleVector dnu(4), dnd(4), cn(4);
  assignHiggsSfermions(higgsm, higgsc, dnu, dnd, cn, beta);

  DoubleMatrix lsStauLStauLR(4, 2), lsStauLStau12(4, 2);
  DoubleMatrix lsStauRStauLR(4, 2), lsStauRStau12(4, 2);
  /// Order (s1 s2 G A, L R)
  lsStauLStauLR(1, 1) = g * mz * geL * cosb / costhDrbar + 
    root2 * htau * mtau;
  
  lsStauLStauLR(1, 2) = forLoops.utau / root2;
  lsStauLStauLR(3, 2) = -1.0 / root2 * 
    (smu * sinb * htau + forLoops.utau * cosb);
  lsStauLStauLR(4, 2) = -1.0 / root2 * 
    (smu * cosb * htau - forLoops.utau * sinb);
  lsStauLStauLR(2, 1) = -g * mz * geL * sinb / costhDrbar;
  lsStauLStauLR(2, 2) = htau * smu / root2;

  lsStauRStauLR(1, 1) = lsStauLStauLR(1, 2);
  lsStauRStauLR(1, 2) = g * mz * geR * cosb / costhDrbar +
    root2 * htau * mtau;
  lsStauRStauLR(2, 1) = lsStauLStauLR(2, 2);
  lsStauRStauLR(2, 2) = -g * mz * geR * sinb / costhDrbar;
  lsStauRStauLR(3, 1) = -lsStauLStauLR(3, 2);
  lsStauRStauLR(4, 1) = -lsStauLStauLR(4, 2);

  /// Mix staus up
  int i;
  DoubleVector temp(2), temp2(2);
  for (i=1; i<=4; i++) {
    temp(1) = lsStauLStauLR(i, 1);
    temp(2) = lsStauLStauLR(i, 2);
    temp2 = rot2d(thetatau) * temp;
    lsStauLStau12(i, 1) = temp2(1);
    lsStauLStau12(i, 2) = temp2(2);
    temp(1) = lsStauRStauLR(i, 1);
    temp(2) = lsStauRStauLR(i, 2);
    temp2 = rot2d(thetatau) * temp;
    lsStauRStau12(i, 1) = temp2(1);
    lsStauRStau12(i, 2) = temp2(2);
 }

  DoubleMatrix lHStauLStau12(lsStauLStau12), lHStauRStau12(lsStauRStau12);
  /// Mix CP-even Higgses up
  for (i=1; i<=2; i++) { /// i is the s1/s2 label
    temp(1) = lsStauLStau12(1, i);
    temp(2) = lsStauLStau12(2, i);
    temp2 = rot2d(alpha) * temp;
    lHStauLStau12(1, i) = temp2(1);
    lHStauLStau12(2, i) = temp2(2);
    temp(1) = lsStauRStau12(1, i);
    temp(2) = lsStauRStau12(2, i);
    temp2 = rot2d(alpha) * temp;
    lHStauRStau12(1, i) = temp2(1);
    lHStauRStau12(2, i) = temp2(2);
  }

  /// Charged Higgs Feynman rules
  DoubleMatrix lChHstauLsnu12(2, 2); /// (H+ G+, L R) basis
  lChHstauLsnu12(2, 1) = -g * displayMwRun() * cos(2.0 * beta) / root2 
    + htau * mtau * cosb;
  lChHstauLsnu12(1, 1) = g * displayMwRun() * sin(2.0 * beta) / root2 
    - htau * mtau * sinb;

  DoubleMatrix lChHstauRsnu12(2, 2); /// (H+ G+, L R) basis
  lChHstauRsnu12(2, 1) = htau * smu * sinb + forLoops.utau * cosb;
  lChHstauRsnu12(1, 1) = htau * smu * cosb - forLoops.utau * sinb;

 for (i=1; i<=4; i++) {
    higgs(1, 1) = higgs(1, 1) + 
      0.5 * (htausq * dnd(i) - sqr(g) * geL * cn(i) 
	     / (2.0 * sqr(costhDrbar))) * a0(higgsm(i), q);
    higgs(2, 2) = higgs(2, 2) + 
      0.5 * (htausq * dnd(i) - sqr(g) * geR * cn(i) 
	     / (2.0 * sqr(costhDrbar))) * a0(higgsm(i), q);
  }
  for (i=3; i<=4; i++) {
    double a0p = a0(higgsc(i - 2), q);
    higgs(1, 1) = higgs(1, 1) + a0p * 
	  (sqr(g) * (geL / (2.0 * sqr(costhDrbar)) + 0.5) * cn(i));
    higgs(2, 2) = higgs(2, 2) + a0p * 
      (htausq * dnu(i) + sqr(g) * geR / (2.0 * sqr(costhDrbar)) * cn(i));
  }
  int j; for(i=1; i<=4; i++)
    for (j=1; j<=2; j++) {
      double b0p = b0(p, higgsm(i), mstau(j), q);
      higgs(1, 1) = higgs(1, 1) + sqr(lHStauLStau12(i, j)) * b0p;
      higgs(1, 2) = higgs(1, 2) + 
	lHStauLStau12(i, j) * lHStauRStau12(i, j) * b0p;
      higgs(2, 2) = higgs(2, 2) + sqr(lHStauRStau12(i, j)) * b0p;
    }
  for(i=1; i<=2; i++) {
    double b0p = b0(p, higgsc(i), msnu(3), q); 
    higgs(1, 1) = higgs(1, 1) + sqr(lChHstauLsnu12(i, 1)) * b0p;
    higgs(1, 2) = higgs(1, 2) + 
      lChHstauLsnu12(i, 1) * lChHstauRsnu12(i, 1) * b0p;
    higgs(2, 2) = higgs(2, 2) + sqr(lChHstauRsnu12(i, 1)) * b0p;
  }

  return higgs;
}

template<class SoftPars>
DoubleMatrix Softsusy<SoftPars>::addStauEweak(double p, double /* mtau */, DoubleMatrix & electroweak) {
  double    mw         = displayMwRun();
  double    mz         = displayMzRun();
  double    sinthDrbar = calcSinthdrbar();
  double    costhDrbar = sqrt(1.0 - sqr(sinthDrbar));
  double    g	       = displayGaugeCoupling(2);
  double    gp         = displayGaugeCoupling(1) * sqrt(0.6);
  double    e          = g * sinthDrbar;
  double    thetat     = forLoops.thetat;
  double    thetab     = forLoops.thetab;
  double    thetatau   = forLoops.thetatau;
  double    ct         = cos(thetat);
  double    st         = sin(thetat);
  double    cb         = cos(thetab);
  double    sb         = sin(thetab);
  double    ctau       = cos(thetatau);
  double    stau       = sin(thetatau);
  double    q          = displayMu();
  DoubleVector msbot(2);
  msbot(1)             = forLoops.md(1, 3);
  msbot(2)             = forLoops.md(2, 3);
  DoubleVector mstop(2);
  mstop(1)             = forLoops.mu(1, 3);
  mstop(2)             = forLoops.mu(2, 3);
  DoubleVector mstau(2);
  mstau(1)             = forLoops.me(1, 3);
  mstau(2)             = forLoops.me(2, 3);
  DoubleVector msnu(3);
  msnu(1)              = forLoops.msnu(1);
  msnu(2)              = forLoops.msnu(2);
  msnu(3)              = forLoops.msnu(3);

  /// LCT: EW contributions
  electroweak(1, 1) = electroweak(1, 1) + 
    4.0 * sqr(g) / sqr(costhDrbar) * sqr(geL) * a0(mz, q) +
    2.0 * sqr(g) * a0(mw, q) + sqr(e) * 
    (sqr(ctau) * ffn(p, mstau(1), 0., q) + 
     sqr(stau) * ffn(p, mstau(2), 0., q));
  electroweak(1, 1) = electroweak(1, 1) +   
    sqr(g) / sqr(costhDrbar) * sqr(geL) * 
    (sqr(ctau) * ffn(p, mstau(1), mz, q) + 
     sqr(stau) * ffn(p, mstau(2), mz, q)) +
    sqr(g) * 0.5 * ffn(p, msnu(3), mw, q);
  electroweak(1, 1) = electroweak(1, 1) +   
    sqr(g) * 0.25 * (sqr(ctau) * a0(mstau(1), q) + 
		     sqr(stau) * a0(mstau(2), q) +
		     2.0 * a0(msnu(3), q));
  electroweak(1, 1) = electroweak(1, 1) + sqr(g) * 
    (-3.0 * 0.25 * (sqr(ct) * a0(mstop(1), q) + sqr(st) * a0(mstop(2), q))
     +3.0 * 0.25 * (sqr(cb) * a0(msbot(1), q) + sqr(sb) * a0(msbot(2), q))
     +      0.25 * (sqr(ctau) * a0(mstau(1), q) + sqr(stau) * a0(mstau(2), q))
     -      0.25 * a0(msnu(3), q)
     -3.0 * 0.25 * (a0(forLoops.mu(1, 2), q) + a0(forLoops.mu(1, 1), q))
     +3.0 * 0.25 * (a0(forLoops.md(1, 2), q) + a0(forLoops.md(1, 1), q))
     +      0.25 * (a0(forLoops.me(1, 2), q) + a0(forLoops.me(1, 1), q))
     -      0.25 * (a0(msnu(2), q) + a0(msnu(1), q)))
    + sqr(gp) * 0.25 * sqr(yeL) * 
    (sqr(ctau) * a0(mstau(1), q) + sqr(stau) * a0(mstau(2), q));
  electroweak(1, 1) = electroweak(1, 1) + sqr(gp) * 0.25 * yeL *
    (+3.0 * yuL * (sqr(ct) * a0(mstop(1), q) + sqr(st) * a0(mstop(2), q))
     +3.0 * ydL * (sqr(cb) * a0(msbot(1), q) + sqr(sb) * a0(msbot(2), q))
     +      yeL * (sqr(ctau) * a0(mstau(1), q) + sqr(stau) * a0(mstau(2), q))
     +      ynuL * a0(msnu(3), q)
     +3.0 * yuL * (a0(forLoops.mu(1, 2), q) + a0(forLoops.mu(1, 1), q))
     +3.0 * ydL * (a0(forLoops.md(1, 2), q) + a0(forLoops.md(1, 1), q))
     +      yeL * (a0(forLoops.me(1, 2), q) + a0(forLoops.me(1, 1), q))
     +      ynuL * (a0(msnu(1), q) + a0(msnu(2), q))
     //
     +3.0 * yuR * (sqr(st) * a0(mstop(1), q) + sqr(ct) * a0(mstop(2), q))
     +3.0 * ydR * (sqr(sb) * a0(msbot(1), q) + sqr(cb) * a0(msbot(2), q))
     +yeR * (sqr(stau) * a0(mstau(1), q) + sqr(ctau) * a0(mstau(2), q))
     +3.0 * yuR * (a0(forLoops.mu(2, 2), q) + a0(forLoops.mu(2, 1), q))
     +3.0 * ydR * (a0(forLoops.md(2, 2), q) + a0(forLoops.md(2, 1), q))
     +      yeR * (a0(forLoops.me(2, 2), q) + a0(forLoops.me(2, 1), q)));


  electroweak(2, 2) = electroweak(2, 2) + 
    4.0 * sqr(g) / sqr(costhDrbar) * sqr(geR) * a0(mz, q) +
    sqr(e) * 
    (sqr(stau) * ffn(p, mstau(1), 0., q) + 
     sqr(ctau) * ffn(p, mstau(2), 0., q));
  electroweak(2, 2) = electroweak(2, 2) + 
    sqr(g) / sqr(costhDrbar) * sqr(geR) * 
    (sqr(stau) * ffn(p, mstau(1), mz, q) + 
     sqr(ctau) * ffn(p, mstau(2), mz, q));
  electroweak(2, 2) = electroweak(2, 2) + 
    sqr(gp) * 0.25 * sqr(yeR) * 
    (sqr(stau) * a0(mstau(1), q) + sqr(ctau) * a0(mstau(2), q));
  electroweak(2, 2) = electroweak(2, 2) + sqr(gp) * 0.25 * yeR *
    (+3.0 * yuL * (sqr(ct) * a0(mstop(1), q) + sqr(st) * a0(mstop(2), q))
     +3.0 * ydL * (sqr(cb) * a0(msbot(1), q) + sqr(sb) * a0(msbot(2), q))
     +      yeL * (sqr(ctau) * a0(mstau(1), q) + sqr(stau) * a0(mstau(2), q))
     +      ynuL * a0(msnu(3), q) 
     +3.0 * yuL * (a0(forLoops.mu(1, 2), q) + a0(forLoops.mu(1, 1), q))
     +3.0 * ydL * (a0(forLoops.md(1, 2), q) + a0(forLoops.md(1, 1), q))
     +      yeL * (a0(forLoops.me(1, 2), q) + a0(forLoops.me(1, 1), q))
     +      ynuL * (a0(msnu(1), q) + a0(msnu(2), q))
     //
     +3.0 * yuR * (sqr(st) * a0(mstop(1), q) + sqr(ct) * a0(mstop(2), q))
     +3.0 * ydR * (sqr(sb) * a0(msbot(1), q) + sqr(cb) * a0(msbot(2), q))
     +      yeR * (sqr(stau) * a0(mstau(1), q) + sqr(ctau) * a0(mstau(2), q))
     +3.0 * yuR * (a0(forLoops.mu(2, 2), q) + a0(forLoops.mu(2, 1), q))
     +3.0 * ydR * (a0(forLoops.md(2, 2), q) + a0(forLoops.md(2, 1), q))
     +      yeR * (a0(forLoops.me(2, 2), q) + a0(forLoops.me(2, 1), q)));

  electroweak(1, 2) = electroweak(1, 2) + 
    sqr(gp) * stau * ctau * 0.25 * yeL * yeR * 
    (a0(mstau(1), q) - a0(mstau(2), q)) +
    sqr(e) * stau * ctau * 
    (ffn(p, mstau(1), 0., q) - ffn(p, mstau(2), 0., q)) -
    sqr(g) / sqr(costhDrbar) * geL * geR * stau * ctau * 
    (ffn(p, mstau(1), mz, q) - ffn(p, mstau(2), mz, q));

  return electroweak;
}

template<class SoftPars>
void Softsusy<SoftPars>::addStauGaugino(double p, double mtau, DoubleMatrix & chargino, DoubleMatrix & neutralino) {
  double    g	    = displayGaugeCoupling(2);
  double    gp      = displayGaugeCoupling(1) * sqrt(0.6); 
  double    q       = displayMu();
  double    htau    = forLoops.htau;

  ComplexMatrix n(forLoops.nBpmz);
  DoubleVector mneut(forLoops.mnBpmz);
  ComplexMatrix u(forLoops.uBpmz), v(forLoops.vBpmz); 
  DoubleVector mch(forLoops.mchBpmz); 

  const int rank = mneut.displayEnd();

  /// Neutralino Feynman rules
  DoubleVector aPsi0TauStaur(rank), bPsi0TauStaur(rank), aPsi0TauStaul(rank),
    bPsi0TauStaul(rank); 
  aPsi0TauStaur(1) = gp * root2;
  bPsi0TauStaul(1) = -gp / root2;
  bPsi0TauStaul(2) = -g / root2;
  aPsi0TauStaul(3) = htau;
  bPsi0TauStaur(3) = htau;

  ComplexVector aChi0TauStaul(rank), bChi0TauStaul(rank), aChi0TauStaur(rank),
    bChi0TauStaur(rank);

  aChi0TauStaul = n.complexConjugate() * aPsi0TauStaul;
  bChi0TauStaul = n * bPsi0TauStaul;
  aChi0TauStaur = n.complexConjugate() * aPsi0TauStaur;
  bChi0TauStaur = n * bPsi0TauStaur;

  DoubleVector gChi0TauStauLL(rank), fChi0TauStauLL(rank);
  DoubleVector gChi0TauStauLR(rank), fChi0TauStauLR(rank);
  DoubleVector gChi0TauStauRR(rank), fChi0TauStauRR(rank);
  int i;
  for (i=1; i<=rank; i++) {
    fChi0TauStauLL(i) = (aChi0TauStaul(i) * aChi0TauStaul(i).conj() + 
      bChi0TauStaul(i) * bChi0TauStaul(i).conj()).real();
    gChi0TauStauLL(i) = (bChi0TauStaul(i).conj() * aChi0TauStaul(i) + 
      bChi0TauStaul(i) * aChi0TauStaul(i).conj()).real();
    fChi0TauStauRR(i) = (aChi0TauStaur(i) * aChi0TauStaur(i).conj() + 
      bChi0TauStaur(i) * bChi0TauStaur(i).conj()).real();
    gChi0TauStauRR(i) = (bChi0TauStaur(i).conj() * aChi0TauStaur(i) + 
      bChi0TauStaur(i) * aChi0TauStaur(i).conj()).real();
    fChi0TauStauLR(i) = (aChi0TauStaur(i) * aChi0TauStaul(i).conj() + 
      bChi0TauStaur(i) * bChi0TauStaul(i).conj()).real();
    gChi0TauStauLR(i) = (bChi0TauStaul(i).conj() * aChi0TauStaur(i) + 
      bChi0TauStaur(i) * aChi0TauStaul(i).conj()).real();
  }

  /// Chargino Feynman Rules
  DoubleVector bPsicNuStaul(2), bPsicNuStaur(2);

  bPsicNuStaul(1) = g;
  bPsicNuStaur(2) = -htau;
  
  ComplexVector bChicNuStaul(2), bChicNuStaur(2);
  ComplexMatrix bChNuStau(2, 2);

  bChicNuStaul = u * bPsicNuStaul;
  bChicNuStaur = u * bPsicNuStaur;

  DoubleVector fChNuStauLL(2) ;
  DoubleVector fChNuStauLR(2); 
  DoubleVector fChNuStauRR(2); 
  for (i=1; i<=2; i++) {
    fChNuStauLL(i) = (bChicNuStaul(i).conj() * bChicNuStaul(i)).real();
    fChNuStauLR(i) = (bChicNuStaul(i).conj() * bChicNuStaur(i)).real();
    fChNuStauRR(i) = (bChicNuStaur(i).conj() * bChicNuStaur(i)).real();
  }

  for (i=1; i<=2; i++) {
    double one = gfn(p, mch(i), 0., q);
    chargino(1, 1) = chargino(1, 1) + fChNuStauLL(i) * one;
    chargino(1, 2) = chargino(1, 2) + fChNuStauLR(i) * one;
    chargino(2, 2) = chargino(2, 2) + fChNuStauRR(i) * one;
  }

  for (i=1; i<=rank; i++) {
    double one = gfn(p, mneut(i), mtau, q);
    double two = 2.0 * mneut(i) * mtau * b0(p, mneut(i), mtau, q);
    neutralino(1, 1) = neutralino(1, 1) + fChi0TauStauLL(i) * one
       - gChi0TauStauLL(i) * two;
    neutralino(2, 2) = neutralino(2, 2) + fChi0TauStauRR(i) * one
       - gChi0TauStauRR(i) * two;
    neutralino(1, 2) = neutralino(1, 2) + fChi0TauStauLR(i) * one
       - gChi0TauStauLR(i) * two;
  }

}

template<class SoftPars>
void Softsusy<SoftPars>::addStauCorrection(double p, DoubleMatrix & mass, 
				     double mtau) {

/// No point adding radiative corrections to tachyonic particles
  if (mass(1, 1) < 0.0 || mass(2, 2) < 0.0) { 
    flagTachyon(stau);
    if (mass(1, 1) < 0.) mass(1, 1) = EPSTOL;
    else mass(2, 2) = EPSTOL;
    return;
  }

  /// one-loop correction matrix
  DoubleMatrix piSq(2, 2); 
	
  DoubleMatrix stop(2, 2), sbottom(2, 2), 
    higgs(2, 2), electroweak(2, 2), chargino(2, 2), neutralino(2, 2);

  /// LCT: Sfermion contribution
  addStauSfermion(p, mtau, stop, sbottom);
  /// LCT: Higgs contributions
  addStauHiggs(p, mtau, higgs);
  /// LCT: Electroweak contribution
  addStauEweak(p, mtau, electroweak);
  /// LCT: Chargino and neutralino contributions
  addStauGaugino(p, mtau, chargino, neutralino);

  piSq = 1.0 / (16.0 * sqr(PI)) * 
    (stop + sbottom + higgs + chargino + neutralino + electroweak);

  piSq(2, 1) = piSq(1, 2);

  mass = mass - piSq;	
}

template<class SoftPars>
DoubleMatrix Softsusy<SoftPars>::addSdownQCD(double p1, double p2, int family, DoubleMatrix & strong) {
  double mg   = forLoops.mGluino;
  double q    = displayMu();
  double g3sq = sqr(displayGaugeCoupling(3));
  DoubleVector msd(2);
  msd(1)      = forLoops.md(1, family);
  msd(2)      = forLoops.md(2, family);
  double a0t1 = a0(msd(1), q), a0t2 = a0(msd(2), q);
  double ft1  = ffn(p1, msd(1), 0.0, q), ft2 = ffn(p2, msd(2), 0.0, q);
  double ggt1 = gfn(p1, mg, 0., q), ggt2 = gfn(p2, mg, 0., q);

  strong(1, 1) = 4.0 * g3sq / 3.0 * (2.0 * ggt1 + ft1 + a0t1);
  strong(2, 2) = 4.0 * g3sq / 3.0 * (2.0 * ggt2 + ft2 + a0t2);

  return strong;
}

template<class SoftPars>
DoubleMatrix Softsusy<SoftPars>::addSdownHiggs(double p1, double p2, int family, DoubleMatrix & higgs) {
  double mz         = displayMzRun();
  double sinthDrbar = calcSinthdrbar();
  double costhDrbar = sqrt(1.0 - sqr(sinthDrbar));
  double alpha      = forLoops.thetaH;
  double g          = displayGaugeCoupling(2);
  double beta       = atan(displayTanb());
  double sinb       = sin(beta), cosb = cos(beta);
  double q          = displayMu();
  DoubleVector msd(2);
  msd(1)            = forLoops.md(1, family);
  msd(2)            = forLoops.md(2, family);
  DoubleVector msup(2);
  msup(1)           = forLoops.mu(1, family);
  msup(2)           = forLoops.mu(2, family);

  /// define Higgs vector in 't-Hooft Feynman gauge, and couplings:
  DoubleVector higgsm(4), higgsc(2);
  DoubleVector dnu(4), dnd(4), cn(4);
  assignHiggsSfermions(higgsm, higgsc, dnu, dnd, cn, beta);

  DoubleMatrix lsSbotLSbotLR(4, 2);
  DoubleMatrix lsSbotRSbotLR(4, 2);
  /// Order (s1 s2 G A, L R)
  lsSbotLSbotLR(1, 1) = g * mz * gdL * cosb / costhDrbar;
  lsSbotLSbotLR(2, 1) = -g * mz * gdL * sinb / costhDrbar;

  lsSbotRSbotLR(1, 1) = lsSbotLSbotLR(1, 2);
  lsSbotRSbotLR(1, 2) = g * mz * gdR * cosb / costhDrbar;
  lsSbotRSbotLR(2, 1) = lsSbotLSbotLR(2, 2);
  lsSbotRSbotLR(2, 2) = -g * mz * gdR * sinb / costhDrbar;
  lsSbotRSbotLR(3, 1) = -lsSbotLSbotLR(3, 2);
  lsSbotRSbotLR(4, 1) = -lsSbotLSbotLR(4, 2);

  /// Mix sbots up
  int i;
  DoubleVector temp(2), temp2(2);
  DoubleMatrix lHSbotLSbot12(lsSbotLSbotLR), lHSbotRSbot12(lsSbotRSbotLR);
  /// Mix CP-even Higgses up
  for (i=1; i<=2; i++) { /// i is the L/R label
    temp(1) = lsSbotLSbotLR(1, i);
    temp(2) = lsSbotLSbotLR(2, i);
    temp2 = rot2d(alpha) * temp;
    lHSbotLSbot12(1, i) = temp2(1);
    lHSbotLSbot12(2, i) = temp2(2);
    temp(1) = lsSbotRSbotLR(1, i);
    temp(2) = lsSbotRSbotLR(2, i);
    temp2 = rot2d(alpha) * temp;
    lHSbotRSbot12(1, i) = temp2(1);
    lHSbotRSbot12(2, i) = temp2(2);
  }

  /// Charged Higgs Feynman rules
  DoubleMatrix lChHsbotLstopLR(2, 2); /// (H+ G+, L R) basis
  lChHsbotLstopLR(2, 1) = -g * displayMwRun() * cos(2.0 * beta) / root2;
  lChHsbotLstopLR(1, 1) = g * displayMwRun() * sin(2.0 * beta) / root2;
  DoubleMatrix lChHsbotLstop12(2, 2);
  temp(1) = lChHsbotLstopLR(1, 1);
  temp(2) = lChHsbotLstopLR(1, 2);
  temp2 = temp;
  lChHsbotLstop12(1, 1) = temp2(1);
  lChHsbotLstop12(1, 2) = temp2(2);
  temp(1) = lChHsbotLstopLR(2, 1);
  temp(2) = lChHsbotLstopLR(2, 2);
  temp2 = temp;
  lChHsbotLstop12(2, 1) = temp2(1);
  lChHsbotLstop12(2, 2) = temp2(2);
  /// (none for sdownL since they are all Yukawa suppressed)

  for (i=1; i<=4; i++) {
    higgs(1, 1) += 
      0.5 * (- sqr(g) * gdL * cn(i) 
	     / (2.0 * sqr(costhDrbar))) * a0(higgsm(i), q);
    higgs(2, 2) += 
      0.5 * (- sqr(g) * gdR * cn(i) 
	     / (2.0 * sqr(costhDrbar))) * a0(higgsm(i), q);
  }
  for (i=3; i<=4; i++) {
    double a0p = a0(higgsc(i - 2), q);
    higgs(1, 1) += a0p * 
	  (sqr(g) * (gdL / (2.0 * sqr(costhDrbar)) + 0.5) * cn(i));
    higgs(2, 2) += a0p *
      (sqr(g) * gdR / (2.0 * sqr(costhDrbar)) * cn(i));
  }
  int j; for(i=1; i<=4; i++)
    for (j=1; j<=2; j++) {
      double b0p = b0(p1, higgsm(i), msd(j), q);
      higgs(1, 1) += sqr(lHSbotLSbot12(i, j)) * b0p;
      b0p = b0(p2, higgsm(i), msd(j), q);
      higgs(2, 2) += sqr(lHSbotRSbot12(i, j)) * b0p;
    }
  for(i=1; i<=2; i++)
    for (j=1; j<=2; j++) {
      double b0p = b0(p1, higgsc(i), msup(j), q); 
      higgs(1, 1) += sqr(lChHsbotLstop12(i, j)) * b0p;
    }

  return higgs; 
}

template<class SoftPars>
DoubleMatrix Softsusy<SoftPars>::addSdownEweak(double p1, double p2, int family, DoubleMatrix & electroweak) {
  double mw         = displayMwRun();
  double mz         = displayMzRun();
  double sinthDrbar = calcSinthdrbar();
  double costhDrbar = sqrt(1.0 - sqr(sinthDrbar));
  double g          = displayGaugeCoupling(2);
  double gp         = displayGaugeCoupling(1) * sqrt(0.6);
  double e          = g * sinthDrbar;
  double thetat     = forLoops.thetat;
  double thetab     = forLoops.thetab;
  double thetatau   = forLoops.thetatau;
  double ct         = cos(thetat);
  double st         = sin(thetat);
  double cb         = cos(thetab);
  double sb         = sin(thetab);
  double ctau       = cos(thetatau);
  double stau       = sin(thetatau);
  double q          = displayMu();
  DoubleVector msbot(2);
  msbot(1)          = forLoops.md(1, 3);
  msbot(2)          = forLoops.md(2, 3);
  DoubleVector msd(2);
  msd(1)            = forLoops.md(1, family);
  msd(2)            = forLoops.md(2, family);
  DoubleVector msup(2);
  msup(1)           = forLoops.mu(1, family);
  msup(2)           = forLoops.mu(2, family);
  DoubleVector mstop(2);
  mstop(1)          = forLoops.mu(1, 3);
  mstop(2)          = forLoops.mu(2, 3);
  DoubleVector mstau(2);
  mstau(1)          = forLoops.me(1, 3);
  mstau(2)          = forLoops.me(2, 3);
  DoubleVector msnu(3);
  msnu(1)           = forLoops.msnu(1);
  msnu(2)           = forLoops.msnu(2);
  msnu(3)           = forLoops.msnu(3);

  electroweak(1, 1) = electroweak(1, 1) + 
    4.0 * sqr(g) / sqr(costhDrbar) * sqr(gdL) * a0(mz, q) +
    2.0 * sqr(g) * a0(mw, q) + sqr(e / 3.0) * 
    ffn(p1, msd(1), 0., q);
  electroweak(1, 1) = electroweak(1, 1) +   
    sqr(g) / sqr(costhDrbar) * sqr(gdL) * ffn(p1, msd(1), mz, q) +
    sqr(g) * 0.5 * ffn(p1, msup(1), mw, q);
  electroweak(1, 1) = electroweak(1, 1) +   
    sqr(g) * 0.25 * (a0(msd(1), q) + 2.0 * (a0(msup(1), q)));
  electroweak(1, 1) = electroweak(1, 1) + sqr(g) * 
    (-3.0 * 0.25 * (sqr(ct) * a0(mstop(1), q) + sqr(st) * a0(mstop(2), q))
     +3.0 * 0.25 * (sqr(cb) * a0(msbot(1), q) + sqr(sb) * a0(msbot(2), q))
     +      0.25 * (sqr(ctau) * a0(mstau(1), q) + sqr(stau) * a0(mstau(2), q))
     -      0.25 * a0(msnu(3), q)
     -3.0 * 0.25 * (a0(forLoops.mu(1, 2), q) + a0(forLoops.mu(1, 1), q))
     +3.0 * 0.25 * (a0(forLoops.md(1, 2), q) + a0(forLoops.md(1, 1), q))
     +      0.25 * (a0(forLoops.me(1, 2), q) + a0(forLoops.me(1, 1), q))
     -      0.25 * (a0(msnu(2), q) + a0(msnu(1), q)))
    + sqr(gp) * 0.25 * sqr(ydL) * a0(msd(1), q);
  electroweak(1, 1) = electroweak(1, 1) + sqr(gp) * 0.25 * ydL *
    (+3.0 * yuL * (sqr(ct) * a0(mstop(1), q) + sqr(st) * a0(mstop(2), q))
     +3.0 * ydL * (sqr(cb) * a0(msbot(1), q) + sqr(sb) * a0(msbot(2), q))
     +      yeL * (sqr(ctau) * a0(mstau(1), q) + sqr(stau) * a0(mstau(2), q))
     +      ynuL * a0(msnu(3), q)
     +3.0 * yuL * (a0(forLoops.mu(1, 2), q) + a0(forLoops.mu(1, 1), q))
     +3.0 * ydL * (a0(forLoops.md(1, 2), q) + a0(forLoops.md(1, 1), q))
     +      yeL * (a0(forLoops.me(1, 2), q) + a0(forLoops.me(1, 1), q))
     +      ynuL * (a0(msnu(1), q) + a0(msnu(2), q))
     //
     +3.0 * yuR * (sqr(st) * a0(mstop(1), q) + sqr(ct) * a0(mstop(2), q))
     +3.0 * ydR * (sqr(sb) * a0(msbot(1), q) + sqr(cb) * a0(msbot(2), q))
     +yeR * (sqr(stau) * a0(mstau(1), q) + sqr(ctau) * a0(mstau(2), q))
     +3.0 * yuR * (a0(forLoops.mu(2, 2), q) + a0(forLoops.mu(2, 1), q))
     +3.0 * ydR * (a0(forLoops.md(2, 2), q) + a0(forLoops.md(2, 1), q))
     +      yeR * (a0(forLoops.me(2, 2), q) + a0(forLoops.me(2, 1), q)));

  electroweak(2, 2) = electroweak(2, 2) + 
    4.0 * sqr(g) / sqr(costhDrbar) * sqr(gdR) * a0(mz, q) +
    sqr(e / 3.0) * 
    (ffn(p2, msd(2), 0., q));
  electroweak(2, 2) = electroweak(2, 2) + 
    sqr(g) / sqr(costhDrbar) * sqr(gdR) * 
    (ffn(p2, msd(2), mz, q));
  electroweak(2, 2) = electroweak(2, 2) + 
    sqr(gp) * 0.25 * sqr(ydR) * 
    (a0(msd(2), q));
  electroweak(2, 2) = electroweak(2, 2) + sqr(gp) * 0.25 * ydR *
    (+3.0 * yuL * (sqr(ct) * a0(mstop(1), q) + sqr(st) * a0(mstop(2), q))
     +3.0 * ydL * (sqr(cb) * a0(msbot(1), q) + sqr(sb) * a0(msbot(2), q))
     +      yeL * (sqr(ctau) * a0(mstau(1), q) + sqr(stau) * a0(mstau(2), q))
     +      ynuL * a0(msnu(3), q) 
     +3.0 * yuL * (a0(forLoops.mu(1, 2), q) + a0(forLoops.mu(1, 1), q))
     +3.0 * ydL * (a0(forLoops.md(1, 2), q) + a0(forLoops.md(1, 1), q))
     +      yeL * (a0(forLoops.me(1, 2), q) + a0(forLoops.me(1, 1), q))
     +      ynuL * (a0(msnu(1), q) + a0(msnu(2), q))
     //
     +3.0 * yuR * (sqr(st) * a0(mstop(1), q) + sqr(ct) * a0(mstop(2), q))
     +3.0 * ydR * (sqr(sb) * a0(msbot(1), q) + sqr(cb) * a0(msbot(2), q))
     +      yeR * (sqr(stau) * a0(mstau(1), q) + sqr(ctau) * a0(mstau(2), q))
     +3.0 * yuR * (a0(forLoops.mu(2, 2), q) + a0(forLoops.mu(2, 1), q))
     +3.0 * ydR * (a0(forLoops.md(2, 2), q) + a0(forLoops.md(2, 1), q))
     +      yeR * (a0(forLoops.me(2, 2), q) + a0(forLoops.me(2, 1), q)));

  return electroweak;
}

template<class SoftPars>
DoubleMatrix Softsusy<SoftPars>::addSdownNeutralino(double p1, double p2, int /* family */, DoubleMatrix & neutralino) {
  double g  = displayGaugeCoupling(2);
  double gp = displayGaugeCoupling(1) * sqrt(0.6);
  double q  = displayMu();

  ComplexMatrix n(forLoops.nBpmz);
  DoubleVector mneut(forLoops.mnBpmz);

  int rank = mneut.displayEnd();

  /// Neutralino Feynman rules
  DoubleVector aPsi0BSbotr(rank), bPsi0BSbotr(rank), 
    aPsi0BSbotl(rank), bPsi0BSbotl(rank); 
  aPsi0BSbotr(1) = gp * ydR / root2;
  bPsi0BSbotl(1) = gp * ydL / root2;
  bPsi0BSbotl(2) = -0.5 * g * root2;

  ComplexVector aChi0BSbotl(rank), bChi0BSbotl(rank), 
    aChi0BSbotr(rank), bChi0BSbotr(rank);
  aChi0BSbotl = n.complexConjugate() * aPsi0BSbotl;
  bChi0BSbotl = n * bPsi0BSbotl;
  aChi0BSbotr = n.complexConjugate() * aPsi0BSbotr;
  bChi0BSbotr = n * bPsi0BSbotr;

  DoubleVector gChi0BotSbotLL(rank), fChi0BotSbotLL(rank);
  DoubleVector gChi0BotSbotLR(rank), fChi0BotSbotLR(rank);
  DoubleVector gChi0BotSbotRR(rank), fChi0BotSbotRR(rank);

  for (int i=1; i<=rank; i++) {
    fChi0BotSbotLL(i) = (aChi0BSbotl(i) * aChi0BSbotl(i).conj() + 
      bChi0BSbotl(i) * bChi0BSbotl(i).conj()).real();
    gChi0BotSbotLL(i) = (bChi0BSbotl(i).conj() * aChi0BSbotl(i) + 
      bChi0BSbotl(i) * aChi0BSbotl(i).conj()).real();
    fChi0BotSbotRR(i) = (aChi0BSbotr(i) * aChi0BSbotr(i).conj() + 
      bChi0BSbotr(i) * bChi0BSbotr(i).conj()).real();
    gChi0BotSbotRR(i) = (bChi0BSbotr(i).conj() * aChi0BSbotr(i) + 
      bChi0BSbotr(i) * aChi0BSbotr(i).conj()).real();
    fChi0BotSbotLR(i) = (aChi0BSbotr(i) * aChi0BSbotl(i).conj() + 
      bChi0BSbotr(i) * bChi0BSbotl(i).conj()).real();
    gChi0BotSbotLR(i) = (bChi0BSbotl(i).conj() * aChi0BSbotr(i) + 
      bChi0BSbotr(i) * aChi0BSbotl(i).conj()).real();
  }

  for (int i=1; i<=rank; i++) {
    double one = gfn(p1, mneut(i), 0., q);
    neutralino(1, 1) = neutralino(1, 1) +
      fChi0BotSbotLL(i) * one;
    one = gfn(p2, mneut(i), 0., q);
    neutralino(2, 2) = neutralino(2, 2) +
      fChi0BotSbotRR(i) * one;
  }

  return neutralino;
}

template<class SoftPars>
DoubleMatrix Softsusy<SoftPars>::addSdownChargino(double p1, double p2, int /* family */, DoubleMatrix & chargino) {
  double g = displayGaugeCoupling(2);
  double q = displayMu();

  ComplexMatrix u(forLoops.uBpmz), v(forLoops.vBpmz); 
  DoubleVector mch(forLoops.mchBpmz); 

  /// Chargino Feynman Rules
  DoubleVector bPsicTSbotl(2), bPsicTSbotr(2), aPsicTSbotl(2), aPsicTSbotr(2); 
  bPsicTSbotl(1) = g;
  
  ComplexVector aChicTSbotr(2), aChicTSbotl(2), bChicTSbotl(2), bChicTSbotr(2);
  aChicTSbotl = v.complexConjugate() * aPsicTSbotl;
  bChicTSbotl = u * bPsicTSbotl;
  aChicTSbotr = v.complexConjugate() * aPsicTSbotr;
  bChicTSbotr = u * bPsicTSbotr;

  DoubleVector fChTSbotLL(2), gChTSbotLL(2);
  DoubleVector fChTSbotLR(2), gChTSbotLR(2); 
  DoubleVector fChTSbotRR(2), gChTSbotRR(2); 

  for (int i=1; i<=2; i++) {
    fChTSbotLL(i) = (aChicTSbotl(i).conj() * aChicTSbotl(i) +
		      bChicTSbotl(i).conj() * bChicTSbotl(i)).real();
    gChTSbotLL(i) = (bChicTSbotl(i).conj() * aChicTSbotl(i) +
		      aChicTSbotl(i).conj() * bChicTSbotl(i)).real();
    fChTSbotLR(i) = (aChicTSbotl(i).conj() * aChicTSbotr(i) +
		      bChicTSbotl(i).conj() * bChicTSbotr(i)).real();
    gChTSbotLR(i) = (bChicTSbotl(i).conj() * aChicTSbotr(i) +
		      aChicTSbotl(i).conj() * bChicTSbotr(i)).real();
    fChTSbotRR(i) = (aChicTSbotr(i).conj() * aChicTSbotr(i) +
		      bChicTSbotr(i).conj() * bChicTSbotr(i)).real();
    gChTSbotRR(i) = (bChicTSbotr(i).conj() * aChicTSbotr(i) +
		      aChicTSbotr(i).conj() * bChicTSbotr(i)).real();
  }

  for (int i=1; i<=2; i++) {
    double one = gfn(p1, mch(i), 0., q);
    chargino(1, 1) = chargino(1, 1) + fChTSbotLL(i) * one;
    one = gfn(p2, mch(i), 0., q);
    chargino(2, 2) = chargino(2, 2) + fChTSbotRR(i) * one;
  }

  return chargino;
}

template<class SoftPars>
void Softsusy<SoftPars>::addSdownCorrection(DoubleMatrix & mass, int family) {

  /// No point adding radiative corrections to tachyonic particles
  if (mass(1, 1) < 0.0 || mass(2, 2) < 0.0) { 
    if (family == 1) flagTachyon(sdown); 
    else if (family == 2) flagTachyon(sstrange);
    return;
  }

  DoubleMatrix piSq(2, 2); /// Self-energy matrix
  DoubleVector msd(2);
  msd(1) = forLoops.md(1, family);
  msd(2) = forLoops.md(2, family);
  double    p1 = msd(1), p2 = msd(2);
  ///  p1 = p2 = sqrt(msd(1) * msd(2)); 

  DoubleMatrix strong(2, 2), stop(2, 2), sbottom(2, 2), 
    higgs(2, 2), chargino(2, 2), neutralino(2, 2), electroweak(2, 2);

  /// LCT: QCD contribution
  addSdownQCD(p1, p2, family, strong);
  /// LCT: Higgs contribution
  addSdownHiggs(p1, p2, family, higgs);
  /// LCT: Electroweak contribution
  addSdownEweak(p1, p2, family, electroweak);
  /// LCT: Chargino contribution
  addSdownChargino(p1, p2, family, chargino);
  /// LCT: Neutralino contribution
  addSdownNeutralino(p1, p2, family, neutralino);
 
  piSq = 1.0 / (16.0 * sqr(PI)) * 
    (strong + higgs + chargino + neutralino + electroweak);

  mass = mass - piSq;	
}

template<class SoftPars>
DoubleMatrix Softsusy<SoftPars>::addSbotQCD(double p, double /* mt */, DoubleMatrix & strong) {
  double    thetab  = forLoops.thetab;
  double    cb      = cos(thetab);
  double    sb      = sin(thetab);
  double    mg      = forLoops.mGluino;
  double    q       = displayMu();
  double    g3sq    = sqr(displayGaugeCoupling(3));
  double    mb      = forLoops.mb;
  DoubleVector msbot(2);
  msbot(1)          = forLoops.md(1, 3);
  msbot(2)          = forLoops.md(2, 3);
  double    a0t1    = a0(msbot(1), q), a0t2 = a0(msbot(2), q);
  double    ft1     = ffn(p, msbot(1), 0.0, q), ft2 = ffn(p, msbot(2), 0.0, q);
  double    ggt     = gfn(p, mg, mb, q);

  strong(1, 1) = 4.0 * g3sq / 3.0 *
    (2.0 * ggt + sqr(cb) * (ft1 + a0t1) + sqr(sb) * (ft2 + a0t2));
  strong(2, 2) = 4.0 * g3sq / 3.0 *
    (2.0 * ggt + sqr(sb) * (ft1 + a0t1) + sqr(cb) * (ft2 + a0t2));
  strong(1, 2) = 4.0 * g3sq / 3.0 *
    (4.0 * mg * mb * b0(p, mg, mb, q) + sb * cb * (ft1 - a0t1 - ft2 + a0t2));

  return strong;
}

template<class SoftPars>
void Softsusy<SoftPars>::addSbotSfermion(double /* p */, double /* mt */, DoubleMatrix & stop, DoubleMatrix & sbottom) {
  double    thetat   = forLoops.thetat;
  double    thetab   = forLoops.thetab;
  double    thetatau = forLoops.thetatau;
  double    ct       = cos(thetat);
  double    st       = sin(thetat);
  double    cb       = cos(thetab);
  double    sb       = sin(thetab);
  double    ctau     = cos(thetatau);
  double    stau     = sin(thetatau);
  double    q        = displayMu();
  double    ht       = forLoops.ht, htsq = sqr(ht);
  double    hb       = forLoops.hb, hbsq = sqr(hb);
  double    htau     = forLoops.htau;
  DoubleVector msbot(2);
  msbot(1)           = forLoops.md(1, 3);
  msbot(2)           = forLoops.md(2, 3);
  DoubleVector mstop(2);
  mstop(1)           = forLoops.mu(1, 3);
  mstop(2)           = forLoops.mu(2, 3);
  DoubleVector mstau(2);
  mstau(1)           = forLoops.me(1, 3);
  mstau(2)           = forLoops.me(2, 3);

  double    a0t1     = a0(msbot(1), q), a0t2 = a0(msbot(2), q);

  stop(1, 1) = hbsq * (sqr(sb) * a0t1 + sqr(cb) * a0t2);
  stop(2, 2) = hbsq * (sqr(cb) * a0t1 + sqr(sb) * a0t2);
  stop(1, 2) = hbsq * cb * sb * 3.0 * (a0t1 - a0t2)
     + hb * htau * ctau * stau * (a0(mstau(1), q) - a0(mstau(2), q));

  sbottom(1, 1) = 
    htsq * (sqr(st) * a0(mstop(1), q) + sqr(ct) * a0(mstop(2), q));
  sbottom(2, 2) = 
    hbsq * (sqr(ct) * a0(mstop(1), q) + sqr(st) * a0(mstop(2), q));
}

template<class SoftPars>
DoubleMatrix Softsusy<SoftPars>::addSbotHiggs(double p, double mt, DoubleMatrix & higgs) {
  double    mz         = displayMzRun();
  double    sinthDrbar = calcSinthdrbar();
  double    costhDrbar = sqrt(1.0 - sqr(sinthDrbar));
  double    alpha      = forLoops.thetaH;
  double    g	       = displayGaugeCoupling(2);
  double    beta       = atan(displayTanb());
  double    sinb       = sin(beta), cosb = cos(beta);
  double    thetat     = forLoops.thetat;
  double    thetab     = forLoops.thetab;
  double    smu        = -displaySusyMu();
  double    q          = displayMu();
  double    ht         = forLoops.ht, htsq = sqr(ht);
  double    hb         = forLoops.hb, hbsq = sqr(hb);
  double    mb         = forLoops.mb;
  DoubleVector msbot(2);
  msbot(1)             = forLoops.md(1, 3);
  msbot(2)             = forLoops.md(2, 3);
  DoubleVector mstop(2);
  mstop(1)             = forLoops.mu(1, 3);
  mstop(2)             = forLoops.mu(2, 3); 

  /// define Higgs vector in 't-Hooft Feynman gauge, and couplings:
  DoubleVector higgsm(4), higgsc(2);
  DoubleVector dnu(4), dnd(4), cn(4);
  assignHiggsSfermions(higgsm, higgsc, dnu, dnd, cn, beta);

  DoubleMatrix lsSbotLSbotLR(4, 2), lsSbotLSbot12(4, 2);
  DoubleMatrix lsSbotRSbotLR(4, 2), lsSbotRSbot12(4, 2);
  /// Order (s1 s2 G A, L R)
  lsSbotLSbotLR(1, 1) = g * mz * gdL * cosb / costhDrbar + 
    root2 * hb * mb;
  lsSbotLSbotLR(1, 2) = forLoops.ub / root2;
  lsSbotLSbotLR(3, 2) = -1.0 / root2 * 
    (smu * sinb * hb + forLoops.ub * cosb);
  lsSbotLSbotLR(4, 2) = -1.0 / root2 * 
    (smu * cosb * hb - forLoops.ub * sinb);
  lsSbotLSbotLR(2, 1) = -g * mz * gdL * sinb / costhDrbar;
  lsSbotLSbotLR(2, 2) = hb * smu / root2;

  lsSbotRSbotLR(1, 1) = lsSbotLSbotLR(1, 2);
  lsSbotRSbotLR(1, 2) = g * mz * gdR * cosb / costhDrbar +
    root2 * hb * mb;
  lsSbotRSbotLR(2, 1) = lsSbotLSbotLR(2, 2);
  lsSbotRSbotLR(2, 2) = -g * mz * gdR * sinb / costhDrbar;
  lsSbotRSbotLR(3, 1) = -lsSbotLSbotLR(3, 2);
  lsSbotRSbotLR(4, 1) = -lsSbotLSbotLR(4, 2);

  /// Mix sbots up
  int i;
  DoubleVector temp(2), temp2(2);
  for (i=1; i<=4; i++) {
    temp(1) = lsSbotLSbotLR(i, 1);
    temp(2) = lsSbotLSbotLR(i, 2);
    temp2 = rot2d(thetab) * temp;
    lsSbotLSbot12(i, 1) = temp2(1);
    lsSbotLSbot12(i, 2) = temp2(2);
    temp(1) = lsSbotRSbotLR(i, 1);
    temp(2) = lsSbotRSbotLR(i, 2);
    temp2 = rot2d(thetab) * temp;
    lsSbotRSbot12(i, 1) = temp2(1);
    lsSbotRSbot12(i, 2) = temp2(2);
 }

  DoubleMatrix lHSbotLSbot12(lsSbotLSbot12), lHSbotRSbot12(lsSbotRSbot12);
  /// Mix CP-even Higgses up
  for (i=1; i<=2; i++) { /// i is the s1/s2 label
    temp(1) = lsSbotLSbot12(1, i);
    temp(2) = lsSbotLSbot12(2, i);
    temp2 = rot2d(alpha) * temp;
    lHSbotLSbot12(1, i) = temp2(1);
    lHSbotLSbot12(2, i) = temp2(2);
    temp(1) = lsSbotRSbot12(1, i);
    temp(2) = lsSbotRSbot12(2, i);
    temp2 = rot2d(alpha) * temp;
    lHSbotRSbot12(1, i) = temp2(1);
    lHSbotRSbot12(2, i) = temp2(2);
  }


  /// Charged Higgs Feynman rules (H+ G+, L R) basis
  DoubleMatrix lChHsbotLstopLR(2, 2); 
  lChHsbotLstopLR(2, 1) = -g * displayMwRun() * cos(2.0 * beta) / root2 
    - ht * mt * sinb + hb * mb * cosb;
  lChHsbotLstopLR(2, 2) = -ht * smu * cosb - forLoops.ut * sinb;
  lChHsbotLstopLR(1, 1) = g * displayMwRun() * sin(2.0 * beta) / root2 
    - ht * mt * cosb - hb * mb * sinb;
  lChHsbotLstopLR(1, 2) = ht * smu * sinb - forLoops.ut * cosb;

  DoubleMatrix lChHsbotLstop12(2, 2);
  temp(1) = lChHsbotLstopLR(1, 1);
  temp(2) = lChHsbotLstopLR(1, 2);
  temp2 = rot2d(thetat) * temp;
  lChHsbotLstop12(1, 1) = temp2(1);
  lChHsbotLstop12(1, 2) = temp2(2);
  temp(1) = lChHsbotLstopLR(2, 1);
  temp(2) = lChHsbotLstopLR(2, 2);
  temp2 = rot2d(thetat) * temp;
  lChHsbotLstop12(2, 1) = temp2(1);
  lChHsbotLstop12(2, 2) = temp2(2);

  DoubleMatrix lChHsbotRstopLR(2, 2); /// (H+ G+, L R) basis
  lChHsbotRstopLR(1, 1) = hb * smu * cosb - forLoops.ub * sinb;
  lChHsbotRstopLR(1, 2) = -ht * mb * cosb - hb * mt * sinb;
  lChHsbotRstopLR(2, 1) = hb * smu * sinb + forLoops.ub * cosb;
  DoubleMatrix lChHsbotRstop12(2, 2);
  temp(1) = lChHsbotRstopLR(1, 1);
  temp(2) = lChHsbotRstopLR(1, 2);
  temp2 = rot2d(thetat) * temp;
  lChHsbotRstop12(1, 1) = temp2(1);
  lChHsbotRstop12(1, 2) = temp2(2);
  temp(1) = lChHsbotRstopLR(2, 1);
  temp(2) = lChHsbotRstopLR(2, 2);
  temp2 = rot2d(thetat) * temp;
  lChHsbotRstop12(2, 1) = temp2(1);
  lChHsbotRstop12(2, 2) = temp2(2);

  for (i=1; i<=4; i++) {
    higgs(1, 1) = higgs(1, 1) + 
      0.5 * (hbsq * dnd(i) - sqr(g) * gdL * cn(i) 
	     / (2.0 * sqr(costhDrbar))) * a0(higgsm(i), q);
    higgs(2, 2) = higgs(2, 2) + 
      0.5 * (hbsq * dnd(i) - sqr(g) * gdR * cn(i) 
	     / (2.0 * sqr(costhDrbar))) * a0(higgsm(i), q);
  }
  for (i=3; i<=4; i++) {
    double a0p = a0(higgsc(i - 2), q);
    higgs(1, 1) = higgs(1, 1) + a0p * 
	  (htsq * dnd(i) +
	   sqr(g) * (gdL / (2.0 * sqr(costhDrbar)) + 0.5) * cn(i));
    higgs(2, 2) = higgs(2, 2) + a0p *
      (hbsq * dnu(i) + sqr(g) * gdR / (2.0 * sqr(costhDrbar)) * cn(i));
  }
  
  int j; for(i=1; i<=4; i++)
    for (j=1; j<=2; j++) {
      double b0p = b0(p, higgsm(i), msbot(j), q);
      higgs(1, 1) = higgs(1, 1) + sqr(lHSbotLSbot12(i, j)) * b0p;
      higgs(1, 2) = higgs(1, 2) + 
	lHSbotLSbot12(i, j) * lHSbotRSbot12(i, j) * b0p;
      higgs(2, 2) = higgs(2, 2) + sqr(lHSbotRSbot12(i, j)) * b0p;
    }
  
  for(i=1; i<=2; i++)
    for (j=1; j<=2; j++) {
      double b0p = b0(p, higgsc(i), mstop(j), q); 
      higgs(1, 1) = higgs(1, 1) + sqr(lChHsbotLstop12(i, j)) * b0p;
      higgs(1, 2) = higgs(1, 2) + 
	lChHsbotLstop12(i, j) * lChHsbotRstop12(i, j) * b0p;
      higgs(2, 2) = higgs(2, 2) + sqr(lChHsbotRstop12(i, j)) * b0p;
    }

  return higgs;
}

template<class SoftPars>
DoubleMatrix Softsusy<SoftPars>::addSbotEweak(double p, DoubleMatrix & electroweak) {
  double    mz         = displayMzRun();
  double    mw         = displayMwRun();
  double    sinthDrbar = calcSinthdrbar();
  double    costhDrbar = sqrt(1.0 - sqr(sinthDrbar));
  double    g  	       = displayGaugeCoupling(2);
  double    gp         = displayGaugeCoupling(1) * sqrt(0.6);
  double    e          = g * sinthDrbar;
  double    thetat     = forLoops.thetat;
  double    thetab     = forLoops.thetab;
  double    thetatau   = forLoops.thetatau;
  double    ct         = cos(thetat);
  double    st         = sin(thetat);
  double    cb         = cos(thetab);
  double    sb         = sin(thetab);
  double    ctau       = cos(thetatau);
  double    stau       = sin(thetatau);
  double    q          = displayMu();
  DoubleVector msbot(2);
  msbot(1)             = forLoops.md(1, 3);
  msbot(2)             = forLoops.md(2, 3);
  DoubleVector mstop(2);
  mstop(1)             = forLoops.mu(1, 3);
  mstop(2)             = forLoops.mu(2, 3);
  DoubleVector mstau(2);
  mstau(1)             = forLoops.me(1, 3);
  mstau(2)             = forLoops.me(2, 3);
  DoubleVector msnu(3);
  msnu(1)              = forLoops.msnu(1);
  msnu(2)              = forLoops.msnu(2);
  msnu(3)              = forLoops.msnu(3);

  electroweak(1, 1) = electroweak(1, 1) + 
    4.0 * sqr(g) / sqr(costhDrbar) * sqr(gdL) * a0(mz, q) +
    2.0 * sqr(g) * a0(mw, q) + sqr(e / 3.0) * 
    (sqr(cb) * ffn(p, msbot(1), 0., q) + 
     sqr(sb) * ffn(p, msbot(2), 0., q));
  electroweak(1, 1) = electroweak(1, 1) +   
    sqr(g) / sqr(costhDrbar) * sqr(gdL) * 
    (sqr(cb) * ffn(p, msbot(1), mz, q) + 
     sqr(sb) * ffn(p, msbot(2), mz, q)) +
    sqr(g) * 0.5 * (sqr(ct) * ffn(p, mstop(1), mw, q) + 
		    sqr(st) * ffn(p, mstop(2), mw, q));
  electroweak(1, 1) = electroweak(1, 1) +   
    sqr(g) * 0.25 * (sqr(cb) * a0(msbot(1), q) + sqr(sb) * a0(msbot(2), q) +
		     2.0 * (sqr(ct) * a0(mstop(1), q) + 
			    sqr(st) * a0(mstop(2), q)));
  electroweak(1, 1) = electroweak(1, 1) + sqr(g) * 
    (-3.0 * 0.25 * (sqr(ct) * a0(mstop(1), q) + sqr(st) * a0(mstop(2), q))
     +3.0 * 0.25 * (sqr(cb) * a0(msbot(1), q) + sqr(sb) * a0(msbot(2), q))
     +      0.25 * (sqr(ctau) * a0(mstau(1), q) + sqr(stau) * a0(mstau(2), q))
     -      0.25 * a0(msnu(3), q)
     -3.0 * 0.25 * (a0(forLoops.mu(1, 2), q) + a0(forLoops.mu(1, 1), q))
     +3.0 * 0.25 * (a0(forLoops.md(1, 2), q) + a0(forLoops.md(1, 1), q))
     +      0.25 * (a0(forLoops.me(1, 2), q) + a0(forLoops.me(1, 1), q))
     -      0.25 * (a0(msnu(2), q) + a0(msnu(1), q)))
    + sqr(gp) * 0.25 * sqr(ydL) * 
    (sqr(cb) * a0(msbot(1), q) + sqr(sb) * a0(msbot(2), q));
  electroweak(1, 1) = electroweak(1, 1) + sqr(gp) * 0.25 * ydL *
    (+3.0 * yuL * (sqr(ct) * a0(mstop(1), q) + sqr(st) * a0(mstop(2), q))
     +3.0 * ydL * (sqr(cb) * a0(msbot(1), q) + sqr(sb) * a0(msbot(2), q))
     +      yeL * (sqr(ctau) * a0(mstau(1), q) + sqr(stau) * a0(mstau(2), q))
     +      ynuL * a0(msnu(3), q)
     +3.0 * yuL * (a0(forLoops.mu(1, 2), q) + a0(forLoops.mu(1, 1), q))
     +3.0 * ydL * (a0(forLoops.md(1, 2), q) + a0(forLoops.md(1, 1), q))
     +      yeL * (a0(forLoops.me(1, 2), q) + a0(forLoops.me(1, 1), q))
     +      ynuL * (a0(msnu(1), q) + a0(msnu(2), q))
     //
     +3.0 * yuR * (sqr(st) * a0(mstop(1), q) + sqr(ct) * a0(mstop(2), q))
     +3.0 * ydR * (sqr(sb) * a0(msbot(1), q) + sqr(cb) * a0(msbot(2), q))
     +yeR * (sqr(stau) * a0(mstau(1), q) + sqr(ctau) * a0(mstau(2), q))
     +3.0 * yuR * (a0(forLoops.mu(2, 2), q) + a0(forLoops.mu(2, 1), q))
     +3.0 * ydR * (a0(forLoops.md(2, 2), q) + a0(forLoops.md(2, 1), q))
     +      yeR * (a0(forLoops.me(2, 2), q) + a0(forLoops.me(2, 1), q)));

  electroweak(2, 2) = electroweak(2, 2) + 
    4.0 * sqr(g) / sqr(costhDrbar) * sqr(gdR) * a0(mz, q) +
    sqr(e / 3.0) * 
    (sqr(sb) * ffn(p, msbot(1), 0., q) + 
     sqr(cb) * ffn(p, msbot(2), 0., q));
  electroweak(2, 2) = electroweak(2, 2) + 
    sqr(g) / sqr(costhDrbar) * sqr(gdR) * 
    (sqr(sb) * ffn(p, msbot(1), mz, q) + 
     sqr(cb) * ffn(p, msbot(2), mz, q));
  electroweak(2, 2) = electroweak(2, 2) + 
    sqr(gp) * 0.25 * sqr(ydR) * 
    (sqr(sb) * a0(msbot(1), q) + sqr(cb) * a0(msbot(2), q));
  electroweak(2, 2) = electroweak(2, 2) + sqr(gp) * 0.25 * ydR *
    (+3.0 * yuL * (sqr(ct) * a0(mstop(1), q) + sqr(st) * a0(mstop(2), q))
     +3.0 * ydL * (sqr(cb) * a0(msbot(1), q) + sqr(sb) * a0(msbot(2), q))
     +      yeL * (sqr(ctau) * a0(mstau(1), q) + sqr(stau) * a0(mstau(2), q))
     +      ynuL * a0(msnu(3), q) 
     +3.0 * yuL * (a0(forLoops.mu(1, 2), q) + a0(forLoops.mu(1, 1), q))
     +3.0 * ydL * (a0(forLoops.md(1, 2), q) + a0(forLoops.md(1, 1), q))
     +      yeL * (a0(forLoops.me(1, 2), q) + a0(forLoops.me(1, 1), q))
     +      ynuL * (a0(msnu(1), q) + a0(msnu(2), q))
     +3.0 * yuR * (sqr(st) * a0(mstop(1), q) + sqr(ct) * a0(mstop(2), q))
     +3.0 * ydR * (sqr(sb) * a0(msbot(1), q) + sqr(cb) * a0(msbot(2), q))
     +      yeR * (sqr(stau) * a0(mstau(1), q) + sqr(ctau) * a0(mstau(2), q))
     +3.0 * yuR * (a0(forLoops.mu(2, 2), q) + a0(forLoops.mu(2, 1), q))
     +3.0 * ydR * (a0(forLoops.md(2, 2), q) + a0(forLoops.md(2, 1), q))
     +      yeR * (a0(forLoops.me(2, 2), q) + a0(forLoops.me(2, 1), q)));

  electroweak(1, 2) = electroweak(1, 2) + 
    sqr(gp) * sb * cb * 0.25 * ydL * ydR * 
    (a0(msbot(1), q) - a0(msbot(2), q)) +
    sqr(e / 3.0) * sb * cb * 
    (ffn(p, msbot(1), 0., q) - ffn(p, msbot(2), 0., q)) -
    sqr(g) / sqr(costhDrbar) * gdL * gdR * sb * cb * 
    (ffn(p, msbot(1), mz, q) - ffn(p, msbot(2), mz, q));

  return electroweak;
}

template<class SoftPars>
DoubleMatrix Softsusy<SoftPars>::addSbotNeutralino(double p, double /* mt */, DoubleMatrix & neutralino) {
  double g  = displayGaugeCoupling(2);
  double gp = displayGaugeCoupling(1) * sqrt(0.6);
  double hb = forLoops.hb;
  double q  = displayMu();
  double mb = forLoops.mb;

  ComplexMatrix n(forLoops.nBpmz);
  DoubleVector mneut(forLoops.mnBpmz);

  int rank = mneut.displayEnd();

  /// Neutralino Feynman rules 
  DoubleVector aPsi0BSbotr(rank), bPsi0BSbotr(rank), 
    aPsi0BSbotl(rank), bPsi0BSbotl(rank);
  aPsi0BSbotr(1) = ydR * gp / root2;
  bPsi0BSbotl(1) = gp / (3.0 * root2);
  bPsi0BSbotl(2) = -g / root2;
  aPsi0BSbotl(3) = hb;
  bPsi0BSbotr(3) = hb;

  ComplexVector aChi0BSbotl(rank), bChi0BSbotl(rank), aChi0BSbotr(rank),
    bChi0BSbotr(rank);
  aChi0BSbotl = n.complexConjugate() * aPsi0BSbotl;
  bChi0BSbotl = n * bPsi0BSbotl;
  aChi0BSbotr = n.complexConjugate() * aPsi0BSbotr;
  bChi0BSbotr = n * bPsi0BSbotr;

  DoubleVector gChi0BotSbotLL(rank), fChi0BotSbotLL(rank);
  DoubleVector gChi0BotSbotLR(rank), fChi0BotSbotLR(rank);
  DoubleVector gChi0BotSbotRR(rank), fChi0BotSbotRR(rank);

  int i;
  for (i=1; i<=rank; i++) {
    fChi0BotSbotLL(i) = (aChi0BSbotl(i) * aChi0BSbotl(i).conj() + 
      bChi0BSbotl(i) * bChi0BSbotl(i).conj()).real();
    gChi0BotSbotLL(i) = (bChi0BSbotl(i).conj() * aChi0BSbotl(i) + 
      bChi0BSbotl(i) * aChi0BSbotl(i).conj()).real();
    fChi0BotSbotRR(i) = (aChi0BSbotr(i) * aChi0BSbotr(i).conj() + 
      bChi0BSbotr(i) * bChi0BSbotr(i).conj()).real();
    gChi0BotSbotRR(i) = (bChi0BSbotr(i).conj() * aChi0BSbotr(i) + 
      bChi0BSbotr(i) * aChi0BSbotr(i).conj()).real();
    fChi0BotSbotLR(i) = (aChi0BSbotr(i) * aChi0BSbotl(i).conj() + 
      bChi0BSbotr(i) * bChi0BSbotl(i).conj()).real();
    gChi0BotSbotLR(i) = (bChi0BSbotl(i).conj() * aChi0BSbotr(i) + 
      bChi0BSbotr(i) * aChi0BSbotl(i).conj()).real();
  }

  for (i=1; i<=rank; i++) {
    double one = gfn(p, mneut(i), mb, q);
    double two = 2.0 * mneut(i) * mb * b0(p, mneut(i), mb, q);
    neutralino(1, 1) = neutralino(1, 1) +
      fChi0BotSbotLL(i) * one - gChi0BotSbotLL(i) * two;
    neutralino(2, 2) = neutralino(2, 2) +
      fChi0BotSbotRR(i) * one - gChi0BotSbotRR(i) * two;
    neutralino(1, 2) = neutralino(1, 2) +
      fChi0BotSbotLR(i) * one - gChi0BotSbotLR(i) * two;
  }

  return neutralino;
}

template<class SoftPars>
DoubleMatrix Softsusy<SoftPars>::addSbotChargino(double p, double mt, DoubleMatrix & chargino) {
  double g  = displayGaugeCoupling(2);
  double ht = forLoops.ht;
  double hb = forLoops.hb;
  double q  = displayMu();

  ComplexMatrix n(forLoops.nBpmz);
  ComplexMatrix u(forLoops.uBpmz), v(forLoops.vBpmz); 
  DoubleVector mch(forLoops.mchBpmz);

  /// Chargino Feynman Rules
  DoubleVector bPsicTSbotl(2), bPsicTSbotr(2), aPsicTSbotl(2), aPsicTSbotr(2); 
  bPsicTSbotl(1) = g;
  bPsicTSbotr(2) = -hb;
  aPsicTSbotl(2) = -ht;
  
  ComplexVector aChicTSbotr(2), aChicTSbotl(2), bChicTSbotl(2), bChicTSbotr(2);

  aChicTSbotl = v.complexConjugate() * aPsicTSbotl;
  bChicTSbotl = u * bPsicTSbotl;
  aChicTSbotr = v.complexConjugate() * aPsicTSbotr;
  bChicTSbotr = u * bPsicTSbotr;

  DoubleVector fChTSbotLL(2), gChTSbotLL(2) ;
  DoubleVector fChTSbotLR(2), gChTSbotLR(2); 
  DoubleVector fChTSbotRR(2), gChTSbotRR(2); 

  int i;
  for (i=1; i<=2; i++) {
    fChTSbotLL(i) = (aChicTSbotl(i).conj() * aChicTSbotl(i) +
		      bChicTSbotl(i).conj() * bChicTSbotl(i)).real();
    gChTSbotLL(i) = (bChicTSbotl(i).conj() * aChicTSbotl(i) +
		      aChicTSbotl(i).conj() * bChicTSbotl(i)).real();
    fChTSbotLR(i) = (aChicTSbotl(i).conj() * aChicTSbotr(i) +
		      bChicTSbotl(i).conj() * bChicTSbotr(i)).real();
    gChTSbotLR(i) = (bChicTSbotl(i).conj() * aChicTSbotr(i) +
		      aChicTSbotl(i).conj() * bChicTSbotr(i)).real();
    fChTSbotRR(i) = (aChicTSbotr(i).conj() * aChicTSbotr(i) +
		      bChicTSbotr(i).conj() * bChicTSbotr(i)).real();
    gChTSbotRR(i) = (bChicTSbotr(i).conj() * aChicTSbotr(i) +
		      aChicTSbotr(i).conj() * bChicTSbotr(i)).real();
  }

  for (i=1; i<=2; i++) {
    double one = gfn(p, mch(i), mt, q);
    double two = mch(i) * mt * b0(p, mch(i), mt, q) * 2.0;
    chargino(1, 1) = chargino(1, 1) + fChTSbotLL(i) * one -
      gChTSbotLL(i) * two;
    chargino(1, 2) = chargino(1, 2) + fChTSbotLR(i) * one -
      gChTSbotLR(i) * two;
    chargino(2, 2) = chargino(2, 2) + fChTSbotRR(i) * one - 
      gChTSbotRR(i) * two;
  }

  return chargino;
}

/// Checked 16.08.2005
template<class SoftPars>
void Softsusy<SoftPars>::addSbotCorrection(double p, DoubleMatrix & mass, double mt) {

/// No point adding radiative corrections to tachyonic particles
  if (mass(1, 1) < 0.0 || mass(2, 2) < 0.0) { 
    flagTachyon(sbottom);
    if (mass(1, 1) < 0.) mass(1, 1) = EPSTOL;
    else mass(2, 2) = EPSTOL;
    return;
  }

  /// one-loop correction matrix
  DoubleMatrix piSq(2, 2); /// Self-energy matrix

  DoubleMatrix strong(2, 2), stop(2, 2), sbottom(2, 2), 
    higgs(2, 2), electroweak(2, 2), chargino(2, 2), neutralino(2, 2);

  /// LCT: Strong contributions
  addSbotQCD(p, mt, strong);
  /// LCT: Stop, sbottom and stau contributions
  addSbotSfermion(p, mt, stop, sbottom);
  /// LCT: Higgs contribution
  addSbotHiggs(p, mt, higgs);
  /// LCT: Electroweak contribution
  addSbotEweak(p, electroweak);
  /// LCT: Chargino contribution
  addSbotChargino(p, mt, chargino);
  /// LCT: Neutralino contribution
  addSbotNeutralino(p, mt, neutralino);
 
  piSq = 1.0 / (16.0 * sqr(PI)) * 
    (strong + stop + sbottom + higgs + electroweak + chargino + neutralino);

  piSq(2, 1) = piSq(1, 2);

  mass = mass - piSq;	
}

template<class SoftPars>
DoubleMatrix Softsusy<SoftPars>::addSupQCD(double p1, double p2, int family, DoubleMatrix & strong) {
  double g3sq = sqr(displayGaugeCoupling(3));
  double mg   = forLoops.mGluino;
  double q    = displayMu(); 
  DoubleVector msup(2); 
  msup(1)     = forLoops.mu(1, family);
  msup(2)     = forLoops.mu(2, family);

  double a0t1 = a0(msup(1), q), a0t2 = a0(msup(2), q);
  double ft1  = ffn(p1, msup(1), 0.0, q), ft2 = ffn(p2, msup(2), 0.0, q);
  double ggt1 = gfn(p1, mg, 0., q), ggt2 = gfn(p2, mg, 0., q);

  strong(1, 1) = 4.0 * g3sq / 3.0 * (2.0 * ggt1 + (ft1 + a0t1));
  strong(2, 2) = 4.0 * g3sq / 3.0 * (2.0 * ggt2 + (ft2 + a0t2));

  return strong;
}

template<class SoftPars>
DoubleMatrix Softsusy<SoftPars>::addSupHiggs(double p1, double p2, int family, DoubleMatrix & higgs) {
  double    sinthDrbar  = calcSinthdrbar();
  double    costhDrbar  = sqrt(1.0 - sqr(sinthDrbar));
  double    costhDrbar2 = 1.0 - sqr(sinthDrbar);
  double    alpha       = forLoops.thetaH;
  double    g	        = displayGaugeCoupling(2);
  double    beta        = atan(displayTanb());
  double    q           = displayMu();
  double    sinb        = sin(beta), cosb = cos(beta);
  double    mz          = displayMzRun();
  DoubleVector msup(2);
  msup(1)               = forLoops.mu(1, family);
  msup(2)               = forLoops.mu(2, family);
  DoubleVector msd(2);
  msd(1)                = forLoops.md(1, family);
  msd(2)                = forLoops.md(2, family);

  /// define Higgs vector in 't-Hooft Feynman gauge, and couplings:
  DoubleVector higgsm(4), higgsc(2);
  DoubleVector dnu(4), dnd(4), cn(4);
  assignHiggsSfermions(higgsm, higgsc, dnu, dnd, cn, beta);

  /// top now stands for up
  DoubleMatrix lsStopLStopLR(4, 2), lsStopLStop12(4, 2);
  DoubleMatrix lsStopRStopLR(4, 2), lsStopRStop12(4, 2);
  /// Order (s1 s2 G A, L R)
  lsStopLStopLR(1, 1) = g * mz * guL * cosb / costhDrbar;
  lsStopLStopLR(2, 1) = -g * mz * guL * sinb / costhDrbar;

  lsStopRStopLR(1, 1) = lsStopLStopLR(1, 2);
  lsStopRStopLR(1, 2) = g * mz * guR * cosb / costhDrbar;
  lsStopRStopLR(2, 1) = lsStopLStopLR(2, 2);
  lsStopRStopLR(2, 2) = -g * mz * guR * sinb / costhDrbar;
  lsStopRStopLR(3, 1) = -lsStopLStopLR(3, 2);
  lsStopRStopLR(4, 1) = -lsStopLStopLR(4, 2);

  /// Mix stops up
  int i;
  DoubleVector temp(2), temp2(2);
  for (i=1; i<=4; i++) {
    temp(1) = lsStopLStopLR(i, 1);
    temp(2) = lsStopLStopLR(i, 2);
    temp2 = temp;
    lsStopLStop12(i, 1) = temp2(1);
    lsStopLStop12(i, 2) = temp2(2);
    temp(1) = lsStopRStopLR(i, 1);
    temp(2) = lsStopRStopLR(i, 2);
    temp2 = temp;
    lsStopRStop12(i, 1) = temp2(1);
    lsStopRStop12(i, 2) = temp2(2);
 }

  DoubleMatrix lHStopLStop12(lsStopLStop12), lHStopRStop12(lsStopRStop12);
  /// Mix CP-even Higgses up
  for (i=1; i<=2; i++) { /// i is the s1/s2 label
    temp(1) = lsStopLStop12(1, i);
    temp(2) = lsStopLStop12(2, i);
    temp2 = rot2d(alpha) * temp;
    lHStopLStop12(1, i) = temp2(1);
    lHStopLStop12(2, i) = temp2(2);
    temp(1) = lsStopRStop12(1, i);
    temp(2) = lsStopRStop12(2, i);
    temp2 = rot2d(alpha) * temp;
    lHStopRStop12(1, i) = temp2(1);
    lHStopRStop12(2, i) = temp2(2);
  }

  /// Charged Higgs Feynman rules (H+ G+, L R) basis
  DoubleMatrix lChHstopLsbotLR(2, 2); 
  lChHstopLsbotLR(1, 1) = g * displayMwRun() * sin(2.0 * beta) / root2;
  lChHstopLsbotLR(2, 1) =-g * displayMwRun() * cos(2.0 * beta) / root2;

  DoubleMatrix lChHstopLsbot12(2, 2);
  temp(1) = lChHstopLsbotLR(1, 1);
  temp(2) = lChHstopLsbotLR(1, 2);
  temp2 = temp;
  lChHstopLsbot12(1, 1) = temp2(1);
  lChHstopLsbot12(1, 2) = temp2(2);
  temp(1) = lChHstopLsbotLR(2, 1);
  temp(2) = lChHstopLsbotLR(2, 2);
  temp2 = temp;
  lChHstopLsbot12(2, 1) = temp2(1);
  lChHstopLsbot12(2, 2) = temp2(2);

  DoubleMatrix lChHstopRsbotLR(2, 2); /// (H+ G+, L R) basis
  DoubleMatrix lChHstopRsbot12(2, 2);
  temp(1) = lChHstopRsbotLR(1, 1);
  temp(2) = lChHstopRsbotLR(1, 2);
  temp2 = temp;
  lChHstopRsbot12(1, 1) = temp2(1);
  lChHstopRsbot12(1, 2) = temp2(2);
  temp(1) = lChHstopRsbotLR(2, 1);
  temp(2) = lChHstopRsbotLR(2, 2);
  temp2 = temp;
  lChHstopRsbot12(2, 1) = temp2(1);
  lChHstopRsbot12(2, 2) = temp2(2);

  /// LCT: Contributions start here
 for (i=1; i<=4; i++) {
    higgs(1, 1) = higgs(1, 1) + 
      0.5 * (- sqr(g) * guL * 0.5 / sqr(costhDrbar) * cn(i)) 
      * a0(higgsm(i), q);
    higgs(2, 2) = higgs(2, 2) + 
      0.5 * (- sqr(g) * guR * 0.5 / sqr(costhDrbar) * cn(i)) 
      * a0(higgsm(i), q);
  }
  for (i=3; i<=4; i++) { 
    double a0p = a0(higgsc(i - 2), q);
    higgs(1, 1) += (sqr(g) * (guL * 0.5 / costhDrbar2 - 0.5) * cn(i))* a0p;
    higgs(2, 2) += (sqr(g) * guR * 0.5 / costhDrbar2 * cn(i))* a0p;    
  }
  int j; for(i=1; i<=4; i++)
    for (j=1; j<=2; j++) {
      double b0p = b0(p1, higgsm(i), msup(j), q);
      higgs(1, 1) += sqr(lHStopLStop12(i, j)) * b0p;
      b0p = b0(p2, higgsm(i), msup(j), q);
      higgs(2, 2) += sqr(lHStopRStop12(i, j)) * b0p;
    }
  for(i=1; i<=2; i++)
    for (j=1; j<=2; j++) {
      double b0p = b0(p1, higgsc(i), msd(j), q); 
      higgs(1, 1) += sqr(lChHstopLsbot12(i, j)) * b0p;
      b0p = b0(p2, higgsc(i), msd(j), q); 
      higgs(2, 2) += sqr(lChHstopRsbot12(i, j)) * b0p;
    }

  return higgs;
}

template<class SoftPars>
DoubleMatrix Softsusy<SoftPars>::addSupEweak(double p1, double p2, int family, DoubleMatrix & electroweak) {
  double    sinthDrbar  = calcSinthdrbar();
  double    costhDrbar  = sqrt(1.0 - sqr(sinthDrbar));
  double    costhDrbar2 = 1.0 - sqr(sinthDrbar);
  double    g	        = displayGaugeCoupling(2);
  double    gp          = displayGaugeCoupling(1) * sqrt(0.6);
  double    e           = g * sinthDrbar; /// DRbar value of e
  double    thetat      = forLoops.thetat;
  double    thetab      = forLoops.thetab;
  double    thetatau    = forLoops.thetatau;
  double    ct          = cos(thetat);
  double    st          = sin(thetat);
  double    cb          = cos(thetab);
  double    sb          = sin(thetab);
  double    ctau        = cos(thetatau);
  double    stau        = sin(thetatau);
  double    q           = displayMu();
  double    mz          = displayMzRun();
  DoubleVector msup(2);
  msup(1)               = forLoops.mu(1, family);
  msup(2)               = forLoops.mu(2, family);
  DoubleVector msd(2);
  msd(1)                = forLoops.md(1, family);
  msd(2)                = forLoops.md(2, family);

  /// EW bosons
  electroweak(1, 1) += 
    4.0 * sqr(g) / costhDrbar2 * sqr(guL) * a0(mz, q) + 
    2.0 * sqr(g) * a0(displayMwRun(), q) + sqr(2.0 / 3.0 * e) * 
    ffn(p1, msup(1), 0.0, q) 
    + sqr(g * guL / costhDrbar) * ffn(p1, msup(1), mz, q) +
    sqr(g) * 0.5 * ffn(p1, msd(1), displayMwRun(), q) +
    sqr(g) * 0.25 * 
    (a0(msup(1), q) + 2.0 * a0(msd(1), q)) +
    sqr(g) * 0.5 * 
    (1.5 * a0(forLoops.mu(1, 1), q) + 
     1.5 * a0(forLoops.mu(1, 2), q) +
     1.5 * (sqr(ct) * a0(forLoops.mu(1, 3), q) + 
	    sqr(st) * a0(forLoops.mu(2, 3), q)) -
     1.5 * a0(forLoops.md(1, 1), q) -
     1.5 * a0(forLoops.md(1, 2), q) -
     1.5 * (sqr(cb) * a0(forLoops.md(1, 3), q) +
	    sqr(sb) * a0(forLoops.md(2, 3), q)) +
     0.5 * (a0(forLoops.msnu(1), q) + 
	    a0(forLoops.msnu(2), q) +
	    a0(forLoops.msnu(3), q)) -
     0.5 * (a0(forLoops.me(1, 1), q) + a0(forLoops.me(1, 2), q) +
	    sqr(ctau) * a0(forLoops.me(1, 3), q) +
	    sqr(stau) * a0(forLoops.me(2, 3), q))) +
    sqr(gp) * 0.25 * sqr(yuL) * a0(msup(1), q) +
    sqr(gp) * 0.25 * yuL *  
    (3.0 * yuL * (a0(forLoops.mu(1, 1), q) + 
		  a0(forLoops.mu(1, 2), q) + 
		  sqr(ct) * a0(forLoops.mu(1, 3), q) + 
		  sqr(st) * a0(forLoops.mu(2, 3), q)) +
     3.0 * yuR * (a0(forLoops.mu(2, 1), q) + 
		  a0(forLoops.mu(2, 2), q) + 
		  sqr(st) * a0(forLoops.mu(1, 3), q) + 
		  sqr(ct) * a0(forLoops.mu(2, 3), q)) +
     3.0 * ydL * (a0(forLoops.md(1, 1), q) + 
		  a0(forLoops.md(1, 2), q) + 
		  sqr(cb) * a0(forLoops.md(1, 3), q) + 
		  sqr(sb) * a0(forLoops.md(2, 3), q)) +
     3.0 * ydR * (a0(forLoops.md(2, 1), q) + 
		  a0(forLoops.md(2, 2), q) + 
		  sqr(sb) * a0(forLoops.md(1, 3), q) + 
		  sqr(cb) * a0(forLoops.md(2, 3), q)) +
     yeL * (a0(forLoops.me(1, 1), q) + 
	    a0(forLoops.me(1, 2), q) + 
	    sqr(ctau) * a0(forLoops.me(1, 3), q) + 
	    sqr(stau) * a0(forLoops.me(2, 3), q)) +
     yeR * (a0(forLoops.me(2, 1), q) + 
	    a0(forLoops.me(2, 2), q) + 
	    sqr(stau) * a0(forLoops.me(1, 3), q) + 
	    sqr(ctau) * a0(forLoops.me(2, 3), q)) +
     ynuL * (a0(forLoops.msnu(1), q) + 
	     a0(forLoops.msnu(2), q) +
	     a0(forLoops.msnu(3), q)));
     
  electroweak(2, 2) += 
    4.0 * sqr(g) / costhDrbar2 * sqr(guR) * a0(mz, q) + 
    sqr(2.0 / 3.0 * e) * 
    (ffn(p2, msup(2), 0.0, q))
    + sqr(g * guR / costhDrbar) * 
    (ffn(p2, msup(2), mz, q)) +
    sqr(gp) * 0.25 * sqr(yuR) * 
    (a0(msup(2), q)) +
    sqr(gp) * 0.25 * yuR * 
    (3.0 * yuL * (a0(forLoops.mu(1, 1), q) + 
		  a0(forLoops.mu(1, 2), q) + 
		  sqr(ct) * a0(forLoops.mu(1, 3), q) + 
		  sqr(st) * a0(forLoops.mu(2, 3), q)) +
     3.0 * yuR * (a0(forLoops.mu(2, 1), q) + 
		  a0(forLoops.mu(2, 2), q) + 
		  sqr(st) * a0(forLoops.mu(1, 3), q) + 
		  sqr(ct) * a0(forLoops.mu(2, 3), q)) +
     3.0 * ydL * (a0(forLoops.md(1, 1), q) + 
		  a0(forLoops.md(1, 2), q) + 
		  sqr(cb) * a0(forLoops.md(1, 3), q) + 
		  sqr(sb) * a0(forLoops.md(2, 3), q)) +
     3.0 * ydR * (a0(forLoops.md(2, 1), q) + 
		  a0(forLoops.md(2, 2), q) + 
		  sqr(sb) * a0(forLoops.md(1, 3), q) + 
		  sqr(cb) * a0(forLoops.md(2, 3), q)) +
     yeL * (a0(forLoops.me(1, 1), q) + 
	    a0(forLoops.me(1, 2), q) + 
	    sqr(ctau) * a0(forLoops.me(1, 3), q) + 
	    sqr(stau) * a0(forLoops.me(2, 3), q)) +
     yeR * (a0(forLoops.me(2, 1), q) + 
	    a0(forLoops.me(2, 2), q) + 
	    sqr(stau) * a0(forLoops.me(1, 3), q) + 
	    sqr(ctau) * a0(forLoops.me(2, 3), q)) +
     ynuL * (a0(forLoops.msnu(1), q) + 
	     a0(forLoops.msnu(2), q) +
	     a0(forLoops.msnu(3), q)));

  return electroweak;
}

template<class SoftPars>
DoubleMatrix Softsusy<SoftPars>::addSupNeutralino(double p1, double p2, int /* family */, DoubleMatrix & neutralino) {
  double q  = displayMu();
  double g  = displayGaugeCoupling(2);
  double gp = displayGaugeCoupling(1) * sqrt(0.6);

  ComplexMatrix n(forLoops.nBpmz);
  DoubleVector mneut(forLoops.mnBpmz);

  int rank = mneut.displayEnd();

  /// Neutralino Feynman rules
  DoubleVector aPsi0TStopr(rank), bPsi0TStopr(rank), 
    aPsi0TStopl(rank), bPsi0TStopl(rank); 
  aPsi0TStopr(1) = - 4.0 * gp / (3.0 * root2);
  bPsi0TStopl(1) = gp / (3.0 * root2);
  bPsi0TStopl(2) = g / root2;

 ComplexVector aChi0TStopl(rank), bChi0TStopl(rank), 
   aChi0TStopr(rank), bChi0TStopr(rank);
  aChi0TStopl = n.complexConjugate() * aPsi0TStopl;
  bChi0TStopl = n * bPsi0TStopl;
  aChi0TStopr = n.complexConjugate() * aPsi0TStopr;
  bChi0TStopr = n * bPsi0TStopr;

  DoubleVector gChi0TopStopLL(rank), fChi0TopStopLL(rank);
  DoubleVector gChi0TopStopLR(rank), fChi0TopStopLR(rank);
  DoubleVector gChi0TopStopRR(rank), fChi0TopStopRR(rank);

  int i;
  for (i=1; i<=rank; i++) {
    fChi0TopStopLL(i) = (aChi0TStopl(i) * aChi0TStopl(i).conj() + 
      bChi0TStopl(i) * bChi0TStopl(i).conj()).real();
    gChi0TopStopLL(i) = (bChi0TStopl(i).conj() * aChi0TStopl(i) + 
      bChi0TStopl(i) * aChi0TStopl(i).conj()).real();
    fChi0TopStopRR(i) = (aChi0TStopr(i) * aChi0TStopr(i).conj() + 
      bChi0TStopr(i) * bChi0TStopr(i).conj()).real();
    gChi0TopStopRR(i) = (bChi0TStopr(i).conj() * aChi0TStopr(i) + 
      bChi0TStopr(i) * aChi0TStopr(i).conj()).real();
    fChi0TopStopLR(i) = (aChi0TStopr(i) * aChi0TStopl(i).conj() + 
      bChi0TStopr(i) * bChi0TStopl(i).conj()).real();
    gChi0TopStopLR(i) = (bChi0TStopl(i).conj() * aChi0TStopr(i) + 
      bChi0TStopr(i) * aChi0TStopl(i).conj()).real();
  }

 for (i=1; i<=rank; i++) {
    double one = gfn(p1, mneut(i), 0., q);
    double two = 0.;
    neutralino(1, 1) = neutralino(1, 1) +
      fChi0TopStopLL(i) * one - gChi0TopStopLL(i) * two;
    one = gfn(p2, mneut(i), 0., q);
    neutralino(2, 2) = neutralino(2, 2) +
      fChi0TopStopRR(i) * one - gChi0TopStopRR(i) * two;
  }
  return neutralino;
}

template<class SoftPars>
DoubleMatrix Softsusy<SoftPars>::addSupChargino(double p1, double p2, int /* family */, DoubleMatrix & chargino) {
  double g  = displayGaugeCoupling(2);
  double q  = displayMu();

  ComplexMatrix n(forLoops.nBpmz);
  ComplexMatrix u(forLoops.uBpmz), v(forLoops.vBpmz); 
  DoubleVector mch(forLoops.mchBpmz);

  /// Chargino Feynman Rules
  /// LCT: Here stop <==> sup and sbot <==> sdown
  DoubleVector bPsicBStopl(2), bPsicBStopr(2), aPsicBStopl(2), aPsicBStopr(2); 
  aPsicBStopl(1) = g;
  
  ComplexVector aChicBStopr(2), aChicBStopl(2), bChicBStopl(2), bChicBStopr(2);

  aChicBStopl = v.complexConjugate() * aPsicBStopl;
  bChicBStopl = u * bPsicBStopl;
  aChicBStopr = v.complexConjugate() * aPsicBStopr;
  bChicBStopr = u * bPsicBStopr;

  DoubleVector fChBStopLL(2), gChBStopLL(2) ;
  DoubleVector fChBStopLR(2), gChBStopLR(2); 
  DoubleVector fChBStopRR(2), gChBStopRR(2); 

  int i;
  for (i=1; i<=2; i++) {
    fChBStopLL(i) = (aChicBStopl(i).conj() * aChicBStopl(i) +
		      bChicBStopl(i).conj() * bChicBStopl(i)).real();
    gChBStopLL(i) = (bChicBStopl(i).conj() * aChicBStopl(i) +
		      aChicBStopl(i).conj() * bChicBStopl(i)).real();
    fChBStopLR(i) = (aChicBStopl(i).conj() * aChicBStopr(i) +
		      bChicBStopl(i).conj() * bChicBStopr(i)).real();
    gChBStopLR(i) = (bChicBStopl(i).conj() * aChicBStopr(i) +
		      aChicBStopl(i).conj() * bChicBStopr(i)).real();
    fChBStopRR(i) = (aChicBStopr(i).conj() * aChicBStopr(i) +
		      bChicBStopr(i).conj() * bChicBStopr(i)).real();
    gChBStopRR(i) = (bChicBStopr(i).conj() * aChicBStopr(i) +
		      aChicBStopr(i).conj() * bChicBStopr(i)).real();
  }

  for (i=1; i<=2; i++) {
    double one = gfn(p1, mch(i), 0., q);
    double two = 0.;
    chargino(1, 1) = chargino(1, 1) + fChBStopLL(i) * one -
      gChBStopLL(i) * two;
    one = gfn(p2, mch(i), 0., q);
    chargino(2, 2) = chargino(2, 2) + fChBStopRR(i) * one - 
      gChBStopRR(i) * two;
  }
  return chargino;
}

/// As in BPMZ appendix, INCLUDING weak boson loops.
template<class SoftPars>
void Softsusy<SoftPars>::addSupCorrection(DoubleMatrix & mass, int family) {

/// No point adding radiative corrections to tachyonic particles
  if (mass(1, 1) < 0.0 || mass(2, 2) < 0.0) { 
    if (mass(1, 1) < 0.) mass(1, 1) = EPSTOL;
    else mass(2, 2) = EPSTOL;
    if (family == 1) flagTachyon(sup);
    else if (family == 2) flagTachyon(scharm);
    return;
  }
  /// one-loop correction matrix
  DoubleMatrix piSq(2, 2); /// Self-energy matrix
  DoubleVector msup(2);
  msup(1)          = forLoops.mu(1, family);
  msup(2)          = forLoops.mu(2, family);

  DoubleMatrix strong(2, 2), stop(2, 2), sbottom(2, 2), 
    higgs(2, 2), electroweak(2, 2), chargino(2, 2), neutralino(2, 2);
  double p1 = msup(1), p2 = msup(2);

  /// LCT: QCD contribution
  addSupQCD(p1, p2, family, strong);
  /// LCT: Higgs contribution
  addSupHiggs(p1, p2, family, higgs);
  /// LCT: Weak contribution
  addSupEweak(p1, p2, family, electroweak);
   /// LCT: Chargino contribution
  addSupChargino(p1, p2, family, chargino);
  /// LCT: Neutralino contribution
  addSupNeutralino(p1, p2, family, neutralino);

  piSq = 1.0 / (16.0 * sqr(PI)) * 
    (strong + stop + sbottom + higgs + electroweak + chargino + neutralino);

  piSq(2, 1) = piSq(1, 2);
  mass = mass - piSq;	  
}

template<class SoftPars>
void Softsusy<SoftPars>::doUpSquarks(double mt, double pizztMS, double sinthDRbarMS, int accuracy) { 

  /// first two families are simpler
  int family; for (family = 1; family <= 2; family++) {
    
    DoubleMatrix mStopSquared(2, 2);
    treeUpSquark(mStopSquared, mt, pizztMS, sinthDRbarMS, family);
    if (accuracy > 0) addSupCorrection(mStopSquared, family); 
    double theta;
    DoubleVector 
      physicalStopMassesSquared(mStopSquared.sym2by2(theta));
    if (physicalStopMassesSquared(1) < 0.0 ||
	physicalStopMassesSquared(2) < 0.0) {
      if (family == 1) flagTachyon(sup);
      else if (family == 2) flagTachyon(scharm);
    }

    DoubleVector physicalStopMasses(physicalStopMassesSquared.apply(ccbSqrt));

    physpars.mu(1, family) = physicalStopMasses(1);
    physpars.mu(2, family) = physicalStopMasses(2);
  }
  
  family = 3;
  
  DoubleMatrix mStopSquared(2, 2), mStopSquared2(2, 2);
  treeUpSquark(mStopSquared, mt, pizztMS, sinthDRbarMS, family);
  mStopSquared2 = mStopSquared; /// StopSquared2 is now tree-level
  
  /// one loop corrections 
  if (accuracy > 0) {
    double pLight = minimum(forLoops.mu(1, 3), forLoops.mu(2, 3));
    double pHeavy = maximum(forLoops.mu(1, 3), forLoops.mu(2, 3));

    addStopCorrection(pLight, mStopSquared, mt); 
    addStopCorrection(pHeavy, mStopSquared2, mt);
  }

  double theta;
  DoubleVector physicalStopMassesSquared2(mStopSquared2.sym2by2(theta));
  DoubleVector physicalStopMassesSquared(mStopSquared.sym2by2(theta));
  physpars.thetat = theta; /// theta is calculated at p=mstop1

  if (physicalStopMassesSquared(1) < 0.0 || physicalStopMassesSquared(2) < 0.0
      ||  physicalStopMassesSquared2(1) < 0.0 || 
      physicalStopMassesSquared2(2) < 0.0) {
    flagTachyon(stop);
  }
  
  DoubleVector physicalStopMasses(physicalStopMassesSquared.apply(ccbSqrt));
  DoubleVector 
      physicalStopMasses2(physicalStopMassesSquared2.apply(ccbSqrt));
  
  double lightStopMass = minimum(physicalStopMasses(1), physicalStopMasses(2));
  double heavyStopMass = maximum(physicalStopMasses2(1), 
				 physicalStopMasses2(2));
  /// twisted measures the ordering of the sbottom masses. If mstop1 > mstop2,
  /// twisted is defined to be true (mstop2 > mstop1 is defined "untwisted").
  bool twisted = false;
  if (physicalStopMassesSquared(1) > physicalStopMassesSquared(2)) 
    twisted = true;

  if (twisted) {
    physpars.mu(1, family) = heavyStopMass; 
    physpars.mu(2, family) = lightStopMass;
  } else {
    physpars.mu(1, family) = lightStopMass;
    physpars.mu(2, family) = heavyStopMass;
  }
}

template<class SoftPars>
void Softsusy<SoftPars>::treeDownSquark(DoubleMatrix & mass, double mbrun, 
                                        double /* pizztMS */, double sinthDRbarMS,
                                        int family) {
  const double cd = 1.0 / 3.0;
  double mz2 = sqr(displayMzRun()), mb2 = sqr(mbrun);
  double beta = atan(displayTanb()), mu = displaySusyMu(),
    c2b = cos(2 * beta), tanb = displayTanb(), sinth2 = sqr(sinthDRbarMS);

  mass(1, 1) = displaySoftMassSquared(mQl, family, family) +
      mz2 * (-0.5 + cd * sinth2) * c2b;
  mass(2, 2) = displaySoftMassSquared(mDr, family, family) -
      mz2 * cd * sinth2 * c2b;

  if (family != 3) mass(1, 2) = 0.0;
  else {
    if (fabs(forLoops.hb) < EPSTOL) 
      mass(1, 2) = mbrun * (displaySoftA(DA, 3, 3) - mu * tanb);
    else
      mass(1, 2) = mbrun * (forLoops.ub / forLoops.hb - mu * tanb);

    mass(1, 1) = mass(1, 1) + mb2;
    mass(2, 2) = mass(2, 2) + mb2;
  }

  mass(2, 1) = mass(1, 2);
}


template<class SoftPars>
void Softsusy<SoftPars>::doDownSquarks(double mb, double pizztMS, double
			    sinthDRbarMS, int accuracy, double mt) {
  int family; for (family = 1; family <= 2; family++) {
    
    DoubleMatrix mSbotSquared(2, 2);
    treeDownSquark(mSbotSquared, mb, pizztMS, sinthDRbarMS, family);

    if (accuracy > 0) addSdownCorrection(mSbotSquared, family);

    double theta;
    DoubleVector 
      physicalSbotMassesSquared(mSbotSquared.sym2by2(theta));
    if (physicalSbotMassesSquared(1) < 0.0 || 
	physicalSbotMassesSquared(2) < 0.0) {
      if (family == 1) flagTachyon(sdown);
      else if (family == 2) flagTachyon(sstrange);
    }

    DoubleVector physicalSbotMasses(physicalSbotMassesSquared.apply(ccbSqrt));

    physpars.md(1, family) = physicalSbotMasses(1);
    physpars.md(2, family) = physicalSbotMasses(2);
  }
    
  /// start 3rd family computation: more complicated!
    family = 3;
    DoubleMatrix mSbotSquared(2, 2), mSbotSquared2(2, 2);
    treeDownSquark(mSbotSquared, mb, pizztMS, sinthDRbarMS, family);
    mSbotSquared2 = mSbotSquared; /// tree-level matrices

  if (accuracy > 0) {
    double pLight = minimum(forLoops.md(1, 3), forLoops.md(2, 3));
    double pHeavy = maximum(forLoops.md(1, 3), forLoops.md(2, 3));

    addSbotCorrection(pLight, mSbotSquared, mt);
    addSbotCorrection(pHeavy, mSbotSquared2, mt);
  }
    
  double theta;
  DoubleVector physicalSbotMassesSquared2(mSbotSquared2.sym2by2(theta));
  /// The ordering and the mixing angle is defined by the following matrix,
  /// with the LIGHT sbottom as the external momentum
  DoubleVector physicalSbotMassesSquared(mSbotSquared.sym2by2(theta));
  physpars.thetab = theta;

  /// Check for tachyonic sbottoms
  if (minimum(physicalSbotMassesSquared(1), physicalSbotMassesSquared(2))
      < 0.0 || minimum(physicalSbotMassesSquared2(1), 
		       physicalSbotMassesSquared2(2)) < 0.0) {
    flagTachyon(sbottom);
  }

  DoubleVector physicalSbotMasses(physicalSbotMassesSquared.apply(ccbSqrt));
  DoubleVector physicalSbotMasses2(physicalSbotMassesSquared2.apply(ccbSqrt));

  /// twisted measures the ordering of the sbottom masses. If msbot1 > msbot2,
  /// twisted is defined to be true (msbot2 > msbot1 is defined "untwisted").
  bool twisted = false;
  if (physicalSbotMassesSquared(1) > physicalSbotMassesSquared(2)) 
    twisted = true;

  double lightSbotMass = minimum(physicalSbotMasses(1), physicalSbotMasses(2));
  double heavySbotMass = maximum(physicalSbotMasses2(1), 
				 physicalSbotMasses2(2));

  if (twisted) {
    physpars.md(1, family) = heavySbotMass; 
    physpars.md(2, family) = lightSbotMass;
  } else {
    physpars.md(1, family) = lightSbotMass;
    physpars.md(2, family) = heavySbotMass;
  }
}
template<class SoftPars>
void Softsusy<SoftPars>::treeChargedSlepton(DoubleMatrix & mass, double mtaurun, 
                                            double /* pizztMS */, double sinthDRbarMS,
                                            int family) {
  double mz2 = sqr(displayMzRun()), mtau2 = sqr(mtaurun);
  double beta = atan(displayTanb()), mu = displaySusyMu(),
    c2b = cos(2. * beta), tanb = displayTanb(), sinth2 = sqr(sinthDRbarMS);

  mass(1, 1) = displaySoftMassSquared(mLl, family, family) -
    mz2 * (0.5 - sinth2) * c2b;
  mass(2, 2) = displaySoftMassSquared(mEr, family, family) -
    mz2 * sinth2 * c2b;

  if (family != 3) mass(1, 2) = 0.0;
  else {
    if (fabs(forLoops.htau) < EPSTOL)
      mass(1, 2) = mtaurun * (displaySoftA(EA, 3, 3) - mu * tanb);
    else 
      mass(1, 2) = mtaurun * (forLoops.utau / forLoops.htau - mu * tanb);

    mass(1, 1) = mass(1, 1) + mtau2;
    mass(2, 2) = mass(2, 2) + mtau2;
  }

  mass(2, 1) = mass(1, 2);
}

template<class SoftPars>
void Softsusy<SoftPars>::doChargedSleptons(double mtau, double pizztMS, double
			    sinthDRbarMS, int accuracy) {
  DoubleMatrix mSlepSquared(2, 2);

  /// no mixing assuming for first two families
  int family; for (family = 1; family <= 2; family++) {

    treeChargedSlepton(mSlepSquared, mtau, pizztMS, sinthDRbarMS, family);
    if (accuracy > 0) addSlepCorrection(mSlepSquared, family);

    if (mSlepSquared(1, 1) < 0.0 || mSlepSquared(2, 2) < 0.0) {
      if (family == 1) flagTachyon(selectron);
      else if (family == 2) flagTachyon(smuon);
    }
      
    physpars.me(1, family) = ccbSqrt(mSlepSquared(1, 1));
    physpars.me(2, family) = ccbSqrt(mSlepSquared(2, 2));
  }

  /// do third family
  family = 3; 
  treeChargedSlepton(mSlepSquared, mtau, pizztMS, sinthDRbarMS, family);
  DoubleMatrix mSlepSquared2(mSlepSquared);

  if (accuracy > 0) {
    double pLight = minimum(forLoops.me(1, 3), forLoops.me(2, 3));
    double pHeavy = maximum(forLoops.me(1, 3), forLoops.me(2, 3));

    addStauCorrection(pLight, mSlepSquared, mtau);
    addStauCorrection(pHeavy, mSlepSquared2, mtau);  
  }
   
  double theta;
  DoubleVector physicalStauMassesSquared2(mSlepSquared2.sym2by2(theta));
  DoubleVector physicalStauMassesSquared(mSlepSquared.sym2by2(theta));
  physpars.thetatau = theta; /// theta is calculated at p=mstau1

  if (physicalStauMassesSquared(1) < 0.0 || physicalStauMassesSquared(2) < 0.0
      ||  physicalStauMassesSquared2(1) < 0.0 || 
      physicalStauMassesSquared2(2) < 0.0) {
    flagTachyon(stau);
  }
  
  DoubleVector physicalStauMasses(physicalStauMassesSquared.apply(ccbSqrt));
  DoubleVector 
      physicalStauMasses2(physicalStauMassesSquared2.apply(ccbSqrt));
  
  double lightStauMass = minimum(physicalStauMasses(1), physicalStauMasses(2));
  double heavyStauMass = maximum(physicalStauMasses2(1), 
				 physicalStauMasses2(2));
  /// twisted measures the ordering of the sbottom masses. If mstop1 > mstop2,
  /// twisted is defined to be true (mstop2 > mstop1 is defined "untwisted").
  bool twisted = false;
  if (physicalStauMassesSquared(1) > physicalStauMassesSquared(2)) 
    twisted = true;

  if (twisted) {
    physpars.me(1, family) = heavyStauMass; 
    physpars.me(2, family) = lightStauMass;
  } else {
    physpars.me(1, family) = lightStauMass;
    physpars.me(2, family) = heavyStauMass;
  }
}

template<class SoftPars>
void Softsusy<SoftPars>::doSnu(double pizztMS, int accuracy) {
  double mSnuSquared;
  int family; for (family = 1; family <= 3; family++) {
    treeSnu(mSnuSquared, pizztMS, family);
    if (family == 3 && accuracy > 0) addSnuTauCorrection(mSnuSquared);
    else if (accuracy > 0) {
      addSnuCorrection(mSnuSquared, family);
    }

    physpars.msnu(family) = ccbSqrt(mSnuSquared);
  }
}

template<class SoftPars>
void Softsusy<SoftPars>::treeSnu(double & mSnuSquared,
                                 double /* pizztMS */, int family) {
  double mz2 = sqr(displayMzRun());
  double beta = atan(displayTanb());
  double c2b = cos(2.0 * beta);

  mSnuSquared = displaySoftMassSquared(mLl, family, family) + 0.5 * mz2
      * c2b;
}

/// Organises calculation of physical quantities such as sparticle masses etc
/// Call AT MSusy
template<class SoftPars>
void Softsusy<SoftPars>::physical(int accuracy) {
  double sinthDRbarMS, piwwtMS, pizztMS;

  calcDrBarPars();

  if (accuracy == 0) {
    sinthDRbarMS = 0.0;
    piwwtMS = 0.0;
    pizztMS = 0.0;
  }
  else {
    sinthDRbarMS = calcSinthdrbar();
    piwwtMS = sqr(displayMwRun()) - sqr(displayMw());
    pizztMS = sqr(displayMzRun()) - sqr(displayMz());
  }

  /// Running masses at MSUSY
  double mt = forLoops.mt;
  double mb = forLoops.mb;
  double mtau = forLoops.mtau;

  /// Re-calculate the 1-loop tadpoles for the calculation
  calcTadpole1Ms1loop(mt, sinthDRbarMS);  
  calcTadpole2Ms1loop(mt, sinthDRbarMS); 

  doUpSquarks(mt, pizztMS, sinthDRbarMS, accuracy); 
  doDownSquarks(mb, pizztMS, sinthDRbarMS, accuracy, mt); 
  doChargedSleptons(mtau, pizztMS, sinthDRbarMS, accuracy); 
  doSnu(pizztMS, accuracy);
  
  /// Charginos/neutralinos/higgs
  Softsusy<SoftPars> * ppp;
  ppp = this;
  ppp->higgs(accuracy, piwwtMS, pizztMS); 

  const int maxHiggsIterations = 20;
  double currentAccuracy = 1.0;
  DoubleVector oldHiggsMasses(4);
  oldHiggsMasses(1) = ppp->displayPhys().mh0(1);   
  oldHiggsMasses(2) = ppp->displayPhys().mA0(1);
  oldHiggsMasses(3) = ppp->displayPhys().mh0(2);
  oldHiggsMasses(4) = ppp->displayPhys().mHpm;
  bool higgsTachyon = false;
  /// Iterate Higgs calculation (unless accuracy=0, in which case we just need
  /// a rough calculation) until the Higgs masses all converge to better than
  /// TOLERANCE fractional accuracy
  int i = 1; while (i < maxHiggsIterations && accuracy > 0 && 
		    currentAccuracy > TOLERANCE) {
    higgsTachyon = ppp->higgs(accuracy, piwwtMS, pizztMS); /// iterate 

    DoubleVector newHiggsMasses(4);
    newHiggsMasses(1) = ppp->displayPhys().mh0(1);
    newHiggsMasses(2) = ppp->displayPhys().mA0(1);
    newHiggsMasses(3) = ppp->displayPhys().mh0(2);
    newHiggsMasses(4) = ppp->displayPhys().mHpm;

    currentAccuracy = oldHiggsMasses.compare(newHiggsMasses);

    oldHiggsMasses = newHiggsMasses;

    i++;
  }

  if (higgsTachyon) { flagTachyon(h0); flagTachyon(softsusy::A0); 
    flagTachyon(hpm); }
  physpars.mh0(1) = ppp->displayPhys().mh0(1);
  physpars.mA0(1) = ppp->displayPhys().mA0(1);
  physpars.mh0(2) = ppp->displayPhys().mh0(2);
  physpars.mHpm = ppp->displayPhys().mHpm;
 
  gluino(accuracy); 
  charginos(accuracy, piwwtMS); 
  neutralinos(accuracy, piwwtMS, pizztMS);
}

/// For a given trial value of the log of field H2, gives the value of the
/// potential at the minimum. The following global variables must be set before
/// it is called:
static double unificationScale, minTol;
inline double minimufb3(double lnH2) {

  /// Save initial parameters
  double initialMu = tempSoft1->displayMu();
  DoubleVector savePars(tempSoft1->display());

  MssmSoftsusy temp(*tempSoft1);

  double mx = unificationScale; 
  double mu = -temp.displaySusyMu();
  
  int family = 1; /// family could be 2 as well if they're very different
  
  /// field VEVS
  double htau = temp.displayYukawaElement(YE, 3, 3);
  double h2 =  - exp(lnH2), 
    eR = sqrt(fabs(mu * h2 / htau)),
    Lisq = -4.0 * tempSoft1->displaySoftMassSquared(mLl, family, family) / 
    (0.6 * sqr(tempSoft1->displayGaugeCoupling(1)) +
     sqr(tempSoft1->displayGaugeCoupling(2))) + sqr(h2) + sqr(eR);
  
  /// qhat is scale at which one-loop potential corrections are smallest:
  /// iterate
  double qhat = getQhat(minTol,eR, h2, Lisq, mx, temp);
  if (qhat < 0.0) {
    std::ostringstream ii;
    ii << "Error: qhat returned negative in ufb3min" << endl;
    throw ii.str();
  }
  
  double vufb3 = ufb3fn(mu, htau, h2, family, temp);
  if (PRINTOUT > 1) cout << vufb3 << endl;

  /// Restore initial parameters
  tempSoft1->setMu(initialMu);
  tempSoft1->set(savePars);

  return vufb3;  
}

/// Input mx the scale up to which you search for minima
/// Returns minimum value of potential along that direction
/// Does ufbs truly properly but takes ages.
template<>
inline double Softsusy<SoftParsMssm>::ufb3sl(double mx) {

  tempSoft1 = this;

  /// Save initial parameters
  DoubleVector parSave(display());

  unificationScale = mx;
  /// Running should be faster from high scales on the whole
  tempSoft1 -> runto(mx);
  
  double ax = log(1.0e3), bx = 0.5*(log(mx) + log(1.0e3)), cx = log(mx);
  
  double tol = 1.0e-3, lnx; minTol = tol * 10.0;
  /// Numerical recipes routine to determine minimum of potential specified in
  /// minimufb3
  double Vmin = findMinimum(ax, bx, cx, minimufb3, tol, &lnx);

  /// Restore initial parameters
  ///  setMu(initialMu);
  ///  set(parSave);
  
  return Vmin;
}

template<class SoftPars>
double Softsusy<SoftPars>::ufb3sl(double mx) {
   throw "ufb3sl not implemented for non-MSSM models";
}

/// Does SUSY (and other) threshold corrections to alphaS
/// Input alphas in MSbar and it returns it in DRbar scheme. 
/// From hep-ph/9606211
template<class SoftPars>
double Softsusy<SoftPars>::qcdSusythresh(double alphasMSbar, double q) const {
  drBarPars tree(displayDrBarPars());
  double mt = tree.mt;
  double deltaAlphas = alphasMSbar / (2.0 * PI) *
    (0.5 - 2.0 / 3.0 * log(mt / q) - 
     2.0 * log(fabs(tree.mGluino) / q));
  
  int i,j;
  for (i=1; i<=2; i++)
    for (j=1; j<=3; j++)
      deltaAlphas = deltaAlphas - alphasMSbar / (12.0 * PI) * 
	(log(tree.mu(i, j) / q) + 
	 log(tree.md(i, j) / q));
  return alphasMSbar / (1.0 - deltaAlphas); 
}

/// Does SUSY (and other) threshold corrections to alphaEm - returns alpha in
/// DRbar scheme at scale Q. From hep-ph/9606211. Input empirical value of
/// alpha at MZ external momentum....
template<class SoftPars>
double Softsusy<SoftPars>::qedSusythresh(double alphaEm, double q) const {

  drBarPars tree(displayDrBarPars());

  if (tree.mHpm < TOLERANCE) return 0.0;

  double deltaASM = -1.0 / 3.0 + 
    16.0 / 9.0 * log(tree.mt / q);

  double deltaASusy = 
    /// commented out since alpha(MZ) includes it!
    log(tree.mHpm / q) / 3.0 + 4.0 / 9.0 * 
     (log(tree.mu(1,1) / q) + 
      log(tree.mu(2,1) / q) + 
      log(tree.mu(1,2) / q) + 
      log(tree.mu(2,2) / q) + 
      log(tree.mu(1,3) / q) + 
      log(tree.mu(2,3) / q)) + 
     (log(tree.md(1,1) / q) + 
      log(tree.md(2,1) / q) + 
      log(tree.md(1,2) / q) + 
      log(tree.md(2,2) / q) + 
      log(tree.md(1,3) / q) + 
      log(tree.md(2,3) / q)) / 9.0 +
     (log(tree.me(1,1) / q) + 
      log(tree.me(2,1) / q) + 
      log(tree.me(1,2) / q) + 
      log(tree.me(2,2) / q) + 
      log(tree.me(1,3) / q) + 
      log(tree.me(2,3) / q)) / 3.0 
    + (log(fabs(tree.mch(1)) / q) 
       + log(fabs(tree.mch(2)) / q)) * 4.0 / 3.0;
  
  double deltaAlpha;
  deltaAlpha = -alphaEm / (2.0 * PI) * (deltaASM + deltaASusy);

  return alphaEm / (1.0 - deltaAlpha);
}

/// Returns lsp mass in mass and function return labels which particle is lsp:
/// 0 is neutralino posi = #, posj = 0
template<class SoftPars>
int Softsusy<SoftPars>::lsp(double & mass, int & posi, int & posj) const {
  int temp = 0, pos1 = 0, pos2 = 0;
  sPhysical s(displayPhys());
  
  double minmass = fabs(s.mneut(1)); posi = pos1; posj = 0;
  
  /// up squarks 1
  double lightest = s.mu.apply(fabs).min(pos1, pos2);
  if (lightest < minmass) { 
    minmass = lightest; posi = pos1; posj = pos2; temp = 1; }
  
  /// down squarks 2
  lightest = s.md.apply(fabs).min(pos1, pos2);
  if (lightest < minmass) { 
    minmass = lightest; posi = pos1; posj = pos2; temp = 2; }
  
  /// sleptons 3
  lightest = s.me.apply(fabs).min(pos1, pos2);
  if (lightest < minmass) { 
    minmass = lightest; posi = pos1; posj = pos2; temp = 3; }
  
  /// charginos 4
  lightest = fabs(s.mch(1));
  if (lightest < minmass) { 
    minmass = lightest; posi = pos1; posj = 0; temp = 4; }
  
  /// sneutrinos 5
  lightest = s.msnu.apply(fabs).min(pos1);
  if (lightest < minmass) { 
    minmass = lightest; posi = pos1; posj = 0; temp = 5; }
  
  /// gluino 6
  lightest = fabs(s.mGluino);
  if (lightest < minmass) { 
    minmass = lightest; posi = 0; posj = 0; temp = 6; }

  /// gravitino -1 
  lightest = displayGravitino();
  if (lightest < minmass) {
    minmass = lightest; posi = 0; posj = 0; temp = -1; }  

  mass = minmass;
  return temp;
}

template<class SoftPars>
double Softsusy<SoftPars>::maxMass() const {
  int pos1, pos2;
  sPhysical s(displayPhys());
  
  double maxmass = s.mneut.apply(fabs).min(pos1); 
  
  /// up squarks 1
  double heaviest = s.mu.apply(fabs).max(pos1, pos2);
  if (heaviest > maxmass) maxmass = heaviest; 
  
  /// down squarks 2
  heaviest = s.md.apply(fabs).max(pos1, pos2);
  if (heaviest > maxmass) maxmass = heaviest; 
  
  /// sleptons 3
  heaviest = s.me.apply(fabs).max(pos1, pos2);
  if (heaviest > maxmass) maxmass = heaviest; 
  
  /// charginos 4
  heaviest = s.mch.apply(fabs).max();
  if (heaviest > maxmass) maxmass = heaviest; 
  
  /// sneutrinos 5
  heaviest = s.msnu.apply(fabs).max();
  if (heaviest > maxmass) maxmass = heaviest; 
  
  /// gluino 6
  heaviest = s.mGluino;
  if (heaviest > maxmass) maxmass = heaviest; 
  
  return maxmass;
}

/// DRbar pars should be defined already for this
template<class SoftPars>
double Softsusy<SoftPars>::calcMs() const {

  drBarPars tree(displayDrBarPars());

  if (QEWSB < EPSTOL) throw("QEWSB Probably too low\n");

  double msf1, msf2;

  if (alternativeMs) {
     const DoubleVector sups(tree.mu.flatten());
     msf1 = sups.min();
     msf2 = sups.max();
  } else {
     msf1 = tree.mu(1,3);
     msf2 = tree.mu(2,3);
  }

  if (QEWSB < MZ) 
    return maximum(QEWSB * sqrt(msf1 * msf2), displayMz());
  else return QEWSB;
}

/// Provides the first guess at a SUSY object at mt, inputting tanb and oneset
/// (should be at MZ) - it's very crude, doesn't take radiative corrections
/// into account etc. 
template<class SoftPars>
MssmSusy Softsusy<SoftPars>::guessAtSusyMt(double tanb, const QedQcd_legacy & oneset) {
  /// This bit gives a guess at a SUSY object
  QedQcd_legacy leAtMt(oneset);

  DoubleVector a(3), g(3);
  double sinth2 = 1.0 - sqr(MW / MZ);
  /// Gauge couplings at mt
  a = leAtMt.getGaugeMu(oneset.displayPoleMt(), sinth2);
  
  MssmSusy t; 
  t.setTanb(tanb);
  double beta = atan(tanb);
  
  /// Yukawa couplings -- at roughly tree level
  double vev = 246.22;
  double ht = (leAtMt.displayMass(legacy::mTop) - 30.0) * root2 
    / (vev * sin(beta));
  double hb =  ht * tanb * leAtMt.displayMass(legacy::mBottom) /
    (leAtMt.displayMass(legacy::mTop) - 30.0);
  double htau =  hb * leAtMt.displayMass(legacy::mTau) /
    leAtMt.displayMass(legacy::mBottom); 
  t.setYukawaElement(YU, 3, 3, ht);
  t.setYukawaElement(YD, 3, 3, hb);
  t.setYukawaElement(YE, 3, 3, htau);

  double hc = ht * leAtMt.displayMass(legacy::mCharm) / 
    (leAtMt.displayMass(legacy::mTop) - 30.0);
  double hs = hb * leAtMt.displayMass(legacy::mStrange) / 
    leAtMt.displayMass(legacy::mBottom);
  double hmu = htau * leAtMt.displayMass(legacy::mMuon) / 
    leAtMt.displayMass(legacy::mTau);
  double hu = ht * leAtMt.displayMass(legacy::mUp) / 
    (leAtMt.displayMass(legacy::mTop) - 30.0);
  double hd = hb * leAtMt.displayMass(legacy::mDown) / 
    leAtMt.displayMass(legacy::mBottom);
  double he = htau * leAtMt.displayMass(legacy::mElectron) / 
    leAtMt.displayMass(legacy::mTau);
  t.setYukawaElement(YU, 2, 2, hc);
  t.setYukawaElement(YD, 2, 2, hs);
  t.setYukawaElement(YE, 2, 2, hmu);
  t.setYukawaElement(YU, 1, 1, hu);
  t.setYukawaElement(YD, 1, 1, hd);
  t.setYukawaElement(YE, 1, 1, he);

  const double pi = 3.141592653589793;
  int i; for(i=1; i<=3; i++) g(i) = sqrt(4.0 * pi * a(i));
  t.setAllGauge(g);
  
  t.setHvev(vev);

  t.setMu(oneset.displayPoleMt());   
  
  return t;
}

/// Returns low energy softsusy object consistent with BC's m0 etc at MGUT.
/// oneset should be at MZ and contains the SM data to fit the model to.
/// If the running comes into difficulty, eg if a Landau pole is reached, it
/// returns a ZERO object: no result is possible!
/// Boundary condition is the theoretical condition on parameters at the high
/// energy scale mx: the parameters themselves are contained within the vector.
template<class SoftPars>
void Softsusy<SoftPars>::fixedPointIteration
(void (*boundaryCondition)(Softsusy<SoftPars> &, const DoubleVector &),
 double mxGuess, 
 const DoubleVector & pars, int sgnMu, double tanb, const QedQcd_legacy &
 oneset, bool gaugeUnification, bool ewsbBCscale) {

  try {

    const static Softsusy<SoftPars> empty;
  
    double muFirst = displaySusyMu(); /// Remember initial values

    bool setTbAtMXflag = displaySetTbAtMX(); 
    bool altFlag = displayAltEwsb();
    double m32 = displayGravitino();
    double muCondFirst = displayMuCond();
    double maCondFirst = displayMaCond();
    bool altMs = alternativeMs;

    setSoftsusy(empty); /// Always starts from an empty object
    /// These are things that are re-written by the new initialisation
    setSetTbAtMX(setTbAtMXflag); 
    if (altFlag) useAlternativeEwsb();
    setData(oneset); 
    setMw(MW); 
    setM32(m32);
    setMuCond(muCondFirst);
    setMaCond(maCondFirst);
    setAlternativeMs(altMs);

    double mz = displayMz();

    /// Here all was same
    if (mxGuess > 0.0) 
      mxBC = mxGuess; 
    else {
      std::string ii("Trying to use negative mx in Softsusy<SoftPars>::fixedPointIteration.\n");
      ii = ii + "Now illegal! Use positive mx for first guess of mx.\n";
      throw ii;
    }
    
    if (oneset.displayMu() != mz) {
      cout << "WARNING: fixedPointIteration in softsusy.cpp called with oneset at scale\n" 
	   << oneset.displayMu() << "\ninstead of " << mz << endl;
    }
    
    int maxtries = 100; 
    double tol = TOLERANCE;
    
    MssmSusy t(guessAtSusyMt(tanb, oneset));

    // default SoftSusy loop number
    int lpnum = 2;

    if (USE_THREE_LOOP_RGE) lpnum = 3;

    t.setLoops(lpnum); /// >= 2 loops should protect against ht Landau pole 
    t.runto(mxBC); 
   
    setSusy(t);
    /// Initial guess: B=0, 
    boundaryCondition(*this, pars);

    if ((sgnMu == 1 || sgnMu == -1) && !ewsbBCscale) {
      setSusyMu(sgnMu * 1.0);
      setM3Squared(0.);
    }
    else {
      if (altEwsb) {
	setSusyMu(displayMuCond());
	setM3Squared(displayMaCond());
      }
      else {
	setSusyMu(muFirst);
	setM3Squared(muFirst); 
      }
    }

    run(mxBC, mz);
    
    if (sgnMu == 1 || sgnMu == -1) rewsbTreeLevel(sgnMu); 
  
    physical(0);
  
    setThresholds(3); setLoops(lpnum);
    
    itLowsoft(maxtries, sgnMu, tol, tanb, boundaryCondition, pars, 
		gaugeUnification, ewsbBCscale);
    
    if (displayProblem().nonperturbative 
	|| displayProblem().higgsUfb || displayProblem().tachyon 
	|| displayProblem().noRhoConvergence)
      return;
    
    runto(maximum(displayMsusy(), mz));
    if (ewsbBCscale) boundaryCondition(*this, pars); 

    physical(3);

    runto(mz);
    
    if (PRINTOUT > 1) cout << " end of iteration" << endl;
  }
    catch(const char *a) {
    std::ostringstream ii;
    ii << "SOFTSUSY problem: " << a << " pars=" << pars << " tanb=" << tanb 
       << " oneset=" << oneset << endl;
    flagProblemThrown(true);
    throw(ii.str());
  }
  catch(const std::string & a) {
    std::ostringstream ii;
    ii << "SOFTSUSY problem: " << a << " pars=" << pars << " tanb=" << tanb 
	 << " oneset=" << oneset << endl;
    flagProblemThrown(true);
    throw ii.str();
  }
  catch(...) {
    std::ostringstream ii;
    ii << "SOFTSUSY problem: " << endl;
    ii << "pars=" << pars << " tanb=" << tanb
       << " oneset=" << oneset << endl;
    flagProblemThrown(true);
    throw ii.str();
  }
}


/// You should evaluate this at a scale MSusy average of stops.
/// Depth of electroweak minimum: hep-ph/9507294. Bug-fixed 19/11/04
template<class SoftPars>
double Softsusy<SoftPars>::realMinMs() const {
  MssmSusy temp(displaySusy());
  temp.runto(calcMs(), TOLERANCE); 
  
  double beta = atan(temp.displayTanb());
  return - 1.0 / 32.0 * 
    (0.6 * sqr(temp.displayGaugeCoupling(1)) + 
     sqr(temp.displayGaugeCoupling(2))) * 
    sqr(sqr(displayHvev()) * cos(2.0 * beta));
}

/// Calculates sin theta at the current scale
template<class SoftPars>
double Softsusy<SoftPars>::calcSinthdrbar() const {

  double sinth = displayGaugeCoupling(1) /
    sqrt(sqr(displayGaugeCoupling(1)) +
	 sqr(displayGaugeCoupling(2)) * 5. / 3.);

  return sinth;   
}

//VEV at current scale, using an input value of Z self-energy
template<class SoftPars>
double Softsusy<SoftPars>::getVev(double pizzt) {

  double vsquared = 4.0 * (sqr(displayMz()) + pizzt) /
    (sqr(displayGaugeCoupling(2)) +
     sqr(displayGaugeCoupling(1)) * 0.6); 

  if (vsquared < 0.0 || testNan(vsquared)) {
    flagTachyon(Z);
    return 246.22;
  }

  return sqrt(vsquared);
}

//VEV at current scale: calculates Z self energy first
template<class SoftPars>
double Softsusy<SoftPars>::getVev() {
  double pizzt = piZZT(displayMz(), displayMu());
  if (pizzt + sqr(displayMz()) < 0.) {
    pizzt = -displayMz() + EPSTOL;
    flagTachyon(Z);
  }

  return getVev(pizzt);
}

/// It'll set the important SUSY couplings: supposed to be applied at MZ
/// You should set up an iteration here since Yuk's depend on top mass which
/// depends on Yuk's etc. 
template<class SoftPars>
void Softsusy<SoftPars>::sparticleThresholdCorrections(double tb) {
  double mz = displayMz();  if (displayMu() != mz) {
    std::ostringstream ii;
    ii << "Called Softsusy<SoftPars>::sparticleThresholdCorrections "
       << "with scale" << displayMu() << endl;
    throw ii.str();
  }
  
  if (!setTbAtMX) setTanb(tb);
  calcDrBarPars(); /// for the up-coming loops
  double alphaMsbar = dataSet.displayAlpha(legacy::ALPHA);
  double alphaDrbar = qedSusythresh(alphaMsbar, displayMu());

  double alphasMZDRbar =
    qcdSusythresh(displayDataSet().displayAlpha(legacy::ALPHAS), displayMu());
  
  /// Do gauge couplings
  double outrho = 1.0, outsin = 0.48, tol = TOLERANCE * 1.0e-8; 
  int maxTries = 20;
  double pizztMZ = piZZT(displayMz(), displayMu(), true);
  double piwwt0  = piWWT(0., displayMu(), true);
  double piwwtMW = piWWT(displayMw(), displayMu(), true);
  
  if (piwwt0 + sqr(displayMw()) < 0.) {
    flagTachyon(W);
    piwwt0 = -sqr(displayMw()) + EPSTOL;
  }
  if (piwwtMW + sqr(displayMw()) < 0.) {
    flagTachyon(W);
    piwwtMW = -sqr(displayMw()) + EPSTOL;
  }
  if (pizztMZ + sqr(displayMz()) < 0.) {
    flagTachyon(Z);
    pizztMZ = -sqr(displayMz()) + EPSTOL;
  }

  rhohat(outrho, outsin, alphaDrbar, pizztMZ, piwwt0, piwwtMW, tol, maxTries);

  //  if (problem.noRhoConvergence) 
  //    outsin = sqrt(1.0 - sqr(displayMw() / displayMz())); 
  
  double eDR = sqrt(4.0 * PI * alphaDrbar), costhDR = cos(asin(outsin));

  DoubleVector newGauge(3);
  newGauge(1) = eDR / costhDR * sqrt(5.0 / 3.0);
  newGauge(2) = eDR / outsin;
  newGauge(3) = sqrt(4.0 * PI * alphasMZDRbar);

  double vev = getVev();

  double beta = atan(displayTanb());
  DoubleMatrix mUq(3, 3), mDq(3, 3), mLep(3, 3);
  massFermions(displayDataSet(), mDq, mUq, mLep);

  /// replace third family quark masses with their loop-corrected values
  mDq(3, 3)  = calcRunningMb();
  mUq(3, 3)  = calcRunningMt();
  mLep(3, 3) = calcRunningMtau();
  /// to do: for a better approximation, the lighter quarks/leptons should
  /// also be corrected with MSSM/DRbar corrections

  /// 3-family mixed-up Yukawa couplings: From PDG 2000
  doQuarkMixing(mDq, mUq); 

  if (MIXING == -1) {
    mDq(1, 1) = 0.; mDq(2, 2) = 0.; mUq(1, 1) = 0.; mUq(2, 2) = 0.;
    mLep(1, 1) = 0.; mLep(2, 2) = 0.;
  }

  double poleMwSq = 0.25 * sqr(newGauge(2)) * sqr(vev) - 
    piWWT(displayMw(), displayMu());
  if (poleMwSq < 0.) flagTachyon(W);
  setMw(ccbSqrt(poleMwSq)); 
  setGaugeCoupling(1, newGauge(1));
  setGaugeCoupling(2, newGauge(2));
  setGaugeCoupling(3, newGauge(3));
  setHvev(vev); 
  setYukawaMatrix(YU, mUq * (root2 / (vev * sin(beta))));
  setYukawaMatrix(YD, mDq * (root2 / (vev * cos(beta)))); 
  setYukawaMatrix(YE, mLep * (root2 / (vev * cos(beta)))); 
}

template<class SoftPars>
double Softsusy<SoftPars>::displayMzRun() const { 
  return displayHvev() * 0.5 * 
    sqrt(sqr(displayGaugeCoupling(2)) + 0.6 * sqr(displayGaugeCoupling(1)));
} 

template<class SoftPars>
void Softsusy<SoftPars>::treeCharginos(DoubleMatrix & mCh, double beta, double mw) {
  mCh(1, 1) = displayGaugino(2);
  mCh(2, 1) = root2 * mw * cos(beta); 
  mCh(1, 2) = mCh(2, 1) * displayTanb();
  mCh(2, 2) = displaySusyMu();
}

template<class SoftPars>
void Softsusy<SoftPars>::treeNeutralinos(DoubleMatrix & mNeut, double beta, double mz, double mw, double sinthDRbar) {
  
  mNeut(1, 1) = displayGaugino(1);
  mNeut(2, 2) = displayGaugino(2);
  mNeut(1, 3) = - mz * cos(beta) * sinthDRbar;
  mNeut(1, 4) = - mNeut(1, 3) * displayTanb();
  mNeut(2, 3) = mw * cos(beta);
  mNeut(2, 4) = - mNeut(2, 3) * displayTanb();
  mNeut(3, 4) = - displaySusyMu();
  mNeut.symmetrise();
  
}

template<class SoftPars>
void Softsusy<SoftPars>::calcDrBarGauginos(double beta, double mw, double mz, double sinth, drBarPars & eg) {
  DoubleMatrix mCh(2, 2);   
  treeCharginos(mCh, beta, mw);
  eg.mch = mCh.asy2by2(eg.thetaL, eg.thetaR);
  eg.mpzCharginos();
  DoubleMatrix mNeut(4, 4);
  treeNeutralinos(mNeut, beta, mz, mw, sinth);
  if (mNeut.diagonaliseSym(eg.mixNeut, eg.mneut) > TOLERANCE *
      1.0e-3) { 
    std::ostringstream ii;
    ii << "accuracy bad in neutralino diagonalisation"<< std::flush;
    ii << mNeut;
    throw ii.str(); 
    }

  eg.mpzNeutralinos();
}
template<class SoftPars>
void Softsusy<SoftPars>::calcDrBarHiggs(double beta, double mz2, double mw2, 
                                        double /* sinthDRbar */, drBarPars & eg) {
   if (eg.mt > 200. || eg.mt < 50.) {
    /// Gone badly off-track
    flagProblemThrown(true);
    if (eg.mt > 200.) eg.mt = 200.;
    if (eg.mt < 50.) eg.mt = 50.;
  }
  /// You could instead do like sPHENO, choose what you'd get from minimising
  /// the potential at tree level, ie (mH2^2-mH1^2)/cos(2 beta)-mz^2. This
  /// *may* be less sensitive to becoming a tachyon at MZ. 
  double mAsq;
  if (sphenoMassConv)
    mAsq = (displayMh2Squared() - displayMh1Squared()) 
        / (cos(2. * beta)) - mz2; 
  else mAsq = displayM3Squared() / (sin(beta) * cos(beta));

  if (mAsq < 0.) {
    /* Previous solution: if we're at MZ, use the pole mA^2
       if (close(displayMu(), MZ, tol)) {
      double mApole = physpars.mA0(1); /// physical value
      setDrBarPars(eg);
      
      double piaa = piAA(mApole, displayMu()); 
      double t1Ov1 = doCalcTadpole1oneLoop(eg.mt, sinthDRbar), 
      t2Ov2 = doCalcTadpole2oneLoop(eg.mt, sinthDRbar); 
      double poleMasq = 
      (displayMh2Squared() - t2Ov2 - 
      displayMh1Squared() + t1Ov1) / 
      cos(2.0 * beta) - mz2 - piaa +
      sqr(sin(beta)) * t1Ov1 + sqr(cos(beta)) * t2Ov2;
      
      mAsq = poleMasq;
      
      if (mAsq < 0.) { flagTachyon(A0); mAsq = fabs(poleMasq); }
      }
     */
    flagTachyon(softsusy::A0); 
    if (mAFlag == false) mAsq = ccbSqrt(mAsq); 
    /// This may be  a bad idea in terms of convergence
    else mAsq = fabs(mAsq);
    
    if (PRINTOUT > 1) cout << " mA^2(tree)=" << mAsq << " since m3sq=" 
			   << displayM3Squared() << " @ "<< displayMu() 
			   << " " << endl; 
  }
    
  DoubleMatrix mH(2, 2); 
  mH(1, 1) = mAsq * sqr(sin(beta)) + mz2 * sqr(cos(beta));
  mH(1, 2) = - sin(beta) * cos(beta) * (mAsq + mz2); 
  mH(2, 2) = mAsq * sqr(cos(beta)) + mz2 * sqr(sin(beta)); 
  mH(2, 1) = mH(1 ,2); 

  DoubleVector mSq(2);
  mSq = mH.sym2by2(eg.thetaH);
  if (mSq(1) < 0. || mSq(2) < 0.) {
    flagTachyon(h0);
  }
  DoubleVector temp(mSq.apply(ccbSqrt));
  if (temp(2) > temp(1)) eg.thetaH = eg.thetaH + PI * 0.5; 

  int pos;
  eg.mh0(1) = temp.min(pos); eg.mh0(2) = temp.max(); 
  eg.mA0(1) = sqrt(mAsq); eg.mHpm = sqrt(mAsq + mw2);  
}

//PA: sets the neutral current couplings
template<class SoftPars>
void Softsusy<SoftPars>::setNeutCurrCouplings(double sinthDRbar, double & sw2, double & guL, double & gdL, double & geL, double & guR, double & gdR, double & geR ) {
  sw2 = sqr(sinthDRbar); 
  guL = 0.5 - 2.0 * sw2 / 3.0;
  gdL = -0.5 + sw2 / 3.0;
  geL = -0.5 + sw2;
  guR = 2.0 * sw2 / 3.0;
  gdR = -sw2 / 3.0;
  geR = -sw2;
}
//PA: sets the Yukawas and Trilinears
template<class SoftPars>
void Softsusy<SoftPars>::calcDRTrilinears(drBarPars & eg, double vev, double beta) {
   if (MIXING > 0) {    
    DoubleMatrix diagUp(displayYukawaMatrix(YU)),
      diagDown(displayYukawaMatrix(YD)),
      diagLep(displayYukawaMatrix(YE));
    
    DoubleMatrix diagTriUp(displayTrilinear(UA)),
      diagTriDown(displayTrilinear(DA)),
      diagTriLep(displayTrilinear(EA));

    DoubleMatrix u(3, 3), v(3, 3); 
    DoubleVector w(3);

    diagUp.diagonalise(u, v, w);
    diagTriUp = u.transpose() * diagTriUp * v;
    eg.ht = w(3);
    eg.ut = diagTriUp(3, 3);

    diagDown.diagonalise(u, v, w);
    diagTriDown = u.transpose() * diagTriDown * v;
    eg.hb = w(3);
    eg.ub = diagTriDown(3, 3);
    
    diagLep.diagonalise(u, v, w);
    diagTriLep = u.transpose() * diagTriLep * v;
    eg.htau = w(3);
    eg.utau = diagTriLep(3, 3);
  } else {
    eg.ht   = displayYukawaElement(YU, 3, 3);
    eg.hb   = displayYukawaElement(YD, 3, 3);
    eg.htau = displayYukawaElement(YE, 3, 3);
    eg.ut   = displayTrilinear(UA, 3, 3);
    eg.ub   = displayTrilinear(DA, 3, 3);
    eg.utau = displayTrilinear(EA, 3, 3);
  }
 
  eg.mt   = eg.ht   * vev * sin(beta) / root2;
  eg.mb   = eg.hb   * vev * cos(beta) / root2;
  eg.mtau = eg.htau * vev * cos(beta) / root2; 

  forLoops.ht = eg.ht; forLoops.mt = eg.mt; forLoops.ut = eg.ut;
  forLoops.hb = eg.hb; forLoops.mb = eg.mb; forLoops.ub = eg.ub;
  forLoops.htau = eg.htau; forLoops.mtau = eg.mtau; forLoops.utau = eg.utau;
 }

/// calculates masses all at tree-level in the DRbar scheme, useful for
/// radiative corrections. 
template<class SoftPars>
void Softsusy<SoftPars>::calcDrBarPars() {
  /// We want to set mu to be the one obtained from the tree-level Higgs
  /// potential for these purposes
  double savedMu = displaySusyMu();
  if (sphenoMassConv) {
     double muSq = treeLevelMuSq();
     double muForNow = ccbSqrt(muSq);
     setSusyMu(muForNow);
  }

  drBarPars eg(displayDrBarPars());
  /// First, must define mstop,sbot,stau and mixing angles in DRbar scheme
  double beta = atan(displayTanb()), mzPole = displayMz();
  double sinthDRbar = calcSinthdrbar();
  double mz = displayMzRun(), mz2 = sqr(mz);
  double pizzt = sqr(mz) - sqr(mzPole);

  setNeutCurrCouplings(sinthDRbar, sw2, guL, gdL, geL, guR, gdR, geR);
  
  double vev = displayHvev();
  calcDRTrilinears(eg, vev, beta);
  
  eg.mGluino = displayGaugino(3);
  
  DoubleVector mSq(2);
  int family; for(family = 1; family <= 3; family++) {
    
    DoubleMatrix mSquared(2, 2); 
    treeUpSquark(mSquared, eg.mt, pizzt, sinthDRbar, family);
    mSq = mSquared.sym2by2(eg.thetat);
    if (mSq(1) < 0. || mSq(2) < 0.) {
      switch(family) {
      case 1: flagTachyon(sup); break;
      case 2: flagTachyon(scharm); break;
      case 3: flagTachyon(stop); break;	
      default: throw("Bad family number in calcDrBarPars\n");
      }
      if (PRINTOUT > 2) cout << " tree sup(" << family << ") tachyon ";
    }
    DoubleVector mstopDRbar(mSq.apply(ccbSqrt));   
    
    treeDownSquark(mSquared, eg.mb, pizzt, sinthDRbar, family);
    mSq = mSquared.sym2by2(eg.thetab);
    if (mSq(1) < 0. || mSq(2) < 0.) {
      switch(family) {
      case 1: flagTachyon(sdown); break;
      case 2: flagTachyon(sstrange); break;
      case 3: flagTachyon(sbottom); break;	
      default: throw("Bad family number in calcDrBarPars\n");
      }
    if (PRINTOUT > 1) cout << " tree sdown(" << family << ") tachyon ";
    }
    DoubleVector msbotDRbar(mSq.apply(ccbSqrt));   
    
    treeChargedSlepton(mSquared, eg.mtau, pizzt, sinthDRbar, family);
    mSq = mSquared.sym2by2(eg.thetatau);
    if (mSq(1) < 0. || mSq(2) < 0.) {
      switch(family) {
      case 1: flagTachyon(selectron); break;
      case 2: flagTachyon(smuon); break;
      case 3: flagTachyon(stau); break;	
      default: throw("Bad family number in calcDrBarPars\n");
      }
    if (PRINTOUT > 1) cout << " tree selectron(" << family << ") tachyon ";
    }
    DoubleVector mstauDRbar(mSq.apply(ccbSqrt));   
    
    int i; for (i=1; i<=2; i++) {
      eg.mu(i, family) = mstopDRbar(i);    eg.md(i, family) = msbotDRbar(i);
      eg.me(i, family) = mstauDRbar(i);    
    }
    double mSnuSquared;
    treeSnu(mSnuSquared, pizzt, family);
    if (mSnuSquared < 0.) {
      switch(family) {
      case 1: flagTachyon(snue); break;
      case 2: flagTachyon(snumu); break;
      case 3: flagTachyon(snutau); break;	
      default: throw("Bad family number in calcDrBarPars\n");
      }
    if (PRINTOUT > 1) cout << " tree sneutrino(" << family << ") tachyon@"
			   << displayMu() << " ";
    }
    eg.msnu(family) = ccbSqrt(mSnuSquared);
  }

  double mw = displayMwRun();
  double mw2 = sqr(mw);
  calcDrBarGauginos(beta, mw, mz, sinthDRbar, eg);
  eg.mw = mw;
  eg.mz = mz; 
  calcDrBarHiggs(beta, mz2, mw2, sinthDRbar, eg); 
  setDrBarPars(eg);

  /// Restore the proper loop corrected value for mu
  setSusyMu(savedMu);
}

template<class SoftPars>
void Softsusy<SoftPars>::itLowsoft
(int maxTries, int sgnMu, double tol, double tanb, 
 void (*boundaryCondition)(Softsusy<SoftPars> &, const DoubleVector &), 
 const DoubleVector & pars, bool gaugeUnification, bool ewsbBCscale) {

  static Softsusy<SoftPars> old;
  static double oldMu = 0.;
  static int numTries = 0;
  double mz = displayMz();

  if (numTries != 0 && sqr(displayMu() / mz - 1.0) > TOLERANCE) {
    cout << "WARNING: itLowsoft called at inappropriate";
    cout << " scale:" << displayMu() << endl; 
    cout << "whereas it should be " << mz << endl; 
  }
  
  if (numTries - 1 > maxTries) {/// Iterating too long: bail out
    flagNoConvergence(true);
    if (PRINTOUT) cout << "itLowsoft reached maxtries\n"; 
    numTries = 0; 
    return;
  }

  if (PRINTOUT > 1) cout << displayProblem(); 

  double mtrun;
  
  /// On first iteration, don't bother with finite corrections  
  numTries = numTries + 1;
  try {
    sparticleThresholdCorrections(tanb); 

    if (problem.noRhoConvergence && PRINTOUT) 
      cout << "No convergence in rhohat\n"; 

    
    /// precision of running/RGE integration: start off low and increase
    double eps = maximum(exp(double(- numTries) * log(10.0)), tol * 0.01); 
    
    /// first stab at MSUSY: root(mstop1(MZ) mstop2(MZ))
    if (numTries == 1) setMsusy(calcMs()); 

    int err = 0;

    err = runto(displayMsusy(), eps);
    double tbIn; double predictedMzSq = 0.;
    predictedMzSq = predMzsq(tbIn);
    setPredMzSq(predictedMzSq);  
    if (!ewsbBCscale) err = runto(mxBC, eps);


    /// Guard against the top Yukawa fixed point
    if (displayYukawaElement(YU, 3, 3) > 3.0 
	|| displayYukawaElement(YD, 3, 3) > 3.0 
	|| displayYukawaElement(YE, 3, 3) > 3.0) {
      setYukawaElement(YU, 3, 3, minimum(3.0, displayYukawaElement(YU, 3, 3)));
      setYukawaElement(YD, 3, 3, minimum(3.0, displayYukawaElement(YD, 3, 3)));
      setYukawaElement(YE, 3, 3, minimum(3.0, displayYukawaElement(YE, 3, 3)));
      flagIrqfp(true); 
    }
    
    if (err) {
      /// problem with running: bail out 
      flagProblemThrown(true);
      if (PRINTOUT) 
	cout << "itLowsoft can't run more approaching mgut\n"; 
      if (PRINTOUT > 1) printObj();
      numTries = 0; 
      return;
    }
  
    if (gaugeUnification) {
      sBrevity (dummy);
      MssmSusy a(this -> MssmSusy::beta(dummy));
      
      /// Equal gauge couplings: let them and their derivatives set the boundary
      /// condition scale -- linear approximation
      mxBC = mxBC * exp((displayGaugeCoupling(2) - displayGaugeCoupling(1))
		    / (a.displayGaugeCoupling(1) - a.displayGaugeCoupling(2)));

      /// if mx is too high/low, will likely get non-perturbative problems
      if (mxBC < 1.0e4) {
	mxBC = 1.0e4;
	if (PRINTOUT > 2) cout << " mxBC too low ";
	flagMgutOutOfBounds(true);
      }
      if (mxBC > 5.0e17) {
	if (PRINTOUT > 2) cout << " mxBC =" << mxBC <<" too high ";
	mxBC = 5.0e17;
	flagMgutOutOfBounds(true);
      }
    }
    
    boundaryCondition(*this, pars); 
    
    if (!ewsbBCscale) err = runto(displayMsusy(), eps);
    calcDrBarPars();

    if (err) {
      // problem with running: bail out 
      flagProblemThrown(true);
      if (PRINTOUT) cout << "itLowsoft gone non-perturbative on way to MZ\n"; 
      if (PRINTOUT > 1) printObj();
      numTries = 0;
      return;
    }

    setMsusy(calcMs());
    if (ewsbBCscale) mxBC = displayMsusy();
    if (PRINTOUT > 0) cout << " mgut=" << mxBC << std::flush;
    
    mtrun = forLoops.mt; ///< This will be at MSUSY
    //    double tbIn; double predictedMzSq = 0.;
    if (numTries < 11) {
      rewsb(sgnMu, mtrun, pars);    
    }
    else { ///< After 11 tries, we start averaging old/new mu values
      double epsi = 0.5;
      if (numTries > 20) epsi = 0.2;
      if (numTries > 30) epsi = 0.1;
      rewsb(sgnMu, mtrun, pars, oldMu, epsi);    
      } 

    oldMu = displaySusyMu();

    fracDiff = sumTol(*this, old, numTries);    
    
    if (numTries !=0 && fracDiff < tol) {///< Accuracy achieved: bail out
      numTries = 0; ///< Reset the number of iterations for the next time
      if (PRINTOUT > 1) cout << " sT=" << fracDiff << " " << std::flush; 
      if (displayProblem().test() && PRINTOUT > 0) 
	cout << " ***problem point***: " << displayProblem() << ".";

      return; 
    } else if (numTries != 0 && fabs(displayM3Squared()) > 1.0e20) {
      /// guards against nasty fatal problems 
      numTries = 0; flagM3sq(true); return;
    }

    // All problems should be reset since only the ones of the final iteration
    // should count (sometimes problems disappear). This can mean that problems
    // only show up as no rho convergence....
    tachyonType nil(none);
    flagAllProblems(false, nil);

    /// Old iteration is 'old': these are the parameters by which convergence is
    /// measured. 
    old.setDrBarPars(displayDrBarPars());
    /// If a print out is desired, print respectively, the difference with the
    /// last iteration (sum tol or sT), the mu parameter and m3^2 from EWSB, and
    /// the predicted MW and MZ boson masses
    if (PRINTOUT > 1) {
      cout << "\n" << numTries << ". sT=" << fracDiff << " mu=" 
	   << displaySusyMu() <<  " m3sq=" <<
	displayM3Squared() << " MWp=" << displayMw() << " Mzp=" 
	   << sqrt(displayPredMzSq()) << std::flush;
    }

    if (problem.noMuConvergence) {
      if (PRINTOUT) 
	cout << "itLowsoft doesn't break EWSB\n"; 
      if (PRINTOUT > 1) printObj();
    }
    
    err = runto(mz, eps);
    if (err) {
      /// problem with running: bail out 
      flagProblemThrown(true);
      if (PRINTOUT) cout << "itLowsoft gone non-perturbative on way to MZ\n"; 
      if (PRINTOUT > 1) printObj();
      ///    old = Softsusy<SoftPars>();
      numTries = 0;
      return;
    }
    
    itLowsoft(maxTries, sgnMu, tol, tanb, boundaryCondition, pars, 
	      gaugeUnification, ewsbBCscale);
  }
  catch(const char *a) {
    numTries = 0;
    throw a;
  }
  catch(const std::string &a) {
    numTries = 0;
    throw a;
  }
}


/// Higgs contribution to the Transverse part of Z self-energy
template<class SoftPars>
double Softsusy<SoftPars>::piZZTHiggs(double p, double q, 
				      double thetaWDRbar) const {
  double    alpha   = displayDrBarPars().thetaH ;
  double    beta    = atan(displayTanb());
  double smHiggs = 0.0, susyHiggs = 0.0;
  double    mz      = displayMzRun();
  double    cw2DRbar    = sqr(cos(thetaWDRbar));
  double    sw2DRbar    = 1.0 - cw2DRbar;

  smHiggs = 
    - sqr(sin(alpha - beta)) *
    (b22bar(p, mz, displayDrBarPars().mh0(1), q) - 
     sqr(mz) * b0(p, mz, displayDrBarPars().mh0(1), q));

  susyHiggs = - sqr(sin(alpha - beta)) *
    b22bar(p, displayDrBarPars().mA0(1), displayDrBarPars().mh0(2), q);
 
  susyHiggs = susyHiggs
    - sqr(cos(alpha - beta)) * 
    (b22bar(p, mz, displayDrBarPars().mh0(2), q) +
     b22bar(p, displayDrBarPars().mA0(1), displayDrBarPars().mh0(1), q) -
     sqr(mz) * b0(p, mz, displayDrBarPars().mh0(2), q));
  
  smHiggs = smHiggs
    - 2.0 * sqr(cw2DRbar) * (2 * sqr(p) + sqr(displayMwRun()) - sqr(mz) *
			     sqr(sw2DRbar) / cw2DRbar)
    * b0(p, displayMwRun(), displayMwRun(), q)
    - (8.0 * sqr(cw2DRbar) + sqr(cos(2.0 * thetaWDRbar))) * 
    b22bar(p, displayMwRun(), displayMwRun(), q);

  susyHiggs = susyHiggs - 
    sqr(cos(2.0 * thetaWDRbar)) 
    * b22bar(p, displayDrBarPars().mHpm, displayDrBarPars().mHpm, q);
  
  double higgs = smHiggs + susyHiggs;

  return higgs;
}

/// sfermion contribution to the Transverse part of Z self-energy
template<class SoftPars>
double Softsusy<SoftPars>::piZZTsfermions(double p, double q) const {
  double    thetat = displayDrBarPars().thetat ;
  double    thetab = displayDrBarPars().thetab;
  double    thetatau= displayDrBarPars().thetatau ;
  double    st      = sin(thetat) ;
  double    sb      = sin(thetab) ;
  double    stau    = sin(thetatau) ;
  double    ct      = cos(thetat) ;
  double    cb      = cos(thetab) ;
  double    ctau    = cos(thetatau);
  double    ct2     = sqr(ct) ;
  double    cb2     = sqr(cb) ;
  double    ctau2   = sqr(ctau) ;
  double    st2     = (1.0 - ct2);
  double    sb2     = (1.0 - cb2);
  double    stau2   = (1.0 - ctau2);
  double  squarks = 0.0, thirdFamily = 0.0, sneutrinos = 0.0, sleptons = 0.0;
  //static 
  DoubleVector vu(2), vd(2), ve(2);
  //static 
  DoubleMatrix vt(2, 2), vb(2, 2), vtau(2, 2);
    
  vu(1) = guL;
  vu(2) = guR;
  
  vt(1, 1) = guL * ct2 - guR * st2;
  vt(1, 2) = (guL + guR) * ct * st;
  vt(2, 2) = guR * ct2 - guL * st2;
  vt(2, 1) = vt(1, 2);
  
  vd(1) = gdL;
  vd(2) = gdR;
  
  vb(1, 1) = gdL * cb2 - gdR * sb2;
  vb(1, 2) = (gdL + gdR) * cb * sb;
  vb(2, 2) = gdR * cb2 - gdL * sb2;
  vb(2, 1) = vb(1, 2);
  
  ve(1) = geL;
  ve(2) = geR;
  
  vtau(1, 1) = geL * ctau2 - geR * stau2;
  vtau(1, 2) = (geL + geR) * ctau * stau;
  vtau(2, 2) = geR * ctau2 - geL * stau2;
  vtau(2, 1) = vtau(1, 2);
  
  /// first two families of sfermions
  int i, j, family;
  for (family = 1; family<=2; family++) {      
    /// sneutrinos
    sneutrinos = sneutrinos -  
      b22bar(p, displayDrBarPars().msnu(family), displayDrBarPars().msnu(family), q);
    for (i=1; i<=2; i++) {
      /// up squarks
      squarks = squarks - 12.0 * sqr(vu(i)) * 
	b22bar(p, displayDrBarPars().mu(i, family), 
	       displayDrBarPars().mu(i, family), q); 
      /// down squarks
      squarks = squarks - 12.0 * sqr(vd(i)) * 
	b22bar(p, displayDrBarPars().md(i, family), 
	       displayDrBarPars().md(i, family), q);
      /// sleptons
      sleptons = sleptons - 4.0 * sqr(ve(i)) * 
	b22bar(p, displayDrBarPars().me(i, family), 
	       displayDrBarPars().me(i, family), q); 
    }
  }

  family = 3;
  
  /// THIRD FAMILY
  /// sneutrinos
  thirdFamily = thirdFamily - 
    b22bar(p, displayDrBarPars().msnu(family), 
	   displayDrBarPars().msnu(family), q);
  for (i=1; i<=2; i++)
    for (j=1; j<=2; j++) {
      /// up squarks
      thirdFamily = thirdFamily - 12.0 * sqr(vt(i, j)) * 
	b22bar(p, displayDrBarPars().mu(i, family), 
	       displayDrBarPars().mu(j, family), q); 
      
      /// down squarks
      thirdFamily = thirdFamily - 12.0 * sqr(vb(i, j)) * 
	b22bar(p, displayDrBarPars().md(i, family), 
	       displayDrBarPars().md(j, family), q);
      /// selectrons
       thirdFamily = thirdFamily - 4.0 * sqr(vtau(i, j)) * 
	b22bar(p, displayDrBarPars().me(i, family), 
	       displayDrBarPars().me(j, family), q); 
    }
  
  double sfermions = squarks + sleptons + sneutrinos  + thirdFamily;

 return sfermions;
}

/// fermion contribution to the Transverse part of Z self-energy: 
template<class SoftPars>
double Softsusy<SoftPars>::piZZTfermions(double p, double q, bool usePoleMt) const {
  /// fermions: these are valid at MZ
  double    mtop =  displayDrBarPars().mt;
  /// We utilise pole mt for these corrections (which the 2-loop Standard Model
  /// pieces assume)
  if (usePoleMt) mtop = displayDataSet().displayPoleMt();

  double    mb   =  displayDrBarPars().mb;
  double    mtau =  displayDrBarPars().mtau;
  double    ms   =  displayDataSet().displayMass(legacy::mStrange) ;
  double    mc   =  displayDataSet().displayMass(legacy::mCharm) ;
  double    mmu  =  displayDataSet().displayMass(legacy::mMuon) ;
  double    mE  =   displayDataSet().displayMass(legacy::mElectron) ;
  double    mD  =   displayDataSet().displayMass(legacy::mDown) ;
  double    mU  =   displayDataSet().displayMass(legacy::mUp);

  double quarks = 0.0;

  quarks = quarks + 3.0 * hfn(p, mU, mU, q) * (sqr(guL) + sqr(guR)); 
  quarks = quarks + 3.0 * hfn(p, mc, mc, q) * (sqr(guL) + sqr(guR)); 
  quarks = quarks + 3.0 * (hfn(p, mtop, mtop, q) * 
		     (sqr(guL) + sqr(guR)) - 4.0 * guL * guR * sqr(mtop) * 
		     b0(p, mtop, mtop, q));

  quarks = quarks + 3.0 * hfn(p, mD, mD, q) * (sqr(gdL) + sqr(gdR)); 
  quarks = quarks + 3.0 * hfn(p, ms, ms, q) * (sqr(gdL) + sqr(gdR)); 
  quarks = quarks + 3.0 * (hfn(p, mb, mb, q) * (sqr(gdL) + sqr(gdR)) 
			   - 4.0 * gdL * gdR * sqr(mb) * b0(p, mb, mb, q)); 

  quarks = quarks + hfn(p, mE, mE, q) * (sqr(geL) + sqr(geR)); 
  quarks = quarks + hfn(p, mmu, mmu, q) * (sqr(geL) + sqr(geR)); 
  quarks = quarks + hfn(p, mtau, mtau, q) * 
    (sqr(geL) + sqr(geR)) - 4.0 * geL * geR 
    * sqr(mtau) * b0(p, mtau, mtau, q); 

  quarks = quarks + 3.0 * hfn(p, 0., 0., q) * 0.25;
  
  return quarks;
}

template<class SoftPars>
double Softsusy<SoftPars>::piZZTNeutralinos(double p, double q, 
					    double thetaWDRbar) const {
  
  double    cw2DRbar    = sqr(cos(thetaWDRbar));
  double    g       = displayGaugeCoupling(2);
  /// Neutralinos
  //static 
  double neutralinos = 0.0;
  ComplexMatrix aPsi(4, 4), bPsi(4, 4), aChi(4, 4), bChi(4, 4);
  ComplexMatrix n(displayDrBarPars().nBpmz);
  DoubleVector mneut(displayDrBarPars().mnBpmz);

  aPsi(3, 3) = g / (2.0 * cos(thetaWDRbar)); aPsi(4, 4) = -1. * aPsi(3, 3);
  bPsi = -1. * aPsi;
  
  aChi = n.complexConjugate() * aPsi * n.transpose();
  bChi = n * bPsi * n.hermitianConjugate();
  
  for (int i=1; i<=4; i++)
    for (int j=1; j<=4; j++) {
      neutralinos = neutralinos + cw2DRbar / (2.0 * sqr(g)) * 
	((sqr(aChi(i, j).mod()) + sqr(bChi(i, j).mod())) * 
	 hfn(p, mneut(i), mneut(j), q)
	 + 4.0 * (bChi(i, j).conj() * aChi(i, j)).real() *
	 mneut(i) * mneut(j) * b0(p, mneut(i), mneut(j), q)); 
    }
  
  return neutralinos;
}

template<class SoftPars>
double Softsusy<SoftPars>::piZZTCharginos(double p, double q, double thetaWDRbar) const {
 
  double    cw2DRbar    = sqr(cos(thetaWDRbar));
  double    g       = displayGaugeCoupling(2);
  double  charginos = 0.0;

  /// Charginos
  ///  static 
  ComplexMatrix u(displayDrBarPars().uBpmz), v(displayDrBarPars().vBpmz); 
  DoubleVector mch(displayDrBarPars().mchBpmz); 
  ComplexMatrix aPsiCh(2, 2), aCh(2, 2), bCh(2, 2);
  aPsiCh(1, 1) = g * cos(thetaWDRbar);
  aPsiCh(2, 2) = g * cos(2.0 * thetaWDRbar) / (2.0 * cos(thetaWDRbar));

  aCh = v.complexConjugate() * aPsiCh * v.transpose();
  bCh = u * aPsiCh * u.hermitianConjugate();
  
  for (int i=1; i<=2; i++)
    for(int j=1; j<=2; j++) {	
      charginos = charginos + cw2DRbar / sqr(g) * 
	((sqr(aCh(i, j).mod()) + sqr(bCh(i, j).mod())) * 
	 hfn(p, mch(i), mch(j), q) 
	 + 4.0 * (bCh(i, j).conj() * aCh(i, j)).real() * mch(i) *
	 mch(j) * b0(p, mch(i), mch(j), q));
    }

  return charginos;

}
/// Transverse part of Z self-energy: has been checked
template<class SoftPars>
double Softsusy<SoftPars>::piZZT(double p, double q, bool usePoleMt) const {
  
  double    thetaWDRbar = asin(calcSinthdrbar());
  double    cw2DRbar    = sqr(cos(thetaWDRbar));
  double    g       = displayGaugeCoupling(2);
  double rhs = 0.0;
 
  //PA: obtain Higgs contributions in separate method
  double higgs = piZZTHiggs(p, q, thetaWDRbar);
  //PA: obtain sfermion contributions in separate method
  double sfermions = piZZTsfermions(p, q);
  //PA: obtain fermion contributions in separate method
  double fermions = piZZTfermions(p, q, usePoleMt);
   //PA: obtain neutralino contributions in separate method
  double neutralinos = piZZTNeutralinos(p, q, thetaWDRbar);
   //PA: obtain neutralino contributions in separate method
  double charginos = piZZTCharginos(p, q, thetaWDRbar);
 
  rhs = higgs + charginos + neutralinos + fermions + sfermions ;
  double pi = rhs * sqr(g) / (cw2DRbar * 16.0 * sqr(PI));

  return pi;
}


template<class SoftPars>
double Softsusy<SoftPars>::piWWTHiggs(double p, double q, double thetaWDRbar) const {
  double    beta      = atan(displayTanb());
  double    alpha     = displayDrBarPars().thetaH ;
  double    cw2DRbar  = sqr(cos(thetaWDRbar));
  double    sw2DRbar  = 1.0 - cw2DRbar;
  double    mH = displayDrBarPars().mh0(2); 
  double    mh0 = displayDrBarPars().mh0(1);
  double    mHc = displayDrBarPars().mHpm;
  double    mA = displayDrBarPars().mA0(1);
  double smHiggs = 0.0, susyHiggs = 0.;
  double    mz      = displayMzRun();
  smHiggs = - sqr(sin(alpha - beta)) * 
    (b22bar(p, mh0, displayMwRun(), q) 
     - sqr(displayMwRun()) * b0(p, mh0, displayMwRun(), q));

  susyHiggs = - sqr(sin(alpha - beta)) * b22bar(p, mH, mHc, q);

  susyHiggs = susyHiggs -
    sqr(cos(alpha - beta)) *
    (b22bar(p, mh0, mHc, q) + b22bar(p, mH, displayMwRun(), q) 
     - sqr(displayMwRun()) * b0(p, mH, displayMwRun(), q)) -
    b22bar(p, mA, mHc, q);

  smHiggs = smHiggs  - 
    (1.0 + 8.0 * cw2DRbar) * b22bar(p, mz, displayMwRun(), q)
    - sw2DRbar * (8.0 * b22bar(p, displayMwRun(), 0.0, q) + 4.0 * sqr(p) * 
		  b0(p, displayMwRun(), 0.0, q))  
    - ((4.0 * sqr(p) + sqr(mz) + sqr(displayMwRun())) * cw2DRbar - sqr(mz)  *
       sqr(sw2DRbar)) * b0(p, mz, displayMwRun(), q);

  double higgs = smHiggs + susyHiggs;

  return higgs;
}

template<class SoftPars>
double Softsusy<SoftPars>::piWWTfermions(double p, double q, bool usePoleMt) const {
 double    mtop =  displayDrBarPars().mt;
  /// We utilise pole mt for these corrections (which the 2-loop Standard Model
  /// pieces assume)
  if (usePoleMt) mtop = displayDataSet().displayPoleMt();
  double    mb   =  displayDrBarPars().mb;
  double    mtau =  displayDrBarPars().mtau;
 
  /// fermions: these are valid at MZ
  double    ms   =  displayDataSet().displayMass(legacy::mStrange) ;
  double    mc   =  displayDataSet().displayMass(legacy::mCharm) ;
  double    mmu  =  displayDataSet().displayMass(legacy::mMuon) ;
  double    mE  =  displayDataSet().displayMass(legacy::mElectron) ;
  double    mD  =  displayDataSet().displayMass(legacy::mDown) ;
  double    mU  =  displayDataSet().displayMass(legacy::mUp);
  
 double fermions =
    1.5 * (hfn(p, mU, mD, q) + hfn(p, mc, ms, q) + hfn(p, mtop, mb, q)) + 0.5
    * (hfn(p, 0.0, mE, q) + hfn(p, 0.0, mmu, q) + hfn(p, 0.0, mtau, q));     

 return fermions;
}

template<class SoftPars>
double Softsusy<SoftPars>::piWWTsfermions(double p, double q) const {
  
  double    thetat = displayDrBarPars().thetat ;
  double    thetab = displayDrBarPars().thetab;
  double    thetatau= displayDrBarPars().thetatau ;
  double    st      = sin(thetat) ;
  double    sb      = sin(thetab) ;
  double    stau    = sin(thetatau) ;
  double    ct      = cos(thetat) ;
  double    cb      = cos(thetab) ;
  double    ctau    = cos(thetatau);

  /// sfermions
  double sfermions = - 6.0 *
    (b22bar(p, displayDrBarPars().mu(1, 1), displayDrBarPars().md(1, 1), q) +
     b22bar(p, displayDrBarPars().mu(1, 2), displayDrBarPars().md(1, 2), q)) -
    2.0 * (b22bar(p, displayDrBarPars().msnu(1), displayDrBarPars().me(1, 1), q) +
	   b22bar(p, displayDrBarPars().msnu(2), displayDrBarPars().me(1, 2), q));
  
  /// stop/sbottom
  DoubleMatrix w(2, 2);
  double stopBot = 0.0;
  w(1, 1) = ct * cb; w(1, 2) = ct * sb; w(2, 1) = st * cb; w(2, 2) = st * sb;
  int i, j; 
  for(i=1; i<=2; i++)
    for(j=1; j<=2; j++)
      stopBot = stopBot - 6.0 * sqr(w(i, j)) * 
	b22bar(p, displayDrBarPars().mu(i, 3), displayDrBarPars().md(j, 3), q);
  
  /// LH slepton
  double slepton = - 2.0 * 
    (sqr(ctau) * 
     b22bar(p, displayDrBarPars().msnu(3), displayDrBarPars().me(1, 3), q)
     + sqr(stau) * 
     b22bar(p, displayDrBarPars().msnu(3), displayDrBarPars().me(2, 3), q));
  
  sfermions += stopBot + slepton; 
  
  return sfermions;
}

template<class SoftPars>
double Softsusy<SoftPars>::piWWTgauginos(double p, double q,
                                         double /* thetaWDRbar */) const {
  double    g       = displayGaugeCoupling(2);
  ComplexMatrix aPsi0PsicW(4, 2), bPsi0PsicW(4, 2), aChi0ChicW(4, 2),
    bChi0ChicW(4, 2);
  DoubleMatrix fW(4, 2), gW(4, 2);
  
  aPsi0PsicW(2, 1) = - g;
  bPsi0PsicW(2, 1) = - g;
  aPsi0PsicW(4, 2) = g / root2;		     
  bPsi0PsicW(3, 2) = -g / root2;		     
  
  ComplexMatrix aPsi(4, 4), bPsi(4, 4), aChi(4, 4), bChi(4, 4);
  ComplexMatrix n(displayDrBarPars().nBpmz);
  DoubleVector mneut(displayDrBarPars().mnBpmz);
  ComplexMatrix u(displayDrBarPars().uBpmz), v(displayDrBarPars().vBpmz); 
  DoubleVector mch(displayDrBarPars().mchBpmz); 

  /// These ought to be in physpars
  aChi0ChicW = n.complexConjugate() * aPsi0PsicW * v.transpose();
  bChi0ChicW = n * bPsi0PsicW * u.hermitianConjugate();

  double gauginos = 0.0;

  for(int i=1;i<=4;i++)
    for(int j=1;j<=2;j++) {
      fW(i, j) = sqr(aChi0ChicW(i, j).mod()) + sqr(bChi0ChicW(i, j).mod());
      gW(i, j) = 2.0 * (bChi0ChicW(i, j).conj() * aChi0ChicW(i, j)).real(); 
      gauginos = gauginos + 
	(fW(i, j) * hfn(p, mneut(i), mch(j), q)
	 + 2.0 * gW(i, j) * mneut(i) * mch(j) * b0(p, mneut(i), mch(j), q)) 
	/ sqr(g);
    }
return gauginos;
}

/// W propagator to 1 loop in MSSM 
/// It's all been checked
template<class SoftPars>
double Softsusy<SoftPars>::piWWT(double p, double q, bool usePoleMt) const {

  double    thetaWDRbar = asin(calcSinthdrbar());
  double    g       = displayGaugeCoupling(2);

  double ans = 0.0;
  double higgs = piWWTHiggs(p, q, thetaWDRbar);
  double fermions = piWWTfermions(p, q, usePoleMt);   
  double sfermions = piWWTsfermions(p, q);   
  double gauginos = piWWTgauginos(p, q, thetaWDRbar);
  ans = higgs + sfermions + fermions + gauginos;

  double pi = ans * sqr(g) / (16.0 * sqr(PI));

  return pi;
}



template<class SoftPars>
double Softsusy<SoftPars>::pis1s1Sfermions(double p, double q,  DoubleMatrix ls1tt,  DoubleMatrix ls1bb,  DoubleMatrix ls1tautau) const {
  drBarPars tree(displayDrBarPars()); 
  double hb = tree.hb, htau = tree.htau;
  double    thetat  = tree.thetat ;
  double    thetab  = tree.thetab;
  double    thetatau= tree.thetatau;
  double    st      = sin(thetat) ;
  double    sb      = sin(thetab) ;
  double    stau    = sin(thetatau);
  double    ct      = cos(thetat) ;
  double    cb      = cos(thetab) ;
  double    ctau    = cos(thetatau);
  double    msbot1  = tree.md(1, 3);
  double    msbot2  = tree.md(2, 3);
  double    mstau1  = tree.me(1, 3);
  double    mstau2  = tree.me(2, 3);
  double    mstop1  = tree.mu(1, 3);
  double    mstop2  = tree.mu(2, 3);
  double    mz      = displayMzRun();
  double    g       = displayGaugeCoupling(2);
  double    thetaWDRbar = asin(calcSinthdrbar());
  double    costhDrbar  = cos(thetaWDRbar);
  double    cw2DRbar    = sqr(cos(thetaWDRbar));
  double    cosb        = cos(atan(displayTanb()));

  double sbots = 3.0 * sqr(hb) * (a0(msbot1, q) + a0(msbot2, q));

  double staus =  sqr(htau) *
    (a0(mstau1, q) + a0(mstau2, q));

  double stops = 3.0 * sqr(g) / (2.0 * cw2DRbar) *
    (guL * (sqr(ct) * a0(mstop1, q) + sqr(st) * a0(mstop2, q)) +
     guR * (sqr(st) * a0(mstop1, q) + sqr(ct) * a0(mstop2, q)));

  sbots = sbots + 3.0 * sqr(g) / (2.0 * cw2DRbar) *
    (gdL * (sqr(cb) * a0(msbot1, q) + sqr(sb) * a0(msbot2, q)) +
     gdR * (sqr(sb) * a0(msbot1, q) + sqr(cb) * a0(msbot2, q)));

  staus = staus + sqr(g) / (2.0 * cw2DRbar) *
    (geL * (sqr(ctau) * a0(mstau1, q) + sqr(stau) * a0(mstau2, q)) +
     geR * (sqr(stau) * a0(mstau1, q) + sqr(ctau) * a0(mstau2, q)));

  double sups = 0.0, sdowns = 0.0, sneutrinos = 0.0, sleps = 0.0;
  int fam; for(fam=1; fam<=2; fam++) {
    sups = sups + 3.0 * sqr(g) / (2.0 * cw2DRbar) *
      (guL * a0(tree.mu(1, fam), q) + 
       guR * a0(tree.mu(2, fam), q));
    sdowns = sdowns + 3.0 * sqr(g) / (2.0 * cw2DRbar) *
      (gdL * a0(tree.md(1, fam), q) + 
       gdR * a0(tree.md(2, fam), q));
    sleps = sleps + sqr(g) / (2.0 * cw2DRbar) * 
      (geL * a0(tree.me(1, fam), q) + 
       geR * a0(tree.me(2, fam), q));
  }
  

/// Mix 3rd family up
  ls1tt = rot2d(thetat) * ls1tt * rot2d(-thetat);
  ls1bb = rot2d(thetab) * ls1bb * rot2d(-thetab);
  ls1tautau = rot2d(thetatau) * ls1tautau * rot2d(-thetatau);

  int i, j; for (i=1; i<=2; i++) {
    for (j=1; j<=2; j++) {
      stops = stops + 3.0 * sqr(ls1tt(i, j)) * 
	b0(p, tree.mu(i, 3), tree.mu(j, 3), q);
      sbots = sbots + 3.0 * sqr(ls1bb(i, j)) * 
	b0(p, tree.md(i, 3), tree.md(j, 3), q);
      staus = staus +  sqr(ls1tautau(i, j)) * 
	b0(p, tree.me(i, 3), tree.me(j, 3), q);
    }}

  /// selectron couplings to s1 Higgs state: neglect Yukawas + mixing
  double ls1eeLL = g * mz * geL * cosb / costhDrbar;
  double ls1eeRR = g * mz * geR * cosb / costhDrbar;
  double ls1ddLL = g * mz * gdL * cosb / costhDrbar;
  double ls1ddRR = g * mz * gdR * cosb / costhDrbar;
  double ls1uuLL = g * mz * guL * cosb / costhDrbar;
  double ls1uuRR = g * mz * guR * cosb / costhDrbar;

  int k; 
  for (k=1; k<=2; k++) {
    sups = sups + 
      3.0 * sqr(ls1uuLL) * 
      b0(p, tree.mu(1, k), tree.mu(1, k), q) +
      + 3.0 * sqr(ls1uuRR) * 
      b0(p, tree.mu(2, k), tree.mu(2, k), q);
    sdowns = sdowns + 
      3.0 * sqr(ls1ddLL) * 
      b0(p, tree.md(1, k), tree.md(1, k), q) +
      + 3.0 * sqr(ls1ddRR) * 
      b0(p, tree.md(2, k), tree.md(2, k), q);
    sleps = sleps + 
      sqr(ls1eeLL) * 
      b0(p, tree.me(1, k), tree.me(1, k), q) +
      + sqr(ls1eeRR) * 
      b0(p, tree.me(2, k), tree.me(2, k), q);
    }  

  sneutrinos = sneutrinos + 
    sqr(g) / (2.0 * cw2DRbar) * gnuL * 
    (a0(tree.msnu(1), q) + a0(tree.msnu(2), q) + 
     a0(tree.msnu(3), q)) +
    sqr(g * mz / sqrt(cw2DRbar) * gnuL * cosb) * 
    (b0(p, tree.msnu(1), tree.msnu(1), q) +
     b0(p, tree.msnu(2), tree.msnu(2), q) +
     b0(p, tree.msnu(3), tree.msnu(3), q));

  double sfermions = sups  + sdowns  + sleps  + stops  + sbots  + staus  + sneutrinos;
   return sfermions;
}

template<class SoftPars>
double Softsusy<SoftPars>::pis1s2Sfermions(double p, double q,  DoubleMatrix ls1tt,  DoubleMatrix ls1bb,  DoubleMatrix ls1tautau, DoubleMatrix ls2tt,  DoubleMatrix ls2bb,  DoubleMatrix ls2tautau) const {
drBarPars tree(displayDrBarPars()); 
  double    thetat  = tree.thetat ;
  double    thetab  = tree.thetab;
  double    thetatau= tree.thetatau;
  double    mz      = displayMzRun();
  double    g       = displayGaugeCoupling(2);
  double    thetaWDRbar = asin(calcSinthdrbar());
  double    costhDrbar  = cos(thetaWDRbar);
  double cosb = cos(atan(displayTanb())), sinb = sin(atan(displayTanb()));

   /// Mix 3rd family up
  ls1tt = rot2d(thetat) * ls1tt * rot2d(-thetat);
  ls1bb = rot2d(thetab) * ls1bb * rot2d(-thetab);
  ls1tautau = rot2d(thetatau) * ls1tautau * rot2d(-thetatau);
 
  ls2tt = rot2d(thetat) * ls2tt * rot2d(-thetat);
  ls2bb = rot2d(thetab) * ls2bb * rot2d(-thetab);
  ls2tautau = rot2d(thetatau) * ls2tautau * rot2d(-thetatau);
  
  double sfermions = 0.0;
  int i, j; for (i=1; i<=2; i++)
    for (j=1; j<=2; j++) {
      sfermions = sfermions + 3.0 * ls1tt(i, j) * ls2tt(i, j) *
	b0(p, tree.mu(i, 3), tree.mu(j, 3), q);
      sfermions = sfermions + 3.0 * ls1bb(i, j) * ls2bb(i, j) * 
	b0(p, tree.md(i, 3), tree.md(j, 3), q);
      sfermions = sfermions + ls1tautau(i, j) * ls2tautau(i, j) * 
	b0(p, tree.me(i, 3), tree.me(j, 3), q);
    }
    
  /// sneutrinos
  double l1 = g * mz / costhDrbar * cosb * gnuL, 
    l2 = -g * mz / costhDrbar * gnuL * sinb; 
  sfermions = sfermions +
    l1 * l2 * (b0(p, tree.msnu(1), tree.msnu(1), q) +
	       b0(p, tree.msnu(2), tree.msnu(2), q) +
	       b0(p, tree.msnu(3), tree.msnu(3), q));

  /// selectron couplings to s1 Higgs state: neglect Yukawas + mixing
  double ls1uuLL = g * mz * guL * cosb / costhDrbar;
  double ls1uuRR = g * mz * guR * cosb / costhDrbar;
  double ls1eeLL = g * mz * geL * cosb / costhDrbar;
  double ls1eeRR = g * mz * geR * cosb / costhDrbar;
  double ls1ddLL = g * mz * gdL * cosb / costhDrbar;
  double ls1ddRR = g * mz * gdR * cosb / costhDrbar;
  /// couplings to s2 Higgs state: neglect Yukawas + mixing
  double ls2uuLL = -g * mz * guL * sinb / costhDrbar;
  double ls2uuRR = -g * mz * guR * sinb / costhDrbar;
  double ls2eeLL = -g * mz * geL * sinb / costhDrbar;
  double ls2eeRR = -g * mz * geR * sinb / costhDrbar;
  double ls2ddLL = -g * mz * gdL * sinb / costhDrbar;
  double ls2ddRR = -g * mz * gdR * sinb / costhDrbar;

  int k; 
  for (k=1; k<=2; k++) {
    sfermions = sfermions + 3.0 * ls1uuLL * ls2uuLL *
      b0(p, tree.mu(1, k), tree.mu(1, k), q) +
      + 3.0 * ls1uuRR * ls2uuRR *
      b0(p, tree.mu(2, k), tree.mu(2, k), q);
    sfermions = sfermions + 3.0 * ls1ddLL * ls2ddLL *
      b0(p, tree.md(1, k), tree.md(1, k), q) +
      + 3.0 * ls1ddRR * ls2ddRR *
      b0(p, tree.md(2, k), tree.md(2, k), q);
    sfermions = sfermions + ls1eeLL * ls2eeLL * 
      b0(p, tree.me(1, k), tree.me(1, k), q) +
      + ls1eeRR * ls2eeRR * 
      b0(p, tree.me(2, k), tree.me(2, k), q);
    }  

  return sfermions;
}

template<class SoftPars>
double Softsusy<SoftPars>::pis2s2Sfermions(double p, double q, DoubleMatrix ls2tt,  DoubleMatrix ls2bb,  DoubleMatrix ls2tautau) const {
  drBarPars tree(displayDrBarPars()); 
  double    thetat  = tree.thetat ;
  double    thetab  = tree.thetab;
  double    thetatau= tree.thetatau;
  double    st      = sin(thetat) ;
  double    sb      = sin(thetab) ;
  double    stau    = sin(thetatau);
  double    ct      = cos(thetat) ;
  double    cb      = cos(thetab) ;
  double    ctau    = cos(thetatau);
  double    msbot1  = tree.md(1, 3);
  double    msbot2  = tree.md(2, 3);
  double    mstau1  = tree.me(1, 3);
  double    mstau2  = tree.me(2, 3);
  double    mstop1  = tree.mu(1, 3);
  double    mstop2  = tree.mu(2, 3);
  double    mz      = displayMzRun();
  double    g       = displayGaugeCoupling(2);
  double    thetaWDRbar = asin(calcSinthdrbar());
  double    costhDrbar  = cos(thetaWDRbar);
  double    cw2DRbar    = sqr(cos(thetaWDRbar));
  double  ht = tree.ht;
  double sinb = sin(atan(displayTanb()));


   /// stop contribution
  double sfermions = 3.0 * sqr(ht) * (a0(mstop1, q) + a0(mstop2, q));
  sfermions = sfermions - 3.0 * sqr(g) / (2.0 * cw2DRbar) *
    (guL * (sqr(ct) * a0(mstop1, q) + sqr(st) * a0(mstop2, q)) +
     guR * (sqr(st) * a0(mstop1, q) + sqr(ct) * a0(mstop2, q)));
 
  /// sbottom contribution
  sfermions = sfermions - 3.0 * sqr(g) / (2.0 * cw2DRbar) *
    (gdL * (sqr(cb) * a0(msbot1, q) + sqr(sb) * a0(msbot2, q)) +
     gdR * (sqr(sb) * a0(msbot1, q) + sqr(cb) * a0(msbot2, q)));

  //stau
  sfermions = sfermions - sqr(g) / (2.0 * cw2DRbar) *
    (geL * (sqr(ctau) * a0(mstau1, q) + sqr(stau) * a0(mstau2, q)) +
     geR * (sqr(stau) * a0(mstau1, q) + sqr(ctau) * a0(mstau2, q)));

  /// first two families of sparticles
  int fam; for(fam=1; fam<=2; fam++) 
    sfermions = sfermions - 3.0 * sqr(g) / (2.0 * cw2DRbar) *
      (guL * a0(tree.mu(1, fam), q) + 
       guR * a0(tree.mu(2, fam), q) +
       gdL * a0(tree.md(1, fam), q) + 
       gdR * a0(tree.md(2, fam), q)) - sqr(g) / (2.0 * cw2DRbar) * 
      (geL * a0(tree.me(1, fam), q) + 
       geR * a0(tree.me(2, fam), q)) -
      sqr(g)  / (2.0 * cw2DRbar) * gnuL *
      (a0(tree.msnu(fam), q));
  sfermions = sfermions -
     sqr(g)  / (2.0 * cw2DRbar) * gnuL *
      (a0(tree.msnu(3), q));

  /// Mix 3rd family up
  ls2tt = rot2d(thetat) * ls2tt * rot2d(-thetat);
  ls2bb = rot2d(thetab) * ls2bb * rot2d(-thetab);
  ls2tautau = rot2d(thetatau) * ls2tautau * rot2d(-thetatau);

  for (int i=1; i<=2; i++)
    for (int j=1; j<=2; j++) {
      /// stop 
      sfermions = sfermions + 3.0 * sqr(ls2tt(i, j)) * 
	b0(p, tree.mu(i, 3), tree.mu(j, 3), q);
      /// sbottom
      sfermions = sfermions + 3.0 * sqr(ls2bb(i, j)) * 
	b0(p, tree.md(i, 3), tree.md(j, 3), q);
      /// stay
      sfermions = sfermions + sqr(ls2tautau(i, j)) * 
	b0(p, tree.me(i, 3), tree.me(j, 3), q);
    }

  /// couplings to s2 Higgs state: neglect Yukawas + mixing
  double ls2uuLL = -g * mz * guL * sinb / costhDrbar;
  double ls2uuRR = -g * mz * guR * sinb / costhDrbar;
  double ls2eeLL = -g * mz * geL * sinb / costhDrbar;
  double ls2eeRR = -g * mz * geR * sinb / costhDrbar;
  double ls2ddLL = -g * mz * gdL * sinb / costhDrbar;
  double ls2ddRR = -g * mz * gdR * sinb / costhDrbar;

  int k; 
  for (k=1; k<=2; k++) {
    sfermions = sfermions + 3.0 * sqr(ls2uuLL) * 
      b0(p, tree.mu(1, k), tree.mu(1, k), q) +
      + 3.0 * sqr(ls2uuRR) * 
      b0(p, tree.mu(2, k), tree.mu(2, k), q);
    sfermions = sfermions + 3.0 * sqr(ls2ddLL) * 
      b0(p, tree.md(1, k), tree.md(1, k), q) +
      + 3.0 * sqr(ls2ddRR) * 
      b0(p, tree.md(2, k), tree.md(2, k), q);
    sfermions = sfermions + sqr(ls2eeLL) * 
      b0(p, tree.me(1, k), tree.me(1, k), q) +
      + sqr(ls2eeRR) * 
      b0(p, tree.me(2, k), tree.me(2, k), q);
    }  

  sfermions = sfermions +
    sqr(g * mz * gnuL * sinb) / cw2DRbar *
    (b0(p, tree.msnu(1), tree.msnu(1), q) +
     b0(p, tree.msnu(2), tree.msnu(2), q) +
     b0(p, tree.msnu(3), tree.msnu(3), q));

return sfermions;
}

template<class SoftPars>
double Softsusy<SoftPars>::pis1s1Fermions(double p, double q) const {
 /// fermions: 3rd family only for now
   double hb = displayDrBarPars().hb, htau = displayDrBarPars().htau;
   double    mtau    = displayDrBarPars().mtau;
   double    mb      = displayDrBarPars().mb;
   double fermions = 3.0 * sqr(hb) *
    ((sqr(p) - 4.0 * sqr(mb)) * 
     b0(p, mb, mb, q) - 2.0 * a0(mb, q));
  
  fermions = fermions + sqr(htau) *
    ((sqr(p) - 4.0 * sqr(mtau)) * b0(p, mtau, mtau, q) - 
     2.0 * a0(mtau, q));

  return fermions;
}
template<class SoftPars>
double Softsusy<SoftPars>::pis2s2Fermions(double p, double q) const {
 /// fermions: 3rd family only for now
   double mt =  displayDrBarPars().mt,  ht = displayDrBarPars().ht;
   double fermions = 3.0 * sqr(ht) *
    ((sqr(p) - 4.0 * sqr(mt)) * 
     b0(p, mt, mt, q) - 2.0 * a0(mt, q));
 return fermions;
}

template<class SoftPars>
double Softsusy<SoftPars>::pis1s1Higgs(double p, double q) const {
  
  double    thetaWDRbar = asin(calcSinthdrbar());
  double    costhDrbar  = cos(thetaWDRbar);
  double    cw2DRbar    = sqr(cos(thetaWDRbar));
  double    sw2DRbar    = 1.0 - cw2DRbar;
  double    alpha   = displayDrBarPars().thetaH;
  double    calpha2 = sqr(cos(alpha)), salpha2 = sqr(sin(alpha)), 
     s2alpha = sin(2.0 * alpha), c2alpha = cos(2.0 * alpha);
  double    mz      = displayMzRun();
  double    g       = displayGaugeCoupling(2);
  double    mHc     = displayDrBarPars().mHpm;
  double    mA      = displayDrBarPars().mA0(1);
  double beta = atan(displayTanb());
  double cosb = cos(beta), cosb2 = sqr(cosb), cos2b = cos(2.0 * beta), 
     sinb = sin(beta), sinb2 = sqr(sinb), sin2b = sin(2.0 * beta);

  double higgs = sqr(g) * 0.25 *
    (sinb2 * (2.0 * ffn(p, mHc, displayMwRun(), q) + 
	      ffn(p, mA, mz, q) / cw2DRbar) +
     cosb2 * (2.0 * ffn(p, displayMwRun(), displayMwRun(), q)  + 
	      ffn(p, mz, mz, q) / cw2DRbar)) +
    1.75 * sqr(g) * cosb2 * 
    (2.0 * sqr(displayMwRun()) * b0(p, displayMwRun(), displayMwRun(), q) + 
     sqr(mz) * b0(p, mz, mz, q) / cw2DRbar) +
    sqr(g) * (2.0 * a0(displayMwRun(), q) + a0(mz, q) / cw2DRbar);

  /// Trilinear Higgs couplings in basis H h G A: have assumed the couplings
  /// are symmetric (ie hHs1 = Hhs1)
  DoubleMatrix hhs1(4, 4);
  hhs1(1, 1) = cosb * (3.0 * calpha2 - salpha2) - sinb * s2alpha;
  hhs1(2, 2) = cosb * (3.0 * salpha2 - calpha2) + sinb * s2alpha;
  hhs1(1, 2) = -2.0 * cosb * s2alpha - sinb * c2alpha;
  hhs1(2, 1) = hhs1(1, 2);
  hhs1(3, 3) = cos2b * cosb;
  hhs1(4, 4) = -cos2b * cosb;
  hhs1(3, 4) = -sin2b * cosb; hhs1(4, 3) = hhs1(3, 4);
  hhs1 = hhs1 * (g * mz / (2.0 * costhDrbar));

  /// Quadrilinear Higgs couplings
  DoubleVector hhs1s1(4);
  hhs1s1(1) = 3.0 * calpha2 - salpha2;
  hhs1s1(2) = 3.0 * salpha2 - calpha2;
  hhs1s1(3) = cos2b; hhs1s1(4) = -cos2b;
  hhs1s1 = hhs1s1 * (sqr(g) * 0.25 / (sqr(costhDrbar)));

  /// define Higgs vector in 't-Hooft Feynman gauge, and couplings:
  DoubleVector higgsm(4), higgsc(2);
  DoubleVector dnu(4), dnd(4), cn(4);
  assignHiggs(higgsm, higgsc, dnu, dnd, cn, beta);

  for (int i=1; i<=4; i++) {
    for (int j=1; j<=4; j++) {
      higgs = higgs + 0.5 * sqr(hhs1(i, j)) * b0(p, higgsm(i), higgsm(j), q);
      ///      cout << "higgs(" << i << "," << j << ")=" << higgs;
    }
    higgs = higgs + 0.5 * hhs1s1(i) * a0(higgsm(i), q);
    ///      cout << "higgs(" << i << "," << j << ")=" << higgs;
  }

  ///  cout << hhs1 << hhs1s1 << higgsm;

  /// Basis (G+ H+, G- H-)
  DoubleMatrix hphps1(2, 2);
  hphps1(1, 1) = cos2b * cosb;
  hphps1(2, 2) = -cos2b * cosb + 2.0 * cw2DRbar * cosb;
  hphps1(1, 2) = -sin2b * cosb + cw2DRbar * sinb; 
  hphps1(2, 1) = hphps1(1, 2);
  hphps1 = hphps1 * (g * mz * 0.5 / costhDrbar);
 
  /// (G+ H+)
  DoubleVector hphps1s1(2);
  hphps1s1(1) = cw2DRbar + sw2DRbar * cos2b;
  hphps1s1(2) = cw2DRbar - sw2DRbar * cos2b;
  hphps1s1 = hphps1s1 * (sqr(g) * 0.25 / cw2DRbar);

  for (int i=1; i<=2; i++) {
    for (int j=1; j<=2; j++) 
      higgs = higgs + sqr(hphps1(i, j)) * b0(p, higgsc(i), higgsc(j), q);
    higgs = higgs + hphps1s1(i) * a0(higgsc(i), q);
  }

  return higgs;
}

template<class SoftPars>
double Softsusy<SoftPars>::pis1s2Higgs(double p, double q) const {

 double    thetaWDRbar = asin(calcSinthdrbar());
  double    costhDrbar  = cos(thetaWDRbar);
  double    cw2DRbar    = sqr(cos(thetaWDRbar));
  double    alpha   = displayDrBarPars().thetaH;
  double    calpha2 = sqr(cos(alpha)), salpha2 = sqr(sin(alpha)), 
     s2alpha = sin(2.0 * alpha), c2alpha = cos(2.0 * alpha);
  double    mz      = displayMzRun();
  double    g       = displayGaugeCoupling(2);
  double    mHc     = displayDrBarPars().mHpm;
  double    mA      = displayDrBarPars().mA0(1);
  double beta = atan(displayTanb());
  double cosb = cos(beta), cos2b = cos(2.0 * beta), 
     sinb = sin(beta), sin2b = sin(2.0 * beta);

 double higgs = sqr(g) * 0.25 * sinb * cosb *
    (2.0 * ffn(p, displayMwRun(), displayMwRun(), q) -2.0 * ffn(p, mHc, displayMwRun(), q) +
     (ffn(p, mz, mz, q) - ffn(p, mA, mz, q)) / cw2DRbar +
     7.0 * (2.0 * sqr(displayMwRun()) * b0(p, displayMwRun(), displayMwRun(), q) + 
	    sqr(mz) * b0(p, mz, mz, q) / cw2DRbar)); 
  
  /// Trilinear Higgs couplings in basis H h G A: have assumed the couplings
  /// are symmetric (ie hHs1 = Hhs1)
  DoubleMatrix hhs1(4, 4);
  hhs1(1, 1) = cosb * (3.0 * calpha2 - salpha2) - sinb * s2alpha;
  hhs1(2, 2) = cosb * (3.0 * salpha2 - calpha2) + sinb * s2alpha;
  hhs1(1, 2) = -2.0 * cosb * s2alpha - sinb * c2alpha;
  hhs1(2, 1) = hhs1(1, 2);
  hhs1(3, 3) = cos2b * cosb;
  hhs1(4, 4) = -cos2b * cosb;
  hhs1(3, 4) = -sin2b * cosb; hhs1(4, 3) = hhs1(3, 4);
  hhs1 = hhs1 * (g * mz / (2.0 * costhDrbar));
  /// Trilinear Higgs couplings in basis H h G A: have assumed the couplings
  /// are symmetric (ie hHs2 = Hhs2)
  DoubleMatrix hhs2(4, 4);
  hhs2(1, 1) = sinb * (3.0 * salpha2 - calpha2) - cosb * s2alpha;
  hhs2(2, 2) = sinb * (3.0 * calpha2 - salpha2) + cosb * s2alpha;
  hhs2(1, 2) = 2.0 * sinb * s2alpha - cosb * c2alpha;
  hhs2(2, 1) = hhs2(1, 2);
  hhs2(3, 3) = -cos2b * sinb;
  hhs2(4, 4) = cos2b * sinb;
  hhs2(3, 4) = sin2b * sinb; hhs2(4, 3) = hhs2(3, 4);
  hhs2 = hhs2 * (g * mz / (2.0 * costhDrbar));

  /// Quadrilinear Higgs couplings
  DoubleVector hhs1s2(4);
  hhs1s2(1) = -s2alpha;
  hhs1s2(2) = s2alpha;
  hhs1s2 = hhs1s2 * (sqr(g) * 0.25 / (sqr(costhDrbar)));

  ///  cout << "alpha=" << alpha << " g=" << g << " " << gp << " cw" << costhDrbar << " " << displayTanb();

  /// define Higgs vector in 't-Hooft Feynman gauge, and couplings:
  DoubleVector higgsm(4), higgsc(2);
  assignHiggs(higgsm, higgsc);

  double trilinear = 0.;
  double quartic = 0., quarticNeut = 0.;

  for (int i=1; i<=4; i++) {
    for (int j=1; j<=4; j++) {
      higgs = higgs + 0.5 * hhs1(i, j) * hhs2(i, j) * 
	b0(p, higgsm(i), higgsm(j), q);
        trilinear = trilinear + 0.5 * hhs1(i, j) * hhs2(i, j) * 
	  b0(p, higgsm(i), higgsm(j), q);
    }
    higgs = higgs + 0.5 * hhs1s2(i) * a0(higgsm(i), q);
    quarticNeut = quarticNeut + 0.5 * hhs1s2(i) * a0(higgsm(i), q);
  }

  DoubleMatrix hphps2(2, 2);
  hphps2(1, 1) = -cos2b * sinb;
  hphps2(2, 2) = cos2b * sinb + 2.0 * cw2DRbar * sinb;
  hphps2(1, 2) = sin2b * sinb - cw2DRbar * cosb; 
  hphps2(2, 1) = hphps2(1, 2);
  hphps2 = hphps2 * (g * mz * 0.5 / costhDrbar);
  DoubleMatrix hphps1(2, 2);
  hphps1(1, 1) = cos2b * cosb;
  hphps1(2, 2) = -cos2b * cosb + 2.0 * cw2DRbar * cosb;
  hphps1(1, 2) = -sin2b * cosb + cw2DRbar * sinb; 
  hphps1(2, 1) = hphps1(1, 2);
  hphps1 = hphps1 * (g * mz * 0.5 / costhDrbar);

  DoubleVector hphps1s2(2);
  hphps1s2(1) = - cw2DRbar * sin2b;
  hphps1s2(2) = cw2DRbar * sin2b;
  hphps1s2 = hphps1s2 * (sqr(g) * 0.25 / cw2DRbar);

  for (int i=1; i<=2; i++) {
    for (int j=1; j<=2; j++) {
      higgs = higgs + hphps1(i, j) * hphps2(i, j) 
	* b0(p, higgsc(i), higgsc(j), q);
        trilinear = trilinear + 
      hphps1(i, j) * hphps2(i, j) 
	  * b0(p, higgsc(i), higgsc(j), q);
    }
    higgs = higgs + hphps1s2(i) * a0(higgsc(i), q);
    quartic = quartic + hphps1s2(i) * a0(higgsc(i), q);
  }

  return higgs;
}

template<class SoftPars>
double Softsusy<SoftPars>::pis2s2Higgs(double p, double q) const {
double    thetaWDRbar = asin(calcSinthdrbar());
  double    costhDrbar  = cos(thetaWDRbar);
  double    cw2DRbar    = sqr(cos(thetaWDRbar));
  double    sw2DRbar    = 1.0 - cw2DRbar;
  double    alpha   = displayDrBarPars().thetaH;
  double    calpha2 = sqr(cos(alpha)), salpha2 = sqr(sin(alpha)), 
     s2alpha = sin(2.0 * alpha), c2alpha = cos(2.0 * alpha);
  double    mz      = displayMzRun();
  double    g       = displayGaugeCoupling(2);
  double    mHc     = displayDrBarPars().mHpm;
  double    mA      = displayDrBarPars().mA0(1);
  double beta = atan(displayTanb());
  double cosb = cos(beta), cosb2 = sqr(cosb), cos2b = cos(2.0 * beta), 
     sinb = sin(beta), sinb2 = sqr(sinb), sin2b = sin(2.0 * beta);
 double higgs = sqr(g) * 0.25 *
    (cosb2 * (2.0 * ffn(p, mHc, displayMwRun(), q) + 
	      ffn(p, mA, mz, q) / cw2DRbar) +
     sinb2 * (2.0 * ffn(p, displayMwRun(), displayMwRun(), q) + 
	      ffn(p, mz, mz, q) / cw2DRbar)) +
    1.75 * sqr(g) * sinb2 * 
    (2.0 * sqr(displayMwRun()) * b0(p, displayMwRun(), displayMwRun(), q) + 
     sqr(mz) * b0(p, mz, mz, q) / cw2DRbar) +
    sqr(g) * (2.0 * a0(displayMwRun(), q) + a0(mz, q) / cw2DRbar);

  double quartic = 0., trilinear = 0.;
 
  /// Trilinear Higgs couplings in basis H h G A: have assumed the couplings
  /// are symmetric (ie hHs2 = Hhs2)
  DoubleMatrix hhs2(4, 4);
  hhs2(1, 1) = sinb * (3.0 * salpha2 - calpha2) - cosb * s2alpha;
  hhs2(2, 2) = sinb * (3.0 * calpha2 - salpha2) + cosb * s2alpha;
  hhs2(1, 2) = 2.0 * sinb * s2alpha - cosb * c2alpha;
  hhs2(2, 1) = hhs2(1, 2);
  hhs2(3, 3) = -cos2b * sinb;
  hhs2(4, 4) = cos2b * sinb;
  hhs2(3, 4) = sin2b * sinb; hhs2(4, 3) = hhs2(3, 4);
  hhs2 = hhs2 * (g * mz / (2.0 * costhDrbar));

  /// Quadrilinear Higgs couplings
  DoubleVector hhs2s2(4);
  hhs2s2(1) = 3.0 * salpha2 - calpha2;
  hhs2s2(2) = 3.0 * calpha2 - salpha2;
  hhs2s2(3) = -cos2b; hhs2s2(4) = cos2b;
  hhs2s2 = hhs2s2 * (sqr(g) * 0.25 / (sqr(costhDrbar)));

  /// define Higgs vector in 't-Hooft Feynman gauge, and couplings:
  DoubleVector higgsm(4), higgsc(2);
  assignHiggs(higgsm, higgsc);

  for (int i=1; i<=4; i++) {
    for (int j=1; j<=4; j++) {
      higgs = higgs + 0.5 * sqr(hhs2(i, j)) * b0(p, higgsm(i), higgsm(j), q);
      trilinear = trilinear + 
	0.5 * sqr(hhs2(i, j)) * b0(p, higgsm(i), higgsm(j), q);
    }
    higgs = higgs + 0.5 * hhs2s2(i) * a0(higgsm(i), q);
    quartic = quartic + 0.5 * hhs2s2(i) * a0(higgsm(i), q);
  }

  DoubleMatrix hphps2(2, 2);
  hphps2(1, 1) = -cos2b * sinb;
  hphps2(2, 2) = cos2b * sinb + 2.0 * cw2DRbar * sinb;
  hphps2(1, 2) = sin2b * sinb - cw2DRbar * cosb; 
  hphps2(2, 1) = hphps2(1, 2);
  hphps2 = hphps2 * (g * mz * 0.5 / costhDrbar);

  DoubleVector hphps2s2(2);
  hphps2s2(1) = cw2DRbar - sw2DRbar * cos2b;
  hphps2s2(2) = cw2DRbar + sw2DRbar * cos2b;
  hphps2s2 = hphps2s2 * (sqr(g) * 0.25 / cw2DRbar);

  for (int i=1; i<=2; i++) {
    for (int j=1; j<=2; j++) {
      higgs = higgs + sqr(hphps2(i, j)) * b0(p, higgsc(i), higgsc(j), q);
      trilinear = trilinear + sqr(hphps2(i, j)) * 
	b0(p, higgsc(i), higgsc(j), q);
    }
    higgs = higgs + hphps2s2(i) * a0(higgsc(i), q);
    quartic = quartic + hphps2s2(i) * a0(higgsc(i), q);
  }
return higgs;
}

template<class SoftPars>
double Softsusy<SoftPars>::pis1s1Neutralinos(double p, double q) const {
  double    g       = displayGaugeCoupling(2);
  double    gp      = displayGaugeCoupling(1) * sqrt(0.6);
  double neutralinos = 0.0;

  DoubleMatrix aPsi(4, 4);
  ComplexMatrix aChi(4, 4), bChi(4, 4);
  ComplexMatrix n(displayDrBarPars().nBpmz);
  DoubleVector mneut(displayDrBarPars().mnBpmz);
 

  aPsi(1, 3) = -gp * 0.5; 
  aPsi(2, 3) = g * 0.5; 
  aPsi.symmetrise();
  aChi = n.complexConjugate() * aPsi * n.hermitianConjugate();
  bChi = n * aPsi * n.transpose();

 
DoubleMatrix fChiChis1s1(4, 4), gChiChis1s1(4, 4);
  for(int i=1; i<=4; i++)
    for (int j=1; j<=4; j++) {
      fChiChis1s1(i, j) = sqr(aChi(i, j).mod()) + sqr(bChi(i, j).mod());
      gChiChis1s1(i, j) = (bChi(i, j).conj() * aChi(i, j) + 
	aChi(i, j).conj() * bChi(i, j)).real();
      neutralinos = neutralinos + 0.5 * 
	(fChiChis1s1(i, j) * gfn(p, mneut(i), mneut(j), q) - 2.0 *
	 gChiChis1s1(i, j) * mneut(i) * mneut(j) * 
	 b0(p, mneut(i), mneut(j), q));
    }

  return neutralinos;
}

template<class SoftPars>
double Softsusy<SoftPars>::pis1s2Neutralinos(double p, double q) const {
  double    g       = displayGaugeCoupling(2);
  double    gp      = displayGaugeCoupling(1) * sqrt(0.6);
  double neutralinos = 0.0;

  DoubleMatrix aPsi1(4, 4);
  ComplexMatrix aChi1(4, 4), bChi1(4, 4);
  ComplexMatrix n(displayDrBarPars().nBpmz);
  DoubleVector mneut(displayDrBarPars().mnBpmz);
 
  aPsi1(1, 3) = -gp * 0.5; 
  aPsi1(2, 3) = g * 0.5; 
  aPsi1.symmetrise();
  aChi1 = n.complexConjugate() * aPsi1 * n.hermitianConjugate();
  bChi1 = n * aPsi1 * n.transpose();
  DoubleMatrix aPsi2(4, 4);
  ComplexMatrix aChi2(4, 4), bChi2(4, 4);
  aPsi2(1, 4) = gp * 0.5; 
  aPsi2(2, 4) = -g * 0.5; 
  aPsi2.symmetrise();
  aChi2 = n.complexConjugate() * aPsi2 * n.hermitianConjugate();
  bChi2 = n * aPsi2 * n.transpose();

  DoubleMatrix fChiChis1s2(4, 4), gChiChis1s2(4, 4);
  for(int i=1; i<=4; i++)
    for (int j=1; j<=4; j++) {
      fChiChis1s2(i, j) = (aChi1(i, j).conj() * aChi2(i, j) + 
	bChi1(i, j).conj() * bChi2(i, j)).real();
      gChiChis1s2(i, j) = (bChi1(i, j).conj() * aChi2(i, j) + 
	aChi1(i, j).conj() * bChi2(i, j)).real();
      neutralinos = neutralinos + 0.5 * 
	(fChiChis1s2(i, j) * gfn(p, mneut(i), mneut(j), q) - 2.0 *
	 gChiChis1s2(i, j) * mneut(i) * mneut(j) * 
	 b0(p, mneut(i), mneut(j), q));
    }

 return neutralinos;
}

template<class SoftPars>
double Softsusy<SoftPars>::pis2s2Neutralinos(double p, double q) const {
  double    g       = displayGaugeCoupling(2);
  double    gp      = displayGaugeCoupling(1) * sqrt(0.6);
/// Neutralino contribution
  double neutralinos = 0.0;

  DoubleMatrix aPsi(4, 4);
  ComplexMatrix aChi(4, 4), bChi(4, 4);
  ComplexMatrix n(displayDrBarPars().nBpmz);
  DoubleVector mneut(displayDrBarPars().mnBpmz);
  
  aPsi(1, 4) = gp * 0.5; 
  aPsi(2, 4) = -g * 0.5; 
  aPsi.symmetrise();
  aChi = n.complexConjugate() * aPsi * n.hermitianConjugate();
  bChi = n * aPsi * n.transpose();

  DoubleMatrix fChiChis2s2(4, 4), gChiChis2s2(4, 4);
  for(int i=1; i<=4; i++)
    for (int j=1; j<=4; j++) {
      fChiChis2s2(i, j) = sqr(aChi(i, j).mod()) + sqr(bChi(i, j).mod());
      gChiChis2s2(i, j) = (bChi(i, j).conj() * aChi(i, j) + 
	aChi(i, j).conj() * bChi(i, j)).real();
      neutralinos = neutralinos + 0.5 * 
	(fChiChis2s2(i, j) * gfn(p, mneut(i), mneut(j), q) - 2.0 *
	 gChiChis2s2(i, j) * mneut(i) * mneut(j) * 
	 b0(p, mneut(i), mneut(j), q));
    }

return neutralinos;
}

template<class SoftPars>
double Softsusy<SoftPars>::pis1s1Charginos(double p, double q) const {
  double chargino = 0.0;
  double g = displayGaugeCoupling(2);
  ComplexMatrix u(displayDrBarPars().uBpmz), v(displayDrBarPars().vBpmz); 
  DoubleVector mch(displayDrBarPars().mchBpmz); 
  DoubleMatrix aPsic(2, 2);
  aPsic(1, 2) = g / root2; 
  ComplexMatrix aChic(2, 2), bChic(2, 2);
  DoubleMatrix fChiChis1s1(2, 2), gChiChis1s1(2, 2);
  aChic = v.complexConjugate() * aPsic * u.hermitianConjugate();
  bChic = u * aPsic.transpose() * v.transpose();
  for(int i=1; i<=2; i++)
    for (int j=1; j<=2; j++) {
      fChiChis1s1(i, j) = sqr(aChic(i, j).mod()) + sqr(bChic(i, j).mod());
      gChiChis1s1(i, j) = (bChic(i, j).conj() * aChic(i, j) + 
	aChic(i, j).conj() * bChic(i, j)).real();
      chargino = chargino + 
	(fChiChis1s1(i, j) * gfn(p, mch(i), mch(j), q) - 2.0 *
	 gChiChis1s1(i, j) * mch(i) * mch(j) * 
	 b0(p, mch(i), mch(j), q));
    }
  
  return chargino;
}

template<class SoftPars>
double Softsusy<SoftPars>::pis1s2Charginos(double p, double q) const {
  double g = displayGaugeCoupling(2);
 double chargino = 0.0;
  ComplexMatrix u(displayDrBarPars().uBpmz), v(displayDrBarPars().vBpmz); 
  DoubleVector mch(displayDrBarPars().mchBpmz); 

  DoubleMatrix aPsic1(2, 2), aPsic2(2, 2);
  aPsic1(1, 2) = g / root2; 
  ComplexMatrix aChic1(2, 2), bChic1(2, 2);
  ComplexMatrix aChic2(2, 2), bChic2(2, 2);
  aChic1 = v.complexConjugate() * aPsic1 * u.hermitianConjugate();
  bChic1 = u * aPsic1.transpose() * v.transpose();
  aPsic2(2, 1) = g / root2;
  aChic2 = v.complexConjugate() * aPsic2 * u.hermitianConjugate();
  bChic2 = u * aPsic2.transpose() * v.transpose();
  DoubleMatrix fChiChis1s2(2, 2), gChiChis1s2(2, 2);
  for(int i=1; i<=2; i++)
    for (int j=1; j<=2; j++) {
      fChiChis1s2(i, j) = (aChic1(i, j).conj() * aChic2(i, j) + 
	bChic1(i, j).conj() * bChic2(i, j)).real();
      gChiChis1s2(i, j) = (bChic1(i, j).conj() * aChic2(i ,j) + 
	aChic1(i, j).conj() * bChic2(i, j)).real();
      chargino = chargino + 
	(fChiChis1s2(i, j) * gfn(p, mch(i), mch(j), q) - 2.0 *
	 gChiChis1s2(i, j) * mch(i) * mch(j) * 
	 b0(p, mch(i), mch(j), q));
    }

 return chargino;
}
template<class SoftPars>
double Softsusy<SoftPars>::pis2s2Charginos(double p, double q) const {
  double chargino = 0.0;
  double g = displayGaugeCoupling(2);
  ComplexMatrix u(displayDrBarPars().uBpmz), v(displayDrBarPars().vBpmz); 
  DoubleVector mch(displayDrBarPars().mchBpmz); 
  DoubleMatrix aPsic(2, 2);
  aPsic(2, 1) = g / root2;
  ComplexMatrix aChic(2, 2), bChic(2, 2);
  aChic = v.complexConjugate() * aPsic * u.hermitianConjugate();
  bChic = u * aPsic.transpose() * v.transpose();
  DoubleMatrix fChiChis2s2(2, 2), gChiChis2s2(2, 2);
  for(int i=1; i<=2; i++)
    for (int j=1; j<=2; j++) {
      fChiChis2s2(i, j) = sqr(aChic(i, j).mod()) + sqr(bChic(i, j).mod());
      gChiChis2s2(i, j) = (bChic(i, j).conj() * aChic(i, j) + 
	aChic(i, j).conj() * bChic(i, j)).real();
      chargino = chargino + 
	(fChiChis2s2(i, j) * gfn(p, mch(i), mch(j), q) - 2.0 *
	 gChiChis2s2(i, j) * mch(i) * mch(j) * 
	 b0(p, mch(i), mch(j), q));
    }

 return chargino;
}


template<class SoftPars>
double Softsusy<SoftPars>::pis1s1(double p, double q) const {
  drBarPars tree(displayDrBarPars());
  double    beta    = atan(displayTanb());
  double    mb      = tree.mb;
  double    hb      = tree.hb;
  double    thetaWDRbar = asin(calcSinthdrbar());
  double    costhDrbar  = cos(thetaWDRbar);
  double    smu     = -displaySusyMu(); /// minus sign taken into acct here!
  double    mz      = displayMzRun();
  double    g       = displayGaugeCoupling(2);
  double cosb = cos(beta);
  double gmzOcthW = g * mz / costhDrbar;
  //PA: get fermion contribution
  double fermions = pis1s1Fermions(p, q);
  // sfermion couplings to s1 Higgs state
  DoubleMatrix ls1tt(2, 2), ls1bb(2, 2), ls1tautau(2, 2);
  H1SfSfCouplings(ls1tt, ls1bb, ls1tautau, gmzOcthW, smu, cosb, root2*mb/hb);
  //PA: get sfermion contribution
  double sfermions = pis1s1Sfermions(p, q, ls1tt, ls1bb, ls1tautau);
  //PA: get Higgs contribution
  double higgs = pis1s1Higgs(p, q);
  /// Neutralino contribution
  double neutralinos = pis1s1Neutralinos(p, q);
  /// Chargino contribution
  double chargino = pis1s1Charginos(p, q);  

  return 
    (sfermions + 
     fermions + higgs + neutralinos + chargino) / (16.0 * sqr(PI));
}

template<class SoftPars>
double Softsusy<SoftPars>::pis1s2(double p, double q) const {
  drBarPars tree(displayDrBarPars());

  double    beta    = atan(displayTanb());
  double    mb   =  tree.mb, hb   =  tree.hb; 
  double    thetaWDRbar = asin(calcSinthdrbar());
  double    costhDrbar  = cos(thetaWDRbar);
  double    smu     = -displaySusyMu(); /// minus sign taken into acct here!
  double    g       = displayGaugeCoupling(2);
  double cosb = cos(beta), sinb = sin(beta);
  double  mz = displayMzRun();
  // sfermion couplings to s1 Higgs state
  DoubleMatrix ls1tt(2, 2), ls1bb(2, 2), ls1tautau(2, 2);
  double gmzOcthW = g * mz / costhDrbar;
  H1SfSfCouplings(ls1tt, ls1bb, ls1tautau, gmzOcthW, smu, cosb, root2*mb/hb);
  /// sfermion couplings to s2 Higgs state
  DoubleMatrix ls2tt(2, 2), ls2bb(2, 2), ls2tautau(2, 2);
  H2SfSfCouplings(ls2tt, ls2bb, ls2tautau, gmzOcthW, smu, sinb);
  //PA: get sfermion contribution
  double sfermions = pis1s2Sfermions(p, q, ls1tt, ls1bb, ls1tautau, ls2tt, ls2bb, ls2tautau);
  //PA: get Higgs contribution
  double higgs = pis1s2Higgs(p, q);
  /// Neutralino contribution
  double neutralinos = pis1s2Neutralinos(p, q); 
  /// Chargino contribution
  double chargino = pis1s2Charginos(p, q);  

  return (sfermions + higgs + neutralinos + chargino) 
    / (16.0 * sqr(PI));
}

/// checked 28.10.02
template<class SoftPars>
double Softsusy<SoftPars>::pis2s2(double p, double q) const {
  drBarPars tree(displayDrBarPars());
  double beta = atan(displayTanb());
  double thetaWDRbar = asin(calcSinthdrbar());
  double costhDrbar = cos(thetaWDRbar);
  double smu = -displaySusyMu(); /// minus sign taken into acct here!
  double g = displayGaugeCoupling(2);
  double sinb = sin(beta);
  double mz = displayMzRun();
  double gmzOcthW = g * mz / costhDrbar;
  
  double fermions = pis2s2Fermions(p, q);
  /// sfermion couplings to s2 Higgs state
  DoubleMatrix ls2tt(2, 2), ls2bb(2, 2), ls2tautau(2, 2);
  H2SfSfCouplings(ls2tt, ls2bb, ls2tautau, gmzOcthW, smu, sinb);
  double sfermions = pis2s2Sfermions(p, q, ls2tt, ls2bb, ls2tautau);
  double higgs = pis2s2Higgs(p, q);
  double neutralinos = pis2s2Neutralinos(p, q); 
  double chargino = pis2s2Charginos(p, q);   

  return (fermions + sfermions + higgs + neutralinos + chargino) 
    / (16.0 * sqr(PI));
}

/// LCT: Returns trilinear neutralino-chargino-hpm coupling in unrotated basis
template<class SoftPars>
void Softsusy<SoftPars>::getNeutralinoCharginoHpmCoup(ComplexMatrix & apph1, ComplexMatrix & /* apph2 */, ComplexMatrix & /* bpph1 */, ComplexMatrix & bpph2) const {
  double gp = sqrt(0.6) * displayGaugeCoupling(1);
  double g = displayGaugeCoupling(2);

  apph1(1, 2) = gp / root2;
  bpph2(1, 2) = gp / root2;
  apph1(2, 2) = g / root2;
  bpph2(2, 2) = g / root2;
  apph1(3, 1) = -g;
  bpph2(4, 1) = g;
}

/// LCT: Returns fermion contributions to charged Higgs self-energy
template<class SoftPars>
double Softsusy<SoftPars>::piHpHmFermions(double p, double q) const {
  /// fermions: 3rd family only for now
  drBarPars tree(displayDrBarPars());
  double    ht      = tree.ht;
  double    hb      = tree.hb;
  double    htau    = tree.htau;
  double    mt      = tree.mt;
  double    mb      = tree.mb;
  double    mtau    = tree.mtau;
  double    beta    = atan(displayTanb());
  double    cosb    = cos(beta), cosb2 = sqr(cosb);
  double    sinb    = sin(beta), sinb2 = sqr(sinb), sin2b = sin(2.0 * beta);

  double fermions = 3.0 * (sqr(ht) * cosb2 + sqr(hb) * sinb2) * 
    gfn(p, mt, mb, q) + sinb2 * sqr(htau) * gfn(p, 0.0, mtau, q) -
    6.0 * hb * ht * mt * mb * sin2b * b0(p, mt, mb, q);

  return fermions;
}

/// LCT: Returns sfermion contributions to charged Higgs self-energy
template<class SoftPars>
double Softsusy<SoftPars>::piHpHmSfermions(double p, double q, double mu) const {
  drBarPars tree(displayDrBarPars());
  double    ht          = tree.ht;
  double    hb          = tree.hb;
  double    htau        = tree.htau;
  double    beta        = atan(displayTanb());
  double    mt          = tree.mt;
  double    mb          = tree.mb;
  double    mtau        = tree.mtau;
  double    thetaWDRbar = asin(calcSinthdrbar());
  double    cwDRbar     = cos(thetaWDRbar);
  double    cw2DRbar    = sqr(cwDRbar);
  double    thetat      = tree.thetat;
  double    thetab      = tree.thetab;
  double    thetatau    = tree.thetatau;
  double    st          = sin(thetat);
  double    sb          = sin(thetab);
  double    ct          = cos(thetat);
  double    cb          = cos(thetab);
  double    stau        = sin(thetatau);
  double    ctau        = cos(thetatau);
  double    g           = displayGaugeCoupling(2);
  double    mwRun       = displayMwRun();
  double    cosb        = cos(beta), cosb2 = sqr(cosb), cos2b = cos(2.0 * beta);
  double    sinb        = sin(beta), sinb2 = sqr(sinb), sin2b = sin(2.0 * beta);

  /// first two generations: forget lighter Yukawas
  DoubleMatrix lHpud(2, 2), lHpen(2, 2);
  lHpud(1, 1) = g * mwRun * sin2b / root2;

  double sfermions = 3.0 * sqr(lHpud(1, 1)) * 
    (b0(p, tree.mu(1, 1), tree.md(1, 1), q) + 
     b0(p, tree.mu(1, 2), tree.md(1, 2), q)) + sqr(lHpud(1, 1)) *
    (b0(p, tree.msnu(1), tree.me(1, 1), q) + 
     b0(p, tree.msnu(2), tree.me(1, 2), q));

  /// 3rd family
  lHpud(1, 1) = g * mwRun * sin2b / root2 
    - ht * mt * cosb - hb * mb * sinb;
  lHpud(1, 2) = hb * mu * cosb - tree.ub * sinb;
  lHpud(2, 1) = ht * mu * sinb - tree.ut * cosb;
  lHpud(2, 2) = -ht * mb * cosb - hb * mt * sinb;
  lHpud = rot2d(tree.thetat) * lHpud * rot2d(-tree.thetab);

  lHpen(1, 1) = g * mwRun * sin2b / root2 
     - htau * mtau * sinb;
  lHpen(1, 2) = htau * mu * cosb - tree.utau * sinb;
  lHpen(2, 1) = 0.;
  lHpen(2, 2) = 0.;
  lHpen = lHpen * rot2d(-tree.thetatau);

  double thirdFamSfermions = 0.;
  int i, j; for(i=1; i<=2; i++) 
    for (j=1; j<=2; j++) 
      thirdFamSfermions = thirdFamSfermions +
	3.0 * sqr(lHpud(i, j)) * b0(p, tree.mu(i, 3), tree.md(j, 3), q);
    
  for (j=1; j<=2; j++) 
    thirdFamSfermions = thirdFamSfermions +
      sqr(lHpen(1, j)) * b0(p, tree.msnu(3), tree.me(j, 3), q);

  /// ups
  sfermions = sfermions + 
    3.0 * sqr(g) * cos2b * (-0.5 * guL / cw2DRbar + 0.5) * 
    (a0(tree.mu(1, 1), q) + a0(tree.mu(1, 2), q)) +
    3.0 * sqr(g) * cos2b * (-0.5 * guR / cw2DRbar) * 
    (a0(tree.mu(2, 1), q) + a0(tree.mu(2, 2), q));
  /// downs
  sfermions = sfermions + 
    3.0 * sqr(g) * cos2b * (-0.5 * gdL / cw2DRbar - 0.5) * 
    (a0(tree.md(1, 1), q) + a0(tree.md(1, 2), q)) +
    3.0 * sqr(g) * cos2b * (-0.5 * gdR / cw2DRbar) * 
    (a0(tree.md(2, 1), q) + a0(tree.md(2, 2), q));
  /// sneutrinos/selectrons
  sfermions = sfermions + 
    sqr(g) * cos2b * (-0.5 * gnuL / cw2DRbar + 0.5) * 
    (a0(tree.msnu(1), q) + a0(tree.msnu(2), q));
  /// downs
  sfermions = sfermions + 
    sqr(g) * cos2b * (-0.5 * geL / cw2DRbar - 0.5) * 
    (a0(tree.me(1, 1), q) + a0(tree.me(1, 2), q)) +
    sqr(g) * cos2b * (-0.5 * geR / cw2DRbar) * 
    (a0(tree.me(2, 1), q) + a0(tree.me(2, 2), q));  
  /// 3rd fam ups
  thirdFamSfermions = thirdFamSfermions +
    3.0 * (sqr(hb) * sinb2 - sqr(g) * cos2b * 0.5 * guL / cw2DRbar + 
	   0.5 * cos2b * sqr(g)) * 
    (sqr(ct) * a0(tree.mu(1, 3), q) + sqr(st) * a0(tree.mu(2, 3), q)) +
    3.0 * (sqr(ht) * cosb2 - sqr(g) * cos2b * 0.5 * guR / cw2DRbar) *
    (sqr(st) * a0(tree.mu(1, 3), q) + sqr(ct) * a0(tree.mu(2, 3), q));
  /// 3rd fam downs
  thirdFamSfermions = thirdFamSfermions +
    3.0 * (sqr(ht) * cosb2 - sqr(g) * cos2b * 0.5 * gdL / cw2DRbar -
	   0.5 * cos2b * sqr(g)) * 
    (sqr(cb) * a0(tree.md(1, 3), q) + sqr(sb) * a0(tree.md(2, 3), q)) +
    3.0 * (sqr(hb) * sinb2 - sqr(g) * cos2b * 0.5 * gdR / cw2DRbar) *
    (sqr(sb) * a0(tree.md(1, 3), q) + sqr(cb) * a0(tree.md(2, 3), q));
  /// 3rd fam snus
  thirdFamSfermions = thirdFamSfermions +
    (sqr(htau) * sinb2 - sqr(g) * cos2b * 0.5 * gnuL / cw2DRbar + 
     0.5 * cos2b * sqr(g)) * a0(tree.msnu(3), q);
  /// 3rd fam es
  thirdFamSfermions = thirdFamSfermions +
    (-sqr(g) * cos2b * 0.5 * geL / cw2DRbar - 0.5 * cos2b * sqr(g)) * 
    (sqr(ctau) * a0(tree.me(1, 3), q) + sqr(stau) * a0(tree.me(2, 3), q)) +
    (sqr(htau) * sinb2 - sqr(g) * cos2b * 0.5 * geR / cw2DRbar) *
    (sqr(stau) * a0(tree.me(1, 3), q) + sqr(ctau) * a0(tree.me(2, 3), q));

  sfermions = sfermions + thirdFamSfermions;

  return sfermions;
}

/// LCT: Returns pure gauge contributions to charged Higgs self-energy
template<class SoftPars>
double Softsusy<SoftPars>::piHpHmGauge(double /* p */, double q) const {
  drBarPars tree(displayDrBarPars());
  double    mz          = displayMzRun();
  double    thetaWDRbar = asin(calcSinthdrbar());
  double    cwDRbar     = cos(thetaWDRbar);
  double    cw2DRbar    = sqr(cwDRbar);
  double    g           = displayGaugeCoupling(2);
  double    mwRun       = displayMwRun();

  double gauge = 2.0 * sqr(g) * a0(mwRun, q) +
  sqr(g) * sqr(cos(2.0 * thetaWDRbar)) / cw2DRbar * a0(mz, q);

 return gauge;
}

/// LCT: Returns Higgs contributions to charged Higgs self-energy
template<class SoftPars>
double Softsusy<SoftPars>::piHpHmHiggs(double p, double q) const {
  drBarPars tree(displayDrBarPars());
  double    alpha       = tree.thetaH;
  double    mz          = displayMzRun();
  double    thetaWDRbar = asin(calcSinthdrbar());
  double    cwDRbar     = cos(thetaWDRbar);
  double    cw2DRbar    = sqr(cwDRbar);
  double    sw2DRbar    = 1.0 - cw2DRbar;
  double    g           = displayGaugeCoupling(2);
  double    mh0         = tree.mh0(1);
  double    mHc         = tree.mHpm;
  double    mH          = tree.mh0(2);
  double    mA          = tree.mA0(1);
  double    mwRun       = displayMwRun();
  double    beta        = atan(displayTanb());
  double    cosb        = cos(beta), cos2b = cos(2.0 * beta);
  double    sinb        = sin(beta), sin2b = sin(2.0 * beta);

  double higgs = sqr(g) * 0.25 * 
    (sqr(sin(alpha - beta)) * ffn(p, mH, mwRun, q) +
     sqr(cos(alpha - beta)) * ffn(p, mh0, mwRun, q) +
     ffn(p, mA, mwRun, q) + sqr(cos(2.0 * thetaWDRbar)) / cw2DRbar * 
     ffn(p, mHc, mz, q)) +
    sqr(g) * sw2 * ffn(p, mHc, 0., q)  +
    sqr(g) * sqr(mwRun) * 0.25 * b0(p, mwRun, mA, q);

  DoubleMatrix lHpH0Hm(2, 2);

  /// BPMZ not so clear: assuming (D.67) that s_i G+ H- coupling is same as 
  /// s_i G- H+ (by C conservation): basis is (s1 s2, G- H-)
  lHpH0Hm(1, 1) = -sin2b * cosb + cw2DRbar * sinb;
  lHpH0Hm(2, 1) =  sin2b * sinb - cw2DRbar * cosb;
  lHpH0Hm(1, 2) = -cos2b * cosb + 2.0 * cw2DRbar * cosb;
  lHpH0Hm(2, 2) =  cos2b * sinb + 2.0 * cw2DRbar * sinb;

  /// in (H h) basis
  lHpH0Hm = rot2d(alpha) * lHpH0Hm * (g * mz * 0.5 / cwDRbar); 
  DoubleVector spech(2), specHm(2); 
  spech(1) = mH; spech(2) = mh0;
  specHm(1) = mwRun; specHm(2) = mHc;
  int i,j;
  for(i=1; i<=2; i++) 
    for(j=1; j<=2; j++)
      higgs = higgs + sqr(lHpH0Hm(i, j)) * b0(p, spech(i), specHm(j), q);

  double l1 = (2.0 * sqr(sin2b) - 1.0) * sqr(g) * 0.25 / cw2DRbar;
  double l2 = 2.0 * sqr(cos2b) * sqr(g) * 0.25 / cw2DRbar;
  higgs = higgs + l1 * a0(specHm(1), q) + l2 * a0(specHm(2), q);
  DoubleMatrix lHpHmH0H0(2, 2);
  lHpHmH0H0(1, 1) = cw2DRbar - sw2DRbar * cos2b;
  lHpHmH0H0(1, 2) = cw2DRbar * sin2b;
  lHpHmH0H0(2, 1) = lHpHmH0H0(1, 2); 
  lHpHmH0H0(2, 2) = cw2DRbar + sw2DRbar * cos2b;
  lHpHmH0H0 = rot2d(alpha) * lHpHmH0H0 * rot2d(-alpha) *
    (0.25 * sqr(g) / cw2DRbar);
  double lHpHmG0G0 = (cw2DRbar * (1.0 + sqr(sin2b)) - sw2DRbar * sqr(cos2b)) *
    0.25 * sqr(g) / cw2DRbar;
  double lHpHmA0A0 = sqr(cos2b) * 0.25 * sqr(g) / cw2DRbar;
  higgs = higgs + 0.5 * (lHpHmH0H0(1, 1) * a0(spech(1), q) +
    lHpHmH0H0(2, 2) * a0(spech(2), q) + lHpHmG0G0 * a0(mz, q) + 
    lHpHmA0A0 * a0(mA, q)); 

  return higgs;
}

/// LCT: Returns neutralino-chargino contributions to charged Higgs 
/// self-energy
template<class SoftPars>
double Softsusy<SoftPars>::piHpHmGauginos(double p, double q) const {
  drBarPars tree(displayDrBarPars());
  double    beta    = atan(displayTanb());
  double    cosb    = cos(beta);
  double    sinb    = sin(beta);

  double gauginos = 0.;

  ComplexMatrix apph1(4, 2), apph2(4, 2), bpph1(4, 2), bpph2(4, 2);

  getNeutralinoCharginoHpmCoup(apph1, apph2, bpph1, bpph2);

  /// Get to physical gaugino estates, and G+H+
  apph1 = -sinb * tree.nBpmz.complexConjugate() * 
    apph1 * tree.uBpmz.hermitianConjugate();
  bpph2 = cosb * tree.nBpmz * bpph2 * tree.vBpmz.transpose();
  
  int i,j;
   for (i=1; i<=4; i++)
    for(j=1; j<=2; j++) {
      double f = sqr(apph1(i, j).mod()) + sqr(bpph2(i, j).mod());
      double ga = 2.0 * (bpph2(i, j).conj() * apph1(i, j)).real();
      gauginos = gauginos + f * gfn(p, tree.mchBpmz(j), tree.mnBpmz(i), q) -
	2.0 * ga * tree.mchBpmz(j) * tree.mnBpmz(i) * 
	b0(p, tree.mchBpmz(j), tree.mnBpmz(i), q);
    }

   return gauginos;
}

/// LCT: Charged Higgs self-energy 
template<class SoftPars>
double Softsusy<SoftPars>::piHpHm(double p, double q) const {
  double mu = -displaySusyMu(); ///<< LCT: Note minus sign. Consistent with 
  /// SOFTSUSY

  /// LCT: fermion contribution
  double fermions = piHpHmFermions(p,q);
  /// LCT: sfermion contribution
  double sfermions = piHpHmSfermions(p, q, mu);
  /// LCT: pure gauge contribution
  double gauge = piHpHmGauge(p, q);
  /// LCT: Higgs contribution
  double higgs = piHpHmHiggs(p, q);
  /// LCT: neutralino-chargino contribution
  double gauginos = piHpHmGauginos(p, q);

  const static double loopFac = 1.0 / (16 * sqr(PI));
  double pihh = fermions + sfermions + gauge + higgs + gauginos;

  return loopFac * pihh;
}


template<class SoftPars>
double Softsusy<SoftPars>::piAA(double p, double q) const {/// checked 30.07.03
  drBarPars tree(displayDrBarPars());

 if (tree.mu(1, 3) == 0.0 || tree.mu(2, 3) == 0.0) {
   if (PRINTOUT > 1)
    cout << "Trying to calculate piAA without having first calculated"
	 << " the DRbar masses.\n";
   return 0.0; 
  }

  double    beta    = atan(displayTanb());
  double    mt     =  tree.mt;
  double    mb    =  tree.mb;
  double    mtau =    tree.mtau;
  double    mmu =  displayDataSet().displayMass(legacy::mMuon);
  double    mE =  displayDataSet().displayMass(legacy::mElectron);
  double    mz = displayMzRun();
  double    thetaWDRbar = asin(calcSinthdrbar());
  double    cw2DRbar    = sqr(cos(thetaWDRbar));
  double    sw2DRbar    = 1.0 - cw2DRbar;
  double    thetat  = tree.thetat ;
  double    thetab  = tree.thetab;
  double    thetatau= tree.thetatau;
  double    st      = sin(thetat) ;
  double    sb      = sin(thetab) ;
  double    ct      = cos(thetat) ;
  double    cb      = cos(thetab) ;
  double    stau    = sin(thetatau);
  double    ctau    = cos(thetatau);
  double    g       = displayGaugeCoupling(2);
  double    mh0     = maximum(tree.mh0(1), EPSTOL); ///< protects vs zeros
  double    mHc     = maximum(tree.mHpm, EPSTOL);
  double    mH     = maximum(tree.mh0(2), EPSTOL);
  double    mA      = maximum(tree.mA0(1), EPSTOL);
  double    gp      = displayGaugeCoupling(1) * sqrt(0.6);
  double cosb = cos(beta), cosb2 = sqr(cosb), cos2b = cos(2.0 * beta),
    sinb = sin(beta), sinb2 = sqr(sinb), sin2b = sin(2.0 * beta);

  /// Up fermions/sfermions
  double fermions = 0., upsfermions = 0.0;
  fermions = fermions + cosb2 * 
    (3.0 * sqr(tree.ht) * 
     (sqr(p) * b0(p, mt, mt, q) - 2.0 * a0(mt, q))); /// ignore 1st 2 families
 
  /// LH gens 1-2 
  upsfermions = upsfermions +
    3.0 * (-sqr(g) / (cw2DRbar) * cos2b * 0.5 * guL) *
    (a0(tree.mu(1, 1), q) + a0(tree.mu(1, 2), q));
 
  /// RH gen 1-2
  upsfermions = upsfermions +
    3.0 * (-sqr(g) / (cw2DRbar) * cos2b * 0.5 * guR) *
    (a0(tree.mu(2, 1), q) + a0(tree.mu(2, 2), q));
  
  upsfermions = upsfermions +
    3.0 * (sqr(displayYukawaElement(YU, 1, 1) * 
	       (-displaySusyMu() * sinb - 
		displaySoftA(UA, 1, 1) * cosb))) * 0.5 * 
    (b0(p, tree.mu(1, 1), tree.mu(2, 1), q) +
     b0(p, tree.mu(2, 1), tree.mu(1, 1), q)) + 
    3.0 * (sqr(displayYukawaElement(YU, 2, 2) * 
	       (-displaySusyMu() * sinb - 
		displaySoftA(UA, 2, 2) * cosb))) * 0.5 * 
    (b0(p, tree.mu(1, 2), tree.mu(2, 2), q) +
     b0(p, tree.mu(2, 2), tree.mu(1, 2), q));

  double stops = 0.0;
  stops = stops + 
    3.0 * (sqr(tree.ht) * cosb2 - 
	   sqr(g) / (cw2DRbar) * cos2b * 0.5 * guL) *
    (sqr(ct) * a0(tree.mu(1, 3), q) + 
     sqr(st) * a0(tree.mu(2, 3), q));

  stops = stops + 
    3.0 * (sqr(tree.ht) * cosb2 - 
	   sqr(g) / (cw2DRbar) * cos2b * 0.5 * guR) *
    (sqr(st) * a0(tree.mu(1, 3), q) + 
     sqr(ct) * a0(tree.mu(2, 3), q));

  double auu12 =  
    (tree.ht *
     displaySusyMu() * sinb + tree.ut * cosb) / root2;

  int i, j;
  stops = stops + 
    3.0 * sqr(auu12) * 
    (b0(p, tree.mu(1, 3), tree.mu(2, 3), q) +
     b0(p, tree.mu(2, 3), tree.mu(1, 3), q));	

  /// Other Higgs'
  double higgs = 0.;
  higgs = higgs +
    sqr(g) * 0.25 * 
    (2.0 * ffn(p, tree.mHpm, displayMwRun(), q) + 
     sqr(sin(tree.thetaH - beta)) / cw2DRbar *
     ffn(p, tree.mh0(2), mz, q) + 
     sqr(cos(tree.thetaH - beta)) / cw2DRbar * 
     ffn(p, tree.mh0(1), mz, q));
  
  /// trilinear Higgs coupling Feynman rules
  DoubleVector lAas(2), lAah(2);
  lAas(1) = -cos2b * cosb; lAas(2) = cos2b * sinb;
  lAas = lAas * (g * mz / (2.0 * cos(thetaWDRbar)));
  lAah = rot2d(tree.thetaH) * lAas;
  /// trilinear Higgs/Goldstone Feynman rules
  DoubleVector lAgs(2), lAgh(2);
  lAgs(1) = -sin2b * cosb; lAgs(2) = sin2b * sinb;
  lAgs = lAgs * (g * mz / (2.0 * cos(thetaWDRbar)));
  lAgh = rot2d(tree.thetaH) * lAgs;

  /// quartic Higgs coupling Ferynman rules: check
  DoubleMatrix lAa12(2, 2), lAahh(2, 2);
  lAa12(1, 1) = -cos(2.0 * beta);
  lAa12(2, 2) = cos(2.0 * beta);
  lAahh = rot2d(tree.thetaH) * lAa12 *
    rot2d(-tree.thetaH);
  lAahh = lAahh * sqr(g) / (4.0 * cw2DRbar);
  double lAaaa = sqr(g) / (4.0 * cw2DRbar) * 3.0 * sqr(cos2b),
    lAagg = (3.0 * sqr(sin2b) - 1.0) * sqr(g) / (4.0 * cw2DRbar);
  
  higgs = higgs +
    0.5 * (sqr(lAah(2)) * b0(p, mA, mh0, q) * 2 + 
	   sqr(lAah(1)) * b0(p, mA, mH, q) * 2 + 
	   sqr(lAgh(2)) * b0(p, mz, mh0, q) * 2 +
	   sqr(lAgh(1)) * b0(p, mz, mH, q) * 2 +
	   lAahh(1, 1) * a0(mH, q) + lAahh(2, 2) * a0(mh0, q) + 
	   lAaaa * a0(mA, q) + lAagg * a0(mz, q)) +             
    sqr(g) * sqr(displayMwRun()) * 0.5 * b0(p, displayMwRun(), mHc, q) + 
    sqr(g) / (4.0 * cw2DRbar)
    * ((cw2DRbar * (1.0 + sqr(sin2b)) - sw2DRbar * sqr(cos2b)) * 
       a0(displayMwRun(), q) + sqr(cos2b) * a0(mHc, q));    
  
  higgs = higgs + sqr(g) * 
    (2.0 * a0(displayMwRun(), q) + a0(mz, q) / cw2DRbar); 
  
  DoubleMatrix aPsiP1(4, 4), aPsiP2(4, 4), bPsiP1(4, 4), 
    bPsiP2(4, 4);
  aPsiP1(1, 3) = - gp * 0.5;   aPsiP1(3, 1) = - gp * 0.5; 
  aPsiP1(2, 3) = g * 0.5;      aPsiP1(3, 2) = g * 0.5; 
  aPsiP2(1, 4) = - gp * 0.5;   aPsiP2(4, 1) = - gp * 0.5; 
  aPsiP2(2, 4) = g * 0.5;      aPsiP2(4, 2) = g * 0.5; 
  
  bPsiP1 = -1.0 * aPsiP1; bPsiP2 = -1.0 * aPsiP2;

  ComplexMatrix aPsi(4, 4), bPsi(4, 4), aChi(4, 4), bChi(4, 4);
  ComplexMatrix n(tree.nBpmz);
  DoubleVector mneut(tree.mnBpmz);
  ComplexMatrix u(tree.uBpmz), v(tree.vBpmz); 
  DoubleVector mch(tree.mchBpmz); 
  
  /// Rotate to A neutralino basis
  ComplexMatrix aChiChiA(4, 4), bChiChiA(4, 4);
  aChiChiA = -sinb * n.complexConjugate() * aPsiP1 * n.hermitianConjugate() +
    cosb * n.complexConjugate() * aPsiP2 * n.hermitianConjugate();
  bChiChiA = -sinb * n * bPsiP1 * n.transpose() + cosb * n * bPsiP2 *
    n.transpose();
  
  double neutralinos = 0.0; 
  DoubleMatrix fChiChiA(4, 4), gChiChiA(4, 4);
  for (i=1; i<=4; i++)
    for (j=1; j<=4; j++) {
      fChiChiA(i, j) = sqr(aChiChiA(i, j).mod()) + sqr(bChiChiA(i, j).mod());
      gChiChiA(i, j) = 2.0 * (bChiChiA(i, j).conj() * aChiChiA(i, j)).real();
      neutralinos = neutralinos + 
	0.5 * fChiChiA(i, j) * 
	gfn(p, mneut(i), mneut(j), q) -
	gChiChiA(i, j) * mneut(i) * mneut(j) * 
	b0(p, mneut(i), mneut(j), q);
    }
  
  ComplexMatrix aChPsiP1(2, 2), aChPsiP2(2, 2), bChPsiP1(2, 2),
    bChPsiP2(2, 2), aChPsiA(2, 2), bChPsiA(2, 2); 
  aChPsiP1(1, 2) = g / root2;   aChPsiP2(2, 1) = - g / root2;
  bChPsiP1 = - 1.0 * aChPsiP1.transpose(); 
  bChPsiP2 = - 1.0 * aChPsiP2.transpose();
  
  /// To A-mass basis
  aChPsiA = -sinb * aChPsiP1 + cosb * aChPsiP2;
  bChPsiA = -sinb * bChPsiP1 + cosb * bChPsiP2;
  
  ComplexMatrix aChChiA(2, 2), bChChiA(2, 2);
  DoubleMatrix  fChChiA(2, 2), gChChiA(2, 2);
  aChChiA = v.complexConjugate() * aChPsiA * u.hermitianConjugate();
  bChChiA = u * bChPsiA * v.transpose();

  double charginos = 0.0;
  for (i=1; i<=2; i++)
    for (j=1; j<=2; j++) {
      fChChiA(i, j) = sqr(aChChiA(i, j).mod()) + sqr(bChChiA(i, j).mod());
      gChChiA(i, j) = 2.0 * (bChChiA(i, j).conj() * aChChiA(i, j)).real();
    
      charginos = charginos +
	fChChiA(i, j) * gfn(p, mch(i), mch(j), q) -
	2.0 * gChChiA(i, j) * mch(i) * mch(j) * b0(p, mch(i), mch(j), q);
    }
  ///  cout << "f=" << fChChiA << " g=" << gChChiA;

  fermions = fermions + sinb2 * 
    (sqr(forLoops.htau) * 
     (sqr(p) * b0(p, mtau, mtau, q) - 2.0 * a0(mtau, q)) +
     sqr(displayYukawaElement(YE, 2, 2)) * 
     (sqr(p) * b0(p, mmu, mmu, q) - 2.0 * a0(mmu, q)) +
     sqr(displayYukawaElement(YE, 1, 1)) * 
     (sqr(p) * b0(p, mE, mE, q) - 2.0 * a0(mE, q)));

  DoubleMatrix aee(2, 2); 
  aee(1, 2) = - forLoops.htau / root2 * 
    (-displaySusyMu() * sinb - forLoops.utau / forLoops.htau * cosb);
  aee(2, 1) = - aee(1, 2);
  aee = rot2d(thetatau) * aee * rot2d(-thetatau);

  /// Down fermions 
  fermions = fermions + sinb2 * 
    (3.0 * sqr(tree.hb) * 
     (sqr(p) * b0(p, mb, mb, q) - 2.0 * a0(mb, q)));
  
  DoubleMatrix add(2, 2); 
  add(1, 2) = - tree.hb / root2 * 
    (-displaySusyMu() * sinb - tree.ub / tree.hb * cosb);
  add(2, 1) = - add(1, 2);
  add = rot2d(thetab) * add * rot2d(-thetab);

  /// sneutrinos
  double sneutrinos = -
    sqr(g) / cw2DRbar * 0.5 * gnuL * cos2b * 
    (a0(tree.msnu(1), q) + a0(tree.msnu(2), q) +
     a0(tree.msnu(3), q));

  /// LH gens 1-2
  double sleptons = 0.;
  sleptons = -(sqr(g) / (cw2DRbar) * cos2b * 0.5 * geL) *
    (a0(tree.me(1, 1), q) + a0(tree.me(1, 2), q));

  /// RH gens 1-2
  sleptons = sleptons -
    (sqr(g) / (cw2DRbar) * cos2b * 0.5 * geR) *
    (a0(tree.me(2, 1), q) + a0(tree.me(2, 2), q));

  double staus = 
    (sqr(tree.htau) * sinb2 -
	   sqr(g) / (cw2DRbar) * cos2b * 0.5 * geL) *
    (sqr(ctau) * a0(tree.me(1, 3), q) + 
     sqr(stau) * a0(tree.me(2, 3), q));
  
  staus = staus +
    (sqr(tree.htau) * sinb2 -
	   sqr(g) / (cw2DRbar) * cos2b * 0.5 * geR) *
    (sqr(stau) * a0(tree.me(1, 3), q) + 
     sqr(ctau) * a0(tree.me(2, 3), q));
  
  sleptons = sleptons +
    (sqr(displayYukawaElement(YE, 1, 1) * 
	       (-displaySusyMu() * cosb - 
		displaySoftA(EA, 1, 1) * sinb))) * 0.5 * 
    (b0(p, tree.me(1, 1), tree.me(2, 1), q) +
     b0(p, tree.me(2, 1), tree.me(1, 1), q)) + 
    (sqr(displayYukawaElement(YE, 2, 2) * 
	       (-displaySusyMu() * cosb - 
		displaySoftA(EA, 2, 2) * sinb))) * 0.5 * 
    (b0(p, tree.me(1, 2), tree.me(2, 2), q) +
     b0(p, tree.me(2, 2), tree.me(1, 2), q));


  double aee12 = (tree.htau * 
    displaySusyMu() * cosb + tree.utau * sinb) / root2;

  staus = staus +  
    sqr(aee12) * 
    (b0(p, tree.me(1, 3), tree.me(2, 3), q) +
     b0(p, tree.me(2, 3), tree.me(1, 3), q));
  
  /// Sign of certain contributions bug-fixed 30-07-03
  double dsfermions = 0.0;
  dsfermions = dsfermions -
    3.0 * (sqr(g) / (cw2DRbar) * cos2b * 0.5 * gdL) *
    (a0(tree.md(1, 1), q) + a0(tree.md(1, 2), q));

  double sbots =  
    3.0 * (sqr(tree.hb) * sinb2 -
	   sqr(g) / (cw2DRbar) * cos2b * 0.5 * gdL) *
    (sqr(cb) * a0(tree.md(1, 3), q) + 
     sqr(sb) * a0(tree.md(2, 3), q));

  dsfermions = dsfermions -
    3.0 * (sqr(g) / (cw2DRbar) * cos2b * 0.5 * gdR) *
    (a0(tree.md(2, 1), q) + a0(tree.md(2, 2), q));
  
  sbots = sbots +
    3.0 * (sqr(tree.hb) * sinb2 -
	   sqr(g) / (cw2DRbar) * cos2b * 0.5 * gdR) *
    (sqr(sb) * a0(tree.md(1, 3), q) + 
     sqr(cb) * a0(tree.md(2, 3), q));

  dsfermions = dsfermions +
    3.0 * (sqr(displayYukawaElement(YD, 1, 1) * 
	       (-displaySusyMu() * cosb - 
		displaySoftA(DA, 1, 1) * sinb))) * 0.5 * 
    (b0(p, tree.md(1, 1), tree.md(2, 1), q) +
     b0(p, tree.md(2, 1), tree.md(1, 1), q)) + 
    3.0 * (sqr(displayYukawaElement(YD, 2, 2) * 
	       (-displaySusyMu() * cosb - 
		displaySoftA(DA, 2, 2) * sinb))) * 0.5 * 
    (b0(p, tree.md(1, 2), tree.md(2, 2), q) +
     b0(p, tree.md(2, 2), tree.md(1, 2), q));

  double add12 = (tree.hb *
    displaySusyMu() * cosb + tree.ub * sinb) / root2;

  sbots = sbots +  
	3.0 * sqr(add12) * 
	(b0(p, tree.md(1, 3), tree.md(2, 3), q) +
	 b0(p, tree.md(2, 3), tree.md(1, 3), q));

  double sfermions = sleptons + staus + upsfermions + dsfermions + 
    sneutrinos + stops + sbots; 

  return (fermions + sfermions + higgs + neutralinos + charginos) 
    / (16.0 * sqr(PI));
}


template<class SoftPars>
double Softsusy<SoftPars>::piZGT(double p, double q) const { ///! checked 7/6/6
  drBarPars tree(displayDrBarPars());
  double alphaMsbar = dataSet.displayAlpha(legacy::ALPHA);
  double alphaDrbar = qedSusythresh(alphaMsbar, displayMu());

  double    thetaWDRbar = asin(calcSinthdrbar());
  double    cw2DRbar    = sqr(cos(thetaWDRbar));
  double    sw2DRbar    = 1.0 - cw2DRbar;
  double    mHc = tree.mHpm;
  double    mtop = displayDataSet().displayPoleMt();
  double    mb   =  tree.mb;
  double    mtau =  tree.mtau;
  double    ms   =  displayDataSet().displayMass(legacy::mStrange) ;
  double    mc   =  displayDataSet().displayMass(legacy::mCharm) ;
  double    mmu  =  displayDataSet().displayMass(legacy::mMuon) ;
  double    mE  =  displayDataSet().displayMass(legacy::mElectron) ;
  double    mD  =  displayDataSet().displayMass(legacy::mDown) ;
  double    mU  =  displayDataSet().displayMass(legacy::mUp);
  double    thetat = tree.thetat ;
  double    thetab = tree.thetab;
  double    thetatau= tree.thetatau ;
  double    ct      = cos(thetat) ;
  double    cb      = cos(thetab) ;
  double    ctau    = cos(thetatau);
  double    ct2     = sqr(ct) ;
  double    cb2     = sqr(cb) ;
  double    ctau2   = sqr(ctau) ;
  double    st2     = (1.0 - ct2);
  double    sb2     = (1.0 - cb2);
  double    stau2   = (1.0 - ctau2);
  double    msbot1  = tree.md(1, 3);
  double    msbot2  = tree.md(2, 3);
  double    mstau1  = tree.me(1, 3);
  double    mstau2  = tree.me(2, 3);
  double    mstop1  = tree.mu(1, 3);
  double    mstop2  = tree.mu(2, 3);
  double    g       = displayGaugeCoupling(2);

  double ans = 0.0;
  
  ans = ans + 
    (12.0 * sw2DRbar - 10.0) * b22bar(p, displayMwRun(), displayMwRun(), q) - 
    2.0 * (sqr(displayMwRun()) + 
	   2.0 * cw2DRbar * sqr(p)) * b0(p, displayMwRun(), displayMwRun(), q);
  
  ans = ans + 2.0 * (guL - guR) * (4.0 * b22bar(p, mU, mU, q) +
				   sqr(p) * b0(p, mU, mU, q));
  ans = ans + 2.0 * (guL - guR) * (4.0 * b22bar(p, mc, mc, q) +
				   sqr(p) * b0(p, mc, mc, q));
  ans = ans + 2.0 * (guL - guR) * (4.0 * b22bar(p, mtop, mtop, q) +
				   sqr(p) * b0(p, mtop, mtop, q));
  ans = ans - (gdL - gdR) * (4.0 * b22bar(p, mD, mD, q) +
			     sqr(p) * b0(p, mD, mD, q));
  ans = ans - (gdL - gdR) * (4.0 * b22bar(p, ms, ms, q) +
			     sqr(p) * b0(p, ms, ms, q));
  ans = ans - (gdL - gdR) * (4.0 * b22bar(p, mb, mb, q) +
			     sqr(p) * b0(p, mb, mb, q));
  ans = ans - (geL - geR) * (4.0 * b22bar(p, mE, mE, q) +
			     sqr(p) * b0(p, mE, mE, q));
  ans = ans - (geL - geR) * (4.0 * b22bar(p, mmu, mmu, q) +
			     sqr(p) * b0(p, mmu, mmu, q));
  ans = ans - (geL - geR) * (4.0 * b22bar(p, mtau, mtau, q) +
			     sqr(p) * b0(p, mtau, mtau, q));
  
  ans = ans - 2.0 * cos(2.0 * thetaWDRbar) * b22bar(p, mHc, mHc, q);
  
  ComplexMatrix u(tree.uBpmz), v(tree.vBpmz); 
  DoubleVector mch(tree.mchBpmz); 
  
  int i; for (i=1; i<=2; i++) {
    ans = ans + 0.5 * (sqr(v(i, 1).mod()) + sqr(u(i, 1).mod()) 
		       + 2.0 * cos(2.0 * thetaWDRbar))
      * (4.0 * b22bar(p, mch(i), mch(i), q) +
	 sqr(p) * b0(p, mch(i), mch(i), q));
    
    ans = ans - 
      8.0 * (guL * b22bar(p, tree.mu(1, i), tree.mu(1, i), q) -
	     guR * b22bar(p, tree.mu(2, i), tree.mu(2, i), q)) +
      4.0 * (gdL * b22bar(p, tree.md(1, i), tree.md(1, i), q) -
	     gdR * b22bar(p, tree.md(2, i), tree.md(2, i), q)) +
      4.0 * (geL * b22bar(p, tree.me(1, i), tree.me(1, i), q) -
	     geR * b22bar(p, tree.me(2, i), tree.me(2, i), q));
  }
  
  ans = ans -  8.0 * 
    ((guL * ct2 - guR * st2) * b22bar(p, mstop1, mstop1, q) + 
     (guL * st2 - guR * ct2) * b22bar(p, mstop2, mstop2, q));
  ans = ans +  4.0 * 
    ((gdL * cb2 - gdR * sb2) * b22bar(p, msbot1, msbot1, q) + 
     (gdL * sb2 - gdR * cb2) * b22bar(p, msbot2, msbot2, q));
  ans = ans +  4.0 * 
    ((geL * ctau2 - geR * stau2) * b22bar(p, mstau1, mstau1, q) + 
     (geL * stau2 - geR * ctau2) * b22bar(p, mstau2, mstau2, q));
  
  double e = sqrt(alphaDrbar * 4.0 * PI);
  
  return ans * e * g / (16.0 * sqr(PI) * cos(thetaWDRbar)); 
} 

inline double rho2(double r) {/// checked 
  if (r <= 1.9)
    return 19.0 - 16.5 * r + 43.0 * sqr(r) / 12.0 + 7.0 / 120.0 * sqr(r) * r -
      PI * sqrt(r) * (4.0 - 1.5 * r + 3.0 / 32.0 * sqr(r) + sqr(r) * r /
		      256.0) - sqr(PI) * (2.0 - 2.0 * r + 0.5 * sqr(r)) -
      log(r) * (3.0 * r - 0.5 * sqr(r));
  else {
    double rm1 = 1.0 / r, rm2 = sqr(rm1), rm3 = rm2 * rm1, rm4 = rm3 * rm1,
      rm5 = rm4 * rm1;
    return sqr(log(r)) * (1.5 - 9.0 * rm1 - 15.0 * rm2 - 48.0 * rm3 - 168.0
			  * rm4 - 612.0 * rm5) -
      log(r) * (13.5 + 4.0 * rm1 - 125.0 / 4.0 * rm2 - 558.0 / 5.0 * rm3 -
		8307.0 / 20.0 * rm4 - 109321.0 / 70.0 * rm5)
      + sqr(PI) * (1.0 - 4.0 * rm1 - 5.0 * rm2 - 16.0 * rm3 -
		   56.0 * rm4 - 204.0 * rm5)
      + 49.0 / 4.0 + 2.0 / 3.0 * rm1 + 1613.0 / 48.0 * rm2 + 87.57 * rm3 +
      341959.0 / 1200.0 * rm4 + 9737663.0 / 9800.0 * rm5;
  }
}

inline double fEff(double x) {
  double arg = 1.0 / (1.0 + x);
  return 2.0 / x + 3.5 - (3.0 + 2.0 / x) * log(x) +
    sqr(1.0 + 1.0 / x) * 
    (2.0 * dilog(arg) - sqr(PI) / 3.0 + sqr(log(1.0 + x)));
}

inline double gEff(double x) {
  double y = sqrt(x / (4.0 - x));
  return (1.0 / x + 0.5) * (atan(y) / y - 1.0) + 9.0 / 8.0 + 0.5 / x -
    (1.0 + 0.5 / x) * 4.0 / x * sqr(atan(y));
}


template<class SoftPars>
double Softsusy<SoftPars>::sinSqThetaEff() {
  if (displayMu() != MZ) {
    throw("Should call Softsusy<SoftPars>::sinSqThetaEff() at MZ only\n");
  }
  double kl = 0.;
 
  calcDrBarPars();

  double alphaMsbar = dataSet.displayAlpha(legacy::ALPHA);
  double alphaDrbar = qedSusythresh(alphaMsbar, displayMu());
  double sinthDrbar = calcSinthdrbar();
  double costhDrbar = sqrt(1.0 - sqr(sinthDrbar));
  double costhW     = MW / MZ;

  double vLmzSq = 0.5 * fEff(1.0 / sqr(costhW)) + 
    4.0 * sqr(costhDrbar) * gEff(1.0 / sqr(costhW)) -
    (1.0 - 6.0 * sqr(sinthDrbar) + 8.0 * sqr(sinthDrbar) * sqr(sinthDrbar)) /
    (4.0 * sqr(costhDrbar)) * fEff(1.0);

  kl = 1.0 + costhDrbar / sinthDrbar * 
    (piZGT(MZ, displayMu()) - piZGT(0., displayMu())) / sqr(MZ) +
     alphaDrbar / PI * sqr(costhDrbar) / sqr(sinthDrbar) * 2.0 * log(costhW)
     - alphaDrbar / (4.0 * PI * sqr(sinthDrbar)) * vLmzSq;

  return sqr(sinthDrbar) * kl;
}

/// outrho, outsin represent the DRbar values
template<class SoftPars>
double Softsusy<SoftPars>::deltaVb(double outrho, double outsin, 
                                   double alphaDRbar, double /* pizztMZ */) const {
  drBarPars tree(displayDrBarPars());

  double g       = displayGaugeCoupling(2);
  double gp      = displayGaugeCoupling(1) * sqrt(0.6);
  double costh   = (displayMw() / displayMz());
  double cw2     = sqr(costh) ;
  double sw2     = (1.0 - cw2);
  double outcos  = sqrt(1.0 - sqr(outsin));
  double q       = displayMu();

  ComplexMatrix n(tree.nBpmz);
  DoubleVector mneut(tree.mnBpmz);
  //PA: get the dimension of menut
  const int dimN =  mneut.displayEnd();
  ComplexMatrix u(tree.uBpmz), v(tree.vBpmz); 
  DoubleVector mch(tree.mchBpmz); 

  double deltaVbSm = outrho * alphaDRbar / (4.0 * PI * sqr(outsin)) *
    (6.0 + log(cw2) / sw2 * 
     (3.5 - 2.5 * sw2 - sqr(outsin) * (5.0 - 1.5 * cw2 / sqr(outcos)))
     - 4. * log(sqr(displayMz()/q)));
  
  DoubleVector bPsi0NuNul(dimN), bPsicNuSell(2);
  DoubleVector bPsi0ESell(dimN), aPsicESnul(2);
  ComplexVector bChi0NuNul(dimN), bChicNuSell(2);
  ComplexVector bChi0ESell(dimN), aChicESnul(2);
  
  bPsicNuSell(1) = g;
  bPsi0NuNul(2) = root2 * g * 0.5;
  bPsi0NuNul(1) = - gp / root2;
  aPsicESnul(1) = g;
  bPsi0ESell(1) = -gp / root2;
  bPsi0ESell(2) = -g * root2 * 0.5;
  
  bChicNuSell = u * bPsicNuSell;
  bChi0ESell =  n * bPsi0ESell;
  bChi0NuNul = n * bPsi0NuNul;

  aChicESnul = v.complexConjugate() * aPsicESnul;
  
  double deltaZnue = 0.0, deltaZe = 0.0;
  double mselL = tree.me(1, 1), 
    msnue = tree.msnu(1);
  int i; for(i=1; i<=dimN; i++) {
   if (i < 3) {
      deltaZnue = deltaZnue -
	sqr(bChicNuSell(i).mod()) * b1(0.0, mch(i), mselL, q);
      deltaZe = deltaZe -
	sqr(aChicESnul(i).mod()) * b1(0.0, mch(i), msnue, q);
    }
    deltaZnue = deltaZnue -
      sqr(bChi0NuNul(i).mod()) * b1(0.0, mneut(i), msnue, q);
    deltaZe = deltaZe -
      sqr(bChi0ESell(i).mod()) * b1(0.0, mneut(i), mselL, q);
  }
  
  DoubleVector bPsicNuSmul(2);
  DoubleVector bPsi0MuSmul(dimN), aPsicMuSnul(2);
  ComplexVector bChicNuSmul(2);
  ComplexVector bChi0MuSmul(dimN), aChicMuSnul(2);
  
  double hmu = displayYukawaElement(YE, 2, 2);
  bPsicNuSmul(1) = g;
  bPsicNuSmul(2) = 0.;
  aPsicMuSnul(1) = g;
  aPsicMuSnul(2) = 0.;
  bPsi0MuSmul(1) = -gp / root2;
  bPsi0MuSmul(2) = -g * root2 * 0.5;

  bChicNuSmul = u * bPsicNuSmul;
  bChi0MuSmul =  n * bPsi0MuSmul;
  bChi0NuNul = n * bPsi0NuNul;
  aChicMuSnul = v.complexConjugate() * aPsicMuSnul;
  
  double deltaZnumu = 0.0, deltaZmu = 0.0;
  double msnumu = tree.msnu(2),
    msmuL = tree.me(1, 2);
  for(i=1; i<=dimN; i++) {
    if (i < 3) {
      deltaZnumu = deltaZnumu -
	sqr(bChicNuSmul(i).mod()) * b1(0.0, mch(i), msmuL, q);
      deltaZmu = deltaZmu -
	sqr(aChicMuSnul(i).mod()) * b1(0.0, mch(i), msnumu, q);
    }
    deltaZnumu = deltaZnumu -
      sqr(bChi0NuNul(i).mod()) * b1(0.0, mneut(i), msnumu, q);
    deltaZmu = deltaZmu -
      sqr(bChi0MuSmul(i).mod()) * b1(0.0, mneut(i), msmuL, q);
  }
  
  DoubleMatrix aPsi0PsicW(dimN, 2), bPsi0PsicW(dimN, 2), fW(dimN, 2), gW(dimN, 2);
  ComplexMatrix aChi0ChicW(dimN, 2), bChi0ChicW(dimN, 2);
  
  aPsi0PsicW(2, 1) = - g;
  bPsi0PsicW(2, 1) = - g;
  aPsi0PsicW(4, 2) = g / root2;		     
  bPsi0PsicW(3, 2) = -g / root2;		     
  
  /// These ought to be in physpars
  aChi0ChicW = n.complexConjugate() * aPsi0PsicW * v.transpose();
  bChi0ChicW = n * bPsi0PsicW * u.hermitianConjugate();
  
  Complex deltaVE = 0.0;
  int j; for(i=1; i<=2; i++)
    for(j=1; j<=dimN; j++) {
      deltaVE = deltaVE + bChicNuSell(i) * bChi0ESell(j).conj() *
	(- root2 / g * aChi0ChicW(j, i) * mch(i) * mneut(j) *
	 c0(mselL, mch(i), mneut(j)) + 1.0 / (root2 * g) *
	 bChi0ChicW(j, i) *
	 (b0(0.0, mch(i), mneut(j), q) + sqr(mselL) * 
	  c0(mselL, mch(i), mneut(j)) - 0.5));
      deltaVE = deltaVE - aChicESnul(i) * bChi0NuNul(j) *
	(- root2 / g * bChi0ChicW(j, i) * mch(i) * mneut(j) *
	 c0(msnue, mch(i), mneut(j)) + 1.0 / (root2 * g) *
	 aChi0ChicW(j, i) * 
	 (b0(0.0, mch(i), mneut(j), q) + sqr(msnue) * 
	  c0(msnue, mch(i), mneut(j)) - 0.5));
      if (i == 1)
	deltaVE = deltaVE +
	  0.5 * bChi0ESell(j).conj() * bChi0NuNul(j) *
	  (b0(0.0, mselL, msnue, q) + sqr(mneut(j)) * 
	   c0(mneut(j), mselL, msnue) + 0.5);
    }
  
  Complex deltaVMu = 0.0;
  for(i=1; i<=2; i++)
    for(j=1; j<=dimN; j++) {
      deltaVMu = deltaVMu + bChicNuSmul(i) * bChi0MuSmul(j).conj() *
	(- root2 / g * aChi0ChicW(j, i) * mch(i) * mneut(j) *
	 c0(msmuL, mch(i), mneut(j)) + 1.0 / (root2 * g) *
	 bChi0ChicW(j, i) *
	 (b0(0.0, mch(i), mneut(j), q) + sqr(msmuL) * 
	  c0(msmuL, mch(i), mneut(j)) - 0.5));
      deltaVMu = deltaVMu - aChicMuSnul(i) * bChi0NuNul(j) *
	(- root2 / g * bChi0ChicW(j, i) * mch(i) * mneut(j) *
	 c0(msnumu, mch(i), mneut(j)) + 1.0 / (root2 * g) *
	 aChi0ChicW(j, i) * 
	 (b0(0.0, mch(i), mneut(j), q) + sqr(msnumu) * 
	  c0(msnumu, mch(i), mneut(j)) - 0.5));
      if (i == 1)
	deltaVMu = deltaVMu +
	  0.5 * bChi0MuSmul(j).conj() * bChi0NuNul(j) *
	  (b0(0.0, msmuL, msnumu, q) + sqr(mneut(j)) * c0(mneut(j), msmuL,
							  msnumu) + 0.5);
    }
  
  Complex a1(0.0, 0.0);
  for(i=1; i<=2; i++)
    for(j=1; j<=dimN; j++) {
      a1 = a1 + 0.5 * aChicMuSnul(i) * bChicNuSell(i).conj() *
	bChi0NuNul(j) * bChi0ESell(j) * mch(i) * mneut(j) * 
	d0(mselL, msnumu, mch(i), mneut(j));
      a1 = a1 + 0.5 * aChicESnul(i).conj() * bChicNuSmul(i) *
	bChi0NuNul(j).conj() * bChi0MuSmul(j).conj() * mch(i) * mneut(j) * 
	d0(msmuL, msnue, mch(i), mneut(j));
      a1 = a1 + bChicNuSmul(i) * bChicNuSell(i).conj() *
	bChi0MuSmul(j).conj() * bChi0ESell(j) * 
	d27(msmuL, mselL, mch(i), mneut(j));
      a1 = a1 + aChicMuSnul(i).conj() * aChicESnul(i) *
	bChi0NuNul(j) * bChi0NuNul(j).conj() * 
	d27(msnumu, msnue, mch(i), mneut(j));
    } 

  double deltaVbSusy = 
    (-sqr(outsin) * sqr(outcos) / (2.0 * PI * alphaDRbar) * sqr(displayMz())
     * a1.real() + deltaVE.real() + deltaVMu.real() + 
     0.5 * (deltaZe + deltaZnue + deltaZmu + deltaZnumu) ) /
    (sqr(PI) * 16.0); 

  double deltaVb;
  deltaVb = deltaVbSusy + deltaVbSm;

  return deltaVb;
}

//PA: returns the mixing of Hu into between h1 
template<class SoftPars>
double Softsusy<SoftPars>::h1s2Mix(){					    
   return cos(displayDrBarPars().thetaH);
}

template<class SoftPars>
double Softsusy<SoftPars>::dRho(double outrho, double outsin, double alphaDRbar, 
			     double pizztMZ, double piwwtMW) {
  double mz = displayMz();
  drBarPars tree(displayDrBarPars());

  /// 2 loop SM contribution
  double mt   = dataSet.displayPoleMt(); 
  double sinb = sin(atan(displayTanb()));
  
  double xt = 3.0 * GMU * sqr(mt) / (8.0 * sqr(PI) * root2);
  
  double deltaRho2LoopSm = alphaDRbar * sqr(displayGaugeCoupling(3)) / 
    (16.0 * PI * sqr(PI) * sqr(outsin)) * /// bug-fixed 24.08.2002
    (-2.145 * sqr(mt) / sqr(displayMw()) + 1.262 * log(mt / mz) - 2.24 
     - 0.85 * sqr(mz)
     / sqr(mt)) + sqr(xt) * sqr(h1s2Mix()) / sqr(sinb) *
    rho2(tree.mh0(1) / mt) / 3.0;

  double deltaRhoOneLoop = pizztMZ / (outrho * sqr(mz))
    - piwwtMW / sqr(displayMw());
  
  double deltaRho = deltaRhoOneLoop + deltaRho2LoopSm;

  return deltaRho;
}

template<class SoftPars>
double Softsusy<SoftPars>::dR(double outrho, double outsin, double alphaDRbar,
			   double pizztMZ, double piwwt0) {
  drBarPars tree(displayDrBarPars());

  double outcos = cos(asin(outsin));
  /// 2 loop SM contribution
  double mt   = dataSet.displayPoleMt();

  double sinb = sin(atan(displayTanb()));
  
  double xt = 3.0 * GMU * sqr(mt) / (8.0 * sqr(PI) * root2);
  
  double dvb = deltaVb(outrho, outsin, alphaDRbar, pizztMZ);

  double mz = displayMz();

  double deltaR =  outrho * piwwt0 / sqr(displayMw()) - 
    pizztMZ / sqr(mz) + dvb;

  /// Dominant two-loop SM term
  double deltaR2LoopSm = alphaDRbar * sqr(displayGaugeCoupling(3)) / 
    (16.0 * sqr(PI) * PI * sqr(outsin) * sqr(outcos)) *
    (2.145 * sqr(mt) / sqr(mz) + 0.575 * log(mt / mz) - 0.224 
     - 0.144 * sqr(mz) / sqr(mt)) - 
    sqr(xt) * sqr(h1s2Mix()) / sqr(sinb) *
    rho2(tree.mh0(1) / mt) * (1.0 - deltaR) * outrho / 3.0;

  deltaR = deltaR + deltaR2LoopSm; 

  return deltaR;
}

/// Checked 20.11.00
/// Flags noconvergence if there's trouble...then don't believe outrho and
/// outsin produced - they are fudged!
template<class SoftPars>
void Softsusy<SoftPars>::rhohat(double & outrho, double & outsin, double alphaDRbar,
			  double pizztMZ, double piwwt0, double piwwtMW, 
			  double tol, int maxTries) {

  static double oldrho = 0.23, oldsin = 0.8;

  double mz = displayMz();
  if (displayMu() != mz) {
    std::ostringstream ii;   
    ii << "Called Softsusy<SoftPars>::rhohat "
       << "with scale" << displayMu() << endl;
    throw ii.str();
  }
  
  static int numTries = 0;
  
  if ((outrho < TOLERANCE || outsin < TOLERANCE) || fabs(outsin) > 1. ||
      (numTries - 1 > maxTries)) {  
    oldrho = 0.23; oldsin = 0.8;
    numTries = 0;
    outrho = 0.23; outsin = 0.8; 
    flagNoRhoConvergence(true);
    if (PRINTOUT) cout << std::flush << "rhohat reached maxtries\n"; 
    return;
  }
  
  /// Difference to last iteration
  double sumTol;
  sumTol = fabs(oldrho / outrho - 1.0) + fabs(oldsin / outsin - 1.0);
  
  if (numTries != 0 && sumTol < tol) {
    numTries = 0;
    oldrho = 0.23, oldsin = 0.8;
    if (PRINTOUT > 2) cout << "sin rho converged\n";
    return;
  }

  numTries = numTries + 1;
  
  oldrho = outrho; oldsin = outsin; 
    
  double deltaR = dR(outrho, outsin, alphaDRbar, pizztMZ, piwwt0);

  if (PRINTOUT > 2) cout << " st2=" << sumTol << " dr=" << deltaR 
			 << " outrho=" << outrho << " outsin=" << outsin 
			 << " aDRbar=" << alphaDRbar << " piZ=" << pizztMZ 
			 << " piwwt0=" << piwwt0 << endl;
  
  double sin2thetasqO4 = PI * alphaDRbar / 
    (root2 * sqr(mz) * GMU * (1.0 - deltaR)); 

  if (sin2thetasqO4 >= 0.25) sin2thetasqO4 = 0.25;
  if (sin2thetasqO4 < 0.0) sin2thetasqO4 = 0.0;

  double sin2theta = sqrt(4.0 * sin2thetasqO4);

  double theta = 0.5 * asin(sin2theta);
  
  outsin = sin(theta); 

  double deltaRho = dRho(outrho, outsin, alphaDRbar, pizztMZ, piwwtMW);

  if (fabs(deltaRho) < 1.0) outrho = 1.0 / (1.0 - deltaRho);
  else outrho = 1.0;

  if (PRINTOUT > 2) cout << " drho=" << deltaRho << " sw=" << outsin << endl; 

  rhohat(outrho, outsin, alphaDRbar, pizztMZ, piwwt0, piwwtMW, tol, maxTries);
}

template<class SoftPars>
void Softsusy<SoftPars>::methodBoundaryCondition(const DoubleVector & /* pars */) {
  std::ostringstream ii;
  ii << "Should only use Softsusy<SoftPars>::methodBoundaryCondition in derived"
     << " objects.\n";
  throw ii.str();
}

template<class SoftPars>
void Softsusy<SoftPars>::rpvSet(const DoubleVector & /* parameters */) {
  std::ostringstream ii;
  ii << "Should only use Softsusy<SoftPars>::rpvSet in derived"
     << " objects.\n";
  throw ii.str();
}

/// isajet routine for getting alpha_s - only used here for getting numbers
/// from isawig interface
inline double sualfs(double qsq, double alam4, double tmass) {
  double anf = 6.0, bmass = 5.0, alam = 0.1, alam5 = 0.1, b0 = 0.;
  double alamsq = 0., b1 = 0., b2 = 0., x = 0., tt = 0., t = 0., sualfs = 0.;
  if (qsq < 4.0 * sqr(bmass)) {
    anf   = 4.0;
    alam  = alam4;
  }
  else if (qsq < 4.0 * sqr(tmass)) {
    anf   = 5.0;
    alam  = alam4 * pow(alam4 / (2.0*bmass), 2.0 / 23.0)
      * pow(log(4.0 * sqr(bmass) / sqr(alam4)), -963.0 / 13225.0);
  }
  else {
    anf   = 6.0;
        alam5 = alam4 * pow(alam4/(2.0*bmass), 2.0/23.0)
	  * pow(log(4.0*sqr(bmass)/sqr(alam4)), -963.0/13225.0);
	alam  = alam5 * pow(alam5/(2.0*tmass), 2.0/21.0)
	  * pow(log(4.0*sqr(tmass)/sqr(alam5)), -107.0/1127.0);
  }
  b0       = 11.0-2.0/3.0*anf;
  alamsq   = sqr(alam);
  t        = log(qsq/alamsq);
  if (t <= 1.0) t = log(4.0/alamsq);
  double alphas   = 4*PI/b0/t;
  b1 = 102.0-38.0/3.0*anf;
  b2 = 0.5*(2857.0-5033.0/9.0*anf+325.0/27.0*sqr(anf));
  x  = b1/(sqr(b0)*t);
  tt = log(t);
  sualfs = alphas*(1.0-x*tt+sqr(x)*(sqr(tt-0.5)
				  +b2*b0/sqr(b1)-1.25));
  return sualfs;
}

/// another function which mimics Isajet's handling of quark masses
inline double ssmqcd(double dm, double dq, double mtopPole) {      
  double dlam4=0.177;
  double dqbt=10.;
  double dqtp=2*mtopPole;
  double ssmqcd=0., renorm = 0.;
  int dneff=4;
  double po=12.0/(33.0-2.*dneff);
  if(dq < dqbt) {
    renorm=pow(log(2*dm/dlam4)/log(dq/dlam4), po);
    ssmqcd=renorm*dm;
    return ssmqcd;
  }
  else
    renorm=pow(log(2*dm/dlam4)/log(dqbt/dlam4), po);
  
  dneff=5;
  po=12.0/(33.0-2.*dneff);
  double dlam5=exp((25.0*log(dlam4)-log(sqr(dqbt)))/23.0);
  if(dq >= dqbt && dq < dqtp) {
    renorm=renorm *pow(log(dqbt/dlam5)/log(dq/dlam5), po);
    ssmqcd=renorm*dm;
    return ssmqcd;
  }
  else {
    renorm=renorm
      *pow(log(dqbt/dlam5)/log(dqtp/dlam5), po);
  }

  dneff=6;
  po=12.0/(33.0-2.*dneff);
  double dlam6=exp((25.0*log(dlam4)-log(sqr(dqbt))
	      -log(4*sqr(mtopPole)))/21.0);
  renorm=renorm
    *pow(log(dqtp/dlam6)/log(dq/dlam6), po);
  ssmqcd=renorm*dm;
  return ssmqcd;
}

/// Works out how best to fit the isajet numbers to the spectrum.
/// There are problems with the Higgs and sbottoms because ISAJET assumes
/// certain tree-level relations between masses that are broken by SOFTSUSY's
/// higher accuracy. The differences get large for high tan beta around 50, at
/// around 10 they're typically only a percent.
template<class SoftPars>
void Softsusy<SoftPars>::isajetNumbers764 
(double & mtopPole, double & mGPole, double & smu, double & mA, double & tanb, 
 double & mq1l, double & mdr, double & mur, double & meL, double & meR, 
 double & mql3, double & mdr3, double & mur3, double &  mtauL, 
 double & mtauR, double & at, double & ab, double & atau, double & mq2l, 
 double & msr, double & mcr, double & mmuL, double & mmuR, double & m1, 
 double & m2) 
  const {

  /// Store a copy of the current object
  Softsusy<SoftPars> store(*this);

  /// Run to MSUSY 
  store.runto(displayMsusy(), TOLERANCE * 0.01);

  DoubleMatrix mNeut(4, 4);
  
  /// Store current object in vector
  DoubleVector storeObject(store.display());
  double muInitial = store.displayMu();

  m1 = store.displayGaugino(1);
  m2 = store.displayGaugino(2); 
  smu = store.displaySusyMu();

  /// restore object at MSUSY
  store.setMu(muInitial);
  store.set(storeObject);

  /// tree level
  mNeut(1, 1) = m1;
  mNeut(2, 2) = m2;
  mNeut(3, 4) = - smu;
  store.addNeutralinoLoop(fabs(m1), mNeut);

  m1 = fabs(mNeut(1, 1));
  m2 = fabs(mNeut(2, 2));
  smu = -mNeut(3, 4);

  double beta = atan(store.displayTanb());

  /// define ISAJET constants
  double amdn=0.0099;
  double  amup=0.0056;
  double  amst=0.199;
  double  amch=1.35;
  double  ambt=5.0;
  double  ame=0.511e-3;
  double  ammu=0.105;
  double  amz=91.17;
  double  amw = 80.0;
  double  sn2thw=0.232;

  mq1l = sqrt(sqr(displayPhys().mu.display(1, 1)) - sqr(amup) -
	      (0.5 - 2./3. * sn2thw) * sqr(amz) * cos(2.0 * beta));
  mq2l = sqrt(sqr(displayPhys().mu.display(1, 2)) - sqr(amch) -
	      (0.5 - 2./3. * sn2thw) * sqr(amz) * cos(2.0 * beta));

  mur = sqrt(sqr(displayPhys().mu.display(2, 1)) - sqr(amup) -
	      (2./3. * sn2thw) * sqr(amz) * cos(2.0 * beta));
  mcr = sqrt(sqr(displayPhys().mu.display(2, 2)) - sqr(amch) -
	      (2./3. * sn2thw) * sqr(amz) * cos(2.0 * beta));

  mdr = sqrt(sqr(displayPhys().md.display(2, 1)) - sqr(amdn) -
	      (- 1./3. * sn2thw) * sqr(amz) * cos(2.0 * beta));
  msr = sqrt(sqr(displayPhys().md.display(2, 2)) - sqr(amst) -
	      (- 1./3. * sn2thw) * sqr(amz) * cos(2.0 * beta));

  meL = sqrt(sqr(displayPhys().me.display(1, 1)) - sqr(ame) +
	      (0.5 - sn2thw) * sqr(amz) * cos(2.0 * beta));
  mmuL = sqrt(sqr(displayPhys().me.display(1, 2)) - sqr(ammu) +
	      (0.5 - sn2thw) * sqr(amz) * cos(2.0 * beta));

  meR = sqrt(sqr(displayPhys().me.display(2, 1)) - sqr(ame) -
	      (- sn2thw) * sqr(amz) * cos(2.0 * beta));
  mmuR = sqrt(sqr(displayPhys().me.display(2, 2)) - sqr(ammu) -
	      (- sn2thw) * sqr(amz) * cos(2.0 * beta));

  mtopPole = store.displayDataSet().displayPoleMt();

  double asmt = sualfs(sqr(mtopPole), 0.36, mtopPole);
  double mtmt = mtopPole/(1.+4*asmt/3./PI+(16.11-1.04*(5.-6.63/mtopPole))*
			  sqr(asmt/PI));

  double qsusy = displayMsusy();
  DoubleMatrix mStopSquared(2, 2);
  DoubleMatrix m(2, 2); 
  m(1, 1) = sqr(displayPhys().mu.display(1, 3));
  m(2, 2) = sqr(displayPhys().mu.display(2, 3));
  mStopSquared = rot2d(-displayPhys().thetat) * m * 
    rot2d(displayPhys().thetat);

  double mtrun;

  /// Loop is for iterative solution
  int i; for (i=1; i<=3; i++) {
    mtrun = ssmqcd(mtmt, qsusy, mtopPole);
    
    mql3 = sqrt(mStopSquared(1, 1) - sqr(mtrun) - 
		(0.5 - 2.0 / 3.0 * sn2thw) * sqr(amz) * cos(2.0 * beta));
    mur3 = sqrt(mStopSquared(2, 2) - sqr(mtrun) -
		2.0 / 3.0 * sn2thw * sqr(amz) * cos(2.0 * beta));
    qsusy = sqrt(mql3 * mur3);
  }
  
  DoubleMatrix mSbotSquared(2, 2);
  m(1, 1) = sqr(displayPhys().md.display(1, 3));
  m(2, 2) = sqr(displayPhys().md.display(2, 3));
  mSbotSquared = rot2d(-displayPhys().thetab) * m * 
    rot2d(displayPhys().thetab);

  double asmb = sualfs(sqr(ambt), 0.36, mtopPole);
  double mbmb = ambt * (1.0 - 4.0 * asmb / (3.0 * PI));
  double mbq = ssmqcd(mbmb, qsusy, mtopPole);

  /// We're fiddling these variables to reproduce inside isajet the squark mass
  /// matrices that we've found in softsusy
  mdr3 = sqrt(mSbotSquared(2, 2) - sqr(mbq) 
	      + 1.0 / 3.0 * sn2thw * sqr(amz) * cos(2.0 * beta));

  DoubleMatrix mStauSquared(2, 2);
  m(1, 1) = sqr(displayPhys().me.display(1, 3));
  m(2, 2) = sqr(displayPhys().me.display(2, 3));
  mStauSquared = rot2d(-displayPhys().thetatau) * m * 
    rot2d(displayPhys().thetatau);

  double mlq = 1.7463;

  mtauL = sqrt(mStauSquared(1, 1) - sqr(mlq) + (sqr(amw) -  sqr(amz) * 0.5) *
	       cos(2.0 * beta));
  mtauR = sqrt(mStauSquared(2, 2) - sqr(mlq) - (sqr(amw) - sqr(amz)) * 
	       cos(2.0 * beta));

  at  = mStopSquared(1, 2) / mtrun + smu / tan(beta);
  ab  = mSbotSquared(1, 2) / mbq + smu * tan(beta);  
  atau = mStauSquared(1, 2) / mlq + smu * tan(beta);

  mGPole   = store.displayPhys().mGluino;
  mA       = store.displayPhys().mA0(1);
  
  tanb = store.displayTanb();
}

/// First name input is the name of an OUTPUT file from ssrun, the second name
/// is the name of the interface file for INPUT to ssrun
template<class SoftPars>
void Softsusy<SoftPars>::ssrunInterface764(const char fname [80], 
					   const char softfname [80]) 
  const {
  std::fstream softOutput(softfname, std::ios::out);

  ssrunInterface764Inside(fname, softOutput);

  softOutput.close();
}

template<class SoftPars>
void Softsusy<SoftPars>::ssrunInterface764Inside(const char fname [80], 
					   std::fstream & softOutput) 
  const { 

  double mtopPole, mGPole, smu, mA, tanb, mq1l, mdr, mur, 
    meL, meR, mql3, mdr3, mur3,  mtauL, mtauR, at, ab, atau, 
    mq2l, msr, mcr, mmuL, mmuR, m1, m2;

  isajetNumbers764(mtopPole, mGPole, smu, mA, tanb, mq1l, mdr, mur, 
		   meL, meR, mql3, mdr3, mur3,  mtauL, mtauR, at, ab, atau, 
		   mq2l, msr, mcr, mmuL, mmuR, m1, m2);

  softOutput << "'" << fname << "'" << endl;

  softOutput << mtopPole << endl;
  
  softOutput << mGPole  << "," << smu << "," << mA << "," << tanb
       << endl;

  softOutput << mq1l    << "," << mdr << "," << mur << "," 
       << meL     << "," << meR << endl;

  softOutput << mql3    << "," << mdr3  << "," << mur3 << "," 
       << mtauL   << "," << mtauR << "," 
       << at      << "," << ab    << "," << atau
       << endl;

  softOutput << mq2l << "," << msr << "," << mcr << ","  
       << mmuL << "," << mmuR
       << endl;

  softOutput << m1 << "," << m2 << endl;

  if (displayGravitino() > 6.66e-66 && displayGravitino() < 1.e18) 
    softOutput << displayGravitino() << endl;
  else softOutput << "/" << endl;
}

template<class SoftPars>
void Softsusy<SoftPars>::isawigInterface764(const char herwigInputFile [80], 
				      const char isajetOutputFile [80],
				      const char softOutputFile [80])  
  const {

  std::fstream softOutput(softOutputFile, std::ios::out);

  ssrunInterface764Inside(isajetOutputFile, softOutput);

  softOutput << "n" << endl << "n" << endl;

  softOutput << "'" << herwigInputFile << "'" << endl;

  softOutput.close();
}

/// For input into isajet parameter file called fname
template<class SoftPars>
void Softsusy<SoftPars>::isajetInterface764(const char fname[80]) const {

  std::fstream softOutput(fname, std::ios::out);

  double mtopPole, mGPole, smu, mA, tanb, mq1l, mdr, mur, 
    meL, meR, mql3, mdr3, mur3,  mtauL, mtauR, at, ab, atau, 
    mq2l, msr, mcr, mmuL, mmuR, m1, m2;

  isajetNumbers764(mtopPole, mGPole, smu, mA, tanb, mq1l, mdr, mur, 
		   meL, meR, mql3, mdr3, mur3,  mtauL, mtauR, at, ab, atau, 
		   mq2l, msr, mcr, mmuL, mmuR, m1, m2);

  /// PRINT out the fiddle in ISAJET PARAMETER file format
  softOutput << "TMASS"  << endl;
  softOutput << mtopPole << "/" << endl;

  softOutput << "MSSMA" << endl;
  softOutput << mGPole  << "," << smu << "," << mA << "," << tanb
       << "/"     << endl;

  softOutput << "MSSMB" << endl;
  softOutput << mq1l    << "," << mdr << "," << mur << "," 
       << meL     << "," << meR << "/" << endl;

  softOutput << "MSSMC" << endl;
  softOutput << mql3    << "," << mdr3  << "," << mur3 << "," 
       << mtauL   << "," << mtauR << "," 
       << at      << "," << ab    << "," << atau
       << "/" << endl;

  softOutput << "MSSMD" << endl;
  softOutput << mq2l << "," << msr << "," << mcr << ","  
       << mmuL << "," << mmuR
       << "/" << endl;

  softOutput << "MSSME" << endl;
  softOutput << m1 << "," << m2 << "/" << endl;

  softOutput.close();
}

template<class SoftPars>
void Softsusy<SoftPars>::headerSLHA(std::ostream & out) {

  out.setf(std::ios::scientific, std::ios::floatfield);
  out.precision(8);

  out << "# SOFTSUSY" << SOFTSUSY_VERSION << " SLHA compliant output" << endl;
  out << "# B.C. Allanach, Comput. Phys. Commun. 143 (2002) 305-331,";
  out << " hep-ph/0104145\n";
}

template<class SoftPars>
void Softsusy<SoftPars>::spinfoSLHA(std::ostream & out) {
  out << "Block SPINFO          # Program information\n"
      << "     1    SOFTSUSY    # spectrum calculator\n";
  out << "     2    " << SOFTSUSY_VERSION << "       # version number\n";
  if (displayProblem().noConvergence)
    out << "     3   Possible problem: Not achieved desired accuracy of "
	<< TOLERANCE << "- got " 
	<< fracDiff << endl;
  if (displayProblem().inaccurateHiggsMass)
    out << "     3   # Warning: Higgs masses are very inaccurate at this point.\n";
  int posj = 0, posi = 0; double mass = 0.;
  int temp = lsp(mass, posi, posj);
  if (temp != 0 && temp != -1) {
    out << "     3   # Warning: " << recogLsp(temp, posj);
    out << " LSP" << endl;
  }
  if (displayProblem().testSeriousProblem()) 
    out << "     4   Point invalid: " << displayProblem() << endl;
}

template<class SoftPars>
void Softsusy<SoftPars>::softsusySLHA(std::ostream & out) {
  out << "# SOFTSUSY-specific non SLHA information:\n";
  out << "# MIXING=" << MIXING << " Desired accuracy=" << TOLERANCE << " Achieved accuracy=" << displayFracDiff() << endl;
}

template<class SoftPars>
void Softsusy<SoftPars>::higgsMSLHA(std::ostream & out) {
  out << "        25    "; printRow(out, displayPhys().mh0(1)); out << "   # h0\n";
  out << "        35    "; printRow(out, displayPhys().mh0(2)); out << "   # H0\n";
  out << "        36    "; printRow(out, displayPhys().mA0(1)); out << "   # A0\n";
  out << "        37    "; printRow(out, displayPhys().mHpm); out << "   # H+\n";
}

template<class SoftPars>
void Softsusy<SoftPars>::neutralinoCharginoMSLHA(std::ostream & out) {
  sPhysical s(displayPhys());

  out << "   1000022    "; printRow(out, s.mneut(1)); 
  out << "   # ~neutralino(1)\n";
  out << "   1000023    "; printRow(out, s.mneut(2)); 
  out << "   # ~neutralino(2)\n";
  out << "   1000024    "; printRow(out, fabs(s.mch(1))); out << "   # ~chargino(1)\n";
  out << "   1000025    "; printRow(out, s.mneut(3));
  out << "   # ~neutralino(3)\n";
  out << "   1000035    "; printRow(out, s.mneut(4));
  out << "   # ~neutralino(4)\n";
  out << "   1000037    "; printRow(out, fabs(s.mch(2))); out << "   # ~chargino(2)\n";
}

template<class SoftPars>
void Softsusy<SoftPars>::massSLHA(std::ostream & out) {
  sPhysical s(displayPhys());

  out << "Block MASS                      # Mass spectrum\n";
  out << "# PDG code     mass             particle\n";
  /// out << "          6   "; printRow(out, displayDataSet().displayPoleMt()); 
  /// out << "   # top\n";
  out << "        24    "; printRow(out, displayMw()); out << "   # MW\n";
  higgsMSLHA(out);
  out << "   1000021    "; printRow(out, s.mGluino); out << "   # ~g\n";
  neutralinoCharginoMSLHA(out);
  const double underflow = 1.0e-120, defaultG = 1.0e18;
  if (fabs(displayGravitino()) > underflow && 
      fabs(displayGravitino()) < defaultG) 
    out << "   1000039     " << displayGravitino() << "   # ~gravitino\n";
  sfermionsSLHA(out);
}

template<class SoftPars>
void Softsusy<SoftPars>::sfermionsSLHA(std::ostream & out) {
  sPhysical s(displayPhys());

  out << "   1000001    "; printRow(out, s.md(1, 1)); out << "   # ~d_L\n";
  out << "   1000002    "; printRow(out, s.mu(1, 1)); out << "   # ~u_L\n";
  out << "   1000003    "; printRow(out, s.md(1, 2)); out << "   # ~s_L\n";
  out << "   1000004    "; printRow(out, s.mu(1, 2)); out << "   # ~c_L\n";
  out << "   1000005    "; printRow(out, minimum(s.md(1, 3), s.md(2, 3)));
  out << "   # ~b_1\n";   out << "   1000006    "; 
  printRow(out, minimum(s.mu(1, 3), s.mu(2, 3)));
  out << "   # ~t_1\n";
  out << "   1000011    "; printRow(out, s.me(1, 1)); out << "   # ~e_L\n";
  out << "   1000012    "; printRow(out, s.msnu(1)); out << "   # ~nue_L\n";
  out << "   1000013    "; printRow(out, s.me(1, 2)); out << "   # ~mu_L\n";
  out << "   1000014    "; printRow(out, s.msnu(2)); out << "   # ~numu_L\n";
  out << "   1000015    "; printRow(out, minimum(s.me(1, 3), s.me(2, 3)));
  out << "   # ~stau_1\n";
  out << "   1000016    "; printRow(out, s.msnu(3)); out << "   # ~nu_tau_L\n";
  out << "   2000001    "; printRow(out, s.md(2, 1)); out << "   # ~d_R\n";
  out << "   2000002    "; printRow(out, s.mu(2, 1)); out << "   # ~u_R\n";
  out << "   2000003    "; printRow(out, s.md(2, 2)); out << "   # ~s_R\n";
  out << "   2000004    "; printRow(out, s.mu(2, 2)); out << "   # ~c_R\n";
  out << "   2000005    "; printRow(out, maximum(s.md(1, 3), s.md(2, 3)));
  out << "   # ~b_2\n";
  out << "   2000006    "; printRow(out, maximum(s.mu(1, 3), s.mu(2, 3)));
  out << "   # ~t_2\n";
  out << "   2000011    "; printRow(out, s.me(2, 1)); out << "   # ~e_R\n";
  out << "   2000013    "; printRow(out, s.me(2, 2)); out << "   # ~mu_R\n";
  out << "   2000015    "; printRow(out, maximum(s.me(1, 3), s.me(2, 3)));
  out << "   # ~stau_2\n";
}

template<class SoftPars>
void Softsusy<SoftPars>::alphaSLHA(std::ostream & out) {
  out << "Block alpha                   " << 
    "  # Effective Higgs mixing parameter\n";
  out << "          "; printRow(out, displayPhys().thetaH);        
  out << "       # alpha\n";
}

template<class SoftPars>
void Softsusy<SoftPars>::neutralinoMixingSLHA(std::ostream & out) {
  const sPhysical s(displayPhys());

  out << "Block nmix                  # neutralino mixing matrix\n";
  const int rank = s.mneut.displayEnd();
  for (int i = 1; i <= rank; i++) {
    for (int j = 1; j <= rank; j++) {
      out << "  " << i << "  " << j << "    ";
      printRow(out, s.mixNeut(j, i));
      out << "   # N_{" << i << "," << j << "}\n";
    }
  }
}

template<class SoftPars>
void Softsusy<SoftPars>::inomixingSLHA(std::ostream & out) {
  sPhysical s(displayPhys());
  int i, j;

  neutralinoMixingSLHA(out);
  
  DoubleMatrix u(rot2d(s.thetaL)), v(rot2d(s.thetaR)); 

  /// implementing positive-only chargino masses in SLHA
  if (s.mch(1) < 0.) {
    v(1, 1) = -v(1, 1); v(1, 2) = -v(1, 2);
  }
  if (s.mch(2) < 0.) {
    v(2, 1) = -v(2, 1); v(2, 2) = -v(2, 2);
  }

  i = 1;
  out << "Block Umix                  # chargino U mixing matrix \n";
  for (i=1; i<=2; i++)
    for (j=1; j<=2; j++) {
      out << "  " << i << "  " << j << "    "; printRow(out, u(i, j));
      out << "   # U_{" << i << "," << j << "}\n";      
    }

  out << "Block Vmix                  # chargino V mixing matrix \n";
  for (i=1; i<=2; i++)
    for (j=1; j<=2; j++) {
      out << "  " << i << "  " << j << "    "; printRow(out, v(i, j));
      out << "   # V_{" << i << "," << j << "}\n";      
    }  

}

template<class SoftPars>
void Softsusy<SoftPars>::modselSLHA(std::ostream & out, const char model[]) {
  out << "Block MODSEL  # Select model\n";
  int modsel = 0;
  if (!strcmp(model, "sugra")) modsel = 1;
  if (!strcmp(model, "gmsb")) modsel = 2;
  if (!strcmp(model, "amsb")) modsel = 3;
  if (!strcmp(model, "splitgmsb")) modsel = 4;
  out << "     1    " << modsel << "   # " << model << "\n"; /// Les Houches
							    /// accord codes
}

template<class SoftPars>
void Softsusy<SoftPars>::sfermionmixSLHA(std::ostream & out) {
    sPhysical s(displayPhys());
    DoubleMatrix m(2, 2);
    out << "Block stopmix               # stop mixing matrix\n";
    if (s.mu(1, 3) < s.mu(2, 3)) m = rot2d(s.thetat);
    else m = rot2dTwist(s.thetat);
    int i, j; 
    for (i=1; i<=2; i++)
      for (j=1; j<=2; j++) {
	out << "  " << i << "  " << j << "    ";  printRow(out, m(i, j));
	out << "   # F_{" << i << j << "}" << endl;
      }

    out << "Block sbotmix               # sbottom mixing matrix\n";
    if (s.md(1, 3) < s.md(2, 3)) m = rot2d(s.thetab);
    else m = rot2dTwist(s.thetab);
    for (i=1; i<=2; i++)
      for (j=1; j<=2; j++) {
	out << "  " << i << "  " << j << "    "; printRow(out, m(i, j));
	out << "   # F_{" << i << j << "}" << endl;
      }
    
    out << "Block staumix               # stau mixing matrix\n";
    if (s.me(1, 3) < s.me(2, 3)) m = rot2d(s.thetatau);
    else m = rot2dTwist(s.thetatau);
    for (i=1; i<=2; i++)
      for (j=1; j<=2; j++) {
	out << "  " << i << "  " << j << "    "; printRow(out, m(i, j));
	out << "   # F_{" << i << j << "}" << endl;
      }
}

template<class SoftPars>
void Softsusy<SoftPars>::gaugeSLHA(std::ostream & out) {
  double gp = displayGaugeCoupling(1) * sqrt(0.6);
  out << "Block gauge Q= " << displayMu() << "  # SM gauge couplings\n";
  out << "     1     " << gp << "   # g'(Q)MSSM DRbar"
      << endl;   
  out << "     2     " << displayGaugeCoupling(2) 
      << "   # g(Q)MSSM DRbar" << endl;   
  out << "     3     " << displayGaugeCoupling(3) << "   # g3(Q)MSSM DRbar" 
      << endl;   
}

template<class SoftPars>
void Softsusy<SoftPars>::yukawasSLHA(std::ostream & out) {
      out << "Block yu Q= " << displayMu() << "  \n"
	  << "  3  3     " << displayYukawaElement(YU, 3, 3) 
	  << "   # Yt(Q)MSSM DRbar" << endl;
      out << "Block yd Q= " << displayMu() << "  \n"
	  << "  3  3     " << displayYukawaElement(YD, 3, 3) 
	  << "   # Yb(Q)MSSM DRbar" << endl;
      out << "Block ye Q= " << displayMu() << "  \n"
	  << "  3  3     " << displayYukawaElement(YE, 3, 3) 
	  << "   # Ytau(Q)MSSM DRbar" << endl;
}

template<class SoftPars>
void Softsusy<SoftPars>::hmixSLHA(std::ostream & out) {
  out << "Block hmix Q= " << displayMu() << 
    " # Higgs mixing parameters\n";
  out << "     1    "; printRow(out, displaySusyMu()); 
  out << "    # mu(Q)MSSM DRbar\n";
  out << "     2    "; printRow(out, displayTanb()); 
  out << "    # tan beta(Q)MSSM DRbar Feynman gauge\n";
  out << "     3    "; printRow(out, displayHvev()); 
  out << "    # higgs vev(Q)MSSM DRbar Feynman gauge\n";
  out << "     4    "; 
  printRow(out, displayM3Squared() / 
	   (sin(atan(displayTanb())) * cos(atan(displayTanb())))); 
  out << "    # mA^2(Q)MSSM DRbar\n";
}

template<class SoftPars>
void Softsusy<SoftPars>::msoftSLHA(std::ostream & out) {
      out << "Block msoft Q= " << displayMu() 
	  << "  # MSSM DRbar SUSY breaking parameters\n"; 
      int i;
      for (i=1; i<=3; i++) {
	out << "     " << i << "    "; 
	printRow(out, displayGaugino(i)); 
	out << "      # M_" << i << "(Q)" << endl;      
      }
      
      out << "    21    "; printRow(out, displayMh1Squared()); 
      out << "      # mH1^2(Q)" << endl;    
      out << "    22    "; printRow(out, displayMh2Squared()); 
      out << "      # mH2^2(Q)" << endl;    
      
      out << "    31    "; printRow(out, ccbSqrt(displaySoftMassSquared(mLl, 1, 1))); 
      out << "      # meL(Q)" << endl;    
      out << "    32    "; printRow(out, ccbSqrt(displaySoftMassSquared(mLl, 2, 2))); 
      out << "      # mmuL(Q)" << endl;    
      out << "    33    "; printRow(out, ccbSqrt(displaySoftMassSquared(mLl, 3, 3))); 
      out << "      # mtauL(Q)" << endl;    
      out << "    34    "; printRow(out, ccbSqrt(displaySoftMassSquared(mEr, 1, 1))); 
      out << "      # meR(Q)" << endl;    
      out << "    35    "; printRow(out, ccbSqrt(displaySoftMassSquared(mEr, 2, 2))); 
      out << "      # mmuR(Q)" << endl;    
      out << "    36    "; printRow(out, ccbSqrt(displaySoftMassSquared(mEr, 3, 3))); 
      out << "      # mtauR(Q)" << endl;    
      out << "    41    "; printRow(out, ccbSqrt(displaySoftMassSquared(mQl, 1, 1))); 
      out << "      # mqL1(Q)" << endl;    
      out << "    42    "; printRow(out, ccbSqrt(displaySoftMassSquared(mQl, 2, 2))); 
      out << "      # mqL2(Q)" << endl;    
      out << "    43    "; printRow(out, ccbSqrt(displaySoftMassSquared(mQl, 3, 3))); 
      out << "      # mqL3(Q)" << endl;    
      out << "    44    "; printRow(out, ccbSqrt(displaySoftMassSquared(mUr, 1, 1))); 
      out << "      # muR(Q)" << endl;    
      out << "    45    "; printRow(out, ccbSqrt(displaySoftMassSquared(mUr, 2, 2))); 
      out << "      # mcR(Q)" << endl;    
      out << "    46    "; printRow(out, ccbSqrt(displaySoftMassSquared(mUr, 3, 3))); 
      out << "      # mtR(Q)" << endl;    
      out << "    47    "; printRow(out, ccbSqrt(displaySoftMassSquared(mDr, 1, 1))); 
      out << "      # mdR(Q)" << endl;    
      out << "    48    "; printRow(out, ccbSqrt(displaySoftMassSquared(mDr, 2, 2))); 
      out << "      # msR(Q)" << endl;    
      out << "    49    "; printRow(out, ccbSqrt(displaySoftMassSquared(mDr, 3, 3))); 
      out << "      # mbR(Q)" << endl;    
      
      out << "Block au Q= " << displayMu() << "  \n" 
	  << "  1  1    "; printRow(out, displaySoftA(UA, 1, 1));
      out  << "      # Au(Q)MSSM DRbar" << endl   
	   << "  2  2    "; printRow(out, displaySoftA(UA, 2, 2));
      out  << "      # Ac(Q)MSSM DRbar" << endl   
	   << "  3  3    "; printRow(out, displaySoftA(UA, 3, 3));
      out  << "      # At(Q)MSSM DRbar" << endl;   
      out << "Block ad Q= " << displayMu() << "  \n"
	 << "  1  1    "; printRow(out, displaySoftA(DA, 1, 1));
      out  << "      # Ad(Q)MSSM DRbar" << endl   
	   << "  2  2    "; printRow(out, displaySoftA(DA, 2, 2));
      out  << "      # As(Q)MSSM DRbar" << endl   
	   << "  3  3    "; printRow(out, displaySoftA(DA, 3, 3));
      out << "      # Ab(Q)MSSM DRbar" << endl;   
      out << "Block ae Q= " << displayMu() << "  \n"
	  << "  1  1    "; printRow(out, displaySoftA(EA, 1, 1));
      out  << "      # Ae(Q)MSSM DRbar" << endl   
	   << "  2  2    "; printRow(out, displaySoftA(EA, 2, 2));
      out  << "      # Amu(Q)MSSM DRbar" << endl   
	   << "  3  3    ";   printRow(out, displaySoftA(EA, 3, 3));
      out << "      # Atau(Q)MSSM DRbar" << endl;   
}

template<class SoftPars>
void Softsusy<SoftPars>::drbarSLHA(std::ostream & out, int numPoints, double qMax, int n) {
  /// Starting non-essential information. The following decides what scale to
  /// output the running parameters at. It depends upon what qMax is and how
  /// many points the user has requested.
  /// For qMax = 0 and 1 point (defaults), Q=MSUSY is printed out.
  /// For qMax = 0 and n points, points are spaced logarithmically between MZ
  /// and MSUSY.
  /// For qMax != 0 and 1 point, Q=qMax is printed.
  /// For qMax != 0 and n points, points are log spaced between MZ and qMax.
  
  double ms = displayMsusy();
  double q = minimum(ms, MZ);
  
  if (numPoints == 1 && qMax < EPSTOL) q = ms;
  else if (numPoints == 1 && qMax > EPSTOL) q = qMax;
  else if (numPoints > 1 && qMax < EPSTOL) qMax = ms;

  if (numPoints > 1) { 
    
    if (n > 1) {
      double logq = (log(qMax) - log(MZ)) * double(n-1) / 
	double(numPoints-1) + log(MZ);
      q = exp(logq);
    }
	else q = MZ;
  }
  
  runto(q);
  gaugeSLHA(out);
  yukawasSLHA(out);
  hmixSLHA(out);
  msoftSLHA(out);
}

template<class SoftPars>
void Softsusy<SoftPars>::sminputsSLHA(std::ostream & out) {
  QedQcd_legacy d(displayDataSet());
  out << "Block SMINPUTS             # Standard Model inputs\n";
  out << "     1   "; printRow(out, 1.0 / d.displayAlpha(legacy::ALPHA)); 
  out << "   # alpha_em^(-1)(MZ) SM MSbar\n";
  out << "     2   "; printRow(out, GMU); out << "   # G_Fermi\n";
  out << "     3   "; printRow(out, d.displayAlpha(legacy::ALPHAS)); 
  out << "   # alpha_s(MZ)MSbar\n";
  out << "     4   "; printRow(out, displayMz()); out << "   # MZ(pole)\n";
  out << "     5   "; printRow(out, d.displayMbMb()); out << "   # mb(mb)\n";
  out << "     6   "; printRow(out, d.displayPoleMt()); 
  out << "   # Mtop(pole)\n";
  out << "     7   "; printRow(out, d.displayPoleMtau()); 
  out << "   # Mtau(pole)\n";
}

template<class SoftPars>
void Softsusy<SoftPars>::extparSLHA(std::ostream & out, 
                                    const DoubleVector & pars,
                                    bool ewsbBCscale) {
  out << "Block EXTPAR               # non-universal SUSY breaking parameters\n";
  if (ewsbBCscale) 
    out << "     0    -1.00000000e+00  # Set MX=MSUSY\n";
  else {
    out << "     0    "; printRow(out, mxBC); out << "  # MX scale\n";
  }
  
  int i;
  for (i=1; i<=3; i++) {
    out << "     " << i << "    "; 
    printRow(out, pars.display(i)); 
    out << "  # M_" << i << "(MX)" << endl;      
  }
  out << "     11   "; printRow(out, pars.display(11)) ; 
  out << "  # At(MX)" << endl;    
  out << "     12   "; printRow(out, pars.display(12)) ; 
  out << "  # Ab(MX)" << endl;    
  out << "     13   "; printRow(out, pars.display(13)) ; 
  out << "  # Atau(MX)" << endl;    
  if (!altEwsb) {
    out << "     21   "; printRow(out, pars.display(21)) ; 
    out << "  # mHd^2(MX)" << endl;    
    out << "     22   "; printRow(out, pars.display(22)) ; 
    out << "  # mHu^2(MX)" << endl;    
  } else {
    out << "     23   "; printRow(out, displayMuCond()) ; 
    out << "  # mu(MX)" << endl;    
    out << "     26   "; printRow(out, displayMaCond()) ; 
    out << "  # mA(pole)" << endl;    
  }
  if (displaySetTbAtMX()) {
    out << "     25   "; printRow(out, pars.display(25)) ; 
    out << "  # tan beta(MX)" << endl;    
  }
  out << "     31   "; printRow(out, pars.display(31)) ; 
  out << "  # meL(MX)" << endl;    
  out << "     32   "; printRow(out, pars.display(32)) ; 
  out << "  # mmuL(MX)" << endl;    
  out << "     33   "; printRow(out, pars.display(33)) ; 
  out << "  # mtauL(MX)" << endl;    
  out << "     34   "; printRow(out, pars.display(34)) ; 
  out << "  # meR(MX)" << endl;    
  out << "     35   "; printRow(out, pars.display(35)) ; 
  out << "  # mmuR(MX)" << endl;    
  out << "     36   "; printRow(out, pars.display(36)) ; 
  out << "  # mtauR(MX)" << endl;    
  out << "     41   "; printRow(out, pars.display(41)) ; 
  out << "  # mqL1(MX)" << endl;    
  out << "     42   "; printRow(out, pars.display(42)) ; 
  out << "  # mqL2(MX)" << endl;    
  out << "     43   "; printRow(out, pars.display(43)) ; 
  out << "  # mqL3(MX)" << endl;    
  out << "     44   "; printRow(out, pars.display(44)) ; 
  out << "  # muR(MX)" << endl;    
  out << "     45   "; printRow(out, pars.display(45)) ; 
  out << "  # mcR(MX)" << endl;    
  out << "     46   "; printRow(out, pars.display(46)) ; 
  out << "  # mtR(MX)" << endl;    
  out << "     47   "; printRow(out, pars.display(47)) ; 
  out << "  # mdR(MX)" << endl;    
  out << "     48   "; printRow(out, pars.display(48)) ; 
  out << "  # msR(MX)" << endl;    
  out << "     49   "; printRow(out, pars.display(49)) ; 
  out << "  # mbR(MX)" << endl;    
}

template<class SoftPars>
void Softsusy<SoftPars>::minparSLHA(std::ostream & out, const char model [], 
			      const DoubleVector & pars, double tanb, 
			      int sgnMu, 
			      bool ewsbBCscale) {
  /// For universal models, users still want to know MX and it has to be
  /// specially printed out as EXTPAR 0
  bool printMX = false;

  out << "Block MINPAR               # SUSY breaking input parameters\n";
  out << "     3   "; printRow(out, tanb)            ; out << "   # tanb, DRbar, Feynman gauge" << endl;
  if (!altEwsb) {
    out << "     4   "; 
    printRow(out, double(sgnMu)); 
    out << "   # sign(mu)"<< endl;
  }
  if (!strcmp(model, "sugra")) {
    out << "     1   "; printRow(out, pars.display(1)); out << "   # m0" << endl;
    out << "     2   "; printRow(out, pars.display(2)) ; out << "   # m12" << endl;
    out << "     5   "; printRow(out, pars.display(3)) ; out << "   # A0" << endl;
    printMX = true;
  }
  else if (!strcmp(model, "gmsb")) {
    out << "     1   "; printRow(out, pars.display(3)); out << "   # lambda" << endl;
    out << "     2   "; printRow(out, pars.display(2)) ; out << "   # M_mess" 
	 << endl;
    out << "     5   "; printRow(out, pars.display(1)) ; out << "   # N5" << endl;
    out << "     6   "; printRow(out, pars.display(4)) ; out << "   # cgrav" << endl;
  }
  else if (!strcmp(model, "splitgmsb")) {
    out << "     1   "; printRow(out, pars.display(2)); out << "   # lambdaL" << endl;
    out << "     2   "; printRow(out, pars.display(3)) ; out << "   # lambdaD" 
	 << endl;
    out << "     5   "; printRow(out, pars.display(1)) ; out << "   # N5" << endl; 
    out << "     6   "; printRow(out, 1.0) ; out << "   # cgrav" << endl;
    out << "     7   "; printRow(out, pars.display(4)) ; out << "   # mMess" << endl;
    out << "     8   "; printRow(out, pars.display(5)) ; out << "   # mu/M2" << endl;
    out << "     9   "; printRow(out, pars.display(6)) ; out << "   # mA(pole)/M2" << endl;
    out << "    10   "; printRow(out, pars.display(7)) ; out << "   # desired mh^0" << endl;
}
  else if (!strcmp(model, "amsb")) {
    out << "     1   "; printRow(out, pars.display(2)) ; out << "   # m0" 
	<< endl;
    out << "     2   "; printRow(out, pars.display(1)); out << "   # m3/2" 
	<< endl;
    printMX = true;
  }
  else 
    if (!strcmp(model, "nonUniversal")) 
      extparSLHA(out, pars, ewsbBCscale);
  else {
    std::ostringstream ii;
    ii << "Attempting to use SUSY Les Houches Accord for model " 
       << model << " - cannot do at present\n";
    throw ii.str();
  }  
  if (printMX) {
  out << "Block EXTPAR               # scale of SUSY breaking BCs\n";
  out << "     0   "; printRow(out, mxBC); out << "   # MX scale\n";
  }
}
 
template<class SoftPars>
void Softsusy<SoftPars>::slha1(std::ostream & out, const char model[], 
			 const DoubleVector & pars, 
			 int sgnMu, double tanb, 
			 double qMax, 
			 int numPoints, 
			 bool ewsbBCscale) {
  lesHouchesAccordOutput(out, model, pars, sgnMu, tanb, qMax, numPoints, 
			 ewsbBCscale);
}

/// SUSY Les Houches accord for interfacing to Monte-Carlos, decay programs etc.
template<class SoftPars>
void Softsusy<SoftPars>::lesHouchesAccordOutput(std::ostream & out, const char model[], 
					  const DoubleVector & pars, 
					  int sgnMu, double tanb, 
					  double qMax, 
					  int numPoints, 
					  bool ewsbBCscale) {
  if (forceSlha1 == true) {
    slha1(out, model, pars, sgnMu, tanb, qMax, numPoints, 
	  ewsbBCscale);
    return;
  }
  int nn = out.precision();
  headerSLHA(out);
  spinfoSLHA(out);
  modselSLHA(out, model);
  sminputsSLHA(out); 
  minparSLHA(out, model, pars, tanb, sgnMu, ewsbBCscale);
  softsusySLHA(out);

  if (!displayProblem().testSeriousProblem() || printRuledOutSpectra) {
    massSLHA(out);
    alphaSLHA(out);
    inomixingSLHA(out);
    sfermionmixSLHA(out);

    int n = 0; while (n < numPoints) {
      n++; drbarSLHA(out, numPoints, qMax, n);
    }
  } else {
    out << "# Declining to write spectrum because of serious problem"
	<< " with point" << endl;
  }  
  out.precision(nn);
}


/// Returns nlsp mass in mass and function return labels which particle is nlsp:
/// 0 is neutralino posi = #, posj = 0
template<class SoftPars>
int Softsusy<SoftPars>::nlsp(double & mass, int & posi, int & posj) const {
  int temp = 0, ntemp, lsppos1, lsppos2, nlsppos1, nlsppos2, pos1, pos2;
  sPhysical s(displayPhys());
  
  double minmass = s.mneut.apply(fabs).min(pos1); lsppos1 = pos1; lsppos2 = 0;
  double nminmass;

  /// up squarks 1
  double lightest = s.mu.apply(fabs).min(pos1, pos2);
  if (lightest < minmass) { 
    nminmass = minmass; nlsppos1 = lsppos1; nlsppos2 = lsppos2; ntemp=0;
    minmass = lightest; lsppos1 = pos1; lsppos2 = pos2; temp = 1; }
  else { nminmass = lightest; nlsppos1 = pos1; nlsppos2 = pos2; ntemp = 1;}

  /// down squarks 2
  lightest = s.md.apply(fabs).min(pos1, pos2);
  if (lightest < minmass) { 
    nminmass = minmass; nlsppos1 = lsppos1; nlsppos2 = lsppos2; ntemp=temp;
    minmass = lightest; lsppos1 = pos1; lsppos2 = pos2; temp = 2;}
  else if (lightest < nminmass){
     nminmass = lightest; nlsppos1 = pos1; nlsppos2 = pos2; ntemp = 2;}

  /// sleptons 3
  lightest = s.me.apply(fabs).min(pos1, pos2);
  if (lightest < minmass) { 
    nminmass = minmass; nlsppos1 = lsppos1; nlsppos2 = lsppos2; ntemp=temp;
    minmass = lightest; lsppos1 = pos1; lsppos2 = pos2; temp = 3;}
  else if (lightest < nminmass){
    nminmass = lightest; nlsppos1 = pos1; nlsppos2 = pos2; ntemp = 3;}

  /// charginos 4
  lightest = s.mch.apply(fabs).min(pos1);
  if (lightest < minmass) { 
    nminmass = minmass; nlsppos1 = lsppos1; nlsppos2 = lsppos2; ntemp=temp;
    minmass = lightest; lsppos1 = pos1; lsppos2 = 0; temp = 4;}
  else if (lightest < nminmass){
    nminmass = lightest; nlsppos1 = pos1; nlsppos2 = 0; ntemp = 4;}

  /// sneutrinos 5
  lightest = s.msnu.apply(fabs).min(pos1);
  if (lightest < minmass) { 
    nminmass = minmass; nlsppos1 = lsppos1; nlsppos2 = lsppos2; ntemp=temp;
    minmass = lightest; lsppos1 = pos1; lsppos2 = 0; temp = 5;}
  else if (lightest < nminmass){
    nminmass = lightest; nlsppos1 = pos1; nlsppos2 = 0; ntemp = 5;}
  
  /// gluino 6
  lightest = s.mGluino;
  if (lightest < minmass) { 
    nminmass = minmass; nlsppos1 = lsppos1; nlsppos2 = lsppos2; ntemp=temp;
    minmass = lightest; lsppos1 = pos1; lsppos2 = 0; temp = 6;}
  else if (lightest < nminmass){
    nminmass = lightest; nlsppos1 = pos1; nlsppos2 = 0; ntemp = 6;}
  
  
  /// check that nlsp is not the next-to-lightest in the group of the lsp
  double nlightest;

  switch (temp){
    /// neutralinos
  case 0: {
    nlightest = s.mneut.apply(fabs).nmin(pos1);
    if (nlightest < nminmass) { 
      nminmass = nlightest; nlsppos1 = pos1; nlsppos2 = 0; ntemp = 0;}
  }
    /// up squarks 1
  case 1: {
    nlightest = s.mu.apply(fabs).nmin(pos1, pos2);
    if (nlightest < nminmass) { 
      nminmass = nlightest; nlsppos1 = pos1; nlsppos2 = pos2; ntemp = 1;}
  }
    /// down squarks 2
  case 2: {
    nlightest = s.md.apply(fabs).nmin(pos1, pos2);
    if (nlightest < nminmass) { 
      nminmass = nlightest; nlsppos1 = pos1; nlsppos2 = pos2; ntemp = 2;}
  }
    /// sleptons 3
  case 3: {
    nlightest = s.me.apply(fabs).nmin(pos1, pos2);
    if (nlightest < nminmass){
      nminmass = nlightest; nlsppos1 = pos1; nlsppos2 = pos2; ntemp = 3;}
  }
    /// charginos 4
  case 4: {
    nlightest = s.mch.apply(fabs).nmin(pos1);
    if (nlightest < nminmass){
      nminmass = nlightest; nlsppos1 = pos1; nlsppos2 = 0; ntemp = 4;}
  }
    /// sneutrinos 5
  case 5: {
    nlightest = s.msnu.apply(fabs).nmin(pos1);
    if (nlightest < nminmass){
      nminmass = nlightest; nlsppos1 = pos1; nlsppos2 = 0; ntemp = 5;}
  }
    /// gluino 6
    /// there is only one gluino mass
  }
  
  mass = nminmass;
  posi = nlsppos1; posj = nlsppos2;

  return ntemp;


}

static double mhTrue = 0.;
const static double sigmaMh = 2.0;

/// Fit to LEP2 Standard Model results
inline double lep2Likelihood(double mh) {
  double minusTwoLnQ = 0.;
  /// the approximation to the LEP2 results in hep-ex/0508037 follows
  if (mh < 114.9) minusTwoLnQ = 718.12 - 6.25 * mh;
  else if (mh < 120.) minusTwoLnQ = 3604.71 - 61.4118 * mh + 
			0.261438 * sqr(mh);
  return exp(-0.5 * minusTwoLnQ);
}

/// smears the likelihood curve for a Standard Model Higgs mass with a 3 GeV
/// Gaussian theoretical error
inline DoubleVector mhIntegrand(double mh, const DoubleVector & /* y */) {
  DoubleVector dydx(1);
  dydx(1) = lep2Likelihood(mh) * 
    exp(-sqr(mhTrue - mh) / (2.0 * sqr(sigmaMh))) ;

  return dydx;
}

/// returns the smeared log likelihood coming from LEP2 Higgs mass bounds
inline double lnLHiggs(double mh) {
  if (mh > 130.) return 0.;
  /// error code
  if (mh < EPSTOL) return -numberOfTheBeast;

  double from = mh - 4.0 * sigmaMh, 
    to = mh + 4.0 * sigmaMh, 
    guess = 1.0, 
    hmin = TOLERANCE * 1.0e-5;
  mhTrue = mh;

  DoubleVector v(1);
  double eps = TOLERANCE * 1.0e-5;
  v(1) = 1.0; 

  /// Runge-Kutta, f(b) = int^b0 I(x) dx, I is integrand => d f / db = I(b)
  /// odeint has a problem at f(0): therefore, define f'(b)=f(b)+1
  integrateOdes(v, from, to, eps, guess, hmin, mhIntegrand, odeStepper); 
  
  if (v(1) < EPSTOL || fabs(v(1) - 1.0) < EPSTOL) return -numberOfTheBeast; 
  else return log((v(1) - 1.0) / (sqrt(2.0 * PI) * sigmaMh));
}


/// checked 22/04/06
template<class SoftPars>
double Softsusy<SoftPars>::smPredictionMW() const {
  double mh = displayPhys().mh0(1);

  double dH = log(mh / 100.);
  double dh = sqr(mh / 100.);
  double mt = displayDataSet().displayPoleMt();
  double dt = sqr(mt / 174.3) - 1.;

  double dZ = MZ / 91.1875 - 1.;
  double deltaAlphaHad = 0.027572;
  /// first contribution is leptonic: from hep-ph/9803313
  /// second from hep-ph/0104304
  double deltaAlpha = 0.0314977 + deltaAlphaHad;
  double dAlpha = deltaAlpha / 0.05907 - 1.;
  double dAlphas = displayDataSet().displayAlpha(legacy::ALPHAS) / 0.119 - 1.;

  double mw0 = 80.38;
  double c1 = 0.05253, c2 = 0.010345, c3 = 0.001021, c4 = -0.000070, 
    c5 = 1.077, c6 = 0.5270, c7 = 0.0698, c8 = 0.004055, c9 = 0.000110, 
    c10 = 0.0716, c11 = 115.0; 

  double ans = mw0 - c1 * dH - c2 * sqr(dH) + c3 * sqr(dH) * sqr(dH) + 
    c4 * (dh - 1.) - c5 * dAlpha + c6 * dt - c7 * sqr(dt) - c8 * dH * dt + 
    c9 * dh * dt - c10 * dAlphas + c11 * dZ;

  return ans;
}

template<class SoftPars>
double Softsusy<SoftPars>::twoLoopGm2(double amu1Loop) const {

  double alpha = displayDataSet().displayAlpha(legacy::ALPHA);
  double mMu   = displayDataSet().displayMass(legacy::mMuon);
  double mSusy = displayMsusy();
  double mu    = displaySusyMu();
  double mtau  = displayDataSet().displayMass(legacy::mTau);
  double cosb  = cos(atan(displayTanb()));
  double sinb  = sin(atan(displayTanb()));
  double sinA  = sin(displayDrBarPars().thetaH);
  double cosA  = cos(displayDrBarPars().thetaH);
  double Atau  = displaySoftA(EA, 3, 3);
  double Atop  = displaySoftA(UA, 3, 3);
  double Abot  = displaySoftA(DA, 3, 3);
  double cosTau = cos(displayDrBarPars().thetatau);
  double sinTau = sin(displayDrBarPars().thetatau);
  double cosTop = cos(displayDrBarPars().thetat);
  double sinTop = sin(displayDrBarPars().thetat);
  double cosBot = cos(displayDrBarPars().thetab);
  double sinBot = sin(displayDrBarPars().thetab);
  double mbot   = displayDataSet().displayMass(legacy::mBottom);
  double mtop   = displayDataSet().displayMass(legacy::mTop);
  double tanb   = displayTanb();
  double mA0    = displayDrBarPars().mA0(1);
  double mh0    = displayDrBarPars().mh0(1);
  double mH0    = displayDrBarPars().mh0(2);
  double sw     = sqrt(1.0 - sqr(MW / MZ));

  DoubleVector mstau(2), msbot(2), mstop(2);
  mstau(1) = displayDrBarPars().me.display(1, 3);   
  mstau(2) = displayDrBarPars().me.display(2, 3);
  mstop(1) = displayDrBarPars().mu.display(1, 3);   
  mstop(2) = displayDrBarPars().mu.display(2, 3);
  msbot(1) = displayDrBarPars().md.display(1, 3);   
  msbot(2) = displayDrBarPars().md.display(2, 3);

  /// Neutralinos
  ComplexMatrix n(forLoops.nBpmz);
  DoubleVector mneut(forLoops.mnBpmz);
  ComplexMatrix u(forLoops.uBpmz), v(forLoops.vBpmz); 
  DoubleVector mch(forLoops.mchBpmz); 

  double logPiece = -amu1Loop * (4.0 * alpha / PI * log(mSusy / mMu));

  /// Higgs-sfermion couplings
  DoubleVector lStauh0(2), lStauH0(2);
  lStauh0(1) = 2.0 * mtau / (sqr(mstau(1)) * cosb) * (-mu * cosA - Atau * sinA)
    * cosTau * sinTau;
  lStauH0(1) = 2.0 * mtau / (sqr(mstau(1)) * cosb) * (-mu * sinA + Atau * cosA)
    * cosTau * sinTau;
  lStauh0(2) = -2.0 * mtau / (sqr(mstau(2)) * cosb) * (-mu * cosA - Atau * sinA)
    * sinTau * cosTau;
  lStauH0(2) = -2.0 * mtau / (sqr(mstau(2)) * cosb) * (-mu * sinA + Atau * cosA)
    * cosTau * sinTau;

  DoubleVector lSboth0(2), lSbotH0(2);
  lSboth0(1) = 2.0 * mbot / (sqr(msbot(1)) * cosb) * (-mu * cosA - Abot * sinA)
    * cosBot * sinBot;
  lSbotH0(1) = 2.0 * mbot / (sqr(msbot(1)) * cosb) * (-mu * sinA + Abot * cosA)
    * cosBot * sinBot;
  lSboth0(2) = -2.0 * mbot / (sqr(msbot(2)) * cosb) * 
    (-mu * cosA - Abot * sinA)
    * sinBot * cosBot;
  lSbotH0(2) = -2.0 * mbot / (sqr(msbot(2)) * cosb) * 
    (-mu * sinA + Abot * cosA)
    * cosBot * sinBot;

  DoubleVector lStoph0(2), lStopH0(2);
  lStoph0(1) = 2.0 * mtop / (sqr(mstop(1)) * sinb) * (mu * sinA + Atop * cosA)
    * cosTop * sinTop;
  lStopH0(1) = 2.0 * mtop / (sqr(mstop(1)) * sinb) * (-mu * cosA + Atop * sinA)
    * cosTop * sinTop;
  lStoph0(2) = -2.0 * mtop / (sqr(mstop(2)) * sinb) * (mu * sinA + Atop * cosA)
    * sinTop * cosTop;
  lStopH0(2) = -2.0 * mtop / (sqr(mstop(2)) * sinb) * 
    (-mu * cosA + Atop * sinA)    * cosTop * sinTop;

  DoubleVector lChiCh0(2), lChiCH0(2), lChiCA0(2);
  int k; for (k=1; k<=2; k++) {
    lChiCh0(k) = root2 * MW / mch(k) *
      ((u(k, 1) * v(k, 2)).real() * cosA - sinA * (u(k, 2) * v(k, 1)).real());
    lChiCH0(k) = root2 * MW / mch(k) *
      ((u(k, 1) * v(k, 2)).real() * sinA + cosA * (u(k, 2) * v(k, 1)).real());
    lChiCA0(k) = root2 * MW / mch(k) *
      (-(u(k, 1) * v(k, 2)).real() * cosb - sinb * (u(k, 2) * v(k, 1)).real());
  }
    
  double lMuh0 = -sinA / cosb;
  double lMuH0 = cosA / cosb;
  double lMuA0 = tanb;

  double achigh = 0.;
  for (k=1; k<=2; k++) {
    achigh = achigh +
      (lMuA0 * lChiCA0(k) * fps(sqr(mch(k) / mA0))) +
      (lMuh0 * lChiCh0(k) * fs(sqr(mch(k) / mh0))) +
      (lMuH0 * lChiCH0(k) * fs(sqr(mch(k) / mH0)));
  }

  double asfgh = 0.;
  for (k=1; k<=2; k++) {
    asfgh = asfgh +
      (4.0 / 3.0 * lMuh0 * lStoph0(k) * ffbar(sqr(mstop(k) / mh0))) +
      (4.0 / 3.0 * lMuH0 * lStopH0(k) * ffbar(sqr(mstop(k) / mH0))) +
      (1.0 / 3.0 * lMuh0 * lSboth0(k) * ffbar(sqr(msbot(k) / mh0))) +
      (1.0 / 3.0 * lMuH0 * lSbotH0(k) * ffbar(sqr(msbot(k) / mH0))) +
      (lMuh0 * lStauh0(k) * ffbar(sqr(mstau(k) / mh0))) +
      (lMuH0 * lStauH0(k) * ffbar(sqr(mstau(k) / mH0)));
  }

  double c = sqr(alpha * mMu / (PI * MW * sw)) / 8.0;
  asfgh  = asfgh * c;
  achigh = achigh * c;

  return logPiece + achigh + asfgh;
}

/// Again, another dummy - useful in alternative EWSB conditions sometimes
template<class SoftPars>
 void Softsusy<SoftPars>::setEwsbConditions(const DoubleVector & inputs) {
   setMuCond(inputs.display(1));
   setMaCond(inputs.display(2));
 }


/// input diagonal matrices and it'll give you back mixed ones
template<class SoftPars>
void Softsusy<SoftPars>::doQuarkMixing(DoubleMatrix & /* mDon */,
                                       DoubleMatrix & /* mUpq */) {
  /// This is a dummy routine - MIXING is ignored in this object (it's all
  /// done in FLAVOURMSSMSOFTSUSY these days).
}


/*
template<class SoftPars>
double MssmSoftsusy::twoLpMt() const {
  const double zt2 = sqr(PI) / 6.;
  double mmsb1 = sqr(displayDrBarPars().md(1, 3));
  double mmsb2 = sqr(displayDrBarPars().md(2, 3));
  double mmst1 = sqr(displayDrBarPars().mu(1, 3));
  double mmst2 = sqr(displayDrBarPars().mu(2, 3));
  double mgl = displayGaugino(3);
  double mmgl = sqr(mgl);
  double mt = displayDrBarPars().mt;
  double mmt = sqr(mt);
  double mb = displayDrBarPars().mb;
  double mmb = sqr(mb);
  double csb = cos(displayDrBarPars().thetab), 
    cs2b = cos(displayDrBarPars().thetab * 2.), 
    cs4b = cos(4 * displayDrBarPars().thetab);
  double snb = sin(displayDrBarPars().thetab), 
    sn2b = sin(displayDrBarPars().thetab * 2.), 
    sn4b = sin(4 * displayDrBarPars().thetab);
  double cst = cos(displayDrBarPars().thetat), 
    cs2t = cos(displayDrBarPars().thetat * 2.), 
    cs4t = cos(4 * displayDrBarPars().thetat);
  double snt = sin(displayDrBarPars().thetat), 
    sn2t = sin(displayDrBarPars().thetat * 2.), 
    sn4t = sin(4 * displayDrBarPars().thetat);
  double mmu = sqr(displayMu());

  /// average of first 2 generations squark mass
  double msq = 0.125 * (displayDrBarPars().mu(1, 1) + 
			displayDrBarPars().mu(2, 1) + 
			displayDrBarPars().md(1, 1) + 
			displayDrBarPars().md(2, 1) + 		       
			displayDrBarPars().mu(1, 2) + 
			displayDrBarPars().mu(2, 2) + 
			displayDrBarPars().md(1, 2) + 
			displayDrBarPars().md(2, 2));
  double mmsusy = sqr(msq);

  double lnMglSq = log(mmgl);
  double lnMsbSq = log(mmsb1);
  double lnMsb2Sq = log(mmsb2);
  double lnMst1Sq = log(mmst1);
  double lnMst2Sq = log(mmst2);
  double lnMmsusy = log(mmsusy);
  double lnMmt = log(mmt);
  double lnMmu = log(mmu);
  
  double resmt =

       + sqr(cs2t) * (
          - 640/9
          - 128/9*zt2
          )

       + fin(mmgl,mmsusy)*den(mmgl - mmst1,1)*sn2t * (
          + 32/3*mmsusy/mt*mgl
          )

       + fin(mmgl,mmsusy)*den(mmgl - mmst1,1) * (
          + 16/3*mmst1
          - 8*mmsusy
          )

       + fin(mmgl,mmsusy)*den(mmgl - mmst1,2)*sn2t * (
          + 32/3*mmst1*mmsusy/mt*mgl
          - 32/3*sqr(mmst1)/mt*mgl
          )

       + fin(mmgl,mmsusy)*den(mmgl - mmst1,2) * (
          - 56/3*mmst1*mmsusy
          + 56/3*sqr(mmst1)
          )

       + fin(mmgl,mmsusy)*den(mmgl - mmst1,3) * (
          - 32/3*sqr(mmst1)*mmsusy
          + 32/3*pow(mmst1,3)
          )

       + fin(mmgl,mmsusy)*den(mmgl - mmst2,1)*sn2t * (
          - 32/3*mmsusy/mt*mgl
          )

       + fin(mmgl,mmsusy)*den(mmgl - mmst2,1) * (
          + 16/3*mmst2
          - 8*mmsusy
          )

       + fin(mmgl,mmsusy)*den(mmgl - mmst2,2)*sn2t * (
          - 32/3*mmst2*mmsusy/mt*mgl
          + 32/3*sqr(mmst2)/mt*mgl
          )

       + fin(mmgl,mmsusy)*den(mmgl - mmst2,2) * (
          - 56/3*mmst2*mmsusy
          + 56/3*sqr(mmst2)
          )

       + fin(mmgl,mmsusy)*den(mmgl - mmst2,3) * (
          - 32/3*sqr(mmst2)*mmsusy
          + 32/3*pow(mmst2,3)
          )

       + fin(mmgl,mmsusy) * (
          - 16/3
          )

       + fin(mmsb1,mmgl)*den(mmgl - mmst1,1)*sn2t * (
          + 4/3*mmsb1/mt*mgl
          )

       + fin(mmsb1,mmgl)*den(mmgl - mmst1,1) * (
          - 1/3*mmsb1
          )

       + fin(mmsb1,mmgl)*den(mmgl - mmst1,2)*sn2t * (
          - 4/3*mmsb1*mmst1/mt*mgl
          + 4/3*sqr(mmsb1)/mt*mgl
          )

       + fin(mmsb1,mmgl)*den(mmgl - mmst1,2) * (
          + mmsb1*sqr(mmst1
          - mmsb1)
          )

       + fin(mmsb1,mmgl)*den(mmgl - mmst1,3) * (
          + 4/3*mmsb1*sqr(mmst1)
          - 4/3*sqr(mmsb1)*mmst1
          )

       + fin(mmsb1,mmgl)*den(mmgl - mmst2,1)*sn2t * (
          - 4/3*mmsb1/mt*mgl
          )

       + fin(mmsb1,mmgl)*den(mmgl - mmst2,1) * (
          - 1/3*mmsb1
          )

       + fin(mmsb1,mmgl)*den(mmgl - mmst2,2)*sn2t * (
          + 4/3*mmsb1*mmst2/mt*mgl
          - 4/3*sqr(mmsb1)/mt*mgl
          )

       + fin(mmsb1,mmgl)*den(mmgl - mmst2,2) * (
          + mmsb1*sqr(mmst2
          - mmsb1)
          )

       + fin(mmsb1,mmgl)*den(mmgl - mmst2,3) * (
          + 4/3*mmsb1*sqr(mmst2)
          - 4/3*sqr(mmsb1)*mmst2
          )

       + fin(mmsb2,mmgl)*den(mmgl - mmst1,1)*sn2t * (
          + 4/3*mmsb2/mt*mgl
          )

       + fin(mmsb2,mmgl)*den(mmgl - mmst1,1) * (
          - 1/3*mmsb2
          )

       + fin(mmsb2,mmgl)*den(mmgl - mmst1,2)*sn2t * (
          - 4/3*mmsb2*mmst1/mt*mgl
          + 4/3*sqr(mmsb2)/mt*mgl
          )

       + fin(mmsb2,mmgl)*den(mmgl - mmst1,2) * (
          + mmsb2*sqr(mmst1
          - mmsb2)
          )

       + fin(mmsb2,mmgl)*den(mmgl - mmst1,3) * (
          + 4/3*mmsb2*sqr(mmst1)
          - 4/3*sqr(mmsb2)*mmst1
          )

       + fin(mmsb2,mmgl)*den(mmgl - mmst2,1)*sn2t * (
          - 4/3*mmsb2/mt*mgl
          )

       + fin(mmsb2,mmgl)*den(mmgl - mmst2,1) * (
          - 1/3*mmsb2
          )

       + fin(mmsb2,mmgl)*den(mmgl - mmst2,2)*sn2t * (
          + 4/3*mmsb2*mmst2/mt*mgl
          - 4/3*sqr(mmsb2)/mt*mgl
          )

       + fin(mmsb2,mmgl)*den(mmgl - mmst2,2) * (
          + mmsb2*sqr(mmst2
          - mmsb2)
          )

       + fin(mmsb2,mmgl)*den(mmgl - mmst2,3) * (
          + 4/3*mmsb2*sqr(mmst2)
          - 4/3*sqr(mmsb2)*mmst2
          )

       + fin(mmst1,mmgl)*den(mmgl - mmst1,1)*sn2t * (
          + 88/9*mmst1/mt*mgl
          )

       + fin(mmst1,mmgl)*den(mmgl - mmst1,1)*sqr(cs2t) * (
          + 5/3*mmst1
          )

       + fin(mmst1,mmgl)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          - 154/9*sqr(mmst1)
          )

       + fin(mmst1,mmgl)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1) * (
          + 34/9*sqr(mmst1)
          )

       + fin(mmst1,mmgl)*den(mmgl - mmst1,1) * (
          + 22/9*mmst1
          )

       + fin(mmst1,mmgl)*den(mmgl - mmst1,2) * (
          + 12*sqr(mmst1)
          )

       + fin(mmst1,mmgl)*den(mmgl - mmst1,3) * (
          + 16/3*pow(mmst1,3)
          )

       + fin(mmst1,mmgl)*den(mmgl - mmst2,1)*sn2t * (
          - 16/9*mmst1/mt*mgl
          )

       + fin(mmst1,mmgl)*den(mmgl - mmst2,1)*sqr(cs2t) * (
          - 5/3*mmst1
          )

       + fin(mmst1,mmgl)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          + 26/9*sqr(mmst1)
          )

       + fin(mmst1,mmgl)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1) * (
          - 34/9*sqr(mmst1)
          )

       + fin(mmst1,mmgl)*den(mmgl - mmst2,1) * (
          + 16/9*mmst1
          )

       + fin(mmst1,mmgl)*den(mmgl - mmst2,2)*sn2t * (
          + 4/3*mmst1*mmst2/mt*mgl
          - 4/3*sqr(mmst1)/mt*mgl
          )

       + fin(mmst1,mmgl)*den(mmgl - mmst2,2)*sqr(cs2t) * (
          + 26/9*mmst1*mmst2
          )

       + fin(mmst1,mmgl)*den(mmgl - mmst2,2) * (
          - 17/9*mmst1*sqr(mmst2
          - mmst1)
          )

       + fin(mmst1,mmgl)*den(mmgl - mmst2,3) * (
          + 4/3*mmst1*sqr(mmst2)
          - 4/3*sqr(mmst1)*mmst2
          )

       + fin(mmst1,mmgl)*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          - 128/9*mmst1
          )

       + fin(mmst1,mmsb1)*den(mmgl - mmst1,1) * (
          - 2/3*mmst1
          )

       + fin(mmst1,mmsb1)*den(mmgl - mmst1,2)*sn2t * (
          - 4/3*mmsb1*mmst1/mt*mgl
          + 4/3*sqr(mmst1)/mt*mgl
          )

       + fin(mmst1,mmsb1)*den(mmgl - mmst1,2) * (
          + mmsb1*mmst1
          - 7/3*sqr(mmst1)
          )

       + fin(mmst1,mmsb1)*den(mmgl - mmst1,3) * (
          + 4/3*mmsb1*sqr(mmst1)
          - 4/3*pow(mmst1,3)
          )

       + fin(mmst1,mmsb2)*den(mmgl - mmst1,1) * (
          - 2/3*mmst1
          )

       + fin(mmst1,mmsb2)*den(mmgl - mmst1,2)*sn2t * (
          - 4/3*mmsb2*mmst1/mt*mgl
          + 4/3*sqr(mmst1)/mt*mgl
          )

       + fin(mmst1,mmsb2)*den(mmgl - mmst1,2) * (
          + mmsb2*mmst1
          - 7/3*sqr(mmst1)
          )

       + fin(mmst1,mmsb2)*den(mmgl - mmst1,3) * (
          + 4/3*mmsb2*sqr(mmst1)
          - 4/3*pow(mmst1,3)
          )

       + fin(mmst1,mmst2)*den(mmgl - mmst1,1)*sn2t * (
          - 4/9*mmst1/mt*mgl
          )

       + fin(mmst1,mmst2)*den(mmgl - mmst1,1)*sqr(cs2t) * (
          - 11/9*mmst1
          )

       + fin(mmst1,mmst2)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1) * (
          - 8/9*sqr(mmst1)
          )

       + fin(mmst1,mmst2)*den(mmgl - mmst1,1) * (
          + mmst1
          )

       + fin(mmst1,mmst2)*den(mmgl - mmst1,2)*sn2t * (
          - 4/3*mmst1*mmst2/mt*mgl
          + 4/3*sqr(mmst1)/mt*mgl
          )

       + fin(mmst1,mmst2)*den(mmgl - mmst1,2)*sqr(cs2t) * (
          - 26/9*sqr(mmst1)
          )

       + fin(mmst1,mmst2)*den(mmgl - mmst1,2) * (
          + mmst1*mmst2
          + 5/9*sqr(mmst1)
          )

       + fin(mmst1,mmst2)*den(mmgl - mmst1,3) * (
          + 4/3*sqr(mmst1)*mmst2
          - 4/3*pow(mmst1,3)
          )

       + fin(mmst1,mmst2)*den(mmgl - mmst2,1)*sn2t * (
          + 4/9*mmst1/mt*mgl
          )

       + fin(mmst1,mmst2)*den(mmgl - mmst2,1)*sqr(cs2t) * (
          - 11/9*mmst1
          )

       + fin(mmst1,mmst2)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1) * (
          + 8/9*sqr(mmst1)
          )

       + fin(mmst1,mmst2)*den(mmgl - mmst2,1) * (
          + 1/9*mmst1
          )

       + fin(mmst1,mmst2)*den(mmgl - mmst2,2)*sn2t * (
          - 4/3*mmst1*mmst2/mt*mgl
          + 4/3*sqr(mmst1)/mt*mgl
          )

       + fin(mmst1,mmst2)*den(mmgl - mmst2,2)*sqr(cs2t) * (
          - 26/9*mmst1*mmst2
          )

       + fin(mmst1,mmst2)*den(mmgl - mmst2,2) * (
          + 5/9*mmst1*sqr(mmst2
          + mmst1)
          )

       + fin(mmst1,mmst2)*den(mmgl - mmst2,3) * (
          - 4/3*mmst1*sqr(mmst2)
          + 4/3*sqr(mmst1)*mmst2
          )

       + fin(mmst1,mmsusy)*den(mmgl - mmst1,1) * (
          - 16/3*mmst1
          )

       + fin(mmst1,mmsusy)*den(mmgl - mmst1,2)*sn2t * (
          - 32/3*mmst1*mmsusy/mt*mgl
          + 32/3*sqr(mmst1)/mt*mgl
          )

       + fin(mmst1,mmsusy)*den(mmgl - mmst1,2) * (
          + 8*mmst1*mmsusy
          - 56/3*sqr(mmst1)
          )

       + fin(mmst1,mmsusy)*den(mmgl - mmst1,3) * (
          + 32/3*sqr(mmst1)*mmsusy
          - 32/3*pow(mmst1,3)
          )

       + fin(mmst2,mmgl)*sqr(cs2t) * (
          - 128/9
          )

       + fin(mmst2,mmgl)*den(mmgl - mmst1,1)*sn2t * (
          + 16/9*mmst2/mt*mgl
          )

       + fin(mmst2,mmgl)*den(mmgl - mmst1,1)*sqr(cs2t) * (
          + 26/9*mmst1
          + 11/9*mmst2
          )

       + fin(mmst2,mmgl)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          - 26/9*sqr(mmst1)
          )

       + fin(mmst2,mmgl)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1) * (
          + 34/9*sqr(mmst1)
          )

       + fin(mmst2,mmgl)*den(mmgl - mmst1,1) * (
          - 34/9*mmst1
          - 2*mmst2
          )

       + fin(mmst2,mmgl)*den(mmgl - mmst1,2)*sn2t * (
          - 4/3*mmst1*mmst2/mt*mgl
          + 4/3*sqr(mmst2)/mt*mgl
          )

       + fin(mmst2,mmgl)*den(mmgl - mmst1,2)*sqr(cs2t) * (
          + 26/9*mmst1*mmst2
          )

       + fin(mmst2,mmgl)*den(mmgl - mmst1,2) * (
          - 17/9*mmst1*sqr(mmst2
          - mmst2)
          )

       + fin(mmst2,mmgl)*den(mmgl - mmst1,3) * (
          - 4/3*mmst1*sqr(mmst2)
          + 4/3*sqr(mmst1)*mmst2
          )

       + fin(mmst2,mmgl)*den(mmgl - mmst2,1)*sn2t * (
          - 88/9*mmst2/mt*mgl
          )

       + fin(mmst2,mmgl)*den(mmgl - mmst2,1)*sqr(cs2t) * (
          - 154/9*mmst1
          - 139/9*mmst2
          )

       + fin(mmst2,mmgl)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          + 154/9*sqr(mmst1)
          )

       + fin(mmst2,mmgl)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1) * (
          - 34/9*sqr(mmst1)
          )

       + fin(mmst2,mmgl)*den(mmgl - mmst2,1) * (
          + 34/9*mmst1
          + 56/9*mmst2
          )

       + fin(mmst2,mmgl)*den(mmgl - mmst2,2) * (
          + 12*sqr(mmst2)
          )

       + fin(mmst2,mmgl)*den(mmgl - mmst2,3) * (
          + 16/3*pow(mmst2,3)
          )

       + fin(mmst2,mmgl)*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          + 128/9*mmst1
          )

       + fin(mmst2,mmsb1)*den(mmgl - mmst2,1) * (
          - 2/3*mmst2
          )

       + fin(mmst2,mmsb1)*den(mmgl - mmst2,2)*sn2t * (
          + 4/3*mmsb1*mmst2/mt*mgl
          - 4/3*sqr(mmst2)/mt*mgl
          )

       + fin(mmst2,mmsb1)*den(mmgl - mmst2,2) * (
          + mmsb1*mmst2
          - 7/3*sqr(mmst2)
          )

       + fin(mmst2,mmsb1)*den(mmgl - mmst2,3) * (
          + 4/3*mmsb1*sqr(mmst2)
          - 4/3*pow(mmst2,3)
          )

       + fin(mmst2,mmsb2)*den(mmgl - mmst2,1) * (
          - 2/3*mmst2
          )

       + fin(mmst2,mmsb2)*den(mmgl - mmst2,2)*sn2t * (
          + 4/3*mmsb2*mmst2/mt*mgl
          - 4/3*sqr(mmst2)/mt*mgl
          )

       + fin(mmst2,mmsb2)*den(mmgl - mmst2,2) * (
          + mmsb2*mmst2
          - 7/3*sqr(mmst2)
          )

       + fin(mmst2,mmsb2)*den(mmgl - mmst2,3) * (
          + 4/3*mmsb2*sqr(mmst2)
          - 4/3*pow(mmst2,3)
          )

       + fin(mmst2,mmsusy)*den(mmgl - mmst2,1) * (
          - 16/3*mmst2
          )

       + fin(mmst2,mmsusy)*den(mmgl - mmst2,2)*sn2t * (
          + 32/3*mmst2*mmsusy/mt*mgl
          - 32/3*sqr(mmst2)/mt*mgl
          )

       + fin(mmst2,mmsusy)*den(mmgl - mmst2,2) * (
          + 8*mmst2*mmsusy
          - 56/3*sqr(mmst2)
          )

       + fin(mmst2,mmsusy)*den(mmgl - mmst2,3) * (
          + 32/3*sqr(mmst2)*mmsusy
          - 32/3*pow(mmst2,3)
          )

       + lnMglSq*sqr(cs2t) * (
          + 128/3
          )

       + sqr(lnMglSq)*sqr(cs2t) * (
          - 64/9
          )

       + sqr(lnMglSq)*den(mmgl - mmst1,1)*sn2t * (
          + 2/3*mmsb1/mt*mgl
          + 2/3*mmsb2/mt*mgl
          + 10*mmst1/mt*mgl
          + 2/3*mmst2/mt*mgl
          + 16/3*mmsusy/mt*mgl
          )

       + sqr(lnMglSq)*den(mmgl - mmst1,1)*sqr(cs2t) * (
          + 53/3*mmst1
          )

       + sqr(lnMglSq)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sn2t * (
          + 8/3*sqr(mmst1)/mt*mgl
          )

       + sqr(lnMglSq)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          - 223/9*sqr(mmst1)
          )

       + sqr(lnMglSq)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1) * (
          + 43/9*sqr(mmst1)
          )

       + sqr(lnMglSq)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,2)*sn2t * (
          - 16/9*pow(mmst1,3)/mt*mgl
          )

       + sqr(lnMglSq)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,2) * (
          + 32/9*pow(mmst1,3)
          )

       + sqr(lnMglSq)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,3) * (
          - 16/9*pow(mmst1,4)
          )

       + sqr(lnMglSq)*den(mmgl - mmst1,1) * (
          - 1/2*mmsb1
          - 1/2*mmsb2
          - 73/6*mmst1
          - 1/2*mmst2
          + 4/3*mmsusy
          )

       + sqr(lnMglSq)*den(mmgl - mmst1,2)*sn2t * (
          + 2/3*mmsb1*mmst1/mt*mgl
          + 2/3*mmsb2*mmst1/mt*mgl
          + 2/3*mmst1*mmst2/mt*mgl
          - 16*mmst1*mmsusy/mt*mgl
          - 122/9*sqr(mmst1)/mt*mgl
          + 32/3*sqr(mmsusy)/mt*mgl
          )

       + sqr(lnMglSq)*den(mmgl - mmst1,2)*sqr(cs2t) * (
          + 53/6*sqr(mmst1)
          )

       + sqr(lnMglSq)*den(mmgl - mmst1,2)*den(mmst1 - mmst2,1)*sn2t * (
          + 8/9*pow(mmst1,3)/mt*mgl
          )

       + sqr(lnMglSq)*den(mmgl - mmst1,2)*den(mmst1 - mmst2,2) * (
          + 8/9*pow(mmst1,4)
          )

       + sqr(lnMglSq)*den(mmgl - mmst1,2) * (
          - 7/6*mmsb1*mmst1
          - 7/6*mmsb2*mmst1
          - 7/6*mmst1*mmst2
          + 52/3*mmst1*mmsusy
          + 151/9*sqr(mmst1)
          - 8*sqr(mmsusy)
          )

       + sqr(lnMglSq)*den(mmgl - mmst1,3)*sn2t * (
          - 8/9*pow(mmst1,3)/mt*mgl
          )

       + sqr(lnMglSq)*den(mmgl - mmst1,3) * (
          - 2/3*mmsb1*sqr(mmst1)
          - 2/3*mmsb2*sqr(mmst1)
          - 32/3*mmst1*sqr(mmsusy)
          - 2/3*sqr(mmst1)*mmst2
          + 16*sqr(mmst1)*mmsusy
          + 46/3*pow(mmst1,3)
          )

       + sqr(lnMglSq)*den(mmgl - mmst1,4) * (
          + 4/9*pow(mmst1,4)
          )

       + sqr(lnMglSq)*den(mmgl - mmst2,1)*sn2t * (
          - 2/3*mmsb1/mt*mgl
          - 2/3*mmsb2/mt*mgl
          + 2/9*mmst1/mt*mgl
          - 98/9*mmst2/mt*mgl
          - 16/3*mmsusy/mt*mgl
          )

       + sqr(lnMglSq)*den(mmgl - mmst2,1)*sqr(cs2t) * (
          - 223/9*mmst1
          - 64/9*mmst2
          )

       + sqr(lnMglSq)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sn2t * (
          - 8/3*sqr(mmst1)/mt*mgl
          )

       + sqr(lnMglSq)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          + 223/9*sqr(mmst1)
          )

       + sqr(lnMglSq)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1) * (
          - 43/9*sqr(mmst1)
          )

       + sqr(lnMglSq)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,2)*sn2t * (
          + 16/9*pow(mmst1,3)/mt*mgl
          )

       + sqr(lnMglSq)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,2) * (
          - 32/9*pow(mmst1,3)
          )

       + sqr(lnMglSq)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,3) * (
          + 16/9*pow(mmst1,4)
          )

       + sqr(lnMglSq)*den(mmgl - mmst2,1) * (
          - 1/2*mmsb1
          - 1/2*mmsb2
          + 109/18*mmst1
          - 101/18*mmst2
          + 4/3*mmsusy
          )

       + sqr(lnMglSq)*den(mmgl - mmst2,2)*sn2t * (
          - 2/3*mmsb1*mmst2/mt*mgl
          - 2/3*mmsb2*mmst2/mt*mgl
          - 14/9*mmst1*mmst2/mt*mgl
          - 8/9*sqr(mmst1)/mt*mgl
          + 16*mmst2*mmsusy/mt*mgl
          + 38/3*sqr(mmst2)/mt*mgl
          - 32/3*sqr(mmsusy)/mt*mgl
          )

       + sqr(lnMglSq)*den(mmgl - mmst2,2)*sqr(cs2t) * (
          + 53/6*sqr(mmst2)
          )

       + sqr(lnMglSq)*den(mmgl - mmst2,2)*den(mmst1 - mmst2,1)*sn2t * (
          + 8/9*pow(mmst1,3)/mt*mgl
          )

       + sqr(lnMglSq)*den(mmgl - mmst2,2)*den(mmst1 - mmst2,1) * (
          - 32/9*pow(mmst1,3)
          )

       + sqr(lnMglSq)*den(mmgl - mmst2,2)*den(mmst1 - mmst2,2) * (
          + 8/9*pow(mmst1,4)
          )

       + sqr(lnMglSq)*den(mmgl - mmst2,2) * (
          - 7/6*mmsb1*mmst2
          - 7/6*mmsb2*mmst2
          + 11/18*mmst1*mmst2
          + 8/3*sqr(mmst1)
          + 52/3*mmst2*mmsusy
          + 53/3*sqr(mmst2)
          - 8*sqr(mmsusy)
          )

       + sqr(lnMglSq)*den(mmgl - mmst2,3)*sn2t * (
          + 8/9*pow(mmst2,3)/mt*mgl
          )

       + sqr(lnMglSq)*den(mmgl - mmst2,3) * (
          - 2/3*mmsb1*sqr(mmst2)
          - 2/3*mmsb2*sqr(mmst2)
          - 2/3*mmst1*sqr(mmst2)
          - 32/3*mmst2*sqr(mmsusy)
          + 16*sqr(mmst2)*mmsusy
          + 46/3*pow(mmst2,3)
          )

       + sqr(lnMglSq)*den(mmgl - mmst2,4) * (
          + 4/9*pow(mmst2,4)
          )

       + sqr(lnMglSq) * (
          - 166/9
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmst1,1)*sn2t * (
          - 4/3*mmsb1/mt*mgl
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmst1,1) * (
          + mmsb1
          - 4/3*mmst1
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmst1,2)*sn2t * (
          - 4/3*mmsb1*mmst1/mt*mgl
          + 8/3*sqr(mmst1)/mt*mgl
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmst1,2) * (
          + 7/3*mmsb1*mmst1
          - 14/3*sqr(mmst1)
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmst1,3) * (
          + 4/3*mmsb1*sqr(mmst1)
          - 8/3*pow(mmst1,3)
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmst2,1)*sn2t * (
          + 4/3*mmsb1/mt*mgl
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmst2,1) * (
          + mmsb1
          - 4/3*mmst2
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmst2,2)*sn2t * (
          + 4/3*mmsb1*mmst2/mt*mgl
          - 8/3*sqr(mmst2)/mt*mgl
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmst2,2) * (
          + 7/3*mmsb1*mmst2
          - 14/3*sqr(mmst2)
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmst2,3) * (
          + 4/3*mmsb1*sqr(mmst2)
          - 8/3*pow(mmst2,3)
          )

       + lnMglSq*lnMsbSq * (
          + 4/3
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmst1,1)*sn2t * (
          - 4/3*mmsb2/mt*mgl
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmst1,1) * (
          + mmsb2
          - 4/3*mmst1
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmst1,2)*sn2t * (
          - 4/3*mmsb2*mmst1/mt*mgl
          + 8/3*sqr(mmst1)/mt*mgl
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmst1,2) * (
          + 7/3*mmsb2*mmst1
          - 14/3*sqr(mmst1)
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmst1,3) * (
          + 4/3*mmsb2*sqr(mmst1)
          - 8/3*pow(mmst1,3)
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmst2,1)*sn2t * (
          + 4/3*mmsb2/mt*mgl
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmst2,1) * (
          + mmsb2
          - 4/3*mmst2
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmst2,2)*sn2t * (
          + 4/3*mmsb2*mmst2/mt*mgl
          - 8/3*sqr(mmst2)/mt*mgl
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmst2,2) * (
          + 7/3*mmsb2*mmst2
          - 14/3*sqr(mmst2)
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmst2,3) * (
          + 4/3*mmsb2*sqr(mmst2)
          - 8/3*pow(mmst2,3)
          )

       + lnMglSq*lnMsb2Sq * (
          + 4/3
          )

       + lnMglSq*lnMst1Sq*sn2t * (
          - 208/9/mt*mgl
          )

       + lnMglSq*lnMst1Sq*sqr(cs2t) * (
          - 64/9
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst1,1)*sn2t * (
          - 28/3*mmst1/mt*mgl
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst1,1)*sn4t*cs2t * (
          - 8/9*mmst1/mt*mgl
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst1,1)*sqr(cs2t) * (
          - 121/9*mmst1
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sn2t
       * (
          - 8/3*sqr(mmst1)/mt*mgl
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sn4t*
      cs2t * (
          + 16/9*sqr(mmst1)/mt*mgl
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sqr(cs2t)
       * (
          + 287/9*sqr(mmst1)
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1) * (
          - 43/9*sqr(mmst1)
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,2)*sn2t
       * (
          + 16/9*pow(mmst1,3)/mt*mgl
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,2) * (
          - 32/9*pow(mmst1,3)
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,3) * (
          + 16/9*pow(mmst1,4)
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst1,1) * (
          - 34/3*mmst1
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst1,2)*sn2t * (
          + 172/9*sqr(mmst1)/mt*mgl
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst1,2)*sn4t*cs2t * (
          - 8/9*sqr(mmst1)/mt*mgl
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst1,2)*sqr(cs2t) * (
          - 11/9*sqr(mmst1)
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst1,2)*den(mmst1 - mmst2,1)*sn2t
       * (
          - 8/9*pow(mmst1,3)/mt*mgl
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst1,2)*den(mmst1 - mmst2,2) * (
          - 8/9*pow(mmst1,4)
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst1,2) * (
          - 130/3*sqr(mmst1)
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst1,3)*sn2t * (
          + 16/9*pow(mmst1,3)/mt*mgl
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst1,3)*sqr(cs2t) * (
          + 16/9*pow(mmst1,3)
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst1,3) * (
          - 212/9*pow(mmst1,3)
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst1,4) * (
          - 8/9*pow(mmst1,4)
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst2,1)*sn2t * (
          + 20/9*mmst1/mt*mgl
          + 8/9*mmst2/mt*mgl
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst2,1)*sn4t*cs2t * (
          + 8/9*mmst1/mt*mgl
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst2,1)*sqr(cs2t) * (
          + 5/3*mmst1
          - 22/9*mmst2
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sn2t
       * (
          + 8/3*sqr(mmst1)/mt*mgl
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sn4t*
      cs2t * (
          - 16/9*sqr(mmst1)/mt*mgl
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sqr(cs2t)
       * (
          - 31/9*sqr(mmst1)
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1) * (
          + 43/9*sqr(mmst1)
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,2)*sn2t
       * (
          - 16/9*pow(mmst1,3)/mt*mgl
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,2) * (
          + 32/9*pow(mmst1,3)
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,3) * (
          - 16/9*pow(mmst1,4)
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst2,1) * (
          - 34/9*mmst1
          + 2/9*mmst2
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst2,2)*sn2t * (
          + 4*mmst1*mmst2/mt*mgl
          + 8/9*sqr(mmst1)/mt*mgl
          - 8/3*sqr(mmst2)/mt*mgl
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst2,2)*sn4t*cs2t * (
          - 8/9*mmst1*mmst2/mt*mgl
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst2,2)*sqr(cs2t) * (
          - 32/9*mmst1*mmst2
          - 52/9*sqr(mmst2)
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst2,2)*den(mmst1 - mmst2,1)*sn2t
       * (
          - 8/9*pow(mmst1,3)/mt*mgl
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst2,2)*den(mmst1 - mmst2,1) * (
          + 32/9*pow(mmst1,3)
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst2,2)*den(mmst1 - mmst2,2) * (
          - 8/9*pow(mmst1,4)
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst2,2) * (
          + 37/9*mmst1*mmst2
          - 8/3*sqr(mmst1)
          + 10/9*sqr(mmst2)
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst2,3)*sqr(cs2t) * (
          - 16/9*mmst1*sqr(mmst2)
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmst2,3) * (
          + 28/9*mmst1*sqr(mmst2)
          - 8/3*pow(mmst2,3)
          )

       + lnMglSq*lnMst1Sq*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          + 256/9*mmgl
          + 256/9*mmst1
          )

       + lnMglSq*lnMst1Sq * (
          + 52/9
          )

       + lnMglSq*lnMst2Sq*sn2t * (
          + 208/9/mt*mgl
          )

       + lnMglSq*lnMst2Sq*sqr(cs2t) * (
          + 64/3
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst1,1)*sn2t * (
          - 28/9*mmst2/mt*mgl
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst1,1)*sn4t*cs2t * (
          + 16/9*mmst1/mt*mgl
          + 8/9*mmst2/mt*mgl
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst1,1)*sqr(cs2t) * (
          - 53/9*mmst1
          - 16/9*mmst2
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sn2t
       * (
          - 8/3*sqr(mmst1)/mt*mgl
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sn4t*
      cs2t * (
          - 16/9*sqr(mmst1)/mt*mgl
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sqr(cs2t)
       * (
          + 31/9*sqr(mmst1)
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1) * (
          - 43/9*sqr(mmst1)
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,2)*sn2t
       * (
          + 16/9*pow(mmst1,3)/mt*mgl
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,2) * (
          - 32/9*pow(mmst1,3)
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,3) * (
          + 16/9*pow(mmst1,4)
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst1,1) * (
          + 61/9*mmst1
          + 25/9*mmst2
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst1,2)*sn2t * (
          - 28/9*mmst1*mmst2/mt*mgl
          + 32/9*sqr(mmst1)/mt*mgl
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst1,2)*sn4t*cs2t * (
          + 8/9*mmst1*mmst2/mt*mgl
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst1,2)*sqr(cs2t) * (
          - 32/9*mmst1*mmst2
          - 52/9*sqr(mmst1)
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst1,2)*den(mmst1 - mmst2,1)*sn2t
       * (
          - 8/9*pow(mmst1,3)/mt*mgl
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst1,2)*den(mmst1 - mmst2,2) * (
          - 8/9*pow(mmst1,4)
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst1,2) * (
          + 53/9*mmst1*mmst2
          + 2*sqr(mmst1)
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst1,3)*sqr(cs2t) * (
          - 16/9*sqr(mmst1)*mmst2
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst1,3) * (
          + 28/9*sqr(mmst1)*mmst2
          - 8/3*pow(mmst1,3)
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst2,1)*sn2t * (
          - 8/9*mmst1/mt*mgl
          + 92/9*mmst2/mt*mgl
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst2,1)*sn4t*cs2t * (
          - 16/9*mmst1/mt*mgl
          - 8/9*mmst2/mt*mgl
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst2,1)*sqr(cs2t) * (
          + 287/9*mmst1
          + 166/9*mmst2
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sn2t
       * (
          + 8/3*sqr(mmst1)/mt*mgl
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sn4t*
      cs2t * (
          + 16/9*sqr(mmst1)/mt*mgl
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sqr(cs2t)
       * (
          - 287/9*sqr(mmst1)
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1) * (
          + 43/9*sqr(mmst1)
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,2)*sn2t
       * (
          - 16/9*pow(mmst1,3)/mt*mgl
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,2) * (
          + 32/9*pow(mmst1,3)
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,3) * (
          - 16/9*pow(mmst1,4)
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst2,1) * (
          - 59/9*mmst1
          - 161/9*mmst2
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst2,2)*sn2t * (
          + 8/9*mmst1*mmst2/mt*mgl
          + 8/9*sqr(mmst1)/mt*mgl
          - 164/9*sqr(mmst2)/mt*mgl
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst2,2)*sn4t*cs2t * (
          + 8/9*sqr(mmst2)/mt*mgl
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst2,2)*sqr(cs2t) * (
          - 11/9*sqr(mmst2)
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst2,2)*den(mmst1 - mmst2,1)*sn2t
       * (
          - 8/9*pow(mmst1,3)/mt*mgl
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst2,2)*den(mmst1 - mmst2,1) * (
          + 32/9*pow(mmst1,3)
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst2,2)*den(mmst1 - mmst2,2) * (
          - 8/9*pow(mmst1,4)
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst2,2) * (
          - 16/9*mmst1*mmst2
          - 8/3*sqr(mmst1)
          - 398/9*sqr(mmst2)
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst2,3)*sn2t * (
          - 16/9*pow(mmst2,3)/mt*mgl
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst2,3)*sqr(cs2t) * (
          + 16/9*pow(mmst2,3)
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst2,3) * (
          - 212/9*pow(mmst2,3)
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmst2,4) * (
          - 8/9*pow(mmst2,4)
          )

       + lnMglSq*lnMst2Sq*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          - 256/9*mmgl
          - 256/9*mmst1
          )

       + lnMglSq*lnMst2Sq * (
          + 52/9
          )

       + lnMglSq*lnMmsusy*den(mmgl - mmst1,1)*sn2t * (
          - 32/3*mmsusy/mt*mgl
          )

       + lnMglSq*lnMmsusy*den(mmgl - mmst1,1) * (
          - 8/3*mmsusy
          )

       + lnMglSq*lnMmsusy*den(mmgl - mmst1,2)*sn2t * (
          + 32*mmst1*mmsusy/mt*mgl
          - 64/3*sqr(mmsusy)/mt*mgl
          )

       + lnMglSq*lnMmsusy*den(mmgl - mmst1,2) * (
          - 104/3*mmst1*mmsusy
          + 16*sqr(mmsusy)
          )

       + lnMglSq*lnMmsusy*den(mmgl - mmst1,3) * (
          + 64/3*mmst1*sqr(mmsusy)
          - 32*sqr(mmst1)*mmsusy
          )

       + lnMglSq*lnMmsusy*den(mmgl - mmst2,1)*sn2t * (
          + 32/3*mmsusy/mt*mgl
          )

       + lnMglSq*lnMmsusy*den(mmgl - mmst2,1) * (
          - 8/3*mmsusy
          )

       + lnMglSq*lnMmsusy*den(mmgl - mmst2,2)*sn2t * (
          - 32*mmst2*mmsusy/mt*mgl
          + 64/3*sqr(mmsusy)/mt*mgl
          )

       + lnMglSq*lnMmsusy*den(mmgl - mmst2,2) * (
          - 104/3*mmst2*mmsusy
          + 16*sqr(mmsusy)
          )

       + lnMglSq*lnMmsusy*den(mmgl - mmst2,3) * (
          + 64/3*mmst2*sqr(mmsusy)
          - 32*sqr(mmst2)*mmsusy
          )

    + lnMglSq*lnMmt*den(mmgl - mmst1,1)*sn2t * (
          + 16/3*mmst1/mt*mgl
          )

       + lnMglSq*lnMmt*den(mmgl - mmst1,1) * (
          - 16/3*mmst1
          )

       + lnMglSq*lnMmt*den(mmgl - mmst1,2) * (
          - 8/3*sqr(mmst1)
          )

       + lnMglSq*lnMmt*den(mmgl - mmst2,1)*sn2t * (
          - 16/3*mmst2/mt*mgl
          )

       + lnMglSq*lnMmt*den(mmgl - mmst2,1) * (
          - 16/3*mmst2
          )

       + lnMglSq*lnMmt*den(mmgl - mmst2,2) * (
          - 8/3*sqr(mmst2)
          )

       + lnMglSq*lnMmt * (
          - 40/3
          )

       + lnMglSq*lnMmu*den(mmgl - mmst1,1)*sn2t * (
          - 16*mmst1/mt*mgl
          + 16/9*mmst2/mt*mgl
          )

       + lnMglSq*lnMmu*den(mmgl - mmst1,1)*sn4t*cs2t * (
          - 8/9*mmst1/mt*mgl
          - 8/9*mmst2/mt*mgl
          )

       + lnMglSq*lnMmu*den(mmgl - mmst1,1)*sqr(cs2t) * (
          - 16*mmst1
          + 16/9*mmst2
          )

       + lnMglSq*lnMmu*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sqr(cs2t)
       * (
          + 128/9*sqr(mmst1)
          )

       + lnMglSq*lnMmu*den(mmgl - mmst1,1) * (
          + 332/9*mmst1
          - 16/9*mmst2
          )

       + lnMglSq*lnMmu*den(mmgl - mmst1,2)*sn2t * (
          + 16/9*mmst1*mmst2/mt*mgl
          - 8/9*sqr(mmst1)/mt*mgl
          )

       + lnMglSq*lnMmu*den(mmgl - mmst1,2)*sn4t*cs2t * (
          - 8/9*mmst1*mmst2/mt*mgl
          + 8/9*sqr(mmst1)/mt*mgl
          )

       + lnMglSq*lnMmu*den(mmgl - mmst1,2)*sqr(cs2t) * (
          + 32/9*mmst1*mmst2
          - 32/3*sqr(mmst1)
          )

       + lnMglSq*lnMmu*den(mmgl - mmst1,2) * (
          - 32/9*mmst1*mmst2
          + 178/9*sqr(mmst1)
          )

       + lnMglSq*lnMmu*den(mmgl - mmst1,3)*sqr(cs2t) * (
          + 16/9*sqr(mmst1)*mmst2
          - 16/9*pow(mmst1,3)
          )

       + lnMglSq*lnMmu*den(mmgl - mmst1,3) * (
          - 16/9*sqr(mmst1)*mmst2
          + 8/9*pow(mmst1,3)
          )

       + lnMglSq*lnMmu*den(mmgl - mmst2,1)*sn2t * (
          - 16/9*mmst1/mt*mgl
          + 16*mmst2/mt*mgl
          )

       + lnMglSq*lnMmu*den(mmgl - mmst2,1)*sn4t*cs2t * (
          + 8/9*mmst1/mt*mgl
          + 8/9*mmst2/mt*mgl
          )

       + lnMglSq*lnMmu*den(mmgl - mmst2,1)*sqr(cs2t) * (
          + 16*mmst1
          - 16/9*mmst2
          )

       + lnMglSq*lnMmu*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sqr(cs2t)
       * (
          - 128/9*sqr(mmst1)
          )

       + lnMglSq*lnMmu*den(mmgl - mmst2,1) * (
          - 16/9*mmst1
          + 332/9*mmst2
          )

       + lnMglSq*lnMmu*den(mmgl - mmst2,2)*sn2t * (
          - 16/9*mmst1*mmst2/mt*mgl
          + 8/9*sqr(mmst2)/mt*mgl
          )

       + lnMglSq*lnMmu*den(mmgl - mmst2,2)*sn4t*cs2t * (
          + 8/9*mmst1*mmst2/mt*mgl
          - 8/9*sqr(mmst2)/mt*mgl
          )

       + lnMglSq*lnMmu*den(mmgl - mmst2,2)*sqr(cs2t) * (
          + 32/9*mmst1*mmst2
          - 32/3*sqr(mmst2)
          )

       + lnMglSq*lnMmu*den(mmgl - mmst2,2) * (
          - 32/9*mmst1*mmst2
          + 178/9*sqr(mmst2)
          )

       + lnMglSq*lnMmu*den(mmgl - mmst2,3)*sqr(cs2t) * (
          + 16/9*mmst1*sqr(mmst2)
          - 16/9*pow(mmst2,3)
          )

       + lnMglSq*lnMmu*den(mmgl - mmst2,3) * (
          - 16/9*mmst1*sqr(mmst2)
          + 8/9*pow(mmst2,3)
          )

       + lnMglSq*lnMmu * (
          + 36
          )

       + lnMglSq*den(mmgl - mmst1,1)*sn2t * (
          - 8/3*mmsb1/mt*mgl
          - 8/3*mmsb2/mt*mgl
          - 176/3*mmst1/mt*mgl
          - 8/9*mmst2/mt*mgl
          + 128/3*mmsusy/mt*mgl
          )

       + lnMglSq*den(mmgl - mmst1,1)*sn4t*cs2t * (
          - 8/9*mmst1/mt*mgl
          - 8/9*mmst2/mt*mgl
          )

       + lnMglSq*den(mmgl - mmst1,1)*sqr(cs2t) * (
          - 428/9*mmst1
          + 16/9*mmst2
          )

       + lnMglSq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          + 796/9*sqr(mmst1)
          )

       + lnMglSq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1) * (
          - 76/3*sqr(mmst1)
          )

       + lnMglSq*den(mmgl - mmst1,1) * (
          + 2*mmsb1
          + 2*mmsb2
          + 290/9*mmst1
          + 2/9*mmst2
          - 16*mmsusy
          )

       + lnMglSq*den(mmgl - mmst1,2)*sn2t * (
          - 8/3*mmsb1*mmst1/mt*mgl
          - 8/3*mmsb2*mmst1/mt*mgl
          - 8/9*mmst1*mmst2/mt*mgl
          - 64/3*mmst1*mmsusy/mt*mgl
          + 56/9*sqr(mmst1)/mt*mgl
          + 32*sqr(mmsusy)/mt*mgl
          )

       + lnMglSq*den(mmgl - mmst1,2)*sn4t*cs2t * (
          - 8/9*mmst1*mmst2/mt*mgl
          + 8/9*sqr(mmst1)/mt*mgl
          )

       + lnMglSq*den(mmgl - mmst1,2)*sqr(cs2t) * (
          + 32/9*mmst1*mmst2
          - 238/9*sqr(mmst1)
          )

       + lnMglSq*den(mmgl - mmst1,2)*den(mmst1 - mmst2,1) * (
          - 8/9*pow(mmst1,3)
          )

       + lnMglSq*den(mmgl - mmst1,2) * (
          + 14/3*mmsb1*mmst1
          + 14/3*mmsb2*mmst1
          + 10/9*mmst1*mmst2
          + 16/3*mmst1*mmsusy
          - 290/9*sqr(mmst1)
          - 24*sqr(mmsusy)
          )

       + lnMglSq*den(mmgl - mmst1,3)*sqr(cs2t) * (
          + 16/9*sqr(mmst1)*mmst2
          - 16/9*pow(mmst1,3)
          )

       + lnMglSq*den(mmgl - mmst1,3) * (
          + 8/3*mmsb1*sqr(mmst1)
          + 8/3*mmsb2*sqr(mmst1)
          - 32*mmst1*sqr(mmsusy)
          + 8/9*sqr(mmst1)*mmst2
          + 64/3*sqr(mmst1)*mmsusy
          - 200/9*pow(mmst1,3)
          )

       + lnMglSq*den(mmgl - mmst2,1)*sn2t * (
          + 8/3*mmsb1/mt*mgl
          + 8/3*mmsb2/mt*mgl
          + 8/9*mmst1/mt*mgl
          + 176/3*mmst2/mt*mgl
          - 128/3*mmsusy/mt*mgl
          )

       + lnMglSq*den(mmgl - mmst2,1)*sn4t*cs2t * (
          + 8/9*mmst1/mt*mgl
          + 8/9*mmst2/mt*mgl
          )

       + lnMglSq*den(mmgl - mmst2,1)*sqr(cs2t) * (
          + 812/9*mmst1
          + 368/9*mmst2
          )

       + lnMglSq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          - 796/9*sqr(mmst1)
          )

       + lnMglSq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1) * (
          + 76/3*sqr(mmst1)
          )

       + lnMglSq*den(mmgl - mmst2,1) * (
          + 2*mmsb1
          + 2*mmsb2
          - 226/9*mmst1
          + 62/9*mmst2
          - 16*mmsusy
          )

       + lnMglSq*den(mmgl - mmst2,2)*sn2t * (
          + 8/3*mmsb1*mmst2/mt*mgl
          + 8/3*mmsb2*mmst2/mt*mgl
          + 8/9*mmst1*mmst2/mt*mgl
          + 64/3*mmst2*mmsusy/mt*mgl
          - 56/9*sqr(mmst2)/mt*mgl
          - 32*sqr(mmsusy)/mt*mgl
          )

       + lnMglSq*den(mmgl - mmst2,2)*sn4t*cs2t * (
          + 8/9*mmst1*mmst2/mt*mgl
          - 8/9*sqr(mmst2)/mt*mgl
          )

       + lnMglSq*den(mmgl - mmst2,2)*sqr(cs2t) * (
          + 32/9*mmst1*mmst2
          - 238/9*sqr(mmst2)
          )

       + lnMglSq*den(mmgl - mmst2,2)*den(mmst1 - mmst2,1) * (
          + 8/9*pow(mmst1,3)
          )

       + lnMglSq*den(mmgl - mmst2,2) * (
          + 14/3*mmsb1*mmst2
          + 14/3*mmsb2*mmst2
          + 2/9*mmst1*mmst2
          - 8/9*sqr(mmst1)
          + 16/3*mmst2*mmsusy
          - 298/9*sqr(mmst2)
          - 24*sqr(mmsusy)
          )

       + lnMglSq*den(mmgl - mmst2,3)*sqr(cs2t) * (
          + 16/9*mmst1*sqr(mmst2)
          - 16/9*pow(mmst2,3)
          )

       + lnMglSq*den(mmgl - mmst2,3) * (
          + 8/3*mmsb1*sqr(mmst2)
          + 8/3*mmsb2*sqr(mmst2)
          + 8/9*mmst1*sqr(mmst2)
          - 32*mmst2*sqr(mmsusy)
          + 64/3*sqr(mmst2)*mmsusy
          - 200/9*pow(mmst2,3)
          )

       + lnMglSq * (
          + 232/3
          )

       + sqr(lnMsbSq)*den(mmgl - mmst1,1)*sn2t * (
          + 2/3*mmsb1/mt*mgl
          )

       + sqr(lnMsbSq)*den(mmgl - mmst1,1) * (
          - 7/6*mmsb1
          + 2/3*mmst1
          )

       + sqr(lnMsbSq)*den(mmgl - mmst1,2)*sn2t * (
          + 8/3*mmsb1*mmst1/mt*mgl
          - 4/3*sqr(mmsb1)/mt*mgl
          - 4/3*sqr(mmst1)/mt*mgl
          )

       + sqr(lnMsbSq)*den(mmgl - mmst1,2) * (
          - 4*mmsb1*sqr(mmst1
          + mmsb1)
          + 7/3*sqr(mmst1)
          )

       + sqr(lnMsbSq)*den(mmgl - mmst1,3) * (
          - 8/3*mmsb1*sqr(mmst1)
          + 4/3*sqr(mmsb1)*mmst1
          + 4/3*pow(mmst1,3)
          )

       + sqr(lnMsbSq)*den(mmgl - mmst2,1)*sn2t * (
          - 2/3*mmsb1/mt*mgl
          )

       + sqr(lnMsbSq)*den(mmgl - mmst2,1) * (
          - 7/6*mmsb1
          + 2/3*mmst2
          )

       + sqr(lnMsbSq)*den(mmgl - mmst2,2)*sn2t * (
          - 8/3*mmsb1*mmst2/mt*mgl
          + 4/3*sqr(mmsb1)/mt*mgl
          + 4/3*sqr(mmst2)/mt*mgl
          )

       + sqr(lnMsbSq)*den(mmgl - mmst2,2) * (
          - 4*mmsb1*sqr(mmst2
          + mmsb1)
          + 7/3*sqr(mmst2)
          )

       + sqr(lnMsbSq)*den(mmgl - mmst2,3) * (
          - 8/3*mmsb1*sqr(mmst2)
          + 4/3*sqr(mmsb1)*mmst2
          + 4/3*pow(mmst2,3)
          )

       + sqr(lnMsbSq) * (
          - 1/3
          )

       + lnMsbSq*lnMst1Sq*den(mmgl - mmst1,1) * (
          + 4/3*mmsb1
          )

       + lnMsbSq*lnMst1Sq*den(mmgl - mmst1,2)*sn2t * (
          - 4*mmsb1*mmst1/mt*mgl
          + 8/3*sqr(mmsb1)/mt*mgl
          )

       + lnMsbSq*lnMst1Sq*den(mmgl - mmst1,2) * (
          + 17/3*mmsb1*mmst1
          - 2*sqr(mmsb1)
          )

       + lnMsbSq*lnMst1Sq*den(mmgl - mmst1,3) * (
          + 4*mmsb1*sqr(mmst1)
          - 8/3*sqr(mmsb1)*mmst1
          )

       + lnMsbSq*lnMst2Sq*den(mmgl - mmst2,1) * (
          + 4/3*mmsb1
          )

       + lnMsbSq*lnMst2Sq*den(mmgl - mmst2,2)*sn2t * (
          + 4*mmsb1*mmst2/mt*mgl
          - 8/3*sqr(mmsb1)/mt*mgl
          )

       + lnMsbSq*lnMst2Sq*den(mmgl - mmst2,2) * (
          + 17/3*mmsb1*mmst2
          - 2*sqr(mmsb1)
          )

       + lnMsbSq*lnMst2Sq*den(mmgl - mmst2,3) * (
          + 4*mmsb1*sqr(mmst2)
          - 8/3*sqr(mmsb1)*mmst2
          )

       + lnMsbSq*lnMmt * (
          - 2/3
          )

       + lnMsbSq*den(mmgl - mmst1,1)*sn2t * (
          + 8/3*mmsb1/mt*mgl
          )

       + lnMsbSq*den(mmgl - mmst1,1) * (
          + 2*mmst1
          )

       + lnMsbSq*den(mmgl - mmst1,2)*sn2t * (
          + 4*sqr(mmsb1)/mt*mgl
          - 4*sqr(mmst1)/mt*mgl
          )

       + lnMsbSq*den(mmgl - mmst1,2) * (
          + 4/3*mmsb1*mmst1
          - 3*sqr(mmsb1)
          + 7*sqr(mmst1)
          )

       + lnMsbSq*den(mmgl - mmst1,3) * (
          - 4*sqr(mmsb1)*mmst1
          + 4*pow(mmst1,3)
          )

       + lnMsbSq*den(mmgl - mmst2,1)*sn2t * (
          - 8/3*mmsb1/mt*mgl
          )

       + lnMsbSq*den(mmgl - mmst2,1) * (
          + 2*mmst2
          )

       + lnMsbSq*den(mmgl - mmst2,2)*sn2t * (
          - 4*sqr(mmsb1)/mt*mgl
          + 4*sqr(mmst2)/mt*mgl
          )

       + lnMsbSq*den(mmgl - mmst2,2) * (
          + 4/3*mmsb1*mmst2
          - 3*sqr(mmsb1)
          + 7*sqr(mmst2)
          )

       + lnMsbSq*den(mmgl - mmst2,3) * (
          - 4*sqr(mmsb1)*mmst2
          + 4*pow(mmst2,3)
          )

       + lnMsbSq * (
          + 1/9
          )

       + sqr(lnMsb2Sq)*den(mmgl - mmst1,1)*sn2t * (
          + 2/3*mmsb2/mt*mgl
          )

       + sqr(lnMsb2Sq)*den(mmgl - mmst1,1) * (
          - 7/6*mmsb2
          + 2/3*mmst1
          )

       + sqr(lnMsb2Sq)*den(mmgl - mmst1,2)*sn2t * (
          + 8/3*mmsb2*mmst1/mt*mgl
          - 4/3*sqr(mmsb2)/mt*mgl
          - 4/3*sqr(mmst1)/mt*mgl
          )

       + sqr(lnMsb2Sq)*den(mmgl - mmst1,2) * (
          - 4*mmsb2*sqr(mmst1
          + mmsb2)
          + 7/3*sqr(mmst1)
          )

       + sqr(lnMsb2Sq)*den(mmgl - mmst1,3) * (
          - 8/3*mmsb2*sqr(mmst1)
          + 4/3*sqr(mmsb2)*mmst1
          + 4/3*pow(mmst1,3)
          )

       + sqr(lnMsb2Sq)*den(mmgl - mmst2,1)*sn2t * (
          - 2/3*mmsb2/mt*mgl
          )

       + sqr(lnMsb2Sq)*den(mmgl - mmst2,1) * (
          - 7/6*mmsb2
          + 2/3*mmst2
          )

       + sqr(lnMsb2Sq)*den(mmgl - mmst2,2)*sn2t * (
          - 8/3*mmsb2*mmst2/mt*mgl
          + 4/3*sqr(mmsb2)/mt*mgl
          + 4/3*sqr(mmst2)/mt*mgl
          )

       + sqr(lnMsb2Sq)*den(mmgl - mmst2,2) * (
          - 4*mmsb2*sqr(mmst2
          + mmsb2)
          + 7/3*sqr(mmst2)
          )

       + sqr(lnMsb2Sq)*den(mmgl - mmst2,3) * (
          - 8/3*mmsb2*sqr(mmst2)
          + 4/3*sqr(mmsb2)*mmst2
          + 4/3*pow(mmst2,3)
          )

       + sqr(lnMsb2Sq) * (
          - 1/3
          )

       + lnMsb2Sq*lnMst1Sq*den(mmgl - mmst1,1) * (
          + 4/3*mmsb2
          )

       + lnMsb2Sq*lnMst1Sq*den(mmgl - mmst1,2)*sn2t * (
          - 4*mmsb2*mmst1/mt*mgl
          + 8/3*sqr(mmsb2)/mt*mgl
          )

       + lnMsb2Sq*lnMst1Sq*den(mmgl - mmst1,2) * (
          + 17/3*mmsb2*mmst1
          - 2*sqr(mmsb2)
          )

       + lnMsb2Sq*lnMst1Sq*den(mmgl - mmst1,3) * (
          + 4*mmsb2*sqr(mmst1)
          - 8/3*sqr(mmsb2)*mmst1
          )

       + lnMsb2Sq*lnMst2Sq*den(mmgl - mmst2,1) * (
          + 4/3*mmsb2
          )

       + lnMsb2Sq*lnMst2Sq*den(mmgl - mmst2,2)*sn2t * (
          + 4*mmsb2*mmst2/mt*mgl
          - 8/3*sqr(mmsb2)/mt*mgl
          )

       + lnMsb2Sq*lnMst2Sq*den(mmgl - mmst2,2) * (
          + 17/3*mmsb2*mmst2
          - 2*sqr(mmsb2)
          )

       + lnMsb2Sq*lnMst2Sq*den(mmgl - mmst2,3) * (
          + 4*mmsb2*sqr(mmst2)
          - 8/3*sqr(mmsb2)*mmst2
          )

       + lnMsb2Sq*lnMmt * (
          - 2/3
          )

       + lnMsb2Sq*den(mmgl - mmst1,1)*sn2t * (
          + 8/3*mmsb2/mt*mgl
          )

       + lnMsb2Sq*den(mmgl - mmst1,1) * (
          + 2*mmst1
          )

       + lnMsb2Sq*den(mmgl - mmst1,2)*sn2t * (
          + 4*sqr(mmsb2)/mt*mgl
          - 4*sqr(mmst1)/mt*mgl
          )

       + lnMsb2Sq*den(mmgl - mmst1,2) * (
          + 4/3*mmsb2*mmst1
          - 3*sqr(mmsb2)
          + 7*sqr(mmst1)
          )

       + lnMsb2Sq*den(mmgl - mmst1,3) * (
          - 4*sqr(mmsb2)*mmst1
          + 4*pow(mmst1,3)
          )

       + lnMsb2Sq*den(mmgl - mmst2,1)*sn2t * (
          - 8/3*mmsb2/mt*mgl
          )

       + lnMsb2Sq*den(mmgl - mmst2,1) * (
          + 2*mmst2
          )

       + lnMsb2Sq*den(mmgl - mmst2,2)*sn2t * (
          - 4*sqr(mmsb2)/mt*mgl
          + 4*sqr(mmst2)/mt*mgl
          )

       + lnMsb2Sq*den(mmgl - mmst2,2) * (
          + 4/3*mmsb2*mmst2
          - 3*sqr(mmsb2)
          + 7*sqr(mmst2)
          )

       + lnMsb2Sq*den(mmgl - mmst2,3) * (
          - 4*sqr(mmsb2)*mmst2
          + 4*pow(mmst2,3)
          )

       + lnMsb2Sq * (
          + 1/9
          )

       + lnMst1Sq*sn2t * (
          + 280/9/mt*mgl
          )

       + lnMst1Sq*sqr(cs2t) * (
          + 64/9
          )

       + sqr(lnMst1Sq)*sn2t * (
          + 8/mt*mgl
          )

       + sqr(lnMst1Sq)*den(mmgl - mmst1,1)*sn2t * (
          - 2/9*mmst1/mt*mgl
          - 4/9*mmst2/mt*mgl
          )

       + sqr(lnMst1Sq)*den(mmgl - mmst1,1)*sn4t*cs2t * (
          + 8/9*mmst1/mt*mgl
          )

       + sqr(lnMst1Sq)*den(mmgl - mmst1,1)*sqr(cs2t) * (
          - 14/9*mmst1
          - 11/9*mmst2
          )

       + sqr(lnMst1Sq)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sn4t*cs2t * (
          - 16/9*sqr(mmst1)/mt*mgl
          )

       + sqr(lnMst1Sq)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          - 77/9*sqr(mmst1)
          )

       + sqr(lnMst1Sq)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1) * (
          + 13/9*sqr(mmst1)
          )

       + sqr(lnMst1Sq)*den(mmgl - mmst1,1) * (
          - 2/3*mmsb1
          - 2/3*mmsb2
          + 419/18*mmst1
          + mmst2
          - 16/3*mmsusy
          )

       + sqr(lnMst1Sq)*den(mmgl - mmst1,2)*sn2t * (
          + 2*mmsb1*mmst1/mt*mgl
          - 4/3*sqr(mmsb1)/mt*mgl
          + 2*mmsb2*mmst1/mt*mgl
          - 4/3*sqr(mmsb2)/mt*mgl
          + 2*mmst1*mmst2/mt*mgl
          + 16*mmst1*mmsusy/mt*mgl
          - 86/9*sqr(mmst1)/mt*mgl
          - 4/3*sqr(mmst2)/mt*mgl
          - 32/3*sqr(mmsusy)/mt*mgl
          )

       + sqr(lnMst1Sq)*den(mmgl - mmst1,2)*sn4t*cs2t * (
          + 8/9*sqr(mmst1)/mt*mgl
          )

       + sqr(lnMst1Sq)*den(mmgl - mmst1,2)*sqr(cs2t) * (
          - 26/9*mmst1*mmst2
          - 85/18*sqr(mmst1)
          )

       + sqr(lnMst1Sq)*den(mmgl - mmst1,2) * (
          - 17/6*mmsb1*sqr(mmst1
          + mmsb1)
          - 17/6*mmsb2*sqr(mmst1
          + mmsb2)
          + 1/18*mmst1*mmst2
          - 68/3*mmst1*mmsusy
          + 92/3*sqr(sqr(mmst1)
          + mmst2)
          + 8*sqr(mmsusy)
          )

       + sqr(lnMst1Sq)*den(mmgl - mmst1,3)*sn2t * (
          - 8/9*pow(mmst1,3)/mt*mgl
          )

       + sqr(lnMst1Sq)*den(mmgl - mmst1,3)*sqr(cs2t) * (
          - 16/9*pow(mmst1,3)
          )

       + sqr(lnMst1Sq)*den(mmgl - mmst1,3) * (
          - 2*mmsb1*sqr(mmst1)
          + 4/3*sqr(mmsb1)*mmst1
          - 2*mmsb2*sqr(mmst1)
          + 4/3*sqr(mmsb2)*mmst1
          + 4/3*mmst1*sqr(mmst2)
          + 32/3*mmst1*sqr(mmsusy)
          - 2*sqr(mmst1)*mmst2
          - 16*sqr(mmst1)*mmsusy
          + 110/9*pow(mmst1,3)
          )

       + sqr(lnMst1Sq)*den(mmgl - mmst1,4) * (
          + 4/9*pow(mmst1,4)
          )

       + sqr(lnMst1Sq)*den(mmgl - mmst2,1)*sn2t * (
          - 2/3*mmst1/mt*mgl
          )

       + sqr(lnMst1Sq)*den(mmgl - mmst2,1)*sqr(cs2t) * (
          - 13/9*mmst1
          )

       + sqr(lnMst1Sq)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          + 13/9*sqr(mmst1)
          )

       + sqr(lnMst1Sq)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1) * (
          - 13/9*sqr(mmst1)
          )

       + sqr(lnMst1Sq)*den(mmgl - mmst2,1) * (
          + 17/18*mmst1
          )

       + sqr(lnMst1Sq)*den(mmgl - mmst2,2) * (
          - 2/3*mmst1*mmst2
          )

       + sqr(lnMst1Sq)*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          - 128/9*mmgl
          - 64/9*mmst1
          )

       + sqr(lnMst1Sq) * (
          + 41/9
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst1,1)*sn2t * (
          - 8/9*mmst1/mt*mgl
          + 8/3*mmst2/mt*mgl
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst1,1)*sn4t*cs2t * (
          - 16/9*mmst1/mt*mgl
          - 8/9*mmst2/mt*mgl
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst1,1)*sqr(cs2t) * (
          + 5/9*mmst1
          + 38/9*mmst2
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sn2t
       * (
          + 8/3*sqr(mmst1)/mt*mgl
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sn4t*
      cs2t * (
          + 16/9*sqr(mmst1)/mt*mgl
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sqr(cs2t)
       * (
          - 5/9*sqr(mmst1)
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1) * (
          + 17/9*sqr(mmst1)
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,2)*sn2t
       * (
          - 16/9*pow(mmst1,3)/mt*mgl
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,2) * (
          + 32/9*pow(mmst1,3)
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,3) * (
          - 16/9*pow(mmst1,4)
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst1,1) * (
          - 11/3*mmst1
          - 34/9*mmst2
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst1,2)*sn2t * (
          - 20/9*mmst1*mmst2/mt*mgl
          - 8/9*sqr(mmst1)/mt*mgl
          + 8/3*sqr(mmst2)/mt*mgl
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst1,2)*sn4t*cs2t * (
          - 8/9*mmst1*mmst2/mt*mgl
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst1,2)*sqr(cs2t) * (
          + 28/3*mmst1*mmst2
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst1,2)*den(mmst1 - mmst2,1)*sn2t
       * (
          + 8/9*pow(mmst1,3)/mt*mgl
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst1,2)*den(mmst1 - mmst2,2) * (
          + 8/9*pow(mmst1,4)
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst1,2) * (
          - 11/3*mmst1*mmst2
          - 8/9*sqr(mmst1)
          - 2*sqr(mmst2)
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst1,3)*sqr(cs2t) * (
          + 16/9*sqr(mmst1)*mmst2
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst1,3) * (
          - 8/3*mmst1*sqr(mmst2)
          + 20/9*sqr(mmst1)*mmst2
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst2,1)*sn2t * (
          - 8/9*mmst1/mt*mgl
          - 8/9*mmst2/mt*mgl
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst2,1)*sn4t*cs2t * (
          - 8/9*mmst1/mt*mgl
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst2,1)*sqr(cs2t) * (
          + 11/9*mmst1
          + 22/9*mmst2
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sn2t
       * (
          - 8/3*sqr(mmst1)/mt*mgl
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sn4t*
      cs2t * (
          + 16/9*sqr(mmst1)/mt*mgl
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sqr(cs2t)
       * (
          + 5/9*sqr(mmst1)
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1) * (
          - 17/9*sqr(mmst1)
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,2)*sn2t
       * (
          + 16/9*pow(mmst1,3)/mt*mgl
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,2) * (
          - 32/9*pow(mmst1,3)
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,3) * (
          + 16/9*pow(mmst1,4)
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst2,1) * (
          + 17/9*mmst1
          - 2/9*mmst2
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst2,2)*sn2t * (
          - 4*mmst1*mmst2/mt*mgl
          - 8/9*sqr(mmst1)/mt*mgl
          + 8/3*sqr(mmst2)/mt*mgl
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst2,2)*sn4t*cs2t * (
          + 8/9*mmst1*mmst2/mt*mgl
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst2,2)*sqr(cs2t) * (
          + 32/9*mmst1*mmst2
          + 52/9*sqr(mmst2)
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst2,2)*den(mmst1 - mmst2,1)*sn2t
       * (
          + 8/9*pow(mmst1,3)/mt*mgl
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst2,2)*den(mmst1 - mmst2,1) * (
          - 32/9*pow(mmst1,3)
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst2,2)*den(mmst1 - mmst2,2) * (
          + 8/9*pow(mmst1,4)
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst2,2) * (
          - 25/9*mmst1*mmst2
          + 8/3*sqr(mmst1)
          - 10/9*sqr(mmst2)
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst2,3)*sqr(cs2t) * (
          + 16/9*mmst1*sqr(mmst2)
          )

       + lnMst1Sq*lnMst2Sq*den(mmgl - mmst2,3) * (
          - 28/9*mmst1*sqr(mmst2)
          + 8/3*pow(mmst2,3)
          )

       + lnMst1Sq*lnMmsusy*den(mmgl - mmst1,1) * (
          + 32/3*mmsusy
          )

       + lnMst1Sq*lnMmsusy*den(mmgl - mmst1,2)*sn2t * (
          - 32*mmst1*mmsusy/mt*mgl
          + 64/3*sqr(mmsusy)/mt*mgl
          )

       + lnMst1Sq*lnMmsusy*den(mmgl - mmst1,2) * (
          + 136/3*mmst1*mmsusy
          - 16*sqr(mmsusy)
          )

       + lnMst1Sq*lnMmsusy*den(mmgl - mmst1,3) * (
          - 64/3*mmst1*sqr(mmsusy)
          + 32*sqr(mmst1)*mmsusy
          )

       + lnMst1Sq*lnMmt*den(mmgl - mmst1,1)*sn2t * (
          - 16/3*mmst1/mt*mgl
          )

       + lnMst1Sq*lnMmt*den(mmgl - mmst1,1) * (
          + 16/3*mmst1
          )

       + lnMst1Sq*lnMmt*den(mmgl - mmst1,2) * (
          + 8/3*sqr(mmst1)
          )

       + lnMst1Sq*lnMmt * (
          - 2/3
          )

       + lnMst1Sq*lnMmu*sn2t * (
          + 64/9/mt*mgl
          )

       + lnMst1Sq*lnMmu*sqr(cs2t) * (
          + 64/9
          )

       + lnMst1Sq*lnMmu*den(mmgl - mmst1,1)*sn2t * (
          + 16*mmst1/mt*mgl
          - 16/9*mmst2/mt*mgl
          )

       + lnMst1Sq*lnMmu*den(mmgl - mmst1,1)*sn4t*cs2t * (
          + 8/9*mmst1/mt*mgl
          + 8/9*mmst2/mt*mgl
          )

       + lnMst1Sq*lnMmu*den(mmgl - mmst1,1)*sqr(cs2t) * (
          + 16*mmst1
          - 16/9*mmst2
          )

       + lnMst1Sq*lnMmu*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sqr(cs2t)
       * (
          - 128/9*sqr(mmst1)
          )

       + lnMst1Sq*lnMmu*den(mmgl - mmst1,1) * (
          - 332/9*mmst1
          + 16/9*mmst2
          )

       + lnMst1Sq*lnMmu*den(mmgl - mmst1,2)*sn2t * (
          - 16/9*mmst1*mmst2/mt*mgl
          + 8/9*sqr(mmst1)/mt*mgl
          )

       + lnMst1Sq*lnMmu*den(mmgl - mmst1,2)*sn4t*cs2t * (
          + 8/9*mmst1*mmst2/mt*mgl
          - 8/9*sqr(mmst1)/mt*mgl
          )

       + lnMst1Sq*lnMmu*den(mmgl - mmst1,2)*sqr(cs2t) * (
          - 32/9*mmst1*mmst2
          + 32/3*sqr(mmst1)
          )

       + lnMst1Sq*lnMmu*den(mmgl - mmst1,2) * (
          + 32/9*mmst1*mmst2
          - 178/9*sqr(mmst1)
          )

       + lnMst1Sq*lnMmu*den(mmgl - mmst1,3)*sqr(cs2t) * (
          - 16/9*sqr(mmst1)*mmst2
          + 16/9*pow(mmst1,3)
          )

       + lnMst1Sq*lnMmu*den(mmgl - mmst1,3) * (
          + 16/9*sqr(mmst1)*mmst2
          - 8/9*pow(mmst1,3)
          )

       + lnMst1Sq*lnMmu*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          - 128/9*mmst1
          )

       + lnMst1Sq*lnMmu * (
          - 128/9
          )

       + lnMst1Sq*den(mmgl - mmst1,1)*sn2t * (
          + 172/3*mmst1/mt*mgl
          - 28/9*mmst2/mt*mgl
          )

       + lnMst1Sq*den(mmgl - mmst1,1)*sn4t*cs2t * (
          + 16/9*mmst1/mt*mgl
          + 8/9*mmst2/mt*mgl
          )

       + lnMst1Sq*den(mmgl - mmst1,1)*sqr(cs2t) * (
          + 229/9*mmst1
          - 49/9*mmst2
          )

       + lnMst1Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sn2t * (
          - 8/9*sqr(mmst1)/mt*mgl
          )

       + lnMst1Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          - 718/9*sqr(mmst1)
          )

       + lnMst1Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1) * (
          + 34/3*sqr(mmst1)
          )

       + lnMst1Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,2) * (
          - 8/9*pow(mmst1,3)
          )

       + lnMst1Sq*den(mmgl - mmst1,1) * (
          - 2*mmsb1
          - 2*mmsb2
          - 13/3*mmst1
          + 43/9*mmst2
          - 16*mmsusy
          )

       + lnMst1Sq*den(mmgl - mmst1,2)*sn2t * (
          + 8/3*mmsb1*mmst1/mt*mgl
          - 4*sqr(mmsb1)/mt*mgl
          + 8/3*mmsb2*mmst1/mt*mgl
          - 4*sqr(mmsb2)/mt*mgl
          + 8/9*mmst1*mmst2/mt*mgl
          + 64/3*mmst1*mmsusy/mt*mgl
          + 52/9*sqr(mmst1)/mt*mgl
          - 4*sqr(mmst2)/mt*mgl
          - 32*sqr(mmsusy)/mt*mgl
          )

       + lnMst1Sq*den(mmgl - mmst1,2)*sn4t*cs2t * (
          + 8/9*mmst1*mmst2/mt*mgl
          - 8/9*sqr(mmst1)/mt*mgl
          )

       + lnMst1Sq*den(mmgl - mmst1,2)*sqr(cs2t) * (
          - 110/9*mmst1*mmst2
          + 16*sqr(mmst1)
          )

       + lnMst1Sq*den(mmgl - mmst1,2)*den(mmst1 - mmst2,1) * (
          + 8/9*pow(mmst1,3)
          )

       + lnMst1Sq*den(mmgl - mmst1,2) * (
          - 6*mmsb1*mmst1
          + 3*sqr(mmsb1)
          - 6*mmsb2*mmst1
          + 3*sqr(mmsb2)
          + 56/9*mmst1*mmst2
          - 48*mmst1*mmsusy
          + 187/9*sqr(mmst1)
          + 3*sqr(mmst2)
          + 24*sqr(mmsusy)
          )

       + lnMst1Sq*den(mmgl - mmst1,3)*sqr(cs2t) * (
          - 16/9*sqr(mmst1)*mmst2
          + 16/9*pow(mmst1,3)
          )

       + lnMst1Sq*den(mmgl - mmst1,3) * (
          - 8/3*mmsb1*sqr(mmst1)
          + 4*sqr(mmsb1)*mmst1
          - 8/3*mmsb2*sqr(mmst1)
          + 4*sqr(mmsb2)*mmst1
          + 4*mmst1*sqr(mmst2)
          + 32*mmst1*sqr(mmsusy)
          - 8/9*sqr(mmst1)*mmst2
          - 64/3*sqr(mmst1)*mmsusy
          + 92/9*pow(mmst1,3)
          )

       + lnMst1Sq*den(mmgl - mmst2,1)*sn2t * (
          - 16/3*mmst1/mt*mgl
          )

       + lnMst1Sq*den(mmgl - mmst2,1)*sn4t*cs2t * (
          + 8/9*mmst1/mt*mgl
          )

       + lnMst1Sq*den(mmgl - mmst2,1)*sqr(cs2t) * (
          - 6*mmst1
          )

       + lnMst1Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sn2t * (
          + 8/9*sqr(mmst1)/mt*mgl
          )

       + lnMst1Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          + 26/3*sqr(mmst1)
          )

       + lnMst1Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1) * (
          - 34/3*sqr(mmst1)
          )

       + lnMst1Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,2) * (
          + 8/9*pow(mmst1,3)
          )

       + lnMst1Sq*den(mmgl - mmst2,1) * (
          + 52/9*mmst1
          )

       + lnMst1Sq*den(mmgl - mmst2,2)*sqr(cs2t) * (
          + 16/9*mmst1*mmst2
          )

       + lnMst1Sq*den(mmgl - mmst2,2) * (
          - 40/9*mmst1*mmst2
          )

       + lnMst1Sq*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          - 128/3*mmgl
          - 640/9*mmst1
          )

       + lnMst1Sq * (
          + 5/9
          )

       + lnMst2Sq*sn2t * (
          - 280/9/mt*mgl
          )

       + lnMst2Sq*sqr(cs2t) * (
          - 64
          )

       + sqr(lnMst2Sq)*sn2t * (
          - 8/mt*mgl
          )

       + sqr(lnMst2Sq)*sqr(cs2t) * (
          - 64/9
          )

       + sqr(lnMst2Sq)*den(mmgl - mmst1,1)*sn2t * (
          + 4/9*mmst1/mt*mgl
          + 2/9*mmst2/mt*mgl
          )

       + sqr(lnMst2Sq)*den(mmgl - mmst1,1)*sqr(cs2t) * (
          + 8/3*mmst1
          - 11/9*mmst2
          )

       + sqr(lnMst2Sq)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          - 13/9*sqr(mmst1)
          )

       + sqr(lnMst2Sq)*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1) * (
          + 13/9*sqr(mmst1)
          )

       + sqr(lnMst2Sq)*den(mmgl - mmst1,1) * (
          - 14/9*mmst1
          + 1/2*mmst2
          )

       + sqr(lnMst2Sq)*den(mmgl - mmst1,2)*sn2t * (
          + 8/3*mmst1*mmst2/mt*mgl
          - 4/3*sqr(mmst1)/mt*mgl
          - 4/3*sqr(mmst2)/mt*mgl
          )

       + sqr(lnMst2Sq)*den(mmgl - mmst1,2)*sqr(cs2t) * (
          - 26/9*mmst1*mmst2
          + 26/9*sqr(mmst1)
          )

       + sqr(lnMst2Sq)*den(mmgl - mmst1,2) * (
          - 10/9*mmst1*mmst2
          - 5/9*sqr(sqr(mmst1)
          + mmst2)
          )

       + sqr(lnMst2Sq)*den(mmgl - mmst1,3) * (
          + 4/3*mmst1*sqr(mmst2)
          - 8/3*sqr(mmst1)*mmst2
          + 4/3*pow(mmst1,3)
          )

       + sqr(lnMst2Sq)*den(mmgl - mmst2,1)*sn2t * (
          + 2/3*mmst2/mt*mgl
          )

       + sqr(lnMst2Sq)*den(mmgl - mmst2,1)*sn4t*cs2t * (
          + 16/9*mmst1/mt*mgl
          + 8/9*mmst2/mt*mgl
          )

       + sqr(lnMst2Sq)*den(mmgl - mmst2,1)*sqr(cs2t) * (
          - 77/9*mmst1
          - 34/3*mmst2
          )

       + sqr(lnMst2Sq)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sn4t*cs2t * (
          - 16/9*sqr(mmst1)/mt*mgl
          )

       + sqr(lnMst2Sq)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          + 77/9*sqr(mmst1)
          )

       + sqr(lnMst2Sq)*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1) * (
          - 13/9*sqr(mmst1)
          )

       + sqr(lnMst2Sq)*den(mmgl - mmst2,1) * (
          - 2/3*mmsb1
          - 2/3*mmsb2
          + 13/9*mmst1
          + 149/6*mmst2
          - 16/3*mmsusy
          )

       + sqr(lnMst2Sq)*den(mmgl - mmst2,2)*sn2t * (
          - 2*mmsb1*mmst2/mt*mgl
          + 4/3*sqr(mmsb1)/mt*mgl
          - 2*mmsb2*mmst2/mt*mgl
          + 4/3*sqr(mmsb2)/mt*mgl
          + 2/3*mmst1*mmst2/mt*mgl
          - 16*mmst2*mmsusy/mt*mgl
          + 74/9*sqr(mmst2)/mt*mgl
          + 32/3*sqr(mmsusy)/mt*mgl
          )

       + sqr(lnMst2Sq)*den(mmgl - mmst2,2)*sn4t*cs2t * (
          - 8/9*sqr(mmst2)/mt*mgl
          )

       + sqr(lnMst2Sq)*den(mmgl - mmst2,2)*sqr(cs2t) * (
          - 137/18*sqr(mmst2)
          )

       + sqr(lnMst2Sq)*den(mmgl - mmst2,2) * (
          - 17/6*mmsb1*sqr(mmst2
          + mmsb1)
          - 17/6*mmsb2*sqr(mmst2
          + mmsb2)
          + 1/2*mmst1*mmst2
          - 68/3*mmst2*mmsusy
          + 281/9*sqr(mmst2)
          + 8*sqr(mmsusy)
          )

       + sqr(lnMst2Sq)*den(mmgl - mmst2,3)*sn2t * (
          + 8/9*pow(mmst2,3)/mt*mgl
          )

       + sqr(lnMst2Sq)*den(mmgl - mmst2,3)*sqr(cs2t) * (
          - 16/9*pow(mmst2,3)
          )

       + sqr(lnMst2Sq)*den(mmgl - mmst2,3) * (
          - 2*mmsb1*sqr(mmst2)
          + 4/3*sqr(mmsb1)*mmst2
          - 2*mmsb2*sqr(mmst2)
          + 4/3*sqr(mmsb2)*mmst2
          + 2/3*mmst1*sqr(mmst2)
          + 32/3*mmst2*sqr(mmsusy)
          - 16*sqr(mmst2)*mmsusy
          + 98/9*pow(mmst2,3)
          )

       + sqr(lnMst2Sq)*den(mmgl - mmst2,4) * (
          + 4/9*pow(mmst2,4)
          )

       + sqr(lnMst2Sq)*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          + 128/9*mmgl
          + 64/9*mmst1
          )

       + sqr(lnMst2Sq) * (
          + 41/9
          )

       + lnMst2Sq*lnMmsusy*den(mmgl - mmst2,1) * (
          + 32/3*mmsusy
          )

       + lnMst2Sq*lnMmsusy*den(mmgl - mmst2,2)*sn2t * (
          + 32*mmst2*mmsusy/mt*mgl
          - 64/3*sqr(mmsusy)/mt*mgl
          )

       + lnMst2Sq*lnMmsusy*den(mmgl - mmst2,2) * (
          + 136/3*mmst2*mmsusy
          - 16*sqr(mmsusy)
          )

       + lnMst2Sq*lnMmsusy*den(mmgl - mmst2,3) * (
          - 64/3*mmst2*sqr(mmsusy)
          + 32*sqr(mmst2)*mmsusy
          )

       + lnMst2Sq*lnMmt*den(mmgl - mmst2,1)*sn2t * (
          + 16/3*mmst2/mt*mgl
          )

       + lnMst2Sq*lnMmt*den(mmgl - mmst2,1) * (
          + 16/3*mmst2
          )

       + lnMst2Sq*lnMmt*den(mmgl - mmst2,2) * (
          + 8/3*sqr(mmst2)
          )

       + lnMst2Sq*lnMmt * (
          - 2/3
          )

       + lnMst2Sq*lnMmu*sn2t * (
          - 64/9/mt*mgl
          )

       + lnMst2Sq*lnMmu*sqr(cs2t) * (
          - 64/9
          )

       + lnMst2Sq*lnMmu*den(mmgl - mmst2,1)*sn2t * (
          + 16/9*mmst1/mt*mgl
          - 16*mmst2/mt*mgl
          )

       + lnMst2Sq*lnMmu*den(mmgl - mmst2,1)*sn4t*cs2t * (
          - 8/9*mmst1/mt*mgl
          - 8/9*mmst2/mt*mgl
          )

       + lnMst2Sq*lnMmu*den(mmgl - mmst2,1)*sqr(cs2t) * (
          - 16*mmst1
          + 16/9*mmst2
          )

       + lnMst2Sq*lnMmu*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sqr(cs2t)
       * (
          + 128/9*sqr(mmst1)
          )

       + lnMst2Sq*lnMmu*den(mmgl - mmst2,1) * (
          + 16/9*mmst1
          - 332/9*mmst2
          )

       + lnMst2Sq*lnMmu*den(mmgl - mmst2,2)*sn2t * (
          + 16/9*mmst1*mmst2/mt*mgl
          - 8/9*sqr(mmst2)/mt*mgl
          )

       + lnMst2Sq*lnMmu*den(mmgl - mmst2,2)*sn4t*cs2t * (
          - 8/9*mmst1*mmst2/mt*mgl
          + 8/9*sqr(mmst2)/mt*mgl
          )

       + lnMst2Sq*lnMmu*den(mmgl - mmst2,2)*sqr(cs2t) * (
          - 32/9*mmst1*mmst2
          + 32/3*sqr(mmst2)
          )

       + lnMst2Sq*lnMmu*den(mmgl - mmst2,2) * (
          + 32/9*mmst1*mmst2
          - 178/9*sqr(mmst2)
          )

       + lnMst2Sq*lnMmu*den(mmgl - mmst2,3)*sqr(cs2t) * (
          - 16/9*mmst1*sqr(mmst2)
          + 16/9*pow(mmst2,3)
          )

       + lnMst2Sq*lnMmu*den(mmgl - mmst2,3) * (
          + 16/9*mmst1*sqr(mmst2)
          - 8/9*pow(mmst2,3)
          )

       + lnMst2Sq*lnMmu*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          + 128/9*mmst1
          )

       + lnMst2Sq*lnMmu * (
          - 128/9
          )

       + lnMst2Sq*den(mmgl - mmst1,1)*sn2t * (
          + 4/9*mmst1/mt*mgl
          + 52/9*mmst2/mt*mgl
          )

       + lnMst2Sq*den(mmgl - mmst1,1)*sn4t*cs2t * (
          - 8/9*mmst2/mt*mgl
          )

       + lnMst2Sq*den(mmgl - mmst1,1)*sqr(cs2t) * (
          + 37/3*mmst1
          + 19/3*mmst2
          )

       + lnMst2Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sn2t * (
          + 8/9*sqr(mmst1)/mt*mgl
          )

       + lnMst2Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          - 26/3*sqr(mmst1)
          )

       + lnMst2Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,1) * (
          + 14*sqr(mmst1)
          )

       + lnMst2Sq*den(mmgl - mmst1,1)*den(mmst1 - mmst2,2) * (
          + 8/9*pow(mmst1,3)
          )

       + lnMst2Sq*den(mmgl - mmst1,1) * (
          - 137/9*mmst1
          - 23/3*mmst2
          )

       + lnMst2Sq*den(mmgl - mmst1,2)*sn2t * (
          - 4*sqr(mmst1)/mt*mgl
          + 4*sqr(mmst2)/mt*mgl
          )

       + lnMst2Sq*den(mmgl - mmst1,2)*sqr(cs2t) * (
          + 94/9*mmst1*mmst2
          + 26/3*sqr(mmst1)
          )

       + lnMst2Sq*den(mmgl - mmst1,2) * (
          - 82/9*mmst1*mmst2
          - 5/3*sqr(mmst1)
          - 3*sqr(mmst2)
          )

       + lnMst2Sq*den(mmgl - mmst1,3) * (
          - 4*mmst1*sqr(mmst2)
          + 4*pow(mmst1,3)
          )

       + lnMst2Sq*den(mmgl - mmst2,1)*sn2t * (
          + 8/3*mmst1/mt*mgl
          - 520/9*mmst2/mt*mgl
          )

       + lnMst2Sq*den(mmgl - mmst2,1)*sn4t*cs2t * (
          - 8/9*mmst1/mt*mgl
          - 16/9*mmst2/mt*mgl
          )

       + lnMst2Sq*den(mmgl - mmst2,1)*sqr(cs2t) * (
          - 734/9*mmst1
          - 152/3*mmst2
          )

       + lnMst2Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sn2t * (
          - 8/9*sqr(mmst1)/mt*mgl
          )

       + lnMst2Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          + 718/9*sqr(mmst1)
          )

       + lnMst2Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,1) * (
          - 14*sqr(mmst1)
          )

       + lnMst2Sq*den(mmgl - mmst2,1)*den(mmst1 - mmst2,2) * (
          - 8/9*pow(mmst1,3)
          )

       + lnMst2Sq*den(mmgl - mmst2,1) * (
          - 2*mmsb1
          - 2*mmsb2
          + 50/3*mmst1
          + 52/9*mmst2
          - 16*mmsusy
          )

       + lnMst2Sq*den(mmgl - mmst2,2)*sn2t * (
          - 8/3*mmsb1*mmst2/mt*mgl
          + 4*sqr(mmsb1)/mt*mgl
          - 8/3*mmsb2*mmst2/mt*mgl
          + 4*sqr(mmsb2)/mt*mgl
          - 8/9*mmst1*mmst2/mt*mgl
          - 64/3*mmst2*mmsusy/mt*mgl
          - 16/9*sqr(mmst2)/mt*mgl
          + 32*sqr(mmsusy)/mt*mgl
          )

       + lnMst2Sq*den(mmgl - mmst2,2)*sn4t*cs2t * (
          - 8/9*mmst1*mmst2/mt*mgl
          + 8/9*sqr(mmst2)/mt*mgl
          )

       + lnMst2Sq*den(mmgl - mmst2,2)*sqr(cs2t) * (
          - 32/9*mmst1*mmst2
          + 74/3*sqr(mmst2)
          )

       + lnMst2Sq*den(mmgl - mmst2,2)*den(mmst1 - mmst2,1) * (
          - 8/9*pow(mmst1,3)
          )

       + lnMst2Sq*den(mmgl - mmst2,2) * (
          - 6*mmsb1*mmst2
          + 3*sqr(mmsb1)
          - 6*mmsb2*mmst2
          + 3*sqr(mmsb2)
          + 22/9*mmst1*mmst2
          + 8/9*sqr(mmst1)
          - 48*mmst2*mmsusy
          + 20*sqr(mmst2)
          + 24*sqr(mmsusy)
          )

       + lnMst2Sq*den(mmgl - mmst2,3)*sqr(cs2t) * (
          - 16/9*mmst1*sqr(mmst2)
          + 16/9*pow(mmst2,3)
          )

       + lnMst2Sq*den(mmgl - mmst2,3) * (
          - 8/3*mmsb1*sqr(mmst2)
          + 4*sqr(mmsb1)*mmst2
          - 8/3*mmsb2*sqr(mmst2)
          + 4*sqr(mmsb2)*mmst2
          - 8/9*mmst1*sqr(mmst2)
          + 32*mmst2*sqr(mmsusy)
          - 64/3*sqr(mmst2)*mmsusy
          + 128/9*pow(mmst2,3)
          )

       + lnMst2Sq*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          + 128/3*mmgl
          + 640/9*mmst1
          )

       + lnMst2Sq * (
          + 5/9
          )

       + sqr(lnMmsusy)*den(mmgl - mmst1,1)*sn2t * (
          + 16/3*mmsusy/mt*mgl
          )

       + sqr(lnMmsusy)*den(mmgl - mmst1,1) * (
          - 4*mmsusy
          )

       + sqr(lnMmsusy)*den(mmgl - mmst1,2) * (
          - 16/3*mmst1*mmsusy
          )

       + sqr(lnMmsusy)*den(mmgl - mmst2,1)*sn2t * (
          - 16/3*mmsusy/mt*mgl
          )

       + sqr(lnMmsusy)*den(mmgl - mmst2,1) * (
          - 4*mmsusy
          )

       + sqr(lnMmsusy)*den(mmgl - mmst2,2) * (
          - 16/3*mmst2*mmsusy
          )

       + sqr(lnMmsusy) * (
          + 8/3
          )

       + lnMmsusy*lnMmt * (
          - 16/3
          )

       + lnMmsusy*den(mmgl - mmst1,1)*sn2t * (
          - 128/3*mmsusy/mt*mgl
          )

       + lnMmsusy*den(mmgl - mmst1,1) * (
          + 32*mmsusy
          )

       + lnMmsusy*den(mmgl - mmst1,2) * (
          + 128/3*mmst1*mmsusy
          )

       + lnMmsusy*den(mmgl - mmst2,1)*sn2t * (
          + 128/3*mmsusy/mt*mgl
          )

       + lnMmsusy*den(mmgl - mmst2,1) * (
          + 32*mmsusy
          )

       + lnMmsusy*den(mmgl - mmst2,2) * (
          + 128/3*mmst2*mmsusy
          )

       + lnMmsusy * (
          + 152/9
          )

       + lnMmt*lnMmu * (
          + 64/3
          )

       + lnMmt*den(mmgl - mmst1,1) * (
          + 8/3*mmst1
          )

       + lnMmt*den(mmgl - mmst2,1) * (
          + 8/3*mmst2
          )

       + lnMmt * (
          + 8
          )

       + lnMmu*sqr(cs2t) * (
          + 128/9
          )

       + sqr(lnMmu) * (
          - 130/9
          )

       + lnMmu*den(mmgl - mmst1,1)*sn2t * (
          + 8/9*mmst1/mt*mgl
          - 16/9*mmst2/mt*mgl
          )

       + lnMmu*den(mmgl - mmst1,1)*sn4t*cs2t * (
          - 8/9*mmst1/mt*mgl
          + 8/9*mmst2/mt*mgl
          )

       + lnMmu*den(mmgl - mmst1,1)*sqr(cs2t) * (
          + 88/9*mmst1
          - 8/3*mmst2
          )

       + lnMmu*den(mmgl - mmst1,1) * (
          - 58/3*mmst1
          + 8/3*mmst2
          )

       + lnMmu*den(mmgl - mmst1,2)*sqr(cs2t) * (
          - 16/9*mmst1*mmst2
          + 16/9*sqr(mmst1)
          )

       + lnMmu*den(mmgl - mmst1,2) * (
          + 16/9*mmst1*mmst2
          - 8/9*sqr(mmst1)
          )

       + lnMmu*den(mmgl - mmst2,1)*sn2t * (
          + 16/9*mmst1/mt*mgl
          - 8/9*mmst2/mt*mgl
          )

       + lnMmu*den(mmgl - mmst2,1)*sn4t*cs2t * (
          - 8/9*mmst1/mt*mgl
          + 8/9*mmst2/mt*mgl
          )

       + lnMmu*den(mmgl - mmst2,1)*sqr(cs2t) * (
          - 8/3*mmst1
          + 88/9*mmst2
          )

       + lnMmu*den(mmgl - mmst2,1) * (
          + 8/3*mmst1
          - 58/3*mmst2
          )

       + lnMmu*den(mmgl - mmst2,2)*sqr(cs2t) * (
          - 16/9*mmst1*mmst2
          + 16/9*sqr(mmst2)
          )

       + lnMmu*den(mmgl - mmst2,2) * (
          + 16/9*mmst1*mmst2
          - 8/9*sqr(mmst2)
          )

       + lnMmu * (
          - 932/9
          )

       + den(mmgl - mmst1,1)*sn2t * (
          + 4/3*mmsb1/mt*mgl*zt2
          + 8*mmsb1/mt*mgl
          + 4/3*mmsb2/mt*mgl*zt2
          + 8*mmsb2/mt*mgl
          + 88/9*mmst1/mt*mgl*zt2
          + 676/9*mmst1/mt*mgl
          + 4/3*mmst2/mt*mgl*zt2
          + 56/9*mmst2/mt*mgl
          + 32/3*mmsusy/mt*mgl*zt2
          + 64*mmsusy/mt*mgl
          )

       + den(mmgl - mmst1,1)*sn4t*cs2t * (
          - 8/9*mmst1/mt*mgl
          + 8/9*mmst2/mt*mgl
          )

       + den(mmgl - mmst1,1)*sqr(cs2t) * (
          + 41/9*mmst1*zt2
          + 439/9*mmst1
          - 8/3*mmst2
          )

       + den(mmgl - mmst1,1)*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          - 20*sqr(mmst1)*zt2
          - 140*sqr(mmst1)
          )

       + den(mmgl - mmst1,1)*den(mmst1 - mmst2,1) * (
          + 20/3*sqr(mmst1)*zt2
          + 428/9*sqr(mmst1)
          )

       + den(mmgl - mmst1,1) * (
          - mmsb1*zt2
          - 6*mmsb1
          - mmsb2*zt2
          - 6*mmsb2
          + 50/9*mmst1*zt2
          + 61/9*mmst1
          - mmst2*zt2
          - 10/3*mmst2
          - 8*mmsusy*zt2
          - 48*mmsusy
          )

       + den(mmgl - mmst1,2)*sqr(cs2t) * (
          - 16/9*mmst1*mmst2
          + 16/9*sqr(mmst1)
          )

       + den(mmgl - mmst1,2) * (
          - 4/3*mmsb1*mmst1*zt2
          - 8*mmsb1*mmst1
          - 4/3*mmsb2*mmst1*zt2
          - 8*mmsb2*mmst1
          - 4/3*mmst1*mmst2*zt2
          - 56/9*mmst1*mmst2
          - 32/3*mmst1*mmsusy*zt2
          - 64*mmst1*mmsusy
          + 44/3*sqr(mmst1)*zt2
          + 868/9*sqr(mmst1)
          )

       + den(mmgl - mmst1,3) * (
          + 16/3*pow(mmst1,3)*zt2
          + 112/3*pow(mmst1,3)
          )

       + den(mmgl - mmst2,1)*sn2t * (
          - 4/3*mmsb1/mt*mgl*zt2
          - 8*mmsb1/mt*mgl
          - 4/3*mmsb2/mt*mgl*zt2
          - 8*mmsb2/mt*mgl
          - 4/3*mmst1/mt*mgl*zt2
          - 56/9*mmst1/mt*mgl
          - 88/9*mmst2/mt*mgl*zt2
          - 676/9*mmst2/mt*mgl
          - 32/3*mmsusy/mt*mgl*zt2
          - 64*mmsusy/mt*mgl
          )

       + den(mmgl - mmst2,1)*sn4t*cs2t * (
          - 8/9*mmst1/mt*mgl
          + 8/9*mmst2/mt*mgl
          )

       + den(mmgl - mmst2,1)*sqr(cs2t) * (
          - 20*mmst1*zt2
          - 428/3*mmst1
          - 139/9*mmst2*zt2
          - 821/9*mmst2
          )

       + den(mmgl - mmst2,1)*den(mmst1 - mmst2,1)*sqr(cs2t) * (
          + 20*sqr(mmst1)*zt2
          + 140*sqr(mmst1)
          )

       + den(mmgl - mmst2,1)*den(mmst1 - mmst2,1) * (
          - 20/3*sqr(mmst1)*zt2
          - 428/9*sqr(mmst1)
          )

       + den(mmgl - mmst2,1) * (
          - mmsb1*zt2
          - 6*mmsb1
          - mmsb2*zt2
          - 6*mmsb2
          + 17/3*mmst1*zt2
          + 398/9*mmst1
          + 110/9*mmst2*zt2
          + 163/3*mmst2
          - 8*mmsusy*zt2
          - 48*mmsusy
          )

       + den(mmgl - mmst2,2)*sqr(cs2t) * (
          - 16/9*mmst1*mmst2
          + 16/9*sqr(mmst2)
          )

       + den(mmgl - mmst2,2) * (
          - 4/3*mmsb1*mmst2*zt2
          - 8*mmsb1*mmst2
          - 4/3*mmsb2*mmst2*zt2
          - 8*mmsb2*mmst2
          - 4/3*mmst1*mmst2*zt2
          - 56/9*mmst1*mmst2
          - 32/3*mmst2*mmsusy*zt2
          - 64*mmst2*mmsusy
          + 44/3*sqr(mmst2)*zt2
          + 868/9*sqr(mmst2)
          )

       + den(mmgl - mmst2,3) * (
          + 16/3*pow(mmst2,3)*zt2
          + 112/3*pow(mmst2,3)
          )

       - 77/2
          + 8/9*zt2
         ;  

  return resmt;
}


template<class SoftPars>
double Softsusy<SoftPars>::twoLpMb() const {
  const double zt2 = sqr(PI) / 6.;
  double mmsb1 = sqr(displayDrBarPars().md(1, 3));
  double mmsb2 = sqr(displayDrBarPars().md(2, 3));
  double mmst1 = sqr(displayDrBarPars().mu(1, 3));
  double mmst2 = sqr(displayDrBarPars().mu(2, 3));
  double mgl = displayGaugino(3);
  double mmgl = sqr(mgl);
  double mt = displayDrBarPars().mt;
  double mmt = sqr(mt);
  double mb = displayDrBarPars().mb;
  double mmb = sqr(mb);
  double csb = cos(displayDrBarPars().thetab), 
    cs2b = cos(displayDrBarPars().thetab * 2.), 
    cs4b = cos(4 * displayDrBarPars().thetab);
  double snb = sin(displayDrBarPars().thetab), 
    sn2b = sin(displayDrBarPars().thetab * 2.), 
    sn4b = sin(4 * displayDrBarPars().thetab);
  double cst = cos(displayDrBarPars().thetat), 
    cs2t = cos(displayDrBarPars().thetat * 2.), 
    cs4t = cos(4 * displayDrBarPars().thetat);
  double snt = sin(displayDrBarPars().thetat), 
    sn2t = sin(displayDrBarPars().thetat * 2.), 
    sn4t = sin(4 * displayDrBarPars().thetat);
  double mmu = sqr(displayMu());

  /// average of first 2 generations squark mass
  double msq = 0.125 * (displayDrBarPars().mu(1, 1) + 
			displayDrBarPars().mu(2, 1) + 
			displayDrBarPars().md(1, 1) + 
			displayDrBarPars().md(2, 1) + 		       
			displayDrBarPars().mu(1, 2) + 
			displayDrBarPars().mu(2, 2) + 
			displayDrBarPars().md(1, 2) + 
			displayDrBarPars().md(2, 2));
  double mmsusy = sqr(msq);

  double lnMglSq = log(mmgl);
  double lnMsbSq = log(mmsb1);
  double lnMsb2Sq = log(mmsb2);
  double lnMst1Sq = log(mmst1);
  double lnMst2Sq = log(mmst2);
  double lnMmsusy = log(mmsusy);
  double lnMmt = log(mmt);
  double lnMmb = log(mmb);
  double lnMmu = log(mmu);

  cout << "# " << den(mmgl - mmsb1,1) << " " <<        + den(mmgl - mmsb2,3) * (
          + 16/3*pow(mmsb2,3)*zt2
          + 112/3*pow(mmsb2,3)
          )
       << " ";

  double resmb =

       + sqr(cs2b) * (
          - 640/9
          - 128/9*zt2
          )

       + fin(mmgl,mmsusy)*den(mmgl - mmsb1,1)*sn2b * (
          + 32/3/mb*mmsusy*mgl
          )

       + fin(mmgl,mmsusy)*den(mmgl - mmsb1,1) * (
          + 16/3*mmsb1
          - 8*mmsusy
          )

       + fin(mmgl,mmsusy)*den(mmgl - mmsb1,2)*sn2b * (
          + 32/3/mb*mmsb1*mmsusy*mgl
          - 32/3/mb*sqr(mmsb1)*mgl
          )

       + fin(mmgl,mmsusy)*den(mmgl - mmsb1,2) * (
          - 56/3*mmsb1*mmsusy
          + 56/3*sqr(mmsb1)
          )

       + fin(mmgl,mmsusy)*den(mmgl - mmsb1,3) * (
          - 32/3*sqr(mmsb1)*mmsusy
          + 32/3*pow(mmsb1,3)
          )

       + fin(mmgl,mmsusy)*den(mmgl - mmsb2,1)*sn2b * (
          - 32/3/mb*mmsusy*mgl
          )

       + fin(mmgl,mmsusy)*den(mmgl - mmsb2,1) * (
          + 16/3*mmsb2
          - 8*mmsusy
          )

       + fin(mmgl,mmsusy)*den(mmgl - mmsb2,2)*sn2b * (
          - 32/3/mb*mmsb2*mmsusy*mgl
          + 32/3/mb*sqr(mmsb2)*mgl
          )

       + fin(mmgl,mmsusy)*den(mmgl - mmsb2,2) * (
          - 56/3*mmsb2*mmsusy
          + 56/3*sqr(mmsb2)
          )

       + fin(mmgl,mmsusy)*den(mmgl - mmsb2,3) * (
          - 32/3*sqr(mmsb2)*mmsusy
          + 32/3*pow(mmsb2,3)
          )

       + fin(mmgl,mmsusy) * (
          - 16/3
          )

       + fin(mmsb1,mmgl)*den(mmgl - mmsb1,1)*sn2b * (
          + 88/9/mb*mmsb1*mgl
          )

       + fin(mmsb1,mmgl)*den(mmgl - mmsb1,1)*sqr(cs2b) * (
          + 5/3*mmsb1
          )

       + fin(mmsb1,mmgl)*den(mmgl - mmsb1,1)*den(mmsb1 - mmsb2,1)*sqr(cs2b) * (
          - 154/9*sqr(mmsb1)
          )

       + fin(mmsb1,mmgl)*den(mmgl - mmsb1,1)*den(mmsb1 - mmsb2,1) * (
          + 34/9*sqr(mmsb1)
          )

       + fin(mmsb1,mmgl)*den(mmgl - mmsb1,1) * (
          + 22/9*mmsb1
          )

       + fin(mmsb1,mmgl)*den(mmgl - mmsb1,2) * (
          + 12*sqr(mmsb1)
          )

       + fin(mmsb1,mmgl)*den(mmgl - mmsb1,3) * (
          + 16/3*pow(mmsb1,3)
          )

       + fin(mmsb1,mmgl)*den(mmgl - mmsb2,1)*sn2b * (
          - 16/9/mb*mmsb1*mgl
          )

       + fin(mmsb1,mmgl)*den(mmgl - mmsb2,1)*sqr(cs2b) * (
          - 5/3*mmsb1
          )

       + fin(mmsb1,mmgl)*den(mmgl - mmsb2,1)*den(mmsb1 - mmsb2,1)*sqr(cs2b) * (
          + 26/9*sqr(mmsb1)
          )

       + fin(mmsb1,mmgl)*den(mmgl - mmsb2,1)*den(mmsb1 - mmsb2,1) * (
          - 34/9*sqr(mmsb1)
          )

       + fin(mmsb1,mmgl)*den(mmgl - mmsb2,1) * (
          + 16/9*mmsb1
          )

       + fin(mmsb1,mmgl)*den(mmgl - mmsb2,2)*sn2b * (
          + 4/3/mb*mmsb1*mmsb2*mgl
          - 4/3/mb*sqr(mmsb1)*mgl
          )

       + fin(mmsb1,mmgl)*den(mmgl - mmsb2,2)*sqr(cs2b) * (
          + 26/9*mmsb1*mmsb2
          )

       + fin(mmsb1,mmgl)*den(mmgl - mmsb2,2) * (
          - 17/9*mmsb1*sqr(mmsb2
          - mmsb1)
          )

       + fin(mmsb1,mmgl)*den(mmgl - mmsb2,3) * (
          + 4/3*mmsb1*sqr(mmsb2)
          - 4/3*sqr(mmsb1)*mmsb2
          )

       + fin(mmsb1,mmgl)*den(mmsb1 - mmsb2,1)*sqr(cs2b) * (
          - 128/9*mmsb1
          )

       + fin(mmsb1,mmsb2)*den(mmgl - mmsb1,1)*sn2b * (
          - 4/9/mb*mmsb1*mgl
          )

       + fin(mmsb1,mmsb2)*den(mmgl - mmsb1,1)*sqr(cs2b) * (
          - 11/9*mmsb1
          )

       + fin(mmsb1,mmsb2)*den(mmgl - mmsb1,1)*den(mmsb1 - mmsb2,1) * (
          - 8/9*sqr(mmsb1)
          )

       + fin(mmsb1,mmsb2)*den(mmgl - mmsb1,1) * (
          + mmsb1
          )

       + fin(mmsb1,mmsb2)*den(mmgl - mmsb1,2)*sn2b * (
          - 4/3/mb*mmsb1*mmsb2*mgl
          + 4/3/mb*sqr(mmsb1)*mgl
          )

       + fin(mmsb1,mmsb2)*den(mmgl - mmsb1,2)*sqr(cs2b) * (
          - 26/9*sqr(mmsb1)
          )

       + fin(mmsb1,mmsb2)*den(mmgl - mmsb1,2) * (
          + mmsb1*mmsb2
          + 5/9*sqr(mmsb1)
          )

       + fin(mmsb1,mmsb2)*den(mmgl - mmsb1,3) * (
          + 4/3*sqr(mmsb1)*mmsb2
          - 4/3*pow(mmsb1,3)
          )

       + fin(mmsb1,mmsb2)*den(mmgl - mmsb2,1)*sn2b * (
          + 4/9/mb*mmsb1*mgl
          )

       + fin(mmsb1,mmsb2)*den(mmgl - mmsb2,1)*sqr(cs2b) * (
          - 11/9*mmsb1
          )

       + fin(mmsb1,mmsb2)*den(mmgl - mmsb2,1)*den(mmsb1 - mmsb2,1) * (
          + 8/9*sqr(mmsb1)
          )

       + fin(mmsb1,mmsb2)*den(mmgl - mmsb2,1) * (
          + 1/9*mmsb1
          )

       + fin(mmsb1,mmsb2)*den(mmgl - mmsb2,2)*sn2b * (
          - 4/3/mb*mmsb1*mmsb2*mgl
          + 4/3/mb*sqr(mmsb1)*mgl
          )

       + fin(mmsb1,mmsb2)*den(mmgl - mmsb2,2)*sqr(cs2b) * (
          - 26/9*mmsb1*mmsb2
          )

       + fin(mmsb1,mmsb2)*den(mmgl - mmsb2,2) * (
          + 5/9*mmsb1*sqr(mmsb2
          + mmsb1)
          )

       + fin(mmsb1,mmsb2)*den(mmgl - mmsb2,3) * (
          - 4/3*mmsb1*sqr(mmsb2)
          + 4/3*sqr(mmsb1)*mmsb2
          )

       + fin(mmsb1,mmst1)*den(mmgl - mmsb1,1) * (
          - 2/3*mmsb1
          )

       + fin(mmsb1,mmst1)*den(mmgl - mmsb1,2)*sn2b * (
          - 4/3/mb*mmsb1*mmst1*mgl
          + 4/3/mb*sqr(mmsb1)*mgl
          )

       + fin(mmsb1,mmst1)*den(mmgl - mmsb1,2) * (
          + mmsb1*mmst1
          - 7/3*sqr(mmsb1)
          )

       + fin(mmsb1,mmst1)*den(mmgl - mmsb1,3) * (
          + 4/3*sqr(mmsb1)*mmst1
          - 4/3*pow(mmsb1,3)
          )

       + fin(mmsb1,mmst2)*den(mmgl - mmsb1,1) * (
          - 2/3*mmsb1
          )

       + fin(mmsb1,mmst2)*den(mmgl - mmsb1,2)*sn2b * (
          - 4/3/mb*mmsb1*mmst2*mgl
          + 4/3/mb*sqr(mmsb1)*mgl
          )

       + fin(mmsb1,mmst2)*den(mmgl - mmsb1,2) * (
          + mmsb1*mmst2
          - 7/3*sqr(mmsb1)
          )

       + fin(mmsb1,mmst2)*den(mmgl - mmsb1,3) * (
          + 4/3*sqr(mmsb1)*mmst2
          - 4/3*pow(mmsb1,3)
          )

       + fin(mmsb1,mmsusy)*den(mmgl - mmsb1,1) * (
          - 16/3*mmsb1
          )

       + fin(mmsb1,mmsusy)*den(mmgl - mmsb1,2)*sn2b * (
          - 32/3/mb*mmsb1*mmsusy*mgl
          + 32/3/mb*sqr(mmsb1)*mgl
          )

       + fin(mmsb1,mmsusy)*den(mmgl - mmsb1,2) * (
          + 8*mmsb1*mmsusy
          - 56/3*sqr(mmsb1)
          )

       + fin(mmsb1,mmsusy)*den(mmgl - mmsb1,3) * (
          + 32/3*sqr(mmsb1)*mmsusy
          - 32/3*pow(mmsb1,3)
          )

       + fin(mmsb2,mmgl)*sqr(cs2b) * (
          - 128/9
          )

       + fin(mmsb2,mmgl)*den(mmgl - mmsb1,1)*sn2b * (
          + 16/9/mb*mmsb2*mgl
          )

       + fin(mmsb2,mmgl)*den(mmgl - mmsb1,1)*sqr(cs2b) * (
          + 26/9*mmsb1
          + 11/9*mmsb2
          )

       + fin(mmsb2,mmgl)*den(mmgl - mmsb1,1)*den(mmsb1 - mmsb2,1)*sqr(cs2b) * (
          - 26/9*sqr(mmsb1)
          )

       + fin(mmsb2,mmgl)*den(mmgl - mmsb1,1)*den(mmsb1 - mmsb2,1) * (
          + 34/9*sqr(mmsb1)
          )

       + fin(mmsb2,mmgl)*den(mmgl - mmsb1,1) * (
          - 34/9*mmsb1
          - 2*mmsb2
          )

       + fin(mmsb2,mmgl)*den(mmgl - mmsb1,2)*sn2b * (
          - 4/3/mb*mmsb1*mmsb2*mgl
          + 4/3/mb*sqr(mmsb2)*mgl
          )

       + fin(mmsb2,mmgl)*den(mmgl - mmsb1,2)*sqr(cs2b) * (
          + 26/9*mmsb1*mmsb2
          )

       + fin(mmsb2,mmgl)*den(mmgl - mmsb1,2) * (
          - 17/9*mmsb1*sqr(mmsb2
          - mmsb2)
          )

       + fin(mmsb2,mmgl)*den(mmgl - mmsb1,3) * (
          - 4/3*mmsb1*sqr(mmsb2)
          + 4/3*sqr(mmsb1)*mmsb2
          )

       + fin(mmsb2,mmgl)*den(mmgl - mmsb2,1)*sn2b * (
          - 88/9/mb*mmsb2*mgl
          )

       + fin(mmsb2,mmgl)*den(mmgl - mmsb2,1)*sqr(cs2b) * (
          - 154/9*mmsb1
          - 139/9*mmsb2
          )

       + fin(mmsb2,mmgl)*den(mmgl - mmsb2,1)*den(mmsb1 - mmsb2,1)*sqr(cs2b) * (
          + 154/9*sqr(mmsb1)
          )

       + fin(mmsb2,mmgl)*den(mmgl - mmsb2,1)*den(mmsb1 - mmsb2,1) * (
          - 34/9*sqr(mmsb1)
          )

       + fin(mmsb2,mmgl)*den(mmgl - mmsb2,1) * (
          + 34/9*mmsb1
          + 56/9*mmsb2
          )

       + fin(mmsb2,mmgl)*den(mmgl - mmsb2,2) * (
          + 12*sqr(mmsb2)
          )

       + fin(mmsb2,mmgl)*den(mmgl - mmsb2,3) * (
          + 16/3*pow(mmsb2,3)
          )

       + fin(mmsb2,mmgl)*den(mmsb1 - mmsb2,1)*sqr(cs2b) * (
          + 128/9*mmsb1
          )

       + fin(mmsb2,mmst1)*den(mmgl - mmsb2,1) * (
          - 2/3*mmsb2
          )

       + fin(mmsb2,mmst1)*den(mmgl - mmsb2,2)*sn2b * (
          + 4/3/mb*mmsb2*mmst1*mgl
          - 4/3/mb*sqr(mmsb2)*mgl
          )

       + fin(mmsb2,mmst1)*den(mmgl - mmsb2,2) * (
          + mmsb2*mmst1
          - 7/3*sqr(mmsb2)
          )

       + fin(mmsb2,mmst1)*den(mmgl - mmsb2,3) * (
          + 4/3*sqr(mmsb2)*mmst1
          - 4/3*pow(mmsb2,3)
          )

       + fin(mmsb2,mmst2)*den(mmgl - mmsb2,1) * (
          - 2/3*mmsb2
          )

       + fin(mmsb2,mmst2)*den(mmgl - mmsb2,2)*sn2b * (
          + 4/3/mb*mmsb2*mmst2*mgl
          - 4/3/mb*sqr(mmsb2)*mgl
          )

       + fin(mmsb2,mmst2)*den(mmgl - mmsb2,2) * (
          + mmsb2*mmst2
          - 7/3*sqr(mmsb2)
          )

       + fin(mmsb2,mmst2)*den(mmgl - mmsb2,3) * (
          + 4/3*sqr(mmsb2)*mmst2
          - 4/3*pow(mmsb2,3)
          )

       + fin(mmsb2,mmsusy)*den(mmgl - mmsb2,1) * (
          - 16/3*mmsb2
          )

       + fin(mmsb2,mmsusy)*den(mmgl - mmsb2,2)*sn2b * (
          + 32/3/mb*mmsb2*mmsusy*mgl
          - 32/3/mb*sqr(mmsb2)*mgl
          )

       + fin(mmsb2,mmsusy)*den(mmgl - mmsb2,2) * (
          + 8*mmsb2*mmsusy
          - 56/3*sqr(mmsb2)
          )

       + fin(mmsb2,mmsusy)*den(mmgl - mmsb2,3) * (
          + 32/3*sqr(mmsb2)*mmsusy
          - 32/3*pow(mmsb2,3)
          )

       + fin(mmst1,mmgl)*den(mmgl - mmsb1,1)*sn2b * (
          + 4/3/mb*mmst1*mgl
          )

       + fin(mmst1,mmgl)*den(mmgl - mmsb1,1) * (
          - 1/3*mmst1
          )

       + fin(mmst1,mmgl)*den(mmgl - mmsb1,2)*sn2b * (
          - 4/3/mb*mmsb1*mmst1*mgl
          + 4/3/mb*sqr(mmst1)*mgl
          )

       + fin(mmst1,mmgl)*den(mmgl - mmsb1,2) * (
          + mmsb1*sqr(mmst1
          - mmst1)
          )

       + fin(mmst1,mmgl)*den(mmgl - mmsb1,3) * (
          - 4/3*mmsb1*sqr(mmst1)
          + 4/3*sqr(mmsb1)*mmst1
          )

       + fin(mmst1,mmgl)*den(mmgl - mmsb2,1)*sn2b * (
          - 4/3/mb*mmst1*mgl
          )

       + fin(mmst1,mmgl)*den(mmgl - mmsb2,1) * (
          - 1/3*mmst1
          )

       + fin(mmst1,mmgl)*den(mmgl - mmsb2,2)*sn2b * (
          + 4/3/mb*mmsb2*mmst1*mgl
          - 4/3/mb*sqr(mmst1)*mgl
          )

       + fin(mmst1,mmgl)*den(mmgl - mmsb2,2) * (
          + mmsb2*sqr(mmst1
          - mmst1)
          )

       + fin(mmst1,mmgl)*den(mmgl - mmsb2,3) * (
          - 4/3*mmsb2*sqr(mmst1)
          + 4/3*sqr(mmsb2)*mmst1
          )

       + fin(mmst2,mmgl)*den(mmgl - mmsb1,1)*sn2b * (
          + 4/3/mb*mmst2*mgl
          )

       + fin(mmst2,mmgl)*den(mmgl - mmsb1,1) * (
          - 1/3*mmst2
          )

       + fin(mmst2,mmgl)*den(mmgl - mmsb1,2)*sn2b * (
          - 4/3/mb*mmsb1*mmst2*mgl
          + 4/3/mb*sqr(mmst2)*mgl
          )

       + fin(mmst2,mmgl)*den(mmgl - mmsb1,2) * (
          + mmsb1*sqr(mmst2
          - mmst2)
          )

       + fin(mmst2,mmgl)*den(mmgl - mmsb1,3) * (
          - 4/3*mmsb1*sqr(mmst2)
          + 4/3*sqr(mmsb1)*mmst2
          )

       + fin(mmst2,mmgl)*den(mmgl - mmsb2,1)*sn2b * (
          - 4/3/mb*mmst2*mgl
          )

       + fin(mmst2,mmgl)*den(mmgl - mmsb2,1) * (
          - 1/3*mmst2
          )

       + fin(mmst2,mmgl)*den(mmgl - mmsb2,2)*sn2b * (
          + 4/3/mb*mmsb2*mmst2*mgl
          - 4/3/mb*sqr(mmst2)*mgl
          )

       + fin(mmst2,mmgl)*den(mmgl - mmsb2,2) * (
          + mmsb2*sqr(mmst2
          - mmst2)
          )

       + fin(mmst2,mmgl)*den(mmgl - mmsb2,3) * (
          - 4/3*mmsb2*sqr(mmst2)
          + 4/3*sqr(mmsb2)*mmst2
          )

       + lnMmb*lnMglSq*den(mmgl - mmsb1,1)*sn2b * (
          + 16/3/mb*mmsb1*mgl
          )

       + lnMmb*lnMglSq*den(mmgl - mmsb1,1) * (
          - 16/3*mmsb1
          )

       + lnMmb*lnMglSq*den(mmgl - mmsb1,2) * (
          - 8/3*sqr(mmsb1)
          )

       + lnMmb*lnMglSq*den(mmgl - mmsb2,1)*sn2b * (
          - 16/3/mb*mmsb2*mgl
          )

       + lnMmb*lnMglSq*den(mmgl - mmsb2,1) * (
          - 16/3*mmsb2
          )

       + lnMmb*lnMglSq*den(mmgl - mmsb2,2) * (
          - 8/3*sqr(mmsb2)
          )

       + lnMmb*lnMglSq * (
          - 40/3
          )

       + lnMmb*lnMsbSq*den(mmgl - mmsb1,1)*sn2b * (
          - 16/3/mb*mmsb1*mgl
          )

       + lnMmb*lnMsbSq*den(mmgl - mmsb1,1) * (
          + 16/3*mmsb1
          )

       + lnMmb*lnMsbSq*den(mmgl - mmsb1,2) * (
          + 8/3*sqr(mmsb1)
          )

       + lnMmb*lnMsbSq * (
          - 2/3
          )

       + lnMmb*lnMsb2Sq*den(mmgl - mmsb2,1)*sn2b * (
          + 16/3/mb*mmsb2*mgl
          )

       + lnMmb*lnMsb2Sq*den(mmgl - mmsb2,1) * (
          + 16/3*mmsb2
          )

       + lnMmb*lnMsb2Sq*den(mmgl - mmsb2,2) * (
          + 8/3*sqr(mmsb2)
          )

       + lnMmb*lnMsb2Sq * (
          - 2/3
          )

       + lnMmb*lnMst1Sq * (
          - 2/3
          )

       + lnMmb*lnMst2Sq * (
          - 2/3
          )

       + lnMmb*lnMmsusy * (
          - 16/3
          )

       + lnMmb*lnMmu * (
          + 64/3
          )

       + lnMmb*den(mmgl - mmsb1,1) * (
          + 8/3*mmsb1
          )

       + lnMmb*den(mmgl - mmsb2,1) * (
          + 8/3*mmsb2
          )

       + lnMmb * (
          + 8
          )

       + lnMglSq*sqr(cs2b) * (
          + 128/3
          )

       + sqr(lnMglSq)*sqr(cs2b) * (
          - 64/9
          )

       + sqr(lnMglSq)*den(mmgl - mmsb1,1)*sn2b * (
          + 10/mb*mmsb1*mgl
          + 2/3/mb*mmsb2*mgl
          + 2/3/mb*mmst1*mgl
          + 2/3/mb*mmst2*mgl
          + 16/3/mb*mmsusy*mgl
          )

       + sqr(lnMglSq)*den(mmgl - mmsb1,1)*sqr(cs2b) * (
          + 53/3*mmsb1
          )

       + sqr(lnMglSq)*den(mmgl - mmsb1,1)*den(mmsb1 - mmsb2,1)*sn2b * (
          + 8/3/mb*sqr(mmsb1)*mgl
          )

       + sqr(lnMglSq)*den(mmgl - mmsb1,1)*den(mmsb1 - mmsb2,1)*sqr(cs2b) * (
          - 223/9*sqr(mmsb1)
          )

       + sqr(lnMglSq)*den(mmgl - mmsb1,1)*den(mmsb1 - mmsb2,1) * (
          + 43/9*sqr(mmsb1)
          )

       + sqr(lnMglSq)*den(mmgl - mmsb1,1)*den(mmsb1 - mmsb2,2)*sn2b * (
          - 16/9/mb*pow(mmsb1,3)*mgl
          )

       + sqr(lnMglSq)*den(mmgl - mmsb1,1)*den(mmsb1 - mmsb2,2) * (
          + 32/9*pow(mmsb1,3)
          )

       + sqr(lnMglSq)*den(mmgl - mmsb1,1)*den(mmsb1 - mmsb2,3) * (
          - 16/9*pow(mmsb1,4)
          )

       + sqr(lnMglSq)*den(mmgl - mmsb1,1) * (
          - 73/6*mmsb1
          - 1/2*mmsb2
          - 1/2*mmst1
          - 1/2*mmst2
          + 4/3*mmsusy
          )

       + sqr(lnMglSq)*den(mmgl - mmsb1,2)*sn2b * (
          + 2/3/mb*mmsb1*mmsb2*mgl
          + 2/3/mb*mmsb1*mmst1*mgl
          + 2/3/mb*mmsb1*mmst2*mgl
          - 16/mb*mmsb1*mmsusy*mgl
          - 122/9/mb*sqr(mmsb1)*mgl
          + 32/3/mb*sqr(mmsusy)*mgl
          )

       + sqr(lnMglSq)*den(mmgl - mmsb1,2)*sqr(cs2b) * (
          + 53/6*sqr(mmsb1)
          )

       + sqr(lnMglSq)*den(mmgl - mmsb1,2)*den(mmsb1 - mmsb2,1)*sn2b * (
          + 8/9/mb*pow(mmsb1,3)*mgl
          )

       + sqr(lnMglSq)*den(mmgl - mmsb1,2)*den(mmsb1 - mmsb2,2) * (
          + 8/9*pow(mmsb1,4)
          )

       + sqr(lnMglSq)*den(mmgl - mmsb1,2) * (
          - 7/6*mmsb1*mmsb2
          - 7/6*mmsb1*mmst1
          - 7/6*mmsb1*mmst2
          + 52/3*mmsb1*mmsusy
          + 151/9*sqr(mmsb1)
          - 8*sqr(mmsusy)
          )

       + sqr(lnMglSq)*den(mmgl - mmsb1,3)*sn2b * (
          - 8/9/mb*pow(mmsb1,3)*mgl
          )

       + sqr(lnMglSq)*den(mmgl - mmsb1,3) * (
          - 32/3*mmsb1*sqr(mmsusy)
          - 2/3*sqr(mmsb1)*mmsb2
          - 2/3*sqr(mmsb1)*mmst1
          - 2/3*sqr(mmsb1)*mmst2
          + 16*sqr(mmsb1)*mmsusy
          + 46/3*pow(mmsb1,3)
          )

       + sqr(lnMglSq)*den(mmgl - mmsb1,4) * (
          + 4/9*pow(mmsb1,4)
          )

       + sqr(lnMglSq)*den(mmgl - mmsb2,1)*sn2b * (
          + 2/9/mb*mmsb1*mgl
          - 98/9/mb*mmsb2*mgl
          - 2/3/mb*mmst1*mgl
          - 2/3/mb*mmst2*mgl
          - 16/3/mb*mmsusy*mgl
          )

       + sqr(lnMglSq)*den(mmgl - mmsb2,1)*sqr(cs2b) * (
          - 223/9*mmsb1
          - 64/9*mmsb2
          )

       + sqr(lnMglSq)*den(mmgl - mmsb2,1)*den(mmsb1 - mmsb2,1)*sn2b * (
          - 8/3/mb*sqr(mmsb1)*mgl
          )

       + sqr(lnMglSq)*den(mmgl - mmsb2,1)*den(mmsb1 - mmsb2,1)*sqr(cs2b) * (
          + 223/9*sqr(mmsb1)
          )

       + sqr(lnMglSq)*den(mmgl - mmsb2,1)*den(mmsb1 - mmsb2,1) * (
          - 43/9*sqr(mmsb1)
          )

       + sqr(lnMglSq)*den(mmgl - mmsb2,1)*den(mmsb1 - mmsb2,2)*sn2b * (
          + 16/9/mb*pow(mmsb1,3)*mgl
          )

       + sqr(lnMglSq)*den(mmgl - mmsb2,1)*den(mmsb1 - mmsb2,2) * (
          - 32/9*pow(mmsb1,3)
          )

       + sqr(lnMglSq)*den(mmgl - mmsb2,1)*den(mmsb1 - mmsb2,3) * (
          + 16/9*pow(mmsb1,4)
          )

       + sqr(lnMglSq)*den(mmgl - mmsb2,1) * (
          + 109/18*mmsb1
          - 101/18*mmsb2
          - 1/2*mmst1
          - 1/2*mmst2
          + 4/3*mmsusy
          )

       + sqr(lnMglSq)*den(mmgl - mmsb2,2)*sn2b * (
          - 14/9/mb*mmsb1*mmsb2*mgl
          - 8/9/mb*sqr(mmsb1)*mgl
          - 2/3/mb*mmsb2*mmst1*mgl
          - 2/3/mb*mmsb2*mmst2*mgl
          + 16/mb*mmsb2*mmsusy*mgl
          + 38/3/mb*sqr(mmsb2)*mgl
          - 32/3/mb*sqr(mmsusy)*mgl
          )

       + sqr(lnMglSq)*den(mmgl - mmsb2,2)*sqr(cs2b) * (
          + 53/6*sqr(mmsb2)
          )

       + sqr(lnMglSq)*den(mmgl - mmsb2,2)*den(mmsb1 - mmsb2,1)*sn2b * (
          + 8/9/mb*pow(mmsb1,3)*mgl
          )

       + sqr(lnMglSq)*den(mmgl - mmsb2,2)*den(mmsb1 - mmsb2,1) * (
          - 32/9*pow(mmsb1,3)
          )

       + sqr(lnMglSq)*den(mmgl - mmsb2,2)*den(mmsb1 - mmsb2,2) * (
          + 8/9*pow(mmsb1,4)
          )

       + sqr(lnMglSq)*den(mmgl - mmsb2,2) * (
          + 11/18*mmsb1*mmsb2
          + 8/3*sqr(mmsb1)
          - 7/6*mmsb2*mmst1
          - 7/6*mmsb2*mmst2
          + 52/3*mmsb2*mmsusy
          + 53/3*sqr(mmsb2)
          - 8*sqr(mmsusy)
          )

       + sqr(lnMglSq)*den(mmgl - mmsb2,3)*sn2b * (
          + 8/9/mb*pow(mmsb2,3)*mgl
          )

       + sqr(lnMglSq)*den(mmgl - mmsb2,3) * (
          - 2/3*mmsb1*sqr(mmsb2)
          - 32/3*mmsb2*sqr(mmsusy)
          - 2/3*sqr(mmsb2)*mmst1
          - 2/3*sqr(mmsb2)*mmst2
          + 16*sqr(mmsb2)*mmsusy
          + 46/3*pow(mmsb2,3)
          )

       + sqr(lnMglSq)*den(mmgl - mmsb2,4) * (
          + 4/9*pow(mmsb2,4)
          )

       + sqr(lnMglSq) * (
          - 166/9
          )

       + lnMglSq*lnMsbSq*sn2b * (
          - 208/9/mb*mgl
          )

       + lnMglSq*lnMsbSq*sqr(cs2b) * (
          - 64/9
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmsb1,1)*sn2b * (
          - 28/3/mb*mmsb1*mgl
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmsb1,1)*sn4b*cs2b * (
          - 8/9/mb*mmsb1*mgl
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmsb1,1)*sqr(cs2b) * (
          - 121/9*mmsb1
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmsb1,1)*den(mmsb1 - mmsb2,1)*sn2b
       * (
          - 8/3/mb*sqr(mmsb1)*mgl
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmsb1,1)*den(mmsb1 - mmsb2,1)*sn4b*
      cs2b * (
          + 16/9/mb*sqr(mmsb1)*mgl
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmsb1,1)*den(mmsb1 - mmsb2,1)*sqr(cs2b)
       * (
          + 287/9*sqr(mmsb1)
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmsb1,1)*den(mmsb1 - mmsb2,1) * (
          - 43/9*sqr(mmsb1)
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmsb1,1)*den(mmsb1 - mmsb2,2)*sn2b
       * (
          + 16/9/mb*pow(mmsb1,3)*mgl
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmsb1,1)*den(mmsb1 - mmsb2,2) * (
          - 32/9*pow(mmsb1,3)
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmsb1,1)*den(mmsb1 - mmsb2,3) * (
          + 16/9*pow(mmsb1,4)
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmsb1,1) * (
          - 34/3*mmsb1
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmsb1,2)*sn2b * (
          + 172/9/mb*sqr(mmsb1)*mgl
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmsb1,2)*sn4b*cs2b * (
          - 8/9/mb*sqr(mmsb1)*mgl
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmsb1,2)*sqr(cs2b) * (
          - 11/9*sqr(mmsb1)
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmsb1,2)*den(mmsb1 - mmsb2,1)*sn2b
       * (
          - 8/9/mb*pow(mmsb1,3)*mgl
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmsb1,2)*den(mmsb1 - mmsb2,2) * (
          - 8/9*pow(mmsb1,4)
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmsb1,2) * (
          - 130/3*sqr(mmsb1)
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmsb1,3)*sn2b * (
          + 16/9/mb*pow(mmsb1,3)*mgl
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmsb1,3)*sqr(cs2b) * (
          + 16/9*pow(mmsb1,3)
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmsb1,3) * (
          - 212/9*pow(mmsb1,3)
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmsb1,4) * (
          - 8/9*pow(mmsb1,4)
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmsb2,1)*sn2b * (
          + 20/9/mb*mmsb1*mgl
          + 8/9/mb*mmsb2*mgl
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmsb2,1)*sn4b*cs2b * (
          + 8/9/mb*mmsb1*mgl
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmsb2,1)*sqr(cs2b) * (
          + 5/3*mmsb1
          - 22/9*mmsb2
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmsb2,1)*den(mmsb1 - mmsb2,1)*sn2b
       * (
          + 8/3/mb*sqr(mmsb1)*mgl
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmsb2,1)*den(mmsb1 - mmsb2,1)*sn4b*
      cs2b * (
          - 16/9/mb*sqr(mmsb1)*mgl
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmsb2,1)*den(mmsb1 - mmsb2,1)*sqr(cs2b)
       * (
          - 31/9*sqr(mmsb1)
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmsb2,1)*den(mmsb1 - mmsb2,1) * (
          + 43/9*sqr(mmsb1)
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmsb2,1)*den(mmsb1 - mmsb2,2)*sn2b
       * (
          - 16/9/mb*pow(mmsb1,3)*mgl
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmsb2,1)*den(mmsb1 - mmsb2,2) * (
          + 32/9*pow(mmsb1,3)
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmsb2,1)*den(mmsb1 - mmsb2,3) * (
          - 16/9*pow(mmsb1,4)
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmsb2,1) * (
          - 34/9*mmsb1
          + 2/9*mmsb2
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmsb2,2)*sn2b * (
          + 4/mb*mmsb1*mmsb2*mgl
          + 8/9/mb*sqr(mmsb1)*mgl
          - 8/3/mb*sqr(mmsb2)*mgl
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmsb2,2)*sn4b*cs2b * (
          - 8/9/mb*mmsb1*mmsb2*mgl
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmsb2,2)*sqr(cs2b) * (
          - 32/9*mmsb1*mmsb2
          - 52/9*sqr(mmsb2)
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmsb2,2)*den(mmsb1 - mmsb2,1)*sn2b
       * (
          - 8/9/mb*pow(mmsb1,3)*mgl
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmsb2,2)*den(mmsb1 - mmsb2,1) * (
          + 32/9*pow(mmsb1,3)
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmsb2,2)*den(mmsb1 - mmsb2,2) * (
          - 8/9*pow(mmsb1,4)
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmsb2,2) * (
          + 37/9*mmsb1*mmsb2
          - 8/3*sqr(mmsb1)
          + 10/9*sqr(mmsb2)
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmsb2,3)*sqr(cs2b) * (
          - 16/9*mmsb1*sqr(mmsb2)
          )

       + lnMglSq*lnMsbSq*den(mmgl - mmsb2,3) * (
          + 28/9*mmsb1*sqr(mmsb2)
          - 8/3*pow(mmsb2,3)
          )

       + lnMglSq*lnMsbSq*den(mmsb1 - mmsb2,1)*sqr(cs2b) * (
          + 256/9*mmgl
          + 256/9*mmsb1
          )

       + lnMglSq*lnMsbSq * (
          + 52/9
          )

       + lnMglSq*lnMsb2Sq*sn2b * (
          + 208/9/mb*mgl
          )

       + lnMglSq*lnMsb2Sq*sqr(cs2b) * (
          + 64/3
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmsb1,1)*sn2b * (
          - 28/9/mb*mmsb2*mgl
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmsb1,1)*sn4b*cs2b * (
          + 16/9/mb*mmsb1*mgl
          + 8/9/mb*mmsb2*mgl
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmsb1,1)*sqr(cs2b) * (
          - 53/9*mmsb1
          - 16/9*mmsb2
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmsb1,1)*den(mmsb1 - mmsb2,1)*sn2b
       * (
          - 8/3/mb*sqr(mmsb1)*mgl
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmsb1,1)*den(mmsb1 - mmsb2,1)*sn4b*
      cs2b * (
          - 16/9/mb*sqr(mmsb1)*mgl
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmsb1,1)*den(mmsb1 - mmsb2,1)*sqr(cs2b)
       * (
          + 31/9*sqr(mmsb1)
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmsb1,1)*den(mmsb1 - mmsb2,1) * (
          - 43/9*sqr(mmsb1)
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmsb1,1)*den(mmsb1 - mmsb2,2)*sn2b
       * (
          + 16/9/mb*pow(mmsb1,3)*mgl
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmsb1,1)*den(mmsb1 - mmsb2,2) * (
          - 32/9*pow(mmsb1,3)
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmsb1,1)*den(mmsb1 - mmsb2,3) * (
          + 16/9*pow(mmsb1,4)
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmsb1,1) * (
          + 61/9*mmsb1
          + 25/9*mmsb2
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmsb1,2)*sn2b * (
          - 28/9/mb*mmsb1*mmsb2*mgl
          + 32/9/mb*sqr(mmsb1)*mgl
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmsb1,2)*sn4b*cs2b * (
          + 8/9/mb*mmsb1*mmsb2*mgl
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmsb1,2)*sqr(cs2b) * (
          - 32/9*mmsb1*mmsb2
          - 52/9*sqr(mmsb1)
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmsb1,2)*den(mmsb1 - mmsb2,1)*sn2b
       * (
          - 8/9/mb*pow(mmsb1,3)*mgl
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmsb1,2)*den(mmsb1 - mmsb2,2) * (
          - 8/9*pow(mmsb1,4)
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmsb1,2) * (
          + 53/9*mmsb1*mmsb2
          + 2*sqr(mmsb1)
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmsb1,3)*sqr(cs2b) * (
          - 16/9*sqr(mmsb1)*mmsb2
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmsb1,3) * (
          + 28/9*sqr(mmsb1)*mmsb2
          - 8/3*pow(mmsb1,3)
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmsb2,1)*sn2b * (
          - 8/9/mb*mmsb1*mgl
          + 92/9/mb*mmsb2*mgl
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmsb2,1)*sn4b*cs2b * (
          - 16/9/mb*mmsb1*mgl
          - 8/9/mb*mmsb2*mgl
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmsb2,1)*sqr(cs2b) * (
          + 287/9*mmsb1
          + 166/9*mmsb2
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmsb2,1)*den(mmsb1 - mmsb2,1)*sn2b
       * (
          + 8/3/mb*sqr(mmsb1)*mgl
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmsb2,1)*den(mmsb1 - mmsb2,1)*sn4b*
      cs2b * (
          + 16/9/mb*sqr(mmsb1)*mgl
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmsb2,1)*den(mmsb1 - mmsb2,1)*sqr(cs2b)
       * (
          - 287/9*sqr(mmsb1)
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmsb2,1)*den(mmsb1 - mmsb2,1) * (
          + 43/9*sqr(mmsb1)
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmsb2,1)*den(mmsb1 - mmsb2,2)*sn2b
       * (
          - 16/9/mb*pow(mmsb1,3)*mgl
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmsb2,1)*den(mmsb1 - mmsb2,2) * (
          + 32/9*pow(mmsb1,3)
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmsb2,1)*den(mmsb1 - mmsb2,3) * (
          - 16/9*pow(mmsb1,4)
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmsb2,1) * (
          - 59/9*mmsb1
          - 161/9*mmsb2
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmsb2,2)*sn2b * (
          + 8/9/mb*mmsb1*mmsb2*mgl
          + 8/9/mb*sqr(mmsb1)*mgl
          - 164/9/mb*sqr(mmsb2)*mgl
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmsb2,2)*sn4b*cs2b * (
          + 8/9/mb*sqr(mmsb2)*mgl
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmsb2,2)*sqr(cs2b) * (
          - 11/9*sqr(mmsb2)
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmsb2,2)*den(mmsb1 - mmsb2,1)*sn2b
       * (
          - 8/9/mb*pow(mmsb1,3)*mgl
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmsb2,2)*den(mmsb1 - mmsb2,1) * (
          + 32/9*pow(mmsb1,3)
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmsb2,2)*den(mmsb1 - mmsb2,2) * (
          - 8/9*pow(mmsb1,4)
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmsb2,2) * (
          - 16/9*mmsb1*mmsb2
          - 8/3*sqr(mmsb1)
          - 398/9*sqr(mmsb2)
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmsb2,3)*sn2b * (
          - 16/9/mb*pow(mmsb2,3)*mgl
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmsb2,3)*sqr(cs2b) * (
          + 16/9*pow(mmsb2,3)
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmsb2,3) * (
          - 212/9*pow(mmsb2,3)
          )

       + lnMglSq*lnMsb2Sq*den(mmgl - mmsb2,4) * (
          - 8/9*pow(mmsb2,4)
          )

       + lnMglSq*lnMsb2Sq*den(mmsb1 - mmsb2,1)*sqr(cs2b) * (
          - 256/9*mmgl
          - 256/9*mmsb1
          )

       + lnMglSq*lnMsb2Sq * (
          + 52/9
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmsb1,1)*sn2b * (
          - 4/3/mb*mmst1*mgl
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmsb1,1) * (
          - 4/3*mmsb1
          + mmst1
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmsb1,2)*sn2b * (
          - 4/3/mb*mmsb1*mmst1*mgl
          + 8/3/mb*sqr(mmsb1)*mgl
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmsb1,2) * (
          + 7/3*mmsb1*mmst1
          - 14/3*sqr(mmsb1)
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmsb1,3) * (
          + 4/3*sqr(mmsb1)*mmst1
          - 8/3*pow(mmsb1,3)
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmsb2,1)*sn2b * (
          + 4/3/mb*mmst1*mgl
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmsb2,1) * (
          - 4/3*mmsb2
          + mmst1
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmsb2,2)*sn2b * (
          + 4/3/mb*mmsb2*mmst1*mgl
          - 8/3/mb*sqr(mmsb2)*mgl
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmsb2,2) * (
          + 7/3*mmsb2*mmst1
          - 14/3*sqr(mmsb2)
          )

       + lnMglSq*lnMst1Sq*den(mmgl - mmsb2,3) * (
          + 4/3*sqr(mmsb2)*mmst1
          - 8/3*pow(mmsb2,3)
          )

       + lnMglSq*lnMst1Sq * (
          + 4/3
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmsb1,1)*sn2b * (
          - 4/3/mb*mmst2*mgl
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmsb1,1) * (
          - 4/3*mmsb1
          + mmst2
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmsb1,2)*sn2b * (
          - 4/3/mb*mmsb1*mmst2*mgl
          + 8/3/mb*sqr(mmsb1)*mgl
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmsb1,2) * (
          + 7/3*mmsb1*mmst2
          - 14/3*sqr(mmsb1)
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmsb1,3) * (
          + 4/3*sqr(mmsb1)*mmst2
          - 8/3*pow(mmsb1,3)
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmsb2,1)*sn2b * (
          + 4/3/mb*mmst2*mgl
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmsb2,1) * (
          - 4/3*mmsb2
          + mmst2
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmsb2,2)*sn2b * (
          + 4/3/mb*mmsb2*mmst2*mgl
          - 8/3/mb*sqr(mmsb2)*mgl
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmsb2,2) * (
          + 7/3*mmsb2*mmst2
          - 14/3*sqr(mmsb2)
          )

       + lnMglSq*lnMst2Sq*den(mmgl - mmsb2,3) * (
          + 4/3*sqr(mmsb2)*mmst2
          - 8/3*pow(mmsb2,3)
          )

       + lnMglSq*lnMst2Sq * (
          + 4/3
          )

       + lnMglSq*lnMmsusy*den(mmgl - mmsb1,1)*sn2b * (
          - 32/3/mb*mmsusy*mgl
          )

       + lnMglSq*lnMmsusy*den(mmgl - mmsb1,1) * (
          - 8/3*mmsusy
          )

       + lnMglSq*lnMmsusy*den(mmgl - mmsb1,2)*sn2b * (
          + 32/mb*mmsb1*mmsusy*mgl
          - 64/3/mb*sqr(mmsusy)*mgl
          )

       + lnMglSq*lnMmsusy*den(mmgl - mmsb1,2) * (
          - 104/3*mmsb1*mmsusy
          + 16*sqr(mmsusy)
          )

       + lnMglSq*lnMmsusy*den(mmgl - mmsb1,3) * (
          + 64/3*mmsb1*sqr(mmsusy)
          - 32*sqr(mmsb1)*mmsusy
          )

       + lnMglSq*lnMmsusy*den(mmgl - mmsb2,1)*sn2b * (
          + 32/3/mb*mmsusy*mgl
          )

       + lnMglSq*lnMmsusy*den(mmgl - mmsb2,1) * (
          - 8/3*mmsusy
          )

       + lnMglSq*lnMmsusy*den(mmgl - mmsb2,2)*sn2b * (
          - 32/mb*mmsb2*mmsusy*mgl
          + 64/3/mb*sqr(mmsusy)*mgl
          )

       + lnMglSq*lnMmsusy*den(mmgl - mmsb2,2) * (
          - 104/3*mmsb2*mmsusy
          + 16*sqr(mmsusy)
          )

       + lnMglSq*lnMmsusy*den(mmgl - mmsb2,3) * (
          + 64/3*mmsb2*sqr(mmsusy)
          - 32*sqr(mmsb2)*mmsusy
          )

       + lnMglSq*lnMmu*den(mmgl - mmsb1,1)*sn2b * (
          - 16/mb*mmsb1*mgl
          + 16/9/mb*mmsb2*mgl
          )

       + lnMglSq*lnMmu*den(mmgl - mmsb1,1)*sn4b*cs2b * (
          - 8/9/mb*mmsb1*mgl
          - 8/9/mb*mmsb2*mgl
          )

       + lnMglSq*lnMmu*den(mmgl - mmsb1,1)*sqr(cs2b) * (
          - 16*mmsb1
          + 16/9*mmsb2
          )

       + lnMglSq*lnMmu*den(mmgl - mmsb1,1)*den(mmsb1 - mmsb2,1)*sqr(cs2b)
       * (
          + 128/9*sqr(mmsb1)
          )

       + lnMglSq*lnMmu*den(mmgl - mmsb1,1) * (
          + 332/9*mmsb1
          - 16/9*mmsb2
          )

       + lnMglSq*lnMmu*den(mmgl - mmsb1,2)*sn2b * (
          + 16/9/mb*mmsb1*mmsb2*mgl
          - 8/9/mb*sqr(mmsb1)*mgl
          )

       + lnMglSq*lnMmu*den(mmgl - mmsb1,2)*sn4b*cs2b * (
          - 8/9/mb*mmsb1*mmsb2*mgl
          + 8/9/mb*sqr(mmsb1)*mgl
          )

       + lnMglSq*lnMmu*den(mmgl - mmsb1,2)*sqr(cs2b) * (
          + 32/9*mmsb1*mmsb2
          - 32/3*sqr(mmsb1)
          )

       + lnMglSq*lnMmu*den(mmgl - mmsb1,2) * (
          - 32/9*mmsb1*mmsb2
          + 178/9*sqr(mmsb1)
          )

       + lnMglSq*lnMmu*den(mmgl - mmsb1,3)*sqr(cs2b) * (
          + 16/9*sqr(mmsb1)*mmsb2
          - 16/9*pow(mmsb1,3)
          )

       + lnMglSq*lnMmu*den(mmgl - mmsb1,3) * (
          - 16/9*sqr(mmsb1)*mmsb2
          + 8/9*pow(mmsb1,3)
          )

       + lnMglSq*lnMmu*den(mmgl - mmsb2,1)*sn2b * (
          - 16/9/mb*mmsb1*mgl
          + 16/mb*mmsb2*mgl
          )

       + lnMglSq*lnMmu*den(mmgl - mmsb2,1)*sn4b*cs2b * (
          + 8/9/mb*mmsb1*mgl
          + 8/9/mb*mmsb2*mgl
          )

       + lnMglSq*lnMmu*den(mmgl - mmsb2,1)*sqr(cs2b) * (
          + 16*mmsb1
          - 16/9*mmsb2
          )

       + lnMglSq*lnMmu*den(mmgl - mmsb2,1)*den(mmsb1 - mmsb2,1)*sqr(cs2b)
       * (
          - 128/9*sqr(mmsb1)
          )

       + lnMglSq*lnMmu*den(mmgl - mmsb2,1) * (
          - 16/9*mmsb1
          + 332/9*mmsb2
          )

       + lnMglSq*lnMmu*den(mmgl - mmsb2,2)*sn2b * (
          - 16/9/mb*mmsb1*mmsb2*mgl
          + 8/9/mb*sqr(mmsb2)*mgl
          )

       + lnMglSq*lnMmu*den(mmgl - mmsb2,2)*sn4b*cs2b * (
          + 8/9/mb*mmsb1*mmsb2*mgl
          - 8/9/mb*sqr(mmsb2)*mgl
          )

       + lnMglSq*lnMmu*den(mmgl - mmsb2,2)*sqr(cs2b) * (
          + 32/9*mmsb1*mmsb2
          - 32/3*sqr(mmsb2)
          )

       + lnMglSq*lnMmu*den(mmgl - mmsb2,2) * (
          - 32/9*mmsb1*mmsb2
          + 178/9*sqr(mmsb2)
          )

       + lnMglSq*lnMmu*den(mmgl - mmsb2,3)*sqr(cs2b) * (
          + 16/9*mmsb1*sqr(mmsb2)
          - 16/9*pow(mmsb2,3)
          )

       + lnMglSq*lnMmu*den(mmgl - mmsb2,3) * (
          - 16/9*mmsb1*sqr(mmsb2)
          + 8/9*pow(mmsb2,3)
          )

       + lnMglSq*lnMmu * (
          + 36
          )

       + lnMglSq*den(mmgl - mmsb1,1)*sn2b * (
          - 176/3/mb*mmsb1*mgl
          - 8/9/mb*mmsb2*mgl
          - 8/3/mb*mmst1*mgl
          - 8/3/mb*mmst2*mgl
          + 128/3/mb*mmsusy*mgl
          )

       + lnMglSq*den(mmgl - mmsb1,1)*sn4b*cs2b * (
          - 8/9/mb*mmsb1*mgl
          - 8/9/mb*mmsb2*mgl
          )

       + lnMglSq*den(mmgl - mmsb1,1)*sqr(cs2b) * (
          - 428/9*mmsb1
          + 16/9*mmsb2
          )

       + lnMglSq*den(mmgl - mmsb1,1)*den(mmsb1 - mmsb2,1)*sqr(cs2b) * (
          + 796/9*sqr(mmsb1)
          )

       + lnMglSq*den(mmgl - mmsb1,1)*den(mmsb1 - mmsb2,1) * (
          - 76/3*sqr(mmsb1)
          )

       + lnMglSq*den(mmgl - mmsb1,1) * (
          + 290/9*mmsb1
          + 2/9*mmsb2
          + 2*mmst1
          + 2*mmst2
          - 16*mmsusy
          )

       + lnMglSq*den(mmgl - mmsb1,2)*sn2b * (
          - 8/9/mb*mmsb1*mmsb2*mgl
          - 8/3/mb*mmsb1*mmst1*mgl
          - 8/3/mb*mmsb1*mmst2*mgl
          - 64/3/mb*mmsb1*mmsusy*mgl
          + 56/9/mb*sqr(mmsb1)*mgl
          + 32/mb*sqr(mmsusy)*mgl
          )

       + lnMglSq*den(mmgl - mmsb1,2)*sn4b*cs2b * (
          - 8/9/mb*mmsb1*mmsb2*mgl
          + 8/9/mb*sqr(mmsb1)*mgl
          )

       + lnMglSq*den(mmgl - mmsb1,2)*sqr(cs2b) * (
          + 32/9*mmsb1*mmsb2
          - 238/9*sqr(mmsb1)
          )

       + lnMglSq*den(mmgl - mmsb1,2)*den(mmsb1 - mmsb2,1) * (
          - 8/9*pow(mmsb1,3)
          )

       + lnMglSq*den(mmgl - mmsb1,2) * (
          + 10/9*mmsb1*mmsb2
          + 14/3*mmsb1*mmst1
          + 14/3*mmsb1*mmst2
          + 16/3*mmsb1*mmsusy
          - 290/9*sqr(mmsb1)
          - 24*sqr(mmsusy)
          )

       + lnMglSq*den(mmgl - mmsb1,3)*sqr(cs2b) * (
          + 16/9*sqr(mmsb1)*mmsb2
          - 16/9*pow(mmsb1,3)
          )

       + lnMglSq*den(mmgl - mmsb1,3) * (
          - 32*mmsb1*sqr(mmsusy)
          + 8/9*sqr(mmsb1)*mmsb2
          + 8/3*sqr(mmsb1)*mmst1
          + 8/3*sqr(mmsb1)*mmst2
          + 64/3*sqr(mmsb1)*mmsusy
          - 200/9*pow(mmsb1,3)
          )

       + lnMglSq*den(mmgl - mmsb2,1)*sn2b * (
          + 8/9/mb*mmsb1*mgl
          + 176/3/mb*mmsb2*mgl
          + 8/3/mb*mmst1*mgl
          + 8/3/mb*mmst2*mgl
          - 128/3/mb*mmsusy*mgl
          )

       + lnMglSq*den(mmgl - mmsb2,1)*sn4b*cs2b * (
          + 8/9/mb*mmsb1*mgl
          + 8/9/mb*mmsb2*mgl
          )

       + lnMglSq*den(mmgl - mmsb2,1)*sqr(cs2b) * (
          + 812/9*mmsb1
          + 368/9*mmsb2
          )

       + lnMglSq*den(mmgl - mmsb2,1)*den(mmsb1 - mmsb2,1)*sqr(cs2b) * (
          - 796/9*sqr(mmsb1)
          )

       + lnMglSq*den(mmgl - mmsb2,1)*den(mmsb1 - mmsb2,1) * (
          + 76/3*sqr(mmsb1)
          )

       + lnMglSq*den(mmgl - mmsb2,1) * (
          - 226/9*mmsb1
          + 62/9*mmsb2
          + 2*mmst1
          + 2*mmst2
          - 16*mmsusy
          )

       + lnMglSq*den(mmgl - mmsb2,2)*sn2b * (
          + 8/9/mb*mmsb1*mmsb2*mgl
          + 8/3/mb*mmsb2*mmst1*mgl
          + 8/3/mb*mmsb2*mmst2*mgl
          + 64/3/mb*mmsb2*mmsusy*mgl
          - 56/9/mb*sqr(mmsb2)*mgl
          - 32/mb*sqr(mmsusy)*mgl
          )

       + lnMglSq*den(mmgl - mmsb2,2)*sn4b*cs2b * (
          + 8/9/mb*mmsb1*mmsb2*mgl
          - 8/9/mb*sqr(mmsb2)*mgl
          )

       + lnMglSq*den(mmgl - mmsb2,2)*sqr(cs2b) * (
          + 32/9*mmsb1*mmsb2
          - 238/9*sqr(mmsb2)
          )

       + lnMglSq*den(mmgl - mmsb2,2)*den(mmsb1 - mmsb2,1) * (
          + 8/9*pow(mmsb1,3)
          )

       + lnMglSq*den(mmgl - mmsb2,2) * (
          + 2/9*mmsb1*mmsb2
          - 8/9*sqr(mmsb1)
          + 14/3*mmsb2*mmst1
          + 14/3*mmsb2*mmst2
          + 16/3*mmsb2*mmsusy
          - 298/9*sqr(mmsb2)
          - 24*sqr(mmsusy)
          )

       + lnMglSq*den(mmgl - mmsb2,3)*sqr(cs2b) * (
          + 16/9*mmsb1*sqr(mmsb2)
          - 16/9*pow(mmsb2,3)
          )

       + lnMglSq*den(mmgl - mmsb2,3) * (
          + 8/9*mmsb1*sqr(mmsb2)
          - 32*mmsb2*sqr(mmsusy)
          + 8/3*sqr(mmsb2)*mmst1
          + 8/3*sqr(mmsb2)*mmst2
          + 64/3*sqr(mmsb2)*mmsusy
          - 200/9*pow(mmsb2,3)
          )

       + lnMglSq * (
          + 232/3
          )

       + lnMsbSq*sn2b * (
          + 280/9/mb*mgl
          )

       + lnMsbSq*sqr(cs2b) * (
          + 64/9
          )

       + sqr(lnMsbSq)*sn2b * (
          + 8/mb*mgl
          )

       + sqr(lnMsbSq)*den(mmgl - mmsb1,1)*sn2b * (
          - 2/9/mb*mmsb1*mgl
          - 4/9/mb*mmsb2*mgl
          )

       + sqr(lnMsbSq)*den(mmgl - mmsb1,1)*sn4b*cs2b * (
          + 8/9/mb*mmsb1*mgl
          )

       + sqr(lnMsbSq)*den(mmgl - mmsb1,1)*sqr(cs2b) * (
          - 14/9*mmsb1
          - 11/9*mmsb2
          )

       + sqr(lnMsbSq)*den(mmgl - mmsb1,1)*den(mmsb1 - mmsb2,1)*sn4b*cs2b * (
          - 16/9/mb*sqr(mmsb1)*mgl
          )

       + sqr(lnMsbSq)*den(mmgl - mmsb1,1)*den(mmsb1 - mmsb2,1)*sqr(cs2b) * (
          - 77/9*sqr(mmsb1)
          )

       + sqr(lnMsbSq)*den(mmgl - mmsb1,1)*den(mmsb1 - mmsb2,1) * (
          + 13/9*sqr(mmsb1)
          )

       + sqr(lnMsbSq)*den(mmgl - mmsb1,1) * (
          + 419/18*mmsb1
          + mmsb2
          - 2/3*mmst1
          - 2/3*mmst2
          - 16/3*mmsusy
          )

       + sqr(lnMsbSq)*den(mmgl - mmsb1,2)*sn2b * (
          + 2/mb*mmsb1*mmsb2*mgl
          + 2/mb*mmsb1*mmst1*mgl
          + 2/mb*mmsb1*mmst2*mgl
          + 16/mb*mmsb1*mmsusy*mgl
          - 86/9/mb*sqr(mmsb1)*mgl
          - 4/3/mb*sqr(mmsb2)*mgl
          - 4/3/mb*sqr(mmst1)*mgl
          - 4/3/mb*sqr(mmst2)*mgl
          - 32/3/mb*sqr(mmsusy)*mgl
          )

       + sqr(lnMsbSq)*den(mmgl - mmsb1,2)*sn4b*cs2b * (
          + 8/9/mb*sqr(mmsb1)*mgl
          )

       + sqr(lnMsbSq)*den(mmgl - mmsb1,2)*sqr(cs2b) * (
          - 26/9*mmsb1*mmsb2
          - 85/18*sqr(mmsb1)
          )

       + sqr(lnMsbSq)*den(mmgl - mmsb1,2) * (
          + 1/18*mmsb1*mmsb2
          - 17/6*mmsb1*mmst1
          - 17/6*mmsb1*mmst2
          - 68/3*mmsb1*mmsusy
          + 92/3*sqr(sqr(sqr(sqr(mmsb1)
          + mmsb2)
          + mmst1)
          + mmst2)
          + 8*sqr(mmsusy)
          )

       + sqr(lnMsbSq)*den(mmgl - mmsb1,3)*sn2b * (
          - 8/9/mb*pow(mmsb1,3)*mgl
          )

       + sqr(lnMsbSq)*den(mmgl - mmsb1,3)*sqr(cs2b) * (
          - 16/9*pow(mmsb1,3)
          )

       + sqr(lnMsbSq)*den(mmgl - mmsb1,3) * (
          + 4/3*mmsb1*sqr(mmsb2)
          + 4/3*mmsb1*sqr(mmst1)
          + 4/3*mmsb1*sqr(mmst2)
          + 32/3*mmsb1*sqr(mmsusy)
          - 2*sqr(mmsb1)*mmsb2
          - 2*sqr(mmsb1)*mmst1
          - 2*sqr(mmsb1)*mmst2
          - 16*sqr(mmsb1)*mmsusy
          + 110/9*pow(mmsb1,3)
          )

       + sqr(lnMsbSq)*den(mmgl - mmsb1,4) * (
          + 4/9*pow(mmsb1,4)
          )

       + sqr(lnMsbSq)*den(mmgl - mmsb2,1)*sn2b * (
          - 2/3/mb*mmsb1*mgl
          )

       + sqr(lnMsbSq)*den(mmgl - mmsb2,1)*sqr(cs2b) * (
          - 13/9*mmsb1
          )

       + sqr(lnMsbSq)*den(mmgl - mmsb2,1)*den(mmsb1 - mmsb2,1)*sqr(cs2b) * (
          + 13/9*sqr(mmsb1)
          )

       + sqr(lnMsbSq)*den(mmgl - mmsb2,1)*den(mmsb1 - mmsb2,1) * (
          - 13/9*sqr(mmsb1)
          )

       + sqr(lnMsbSq)*den(mmgl - mmsb2,1) * (
          + 17/18*mmsb1
          )

       + sqr(lnMsbSq)*den(mmgl - mmsb2,2) * (
          - 2/3*mmsb1*mmsb2
          )

       + sqr(lnMsbSq)*den(mmsb1 - mmsb2,1)*sqr(cs2b) * (
          - 128/9*mmgl
          - 64/9*mmsb1
          )

       + sqr(lnMsbSq) * (
          + 41/9
          )

       + lnMsbSq*lnMsb2Sq*den(mmgl - mmsb1,1)*sn2b * (
          - 8/9/mb*mmsb1*mgl
          + 8/3/mb*mmsb2*mgl
          )

       + lnMsbSq*lnMsb2Sq*den(mmgl - mmsb1,1)*sn4b*cs2b * (
          - 16/9/mb*mmsb1*mgl
          - 8/9/mb*mmsb2*mgl
          )

       + lnMsbSq*lnMsb2Sq*den(mmgl - mmsb1,1)*sqr(cs2b) * (
          + 5/9*mmsb1
          + 38/9*mmsb2
          )

       + lnMsbSq*lnMsb2Sq*den(mmgl - mmsb1,1)*den(mmsb1 - mmsb2,1)*sn2b
       * (
          + 8/3/mb*sqr(mmsb1)*mgl
          )

       + lnMsbSq*lnMsb2Sq*den(mmgl - mmsb1,1)*den(mmsb1 - mmsb2,1)*sn4b*
      cs2b * (
          + 16/9/mb*sqr(mmsb1)*mgl
          )

       + lnMsbSq*lnMsb2Sq*den(mmgl - mmsb1,1)*den(mmsb1 - mmsb2,1)*sqr(cs2b)
       * (
          - 5/9*sqr(mmsb1)
          )

       + lnMsbSq*lnMsb2Sq*den(mmgl - mmsb1,1)*den(mmsb1 - mmsb2,1) * (
          + 17/9*sqr(mmsb1)
          )

       + lnMsbSq*lnMsb2Sq*den(mmgl - mmsb1,1)*den(mmsb1 - mmsb2,2)*sn2b
       * (
          - 16/9/mb*pow(mmsb1,3)*mgl
          )

       + lnMsbSq*lnMsb2Sq*den(mmgl - mmsb1,1)*den(mmsb1 - mmsb2,2) * (
          + 32/9*pow(mmsb1,3)
          )

       + lnMsbSq*lnMsb2Sq*den(mmgl - mmsb1,1)*den(mmsb1 - mmsb2,3) * (
          - 16/9*pow(mmsb1,4)
          )

       + lnMsbSq*lnMsb2Sq*den(mmgl - mmsb1,1) * (
          - 11/3*mmsb1
          - 34/9*mmsb2
          )

       + lnMsbSq*lnMsb2Sq*den(mmgl - mmsb1,2)*sn2b * (
          - 20/9/mb*mmsb1*mmsb2*mgl
          - 8/9/mb*sqr(mmsb1)*mgl
          + 8/3/mb*sqr(mmsb2)*mgl
          )

       + lnMsbSq*lnMsb2Sq*den(mmgl - mmsb1,2)*sn4b*cs2b * (
          - 8/9/mb*mmsb1*mmsb2*mgl
          )

       + lnMsbSq*lnMsb2Sq*den(mmgl - mmsb1,2)*sqr(cs2b) * (
          + 28/3*mmsb1*mmsb2
          )

       + lnMsbSq*lnMsb2Sq*den(mmgl - mmsb1,2)*den(mmsb1 - mmsb2,1)*sn2b
       * (
          + 8/9/mb*pow(mmsb1,3)*mgl
          )

       + lnMsbSq*lnMsb2Sq*den(mmgl - mmsb1,2)*den(mmsb1 - mmsb2,2) * (
          + 8/9*pow(mmsb1,4)
          )

       + lnMsbSq*lnMsb2Sq*den(mmgl - mmsb1,2) * (
          - 11/3*mmsb1*mmsb2
          - 8/9*sqr(mmsb1)
          - 2*sqr(mmsb2)
          )

       + lnMsbSq*lnMsb2Sq*den(mmgl - mmsb1,3)*sqr(cs2b) * (
          + 16/9*sqr(mmsb1)*mmsb2
          )

       + lnMsbSq*lnMsb2Sq*den(mmgl - mmsb1,3) * (
          - 8/3*mmsb1*sqr(mmsb2)
          + 20/9*sqr(mmsb1)*mmsb2
          )

       + lnMsbSq*lnMsb2Sq*den(mmgl - mmsb2,1)*sn2b * (
          - 8/9/mb*mmsb1*mgl
          - 8/9/mb*mmsb2*mgl
          )

       + lnMsbSq*lnMsb2Sq*den(mmgl - mmsb2,1)*sn4b*cs2b * (
          - 8/9/mb*mmsb1*mgl
          )

       + lnMsbSq*lnMsb2Sq*den(mmgl - mmsb2,1)*sqr(cs2b) * (
          + 11/9*mmsb1
          + 22/9*mmsb2
          )

       + lnMsbSq*lnMsb2Sq*den(mmgl - mmsb2,1)*den(mmsb1 - mmsb2,1)*sn2b
       * (
          - 8/3/mb*sqr(mmsb1)*mgl
          )

       + lnMsbSq*lnMsb2Sq*den(mmgl - mmsb2,1)*den(mmsb1 - mmsb2,1)*sn4b*
      cs2b * (
          + 16/9/mb*sqr(mmsb1)*mgl
          )

       + lnMsbSq*lnMsb2Sq*den(mmgl - mmsb2,1)*den(mmsb1 - mmsb2,1)*sqr(cs2b)
       * (
          + 5/9*sqr(mmsb1)
          )

       + lnMsbSq*lnMsb2Sq*den(mmgl - mmsb2,1)*den(mmsb1 - mmsb2,1) * (
          - 17/9*sqr(mmsb1)
          )

       + lnMsbSq*lnMsb2Sq*den(mmgl - mmsb2,1)*den(mmsb1 - mmsb2,2)*sn2b
       * (
          + 16/9/mb*pow(mmsb1,3)*mgl
          )

       + lnMsbSq*lnMsb2Sq*den(mmgl - mmsb2,1)*den(mmsb1 - mmsb2,2) * (
          - 32/9*pow(mmsb1,3)
          )

       + lnMsbSq*lnMsb2Sq*den(mmgl - mmsb2,1)*den(mmsb1 - mmsb2,3) * (
          + 16/9*pow(mmsb1,4)
          )

       + lnMsbSq*lnMsb2Sq*den(mmgl - mmsb2,1) * (
          + 17/9*mmsb1
          - 2/9*mmsb2
          )

       + lnMsbSq*lnMsb2Sq*den(mmgl - mmsb2,2)*sn2b * (
          - 4/mb*mmsb1*mmsb2*mgl
          - 8/9/mb*sqr(mmsb1)*mgl
          + 8/3/mb*sqr(mmsb2)*mgl
          )

       + lnMsbSq*lnMsb2Sq*den(mmgl - mmsb2,2)*sn4b*cs2b * (
          + 8/9/mb*mmsb1*mmsb2*mgl
          )

       + lnMsbSq*lnMsb2Sq*den(mmgl - mmsb2,2)*sqr(cs2b) * (
          + 32/9*mmsb1*mmsb2
          + 52/9*sqr(mmsb2)
          )

       + lnMsbSq*lnMsb2Sq*den(mmgl - mmsb2,2)*den(mmsb1 - mmsb2,1)*sn2b
       * (
          + 8/9/mb*pow(mmsb1,3)*mgl
          )

       + lnMsbSq*lnMsb2Sq*den(mmgl - mmsb2,2)*den(mmsb1 - mmsb2,1) * (
          - 32/9*pow(mmsb1,3)
          )

       + lnMsbSq*lnMsb2Sq*den(mmgl - mmsb2,2)*den(mmsb1 - mmsb2,2) * (
          + 8/9*pow(mmsb1,4)
          )

       + lnMsbSq*lnMsb2Sq*den(mmgl - mmsb2,2) * (
          - 25/9*mmsb1*mmsb2
          + 8/3*sqr(mmsb1)
          - 10/9*sqr(mmsb2)
          )

       + lnMsbSq*lnMsb2Sq*den(mmgl - mmsb2,3)*sqr(cs2b) * (
          + 16/9*mmsb1*sqr(mmsb2)
          )

       + lnMsbSq*lnMsb2Sq*den(mmgl - mmsb2,3) * (
          - 28/9*mmsb1*sqr(mmsb2)
          + 8/3*pow(mmsb2,3)
          )

       + lnMsbSq*lnMst1Sq*den(mmgl - mmsb1,1) * (
          + 4/3*mmst1
          )

       + lnMsbSq*lnMst1Sq*den(mmgl - mmsb1,2)*sn2b * (
          - 4/mb*mmsb1*mmst1*mgl
          + 8/3/mb*sqr(mmst1)*mgl
          )

       + lnMsbSq*lnMst1Sq*den(mmgl - mmsb1,2) * (
          + 17/3*mmsb1*mmst1
          - 2*sqr(mmst1)
          )

       + lnMsbSq*lnMst1Sq*den(mmgl - mmsb1,3) * (
          - 8/3*mmsb1*sqr(mmst1)
          + 4*sqr(mmsb1)*mmst1
          )

       + lnMsbSq*lnMst2Sq*den(mmgl - mmsb1,1) * (
          + 4/3*mmst2
          )

       + lnMsbSq*lnMst2Sq*den(mmgl - mmsb1,2)*sn2b * (
          - 4/mb*mmsb1*mmst2*mgl
          + 8/3/mb*sqr(mmst2)*mgl
          )

       + lnMsbSq*lnMst2Sq*den(mmgl - mmsb1,2) * (
          + 17/3*mmsb1*mmst2
          - 2*sqr(mmst2)
          )

       + lnMsbSq*lnMst2Sq*den(mmgl - mmsb1,3) * (
          - 8/3*mmsb1*sqr(mmst2)
          + 4*sqr(mmsb1)*mmst2
          )

       + lnMsbSq*lnMmsusy*den(mmgl - mmsb1,1) * (
          + 32/3*mmsusy
          )

       + lnMsbSq*lnMmsusy*den(mmgl - mmsb1,2)*sn2b * (
          - 32/mb*mmsb1*mmsusy*mgl
          + 64/3/mb*sqr(mmsusy)*mgl
          )

       + lnMsbSq*lnMmsusy*den(mmgl - mmsb1,2) * (
          + 136/3*mmsb1*mmsusy
          - 16*sqr(mmsusy)
          )

       + lnMsbSq*lnMmsusy*den(mmgl - mmsb1,3) * (
          - 64/3*mmsb1*sqr(mmsusy)
          + 32*sqr(mmsb1)*mmsusy
          )

       + lnMsbSq*lnMmu*sn2b * (
          + 64/9/mb*mgl
          )

       + lnMsbSq*lnMmu*sqr(cs2b) * (
          + 64/9
          )

       + lnMsbSq*lnMmu*den(mmgl - mmsb1,1)*sn2b * (
          + 16/mb*mmsb1*mgl
          - 16/9/mb*mmsb2*mgl
          )

       + lnMsbSq*lnMmu*den(mmgl - mmsb1,1)*sn4b*cs2b * (
          + 8/9/mb*mmsb1*mgl
          + 8/9/mb*mmsb2*mgl
          )

       + lnMsbSq*lnMmu*den(mmgl - mmsb1,1)*sqr(cs2b) * (
          + 16*mmsb1
          - 16/9*mmsb2
          )

       + lnMsbSq*lnMmu*den(mmgl - mmsb1,1)*den(mmsb1 - mmsb2,1)*sqr(cs2b)
       * (
          - 128/9*sqr(mmsb1)
          )

       + lnMsbSq*lnMmu*den(mmgl - mmsb1,1) * (
          - 332/9*mmsb1
          + 16/9*mmsb2
          )

       + lnMsbSq*lnMmu*den(mmgl - mmsb1,2)*sn2b * (
          - 16/9/mb*mmsb1*mmsb2*mgl
          + 8/9/mb*sqr(mmsb1)*mgl
          )

       + lnMsbSq*lnMmu*den(mmgl - mmsb1,2)*sn4b*cs2b * (
          + 8/9/mb*mmsb1*mmsb2*mgl
          - 8/9/mb*sqr(mmsb1)*mgl
          )

       + lnMsbSq*lnMmu*den(mmgl - mmsb1,2)*sqr(cs2b) * (
          - 32/9*mmsb1*mmsb2
          + 32/3*sqr(mmsb1)
          )

       + lnMsbSq*lnMmu*den(mmgl - mmsb1,2) * (
          + 32/9*mmsb1*mmsb2
          - 178/9*sqr(mmsb1)
          )

       + lnMsbSq*lnMmu*den(mmgl - mmsb1,3)*sqr(cs2b) * (
          - 16/9*sqr(mmsb1)*mmsb2
          + 16/9*pow(mmsb1,3)
          )

       + lnMsbSq*lnMmu*den(mmgl - mmsb1,3) * (
          + 16/9*sqr(mmsb1)*mmsb2
          - 8/9*pow(mmsb1,3)
          )

       + lnMsbSq*lnMmu*den(mmsb1 - mmsb2,1)*sqr(cs2b) * (
          - 128/9*mmsb1
          )

       + lnMsbSq*lnMmu * (
          - 128/9
          )

       + lnMsbSq*den(mmgl - mmsb1,1)*sn2b * (
          + 172/3/mb*mmsb1*mgl
          - 28/9/mb*mmsb2*mgl
          )

       + lnMsbSq*den(mmgl - mmsb1,1)*sn4b*cs2b * (
          + 16/9/mb*mmsb1*mgl
          + 8/9/mb*mmsb2*mgl
          )

       + lnMsbSq*den(mmgl - mmsb1,1)*sqr(cs2b) * (
          + 229/9*mmsb1
          - 49/9*mmsb2
          )

       + lnMsbSq*den(mmgl - mmsb1,1)*den(mmsb1 - mmsb2,1)*sn2b * (
          - 8/9/mb*sqr(mmsb1)*mgl
          )

       + lnMsbSq*den(mmgl - mmsb1,1)*den(mmsb1 - mmsb2,1)*sqr(cs2b) * (
          - 718/9*sqr(mmsb1)
          )

       + lnMsbSq*den(mmgl - mmsb1,1)*den(mmsb1 - mmsb2,1) * (
          + 34/3*sqr(mmsb1)
          )

       + lnMsbSq*den(mmgl - mmsb1,1)*den(mmsb1 - mmsb2,2) * (
          - 8/9*pow(mmsb1,3)
          )

       + lnMsbSq*den(mmgl - mmsb1,1) * (
          - 13/3*mmsb1
          + 43/9*mmsb2
          - 2*mmst1
          - 2*mmst2
          - 16*mmsusy
          )

       + lnMsbSq*den(mmgl - mmsb1,2)*sn2b * (
          + 8/9/mb*mmsb1*mmsb2*mgl
          + 8/3/mb*mmsb1*mmst1*mgl
          + 8/3/mb*mmsb1*mmst2*mgl
          + 64/3/mb*mmsb1*mmsusy*mgl
          + 52/9/mb*sqr(mmsb1)*mgl
          - 4/mb*sqr(mmsb2)*mgl
          - 4/mb*sqr(mmst1)*mgl
          - 4/mb*sqr(mmst2)*mgl
          - 32/mb*sqr(mmsusy)*mgl
          )

       + lnMsbSq*den(mmgl - mmsb1,2)*sn4b*cs2b * (
          + 8/9/mb*mmsb1*mmsb2*mgl
          - 8/9/mb*sqr(mmsb1)*mgl
          )

       + lnMsbSq*den(mmgl - mmsb1,2)*sqr(cs2b) * (
          - 110/9*mmsb1*mmsb2
          + 16*sqr(mmsb1)
          )

       + lnMsbSq*den(mmgl - mmsb1,2)*den(mmsb1 - mmsb2,1) * (
          + 8/9*pow(mmsb1,3)
          )

       + lnMsbSq*den(mmgl - mmsb1,2) * (
          + 56/9*mmsb1*mmsb2
          - 6*mmsb1*mmst1
          - 6*mmsb1*mmst2
          - 48*mmsb1*mmsusy
          + 187/9*sqr(mmsb1)
          + 3*sqr(mmsb2)
          + 3*sqr(mmst1)
          + 3*sqr(mmst2)
          + 24*sqr(mmsusy)
          )

       + lnMsbSq*den(mmgl - mmsb1,3)*sqr(cs2b) * (
          - 16/9*sqr(mmsb1)*mmsb2
          + 16/9*pow(mmsb1,3)
          )

       + lnMsbSq*den(mmgl - mmsb1,3) * (
          + 4*mmsb1*sqr(mmsb2)
          + 4*mmsb1*sqr(mmst1)
          + 4*mmsb1*sqr(mmst2)
          + 32*mmsb1*sqr(mmsusy)
          - 8/9*sqr(mmsb1)*mmsb2
          - 8/3*sqr(mmsb1)*mmst1
          - 8/3*sqr(mmsb1)*mmst2
          - 64/3*sqr(mmsb1)*mmsusy
          + 92/9*pow(mmsb1,3)
          )

       + lnMsbSq*den(mmgl - mmsb2,1)*sn2b * (
          - 16/3/mb*mmsb1*mgl
          )

       + lnMsbSq*den(mmgl - mmsb2,1)*sn4b*cs2b * (
          + 8/9/mb*mmsb1*mgl
          )

       + lnMsbSq*den(mmgl - mmsb2,1)*sqr(cs2b) * (
          - 6*mmsb1
          )

       + lnMsbSq*den(mmgl - mmsb2,1)*den(mmsb1 - mmsb2,1)*sn2b * (
          + 8/9/mb*sqr(mmsb1)*mgl
          )

       + lnMsbSq*den(mmgl - mmsb2,1)*den(mmsb1 - mmsb2,1)*sqr(cs2b) * (
          + 26/3*sqr(mmsb1)
          )

       + lnMsbSq*den(mmgl - mmsb2,1)*den(mmsb1 - mmsb2,1) * (
          - 34/3*sqr(mmsb1)
          )

       + lnMsbSq*den(mmgl - mmsb2,1)*den(mmsb1 - mmsb2,2) * (
          + 8/9*pow(mmsb1,3)
          )

       + lnMsbSq*den(mmgl - mmsb2,1) * (
          + 52/9*mmsb1
          )

       + lnMsbSq*den(mmgl - mmsb2,2)*sqr(cs2b) * (
          + 16/9*mmsb1*mmsb2
          )

       + lnMsbSq*den(mmgl - mmsb2,2) * (
          - 40/9*mmsb1*mmsb2
          )

       + lnMsbSq*den(mmsb1 - mmsb2,1)*sqr(cs2b) * (
          - 128/3*mmgl
          - 640/9*mmsb1
          )

       + lnMsbSq * (
          + 5/9
          )

       + lnMsb2Sq*sn2b * (
          - 280/9/mb*mgl
          )

       + lnMsb2Sq*sqr(cs2b) * (
          - 64
          )

       + sqr(lnMsb2Sq)*sn2b * (
          - 8/mb*mgl
          )

       + sqr(lnMsb2Sq)*sqr(cs2b) * (
          - 64/9
          )

       + sqr(lnMsb2Sq)*den(mmgl - mmsb1,1)*sn2b * (
          + 4/9/mb*mmsb1*mgl
          + 2/9/mb*mmsb2*mgl
          )

       + sqr(lnMsb2Sq)*den(mmgl - mmsb1,1)*sqr(cs2b) * (
          + 8/3*mmsb1
          - 11/9*mmsb2
          )

       + sqr(lnMsb2Sq)*den(mmgl - mmsb1,1)*den(mmsb1 - mmsb2,1)*sqr(cs2b) * (
          - 13/9*sqr(mmsb1)
          )

       + sqr(lnMsb2Sq)*den(mmgl - mmsb1,1)*den(mmsb1 - mmsb2,1) * (
          + 13/9*sqr(mmsb1)
          )

       + sqr(lnMsb2Sq)*den(mmgl - mmsb1,1) * (
          - 14/9*mmsb1
          + 1/2*mmsb2
          )

       + sqr(lnMsb2Sq)*den(mmgl - mmsb1,2)*sn2b * (
          + 8/3/mb*mmsb1*mmsb2*mgl
          - 4/3/mb*sqr(mmsb1)*mgl
          - 4/3/mb*sqr(mmsb2)*mgl
          )

       + sqr(lnMsb2Sq)*den(mmgl - mmsb1,2)*sqr(cs2b) * (
          - 26/9*mmsb1*mmsb2
          + 26/9*sqr(mmsb1)
          )

       + sqr(lnMsb2Sq)*den(mmgl - mmsb1,2) * (
          - 10/9*mmsb1*mmsb2
          - 5/9*sqr(sqr(mmsb1)
          + mmsb2)
          )

       + sqr(lnMsb2Sq)*den(mmgl - mmsb1,3) * (
          + 4/3*mmsb1*sqr(mmsb2)
          - 8/3*sqr(mmsb1)*mmsb2
          + 4/3*pow(mmsb1,3)
          )

       + sqr(lnMsb2Sq)*den(mmgl - mmsb2,1)*sn2b * (
          + 2/3/mb*mmsb2*mgl
          )

       + sqr(lnMsb2Sq)*den(mmgl - mmsb2,1)*sn4b*cs2b * (
          + 16/9/mb*mmsb1*mgl
          + 8/9/mb*mmsb2*mgl
          )

       + sqr(lnMsb2Sq)*den(mmgl - mmsb2,1)*sqr(cs2b) * (
          - 77/9*mmsb1
          - 34/3*mmsb2
          )

       + sqr(lnMsb2Sq)*den(mmgl - mmsb2,1)*den(mmsb1 - mmsb2,1)*sn4b*cs2b * (
          - 16/9/mb*sqr(mmsb1)*mgl
          )

       + sqr(lnMsb2Sq)*den(mmgl - mmsb2,1)*den(mmsb1 - mmsb2,1)*sqr(cs2b) * (
          + 77/9*sqr(mmsb1)
          )

       + sqr(lnMsb2Sq)*den(mmgl - mmsb2,1)*den(mmsb1 - mmsb2,1) * (
          - 13/9*sqr(mmsb1)
          )

       + sqr(lnMsb2Sq)*den(mmgl - mmsb2,1) * (
          + 13/9*mmsb1
          + 149/6*mmsb2
          - 2/3*mmst1
          - 2/3*mmst2
          - 16/3*mmsusy
          )

       + sqr(lnMsb2Sq)*den(mmgl - mmsb2,2)*sn2b * (
          + 2/3/mb*mmsb1*mmsb2*mgl
          - 2/mb*mmsb2*mmst1*mgl
          - 2/mb*mmsb2*mmst2*mgl
          - 16/mb*mmsb2*mmsusy*mgl
          + 74/9/mb*sqr(mmsb2)*mgl
          + 4/3/mb*sqr(mmst1)*mgl
          + 4/3/mb*sqr(mmst2)*mgl
          + 32/3/mb*sqr(mmsusy)*mgl
          )

       + sqr(lnMsb2Sq)*den(mmgl - mmsb2,2)*sn4b*cs2b * (
          - 8/9/mb*sqr(mmsb2)*mgl
          )

       + sqr(lnMsb2Sq)*den(mmgl - mmsb2,2)*sqr(cs2b) * (
          - 137/18*sqr(mmsb2)
          )

       + sqr(lnMsb2Sq)*den(mmgl - mmsb2,2) * (
          + 1/2*mmsb1*mmsb2
          - 17/6*mmsb2*mmst1
          - 17/6*mmsb2*mmst2
          - 68/3*mmsb2*mmsusy
          + 281/9*sqr(sqr(sqr(mmsb2)
          + mmst1)
          + mmst2)
          + 8*sqr(mmsusy)
          )

       + sqr(lnMsb2Sq)*den(mmgl - mmsb2,3)*sn2b * (
          + 8/9/mb*pow(mmsb2,3)*mgl
          )

       + sqr(lnMsb2Sq)*den(mmgl - mmsb2,3)*sqr(cs2b) * (
          - 16/9*pow(mmsb2,3)
          )

       + sqr(lnMsb2Sq)*den(mmgl - mmsb2,3) * (
          + 2/3*mmsb1*sqr(mmsb2)
          + 4/3*mmsb2*sqr(mmst1)
          + 4/3*mmsb2*sqr(mmst2)
          + 32/3*mmsb2*sqr(mmsusy)
          - 2*sqr(mmsb2)*mmst1
          - 2*sqr(mmsb2)*mmst2
          - 16*sqr(mmsb2)*mmsusy
          + 98/9*pow(mmsb2,3)
          )

       + sqr(lnMsb2Sq)*den(mmgl - mmsb2,4) * (
          + 4/9*pow(mmsb2,4)
          )

       + sqr(lnMsb2Sq)*den(mmsb1 - mmsb2,1)*sqr(cs2b) * (
          + 128/9*mmgl
          + 64/9*mmsb1
          )

       + sqr(lnMsb2Sq) * (
          + 41/9
          )

       + lnMsb2Sq*lnMst1Sq*den(mmgl - mmsb2,1) * (
          + 4/3*mmst1
          )

       + lnMsb2Sq*lnMst1Sq*den(mmgl - mmsb2,2)*sn2b * (
          + 4/mb*mmsb2*mmst1*mgl
          - 8/3/mb*sqr(mmst1)*mgl
          )

       + lnMsb2Sq*lnMst1Sq*den(mmgl - mmsb2,2) * (
          + 17/3*mmsb2*mmst1
          - 2*sqr(mmst1)
          )

       + lnMsb2Sq*lnMst1Sq*den(mmgl - mmsb2,3) * (
          - 8/3*mmsb2*sqr(mmst1)
          + 4*sqr(mmsb2)*mmst1
          )

       + lnMsb2Sq*lnMst2Sq*den(mmgl - mmsb2,1) * (
          + 4/3*mmst2
          )

       + lnMsb2Sq*lnMst2Sq*den(mmgl - mmsb2,2)*sn2b * (
          + 4/mb*mmsb2*mmst2*mgl
          - 8/3/mb*sqr(mmst2)*mgl
          )

       + lnMsb2Sq*lnMst2Sq*den(mmgl - mmsb2,2) * (
          + 17/3*mmsb2*mmst2
          - 2*sqr(mmst2)
          )

       + lnMsb2Sq*lnMst2Sq*den(mmgl - mmsb2,3) * (
          - 8/3*mmsb2*sqr(mmst2)
          + 4*sqr(mmsb2)*mmst2
          )

       + lnMsb2Sq*lnMmsusy*den(mmgl - mmsb2,1) * (
          + 32/3*mmsusy
          )

       + lnMsb2Sq*lnMmsusy*den(mmgl - mmsb2,2)*sn2b * (
          + 32/mb*mmsb2*mmsusy*mgl
          - 64/3/mb*sqr(mmsusy)*mgl
          )

       + lnMsb2Sq*lnMmsusy*den(mmgl - mmsb2,2) * (
          + 136/3*mmsb2*mmsusy
          - 16*sqr(mmsusy)
          )

       + lnMsb2Sq*lnMmsusy*den(mmgl - mmsb2,3) * (
          - 64/3*mmsb2*sqr(mmsusy)
          + 32*sqr(mmsb2)*mmsusy
          )

       + lnMsb2Sq*lnMmu*sn2b * (
          - 64/9/mb*mgl
          )

       + lnMsb2Sq*lnMmu*sqr(cs2b) * (
          - 64/9
          )

       + lnMsb2Sq*lnMmu*den(mmgl - mmsb2,1)*sn2b * (
          + 16/9/mb*mmsb1*mgl
          - 16/mb*mmsb2*mgl
          )

       + lnMsb2Sq*lnMmu*den(mmgl - mmsb2,1)*sn4b*cs2b * (
          - 8/9/mb*mmsb1*mgl
          - 8/9/mb*mmsb2*mgl
          )

       + lnMsb2Sq*lnMmu*den(mmgl - mmsb2,1)*sqr(cs2b) * (
          - 16*mmsb1
          + 16/9*mmsb2
          )

       + lnMsb2Sq*lnMmu*den(mmgl - mmsb2,1)*den(mmsb1 - mmsb2,1)*sqr(cs2b)
       * (
          + 128/9*sqr(mmsb1)
          )

       + lnMsb2Sq*lnMmu*den(mmgl - mmsb2,1) * (
          + 16/9*mmsb1
          - 332/9*mmsb2
          )

       + lnMsb2Sq*lnMmu*den(mmgl - mmsb2,2)*sn2b * (
          + 16/9/mb*mmsb1*mmsb2*mgl
          - 8/9/mb*sqr(mmsb2)*mgl
          )

       + lnMsb2Sq*lnMmu*den(mmgl - mmsb2,2)*sn4b*cs2b * (
          - 8/9/mb*mmsb1*mmsb2*mgl
          + 8/9/mb*sqr(mmsb2)*mgl
          )

       + lnMsb2Sq*lnMmu*den(mmgl - mmsb2,2)*sqr(cs2b) * (
          - 32/9*mmsb1*mmsb2
          + 32/3*sqr(mmsb2)
          )

       + lnMsb2Sq*lnMmu*den(mmgl - mmsb2,2) * (
          + 32/9*mmsb1*mmsb2
          - 178/9*sqr(mmsb2)
          )

       + lnMsb2Sq*lnMmu*den(mmgl - mmsb2,3)*sqr(cs2b) * (
          - 16/9*mmsb1*sqr(mmsb2)
          + 16/9*pow(mmsb2,3)
          )

       + lnMsb2Sq*lnMmu*den(mmgl - mmsb2,3) * (
          + 16/9*mmsb1*sqr(mmsb2)
          - 8/9*pow(mmsb2,3)
          )

       + lnMsb2Sq*lnMmu*den(mmsb1 - mmsb2,1)*sqr(cs2b) * (
          + 128/9*mmsb1
          )

       + lnMsb2Sq*lnMmu * (
          - 128/9
          )

       + lnMsb2Sq*den(mmgl - mmsb1,1)*sn2b * (
          + 4/9/mb*mmsb1*mgl
          + 52/9/mb*mmsb2*mgl
          )

       + lnMsb2Sq*den(mmgl - mmsb1,1)*sn4b*cs2b * (
          - 8/9/mb*mmsb2*mgl
          )

       + lnMsb2Sq*den(mmgl - mmsb1,1)*sqr(cs2b) * (
          + 37/3*mmsb1
          + 19/3*mmsb2
          )

       + lnMsb2Sq*den(mmgl - mmsb1,1)*den(mmsb1 - mmsb2,1)*sn2b * (
          + 8/9/mb*sqr(mmsb1)*mgl
          )

       + lnMsb2Sq*den(mmgl - mmsb1,1)*den(mmsb1 - mmsb2,1)*sqr(cs2b) * (
          - 26/3*sqr(mmsb1)
          )

       + lnMsb2Sq*den(mmgl - mmsb1,1)*den(mmsb1 - mmsb2,1) * (
          + 14*sqr(mmsb1)
          )

       + lnMsb2Sq*den(mmgl - mmsb1,1)*den(mmsb1 - mmsb2,2) * (
          + 8/9*pow(mmsb1,3)
          )

       + lnMsb2Sq*den(mmgl - mmsb1,1) * (
          - 137/9*mmsb1
          - 23/3*mmsb2
          )

       + lnMsb2Sq*den(mmgl - mmsb1,2)*sn2b * (
          - 4/mb*sqr(mmsb1)*mgl
          + 4/mb*sqr(mmsb2)*mgl
          )

       + lnMsb2Sq*den(mmgl - mmsb1,2)*sqr(cs2b) * (
          + 94/9*mmsb1*mmsb2
          + 26/3*sqr(mmsb1)
          )

       + lnMsb2Sq*den(mmgl - mmsb1,2) * (
          - 82/9*mmsb1*mmsb2
          - 5/3*sqr(mmsb1)
          - 3*sqr(mmsb2)
          )

       + lnMsb2Sq*den(mmgl - mmsb1,3) * (
          - 4*mmsb1*sqr(mmsb2)
          + 4*pow(mmsb1,3)
          )

       + lnMsb2Sq*den(mmgl - mmsb2,1)*sn2b * (
          + 8/3/mb*mmsb1*mgl
          - 520/9/mb*mmsb2*mgl
          )

       + lnMsb2Sq*den(mmgl - mmsb2,1)*sn4b*cs2b * (
          - 8/9/mb*mmsb1*mgl
          - 16/9/mb*mmsb2*mgl
          )

       + lnMsb2Sq*den(mmgl - mmsb2,1)*sqr(cs2b) * (
          - 734/9*mmsb1
          - 152/3*mmsb2
          )

       + lnMsb2Sq*den(mmgl - mmsb2,1)*den(mmsb1 - mmsb2,1)*sn2b * (
          - 8/9/mb*sqr(mmsb1)*mgl
          )

       + lnMsb2Sq*den(mmgl - mmsb2,1)*den(mmsb1 - mmsb2,1)*sqr(cs2b) * (
          + 718/9*sqr(mmsb1)
          )

       + lnMsb2Sq*den(mmgl - mmsb2,1)*den(mmsb1 - mmsb2,1) * (
          - 14*sqr(mmsb1)
          )

       + lnMsb2Sq*den(mmgl - mmsb2,1)*den(mmsb1 - mmsb2,2) * (
          - 8/9*pow(mmsb1,3)
          )

       + lnMsb2Sq*den(mmgl - mmsb2,1) * (
          + 50/3*mmsb1
          + 52/9*mmsb2
          - 2*mmst1
          - 2*mmst2
          - 16*mmsusy
          )

       + lnMsb2Sq*den(mmgl - mmsb2,2)*sn2b * (
          - 8/9/mb*mmsb1*mmsb2*mgl
          - 8/3/mb*mmsb2*mmst1*mgl
          - 8/3/mb*mmsb2*mmst2*mgl
          - 64/3/mb*mmsb2*mmsusy*mgl
          - 16/9/mb*sqr(mmsb2)*mgl
          + 4/mb*sqr(mmst1)*mgl
          + 4/mb*sqr(mmst2)*mgl
          + 32/mb*sqr(mmsusy)*mgl
          )

       + lnMsb2Sq*den(mmgl - mmsb2,2)*sn4b*cs2b * (
          - 8/9/mb*mmsb1*mmsb2*mgl
          + 8/9/mb*sqr(mmsb2)*mgl
          )

       + lnMsb2Sq*den(mmgl - mmsb2,2)*sqr(cs2b) * (
          - 32/9*mmsb1*mmsb2
          + 74/3*sqr(mmsb2)
          )

       + lnMsb2Sq*den(mmgl - mmsb2,2)*den(mmsb1 - mmsb2,1) * (
          - 8/9*pow(mmsb1,3)
          )

       + lnMsb2Sq*den(mmgl - mmsb2,2) * (
          + 22/9*mmsb1*mmsb2
          + 8/9*sqr(mmsb1)
          - 6*mmsb2*mmst1
          - 6*mmsb2*mmst2
          - 48*mmsb2*mmsusy
          + 20*sqr(mmsb2)
          + 3*sqr(mmst1)
          + 3*sqr(mmst2)
          + 24*sqr(mmsusy)
          )

       + lnMsb2Sq*den(mmgl - mmsb2,3)*sqr(cs2b) * (
          - 16/9*mmsb1*sqr(mmsb2)
          + 16/9*pow(mmsb2,3)
          )

       + lnMsb2Sq*den(mmgl - mmsb2,3) * (
          - 8/9*mmsb1*sqr(mmsb2)
          + 4*mmsb2*sqr(mmst1)
          + 4*mmsb2*sqr(mmst2)
          + 32*mmsb2*sqr(mmsusy)
          - 8/3*sqr(mmsb2)*mmst1
          - 8/3*sqr(mmsb2)*mmst2
          - 64/3*sqr(mmsb2)*mmsusy
          + 128/9*pow(mmsb2,3)
          )

       + lnMsb2Sq*den(mmsb1 - mmsb2,1)*sqr(cs2b) * (
          + 128/3*mmgl
          + 640/9*mmsb1
          )

       + lnMsb2Sq * (
          + 5/9
          )

       + sqr(lnMst1Sq)*den(mmgl - mmsb1,1)*sn2b * (
          + 2/3/mb*mmst1*mgl
          )

       + sqr(lnMst1Sq)*den(mmgl - mmsb1,1) * (
          + 2/3*mmsb1
          - 7/6*mmst1
          )

       + sqr(lnMst1Sq)*den(mmgl - mmsb1,2)*sn2b * (
          + 8/3/mb*mmsb1*mmst1*mgl
          - 4/3/mb*sqr(mmsb1)*mgl
          - 4/3/mb*sqr(mmst1)*mgl
          )

       + sqr(lnMst1Sq)*den(mmgl - mmsb1,2) * (
          - 4*mmsb1*mmst1
          + 7/3*sqr(sqr(mmsb1)
          + mmst1)
          )

       + sqr(lnMst1Sq)*den(mmgl - mmsb1,3) * (
          + 4/3*mmsb1*sqr(mmst1)
          - 8/3*sqr(mmsb1)*mmst1
          + 4/3*pow(mmsb1,3)
          )

       + sqr(lnMst1Sq)*den(mmgl - mmsb2,1)*sn2b * (
          - 2/3/mb*mmst1*mgl
          )

       + sqr(lnMst1Sq)*den(mmgl - mmsb2,1) * (
          + 2/3*mmsb2
          - 7/6*mmst1
          )

       + sqr(lnMst1Sq)*den(mmgl - mmsb2,2)*sn2b * (
          - 8/3/mb*mmsb2*mmst1*mgl
          + 4/3/mb*sqr(mmsb2)*mgl
          + 4/3/mb*sqr(mmst1)*mgl
          )

       + sqr(lnMst1Sq)*den(mmgl - mmsb2,2) * (
          - 4*mmsb2*mmst1
          + 7/3*sqr(sqr(mmsb2)
          + mmst1)
          )

       + sqr(lnMst1Sq)*den(mmgl - mmsb2,3) * (
          + 4/3*mmsb2*sqr(mmst1)
          - 8/3*sqr(mmsb2)*mmst1
          + 4/3*pow(mmsb2,3)
          )

       + sqr(lnMst1Sq) * (
          - 1/3
          )

       + lnMst1Sq*den(mmgl - mmsb1,1)*sn2b * (
          + 8/3/mb*mmst1*mgl
          )

       + lnMst1Sq*den(mmgl - mmsb1,1) * (
          + 2*mmsb1
          )

       + lnMst1Sq*den(mmgl - mmsb1,2)*sn2b * (
          - 4/mb*sqr(mmsb1)*mgl
          + 4/mb*sqr(mmst1)*mgl
          )

       + lnMst1Sq*den(mmgl - mmsb1,2) * (
          + 4/3*mmsb1*mmst1
          + 7*sqr(mmsb1)
          - 3*sqr(mmst1)
          )

       + lnMst1Sq*den(mmgl - mmsb1,3) * (
          - 4*mmsb1*sqr(mmst1)
          + 4*pow(mmsb1,3)
          )

       + lnMst1Sq*den(mmgl - mmsb2,1)*sn2b * (
          - 8/3/mb*mmst1*mgl
          )

       + lnMst1Sq*den(mmgl - mmsb2,1) * (
          + 2*mmsb2
          )

       + lnMst1Sq*den(mmgl - mmsb2,2)*sn2b * (
          + 4/mb*sqr(mmsb2)*mgl
          - 4/mb*sqr(mmst1)*mgl
          )

       + lnMst1Sq*den(mmgl - mmsb2,2) * (
          + 4/3*mmsb2*mmst1
          + 7*sqr(mmsb2)
          - 3*sqr(mmst1)
          )

       + lnMst1Sq*den(mmgl - mmsb2,3) * (
          - 4*mmsb2*sqr(mmst1)
          + 4*pow(mmsb2,3)
          )

       + lnMst1Sq * (
          + 1/9
          )

       + sqr(lnMst2Sq)*den(mmgl - mmsb1,1)*sn2b * (
          + 2/3/mb*mmst2*mgl
          )

       + sqr(lnMst2Sq)*den(mmgl - mmsb1,1) * (
          + 2/3*mmsb1
          - 7/6*mmst2
          )

       + sqr(lnMst2Sq)*den(mmgl - mmsb1,2)*sn2b * (
          + 8/3/mb*mmsb1*mmst2*mgl
          - 4/3/mb*sqr(mmsb1)*mgl
          - 4/3/mb*sqr(mmst2)*mgl
          )

       + sqr(lnMst2Sq)*den(mmgl - mmsb1,2) * (
          - 4*mmsb1*mmst2
          + 7/3*sqr(sqr(mmsb1)
          + mmst2)
          )

       + sqr(lnMst2Sq)*den(mmgl - mmsb1,3) * (
          + 4/3*mmsb1*sqr(mmst2)
          - 8/3*sqr(mmsb1)*mmst2
          + 4/3*pow(mmsb1,3)
          )

       + sqr(lnMst2Sq)*den(mmgl - mmsb2,1)*sn2b * (
          - 2/3/mb*mmst2*mgl
          )

       + sqr(lnMst2Sq)*den(mmgl - mmsb2,1) * (
          + 2/3*mmsb2
          - 7/6*mmst2
          )

       + sqr(lnMst2Sq)*den(mmgl - mmsb2,2)*sn2b * (
          - 8/3/mb*mmsb2*mmst2*mgl
          + 4/3/mb*sqr(mmsb2)*mgl
          + 4/3/mb*sqr(mmst2)*mgl
          )

       + sqr(lnMst2Sq)*den(mmgl - mmsb2,2) * (
          - 4*mmsb2*mmst2
          + 7/3*sqr(sqr(mmsb2)
          + mmst2)
          )

       + sqr(lnMst2Sq)*den(mmgl - mmsb2,3) * (
          + 4/3*mmsb2*sqr(mmst2)
          - 8/3*sqr(mmsb2)*mmst2
          + 4/3*pow(mmsb2,3)
          )

       + sqr(lnMst2Sq) * (
          - 1/3
          )

       + lnMst2Sq*den(mmgl - mmsb1,1)*sn2b * (
          + 8/3/mb*mmst2*mgl
          )

       + lnMst2Sq*den(mmgl - mmsb1,1) * (
          + 2*mmsb1
          )

       + lnMst2Sq*den(mmgl - mmsb1,2)*sn2b * (
          - 4/mb*sqr(mmsb1)*mgl
          + 4/mb*sqr(mmst2)*mgl
          )

       + lnMst2Sq*den(mmgl - mmsb1,2) * (
          + 4/3*mmsb1*mmst2
          + 7*sqr(mmsb1)
          - 3*sqr(mmst2)
          )

       + lnMst2Sq*den(mmgl - mmsb1,3) * (
          - 4*mmsb1*sqr(mmst2)
          + 4*pow(mmsb1,3)
          )

       + lnMst2Sq*den(mmgl - mmsb2,1)*sn2b * (
          - 8/3/mb*mmst2*mgl
          )

       + lnMst2Sq*den(mmgl - mmsb2,1) * (
          + 2*mmsb2
          )

       + lnMst2Sq*den(mmgl - mmsb2,2)*sn2b * (
          + 4/mb*sqr(mmsb2)*mgl
          - 4/mb*sqr(mmst2)*mgl
          )

       + lnMst2Sq*den(mmgl - mmsb2,2) * (
          + 4/3*mmsb2*mmst2
          + 7*sqr(mmsb2)
          - 3*sqr(mmst2)
          )

       + lnMst2Sq*den(mmgl - mmsb2,3) * (
          - 4*mmsb2*sqr(mmst2)
          + 4*pow(mmsb2,3)
          )

       + lnMst2Sq * (
          + 1/9
          )

       + sqr(lnMmsusy)*den(mmgl - mmsb1,1)*sn2b * (
          + 16/3/mb*mmsusy*mgl
          )

       + sqr(lnMmsusy)*den(mmgl - mmsb1,1) * (
          - 4*mmsusy
          )

       + sqr(lnMmsusy)*den(mmgl - mmsb1,2) * (
          - 16/3*mmsb1*mmsusy
          )

       + sqr(lnMmsusy)*den(mmgl - mmsb2,1)*sn2b * (
          - 16/3/mb*mmsusy*mgl
          )

       + sqr(lnMmsusy)*den(mmgl - mmsb2,1) * (
          - 4*mmsusy
          )

       + sqr(lnMmsusy)*den(mmgl - mmsb2,2) * (
          - 16/3*mmsb2*mmsusy
          )

       + sqr(lnMmsusy) * (
          + 8/3
          )

       + lnMmsusy*den(mmgl - mmsb1,1)*sn2b * (
          - 128/3/mb*mmsusy*mgl
          )

       + lnMmsusy*den(mmgl - mmsb1,1) * (
          + 32*mmsusy
          )

       + lnMmsusy*den(mmgl - mmsb1,2) * (
          + 128/3*mmsb1*mmsusy
          )

       + lnMmsusy*den(mmgl - mmsb2,1)*sn2b * (
          + 128/3/mb*mmsusy*mgl
          )

       + lnMmsusy*den(mmgl - mmsb2,1) * (
          + 32*mmsusy
          )

       + lnMmsusy*den(mmgl - mmsb2,2) * (
          + 128/3*mmsb2*mmsusy
          )

       + lnMmsusy * (
          + 152/9
          )

       + lnMmu*sqr(cs2b) * (
          + 128/9
          )

       + sqr(lnMmu) * (
          - 130/9
          )

       + lnMmu*den(mmgl - mmsb1,1)*sn2b * (
          + 8/9/mb*mmsb1*mgl
          - 16/9/mb*mmsb2*mgl
          )

       + lnMmu*den(mmgl - mmsb1,1)*sn4b*cs2b * (
          - 8/9/mb*mmsb1*mgl
          + 8/9/mb*mmsb2*mgl
          )

       + lnMmu*den(mmgl - mmsb1,1)*sqr(cs2b) * (
          + 88/9*mmsb1
          - 8/3*mmsb2
          )

       + lnMmu*den(mmgl - mmsb1,1) * (
          - 58/3*mmsb1
          + 8/3*mmsb2
          )

       + lnMmu*den(mmgl - mmsb1,2)*sqr(cs2b) * (
          - 16/9*mmsb1*mmsb2
          + 16/9*sqr(mmsb1)
          )

       + lnMmu*den(mmgl - mmsb1,2) * (
          + 16/9*mmsb1*mmsb2
          - 8/9*sqr(mmsb1)
          )

       + lnMmu*den(mmgl - mmsb2,1)*sn2b * (
          + 16/9/mb*mmsb1*mgl
          - 8/9/mb*mmsb2*mgl
          )

       + lnMmu*den(mmgl - mmsb2,1)*sn4b*cs2b * (
          - 8/9/mb*mmsb1*mgl
          + 8/9/mb*mmsb2*mgl
          )

       + lnMmu*den(mmgl - mmsb2,1)*sqr(cs2b) * (
          - 8/3*mmsb1
          + 88/9*mmsb2
          )

       + lnMmu*den(mmgl - mmsb2,1) * (
          + 8/3*mmsb1
          - 58/3*mmsb2
          )

       + lnMmu*den(mmgl - mmsb2,2)*sqr(cs2b) * (
          - 16/9*mmsb1*mmsb2
          + 16/9*sqr(mmsb2)
          )

       + lnMmu*den(mmgl - mmsb2,2) * (
          + 16/9*mmsb1*mmsb2
          - 8/9*sqr(mmsb2)
          )

       + lnMmu * (
          - 932/9
          )

       + den(mmgl - mmsb1,1)*sn2b * (
          + 88/9/mb*mmsb1*mgl*zt2
          + 676/9/mb*mmsb1*mgl
          + 4/3/mb*mmsb2*mgl*zt2
          + 56/9/mb*mmsb2*mgl
          + 4/3/mb*mmst1*mgl*zt2
          + 8/mb*mmst1*mgl
          + 4/3/mb*mmst2*mgl*zt2
          + 8/mb*mmst2*mgl
          + 32/3/mb*mmsusy*mgl*zt2
          + 64/mb*mmsusy*mgl
          )

       + den(mmgl - mmsb1,1)*sn4b*cs2b * (
          - 8/9/mb*mmsb1*mgl
          + 8/9/mb*mmsb2*mgl
          )

       + den(mmgl - mmsb1,1)*sqr(cs2b) * (
          + 41/9*mmsb1*zt2
          + 439/9*mmsb1
          - 8/3*mmsb2
          )

       + den(mmgl - mmsb1,1)*den(mmsb1 - mmsb2,1)*sqr(cs2b) * (
          - 20*sqr(mmsb1)*zt2
          - 140*sqr(mmsb1)
          )

       + den(mmgl - mmsb1,1)*den(mmsb1 - mmsb2,1) * (
          + 20/3*sqr(mmsb1)*zt2
          + 428/9*sqr(mmsb1)
          )

       + den(mmgl - mmsb1,1) * (
          + 50/9*mmsb1*zt2
          + 61/9*mmsb1
          - mmsb2*zt2
          - 10/3*mmsb2
          - mmst1*zt2
          - 6*mmst1
          - mmst2*zt2
          - 6*mmst2
          - 8*mmsusy*zt2
          - 48*mmsusy
          )

       + den(mmgl - mmsb1,2)*sqr(cs2b) * (
          - 16/9*mmsb1*mmsb2
          + 16/9*sqr(mmsb1)
          )

       + den(mmgl - mmsb1,2) * (
          - 4/3*mmsb1*mmsb2*zt2
          - 56/9*mmsb1*mmsb2
          - 4/3*mmsb1*mmst1*zt2
          - 8*mmsb1*mmst1
          - 4/3*mmsb1*mmst2*zt2
          - 8*mmsb1*mmst2
          - 32/3*mmsb1*mmsusy*zt2
          - 64*mmsb1*mmsusy
          + 44/3*sqr(mmsb1)*zt2
          + 868/9*sqr(mmsb1)
          )

       + den(mmgl - mmsb1,3) * (
          + 16/3*pow(mmsb1,3)*zt2
          + 112/3*pow(mmsb1,3)
          )

       + den(mmgl - mmsb2,1)*sn2b * (
          - 4/3/mb*mmsb1*mgl*zt2
          - 56/9/mb*mmsb1*mgl
          - 88/9/mb*mmsb2*mgl*zt2
          - 676/9/mb*mmsb2*mgl
          - 4/3/mb*mmst1*mgl*zt2
          - 8/mb*mmst1*mgl
          - 4/3/mb*mmst2*mgl*zt2
          - 8/mb*mmst2*mgl
          - 32/3/mb*mmsusy*mgl*zt2
          - 64/mb*mmsusy*mgl
          )

       + den(mmgl - mmsb2,1)*sn4b*cs2b * (
          - 8/9/mb*mmsb1*mgl
          + 8/9/mb*mmsb2*mgl
          )

       + den(mmgl - mmsb2,1)*sqr(cs2b) * (
          - 20*mmsb1*zt2
          - 428/3*mmsb1
          - 139/9*mmsb2*zt2
          - 821/9*mmsb2
          )

       + den(mmgl - mmsb2,1)*den(mmsb1 - mmsb2,1)*sqr(cs2b) * (
          + 20*sqr(mmsb1)*zt2
          + 140*sqr(mmsb1)
          )

       + den(mmgl - mmsb2,1)*den(mmsb1 - mmsb2,1) * (
          - 20/3*sqr(mmsb1)*zt2
          - 428/9*sqr(mmsb1)
          )

       + den(mmgl - mmsb2,1) * (
          + 17/3*mmsb1*zt2
          + 398/9*mmsb1
          + 110/9*mmsb2*zt2
          + 163/3*mmsb2
          - mmst1*zt2
          - 6*mmst1
          - mmst2*zt2
          - 6*mmst2
          - 8*mmsusy*zt2
          - 48*mmsusy
          )

       + den(mmgl - mmsb2,2)*sqr(cs2b) * (
          - 16/9*mmsb1*mmsb2
          + 16/9*sqr(mmsb2)
          )

       + den(mmgl - mmsb2,2) * (
          - 4/3*mmsb1*mmsb2*zt2
          - 56/9*mmsb1*mmsb2
          - 4/3*mmsb2*mmst1*zt2
          - 8*mmsb2*mmst1
          - 4/3*mmsb2*mmst2*zt2
          - 8*mmsb2*mmst2
          - 32/3*mmsb2*mmsusy*zt2
          - 64*mmsb2*mmsusy
          + 44/3*sqr(mmsb2)*zt2
          + 868/9*sqr(mmsb2)
          )

       + den(mmgl - mmsb2,3) * (
          + 16/3*pow(mmsb2,3)*zt2
          + 112/3*pow(mmsb2,3)
          )

       - 77/2
          + 8/9*zt2
         ;

  return resmb;
  }*/

#endif
