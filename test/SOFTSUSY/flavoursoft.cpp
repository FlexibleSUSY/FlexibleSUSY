/** \file flavoursoft.cpp
   Project:     SOFTSUSY 
   Author:      Ben Allanach 
   Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   Webpage:     http://hepforge.cedar.ac.uk/softsusy/

*/

#include "flavoursoft.h"
#include <utility>

namespace softsusy {

extern template class Softsusy<SoftParsMssm>;

extern double sw2, gnuL, guL, gdL, geL, guR, gdR, geR, yuL, yuR, ydL, 
  ydR, yeL, yeR, ynuL;

const FlavourMssmSoftsusy & 
FlavourMssmSoftsusy::operator=(const FlavourMssmSoftsusy & s) {
    if (this == &s) return *this;
    setSoftsusy(s.displaySoftsusy());
    fv = s.fv;
    theta12 = s.theta12;
    theta23 = s.theta23;
    theta13 = s.theta13;
    deltaCkm = s.deltaCkm;
    thetaB12 = s.thetaB12;
    thetaB23 = s.thetaB23;
    thetaB13 = s.thetaB13;
    mNuE   = s.mNuE;
    mNuMu  = s.mNuMu;
    mNuTau = s.mNuTau;
    md2GeV = s.md2GeV;
    mu2GeV = s.mu2GeV;
    ms2GeV = s.ms2GeV;
    mcMc   = s.mcMc;
    mePole = s.mePole;
    mmuPole = s.mmuPole;
    return *this;
}

void convertFromWolfenstein(double lambdaW, double A, double rhobar, 
			    double etabar, double & sin12, double & sin23, 
			    double & sin13) {
  sin12 = lambdaW;
  sin23 = A * sqr(lambdaW);
  Complex x; ///< sin theta13 e^i delta
  Complex i(0., 1.);
  double a2LambdaW4 = sqr(A) * sqr(sqr(lambdaW));
  x = A * lambdaW * sqr(lambdaW) * (rhobar + i * etabar) * 
    sqrt(1.0 - a2LambdaW4) / sqrt(1.0 - sqr(lambdaW)) /
    (1.0 - a2LambdaW4 * (rhobar + i * etabar));
  sin13 = x.mod();
}

int positionOfSym(int i, int j) {
  if (i == 1 && j == 1) return 1;
  else if (i == 2 && j == 2) return 4;
  else if (i == 3 && j == 3) return 6;
  else if ((i == 1 && j == 2) || (i == 2 && j == 1)) return 2;
  else if ((i == 1 && j == 3) || (i == 3 && j == 1)) return 3;
  else if ((i == 2 && j == 3) || (i == 3 && j == 2)) return 5;
  else {
    std::ostringstream ii;
    ii << " # WARNING: can't convert position in 3 by 3 matrix" 
       << i << j << endl; 
    throw ii.str();
    return 0;
  }
}

void flavourBcs(MssmSoftsusy & m, 
		 const DoubleVector & inputParameters) {
  int i, count = 1;
  for (i=1; i<=3; i++) {
    m.setGauginoMass(i, inputParameters.display(count));
    count++;
  }
 
  /// Find the Yukawa rotation matrices such that the super-CKM basis may be
  /// defined
  DoubleMatrix Ud(3, 3), Vd(3, 3), Uu(3, 3), Vu(3, 3);
  DoubleVector yu(3), yd(3);
  //  m.displayYukawaMatrix(YU).diagonalise(Uut, Vu, yu);
  m.displayYukawaMatrix(YU).diagonalise(Vu, Uu, yu);
  m.displayYukawaMatrix(YD).diagonalise(Vd, Ud, yd);
  ckmNormalise(Vu, Vd, Uu, Ud); 

  DoubleMatrix Ue(3, 3), Ve(3, 3);
  DoubleVector ye(3);
  m.displayYukawaMatrix(YE).diagonalise(Ve, Ue, ye);

  DoubleMatrix msq2in(3, 3), msu2in(3, 3), msd2in(3, 3), 
    msl2in(3, 3), mse2in(3, 3);

  int j; for (i=1; i<=3; i++)
    for (j=i; j<=3; j++) {
      /// NB switched i and j: may be a problem in general for complex masses
      msq2in(i, j) = inputParameters.display(count);
      count++;
    }

   for (i=1; i<=3; i++)
    for (j=i; j<=3; j++) {
      msu2in(i, j) = inputParameters.display(count);
      count++;
    }
   for (i=1; i<=3; i++)
    for (j=i; j<=3; j++) {
      msd2in(i, j) = inputParameters.display(count);
      count++;
    }
   for (i=1; i<=3; i++)
    for (j=i; j<=3; j++) {
      msl2in(i, j) = inputParameters.display(count);
      count++;
    }
   for (i=1; i<=3; i++)
    for (j=i; j<=3; j++) {
      mse2in(i, j) = inputParameters.display(count);
      count++;
    }

   msl2in.symmetrise();
   msq2in.symmetrise();
   msu2in.symmetrise();
   msd2in.symmetrise();
   mse2in.symmetrise();
   
  m.setSoftMassMatrix(mLl, Ve * msl2in * Ve.transpose());
  m.setSoftMassMatrix(mEr, Ue * mse2in.transpose() * Ue.transpose());
  m.setSoftMassMatrix(mQl, Vd * msq2in * Vd.transpose());
  m.setSoftMassMatrix(mUr, Uu * msu2in.transpose() * Uu.transpose());
  m.setSoftMassMatrix(mDr, Ud * msd2in.transpose() * Ud.transpose());

  /// This sets the trilinears with any universal ones set
  m.universalTrilinears(inputParameters.display(62));

  /// This overwrites (in the CKM basis) any of trilinears we need to
  /// First, must transform to SCKM basis
  DoubleMatrix sckm(Uu.transpose() * m.displayTrilinear(UA).transpose() * 
		    Vu);
  DoubleMatrix tuin(3, 3), tdin(3, 3), tein(3, 3);
  int count2 = 0;
  for (i=1; i<=3; i++)
    for (j=1; j<=3; j++) {
      if (slha2setTrilinear[count2]) 
	sckm(i, j) = inputParameters.display(count);
      count++; count2++;
    }
  /// transform back to interaction basis
  m.setTrilinearMatrix(UA, Vu * sckm.transpose() * Uu.transpose());
  
  sckm = Ud.transpose() * m.displayTrilinear(DA).transpose() * Vd; 
  for (i=1; i<=3; i++)
    for (j=1; j<=3; j++) {
      if (slha2setTrilinear[count2]) 
	sckm(i, j) = inputParameters.display(count);
      count++; count2++;
    }
  m.setTrilinearMatrix(DA, Vd * sckm.transpose() * Ud.transpose());

  sckm = Ue.transpose() * m.displayTrilinear(EA).transpose() * 
    Ve; 
  for (i=1; i<=3; i++)
    for (j=1; j<=3; j++) {
      if (slha2setTrilinear[count2]) 
	sckm(i, j) = inputParameters.display(count);
      count++; count2++;
    }
  m.setTrilinearMatrix(EA, Ve * sckm.transpose() * Ue.transpose());  

  if (!m.displayAltEwsb()) {
    m.setMh1Squared(inputParameters.display(63));
    m.setMh2Squared(inputParameters.display(64));
  }
}

const ComplexMatrix FlavourMssmSoftsusy::displayComplexCkm() const {
  double d   = displayDelta();
  Complex eID(cos(d), sin(d));
  double s12 = sin(displayTheta12());
  double s13 = sin(displayTheta13());
  double s23 = sin(displayTheta23());
  
  double c12 = cos(displayTheta12());
  double c13 = cos(displayTheta13());
  double c23 = cos(displayTheta23());

  ComplexMatrix ckmMatrix(3, 3);
  ckmMatrix(1, 1) = c12 * c13;      
  ckmMatrix(1, 2) = s12 * c13; 
  ckmMatrix(1, 3) = (s13 / eID);            
  ckmMatrix(2, 1) = (-s12 * c23 - c12 * s23 * s13 * eID);
  ckmMatrix(2, 2) = (c12 * c23 - s12 * s23 * s13 * eID);
  ckmMatrix(2, 3) = s23 * c13; 
  ckmMatrix(3, 1) = (s12 * s23 - c12 * c23 * s13 * eID); 
  ckmMatrix(3, 2) = (-c12 * s23 - s12 * c23 * s13 * eID); 
  ckmMatrix(3, 3) = c23 * c13; 

  return ckmMatrix;
}

const DoubleMatrix FlavourMssmSoftsusy::displayCkm() const {
  return display3x3RealMixing(displayTheta12(), displayTheta13(), 
			      displayTheta23(), displayDelta());
}

const DoubleMatrix FlavourMssmSoftsusy::displayMns() const {
  double d = 0.;
  return display3x3RealMixing(displayThetaB12(), displayThetaB13(), 
			      displayThetaB23(), d);
}

std::ostream & operator <<(std::ostream & left, const FlavourMssmSoftsusy & m) {
  const std::string HR = "----------------------------------------------------------";
  left << m.displayMssmSoft();
  left << "Flavour violating parameters:\n";
  left << "CKM data: 12=" << m.displayTheta12() 
       << " 23=" << m.displayTheta23() 
       << " 13=" << m.displayTheta13() << endl;
  left << "md(2): " << m.displayMd2GeV();
  left << " mu(2): " << m.displayMu2GeV();
  left << " ms(2): " << m.displayMs2GeV();
  left << " mc(mc): " << m.displayMcMc() << endl;
  left << "mePole: " << m.displayPoleMe() 
       << " mmuPole: " << m.displayPoleMmu() << endl;
  left << "m_nue: " << m.displayMnuE() << "  m_numu: " << m.displayMnuMu()
       << " m_nutau: " << m.displayMnuTau() << endl;
  left << "MNS data: 12=" << m.displayThetaB12() 
       << " 23=" << m.displayThetaB23() 
       << " 13=" << m.displayThetaB13() << endl;
  left << "uSqMix" << m.displayFlavourPhysical().uSqMix;
  left << "dSqMix" << m.displayFlavourPhysical().dSqMix;
  left << "eSqMix" << m.displayFlavourPhysical().eSqMix;
  left << "nuSqMix" << m.displayFlavourPhysical().nuSqMix;

  left << HR;
  return left;
}

void FlavourMssmSoftsusy::modselSLHA(std::ostream & out, const char model []) {
  MssmSoftsusy::modselSLHA(out, model);
  out << "     6    " << 1 << "   # flavour violating MSSM\n"; 
}

void FlavourMssmSoftsusy::sminputsSLHA(std::ostream & out) {
  MssmSoftsusy::sminputsSLHA(out);
  /*  out << "     8   "; printRow(out, displayMnuTau()); 
  out << "   # Mnu3(pole)\n";
  out << "    14   "; printRow(out, displayMnuMu()); 
  out << "   # Mnu2(pole)\n";
  out << "    15   "; printRow(out, displayMnuE()); 
  out << "   # Mnu1(pole)\n";*/

  out << "    21   ";  printRow(out, displayMd2GeV()); 
  out << "   # Mdown(2 GeV) MSbar\n";
  out << "    22   ";  printRow(out, displayMu2GeV()); 
  out << "   # Mup(2 GeV) MSbar\n";
  out << "    23   ";  printRow(out, displayMs2GeV()); 
  out << "   # Mstrange(2 GeV) MSbar\n";
  //  d.runto(maximum(d.displayMass(mCharm), 1.0));
  out << "    24   ";  printRow(out, displayMcMc()); 
  out << "   # Mcharm(Mcharm) MSbar\n";
  out << "    11   "; printRow(out, displayPoleMe()); 
  out << "   # Me(pole)\n";
  out << "    13   "; printRow(out, displayPoleMmu()); 
  out << "   # Mmu(pole)\n";
 
}

void FlavourMssmSoftsusy::sleptonsSLHA(std::ostream & out) {
  DoubleVector m(6); int count = 0, i, j;
  for (i=1; i<=2; i++)
    for (j=1; j<=3; j++) {
      count++;
      m(count) = displayPhys().me(i, j);
    }
  m = m.sort();
  out << "   1000011    "; printRow(out, m(1)); out << "   # ~e_1\n";
  out << "   1000013    "; printRow(out, m(2)); out << "   # ~e_2\n";
  out << "   1000015    "; printRow(out, m(3)); out << "   # ~e_3\n";
  out << "   2000011    "; printRow(out, m(4)); out << "   # ~e_4\n";
  out << "   2000013    "; printRow(out, m(5)); out << "   # ~e_5\n";
  out << "   2000015    "; printRow(out, m(6)); out << "   # ~e_6\n";
  
  m.setEnd(3); 
  for (i=1; i<=3; i++) m(i) = displayPhys().msnu(i);
  m = m.sort();
  out << "   1000012    "; printRow(out, m(1)); out << "   # ~nu_1\n";
  out << "   1000014    "; printRow(out, m(2)); out << "   # ~nu_2\n";
  out << "   1000016    "; printRow(out, m(3)); out << "   # ~nu_3\n";
}
 
void FlavourMssmSoftsusy::sfermionsSLHA(std::ostream & out) {
  /// SLHA2 flavour sector defined in increasing order of mass
  sPhysical s(displayPhys());
  int i, j; 

  sleptonsSLHA(out);
  DoubleVector m(6);
  int count = 0;
  for (i=1; i<=2; i++)
    for (j=1; j<=3; j++) {
      count++;
      m(count) = s.md(i, j);
    }
  m = m.sort();
  out << "   1000001    "; printRow(out, m(1)); out << "   # ~d_1\n";
  out << "   1000003    "; printRow(out, m(2)); out << "   # ~d_2\n";
  out << "   1000005    "; printRow(out, m(3)); out << "   # ~d_3\n";
  out << "   2000001    "; printRow(out, m(4)); out << "   # ~d_4\n";
  out << "   2000003    "; printRow(out, m(5)); out << "   # ~d_5\n";
  out << "   2000005    "; printRow(out, m(6)); out << "   # ~d_6\n";
      
  count = 0;
  for (i=1; i<=2; i++)
    for (j=1; j<=3; j++) {
      count++;
      m(count) = s.mu(i, j);
    }
  m = m.sort();
  out << "   1000002    "; printRow(out, m(1)); out << "   # ~u_1\n";
  out << "   1000004    "; printRow(out, m(2)); out << "   # ~u_2\n";
  out << "   1000006    "; printRow(out, m(3)); out << "   # ~u_3\n";
  out << "   2000002    "; printRow(out, m(4)); out << "   # ~u_4\n";
  out << "   2000004    "; printRow(out, m(5)); out << "   # ~u_5\n";
  out << "   2000006    "; printRow(out, m(6)); out << "   # ~u_6\n";      
}

void FlavourMssmSoftsusy::selmixSLHA(std::ostream & out) {
  out << "Block SELMIX  # super MNS slepton mass^2 matrix\n";
  for (int i=1; i<=6; i++)
    for (int j=1; j<=6; j++) {
      out << "  " << i << "  " << j << "    "; 
      printRow(out, displayFlavourPhysical().eSqMix.display(i, j));
      out << "   # (SELMIX)_{" << i << j << "}" << endl;
    }
}

void FlavourMssmSoftsusy::snumixSLHA(std::ostream & out) {
  out << "Block SNUMIX                # super MNS slepton mass^2 matrix\n";
  for (int i=1; i<=3; i++)
    for (int j=1; j<=3; j++) {
      out << "  " << i << "  " << j << "    "; 
      printRow(out, displayFlavourPhysical().nuSqMix.display(i, j));
      out << "   # (SNUMIX)_{" << i << j << "}" << endl;
    }
}

void FlavourMssmSoftsusy::sfermionmixSLHA(std::ostream & out) {
  out << "Block USQMIX  # super CKM squark mass^2 matrix\n";
  int i, j; for (i=1; i<=6; i++)
    for (j=1; j<=6; j++) {
      out << "  " << i << "  " << j << "    "; 
      printRow(out, displayFlavourPhysical().uSqMix.display(i, j));
      out << "   # (USQMIX)_{" << i << j << "}" << endl;
    }
  
  out << "Block DSQMIX  # super CKM squark mass^2 matrix\n";
  for (i=1; i<=6; i++)
    for (j=1; j<=6; j++) {
      out << "  " << i << "  " << j << "    "; 
      printRow(out, displayFlavourPhysical().dSqMix.display(i, j));
      out << "   # (DSQMIX)_{" << i << j << "}" << endl;
    }

  selmixSLHA(out);
  snumixSLHA(out);
}

void FlavourMssmSoftsusy::vckminSLHA(std::ostream & out) {
  /// these come from 2006 PDG
  double s12 = sin(theta12);
  double s23 = sin(theta23);
  double s13 = sin(theta13);
  
  out << "Block VCKMIN               # input CKM mixing matrix parameters\n";
  double lambdaW = s12; 
  double aCkm;
  if (s23 == 0.) aCkm = 0.; 
    else aCkm = s23 / sqr(lambdaW);
  double c = sqrt((1.0 - sqr(s23)) / (1.0 - sqr(lambdaW)));
  Complex eid(cos(deltaCkm), sin(deltaCkm));
  Complex r;
  if (s13 == 0.) r = Complex(0., 0.);
  else r = (s13 * eid / 
	    (c * aCkm * lambdaW * sqr(lambdaW) + s13 * eid * sqr(s23)));
  
  out << "     1    " << lambdaW << "   # lambda_W\n";
  out << "     2    " << aCkm << "   # A\n";
  out << "     3    " << r.real() << "   # rhobar\n";  
  out << "     4    " << r.imag() 
      << "   # etabar (no phases used in SOFTSUSY yet though)\n";    
}

void FlavourMssmSoftsusy::extparSLHA(std::ostream & out, 
				     const DoubleVector & pars, 
				     bool ewsbBCscale) {
    out << "Block EXTPAR              "
	<< " # non-universal SUSY breaking parameters\n";
    if (ewsbBCscale) 
      out << "     0    -1.00000000e+00  # Set MX=MSUSY\n";
    else {
      out << "     0    "; printRow(out, displayMxBC()); 
      cout << "  # MX scale\n";
    }

    int count = 0, i;
    for (i=1; i<=3; i++) {
      count++;
      out << "     " << i << "   ";
      printRow(out, pars.display(count)); 
      out << "   # M_" << count << "(MX)" << endl;      
    }
    
    if (!displayAltEwsb()) {
      out << "     21  "; printRow(out, pars.display(63)) ; 
      out << "   # mHd^2(MX)" << endl;    
      out << "     22  "; printRow(out, pars.display(64)) ; 
      out << "   # mHu^2(MX)" << endl;    
    } else {
      out << "     23   "; printRow(out, displayMuCond()) ; 
      out << "  # mu(MX)" << endl;    
      out << "     26   "; printRow(out, displayMaCond()) ; 
      out << "  # mA(pole)" << endl;    
    }

    out << "Block MSQ2IN               # GUT-scale input for DRbar mass^2 parameter\n";
    int j; for (i=1; i<=3; i++)
      for (j=i; j<=3; j++) {
	count++;
	out << "   " << i << " " << j << "   ";
	printRow(out, pars.display(count));
	out << "   # input (mhat_Q^2)_{" << i << "," << j << "} in SCKM basis";
	out << endl;
      }
    out << "Block MSU2IN               # GUT-scale input for DRbar mass^2 parameter\n";
    for (i=1; i<=3; i++)
      for (j=i; j<=3; j++) {
	count++;
	out << "   " << i << " " << j << "   ";
	printRow(out, pars.display(count));
	out << "   # input (mhat_u^2)_{" << i << "," << j << "} in SCKM basis";
	out << endl;
      }
    out << "Block MSD2IN               # GUT-scale input for DRbar mass^2 parameter\n";
    for (i=1; i<=3; i++)
      for (j=i; j<=3; j++) {
	count++;
	out << "   " << i << " " << j << "   ";
	printRow(out, pars.display(count));
	out << "   # input (mhat_d^2)_{" << i << "," << j << "} in SCKM basis";
	out << endl;
      }
    out << "Block MSL2IN               # GUT-scale input for DRbar mass^2 parameter\n";
    for (i=1; i<=3; i++)
      for (j=i; j<=3; j++) {
	count++;
	out << "   " << i << " " << j << "   ";
	printRow(out, pars.display(count));
	out << "   # input (mhat_L^2)_{" << i << "," << j << "} in SCKM basis";
	out << endl;
      }
    out << "Block MSE2IN               # GUT-scale input for DRbar mass^2 parameter\n";
    for (i=1; i<=3; i++)
      for (j=i; j<=3; j++) {
	count++;
	out << "   " << i << " " << j << "   ";
	printRow(out, pars.display(count));
	out << "   # input (mhat_L^2)_{" << i << "," << j << "} in SCKM basis";
	out << endl;
      }

    bool anyTus = false;       int count2 = 0;
    for (i=0; i<9; i++) if (slha2setTrilinear[i]) anyTus = true;
    if (anyTus) 
      out << "Block TUIN                 # GUT-scale input for DRbar mass parameter\n";    
    for (i=1; i<=3; i++)
      for (j=1; j<=3; j++) {
	count++; 
	if (slha2setTrilinear[count2]) {	  
	  out << "   " << i << " " << j << "   ";
	  printRow(out, pars.display(count));
	out << "   # input (That_U)_{" << i << "," << j << "} in SCKM basis";
	  out << endl;
	}
	count2++;
      }
  
    anyTus = false; 
    for (i=9; i<18; i++) if (slha2setTrilinear[i]) anyTus = true;
    if (anyTus) 
      out << "Block TDIN                 # GUT-scale input for DRbar mass parameter\n";    
    for (i=1; i<=3; i++)
      for (j=1; j<=3; j++) {
	count++;
	if (slha2setTrilinear[count2]) {	  
	  out << "   " << i << " " << j << "   ";
	  printRow(out, pars.display(count));
	  out << "   # input (That_D)_{" << i << "," << j << "} in SCKM basis";
	  out << endl;
	}
	count2++;
      }

    anyTus = false; 
    for (i=19; i<27; i++) if (slha2setTrilinear[i]) anyTus = true;
    if (anyTus) 
    out << "Block TEIN                   # GUT-scale input for DRbar mass parameter\n";    
    for (i=1; i<=3; i++)
      for (j=1; j<=3; j++) {
	count++;
	if (slha2setTrilinear[count2]) {	  
	  out << "   " << i << " " << j << "   ";
	  printRow(out, pars.display(count));
	  out << "   # input (That_E)_{" << i << "," << j << "} in SCKM basis";
	  out << endl;
	}
	count2++;
    }
    if (fabs(pars.display(62)) > EPSTOL) 
      out << "# CAUTION: a universal A0=" << pars.display(62) 
	  << " was set in addition.\n";
}

void FlavourMssmSoftsusy::yukawasSLHA(std::ostream & out) {
  /// Find the Yukawa rotation matrices such that the super-CKM basis may be
  /// defined
  DoubleMatrix Ud(3, 3), Vd(3, 3), Uu(3, 3), Vu(3, 3);
  DoubleVector yu(3), yd(3);
  displayYukawaMatrix(YU).diagonalise(Vu, Uu, yu);
  displayYukawaMatrix(YD).diagonalise(Vd, Ud, yd);
  ckmNormalise(Vu, Vd, Uu, Ud); 
  
  DoubleMatrix Ue(3, 3), Ve(3, 3);
  DoubleVector ye(3);
  displayYukawaMatrix(YE).diagonalise(Ve, Ue, ye);
  
  out << "Block yu Q= " << displayMu() << "   # diagonal Up Yukawa matrix\n";
  int i; for (i=1; i<=3; i++) {
    out << "  " << i << "  " << i << "     " << yu(i)
	<< "    # YU_{" << i << i << "}(Q)MSSM DRbar" << endl;
  }
  
  out << "Block yd Q= " << displayMu() << "   # diagonal down Yukawa matrix\n";
  for (i=1; i<=3; i++) {
    out << "  " << i << "  " << i << "     " << yd(i)
	<< "    # YD_{" << i << i << "}(Q)MSSM DRbar" << endl;
  }
  
  out << "Block ye Q= " << displayMu() 
      << "   # diagonal lepton Yukawa matrix\n";
  for (i=1; i<=3; i++) {
    out << "  " << i << "  " << i << "     " << ye(i)
	<< "    # YE_{" << i << i << "}(Q)MSSM DRbar" << endl;
  }    
}

void FlavourMssmSoftsusy::msoftSLHA(std::ostream & out) {
    /// Find the Yukawa rotation matrices such that the super-CKM basis may be
    /// defined
    DoubleMatrix Ud(3, 3), Vd(3, 3), Uu(3, 3), Vu(3, 3);
    DoubleVector yu(3), yd(3);
    displayYukawaMatrix(YU).diagonalise(Vu, Uu, yu);
    displayYukawaMatrix(YD).diagonalise(Vd, Ud, yd);
    ckmNormalise(Vu, Vd, Uu, Ud); 
  
    DoubleMatrix Ue(3, 3), Ve(3, 3);
    DoubleVector ye(3);
    displayYukawaMatrix(YE).diagonalise(Ve, Ue, ye);
    
    DoubleMatrix tuSckm(3, 3), muSqSckm(3, 3);
    DoubleMatrix mqSqSckm(3, 3), tdSckm(3, 3), mdSqSckm(3, 3);
    mqSqSckm  = Vd.transpose() * displaySoftMassSquared(mQl) * Vd;
    tdSckm    = Ud.transpose() * displayTrilinear(DA).transpose() * Vd;
    mdSqSckm  = Ud.transpose() * displaySoftMassSquared(mDr).transpose() * 
      Ud; 

    tuSckm    = Uu.transpose() * displayTrilinear(UA).transpose() * Vu;
    muSqSckm  = Uu.transpose() * displaySoftMassSquared(mUr).transpose() * 
      Uu; 
    
    DoubleMatrix mlSqSckm(3, 3), teSckm(3, 3), meSqSckm(3, 3);
    mlSqSckm  = Ve.transpose() * displaySoftMassSquared(mLl) * Ve;
    teSckm    = Ue.transpose() * displayTrilinear(EA).transpose() * Ve;
    meSqSckm  = Ue.transpose() * displaySoftMassSquared(mEr).transpose() * 
      Ue; 

    /// There are often high degrees of cancellation in these matrices,
    /// meaning that at the 5th significant figure, one can have
    /// asymmetries. In order to avoid problems with users' codes, these are
    /// just symmetrised.
    mdSqSckm.symmetrise();
    muSqSckm.symmetrise();
    meSqSckm.symmetrise();   
    mqSqSckm.symmetrise();
    mlSqSckm.symmetrise();

    out << "Block msq2 Q= " << displayMu() << " # super CKM squark mass^2 matrix"
	<< " - DRbar" << endl;
    int i, j; for (i=1; i<=3; i++)
      for (j=1; j<=3; j++) {
	out << "  " << i << "  " << j << "    "; 
	printRow(out, mqSqSckm(i, j));
	out << "    # (m_Q^2)_{" << i << j << "}" << endl;
      }
    
    out << "Block msl2 Q= " << displayMu() << " # super MNS slepton mass^2 matrix"
	<< " - DRbar" << endl;
    for (i=1; i<=3; i++)
      for (j=1; j<=3; j++) {
	out << "  " << i << "  " << j << "    "; 
	printRow(out, mlSqSckm(i, j));
	out << "    # (m_L^2)_{" << i << j << "}" << endl;
      }
    
    out << "Block msd2 Q= " << displayMu() << " # super CKM squark mass^2 matrix"
	<< " - DRbar" << endl;
    for (i=1; i<=3; i++)
      for (j=1; j<=3; j++) {
	out << "  " << i << "  " << j << "    "; 
	printRow(out, mdSqSckm(i, j));
	out << "    # (m_d^2)_{" << i << j << "}" << endl;
      }
    
    out << "Block msu2 Q= " << displayMu() << " # super CKM squark mass^2 matrix"
	<< " - DRbar" << endl;
    for (i=1; i<=3; i++)
      for (j=1; j<=3; j++) {
	out << "  " << i << "  " << j << "    "; 
	printRow(out, muSqSckm(i, j));
	out << "    # (m_u^2)_{" << i << j << "}" << endl;
      }
    
    out << "Block mse2 Q= " << displayMu() << " # super MNS slepton mass^2 matrix"
	<< " - DRbar" << endl;
    for (i=1; i<=3; i++)
      for (j=1; j<=3; j++) {
	out << "  " << i << "  " << j << "    "; 
	printRow(out, meSqSckm(i, j));
	out << "    # (m_e^2)_{" << i << j << "}" << endl;
      }
    
    out << "Block tu Q= " << displayMu() << "   # super CKM trilinear matrix"
	<< " - DRbar" << endl;
    for (i=1; i<=3; i++)
      for (j=1; j<=3; j++) {
	out << "  " << i << "  " << j << "    "; 
	printRow(out, tuSckm(i, j));
	out << "    # (T_u)_{" << i << j << "}" << endl;
      }
    
    out << "Block td Q= " << displayMu() << "   # super CKM trilinear matrix"
	<< " - DRbar" << endl;
    for (i=1; i<=3; i++)
      for (j=1; j<=3; j++) {
	out << "  " << i << "  " << j << "    "; 
	printRow(out, tdSckm(i, j));
	out << "    # (T_d)_{" << i << j << "}" << endl;
      }
    
    out << "Block te Q= " << displayMu() << "   # super CKM trilinear matrix"
	<< " - DRbar" << endl;
    for (i=1; i<=3; i++)
      for (j=1; j<=3; j++) {
	out << "  " << i << "  " << j << "    "; 
	printRow(out, teSckm(i, j));
	out << "    # (T_e)_{" << i << j << "}" << endl;
      }
  DoubleMatrix vCkm(Vu.transpose() * Vd);

  out << "Block VCKM Q= " << displayMu() << " # DRbar CKM mixing matrix\n";
  for (i=1; i<=3; i++)
    for (j=1; j<=3; j++) {
      out << "  " << i << "  " << j << "    "; printRow(out, vCkm(i, j));
      out << "    # CKM_{" << i << j << "} matrix element" << endl;
    }
    
    out << "Block msoft Q= " << displayMu() 
	 << " # MSSM DRbar SUSY breaking parameters\n"; 

    for (i=1; i<=3; i++) {
      out << "     " << i << "    "; 
      printRow(out, displayGaugino(i)); 
      out << "     # M_" << i << "(Q)" << endl;      
    }

    out << "    21    "; printRow(out, displayMh1Squared()); 
    out << "     # mH1^2(Q)" << endl;    
    out << "    22    "; printRow(out, displayMh2Squared()); 
    out << "     # mH2^2(Q)" << endl;    
}

void FlavourMssmSoftsusy::lesHouchesAccordOutput(std::ostream & out, 
						 const char model[], 
						 const DoubleVector & pars, 
						 int sgnMu, double tanb, 
						 double qMax, 
						 int numPoints, 
						 bool ewsbBCscale) {

  int nn = out.precision();
  headerSLHA(out);
  spinfoSLHA(out);
  modselSLHA(out, model);
  sminputsSLHA(out);
  vckminSLHA(out);
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

void FlavourMssmSoftsusy::minparSLHA(std::ostream & out, const char model [], 
				     const DoubleVector & pars, double tanb, 
				     int sgnMu, 
				     bool ewsbBCscale) {
  /// For universal models, users still want to know MX and it has to be
  /// specially printed out as EXTPAR 0
  bool printMX = false;

  out << "Block MINPAR               # SUSY breaking input parameters\n";
  out << "     3   "; printRow(out, tanb)            ; out << "   # tanb" << endl;
  if (!displayAltEwsb()) {
    out << "     4   "; 
    printRow(out, double(sgnMu)); 
    out << "   # sign(mu)"<< endl;
  }
  if (!strcmp(model, "sugra")) {
    out << "     1   "; printRow(out, sqrt(pars.display(4))); out << "   # m0" << endl;
    out << "     2   "; printRow(out, pars.display(2)) ; out << "   # m12" << endl;
    out << "     5   "; printRow(out, pars.display(62)) ; out << "   # A0" << endl;
    printMX = true;
  }
  else if (!strcmp(model, "gmsb")) {
    out << "     1   "; printRow(out, pars.display(3)); out << "   # LAMBDA" << endl;
    out << "     2   "; printRow(out, pars.display(2)) ; out << "   # M_mess" 
	 << endl;
    out << "     5   "; printRow(out, pars.display(1)) ; out << "   # N5" << endl;
    out << "     6   "; printRow(out, pars.display(4)) ; out << "   # cgrav" << endl;
  }
  else if (!strcmp(model, "splitgmsb")) {
    out << "     1   "; printRow(out, pars.display(2)); out << "   # LAMBDAL" << endl;
    out << "     2   "; printRow(out, pars.display(3)) ; out << "   # LAMBDAD" 
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
       << model << " - cannot do at present in flavour violating version\n";
    throw ii.str();
  }  
  if (printMX) {
  out << "Block EXTPAR               # scale of SUSY breaking BCs\n";
  out << "     0   "; printRow(out, displayMxBC()); out << "   # MX scale\n";
  }
}
 


void FlavourMssmSoftsusy::doUpSquarks(double /* mt */, double /* pizztMS */,
				      double sinthDRbarMS, int /* accuracy */) {
  /// now need to re-do 6 by 6 scalar mass matrices
  DoubleMatrix massUsq(6, 6);

  /// super CKM basis mass squared matrices
  DoubleMatrix Ud(3, 3), Vd(3, 3), Uu(3, 3), Vu(3, 3);
  DoubleVector yu(3), yd(3);
  displayYukawaMatrix(YU).diagonalise(Vu, Uu, yu);
  displayYukawaMatrix(YD).diagonalise(Vd, Ud, yd);
  ckmNormalise(Vu, Vd, Uu, Ud); 

  DoubleMatrix mqSqSckm(3, 3), tuSckm(3, 3), muSqSckm(3, 3);
  mqSqSckm  = Vu.transpose() * displaySoftMassSquared(mQl) * Vu;
  tuSckm    = Uu.transpose() * displayTrilinear(UA).transpose() * Vu;
  muSqSckm  = Uu.transpose() * displaySoftMassSquared(mUr).transpose() * 
    Uu; 

  DoubleVector mu(3), md(3);

  double susyMu = displaySusyMu();
  double vev = displayHvev();
  double tanb = displayTanb(), beta = atan(tanb);
  double sinb = sin(beta), cosb = cos(beta);
  /// the following contributions are only calculated at tree-level:
  int i, j;   for (i=1; i<=3; i++) {
    mu(i) = yu(i) * vev / root2 * sinb;
    md(i) = yd(i) * vev / root2 * cosb;

    for (j=i; j<=3; j++) {
      massUsq(i, j) = mqSqSckm(i, j);
      if (i == j) {
	massUsq(i, j) = massUsq(i, j) + sqr(mu(i)) + 
	  cos(2.0 * beta) * sqr(displayMzRun()) * 
	  (0.5 - 2. / 3. * sqr(sinthDRbarMS));
      }
    }
  
    for (j=4; j<=6; j++) {
      massUsq(i, j) = vev * sinb / root2 * tuSckm(j-3, i);
      if (i == j-3) massUsq(i, j) = massUsq(i, j) - susyMu * mu(i) / tanb;
    }
  }

  for (i=4; i<=6; i++)
    for (j=i; j<=6; j++) {
      massUsq(i, j) = muSqSckm(i - 3, j - 3);
      if (i == j) massUsq(i, j) = massUsq(i, j) + 
		    cos(2.0 * beta) * sqr(displayMzRun()) * 
		    sqr(sinthDRbarMS) * 2.0 / 3.0 + sqr(mu(i - 3)); 
    }

  /// add loop corrections to the family-universal pieces
  DoubleMatrix temp(2, 2);
  double p = minimum(displayDrBarPars().mu(1, 3), displayDrBarPars().mu(2, 3));
  int family; for (family=1; family<=3; family++) {
    temp(1, 1) = massUsq(family, family);
    temp(1, 2) = massUsq(family, family + 3);
    temp(2, 2) = massUsq(family + 3, family + 3);

    DoubleMatrix aa(temp);
    if (family < 3) addSupCorrection(temp, family);
    else {
      /// DEBUG: take the real part of the stop mass matrix for the time being
      addStopCorrection(p, aa, mu(family));
      for (i=1; i<=2; i++) for (j=1; j<=2; j++) temp(i, j) = aa(i, j);
    }
    massUsq(family, family) = temp(1, 1);
    massUsq(family, family + 3) = temp(1, 2);
    massUsq(family + 3, family + 3) = temp(2, 2);
  }
  
  massUsq.symmetrise();

  DoubleMatrix uSqMixT(6, 6);
  DoubleVector uSqMasses(6);
  if (massUsq.diagonaliseSym(uSqMixT, uSqMasses) > EPSTOL) {
    throw("WARNING: up squark flavour diagonalisation bad accuracy\n");
  }

  fv.uSqMix = uSqMixT.transpose();
  if (fv.msU.min(j) < 0.) {
    flagTachyon(sup);     flagTachyon(scharm);     flagTachyon(stop);
    if (PRINTOUT > 1) cout << " up squark(" << j << ") tachyon";
  }
  fv.msU = uSqMasses.apply(ccbSqrt);

  // fill in the flavour conserving parts with the correct masses
  int t1Pos = 0;
  sPhysical s(displayPhys());
  DoubleMatrix a(fv.uSqMix);

  for(i=3; i>=1; i--) {
    a.displayCol(i).apply(fabs).max(j); 
    s.mu(1, i) = fv.msU(j);
    for (int k=1; k<=6; k++) a(j, k) = 0.;
    if (i == 3) t1Pos = j;
    a.displayCol(i+3).apply(fabs).max(j); 
    for (int k=1; k<=6; k++) a(j, k) = 0.;
    s.mu(2, i) = fv.msU(j);
  }

  s.thetat = asin(fv.uSqMix(t1Pos, 6));
  setPhys(s);

}

void FlavourMssmSoftsusy::doDownSquarks(double /* mb */, double /* pizztMS */,
					double sinthDRbarMS, int accuracy, 
					double /* mt */) {
  /// now need to re-do 6 by 6 scalar mass matrices
  DoubleMatrix massDsq(6, 6);

  /// Find the Yukawa rotation matrices such that the super-CKM basis may be
  /// defined
  DoubleMatrix Ud(3, 3), Vd(3, 3), Uu(3, 3), Vu(3, 3);
  DoubleVector yu(3), yd(3);
  displayYukawaMatrix(YU).diagonalise(Vu, Uu, yu);
  displayYukawaMatrix(YD).diagonalise(Vd, Ud, yd);
  ckmNormalise(Vu, Vd, Uu, Ud); 

  DoubleMatrix mqSqSckm(3, 3), tdSckm(3, 3), mdSqSckm(3, 3);
  mqSqSckm  = Vd.transpose() * displaySoftMassSquared(mQl) * Vd;
  tdSckm    = Ud.transpose() * displayTrilinear(DA).transpose() * Vd;
  mdSqSckm  = Ud.transpose() * displaySoftMassSquared(mDr).transpose() * 
    Ud; 

  DoubleVector mu(3), md(3);

  double susyMu = displaySusyMu();
  double vev = displayHvev();
  double tanb = displayTanb(), beta = atan(tanb);
  double sinb = sin(beta), cosb = cos(beta);
  /// the following contributions are only calculated at tree-level:
  int i, j;   for (i=1; i<=3; i++) {
    mu(i) = yu(i) * vev / root2 * sinb;
    md(i) = yd(i) * vev / root2 * cosb;

    for (j=i; j<=3; j++) {
      massDsq(i, j) = mqSqSckm(i, j);
      if (i == j) {
	massDsq(i, j) = massDsq(i, j) + sqr(md(i)) + 
	  cos(2.0 * beta) * sqr(displayMzRun()) * 
	  (-0.5 + 1. / 3. * sqr(sinthDRbarMS));
      }
    }

    for (j=4; j<=6; j++) {
      massDsq(i, j) = vev * cosb / root2 * tdSckm(j-3, i);
      if (i == j-3) massDsq(i, j) = massDsq(i, j) - susyMu * md(i) * tanb;
    }
  }

  for (i=4; i<=6; i++)
    for (j=i; j<=6; j++) {
      massDsq(i, j) = mdSqSckm(i - 3, j - 3);
	  if (i == j) massDsq(i, j) = massDsq(i, j) -
		    cos(2.0 * beta) * sqr(displayMzRun()) * 
		    sqr(sinthDRbarMS) * 1.0 / 3.0 + sqr(md(i - 3));
    }

  /// add loop corrections to the family-universal pieces
  DoubleMatrix temp(2, 2);

  /// p is for sbottom corrections
  double p = minimum(displayDrBarPars().md.display(1, 3), 
		     displayDrBarPars().md.display(2, 3));
  if (accuracy > 0) {
    int family; for (family=1; family<=3; family++) {
      temp(1, 1) = massDsq(family, family);
      temp(1, 2) = massDsq(family, family + 3);
      temp(2, 2) = massDsq(family + 3, family + 3);
      
      if (family < 3) addSdownCorrection(temp, family);
      else addSbotCorrection(p, temp, md(family));
      massDsq(family, family) = temp(1, 1);
      massDsq(family, family + 3) = temp(1, 2);
      massDsq(family + 3, family + 3) = temp(2, 2);
    }
  }

  massDsq.symmetrise();

  DoubleMatrix dSqMixT(6, 6);
  DoubleVector dSqMasses(6);
  if (massDsq.diagonaliseSym(dSqMixT, dSqMasses) > EPSTOL) {
    throw("WARNING: dn squark flavour diagonalisation bad accuracy\n");
  }

  fv.dSqMix = dSqMixT.transpose();
  if (fv.msD.min(j) < 0.) {
    flagTachyon(sdown); flagTachyon(sstrange); flagTachyon(sbottom);
    if (PRINTOUT > 1) cout << " down squark(" << j << ") tachyon";
  }
  fv.msD = dSqMasses.apply(ccbSqrt);

  // fill in the flavour conserving parts with the correct masses
  sPhysical s(displayPhys());
  int b1Pos = 0;
  DoubleMatrix a(fv.dSqMix);

  for(i=3; i>=1; i--) {
    a.displayCol(i).apply(fabs).max(j); 
    s.md(1, i) = fv.msD(j);
    for (int k=1; k<=6; k++) a(j, k) = 0.;
    if (i == 3) b1Pos = j;
    a.displayCol(i+3).apply(fabs).max(j); 
    for (int k=1; k<=6; k++) a(j, k) = 0.;
    s.md(2, i) = fv.msD(j);
  }

  s.thetab = asin(fv.dSqMix(b1Pos, 6));
  setPhys(s);
}

void FlavourMssmSoftsusy::doChargedSleptons(double /* mtau */, double /* pizztMS */,
					    double sinthDRbarMS, 
					    int accuracy) { 
  //  MssmSoftsusy::doDownSquarks(mb, pizztMS, sinthDRbarMS, accuracy, mt); 
  /// now need to re-do 6 by 6 scalar mass matrices
  DoubleMatrix massEsq(6, 6);

  /// Find the Yukawa rotation matrices such that the super-CKM basis may be
  /// defined
  DoubleMatrix Ue(3, 3), Ve(3, 3);
  DoubleVector ye(3);
  displayYukawaMatrix(YE).diagonalise(Ue, Ve, ye);

  /// super MNS basis mass squared matrices
  DoubleMatrix mlSqSckm(3, 3), teSckm(3, 3), meSqSckm(3, 3);
  mlSqSckm  = Ve.transpose() * displaySoftMassSquared(mLl) * Ve;
  teSckm    = Ue.transpose() * displayTrilinear(EA).transpose() * Ve;
  meSqSckm  = Ue.transpose() * displaySoftMassSquared(mEr).transpose() * 
    Ue; 

  DoubleVector me(3);

  double susyMu = displaySusyMu();
  double vev = displayHvev();
  double tanb = displayTanb(), beta = atan(tanb);
  double cosb = cos(beta);
  /// the following contributions are only calculated at tree-level:
  int i, j;   for (i=1; i<=3; i++) {
    me(i) = ye(i) * vev / root2 * cosb;

    for (j=i; j<=3; j++) {
      massEsq(i, j) = mlSqSckm(i, j);
      if (i == j) {
	massEsq(i, j) = massEsq(i, j) + sqr(me(i)) + 
	  cos(2.0 * beta) * sqr(displayMzRun()) * 
	  (-0.5 + sqr(sinthDRbarMS));
      }
    }
  
    for (j=4; j<=6; j++) {
      massEsq(i, j) = vev * cosb / root2 * teSckm(j-3, i);
      if (i == j-3) massEsq(i, j) = massEsq(i, j) - susyMu * me(i) * tanb;
    }
  }

  for (i=4; i<=6; i++)
    for (j=i; j<=6; j++) {
      massEsq(i, j) = meSqSckm(i - 3, j - 3);
      if (i == j) massEsq(i, j) = massEsq(i, j) -
		    cos(2.0 * beta) * sqr(displayMzRun()) * 
		    sqr(sinthDRbarMS) + sqr(me(i - 3));
    }

  /// add loop corrections to the family-universal pieces
  if (accuracy > 0) {
  DoubleMatrix temp(2, 2);
  double p = minimum(displayDrBarPars().me(1, 3), displayDrBarPars().me(2, 3));
  int family; for (family=1; family<=3; family++) {
    temp(1, 1) = massEsq(family, family);
    temp(1, 2) = massEsq(family, family + 3);
    temp(2, 2) = massEsq(family + 3, family + 3);

    if (family < 3) addSlepCorrection(temp, family);
    else addStauCorrection(p, temp, me(family));
    massEsq(family, family) = temp(1, 1);
    massEsq(family, family + 3) = temp(1, 2);
    massEsq(family + 3, family + 3) = temp(2, 2);
    }
  }
  
  massEsq.symmetrise();

  DoubleMatrix eSqMixT(6, 6);
  DoubleVector eSqMasses(6);
  if (massEsq.diagonaliseSym(eSqMixT, eSqMasses) > EPSTOL) {
    throw("WARNING: charged slepton flavour diagonalisation bad accuracy\n");
  }

  fv.eSqMix = eSqMixT.transpose();
  if (fv.msE.min(j) < 0.) {
    flagTachyon(selectron); flagTachyon(smuon); flagTachyon(stau);
    if (PRINTOUT > 1) cout << " down slepton(" << j << ") tachyon";
  }
  fv.msE = eSqMasses.apply(ccbSqrt);

  // fill in the flavour conserving parts with the correct masses
  sPhysical s(displayPhys());
  int tau1Pos = 0;
  DoubleMatrix a(fv.eSqMix);

  for(i=3; i>=1; i--) {
    a.displayCol(i).apply(fabs).max(j); 
    s.me(1, i) = fv.msE(j);
    for (int k=1; k<=6; k++) a(j, k) = 0.;
    if (i == 3) tau1Pos = j;
    a.displayCol(i+3).apply(fabs).max(j); 
    for (int k=1; k<=6; k++) a(j, k) = 0.;
    s.me(2, i) = fv.msE(j);
  }

  s.thetatau = asin(fv.eSqMix(tau1Pos, 6));
  setPhys(s);
}

void FlavourMssmSoftsusy::doSnu(double /* pizztMS */, int accuracy) {
  
  /// Find the Yukawa rotation matrices such that the super-CKM basis may be
  /// defined
  DoubleMatrix Ue(3, 3), Ve(3, 3);
  DoubleVector ye(3);
  displayYukawaMatrix(YE).diagonalise(Ve, Ue, ye);

  DoubleMatrix uMns(displayMns());

  DoubleMatrix mNuSq(uMns.transpose() * Ve.transpose() * 
		     displaySoftMassSquared(mLl) * Ve * uMns);
  double beta = atan(displayTanb());
  int i, j; for (i=1; i<=3; i++)
    for (j=1; j<=3; j++) {
      if (i == j) {
	mNuSq(i, j) = mNuSq(i, j) + 0.5 * sqr(displayMzRun()) * 
	  cos(2.0 * beta);
	if (accuracy > 0) addSnuCorrection(mNuSq(i, j), i);	
      }
    }
  
  DoubleMatrix nuSqMix(3, 3);
  DoubleVector nuSqMasses(3);
  if (mNuSq.diagonaliseSym(nuSqMix, nuSqMasses) > EPSTOL) {
    std::ostringstream ii; 
    ii << "WARNING:  sneutrino flavour diagonalisation bad accuracy\n";
    ii.setf(std::ios::scientific, std::ios::floatfield);
    ii.precision(6);
    ii << "mNuSq=" << mNuSq << " on ";
    ii << *this;
    throw ii.str();
  }
  
  fv.nuSqMix = nuSqMix.transpose();
  if (fv.msNu.min(j) < 0.) {
    flagTachyon(snue); flagTachyon(snumu); flagTachyon(snutau);
    if (PRINTOUT > 1) cout << " sneutrino(" << j << ") tachyon";
  }
  fv.msNu = nuSqMasses.apply(ccbSqrt);
  
  sPhysical s(displayPhys());
  /*  for (i=1; i<=3; i++) {
    fv.nuSqMix.displayRow(i).max(j);
    // s.msnu(i) = fv.msNu(j);
    s.msnu(i) = fv.msNu(i);
    }*/
  s.msnu = fv.msNu;
  setPhys(s);  
}

const flavourPhysical & flavourPhysical::operator=(const flavourPhysical &s) {
  if (this == &s) return *this;
  dSqMix = s.dSqMix;
  uSqMix = s.uSqMix;
  eSqMix = s.eSqMix;
  nuSqMix = s.nuSqMix;
  msD = s.msD;
  msU = s.msU;
  msE = s.msE;
  msNu = s.msNu;
  return *this;
}

void FlavourMssmSoftsusy::doQuarkMixing(DoubleMatrix & mDon, 
					DoubleMatrix & mUpq) {
  if (MIXING < 0 || MIXING > 4) {
    std::ostringstream ii;
    ii << "In FlavourMssmSoftsusy::sparticleThresholdCorrections(double tb) ";
    ii << "\n MIXING=" << MIXING << " is out of range (0 -> 4)\n";
    throw ii.str();
  }

  DoubleMatrix vCkm(displayCkm());

  switch(MIXING) {
  case 1:       
    mUpq = vCkm.transpose() * mUpq * vCkm;
    break;
  case 2: 
    mDon = vCkm * mDon * vCkm.transpose();
    break;
  case 3:
    mUpq = vCkm.transpose() * mUpq;
    break;
  case 4:
    mDon = vCkm * mDon;
    break;
  default:
    return;
  }    
}

std::ostream & operator <<(std::ostream & left, const flavourPhysical &s) {
  left << "mu: " << s.msU << "uMix: " << s.uSqMix;
  left << "md: " << s.msD << "dMix: " << s.dSqMix;
  left << "me: " << s.msE << "eMix: " << s.eSqMix;
  left << "mnu: " << s.msNu << "nuMix: " << s.nuSqMix;
  return left;
}

void extractFlavour(
   int gen, DoubleMatrix& Z, const DoubleVector& m,
   double& m1, double& m2, double& theta)
{
   int j1 = gen, j2 = gen;

   Z.displayRow(gen).absmax(j1);
   m1 = m(j1);

   if (gen == 3) theta = asin(Z(gen+3, j1));

   // set j1'th column to zero in order to not pick j1 again
   for (int i = 1; i <= Z.displayRows(); i++)
      Z(i, j1) = 0;

   Z.displayRow(gen+3).absmax(j2);
   m2 = m(j2);

   // set j2'th column to zero in order to not pick j2 again
   for (int i = 1; i <= Z.displayRows(); i++)
      Z(i, j2) = 0;
}

// calculates masses all at tree-level in the DRbar scheme, useful for
// radiative corrections. 
void FlavourMssmSoftsusy::calcDrBarPars() {
  // First, must define mstop,sbot,stau and mixing angles in DRbar scheme
  double beta = atan(displayTanb());
  double sinthDRbar = calcSinthdrbar();
  double mz = displayMzRun(), mz2 = sqr(mz);

  sw2 = sqr(sinthDRbar);
  guL = 0.5 - 2.0 * sw2 / 3.0;
  gdL = -0.5 + sw2 / 3.0;
  geL = -0.5 + sw2;
  guR = 2.0 * sw2 / 3.0;
  gdR = -sw2 / 3.0;
  geR = -sw2;

  double vev = displayHvev();

  DoubleMatrix massUsq(6, 6), massDsq(6, 6), massEsq(6, 6);
  /// Find the Yukawa rotation matrices such that the super-CKM basis may be
  /// defined
  DoubleMatrix Ud(3, 3), Vd(3, 3), Uu(3, 3), Vu(3, 3), Ue(3, 3), Ve(3, 3);
  DoubleVector ye(3), yu(3), yd(3);
  displayYukawaMatrix(YU).diagonalise(Vu, Uu, yu);
  displayYukawaMatrix(YD).diagonalise(Vd, Ud, yd);
  displayYukawaMatrix(YE).diagonalise(Ue, Ve, ye);
    ckmNormalise(Vu, Vd, Uu, Ud); 

  DoubleMatrix mqSqSckm(3, 3), mqSqUckm(3, 3), tdSckm(3, 3), mdSqSckm(3, 3), 
    tuSckm(3, 3), muSqSckm(3, 3), mlSqSckm(3, 3), teSckm(3, 3), meSqSckm(3, 3);
  mqSqSckm  = Vd.transpose() * displaySoftMassSquared(mQl) * Vd;
  mqSqUckm  = Vu.transpose() * displaySoftMassSquared(mQl) * Vu;
  tdSckm    = Ud.transpose() * displayTrilinear(DA).transpose() * Vd;
  mdSqSckm  = Ud.transpose() * displaySoftMassSquared(mDr).transpose() * 
    Ud; 
  tuSckm    = Uu.transpose() * displayTrilinear(UA).transpose() * Vu;
  muSqSckm  = Uu.transpose() * displaySoftMassSquared(mUr).transpose() * 
    Uu; 
  mlSqSckm  = Ve.transpose() * displaySoftMassSquared(mLl) * Ve;
  teSckm    = Ue.transpose() * displayTrilinear(EA).transpose() * Ve;
  meSqSckm  = Ue.transpose() * displaySoftMassSquared(mEr).transpose() * Ue; 

  DoubleVector mu(3), md(3), me(3);

  double susyMu = displaySusyMu();
  double tanb = displayTanb();
  double sinb = sin(beta), cosb = cos(beta);
  /// the following contributions are only calculated at tree-level:
  int i, j;   for (i=1; i<=3; i++) {
    me(i) = ye(i) * vev / root2 * cosb;
    mu(i) = yu(i) * vev / root2 * sinb;
    md(i) = yd(i) * vev / root2 * cosb;

    for (j=i; j<=3; j++) {
      massUsq(i, j) = mqSqUckm(i, j);
      massDsq(i, j) = mqSqSckm(i, j);
      massEsq(i, j) = mlSqSckm(i, j);
      if (i == j) {
	massUsq(i, j) = massUsq(i, j) + sqr(mu(i)) + 
	  cos(2.0 * beta) * sqr(displayMzRun()) * 
	  (0.5 - 2. / 3. * sqr(sinthDRbar));
	massDsq(i, j) = massDsq(i, j) + sqr(md(i)) + 
	  cos(2.0 * beta) * sqr(displayMzRun()) * 
	  (-0.5 + 1. / 3. * sqr(sinthDRbar));
	massEsq(i, j) = massEsq(i, j) + sqr(me(i)) + 
	  cos(2.0 * beta) * sqr(displayMzRun()) * 
	  (-0.5 + sqr(sinthDRbar));
	
      }
    }
  
    for (j=4; j<=6; j++) {
      massUsq(i, j) = vev * sinb / root2 * tuSckm(j-3, i);
      massDsq(i, j) = vev * cosb / root2 * tdSckm(j-3, i);
      massEsq(i, j) = vev * cosb / root2 * teSckm(j-3, i);
      if (i == j-3) {
	massUsq(i, j) = massUsq(i, j) - susyMu * mu(i) / tanb;
	massDsq(i, j) = massDsq(i, j) - susyMu * md(i) * tanb;
	massEsq(i, j) = massEsq(i, j) - susyMu * me(i) * tanb;
      }
    }  
  }
    
  for (i=4; i<=6; i++)
    for (j=i; j<=6; j++) {
      massUsq(i, j) = muSqSckm(i - 3, j - 3);
      massDsq(i, j) = mdSqSckm(i - 3, j - 3);
      massEsq(i, j) = meSqSckm(i - 3, j - 3);
      if (i == j) {
	massUsq(i, j) = massUsq(i, j) + 
	  cos(2.0 * beta) * sqr(displayMzRun()) * 
	  sqr(sinthDRbar) * 2.0 / 3.0 + sqr(mu(i - 3)); 
	massDsq(i, j) = massDsq(i, j) -
	  cos(2.0 * beta) * sqr(displayMzRun()) * 
	  sqr(sinthDRbar) * 1.0 / 3.0 + sqr(md(i - 3));
	massEsq(i, j) = massEsq(i, j) -
	  cos(2.0 * beta) * sqr(displayMzRun()) * 
	  sqr(sinthDRbar) + sqr(me(i - 3));
      }
    }

  massUsq.symmetrise();
  massDsq.symmetrise();
  massEsq.symmetrise();

  /// These are the TRANSPOSES of the matrices defined in SLHA2. 
  DoubleMatrix uSqMixT(6, 6), dSqMixT(6, 6), eSqMixT(6, 6);
  DoubleVector uSqMasses(6), dSqMasses(6), eSqMasses(6);
  massUsq.diagonaliseSym(uSqMixT, uSqMasses);
  massDsq.diagonaliseSym(dSqMixT, dSqMasses);
  massEsq.diagonaliseSym(eSqMixT, eSqMasses);
  for (i=1; i<=6; i++) 
    if (uSqMasses(i) < 0. || eSqMasses(i) < 0. || dSqMasses(i) < 0.) {
      if (uSqMasses(i) < 0.) { 
	flagTachyon(sup); flagTachyon(scharm); flagTachyon(stop); 
      if (PRINTOUT > 1) cout << " tree sup-type tachyon ";
      } else if (eSqMasses(i) < 0.) {
	flagTachyon(selectron); flagTachyon(smuon); flagTachyon(stau); 
      if (PRINTOUT > 1) cout << " tree slepton tachyon ";
      } else if (dSqMasses(i) < 0.) {
	flagTachyon(sdown); flagTachyon(sstrange); flagTachyon(sbottom); 
      if (PRINTOUT > 1) cout << " tree sdown-type tachyon ";
      }
      i=7; ///< bail if you find a tachyon
    }

  if (PRINTOUT > 2) cout << uSqMasses << dSqMasses << eSqMasses;
  
  uSqMasses = uSqMasses.apply(zeroSqrt);
  dSqMasses = dSqMasses.apply(zeroSqrt);
  eSqMasses = eSqMasses.apply(zeroSqrt);

  drBarPars s(displayDrBarPars());
  for(i=1; i<=3; i++) {
     extractFlavour(i, uSqMixT, uSqMasses, s.mu(1,i), s.mu(2,i), s.thetat);
     extractFlavour(i, dSqMixT, dSqMasses, s.md(1,i), s.md(2,i), s.thetab);
     extractFlavour(i, eSqMixT, eSqMasses, s.me(1,i), s.me(2,i), s.thetatau);
  }

  DoubleMatrix uMns(displayMns());
  DoubleMatrix mNuSq(uMns.transpose() * Ve.transpose() * 
		     displaySoftMassSquared(mLl) * Ve * uMns);
  for (i=1; i<=3; i++)
    for (j=1; j<=3; j++) {
      if (i == j) {
	mNuSq(i, j) = mNuSq(i, j) + 0.5 * sqr(displayMzRun()) * 
	  cos(2.0 * beta);
      }
    }
  DoubleMatrix nuSqMix(3, 3);
  DoubleVector nuSqMasses(3);
  mNuSq.diagonaliseSym(nuSqMix, nuSqMasses);
  nuSqMasses = nuSqMasses.apply(zeroSqrt);
  for (i=1; i<=3; i++) {
    if (nuSqMasses(i) < EPSTOL) {
      flagTachyon(snue); flagTachyon(snumu); flagTachyon(snutau);
      if (PRINTOUT > 0) cout << " sneutrino type tachyon ";
      i=4; ///< bail out
    }
    nuSqMix.displayRow(i).max(j);
    s.msnu(i) = nuSqMasses(j);
  }

  s.ht = yu(3);        s.hb = yd(3);        s.htau = ye(3);
  s.ut = tuSckm(3, 3); s.ub = tdSckm(3, 3); s.utau = teSckm(3, 3);

  s.mt   = s.ht   * vev * sin(beta) / root2;
  s.mb   = s.hb   * vev * cos(beta) / root2;
  s.mtau = s.htau * vev * cos(beta) / root2; 

  s.mGluino = displayGaugino(3);

  double mw = displayMwRun();
  double mw2 = sqr(mw);
  calcDrBarGauginos(beta, mw, mz, sinthDRbar, s);
  s.mw = mw;
  s.mz = mz;

  //  setDrBarPars(s); 
  calcDrBarHiggs(beta, mz2, mw2, sinthDRbar, s);  

  setDrBarPars(s);

  return;
}

MssmSusy FlavourMssmSoftsusy::guessAtSusyMt(double tanb, 
					    const QedQcd_legacy & oneset) { 

  // This bit gives a guess at a SUSY object
  QedQcd_legacy leAtMt(oneset);

  DoubleVector a(3), g(3);
  double sinth2 = 1.0 - sqr(MW / MZ);
  // Gauge couplings at mt
  a = leAtMt.getGaugeMu(oneset.displayPoleMt(), sinth2);
  
  MssmSusy t(MssmSoftsusy::guessAtSusyMt(tanb, oneset));
  DoubleMatrix yu(t.displayYukawaMatrix(YU)), yd(t.displayYukawaMatrix(YD)); 
  
  doQuarkMixing(yd, yu);
  t.setYukawaMatrix(YU, yu);
  t.setYukawaMatrix(YD, yd);

  return t;
}

void FlavourMssmSoftsusy::setAngles(double lambdaW, double aCkm, double rhobar, 
				    double etabar) {
  theta12 = asin(lambdaW);
  theta23 = asin(aCkm * sqr(lambdaW));

  double lambdaW3 = lambdaW * sqr(lambdaW);
  double lambdaW4 = lambdaW * lambdaW3;

  Complex rpe(rhobar, etabar);
  Complex lhs = aCkm * lambdaW3 * rpe * sqrt(1.0 - sqr(aCkm) * lambdaW4) / 
    sqrt(1.0 - sqr(lambdaW)) / (1.0 - sqr(aCkm) * lambdaW4 * rpe);

  theta13 = asin(lhs.mod());
  deltaCkm = lhs.arg();
}

void FlavourMssmSoftsusy::sCkm 
(DoubleMatrix & deltaULL, DoubleMatrix & deltaURR, DoubleMatrix & deltaULR, 
 DoubleMatrix & deltaDLL, DoubleMatrix & deltaDRR, DoubleMatrix & deltaDLR) 
  const {

  flavourPhysical fv(displayFlavourPhysical());

  DoubleMatrix mupsq(6, 6), msdsq(6, 6), mudiagsq(6, 6), mddiagsq(6, 6);
  int i; for (i=1; i<=6; i++) {
    mudiagsq(i, i) = sqr(fv.msU(i)); mddiagsq(i, i) = sqr(fv.msD(i));
  }

  mupsq = fv.uSqMix.transpose() * mudiagsq * fv.uSqMix;
  msdsq = fv.dSqMix.transpose() * mddiagsq * fv.dSqMix;

  int j;
  for (i=1; i<=3; i++)
    for (j=1; j<=3; j++) {
      deltaULL(i, j) = mupsq(i, j) / 
	(sqrt(fabs(mupsq(i, i) * mupsq(j, j))));
      deltaURR(i, j) = mupsq(i+3, j+3) / 
	(sqrt(fabs(mupsq(i+3, i+3) * mupsq(j+3, j+3))));
      deltaULR(i, j) = mupsq(i, j+3) / (sqrt(fabs(mupsq(i, i) * 
						 mupsq(j+3, j+3))));
      deltaDLL(i, j) = msdsq(i, j) / 
	  (sqrt(fabs(msdsq(i, i) * msdsq(j, j))));
      deltaDRR(i, j) = msdsq(i+3, j+3) / 
	(sqrt(fabs(msdsq(i+3, i+3) * msdsq(j+3, j+3))));
      deltaDLR(i, j) = msdsq(i, j+3) / 
	(sqrt(fabs(msdsq(i, i) * msdsq(j+3, j+3))));
    }
}


void FlavourMssmSoftsusy::sCkmAnalytic
(double m32, double & deltaULL23, double & deltaULR23,
 double & deltaULL12, double & deltaDLL23, double & deltaDLR23, 
 double & deltaDLL12)  const { 

  drBarPars d(displayDrBarPars());

  double mfac = m32 / (16 * sqr(PI));
  
  double beta = atan(displayTanb());
  double c2b = cos(2 * beta);
  double mt = d.mt;
  double sinth = calcSinthdrbar();
  double sinthw2 = sqr(sinth);
  double lt = d.ht, lb = d.hb, ltau = d.htau;
  //  double ls = displayYukawaElement(YD, 2, 2);
  double gs = displayGaugeCoupling(3);
  double gw = displayGaugeCoupling(2);
  double g1 = displayGaugeCoupling(1);
  double mb = d.mb;


  double vd = displayHvev() * cos(beta);
  double vu = displayHvev() * sin(beta);
  
  DoubleMatrix Vu(3, 3), Uu(3, 3); DoubleVector yu(3);
  displayYukawaMatrix(YU).diagonalise(Vu, Uu, yu);
  DoubleMatrix Vd(3, 3), Ud(3, 3); DoubleVector yd(3);
  displayYukawaMatrix(YD).diagonalise(Vd, Ud, yd);
  ckmNormalise(Vu, Vd, Uu, Ud); 
  DoubleMatrix V(displayCkm());
  // DoubleMatrix V(Vu.transpose() * Vd);

  double ddownL = sqr(MZ) * (-1./2. + 1./3. * sinthw2) * c2b;
  double ddownR = sqr(MZ) * (-1./3.) * sinthw2 * c2b;
  double dupL = sqr(MZ) * (0.5 - 2./3. * sinthw2) * c2b;
  double dupR = sqr(MZ) * 2/3 * sinthw2 * c2b;
  double bett = -13./15. * sqr(g1) - 3. * sqr(gw) - 16/3. * sqr(gs) + 
    6. * sqr(lt) + sqr(lb);
  double betb = -7./15. * sqr(g1) - 3. * sqr(gw) - 16/3. * sqr(gs) + 
    6. * sqr(lb) + sqr(ltau) + sqr(lt);

  double mdLR23 = -vd / sqrt(2.0) * mfac * lb * sqr(lt) * V(3, 2) * V(3, 3);
  double mdL23 = sqr(mfac) * V(3, 2) * V(3, 3) * sqr(lt) * (bett);
  double mdLii = ddownL + sqr(mfac) * 
      (-11./50. * sqr(sqr(g1)) - 1.5 * sqr(sqr(gw)) + 8 * sqr(sqr(gs)));
  double mdLij = sqr(mfac) * V(3, 1) * V(3, 2) * sqr(lt) * (bett - sqr(lb));

  /// Ben's correction factor
  /*  mdLij += sqr(mfac) * V(2, 1) * V(2, 2) * sqr(lc) * 
      (bett - sqr(lb) - 3.0 * sqr(lt));*/

  double mdL33 = sqr(mb) + ddownL + 
    sqr(mfac) * (-11./50. * sqr(sqr(g1)) - 1.5 * sqr(sqr(gw)) + 
		 8 * sqr(sqr(gs)) + sqr(lt) * bett + sqr(lb) * betb);
  double mdR33 = sqr(mb) + ddownR + sqr(mfac) * 
    (-22./25. * sqr(sqr(g1)) + 8 * sqr(sqr(gs)) + 2 * sqr(lb) * betb);

  deltaDLL23 = mdL23 / sqrt(mdL33 * mdLii);
  deltaDLR23 = mdLR23 / sqrt(mdLii * mdR33);
  deltaDLL12 = mdLij / mdLii;

  double muL33 = sqr(mt) + dupL + 
    sqr(mfac) * (-11./50. * 
		 sqr(sqr(g1)) - 1.5 * sqr(sqr(gw)) + 8 * sqr(sqr(gs)) 
		 + sqr(lt) * bett + sqr(lb) * betb);
  double muLi3 = sqr(mfac) * V(2, 3) * V(3, 3) * sqr(lb) * betb;

  double muR33 = sqr(mt) + dupR + sqr(mfac) * (-88./25. * sqr(sqr(g1)) + 
					       8 * sqr(sqr(gs)) + 
					       2 * sqr(lt) * bett);
  
  double muLRi3 = -vu * mfac * lt * sqr(lb) * V(2, 3) * V(3, 3) / sqrt(2.0);

  double muLij = sqr(mfac) * V(1, 3) * V(2, 3) * sqr(lb) * (betb - sqr(lt));

  // NB this is Ben's correction factor
  /*  muLij += sqr(mfac) * V(1, 2) * V(2, 2) * sqr(ls) * 
      (betb - sqr(lt) - 3.0 * sqr(lb));*/

  double muLii = dupL + sqr(mfac) * (-11./50. * sqr(sqr(g1)) - 
				     1.5 * sqr(sqr(gw)) + 8 * sqr(sqr(gs)));
  deltaULL12 = muLij / muLii;
  deltaULL23 = muLi3 / sqrt(muL33 * muLii);
  deltaULR23 = muLRi3 / sqrt(muLii * muR33);
}

void ckmNormalise(DoubleMatrix & Vu, DoubleMatrix & Vd, DoubleMatrix & Uu,
		  DoubleMatrix & Ud) {
  if (Vu.displayRows() != 3 || Vu.displayCols() != 3 || 
      Vd.displayRows() != 3 || Vd.displayCols() != 3) {
    throw("Wrong size of matrices in flavoursoft.cpp::ckmNormalise\n");
  }

  DoubleMatrix vCkm(Vu.transpose() * Vd);
  
  DoubleMatrix signsU(3, 3), signsD(3, 3);

  int i; for (i=1; i<=3; i++) {
    if (vCkm(i, i) < 0.) {
      signsU(i, i) = -1.; 
      int j; for (j=1; j<=3; j++) vCkm(i, j) = -vCkm(i, j); 
    }
    else signsU(i, i) = 1.;
    signsD(i, i) = 1.;
  }

  /// Now, we make 12 and 23 elements positive while keeping diagonal
  /// element positive
  if (vCkm(1, 2) < 0.) {
    signsD(2, 2) = -1.; signsU(2, 2) = -signsU(2, 2);
    int j; for (j=1; j<=3; j++) {
      vCkm(2, j) = -vCkm(2, j);  
      vCkm(j, 2) = -vCkm(j, 2);
    }
  }

  if (vCkm(2, 3) < 0.) {
    signsD(3, 3) = -1.; signsU(3, 3) = -signsU(3, 3);
  }
  /// We've now exhausted the 5 re-signing degrees of freedom that we have to
  /// play with. 

  /// On the other hand, this fixes the phases definitively to avoid any
  /// problems 
  /*  if (signsU(1, 1) * Vu(1, 1) < 0.) {
    for (i = 1; i<=3; i++) {
      signsU(i, i) = -signsU(i, i);
      signsD(i, i) = -signsD(i, i);
    }
    }*/

  Vu = Vu * signsU; Uu = Uu * signsU;
  Vd = Vd * signsD; Ud = Ud * signsD;
  /// latter relations in order to obtain same sign for the fermion mass
}


void FlavourMssmSoftsusy::sCkmRunning
(DoubleMatrix & deltaULL, DoubleMatrix & deltaURR, DoubleMatrix & deltaULR, 
 DoubleMatrix & deltaDLL, DoubleMatrix & deltaDRR, DoubleMatrix & deltaDLR) 
  {

  calcDrBarPars();

  // First, must define mstop,sbot,stau and mixing angles in DRbar scheme
  double beta = atan(displayTanb());
  double sinthDRbar = calcSinthdrbar();

  sw2 = sqr(sinthDRbar);
  guL = 0.5 - 2.0 * sw2 / 3.0;
  gdL = -0.5 + sw2 / 3.0;
  geL = -0.5 + sw2;
  guR = 2.0 * sw2 / 3.0;
  gdR = -sw2 / 3.0;
  geR = -sw2;

  double vev = displayHvev();

  DoubleMatrix massUsq(6, 6), massDsq(6, 6), massEsq(6, 6);
  /// Find the Yukawa rotation matrices such that the super-CKM basis may be
  /// defined
  DoubleMatrix Ud(3, 3), Vd(3, 3), Uu(3, 3), Vu(3, 3), Ue(3, 3), Ve(3, 3);
  DoubleVector ye(3), yu(3), yd(3);
  displayYukawaMatrix(YU).diagonalise(Vu, Uu, yu);
  displayYukawaMatrix(YD).diagonalise(Vd, Ud, yd);
  displayYukawaMatrix(YE).diagonalise(Ue, Ve, ye);
  ckmNormalise(Vu, Vd, Uu, Ud); 

  DoubleMatrix mqSqSckm(3, 3), mqSqUckm(3, 3), tdSckm(3, 3), mdSqSckm(3, 3), 
    tuSckm(3, 3), muSqSckm(3, 3), mlSqSckm(3, 3), teSckm(3, 3), meSqSckm(3, 3);
  mqSqSckm  = Vd.transpose() * displaySoftMassSquared(mQl) * Vd;
  mqSqUckm  = Vu.transpose() * displaySoftMassSquared(mQl) * Vu;
  tdSckm    = Ud.transpose() * displayTrilinear(DA).transpose() * Vd;
  mdSqSckm  = Ud.transpose() * displaySoftMassSquared(mDr).transpose() * 
    Ud; 
  tuSckm    = Uu.transpose() * displayTrilinear(UA).transpose() * Vu;
  muSqSckm  = Uu.transpose() * displaySoftMassSquared(mUr).transpose() * 
    Uu; 
  mlSqSckm  = Ve.transpose() * displaySoftMassSquared(mLl) * Ve;
  teSckm    = Ue.transpose() * displayTrilinear(EA).transpose() * Ve;
  meSqSckm  = Ue.transpose() * displaySoftMassSquared(mEr).transpose() * Ue; 

  DoubleVector mu(3), md(3), me(3);

  double susyMu = displaySusyMu();
  double tanb = displayTanb();
  double sinb = sin(beta), cosb = cos(beta);
  /// the following contributions are only calculated at tree-level:
  int i, j;   for (i=1; i<=3; i++) {
    me(i) = ye(i) * vev / root2 * cosb;
    mu(i) = yu(i) * vev / root2 * sinb;
    md(i) = yd(i) * vev / root2 * cosb;

    for (j=i; j<=3; j++) {
      massUsq(i, j) = mqSqUckm(i, j);
      massDsq(i, j) = mqSqSckm(i, j);
      massEsq(i, j) = mlSqSckm(i, j);
      if (i == j) {
	massUsq(i, j) = massUsq(i, j) + sqr(mu(i)) + 
	  cos(2.0 * beta) * sqr(displayMzRun()) * 
	  (0.5 - 2. / 3. * sqr(sinthDRbar));
	massDsq(i, j) = massDsq(i, j) + sqr(md(i)) + 
	  cos(2.0 * beta) * sqr(displayMzRun()) * 
	  (-0.5 + 1. / 3. * sqr(sinthDRbar));
	massEsq(i, j) = massEsq(i, j) + sqr(me(i)) + 
	  cos(2.0 * beta) * sqr(displayMzRun()) * 
	  (-0.5 + sqr(sinthDRbar));
	
      }
    }
  
    for (j=4; j<=6; j++) {
      massUsq(i, j) = vev * sinb / root2 * tuSckm(j-3, i);
      massDsq(i, j) = vev * cosb / root2 * tdSckm(j-3, i);
      massEsq(i, j) = vev * cosb / root2 * teSckm(j-3, i);
      if (i == j-3) {
	massUsq(i, j) = massUsq(i, j) - susyMu * mu(i) / tanb;
	massDsq(i, j) = massDsq(i, j) - susyMu * md(i) * tanb;
	massEsq(i, j) = massEsq(i, j) - susyMu * me(i) * tanb;
      }
    }  
  }
    
  for (i=4; i<=6; i++)
    for (j=i; j<=6; j++) {
      massUsq(i, j) = muSqSckm(i - 3, j - 3);
      massDsq(i, j) = mdSqSckm(i - 3, j - 3);
      massEsq(i, j) = meSqSckm(i - 3, j - 3);
      if (i == j) {
	massUsq(i, j) = massUsq(i, j) + 
	  cos(2.0 * beta) * sqr(displayMzRun()) * 
	  sqr(sinthDRbar) * 2.0 / 3.0 + sqr(mu(i - 3)); 
	massDsq(i, j) = massDsq(i, j) -
	  cos(2.0 * beta) * sqr(displayMzRun()) * 
	  sqr(sinthDRbar) * 1.0 / 3.0 + sqr(md(i - 3));
	massEsq(i, j) = massEsq(i, j) -
	  cos(2.0 * beta) * sqr(displayMzRun()) * 
	  sqr(sinthDRbar) + sqr(me(i - 3));
      }
    }

  massUsq.symmetrise();
  massDsq.symmetrise();
  massEsq.symmetrise();

  for (i=1; i<=3; i++)
    for (j=1; j<=3; j++) {
      deltaULL(i, j) = massUsq(i, j) / 
	(sqrt(fabs(massUsq(i, i) * massUsq(j, j))));
      deltaURR(i, j) = massUsq(i+3, j+3) / 
	(sqrt(fabs(massUsq(i+3, i+3) * massUsq(j+3, j+3))));
      deltaULR(i, j) = massUsq(i, j+3) / (sqrt(fabs(massUsq(i, i) * 
						 massUsq(j+3, j+3))));
      deltaDLL(i, j) = massDsq(i, j) / 
	  (sqrt(fabs(massDsq(i, i) * massDsq(j, j))));
      deltaDRR(i, j) = massDsq(i+3, j+3) / 
	(sqrt(fabs(massDsq(i+3, i+3) * massDsq(j+3, j+3))));
      deltaDLR(i, j) = massDsq(i, j+3) / 
	(sqrt(fabs(massDsq(i, i) * massDsq(j+3, j+3))));
    }
}

} // namespace softsusy
