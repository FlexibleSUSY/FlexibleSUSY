
/** \file softsusy.h
   - Project:     SOFTSUSY 
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/

   \brief Header file for RP conserving MSSM object including all 
              (real) soft SUSY breaking parameters and (real) SUSY couplings.
	      Also contains MssmSoftsusy2 object, which has alternative EWSB
	      boundary conditions at MSUSY.
*/


#ifndef SOFTSUSY_H
#define SOFTSUSY_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include "def.h"
#include "utils.h"
#include "numerics.h"
#include "physpars.h"
#include "lowe_legacy.h"
#include "softpars.h"

#include "mssm_twoloophiggs.h"
#include "mssmUtils.h"

namespace softsusy {

using std::fabs;
using std::sqrt;

///< default SUSY breaking boundary condition scale
const double mxDefault = 1.9e16; 

/// A different REWSB condition - given by mu and MA instead of mh1,2
class AltEwsbMssm {
private:
  /// Sets which alternative condition to use: 1 is mu(MSUSY) and MApole
  /// 2 is mu(mInput) and tan beta
  double mAcond, muCond; ///< user set conditions on mA and mu at M_SUSY
public:
  inline void setAltEwsb(double ma, double mu) { mAcond = ma; muCond = mu; };
  inline void setAltEwsbMssm(const AltEwsbMssm & s) { *this = s; };

  //  inline int displayConditionStyle() { return conditionStyle; };
  inline double displayMaCond() const { return mAcond; };
  inline double displayMuCond() const { return muCond; };
  inline const AltEwsbMssm & displayAltEwsbMssm() const { return *this; };

  inline void setMaCond(double maInput) { mAcond = maInput;  };
  inline void setMuCond(double muInput) { muCond = muInput;  };
};

/// Contains all supersymmetric MSSM parameters, incorporating R_p MSSM
template<class SoftPars>
class Softsusy: public SoftPars, public AltEwsbMssm {
  /// Includes 
  /// - Soft terms
  /// - DRbar masses
  /// - RGEs
  /// - Higgs VEVs, tadpoles
  /// - Info on any problems with the parameter space point
  /// - main routines for driving the calculation
  /// - loop corrections for sparticle masses and Standard Model couplings
  /// - fine-tuning and flavour changing neutral currents calculation
private:
  sPhysical physpars; ///< Contains pole masses and mixings of sparticles
  drBarPars forLoops; ///< Contains DRbar tree-level masses
  sProblem problem;   ///< Contains problem flags 
  double msusy;       ///< Scale at which Higgs potential is minimised
  double minV;        ///< Value of Higgs potential at minimum
  double mw;          ///< Pole W mass prediction
  QedQcd_legacy dataSet;     ///< contains low energy data on quark masses etc
  double fracDiff;    ///< fractional difference to last iteration
  bool setTbAtMX;     ///< flag: do we set tan beta at the SUSY breaking scale?
  bool altEwsb;       ///< flag: do we set mu, mA at the SUSY breaking scale?
  double predMzSq;    ///< predicted Z mass squared after iteration
  double t1OV1Ms, t2OV2Ms;  ///< DRbar tadpoles(MSusy): incl 2 loops
  double t1OV1Ms1loop, t2OV2Ms1loop; ///< DRbar tadpoles(MSusy): excl 2 loops
  bool alternativeMs;

protected:
  void setT1OV1Ms(double t1) { t1OV1Ms = t1; } 
  void setT2OV2Ms(double t2) { t2OV2Ms = t2; } 
  void setT1OV1Ms1loop(double t1) { t1OV1Ms1loop = t1; }
  void setT2OV2Ms1loop(double t2) { t2OV2Ms1loop = t2; }
  double mxBC;        ///< Scale at which SUSY breaking boundary conditions set

public:
  using SoftPars::set;
  using SoftPars::setPars;
  using SoftPars::setSoftPars;
  using SoftPars::setMu;
  using SoftPars::setLoops;
  using SoftPars::setThresholds;
  using SoftPars::setSusy;
  using SoftPars::setSusyMu;
  using SoftPars::setHvev;
  using SoftPars::setMh1Squared;
  using SoftPars::setMh2Squared;
  using SoftPars::setM3Squared;
  using SoftPars::setM32;
  using SoftPars::setTanb;
  using SoftPars::setGaugeCoupling;
  using SoftPars::setYukawaMatrix;
  using SoftPars::setYukawaElement;
  using SoftPars::display;
  using SoftPars::displayGaugeCoupling;
  using SoftPars::displayYukawaElement;
  using SoftPars::displayYukawaMatrix;
  using SoftPars::displayLoops;
  using SoftPars::displaySusyMu;
  using SoftPars::displaySusy;
  using SoftPars::displayMu;
  using SoftPars::displayThresholds;
  using SoftPars::displayGauge;
  using SoftPars::displayMh1Squared;
  using SoftPars::displayMh2Squared;
  using SoftPars::displayM3Squared;
  using SoftPars::displayTanb;
  using SoftPars::displayHvev;
  using SoftPars::displayGaugino;
  using SoftPars::displaySoftMassSquared;
  using SoftPars::displaySoftA;
  using SoftPars::displayGravitino;
  using SoftPars::displayTrilinear;
  using SoftPars::run;
  using SoftPars::runto;

  typedef typename SoftPars::susy_type Susy;

  //  void (*boundaryCondition)(Softsusy &, const DoubleVector &);
  /// Default constructor fills object with zeroes
  Softsusy();
  /// Constructor sets SUSY parameters only from another object
  Softsusy(const Susy &);
  /// Constructor copies another object

  Softsusy(const Softsusy &);
  /// Sets all parameters from s, sp, mu is the mu superpotential parameter, l
  /// is the number of loops used for RG evolution, t is the thresholds
  /// accuracy parameter, mg is the gravitino mass, hv is the Higgs VEV
  /// parameter. 
  Softsusy(const SoftPars & s, const sPhysical & sp, double mu, int l, 
	       int t, double hv);

  /// Set all data in the object equal to another
  const Softsusy & operator=(const Softsusy & s);
  
  /// Displays whole object as a const
  inline const Softsusy & displayMssmSoft() const;

  /// Displays physical parameters only
  inline const sPhysical & displayPhys() const;

  /// Displays tree-level masses and mixings of sparticles and third
  /// generation fermions
  inline const drBarPars & displayDrBarPars() const;
  /// Returns any problem flags associated with the object
  const sProblem & displayProblem() const {return problem; };
  /// Gives the low energy Standard Model data set used for the object
  inline const QedQcd_legacy & displayDataSet() const;
  /// Displays iteration accuracy attained 
  inline double displayFracDiff() const { return fracDiff; }; 
  double displayMinpot() const;    ///< Returns minimum of Higgs potential
  double displayMsusy() const; ///< Returns Higgs minimisation scale
  double displayMw() const; ///< Returns predicted pole MW
  /// Returns DRbar MW, must be calculated by calcDrBarPars first
  double displayMwRun() const; 
  /// Returns DRbar MZ, must be calculated by calcDrBarPars first
  double displayMzRun() const; 
  double displayTadpole1Ms() const; ///< displays t_1/v_1 tadpole
  double displayTadpole2Ms() const; ///< displays t_2/v_2 tadpole
  double displayTadpole1Ms1loop() const; ///< displays t_1/v_1 tadpole@1 loop
  double displayTadpole2Ms1loop() const; ///< displays t_2/v_2 tadpole@1 loop
  double displayMxBC() const { return mxBC; }; ///< displays M_X value
  /// Returns object as a const
  const Softsusy & displaySoftsusy() const { return *this; }
  /// Returns value of pole MZ being used
  double displayMz() const { return displayDataSet().displayMu(); }
  /// Is tan beta set at the user defined SUSY breaking scale?
  bool displaySetTbAtMX() const { return setTbAtMX; } 
  bool displayAltEwsb() const { return altEwsb; }
  double displayPredMzSq() const { return predMzSq; }

  /// Flags weird mgut-type problems
  void flagMgutOutOfBounds(bool a) { problem.mgutOutOfBounds = a; };
  /// Flags Infra-red quasi fixed point breach problem
  void flagIrqfp(bool a) { problem.irqfp = a; };
  /// Flags non-perturbative RG evolution
  void flagNonperturbative(bool a) { problem.nonperturbative = a; };
  /// Flags a negative-mass squared scalar (really a CCB problem)
  void flagTachyon(tachyonType a) { 
    problem.tachyon = a; 
    if (PRINTOUT > 2) cout << tachyonNames[a] << " tachyon ";
  };
  /// Flags problem with Higgs potential minimum
  void flagM3sq(bool a) { problem.m3sq = a; };
  /// Flags fact that calculation hasn't acheived required accuracy
  void flagNoConvergence(bool a) { problem.noConvergence = a; };
  /// Flags fact that mu sub iteration didn't converge
  void flagNoMuConvergence(bool a) { problem.noMuConvergence = a; };
  /// Flags fact that rho parameter sub iteration didn't converge
  void flagNoRhoConvergence(bool a) { problem.noRhoConvergence = a; };
  /// Flags point inconsistent with electroweak symmetry breaking
  void flagMusqwrongsign(bool a) { problem.muSqWrongSign = a; };
  /// Flags point where higher order corrections to higgs mass are potentially
  /// very large
  void flagInaccurateHiggsMass(bool a) { problem.inaccurateHiggsMass = a; };
  /// Flags an inconsistent Higgs minimum
  void flagHiggsufb(bool a) { problem.higgsUfb = a; };
  /// Sets all problems equal to either true or false (contained in a)
  void flagAllProblems(bool a, tachyonType b) { problem.irqfp = a; 
    problem.tachyon = b; problem.m3sq = a; problem.badConvergence = a;
    problem.noConvergence = a; problem.higgsUfb = a;
    problem.nonperturbative = a; problem.noRhoConvergence = a; 
    problem.noMuConvergence = a; problem.muSqWrongSign = a; 
    problem.inaccurateHiggsMass = b; problem.mgutOutOfBounds = a; }
  /// Flags a numerical exception eg number too big/small
  void flagProblemThrown(bool a) { problem.problemThrown = a; }
  
  /// Sets whole object equal to another  
  void setSoftsusy(const Softsusy & s) { *this = s; };
  /// Sets low energy Standard Model fermion mass and gauge coupling data
  void setData(const QedQcd_legacy & r) { dataSet = r; };
  /// Sets potential value at minimum of Higgs potential
  void setMinpot(double);
  /// Sets scale of Higgs potential minimisation and sparticle mass calculation
  void setMsusy(double);
  /// sets pole MW prediction
  void setMw(double);
  /// Sets all physical parameters
  void setPhys(const sPhysical & s) { physpars = s; };
  /// Sets the scale at which high-scale boundary conditions are applied
  void setMxBC(double mx) { mxBC = mx; };
  /// Sets tree-level DRbar parameters
  void setDrBarPars(const drBarPars & s) { forLoops = s; };
  /// Sets the setTbAtMX flag
  void setSetTbAtMX(bool a) { setTbAtMX = a; }
  /// Use alernative EWSB conditions: set mu and MA(pole)
  void useAlternativeEwsb() { altEwsb = true; }
  /// Set MZ^2 predicted after iteration
  void setPredMzSq(double a) { predMzSq = a; }
  void setAlternativeMs(bool flag) { alternativeMs = flag; }
  void setTadpole1Ms(double t) { t1OV1Ms = t; }
  void setTadpole2Ms(double t) { t2OV2Ms = t; }

  //PA: sets fracDiff, needed for access by NmssmSoftsusy methods
  void setFracDiff(double fD) { fracDiff = fD; };

  //PA: fixes trilnear H1-sfermion-sfermion couplings 
  //for use in doCalcTadpole1oneLoop  
  void H1SfSfCouplings(DoubleMatrix & lTS1Lr, DoubleMatrix & lBS1Lr, DoubleMatrix & lTauS1Lr, double gmzOcthW, double mu, double cosb, double v1) const;
  //PA: fixes trilnear H2-sfermion-sfermion couplings 
  //for use in doCalcTadpole1oneLoop
  void H2SfSfCouplings(DoubleMatrix & lTS2Lr, DoubleMatrix & lBS2Lr, DoubleMatrix & lTauS2Lr, double gmzOcthW, double mu, double sinb) const;
  //PA: routine to calculate sfermiom contributions to H1 tadpole / v1
  double doCalcTad1Sfermions(DoubleMatrix lTS1Lr,  DoubleMatrix lBS1Lr,  DoubleMatrix lTauS1Lr, double costhDRbar) const;
  //PA: routine to calculate sfermiom contributions to (16 \pi^2) t1 / v1
  double doCalcTad2Sfermions(DoubleMatrix lTS2Lr, DoubleMatrix lBS2Lr, DoubleMatrix lTauS2Lr, double costhDRbar) const;
  //PA: fixes trilnear H1-fermion-fermion couplings 
  //for use in doCalcTadpole1oneLoop  
  double doCalcTad1fermions(double q, double v1) const;
  //PA: fixes trilnear H2-fermion-fermion couplings 
  //for use in doCalcTadpole1oneLoop  
  double doCalcTad2fermions(double q) const;
  //PA: one loop H1 tadpole contributions from Higgs bosons in the loops
  // Follwing BPMZ Goldstone bosons are not included in this.
  double doCalcTad1Higgs(double q, double costhDRbar2, double g, double tanb) const;
  //one loop H2 tadpole contributions from Higgs bosons in the loops
  // Follwing BPMZ Goldstone bosons are not included in this.
  double doCalcTad2Higgs(double q, double costhDRbar2, double g, double tanb) const;
  //PA: one loop H1 tadpole contributions from Neutralinos in the loops
  double doCalcTad1Neutralinos(double q, double costhDRbar, double g, 
                             double cosb) const;
 //PA: one loop H2 tadpole contributions from Neutralinos in the loops
  double doCalcTad2Neutralinos(double q, double costhDRbar, double g, double sinb) const;
   //PA: one loop H1 tadpole contributions from Charginos in the loops
  double doCalcTad1Charginos(double q, double g, double cosb) const;

  double doCalcTad2Charginos(double q, double g, double sinb) const;
  
 //PA: one loop H1 tadpole contributions from Charginos in the loops
  double doCalcTad1GaugeBosons(double q, double costhDRbar2, double g, 
                               double tanb) const;

  double doCalcTad2GaugeBosons(double q, double costhDRbar2, 
			       double g, double tanb) const;
  /// Does the full 2-loop calculation of both tadpoles and sets them
  void doTadpoles(double mt, double sinthDRbar);
  /// Does the calculation of one-loop pieces of \f$ t_1 / v_1 \f$ 
  virtual double doCalcTadpole1oneLoop(double mt, double sinthDRbar);
  /// Does the calculation of one-loop pieces of \f$ t_2 / v_2 \f$ 
  virtual double doCalcTadpole2oneLoop(double mt, double sinthDRbar);
  /// Calculates and sets the one-loop pieces of \f$ t_1 / v_1 \f$: sets both
  /// 1-loop and total pieces equal to the one-loop piece
  virtual void calcTadpole1Ms1loop(double mt, double sinthDRbar);
  /// Calculates then sets the one-loop pieces of \f$ t_2 / v_2 \f$: sets both
  /// 1-loop and total pieces equal to the one-loop piece
  virtual void calcTadpole2Ms1loop(double mt, double sinthDRbar);
  /// LCT: Adds one-loop strong corrections to stop mass matrix
  /// IO parameters: p=external momentum, mt=DR bar top mass, 
  /// outputs strong=2x2 matrix of self-energies in LR basis
  DoubleMatrix addStopQCD(double p, double mt, DoubleMatrix & strong);
  /// LCT: Adds one-loop stop corrections to stop mass matrix
  /// IO parameters: p=external momentum, mt=DR bar top mass, 
  /// outputs stop=2x2 matrix of self-energies in LR basis
  DoubleMatrix addStopStop(double p, double mt, DoubleMatrix & stop);
  /// LCT: Adds one-loop sbottom corrections to stop mass matrix
  /// IO parameters: p=external momentum, mt=DR bar top mass, 
  /// outputs sbottom=2x2 matrix of self-energies in LR basis
  DoubleMatrix addStopSbottom(double p, double mt, DoubleMatrix & sbottom);
  /// LCT: Adds one-loop Higgs corrections to stop mass matrix
  /// IO parameters: p=external momentum, mt=DR bar top mass, 
  /// outputs higgs=2x2 matrix of self-energies in LR basis
  DoubleMatrix addStopHiggs(double p, double mt, DoubleMatrix & higgs);
 /// LCT: Adds one-loop electroweak corrections to stop mass matrix
  /// IO parameters: p=external momentum, mt=DR bar top mass, 
  /// outputs electroweak=2x2 matrix of self-energies in LR basis
  DoubleMatrix addStopEweak(double p, DoubleMatrix & electroweak);
  /// LCT: Adds one-loop neutralino corrections to stop mass matrix
  /// IO parameters: p=external momentum, outputs neutralino=2x2 matrix of 
  /// self-energies in LR basis
  DoubleMatrix addStopNeutralino(double p, double mt, DoubleMatrix & neutralino);
  /// LCT: Adds one-loop chargino corrections to stop mass matrix
  /// IO parameters: p=external momentum, outputs chargino=2x2 matrix of 
  /// self-energies in LR basis
  DoubleMatrix addStopChargino(double p, DoubleMatrix & chargino);
  /// Adds one-loop corrections to stop mass matrix
  /// IO parameters: p=external momentum, mass=tree level mass matrix on
  /// input, is returned with radiative corrections added, mt=DR bar top mass
  virtual void addStopCorrection(double p, DoubleMatrix & mass, double mt);
  /// LCT: Adds one-loop strong corrections to sdown mass matrix
  /// IO parameters: {p1,p2}=external momentum, outputs strong=2x2 matrix of 
  /// self-energies in LR basis
  DoubleMatrix addSdownQCD(double p1, double p2, int family, 
  DoubleMatrix & strong); 
  /// LCT: Adds one-loop Higgs corrections to sdown mass matrix
  /// IO parameters: {p1,p2}=external momentum, outputs higgs=2x2 matrix of 
  /// self-energies in LR basis
  DoubleMatrix addSdownHiggs(double p1, double p2, int family, 
  DoubleMatrix & higgs); 
  /// LCT: Adds one-loop electroweak corrections to sdown mass matrix
  /// IO parameters: {p1,p2}=external momentum, outputs electroweak=2x2 matrix 
  /// of self-energies in LR basis
  DoubleMatrix addSdownEweak(double p1, double p2, int family, 
  DoubleMatrix & electroweak); 
  /// LCT: Adds one-loop neutralino corrections to sdown mass matrix
  /// IO parameters: {p1,p2}=external momentum, outputs neutralino=2x2 matrix 
  /// of self-energies in LR basis
  DoubleMatrix addSdownNeutralino(double p1, double p2, int family, 
  DoubleMatrix & neutralino);
  /// LCT: Adds one-loop chargino corrections to sdown mass matrix
  /// IO parameters: {p1,p2}=external momentum, outputs chargino=2x2 matrix of 
  /// self-energies in LR basis
  DoubleMatrix addSdownChargino(double p1, double p2, int family, 
  DoubleMatrix & chargino); 
  /// Adds one-loop corrections to sbottom mass matrix at p=root(mb1 mb2)
  /// IO parameters: mass=tree level mass matrix on
  /// input, is returned with radiative corrections added, mt=DR bar top mass
  virtual void addSdownCorrection(DoubleMatrix & mass, int family);
  /// LCT: Adds one-loop strong corrections to sbottom mass matrix
  /// IO parameters: p=external momentum, mt=DR bar top mass, 
  /// outputs strong=2x2 matrix of self-energies in LR basis
  DoubleMatrix addSbotQCD(double p, double mt, DoubleMatrix & strong);
  /// LCT: Adds one-loop stop, sbottom and stau corrections to sbottom mass 
  /// matrix.  IO parameters: p=external momentum, mt=DR bar top mass, 
  /// outputs {stop,sbottom}=2x2 matrix of self-energies in LR basis
  void addSbotSfermion(double p, double mt, DoubleMatrix & stop, 
  DoubleMatrix & sbottom);
  /// LCT: Adds one-loop Higgs corrections to sbottom mass matrix
  /// IO parameters: p=external momentum, mt=DR bar top mass, 
  /// outputs higgs=2x2 matrix of self-energies in LR basis
  DoubleMatrix addSbotHiggs(double p, double mt, DoubleMatrix & higgs);
  /// LCT: Adds one-loop electroweak corrections to sbottom mass matrix
  /// IO parameters: p=external momentum, outputs electroweak=2x2 matrix of 
  /// self-energies in LR basis
  DoubleMatrix addSbotEweak(double p, DoubleMatrix & electroweak);
  /// LCT: Adds one-loop neutralino corrections to sbottom mass matrix
  /// IO parameters: p=external momentum, outputs neutralino=2x2 matrix of 
  /// self-energies in LR basis
  DoubleMatrix addSbotNeutralino(double p, double mt, DoubleMatrix & neutralino);
  /// LCT: Adds one-loop chargino corrections to sbottom mass matrix
  /// IO parameters: p=external momentum, outputs chargino=2x2 matrix of 
  /// self-energies in LR basis
  DoubleMatrix addSbotChargino(double p, double mt, DoubleMatrix & chargino);
  /// Adds one-loop corrections to sbottom mass matrix at p=root(mb1 mb2)
  /// IO parameters: p=external momentum scale, mass=tree level mass matrix on
  /// input, is returned with radiative corrections added, mt=DR bar top mass
  virtual void addSbotCorrection(double p, DoubleMatrix & mass, double mb);
  /// LCT: Adds one-loop Higgs corrections to slepton self energy
  /// IO parameters: p=external momentum, family=generation, 
  /// outputs higgs=2x2 matrix of self-energies
  DoubleMatrix addSlepHiggs(double p1, double p2, int family, DoubleMatrix & higgs);
  /// LCT: Adds one-loop electroweak corrections to slepton self energy
  /// IO parameters: p=external momentum, family=generation, 
  /// outputs electroweak=2x2 matrix of self-energies
  DoubleMatrix addSlepEweak(double p1, double p2, int family, DoubleMatrix & electroweak);
  /// LCT: Adds one-loop gaugino corrections to slepton self energy
  /// IO parameters: p=external momentum, family=generation, 
  /// outputs {chargino,neutralino}=2x2 matrix of self-energies
  void addSlepGaugino(double p1, double p2, int family, DoubleMatrix & chargino, DoubleMatrix & neutralino);
  /// Adds one-loop corrections to sel_fam mass matrix at p=root(msel1 msel2)
  /// IO parameters: mass=tree level mass on
  /// input, is returned with radiative corrections added, mt=DR bar top mass
  virtual void addSlepCorrection(DoubleMatrix & mass, int family);
  /// This is a special version that I've added to correct lower families'
  /// mixings too
  void addSlepCorrection(double p, DoubleMatrix & mass, int family);
  /// LCT: Adds one-loop sfermion corrections to stau self energy
  /// IO parameters: p=external momentum, mtau=DR bar tau mass,
  /// outputs {stop,sbottom}=2x2 matrix of self-energies
  void addStauSfermion(double p, double mtau, DoubleMatrix & stop, DoubleMatrix & sbottom);
  /// LCT: Adds one-loop Higgs corrections to stau self energy
  /// IO parameters: p=external momentum, mtau=DR bar tau mass,
  /// outputs higgs=2x2 matrix of self-energies
  void addStauGaugino(double p, double mtau, DoubleMatrix & chargino, DoubleMatrix & neutralino);
  /// LCT: Adds one-loop gaugino corrections to stau self energy
  /// IO parameters: p=external momentum, mtau=DR bar tau mass,
  /// outputs {chargino,neutralino}=2x2 matrix of self-energies
  DoubleMatrix addStauEweak(double p, double mtau, DoubleMatrix & electroweak);
  /// LCT: Adds one-loop electroweak corrections to stau self energy
  /// IO parameters: p=external momentum, mtau=DR bar tau mass,
  /// outputs electroweak=2x2 matrix of self-energies
  DoubleMatrix addStauHiggs(double p, double mtau, DoubleMatrix & higgs);
  /// Adds one-loop corrections to stau mass matrix at p=root(mtau1 mtau2)
  /// IO parameters: mass=tree level mass on
  /// input, is returned with radiative corrections added, mt=DR bar top mass
  virtual void addStauCorrection(double p, DoubleMatrix & mass, double mtau);
  /// LCT: Adds one-loop strong corrections to sup mass matrix
  /// IO parameters: {p1,p2}=external momentum, outputs strong=2x2 matrix of 
  /// self-energies in LR basis
  DoubleMatrix addSupQCD(double p1, double p2, int family, 
  DoubleMatrix & strong); 
  /// LCT: Adds one-loop Higgs corrections to sup mass matrix
  /// IO parameters: {p1,p2}=external momentum, outputs higgs=2x2 matrix of 
  /// self-energies in LR basis
  DoubleMatrix addSupHiggs(double p1, double p2, int family, 
  DoubleMatrix & higgs); 
  /// LCT: Adds one-loop electroweak corrections to sup mass matrix
  /// IO parameters: {p1,p2}=external momentum, outputs electroweak=2x2 matrix 
  /// of self-energies in LR basis
  DoubleMatrix addSupEweak(double p1, double p2, int family, 
  DoubleMatrix & electroweak); 
  /// LCT: Adds one-loop neutralino corrections to sup mass matrix
  /// IO parameters: {p1,p2}=external momentum, outputs neutralino=2x2 matrix 
  /// of self-energies in LR basis
  DoubleMatrix addSupNeutralino(double p1, double p2, int family, 
  DoubleMatrix & neutralino); 
  /// LCT: Adds one-loop chargino corrections to sup mass matrix
  /// IO parameters: {p1,p2}=external momentum, outputs chargino=2x2 matrix of 
  /// self-energies in LR basis
  DoubleMatrix addSupChargino(double p1, double p2, int family, 
  DoubleMatrix & chargino); 
  /// Adds one-loop corrections to sup mass matrix at p=root(msup1 msup2)
  /// IO parameters: mass=tree level mass on
  /// input, is returned with radiative corrections added, mt=DR bar top mass
  virtual void addSupCorrection(DoubleMatrix & mass, int family);
  /// LCT: Adds one-loop sfermion corrections to tau sneutrino self energy
  /// IO parameters: p=external momentum, outputs {stop,sbottom}=self-energy
  void addSnuTauSfermion(double p, double & stop, double & sbottom);
  /// LCT: Adds one-loop Higgs corrections to tau sneutrino self energy
  /// IO parameters: p=external momentum, outputs higgs=self-energy
  double addSnuTauHiggs(double p, double & higgs);
  /// LCT: Adds one-loop electroweak corrections to tau sneutrino self energy
  /// IO parameters: p=external momentum, outputs electroweak=self-energy
  double addSnuTauEweak(double p, double & electroweak);
  /// LCT: Adds one-loop gaugino corrections to tau sneutrino self energy
  /// IO parameters: p=external momentum, outputs 
  /// {chargino,neutralino}=self-energy
  void addSnuTauGaugino(double p, double & chargino, double & neutralino);
  /// Adds one-loop corrections to tau sneutrino mass 
  /// IO parameters: p=external momentum, mass=tree level mass matrix on
  /// input, is returned with radiative corrections added, mt=DR bar top mass
  virtual void addSnuTauCorrection(double & mass);
  /// LCT: Adds one-loop Higgs corrections to sneutrino self energy
  /// IO parameters: p=external momentum, family=generation, 
  /// outputs higgs=self-energy
  double addSnuHiggs(double p, int family, double & higgs);
  /// LCT: Adds one-loop electroweak corrections to sneutrino self energy
  /// IO parameters: p=external momentum, family=generation,
  /// outputs electroweak=self-energy
  double addSnuEweak(double p, int family, double & electroweak);
  /// LCT: Adds one-loop gaugino corrections to sneutrino self energy
  /// IO parameters: p=external momentum, outputs 
  /// {chargino,neutralino}=self-energy
  void addSnuGaugino(double p, int family, double & chargino, double & neutralino);
  /// Adds one-loop corrections to sneutrino mass of family "family"
  /// IO parameters: mass=tree level mass matrix on
  /// input, is returned with radiative corrections added
  virtual void addSnuCorrection(double & mass, int family);
  /// Adds approximate one-loop corrections to squark mass matrix for first
  /// two families.
  /// IO parameters: mass=tree level mass matrix on
  /// input, is returned with radiative corrections added
  void addSquarkCorrection(DoubleMatrix & mass);
  /// Organises calculation of all up squark masses.
  /// IO parameters: mt=DRbar top mass, pizztMS=Z self energy at Q=M_SUSY,
  /// sinthDRbarMS=DRbar value of sin theta_w, accuracy=0,1: number of loops
  /// to add to tree-level squark mass matrix
  virtual void doUpSquarks(double mt, double pizztMS, double sinthDRbarMS, int
		      accuracy); 
  /// Organises calculation of all down squark masses.
  /// IO parameters: mb=DRbar bottom mass, pizztMS=Z self energy at Q=M_SUSY,
  /// sinthDRbarMS=DRbar value of sin theta_w, accuracy=0,1: number of loops
  /// to add to tree-level squark mass matrix
  virtual void doDownSquarks(double mb, double pizztMS, double sinthDRbarMS, 
			     int accuracy, double mt);
  /// Organises calculation of all slepton masses.
  /// IO parameters: mT=DRbar tau mass, pizztMS=Z self energy at Q=M_SUSY,
  /// sinthDRbarMS=DRbar value of sin theta_w, accuracy=0,1: number of loops
  /// to add to tree-level squark mass matrix
  virtual void doChargedSleptons(double mT, double pizztMS, 
				 double sinthDRbarMS, int accuracy);
  /// Organises calculation of all sneutrino masses, pizztMS=Z self energy at
  /// Q=M_SUSY
  virtual void doSnu(double pizztMS, int accuracy = 0);
  /// Returns tree-level up squark mass matrix in "mass".
  /// IO parameters: mass=tree level mass matrix on
  /// input, is returned with radiative corrections added, mtrun=DR bar top
  /// mass, family=generation of squark, pizztMS=Z self energy at Q=M_SUSY,
  /// sinthDRbarMS=DRbar value of sin theta_w 
  virtual void treeUpSquark(DoubleMatrix & mass, double mtrun, double pizztMS, 
                            double sinthDRbarMS, int family);
  /// Returns tree-level down squark mass matrix in "mass".
  /// IO parameters: mass=tree level mass matrix on
  /// input, is returned with radiative corrections added, mbrun=DR bar bottom
  /// mass, family=generation of squark, pizztMS=Z self energy at Q=M_SUSY,
  /// sinthDRbarMS=DRbar value of sin theta_w 
  virtual void treeDownSquark(DoubleMatrix & mass, double mbrun, double pizztMS,
                              double sinthDRbarMS, int family);
  /// Returns tree-level down squark mass matrix in "mass".
  /// IO parameters: mass=tree level mass matrix on
  /// input, is returned with radiative corrections added, mTrun=DR bar tau
  /// mass, family=generation of slepton, pizztMS=Z self energy at Q=M_SUSY,
  /// sinthDRbarMS=DRbar value of sin theta_w
  virtual void treeChargedSlepton(DoubleMatrix & mass, double mTrun, 
                                  double pizztMS, double sinthDRbarMS, 
                                  int family);
  /// Organises calculation of all sneutrino masses, pizztMS=Z self energy at
  /// Q=M_SUSY, mSq=mass of sneutrino, family=generation of sneutrino
  void treeSnu(double & mSq, double pizztMS, int family);

  /// Organises tree-level calculation of all sparticle masses and mixings
  virtual void calcDrBarPars();
  /// calculates the higgs DRbar parameters. Make sure mt is set in eg. It
  /// will fill in the Higgs masses with the appropriate values on exit.
  void setNeutCurrCouplings(double sinthDRbar, double & sw2, double & guL, double & gdL, double & geL, double & guR, double & gdR, double & geR );
//PA: sets the Yukawas and Trilinears
  void calcDRTrilinears(drBarPars & eg, double vev, double beta);

  void calcDrBarHiggs(double beta, double mz2, double mw2, 
		      double sinthDRbar, drBarPars & eg);
  /// calculates the chargino DRbar parameters. It will fill in the chargino
  /// masses in eg with the appropriate values on exit. 
  void treeCharginos(DoubleMatrix & mCh, double beta, double mw);
  /// calculates the chargino DRbar parameters. It will fill in the chargino
  /// masses in eg with the appropriate values on exit. 
  void treeNeutralinos(DoubleMatrix & mN, double beta, double mz, 
                            double mw, double sinth);
  // calculates the chargino and neutralino DRbar parameters. 
  //It will fill in the chargino and neutralino
  /// masses in eg with the appropriate values on exit. 
  void calcDrBarGauginos(double beta, double mw, double mz, double sinth, 
                         drBarPars & eg);
  /// For an input tan beta=tb, sets gauge and Yukawa couplings according to
  /// the tree-level spectrum and data set: pars contains the boundary
  /// conditions. They aren't used in R-parity conservation, though. 
  virtual void sparticleThresholdCorrections(double tb);
  //					     const DoubleVector & pars);
  /// Does SUSY (and other) threshold corrections to alphaEm - returns alpha in
  /// DRbar scheme at scale Q. From hep-ph/9606211. Input empirical value of
  /// alpha at 0 momentum....
  double qedSusythresh(double alphaEm, double Q) const;
  /// Does SUSY (and other) threshold corrections to alphaS
  /// Input alphas in MSbar and it returns it in DRbar scheme. 
  /// From hep-ph/9606211
  double qcdSusythresh(double alphasMSbar, double Q) const;
  /// Calculates the best scale at which to do symmetry breaking:
  /// \f$ M_{SUSY}=\sqrt{m_{{\tilde t}_1 {\tilde t}_2}} \f$. Should only be
  /// called after calcDrBarPars.
  double calcMs() const;
  /// Calculates physical sparticle masses to accuracy number of loops. Should
  /// be called at M_{SUSY}.
  virtual void physical(int accuracy);
   /// Obtains (16 Pi^2) / mt  times 1-loop and 2-loop QCD corrections 
  //to the pole mt for use in calcRunningMt 
  virtual double calcRunMtQCD() const;
  /// Obtains (16 Pi^2) / mt  times 1-loop stop-Gluino corrections 
  //to the pole mt for use in calcRunningMt 
  virtual double calcRunMtStopGluino() const;
  /// Obtains (16 Pi^2) / mt  times 1-loop Higgs corrections 
  //to the pole mt for use in calcRunningMt 
  virtual double calcRunMtHiggs() const;

  /// Obtains (16 Pi^2) / mt  times 1-loop Neutralino corrections 
  //to the pole mt for use in calcRunningMt 
  virtual double calcRunMtNeutralinos() const;
  /// Obtains (16 Pi^2) / mt  times 1-loop Charginos corrections 
  //to the pole mt for use in calcRunningMt 
  virtual double calcRunMtCharginos() const;

  /// Applies 1-loop SUSY and 2-loop QCD corrections to pole mt in order to
  /// return the DRbar running value at the current scale
  virtual double calcRunningMt();
  //PA: calculates factor to convert to DrBar for mtau
  virtual double calcRunMtauDrBarConv() const;
  // Obtains (1 / mTAU)  times 1-loop squark-chargino corrections 
  //to mtau for use in calcRunningMtau 
  virtual double calcRunMtauCharginos(double mTauSMMZ) const;
  // Obtains (1 / mTAU)  times 1-loop squark-neutralino corrections 
  //to mtau for use in calcRunningMtau 
  virtual double calcRunMtauNeutralinos(double mTauSMMZ) const; 
  /// Obtains (1 / mtau)  times 1-loop Higgs corrections 
  // to mtau for use in calcRunningMtau 
  virtual double calcRunMtauHiggs() const;
  /// Applies 1-loop SUSY corrections to pole mtau in order to
  /// return the DRbar running value at the current scale
  virtual double calcRunningMtau() const;  
 
  //PA: calculates factor to convert to DrBar for mb
  virtual double calcRunMbDrBarConv() const;
  /// Obtains (1 / mb)  times 1-loop sbottom-Gluino corrections 
  //to mb for use in calcRunningMb 
  virtual double calcRunMbSquarkGluino() const;
  /// Obtains (1 / mb)  times 1-loop squark-chargino corrections 
  //to mb for use in calcRunningMb 
  virtual double calcRunMbChargino() const;
  /// Obtains (1 / mb)  times 1-loop Higgs corrections 
  // to mb for use in calcRunningMb 
  virtual double calcRunMbHiggs() const;
   /// Obtains (1 / mb)  times 1-loop neutralino-suqrk corrections 
  // to mb for use in calcRunningMb 
  virtual double calcRunMbNeutralinos() const;
  /// Applies approximate 1-loop SUSY corrections to mb(MZ) in order to
  /// return the DRbar running value
  virtual double calcRunningMb() const;
  /*
  /// Calculates top Yukawa coupling, supply Higgs vev parameter at current
  /// scale 
  double calcHt(double vev);
  /// Calculates bottom Yukawa coupling, supply Higgs vev parameter at current
  /// scale 
  double calcHb(double vev) const;
  /// Calculates tau Yukawa coupling, supply Higgs vev parameter at current
  /// scale 
  double calcHtau(double vev) const;
  */
  /// Calculates DRbar sin theta_w at the current scale from gauge couplings 
  double calcSinthdrbar() const;
  /// Calculates Higgs VEV parameter from gauge couplings and MZ
  double getVev();
  /// Input for this one (saves time, possibly) is to give the self-energy of
  /// the Z at the current scale
  double getVev(double pizzt);
  /// Calculates pole chargino masses and mixing using approximate 1-loop SUSY
  /// corrections. IO parameters: piwwt is the W self-energy at the current,
  /// accuracy is the number of loops required (0 or 1 currently)
  virtual void charginos(int accuracy, double piwwt);
   ///PA: Fills sigmaL, sigmaR and sigmaS with sfermion parts of the charginos 
  //loop corrections.
  virtual void addChaLoopSfermion(double p, DoubleMatrix & sigmaL, DoubleMatrix & sigmaR, DoubleMatrix & sigmaS) const;
   ///PA: Fills sigmaL, sigmaR and sigmaS with Gauge parts of the charginos 
  //loop corrections.
  virtual void addChaLoopGauge(double p, DoubleMatrix & sigmaL, DoubleMatrix & sigmaR, DoubleMatrix & sigmaS, DoubleMatrix b1pCha, DoubleMatrix b0pCha, DoubleMatrix b1pNeut, DoubleMatrix b0pNeut) const;
  ///PA: Fills sigmaL, sigmaR and sigmaS with Higgs parts of the charginos 
  //loop corrections.
  virtual void addChaLoopHiggs(double p, DoubleMatrix & sigmaL, DoubleMatrix & sigmaR, DoubleMatrix & sigmaS, DoubleMatrix b1pCha, DoubleMatrix b0pCha, DoubleMatrix b1pNeut, DoubleMatrix b0pNeut) const;
  /// Adds the loop corrections on to an input tree-level chargino mass
  virtual void addCharginoLoop(double p, DoubleMatrix &);
  /// Calculates pole neutralino masses and mixingusing approximate 1-loop SUSY
  /// corrections. IO parameters: piwwt is the W self-energy at M_SUSY,
  /// accuracy is the number of loops required (0 or 1 currently), pizzt is
  /// the Z self-energy at M_SUSY
  virtual void neutralinos(int accuracy, double piwwt, double pizzt);
  /// LCT: Adds sfermion loop corrections to neutralino mass matrix.
  /// IO parameters: p=external momentum, {sigmaL,sigmaR,sigmaS}= 4x4 matrices 
  /// of {left,right,scalar} corrections
  void addNeutLoopSfermion(double p, DoubleMatrix & sigmaL, 
  DoubleMatrix & sigmaR, DoubleMatrix & sigmaS);
  /// LCT: Returns matrix of Passarino-Veltman B0 and B1 functions.
  void getNeutPassarinoVeltman(double p, double q, DoubleMatrix & b0fn, 
  DoubleMatrix & b1fn); 
  /// LCT: Adds gauge loop corrections to neutralino mass matrix.
  /// IO parameters: p=external momentum, {sigmaL,sigmaR,sigmaS}= 4x4 matrices 
  /// of {left,right,scalar} corrections
  void addNeutLoopGauge(double p, DoubleMatrix & sigmaL, 
  DoubleMatrix & sigmaR, DoubleMatrix & sigmaS);
  /// LCT: Adds Higgs loop corrections to neutralino mass matrix.
  /// IO parameters: p=external momentum, {sigmaL,sigmaR,sigmaS}= 4x4 matrices 
  /// of {left,right,scalar} corrections
  void addNeutLoopHiggs(double p, DoubleMatrix & sigmaL, 
  DoubleMatrix & sigmaR, DoubleMatrix & sigmaS);
  /// Adds the loop corrections on to an input tree-level neutralino mass
  virtual void addNeutralinoLoop(double p, DoubleMatrix &);
  /// Calculates pole gluino mass to 1-loop SUSY corrections
  virtual void gluino(int accuracy);
  /// Calculates pole Higgs masses and mixings: full 1-loop SUSY corrections
  /// and 2-loop alpha_t (alpha_s + alpha_t) + alpha_s alpha_b effective
  /// potential corrections. 
  /// IO parameters: piwwt is the W self-energy at M_SUSY, accuracy is number
  /// of loops (0 or 1) to use and pizzt is the Z self-energy at M_SUSY
  /// Returns "true" if there's a tachyon problem
  bool higgs(int accuracy, double piwwt, double pizzt);
  /// Calculates pole Higgs masses and mixings: full 1-loop SUSY corrections
  /// and 2-loop alpha_t (alpha_s + alpha_t) + alpha_s alpha_b effective
  /// potential corrections. 
  /// IO parameters: piwwt is the W self-energy at M_SUSY, accuracy is number
  /// of loops (0 or 1) to use and pizzt is the Z self-energy at M_SUSY
  //  virtual void newhiggs(int accuracy, double piwwt, double pizzt);
  /// Tree-level REWSB calculation, returning mu at correct value. sgnMu is the
  /// required sign (+/- 1). Returns 1 if mu^2<0, indicating an inconsistent
  /// minimum 
  virtual int rewsbMu(int sgnMu, double & mu) const;
  /// Tree-level REWSB calculation, returning m3sq at correct value consistent
  /// with mu
  virtual int rewsbM3sq(double, double &) const;
  /// Fixes mH1^2 and mH2^2 in order to get muCond and mAcond correct
  void alternativeEwsb(double mt);
  /// Organises high accuracy rewsb: call it at the low scale M_{SUSY}
  /// IO parameters: sgnMu is +/-1 (desired sign of mu), mt is DRbar top mass
  /// muOld is a value of mu with which the current one is to be averaged, if
  /// it is set to some value above the number of the beast. This can be
  /// useful in attaining convergence in regions where it is difficult.
  /// eps reflects how the old and new values of mu are to be averaged:
  /// it's (eps * mu + (1-eps) * muOld)
  virtual void rewsb(int sgnMu, double mt, const DoubleVector & pars,
		     double muOld = -6.66e66, double eps = 0.5);
  /// Organises tree-level rewsb: call it at the low scale M_{SUSY}
  /// IO parameters: sgnMu is +/-1 (desired sign of mu)
  virtual void rewsbTreeLevel(int sgnMu);
  /// This is the value you'd get from minimising the Higgs potential at
  /// tree-level 
  double treeLevelMuSq();
  /// Obtains solution of one-loop effective potential minimisation via
  /// iteration technique. Currently includes: all 1-loop SUSY tadpoles, plus
  /// 2-loop alpha_t (alpha_t + alpha_s) + alpha_b alpha_s corrections
  /// IO parameters: 
  /// munew=current value of mu for this iteration, sgnMu=desired sign of mu,
  /// mt=DRbar mtop, maxTries=maximum number of iterations before it bails
  /// out, pizzMS=self-energy of MZ at current scale, sinthDRbar=DRbar value
  /// of sin theta_w, tol=desired fractional accuracy on mu, err=error flag:
  /// err=1 if no iteration reached, 2 if incorrect rewsb
  void iterateMu(double & munew, int sgnMu, double mt, 
		  int maxTries, double pizztMS, double sinthDRbar, double tol,
		 int  & err);
  /// This is a check: predicts tan beta from the values of soft parameters
  /// and mu that we have
  double predTanb(double muSusy = -6.66e66) const;
  /// Predicts value of MZ(pole) from values of soft parameters and mu that we
  /// have. tanb=tan beta is also predicted
  double predMzsq(double & tanb, double muOld = -6.66e66, double eps = 0.);
  /// Calculates fine-tuning for soft parameters and mu, m_3^2, top Yukawa. 
  /// IO parameters: bcPars 
  /// should be a vector giving the high-scale SUSY breaking boundary
  /// condition parameters, MX is the high-scale, boundaryCondition is the
  /// user-supplied function that sets the SUSY breaking BCs.
  /// If doTop is true, it also calculates the fine tuning associated with the
  /// top Yukawa coupling. 
  DoubleVector fineTune(void (*boundaryCondition)(Softsusy &, 
						  const DoubleVector &), 
			const DoubleVector & bcPars, double MX, 
			bool doTop = false);
  /// Give it a GUT scale object consistent with rewsb
  /// and it'll return the fine tuning of one parameter specified by numPar
  /// ht, mu  and m3sq at the high 
  double it1par(int numPar, const DoubleVector & bcPars);
  /// Input mx the scale up to which you search for minima.
  /// Returns minimum value of potential along UFB3 direction.
  /// Does ufbs truly properly but takes a long time.
  double ufb3sl(double);
  /// You should evaluate this at a scale MSusy average of stops.
  /// Returns depth of electroweak minimum
  double realMinMs() const;
  
  /// Calculates Higgs contribution to the transverse part of Z self-energy: 
  //for p=external momentum, Q=renormalisation scale
  virtual double piZZTHiggs(double p, double Q, double thetaWDRbar) const;
  /// Calculates sfermion contribution to the transverse part of Z self-energy: 
  //for p=external momentum, Q=renormalisation scale
  virtual double piZZTsfermions(double p, double Q) const;
  /// Calculates fermion contribution to the transverse part of Z self-energy: 
  //for p=external momentum, Q=renormalisation scale
  virtual double piZZTfermions(double p, double Q, bool usePoleMt) const;
 /// Calculates neutralino contrib. to the transverse part of Z self-energy: 
  //for p=external momentum, Q=renormalisation scale
  virtual double piZZTNeutralinos(double p, double Q, double thetaWDRbar) const;
   /// Calculates chargino contrib. to the transverse part of Z self-energy: 
  //for p=external momentum, Q=renormalisation scale
  virtual double piZZTCharginos(double p, double q, double thetaWDRbar) const;
  /// Calculates transverse part of Z self-energy: for p=external momentum,
  /// Q=renormalisation scale
  virtual double piZZT(double p, double Q, bool usePoleMt = false) const;
   /// Calculates Higgs contribution to the transverse part of W self-energy: 
  //for p=external momentum, Q=renormalisation scale
  virtual double piWWTHiggs(double p, double q, double thetaWDRbar) const;
  /// Calculates fermion contribution to the transverse part of W self-energy: 
  //for p=external momentum, Q=renormalisation scale
  virtual double piWWTfermions(double p, double Q, bool usePoleMt) const;
  /// Calculates sfermion contribution to the transverse part of W self-energy: 
  //for p=external momentum, Q=renormalisation scale
  virtual double piWWTsfermions(double p, double Q) const;
  /// Calculates transverse part of W self-energy: for p=external momentum,
  /// Q=renormalisation scale
  virtual double piWWTgauginos(double p, double Q, double thetaWDRbar) const;
  /// Calculates transverse part of W self-energy: for p=external momentum,
  /// Q=renormalisation scale
  virtual double piWWT(double p, double Q, bool usePoleMt = false) const;
  /// LCT: Returns neutralino-chargino-charged-Higgs trilinear couplings in 
  /// weak basis
  virtual void getNeutralinoCharginoHpmCoup(ComplexMatrix & apph1, ComplexMatrix & apph2, ComplexMatrix & bpph1, ComplexMatrix & bpph2) const;
  /// LCT: Calculates (16 Pi^2) times the fermion contribution to H^+H^- 
  /// self-energy for p=external momentum, q=renormalisation scale
  double piHpHmFermions(double p, double q) const;
  /// LCT: Calculates (16 Pi^2) times the sfermion contribution to H^+H^-
  /// self-energy for p=external momentum, q=renormalisation scale, 
  /// and SUSY parameter mu
  double piHpHmSfermions(double p, double q, double mu) const;
  /// LCT: Calculates (16 Pi^2) times the gauge contribution to H^+H^-
  /// self-energy for p=external momentum, q=renormalisation scale
  double piHpHmGauge(double p, double q) const;
  /// LCT: Calculates (16 Pi^2) times the Higgs contribution to H^+H^-
  /// self-energy for p=external momentum, q=renormalisation scale
  double piHpHmHiggs(double p, double q) const;
  /// LCT: Calculates (16 Pi^2) times the gaugino contribution to H^+H^-
  /// self-energy for p=external momentum, q=renormalisation scale
  double piHpHmGauginos(double p, double q) const;
  /// Calculates transverse part of H^+H^- self-energy:  
  /// for p=external momentum, Q=renormalisation scale
  virtual double piHpHm(double p, double Q) const;
  /// Calculates Z gamma self-energy:  
  /// for p=external momentum, Q=renormalisation scale
  double piZGT(double p, double Q) const;
  /// Calculates transverse part of A^0 self-energy: for p=external momentum,
  /// Q=renormalisation scale
  virtual double piAA(double p, double Q) const;
  //PA:Calculates (16 Pi^2) times the sfermion contribution to Higgs 
  //self-energy: for p=external momentum, q=renormalisation scale
  double pis1s1Sfermions(double p, double q,  DoubleMatrix ls1tt,  DoubleMatrix ls1bb,  DoubleMatrix ls1tautau) const;
 //PA:Calculates (16 Pi^2) times the sfermion contribution to Higgs 
  //self-energy: for p=external momentum, q=renormalisation scale
  double pis1s2Sfermions(double p, double q,  DoubleMatrix ls1tt,  DoubleMatrix ls1bb,  DoubleMatrix ls1tautau, DoubleMatrix ls2tt,  DoubleMatrix ls2bb,  DoubleMatrix ls2tautau) const;
//PA:Calculates (16 Pi^2) times the sfermion contribution to Higgs 
  //self-energy: for p=external momentum, q=renormalisation scale
  double pis2s2Sfermions(double p, double q, DoubleMatrix ls2tt,  DoubleMatrix ls2bb,  DoubleMatrix ls2tautau) const;
 //PA:Calculates (16 Pi^2) times the fermion contribution to Higgs self-energy:
  //for p=external momentum, q=renormalisation scale
  double pis1s1Fermions(double p, double q) const;
//PA:Calculates (16 Pi^2) times the fermion contribution to Higgs self-energy:
  //for p=external momentum, q=renormalisation scale
  double pis2s2Fermions(double p, double q) const;
//PA:Calculates (16 Pi^2) times the Higgs contribution to Higgs self-energy: 
//for p=external momentum, q=renormalisation scale
  double pis1s1Higgs(double p, double q) const;
//PA:Calculates (16 Pi^2) times the Higgs contribution to Higgs self-energy: 
//for p=external momentum, q=renormalisation scale
  double pis1s2Higgs(double p, double q) const;
//PA:Calculates (16 Pi^2) times the Higgs contribution to Higgs self-energy: 
//for p=external momentum, q=renormalisation scale
  double pis2s2Higgs(double p, double q) const;
  //PA:Calculates (16 Pi^2) times the Neutralino contribution to the
  // Higgs self-energy: for p=external momentum, q=renormalisation scale
  double pis1s1Neutralinos(double p, double q) const;
  //PA:Calculates (16 Pi^2) times the Neutralino contribution to the
  // Higgs self-energy: for p=external momentum, q=renormalisation scale
  double pis1s2Neutralinos(double p, double q) const;
  //PA:Calculates (16 Pi^2) times the Neutralino contribution to the
  // Higgs self-energy: for p=external momentum, q=renormalisation scale
  double pis2s2Neutralinos(double p, double q) const;
  //PA:Calculates (16 Pi^2) times the Chargino contribution to the
  // Higgs self-energy: for p=external momentum, q=renormalisation scale
  double pis1s1Charginos(double p, double q) const;
   //PA:Calculates (16 Pi^2) times the Chargino contribution to the
  // Higgs self-energy: for p=external momentum, q=renormalisation scale
  double pis1s2Charginos(double p, double q) const;
  //PA:Calculates (16 Pi^2) times the Chargino contribution to the
  // Higgs self-energy: for p=external momentum, q=renormalisation scale
  double pis2s2Charginos(double p, double q) const;
/// Calculates transverse part of Higgs self-energy: for p=external momentum,
  /// Q=renormalisation scale
  double pis1s1(double p, double q) const;
  /// Calculates transverse part of Higgs self-energy: for p=external momentum,
  /// Q=renormalisation scale
  double pis1s2(double p, double q) const;
  /// Calculates transverse part of Higgs self-energy: for p=external momentum,
  /// Q=renormalisation scale
  double pis2s2(double p, double q) const;
  /// Calculates sin^2 theta^l_eff
  double sinSqThetaEff();
  //PA: gets h1 mixing element with Hu.
  virtual double h1s2Mix();
  /// Iterative determination of rho parameter consistent with muon decay
  /// constant, MZ and alpha_0. 
  /// IO parameters: outrho and outsin are the current DRbar values of sin
  /// theta_w and the rho parameter. alphaMZDRbar=already calculated value of
  /// fine-structure constant, pizztMZ=Z self-energy at Q=MZ, tol=desired
  /// piwwt0=self-energy of the W at p=0, piwwtMW=self-energy of the W at p=MW
  /// accuracy and maxTries is the maximum number of iterations the routines
  /// will allow
  virtual void rhohat(double & outrho, double & outsin, double alphaMZDRbar,
		      double pizztMZ, double piwwt0, double piwwtMW, 
		      double tol, int maxTries);
  /// Calculates delta_v corrections for outrho=DRbar rho parameter,
  /// outsin=DRbar sin theta_w, alphaDRbar=alpha(Q) in the DR bar scheme,
  /// pizztMZ=self-energy of the Z at MZ
  virtual double deltaVb(double outrho, double outsin, double alphaDRbar, 
		double pizztMZ) const;
  /// Calculates delta rho corrections for outrho=DRbar rho parameter,
  /// outsin=DRbar sin theta_w, alphaDRbar=alpha(Q) in the DR bar scheme,
  /// pizztMZ=self-energy of the Z at MZ, piwwtMW=self-energy of the W at p=MW
  virtual double dRho(double outrho, double outsin, double alphaDRbar, 
		 double pizztMZ, double piwwtMW);
  /// Calculates delta r corrections for outrho=DRbar rho parameter,
  /// outsin=DRbar sin theta_w, alphaDRbar=alpha(Q) in the DR bar scheme,
  /// pizztMZ=self-energy of the Z at p=MZ, pizzt0=self-energy of the W at p=0
  virtual double dR(double outrho, double outsin, double alphaDRbar, 
                    double pizztMZ, double piwwt0);
  /// Returns the mass of the heaviest SUSY particle, excluding gravitino
  double maxMass() const;
  /// Returns lsp mass in mass and function return labels which particle is 
  /// lsp=0 means LSP is neutralino, 1=up squark, 2=down squark, 3=sleptons,
  /// 4=charginos, 5=sneutrinos, 6=gluino. posi, posj give the 
  /// "handedness" (ie 1 or 2) for scalars and family respectively. 
  int lsp(double & mass, int & posi, int & posj) const;
  /// Returns nlsp mass in mass and function return labels which particle is 
  /// nlsp=0 means NLSP is neutralino, 1=up squark, 2=down squark, 3=sleptons,
  /// 4=charginos, 5=sneutrinos, 6=gluino. posi, posj give the 
  /// "handedness" (ie 1 or 2) for scalars and family respectively.
  int nlsp(double & mass, int & posi, int & posj) const;
  
  /// Prints a list of important sparticle/Higgs masses to standard output
  std::string printShort() const;
  /// Prints a list of all sparticle/Higgs masses to standard output
  std::string printLong();

  /// Prints whols object to standard output
  virtual void printObj() { cout << *this; };
  
  /// log(max(a^2, b^2, c^2) / Q^2)
  double thet(double a, double b, double c = 0.0);
 
  /// Driver calculation to determine all sparticle masses and parameters.
  /// Returns low energy softsusy object consistent with BC's m0 etc at MGUT.
  /// oneset should be at MZ and contains the SM data to fit the model to.
  /// If the running comes into difficulty, eg if a Landau pole is reached, it
  /// returns a ZERO object: no result is possible!
  /// Boundary condition is the theoretical condition on parameters at the high
  /// energy scale mx: the parameters themselves are contained within the
  /// vector. IO parameters:
  /// sgnMu is the desired sign of mu: + or - 1. If mu is 0, mu is set
  /// initially as a boundary condition. tanb = desired value of DR bar tan
  /// beta(MZ). boundaryCondition is the function which sets to SUSY BCs at
  /// the high scale, mxGuess is the GUT scale
  /// gaugeUnification=true if the user requests true gauge unification, in
  /// which case mxGuess is the first guess. Returns actual GUT scale.
  /// ewsbBCscale = false (or omitted) means that the boundary condition on
  /// SUSY breaking is set in the usual way. If it is true, the boundary
  /// condition is set to \f$\sqrt{m_{{\tilde t}_1} m_{{\tilde t}_2}} \f$, ie
  /// like in the "pheno MSSM".
  void fixedPointIteration(void (*boundaryCondition)
			   (Softsusy<SoftPars>&, const DoubleVector &),
			   double mxGuess, 
			   const DoubleVector & pars, int sgnMu, double tanb,
			   const QedQcd_legacy & oneset, bool gaugeUnification, 
			   bool ewsbBCscale =  false); 
  /// legacy wrapper to provide backward compatibility: does the same as the
  /// above 
  void lowOrg(void (*boundaryCondition)
		(Softsusy<SoftPars>&, const DoubleVector &),
		double mxGuess, 
		const DoubleVector & pars, int sgnMu, double tanb,
		const QedQcd_legacy & oneset, bool gaugeUnification, 
		bool ewsbBCscale =  false) {
    fixedPointIteration(boundaryCondition, mxGuess, pars, sgnMu, tanb, oneset,
			gaugeUnification, ewsbBCscale);
  }; 

  /// Main iteration routine: 
  /// Boundary condition is the theoretical condition on parameters at the high
  /// energy scale mx: the parameters themselves are contained within the
  /// vector. IO parameters:  
  /// maxTries is the maximum number of iterations allowed, mx is the GUT
  /// scale (negative if you require gauge unification),
  /// sgnMu is the desired sign of mu: + or - 1. If mu is 0, mu is set
  /// initially as a boundary condition. tanb = desired value of DR bar tan
  /// beta(MZ).

  void itLowsoft(int maxTries, int sgnMu, double tol, 
		 double tanb, void (*boundaryCondition)(Softsusy &, 
							const DoubleVector &), 
		 const DoubleVector & pars, bool gaugeUnification, 
		 bool ewsbBCscale);
  
  /// Dummy function to allow users to re-define it in user supplied
  /// objects. Pars will contain necessary parameters to describe high-scale
  /// boundary conditions on SUSY breaking terms
  virtual void rpvSet(const DoubleVector & parameters);

  /// Dummy function to allow users to re-define it in user supplied
  /// objects. Pars will contain necessary parameters to describe high-scale
  /// boundary conditions on SUSY breaking terms
  virtual void methodBoundaryCondition(const DoubleVector & pars);

  /// Works out how best to fit the isajet numbers to the spectrum.
  /// There are problems with the Higgs and sbottoms because ISAJET assumes
  /// certain tree-level relations between masses that are broken by SOFTSUSY's
  /// higher accuracy. The differences get large for high tan beta around 50,
  /// at around 10 they're typically only a percent.
  /// Output parameters: mtopPole is the pole mass of the top quark, mGPole is
  /// the gluino pole mass, smu is the superpotential Higgs parameter, mA is
  /// the pseudoscalar Higgs mass, tanb is tan beta, mq1l is the mass of the
  /// LH first family squark, mdr is the mass of the RH first family squark,
  /// meL is the LH selectron mass, meR is the RH selectron mass, mql3 is the
  /// LH stop mass, mdr3 is the RH sbottom mass, mur3 is the RH stop mass,
  /// mtauL is the LH stau mass, mtauR is the RH stau mass, at is the stop
  /// trilinear term, ab is the sbottom trilinear term, atau is the stau
  /// trilinear term, mq2l, msr, mcr are 2nd family squark masses, mmuL, mmuR
  /// are the smu masses, m1 and m2 are the gaugino mass parameters.
  void isajetNumbers764
  (double & mtopPole, double & mGPole, double & smu, double & mA, 
   double & tanb, double & mq1l, double & mdr, double & mur, double & meL, 
   double & meR, double & mql3, double & mdr3, double & mur3, double &  mtauL, 
   double & mtauR, double & at, double & ab, double & atau, double & mq2l, 
   double & msr, double & mcr, double & mmuL, double & mmuR, double & m1, 
   double & m2) const; 
  /// prints a file into fname which acts as an input to isajet
  void isajetInterface764(const char fname[80]) const;
  /// prints a file into fname which acts as an input to isassrun: fstream
  /// should be opened before calling
  void ssrunInterface764Inside(const char fname [80], std::fstream & ) const;
  /// prints a file into fname which acts as an input to isajet
  /// First name input is the name of an OUTPUT file from ssrun, the second
  /// name is the name of the interface file for INPUT to ssrun
  void ssrunInterface764(const char fname [80], const char softfname [80]) const;
  /// Prints a file into fnamesoft that can be input into isawig. fnamein
  /// gives the ISASSRUN output file name, fnameout is the filename for the
  /// ISAWIG input file
  void isawigInterface764(const char fnamein [80], const char fnameout [80], 
			  const char fnamesoft[80]) const;
  /// Outputs with Les Houches accord conventions to standard output. 
  /// Inputs:
  /// out should be something like cout or fout depending on whether you want
  /// output in a file or not.
  /// model contains what form of model is used for the SUSY breaking terms
  /// (eg sugra, gmsb, amsb, nonUniversal). qMax is only relevant if you want
  /// a gridded output of running parameters up to some scale qMax. Put
  /// numPoints = 1 if you don't want to use this option - then qMaz is
  /// immaterial. mb is mb(mb) in the MSbar scheme used to produce the output,
  /// whereas mtau is the pole mass used (eg 1.777). mgut is the GUT scale
  /// that has been determined, and altEwsb is true if you specified mu and mA
  /// as input parameters (not tan beta and mH1, mH2).
  virtual void lesHouchesAccordOutput(std::ostream & out, const char model[], 
				      const DoubleVector & pars, 
				      int sgnMu, double tanb, double qMax, 
				      int numPoints, 
				      bool ewsbBCscale);
  void slha1(std::ostream & out, const char model[], const DoubleVector & pars, 
	     int sgnMu, double tanb, double qMax, int numPoints, 
	     bool ewsbBCscale);
  /// Normally, this is just a dummy function that is un-used. But sometimes,
  /// it can be used to re-set electroweak symmetry breaking conditions,
  /// depending on the value of inputs at the GUT scale. 
  virtual void setEwsbConditions(const DoubleVector & inputs);
  /// This does the job of the above method, but performs the header and
  /// CHANGES TO 8 SIG FIGS
  void headerSLHA(std::ostream & out);
  /// This does the SPINFO block of SLHA
  void spinfoSLHA(std::ostream & out);
  /// This does the MODSEL block of SLHA
  void modselSLHA(std::ostream & out, const char model[]);
  /// SMINPUTS block of SLHA
  void sminputsSLHA(std::ostream & out);
  /// MINPAR block of SLHA
  void minparSLHA(std::ostream & out, const char model [], 
		  const DoubleVector & pars, double tanb, int sgnMu, 
		  bool ewsbBCscale);
  /// EXTPAR block of SLHA
  virtual void extparSLHA(std::ostream & out, const DoubleVector & pars, 
			  bool ewsbBCscale);
  /// This does the job of the above method, but outputs the Mass block
  void massSLHA(std::ostream & out);
  /// higgs part of mass block of SLHA
  virtual void higgsMSLHA(std::ostream & out);
  /// neutralino and charigno part of mass block of SLHA
  virtual void neutralinoCharginoMSLHA(std::ostream & out);
  /// sfermions part of mass block of SLHA
  virtual void sfermionsSLHA(std::ostream & out);
  /// hmix block output of SLHA
  void sfermionmixSLHA(std::ostream & out);
  /// nmix block output of SLHA
  virtual void neutralinoMixingSLHA(std::ostream & out);
  /// This does the job of the above method, but outputs the UMIX/VMIX blocks
  void inomixingSLHA(std::ostream & out);
  /// SOFTSUSY comments in SLHA
  void softsusySLHA(std::ostream & out);
  /// sfermionic part of mixing blocks
  void alphaSLHA(std::ostream & out);
  /// higgs part of mixing blocks
  virtual void hmixSLHA(std::ostream & out);
  /// sfermionic part of mixing blocks
  void gaugeSLHA(std::ostream & out);
  /// sfermionic part of mixing blocks
  virtual void yukawasSLHA(std::ostream & out);
  /// MSOFT block of SLHA
  virtual void msoftSLHA(std::ostream & out);
  /// outputs DRbar parameters at scale Q (default is at MSUSY)
  virtual void drbarSLHA(std::ostream & out, int numPoints, double qMax, int n);

  /// Sets the minimum of potential to be the difference between the UFB-3
  /// direction minimum and the standard EW breaking minimum. mgut is
  /// obviously the high-scale at which boundary conditions are employed
  void doUfb3(double mgut);

  /// Utility function: sets Higgs masses of neutral Higgs' (higgsm) and
  /// charged (higgsc). They should be of dimension 4 and 2 respectively.
  /// Also sets couplings dnu(4), dnd(4) and cn(4). beta is from tan beta.
  void assignHiggs(DoubleVector & higgsm, DoubleVector & higgsc, 
		   DoubleVector & dnu, DoubleVector & dnd, 
		   DoubleVector & cn, double beta) const;

  /// Utility function: sets Higgs masses of neutral Higgs' (higgsm) and
  /// charged (higgsc). They should be of dimension 4 and 2 respectively.
  void assignHiggs(DoubleVector & higgsm, DoubleVector & higgsc) const;

  /// Utility function: sets Higgs masses of neutral Higgs' (higgsm) and
  /// charged (higgsc) constributions for sfermions. They should be of
  /// dimension 4 and 2 respectively. 
  /// Also sets couplings dnu(4), dnd(4) and cn(4). beta is from tan beta.
  void assignHiggsSfermions(DoubleVector & higgsm, DoubleVector & higgsc, 
			    DoubleVector & dnu, DoubleVector & dnd, 
			    DoubleVector & cn, double beta) const;

  /// predicts MW in Standard Model limit: parameterisation from
  /// hep-ph/0311148: estimated 4 MeV on the error
  double smPredictionMW() const;

  /// gives two-loop contribution to (g-2)_mu according to the calculation in 
  /// Stockinger et al
  double twoLoopGm2(double amu1Loop) const;

  /// Two-loop O(\alpha_s^2) MSSM corrections to the pole masses of heavy quarks
  /// by A.Bednyakov, A.Onishchenko, V.Velizhanin, O.Veretin  
  double twoLpMt() const;
  double twoLpMb() const;

  /// Input diagonal mass matrices and it'll give you back mixed ones, based on
  /// the CKM quark mixing matrix you supplied in vCkm. 
  /// The idea is that Softsusy objects are UNmixed. Therefore this method
  /// does nothing. Derived objects may have mixing implemented
  virtual void doQuarkMixing(DoubleMatrix & mDon, DoubleMatrix & mUpq);

  /// Provides the first guess at a SUSY object at mt, inputting tanb and
  /// oneset (should be at MZ) - it's very crude, doesn't take radiative
  /// corrections into account etc. oneset provides low energy data and
  /// tanb=tan beta 
  virtual MssmSusy guessAtSusyMt(double tanb, const QedQcd_legacy & oneset);
};

typedef Softsusy<SoftParsMssm> MssmSoftsusy;

template<class SoftPars>
Softsusy<SoftPars>::Softsusy()
  : SoftPars(), AltEwsbMssm(), physpars(), forLoops(), 
    problem(), msusy(0.0), minV(6.66e66), 
    mw(0.0), dataSet(), fracDiff(1.), setTbAtMX(false), altEwsb(false), 
    predMzSq(0.), t1OV1Ms(0.), t2OV2Ms(0.), t1OV1Ms1loop(0.), 
    t2OV2Ms1loop(0.), alternativeMs(false), mxBC(mxDefault) { 
      setPars(110);
      setMu(0.0);
      setLoops(0);
      setThresholds(0);
}


template<class SoftPars>
Softsusy<SoftPars>::Softsusy(const Softsusy & s)
  : SoftPars(s.displaySoftPars()), AltEwsbMssm(s.displayAltEwsbMssm()), 
    physpars(s.displayPhys()), 
    forLoops(s.displayDrBarPars()), 
    problem(s.problem), msusy(s.msusy), minV(s.minV), 
    mw(s.mw), dataSet(s.displayDataSet()), fracDiff(s.displayFracDiff()), 
    setTbAtMX(s.displaySetTbAtMX()), 
    altEwsb(s.displayAltEwsb()), predMzSq(s.displayPredMzSq()), 
    t1OV1Ms(s.displayTadpole1Ms()), t2OV2Ms(s.displayTadpole2Ms()), 
    t1OV1Ms1loop(s.displayTadpole1Ms1loop()), 
    t2OV2Ms1loop(s.displayTadpole2Ms1loop()), alternativeMs(s.alternativeMs),
    mxBC(s.displayMxBC()) {

    setPars(110);
    setMu(s.displayMu()); 
    setLoops(s.displayLoops());
    setThresholds(s.displayThresholds());
}

template<class SoftPars>
Softsusy<SoftPars>::Softsusy(const Susy &s)
  : SoftPars(s), AltEwsbMssm(), 
    physpars(), forLoops(), problem(), 
    msusy(0.0), minV(6.66e66), mw(0.0), dataSet(), fracDiff(1.), 
    setTbAtMX(false), altEwsb(false), predMzSq(0.), t1OV1Ms(0.), 
    t2OV2Ms(0.), t1OV1Ms1loop(0.), t2OV2Ms1loop(0.), alternativeMs(false),
    mxBC(mxDefault) { 
  setPars(110);
  setMu(s.displayMu()); 
  setLoops(s.displayLoops());
  setThresholds(s.displayThresholds());
}


template<class SoftPars>
Softsusy<SoftPars>::Softsusy
(const SoftPars & s, const sPhysical & sp, double mu, int l, int t, 
 double hv) 
  : SoftPars(s), AltEwsbMssm(), physpars(sp), forLoops(), problem(), msusy(0.0),
    minV(6.66e66), mw(0.0), dataSet(), fracDiff(1.), setTbAtMX(false), 
    altEwsb(false), predMzSq(0.), t1OV1Ms(0.), 
    t2OV2Ms(0.), t1OV1Ms1loop(0.), t2OV2Ms1loop(0.), alternativeMs(false),
    mxBC(mxDefault) {
  setHvev(hv);
  setPars(110);
  setMu(mu);
  setLoops(l);
  setThresholds(t);
}



template<class SoftPars>
const Softsusy<SoftPars> & Softsusy<SoftPars>::displayMssmSoft() const { return *this; }

template<class SoftPars>
const QedQcd_legacy & Softsusy<SoftPars>::displayDataSet() const { return dataSet; }

template<class SoftPars>
const sPhysical & Softsusy<SoftPars>::displayPhys() const { return physpars; }

template<class SoftPars>
const drBarPars & Softsusy<SoftPars>::displayDrBarPars() const { 
  return forLoops; 
}

template<class SoftPars>
double Softsusy<SoftPars>::displayMinpot() const { return minV; } 
template<class SoftPars>
double Softsusy<SoftPars>::displayMsusy() const { return msusy; } 
template<class SoftPars>
double Softsusy<SoftPars>::displayMw() const { return mw; } 

template<class SoftPars>
double Softsusy<SoftPars>::displayTadpole1Ms() const {
  return t1OV1Ms; 
}

template<class SoftPars>
double Softsusy<SoftPars>::displayTadpole2Ms() const {
  return t2OV2Ms; 
}

template<class SoftPars>
double Softsusy<SoftPars>::displayTadpole1Ms1loop() const {
  return t1OV1Ms1loop; 
}

template<class SoftPars>
double Softsusy<SoftPars>::displayTadpole2Ms1loop() const {
  return t2OV2Ms1loop; 
}

template<class SoftPars>
void Softsusy<SoftPars>::setMinpot(double f) { minV = f; }
template<class SoftPars>
void Softsusy<SoftPars>::setMsusy(double f) { msusy = f; }
template<class SoftPars>
void Softsusy<SoftPars>::setMw(double f) { mw = f; }
template<class SoftPars>
void Softsusy<SoftPars>::doUfb3(double mgut) { setMinpot(ufb3sl(mgut) -
							  realMinMs()); } 
/// Prints out header line for print-short output
void printShortInitialise();
/// returns the square root of the absolute value of the argument
// returns sqrt(f) for f>0 
inline double ccbSqrt(double f){ return sqrt(fabs(f)); }
/// returns f * f * sign(f)
inline double signedSqr(double f){ if (f > 0.) return sqr(f); 
  else return -sqr(f); }
/// Two-loop Standard Model corrections to rho parameter
double rho2(double r);

/// For a given trial value of the log of field H2, gives the value of the
/// potential at the minimum. The following global variables must be set before
/// it is called: unificationScale=high BC scale, minTol=fractional accuracy
/// with which minimum is found
double minimufb3(double);

/// function used for calculating sin theta_eff
double gEff(double x);
/// function used for calculating sin theta_eff
double fEff(double x);

/// parameterisation of LEP2 likelihood for Standard Model higgs mass
double lep2Likelihood(double mh);
/// integrand for 2 GeV Gaussian theory error for SM higgs likelihood from LEP2
DoubleVector mhIntegrand(double mh, const DoubleVector & y);
/// integrates LEP2 likelihood for SM higgs with 2 GeV Gaussian error, returns
/// the log likelihood
double lnLHiggs(double mh);

//double averageMus(susyMu, muOld);

#include "softsusy.cpp"

} // namespace softsusy

#endif
