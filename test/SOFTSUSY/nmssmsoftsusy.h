/** \file nmssmsoftsusy.h
    - Project: Next-to-Minimal SOFTSUSY 
    - Author: Ben Allanach, Peter Athron, Lewis Tunstall, Alexander Voigt,
      Anthony Williams
    - Manual: http://arxiv.org/abs/1311.7659
    - Webpage: http://hepforge.cedar.ac.uk/softsusy/

    \brief NmssmSoftsusy object contains all the NMSSM functionality: SUSY
    breaking and preserving parameters, and physical parameters etc
*/


#ifndef NMSSMSOFTSUSY_H
#define NMSSMSOFTSUSY_H

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
#include "nmssmsoftpars.h"
#include "softsusy.h"
#include "mssmUtils.h"
#include "nmssm2loop.h"
#include "nmssmUtils.h"

namespace softsusy {

/// Contains all supersymmetric NMSSM parameters, incorporating R_p NMSSM
class NmssmSoftsusy: public Softsusy<SoftParsNmssm> {
private:
  
   double tSOVSMs;  ///< New Nmssm DRbar tadpole(MSusy): incl 2 loops
   double tSOVSMs1loop; ///<New Nmssm DRbar tadpole(MSusy): excl 2 loops

public:
//  void (*boundaryCondition)(NmssmSoftsusy &, const DoubleVector &);
  /// Default constructor fills object with zeroes
  NmssmSoftsusy();
  /// Constructor sets SUSY parameters only from another object
  NmssmSoftsusy(const NmssmSusy &);
  /// Constructor copies another object
  NmssmSoftsusy(const NmssmSoftsusy &);
  /// Constructor copies an MSSM object, needed for test vs mssm
  NmssmSoftsusy(const MssmSoftsusy &);
  /// Sets all parameters from s, sp, mu is the mu superpotential parameter, l
  /// is the number of loops used for RG evolution, t is the thresholds
  /// accuracy parameter, mg is the gravitino mass, hv is the Higgs VEV
  /// parameter.
  NmssmSoftsusy(const SoftParsNmssm & s, const sPhysical & sp, double mu, int l, int t, double hv);
  /// Set all data in the object equal to another
  const NmssmSoftsusy & operator=(const NmssmSoftsusy & s);
  
  double displayTadpoleSMs() const; ///< displays t_s/s tadpole
  double displayTadpoleSMs1loop() const; ///< displays t_2/v_s tadpole@1 loop
  void setTadpoleSMs(double t) { tSOVSMs = t; }
  /// PA: obtains NMSSM P1-sfermion-sfermion couplings 
  //for 3rd generation sfermions
  void  P1SfSfCouplings(DoubleMatrix & lTP1Lr, DoubleMatrix & lBP1Lr, 
			DoubleMatrix  & lTauP1Lr) const;
  /// PA: obtains NMSSM P2-sfermion-sfermion couplings 
  //for 3rd generation sfermions
  void  P2SfSfCouplings(DoubleMatrix & lTP2Lr, DoubleMatrix & lBP2Lr, 
			DoubleMatrix  & lTauP2Lr) const;
  /// PA: obtains NMSSM P3-sfermion-sfermion couplings 
  //for 3rd generation sfermions
  void  P3SfSfCouplings(DoubleMatrix & lTP3Lr, DoubleMatrix & lBP3Lr, 
			DoubleMatrix  & lTauP3Lr) const;

  /// PA: obtains NMSSM H1-sfermion-sfermion couplings 
  //for 3rd generation sfermions
  void  H1SfSfCouplings(DoubleMatrix & lTS1Lr, DoubleMatrix & lBS1Lr, 
			DoubleMatrix  & lTauS1Lr, double gmzOcthW, 
			double mu,  double cosb, double v1) const;
  /// PA: obtains NMSSM H2-sfermion-sfermion couplings 
  //for 3rd generation sfermions
  void H2SfSfCouplings(DoubleMatrix & lTS1Lr, DoubleMatrix & lBS1Lr, 
		       DoubleMatrix  & lTauS1Lr, double gmzOcthW, 
		       double mu,  double sinb) const;
  /// PA: obtains NMSSM S-sfermion-sfermion couplings 
  //for 3rd generation sfermion
  void SSfSfCouplings(DoubleMatrix & lTS3Lr, DoubleMatrix & lBS3Lr, 
		      DoubleMatrix & lTauS3Lr, double lam) const;
  /// PA:routine to calculate sfermiom contributions to (16 \pi^2) ts / s 
  double doCalcTadSSfermions(DoubleMatrix lTS3Lr, DoubleMatrix lBS3Lr, 
			     DoubleMatrix lTauS3Lr, double q, double s) const;
  /// PA: for loop corrections, helps adding Higgs corrections in a tidy way
  void assignHiggs(DoubleVector & higgsm, DoubleVector & higgsa, 
                     DoubleVector & higgsc) const;
  /// PA: NMSSM routine to obtain Higgs loop parts of (16 \pi^2) t1/v1
//Includes goldstone bosons. 
  double doCalcTad1Higgs(double q, double costhDRbar, double g, double tanb) const;
  /// PA: NMSSM routine to obtain Higgs loop parts of (16 \pi^2) t2/v2
  //Includes goldstone bosons.
  double doCalcTad2Higgs(double q, double costhDRbar, double g, double tanb) const;
  /// PA: NMSSM routine to obtain Higgs loop parts of (16 \pi^2) ts/s
  //Includes goldstone bosons.
  double doCalcTadSHiggs(double q, double tanb) const;
  /// PA: NMSSM routine to obtain Neutralino loop parts of (16 \pi^2) t1/v1
  double doCalcTad1Neutralinos(double q, double costhDRbar, double g, 
                               double tanb) const;
  /// PA: NMSSM routine to obtain Neutralino loop parts of (16 \pi^2) t2/v2
  double doCalcTad2Neutralinos(double q, double costhDRbar, 
			     double g, double sinb) const;
  /// PA: NMSSM routine to obtain Neutralino loop parts of (16 \pi^2) tS/s
  double doCalcTadSNeutralinos(double q,  double lam, double kap) const;
  /// PA: NMSSM routine to obtain Chargino loop parts of (16 \pi^2) tS/s
  double doCalcTadSCharginos(double q,  double lam) const;
/// Does the calculation of one-loop pieces of \f$ t_1 / v_1 \f$ 
  virtual double doCalcTadpole1oneLoop(double mt, double sinthDRbar) const;
  /// Does the calculation of one-loop pieces of \f$ t_2 / v_2 \f$ 
  virtual double doCalcTadpole2oneLoop(double mt, double sinthDRbar) const;
  /// Does the calculation of one-loop pieces of \f$ t_s / s \f$ 
  virtual double doCalcTadpoleSoneLoop(double mt, double sinthDRbar) const;

/// Calculates and sets the one-loop pieces of \f$ t_1 / v_1 \f$: sets both
  /// 1-loop and total pieces equal to the one-loop piece
  virtual void calcTadpole1Ms1loop(double mt, double sinthDRbar);
  /// Calculates then sets the one-loop pieces of \f$ t_2 / v_2 \f$: sets both
  /// 1-loop and total pieces equal to the one-loop piece
  virtual void calcTadpole2Ms1loop(double mt, double sinthDRbar);
  /// Calculates then sets the one-loop pieces of \f$ t_s / s \f$: sets both
  /// 1-loop and total pieces equal to the one-loop piece
  /// Calculates then sets the one-loop pieces of \f$ t_s / s \f$: sets both
  /// 1-loop and total pieces equal to the one-loop piece
  virtual void calcTadpoleSMs1loop(double mt, double sinthDRbar);
  /// PA: calls routines to calculate all three tadpoles and sets them.
  // Currently only works at one loop.  
  // Two loop should be added later. 
  void doTadpoles(double mt, double sinthDRbar);
  /// Organises tree-level calculation of all sparticle masses and mixings
  virtual void calcDrBarPars();
 
  /// LCT: Adds one-loop higgs corrections to stop mass matrix
  /// IO parameters: p=external momentum, mt=DR bar top mass, 
  /// outputs higgs=2x2 matrix of self-energies in LR basis
  DoubleMatrix addStopHiggs(double p, double mt, DoubleMatrix & higgs);
  /// Adds one-loop corrections to stop mass matrix
  /// IO parameters: p=external momentum, mass=tree level mass matrix on
  /// input, is returned with radiative corrections added, mt=DR bar top mass
  virtual void addStopCorrection(double p, DoubleMatrix & mass, double mt);
  /// LCT: Adds one-loop higgs corrections to sup mass matrix
  /// IO parameters: {p1,p2}=external momentum, 
  /// outputs higgs=2x2 matrix of self-energies in LR basis
  DoubleMatrix addSupHiggs(double p1, double p2, int family, 
  DoubleMatrix & higgs);
  /// Adds one-loop corrections to sup mass matrix
  /// IO parameters: mass=tree level mass matrix on
  /// input, is returned with radiative corrections added for a given family
  virtual void addSupCorrection(DoubleMatrix & mass, int family);
  /// LCT: Adds one-loop Higgs corrections to sbottom mass matrix
  /// IO parameters: p=external momentum, mt=DR bar top mass, 
  /// outputs higgs=2x2 matrix of self-energies in LR basis
  DoubleMatrix addSbotHiggs(double p, double mt, DoubleMatrix & higgs);
  /// Adds one-loop corrections to sbottom mass matrix
  /// IO parameters: p=external momentum, mass=tree level mass matrix on
  /// input, is returned with radiative corrections added, mt=DR bar top mass
  virtual void addSbotCorrection(double p, DoubleMatrix & mass, double mt);
  /// LCT: Adds one-loop higgs corrections to sdown mass matrix
  /// IO parameters: {p1,p2}=external momentum, 
  /// outputs higgs=2x2 matrix of self-energies in LR basis
  DoubleMatrix addSdownHiggs(double p1, double p2, int family, 
  DoubleMatrix & higgs);
  /// Adds one-loop corrections to sbottom mass matrix at p=root(mb1 mb2)
  /// IO parameters: mass=tree level mass matrix on
  /// input, is returned with radiative corrections added, mt=DR bar top mass
  virtual void addSdownCorrection(DoubleMatrix & mass, int family);
  /// LCT: Adds one-loop Higgs corrections to tau sneutrino self energy
  /// IO parameters: p=external momentum, outputs higgs=self-energy
  double addSnuTauHiggs(double p, double & higgs);
  /// Adds one-loop corrections to tau sneutrino mass 
  /// IO parameters: p=external momentum, mass=tree level mass matrix on
  /// input, is returned with radiative corrections added
  virtual void addSnuTauCorrection(double & mass);
  /// LCT: Adds one-loop Higgs corrections to sneutrino self energy
  /// IO parameters: p=external momentum, family=generation, 
  /// outputs higgs=self-energy
  double addSnuHiggs(double p, int family, double & higgs);
  /// Adds one-loop corrections to sneutrino mass of family "family"
  /// IO parameters: mass=tree level mass matrix on
  /// input, is returned with radiative corrections added
  virtual void addSnuCorrection(double & mass, int family);
  /// LCT: Adds one-loop electroweak corrections to stau self energy
  /// IO parameters: p=external momentum, mtau=DR bar tau mass,
  /// outputs electroweak=2x2 matrix of self-energies
  DoubleMatrix addStauHiggs(double p, double mtau, DoubleMatrix & higgs);
  /// Adds one-loop corrections to stau mass matrix at p=root(mtau1 mtau2)
  /// IO parameters: mass=tree level mass on
  /// input, is returned with radiative corrections added, mtau=DR bar tau mass
  virtual void addStauCorrection(double p, DoubleMatrix & mass, double mtau);
  /// LCT: Adds one-loop Higgs corrections to slepton self energy
  /// IO parameters: p=external momentum, family=generation, 
  /// outputs higgs=2x2 matrix of self-energies
  DoubleMatrix addSlepHiggs(double p1, double p2, int family, DoubleMatrix & higgs);
  /// Adds one-loop corrections to sel_fam mass matrix at p=root(msel1 msel2)
  /// IO parameters: mass=tree level mass on
  /// input, is returned with radiative corrections added, mt=DR bar top mass
  virtual void addSlepCorrection(DoubleMatrix & mass, int family);
  void getNeutPassarinoVeltman(double p, double q, DoubleMatrix & b0fn, 
  DoubleMatrix & b1fn); 
  /// LCT: Adds Higgs loop corrections to neutralino mass matrix.
  /// IO parameters: p=external momentum, {sigmaL,sigmaR,sigmaS}= 4x4 matrices 
  /// of {left,right,scalar} corrections
  void addNeutLoopHiggs(double p, DoubleMatrix & sigmaL, DoubleMatrix & sigmaR, DoubleMatrix & sigmaS);
  /// Adds the loop corrections on to an input tree-level neutralino mass
  virtual void addNeutralinoLoop(double p, DoubleMatrix &);
  
  /// Returns tree-level up squark mass matrix in "mass" for NMSSM.
  /// IO parameters: mass=tree level mass matrix on
  /// input, is returned with radiative corrections added, mtrun=DR bar top
  /// mass, family=generation of squark, pizztMS=Z self energy at Q=M_SUSY,
  /// sinthDRbarMS=DRbar value of sin theta_w 
virtual void treeUpSquark(DoubleMatrix & mass, double mtrun, 
                    double pizztMS, double sinthDRbarMS, 
                    int family);

/// Returns tree-level down squark mass matrix in "mass" for NMSSM.
  /// IO parameters: mass=tree level mass matrix on
  /// input, is returned with radiative corrections added, mbrun=DR bar bottom
  /// mass, family=generation of squark, pizztMS=Z self energy at Q=M_SUSY,
  /// sinthDRbarMS=DRbar value of sin theta_w 
virtual void treeDownSquark(DoubleMatrix & mass, double mbrun, double pizztMS, 
		double sinthDRbarMS, int family);
  /// Returns tree-level down squark mass matrix in "mass" for NMSSM.
  /// IO parameters: mass=tree level mass matrix on
  /// input, is returned with radiative corrections added, mTrun=DR bar tau
  /// mass, family=generation of slepton, pizztMS=Z self energy at Q=M_SUSY,
  /// sinthDRbarMS=DRbar value of sin theta_w
virtual  void treeChargedSlepton(DoubleMatrix & mass, double mTrun, double pizztMS, 
		double sinthDRbarMS, int family);
  /// LCT: new routine to allocate NMSSM chargino masses, 
  //Returns tree-level chargino mass matrix in the NMSSM 
  void treeCharginos(DoubleMatrix & mass, double beta, double mw);
  
/// LCT: new routine for NMSSM neutralino masses, 
  //Returns tree-level Neutralino mass matrix in the NMSSM 
  void treeNeutralinos(DoubleMatrix & mass, double beta, double mz, double mw, 
                       double sinthDRbar);
  /// PA:  calls  treeCharginos and treeNeutralinos, 
  // performs diagonalisation and fills eg with appropriate values.
  void calcDrBarGauginos(double beta, double mw, double mz, double sinth, 
                         drBarPars & eg);
   /// PA: fills tree level CP even and CP odd Higgs mass matrices 
  //and tree level mHPm. Does *not* use EWSB substitution   
  void treeHiggsAlt(DoubleMatrix & mS, DoubleMatrix & mPpr, DoubleMatrix & mP2, 
                  double & mHpm, double beta) const;
  /// PA: fills tree level CP even and CP odd Higgs mass matrices 
  //and tree level mHPm. Uses EWSB substitution. 
  //Called in higgs and calcDrBarParsHiggs   
  void treeHiggs(DoubleMatrix & mS, DoubleMatrix & mPpr, DoubleMatrix & mP2, 
                  double & mHpm, double beta) const;
  //calculates DrBar Higgs masses and sets them    
  void calcDrBarHiggs(double beta, double mz2, double mw2, double sinthDRbar, 
                      drBarPars & eg);
  
  /// Returns mu from rewsb requirement. 
  /// returns 1 if there's a problem. Call at MSusy
  /// PA: To be used in general Z3 violating nmssm 
  virtual int rewsbMu(int sgnMu, double & mu) const;
  // PA: NMssm rewsb routine which fixes imn much the same way as 
  // mu is fixed in the Mssm using mueff = lambda * s / root 
  // For use in Z3 constrained version or when other scenarios 
  // where mu = 0
  virtual int rewsbSvev(int sgnMu, double & svev) const;
  /// returns 1 if mu < 1.0e-9
  /// PA:  nmssm version for use in Z3 violating case.  
  virtual int rewsbM3sq(double mu, double & m3sq) const;
  /// PA:: In case of Z3 invariance EWSB outputs kappa instead.
  virtual int rewsbKap(double & kap) const;
  /// PA: third EWSB condition (for the singlet Higgs field) 
  //new with respect to the MSSM.
  virtual int rewsbXiS(double & xiS) const;
  // PA:For Z3 invariant NMSSM when we solve for s, kappa and mS
  // or for non universal Higgs
  virtual int rewsbmSsq(double & mSsq) const;
  /// PA: for non universal Higgs
  virtual int rewsbmH1sq(double & mH1sq) const;
  /// PA: for non universal Higgs
  virtual int rewsbmH2sq(double & mH211sq) const;
  /// PA: Imposes EWSB at the tree level. 
  virtual void rewsbTreeLevel(int sgnMu);
  /// PA: finds mu iteratively in the casew where we use EWSB to swap 
  //(mu, m3sq, XiS) --> (mZ, tb, s) 
  // Uses the full one loop tadpole from Degrassi and Slavich.  
  // No two loop added yet.
  void iterateMu(double & munew, int sgnMu, double mt, 
		 int maxTries, double pizztMS, double sinthDRbar, double tol,
		 int  & err);
  //Routine for iteratively solving for the singlet vev, s = <S>.
  // where the EWSB is used to swap (kappa, mS) --> (mZ, tb)   
  // and determine s.
  // Uses the full one loop tadpole from Degrassi and Slavich.  
  // No two loop added yet.
  void iterateSvev(double & sold, int sgnMu,
		   double mt, int maxTries, double pizzMS,
		   double sinthDRbar, double tol, int & err);
  /// PA: organises imposition of EWSB conditions. 
  // Currently at full one loop, two loop to be added.
  // Currently works for two cases depending on Z3 switch
  //Z3 = true:  s --> mZ, kappa --> tan beta, mS --> s 
  //ie (kappa, mS) --> (mZ, tb)
  //Z3 = false:  mu --> mZ, m3sq --> tan beta, s --> XiS  (Z3 = false) 
  // ie (mu, m3sq, XiS) --> (mZ, tb, s) 
  void rewsb(int sgnMu, double mt, double muOld = -6.66e66, double eps = 0.);
  /// Calculates physical sparticle masses to accuracy number of loops. Should
  /// be called at M_{SUSY}.
  virtual void physical(int accuracy);
  //// PA: Fills sigmaL, sigmaR and sigmaS with NMSSM Higgs parts 
  ///of the charginos loop corrections.
  virtual void addChaLoopHiggs(double p, DoubleMatrix & sigmaL, DoubleMatrix & sigmaR, DoubleMatrix & sigmaS, DoubleMatrix b1pCha, DoubleMatrix b0pCha, DoubleMatrix b1pNeut, DoubleMatrix b0pNeut) const;
/// Adds the loop corrections on to an input tree-level chargino mass
 virtual void addCharginoLoop(double p, DoubleMatrix &);
  /// Calculates pole chargino masses and mixing.
  // IO parameters: piwwt is the W self-energy at the current,
  /// accuracy is the number of loops required (0 or 1 currently)
  virtual void charginos(int accuracy, double piwwt, sPhysical & phys);
/// Calculates pole neutralino masses and mixing. 
//IO parameters: piwwt is the W self-energy at M_SUSY,
  /// accuracy is the number of loops required (0 or 1 currently), pizzt is
  /// the Z self-energy at M_SUSY
  virtual void neutralinos(int accuracy, double piwwt, double pizzt, 
                           sPhysical & phys);
  /// Calculates pole Higgs masses and mixings: full 1-loop SUSY corrections
  // leading log two loop for general nmssm and full effective potential 
  //corrections at order alpha_s alpha_t for Z_3 case from Degrassi and Slavich.
  /// IO parameters: piwwt is the W self-energy at M_SUSY, accuracy is number
  /// of loops (0 or 1) to use and pizzt is the Z self-energy at M_SUSY
  /// Returns "true" if there's a tachyon problem
  bool higgs(int accuracy, double piwwt, double pizzt, sPhysical & phys);
  /// PA:: fixes The CP odd mixing matrix with the conventions 
  // Degrassi and Slavich arXiv:0907.4682
  void DegrassiSlavicMix(DoubleMatrix & P) const; 
  /// Calculates Higgs contribution to the transverse part of Z self-energy: 
  //for p=external momentum, Q=renormalisation scale
  virtual double piZZTHiggs(double p, double q,	double thetaWDRbar) const;
  /// Calculates neutralino contrib. to the transverse part of Z self-energy: 
  //for p=external momentum, Q=renormalisation scale
  virtual double piZZTNeutralinos(double p, double Q, double thetaWDRbar) const;
  /// Calculates transverse part of Z self-energy: for p=external momentum,
  /// Q=renormalisation scale
  virtual double piZZT(double p, double Q, bool usePoleMt = false) const;
  /// Calculates Higgs contribution to the transverse part of W self-energy: 
  //for p=external momentum, q=renormalisation scale
  virtual double piWWTHiggs(double p, double q, double thetaWDRbar) const;
  /// Calculates gaugino contr. to transverse part of W self-energy: 
  //for p=external momentum, Q=renormalisation scale
  virtual double piWWTgauginos(double p, double Q, double thetaWDRbar) const;
  /// Calculates transverse part of W self-energy: for p=external momentum,
  /// Q=renormalisation scale
  virtual double piWWT(double p, double Q, bool usePoleMt = false) const;
  /// PA: self energy routines for pseudo scalar self energies
  double pip1p1(double p, double q) const;
  double pip1p2(double p, double q) const;
  double pip2p2(double p, double q) const;
  double pip1p3(double p, double q) const;
  double pip2p3(double p, double q) const;
  double pip3p3(double p, double q) const;
  /// PA: Obtains trilnear couplings of P1-Higgs, for use in loop functions
  void getP1HiggsTriCoup(DoubleMatrix & spp1, DoubleMatrix & hphpp1, double cw2DRbar) const;
   /// PA: Obtains trilnear couplings of P2-Higgs for use in loop functions
  void getP2HiggsTriCoup(DoubleMatrix & spp2, DoubleMatrix & hphpp2, double cw2DRbar) const;
   /// PA: Obtains trilnear couplings of P3-Higgs for use in loop functions
  void getP3HiggsTriCoup(DoubleMatrix & spp2, DoubleMatrix & hphpp2) const;
  /// PA: Obtains trilnear couplings of s1-higgs-higgs for use in loop functions
   void getS1HiggsTriCoup(DoubleMatrix & sss1, DoubleMatrix & pps1, DoubleMatrix & hphps1, double thetaWDRbar) const; 
  /// PA: Obtains trilnear couplings of s2-higgs-higgs for use in loop functions
  void getS2HiggsTriCoup(DoubleMatrix & sss2, DoubleMatrix & pps2, DoubleMatrix & hphps2, double thetaWDRbar) const; 
  /// PA: Obtains trilnear couplings of s3-higgs-higgs for use in loop functions
  void getS3HiggsTriCoup(DoubleMatrix & sss3, DoubleMatrix & pps3, DoubleMatrix & hphps3) const; 
  /// PA:Calculates (16 Pi^2) times the Higgs contribution to Higgs self-energy: 
  //for p=external momentum, q=renormalisation scale
  double pis1s1Higgs(double p, double q) const;
  double pis1s2Higgs(double p, double q) const;
  double pis2s2Higgs(double p, double q) const;
  double pis1s3Higgs(double p, double q) const;
  double pis2s3Higgs(double p, double q) const;
  double pis3s3Higgs(double p, double q) const;
  /// PA: obtains CP even Higgs-Neutralino couplings
  void getS1NeutralinoCoup(ComplexMatrix & aChi, ComplexMatrix & bChi) const;
  void getS2NeutralinoCoup(ComplexMatrix & aChi, ComplexMatrix & bChi) const;
  void getS3NeutralinoCoup(ComplexMatrix & aChi, ComplexMatrix & bChi) const;
  /// PA: obtains CP odd Higgs-Neutralino couplings
  void getP1NeutralinoCoup(ComplexMatrix & aChi, ComplexMatrix & bChi) const;
  void getP2NeutralinoCoup(ComplexMatrix & aChi, ComplexMatrix & bChi) const;
  void getP3NeutralinoCoup(ComplexMatrix & aChi, ComplexMatrix & bChi) const;
  /// PA:Calculates (16 Pi^2) times Neutralino contrib. to Higgs self-energy: 
  //for p=external momentum, q=renormalisation scale
  double pis1s1Neutralinos(double p, double q) const;
  double pis1s2Neutralinos(double p, double q) const;
  double pis2s2Neutralinos(double p, double q) const;
  double pis1s3Neutralinos(double p, double q) const;
  double pis2s3Neutralinos(double p, double q) const;
  double pis3s3Neutralinos(double p, double q) const;
  /// PA:Calculates (16 Pi^2) times Chargino contrib. to Higgs self-energy: 
  //for p=external momentum, q=renormalisation scale
  double pis1s3Charginos(double p, double q) const; 
  double pis2s3Charginos(double p, double q) const;
  double pis3s3Charginos(double p, double q) const;
  
   /// PA:Calculates (16 Pi^2) times sfermion contrib. to Higgs self-energy: 
  //for p=external momentum, q=renormalisation scale
  double pis1s3Sfermions(double p, double q, DoubleMatrix ls1tt,  
                         DoubleMatrix ls1bb,  DoubleMatrix ls1tautau, 
                         DoubleMatrix ls3tt,  DoubleMatrix ls3bb,  
                         DoubleMatrix ls3tautau) const;
  double pis2s3Sfermions(double p, double q, DoubleMatrix ls2tt,
                         DoubleMatrix ls2bb,  DoubleMatrix ls2tautau, 
                         DoubleMatrix ls3tt,  DoubleMatrix ls3bb,  
                         DoubleMatrix ls3tautau) const;
  double pis3s3Sfermions(double p, double q, DoubleMatrix ls3tt,  
                         DoubleMatrix ls3bb, DoubleMatrix ls3tautau) const;
     
  /// Calculates the Higgs self-energy: for p=external momentum,
  /// Q=renormalisation scale
  double pis1s1(double p, double q) const;
  double pis1s2(double p, double q) const;
  double pis2s2(double p, double q) const;
  double pis1s3(double p, double q) const;
  double pis2s3(double p, double q) const;
  double pis3s3(double p, double q) const;

  /// PA: gets h1 mixing element with Hu.
  virtual double h1s2Mix();

  /// LCT: Charged Higgs trilinear couplings with CP-odd and CP-even Higgs'
  void getHp1HiggsTriCoup(DoubleMatrix & ahphp1, DoubleMatrix & hhphp1) const;
  void getHp2HiggsTriCoup(DoubleMatrix & ahphp2, DoubleMatrix & hhphp2) const;
  /// LCT: Charged Higgs trilinear couplings with neutralinos/charginos
  void getNeutralinoCharginoHpmCoup(ComplexMatrix & apph1, 
  ComplexMatrix & apph2, ComplexMatrix & bpph1, ComplexMatrix & bpph2) const;

  /// LCT: Returns gauge-Higgs-hpm trilinear couplings with hpm in weak basis
  void getGaugeHiggsHpmCoup(DoubleVector & wmhhp1, DoubleVector & wmahp1, 
  DoubleVector & gamhphp1, DoubleVector & zhphp1, DoubleVector & wmhhp2, 
  DoubleVector & wmahp2, DoubleVector & gamhphp2, DoubleVector & zhphp2) const;

  /// LCT: Calculates (16 Pi^2) times Higgs contribution to (i, j) component of 
  /// the charged Higgs self-energy in the weak basis
  /// for p=external momentum, q=renormalisation scale
  double piHpm11Higgs(double p, double q) const;
  double piHpm12Higgs(double p, double q) const;
  double piHpm22Higgs(double p, double q) const;

  /// LCT: Calculates (16 Pi^2) times gaugino contribution to (i, j) component 
  /// of charged Higgs self-energy in the weak basis
  /// for p=external momentum, q=renormalisation scale
  double piHpm11Gaugino(double p, double q) const;
  double piHpm12Gaugino(double p, double q) const;
  double piHpm22Gaugino(double p, double q) const;

 /// LCT: Calculates the charged Higgs self-energy in mass basis: 
  /// for p=external momentum, q=renormalisation scale
  double piHpHm(double p, double q) const;
  /// PA: routine to calculate the Higgs loop contributions to running Mt
  virtual double calcRunMtHiggs() const;
  /// PA: routine to calculate the Neutralino loop contributions to running Mt
  virtual double calcRunMtNeutralinos() const;
  /// PA: routine to calculate the DrBar running top mass from the pole mass
  // and full one loop self energy, with also two loop qcd included.
  virtual double calcRunningMt();
   /// PA: routine to calculate the Higgs loop contributions to running Mb
  virtual double calcRunMbHiggs() const;
  /// PA: routine to calculate the Neutralino loop contributions to running Mb
  virtual double calcRunMbNeutralinos() const;
  //  PA: routine to calculate the DrBar running bottom mass from the SM MS-bar 
  //bottom m and one loop self energy minus photon and gluon parts.
  virtual double calcRunningMb() const;
  /// PA: routine to calculate the Higgs loop contributions to running Mtau
  virtual double calcRunMtauHiggs() const;
  /// PA: routine to calculate the Neutralino loop contributions to running Mtau
  virtual double calcRunMtauNeutralinos(double mTauSMMZ) const;
  //  PA: routine to calculate the DrBar running tau mass from the SM MS-bar 
  //tau mass and one loop self energy minus photon parts.
  virtual double calcRunningMtau() const;


  NmssmSusy guessAtSusyMt(double tanb, DoubleVector nmpars, const QedQcd_legacy & oneset);
  
  double predTanb(double muSusy = -6.66e66) const;
  double predMzsq(double & tanb, double muOld = -6.66e66, double eps = 0.);

  /// Main iteration routine: 
  /// Boundary condition is the theoretical condition on parameters at the high
  /// energy scale mx: the parameters themselves are contained within the
  /// vector. IO parameters:  
  /// maxTries is the maximum number of iterations allowed, mx is the GUT
  /// scale (negative if you require gauge unification),
  /// sgnMu is the desired sign of mu: + or - 1 when lambda * svev = 0, 
  /// the desired sign of svev when mu = 0 and the sign of in the quadratic 
  /// equation for mu when  both mu and lambda * svev are non-zero. 
  /// tanb = desired value of DR bar tan beta(MZ).
  void itLowsoft(int maxTries, int sgnMu, double tol, 
		 double tanb, void (*boundaryCondition)(NmssmSoftsusy &, 
							const DoubleVector &), 
		 const DoubleVector & pars, const DoubleVector nmpars, 
                 bool gaugeUnification, bool ewsbBCscale);
  

  void lowOrg(void (*boundaryCondition)
		(NmssmSoftsusy &, const DoubleVector &),
		double mxGuess, const DoubleVector & pars, 
                const DoubleVector nmpars, int sgnMu, double tanb,
		const QedQcd_legacy & oneset, bool gaugeUnification, 
		bool ewsbBCscale =  false); 
  
 /// PAq: A print method used in development.  I find it useful and easier to read than couting the normal display function or calling printlong etc.    
  void printall() const;

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
  /// This does the MODSEL block of SLHA
  void modselSLHA(std::ostream & out, const char model[]);
  /// EXTPAR block of SLHA
  virtual void extparSLHA(std::ostream & out, const DoubleVector & pars,
			  bool ewsbBCscale);
  /// higgs part of mass block of SLHA
  virtual void higgsMSLHA(std::ostream & out);
  /// neutralino and charigno part of mass block of SLHA
  virtual void neutralinoCharginoMSLHA(std::ostream & out);
  /// nmnmix block output of SLHA
  virtual void neutralinoMixingSLHA(std::ostream & out);
  /// CP even mixing matrix NMHMIX block of SLHA
  void nmhmixSLHA(std::ostream&);
  /// CP odd mixing matrix NMHMIX block of SLHA
  void nmamixSLHA(std::ostream&);
  /// NMSSMRUN block of SLHA
  void nmssmrunSLHA(std::ostream&);
  /// outputs DRbar parameters at scale Q (default is at MSUSY)
  virtual void drbarSLHA(std::ostream&, int numPoints, double qMax, int n);
};

inline NmssmSoftsusy::NmssmSoftsusy()
                     : Softsusy<SoftParsNmssm>(), tSOVSMs(0.0), tSOVSMs1loop(0.0)  {}


inline NmssmSoftsusy::NmssmSoftsusy(const NmssmSoftsusy & s)
		     : Softsusy<SoftParsNmssm>(s),
                     tSOVSMs(s.tSOVSMs), tSOVSMs1loop(s.tSOVSMs1loop)  {
   
   setPars(121);   
}


inline NmssmSoftsusy::NmssmSoftsusy(const MssmSoftsusy & s)
		     : Softsusy<SoftParsNmssm>(s),
                     tSOVSMs(0), tSOVSMs1loop(0)  {
   setPars(121);   
}


inline NmssmSoftsusy::NmssmSoftsusy(const NmssmSusy &s)
                     : Softsusy<SoftParsNmssm>(s),tSOVSMs(0.0), tSOVSMs1loop(0.0)  {
  setPars(121);
}


inline NmssmSoftsusy::NmssmSoftsusy
(const SoftParsNmssm & s, const sPhysical & sp, double mu, int l, int t,
 double hv): Softsusy<SoftParsNmssm>(s, sp, mu, l, t, hv),tSOVSMs(0.0), tSOVSMs1loop(0.0)  {
 setPars(121);
  
}


inline double NmssmSoftsusy::displayTadpoleSMs() const {
  return tSOVSMs; 
}

inline double NmssmSoftsusy::displayTadpoleSMs1loop() const {
  return tSOVSMs1loop; 
}

} // namespace softsusy

#endif
