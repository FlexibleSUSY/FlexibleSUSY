
/** \file def.cpp
   - Project:     SOFTSUSY
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/

*/

#include "def.h"

/// global variable declaration
/// no quark mixing (dominant third family approx), and no verbose output
namespace softsusy {
  int MIXING = 1, PRINTOUT = 0;
  /// fractional accuracy required - may need to make smaller for accurate
  /// neutrino masses
  double TOLERANCE = 1.0e-4;
  /// decay constant of muon
  double GMU = flexiblesusy::Electroweak_constants::gfermi;
  /// there are two possible conventions: if QEWSB > MZ, its value is assumed
  /// in GeV and used as a constant MSUSY. Otherwise, it MULTIPLIES the usual 
  /// MSUSY value, of root(mstop1 mstop2)
  double QEWSB = 1.0; 
  /// Do we include 2-loop RGEs of *all* scalar masses and A-terms, or only the
  /// scalar mass Higgs parameters? (Other quantities all 2-loop anyway): the
  /// default in SOFTSUSY 3 is to include all 2-loop terms
  bool INCLUDE_2_LOOP_SCALAR_CORRECTIONS = true;
  /// PA: switch for selecting between the Z3 preserving case and the
  // Z3 violating case with mu and m3sq allowed
  // needed fpr selecting which EWSB routines to use.
  bool Z3 = true;
  //If true the input value of lambda is set at the GUT scale
  //If false the input lambda is set at MSUSY.
  bool GUTlambda = false;
  //If true the input value of kappa is set at the GUT scale
  //If false the input kappa is set at MSUSY.
  bool GUTkappa = false;
  bool GUTmuPrime = false;
  bool GUTxiF = false;
  bool GUTsVev = false;
  /// Enable/Disable NMSSMTools compatible SLHA output
  bool NMSSMTools = false;
  //If true then the EWSB conditions will output soft Higgs masses
  //Will be inconsistent with constrained models
  //but can be useful for non-universal Higgs cases 
  bool SoftHiggsOut = false;
  int MICROMEGAS = 0;
  int NMSDECAY = 0;
  /// number of loops used to calculate Higgs mass and tadpoles. They should be
  /// identical for a consistent calculation
  int numHiggsMassLoops = 2, numRewsbLoops = 2;
  /// global pole mass of MZ in GeV - MZCENT is defined in def.h
  double MZ = MZCENT;
  /// default is to fix RPV parameters at the GUT scale
  bool susyRpvBCatMSUSY = false;
  /// default: use SLHA1 conventions for RPV output
  bool forceSlha1 = false;
  /// default is to have no trilinears set by SLHA2 conventions
  bool slha2setTrilinear[] = {false, false, false, false, false, false, false, 
			      false, false, false, false, false, false, false, 
			      false, false, false, false, false, false, false, 
			      false, false, false, false, false, false};
  /// default is to *not* print out theoretically excluded spectra
  bool printRuledOutSpectra = false;
  /// default is to do a quick calculation instead of trying hard to get
  /// convergence 
  bool tryToConvergeHard = false;
  /// default is to set tree-level tachyonic A masses to 0 in loops
  bool mAFlag = false;
  /// default is that developer-level debugging printout flag is OFF
  bool printDEBUG = false;
  /// random number seed
  long idummySave = -22;
  /// Default is to not include 2-loop MSSM corrections in degenerate
  /// gluino/squark mass limit 
  bool includeTwoLoopMssmCorrectionsToMt = false;
  /// Default: use SOFTSUSY conventions for masses of sparticles in loops, ie
  /// tree-level masses computed with the 2-loop Higgs potential
  bool sphenoMassConv = false;

  /// Controls the use of MSSM three-loop RGEs
  bool USE_THREE_LOOP_RGE = false;
  
#ifdef COMPILE_FULL_SUSY_THRESHOLD 

  /// Threshold to prevent the re-evaluation of two-loop leading SUSY thresholds
  /// corrections. By default we set it to 10%. 
  double TWOLOOP_NUM_THRESH = 0.1;
  /// Includes the evaluation of leading two-loop thresholds corrections
  /// to the strong coupling constant and to the third family of fermion masses 
  bool USE_TWO_LOOP_THRESHOLD = false;
  /// just implements decoupling procedure "consistently" for
  /// the case of b-quark mass. It requires the external momentum to be zero. 
  /// However, the difference between the p^2 = 0 and p^2 = mb^2 cases
  /// is expected to be of O((mb/MSUSY)^2), which we can formally neglect.
  bool MB_DECOUPLING = false;

#endif //COMPILE_FULL_SUSY_THRESHOLD
} // namespace softsusy
/// end of global variable declaration
