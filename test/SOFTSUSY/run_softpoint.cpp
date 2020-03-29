
/** \file softpoint.cpp
    - Project:     SOFTSUSY 
   - Authors:     Ben Allanach, Markus Bernhardt 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 

   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/
   - Description: main calling program: command line interface. Reads Les
   - Houches files and command-line inputs and drives the calculation of a point
   - in parameter space.
*/ 

#include "run_softpoint.h"

extern template class softsusy::Softsusy<softsusy::SoftParsNmssm>;
extern template class softsusy::Softsusy<softsusy::SoftParsMssm>;

using std::istringstream;
using std::ostringstream;
using std::string;
using std::cout;
using std::cerr;
using std::cin;
using std::endl;

using namespace softsusy;

// Returns a string with all characters in upper case: very handy
string ToUpper(const string & s) {
        string result;
        unsigned int index;
        for (index = 0; index < s.length(); index++) {
	  char a = s[index];
	  a = toupper(a);
	  result = result + a;
        }
	
	return result;
    }

void errorCall() {
  ostringstream ii;
  ii << "\n\nSOFTSUSY" << SOFTSUSY_VERSION 
     << " called with incorrect arguments. Need to put either:\n";
  ii << "./softpoint.x leshouches < lesHouchesInput\n for SLHA/SLAH2 input, or\n";
  ii << "./softpoint.x sugra [SUGRA parameters] [other options]\n";
  ii << "./softpoint.x amsb [mAMSB parameters] [other options]\n";
  ii << "./softpoint.x gmsb [mGMSB parameters] [other options]\n";
  ii << "./softpoint.x nmssm sugra [NMSSM flags] [NMSSM parameters] [other options]\n\n";
  ii << "[other options]: --mbmb=<value> --mt=<value> --alpha_s=<value> --QEWSB=<value>\n";
  ii << "--alpha_inverse=<value> --tanBeta=<value> --sgnMu=<value> --tol=<value>\n";
  ii << "--higgsUncertainties gives an estimate of Higgs mass uncertainties\n";
#ifdef COMPILE_FULL_SUSY_THRESHOLD
  if (USE_TWO_LOOP_THRESHOLD) ii << "--two-loop-susy-thresholds switches on leading 2-loop SUSY threshold corrections to third generation Yukawa couplings and g3.\n";
#endif //COMPILE_FULL_SUSY_THRESHOLD
  if (USE_THREE_LOOP_RGE) ii << "--three-loop-rges switches on 3-loop RGEs\n";
  ii << "--mgut=unified sets the scale at which SUSY breaking terms are set to the GUT\n";
  ii << "scale where g1=g2. --mgut=<value> sets it to a fixed scale, ";
  ii << "whereas --mgut=msusy\nsets it to MSUSY\n\n";
  ii << "If you want the R-parity violating MSSM calculation, set any of the following:\n";
  ii << "--lambda <i> <j> <k> <coupling>, the word lambda replaceable with lambdaP\nor lambdaPP for LLE, LQD, UDD coupling, respectively.\n\n";
  ii << "[SUGRA parameters]: --m0=<value> --m12=<value> --a0=<value>\n";
  ii << "[mAMSB parameters]: --m0=<value> --m32=<value>\n";
  ii << "[mGMSB parameters]: --n5=<value> --mMess=<value> --LAMBDA=<value> --cgrav=<value>\n\n";
  ii << "Bracketed entries are numerical values, in units of GeV if they are massive.\n";
  ii << "Warning: entries left unspecified will be assumed to be zero for SUSY breaking\nterms, unified (for mgut) or at their default central values for Standard Model parameters\n";
  ii << "\n"
     "[NMSSM flags]:\n"
     "  --lambdaAtMsusy   input lambda at renormalization scale Q = Msusy\n"
     "\n"
     "[NMSSM parameters]:\n"
     "  --m0= , --m12= , --a0= , --tanBeta= , --mHd2= , --mHu2= ,\n"
     "  --mu= , --m3SqrOverCosBetaSinBeta= , --lambda= , --kappa= ,\n"
     "  --Alambda= , --Akappa= , --lambdaS= , --xiF= , --xiS= ,\n"
     "  --muPrime= , --mPrimeS2= , --mS2=\n"
     "\n"
     "  Unset NMSSM parameters are assumed to be zero\n"
     "\n"
     "NMSSM example:\n"
     "  ./softpoint.x nmssm sugra --m0=125 --m12=200 --tanBeta=10 --a0=-300 \\\n"
     "     --lambda=0.1 --lambdaAtMsusy\n";

  throw ii.str();
}

int main(int argc, char *argv[]) {
  /// Sets up exception handling
  signal(SIGFPE, FPE_ExceptionHandler); 

  double lambdaW = 0., aCkm = 0., rhobar = 0., etabar = 0.;
  NMSSM_input nmssm_input; // NMSSM input parameters

  bool flavourViolation = false, gutScaleOutput = false, 
    higgsUncertainties = false;

  int numPoints = 1;

  double qMax = 0.;

  // Sets format of output: 4 decimal places
  outputCharacteristics(6);

  void (*boundaryCondition)(MssmSoftsusy &, const DoubleVector &)=extendedSugraBcs;
  void (*nmssmBoundaryCondition)(NmssmSoftsusy&, const DoubleVector&) = NmssmMsugraBcs;

  QedQcd_legacy oneset;
  MssmSoftsusy m;
  FlavourMssmSoftsusy k;
  NmssmSoftsusy nmssm;
  MssmSoftsusy * r = &m; 
  // RpvNeutrino kw;
  bool RPVflag = false;
  enum Model_t { MSSM, NMSSM } susy_model = MSSM; // susy model (MODSEL entry 3)
  softsusy::GUTlambda = true;
  softsusy::GUTkappa = true;
  softsusy::GUTmuPrime = true;
  softsusy::GUTxiF = true;
  softsusy::GUTsVev = true;

  k.setInitialData(oneset);

  try {
  if (argc !=1 && strcmp(argv[1],"leshouches") != 0) {
    cerr << "SOFTSUSY" << SOFTSUSY_VERSION << endl;
    if (!strcmp(argv[1], "-v") || !strcmp(argv[1], "--version")) exit(0);
    cerr << "B.C. Allanach, Comput. Phys. Commun. 143 (2002) 305-331,";
    cerr << " hep-ph/0104145\n";
    cerr << "For RPV aspects, B.C. Allanach and M.A. Bernhardt, Comput. "
	 << "Phys. Commun. 181 (2010) 232, arXiv:0903.1805.\n\n";
    cerr << "Low energy data in SOFTSUSY: MIXING=" << MIXING << " TOLERANCE=" 
	 << TOLERANCE << endl;
    cerr << "G_F=" << GMU << " GeV^2" << endl;
  }
  
  double mgutGuess = 2.0e16, tanb = 0.;
  int sgnMu = 1;
  bool gaugeUnification = true, ewsbBCscale = false;
  double desiredMh = 0.;

  // If there are no arguments, give error message,
  // or if none of the options are called, then go to error message
  if (argc == 1 || ( (strcmp(argv[1], "sugra") && 
		      strcmp(argv[1], "amsb") &&
		      strcmp(argv[1], "gmsb") && 
		      strcmp(argv[1], "runto") && 
		      strcmp(argv[1], "leshouches") && 
                      strcmp(argv[1], "nmssm") &&
		      strcmp(argv[1], "-v") &&
		      strcmp(argv[1], "--version"))))
    errorCall();
  
  DoubleVector pars(3); 
  
  const char* modelIdent = "";

  bool compilationProblem = false;
  
  /// Non model specific options
  if (strcmp(argv[1], "leshouches")) {
      for (int i = 2; i < argc; i++) {
	if (starts_with(argv[i],"--mbmb=")) {
	  double mbIn = get_value(argv[i], "--mbmb=");
	  oneset.setMbMb(mbIn);
	}
	else if (starts_with(argv[i],"--higgsUncertainties")) 
	  higgsUncertainties = true;
	else if (starts_with(argv[i],"--tol=")) 
	  TOLERANCE = get_value(argv[i], "--tol=");
	else if (starts_with(argv[i],"--mt=")) 
	  oneset.setPoleMt(get_value(argv[i], "--mt="));
	else if (starts_with(argv[i],"--alpha_s="))
	  oneset.setAlpha(legacy::ALPHAS, get_value(argv[i], "--alpha_s="));      
	else if (starts_with(argv[i],"--alpha_inverse="))
	  oneset.setAlpha(legacy::ALPHA, 1.0 / get_value(argv[i],"--alpha_inverse="));
	else if (starts_with(argv[i],"--RPV")) 
	  RPVflag = true;
	else if (starts_with(argv[i], "--tanBeta=")) 
	  tanb = get_value(argv[i], "--tanBeta=");
	else if (starts_with(argv[i], "--sgnMu=")) 
	  sgnMu = get_valuei(argv[i], "--sgnMu=");
	else if (starts_with(argv[i], "--mgut=")) 
	  mgutGuess = mgutCheck(argv[i], gaugeUnification, ewsbBCscale); 
	else if (starts_with(argv[i], "--disable-two-loop-susy-thresholds")) {
#ifdef COMPILE_FULL_SUSY_THRESHOLD
	  USE_TWO_LOOP_THRESHOLD = false;
	  m.setAllTwoLoopThresholds(false);
#else
	  compilationProblem = true;
	  cout << "Two-loop thresholds not compiled.\n";
	  cout << "Please use the --two-loop-susy-thresholds with the configure option.\n";
	  cout << "Make sure you install the CLN and GiNaC packages beforehand.\n";
#endif
	}
	else if (starts_with(argv[i], "--two-loop-susy-thresholds")) {
#ifdef COMPILE_FULL_SUSY_THRESHOLD
	  USE_TWO_LOOP_THRESHOLD = true;
	  m.setAllTwoLoopThresholds(true);
#else
	  compilationProblem = true;
	  cout << "Two-loop thresholds not compiled.\n";
	  cout << "Please use the --enable-two-loop-susy-thresholds with the configure option.\n";
	  cout << "Make sure you install the CLN and GiNaC packages beforehand.\n";
#endif
	}
	else if (starts_with(argv[i], "--disable-three-loop-rges")) {
	  USE_THREE_LOOP_RGE = false;
	  compilationProblem = true;
	  cout << "Two-loop thresholds not compiled.\n";
	  cout << "Please use the --enable-two-loop-susy-thresholds with ./configure\n";
	  cout << "Make sure you install the CLN and GiNaC packages beforehand.\n";
	}
	else if (starts_with(argv[i], "--three-loop-rges")) {
	  USE_THREE_LOOP_RGE = true;
	  compilationProblem = true;
	  cout << "Three-loop RGEs not compiled.\n";
	  cout << "Please use the --enable-three-loop-rges with ./configure\n";
	}
	else if (starts_with(argv[i], "--QEWSB=")) 
	  QEWSB = get_value(argv[i], "--QEWSB=");
      }
      if (compilationProblem) exit(-1);
      if (tanb < 1.5 || tanb > 70.) {
	ostringstream ii; 
	ii << "tanBeta=" << tanb 
	   << " in SUGRA input. The point will not yield a sensible answer\n";
	throw ii.str();
      }
      
      /// Pass through to see if there are any RPV options
      if (strcmp(argv[1], "nmssm")) {
        for (int i = 2; i < argc; i++) {
          if (starts_with(argv[i], "--lambda")) {
            if (i + 4 >= argc) {
              throw "ERROR: three indices and one value need to be provided"
                " after --lambda or --lambdaP or --lambdaPP\n";
            }
            RPVflag = true;
            // int ii= int(atof(argv[i+1]));
            // int j = int(atof(argv[i+2]));
            // int k = int(atof(argv[i+3]));
            // double d = atof(argv[i+4]);
            // if (starts_with(argv[i], "--lambdaPP"))
            //   kw.setLambda(LU, ii, j, k, d);
            // else if (starts_with(argv[i], "--lambdaP"))
            //   kw.setLambda(LD, k, ii, j, d);
            // else if (starts_with(argv[i], "--lambda"))
            //   kw.setLambda(LE, k, ii, j, d);
          }
          if (starts_with(argv[i], "--kappa")) {
            if (i + 2 >= argc) {
              throw "ERROR: one index and one value need to be provided"
                " after --kappa\n";
            }
            // int ii = int(atof(argv[i+1]));
            // double d = atof(argv[i+2]);
            // kw.setKappa(ii, d);
            RPVflag = true;
          }
        }
      }
      
      /// Model specific options
      if (!strcmp(argv[1], "sugra")) {
	cout << "# SOFTSUSY SUGRA calculation" << endl;
	boundaryCondition = &sugraBcs;
	modelIdent = "sugra";
	double m0 = 0., m12 = 0., a0 = 0.;
	for (int i = 2; i < argc; i++) {
	  if (starts_with(argv[i], "--m0=")) m0 = get_value(argv[i], "--m0=");
	  else if (starts_with(argv[i], "--m12=")) 
	    m12 = get_value(argv[i], "--m12=");
	  else if (starts_with(argv[i], "--a0="))
	    a0  = get_value(argv[i], "--a0=");
	}
	pars(1) = m0; pars(2) = m12; pars(3) = a0;      
	if (m12 < MZ) {
	  ostringstream ii; 
	  ii << "m12=" << m12 
	     << " in SUGRA input. The point will not yield a sensible answer\n";
	  throw ii.str();
	}
	r = &m;
      }
      
      if (!strcmp(argv[1], "amsb")) {
	cout << "# SOFTSUSY mAMSB calculation" << endl;
	boundaryCondition = &amsbBcs;
	modelIdent = "amsb";
	double m0 = 0., m32 = 0.;
	for (int i = 2; i < argc; i++) {
	  if (starts_with(argv[i], "--m0=")) 
	    m0 = get_value(argv[i], "--m0=");
	  else if (starts_with(argv[i], "--m32=")) 
	    m32 = get_value(argv[i], "--m32=");
	}
	pars(1) = m0; pars(2) = m32; 
	if (m32 < 1.0e3) {
	  ostringstream ii; 
	  ii << "m32=" << m32 
	     << " in SUGRA input (too low). The point will not yield a sensible answer\n";
	  throw ii.str();
	}
	r = &m;
      }
      
      if (!strcmp(argv[1], "gmsb")) {
	cout << "# SOFTSUSY mGMSB calculation" << endl;
	boundaryCondition = &gmsbBcs;
	modelIdent = "gmsb";
	double n5 = 0., mMess = 0., LAMBDA = 0., cgrav = 1.;
	for (int i = 2; i < argc; i++) {
	  if (starts_with(argv[i], "--n5=")) n5 = get_value(argv[i], "--n5=");
	  else if (starts_with(argv[i], "--mMess=")) { 
	    gaugeUnification = false; 
	    mMess = get_value(argv[i], "--mMess=");
	    mgutGuess = mMess;
	  }
	  else if (starts_with(argv[i], "--LAMBDA=")) 
	    LAMBDA = get_value(argv[i], "--LAMBDA=");
	  else if (starts_with(argv[i], "--cgrav=")) 
	    cgrav = get_value(argv[i], "--cgrav=");
	}
	
	pars.setEnd(4);
	pars(1) = n5; pars(2) = mMess; pars(3) = LAMBDA; pars(4) = cgrav;
	if (mMess < 1.0e3) {
	  ostringstream ii; 
	  ii << " mMess=" << mMess
	     << " in SUGRA input (too low). The point will not yield a sensible answer\n";
	  throw ii.str();
	}
	
	r = &m;
	if (LAMBDA > mMess) {
	  ostringstream ii;
	  ii << "Input LAMBDA=" << LAMBDA << " should be less than mMess="
	     << mMess << endl;
	  throw ii.str();
	}
	if (cgrav > 1.0) {
	  ostringstream ii;
	  ii << "Input cgrav=" << cgrav << " a real number bigger than or "
	     << " equal to 1 (you can use 1 as a default value).\n";
	  throw ii.str();
	}
      }
    }
    if (!strcmp(argv[1], "nmssm")) {
      susy_model = NMSSM;
      NMSSM_command_line_parser nmssm_parser(&nmssm_input);
      nmssm_parser.parse(argc, argv);
      modelIdent = nmssm_parser.get_modelIdent();
      pars = nmssm_parser.get_pars();
    }
    
  if (!strcmp(argv[1], "leshouches")) {
	/// SLHA option "leshouches" used.
	outputCharacteristics(8);
	if (argc == 2) {
	  string line, block;
	  int model;
	  
	  while (getline(cin,line)) {
	    //	  mgutGuess = mgutCheck("unified", gaugeUnification);
	    
	    //	cout << line << endl;
	    istringstream input(line); 
	    string word1, word2;
	    input >> word1;
	    
	    if (word1.find("#") == string::npos &&
                !contains_only_whitespace(word1)) {
	      // read in another word if there's no comment
	      input >> word2; 
	      
	      if (ToUpper(word1) == "BLOCK")  { 
		block = ToUpper(word2);
		
	      } else { // ought to be data
		istringstream kk(line);
		if (block == "MODSEL") {
		  int i; kk >> i; 
		  
		  switch(i) {
		  case 1: kk >> model; 
		    switch(model) {
		    case 0: boundaryCondition = &extendedSugraBcs;
		      modelIdent = "nonUniversal"; r=&m;
		      break;
		    case 1: 
		      if (!flavourViolation) {
			pars.setEnd(3); 
			boundaryCondition = &sugraBcs; 
		      }
		      modelIdent = "sugra";
		      break;
		    case 2: 
		      if (!flavourViolation) {
			boundaryCondition = &gmsbBcs; 
			pars.setEnd(4); 
		      } 
		      modelIdent = "gmsb";
		      break;
		    case 3: 		    
		      boundaryCondition = &amsbBcs; 
		      pars.setEnd(2); 
		      modelIdent = "amsb";
		      break;
		    case 4:
		    boundaryCondition = &splitGmsb;
		    pars.setEnd(7); sgnMu = 0; 
		    modelIdent = "splitgmsb";
		    break;
		    default: 
		      ostringstream ii;
		      ii << "SOFTSUSY" << SOFTSUSY_VERSION 
			 << " cannot yet do model " 
			 << model << ": terminal error\n";
		      throw ii.str();
		    }
		    break;
                    // reading entry 3: susy model (MSSM, NMSSM, ...)
                  case 3: { int i; kk >> i;
                      switch(i) {
                      case 0: susy_model = MSSM; // default
                         break;
                      case 1: susy_model = NMSSM;
                         if (flavourViolation) {
                            flavourViolation = false;
                            cout << "# Warning: flavour violation is currtently"
                               " not supported in the NMSSM\n";
                         }
                         break;
                      default:
                         ostringstream ii;
                         ii << "MODSEL 3 choosing silly model switch\n"
                            << "(" << i << ") not a valid switch" << endl;
                         throw ii.str();
                      }
                    }
                    break;
		  case 4: int i; kk >> i;
		    switch(i) {
		    case 0: RPVflag = false;
		      break;
		    case 1: RPVflag = true;
		      break;
		    default:
		      ostringstream ii;
		      ii << "MODSEL 4 choosing silly RPV switch\n"
			 << "(" << i << ") not a valid switch" << endl;
		      throw ii.str();
		    }
		    break;
		  case 6: int j; kk >> j;
                    switch(j) {
                    case 0: flavourViolation = false; break;
                    default:
                       if (susy_model == NMSSM) {
                          flavourViolation = false;
                          cout << "# Warning: flavour violation is currtently"
                             " not supported in the NMSSM\n";
                       } else {
                          r = &k; 
                          flavourViolation = true;
                          if (boundaryCondition != & amsbBcs && 
			      boundaryCondition != & gmsbBcs) {
                             pars.setEnd(64);
                             boundaryCondition = &flavourBcs;
                          }
                       }
                    }
                    break;
		  case 11: kk >> numPoints;
		    if (numPoints < 1) {
		      ostringstream ii;
		      ii << "MODSEL 11 selecting silly number of points"
			 << "(" << numPoints << ") to output" << endl;
		      throw ii.str();
		    }
		    break;
		  case 12: double d; kk >> d;
		    if (d < MZ && d > 0.) {
		      ostringstream ii;
		      ii << "MODSEL 12 selecting silly scale Qmax"
			 << "(" << d << ") < MZ to output" << endl;
		      throw ii.str();
		    }
		    if (close(d + 1., 0., EPSTOL)) gutScaleOutput = true;
		    qMax = d; break;
		  default:
		    cout << "# WARNING: don't understand first integer " 
			 << word1 << " " << word2 << " in block " << block
			 << ": ignoring it\n";
		    break;
		  }
		}
		else if (block == "MINPAR") {
		  int i; double d; kk >> i >> d; 
		  switch (i) {
		  case 3: tanb = d;
                    nmssm_input.set(NMSSM_input::tanBeta, d);
                    break;
		  case 4: sgnMu = int(d); break;
		  default: 
		    switch(model) {
		    case 0:
		      // SUGRA inputs to fill out the pheno MSSM case
		      switch(i) {
		      case 1: pars(1) = d; break;
		      case 2: pars(2) = d; break;
		      case 5: pars(3) = d; break;
		      default: 
			ostringstream ii;
			ii << "Didn't understand pheno MSSM input " << i << endl;
			break;
		      } break;
		    case 1: // SUGRA inputs
		      switch(i) {
		      case 1: 
			if (flavourViolation) { pars.setEnd(77);
			  double m0 = sqr(d);
			  pars(4) = m0; pars(7) = m0; pars(9) = m0;
			  pars(10) = m0; pars(13) = m0; pars(15) = m0; 
			  pars(16) = m0; pars(19) = m0; pars(21) = m0; 
			  pars(22) = m0; pars(25) = m0; pars(27) = m0; 
			  pars(28) = m0; pars(31) = m0; pars(33) = m0; 
			  pars(63) = m0; pars(64) = m0;
			} else pars(1) = d; 
			break;
		      case 2: 
			if (flavourViolation) {
			  pars(1) = d; pars(2) = d; pars(3) = d;
			} else pars(2) = d; 
			break;
		      case 5: 
			if (flavourViolation) {
			  pars.setEnd(77);
			  pars(62) = d;
			  
			} else pars(3) = d; 
			break;
		      default: 
			ostringstream ii;
			ii << "Didn't understand SUGRA input " << i << endl;
			break;
		      } break;
		    case 2: // GMSB inputs
		      switch(i) {
		      case 1: pars(3) = d; break;
		    case 2: pars(2) = d; mgutGuess = d;
		      gaugeUnification = false; break;
		      case 5: pars(1) = d; break;
		      case 6: pars(4) = d; break;
		      default: 
			ostringstream ii;
			ii << "Didn't understand GMSB input " << i << endl;
			break;
		      } break;
		    case 3: ///< AMSB inputs
		      switch(i) {
		      case 1: pars(2) = d; break;
		      case 2: pars(1) = d; break;
		      default: 
			ostringstream ii;
		      ii << "Didn't understand AMSB input " << i << endl;
		      break;
		      } break;
		    case 4: ///< split GMSB inputs 
		      switch(i) {
		      case 1: pars(2) = d; break;
		      case 2: pars(3) = d; break;
		      case 5: pars(1) = d; break;
		      case 6: pars(7) = d; break;
		      case 7: pars(4) = d; mgutGuess = d; 
			gaugeUnification = false; break;
		      case 8: pars(5) = d; break;
		      case 9: pars(6) = d; m.useAlternativeEwsb(); 
			// kw.useAlternativeEwsb();
			break;
		      case 10: desiredMh = d; break;
		      default: 
			ostringstream ii;
			ii << "Didn't understand GMSB input " << i << endl;
			break;
		      } break;
		    default: 
		      ostringstream ii;
		    ii << "Didn't understand model input " << model << endl;
		    break;
		    }
		    break;
		  }
		}
		// Adding non-minimal options. 
		else if (block == "EXTPAR") {
                  int i; double d; kk >> i >> d;

                  // read extra NMSSM input parameters from EXTPAR
                  // (skipping NMSSM parameters if the MSSM was selected)
                  if (susy_model == MSSM) {
                     switch (i) {
                     case 61:
                     case 62:
                     case 63:
                     case 64:
                     case 65:
                     case 66:
                     case 67:
                     case 68:
                     case 69:
                     case 70:
                        cout << "# Warning: NMSSM parameter EXTPAR " << i
                             << " given but MSSM chosen -- ignoring it.\n";
                        continue;
                     }
                  } else if (susy_model == NMSSM) {
                     // read NMSSM susy parameters only and continue
                     switch (i) {
                     case 23: nmssm_input.set(NMSSM_input::mu     , d); 
		       continue;
                     case 61: nmssm_input.set(NMSSM_input::lambda , d); 
		       continue;
                     case 62: nmssm_input.set(NMSSM_input::kappa  , d); 
		       continue;
                     case 65: nmssm_input.set(NMSSM_input::lambdaS, d); 
		       continue;
                     case 66: nmssm_input.set(NMSSM_input::xiF    , d); 
		       continue;
                     case 68: nmssm_input.set(NMSSM_input::muPrime, d);
		       continue;
                     }
                  }
		  /// First, we want to convert our input to EXTPAR if we have
		  /// mSUGRA already
		  if (!strcmp(modelIdent, "sugra")) {
		    modelIdent = "nonUniversal";
		    if (!flavourViolation) {
		      /// We assume mSUGRA BCs with no flavour violation
		      r=&m; 
		      boundaryCondition = &extendedSugraBcs;
		      double m0 = pars(1), m12 = pars(2), a0 = pars(3);
		      pars.setEnd(49);
		      int i; for (i=1; i<=3; i++) pars(i) = m12;
		      for (i=11; i<=13; i++) pars(i) = a0;
		      pars(21) = m0*m0; pars(22) = m0*m0;
		      for (i=31; i<=36; i++) pars(i) = m0;		    
		      for (i=41; i<=49; i++) pars(i) = m0;		    
		      // kw.setNumRpcBcs(50); 
                      if (susy_model == NMSSM) {
                         pars.setEnd(56);
                         pars(50) = a0; // Alambda
                         pars(51) = a0; // Akappa
                         pars(52) = 0.; // mS'^2
                         pars(53) = m0*m0; // mS^2
                         pars(54) = 0.; // mu
                         pars(55) = 0.; // Bmu
                         pars(56) = 0.; // xiS
                      }
  		    } else {
		      /// This is flavour violation with EXTPAR: mSUGRA BCs
		      /// with flavour violation
		      // boundaryCondition = &flavourBcs;
		      if (pars.displayEnd() == 3) {
			double m0 = pars(1), m12 = pars(2), a0 = pars(3);
			double msq = m0 * m0;
			pars.setEnd(77);
			int i; for (i=1; i<=3; i++) pars(i) = m12;
			/// Fill in scalar mass squareds
			for (i=1; i<=5; i++) {
			  int num = (i-1) * 6 + 4;
			  pars(num)  = msq; 
			  pars(num+3)  = msq; 
			  pars(num+5)  = msq;
			}
			
			pars(62) = a0;
			pars(63) = msq; pars(64) = msq;
		      }
		      // kw.setNumRpcBcs(65);
		    }
		  }
		  
		  if (!strcmp(modelIdent, "nonUniversal")) {
		    /// First, put parameters that depend not on
		    /// flavoured/unflavoured input
		    if (i == 0) { 
		      mgutGuess = d;
		      gaugeUnification = false;
		      // setting Minput=-1 should yield MSSM BCs at MSUSY
		      if (fabs(d + 1.0) < EPSTOL) {
			mgutGuess = 1.0e3;
			ewsbBCscale = true;
			if (gaugeUnification) 
			  cout << "# Gauge unification ignored since pheno MSSM"
			       << " assumes BC set at QEWSB\n"; 
			gaugeUnification = false;
		      }
		    }
		    else if (i == 25) {
		      tanb = d;
		      if (!flavourViolation && pars.displayEnd() != 49) 
			pars.setEnd(49);
		      pars(i) = d;
		      r->setSetTbAtMX(true);
                      nmssm.setSetTbAtMX(true);
		    } 
		    else if (i == 23 || i == 26) {
		      m.useAlternativeEwsb();
		      k.useAlternativeEwsb();
		      if (i == 23) {
                         r->setMuCond(d); r->setSusyMu(d);
                         k.setMuCond(d); k.setSusyMu(d);
                         if (susy_model == NMSSM) {
                           if (pars.displayEnd() < 56) pars.setEnd(56);
                           nmssm_input.set(NMSSM_input::mu, d);
                           pars(54) = d;
                         }
                      }
		      if (i == 26) {
			r->setMaCond(d); k.setMaCond(d);
		      }
		    }
		    else if (!flavourViolation) {
		      if ((i > 0 && i <=  3) || (i >= 11 && i <= 13) || 
			  (i >= 21 && i <= 23) || (i == 26 || i == 25) 
			  || (i >= 31 && i <= 36) || 
			  (i >= 41 && i <= 49)) {
			if (pars.displayEnd() < 49) pars.setEnd(49);
			pars(i) = d;
                        if (susy_model == NMSSM) {
                           if (pars.displayEnd() < 56) pars.setEnd(56);
                           switch (i) {
                           case 21: nmssm_input.set(NMSSM_input::mHd2, d); break;
                           case 22: nmssm_input.set(NMSSM_input::mHu2, d); break;
                           case 23:
                             nmssm_input.set(NMSSM_input::mu, d);
                             pars(54) = d;
                             break;
                           }
                        }
  		      } else if ((61 <= i && i <= 70) || i == 24) {
			if (pars.displayEnd() < 56) pars.setEnd(56);
                        switch (i) {
                        case 24:
                           nmssm_input.set(NMSSM_input::BmuOverCosBetaSinBeta, d);
                           pars(55) = d;
                           break;
                        case 63:
                           nmssm_input.set(NMSSM_input::Alambda, d);
                           pars(50) = d;
                           break;
                        case 64:
                           nmssm_input.set(NMSSM_input::Akappa, d);
                           pars(51) = d;
                           break;
                        case 67:
                           nmssm_input.set(NMSSM_input::xiS, d);
                           pars(56) = d;
                           break;
                        case 69:
                           nmssm_input.set(NMSSM_input::mPrimeS2, d);
                           // setting pars(52) = B' = mS'^2 / mu'
                           if (nmssm_input.is_set(NMSSM_input::muPrime)) {
                              const double muPrime = nmssm_input.get(NMSSM_input::muPrime);
                              if (!close(muPrime, 0.0, EPSTOL))
                                 pars(52) = d / muPrime;
                           }
                           break;
                        case 70:
                           nmssm_input.set(NMSSM_input::mS2, d);
                           pars(53) = d;
                           break;
                        }
                      } else {
                        cout << "WARNING: did not understand parameter " 
                             << i << " in non-flavoured EXTPAR inputs\n";
                      }
		    } else {
		      /// Have to translate the numbers from SLHA to your
		      /// convention with flavour violation
		      if ((i > 0 && i < 4)) pars(i) = d; 
		      else if (i == 31) pars(22) = sqr(d);
		      else if (i == 32) pars(25) = sqr(d);
		      else if (i == 33) pars(27) = sqr(d);
		      else if (i == 34) pars(28) = sqr(d);
		      else if (i == 35) pars(31) = sqr(d);
		      else if (i == 36) pars(33) = sqr(d);
		      else if (i == 41) pars(4) = sqr(d);
		      else if (i == 42) pars(7) = sqr(d);
		      else if (i == 43) pars(9) = sqr(d);
		      else if (i == 44) pars(10) = sqr(d);
		      else if (i == 45) pars(13) = sqr(d);
		      else if (i == 46) pars(15) = sqr(d);
		      else if (i == 47) pars(16) = sqr(d);
		      else if (i == 48) pars(19) = sqr(d);
		      else if (i == 49) pars(21) = sqr(d);
		      else if (i == 21) pars(63) = d;
		      else if (i == 22) pars(64) = d;
		      else if (i > 10 && i < 14) 
			cout << "WARNING: At,Ab,Atau are for SLHA1 only. "
			     << "Setting them to zero.\n"
			     << "Please use blocks TUIN, TDIN, TEIN for "
			     << "flavour violating SLHA input\n";
		      else {
			cout << "WARNING: did not understand parameter " 
			     << i << " in flavoured EXTPAR inputs\n";
		      }
		    }
		  }
		}
                else if (block == "QEXTPAR") {
                  int i; double d; kk >> i >> d;
                  if (susy_model == NMSSM) {
                     switch (i) {
                     case 61: // scale where to input lambda
                        if (fabs(d + 1.0) < EPSTOL) {
                           softsusy::GUTlambda = false;
                        } else {
                           cout << "# WARNING: cannot input NMSSM parameter lambda"
                              " (set in QEXTPAR " << i << ") at a scale "
                              "different from M_susy.  Please set QEXTPAR "
                                << i << " to -1 (M_susy) or remove the entry.\n";
                        }
                        break;
                     case 62: // scale where to input kappa
                        if (fabs(d + 1.0) < EPSTOL) {
                           softsusy::GUTkappa = false;
                        } else {
                           cout << "# WARNING: cannot input NMSSM parameter kappa"
                              " (set in QEXTPAR " << i << ") at a scale "
                              "different from M_susy.  Please set QEXTPAR "
                                << i << " to -1 (M_susy) or remove the entry.\n";
                        }
                        break;
                     case 65: // scale where to input <S>
                        if (fabs(d + 1.0) < EPSTOL) {
                           softsusy::GUTsVev = false;
                        } else {
                           cout << "# WARNING: cannot input NMSSM parameter <S>"
                              " (set in QEXTPAR " << i << ") at a scale "
                              "different from M_susy.  Please set QEXTPAR "
                                << i << " to -1 (M_susy) or remove the entry.\n";
                        }
                        break;
                     case 66: // scale where to input xiF
                        if (fabs(d + 1.0) < EPSTOL) {
                           softsusy::GUTxiF = false;
                        } else {
                           cout << "# WARNING: cannot input NMSSM parameter xiF"
                              " (set in QEXTPAR " << i << ") at a scale "
                              "different from M_susy.  Please set QEXTPAR "
                                << i << " to -1 (M_susy) or remove the entry.\n";
                        }
                        break;
                     case 68: // scale where to input mu'
                        if (fabs(d + 1.0) < EPSTOL) {
                           softsusy::GUTmuPrime = false;
                        } else {
                           cout << "# WARNING: cannot input NMSSM parameter mu'"
                              " (set in QEXTPAR " << i << ") at a scale "
                              "different from M_susy.  Please set QEXTPAR "
                                << i << " to -1 (M_susy) or remove the entry.\n";
                        }
                        break;
                     default:
                        cout << "# WARNING: cannot use parameter " << i <<
                           " (set in QEXTPAR) as input at a different"
                           " scale (in the NMSSM) -- ignoring the scale choice\n";
                     }
                     continue;
                  }
                  cout << "# WARNING: cannot use parameter " << i <<
                     " (set in QEXTPAR) as input at a different"
                     " scale -- ignoring the scale choice\n";
                }
		else if (block == "VCKMIN") {
		  int i; double d; kk >> i >> d;
		  switch(i) {
		  case 1: lambdaW = d; break;
		  case 2: aCkm = d;   break;		  
		  case 3: rhobar = d; break;
		  case 4: etabar = d; break;
		  default:
		    cout << "# WARNING: Don't understand data input " << i 
			 << " " << d << " in block "
			 << block << ": ignoring it\n";
		    break;
		  }
		}
		else if (block == "UMNSIN") {
		  int i; double d; kk >> i >> d;
		  switch(i) {
		  case 1: k.setThetaB12(asin(d)); break;
		  case 2: k.setThetaB23(asin(d)); break;
		  case 3: k.setThetaB13(asin(d)); break;
		  case 4: cout << "# Cannot yet do complex phases: ";
		    cout << "setting it to zero" << endl; break;
		  case 5: cout << "# Cannot yet do complex phases: ";
		    cout << "setting it to zero" << endl; break;
		  case 6: cout << "# Cannot yet do complex phases: ";
		    cout << "setting it to zero" << endl; break;
		  default:
		    cout << "# WARNING: Don't understand data input " << i 
			 << " " << d << " in block "
			 << block << ": ignoring it\n";
		    break;
		  }
		}
		else if (block == "SMINPUTS") {
		  int i; double d; kk >> i >> d; 
		  switch (i) {
		  case 1: oneset.setAlpha(legacy::ALPHA, 1.0 / d); break;
		  case 2: GMU = d; break;
		  case 3: oneset.setAlpha(legacy::ALPHAS, d); break; 
		  case 4: oneset.setMu(d); m.setData(oneset); MZ = d; break;
		  case 5: oneset.setMass(legacy::mBottom, d); 
		    oneset.setMbMb(d); break;
		  case 6: oneset.setPoleMt(d); break;
		  case 7: oneset.setMass(legacy::mTau, d); 
		    oneset.setPoleMtau(d); break;
		  case 8: k.setMnuTau(d);
                     break;
		  case 11: oneset.setMass(legacy::mElectron, d);
                     k.setPoleMe(d);
                     break;
		  case 12: k.setMnuMu(d);
                     break;
		  case 13: oneset.setMass(legacy::mMuon, d); k.setPoleMmu(d);
                     break;
		  case 14: k.setMnuTau(d);
                     break;
		  case 21: oneset.setMass(legacy::mDown, d); k.setMd2GeV(d);
		    break;
		  case 22: oneset.setMass(legacy::mUp, d); k.setMu2GeV(d);
                     break;
		  case 23: oneset.setMass(legacy::mStrange, d); k.setMs2GeV(d);
                     break;
		  case 24: oneset.setMass(legacy::mCharm, d); k.setMcMc(d);
                     break;
		  default: 
		    cout << "# WARNING: Don't understand data input " << i 
			 << " " << d << " in block "
			 << block << ": ignoring it\n"; break;
		  }
		} 
		else if (block == "MSQ2IN") {
		  modelIdent = "nonUniversal";
		  // int i, j; double d; kk >> i >> j >> d;
		  // pars(positionOfSym(i, j) + 3) = d;
		}
		else if (block == "MSU2IN") {
		  modelIdent = "nonUniversal";
		  // int i, j; double d; kk >> i >> j >> d;
		  // pars(positionOfSym(i, j) + 9) = d;
		}
		else if (block == "MSD2IN") {
		  modelIdent = "nonUniversal";
		  // int i, j; double d; kk >> i >> j >> d;
		  // pars(positionOfSym(i, j) + 15) = d;
		}
		else if (block == "MSL2IN") {
		  modelIdent = "nonUniversal";
		  // int i, j; double d; kk >> i >> j >> d;
		  // pars(positionOfSym(i, j) + 21) = d;
	      }
		else if (block == "MSE2IN") {
		  modelIdent = "nonUniversal";
		  // int i, j; double d; kk >> i >> j >> d;
		  // pars(positionOfSym(i, j) + 27) = d;
		}
		else if (block == "TUIN") {
		  modelIdent = "nonUniversal";
		  // int i, j; double d; kk >> i >> j >> d;
		  // pars((i-1) * 3 + j + 33) = d;
		  // slha2setTrilinear[(i-1) * 3 + j - 1] = true;
		}
		else if (block == "TDIN") {
		  modelIdent = "nonUniversal";
		  // int i, j; double d; kk >> i >> j >> d;
		  // pars((i-1) * 3 + j + 42) = d;
		  // slha2setTrilinear[(i-1) * 3 + j + 8] = true;
		}
		else if (block == "TEIN") {
		  modelIdent = "nonUniversal";
		  // int i, j; double d; kk >> i >> j >> d;
		  // pars((i-1) * 3 + j + 51) = d;
		  // slha2setTrilinear[(i-1) * 3 + j + 17] = true;
		}
		else if (block == "RVLAMLLEIN") {
		  // int i,j,k; double d; kk >> i >> j >> k >> d;
		  // if ((i > 0 && i <=  3) || (j > 0 && j <=  3) ||
		  //     (k > 0 && k <=  3)) {
		  //   kw.setLambda(LE, k, i, j, d);
		  // }
		  // else {
		  //   cout << "# WARNING: Don't understand data input " << i 
		  //        << " " << j << " " << k << " " << d << " in block "
		  //        << block << ": ignoring it\n"; break;
		  // }
		}
		else if (block == "RVLAMLQDIN") {
		  // int i,j,k; double d; kk >> i >> j >> k >> d;
		  // if((i > 0 && i <=  3) || (j > 0 && j <=  3) ||
		  //    (k > 0 && k <=  3)) {
		  //   kw.setLambda(LD, k, i, j, d);
		  // }
		  // else {
		  //   cout << "# WARNING: Don't understand data input " << i 
		  //        << " " << j << " " << k << " " << d << " in block "
		  //        << block << ": ignoring it\n"; break;
		  // }
		}
		else if (block == "RVLAMUDDIN") {
		  // int i,j,k; double d; kk >> i >> j >> k >> d;
		  // if((i > 0 && i <=  3) || (j > 0 && j <=  3) || 
		  //    (k > 0 && k <=  3)) {
		  //   kw.setLambda(LU, i, j, k, d);
		  // }
		}
		else if (block == "RVTIN") {
		  // int i,j,k; double d; kk >> i >> j >> k >> d;
		  // if((i > 0 && i <=  3) || (j > 0 && j <=  3) || 
		  //    (k > 0 && k <=  3)) {
		  //   kw.setHr(LE, k, i, j, d);
		  // }
		  // else {
		  //   cout << "# WARNING: Don't understand data input " << i 
		  //        << " " << j << " " << k << " " << d << " in block "
		  //        << block << ": ignoring it\n"; break;
		  // }
		}
		else if (block == "RVTPIN") {
		  // int i,j,k; double d; kk >> i >> j >> k >> d;
		  // if((i > 0 && i <=  3) || (j > 0 && j <=  3) || 
		  //    (k > 0 && k <=  3)) {
		  //   kw.setHr(LD, k, i, j, d);
		  // }
		  // else {
		  //   cout << "# WARNING: Don't understand data input " << i 
		  //        << " " << j << " " << k << " " << d << " in block "
		  //        << block << ": ignoring it\n"; break;
		  // }
		}
		else if (block == "RVTPPIN") {
		  // int i,j,k; double d; kk >> i >> j >> k >> d;
		  // if((i > 0 && i <=  3) || (j > 0 && j <=  3) ||
		  //    (k > 0 && k <=  3)) {
		  //   kw.setHr(LU, i, j, k, d);
		  // }
		  // else {
		  //   cout << "# WARNING: Don't understand data input " << i 
		  //        << " " << j << " " << k << " " << d << " in block "
		  //        << block << ": ignoring it\n"; break;
		  // }
		}
		else if (block == "RVKAPPAIN") {
		  // int i; double d; kk >> i >> d;
		  // if (i > 0 && i <=  3) {
		  //   kw.setKappa(i, d);
                  // }
                  // else {
                  //   cout << "# WARNING: Don't understand data input " << i << d 
                  //        << " in block " << block << ": ignoring it\n"; break;
                  // }
		}
		else if (block == "RVDIN") {
		  // int i; double d; kk >> i >> d;
		  // if (i > 0 && i <=  3) {
		  //   kw.setD(i, d);
		  // }
		  // else {
		  //   cout << "# WARNING: Don't understand data input " 
		  //        << i << d << " in block " << block 
		  //        << ": ignoring it\n"; break;
		  // }
		}
		// input of sneutrino VEVs not supported yet
		else if (block == "RVSNUVEVIN") {
		  // int i; double d; kk >> i >> d;
		  // cout << "# WARNING: in block " << block 
		  //      << ": SOFTSUSY does not support setting of sneutrino VEVs"
		  //      << " yet : ignoring them\n"; break;
		}
		else if (block == "RVMLH1SQIN") {
		  // int i; double d; kk >> i >> d;
		  // if ((i > 0 && i <=  3)) {
		  //   kw.setMh1lSquared(i, d);
		  // }
		  // else {
		  //   cout << "# WARNING: Don't understand data input " 
		  //        << i << d << " in block " << block 
		  //        << ": ignoring it\n"; break;
		  // }
		}
		else if (block == "SOFTSUSY") {
		  int i; double d; kk >> i >> d;
		  switch(i) {
		  case 1: TOLERANCE = d; break;
		  case 2: 
		    MIXING = int(d); 
		    //if (MIXING > 0) flavourViolation = true;
		    break;
		  case 3: PRINTOUT = int(d); break;
		  case 4: QEWSB = d; break;
		  case 5: INCLUDE_2_LOOP_SCALAR_CORRECTIONS = 
		      bool(int(d+EPSTOL)); break;
		  case 6: outputCharacteristics(int(d+EPSTOL)-1); break;  
		  case 7: {
		    int num = int(d+EPSTOL);
		    if (num != 1 && num!= 2) {
		      cout << "# WARNING: Can only set number of loops for"
			   << " higgs masses and REWSB to be 1 or 2 in "
			   << " BLOCK SOFTSUSY parameter 7, not " << num 
			   << ". Ignoring.\n";
		  } else {
		      numHiggsMassLoops = num;
		      numRewsbLoops = num;
		    }
		  }
		    break;
		  case 8: {
		  int num = int(d + EPSTOL);
		  if (num != 0 && num != 1) {
		    cout << "# Incorrect flag in BLOCK SOFTSUSY 8: " 
			 << num << ". Ignoring it.\n";
		  } else {
		    if (num == 1) susyRpvBCatMSUSY = true;
		    else susyRpvBCatMSUSY = false;
		  }
		  }
		    break;
		  case 9: {
		    // int num = int(d + EPSTOL);
		    // if (num != 0 && num != 1) {
		    //   cout << "# Incorrect flag in BLOCK SOFTSUSY 9: "
		    //        << num << ". Ignoring it.\n";
		    // } else {
		    //   if (num == 1) kw.setInvertedOutput();
		    //   else kw.setNormalOutput();
                    // }
		  }
		    break;
		  case 10: {
		    int num = int(d + EPSTOL);
		    if (num != 0 && num != 1) {
		      cout << "# Incorrect flag in BLOCK SOFTSUSY 10: " 
			   << num << ". Ignoring it.\n";
		    } else if (num == 1) forceSlha1 = true;
		    else forceSlha1 = false;
		  }
		    break;
		  case 11: {
		    /// Set gravitino mass
		    r->setM32(d); 
		  }
		    break;
		  case 12: {
		    int num = int(d + EPSTOL);
		    if (num == 1) printRuledOutSpectra = true;
		    else if (num == 0) printRuledOutSpectra = false;
		    else cout << "# WARNING: Don't understand Block SOFTSUSY "
			      << "parameter 12 " << d << ". Ignoring it." 
			      << endl;
		  }
		    break;
		  case 13: {
		    int num = int(d + EPSTOL);
		    if (num == 1) mAFlag = true;		  
		  }
                    break;
                  case 15: {
                    int num = int(d + EPSTOL);
                    softsusy::NMSSMTools = num;
                  }
		    break;
                  case 16: {
                     int num = int(d + EPSTOL);
                     softsusy::MICROMEGAS = num;
                  }
                    break;
                  case 17: {
                     int num = int(d + EPSTOL);
                     softsusy::NMSDECAY = num;
                  }
		    break;
                  case 18: {
                    int num = int(d + EPSTOL);
                    if(num == 1) softsusy::SoftHiggsOut = true;
                  }
                     break;
		  case 19: {
                    int num = int(d + EPSTOL);
		    if (num == 1) USE_THREE_LOOP_RGE = true;
		    else if (num == 0) USE_THREE_LOOP_RGE = false;
		    else cout << "WARNING: incorrect setting for SOFTSUSY Block 19 (should be 0 or 1)\n";
		    break;			     
		  }
#ifdef COMPILE_FULL_SUSY_THRESHOLD
		  case 20: {
                    int num = int(d + EPSTOL);
		    // AVB: can be set to just 1 Turn on all thresholds
		    //      can be set to 1 + 2 * ( flags for included thresholds)
		    //      to have a finer control over included thresholds
		    if (num > 0) {
		      USE_TWO_LOOP_THRESHOLD = true;
		      r->included_thresholds = (num & 
						    (ENABLE_TWO_LOOP_AS_AS_YUK | 
						     ENABLE_TWO_LOOP_MT_AS | 
						     ENABLE_TWO_LOOP_MB_AS | 
						     ENABLE_TWO_LOOP_MB_YUK | 
						     ENABLE_TWO_LOOP_MTAU_YUK));
				   
		    } else if (num == 0) { 
		      USE_TWO_LOOP_THRESHOLD = false;  
		      r->included_thresholds = 0; 
		    } else cout << "WARNING: incorrect setting for SOFTSUSY Block 20 (should be an integer number in range 0,...,31)\n";
		    break;
		  }
#endif
		  case 21: {
		    int num = int(d + EPSTOL);
		    if (num > 0) higgsUncertainties = true;
		    break;
		  }
		  default:
		    cout << "# WARNING: Don't understand data input " << i 
			 << " " << d << " in block "
			 << block << ": ignoring it\n"; break;
		  }
		}
		else {
		  cout << "# WARNING: don't recognise block " << block 
		       << ": ignoring all data in it" << endl;
		}
		// end if blocks
	      
	      } // end of data
	    } // end of no-comment
            
	} // end of file
          
        }
	else errorCall();
      }

      /// prepare CKM angles
      if (flavourViolation || RPVflag) k.setAngles(lambdaW, aCkm, rhobar,
        					   etabar);
      
      if (r->displayAltEwsb()) {
	if (strcmp(modelIdent, "splitgmsb")) {
	  //	boundaryCondition = &extendedSugraBcs2;
	  r->setSusyMu(pars(23)); 
	} else {
	  ostringstream ii;
	  ii << "Split GMSB BCs should not supported with alternative EWSB\n";
	  throw ii.str();
	}
	sgnMu = 0; // Flags different BCs
      }

      // set NMSSM boundary conditions
      if (susy_model == NMSSM) {
         softsusy::Z3 = nmssm_input.is_Z3_symmetric();

         if (flavourViolation) {
            string msg("# Error: flavour violation in the NMSSM is currenty"
                       " not supported\n");
            throw msg;
         }
         if (strcmp(modelIdent, "sugra") == 0) {
           // if we chose the soft Higgs masses as EWSB output, we
           // must not use msugra, because it would overwrite mHu2,
           // mHd2, mS2 at the GUT scale
           if (softsusy::SoftHiggsOut) {
             if (pars.size() != 6)
               pars.setEnd(6);
             pars(4) = nmssm_input.get(NMSSM_input::mu);
             pars(5) = nmssm_input.get(NMSSM_input::BmuOverCosBetaSinBeta);
             pars(6) = nmssm_input.get(NMSSM_input::xiS);
             nmssmBoundaryCondition = &NmssmSugraNoSoftHiggsMassBcs;
           } else {
             if (softsusy::Z3) {
               // Here we must use SemiMsugraBcs to avoid setting mS2
               // at the GUT scale
               if (pars.size() != 5)
                 pars.setEnd(5);
               pars(4) = pars(3); // sets Al to A0
               pars(5) = pars(3); // sets Ak to A0
               nmssmBoundaryCondition = &SemiMsugraBcs;
             } else {
               nmssmBoundaryCondition = &NmssmMsugraBcs;
             }
           }
         } else if (strcmp(modelIdent, "nonUniversal") == 0) {
           nmssmBoundaryCondition = &extendedNMSugraBcs;
           if (pars.size() != 56) {
              string msg("# Error: NMSSM non-minmal sugra boundary condition"
                         " chosen, but pars does not have 56 entries\n");
              throw msg;
           }
           if (nmssm_input.is_set(NMSSM_input::mu))
             pars(54) = nmssm_input.get(NMSSM_input::mu);
           if (nmssm_input.is_set(NMSSM_input::BmuOverCosBetaSinBeta))
             pars(55) = nmssm_input.get(NMSSM_input::BmuOverCosBetaSinBeta);
           if (nmssm_input.is_set(NMSSM_input::xiS))
             pars(56) = nmssm_input.get(NMSSM_input::xiS);
         } else {
            string msg("# Error: non-sugra boundary conditions for the NMSSM"
                       " are currently not supported\n");
            throw msg;
         }
      }

      if (RPVflag) {
	// kw.rpvDisplay(pars);
	// kw.setFlavourSoftsusy(k);
	// r = &kw;
        //
	// if (boundaryCondition == &sugraBcs)
	//   boundaryCondition = &rpvSugraBcs;
	// else if (boundaryCondition == &amsbBcs)
	//   boundaryCondition = &rpvAmsbBcs;
	// else if (boundaryCondition == &gmsbBcs)
	//   boundaryCondition = &rpvGmsbBcs;
	// else if (boundaryCondition == &extendedSugraBcs) {
	//   boundaryCondition = &rpvExtendedSugraBcs;
	// }
	// else if (boundaryCondition == &extendedSugraBcs) {
	//   boundaryCondition = &rpvExtendedSugraBcs;
	//   kw.useAlternativeEwsb();
	// }
	// else {
	//   ostringstream ii;
	//   ii << "# ERROR: there is no RPV version for selected "
	//      << " boundary condition. If you have" << endl
	//      << "# *BOTH* Flavour violating *AND* RPV flags switched on, "
	//      << " remove the flavour" << endl
	//      << " violating one." << endl;
	//   throw ii.str();
	// }
      }
      // intput error checking  
      if (sgnMu != 1 && sgnMu != -1 && sgnMu != 0) {
	ostringstream ii;
	ii << "Incorrect input for sign(mu)=" << sgnMu <<endl;
	throw ii.str();
      }
      if (tanb < 1.0 || tanb > 5.0e2)  {
	ostringstream ii;
	ii << "Incorrect input for tan beta=" << tanb <<endl;
	throw ii.str();
      }
      
      oneset.toMz();
      /// Set quark mixing correctly
      // r->setMixing(MIXING); k.setMixing(MIXING); nmssm.setMixing(MIXING);
      // r->setQewsb(qewsb); k.setQewsb(qewsb); nmssm.setQewsb(qewsb);
      if (USE_THREE_LOOP_RGE) {
	r->setLoops(3); k.setLoops(3); nmssm.setLoops(3);
      }
    switch (susy_model) {
    case MSSM:
      r->fixedPointIteration(boundaryCondition, mgutGuess, pars, sgnMu, tanb, 
			     oneset, gaugeUnification, ewsbBCscale);
      
      /// Fix to mh if additional operators are assumed
      if (desiredMh > 0.1) {
        sPhysical s(r->displayPhys()); s.mh0(1) = desiredMh; r->setPhys(s);
      }
      
      if (gutScaleOutput) qMax = r->displayMxBC();

      r->lesHouchesAccordOutput(cout, modelIdent, pars, sgnMu, tanb, qMax,  
				numPoints, ewsbBCscale);

      // if (higgsUncertainties) {
      //   int numPts = 30;
      //   DoubleVector mh(numPts), mH(numPts), mA(numPts), mHp(numPts);
      //   double dmh = 0., dmH = 0., dmA = 0., dmHp = 0.;
      //   double lnqMin = log(MZ), 
      //     lnqMax = log(2.0 * r->displayMsusy());
      //   for (int i = 0; i< numPts; i++) {
      //     MssmSoftsusy a(r->displaySoftsusy());
      //     double lnq = (lnqMax - lnqMin) * double(i) / double(numPts) + lnqMin;
      //     double q = exp(lnq);
      //     int accuracy = 3; 
      //     double mt = 0., sinth = 0., piww = 0., pizz = 0;
      //     a.calcHiggsAtScale(accuracy, mt, sinth, piww, pizz, q);
      //     if (!a.displayProblem().testSeriousProblem()) {
      //       mh.set(i+1, a.displayPhys().mh0(1));
      //       mH.set(i+1, a.displayPhys().mh0(2));
      //       mA.set(i+1, a.displayPhys().mA0(1));
      //       mHp.set(i+1, a.displayPhys().mHpm);
      //     }
      //    }
      //   /// Try to add some fixed-order uncertainties to this
      //   int p;
      //   dmh  =  mh.max(p) -  mh.min(p);
      //   dmA  =  mA.max(p) -  mA.min(p);
      //   dmH  =  mH.max(p) -  mH.min(p);
      //   dmHp = mHp.max(p) - mHp.min(p);
      //   cout << "# Estimated Higgs uncertainties in GeV:\n# Dmh=" 
      //        << dmh << " DmH=" << dmH << " DmA=" << dmA 
      //        << " DmH+-=" << dmHp << endl;
      // }
      
      if (r->displayProblem().test()) {
	cout << "# SOFTSUSY problem with point: " << r->displayProblem() 
	     << endl;
        return -1;
      }
      break;
    case NMSSM: {
      nmssm_input.check_setup();

      DoubleVector nmpars(nmssm_input.get_nmpars());
      nmssm.lowOrg(nmssmBoundaryCondition, mgutGuess, pars, nmpars, sgnMu,
                   tanb, oneset, gaugeUnification, ewsbBCscale);
      nmssm.lesHouchesAccordOutput(cout, modelIdent, pars, sgnMu, tanb, qMax,
                                   numPoints, ewsbBCscale);
      if (nmssm.displayProblem().test()) {
         cout << "# SOFTSUSY problem with NMSSM point: "
              << nmssm.displayProblem() << endl;
         return -1;
      }
      }
      break;
    default:
      cout << "# Error: unknown susy model " << susy_model
           << ", please check your MODSEL (entry 3) settings" << endl;
      break;
    }
  }
  catch(const string & a) { cout << a; return -1; }
  catch(const char * a) { cout << a; return -1; }
  catch(...) { cout << "Unknown type of exception caught.\n"; return -1; }
  
  return 0;
}

