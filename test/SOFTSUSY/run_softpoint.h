
/** \file softpoint.h
   - Project:     SOFTSUSY 
   - Author:      Ben Allanach, Markus Bernhardt 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/
   - Description: main calling program: command line interface
   \brief Main program: reads input in SLHA or command line format
*/

#include <iostream>
#include <sstream>
#include <cstring>
#include "mycomplex.h"
#include "def.h"
#include "linalg.h"
#include "lowe_legacy.h"
#include "rge.h"
#include "softsusy.h"
#include "flavoursoft.h"
#include "nmssmsoftsusy.h"
#include "softpars.h"
#include "physpars.h"
#include "susy.h"
#include "utils.h"
#include "numerics.h"



/// string routine for options
bool starts_with(const std::string& str,
		 const std::string& prefix) {
  return !str.compare(0, prefix.size(), prefix);
}

double get_value(const std::string& str, const std::string& prefix) {
   return atof(str.substr(prefix.size()).c_str());
}

int get_valuei(const std::string& str, const std::string& prefix) {
   return atoi(str.substr(prefix.size()).c_str());
}

bool contains_only_whitespace(const std::string& str) {
  return str.find_first_not_of(" \t") == std::string::npos;
}

namespace softsusy {
   /// Requested by CMS
   void splitGmsb(MssmSoftsusy & m, const DoubleVector & inputParameters);
}

/// Does the user require gauge unification or not -- gaugeUnification changed
/// to be correct value
inline double mgutCheck(char * a, bool & gaugeUnification, 
			bool & ewsbBCscale) { 
  gaugeUnification = false; ewsbBCscale = false;
  if (!strcmp(a, "--mgut=?") || !strcmp(a,"--mgut=unified")) {
    gaugeUnification = true; 
    return 2.0e16;
  }
  if (!strcmp(a, "--mgut=msusy")) {
    ewsbBCscale = true;
    return 1.0e3;
  }
  else return get_value(a, "--mgut");
}

/// Incorrect input: gives advice on how to supply it
void errorCall();
