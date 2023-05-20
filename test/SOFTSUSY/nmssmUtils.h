/** \file nmssmUtils.h
    - Project: Next-to-minimal SOFTSUSY 
    - Author: Ben Allanach, Peter Athron, Lewis Tunstall, Alexander Voigt,
      Anthony Williams
    - Manual: http://arxiv.org/abs/1311.7659
    - Webpage: http://hepforge.cedar.ac.uk/softsusy/

    \brief A few handy routines for the NMSSM: SUSY breaking conditions etc
*/

#ifndef NMSSM_UTILS_H
#define NMSSM_UTILS_H

#include <iosfwd>
#include <string>

#include "linalg.h"

namespace softsusy {

class NmssmSoftsusy;
class SoftParsNmssm;
template <class SoftPars> class Softsusy;

/// class for NMSSM input parameters
class NMSSM_input {
public:
   enum NMSSM_parameters {
      tanBeta,               // MINPAR entry 3
      mHd2,                  // EXTPAR entry 21
      mHu2,                  // EXTPAR entry 22
      mu,                    // EXTPAR entry 23
      BmuOverCosBetaSinBeta, // m3^2/(cos(Beta)sin(Beta)) EXTPAR entry 24
      lambda,                // EXTPAR entry 61
      kappa,                 // EXTPAR entry 62
      Alambda,               // EXTPAR entry 63
      Akappa,                // EXTPAR entry 64
      lambdaS,               // lambda * <S> EXTPAR entry 65
      xiF,                   // EXTPAR entry 66
      xiS,                   // EXTPAR entry 67
      muPrime,               // EXTPAR entry 68
      mPrimeS2,              // EXTPAR entry 69
      mS2,                   // EXTPAR entry 70
      NUMBER_OF_NMSSM_INPUT_PARAMETERS
   };

   NMSSM_input();
   ~NMSSM_input() {}

   /// set parameter to given value
   void set(NMSSM_parameters, double);
   /// get value of parameter
   double get(NMSSM_parameters) const;
   /// returns vector with supersymmetric NMSSM parameters
   DoubleVector get_nmpars() const;
   /// returns true if parameter was set, false otherwise
   bool is_set(NMSSM_parameters par) const;
   /// returns true if input parameter set defines a Z3 symmetric NMSSM
   bool is_Z3_symmetric() const;
   /// checks the NMSSM parameter setup, throws if not SLHA2 conform
   void check_setup();

   friend std::ostream& operator<<(std::ostream&, const NMSSM_input&);

private:
   double parameter[NUMBER_OF_NMSSM_INPUT_PARAMETERS];  ///< NMSSM parameters
   bool has_been_set[NUMBER_OF_NMSSM_INPUT_PARAMETERS];
   static char const * const parameter_names[NUMBER_OF_NMSSM_INPUT_PARAMETERS];

   /// checks if the unset parameters can be used as EWSB output
   void check_ewsb_output_parameters() const;
};

/// class which parses command line options of the NMSSM
class NMSSM_command_line_parser {
public:
   NMSSM_command_line_parser(NMSSM_input*);
   ~NMSSM_command_line_parser() {}

   /// parse command line options
   void parse(int argc, char* argv[]);
   /// returns model identification string
   const char* get_modelIdent() const;
   /// returns vector of GUT parameters
   DoubleVector get_pars() const;

private:
   NMSSM_input* nmssm_input;  ///< pointer to NMSSM input parameters
   const char* model_ident;   ///< model identification string
   double m0, m12, a0;        ///< GUT parameters

   /// returns floating point number from string like --m0=125.4
   static double get_value(const std::string& str, const std::string& prefix);
   /// checks if a string starts with a certain prefix
   static bool starts_with(const std::string& str, const std::string& prefix);
};

void NmssmMsugraBcs(NmssmSoftsusy & m, const DoubleVector & inputParameters);

void MssmMsugraBcs(NmssmSoftsusy & m, const DoubleVector & inputParameters);

void SemiMsugraBcs(NmssmSoftsusy & m, const DoubleVector & inputParameters);

void NmssmSugraNoSoftHiggsMassBcs(NmssmSoftsusy&, const DoubleVector&);

void NmssmSemiMsugraNoSoftHiggsMassBcs(NmssmSoftsusy&, const DoubleVector&);

void extendedNMSugraBcs(NmssmSoftsusy & m, const DoubleVector & inputParameters);

/// Calculates fractional difference in Drbar masses between in and out
double sumTol(const Softsusy<SoftParsNmssm> & in, const Softsusy<SoftParsNmssm> & out, int numTries);

} // namespace softsusy

#endif
