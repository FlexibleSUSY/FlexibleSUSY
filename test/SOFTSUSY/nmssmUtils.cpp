
#include "nmssmUtils.h"
#include "nmssmsoftsusy.h"
#include "def.h"
#include "utils.h"

#include <cassert>
#include <sstream>

namespace softsusy {

char const * const NMSSM_input::parameter_names[NUMBER_OF_NMSSM_INPUT_PARAMETERS] = {
   "tan(beta)", "mHd^2", "mHu^2", "mu", "m3^2/(cos(beta)sin(beta))", "lambda",
   "kappa", "Alambda", "Akappa", "lambda*S", "xiF", "xiS", "mu'",
   "mS'^2", "mS^2"
};

NMSSM_input::NMSSM_input()
  : parameter()    // sets all parameters to zero
  , has_been_set() // sets all values to zero (false)
{}

void NMSSM_input::set(NMSSM_parameters par, double value) {
   assert(par < NUMBER_OF_NMSSM_INPUT_PARAMETERS);
   parameter[par] = value;
   has_been_set[par] = true;
}

double NMSSM_input::get(NMSSM_parameters par) const {
   assert(par < NUMBER_OF_NMSSM_INPUT_PARAMETERS);
   return parameter[par];
}

DoubleVector NMSSM_input::get_nmpars() const {
   DoubleVector nmpars(5);
   nmpars(1) = get(NMSSM_input::lambda);
   nmpars(2) = get(NMSSM_input::kappa);
   if (is_set(NMSSM_input::lambdaS)) {
      if (!close(get(NMSSM_input::lambda),0.,EPSTOL)) {
         nmpars(3) = get(NMSSM_input::lambdaS) * sqrt(2.)
            / get(NMSSM_input::lambda);
      } else {
         std::string msg =
            "# Error: you set lambda * <S> to a non-zero value"
            ", but lambda is zero.  "
            "Please set lambda (EXTPAR entry 61) to a non-zero value.\n";
         throw msg;
      }
   }
   nmpars(4) = get(NMSSM_input::xiF);
   nmpars(5) = get(NMSSM_input::muPrime);
   return nmpars;
};

bool NMSSM_input::is_set(NMSSM_parameters par) const {
   assert(par < NUMBER_OF_NMSSM_INPUT_PARAMETERS);
   return has_been_set[par];
}

bool NMSSM_input::is_Z3_symmetric() const {
   return (!is_set(mu)      || close(parameter[mu]      , 0., EPSTOL))
      && (!is_set(BmuOverCosBetaSinBeta)
          || close(parameter[BmuOverCosBetaSinBeta], 0., EPSTOL))
      && (!is_set(muPrime)  || close(parameter[muPrime] , 0., EPSTOL))
      && (!is_set(mPrimeS2) || close(parameter[mPrimeS2], 0., EPSTOL))
      && (!is_set(xiF)      || close(parameter[xiF]     , 0., EPSTOL));
}

void NMSSM_input::check_setup() {
   check_ewsb_output_parameters();
}

void NMSSM_input::check_ewsb_output_parameters() const {
   bool supported = false;

   // check supported cases
   const bool Z3_symmetric = is_Z3_symmetric();
   if (SoftHiggsOut) {
      if (!is_set(mHd2) && !is_set(mHu2) && !is_set(mS2))
         supported = true;
      if (!is_set(lambdaS) || close(parameter[lambdaS], 0., EPSTOL))
        throw "# ERROR: <S> is zero!  Since SoftHiggsOut == true, <S> is not"
          " determined by the EWSB conditions, so <S> has to be set to"
          " a non-zero value on the user-side!\n";
   } else {
      if (Z3_symmetric) {
         if (!is_set(lambdaS) && !is_set(kappa) && !is_set(mS2))
            supported = true;
      } else {
         if (!is_set(mu) && !is_set(BmuOverCosBetaSinBeta) && !is_set(xiS))
            supported = true;
         if (!is_set(lambdaS) || close(parameter[lambdaS], 0., EPSTOL))
            throw "# ERROR: <S> is zero!  In the Z3 violating NMSSM <S> is not"
               " determined by the EWSB conditions, so <S> has to be set to"
               " a non-zero value on the user-side!\n";
      }
   }

   if (!supported) {
      std::ostringstream msg;
      msg << "# Error: no combination of the following unset parameters is"
         " currently supported as EWSB output for a Z3 "
          << (Z3_symmetric ? "symmetric" : "violating") << " NMSSM: ";
      for (unsigned i = 0; i < NUMBER_OF_NMSSM_INPUT_PARAMETERS; i++) {
         if (!is_set(static_cast<NMSSM_parameters>(i)))
            msg << parameter_names[i] << ", ";
      }
      msg << "\n" "# Note: supported are: ";
      if (SoftHiggsOut) {
            msg << "{" << parameter_names[mHd2]
                << ", " << parameter_names[mHu2]
                << ", " << parameter_names[mS2] << "}";
      } else {
         if (Z3_symmetric) {
            msg << "{" << parameter_names[lambdaS]
                << ", " << parameter_names[kappa]
                << ", " << parameter_names[mS2] << "}";
         } else {
            msg << "{" << parameter_names[mu]
                << ", " << parameter_names[BmuOverCosBetaSinBeta]
                << ", " << parameter_names[xiS] << "}";
         }
      }
      msg << '\n';
      throw msg.str();
   }
}

std::ostream& operator<<(std::ostream& lhs, const NMSSM_input& rhs) {
   for (unsigned i = 0; i < NMSSM_input::NUMBER_OF_NMSSM_INPUT_PARAMETERS; i++) {
      if (rhs.is_set(static_cast<NMSSM_input::NMSSM_parameters>(i))) {
         lhs << rhs.parameter_names[i] << " = "
             << rhs.get(static_cast<NMSSM_input::NMSSM_parameters>(i)) << ", ";
      }
   }
   return lhs;
}


NMSSM_command_line_parser::NMSSM_command_line_parser(NMSSM_input* nmssm_input_)
   : nmssm_input(nmssm_input_)
   , model_ident("")
   , m0(0.)
   , m12(0.)
   , a0(0.)
{
   assert(nmssm_input_);
}

void NMSSM_command_line_parser::parse(int argc, char* argv[]) {
   if (argc < 3)
      throw "# Error: NMSSM_command_line_parser: not enough command line"
         " arguments given";

   if (strcmp(argv[1], "nmssm"))
      throw "# Error: NMSSM_command_line_parser: first argument is not"
         " nmssm";

   model_ident = argv[2];

   for (int i = 3; i < argc; i++) {
      if (starts_with(argv[i], "--m0="))
         m0  = get_value(argv[i], "--m0=");
      else if (starts_with(argv[i], "--m12="))
         m12 = get_value(argv[i], "--m12=");
      else if (starts_with(argv[i], "--a0="))
         a0  = get_value(argv[i], "--a0=");
      else if (strcmp(argv[i], "--lambdaAtMsusy") == 0)
         softsusy::GUTlambda = false;
      else if (starts_with(argv[i], "--tanBeta="))
         nmssm_input->set(NMSSM_input::tanBeta, get_value(argv[i], "--tanBeta="));
      else if (starts_with(argv[i], "--mHu2="))
         nmssm_input->set(NMSSM_input::mHu2, get_value(argv[i], "--mHu2="));
      else if (starts_with(argv[i], "--mHd2="))
         nmssm_input->set(NMSSM_input::mHd2, get_value(argv[i], "--mHd2="));
      else if (starts_with(argv[i], "--mu="))
         nmssm_input->set(NMSSM_input::mu, get_value(argv[i], "--mu="));
      else if (starts_with(argv[i], "--m3SqrOverCosBetaSinBeta="))
         nmssm_input->set(NMSSM_input::BmuOverCosBetaSinBeta, get_value(argv[i], "--m3SqrOverCosBetaSinBeta="));
      else if (starts_with(argv[i], "--BmuSqrOverCosBetaSinBeta=")) {
         nmssm_input->set(NMSSM_input::BmuOverCosBetaSinBeta, get_value(argv[i], "--BmuSqrOverCosBetaSinBeta="));
         cout << "# Warning: --BmuSqrOverCosBetaSinBeta= is deprecated, "
            "please use --m3SqrOverCosBetaSinBeta= instead.\n";
      }
      else if (starts_with(argv[i], "--lambda="))
         nmssm_input->set(NMSSM_input::lambda, get_value(argv[i], "--lambda="));
      else if (starts_with(argv[i], "--kappa="))
         nmssm_input->set(NMSSM_input::kappa, get_value(argv[i], "--kappa="));
      else if (starts_with(argv[i], "--Alambda="))
         nmssm_input->set(NMSSM_input::Alambda, get_value(argv[i], "--Alambda="));
      else if (starts_with(argv[i], "--Akappa="))
         nmssm_input->set(NMSSM_input::Akappa, get_value(argv[i], "--Akappa="));
      else if (starts_with(argv[i], "--lambdaS="))
         nmssm_input->set(NMSSM_input::lambdaS, get_value(argv[i], "--lambdaS="));
      else if (starts_with(argv[i], "--xiF="))
         nmssm_input->set(NMSSM_input::xiF, get_value(argv[i], "--xiF="));
      else if (starts_with(argv[i], "--xiS="))
         nmssm_input->set(NMSSM_input::xiS, get_value(argv[i], "--xiS="));
      else if (starts_with(argv[i], "--muPrime="))
         nmssm_input->set(NMSSM_input::muPrime, get_value(argv[i], "--muPrime="));
      else if (starts_with(argv[i], "--mPrimeS2="))
         nmssm_input->set(NMSSM_input::mPrimeS2, get_value(argv[i], "--mPrimeS2="));
      else if (starts_with(argv[i], "--mS2="))
         nmssm_input->set(NMSSM_input::mS2, get_value(argv[i], "--mS2="));
      else
         throw std::string("# Error: unknown NMSSM command line option: ")
            + argv[i] + '\n';
   }

   // check universality condition
   if (strcmp(model_ident, "sugra") == 0) {
      // relax sugra condition if one of the following parameters is
      // set
      if (nmssm_input->is_set(NMSSM_input::Alambda) ||
          nmssm_input->is_set(NMSSM_input::Akappa) ||
          nmssm_input->is_set(NMSSM_input::mHd2) ||
          nmssm_input->is_set(NMSSM_input::mHu2) ||
          nmssm_input->is_set(NMSSM_input::mu) ||
          nmssm_input->is_set(NMSSM_input::BmuOverCosBetaSinBeta) ||
          nmssm_input->is_set(NMSSM_input::xiS) ||
          nmssm_input->is_set(NMSSM_input::mPrimeS2) ||
          nmssm_input->is_set(NMSSM_input::mS2)) {
         model_ident = "nonUniversal";
      }
   }
}

const char* NMSSM_command_line_parser::get_modelIdent() const {
   return model_ident;
}

DoubleVector NMSSM_command_line_parser::get_pars() const {
   DoubleVector pars(3);

   if (strcmp(model_ident, "sugra") == 0) {
      pars(1) = m0;
      pars(2) = m12;
      pars(3) = a0;
   } else if (strcmp(model_ident, "nonUniversal") == 0) {
      pars.setEnd(56);
      for (int i = 1; i <= 3; i++) pars(i) = m12;
      for (int i = 11; i <= 13; i++) pars(i) = a0;
      pars(21) = m0*m0;
      pars(22) = m0*m0;
      for (int i = 31; i <= 36; i++) pars(i) = m0;
      for (int i = 41; i <= 49; i++) pars(i) = m0;
      pars(50) = a0; // Alambda
      pars(51) = a0; // Akappa
      pars(52) = 0.; // mS'^2
      pars(53) = m0*m0; // mS^2

      if (nmssm_input->is_set(NMSSM_input::Alambda))
         pars(50) = nmssm_input->get(NMSSM_input::Alambda);
      if (nmssm_input->is_set(NMSSM_input::Akappa))
         pars(51) = nmssm_input->get(NMSSM_input::Akappa);
      if (nmssm_input->is_set(NMSSM_input::mHd2))
         pars(21) = nmssm_input->get(NMSSM_input::mHd2);
      if (nmssm_input->is_set(NMSSM_input::mHu2))
         pars(22) = nmssm_input->get(NMSSM_input::mHu2);
      if (nmssm_input->is_set(NMSSM_input::mu)) {
         pars(23) = nmssm_input->get(NMSSM_input::mu);
         pars(54) = nmssm_input->get(NMSSM_input::mu);
      }
      if (nmssm_input->is_set(NMSSM_input::BmuOverCosBetaSinBeta)) {
         pars(24) = nmssm_input->get(NMSSM_input::BmuOverCosBetaSinBeta);
         pars(55) = nmssm_input->get(NMSSM_input::BmuOverCosBetaSinBeta);
      }
      if (nmssm_input->is_set(NMSSM_input::mPrimeS2) &&
          nmssm_input->is_set(NMSSM_input::muPrime)) {
         // setting pars(52) = B' = mS'^2 / mu'
         const double muPrime  = nmssm_input->get(NMSSM_input::muPrime);
         const double mPrimeS2 = nmssm_input->get(NMSSM_input::mPrimeS2);
         if (!close(muPrime, 0.0, EPSTOL))
            pars(52) = mPrimeS2 / muPrime;
      }
      if (nmssm_input->is_set(NMSSM_input::mS2))
         pars(53) = nmssm_input->get(NMSSM_input::mS2);
      if (nmssm_input->is_set(NMSSM_input::xiS))
         pars(56) = nmssm_input->get(NMSSM_input::xiS);
   } else {
      throw std::string("# Error: NMSSM boundary condition ") + model_ident
         + " currently not supported at the command line\n";
   }

   return pars;
}

double NMSSM_command_line_parser::get_value(const std::string& str,
                                            const std::string& prefix) {
   return atof(str.substr(prefix.size()).c_str());
}

bool NMSSM_command_line_parser::starts_with(const std::string& str,
                                            const std::string& prefix) {
   return !str.compare(0, prefix.size(), prefix);
}

/// User supplied routine. Inputs m at the unification scale, and uses
/// inputParameters vector to output m with high energy soft boundary
/// conditions.
void NmssmMsugraBcs(NmssmSoftsusy & m, const DoubleVector & inputParameters) {
  double m0 = inputParameters.display(1);
  double m12 = inputParameters.display(2);
  double a0 = inputParameters.display(3);

  /// Sets scalar soft masses equal to m0, fermion ones to m12 and sets the
  /// trilinear scalar coupling to be a0
  m.standardSugra(m0, m12, a0);
}

//PA: msugra bcs in the mssm limit of the general nmssm
void MssmMsugraBcs(NmssmSoftsusy & m, const DoubleVector & inputParameters) {
  double m0 = inputParameters.display(1);
  double m12 = inputParameters.display(2);
  double a0 = inputParameters.display(3);

  /// Sets scalar soft masses equal to m0, fermion ones to m12 and sets the
  /// trilinear scalar coupling to be a0
  m.standardsemiSugra(m0, m12, a0, 0.0, 1e-15);
  m.setMspSquared(1e6);
}

//PA: semi-msugra bcs for the nmssm
void SemiMsugraBcs(NmssmSoftsusy & m, const DoubleVector & inputParameters) {
  assert(inputParameters.size() == 5 &&
         "SemiMsugraBcs: input parameter vector is not of length 5");

  double m0 = inputParameters.display(1);
  double m12 = inputParameters.display(2);
  double a0 = inputParameters.display(3);
  double Al = inputParameters.display(4);
  double Ak = inputParameters.display(5);

  /// Sets scalar soft masses equal to m0, fermion ones to m12 and sets the
  /// trilinear scalar coupling to be a0
  m.standardsemiSugra(m0, m12, a0, Al, Ak);
}

/// NMSSM Msugra, without setting the soft Higgs masses
void NmssmSugraNoSoftHiggsMassBcs(NmssmSoftsusy & m, const DoubleVector & inputParameters) {
  assert(inputParameters.size() >= 6 &&
         "NmssmSugraNoSoftHiggsMassBcs: input parameter"
         " vector has size < 6");

  const double mHd2 = m.displayMh1Squared(),
    mHu2 = m.displayMh2Squared(),
    ms2 = m.displayMsSquared();

  NmssmMsugraBcs(m, inputParameters);
  m.setMh1Squared(mHd2);
  m.setMh2Squared(mHu2);
  m.setMsSquared(ms2);

  // If SoftHiggsOut == true, then mu, Bmu and xiS are not fixed by
  // EWSB.  In this case they must be set in the BCS.
  if (!softsusy::Z3 && softsusy::SoftHiggsOut) {
    m.setSusyMu(inputParameters(4));
    m.setM3Squared(inputParameters(5));
    m.setXiS(inputParameters(6));
  }
}

/// NMSSM Msugra, without setting the soft Higgs masses
/// and non-universal Alambda and Akappa
void NmssmSemiMsugraNoSoftHiggsMassBcs(NmssmSoftsusy& m, const DoubleVector& inputParameters)
{
  assert(inputParameters.size() == 8 &&
         "NmssmSemiMsugraNoSoftHiggsMassBcs: input parameter"
         " vector is not of length 8");

  // mSUGRA without setting mHd2, mHu2, ms2
  NmssmSugraNoSoftHiggsMassBcs(m, inputParameters);

  const double Alambda = inputParameters(7);
  const double Akappa  = inputParameters(8);

  m.setTrialambda(Alambda * m.displayLambda());
  m.setTriakappa(Akappa * m.displayKappa());
}

void generalNmssmBcs(NmssmSoftsusy & m, const DoubleVector & inputParameters) {
  NmssmSusy s; SoftParsNmssm r;
  s = m.displaySusy();
  r.set(inputParameters);

  if (Z3 == false) {
    double m3sq = m.displayM3Squared();
    double XiS = m.displayXiS();
    r.setM3Squared(m3sq);
    r.setXiS(XiS);
  } else {
    double mSsq = m.displayMsSquared();
    r.setMsSquared(mSsq);
  }

  m.setSoftPars(r);
  m.setSusy(s);
}


/// This one doesn't overwrite mh1sq or mh2sq at the high scale.  For cases where the soft scalar Higgs masses are set in EWSB.
void generalNmssmBcs2(NmssmSoftsusy & m, const DoubleVector & inputParameters) {
  NmssmSusy s; SoftParsNmssm r;
  double mh1sq = m.displayMh1Squared();
  double mh2sq = m.displayMh2Squared();
  double mSsq = m.displayMsSquared();
  double m3sq = m.displayM3Squared();
  s = m.displaySusy();
  r.set(inputParameters);
  r.setMh1Squared(mh1sq);
  r.setMh2Squared(mh2sq);
  r.setM3Squared(m3sq);
  r.setMsSquared(mSsq);
  m.setSoftPars(r);
  m.setSusy(s);
}

void extendedNMSugraBcs(NmssmSoftsusy & m, const DoubleVector & inputParameters) {
  for (int i=1; i<=3; i++)
     m.setGauginoMass(i, inputParameters.display(i));
  if (inputParameters.display(25) > 1. && m.displaySetTbAtMX())
    m.setTanb(inputParameters.display(25));
  m.setTrilinearElement(UA, 1, 1, m.displayYukawaElement(YU, 1, 1) *
			inputParameters.display(11));
  m.setTrilinearElement(UA, 2, 2, m.displayYukawaElement(YU, 2, 2) *
			inputParameters.display(11));
  m.setTrilinearElement(UA, 3, 3, m.displayYukawaElement(YU, 3, 3) *
			inputParameters.display(11));
  m.setTrilinearElement(DA, 1, 1, m.displayYukawaElement(YD, 1, 1) *
			inputParameters.display(12));
  m.setTrilinearElement(DA, 2, 2, m.displayYukawaElement(YD, 2, 2) *
			inputParameters.display(12));
  m.setTrilinearElement(DA, 3, 3, m.displayYukawaElement(YD, 3, 3) *
			inputParameters.display(12));
  m.setTrilinearElement(EA, 1, 1, m.displayYukawaElement(YE, 1, 1) *
			inputParameters.display(13));
  m.setTrilinearElement(EA, 2, 2, m.displayYukawaElement(YE, 2, 2) *
			inputParameters.display(13));
  m.setTrilinearElement(EA, 3, 3, m.displayYukawaElement(YE, 3, 3) *
			inputParameters.display(13));
  m.setSoftMassElement(mLl, 1, 1, signedSqr(inputParameters.display(31)));
  m.setSoftMassElement(mLl, 2, 2, signedSqr(inputParameters.display(32)));
  m.setSoftMassElement(mLl, 3, 3, signedSqr(inputParameters.display(33)));
  m.setSoftMassElement(mEr, 1, 1, signedSqr(inputParameters.display(34)));
  m.setSoftMassElement(mEr, 2, 2, signedSqr(inputParameters.display(35)));
  m.setSoftMassElement(mEr, 3, 3, signedSqr(inputParameters.display(36)));
  m.setSoftMassElement(mQl, 1, 1, signedSqr(inputParameters.display(41)));
  m.setSoftMassElement(mQl, 2, 2, signedSqr(inputParameters.display(42)));
  m.setSoftMassElement(mQl, 3, 3, signedSqr(inputParameters.display(43)));
  m.setSoftMassElement(mUr, 1, 1, signedSqr(inputParameters.display(44)));
  m.setSoftMassElement(mUr, 2, 2, signedSqr(inputParameters.display(45)));
  m.setSoftMassElement(mUr, 3, 3, signedSqr(inputParameters.display(46)));
  m.setSoftMassElement(mDr, 1, 1, signedSqr(inputParameters.display(47)));
  m.setSoftMassElement(mDr, 2, 2, signedSqr(inputParameters.display(48)));
  m.setSoftMassElement(mDr, 3, 3, signedSqr(inputParameters.display(49)));

  // If SoftHiggsOut == true, then mHu2, mHd2 and mS2 will be EWSB
  // output.  If Z3 == true, mS2 will be EWSB output.
  if (!softsusy::SoftHiggsOut) {
    m.setMh1Squared(inputParameters.display(21));
    m.setMh2Squared(inputParameters.display(22));
    if (!softsusy::Z3)
      m.setMsSquared(inputParameters.display(53));
  }

  m.setTrialambda(m.displayLambda() * inputParameters.display(50));
  m.setTriakappa(m.displayKappa() * inputParameters.display(51));

  if (!softsusy::Z3) {
    m.setMspSquared(inputParameters.display(52) * m.displayMupr());
  }

  // If SoftHiggsOut == true, then mu, Bmu and xiS are not fixed by
  // EWSB.  In this case they must be set in the BCS.
  if (!softsusy::Z3 && softsusy::SoftHiggsOut) {
    m.setSusyMu(inputParameters(54));
    m.setM3Squared(inputParameters(55));
    m.setXiS(inputParameters(56));
  }
}

//This won't work for Z3 == true until we allow mS not to be set EWSB consitions
void nuhmINM(NmssmSoftsusy & m, const DoubleVector & inputParameters) {
  double m0 = inputParameters.display(1);
  double m12 = inputParameters.display(2);
  double mH  = inputParameters.display(3);
  double A0 = inputParameters.display(4);

  double Al = inputParameters.display(5);
  double Ak = inputParameters.display(6);
  /// Sets scalar soft masses equal to m0, fermion ones to m12 and sets the
  /// trilinear scalar coupling to be A0
  ///  if (m0 < 0.0) m.flagTachyon(true); Deleted on request from A Pukhov
  m.standardsemiSugra(m0, m12, A0, Al, Ak);

  m.setMh1Squared(mH * mH); m.setMh2Squared(mH * mH);
  m.setMsSquared(mH * mH);
  m.setTrialambda(m.displayLambda() * Al);
  m.setTriakappa(m.displayKappa() * Ak);
  if(Z3 == false) {
     m.setMspSquared(inputParameters.display(52) * m.displayMupr());
  }
}

void nuhmIINM(NmssmSoftsusy & m, const DoubleVector & inputParameters) {
  double m0 = inputParameters.display(1);
  double m12 = inputParameters.display(2);
  double mH1  = inputParameters.display(3);
  double mH2  = inputParameters.display(4);
  double A0 = inputParameters.display(5);
  double mS = 0.0;
  if (Z3 == false)
     mS = inputParameters.display(6);
  double Al = inputParameters.display(5);
  double Ak = inputParameters.display(6);
  m.standardsemiSugra(m0, m12, A0, Al, Ak);

  m.setMh1Squared(mH1 * mH1); m.setMh2Squared(mH2 * mH2);
  if (Z3 == false)
     m.setMsSquared(mS * mS);
  m.setTrialambda(m.displayLambda() * Al);
  m.setTriakappa(m.displayKappa() * Ak);
  if (Z3 == false) {
    m.setMspSquared(inputParameters.display(52) * m.displayMupr());
  }
}

//PA: anomally mediation Bcs with an extra universal m0 contribution.
void amsbBcs(NmssmSoftsusy & m, const DoubleVector & inputParameters) {
  double m32 = inputParameters.display(1);
  double m0 = inputParameters.display(2);

  m.standardSugra(m0, 0., 0.);
  m.addAmsb(m32);
}

/// LCT: Difference between two NMSSM SOFTSUSY objects in and out: EWSB terms only
double sumTol(const Softsusy<SoftParsNmssm> & in, const Softsusy<SoftParsNmssm> & out, int numTries) {
  DoubleVector sT(37);
  sumTol(in.displayDrBarPars(), out.displayDrBarPars(), sT);
  /// The predicted value of MZ^2 is an absolute measure of how close to a
  /// true solution we are:
  // double tbPred = 0.;
  double predictedMzSq = in.displayPredMzSq();
  /// We allow an extra factor of 10 for the precision in the predicted value
  /// of MZ compared to TOLERANCE if the program is struggling and gone beyond
  /// 10 tries - an extra 2 comes from MZ v MZ^2
  int k = sT.size() - 2;
  if (!in.displayProblem().testSeriousProblem()) {
    sT(k) = 0.5 * sTfn(predictedMzSq, sqr(MZ));
    if (numTries > 10) sT(k) *= 0.1;
  }

  return sT.max();
}
  
  /// explicit template instantiations
  template class Softsusy<SoftParsNmssm>;
  template class SoftPars<NmssmSusy, nmsBrevity>;
  
} ///< namespace softsusy
