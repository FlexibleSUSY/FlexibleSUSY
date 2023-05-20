
#ifndef MSSM_UTILS_H
#define MSSM_UTILS_H

#include <iosfwd>
#include <string>

namespace softsusy {

class DoubleVector;

struct sBrevity;
class MssmSusy;

template<class SoftPars> class Softsusy;
template<class Susy, class Brevity> class SoftPars;

typedef SoftPars<MssmSusy, sBrevity> SoftParsMssm;
typedef Softsusy<SoftParsMssm> MssmSoftsusy;

/// Allows user to specify a boundary condition where ALL SUSY breaking
/// parameters are specified in inputParameters
void generalBcs(MssmSoftsusy & m, const DoubleVector & inputParameters);
/// Sets all soft parameters in m except for mh1sq or mh2sq: it is intended
/// for the case where mu and M_A^0(pole) is specified
void generalBcs2(MssmSoftsusy & m, const DoubleVector & inputParameters);
/// non-universal mSUGRA boundary conditions including mH1^2 and mH2^2
void extendedSugraBcs(MssmSoftsusy & m, const DoubleVector & inputParameters);
/// User supplied routine. Inputs m at the unification scale, and uses
/// inputParameters vector to output m with high energy soft boundary
/// conditions.
void sugraBcs(MssmSoftsusy & m, const DoubleVector & inputParameters);
/// Non-universal higgs mass conditions. Paramaters are, in order: m0,m12,mH,a0
void nuhmI(MssmSoftsusy & m, const DoubleVector & inputParameters);
/// Non-universal higgs mass conditions. Paramaters are, in order:
/// m0,m12,mH1,mH2,a0
void nuhmII(MssmSoftsusy & m, const DoubleVector & inputParameters);
/// Adds 2-loop AMSB boundary conditions onto m
void amsbBcs(MssmSoftsusy & m, const DoubleVector & inputParameters);
/// Large Volume string compactification boundary conditions
void lvsBcs(MssmSoftsusy & m, const DoubleVector & inputParameters);
/// One-loop GMSB boundary conditions
void gmsbBcs(MssmSoftsusy & m, const DoubleVector & inputParameters);
/// For the user to define....
void userDefinedBcs(MssmSoftsusy & m, const DoubleVector & inputParameters);
void nonUniGauginos(MssmSoftsusy & m, const DoubleVector & inputParameters);
void splitGmsb(MssmSoftsusy & m, const DoubleVector & inputParameters);

/// Pietro's fit to the totality of LEP2 data: gives whether the lightest
/// Higgs is allowed or not - depends only upon the mass and beta-alpha.
/// true means allowed, false means ruled out.
bool testLEPHiggs(const MssmSoftsusy & r, double error = 3.0);

/// Given mu paramer, tau Yukawa htau, family number examined, finds height of
/// potential for temp at |H_2|=h2.
double ufb3fn(double mu, double htau, double h2, int family, const MssmSoftsusy
	      & temp);

/// For UFB-3direction, returns scale at which one-loop corrections are
/// smallest. IO parameters: inminTol is fractional accuracy with which
/// minimum is found, eR is value of RH selectron field, h2 is value of H2
/// field, Lisq=|L_i|^2 slepton VEV value, mx=high BC-scale
double getQhat(double inminTol,double eR, double h2, double Lisq, double mx,
	       MssmSoftsusy & temp);

/// Calculates fractional difference in Drbar masses between in and out
double sumTol(const MssmSoftsusy & in, const MssmSoftsusy & out, int numTries);

/// Prints out the identity of the LSP to standard output.
/// temp, j are defined from lsp function
std::string recogLsp(int temp, int j);

/// Formatted output
std::ostream & operator<<(std::ostream&, const MssmSoftsusy&);

} // namespace softsusy

#endif
