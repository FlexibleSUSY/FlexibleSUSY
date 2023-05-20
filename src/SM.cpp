#include "SM.hpp"

namespace flexiblesusy {

using namespace std;


CM33 standard_unitary_matrix(Real s12, Real s13, Real s23, Real delta)
{
    Real c12 = sqrt(1 - sqr(s12));
    Real c23 = sqrt(1 - sqr(s23));
    Real c13 = sqrt(1 - sqr(s13));
    Comp P = polar(1.0, delta);
    CM33 V;
    V <<  c12*c13              ,  s12*c13              , s13*conj(P),
	 -s12*c23-c12*s23*s13*P,  c12*c23-s12*s23*s13*P, s23*c13    ,
	  s12*s23-c12*c23*s13*P, -c12*s23-s12*c23*s13*P, c23*c13    ;
    return V;
}

CM33 standard_VCKM(Real gamma, Real Vub)
{
    const Real Vus = .2257;		// +- .0021 , PDG 2007
    const Real Vud = .97377;		// +- .00027, PDG 2007
    const Real Vcb = 41.6e-3;		// +- .6e-3 , PDG 2007

    const Real s12 = Vus / sqrt(sqr(Vud)+sqr(Vus)); // = lambda
    const Real s23 = s12 * Vcb / Vus;
    const Real s13 = Vub;
    const Real delta = gamma;	// up to 0.000774555 radians of error

    return standard_unitary_matrix(s12, s13, s23, delta);
}

CM33 standard_VCKM(Real gamma)
{
    return standard_VCKM(gamma, modVub_cent);
}

} // namespace flexiblesusy
