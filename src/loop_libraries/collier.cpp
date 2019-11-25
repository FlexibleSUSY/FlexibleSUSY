#include "collier.hpp"
#include <limits>

#include "ukC0.hpp"

#define two_point_impl(NAME)\
   std::complex<double> NAME##_impl(\
      const std::complex<double>*,\
      const std::complex<double>*, const std::complex<double>*,\
      const double*);
#define three_point_impl(NAME)\
   std::complex<double> NAME##_impl(\
      const std::complex<double>*, const std::complex<double>*, const std::complex<double>*,\
      const std::complex<double>*, const std::complex<double>*, const std::complex<double>*,\
      const double*);

/* Non-vanishing imaginary parts of momentum invariants are not yet
 * suppoted by the current version (1.2.4) of COLLIER. */
#define two_point_collier(NAME)\
   std::complex<double> Collier::NAME(\
      std::complex<double> p10_in,\
      std::complex<double> m02_in, std::complex<double> m12_in,\
      double scl2_in) noexcept {\
   \
   const std::complex<double> p10 (p10_in.real(), 0.);\
   const std::complex<double> m02 = m02_in;\
   const std::complex<double> m12 = m12_in;\
   double scl2 = scl2_in;\
   \
   return NAME##_impl(&p10, &m02, &m12, &scl2);\
}
#define three_point_collier(NAME)\
   std::complex<double> Collier::NAME(\
      std::complex<double> p10_in, std::complex<double> p21_in, std::complex<double> p20_in,\
      std::complex<double> m02_in, std::complex<double> m12_in, std::complex<double> m22_in,\
      double scl2_in) noexcept {\
   \
   const std::complex<double> p10 (p10_in.real(), 0.);\
   const std::complex<double> p21 (p21_in.real(), 0.);\
   const std::complex<double> p20 (p20_in.real(), 0.);\
   const std::complex<double> m02 = m02_in;\
   const std::complex<double> m12 = m12_in;\
   const std::complex<double> m22 = m22_in;\
   double scl2 = scl2_in;\
   \
   return NAME##_impl(&p10, &p21, &p20, &m02, &m12, &m22, &scl2);\
}

// Fortran wrapper routines
extern "C" {
   two_point_impl(B0)
   two_point_impl(B1)

   three_point_impl(C0)
   three_point_impl(C1)
   three_point_impl(C2)
   three_point_impl(C00)
   three_point_impl(C11)
   three_point_impl(C12)
   three_point_impl(C22)
}

//namespace collier {
   // C++ wrapper with non-pointer parameter
   two_point_collier(B0)
   two_point_collier(B1)

std::complex<double> Collier::C0(
      std::complex<double> p10_in, std::complex<double> p21_in, std::complex<double> p20_in,
      std::complex<double> m02_in, std::complex<double> m12_in, std::complex<double> m22_in,
      double scl2_in) noexcept {

   const std::complex<double> p10 (p10_in.real(), 0.);
   const std::complex<double> p21 (p21_in.real(), 0.);
   const std::complex<double> p20 (p20_in.real(), 0.);
   const std::complex<double> m02 = m02_in;
   const std::complex<double> m12 = m12_in;
   const std::complex<double> m22 = m22_in;
   double scl2 = scl2_in;

   if( std::min({m02_in.real(),m12_in.real(),m22_in.real()}) / std::max({abs(p10),abs(p21),abs(p20)}) > 10000.) {
      return ukC0(p10_in.real(),p21_in.real(),p20_in.real(),m02_in.real(),m12_in.real(),m22_in.real());
   }

   return C0_impl(&p10, &p21, &p20, &m02, &m12, &m22, &scl2);
}

   three_point_collier(C1)
   three_point_collier(C2)
   three_point_collier(C00)
   three_point_collier(C11)
   three_point_collier(C12)
   three_point_collier(C22)
//}
