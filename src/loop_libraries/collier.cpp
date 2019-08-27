#include "collier.hpp"

extern "C" {
   // Fortran wrapper routine
   std::complex<double> C00_cll(
         const std::complex<double>*, const std::complex<double>*, const std::complex<double>*,
         const std::complex<double>*, const std::complex<double>*, const std::complex<double>*,
         const double*);
}

//namespace collier {
   // C++ wrapper with non-pointer parameter
   std::complex<double> Collier::C00(
         std::complex<double> p10_in, std::complex<double> p21_in, std::complex<double> p20_in,
         std::complex<double> m02_in, std::complex<double> m12_in, std::complex<double> m22_in,
         double scl2_in) noexcept {

      const std::complex<double> p10 = p10_in;
      const std::complex<double> p21 = p21_in;
      const std::complex<double> p20 = p20_in;
      const std::complex<double> m02 = m02_in;
      const std::complex<double> m12 = m12_in;
      const std::complex<double> m22 = m22_in;
      double scl2 = scl2_in;

      return C00_cll(&p10, &p21, &p20, &m02, &m12, &m22, &scl2);
   }
//}
