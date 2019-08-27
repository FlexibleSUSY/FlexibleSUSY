#include <complex>
#include "loop_library_interface.hpp"

class Collier : public Loop_library_interface {
   public:
      std::complex<double> C00(
         std::complex<double> p10, std::complex<double> p21, std::complex<double> p20,
         std::complex<double> m02, std::complex<double> m12, std::complex<double> m22,
         double scl2) noexcept;
};
