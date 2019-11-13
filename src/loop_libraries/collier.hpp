#include <complex>
#include "loop_library_interface.hpp"

#define two_point(NAME)\
   std::complex<double> NAME(\
      std::complex<double> p10,\
      std::complex<double> m02, std::complex<double> m12,\
      double scl2) noexcept;
#define three_point(NAME)\
   std::complex<double> NAME(\
      std::complex<double> p10, std::complex<double> p21, std::complex<double> p20,\
      std::complex<double> m02, std::complex<double> m12, std::complex<double> m22,\
      double scl2) noexcept;

class Collier : public Loop_library_interface {
   public:
      two_point(B0)
      two_point(B1)

      three_point(C0)
      three_point(C1)
      three_point(C2)
      three_point(C00)
      three_point(C11)
      three_point(C12)
      three_point(C22)
};
