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

#define four_point(NAME)\
   std::complex<double> NAME(\
      std::complex<double> p10, std::complex<double> p21, std::complex<double> p32,\
      std::complex<double> p30, std::complex<double> p20, std::complex<double> p31,\
      std::complex<double> m02, std::complex<double> m12, std::complex<double> m22, std::complex<double> m32,\
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

      four_point(D0)
      four_point(D00)
      four_point(D1)
      four_point(D11)
      four_point(D12)
      four_point(D13)
      four_point(D2)
      four_point(D22)
      four_point(D23)
      four_point(D3)
      four_point(D33)
};
