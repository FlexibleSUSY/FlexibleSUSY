#include <complex>
#include "loop_library_interface.hpp"

class Loop_tools : public Loop_library_interface {
   public:
      std::complex<double> C00(
         double p10, double p21, double p20,
         double m02, double m12, double m22,
         double scl2) noexcept;
};

