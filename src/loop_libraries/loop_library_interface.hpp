#include <complex>
#define two_point_virtual(NAME)\
   virtual std::complex<double> NAME(\
      std::complex<double> p10,\
      std::complex<double> m02, std::complex<double> m12,\
      double scl2) = 0;
#define three_point_virtual(NAME)\
   virtual std::complex<double> NAME(\
      std::complex<double> p10, std::complex<double> p21, std::complex<double> p20,\
      std::complex<double> m02, std::complex<double> m12, std::complex<double> m22,\
      double scl2) = 0;

#define four_point_virtual(NAME)\
   virtual std::complex<double> NAME(\
      std::complex<double> p10, std::complex<double> p21, std::complex<double> p32,\
      std::complex<double> p30, std::complex<double> p20, std::complex<double> p31,\
      std::complex<double> m02, std::complex<double> m12, std::complex<double> m22, std::complex<double> m32,\
      double scl2) = 0;

enum class Loop_library {
   LoopTools,
   Collier
};

class Loop_library_interface {
   public:
      /*
      static std::unique_ptr<LoopLibraryInterface> get_library(LoopLibrary type) {
         switch (type) {
            case LoopLibrary::Collier:
               return std::make_unique<Collier>();
            case LoopLibrary::LoopTools:
               return std::make_unique<LoopTools>();
         }
      }
      */
      two_point_virtual(B0)
      two_point_virtual(B1)

      three_point_virtual(C0)
      three_point_virtual(C1)
      three_point_virtual(C2)
      three_point_virtual(C00)
      three_point_virtual(C11)
      three_point_virtual(C12)
      three_point_virtual(C22)

      four_point_virtual(D0)
      four_point_virtual(D00)
      four_point_virtual(D1)
      four_point_virtual(D11)
      four_point_virtual(D12)
      four_point_virtual(D13)
      four_point_virtual(D2)
      four_point_virtual(D22)
      four_point_virtual(D23)
      four_point_virtual(D3)
      four_point_virtual(D33)
};
