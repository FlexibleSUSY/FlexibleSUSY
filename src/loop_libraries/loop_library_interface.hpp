#include <complex>

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
      virtual std::complex<double> B0(
         std::complex<double> p10,
         std::complex<double> m02, std::complex<double> m12,
         double scl2) = 0;
      virtual std::complex<double> C0(
         std::complex<double> p10, std::complex<double> p21, std::complex<double> m20,
         std::complex<double> m02, std::complex<double> m12, std::complex<double> m22,
         double scl2) = 0;
      virtual std::complex<double> C00(
         std::complex<double> p10, std::complex<double> p21, std::complex<double> m20,
         std::complex<double> m02, std::complex<double> m12, std::complex<double> m22,
         double scl2) = 0;
};
