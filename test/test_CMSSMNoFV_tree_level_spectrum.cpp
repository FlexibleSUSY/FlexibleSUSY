
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_CMSSMNoFV_tree_level_spectrum

#include <boost/test/unit_test.hpp>

#include "test.hpp"
#include "test_CMSSMNoFV.hpp"
#include "CMSSM_two_scale_model.hpp"
#include "CMSSMNoFV_two_scale_model.hpp"

#define COMPARE_MASS(p) TEST_EQUALITY(a.get_##p(), b.get_##p());

using namespace flexiblesusy;

void test_spectrum_equality(const CMSSMNoFV<Two_scale>& a, const CMSSM<Two_scale>& b)
{
   COMPARE_MASS(Mhh);
   COMPARE_MASS(MAh);
   COMPARE_MASS(MHpm);

   COMPARE_MASS(MChi);
   COMPARE_MASS(MCha);
   COMPARE_MASS(MGlu);

   COMPARE_MASS(MVP);
   COMPARE_MASS(MVZ);
   COMPARE_MASS(MVWm);
   COMPARE_MASS(MVG);

   // up-type squarks
   const Eigen::Array<double,2,1> MSu(a.get_MSu());
   const Eigen::Array<double,2,1> MSc(a.get_MSc());
   const Eigen::Array<double,2,1> MSt(a.get_MSt());
   const Eigen::Array<double,6,1> MSu_full(b.get_MSu());
   Eigen::Array<double,6,1> MSu_mixed;
   MSu_mixed(0) = MSu(0);
   MSu_mixed(1) = MSu(1);
   MSu_mixed(2) = MSc(0);
   MSu_mixed(3) = MSc(1);
   MSu_mixed(4) = MSt(0);
   MSu_mixed(5) = MSt(1);
   std::sort(MSu_mixed.data(), MSu_mixed.data() + MSu_mixed.size());

   TEST_EQUALITY(MSu_full, MSu_mixed);

   // down-type squarks
   const Eigen::Array<double,2,1> MSd(a.get_MSd());
   const Eigen::Array<double,2,1> MSs(a.get_MSs());
   const Eigen::Array<double,2,1> MSb(a.get_MSb());
   const Eigen::Array<double,6,1> MSd_full(b.get_MSd());
   Eigen::Array<double,6,1> MSd_mixed;
   MSd_mixed(0) = MSd(0);
   MSd_mixed(1) = MSd(1);
   MSd_mixed(2) = MSs(0);
   MSd_mixed(3) = MSs(1);
   MSd_mixed(4) = MSb(0);
   MSd_mixed(5) = MSb(1);
   std::sort(MSd_mixed.data(), MSd_mixed.data() + MSd_mixed.size());

   TEST_EQUALITY(MSd_full, MSd_mixed);

   // up-type sleptons
   const double MSveL(a.get_MSveL());
   const double MSvmL(a.get_MSvmL());
   const double MSvtL(a.get_MSvtL());
   const Eigen::Array<double,3,1> MSv_full(b.get_MSv());
   Eigen::Array<double,3,1> MSv_mixed;
   MSv_mixed(0) = MSveL;
   MSv_mixed(1) = MSvmL;
   MSv_mixed(2) = MSvtL;
   std::sort(MSv_mixed.data(), MSv_mixed.data() + MSv_mixed.size());

   TEST_EQUALITY(MSv_full, MSv_mixed);

   // down-type sleptons
   const Eigen::Array<double,2,1> MSe(a.get_MSe());
   const Eigen::Array<double,2,1> MSm(a.get_MSm());
   const Eigen::Array<double,2,1> MStau(a.get_MStau());
   const Eigen::Array<double,6,1> MSe_full(b.get_MSe());
   Eigen::Array<double,6,1> MSe_mixed;
   MSe_mixed(0) = MSe(0);
   MSe_mixed(1) = MSe(1);
   MSe_mixed(2) = MSm(0);
   MSe_mixed(3) = MSm(1);
   MSe_mixed(4) = MStau(0);
   MSe_mixed(5) = MStau(1);
   std::sort(MSe_mixed.data(), MSe_mixed.data() + MSe_mixed.size());

   TEST_EQUALITY(MSe_full, MSe_mixed);

   // up-type quarks
   const double MFu(a.get_MFu());
   const double MFc(a.get_MFc());
   const double MFt(a.get_MFt());
   const Eigen::Array<double,3,1> MFu_full(b.get_MFu());
   Eigen::Array<double,3,1> MFu_mixed;
   MFu_mixed(0) = MFu;
   MFu_mixed(1) = MFc;
   MFu_mixed(2) = MFt;
   std::sort(MFu_mixed.data(), MFu_mixed.data() + MFu_mixed.size());

   TEST_EQUALITY(MFu_full, MFu_mixed);

   // down-type quarks
   const double MFd(a.get_MFd());
   const double MFs(a.get_MFs());
   const double MFb(a.get_MFb());
   const Eigen::Array<double,3,1> MFd_full(b.get_MFd());
   Eigen::Array<double,3,1> MFd_mixed;
   MFd_mixed(0) = MFd;
   MFd_mixed(1) = MFs;
   MFd_mixed(2) = MFb;
   std::sort(MFd_mixed.data(), MFd_mixed.data() + MFd_mixed.size());

   TEST_EQUALITY(MFd_full, MFd_mixed);

   // up-type leptons
   const double MFve(a.get_MFve());
   const double MFvm(a.get_MFvm());
   const double MFvt(a.get_MFvt());
   const Eigen::Array<double,3,1> MFv_full(b.get_MFv());
   Eigen::Array<double,3,1> MFv_mixed;
   MFv_mixed(0) = MFve;
   MFv_mixed(1) = MFvm;
   MFv_mixed(2) = MFvt;
   std::sort(MFv_mixed.data(), MFv_mixed.data() + MFv_mixed.size());

   TEST_EQUALITY(MFv_full, MFv_mixed);

   // down-type leptons
   const double MFe(a.get_MFe());
   const double MFm(a.get_MFm());
   const double MFtau(a.get_MFtau());
   const Eigen::Array<double,3,1> MFe_full(b.get_MFe());
   Eigen::Array<double,3,1> MFe_mixed;
   MFe_mixed(0) = MFe;
   MFe_mixed(1) = MFm;
   MFe_mixed(2) = MFtau;
   std::sort(MFe_mixed.data(), MFe_mixed.data() + MFe_mixed.size());

   TEST_EQUALITY(MFe_full, MFe_mixed);
}

BOOST_AUTO_TEST_CASE( test_CMSSMNoFV_tree_level_spectrum )
{
   CMSSMNoFV_input_parameters input;
   input.TanBeta = 10.;
   input.m0 = 125.;
   input.m12 = 200.;
   input.SignMu = 1;
   input.Azero = 0.;

   CMSSMNoFV<Two_scale> m1;
   CMSSM<Two_scale> m2;
   setup_CMSSM_models(m1, m2, input);

   test_spectrum_equality(m1, m2);
   BOOST_REQUIRE(gErrors == 0);
   if (gErrors) {
      BOOST_FAIL("spectra are not equal");
      gErrors = 0;
   }
}
