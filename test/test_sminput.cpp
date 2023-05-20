
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_sminput

#include <boost/test/unit_test.hpp>
#include "slhaea.h"

const std::string input_str =
"# Example input in SLHA2 format\n"
"Block SMINPUTS                 # Standard Model inputs\n"
"    1   1.279340000e+02        # alpha^(-1) SM MSbar(MZ)\n"
"    2   1.166370000e-05        # G_Fermi\n"
"    3   1.176000000e-01        # alpha_s(MZ) SM MSbar\n"
"    4   9.118760000e+01        # MZ(pole)\n"
"    5   4.200000000e+00        # mb(mb) SM MSbar\n"
"    6   1.733000000e+02        # mtop(pole)\n"
"    7   1.777000000e+00        # mtau(pole)\n"
"    8   0.000000000e+00        # mnu3(pole)\n"
"   11   5.109989020e-04        # melectron(pole)\n"
"   12   0.000000000e+00        # mnu1(pole)\n"
"   13   1.056583570e-01        # mmuon(pole)\n"
"   14   0.000000000e+00        # mnu2(pole)\n"
"   21   4.750000000e-03        # md(2 GeV)\n"
"   22   2.400000000e-03        # mu(2 GeV)\n"
"   23   1.040000000e-01        # ms(2 GeV)\n"
"   24   1.270000000e+00        # mc(2 GeV)\n";

BOOST_AUTO_TEST_CASE( test_reading_sminput )
{
   const auto input = SLHAea::Coll::from_str(input_str);

   BOOST_REQUIRE(!input.empty());

   const double alpha_em = SLHAea::to<double>(input.at("SMINPUTS").at("1").at(1));
   const double G_mu     = SLHAea::to<double>(input.at("SMINPUTS").at("2").at(1));
   const double alpha_s  = SLHAea::to<double>(input.at("SMINPUTS").at("3").at(1));
   const double m_Z      = SLHAea::to<double>(input.at("SMINPUTS").at("4").at(1));
   const double m_b      = SLHAea::to<double>(input.at("SMINPUTS").at("5").at(1));
   const double m_t      = SLHAea::to<double>(input.at("SMINPUTS").at("6").at(1));
   const double m_tau    = SLHAea::to<double>(input.at("SMINPUTS").at("7").at(1));
   const double m_nu3    = SLHAea::to<double>(input.at("SMINPUTS").at("8").at(1));
   const double m_e      = SLHAea::to<double>(input.at("SMINPUTS").at("11").at(1));
   const double m_nu1    = SLHAea::to<double>(input.at("SMINPUTS").at("12").at(1));
   const double m_mu     = SLHAea::to<double>(input.at("SMINPUTS").at("13").at(1));
   const double m_nu2    = SLHAea::to<double>(input.at("SMINPUTS").at("14").at(1));
   const double m_d      = SLHAea::to<double>(input.at("SMINPUTS").at("21").at(1));
   const double m_u      = SLHAea::to<double>(input.at("SMINPUTS").at("22").at(1));
   const double m_s      = SLHAea::to<double>(input.at("SMINPUTS").at("23").at(1));
   const double m_c      = SLHAea::to<double>(input.at("SMINPUTS").at("24").at(1));

   BOOST_CHECK_EQUAL(alpha_em, 1.279340000e+02);
   BOOST_CHECK_EQUAL(G_mu    , 1.166370000e-05);
   BOOST_CHECK_EQUAL(alpha_s , 1.176000000e-01);
   BOOST_CHECK_EQUAL(m_Z     , 9.118760000e+01);
   BOOST_CHECK_EQUAL(m_b     , 4.200000000e+00);
   BOOST_CHECK_EQUAL(m_t     , 1.733000000e+02);
   BOOST_CHECK_EQUAL(m_tau   , 1.777000000e+00);
   BOOST_CHECK_EQUAL(m_nu3   , 0.000000000e+00);
   BOOST_CHECK_EQUAL(m_e     , 5.109989020e-04);
   BOOST_CHECK_EQUAL(m_nu1   , 0.000000000e+00);
   BOOST_CHECK_EQUAL(m_mu    , 1.056583570e-01);
   BOOST_CHECK_EQUAL(m_nu2   , 0.000000000e+00);
   BOOST_CHECK_EQUAL(m_d     , 4.750000000e-03);
   BOOST_CHECK_EQUAL(m_u     , 2.400000000e-03);
   BOOST_CHECK_EQUAL(m_s     , 1.040000000e-01);
   BOOST_CHECK_EQUAL(m_c     , 1.270000000e+00);
}

BOOST_AUTO_TEST_CASE( test_iteration_sminput )
{
   const auto input = SLHAea::Coll::from_str(input_str);

   BOOST_REQUIRE(!input.empty());

   double alpha_em = 0.0;
   double G_mu     = 0.0;
   double alpha_s  = 0.0;
   double m_Z      = 0.0;
   double m_b      = 0.0;
   double m_t      = 0.0;
   double m_tau    = 0.0;
   double m_nu3    = 0.0;
   double m_e      = 0.0;
   double m_nu1    = 0.0;
   double m_mu     = 0.0;
   double m_nu2    = 0.0;
   double m_d      = 0.0;
   double m_u      = 0.0;
   double m_s      = 0.0;
   double m_c      = 0.0;

   for (SLHAea::Block::const_iterator line = input.at("SMINPUTS").cbegin();
        line != input.at("SMINPUTS").cend(); ++line) {
      if (!line->is_data_line()) continue;

      if (line->size() >= 2) {
         const int key = SLHAea::to<int>((*line)[0]);
         switch (key) {
         case 1:  alpha_em = SLHAea::to<double>((*line)[1]); break;
         case 2:  G_mu     = SLHAea::to<double>((*line)[1]); break;
         case 3:  alpha_s  = SLHAea::to<double>((*line)[1]); break;
         case 4:  m_Z      = SLHAea::to<double>((*line)[1]); break;
         case 5:  m_b      = SLHAea::to<double>((*line)[1]); break;
         case 6:  m_t      = SLHAea::to<double>((*line)[1]); break;
         case 7:  m_tau    = SLHAea::to<double>((*line)[1]); break;
         case 8:  m_nu3    = SLHAea::to<double>((*line)[1]); break;
         case 11: m_e      = SLHAea::to<double>((*line)[1]); break;
         case 12: m_nu1    = SLHAea::to<double>((*line)[1]); break;
         case 13: m_mu     = SLHAea::to<double>((*line)[1]); break;
         case 14: m_nu2    = SLHAea::to<double>((*line)[1]); break;
         case 21: m_d      = SLHAea::to<double>((*line)[1]); break;
         case 22: m_u      = SLHAea::to<double>((*line)[1]); break;
         case 23: m_s      = SLHAea::to<double>((*line)[1]); break;
         case 24: m_c      = SLHAea::to<double>((*line)[1]); break;
         default: BOOST_CHECK(false); break;
         }
      } else {
         BOOST_CHECK(false);
      }
   }

   BOOST_CHECK_EQUAL(alpha_em, 1.279340000e+02);
   BOOST_CHECK_EQUAL(G_mu    , 1.166370000e-05);
   BOOST_CHECK_EQUAL(alpha_s , 1.176000000e-01);
   BOOST_CHECK_EQUAL(m_Z     , 9.118760000e+01);
   BOOST_CHECK_EQUAL(m_b     , 4.200000000e+00);
   BOOST_CHECK_EQUAL(m_t     , 1.733000000e+02);
   BOOST_CHECK_EQUAL(m_tau   , 1.777000000e+00);
   BOOST_CHECK_EQUAL(m_nu3   , 0.000000000e+00);
   BOOST_CHECK_EQUAL(m_e     , 5.109989020e-04);
   BOOST_CHECK_EQUAL(m_nu1   , 0.000000000e+00);
   BOOST_CHECK_EQUAL(m_mu    , 1.056583570e-01);
   BOOST_CHECK_EQUAL(m_nu2   , 0.000000000e+00);
   BOOST_CHECK_EQUAL(m_d     , 4.750000000e-03);
   BOOST_CHECK_EQUAL(m_u     , 2.400000000e-03);
   BOOST_CHECK_EQUAL(m_s     , 1.040000000e-01);
   BOOST_CHECK_EQUAL(m_c     , 1.270000000e+00);
}
