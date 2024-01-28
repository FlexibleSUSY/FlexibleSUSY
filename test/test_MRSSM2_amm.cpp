// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_MRSSM2_amm

#include <boost/test/unit_test.hpp>

#include "test_MRSSM2.hpp"

#include "lowe.h"
#include "wrappers.hpp"
#include "MRSSM2_amm.hpp"
#include "MRSSM2_lepton_amm_wrapper.hpp"
#include "cxx_qft/MRSSM2_qft.hpp"

using DiagonalMatrix3 = Eigen::DiagonalMatrix<double, 3>;

using namespace flexiblesusy;

struct AMM {
   double ae{};
   double amu{};
   double amu_wrapper{};
   double damu{};
   double atau{};
};

struct Data {
   AMM a_1L{};
   AMM a_1L_2LQCD{};
   AMM a_1L_2LQCD_2LBarrZee{};
};


Data calc_amm(const MRSSM2_slha& m, const softsusy::QedQcd& qedqcd, Spectrum_generator_settings settings)
{
   using MRSSM2_cxx_diagrams::fields::Fe;

   Data result;

   // 1L
   settings.set(Spectrum_generator_settings::calculate_amm, 1.0);

   result.a_1L.ae = MRSSM2_amm::calculate_amm<Fe>(m, qedqcd, settings, 0);
   result.a_1L.amu = MRSSM2_amm::calculate_amm<Fe>(m, qedqcd, settings, 1);
   result.a_1L.atau = MRSSM2_amm::calculate_amm<Fe>(m, qedqcd, settings, 2);

   // 1L + 2L QED
   settings.set(Spectrum_generator_settings::calculate_amm, 1.5);

   result.a_1L_2LQCD.ae = MRSSM2_amm::calculate_amm<Fe>(m, qedqcd, settings, 0);
   result.a_1L_2LQCD.amu = MRSSM2_amm::calculate_amm<Fe>(m, qedqcd, settings, 1);
   result.a_1L_2LQCD.amu_wrapper = MRSSM2_lepton_amm_wrapper::calculate_Fe_amm(m, qedqcd, settings, 1);
   result.a_1L_2LQCD.damu = MRSSM2_amm::calculate_amm_uncertainty<Fe>(m, qedqcd, settings, 1);
   result.a_1L_2LQCD.atau = MRSSM2_amm::calculate_amm<Fe>(m, qedqcd, settings, 2);

   // 1L + 2L QED + Barr-Zee
   settings.set(Spectrum_generator_settings::calculate_amm, 2.0);

   result.a_1L_2LQCD_2LBarrZee.ae = MRSSM2_amm::calculate_amm<Fe>(m, qedqcd, settings, 0);
   result.a_1L_2LQCD_2LBarrZee.amu = MRSSM2_amm::calculate_amm<Fe>(m, qedqcd, settings, 1);
   result.a_1L_2LQCD_2LBarrZee.atau = MRSSM2_amm::calculate_amm<Fe>(m, qedqcd, settings, 2);

   return result;
}


BOOST_AUTO_TEST_CASE( test_amu_chargino_dominance )
{
   softsusy::QedQcd qedqcd;
   Spectrum_generator_settings settings;
   MRSSM2_input_parameters input;

   input.TanBeta = 10;
   input.Ms = 1000;
   input.LamTDInput = -1.0;
   input.LamTUInput = -1.0;
   input.LamSDInput = 1.1;
   input.LamSUInput = -1.1;
   input.MuDInput = 400;
   input.MuUInput = 400;
   input.BMuInput = Sqr(300);
   input.mq2Input = DiagonalMatrix3(Sqr(1000), Sqr(1000), Sqr(1000));
   input.ml2Input = DiagonalMatrix3(Sqr( 500), Sqr( 500), Sqr( 500));
   input.md2Input = DiagonalMatrix3(Sqr(1000), Sqr(1000), Sqr(1000));
   input.mu2Input = DiagonalMatrix3(Sqr(1000), Sqr(1000), Sqr(1000));
   input.me2Input = DiagonalMatrix3(Sqr( 500), Sqr( 500), Sqr( 500));
   input.mS2Input = Sqr(2000);
   input.mT2Input = Sqr(3000);
   input.moc2Input = Sqr(1000);
   input.mRd2Input = Sqr(700);
   input.mRu2Input = Sqr(1000);
   input.MDBSInput = 1000;
   input.MDWBTInput = 500;
   input.MDGocInput = 1500;

   const auto m = setup_MRSSM2(input, qedqcd, settings);
   const auto data = calc_amm(m, qedqcd, settings);

   BOOST_CHECK_CLOSE_FRACTION(data.a_1L.ae, -2.082466293438991e-15, 1e-7);
   BOOST_CHECK_CLOSE_FRACTION(data.a_1L.amu, -8.8989415316380609e-11, 1e-7);
   BOOST_CHECK_CLOSE_FRACTION(data.a_1L.atau, -2.3322905904784368e-08, 1e-7);

   BOOST_CHECK_CLOSE_FRACTION(data.a_1L_2LQCD.ae, -1.8019087470353613e-15, 1e-7);
   BOOST_CHECK_CLOSE_FRACTION(data.a_1L_2LQCD.amu, -8.1717853946513499e-11, 1e-7);
   BOOST_CHECK_EQUAL(data.a_1L_2LQCD.amu, data.a_1L_2LQCD.amu_wrapper);
   BOOST_CHECK_CLOSE_FRACTION(data.a_1L_2LQCD.damu, 3.0887729897668767e-11, 1e-7);
   BOOST_CHECK_CLOSE_FRACTION(data.a_1L_2LQCD.atau, -2.2071645085574141e-08, 1e-7);

   BOOST_CHECK_CLOSE_FRACTION(data.a_1L_2LQCD_2LBarrZee.ae, -1.1468646899659171e-15, 1e-7);
   BOOST_CHECK_CLOSE_FRACTION(data.a_1L_2LQCD_2LBarrZee.amu, -5.371267977475807e-11, 1e-7);
   BOOST_CHECK_CLOSE_FRACTION(data.a_1L_2LQCD_2LBarrZee.atau, -1.4150896272117654e-08, 1e-7);
}


BOOST_AUTO_TEST_CASE( test_amu_neutralino_dominance )
{
   softsusy::QedQcd qedqcd;
   Spectrum_generator_settings settings;
   MRSSM2_input_parameters input;

   input.TanBeta = 10;
   input.Ms = 1000;
   input.LamTDInput = -1.0;
   input.LamTUInput = -1.0;
   input.LamSDInput = 1.1;
   input.LamSUInput = -1.1;
   input.MuDInput = 400;
   input.MuUInput = 400;
   input.BMuInput = Sqr(300);
   input.mq2Input = DiagonalMatrix3(Sqr(1000), Sqr(1000), Sqr(1000));
   input.ml2Input = DiagonalMatrix3(Sqr(8000), Sqr(8000), Sqr(8000));
   input.md2Input = DiagonalMatrix3(Sqr(1000), Sqr(1000), Sqr(1000));
   input.mu2Input = DiagonalMatrix3(Sqr(1000), Sqr(1000), Sqr(1000));
   input.me2Input = DiagonalMatrix3(Sqr( 500), Sqr( 500), Sqr( 500));
   input.mS2Input = Sqr(2000);
   input.mT2Input = Sqr(3000);
   input.moc2Input = Sqr(1000);
   input.mRd2Input = Sqr(700);
   input.mRu2Input = Sqr(1000);
   input.MDBSInput = 1000;
   input.MDWBTInput = 500;
   input.MDGocInput = 1500;

   const auto m = setup_MRSSM2(input, qedqcd, settings);
   const auto data = calc_amm(m, qedqcd, settings);

   BOOST_CHECK_CLOSE_FRACTION(data.a_1L.ae, 1.5869774606726533e-16, 1e-7);
   BOOST_CHECK_CLOSE_FRACTION(data.a_1L.amu, 6.830270770956226e-12, 1e-7);
   BOOST_CHECK_CLOSE_FRACTION(data.a_1L.atau, 3.9423844338780157e-09, 1e-7);

   BOOST_CHECK_CLOSE_FRACTION(data.a_1L_2LQCD.ae, 1.3740221631863563e-16, 1e-7);
   BOOST_CHECK_CLOSE_FRACTION(data.a_1L_2LQCD.amu, 6.2743670477551711e-12, 1e-7);
   BOOST_CHECK_EQUAL(data.a_1L_2LQCD.amu, data.a_1L_2LQCD.amu_wrapper);
   BOOST_CHECK_CLOSE_FRACTION(data.a_1L_2LQCD.atau, 3.7317183255799172e-09, 1e-7);
   BOOST_CHECK_CLOSE_FRACTION(data.a_1L_2LQCD.damu, 3.042784947881566e-11, 1e-7);

   BOOST_CHECK_CLOSE_FRACTION(data.a_1L_2LQCD_2LBarrZee.ae, 7.9652674799007061e-16, 1e-7);
   BOOST_CHECK_CLOSE_FRACTION(data.a_1L_2LQCD_2LBarrZee.amu, 3.4453972095834096e-11, 1e-7);
   BOOST_CHECK_CLOSE_FRACTION(data.a_1L_2LQCD_2LBarrZee.atau, 1.170004365157329e-08, 1e-7);
}
