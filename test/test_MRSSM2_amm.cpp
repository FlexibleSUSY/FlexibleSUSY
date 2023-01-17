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
#include "MRSSM2_lepton_gm2_wrapper.hpp"
#include "cxx_qft/MRSSM2_qft.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_amu )
{
   typedef Eigen::DiagonalMatrix<double, 3> DiagonalMatrix3;
   MRSSM2_input_parameters input;

   // chargino dominance
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

   softsusy::QedQcd qedqcd;
   Spectrum_generator_settings settings;

   MRSSM2_slha m = setup_MRSSM2(input, qedqcd, settings);

   using MRSSM2_cxx_diagrams::fields::Fe;

   // 1L
   settings.set(Spectrum_generator_settings::calculate_amm, 1.0);

   auto ae = MRSSM2_amm::calculate_amm<Fe>(m, qedqcd, settings, 0);
   BOOST_CHECK_CLOSE_FRACTION(ae, -2.0824982517457706e-15, 1e-7);

   auto amu = MRSSM2_amm::calculate_amm<Fe>(m, qedqcd, settings, 1);
   BOOST_CHECK_CLOSE_FRACTION(amu, -8.8990780706146349e-11, 1e-7);

   auto atau = MRSSM2_amm::calculate_amm<Fe>(m, qedqcd, settings, 2);
   BOOST_CHECK_CLOSE_FRACTION(atau, -2.3323253290458326e-08, 1e-7);

   // 1L + 2L QED
   settings.set(Spectrum_generator_settings::calculate_amm, 1.5);

   ae = MRSSM2_amm::calculate_amm<Fe>(m, qedqcd, settings, 0);
   BOOST_CHECK_CLOSE_FRACTION(ae, -1.8019363934392491e-15, 1e-7);

   amu = MRSSM2_amm::calculate_amm<Fe>(m, qedqcd, settings, 1);
   BOOST_CHECK_CLOSE_FRACTION(amu, -8.1719107571588341e-11, 1e-7);
   BOOST_CHECK_CLOSE_FRACTION(amu, MRSSM2_lepton_gm2_wrapper::calculate_Fe_gm2(m, qedqcd, settings, 1), 1e-16);
   double damu = MRSSM2_amm::calculate_amm_uncertainty<Fe>(m, qedqcd, settings, 1);
   BOOST_CHECK_CLOSE_FRACTION(damu, 9.0675554049399791e-13, 1e-7);

   atau = MRSSM2_amm::calculate_amm<Fe>(m, qedqcd, settings, 2);
   BOOST_CHECK_CLOSE_FRACTION(atau, -2.2071973783347808e-08, 1e-7);

   // 1L + 2L QED + Barr-Zee
   settings.set(Spectrum_generator_settings::calculate_amm, 2.0);

   ae = MRSSM2_amm::calculate_amm<Fe>(m, qedqcd, settings, 0);
   BOOST_CHECK_CLOSE_FRACTION(ae, -1.1528572931380434e-15, 1e-7);

   amu = MRSSM2_amm::calculate_amm<Fe>(m, qedqcd, settings, 1);
   BOOST_CHECK_CLOSE_FRACTION(amu, -5.3968953868935609e-11, 1e-7);

   atau = MRSSM2_amm::calculate_amm<Fe>(m, qedqcd, settings, 2);
   BOOST_CHECK_CLOSE_FRACTION(atau, -1.4223352839313197e-08, 1e-7);

   // neutralino dominance

   // 1L
   settings.set(Spectrum_generator_settings::calculate_amm, 1.0);

   ae = MRSSM2_amm::calculate_amm<Fe>(m, qedqcd, settings, 0);
   BOOST_CHECK_CLOSE_FRACTION(ae, -2.0824982517457706e-15, 1e-7);

   amu = MRSSM2_amm::calculate_amm<Fe>(m, qedqcd, settings, 1);
   BOOST_CHECK_CLOSE_FRACTION(amu, -8.8990780706146349e-11, 1e-7);

   atau = MRSSM2_amm::calculate_amm<Fe>(m, qedqcd, settings, 2);
   BOOST_CHECK_CLOSE_FRACTION(atau, -2.3323253290458326e-08, 1e-7);

   // 1L + 2L QED
   settings.set(Spectrum_generator_settings::calculate_amm, 1.5);

   input.ml2Input = DiagonalMatrix3(Sqr(8000), Sqr(8000), Sqr(8000));
   m = setup_MRSSM2(input, qedqcd, settings);

   ae = MRSSM2_amm::calculate_amm<Fe>(m, qedqcd, settings, 0);
   BOOST_CHECK_CLOSE_FRACTION(ae, 1.3730107649079216e-16, 1e-7);

   amu = MRSSM2_amm::calculate_amm<Fe>(m, qedqcd, settings, 1);
   BOOST_CHECK_CLOSE_FRACTION(amu, 6.2697795623223974e-12, 1e-7);

   atau = MRSSM2_amm::calculate_amm<Fe>(m, qedqcd, settings, 2);
   BOOST_CHECK_CLOSE_FRACTION(atau, 3.7303931716651099e-09, 1e-7);

   damu = MRSSM2_amm::calculate_amm_uncertainty<Fe>(m, qedqcd, settings, 1);
   BOOST_CHECK_CLOSE_FRACTION(damu, 7.277946979761437e-12, 1e-7);

   // 1L + 2L QED + Barr-Zee
   settings.set(Spectrum_generator_settings::calculate_amm, 2.0);

   ae = MRSSM2_amm::calculate_amm<Fe>(m, qedqcd, settings, 0);
   BOOST_CHECK_CLOSE_FRACTION(ae, 7.9102513947685461e-16, 1e-7);

   amu = MRSSM2_amm::calculate_amm<Fe>(m, qedqcd, settings, 1);
   BOOST_CHECK_CLOSE_FRACTION(amu, 3.4218497896247393e-11, 1e-7);

   atau = MRSSM2_amm::calculate_amm<Fe>(m, qedqcd, settings, 2);
   BOOST_CHECK_CLOSE_FRACTION(atau, 1.1633430821565949e-08, 1e-7);
}
