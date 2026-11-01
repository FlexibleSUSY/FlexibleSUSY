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
#include "MRSSM2_mass_eigenstates_running.hpp"

using DiagonalMatrix3 = Eigen::DiagonalMatrix<double, 3>;

using namespace flexiblesusy;

struct AMM {
   double ae{}, dae{};
   double amu{}, damu{};
   double atau{}, datau{};
};

struct Dataset {
   AMM a_1L{};                 ///< 1L
   AMM a_1L_2LQCD{};           ///< 1L + 2L QED
   AMM a_1L_2LQCD_2LBarrZee{}; ///< 1L + 2L QED + 2L Barr-Zee
};


AMM calc_amm(const MRSSM2_slha& m, const softsusy::QedQcd& qedqcd, Spectrum_generator_settings settings)
{
   using MRSSM2_cxx_diagrams::fields::Fe;

   AMM result{};

   result.ae = MRSSM2_amm::calculate_amm<Fe>(m, qedqcd, settings, 0);
   result.dae = MRSSM2_amm::calculate_amm_uncertainty<Fe>(m, qedqcd, settings, 0);
   result.amu = MRSSM2_amm::calculate_amm<Fe>(m, qedqcd, settings, 1);
   result.damu = MRSSM2_amm::calculate_amm_uncertainty<Fe>(m, qedqcd, settings, 1);
   result.atau = MRSSM2_amm::calculate_amm<Fe>(m, qedqcd, settings, 2);
   result.datau = MRSSM2_amm::calculate_amm_uncertainty<Fe>(m, qedqcd, settings, 2);

   BOOST_CHECK_EQUAL(result.amu, MRSSM2_lepton_amm_wrapper::calculate_Fe_amm(m, qedqcd, settings, 1));

   return result;
}


Dataset calc_amm_sets(const MRSSM2_slha& m, const softsusy::QedQcd& qedqcd, Spectrum_generator_settings settings)
{
   Dataset result;

   // 1L
   settings.set(Spectrum_generator_settings::calculate_amm, 1.0);
   result.a_1L = calc_amm(m, qedqcd, settings);

   // 1L + 2L QED
   settings.set(Spectrum_generator_settings::calculate_amm, 1.5);
   result.a_1L_2LQCD = calc_amm(m, qedqcd, settings);

   // 1L + 2L QED + Barr-Zee
   settings.set(Spectrum_generator_settings::calculate_amm, 2.0);
   result.a_1L_2LQCD_2LBarrZee = calc_amm(m, qedqcd, settings);

   return result;
}


// check consistency of uncertainties at different loop levels
void check_uncertainty(const Dataset& data)
{
   BOOST_CHECK_GT(data.a_1L.damu, data.a_1L_2LQCD.damu);
   BOOST_CHECK_GT(data.a_1L_2LQCD.damu, data.a_1L_2LQCD_2LBarrZee.damu);

   const double oneL = data.a_1L.amu;
   const double twoL_QCD = data.a_1L_2LQCD.amu - data.a_1L.amu;
   const double twoL_BarrZee = data.a_1L_2LQCD_2LBarrZee.amu - data.a_1L_2LQCD.amu;

   BOOST_CHECK_GT(data.a_1L.damu, twoL_QCD + twoL_BarrZee);
   BOOST_CHECK_GT(data.a_1L_2LQCD.damu, twoL_BarrZee);
   BOOST_CHECK_GT(data.a_1L_2LQCD_2LBarrZee.damu, 0.0);
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

   const auto model = setup_MRSSM2(input, qedqcd, settings);
   const auto data = calc_amm_sets(model, qedqcd, settings);

   BOOST_CHECK_CLOSE_FRACTION(data.a_1L.ae, -2.0825031661019949e-15, 1e-7);
   BOOST_CHECK_CLOSE_FRACTION(data.a_1L.amu, -8.8990990810185124e-11, 1e-7);
   BOOST_CHECK_CLOSE_FRACTION(data.a_1L.atau, -2.3323312714524165e-08, 1e-7);

   BOOST_CHECK_CLOSE_FRACTION(data.a_1L_2LQCD.ae, -1.8019406457152465e-15, 1e-7);
   BOOST_CHECK_CLOSE_FRACTION(data.a_1L_2LQCD.amu, -8.1719300481563679e-11, 1e-7);
   BOOST_CHECK_CLOSE_FRACTION(data.a_1L_2LQCD.atau, -2.207203001934572e-08, 1e-7);

   BOOST_CHECK_CLOSE_FRACTION(data.a_1L_2LQCD_2LBarrZee.ae, -1.1468824936542875e-15, 1e-7);
   BOOST_CHECK_CLOSE_FRACTION(data.a_1L_2LQCD_2LBarrZee.amu, -5.3713523714967649e-11, 1e-7);
   BOOST_CHECK_CLOSE_FRACTION(data.a_1L_2LQCD_2LBarrZee.atau, -1.4151110782475472e-08, 1e-7);

   check_uncertainty(data);
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

   const auto model = setup_MRSSM2(input, qedqcd, settings);
   const auto data = calc_amm_sets(model, qedqcd, settings);

   BOOST_CHECK_CLOSE_FRACTION(data.a_1L.ae, 1.5869696382545575e-16, 1e-7);
   BOOST_CHECK_CLOSE_FRACTION(data.a_1L.amu, 6.8302375916443017e-12, 1e-7);
   BOOST_CHECK_CLOSE_FRACTION(data.a_1L.atau, 3.942387224198737e-09, 1e-7);

   BOOST_CHECK_CLOSE_FRACTION(data.a_1L_2LQCD.ae, 1.3740153891433592e-16, 1e-7);
   BOOST_CHECK_CLOSE_FRACTION(data.a_1L_2LQCD.amu, 6.274336563654111e-12, 1e-7);
   BOOST_CHECK_CLOSE_FRACTION(data.a_1L_2LQCD.atau, 3.7317209639323539e-09, 1e-7);

   BOOST_CHECK_CLOSE_FRACTION(data.a_1L_2LQCD_2LBarrZee.ae, 7.9653070906802407e-16, 1e-7);
   BOOST_CHECK_CLOSE_FRACTION(data.a_1L_2LQCD_2LBarrZee.amu, 3.4454139896377721e-11, 1e-7);
   BOOST_CHECK_CLOSE_FRACTION(data.a_1L_2LQCD_2LBarrZee.atau, 1.1700102338909448e-08, 1e-7);

   check_uncertainty(data);
}
