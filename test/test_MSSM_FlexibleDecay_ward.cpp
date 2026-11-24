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

/**
 * Ward identity of the generated 1-loop amplitudes h -> gamma gamma and
 * h -> gamma Z.  For an external photon gauge invariance requires
 *
 *    form_factor_g = -(p1.p2) form_factor_21 ,
 *
 * exactly if all masses (internal, and of the external Higgs and Z) are
 * the tree-level masses after tree-level EWSB.  In Feynman gauge the W
 * loops use triple- and four-vector, momentum-difference, ghost and
 * scalar-vector-vector vertices, and the sfermion and charged Higgs
 * loops scalar-scalar-vector-vector vertices, so the test checks their
 * translation to the CXXDiagrams vertex methods (meta/Decays.m,
 * ConvertCouplingToCPP).
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_MSSM_FlexibleDecay_ward

#include <boost/test/unit_test.hpp>

#include "MSSM_input_parameters.hpp"
#include "MSSM_mass_eigenstates.hpp"
#include "MSSM_slha_io.hpp"
#include "MSSM_two_scale_spectrum_generator.hpp"
#include "cxx_qft/MSSM_context_base.hpp"
#include "cxx_qft/MSSM_fields.hpp"
#include "decays/MSSM_decays.hpp"

#include "decays/flexibledecay_settings.hpp"
#include "loop_libraries/loop_library.hpp"
#include "lowe.h"
#include "physical_input.hpp"
#include "spectrum_generator_settings.hpp"

#include <Eigen/Core>

#include <complex>

using namespace flexiblesusy;
using namespace flexiblesusy::MSSM_cxx_diagrams;

namespace
{

/// |form_factor_g/(-(p1.p2) form_factor_21) - 1|
template <typename Amplitude>
double ward_violation(const Amplitude& amp, double p1p2)
{
   return std::abs(amp.form_factor_g / (-p1p2 * amp.form_factor_21) - 1.);
}

} // anonymous namespace

BOOST_AUTO_TEST_CASE(test_ward_identity_S_to_gamma_gamma_and_gamma_Z)
{
   Loop_library::set(-1);

   MSSM_slha_io slha_io;
   slha_io.read_from_file("model_files/MSSM/LesHouches.in.MSSM");

   softsusy::QedQcd qedqcd;
   Physical_input physical_input;
   Spectrum_generator_settings settings;
   MSSM_input_parameters input;
   slha_io.fill(settings);
   slha_io.fill(qedqcd);
   slha_io.fill(physical_input);
   slha_io.fill(input);

   // heavy SUSY spectrum, light Higgs boson not too close to MZ
   const double m0sq = 1e6;
   input.mHd2IN = m0sq;
   input.mHu2IN = m0sq;
   for (auto* m2 : {&input.mq2Input, &input.ml2Input, &input.md2Input,
                    &input.mu2Input, &input.me2Input}) {
      *m2 = m0sq * Eigen::Matrix<double, 3, 3>::Identity();
   }
   input.MassBInput = 2000.;
   input.MassWBInput = 2000.;
   input.MassGInput = 2000.;
   input.Auij(2, 2) = -4000.;

   MSSM_spectrum_generator<Two_scale> spectrum_generator;
   spectrum_generator.set_settings(settings);
   spectrum_generator.run(qedqcd, input);
   BOOST_REQUIRE(!spectrum_generator.get_problems().have_problem());

   // tree-level EWSB and tree-level masses everywhere
   MSSM_mass_eigenstates model = spectrum_generator.get_model();
   model.solve_ewsb_tree_level();
   model.calculate_DRbar_masses();
   auto physical = model.get_physical();
   physical.Mhh = model.get_Mhh();
   model.set_physical(physical);

   // the external Z of S -> gamma Z has the pole mass of qedqcd
   auto qedqcd_tree = qedqcd;
   qedqcd_tree.setPoleMZ(model.get_MVZ());

   const MSSM_decays decays(model, qedqcd_tree, physical_input,
                            FlexibleDecay_settings());
   const context_base context{&model};

   const double mz = model.get_MVZ();

   for (int i = 0; i < 2; i++) {
      const double mh = model.get_Mhh(i);
      BOOST_TEST_MESSAGE("hh(" << i << "): m = " << mh << ", MZ = " << mz);
      BOOST_REQUIRE_GT(std::abs(mh - mz), 1.);

      const auto aa =
         decays.calculate_amplitude<fields::hh, fields::VP, fields::VP>(
            context, {i}, {}, {});
      const auto az =
         decays.calculate_amplitude<fields::hh, fields::VP, fields::VZ>(
            context, {i}, {}, {});

      BOOST_REQUIRE_GT(std::abs(aa.form_factor_21), 0.);
      BOOST_REQUIRE_GT(std::abs(az.form_factor_21), 0.);

      BOOST_CHECK_SMALL(ward_violation(aa, 0.5 * mh * mh), 1e-8);
      BOOST_CHECK_SMALL(ward_violation(az, 0.5 * (mh * mh - mz * mz)), 1e-6);
   }
}
