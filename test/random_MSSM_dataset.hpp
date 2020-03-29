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

#ifndef TEST_RANDOM_MSSM_DATASET_H
#define TEST_RANDOM_MSSM_DATASET_H

#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

#include "MSSM_two_scale_spectrum_generator.hpp"
#include "MSSM_input_parameters.hpp"
#include "MSSM_mass_eigenstates.hpp"

#include "random_sign_dataset.hpp"

#include "lowe.h"

namespace flexiblesusy
{
namespace boost_test_tools
{
namespace detail
{
static auto random_MSSM_dataset( void ) -> decltype(
	boost::unit_test::data::random( 5.0, 30.0 ) ^
	boost_test_tools::random_sign_dataset{} ^
	boost::unit_test::data::random( 1.94e+16, 1.95e+16 ) ^
	boost::unit_test::data::random( 500000.0, 10000000.0 ) ^
	boost::unit_test::data::random( 500000.0, 10000000.0 ) ^
	boost::unit_test::data::random( 700.0, 3000.0 ) ^
	boost::unit_test::data::random( 700.0, 3000.0 ) ^
	boost::unit_test::data::random( 700.0, 3000.0 )
)
{
	auto TanBeta_random = boost::unit_test::data::random( 5.0, 30.0 );
	auto signMu_random = boost_test_tools::random_sign_dataset{};
	auto Q_random = boost::unit_test::data::random( 1.94e+16, 1.95e+16 );
	auto mHd2IN_random = boost::unit_test::data::random( 500000.0, 10000000.0 );
	auto mHu2IN_random = boost::unit_test::data::random( 500000.0, 10000000.0 );
	auto MassBInput_random = boost::unit_test::data::random( 700.0, 3000.0 );
	auto MassWBInput_random = boost::unit_test::data::random( 700.0, 3000.0 );
	auto MassGInput_random = boost::unit_test::data::random( 700.0, 3000.0 );
	
	return (
		std::move(TanBeta_random) ^ std::move(signMu_random) ^
		std::move(Q_random) ^ std::move(mHd2IN_random) ^
		std::move(mHu2IN_random) ^ std::move(MassBInput_random) ^
		std::move(MassWBInput_random) ^ std::move(MassGInput_random)
	);
}
}

static auto random_MSSM_dataset( int number_of_random_samples ) ->
decltype(
	boost::unit_test::data::xrange( 1, number_of_random_samples + 1 ) ^
	detail::random_MSSM_dataset()
)
{
	auto random_datset =
		boost::unit_test::data::xrange( 1, number_of_random_samples + 1 ) ^
		detail::random_MSSM_dataset();
	
	return std::move(random_datset);
}

#define FS_TEST_MSSM_PARAMETER_SEQUENCE TanBeta, SignMu, Qin, mHd2IN, \
	mHu2IN, MassBInput, MassWBInput, MassGInput

static MSSM_input_parameters wrap_input_parameters(
	double TanBeta, int SignMu, double Qin, double mHd2IN,
	double mHu2IN, double MassBInput, double MassWBInput,
	double MassGInput )
{
	MSSM_input_parameters input;
	
	input.TanBeta = TanBeta;
	input.SignMu = SignMu;
	input.Qin = Qin;
	input.mHd2IN = mHd2IN;
	input.mHu2IN = mHu2IN;
	input.MassBInput = MassBInput;
	input.MassWBInput = MassWBInput;
	input.MassGInput = MassGInput;
	
	// Pick a sensible value for sfermion masses
	double sfermion_massSqr = (input.mHd2IN + input.mHu2IN) / 2.0;
	Eigen::Matrix<double,3,3> sfermion_mass_matrix;
	
	sfermion_mass_matrix <<
		sfermion_massSqr, 0, 0,
		0, sfermion_massSqr, 0,
		0, 0, sfermion_massSqr;
	
	input.mq2Input = sfermion_mass_matrix;
	input.ml2Input = sfermion_mass_matrix;
	input.md2Input = sfermion_mass_matrix;
	input.mu2Input = sfermion_mass_matrix;
	input.me2Input = sfermion_mass_matrix;
	
	return input;
}

static MSSM_mass_eigenstates calculate_spectrum(
	const MSSM_input_parameters &input )
{
	softsusy::QedQcd qedqcd;
	
	Spectrum_generator_settings spectrum_generator_settings;
	spectrum_generator_settings.set(
		Spectrum_generator_settings::calculate_sm_masses, 1.0 );
	
	MSSM_spectrum_generator<Two_scale> spectrum_generator;
	spectrum_generator.set_settings( spectrum_generator_settings );
	spectrum_generator.set_parameter_output_scale( qedqcd.get_scale() );

	spectrum_generator.run( qedqcd, input );
	return spectrum_generator.get_model();
}
}
}

#endif
