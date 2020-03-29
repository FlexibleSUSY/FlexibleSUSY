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

#ifndef TEST_RANDOM_SM_DATASET_H
#define TEST_RANDOM_SM_DATASET_H

#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

#include "SM_two_scale_spectrum_generator.hpp"
#include "SM_input_parameters.hpp"
#include "SM_mass_eigenstates.hpp"

#include "lowe.h"

namespace flexiblesusy
{
namespace boost_test_tools
{
namespace detail
{
static auto random_SM_dataset( void ) -> decltype(
	boost::unit_test::data::random( 0.15, 0.3 ) ^
	boost::unit_test::data::random( 200.0, 3000.0 ) ^
	boost::unit_test::data::random( 100.0, 1000.0 )
)
{
	auto LambdaIN_random = boost::unit_test::data::random( 0.15, 0.3 );
	auto Qin_random = boost::unit_test::data::random( 200.0, 3000.0 );
	auto QEWSB_random = boost::unit_test::data::random( 100.0, 1000.0 );
	
	return (
		std::move(LambdaIN_random) ^ std::move(Qin_random) ^
		std::move(QEWSB_random)
	);
}
}

static auto random_SM_dataset( int number_of_random_samples ) ->
decltype(
	boost::unit_test::data::xrange( 1, number_of_random_samples + 1 ) ^
	detail::random_SM_dataset()
)
{
	if( !(number_of_random_samples >= 2) )
		throw std::invalid_argument( "random_SM::data_points(): \
number_of_random_samples must be >= 2." );
	
	auto random_datset =
		boost::unit_test::data::xrange( 1, number_of_random_samples + 1 ) ^
		detail::random_SM_dataset();
	
	return std::move(random_datset);
}

#define FS_TEST_SM_PARAMETER_SEQUENCE LambdaIN_p, Qin_p, QEWSB_p

static SM_input_parameters wrap_SM_parameters(
	double LambdaIN, int Qin, double QEWSB )
{
	SM_input_parameters input;
	
	input.LambdaIN = LambdaIN;
	input.Qin = Qin;
	input.QEWSB = QEWSB;
	
	return input;
}

static SM_mass_eigenstates calculate_spectrum(
	const SM_input_parameters &input )
{
	softsusy::QedQcd qedqcd;
	
	Spectrum_generator_settings spectrum_generator_settings;
	spectrum_generator_settings.set(
		Spectrum_generator_settings::calculate_sm_masses, 1.0 );
	
	SM_spectrum_generator<Two_scale> spectrum_generator;
	spectrum_generator.set_settings( spectrum_generator_settings );
	spectrum_generator.set_parameter_output_scale( qedqcd.get_scale() );

	spectrum_generator.run( qedqcd, input );
	return spectrum_generator.get_model();
}
}
}

#endif
