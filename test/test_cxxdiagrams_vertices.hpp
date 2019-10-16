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

#ifndef TEST_VERTICES_H
#define TEST_VERTICES_H

#include <string>
#include <iterator>

#include "cxx_qft/SM_qft.hpp"

#include "test_complex_equality.hpp"
#include "SM_mass_eigenstates.hpp"

namespace flexiblesusy
{
namespace boost_test_tools
{
namespace detail
{
template<class Sequence>
std::string sequence_to_string( const Sequence &s )
{
	auto it = std::begin( s );
	auto end = std::end( s );
	
	if( it == end )
		return "{}";
	
	std::stringstream stream{ "{" };
	--end;
	
	while( it != end )
	{
		stream << *it << ", ";
		++it;
	}
	
	stream << *it << "}";
	return stream.str();
}
}

struct test_zero_vertex
{
	SM_cxx_diagrams::context_base context;
	
	test_zero_vertex( const SM_mass_eigenstates &model )
	: context( model ) {}
	
	template<class Vertex> void operator()( Vertex )
	{
		BOOST_TEST_CONTEXT( "Vertex: " <<
			boost::core::demangle( typeid(Vertex).name() ) )
		{
		for( auto &index : index_range<Vertex>() )
		{
			BOOST_TEST_CONTEXT( "indices: " <<
				detail::sequence_to_string( index ) )
			{
      auto vertex = Vertex::evaluate( index, context );
      this->test( vertex );
			}
		}
		}
	}
	
	void test( const SM_cxx_diagrams::ScalarVertex &v )
	{ TEST_COMPLEX_EQUALITY( v.value(), 0 ); }
	
	void test( const SM_cxx_diagrams::MomentumVertex &v )
	{ TEST_COMPLEX_EQUALITY( v.value( v.index() ), 0 ); }
	
	void test( const SM_cxx_diagrams::MomentumDifferenceVertex &v )
	{
		TEST_COMPLEX_EQUALITY(
			v.value( v.incoming_index(), v.outgoing_index() ), 0 );
	}
	
	void test( const SM_cxx_diagrams::InverseMetricVertex &v )
	{ TEST_COMPLEX_EQUALITY( v.value(), 0 ); }
	
	void test( const SM_cxx_diagrams::ChiralVertex &v )
	{
		TEST_COMPLEX_EQUALITY( v.left(), 0 );
		TEST_COMPLEX_EQUALITY( v.right(), 0 );
	}
	
	void test( const SM_cxx_diagrams::TripleVectorVertex &v )
	{ TEST_COMPLEX_EQUALITY( v.value(
			SM_cxx_diagrams::TripleVectorVertex::even_permutation{} ), 0 ); }
	
	void test( const SM_cxx_diagrams::QuadrupleVectorVertex &v )
	{
		TEST_COMPLEX_EQUALITY( v.value1(), 0 );
		TEST_COMPLEX_EQUALITY( v.value2(), 0 );
		TEST_COMPLEX_EQUALITY( v.value3(), 0 );
	}
};

struct test_vertex_equality
{
	SM_cxx_diagrams::context_base context;
	
	test_vertex_equality( const SM_mass_eigenstates &model )
	: context( model ) {}
	
	template<class VertexPair> void operator()( VertexPair )
	{
		BOOST_TEST_CONTEXT( "VertexPair: " <<
			boost::core::demangle( typeid(VertexPair).name() ) )
		{
		using first = typename boost::mpl::first<VertexPair>::type;
		using second = typename boost::mpl::second<VertexPair>::type;
		
		static_assert( std::is_same<
			decltype( index_range<first>() ),
			decltype( index_range<second>() )
		>::value, "Incompatible index types in vertex pair" );
		
		for( auto &index : index_range<first>() )
		{
			BOOST_TEST_CONTEXT( "indices: " <<
				detail::sequence_to_string( index ) )
			{
      auto first_vertex = first::evaluate( index, context );
      auto second_vertex = second::evaluate( index, context );
      
      this->test( first_vertex, second_vertex );
			}
		}
		}
	}
	
	void test( const SM_cxx_diagrams::ScalarVertex &v1,
		const SM_cxx_diagrams::ScalarVertex &v2 )
	{ TEST_COMPLEX_EQUALITY( v1.value(), v2.value() ); }
	
	void test( const SM_cxx_diagrams::MomentumVertex &v1,
		const SM_cxx_diagrams::MomentumVertex &v2 )
	{
		BOOST_TEST( v1.index() == v2.index() );
		TEST_COMPLEX_EQUALITY( v1.value( v1.index() ),
			v2.value( v1.index() ) );
	}
	
	void test( const SM_cxx_diagrams::MomentumDifferenceVertex &v1,
		const SM_cxx_diagrams::MomentumDifferenceVertex &v2 )
	{
		if( v1.incoming_index() == v2.incoming_index() )
		{
			BOOST_TEST( v1.outgoing_index() == v2.outgoing_index() );
		} else
		{
			BOOST_TEST( v1.incoming_index() == v2.outgoing_index() );
			BOOST_TEST( v1.outgoing_index() == v2.incoming_index() );
		}
			
		TEST_COMPLEX_EQUALITY(
			v1.value( v1.incoming_index(), v1.outgoing_index() ),
			v2.value( v1.incoming_index(), v1.outgoing_index() )
		);
	}
	
	void test( const SM_cxx_diagrams::InverseMetricVertex &v1,
		const SM_cxx_diagrams::InverseMetricVertex &v2 )
	{ TEST_COMPLEX_EQUALITY( v1.value(), v2.value() ); }
	
	void test( const SM_cxx_diagrams::ChiralVertex &v1,
		const SM_cxx_diagrams::ChiralVertex &v2 )
	{
		TEST_COMPLEX_EQUALITY( v1.left(), v2.left() );
		TEST_COMPLEX_EQUALITY( v1.right(), v2.right() );
	}
	
	void test( const SM_cxx_diagrams::TripleVectorVertex &v1,
		const SM_cxx_diagrams::TripleVectorVertex &v2 )
	{
		using even = SM_cxx_diagrams::TripleVectorVertex::even_permutation;
		TEST_COMPLEX_EQUALITY( v1.value( even{} ), v2.value( even{} ) );
	}
	
	void test( const SM_cxx_diagrams::QuadrupleVectorVertex &v1,
		const SM_cxx_diagrams::QuadrupleVectorVertex &v2 )
	{
		TEST_COMPLEX_EQUALITY( v1.value1(), v2.value1() );
		TEST_COMPLEX_EQUALITY( v1.value2(), v2.value2() );
		TEST_COMPLEX_EQUALITY( v1.value3(), v2.value3() );
	}
};
}
}

#endif
