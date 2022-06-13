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
 * @file cxx_qft/SM_vertices.cpp
 *
 * This file was generated with FlexibleSUSY 2.7.1 and SARAH 4.14.5 .
 */

#include "standard_model_context_base.hpp"
#include "standard_model_vertices.hpp"

#include "concatenate.hpp"
#include "wrappers.hpp"

#define INPUTPARAMETER(p) context.model.get_input_parameters().p
#define MODELPARAMETER(p) context.model.get_##p()
#define DERIVEDPARAMETER(p) context.model.p()
#define PHASE(p) context.model.get_##p()

namespace flexiblesusy {
namespace standard_model_cxx_diagrams {
namespace detail {

cxx_diagrams::ScalarVertex VertexImpl<standard_model_cxx_diagrams::fields::Ah, standard_model_cxx_diagrams::fields::Ah, standard_model_cxx_diagrams::fields::hh, standard_model_cxx_diagrams::fields::hh>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto Lambdax = MODELPARAMETER(Lambdax);

   const std::complex<double> result = -Lambdax;

   return {result};
}

cxx_diagrams::ScalarVertex VertexImpl<standard_model_cxx_diagrams::fields::Ah, standard_model_cxx_diagrams::fields::Ah, standard_model_cxx_diagrams::fields::hh>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto v = MODELPARAMETER(v);
   const auto Lambdax = MODELPARAMETER(Lambdax);

   const std::complex<double> result = -(v*Lambdax);

   return {result};
}

cxx_diagrams::InverseMetricVertex VertexImpl<standard_model_cxx_diagrams::fields::Ah, standard_model_cxx_diagrams::fields::Ah, standard_model_cxx_diagrams::fields::VZ, standard_model_cxx_diagrams::fields::VZ>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*Sqr(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW));

   return {result};
}

cxx_diagrams::InverseMetricVertex VertexImpl<standard_model_cxx_diagrams::fields::Ah, standard_model_cxx_diagrams::fields::Ah, typename standard_model_cxx_diagrams::fields::conj<standard_model_cxx_diagrams::fields::VWp>::type, standard_model_cxx_diagrams::fields::VWp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);

   const std::complex<double> result = 0.5*Sqr(g2);

   return {result};
}

cxx_diagrams::MomentumDifferenceVertex VertexImpl<standard_model_cxx_diagrams::fields::Ah, standard_model_cxx_diagrams::fields::hh, standard_model_cxx_diagrams::fields::VZ>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = std::complex<double>(0,0.1)*(5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW));

   return {result, minuend_index, subtrahend_index};
}

cxx_diagrams::MomentumDifferenceVertex VertexImpl<standard_model_cxx_diagrams::fields::Ah, standard_model_cxx_diagrams::fields::Hp, typename standard_model_cxx_diagrams::fields::conj<standard_model_cxx_diagrams::fields::VWp>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const auto g2 = MODELPARAMETER(g2);

   const std::complex<double> result = std::complex<double>(0,-0.5)*g2;

   return {result, minuend_index, subtrahend_index};
}

cxx_diagrams::ChiralVertex VertexImpl<standard_model_cxx_diagrams::fields::Ah, typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fd>::type, standard_model_cxx_diagrams::fields::Fd>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Yd = MODELPARAMETER(Yd);
   const auto Vd = MODELPARAMETER(Vd);
   const auto Ud = MODELPARAMETER(Ud);

   const std::complex<double> left = std::complex<double>(0.,0.7071067811865475)*SUM(j2,0,2,Conj(Vd(gt2,j2))*SUM(j1,0,2,Conj(Ud(gt1,j1))*Yd(j1,j2)));

   const std::complex<double> right = std::complex<double>(0.,-0.7071067811865475)*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Ud(gt2,j1))*Vd(gt1,j2));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<standard_model_cxx_diagrams::fields::Ah, typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fe>::type, standard_model_cxx_diagrams::fields::Fe>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Ye = MODELPARAMETER(Ye);
   const auto Ve = MODELPARAMETER(Ve);
   const auto Ue = MODELPARAMETER(Ue);

   const std::complex<double> left = std::complex<double>(0.,0.7071067811865475)*SUM(j2,0,2,Conj(Ve(gt2,j2))*SUM(j1,0,2,Conj(Ue(gt1,j1))*Ye(j1,j2)));

   const std::complex<double> right = std::complex<double>(0.,-0.7071067811865475)*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Ue(gt2,j1))*Ve(gt1,j2));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<standard_model_cxx_diagrams::fields::Ah, typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fu>::type, standard_model_cxx_diagrams::fields::Fu>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Yu = MODELPARAMETER(Yu);
   const auto Vu = MODELPARAMETER(Vu);
   const auto Uu = MODELPARAMETER(Uu);

   const std::complex<double> left = std::complex<double>(0.,-0.7071067811865475)*SUM(j2,0,2,Conj(Vu(gt2,j2))*SUM(j1,0,2,Conj(Uu(gt1,j1))*Yu(j1,j2)));

   const std::complex<double> right = std::complex<double>(0.,0.7071067811865475)*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Uu(gt2,j1))*Vu(gt1,j2));

   return {left, right};
}

cxx_diagrams::MomentumDifferenceVertex VertexImpl<standard_model_cxx_diagrams::fields::Ah, typename standard_model_cxx_diagrams::fields::conj<standard_model_cxx_diagrams::fields::Hp>::type, standard_model_cxx_diagrams::fields::VWp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const auto g2 = MODELPARAMETER(g2);

   const std::complex<double> result = std::complex<double>(0,-0.5)*g2;

   return {result, minuend_index, subtrahend_index};
}

cxx_diagrams::ChiralVertex VertexImpl<standard_model_cxx_diagrams::fields::Fe, typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fe>::type, standard_model_cxx_diagrams::fields::VP>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = -0.7745966692414834*g1*Cos(ThetaW)*KroneckerDelta(gt1,gt2);

   const std::complex<double> right = -0.5*KroneckerDelta(gt1,gt2)*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<standard_model_cxx_diagrams::fields::Fe, typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fe>::type, standard_model_cxx_diagrams::fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.7745966692414834*g1*KroneckerDelta(gt1,gt2)*Sin(ThetaW);

   const std::complex<double> right = 0.1*KroneckerDelta(gt1,gt2)*(-5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW));

   return {left, right};
}

cxx_diagrams::ScalarVertex VertexImpl<standard_model_cxx_diagrams::fields::hh, standard_model_cxx_diagrams::fields::hh, standard_model_cxx_diagrams::fields::hh, standard_model_cxx_diagrams::fields::hh>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto Lambdax = MODELPARAMETER(Lambdax);

   const std::complex<double> result = -3*Lambdax;

   return {result};
}

cxx_diagrams::ScalarVertex VertexImpl<standard_model_cxx_diagrams::fields::hh, standard_model_cxx_diagrams::fields::hh, standard_model_cxx_diagrams::fields::hh>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto v = MODELPARAMETER(v);
   const auto Lambdax = MODELPARAMETER(Lambdax);

   const std::complex<double> result = -3*v*Lambdax;

   return {result};
}

cxx_diagrams::ScalarVertex VertexImpl<standard_model_cxx_diagrams::fields::hh, standard_model_cxx_diagrams::fields::hh, standard_model_cxx_diagrams::fields::Hp, typename standard_model_cxx_diagrams::fields::conj<standard_model_cxx_diagrams::fields::Hp>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto Lambdax = MODELPARAMETER(Lambdax);

   const std::complex<double> result = -Lambdax;

   return {result};
}

cxx_diagrams::InverseMetricVertex VertexImpl<standard_model_cxx_diagrams::fields::hh, standard_model_cxx_diagrams::fields::hh, standard_model_cxx_diagrams::fields::VZ, standard_model_cxx_diagrams::fields::VZ>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*Sqr(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW));

   return {result};
}

cxx_diagrams::InverseMetricVertex VertexImpl<standard_model_cxx_diagrams::fields::hh, standard_model_cxx_diagrams::fields::hh, typename standard_model_cxx_diagrams::fields::conj<standard_model_cxx_diagrams::fields::VWp>::type, standard_model_cxx_diagrams::fields::VWp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);

   const std::complex<double> result = 0.5*Sqr(g2);

   return {result};
}

cxx_diagrams::ScalarVertex VertexImpl<standard_model_cxx_diagrams::fields::hh, standard_model_cxx_diagrams::fields::Hp, typename standard_model_cxx_diagrams::fields::conj<standard_model_cxx_diagrams::fields::Hp>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto v = MODELPARAMETER(v);
   const auto Lambdax = MODELPARAMETER(Lambdax);

   const std::complex<double> result = -(v*Lambdax);

   return {result};
}

cxx_diagrams::InverseMetricVertex VertexImpl<standard_model_cxx_diagrams::fields::hh, standard_model_cxx_diagrams::fields::Hp, typename standard_model_cxx_diagrams::fields::conj<standard_model_cxx_diagrams::fields::VWp>::type, standard_model_cxx_diagrams::fields::VP>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.3872983346207417*g1*g2*Cos(ThetaW);

   return {result};
}

cxx_diagrams::InverseMetricVertex VertexImpl<standard_model_cxx_diagrams::fields::hh, standard_model_cxx_diagrams::fields::Hp, typename standard_model_cxx_diagrams::fields::conj<standard_model_cxx_diagrams::fields::VWp>::type, standard_model_cxx_diagrams::fields::VZ>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.3872983346207417*g1*g2*Sin(ThetaW);

   return {result};
}

cxx_diagrams::MomentumDifferenceVertex VertexImpl<standard_model_cxx_diagrams::fields::hh, standard_model_cxx_diagrams::fields::Hp, typename standard_model_cxx_diagrams::fields::conj<standard_model_cxx_diagrams::fields::VWp>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const auto g2 = MODELPARAMETER(g2);

   const std::complex<double> result = 0.5*g2;

   return {result, minuend_index, subtrahend_index};
}

cxx_diagrams::InverseMetricVertex VertexImpl<standard_model_cxx_diagrams::fields::hh, standard_model_cxx_diagrams::fields::VZ, standard_model_cxx_diagrams::fields::VZ>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto v = MODELPARAMETER(v);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*v*Sqr(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW));

   return {result};
}

cxx_diagrams::ChiralVertex VertexImpl<standard_model_cxx_diagrams::fields::hh, typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fd>::type, standard_model_cxx_diagrams::fields::Fd>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Yd = MODELPARAMETER(Yd);
   const auto Vd = MODELPARAMETER(Vd);
   const auto Ud = MODELPARAMETER(Ud);

   const std::complex<double> left = -0.7071067811865475*SUM(j2,0,2,Conj(Vd(gt2,j2))*SUM(j1,0,2,Conj(Ud(gt1,j1))*Yd(j1,j2)));

   const std::complex<double> right = -0.7071067811865475*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Ud(gt2,j1))*Vd(gt1,j2));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<standard_model_cxx_diagrams::fields::hh, typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fe>::type, standard_model_cxx_diagrams::fields::Fe>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Ye = MODELPARAMETER(Ye);
   const auto Ve = MODELPARAMETER(Ve);
   const auto Ue = MODELPARAMETER(Ue);

   const std::complex<double> left = -0.7071067811865475*SUM(j2,0,2,Conj(Ve(gt2,j2))*SUM(j1,0,2,Conj(Ue(gt1,j1))*Ye(j1,j2)));

   const std::complex<double> right = -0.7071067811865475*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Ue(gt2,j1))*Ve(gt1,j2));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<standard_model_cxx_diagrams::fields::hh, typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fu>::type, standard_model_cxx_diagrams::fields::Fu>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Yu = MODELPARAMETER(Yu);
   const auto Vu = MODELPARAMETER(Vu);
   const auto Uu = MODELPARAMETER(Uu);

   const std::complex<double> left = -0.7071067811865475*SUM(j2,0,2,Conj(Vu(gt2,j2))*SUM(j1,0,2,Conj(Uu(gt1,j1))*Yu(j1,j2)));

   const std::complex<double> right = -0.7071067811865475*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Uu(gt2,j1))*Vu(gt1,j2));

   return {left, right};
}

cxx_diagrams::ScalarVertex VertexImpl<standard_model_cxx_diagrams::fields::hh, typename standard_model_cxx_diagrams::fields::conj<standard_model_cxx_diagrams::fields::Hp>::type, standard_model_cxx_diagrams::fields::Hp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto v = MODELPARAMETER(v);
   const auto Lambdax = MODELPARAMETER(Lambdax);

   const std::complex<double> result = -(v*Lambdax);

   return {result};
}

cxx_diagrams::InverseMetricVertex VertexImpl<standard_model_cxx_diagrams::fields::hh, typename standard_model_cxx_diagrams::fields::conj<standard_model_cxx_diagrams::fields::Hp>::type, standard_model_cxx_diagrams::fields::VP, standard_model_cxx_diagrams::fields::VWp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.3872983346207417*g1*g2*Cos(ThetaW);

   return {result};
}

cxx_diagrams::InverseMetricVertex VertexImpl<standard_model_cxx_diagrams::fields::hh, typename standard_model_cxx_diagrams::fields::conj<standard_model_cxx_diagrams::fields::Hp>::type, standard_model_cxx_diagrams::fields::VWp, standard_model_cxx_diagrams::fields::VZ>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.3872983346207417*g1*g2*Sin(ThetaW);

   return {result};
}

cxx_diagrams::MomentumDifferenceVertex VertexImpl<standard_model_cxx_diagrams::fields::hh, typename standard_model_cxx_diagrams::fields::conj<standard_model_cxx_diagrams::fields::Hp>::type, standard_model_cxx_diagrams::fields::VWp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const auto g2 = MODELPARAMETER(g2);

   const std::complex<double> result = -0.5*g2;

   return {result, minuend_index, subtrahend_index};
}

cxx_diagrams::InverseMetricVertex VertexImpl<standard_model_cxx_diagrams::fields::hh, typename standard_model_cxx_diagrams::fields::conj<standard_model_cxx_diagrams::fields::VWp>::type, standard_model_cxx_diagrams::fields::VWp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto v = MODELPARAMETER(v);

   const std::complex<double> result = 0.5*v*Sqr(g2);

   return {result};
}

cxx_diagrams::InverseMetricVertex VertexImpl<standard_model_cxx_diagrams::fields::Hp, typename standard_model_cxx_diagrams::fields::conj<standard_model_cxx_diagrams::fields::Hp>::type, standard_model_cxx_diagrams::fields::VP, standard_model_cxx_diagrams::fields::VP>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*Sqr(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW));

   return {result};
}

cxx_diagrams::InverseMetricVertex VertexImpl<standard_model_cxx_diagrams::fields::Hp, typename standard_model_cxx_diagrams::fields::conj<standard_model_cxx_diagrams::fields::Hp>::type, standard_model_cxx_diagrams::fields::VP, standard_model_cxx_diagrams::fields::VZ>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.05*(7.745966692414834*g1*g2*Cos(2*ThetaW) + Sin(2*ThetaW)*(-3*Sqr(g1) + 5*Sqr(g2)));

   return {result};
}

cxx_diagrams::MomentumDifferenceVertex VertexImpl<standard_model_cxx_diagrams::fields::Hp, typename standard_model_cxx_diagrams::fields::conj<standard_model_cxx_diagrams::fields::Hp>::type, standard_model_cxx_diagrams::fields::VP>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.1*(-3.872983346207417*g1*Cos(ThetaW) - 5*g2*Sin(ThetaW));

   return {result, minuend_index, subtrahend_index};
}

cxx_diagrams::InverseMetricVertex VertexImpl<standard_model_cxx_diagrams::fields::Hp, typename standard_model_cxx_diagrams::fields::conj<standard_model_cxx_diagrams::fields::Hp>::type, standard_model_cxx_diagrams::fields::VZ, standard_model_cxx_diagrams::fields::VZ>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*Sqr(g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW));

   return {result};
}

cxx_diagrams::MomentumDifferenceVertex VertexImpl<standard_model_cxx_diagrams::fields::Hp, typename standard_model_cxx_diagrams::fields::conj<standard_model_cxx_diagrams::fields::Hp>::type, standard_model_cxx_diagrams::fields::VZ>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.1*(-5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW));

   return {result, minuend_index, subtrahend_index};
}

cxx_diagrams::InverseMetricVertex VertexImpl<standard_model_cxx_diagrams::fields::Hp, typename standard_model_cxx_diagrams::fields::conj<standard_model_cxx_diagrams::fields::Hp>::type, typename standard_model_cxx_diagrams::fields::conj<standard_model_cxx_diagrams::fields::VWp>::type, standard_model_cxx_diagrams::fields::VWp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);

   const std::complex<double> result = 0.5*Sqr(g2);

   return {result};
}

cxx_diagrams::InverseMetricVertex VertexImpl<standard_model_cxx_diagrams::fields::Hp, typename standard_model_cxx_diagrams::fields::conj<standard_model_cxx_diagrams::fields::VWp>::type, standard_model_cxx_diagrams::fields::VP>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto v = MODELPARAMETER(v);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.3872983346207417*g1*g2*v*Cos(ThetaW);

   return {result};
}

cxx_diagrams::InverseMetricVertex VertexImpl<standard_model_cxx_diagrams::fields::Hp, typename standard_model_cxx_diagrams::fields::conj<standard_model_cxx_diagrams::fields::VWp>::type, standard_model_cxx_diagrams::fields::VZ>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto v = MODELPARAMETER(v);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.3872983346207417*g1*g2*v*Sin(ThetaW);

   return {result};
}

cxx_diagrams::ChiralVertex VertexImpl<standard_model_cxx_diagrams::fields::VP, standard_model_cxx_diagrams::fields::Fd, typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fd>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = -0.2581988897471611*g1*Cos(ThetaW)*KroneckerDelta(gt1,gt2);

   const std::complex<double> right = 0.16666666666666666*KroneckerDelta(gt1,gt2)*(0.7745966692414834*g1*Cos(ThetaW) - 3*g2*Sin(ThetaW));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<standard_model_cxx_diagrams::fields::VP, standard_model_cxx_diagrams::fields::Fe, typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fe>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = -0.7745966692414834*g1*Cos(ThetaW)*KroneckerDelta(gt1,gt2);

   const std::complex<double> right = -0.5*KroneckerDelta(gt1,gt2)*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<standard_model_cxx_diagrams::fields::VP, standard_model_cxx_diagrams::fields::Fu, typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fu>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.5163977794943222*g1*Cos(ThetaW)*KroneckerDelta(gt1,gt2);

   const std::complex<double> right = 0.16666666666666666*KroneckerDelta(gt1,gt2)*(0.7745966692414834*g1*Cos(ThetaW) + 3*g2*Sin(ThetaW));

   return {left, right};
}

cxx_diagrams::MomentumDifferenceVertex VertexImpl<standard_model_cxx_diagrams::fields::VP, standard_model_cxx_diagrams::fields::Hp, typename standard_model_cxx_diagrams::fields::conj<standard_model_cxx_diagrams::fields::Hp>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   int minuend_index = 1;
   int subtrahend_index = 2;

   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.1*(-3.872983346207417*g1*Cos(ThetaW) - 5*g2*Sin(ThetaW));

   return {result, minuend_index, subtrahend_index};
}

cxx_diagrams::TripleVectorVertex VertexImpl<standard_model_cxx_diagrams::fields::VP, standard_model_cxx_diagrams::fields::VWp, typename standard_model_cxx_diagrams::fields::conj<standard_model_cxx_diagrams::fields::VWp>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -(g2*Sin(ThetaW));

   return {result, cxx_diagrams::TripleVectorVertex::odd_permutation{}};
}

cxx_diagrams::ChiralVertex VertexImpl<standard_model_cxx_diagrams::fields::VZ, standard_model_cxx_diagrams::fields::Fd, typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fd>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.2581988897471611*g1*KroneckerDelta(gt1,gt2)*Sin(ThetaW);

   const std::complex<double> right = -0.16666666666666666*KroneckerDelta(gt1,gt2)*(3*g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<standard_model_cxx_diagrams::fields::VZ, standard_model_cxx_diagrams::fields::Fe, typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fe>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.7745966692414834*g1*KroneckerDelta(gt1,gt2)*Sin(ThetaW);

   const std::complex<double> right = 0.1*KroneckerDelta(gt1,gt2)*(-5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<standard_model_cxx_diagrams::fields::VZ, standard_model_cxx_diagrams::fields::Fu, typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fu>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt2 = indices[0];
   const int gt1 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = -0.5163977794943222*g1*KroneckerDelta(gt1,gt2)*Sin(ThetaW);

   const std::complex<double> right = 0.16666666666666666*KroneckerDelta(gt1,gt2)*(3*g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fd>::type, standard_model_cxx_diagrams::fields::Fd, standard_model_cxx_diagrams::fields::Ah>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Yd = MODELPARAMETER(Yd);
   const auto Vd = MODELPARAMETER(Vd);
   const auto Ud = MODELPARAMETER(Ud);

   const std::complex<double> left = std::complex<double>(0.,0.7071067811865475)*SUM(j2,0,2,Conj(Vd(gt2,j2))*SUM(j1,0,2,Conj(Ud(gt1,j1))*Yd(j1,j2)));

   const std::complex<double> right = std::complex<double>(0.,-0.7071067811865475)*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Ud(gt2,j1))*Vd(gt1,j2));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fd>::type, standard_model_cxx_diagrams::fields::Fd, standard_model_cxx_diagrams::fields::hh>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Yd = MODELPARAMETER(Yd);
   const auto Vd = MODELPARAMETER(Vd);
   const auto Ud = MODELPARAMETER(Ud);

   const std::complex<double> left = -0.7071067811865475*SUM(j2,0,2,Conj(Vd(gt2,j2))*SUM(j1,0,2,Conj(Ud(gt1,j1))*Yd(j1,j2)));

   const std::complex<double> right = -0.7071067811865475*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Ud(gt2,j1))*Vd(gt1,j2));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fd>::type, standard_model_cxx_diagrams::fields::Fd, standard_model_cxx_diagrams::fields::VG>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g3 = MODELPARAMETER(g3);

   const std::complex<double> left = -0.5*g3*KroneckerDelta(gt1,gt2);

   const std::complex<double> right = -0.5*g3*KroneckerDelta(gt1,gt2);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fd>::type, standard_model_cxx_diagrams::fields::Fd, standard_model_cxx_diagrams::fields::VP>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = -0.16666666666666666*KroneckerDelta(gt1,gt2)*(0.7745966692414834*g1*Cos(ThetaW) - 3*g2*Sin(ThetaW));

   const std::complex<double> right = 0.2581988897471611*g1*Cos(ThetaW)*KroneckerDelta(gt1,gt2);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fd>::type, standard_model_cxx_diagrams::fields::Fd, standard_model_cxx_diagrams::fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.16666666666666666*KroneckerDelta(gt1,gt2)*(3*g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW));

   const std::complex<double> right = -0.2581988897471611*g1*KroneckerDelta(gt1,gt2)*Sin(ThetaW);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fd>::type, standard_model_cxx_diagrams::fields::Fu, typename standard_model_cxx_diagrams::fields::conj<standard_model_cxx_diagrams::fields::Hp>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Yd = MODELPARAMETER(Yd);
   const auto Yu = MODELPARAMETER(Yu);
   const auto Vu = MODELPARAMETER(Vu);
   const auto Ud = MODELPARAMETER(Ud);
   const auto Uu = MODELPARAMETER(Uu);
   const auto Vd = MODELPARAMETER(Vd);

   const std::complex<double> left = -SUM(j2,0,2,Conj(Vu(gt2,j2))*SUM(j1,0,2,Conj(Ud(gt1,j1))*Yd(j1,j2)));

   const std::complex<double> right = SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Uu(gt2,j1))*Vd(gt1,j2));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fd>::type, standard_model_cxx_diagrams::fields::Fu, typename standard_model_cxx_diagrams::fields::conj<standard_model_cxx_diagrams::fields::VWp>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto Vu = MODELPARAMETER(Vu);
   const auto Vd = MODELPARAMETER(Vd);

   const std::complex<double> left = -0.7071067811865475*g2*SUM(j1,0,2,Conj(Vu(gt2,j1))*Vd(gt1,j1));

   const std::complex<double> right = 0;

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fe>::type, standard_model_cxx_diagrams::fields::Fe, standard_model_cxx_diagrams::fields::Ah>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Ye = MODELPARAMETER(Ye);
   const auto Ve = MODELPARAMETER(Ve);
   const auto Ue = MODELPARAMETER(Ue);

   const std::complex<double> left = std::complex<double>(0.,0.7071067811865475)*SUM(j2,0,2,Conj(Ve(gt2,j2))*SUM(j1,0,2,Conj(Ue(gt1,j1))*Ye(j1,j2)));

   const std::complex<double> right = std::complex<double>(0.,-0.7071067811865475)*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Ue(gt2,j1))*Ve(gt1,j2));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fe>::type, standard_model_cxx_diagrams::fields::Fe, standard_model_cxx_diagrams::fields::hh>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Ye = MODELPARAMETER(Ye);
   const auto Ve = MODELPARAMETER(Ve);
   const auto Ue = MODELPARAMETER(Ue);

   const std::complex<double> left = -0.7071067811865475*SUM(j2,0,2,Conj(Ve(gt2,j2))*SUM(j1,0,2,Conj(Ue(gt1,j1))*Ye(j1,j2)));

   const std::complex<double> right = -0.7071067811865475*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Ue(gt2,j1))*Ve(gt1,j2));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fe>::type, standard_model_cxx_diagrams::fields::Fe, standard_model_cxx_diagrams::fields::VP>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.5*KroneckerDelta(gt1,gt2)*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW));

   const std::complex<double> right = 0.7745966692414834*g1*Cos(ThetaW)*KroneckerDelta(gt1,gt2);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fe>::type, standard_model_cxx_diagrams::fields::Fe, standard_model_cxx_diagrams::fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.5*KroneckerDelta(gt1,gt2)*(g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW));

   const std::complex<double> right = -0.7745966692414834*g1*KroneckerDelta(gt1,gt2)*Sin(ThetaW);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fe>::type, standard_model_cxx_diagrams::fields::Fv, typename standard_model_cxx_diagrams::fields::conj<standard_model_cxx_diagrams::fields::Hp>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Ye = MODELPARAMETER(Ye);
   const auto Ue = MODELPARAMETER(Ue);

   const std::complex<double> left = -SUM(j1,0,2,Conj(Ue(gt1,j1))*Ye(j1,gt2));

   const std::complex<double> right = 0;

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fe>::type, standard_model_cxx_diagrams::fields::Fv, typename standard_model_cxx_diagrams::fields::conj<standard_model_cxx_diagrams::fields::VWp>::type>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto Ve = MODELPARAMETER(Ve);

   const std::complex<double> left = IF(gt2 < 3,-0.7071067811865475*g2*Ve(gt1,gt2),0);

   const std::complex<double> right = 0;

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fe>::type, typename standard_model_cxx_diagrams::fields::conj<standard_model_cxx_diagrams::fields::Hp>::type, standard_model_cxx_diagrams::fields::Fv>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Ye = MODELPARAMETER(Ye);
   const auto Ue = MODELPARAMETER(Ue);

   const std::complex<double> left = -SUM(j1,0,2,Conj(Ue(gt1,j1))*Ye(j1,gt2));

   const std::complex<double> right = 0;

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fu>::type, standard_model_cxx_diagrams::fields::Fd, standard_model_cxx_diagrams::fields::Hp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Yu = MODELPARAMETER(Yu);
   const auto Yd = MODELPARAMETER(Yd);
   const auto Vd = MODELPARAMETER(Vd);
   const auto Uu = MODELPARAMETER(Uu);
   const auto Ud = MODELPARAMETER(Ud);
   const auto Vu = MODELPARAMETER(Vu);

   const std::complex<double> left = SUM(j2,0,2,Conj(Vd(gt2,j2))*SUM(j1,0,2,Conj(Uu(gt1,j1))*Yu(j1,j2)));

   const std::complex<double> right = -SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Ud(gt2,j1))*Vu(gt1,j2));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fu>::type, standard_model_cxx_diagrams::fields::Fd, standard_model_cxx_diagrams::fields::VWp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto Vd = MODELPARAMETER(Vd);
   const auto Vu = MODELPARAMETER(Vu);

   const std::complex<double> left = -0.7071067811865475*g2*SUM(j1,0,2,Conj(Vd(gt2,j1))*Vu(gt1,j1));

   const std::complex<double> right = 0;

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fu>::type, standard_model_cxx_diagrams::fields::Fu, standard_model_cxx_diagrams::fields::Ah>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Yu = MODELPARAMETER(Yu);
   const auto Vu = MODELPARAMETER(Vu);
   const auto Uu = MODELPARAMETER(Uu);

   const std::complex<double> left = std::complex<double>(0.,-0.7071067811865475)*SUM(j2,0,2,Conj(Vu(gt2,j2))*SUM(j1,0,2,Conj(Uu(gt1,j1))*Yu(j1,j2)));

   const std::complex<double> right = std::complex<double>(0.,0.7071067811865475)*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Uu(gt2,j1))*Vu(gt1,j2));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fu>::type, standard_model_cxx_diagrams::fields::Fu, standard_model_cxx_diagrams::fields::hh>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Yu = MODELPARAMETER(Yu);
   const auto Vu = MODELPARAMETER(Vu);
   const auto Uu = MODELPARAMETER(Uu);

   const std::complex<double> left = -0.7071067811865475*SUM(j2,0,2,Conj(Vu(gt2,j2))*SUM(j1,0,2,Conj(Uu(gt1,j1))*Yu(j1,j2)));

   const std::complex<double> right = -0.7071067811865475*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Uu(gt2,j1))*Vu(gt1,j2));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fu>::type, standard_model_cxx_diagrams::fields::Fu, standard_model_cxx_diagrams::fields::VG>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g3 = MODELPARAMETER(g3);

   const std::complex<double> left = -0.5*g3*KroneckerDelta(gt1,gt2);

   const std::complex<double> right = -0.5*g3*KroneckerDelta(gt1,gt2);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fu>::type, standard_model_cxx_diagrams::fields::Fu, standard_model_cxx_diagrams::fields::VP>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = -0.16666666666666666*KroneckerDelta(gt1,gt2)*(0.7745966692414834*g1*Cos(ThetaW) + 3*g2*Sin(ThetaW));

   const std::complex<double> right = -0.5163977794943222*g1*Cos(ThetaW)*KroneckerDelta(gt1,gt2);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fu>::type, standard_model_cxx_diagrams::fields::Fu, standard_model_cxx_diagrams::fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = 0.03333333333333333*KroneckerDelta(gt1,gt2)*(-15*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW));

   const std::complex<double> right = 0.5163977794943222*g1*KroneckerDelta(gt1,gt2)*Sin(ThetaW);

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fv>::type, standard_model_cxx_diagrams::fields::Fe, standard_model_cxx_diagrams::fields::Hp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Ye = MODELPARAMETER(Ye);
   const auto Ue = MODELPARAMETER(Ue);

   const std::complex<double> left = 0;

   const std::complex<double> right = -SUM(j1,0,2,Conj(Ye(j1,gt1))*Ue(gt2,j1));

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fv>::type, standard_model_cxx_diagrams::fields::Fe, standard_model_cxx_diagrams::fields::VWp>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto Ve = MODELPARAMETER(Ve);

   const std::complex<double> left = IF(gt1 < 3,-0.7071067811865475*g2*Conj(Ve(gt2,gt1)),0);

   const std::complex<double> right = 0;

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fv>::type, standard_model_cxx_diagrams::fields::Fv, standard_model_cxx_diagrams::fields::VZ>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> left = -0.5*KroneckerDelta(gt1,gt2)*(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW));

   const std::complex<double> right = 0;

   return {left, right};
}

cxx_diagrams::ChiralVertex VertexImpl<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::Fv>::type, standard_model_cxx_diagrams::fields::Hp, standard_model_cxx_diagrams::fields::Fe>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto Ye = MODELPARAMETER(Ye);
   const auto Ue = MODELPARAMETER(Ue);

   const std::complex<double> left = 0;

   const std::complex<double> right = -SUM(j1,0,2,Conj(Ye(j1,gt1))*Ue(gt2,j1));

   return {left, right};
}

cxx_diagrams::MomentumVertex VertexImpl<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::gP>::type, standard_model_cxx_diagrams::fields::gWpC, standard_model_cxx_diagrams::fields::VWp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = g2*Sin(ThetaW);

   return {result, 1};
}

cxx_diagrams::MomentumVertex VertexImpl<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::gP>::type, standard_model_cxx_diagrams::fields::gWp, typename standard_model_cxx_diagrams::fields::conj<standard_model_cxx_diagrams::fields::VWp>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -(g2*Sin(ThetaW));

   return {result, 1};
}

cxx_diagrams::ScalarVertex VertexImpl<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::gWpC>::type, standard_model_cxx_diagrams::fields::gP, typename standard_model_cxx_diagrams::fields::conj<standard_model_cxx_diagrams::fields::Hp>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto v = MODELPARAMETER(v);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.25*g2*v*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW));

   return {result};
}

cxx_diagrams::MomentumVertex VertexImpl<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::gWpC>::type, standard_model_cxx_diagrams::fields::gP, typename standard_model_cxx_diagrams::fields::conj<standard_model_cxx_diagrams::fields::VWp>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = g2*Sin(ThetaW);

   return {result, 1};
}

cxx_diagrams::ScalarVertex VertexImpl<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::gWpC>::type, standard_model_cxx_diagrams::fields::gWpC, standard_model_cxx_diagrams::fields::Ah>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto v = MODELPARAMETER(v);

   const std::complex<double> result = std::complex<double>(0,0.25)*v*Sqr(g2);

   return {result};
}

cxx_diagrams::ScalarVertex VertexImpl<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::gWpC>::type, standard_model_cxx_diagrams::fields::gWpC, standard_model_cxx_diagrams::fields::hh>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto v = MODELPARAMETER(v);

   const std::complex<double> result = -0.25*v*Sqr(g2);

   return {result};
}

cxx_diagrams::MomentumVertex VertexImpl<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::gWpC>::type, standard_model_cxx_diagrams::fields::gWpC, standard_model_cxx_diagrams::fields::VP>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -(g2*Sin(ThetaW));

   return {result, 1};
}

cxx_diagrams::MomentumVertex VertexImpl<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::gWpC>::type, standard_model_cxx_diagrams::fields::gWpC, standard_model_cxx_diagrams::fields::VZ>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -(g2*Cos(ThetaW));

   return {result, 1};
}

cxx_diagrams::ScalarVertex VertexImpl<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::gWpC>::type, standard_model_cxx_diagrams::fields::gZ, typename standard_model_cxx_diagrams::fields::conj<standard_model_cxx_diagrams::fields::Hp>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto v = MODELPARAMETER(v);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.05*g2*v*(-5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW));

   return {result};
}

cxx_diagrams::MomentumVertex VertexImpl<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::gWpC>::type, standard_model_cxx_diagrams::fields::gZ, typename standard_model_cxx_diagrams::fields::conj<standard_model_cxx_diagrams::fields::VWp>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = g2*Cos(ThetaW);

   return {result, 1};
}

cxx_diagrams::ScalarVertex VertexImpl<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::gWp>::type, standard_model_cxx_diagrams::fields::gP, standard_model_cxx_diagrams::fields::Hp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto v = MODELPARAMETER(v);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.25*g2*v*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW));

   return {result};
}

cxx_diagrams::MomentumVertex VertexImpl<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::gWp>::type, standard_model_cxx_diagrams::fields::gP, standard_model_cxx_diagrams::fields::VWp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -(g2*Sin(ThetaW));

   return {result, 1};
}

cxx_diagrams::ScalarVertex VertexImpl<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::gWp>::type, standard_model_cxx_diagrams::fields::gWp, standard_model_cxx_diagrams::fields::Ah>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto v = MODELPARAMETER(v);

   const std::complex<double> result = std::complex<double>(0,-0.25)*v*Sqr(g2);

   return {result};
}

cxx_diagrams::ScalarVertex VertexImpl<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::gWp>::type, standard_model_cxx_diagrams::fields::gWp, standard_model_cxx_diagrams::fields::hh>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto v = MODELPARAMETER(v);

   const std::complex<double> result = -0.25*v*Sqr(g2);

   return {result};
}

cxx_diagrams::MomentumVertex VertexImpl<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::gWp>::type, standard_model_cxx_diagrams::fields::gWp, standard_model_cxx_diagrams::fields::VP>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = g2*Sin(ThetaW);

   return {result, 1};
}

cxx_diagrams::MomentumVertex VertexImpl<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::gWp>::type, standard_model_cxx_diagrams::fields::gWp, standard_model_cxx_diagrams::fields::VZ>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = g2*Cos(ThetaW);

   return {result, 1};
}

cxx_diagrams::ScalarVertex VertexImpl<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::gWp>::type, standard_model_cxx_diagrams::fields::gZ, standard_model_cxx_diagrams::fields::Hp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto v = MODELPARAMETER(v);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.05*g2*v*(-5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW));

   return {result};
}

cxx_diagrams::MomentumVertex VertexImpl<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::gWp>::type, standard_model_cxx_diagrams::fields::gZ, standard_model_cxx_diagrams::fields::VWp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -(g2*Cos(ThetaW));

   return {result, 1};
}

cxx_diagrams::ScalarVertex VertexImpl<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::gZ>::type, standard_model_cxx_diagrams::fields::gP, standard_model_cxx_diagrams::fields::hh>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto v = MODELPARAMETER(v);
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.025*v*(7.745966692414834*g1*g2*Cos(2*ThetaW) + Sin(2*ThetaW)*(3*Sqr(g1) - 5*Sqr(g2)));

   return {result};
}

cxx_diagrams::ScalarVertex VertexImpl<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::gZ>::type, standard_model_cxx_diagrams::fields::gWpC, standard_model_cxx_diagrams::fields::Hp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto v = MODELPARAMETER(v);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.25*g2*v*(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW));

   return {result};
}

cxx_diagrams::MomentumVertex VertexImpl<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::gZ>::type, standard_model_cxx_diagrams::fields::gWpC, standard_model_cxx_diagrams::fields::VWp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = g2*Cos(ThetaW);

   return {result, 1};
}

cxx_diagrams::ScalarVertex VertexImpl<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::gZ>::type, standard_model_cxx_diagrams::fields::gWp, typename standard_model_cxx_diagrams::fields::conj<standard_model_cxx_diagrams::fields::Hp>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto v = MODELPARAMETER(v);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.25*g2*v*(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW));

   return {result};
}

cxx_diagrams::MomentumVertex VertexImpl<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::gZ>::type, standard_model_cxx_diagrams::fields::gWp, typename standard_model_cxx_diagrams::fields::conj<standard_model_cxx_diagrams::fields::VWp>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -(g2*Cos(ThetaW));

   return {result, 1};
}

cxx_diagrams::ScalarVertex VertexImpl<typename standard_model_cxx_diagrams::fields::bar<standard_model_cxx_diagrams::fields::gZ>::type, standard_model_cxx_diagrams::fields::gZ, standard_model_cxx_diagrams::fields::hh>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto v = MODELPARAMETER(v);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.25*v*Sqr(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW));

   return {result};
}

cxx_diagrams::InverseMetricVertex VertexImpl<typename standard_model_cxx_diagrams::fields::conj<standard_model_cxx_diagrams::fields::Hp>::type, standard_model_cxx_diagrams::fields::VP, standard_model_cxx_diagrams::fields::VWp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto v = MODELPARAMETER(v);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.3872983346207417*g1*g2*v*Cos(ThetaW);

   return {result};
}

cxx_diagrams::InverseMetricVertex VertexImpl<typename standard_model_cxx_diagrams::fields::conj<standard_model_cxx_diagrams::fields::Hp>::type, standard_model_cxx_diagrams::fields::VWp, standard_model_cxx_diagrams::fields::VZ>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto v = MODELPARAMETER(v);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.3872983346207417*g1*g2*v*Sin(ThetaW);

   return {result};
}

cxx_diagrams::QuadrupleVectorVertex VertexImpl<typename standard_model_cxx_diagrams::fields::conj<standard_model_cxx_diagrams::fields::VWp>::type, standard_model_cxx_diagrams::fields::VP, standard_model_cxx_diagrams::fields::VP, standard_model_cxx_diagrams::fields::VWp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> part1 = Sqr(g2)*Sqr(Sin(ThetaW));

   const std::complex<double> part2 = Sqr(g2)*Sqr(Sin(ThetaW));

   const std::complex<double> part3 = -2*Sqr(g2)*Sqr(Sin(ThetaW));

   return {part1, part2, part3};
}

cxx_diagrams::QuadrupleVectorVertex VertexImpl<typename standard_model_cxx_diagrams::fields::conj<standard_model_cxx_diagrams::fields::VWp>::type, standard_model_cxx_diagrams::fields::VP, standard_model_cxx_diagrams::fields::VWp, standard_model_cxx_diagrams::fields::VZ>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> part1 = Cos(ThetaW)*Sin(ThetaW)*Sqr(g2);

   const std::complex<double> part2 = -(Sin(2*ThetaW)*Sqr(g2));

   const std::complex<double> part3 = Cos(ThetaW)*Sin(ThetaW)*Sqr(g2);

   return {part1, part2, part3};
}

cxx_diagrams::TripleVectorVertex VertexImpl<typename standard_model_cxx_diagrams::fields::conj<standard_model_cxx_diagrams::fields::VWp>::type, standard_model_cxx_diagrams::fields::VP, standard_model_cxx_diagrams::fields::VWp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -(g2*Sin(ThetaW));

   return {result, cxx_diagrams::TripleVectorVertex::odd_permutation{}};
}

cxx_diagrams::QuadrupleVectorVertex VertexImpl<typename standard_model_cxx_diagrams::fields::conj<standard_model_cxx_diagrams::fields::VWp>::type, standard_model_cxx_diagrams::fields::VWp, standard_model_cxx_diagrams::fields::VZ, standard_model_cxx_diagrams::fields::VZ>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> part1 = -2*Sqr(g2)*Sqr(Cos(ThetaW));

   const std::complex<double> part2 = Sqr(g2)*Sqr(Cos(ThetaW));

   const std::complex<double> part3 = Sqr(g2)*Sqr(Cos(ThetaW));

   return {part1, part2, part3};
}

cxx_diagrams::TripleVectorVertex VertexImpl<typename standard_model_cxx_diagrams::fields::conj<standard_model_cxx_diagrams::fields::VWp>::type, standard_model_cxx_diagrams::fields::VWp, standard_model_cxx_diagrams::fields::VZ>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = g2*Cos(ThetaW);

   return {result, cxx_diagrams::TripleVectorVertex::odd_permutation{}};
}

cxx_diagrams::QuadrupleVectorVertex VertexImpl<typename standard_model_cxx_diagrams::fields::conj<standard_model_cxx_diagrams::fields::VWp>::type, typename standard_model_cxx_diagrams::fields::conj<standard_model_cxx_diagrams::fields::VWp>::type, standard_model_cxx_diagrams::fields::VWp, standard_model_cxx_diagrams::fields::VWp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);

   const std::complex<double> part1 = 2*Sqr(g2);

   const std::complex<double> part2 = -Sqr(g2);

   const std::complex<double> part3 = -Sqr(g2);

   return {part1, part2, part3};
}

} // namespace detail
} // namespace standard_model_cxx_diagrams
} // namespace flexiblesusy
