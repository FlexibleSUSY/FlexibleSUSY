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
 * This file was generated with FlexibleSUSY 2.6.0 and SARAH 4.14.3 .
 */

#include "SM_context_base.hpp"
#include "SM_input_parameters.hpp"
#include "SM_vertices.hpp"

#include "concatenate.hpp"
#include "wrappers.hpp"

#define INPUTPARAMETER(p) context.model.get_input_parameters().p
#define MODELPARAMETER(p) context.model.get_##p()
#define DERIVEDPARAMETER(p) context.model.p()
#define PHASE(p) context.model.get_##p()

namespace flexiblesusy {
namespace SM_cxx_diagrams {
namespace detail {

ScalarVertex VertexImpl<fields::Ah, fields::Ah, fields::hh, fields::hh>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto Lambdax = MODELPARAMETER(Lambdax);

   const std::complex<double> result = -Lambdax;

   return {result};
}

ScalarVertex VertexImpl<fields::Ah, fields::Ah, fields::hh>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto v = MODELPARAMETER(v);
   const auto Lambdax = MODELPARAMETER(Lambdax);

   const std::complex<double> result = -(v*Lambdax);

   return {result};
}

InverseMetricVertex VertexImpl<fields::Ah, fields::Ah, fields::VZ, fields::VZ>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*Sqr(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Ah, fields::Ah, typename fields::conj<fields::VWp>::type, fields::VWp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);

   const std::complex<double> result = 0.5*Sqr(g2);

   return {result};
}

MomentumDifferenceVertex VertexImpl<fields::Ah, fields::hh, fields::VZ>::evaluate(
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

MomentumDifferenceVertex VertexImpl<fields::Ah, fields::Hp, typename fields::conj<fields::VWp>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const auto g2 = MODELPARAMETER(g2);

   const std::complex<double> result = std::complex<double>(0,-0.5)*g2;

   return {result, minuend_index, subtrahend_index};
}

ChiralVertex VertexImpl<fields::Ah, typename fields::bar<fields::Fe>::type, fields::Fe>::evaluate(
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

MomentumDifferenceVertex VertexImpl<fields::Ah, typename fields::conj<fields::Hp>::type, fields::VWp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const auto g2 = MODELPARAMETER(g2);

   const std::complex<double> result = std::complex<double>(0,-0.5)*g2;

   return {result, minuend_index, subtrahend_index};
}

ScalarVertex VertexImpl<fields::hh, fields::hh, fields::hh, fields::hh>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto Lambdax = MODELPARAMETER(Lambdax);

   const std::complex<double> result = -3*Lambdax;

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::hh, fields::hh>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto v = MODELPARAMETER(v);
   const auto Lambdax = MODELPARAMETER(Lambdax);

   const std::complex<double> result = -3*v*Lambdax;

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::hh, fields::Hp, typename fields::conj<fields::Hp>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto Lambdax = MODELPARAMETER(Lambdax);

   const std::complex<double> result = -Lambdax;

   return {result};
}

InverseMetricVertex VertexImpl<fields::hh, fields::hh, fields::VZ, fields::VZ>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*Sqr(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW));

   return {result};
}

InverseMetricVertex VertexImpl<fields::hh, fields::hh, typename fields::conj<fields::VWp>::type, fields::VWp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);

   const std::complex<double> result = 0.5*Sqr(g2);

   return {result};
}

ScalarVertex VertexImpl<fields::hh, fields::Hp, typename fields::conj<fields::Hp>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto v = MODELPARAMETER(v);
   const auto Lambdax = MODELPARAMETER(Lambdax);

   const std::complex<double> result = -(v*Lambdax);

   return {result};
}

InverseMetricVertex VertexImpl<fields::hh, fields::Hp, typename fields::conj<fields::VWp>::type, fields::VP>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.3872983346207417*g1*g2*Cos(ThetaW);

   return {result};
}

InverseMetricVertex VertexImpl<fields::hh, fields::Hp, typename fields::conj<fields::VWp>::type, fields::VZ>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.3872983346207417*g1*g2*Sin(ThetaW);

   return {result};
}

MomentumDifferenceVertex VertexImpl<fields::hh, fields::Hp, typename fields::conj<fields::VWp>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const auto g2 = MODELPARAMETER(g2);

   const std::complex<double> result = 0.5*g2;

   return {result, minuend_index, subtrahend_index};
}

InverseMetricVertex VertexImpl<fields::hh, fields::VZ, fields::VZ>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto v = MODELPARAMETER(v);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*v*Sqr(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW));

   return {result};
}

ChiralVertex VertexImpl<fields::hh, typename fields::bar<fields::Fe>::type, fields::Fe>::evaluate(
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

InverseMetricVertex VertexImpl<fields::hh, typename fields::conj<fields::Hp>::type, fields::VP, fields::VWp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.3872983346207417*g1*g2*Cos(ThetaW);

   return {result};
}

InverseMetricVertex VertexImpl<fields::hh, typename fields::conj<fields::Hp>::type, fields::VWp, fields::VZ>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.3872983346207417*g1*g2*Sin(ThetaW);

   return {result};
}

MomentumDifferenceVertex VertexImpl<fields::hh, typename fields::conj<fields::Hp>::type, fields::VWp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   int minuend_index = 0;
   int subtrahend_index = 1;

   const auto g2 = MODELPARAMETER(g2);

   const std::complex<double> result = -0.5*g2;

   return {result, minuend_index, subtrahend_index};
}

InverseMetricVertex VertexImpl<fields::hh, typename fields::conj<fields::VWp>::type, fields::VWp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto v = MODELPARAMETER(v);

   const std::complex<double> result = 0.5*v*Sqr(g2);

   return {result};
}

InverseMetricVertex VertexImpl<fields::Hp, typename fields::conj<fields::Hp>::type, fields::VP, fields::VP>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*Sqr(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW));

   return {result};
}

InverseMetricVertex VertexImpl<fields::Hp, typename fields::conj<fields::Hp>::type, fields::VP, fields::VZ>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.05*(7.745966692414834*g1*g2*Cos(2*ThetaW) + Sin(2*ThetaW)*(-3*Sqr(g1) + 5*Sqr(g2)));

   return {result};
}

MomentumDifferenceVertex VertexImpl<fields::Hp, typename fields::conj<fields::Hp>::type, fields::VP>::evaluate(
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

InverseMetricVertex VertexImpl<fields::Hp, typename fields::conj<fields::Hp>::type, fields::VZ, fields::VZ>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.5*Sqr(g2*Cos(ThetaW) - 0.7745966692414834*g1*Sin(ThetaW));

   return {result};
}

MomentumDifferenceVertex VertexImpl<fields::Hp, typename fields::conj<fields::Hp>::type, fields::VZ>::evaluate(
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

InverseMetricVertex VertexImpl<fields::Hp, typename fields::conj<fields::Hp>::type, typename fields::conj<fields::VWp>::type, fields::VWp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);

   const std::complex<double> result = 0.5*Sqr(g2);

   return {result};
}

InverseMetricVertex VertexImpl<fields::Hp, typename fields::conj<fields::VWp>::type, fields::VP>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto v = MODELPARAMETER(v);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.3872983346207417*g1*g2*v*Cos(ThetaW);

   return {result};
}

InverseMetricVertex VertexImpl<fields::Hp, typename fields::conj<fields::VWp>::type, fields::VZ>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto v = MODELPARAMETER(v);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.3872983346207417*g1*g2*v*Sin(ThetaW);

   return {result};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fd>::type, fields::Fd, fields::Ah>::evaluate(
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

ChiralVertex VertexImpl<typename fields::bar<fields::Fd>::type, fields::Fd, fields::hh>::evaluate(
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

ChiralVertex VertexImpl<typename fields::bar<fields::Fd>::type, fields::Fd, fields::VG>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g3 = MODELPARAMETER(g3);

   const std::complex<double> left = -0.5*g3*KroneckerDelta(gt1,gt2);

   const std::complex<double> right = -0.5*g3*KroneckerDelta(gt1,gt2);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fd>::type, fields::Fd, fields::VP>::evaluate(
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

ChiralVertex VertexImpl<typename fields::bar<fields::Fd>::type, fields::Fd, fields::VZ>::evaluate(
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

ChiralVertex VertexImpl<typename fields::bar<fields::Fd>::type, fields::Fu, typename fields::conj<fields::Hp>::type>::evaluate(
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

ChiralVertex VertexImpl<typename fields::bar<fields::Fd>::type, fields::Fu, typename fields::conj<fields::VWp>::type>::evaluate(
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

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, fields::Fe, fields::Ah>::evaluate(
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

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, fields::Fe, fields::hh>::evaluate(
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

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, fields::Fe, fields::VP>::evaluate(
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

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, fields::Fe, fields::VZ>::evaluate(
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

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, fields::Fv, typename fields::conj<fields::Hp>::type>::evaluate(
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

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, fields::Fv, typename fields::conj<fields::VWp>::type>::evaluate(
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

ChiralVertex VertexImpl<typename fields::bar<fields::Fe>::type, typename fields::conj<fields::Hp>::type, fields::Fv>::evaluate(
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

ChiralVertex VertexImpl<typename fields::bar<fields::Fu>::type, fields::Fd, fields::Hp>::evaluate(
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

ChiralVertex VertexImpl<typename fields::bar<fields::Fu>::type, fields::Fd, fields::VWp>::evaluate(
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

ChiralVertex VertexImpl<typename fields::bar<fields::Fu>::type, fields::Fu, fields::Ah>::evaluate(
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

ChiralVertex VertexImpl<typename fields::bar<fields::Fu>::type, fields::Fu, fields::hh>::evaluate(
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

ChiralVertex VertexImpl<typename fields::bar<fields::Fu>::type, fields::Fu, fields::VG>::evaluate(
   const std::array<int, 2>& indices, const context_base& context)
{
   const int gt1 = indices[0];
   const int gt2 = indices[1];
   const auto g3 = MODELPARAMETER(g3);

   const std::complex<double> left = -0.5*g3*KroneckerDelta(gt1,gt2);

   const std::complex<double> right = -0.5*g3*KroneckerDelta(gt1,gt2);

   return {left, right};
}

ChiralVertex VertexImpl<typename fields::bar<fields::Fu>::type, fields::Fu, fields::VP>::evaluate(
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

ChiralVertex VertexImpl<typename fields::bar<fields::Fu>::type, fields::Fu, fields::VZ>::evaluate(
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

ChiralVertex VertexImpl<typename fields::bar<fields::Fv>::type, fields::Fe, fields::Hp>::evaluate(
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

ChiralVertex VertexImpl<typename fields::bar<fields::Fv>::type, fields::Fe, fields::VWp>::evaluate(
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

ChiralVertex VertexImpl<typename fields::bar<fields::Fv>::type, fields::Fv, fields::VZ>::evaluate(
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

ChiralVertex VertexImpl<typename fields::bar<fields::Fv>::type, fields::Hp, fields::Fe>::evaluate(
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

MomentumVertex VertexImpl<typename fields::bar<fields::gP>::type, fields::gWpC, fields::VWp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = g2*Sin(ThetaW);

   return {result, 1};
}

MomentumVertex VertexImpl<typename fields::bar<fields::gP>::type, fields::gWp, typename fields::conj<fields::VWp>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -(g2*Sin(ThetaW));

   return {result, 1};
}

ScalarVertex VertexImpl<typename fields::bar<fields::gWpC>::type, fields::gP, typename fields::conj<fields::Hp>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto v = MODELPARAMETER(v);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.25*g2*v*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW));

   return {result};
}

MomentumVertex VertexImpl<typename fields::bar<fields::gWpC>::type, fields::gP, typename fields::conj<fields::VWp>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = g2*Sin(ThetaW);

   return {result, 1};
}

ScalarVertex VertexImpl<typename fields::bar<fields::gWpC>::type, fields::gWpC, fields::Ah>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto v = MODELPARAMETER(v);

   const std::complex<double> result = std::complex<double>(0,0.25)*v*Sqr(g2);

   return {result};
}

ScalarVertex VertexImpl<typename fields::bar<fields::gWpC>::type, fields::gWpC, fields::hh>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto v = MODELPARAMETER(v);

   const std::complex<double> result = -0.25*v*Sqr(g2);

   return {result};
}

MomentumVertex VertexImpl<typename fields::bar<fields::gWpC>::type, fields::gWpC, fields::VP>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -(g2*Sin(ThetaW));

   return {result, 1};
}

MomentumVertex VertexImpl<typename fields::bar<fields::gWpC>::type, fields::gWpC, fields::VZ>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -(g2*Cos(ThetaW));

   return {result, 1};
}

ScalarVertex VertexImpl<typename fields::bar<fields::gWpC>::type, fields::gZ, typename fields::conj<fields::Hp>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto v = MODELPARAMETER(v);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.05*g2*v*(-5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW));

   return {result};
}

MomentumVertex VertexImpl<typename fields::bar<fields::gWpC>::type, fields::gZ, typename fields::conj<fields::VWp>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = g2*Cos(ThetaW);

   return {result, 1};
}

ScalarVertex VertexImpl<typename fields::bar<fields::gWp>::type, fields::gP, fields::Hp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto v = MODELPARAMETER(v);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.25*g2*v*(0.7745966692414834*g1*Cos(ThetaW) + g2*Sin(ThetaW));

   return {result};
}

MomentumVertex VertexImpl<typename fields::bar<fields::gWp>::type, fields::gP, fields::VWp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -(g2*Sin(ThetaW));

   return {result, 1};
}

ScalarVertex VertexImpl<typename fields::bar<fields::gWp>::type, fields::gWp, fields::Ah>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto v = MODELPARAMETER(v);

   const std::complex<double> result = std::complex<double>(0,-0.25)*v*Sqr(g2);

   return {result};
}

ScalarVertex VertexImpl<typename fields::bar<fields::gWp>::type, fields::gWp, fields::hh>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto v = MODELPARAMETER(v);

   const std::complex<double> result = -0.25*v*Sqr(g2);

   return {result};
}

MomentumVertex VertexImpl<typename fields::bar<fields::gWp>::type, fields::gWp, fields::VP>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = g2*Sin(ThetaW);

   return {result, 1};
}

MomentumVertex VertexImpl<typename fields::bar<fields::gWp>::type, fields::gWp, fields::VZ>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = g2*Cos(ThetaW);

   return {result, 1};
}

ScalarVertex VertexImpl<typename fields::bar<fields::gWp>::type, fields::gZ, fields::Hp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto v = MODELPARAMETER(v);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.05*g2*v*(-5*g2*Cos(ThetaW) + 3.872983346207417*g1*Sin(ThetaW));

   return {result};
}

MomentumVertex VertexImpl<typename fields::bar<fields::gWp>::type, fields::gZ, fields::VWp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -(g2*Cos(ThetaW));

   return {result, 1};
}

ScalarVertex VertexImpl<typename fields::bar<fields::gZ>::type, fields::gWpC, fields::Hp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto v = MODELPARAMETER(v);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.25*g2*v*(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW));

   return {result};
}

MomentumVertex VertexImpl<typename fields::bar<fields::gZ>::type, fields::gWpC, fields::VWp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = g2*Cos(ThetaW);

   return {result, 1};
}

ScalarVertex VertexImpl<typename fields::bar<fields::gZ>::type, fields::gWp, typename fields::conj<fields::Hp>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto v = MODELPARAMETER(v);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.25*g2*v*(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW));

   return {result};
}

MomentumVertex VertexImpl<typename fields::bar<fields::gZ>::type, fields::gWp, typename fields::conj<fields::VWp>::type>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -(g2*Cos(ThetaW));

   return {result, 1};
}

ScalarVertex VertexImpl<typename fields::bar<fields::gZ>::type, fields::gZ, fields::hh>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto v = MODELPARAMETER(v);
   const auto g2 = MODELPARAMETER(g2);
   const auto g1 = MODELPARAMETER(g1);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.25*v*Sqr(g2*Cos(ThetaW) + 0.7745966692414834*g1*Sin(ThetaW));

   return {result};
}

InverseMetricVertex VertexImpl<typename fields::conj<fields::Hp>::type, fields::VP, fields::VWp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto v = MODELPARAMETER(v);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = 0.3872983346207417*g1*g2*v*Cos(ThetaW);

   return {result};
}

InverseMetricVertex VertexImpl<typename fields::conj<fields::Hp>::type, fields::VWp, fields::VZ>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);
   const auto v = MODELPARAMETER(v);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -0.3872983346207417*g1*g2*v*Sin(ThetaW);

   return {result};
}

QuadrupleVectorVertex VertexImpl<typename fields::conj<fields::VWp>::type, fields::VP, fields::VP, fields::VWp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> part1 = Sqr(g2)*Sqr(Sin(ThetaW));

   const std::complex<double> part2 = Sqr(g2)*Sqr(Sin(ThetaW));

   const std::complex<double> part3 = -2*Sqr(g2)*Sqr(Sin(ThetaW));

   return {part1, part2, part3};
}

QuadrupleVectorVertex VertexImpl<typename fields::conj<fields::VWp>::type, fields::VP, fields::VWp, fields::VZ>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> part1 = Cos(ThetaW)*Sin(ThetaW)*Sqr(g2);

   const std::complex<double> part2 = -(Sin(2*ThetaW)*Sqr(g2));

   const std::complex<double> part3 = Cos(ThetaW)*Sin(ThetaW)*Sqr(g2);

   return {part1, part2, part3};
}

TripleVectorVertex VertexImpl<typename fields::conj<fields::VWp>::type, fields::VP, fields::VWp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = -(g2*Sin(ThetaW));

   return {result, TripleVectorVertex::odd_permutation{}};
}

QuadrupleVectorVertex VertexImpl<typename fields::conj<fields::VWp>::type, fields::VWp, fields::VZ, fields::VZ>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> part1 = -2*Sqr(g2)*Sqr(Cos(ThetaW));

   const std::complex<double> part2 = Sqr(g2)*Sqr(Cos(ThetaW));

   const std::complex<double> part3 = Sqr(g2)*Sqr(Cos(ThetaW));

   return {part1, part2, part3};
}

TripleVectorVertex VertexImpl<typename fields::conj<fields::VWp>::type, fields::VWp, fields::VZ>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);
   const auto ThetaW = DERIVEDPARAMETER(ThetaW);

   const std::complex<double> result = g2*Cos(ThetaW);

   return {result, TripleVectorVertex::odd_permutation{}};
}

QuadrupleVectorVertex VertexImpl<typename fields::conj<fields::VWp>::type, typename fields::conj<fields::VWp>::type, fields::VWp, fields::VWp>::evaluate(
   const std::array<int, 0>& indices, const context_base& context)
{
   const auto g2 = MODELPARAMETER(g2);

   const std::complex<double> part1 = 2*Sqr(g2);

   const std::complex<double> part2 = -Sqr(g2);

   const std::complex<double> part3 = -Sqr(g2);

   return {part1, part2, part3};
}

} // namespace detail
} // namespace SM_cxx_diagrams
} // namespace flexiblesusy
