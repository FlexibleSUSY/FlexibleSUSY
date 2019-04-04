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

#define BOOST_TEST_MODULE Test BTOSGAMMA
#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include "random_MSSM_dataset.hpp"
#include "test_complex_equality.hpp"
#include "concatenate.hpp"
#include "numerics.h"

#include "cxx_qft/MSSM_qft.hpp"
#include "MSSM_b_to_s_gamma.hpp"

using namespace flexiblesusy;
using namespace MSSM_cxx_diagrams;
using namespace MSSM_cxx_diagrams::fields;
using namespace boost_test_tools;

#define INPUTPARAMETER(p) context.model.get_input().p
#define MODELPARAMETER(p) context.model.get_##p()
#define DERIVEDPARAMETER(p) context.model.p()
#define PHASE(p) context.model.get_##p()

inline double LoopF1 (const double& x)
{
    if (is_zero(1-x)) {
        return -5./12;
    } else if (is_zero(x)) {
        return -7./6;
    } else {
        return (-7+5*x+8*x*x)/(6*pow(1-x, 3)) - x*(2-3*x)*Log(x)/pow(1-x, 4);
    }
}

inline double LoopF2 (const double& x)
{
    if (is_zero(1-x)) {
        return -7./6;
    }
    else if (is_zero(x)) {
        return 0;
    } else {
        return x*(3-5*x)/(2*Sqr(1-x)) + x*(2-3*x)*Log(x)/pow(1-x, 3);
    }
}

inline double LoopF3 (const double& x)
{
    if (is_zero(1-x)) {
        return 1./12;
    } else if (is_zero(x)) {
        1./3;
    } else {
        return (2+5*x-x*x)/(6*pow(1-x, 3)) + x*Log(x)/pow(1-x, 4);
    }
}

inline double LoopF4 (const double& x)
{
    if (is_zero(1-x)) {
        return 1./6;
    } else if (is_zero(x)) {
        return 1./2;
    } else {
        return (1+x)/(2*pow(1-x, 2)) + x*Log(x)/pow(1-x, 3);
    }
}

inline double LoopG1 (const double& x)
{
    return -3*LoopF3(x);
}

inline double LoopG2 (const double& x)
{
    return 3*(LoopF4(x)-1./2);
}

inline double LoopG3 (const double& x)
{
    return -3*LoopF3(x);
}
inline double LoopG4 (const double& x)
{
    return -3*LoopF4(x);
}

inline double LoopG5 (const double& x)
{
    if (is_zero(1-x)) {
        return 5./16;
    } else if (is_zero(x)) {
        return 11./16;
    } else {
        return 1./8 * ((11-40*x-19*x*x)/(2*pow(1-x, 3)) + 3*x*(1-9*x)/pow(1-x, 4)*Log(x));
    }
}

inline double LoopG6 (const double& x)
{
    if (is_zero(1-x)) {
        return 19./16;
    } else if (is_zero(x)) {
        return 15./8;
    } else {
        return 3./8 * ((5-13*x)/(Sqr(1-x)) + x*(1-9*x)/pow(1-x, 3)*Log(x));
    }
}

namespace
{
    static constexpr int number_of_random_samples = 20;
}


BOOST_AUTO_TEST_SUITE(b_to_s_gamma_test_suite,
	* boost::unit_test::tolerance(1.0e-2) )

BOOST_DATA_TEST_CASE( test_b_to_s_gamma,
	random_MSSM_dataset( number_of_random_samples ),
  index, FS_TEST_MSSM_PARAMETER_SEQUENCE )
{
	const auto input_parameters = wrap_input_parameters(
		FS_TEST_MSSM_PARAMETER_SEQUENCE );
	const auto model = calculate_spectrum( input_parameters );

	const auto problems = model.get_problems();
	BOOST_WARN_MESSAGE( !problems.have_problem(),
		"There was an error calculating the spectrum: " <<
		problems << "Skipping data point..." );

    context_base context { model };
    softsusy::QedQcd qedqcd;

	if( problems.have_problem() == false )
	{
        const auto complex_CKM = qedqcd.get_complex_ckm();
        std::complex<double> Vtb (complex_CKM(2, 2));
        std::complex<double> Vts;
        // avoid division by zero
        if (complex_CKM(2, 1) == std::complex<double>(0, 0)) {
            Vts = std::complex<double>(-0.0404063888996, -0.00072107415577);
        } else {
            Vts = complex_CKM(2, 1);
        }
        const auto g1 = MODELPARAMETER(g1);
        const auto g2 = MODELPARAMETER(g2);
        const auto g3 = MODELPARAMETER(g3);
        const auto mW = context.mass<VWm>({});
        const auto mb = context.mass<Fd>({ 2 });
        const auto mt = context.mass<Fu>({ 2 });
        const auto mHpm = context.mass<Hpm>({ 1 });
        const auto mGlu = context.mass<Glu> ({});
        std::array<double, 2> mCha = { context.mass<Cha> ({0}),
                                       context.mass<Cha> ({1}) };
        std::array<double, 4> mChi = { context.mass<Chi> ({0}),
                                       context.mass<Chi> ({1}),
                                       context.mass<Chi> ({2}),
                                       context.mass<Chi> ({3}) };
        std::array<double, 6> mSd = { context.mass<Sd>({0}), context.mass<Sd>({1}),
            context.mass<Sd>({2}), context.mass<Sd>({3}), context.mass<Sd>({4}), context.mass<Sd>({5}) };
        std::array<double, 6> mSu = { context.mass<fields::Su>({0}), context.mass<fields::Su>({1}),
            context.mass<Su>({2}), context.mass<Su>({3}), context.mass<Su>({4}), context.mass<Su>({5}) };

        const auto ratioH = Sqr(mt/mHpm);
        const auto cotbeta = MODELPARAMETER(vd) / MODELPARAMETER(vu);

        const auto GD = MODELPARAMETER(ZD);
        const auto ZDL = MODELPARAMETER(ZDL);
        const auto ZDR = MODELPARAMETER(ZDR);
        const auto GDL = GD.block<6, 3>(0, 0);
        const auto GDR = GD.block<6, 3>(0, 3);
        const auto GU = MODELPARAMETER(ZU);
        const auto GUL = GU.block<6, 3>(0, 0);
        const auto GUR = GU.block<6, 3>(0, 3);
        const auto Yu = MODELPARAMETER(Yu);
        const auto Yd = MODELPARAMETER(Yd);
        const auto UP = MODELPARAMETER(UP);
        const auto UM = MODELPARAMETER(UM);
        const auto ZN = MODELPARAMETER(ZN);
        const auto normC7 = -16*Sqr(Pi)*Sqr(mW)/(unit_charge(context)*Sqr(g2)*Vtb*Conj(Vts));

        // chargino vertices
        const auto XUL1 = -(g2*Conj(UP(0,0))*GUL*ZDL.adjoint() - Conj(UP(0,1))*(GUR*Yu)*ZDL.adjoint());
        const auto XUL2 = -(g2*Conj(UP(1,0))*GUL*ZDL.adjoint() - Conj(UP(1,1))*(GUR*Yu)*ZDL.adjoint());
        const auto XUR1 = GUL * ZDR * Yd.conjugate() * UM(0,1);
        const auto XUR2 = GUL * ZDR * Yd.conjugate() * UM(1,1);

        //neutralino vertices
        const auto XDL1 = -0.18257418583505536*g1*Conj(ZN(0,0))*GDL*ZDL.adjoint() + 0.7071067811865475*g2*Conj(ZN(0,1))*GDL*ZDL.adjoint() - Conj(ZN(0,2))*(GDL*Yd)*ZDL.adjoint();
        const auto XDL2 = -0.18257418583505536*g1*Conj(ZN(1,0))*GDL*ZDL.adjoint() + 0.7071067811865475*g2*Conj(ZN(1,1))*GDL*ZDL.adjoint() - Conj(ZN(1,2))*(GDL*Yd)*ZDL.adjoint();
        const auto XDL3 = -0.18257418583505536*g1*Conj(ZN(2,0))*GDL*ZDL.adjoint() + 0.7071067811865475*g2*Conj(ZN(2,1))*GDL*ZDL.adjoint() - Conj(ZN(2,2))*(GDL*Yd)*ZDL.adjoint();
        const auto XDL4 = -0.18257418583505536*g1*Conj(ZN(3,0))*GDL*ZDL.adjoint() + 0.7071067811865475*g2*Conj(ZN(3,1))*GDL*ZDL.adjoint() - Conj(ZN(3,2))*(GDL*Yd)*ZDL.adjoint();
        const auto XDR1 = -0.3651483716701107*g1*GDR*ZDR.transpose()*ZN(0,0) - GDL * (ZDR * Yd.conjugate()) * ZN(0,2);
        const auto XDR2 = -0.3651483716701107*g1*GDR*ZDR.transpose()*ZN(1,0) - GDL * (ZDR * Yd.conjugate()) * ZN(1,2);
        const auto XDR3 = -0.3651483716701107*g1*GDR*ZDR.transpose()*ZN(2,0) - GDL * (ZDR * Yd.conjugate()) * ZN(2,2);
        const auto XDR4 = -0.3651483716701107*g1*GDR*ZDR.transpose()*ZN(3,0) - GDL * (ZDR * Yd.conjugate()) * ZN(3,2);

        // formula taken from https://arxiv.org/abs/hep-ph/0002089v2
        const auto C7Hpm = 1./12 * (ratioH*Sqr(cotbeta)*LoopF1(ratioH) + 2*LoopF2(ratioH));
        const auto C8Hpm = 1./12 * (ratioH*Sqr(cotbeta)*(-3)*LoopF3(ratioH) + 6*(LoopF4(ratioH) - 1./2.));

        const auto C7Cha = - Sqr(mW)/(6*g2*g2*Vtb*Conj(Vts)) * SUM(j1, 0, 5,
              1/Sqr(mCha.at(0)) * Conj(XUL1(j1, 1)) * XUL1(j1, 2) * LoopF1(Sqr(mSu.at(j1)/mCha.at(0)))
            + 1/Sqr(mCha.at(1)) * Conj(XUL2(j1, 1)) * XUL2(j1, 2) * LoopF1(Sqr(mSu.at(j1)/mCha.at(1)))
            - 1/Sqr(mCha.at(0)) * 2 * Conj(XUL1(j1, 1)) * XUR1(j1, 2) * mCha.at(0)/mb * LoopF2(Sqr(mSu.at(j1)/mCha.at(0)))
            - 1/Sqr(mCha.at(1)) * 2 * Conj(XUL2(j1, 1)) * XUR2(j1, 2) * mCha.at(1)/mb * LoopF2(Sqr(mSu.at(j1)/mCha.at(1))) );

        const auto C8Cha = - Sqr(mW)/(6*g2*g2*Vtb*Conj(Vts)) * SUM(j1, 0, 5,
              1/Sqr(mCha.at(0)) * Conj(XUL1(j1, 1)) * XUL1(j1, 2) * LoopG1(Sqr(mSu.at(j1)/mCha.at(0)))
            + 1/Sqr(mCha.at(1)) * Conj(XUL2(j1, 1)) * XUL2(j1, 2) * LoopG1(Sqr(mSu.at(j1)/mCha.at(1)))
            - 1/Sqr(mCha.at(0)) * 2 * Conj(XUL1(j1, 1)) * XUR1(j1, 2) * mCha.at(0)/mb * LoopG2(Sqr(mSu.at(j1)/mCha.at(0)))
            - 1/Sqr(mCha.at(1)) * 2 * Conj(XUL2(j1, 1)) * XUR2(j1, 2) * mCha.at(1)/mb * LoopG2(Sqr(mSu.at(j1)/mCha.at(1))) );

        const auto C7Chi = - Sqr(mW)/(6*g2*g2*Vtb*Conj(Vts)) * SUM(j1, 0, 5,
                  1/Sqr(mChi.at(0)) * Conj(XDL1(j1, 1)) * XDL1(j1, 2) * LoopF3(Sqr(mSd.at(j1)/context.mass<Chi>({0})))
                + 1/Sqr(mChi.at(1)) * Conj(XDL2(j1, 1)) * XDL2(j1, 2) * LoopF3(Sqr(mSd.at(j1)/context.mass<Chi>({1})))
                + 1/Sqr(mChi.at(2)) * Conj(XDL3(j1, 1)) * XDL3(j1, 2) * LoopF3(Sqr(mSd.at(j1)/context.mass<Chi>({2})))
                + 1/Sqr(mChi.at(3)) * Conj(XDL4(j1, 1)) * XDL4(j1, 2) * LoopF3(Sqr(mSd.at(j1)/context.mass<Chi>({3})))
                + 2/Sqr(mChi.at(0)) * Conj(XDL1(j1, 1)) * XDR1(j1, 2) * mChi.at(0)/mb * LoopF4(Sqr(mSd.at(j1)/context.mass<Chi>({0})))
                + 2/Sqr(mChi.at(1)) * Conj(XDL2(j1, 1)) * XDR2(j1, 2) * mChi.at(1)/mb * LoopF4(Sqr(mSd.at(j1)/context.mass<Chi>({1})))
                + 2/Sqr(mChi.at(2)) * Conj(XDL3(j1, 1)) * XDR3(j1, 2) * mChi.at(2)/mb * LoopF4(Sqr(mSd.at(j1)/context.mass<Chi>({2})))
                + 2/Sqr(mChi.at(3)) * Conj(XDL4(j1, 1)) * XDR4(j1, 2) * mChi.at(3)/mb * LoopF4(Sqr(mSd.at(j1)/context.mass<Chi>({3})))
                );

        const auto C8Chi = - Sqr(mW)/(6*g2*g2*Vtb*Conj(Vts)) * SUM(j1, 0, 5,
                  1/Sqr(mChi.at(0)) * Conj(XDL1(j1, 1)) * XDL1(j1, 2) * LoopG3(Sqr(mSd.at(j1)/context.mass<Chi>({0})))
                + 1/Sqr(mChi.at(1)) * Conj(XDL2(j1, 1)) * XDL2(j1, 2) * LoopG3(Sqr(mSd.at(j1)/context.mass<Chi>({1})))
                + 1/Sqr(mChi.at(2)) * Conj(XDL3(j1, 1)) * XDL3(j1, 2) * LoopG3(Sqr(mSd.at(j1)/context.mass<Chi>({2})))
                + 1/Sqr(mChi.at(3)) * Conj(XDL4(j1, 1)) * XDL4(j1, 2) * LoopG3(Sqr(mSd.at(j1)/context.mass<Chi>({3})))
                + 2/Sqr(mChi.at(0)) * Conj(XDL1(j1, 1)) * XDR1(j1, 2) * mChi.at(0)/mb * LoopG4(Sqr(mSd.at(j1)/context.mass<Chi>({0})))
                + 2/Sqr(mChi.at(1)) * Conj(XDL2(j1, 1)) * XDR2(j1, 2) * mChi.at(1)/mb * LoopG4(Sqr(mSd.at(j1)/context.mass<Chi>({1})))
                + 2/Sqr(mChi.at(2)) * Conj(XDL3(j1, 1)) * XDR3(j1, 2) * mChi.at(2)/mb * LoopG4(Sqr(mSd.at(j1)/context.mass<Chi>({2})))
                + 2/Sqr(mChi.at(3)) * Conj(XDL4(j1, 1)) * XDR4(j1, 2) * mChi.at(3)/mb * LoopG4(Sqr(mSd.at(j1)/context.mass<Chi>({3})))
                );

        const auto C7Glu = - 4 * g3 * g3 / (9*g2*g2*Vtb*Conj(Vts)) * Sqr(mW/mGlu) *
            SUM(j2, 0, 5, SUM(j1,0,2,Conj(GD(j2,j1))*ZDL(1,j1)) * SUM(j1,0,2,Conj(ZDL(2,j1))*GD(j2,j1)) * LoopF3(Sqr(mSd.at(j2)/context.mass<Glu>({})))
                - 2 * SUM(j1,0,2,Conj(GD(j2,j1))*ZDL(1,j1)) * SUM(j1,0,2,GD(j2,3 + j1)*ZDR(2,j1)) * context.mass<Glu>({})/mb *
                LoopF4(Sqr(mSd.at(j2)/context.mass<Glu>({}))));

        const auto C8Glu = - 4 * g3 * g3 / (9*g2*g2*Vtb*Conj(Vts)) * Sqr(mW/mGlu) *
            SUM(j2, 0, 5, SUM(j1,0,2,Conj(GD(j2,j1))*ZDL(1,j1)) * SUM(j1,0,2,Conj(ZDL(2,j1))*GD(j2,j1)) * LoopG5(Sqr(mSd.at(j2)/context.mass<Glu>({})))
                - 2 * SUM(j1,0,2,Conj(GD(j2,j1))*ZDL(1,j1)) * SUM(j1,0,2,GD(j2,3 + j1)*ZDR(2,j1)) * context.mass<Glu>({})/mb *
                LoopG6(Sqr(mSd.at(j2)/context.mass<Glu>({}))));


        //TODO: add Hpm contribution if CKM matrix is non-diagonal
        const auto C7NP = 0*C7Hpm + C7Cha + C7Chi + C7Glu;
        const auto C8NP = 0*C8Hpm + C8Cha + C8Chi + C8Glu;

        const auto C7_C8 = MSSM_b_to_s_gamma::calculate_b_to_s_gamma(model, qedqcd);

        TEST_COMPLEX_EQUALITY(C7_C8[0], C7NP);
        TEST_COMPLEX_EQUALITY(C7_C8[2], C8NP);

	} else
		BOOST_TEST( true );
}

BOOST_AUTO_TEST_SUITE_END()
