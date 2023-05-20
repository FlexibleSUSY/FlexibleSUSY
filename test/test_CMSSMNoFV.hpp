
#ifndef TEST_CMSSMNOFV_H
#define TEST_CMSSMNOFV_H

#include "test.hpp"
#include "wrappers.hpp"
#include "ew_input.hpp"

using namespace flexiblesusy;

template <class T>
void ensure_tree_level_ewsb(T& m)
{
   // ensure that the EWSB eqs. are satisfied (Drees p.222) for an
   // CMSSM-like model
   const double vu = m.get_vu();
   const double vd = m.get_vd();
   const double gY = m.get_g1() * Sqrt(0.6);
   const double g2 = m.get_g2();
   const double Mu = m.get_Mu();
   const double BMu = m.get_BMu();
   const double mHd2 = BMu*vu/vd - (Sqr(gY) + Sqr(g2))*(Sqr(vd) - Sqr(vu))/8. - Sqr(Mu);
   const double mHu2 = BMu*vd/vu + (Sqr(gY) + Sqr(g2))*(Sqr(vd) - Sqr(vu))/8. - Sqr(Mu);
   m.set_mHd2(mHd2);
   m.set_mHu2(mHu2);
}

template <class T>
void ensure_n_loop_ewsb(T& m, int loop_level)
{
   ensure_tree_level_ewsb(m);

   const double precision = m.get_ewsb_iteration_precision();
   m.set_ewsb_loop_order(loop_level);
   m.solve_ewsb();

   if (loop_level == 1) {
      TEST_CLOSE(m.get_ewsb_eq_hh_1() - m.tadpole_hh_1loop(0).real(), 0.0, precision);
      TEST_CLOSE(m.get_ewsb_eq_hh_2() - m.tadpole_hh_1loop(1).real(), 0.0, precision);
   }
}

template <class T1, class T2, class TInput>
void setup_CMSSM_models(T1& m1, T2& m2, const TInput& input)
{
   const double ALPHASMZ = 0.1176;
   const double ALPHAMZ = 1.0 / 127.918;
   const double sinthWsq = 0.23122;
   const double alpha1 = 5 * ALPHAMZ / (3 * (1 - sinthWsq));
   const double alpha2 = ALPHAMZ / sinthWsq;
   const double g1 = Sqrt(4 * Pi * alpha1);
   const double g2 = Sqrt(4 * Pi * alpha2);
   const double g3 = Sqrt(4 * Pi * ALPHASMZ);
   const double tanBeta = input.TanBeta;
   const double sinBeta = sin(atan(tanBeta));
   const double cosBeta = cos(atan(tanBeta));
   const double M12 = input.m12;
   const double m0 = input.m0;
   const double a0 = input.Azero;
   const double root2 = Sqrt(2.0);
   const double vev = 246.0;
   const double vu = vev * sinBeta;
   const double vd = vev * cosBeta;
   const double susyMu = input.SignMu * 120.0;
   const double BMu = Sqr(2.0 * susyMu);
   const double scale = Electroweak_constants::MZ;

   Eigen::Matrix<double,3,3> Yu(Eigen::Matrix<double,3,3>::Zero()),
      Yd(Eigen::Matrix<double,3,3>::Zero()),
      Ye(Eigen::Matrix<double,3,3>::Zero()),
      mm0(Eigen::Matrix<double,3,3>::Zero());
   Yu(2,2) = 165.0   * root2 / (vev * sinBeta);
   Yd(2,2) = 2.9     * root2 / (vev * cosBeta);
   Ye(2,2) = 1.77699 * root2 / (vev * cosBeta);
   mm0 = Sqr(m0) * Eigen::Matrix<double,3,3>::Identity();

   m1.set_scale(scale);
   m1.set_loops(2);
   m1.set_g1(g1);
   m1.set_g2(g2);
   m1.set_g3(g3);
   m1.set_Yu(Yu);
   m1.set_Yd(Yd);
   m1.set_Ye(Ye);
   m1.set_MassB(M12);
   m1.set_MassG(M12);
   m1.set_MassWB(M12);
   m1.set_mq2(mm0);
   m1.set_ml2(mm0);
   m1.set_md2(mm0);
   m1.set_mu2(mm0);
   m1.set_me2(mm0);
   m1.set_mHd2(Sqr(m0));
   m1.set_mHu2(Sqr(m0));
   m1.set_TYu(a0 * Yu);
   m1.set_TYd(a0 * Yd);
   m1.set_TYe(a0 * Ye);
   m1.set_Mu(susyMu);
   m1.set_BMu(BMu);
   m1.set_vu(vu);
   m1.set_vd(vd);

   m2.set_scale(scale);
   m2.set_loops(2);
   m2.set_g1(g1);
   m2.set_g2(g2);
   m2.set_g3(g3);
   m2.set_Yu(Yu);
   m2.set_Yd(Yd);
   m2.set_Ye(Ye);
   m2.set_MassB(M12);
   m2.set_MassG(M12);
   m2.set_MassWB(M12);
   m2.set_mq2(mm0);
   m2.set_ml2(mm0);
   m2.set_md2(mm0);
   m2.set_mu2(mm0);
   m2.set_me2(mm0);
   m2.set_mHd2(Sqr(m0));
   m2.set_mHu2(Sqr(m0));
   m2.set_TYu(a0 * Yu);
   m2.set_TYd(a0 * Yd);
   m2.set_TYe(a0 * Ye);
   m2.set_Mu(susyMu);
   m2.set_BMu(BMu);
   m2.set_vu(vu);
   m2.set_vd(vd);

   ensure_tree_level_ewsb(m1);
   m1.calculate_DRbar_masses();

   ensure_tree_level_ewsb(m2);
   m2.calculate_DRbar_masses();
}

template <class T1, class TInput>
void setup_CMSSM_const(T1& m, const TInput& input)
{
   const double ALPHASMZ = 0.1176;
   const double ALPHAMZ = 1.0 / 127.918;
   const double sinthWsq = 0.23122;
   const double alpha1 = 5 * ALPHAMZ / (3 * (1 - sinthWsq));
   const double alpha2 = ALPHAMZ / sinthWsq;
   const double g1 = Sqrt(4 * Pi * alpha1);
   const double g2 = Sqrt(4 * Pi * alpha2);
   const double g3 = Sqrt(4 * Pi * ALPHASMZ);
   const double tanBeta = input.TanBeta;
   const double sinBeta = sin(atan(tanBeta));
   const double cosBeta = cos(atan(tanBeta));
   const double M12 = input.m12;
   const double m0 = input.m0;
   const double a0 = input.Azero;
   const double root2 = Sqrt(2.0);
   const double vev = 246.0;
   const double vu = vev * sinBeta;
   const double vd = vev * cosBeta;
   const double susyMu = input.SignMu * 120.0;
   const double BMu = Sqr(2.0 * susyMu);
   const double scale = Electroweak_constants::MZ;

   Eigen::Matrix<double,3,3> Yu(Eigen::Matrix<double,3,3>::Zero()),
      Yd(Eigen::Matrix<double,3,3>::Zero()),
      Ye(Eigen::Matrix<double,3,3>::Zero()),
      mm0(Eigen::Matrix<double,3,3>::Zero());
   Yu(2,2) = 165.0   * root2 / (vev * sinBeta);
   Yd(2,2) = 2.9     * root2 / (vev * cosBeta);
   Ye(2,2) = 1.77699 * root2 / (vev * cosBeta);
   mm0 = Sqr(m0) * Eigen::Matrix<double,3,3>::Identity();

   m.set_scale(scale);
   m.set_loops(2);
   m.set_g1(g1);
   m.set_g2(g2);
   m.set_g3(g3);
   m.set_Yu(Yu);
   m.set_Yd(Yd);
   m.set_Ye(Ye);
   m.set_MassB(M12);
   m.set_MassG(M12);
   m.set_MassWB(M12);
   m.set_mq2(mm0);
   m.set_ml2(mm0);
   m.set_md2(mm0);
   m.set_mu2(mm0);
   m.set_me2(mm0);
   m.set_mHd2(Sqr(m0));
   m.set_mHu2(Sqr(m0));
   m.set_TYu(a0 * Yu);
   m.set_TYd(a0 * Yd);
   m.set_TYe(a0 * Ye);
   m.set_Mu(susyMu);
   m.set_BMu(BMu);
   m.set_vu(vu);
   m.set_vd(vd);

   ensure_tree_level_ewsb(m);
   m.calculate_DRbar_masses();
}

#endif
