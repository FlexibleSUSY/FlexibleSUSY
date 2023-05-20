
#include "MSSM_two_scale_model.hpp"
#include "ew_input.hpp"
#include "logger.hpp"
#include "wrappers.hpp"

using namespace flexiblesusy;

void setup(MSSM<Two_scale>& mssm)
{
   Eigen::Matrix<double,3,3> Yu;
   Eigen::Matrix<double,3,3> Yd;
   Eigen::Matrix<double,3,3> Ye;
   double Mu;
   double g1;
   double g2;
   double g3;
   double vd;
   double vu;
   Eigen::Matrix<double,3,3> TYu;
   Eigen::Matrix<double,3,3> TYd;
   Eigen::Matrix<double,3,3> TYe;
   double BMu;
   Eigen::Matrix<double,3,3> mq2;
   Eigen::Matrix<double,3,3> ml2;
   double mHd2;
   double mHu2;
   Eigen::Matrix<double,3,3> md2;
   Eigen::Matrix<double,3,3> mu2;
   Eigen::Matrix<double,3,3> me2;
   double MassB;
   double MassWB;
   double MassG;

   // susy parameters
   Yu << 1.26136e-05, 0, 0,
                   0, 0.00667469, 0,
                   0, 0, 0.857849;

   Yd << 0.000242026, 0, 0,
                   0, 0.00529911, 0,
                   0, 0, 0.193602;

   Ye << 2.84161e-05, 0, 0,
                   0, 0.00587557, 0,
                   0, 0, 0.10199;

   Mu = 627.164;
   g1 = 0.468171;
   g2 = 0.642353;
   g3 = 1.06459;
   vd = 25.0944;
   vu = 242.968;

   // soft parameters
   TYu << -0.0144387, 0, 0,
                   0, -7.64037, 0,
                   0, 0, -759.305;

   TYd << -0.336207, 0, 0,
                  0, -7.36109, 0,
                  0, 0, -250.124;

   TYe << -0.00825134, 0, 0,
                    0, -1.70609, 0,
                    0, 0, -29.4466;

   BMu = 52140.8;

   mq2 << 1.03883e+06, 0, 0,
                    0, 1.03881e+06, 0,
                    0, 0, 879135;

   ml2 << 124856, 0, 0,
               0, 124853, 0,
               0, 0, 124142;

   mHd2 = 92436.9;
   mHu2 = -380337;

   md2 << 954454, 0, 0,
               0, 954439, 0,
               0, 0, 934727;

   mu2 << 963422, 0, 0,
               0, 963400, 0,
               0, 0, 656621;

   me2 << 49215.8, 0, 0,
                0, 49210.9, 0,
                0, 0, 47759.2;

   MassB = 210.328;
   MassWB = 389.189;
   MassG = 1114.45;

   // set parameters
   mssm.set_scale(Electroweak_constants::MZ);
   mssm.set_Yu(Yu);
   mssm.set_Yd(Yd);
   mssm.set_Ye(Ye);
   mssm.set_Mu(Mu);
   mssm.set_g1(g1);
   mssm.set_g2(g2);
   mssm.set_g3(g3);
   mssm.set_vd(vd);
   mssm.set_vu(vu);
   mssm.set_TYu(TYu);
   mssm.set_TYd(TYd);
   mssm.set_TYe(TYe);
   mssm.set_BMu(BMu);
   mssm.set_mq2(mq2);
   mssm.set_ml2(ml2);
   mssm.set_mHd2(mHd2);
   mssm.set_mHu2(mHu2);
   mssm.set_md2(md2);
   mssm.set_mu2(mu2);
   mssm.set_me2(me2);
   mssm.set_MassB(MassB);
   mssm.set_MassWB(MassWB);
   mssm.set_MassG(MassG);
}

void self_energy_example()
{
   MSSM<Two_scale> mssm;

   setup(mssm);
   mssm.calculate_DRbar_masses();
   // INFO(mssm);

   double p = Electroweak_constants::MZ;
   double self_energy_VZ = Re(mssm.self_energy_VZ_1loop(p));

   double vertex_barFe_Fe_VZ_PR = mssm.CpbarFeFeVZPR(0, 0);
   double vertex_barFe_Fe_VZ_PL = mssm.CpbarFeFeVZPL(0, 0);

   INFO("Sigma_VZ(MZ) = " << self_energy_VZ);
   INFO("Vertex(bar(Fe),Fe,VZ)_PR = " << vertex_barFe_Fe_VZ_PR);
   INFO("Vertex(bar(Fe),Fe,VZ)_PL = " << vertex_barFe_Fe_VZ_PL);
}

int main()
{
   INFO("=============================");
   INFO("running self_energy_example()");
   INFO("=============================");

   self_energy_example();

   return 0;
}
