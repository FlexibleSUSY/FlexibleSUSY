// ====================================================================
// This file is part of GM2Calc.
//
// GM2Calc is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// GM2Calc is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with GM2Calc.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

/**
 * @file MSSMNoFV_onshell_mass_eigenstates.hpp
 *
 * @brief contains class for MSSMNoFV with routines needed to solve EWSB
 *        and determine the pole masses and mixings
 *
 * This file was generated at Wed 22 Jul 2015 18:14:31 with FlexibleSUSY
 * 1.2.1 (git commit: v1.2.1-49-gfc4c300) and SARAH 4.5.8 .
 */

#ifndef GM2_MSSMNoFV_ONSHELL_MASS_EIGENSTATES_HPP
#define GM2_MSSMNoFV_ONSHELL_MASS_EIGENSTATES_HPP

#include "MSSMNoFV_onshell_soft_parameters.hpp"
#include "MSSMNoFV_onshell_physical.hpp"
#include "MSSMNoFV_onshell_problems.hpp"

#include <iosfwd>

#include <Eigen/Core>

namespace gm2calc {

/**
 * @class MSSMNoFV_onshell_mass_eigenstates
 * @brief model class with routines determine masses, mixings and EWSB
 */
class MSSMNoFV_onshell_mass_eigenstates : public MSSMNoFV_onshell_soft_parameters {
public:
   virtual ~MSSMNoFV_onshell_mass_eigenstates() {}

   void print(std::ostream&) const override;

   void calculate_DRbar_masses();
   void copy_DRbar_masses_to_pole_masses();
   void do_force_output(bool);
   bool do_force_output() const;
   void reorder_DRbar_masses();
   void reorder_pole_masses();
   void set_physical(const MSSMNoFV_onshell_physical&);
   const MSSMNoFV_onshell_physical& get_physical() const;
   MSSMNoFV_onshell_physical& get_physical();
   const MSSMNoFV_onshell_problems& get_problems() const;
   MSSMNoFV_onshell_problems& get_problems();
   int solve_ewsb_tree_level();
   int solve_ewsb();

   double get_MVG() const { return MVG; }
   double get_MGlu() const { return MGlu; }
   double get_MVP() const { return MVP; }
   double get_MVZ() const { return MVZ; }
   double get_MFd() const { return MFd; }
   double get_MFs() const { return MFs; }
   double get_MFb() const { return MFb; }
   double get_MFu() const { return MFu; }
   double get_MFc() const { return MFc; }
   double get_MFt() const { return MFt; }
   double get_MFve() const { return MFve; }
   double get_MFvm() const { return MFvm; }
   double get_MFvt() const { return MFvt; }
   double get_MFe() const { return MFe; }
   double get_MFm() const { return MFm; }
   double get_MFtau() const { return MFtau; }
   double get_MSveL() const { return MSveL; }
   double get_MSvmL() const { return MSvmL; }
   double get_MSvtL() const { return MSvtL; }
   const Eigen::Array<double,2,1>& get_MSd() const { return MSd; }
   double get_MSd(int i) const { return MSd(i); }
   const Eigen::Array<double,2,1>& get_MSu() const { return MSu; }
   double get_MSu(int i) const { return MSu(i); }
   const Eigen::Array<double,2,1>& get_MSe() const { return MSe; }
   double get_MSe(int i) const { return MSe(i); }
   const Eigen::Array<double,2,1>& get_MSm() const { return MSm; }
   double get_MSm(int i) const { return MSm(i); }
   const Eigen::Array<double,2,1>& get_MStau() const { return MStau; }
   double get_MStau(int i) const { return MStau(i); }
   const Eigen::Array<double,2,1>& get_MSs() const { return MSs; }
   double get_MSs(int i) const { return MSs(i); }
   const Eigen::Array<double,2,1>& get_MSc() const { return MSc; }
   double get_MSc(int i) const { return MSc(i); }
   const Eigen::Array<double,2,1>& get_MSb() const { return MSb; }
   double get_MSb(int i) const { return MSb(i); }
   const Eigen::Array<double,2,1>& get_MSt() const { return MSt; }
   double get_MSt(int i) const { return MSt(i); }
   const Eigen::Array<double,2,1>& get_Mhh() const { return Mhh; }
   double get_Mhh(int i) const { return Mhh(i); }
   const Eigen::Array<double,2,1>& get_MAh() const { return MAh; }
   double get_MAh(int i) const { return MAh(i); }
   const Eigen::Array<double,2,1>& get_MHpm() const { return MHpm; }
   double get_MHpm(int i) const { return MHpm(i); }
   const Eigen::Array<double,4,1>& get_MChi() const { return MChi; }
   double get_MChi(int i) const { return MChi(i); }
   const Eigen::Array<double,2,1>& get_MCha() const { return MCha; }
   double get_MCha(int i) const { return MCha(i); }
   double get_MVWm() const { return MVWm; }


   Eigen::Array<double,1,1> get_MChargedHiggs() const;
   Eigen::Array<double,1,1> get_MPseudoscalarHiggs() const;

   const Eigen::Matrix<double,2,2>& get_ZD() const { return ZD; }
   double get_ZD(int i, int k) const { return ZD(i,k); }
   const Eigen::Matrix<double,2,2>& get_ZU() const { return ZU; }
   double get_ZU(int i, int k) const { return ZU(i,k); }
   const Eigen::Matrix<double,2,2>& get_ZE() const { return ZE; }
   double get_ZE(int i, int k) const { return ZE(i,k); }
   const Eigen::Matrix<double,2,2>& get_ZM() const { return ZM; }
   double get_ZM(int i, int k) const { return ZM(i,k); }
   const Eigen::Matrix<double,2,2>& get_ZTau() const { return ZTau; }
   double get_ZTau(int i, int k) const { return ZTau(i,k); }
   const Eigen::Matrix<double,2,2>& get_ZS() const { return ZS; }
   double get_ZS(int i, int k) const { return ZS(i,k); }
   const Eigen::Matrix<double,2,2>& get_ZC() const { return ZC; }
   double get_ZC(int i, int k) const { return ZC(i,k); }
   const Eigen::Matrix<double,2,2>& get_ZB() const { return ZB; }
   double get_ZB(int i, int k) const { return ZB(i,k); }
   const Eigen::Matrix<double,2,2>& get_ZT() const { return ZT; }
   double get_ZT(int i, int k) const { return ZT(i,k); }
   const Eigen::Matrix<double,2,2>& get_ZH() const { return ZH; }
   double get_ZH(int i, int k) const { return ZH(i,k); }
   const Eigen::Matrix<double,2,2>& get_ZA() const { return ZA; }
   double get_ZA(int i, int k) const { return ZA(i,k); }
   const Eigen::Matrix<double,2,2>& get_ZP() const { return ZP; }
   double get_ZP(int i, int k) const { return ZP(i,k); }
   const Eigen::Matrix<std::complex<double>,4,4>& get_ZN() const { return ZN; }
   const std::complex<double>& get_ZN(int i, int k) const { return ZN(i,k); }
   const Eigen::Matrix<std::complex<double>,2,2>& get_UM() const { return UM; }
   const std::complex<double>& get_UM(int i, int k) const { return UM(i,k); }
   const Eigen::Matrix<std::complex<double>,2,2>& get_UP() const { return UP; }
   const std::complex<double>& get_UP(int i, int k) const { return UP(i,k); }

   void set_PhaseGlu(std::complex<double> PhaseGlu_) { PhaseGlu = PhaseGlu_; }
   std::complex<double> get_PhaseGlu() const { return PhaseGlu; }

   double get_mass_matrix_VG() const;
   void calculate_MVG();
   double get_mass_matrix_Glu() const;
   void calculate_MGlu();
   double get_mass_matrix_VP() const;
   void calculate_MVP();
   double get_mass_matrix_VZ() const;
   void calculate_MVZ();
   double get_mass_matrix_Fd() const;
   void calculate_MFd();
   double get_mass_matrix_Fs() const;
   void calculate_MFs();
   double get_mass_matrix_Fb() const;
   void calculate_MFb();
   double get_mass_matrix_Fu() const;
   void calculate_MFu();
   double get_mass_matrix_Fc() const;
   void calculate_MFc();
   double get_mass_matrix_Ft() const;
   void calculate_MFt();
   double get_mass_matrix_Fve() const;
   void calculate_MFve();
   double get_mass_matrix_Fvm() const;
   void calculate_MFvm();
   double get_mass_matrix_Fvt() const;
   void calculate_MFvt();
   double get_mass_matrix_Fe() const;
   void calculate_MFe();
   double get_mass_matrix_Fm() const;
   void calculate_MFm();
   double get_mass_matrix_Ftau() const;
   void calculate_MFtau();
   double get_mass_matrix_SveL() const;
   void calculate_MSveL();
   double get_mass_matrix_SvmL() const;
   void calculate_MSvmL();
   double get_mass_matrix_SvtL() const;
   void calculate_MSvtL();
   Eigen::Matrix<double,2,2> get_mass_matrix_Sd() const;
   void calculate_MSd();
   Eigen::Matrix<double,2,2> get_mass_matrix_Su() const;
   void calculate_MSu();
   Eigen::Matrix<double,2,2> get_mass_matrix_Se() const;
   void calculate_MSe();
   Eigen::Matrix<double,2,2> get_mass_matrix_Sm() const;
   void calculate_MSm();
   Eigen::Matrix<double,2,2> get_mass_matrix_Stau() const;
   void calculate_MStau();
   Eigen::Matrix<double,2,2> get_mass_matrix_Ss() const;
   void calculate_MSs();
   Eigen::Matrix<double,2,2> get_mass_matrix_Sc() const;
   void calculate_MSc();
   Eigen::Matrix<double,2,2> get_mass_matrix_Sb() const;
   void calculate_MSb();
   Eigen::Matrix<double,2,2> get_mass_matrix_St() const;
   void calculate_MSt();
   Eigen::Matrix<double,2,2> get_mass_matrix_hh() const;
   void calculate_Mhh();
   Eigen::Matrix<double,2,2> get_mass_matrix_Ah() const;
   void calculate_MAh();
   Eigen::Matrix<double,2,2> get_mass_matrix_Hpm() const;
   void calculate_MHpm();
   Eigen::Matrix<double,4,4> get_mass_matrix_Chi() const;
   void calculate_MChi();
   Eigen::Matrix<double,2,2> get_mass_matrix_Cha() const;
   void calculate_MCha();
   double get_mass_matrix_VWm() const;
   void calculate_MVWm();

   double get_ewsb_eq_hh_1() const;
   double get_ewsb_eq_hh_2() const;

   double ThetaW() const;
   double v() const;


private:
   bool force_output{false};           ///< switch to force output of pole masses
   MSSMNoFV_onshell_physical physical; ///< contains the pole masses and mixings
   MSSMNoFV_onshell_problems problems; ///< problems

   int solve_ewsb_tree_level_via_soft_higgs_masses();

   // DR-bar masses
   double MVG{0.0};
   double MGlu{0.0};
   double MVP{0.0};
   double MVZ{0.0};
   double MVWm{0.0};
   double MFd{0.0};
   double MFs{0.0};
   double MFb{0.0};
   double MFu{0.0};
   double MFc{0.0};
   double MFt{0.0};
   double MFve{0.0};
   double MFvm{0.0};
   double MFvt{0.0};
   double MFe{0.0};
   double MFm{0.0};
   double MFtau{0.0};
   double MSveL{0.0};
   double MSvmL{0.0};
   double MSvtL{0.0};
   Eigen::Array<double,2,1> MSd{Eigen::Array<double,2,1>::Zero()};
   Eigen::Array<double,2,1> MSu{Eigen::Array<double,2,1>::Zero()};
   Eigen::Array<double,2,1> MSe{Eigen::Array<double,2,1>::Zero()};
   Eigen::Array<double,2,1> MSm{Eigen::Array<double,2,1>::Zero()};
   Eigen::Array<double,2,1> MStau{Eigen::Array<double,2,1>::Zero()};
   Eigen::Array<double,2,1> MSs{Eigen::Array<double,2,1>::Zero()};
   Eigen::Array<double,2,1> MSc{Eigen::Array<double,2,1>::Zero()};
   Eigen::Array<double,2,1> MSb{Eigen::Array<double,2,1>::Zero()};
   Eigen::Array<double,2,1> MSt{Eigen::Array<double,2,1>::Zero()};
   Eigen::Array<double,2,1> Mhh{Eigen::Array<double,2,1>::Zero()};
   Eigen::Array<double,2,1> MAh{Eigen::Array<double,2,1>::Zero()};
   Eigen::Array<double,2,1> MHpm{Eigen::Array<double,2,1>::Zero()};
   Eigen::Array<double,4,1> MChi{Eigen::Array<double,4,1>::Zero()};
   Eigen::Array<double,2,1> MCha{Eigen::Array<double,2,1>::Zero()};

   // DR-bar mixing matrices
   Eigen::Matrix<double,2,2> ZD{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<double,2,2> ZU{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<double,2,2> ZE{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<double,2,2> ZM{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<double,2,2> ZTau{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<double,2,2> ZS{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<double,2,2> ZC{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<double,2,2> ZB{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<double,2,2> ZT{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<double,2,2> ZH{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<double,2,2> ZA{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<double,2,2> ZP{Eigen::Matrix<double,2,2>::Zero()};
   Eigen::Matrix<std::complex<double>,4,4> ZN{Eigen::Matrix<std::complex<double>,4,4>::Zero()};
   Eigen::Matrix<std::complex<double>,2,2> UM{Eigen::Matrix<std::complex<double>,2,2>::Zero()};
   Eigen::Matrix<std::complex<double>,2,2> UP{Eigen::Matrix<std::complex<double>,2,2>::Zero()};

   // phases
   std::complex<double> PhaseGlu{1.0, 0.0};
};

std::ostream& operator<<(std::ostream&, const MSSMNoFV_onshell_mass_eigenstates&);

} // namespace gm2calc

#endif
