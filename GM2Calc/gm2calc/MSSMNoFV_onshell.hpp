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

#ifndef GM2_MSSMNoFV_ONSHELL_HPP
#define GM2_MSSMNoFV_ONSHELL_HPP

#include "MSSMNoFV_onshell_mass_eigenstates.hpp"

#include <iosfwd>

#include <Eigen/Core>

namespace gm2calc {

/**
 * @class MSSMNoFV_onshell
 * @brief contains the MSSMNoFV parameters in the on-shell scheme
 *
 * In addition, the class stores parameters necessary for the
 * calculation of \f$a_\mu\f$: The electromagnetic coupling at MZ w/o
 * hadronic corrections, \f$\alpha_{\text{em}}(M_Z)\f$, and the
 * electromagnetic coupling in the Thomson limit,
 * \f$\alpha_{\text{em}}(0)\f$.
 *
 * The class also contains helper functions to convert DR-bar
 * parameters to the on-shell scheme.
 */
class MSSMNoFV_onshell : public MSSMNoFV_onshell_mass_eigenstates {
public:
   MSSMNoFV_onshell();
   virtual ~MSSMNoFV_onshell() {}

   /// enable/disable verbose output
   void set_verbose_output(bool flag) { verbose_output = flag; }
   /// tests for verbose output
   bool do_verbose_output() const { return verbose_output; }

   /// set alpha(MZ) w/o hadronic corrections
   void set_alpha_MZ(double);
   /// set alpha in the Thomson limit
   void set_alpha_thompson(double);
   /// soft-breaking trilinear on-shell down-type slepton coupling
   void set_Ae(const Eigen::Matrix<double,3,3>& A) { Ae = A; }
   /// soft-breaking trilinear up-type squark coupling
   void set_Au(const Eigen::Matrix<double,3,3>& A) { Au = A; }
   /// soft-breaking trilinear down-type squark coupling
   void set_Ad(const Eigen::Matrix<double,3,3>& A) { Ad = A; }
   /// soft-breaking trilinear on-shell down-type slepton coupling
   void set_Ae(unsigned i, unsigned k, double a) { Ae(i,k) = a; }
   /// soft-breaking trilinear up-type squark coupling
   void set_Au(unsigned i, unsigned k, double a) { Au(i,k) = a; }
   /// soft-breaking trilinear down-type squark coupling
   void set_Ad(unsigned i, unsigned k, double a) { Ad(i,k) = a; }
   /// set CP-odd Higgs pole mass
   void set_MA0(double m) { get_physical().MAh(1) = m; }
   /// set tan(beta)
   void set_TB(double);

   /// electromagnetic gauge coupling at MZ w/o hadronic corrections
   double get_EL() const { return EL; }
   /// electromagnetic gauge coupling in Thomson limit
   double get_EL0() const { return EL0; }
   /// Hypercharge gauge coupling
   double get_gY() const { return sqrt(0.6) * get_g1(); }
   /// tan(beta) DR-bar
   double get_TB() const;
   /// Vacuum expectation value v
   double get_vev() const;
   /// soft-breaking trilinear on-shell down-type slepton coupling
   const Eigen::Matrix<double,3,3>& get_Ae() const { return Ae; }
   /// soft-breaking trilinear up-type squark coupling
   const Eigen::Matrix<double,3,3>& get_Au() const { return Au; }
   /// soft-breaking trilinear down-type squark coupling
   const Eigen::Matrix<double,3,3>& get_Ad() const { return Ad; }
   /// soft-breaking trilinear on-shell down-type slepton coupling
   double get_Ae(unsigned i, unsigned k) const { return Ae(i,k); }
   /// soft-breaking trilinear up-type squark coupling
   double get_Au(unsigned i, unsigned k) const { return Au(i,k); }
   /// soft-breaking trilinear on-shell down-type slepton coupling
   double get_Ad(unsigned i, unsigned k) const { return Ad(i,k); }
   /// returns W boson pole mass
   double get_MW() const { return get_physical().MVWm; }
   /// returns Z boson pole mass
   double get_MZ() const { return get_physical().MVZ; }
   /// returns electron mass
   double get_ME() const { return get_physical().MFe; }
   /// returns muon pole mass
   double get_MM() const { return get_physical().MFm; }
   /// returns tau mass
   double get_ML() const { return get_physical().MFtau; }
   /// returns up-quark mass
   double get_MU() const { return get_physical().MFu; }
   /// returns charm-quark mass
   double get_MC() const { return get_physical().MFc; }
   /// returns top-quark mass
   double get_MT() const { return get_physical().MFt; }
   /// returns down-quark mass
   double get_MD() const { return get_physical().MFd; }
   /// returns strange-quark mass
   double get_MS() const { return get_physical().MFs; }
   /// returns mb(mb) MS-bar
   double get_MBMB() const { return get_physical().MFb; }
   /// returns mb(MZ) DR-bar
   double get_MB() const { return mb_DRbar_MZ; }
   /// returns CP-odd Higgs mass
   double get_MA0() const { return get_physical().MAh(1); }
   /// returns selectron mixing matrix
   const Eigen::Matrix<double,2,2>& get_USe() const { return get_ZE(); }
   /// returns smuon pole mixing matrix
   const Eigen::Matrix<double,2,2>& get_USm() const { return get_ZM(); }
   /// returns stau mixing matrix
   const Eigen::Matrix<double,2,2>& get_UStau() const { return get_ZTau(); }
   /// returns sup mixing matrix
   const Eigen::Matrix<double,2,2>& get_USu() const { return get_ZU(); }
   /// returns sdown mixing matrix
   const Eigen::Matrix<double,2,2>& get_USd() const { return get_ZD(); }
   /// returns scharm mixing matrix
   const Eigen::Matrix<double,2,2>& get_USc() const { return get_ZC(); }
   /// returns sstrange mixing matrix
   const Eigen::Matrix<double,2,2>& get_USs() const { return get_ZS(); }
   /// returns sbottom mixing matrix
   const Eigen::Matrix<double,2,2>& get_USb() const { return get_ZB(); }
   /// returns stop mixing matrix
   const Eigen::Matrix<double,2,2>& get_USt() const { return get_ZT(); }

   /// calculate SUSY masses in mixed on-shell/DR-bar scheme from given on-shell/DR-bar parameters
   void calculate_masses();

   /// convert given (DR-bar) parameters to mixed on-shell/DR-bar scheme
   void convert_to_onshell(double precision = 1e-8,
                           unsigned max_iterations = 1000);

   /// convert mixed on-shell/DR-bar parameters to non-tan(beta) resummed case
   void convert_to_non_tan_beta_resummed();

private:
   bool verbose_output{false}; ///< verbose output
   double EL{0.0};             ///< electromagnetic gauge coupling at MZ w/o hadronic corrections
   double EL0{0.0};            ///< electromagnetic gauge coupling in the Thomson limit
   double mb_DRbar_MZ{2.8};    ///< mb(MZ) DR-bar
   Eigen::Matrix<double,3,3> Au{Eigen::Matrix<double,3,3>::Zero()}; ///< trilinear couplings
   Eigen::Matrix<double,3,3> Ad{Eigen::Matrix<double,3,3>::Zero()}; ///< trilinear couplings
   Eigen::Matrix<double,3,3> Ae{Eigen::Matrix<double,3,3>::Zero()}; ///< trilinear couplings

   void check_input() const;
   void check_problems() const;
   void calculate_mb_DRbar_MZ();
   void convert_gauge_couplings();
   void convert_yukawa_couplings_treelevel();
   void convert_BMu();
   void convert_ml2();
   void convert_me2(double, unsigned);
   double convert_me2_fpi(double, unsigned);
   double convert_me2_fpi_modify(double, unsigned);
   double convert_me2_root(double, unsigned);
   double convert_me2_root_modify(double, unsigned);
   void convert_Mu_M1_M2(double, unsigned);
   void convert_vev();
   void convert_yukawa_couplings();
   void copy_susy_masses_to_pole();
   unsigned find_bino_like_neutralino();
};

/// streaming operator
std::ostream& operator<<(std::ostream&, const MSSMNoFV_onshell&);

} // namespace gm2calc

#endif
