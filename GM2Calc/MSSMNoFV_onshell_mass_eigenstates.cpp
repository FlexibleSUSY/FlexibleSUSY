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
 * @file MSSMNoFV_onshell_mass_eigenstates.cpp
 * @brief implementation of the MSSMNoFV_onshell model class
 *
 * Contains the definition of the MSSMNoFV_onshell model class methods
 * which solve EWSB and calculate pole masses and mixings from DRbar
 * parameters.
 *
 * This file was generated at Wed 22 Jul 2015 18:14:31 with FlexibleSUSY
 * 1.2.1 (git commit: v1.2.1-49-gfc4c300) and SARAH 4.5.8 .
 */

#include "gm2calc/MSSMNoFV_onshell_mass_eigenstates.hpp"

#include "gm2_linalg.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

namespace {

/**
 * @class RAII_save
 * @brief Saves value of variable and restores it at destruction
 */
template <typename T>
class RAII_save {
public:
   explicit RAII_save(T& var_) noexcept : var(var_), value(var_) {}
   RAII_save(const RAII_save&) = delete;
   RAII_save(RAII_save&&) noexcept = default;
   ~RAII_save() { var = value; }
   RAII_save& operator=(const RAII_save&) = delete;
   RAII_save& operator=(RAII_save&& other) noexcept = default;

private:
   T& var;
   T value{};
};

template <typename T>
constexpr RAII_save<T> make_raii_save(T& var)
{
   return RAII_save<T>(var);
}

template <typename Derived>
void symmetrize(Eigen::MatrixBase<Derived>& m)
{
   static_assert(Eigen::MatrixBase<Derived>::RowsAtCompileTime ==
                 Eigen::MatrixBase<Derived>::ColsAtCompileTime,
                 "symmetrize is only defined for squared matrices");

   for (int i = 0; i < Eigen::MatrixBase<Derived>::RowsAtCompileTime; i++) {
      for (int k = 0; k < i; k++) {
         m(i,k) = m(k,i);
      }
   }
}

template <typename Derived>
void hermitianize(Eigen::MatrixBase<Derived>& m)
{
   symmetrize(m);
}

template <typename T> T sqr(T x) { return x*x; }

template<class Real, int N>
struct Abs_less {
    explicit Abs_less(const Eigen::Array<Real, N, 1>& w_) : w(w_) {}
    bool operator() (int i, int j) { return std::abs(w[i]) < std::abs(w[j]); }
    const Eigen::Array<Real, N, 1>& w;
};

template <class T>
struct Is_not_finite {
   bool operator()(T x) { return !std::isfinite(x); }
};

template <typename Derived>
unsigned closest_index(double mass, const Eigen::ArrayBase<Derived>& v)
{
   unsigned pos;
   typename Derived::PlainObject tmp;
   tmp.setConstant(mass);

   (v - tmp).abs().minCoeff(&pos);

   return pos;
}

/**
 * The element of v, which is closest to mass, is moved to the
 * position idx.
 *
 * @param idx new index of the mass eigenvalue
 * @param mass mass to compare against
 * @param v vector of masses
 * @param z corresponding mixing matrix
 */
template <typename DerivedArray, typename DerivedMatrix>
void move_goldstone_to(int idx, double mass, Eigen::ArrayBase<DerivedArray>& v,
                       Eigen::MatrixBase<DerivedMatrix>& z)
{
   int pos = closest_index(mass, v);
   if (pos == idx) {
      return;
   }

   const int sign = (idx - pos) < 0 ? -1 : 1;
   int steps = std::abs(idx - pos);

   // now we shuffle the states
   while (steps--) {
      const int new_pos = pos + sign;
      v.row(new_pos).swap(v.row(pos));
      z.row(new_pos).swap(z.row(pos));
      pos = new_pos;
   }
}

/**
 * Copies all elements from src to dst which are not close to the
 * elements in cmp.
 *
 * @param src source vector
 * @param cmp vector with elements to compare against
 * @param dst destination vector
 */
template<class Real, int Nsrc, int Ncmp, int Ndst>
void remove_if_equal(const Eigen::Array<Real,Nsrc,1>& src,
                     const Eigen::Array<Real,Ncmp,1>& cmp,
                     Eigen::Array<Real,Ndst,1>& dst)
{
   static_assert(Nsrc == Ncmp + Ndst,
                 "Error: remove_if_equal: vectors have incompatible length!");

   Eigen::Array<Real,Nsrc,1> non_equal(src);

   for (int i = 0; i < Ncmp; i++) {
      const int idx = closest_index(cmp(i), non_equal);
      non_equal(idx) = std::numeric_limits<double>::infinity();
   }

   std::remove_copy_if(non_equal.data(), non_equal.data() + Nsrc,
                       dst.data(), Is_not_finite<Real>());
}

/**
 * @brief reorders vector v according to ordering in vector v2
 * @param v vector with elementes to be reordered
 * @param v2 vector with reference ordering
 */
template<class Real, int N>
void reorder_vector(
   Eigen::Array<Real,N,1>& v,
   const Eigen::Array<Real,N,1>& v2)
{
   Eigen::PermutationMatrix<N> p;
   p.setIdentity();
   std::sort(p.indices().data(), p.indices().data() + p.indices().size(),
             Abs_less<Real, N>(v2));

#if EIGEN_VERSION_AT_LEAST(3,1,4)
   v.matrix().transpose() *= p.inverse();
#else
   Eigen::Map<Eigen::Matrix<Real,N,1> >(v.data()).transpose() *= p.inverse();
#endif
}

/**
 * @brief reorders vector v according to ordering of diagonal elements in mass_matrix
 * @param v vector with elementes to be reordered
 * @param matrix matrix with diagonal elements with reference ordering
 */
template<class Derived>
void reorder_vector(
   Eigen::Array<double,Eigen::MatrixBase<Derived>::RowsAtCompileTime,1>& v,
   const Eigen::MatrixBase<Derived>& matrix)
{
   reorder_vector(v, matrix.diagonal().array().eval());
}

} // anonymous namespace

namespace gm2calc {

#define CLASSNAME MSSMNoFV_onshell_mass_eigenstates
#define PHYSICAL(parameter) physical.parameter
#define MODELPARAMETER(parameter) model->get_##parameter()

void CLASSNAME::do_force_output(bool flag)
{
   force_output = flag;
}

bool CLASSNAME::do_force_output() const
{
   return force_output;
}

const MSSMNoFV_onshell_physical& CLASSNAME::get_physical() const
{
   return physical;
}

MSSMNoFV_onshell_physical& CLASSNAME::get_physical()
{
   return physical;
}

void CLASSNAME::set_physical(const MSSMNoFV_onshell_physical& physical_)
{
   physical = physical_;
}

const MSSMNoFV_onshell_problems& CLASSNAME::get_problems() const
{
   return problems;
}

MSSMNoFV_onshell_problems& CLASSNAME::get_problems()
{
   return problems;
}

int CLASSNAME::solve_ewsb_tree_level()
{
   return solve_ewsb_tree_level_via_soft_higgs_masses();
}

int CLASSNAME::solve_ewsb_tree_level_via_soft_higgs_masses()
{
   int error = 0;

   const double old_mHd2 = mHd2;
   const double old_mHu2 = mHu2;

   mHd2 = (0.025*(-40*vd*sqr(Mu) + 20*vu*BMu + 20*vu*
      BMu - 3*std::pow(vd,3)*sqr(g1) - 5*std::pow(vd,3)*sqr(g2) + 3*vd*sqr(g1)*sqr
      (vu) + 5*vd*sqr(g2)*sqr(vu)))/vd;
   mHu2 = (0.025*(-40*vu*sqr(Mu) + 20*vd*BMu + 20*vd*
      BMu - 3*std::pow(vu,3)*sqr(g1) - 5*std::pow(vu,3)*sqr(g2) + 3*vu*sqr(g1)*sqr
      (vd) + 5*vu*sqr(g2)*sqr(vd)))/vu;

   const bool is_finite = std::isfinite(mHd2) && std::isfinite(mHu2);

   if (!is_finite) {
      mHd2 = old_mHd2;
      mHu2 = old_mHu2;
      error = 1;
   }

   return error;
}

int CLASSNAME::solve_ewsb()
{
   return solve_ewsb_tree_level();
}

void CLASSNAME::print(std::ostream& ostr) const
{
   ostr << "========================================\n"
           "MSSMNoFV_onshell\n"
           "========================================\n";
   MSSMNoFV_onshell_soft_parameters::print(ostr);
   ostr << "----------------------------------------\n"
           "tree-level DRbar masses:\n"
           "----------------------------------------\n";
   ostr << "MVG = " << MVG << '\n';
   ostr << "MGlu = " << MGlu << '\n';
   ostr << "MVP = " << MVP << '\n';
   ostr << "MVZ = " << MVZ << '\n';
   ostr << "MFd = " << MFd << '\n';
   ostr << "MFs = " << MFs << '\n';
   ostr << "MFb = " << MFb << '\n';
   ostr << "MFu = " << MFu << '\n';
   ostr << "MFc = " << MFc << '\n';
   ostr << "MFt = " << MFt << '\n';
   ostr << "MFve = " << MFve << '\n';
   ostr << "MFvm = " << MFvm << '\n';
   ostr << "MFvt = " << MFvt << '\n';
   ostr << "MFe = " << MFe << '\n';
   ostr << "MFm = " << MFm << '\n';
   ostr << "MFtau = " << MFtau << '\n';
   ostr << "MSveL = " << MSveL << '\n';
   ostr << "MSvmL = " << MSvmL << '\n';
   ostr << "MSvtL = " << MSvtL << '\n';
   ostr << "MSd = " << MSd.transpose() << '\n';
   ostr << "MSu = " << MSu.transpose() << '\n';
   ostr << "MSe = " << MSe.transpose() << '\n';
   ostr << "MSm = " << MSm.transpose() << '\n';
   ostr << "MStau = " << MStau.transpose() << '\n';
   ostr << "MSs = " << MSs.transpose() << '\n';
   ostr << "MSc = " << MSc.transpose() << '\n';
   ostr << "MSb = " << MSb.transpose() << '\n';
   ostr << "MSt = " << MSt.transpose() << '\n';
   ostr << "Mhh = " << Mhh.transpose() << '\n';
   ostr << "MAh = " << MAh.transpose() << '\n';
   ostr << "MHpm = " << MHpm.transpose() << '\n';
   ostr << "MChi = " << MChi.transpose() << '\n';
   ostr << "MCha = " << MCha.transpose() << '\n';
   ostr << "MVWm = " << MVWm << '\n';

   ostr << "----------------------------------------\n"
           "tree-level DRbar mixing matrices:\n"
           "----------------------------------------\n";
   ostr << "ZD = " << ZD << '\n';
   ostr << "ZU = " << ZU << '\n';
   ostr << "ZE = " << ZE << '\n';
   ostr << "ZM = " << ZM << '\n';
   ostr << "ZTau = " << ZTau << '\n';
   ostr << "ZS = " << ZS << '\n';
   ostr << "ZC = " << ZC << '\n';
   ostr << "ZB = " << ZB << '\n';
   ostr << "ZT = " << ZT << '\n';
   ostr << "ZH = " << ZH << '\n';
   ostr << "ZA = " << ZA << '\n';
   ostr << "ZP = " << ZP << '\n';
   ostr << "ZN = " << ZN << '\n';
   ostr << "UM = " << UM << '\n';
   ostr << "UP = " << UP << '\n';

   physical.print(ostr);
   problems.print(ostr);
}

/**
 * routine which finds the DRbar mass eigenstates and mixings.
 */
void CLASSNAME::calculate_DRbar_masses()
{
   const auto save_mHd2_raii = make_raii_save(mHd2);
   const auto save_mHu2_raii = make_raii_save(mHu2);

   solve_ewsb_tree_level_via_soft_higgs_masses();

   calculate_MVG();
   calculate_MVP();
   calculate_MVZ();
   calculate_MVWm();
   calculate_MGlu();
   calculate_MFd();
   calculate_MFs();
   calculate_MFb();
   calculate_MFu();
   calculate_MFc();
   calculate_MFt();
   calculate_MFve();
   calculate_MFvm();
   calculate_MFvt();
   calculate_MFe();
   calculate_MFm();
   calculate_MFtau();
   calculate_MSveL();
   calculate_MSvmL();
   calculate_MSvtL();
   calculate_MSd();
   calculate_MSu();
   calculate_MSe();
   calculate_MSm();
   calculate_MStau();
   calculate_MSs();
   calculate_MSc();
   calculate_MSb();
   calculate_MSt();
   calculate_Mhh();
   calculate_MAh();
   calculate_MHpm();
   calculate_MChi();
   calculate_MCha();
}

void CLASSNAME::copy_DRbar_masses_to_pole_masses()
{
   PHYSICAL(MVG) = MVG;
   PHYSICAL(MGlu) = MGlu;
   PHYSICAL(MVP) = MVP;
   PHYSICAL(MVZ) = MVZ;
   PHYSICAL(MFd) = MFd;
   PHYSICAL(MFs) = MFs;
   PHYSICAL(MFb) = MFb;
   PHYSICAL(MFu) = MFu;
   PHYSICAL(MFc) = MFc;
   PHYSICAL(MFt) = MFt;
   PHYSICAL(MFve) = MFve;
   PHYSICAL(MFvm) = MFvm;
   PHYSICAL(MFvt) = MFvt;
   PHYSICAL(MFe) = MFe;
   PHYSICAL(MFm) = MFm;
   PHYSICAL(MFtau) = MFtau;
   PHYSICAL(MSveL) = MSveL;
   PHYSICAL(MSvmL) = MSvmL;
   PHYSICAL(MSvtL) = MSvtL;
   PHYSICAL(MSd) = MSd;
   PHYSICAL(ZD) = ZD;
   PHYSICAL(MSu) = MSu;
   PHYSICAL(ZU) = ZU;
   PHYSICAL(MSe) = MSe;
   PHYSICAL(ZE) = ZE;
   PHYSICAL(MSm) = MSm;
   PHYSICAL(ZM) = ZM;
   PHYSICAL(MStau) = MStau;
   PHYSICAL(ZTau) = ZTau;
   PHYSICAL(MSs) = MSs;
   PHYSICAL(ZS) = ZS;
   PHYSICAL(MSc) = MSc;
   PHYSICAL(ZC) = ZC;
   PHYSICAL(MSb) = MSb;
   PHYSICAL(ZB) = ZB;
   PHYSICAL(MSt) = MSt;
   PHYSICAL(ZT) = ZT;
   PHYSICAL(Mhh) = Mhh;
   PHYSICAL(ZH) = ZH;
   PHYSICAL(MAh) = MAh;
   PHYSICAL(ZA) = ZA;
   PHYSICAL(MHpm) = MHpm;
   PHYSICAL(ZP) = ZP;
   PHYSICAL(MChi) = MChi;
   PHYSICAL(ZN) = ZN;
   PHYSICAL(MCha) = MCha;
   PHYSICAL(UM) = UM;
   PHYSICAL(UP) = UP;
   PHYSICAL(MVWm) = MVWm;

}

/**
 * reorders DRbar masses so that golstones are placed at the index
 * specified in the model files definition of the associated
 * gauge boson (see Z-boson definition in default particles.m file
 * in the Models directory of your SARAH distribution for example)
 */
void CLASSNAME::reorder_DRbar_masses()
{
   move_goldstone_to(0, MVZ, MAh, ZA);
   move_goldstone_to(0, MVWm, MHpm, ZP);
}

/**
 * reorders pole masses so that golstones are placed at the index
 * specified in the model files definition of the associated
 * gauge boson (see Z-boson definition in default particles.m file
 * in the Models directory of your SARAH distribution for example)
 */
void CLASSNAME::reorder_pole_masses()
{
   move_goldstone_to(0, MVZ, PHYSICAL(MAh), PHYSICAL(ZA));
   move_goldstone_to(0, MVWm, PHYSICAL(MHpm), PHYSICAL(ZP));
}

Eigen::Array<double,1,1> CLASSNAME::get_MChargedHiggs() const
{
   Eigen::Array<double,1,1> MHpm_ChargedHiggs;
   Eigen::Array<double,1,1> MHpm_goldstone;

   MHpm_goldstone(0) = MVWm;

   remove_if_equal(MHpm, MHpm_goldstone, MHpm_ChargedHiggs);

   return MHpm_ChargedHiggs;
}

Eigen::Array<double,1,1> CLASSNAME::get_MPseudoscalarHiggs() const
{
   Eigen::Array<double,1,1> MAh_PseudoscalarHiggs;
   Eigen::Array<double,1,1> MAh_goldstone;

   MAh_goldstone(0) = MVZ;

   remove_if_equal(MAh, MAh_goldstone, MAh_PseudoscalarHiggs);

   return MAh_PseudoscalarHiggs;
}




double CLASSNAME::get_mass_matrix_VG() const
{
   return 0;
}

void CLASSNAME::calculate_MVG()
{
   MVG = 0;
}

double CLASSNAME::get_mass_matrix_Glu() const
{
   return MassG;
}

void CLASSNAME::calculate_MGlu()
{
   const auto mass_matrix_Glu = get_mass_matrix_Glu();
   PhaseGlu = std::polar(1., 0.5 * std::arg(std::complex<double>(mass_matrix_Glu)));
   MGlu = std::abs(mass_matrix_Glu);
}

double CLASSNAME::get_mass_matrix_VP() const
{
   return 0;
}

void CLASSNAME::calculate_MVP()
{
   const auto mass_matrix_VP = get_mass_matrix_VP();
   MVP = std::abs(mass_matrix_VP);
}

double CLASSNAME::get_mass_matrix_VZ() const
{
   const double mass_matrix_VZ = 0.25*(sqr(vd) + sqr(vu))*sqr(g2*std::cos(
      ThetaW()) + 0.7745966692414834*g1*std::sin(ThetaW()));

   return mass_matrix_VZ;
}

void CLASSNAME::calculate_MVZ()
{
   const auto mass_matrix_VZ = get_mass_matrix_VZ();
   MVZ = std::sqrt(std::abs(mass_matrix_VZ));
}

double CLASSNAME::get_mass_matrix_Fd() const
{
   const double mass_matrix_Fd = 0.7071067811865475*vd*Yd(0,0);

   return mass_matrix_Fd;
}

void CLASSNAME::calculate_MFd()
{
   const auto mass_matrix_Fd = get_mass_matrix_Fd();
   MFd = std::abs(mass_matrix_Fd);
}

double CLASSNAME::get_mass_matrix_Fs() const
{
   const double mass_matrix_Fs = 0.7071067811865475*vd*Yd(1,1);

   return mass_matrix_Fs;
}

void CLASSNAME::calculate_MFs()
{
   const auto mass_matrix_Fs = get_mass_matrix_Fs();
   MFs = std::abs(mass_matrix_Fs);
}

double CLASSNAME::get_mass_matrix_Fb() const
{
   const double mass_matrix_Fb = 0.7071067811865475*vd*Yd(2,2);

   return mass_matrix_Fb;
}

void CLASSNAME::calculate_MFb()
{
   const auto mass_matrix_Fb = get_mass_matrix_Fb();
   MFb = std::abs(mass_matrix_Fb);
}

double CLASSNAME::get_mass_matrix_Fu() const
{
   const double mass_matrix_Fu = 0.7071067811865475*vu*Yu(0,0);

   return mass_matrix_Fu;
}

void CLASSNAME::calculate_MFu()
{
   const auto mass_matrix_Fu = get_mass_matrix_Fu();
   MFu = std::abs(mass_matrix_Fu);
}

double CLASSNAME::get_mass_matrix_Fc() const
{
   const double mass_matrix_Fc = 0.7071067811865475*vu*Yu(1,1);

   return mass_matrix_Fc;
}

void CLASSNAME::calculate_MFc()
{
   const auto mass_matrix_Fc = get_mass_matrix_Fc();
   MFc = std::abs(mass_matrix_Fc);
}

double CLASSNAME::get_mass_matrix_Ft() const
{
   const double mass_matrix_Ft = 0.7071067811865475*vu*Yu(2,2);

   return mass_matrix_Ft;
}

void CLASSNAME::calculate_MFt()
{
   const auto mass_matrix_Ft = get_mass_matrix_Ft();
   MFt = std::abs(mass_matrix_Ft);
}

double CLASSNAME::get_mass_matrix_Fve() const
{
   return 0;
}

void CLASSNAME::calculate_MFve()
{
   const auto mass_matrix_Fve = get_mass_matrix_Fve();
   MFve = std::abs(mass_matrix_Fve);
}

double CLASSNAME::get_mass_matrix_Fvm() const
{
   return 0;
}

void CLASSNAME::calculate_MFvm()
{
   const auto mass_matrix_Fvm = get_mass_matrix_Fvm();
   MFvm = std::abs(mass_matrix_Fvm);
}

double CLASSNAME::get_mass_matrix_Fvt() const
{
   return 0;
}

void CLASSNAME::calculate_MFvt()
{
   const auto mass_matrix_Fvt = get_mass_matrix_Fvt();
   MFvt = std::abs(mass_matrix_Fvt);
}

double CLASSNAME::get_mass_matrix_Fe() const
{
   const double mass_matrix_Fe = 0.7071067811865475*vd*Ye(0,0);

   return mass_matrix_Fe;
}

void CLASSNAME::calculate_MFe()
{
   const auto mass_matrix_Fe = get_mass_matrix_Fe();
   MFe = std::abs(mass_matrix_Fe);
}

double CLASSNAME::get_mass_matrix_Fm() const
{
   const double mass_matrix_Fm = 0.7071067811865475*vd*Ye(1,1);

   return mass_matrix_Fm;
}

void CLASSNAME::calculate_MFm()
{
   const auto mass_matrix_Fm = get_mass_matrix_Fm();
   MFm = std::abs(mass_matrix_Fm);
}

double CLASSNAME::get_mass_matrix_Ftau() const
{
   const double mass_matrix_Ftau = 0.7071067811865475*vd*Ye(2,2);

   return mass_matrix_Ftau;
}

void CLASSNAME::calculate_MFtau()
{
   const auto mass_matrix_Ftau = get_mass_matrix_Ftau();
   MFtau = std::abs(mass_matrix_Ftau);
}

double CLASSNAME::get_mass_matrix_SveL() const
{
   const double mass_matrix_SveL = 0.125*(8*ml2(0,0) - 0.6*sqr(g1)*(
      -sqr(vd) + sqr(vu)) - sqr(g2)*(-sqr(vd) + sqr(vu)));

   return mass_matrix_SveL;
}

void CLASSNAME::calculate_MSveL()
{
   const auto mass_matrix_SveL = get_mass_matrix_SveL();
   MSveL = std::sqrt(std::abs(mass_matrix_SveL));
}

double CLASSNAME::get_mass_matrix_SvmL() const
{
   const double mass_matrix_SvmL = 0.125*(8*ml2(1,1) - 0.6*sqr(g1)*(
      -sqr(vd) + sqr(vu)) - sqr(g2)*(-sqr(vd) + sqr(vu)));

   return mass_matrix_SvmL;
}

void CLASSNAME::calculate_MSvmL()
{
   const auto mass_matrix_SvmL = get_mass_matrix_SvmL();

   if (mass_matrix_SvmL < 0.) {
      problems.flag_tachyon("SvmL");
   }

   MSvmL = std::sqrt(std::abs(mass_matrix_SvmL));
}

double CLASSNAME::get_mass_matrix_SvtL() const
{
   const double mass_matrix_SvtL = 0.125*(8*ml2(2,2) - 0.6*sqr(g1)*(
      -sqr(vd) + sqr(vu)) - sqr(g2)*(-sqr(vd) + sqr(vu)));

   return mass_matrix_SvtL;
}

void CLASSNAME::calculate_MSvtL()
{
   const auto mass_matrix_SvtL = get_mass_matrix_SvtL();
   MSvtL = std::sqrt(std::abs(mass_matrix_SvtL));
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Sd() const
{
   Eigen::Matrix<double,2,2> mass_matrix_Sd;

   mass_matrix_Sd(0,0) = mq2(0,0) + 0.5*sqr(Yd(0,0))*sqr(vd) - 0.025*
      sqr(g1)*sqr(vd) - 0.125*sqr(g2)*sqr(vd) + 0.025*sqr(g1)*sqr(vu) + 0.125*
      sqr(g2)*sqr(vu);
   mass_matrix_Sd(0,1) = 0.7071067811865475*vd*TYd(0,0) -
      0.7071067811865475*vu*Yd(0,0)*Mu;
   mass_matrix_Sd(1,1) = md2(0,0) + 0.5*sqr(Yd(0,0))*sqr(vd) - 0.05*
      sqr(g1)*sqr(vd) + 0.05*sqr(g1)*sqr(vu);

   hermitianize(mass_matrix_Sd);

   return mass_matrix_Sd;
}

void CLASSNAME::calculate_MSd()
{
   const auto mass_matrix_Sd(get_mass_matrix_Sd());
   fs_diagonalize_hermitian(mass_matrix_Sd, MSd, ZD);
   MSd = sqrt(MSd.cwiseAbs());
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Su() const
{
   Eigen::Matrix<double,2,2> mass_matrix_Su;

   mass_matrix_Su(0,0) = mq2(0,0) - 0.025*sqr(g1)*sqr(vd) + 0.125*sqr(g2)
      *sqr(vd) + 0.5*sqr(Yu(0,0))*sqr(vu) + 0.025*sqr(g1)*sqr(vu) - 0.125*
      sqr(g2)*sqr(vu);
   mass_matrix_Su(0,1) = 0.7071067811865475*vu*TYu(0,0) -
      0.7071067811865475*vd*Yu(0,0)*Mu;
   mass_matrix_Su(1,1) = mu2(0,0) + 0.1*sqr(g1)*sqr(vd) + 0.5*sqr(Yu(0
      ,0))*sqr(vu) - 0.1*sqr(g1)*sqr(vu);

   hermitianize(mass_matrix_Su);

   return mass_matrix_Su;
}

void CLASSNAME::calculate_MSu()
{
   const auto mass_matrix_Su(get_mass_matrix_Su());
   fs_diagonalize_hermitian(mass_matrix_Su, MSu, ZU);
   MSu = sqrt(MSu.cwiseAbs());
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Se() const
{
   Eigen::Matrix<double,2,2> mass_matrix_Se;

   mass_matrix_Se(0,0) = ml2(0,0) + 0.5*sqr(Ye(0,0))*sqr(vd) + 0.075*
      sqr(g1)*sqr(vd) - 0.125*sqr(g2)*sqr(vd) - 0.075*sqr(g1)*sqr(vu) + 0.125*
      sqr(g2)*sqr(vu);
   mass_matrix_Se(0,1) = 0.7071067811865475*vd*TYe(0,0) -
      0.7071067811865475*vu*Ye(0,0)*Mu;
   mass_matrix_Se(1,1) = me2(0,0) + 0.5*sqr(Ye(0,0))*sqr(vd) - 0.15*
      sqr(g1)*sqr(vd) + 0.15*sqr(g1)*sqr(vu);

   hermitianize(mass_matrix_Se);

   return mass_matrix_Se;
}

void CLASSNAME::calculate_MSe()
{
   const auto mass_matrix_Se(get_mass_matrix_Se());
   fs_diagonalize_hermitian(mass_matrix_Se, MSe, ZE);
   MSe = sqrt(MSe.cwiseAbs());
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Sm() const
{
   Eigen::Matrix<double,2,2> mass_matrix_Sm;

   mass_matrix_Sm(0,0) = ml2(1,1) + 0.5*sqr(Ye(1,1))*sqr(vd) + 0.075*
      sqr(g1)*sqr(vd) - 0.125*sqr(g2)*sqr(vd) - 0.075*sqr(g1)*sqr(vu) + 0.125*
      sqr(g2)*sqr(vu);
   mass_matrix_Sm(0,1) = 0.7071067811865475*vd*TYe(1,1) -
      0.7071067811865475*vu*Ye(1,1)*Mu;
   mass_matrix_Sm(1,1) = me2(1,1) + 0.5*sqr(Ye(1,1))*sqr(vd) - 0.15*
      sqr(g1)*sqr(vd) + 0.15*sqr(g1)*sqr(vu);

   hermitianize(mass_matrix_Sm);

   return mass_matrix_Sm;
}

void CLASSNAME::calculate_MSm()
{
   const auto mass_matrix_Sm(get_mass_matrix_Sm());
   fs_diagonalize_hermitian(mass_matrix_Sm, MSm, ZM);

   if (MSm.minCoeff() < 0.) {
      problems.flag_tachyon("Sm");
   }

   MSm = sqrt(MSm.cwiseAbs());
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Stau() const
{
   Eigen::Matrix<double,2,2> mass_matrix_Stau;

   mass_matrix_Stau(0,0) = ml2(2,2) + 0.5*sqr(Ye(2,2))*sqr(vd) + 0.075
      *sqr(g1)*sqr(vd) - 0.125*sqr(g2)*sqr(vd) - 0.075*sqr(g1)*sqr(vu) + 0.125*
      sqr(g2)*sqr(vu);
   mass_matrix_Stau(0,1) = 0.7071067811865475*vd*TYe(2,2) -
      0.7071067811865475*vu*Ye(2,2)*Mu;
   mass_matrix_Stau(1,1) = me2(2,2) + 0.5*sqr(Ye(2,2))*sqr(vd) - 0.15*
      sqr(g1)*sqr(vd) + 0.15*sqr(g1)*sqr(vu);

   hermitianize(mass_matrix_Stau);

   return mass_matrix_Stau;
}

void CLASSNAME::calculate_MStau()
{
   const auto mass_matrix_Stau(get_mass_matrix_Stau());
   fs_diagonalize_hermitian(mass_matrix_Stau, MStau, ZTau);

   if (MStau.minCoeff() < 0.) {
      problems.flag_tachyon("Stau");
   }

   MStau = sqrt(MStau.cwiseAbs());
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Ss() const
{
   Eigen::Matrix<double,2,2> mass_matrix_Ss;

   mass_matrix_Ss(0,0) = mq2(1,1) + 0.5*sqr(Yd(1,1))*sqr(vd) - 0.025*
      sqr(g1)*sqr(vd) - 0.125*sqr(g2)*sqr(vd) + 0.025*sqr(g1)*sqr(vu) + 0.125*
      sqr(g2)*sqr(vu);
   mass_matrix_Ss(0,1) = 0.7071067811865475*vd*TYd(1,1) -
      0.7071067811865475*vu*Yd(1,1)*Mu;
   mass_matrix_Ss(1,1) = md2(1,1) + 0.5*sqr(Yd(1,1))*sqr(vd) - 0.05*
      sqr(g1)*sqr(vd) + 0.05*sqr(g1)*sqr(vu);

   hermitianize(mass_matrix_Ss);

   return mass_matrix_Ss;
}

void CLASSNAME::calculate_MSs()
{
   const auto mass_matrix_Ss(get_mass_matrix_Ss());
   fs_diagonalize_hermitian(mass_matrix_Ss, MSs, ZS);
   MSs = sqrt(MSs.cwiseAbs());
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Sc() const
{
   Eigen::Matrix<double,2,2> mass_matrix_Sc;

   mass_matrix_Sc(0,0) = mq2(1,1) - 0.025*sqr(g1)*sqr(vd) + 0.125*sqr(g2)
      *sqr(vd) + 0.5*sqr(Yu(1,1))*sqr(vu) + 0.025*sqr(g1)*sqr(vu) - 0.125*
      sqr(g2)*sqr(vu);
   mass_matrix_Sc(0,1) = 0.7071067811865475*vu*TYu(1,1) -
      0.7071067811865475*vd*Yu(1,1)*Mu;
   mass_matrix_Sc(1,1) = mu2(1,1) + 0.1*sqr(g1)*sqr(vd) + 0.5*sqr(Yu(1
      ,1))*sqr(vu) - 0.1*sqr(g1)*sqr(vu);

   hermitianize(mass_matrix_Sc);

   return mass_matrix_Sc;
}

void CLASSNAME::calculate_MSc()
{
   const auto mass_matrix_Sc(get_mass_matrix_Sc());
   fs_diagonalize_hermitian(mass_matrix_Sc, MSc, ZC);
   MSc = sqrt(MSc.cwiseAbs());
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Sb() const
{
   Eigen::Matrix<double,2,2> mass_matrix_Sb;

   mass_matrix_Sb(0,0) = mq2(2,2) + 0.5*sqr(Yd(2,2))*sqr(vd) - 0.025*
      sqr(g1)*sqr(vd) - 0.125*sqr(g2)*sqr(vd) + 0.025*sqr(g1)*sqr(vu) + 0.125*
      sqr(g2)*sqr(vu);
   mass_matrix_Sb(0,1) = 0.7071067811865475*vd*TYd(2,2) -
      0.7071067811865475*vu*Yd(2,2)*Mu;
   mass_matrix_Sb(1,1) = md2(2,2) + 0.5*sqr(Yd(2,2))*sqr(vd) - 0.05*
      sqr(g1)*sqr(vd) + 0.05*sqr(g1)*sqr(vu);

   hermitianize(mass_matrix_Sb);

   return mass_matrix_Sb;
}

void CLASSNAME::calculate_MSb()
{
   const auto mass_matrix_Sb(get_mass_matrix_Sb());
   fs_diagonalize_hermitian(mass_matrix_Sb, MSb, ZB);

   if (MSb.minCoeff() < 0.) {
      problems.flag_tachyon("Sb");
   }

   MSb = sqrt(MSb.cwiseAbs());
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_St() const
{
   Eigen::Matrix<double,2,2> mass_matrix_St;

   mass_matrix_St(0,0) = mq2(2,2) - 0.025*sqr(g1)*sqr(vd) + 0.125*sqr(g2)
      *sqr(vd) + 0.5*sqr(Yu(2,2))*sqr(vu) + 0.025*sqr(g1)*sqr(vu) - 0.125*
      sqr(g2)*sqr(vu);
   mass_matrix_St(0,1) = 0.7071067811865475*vu*TYu(2,2) -
      0.7071067811865475*vd*Yu(2,2)*Mu;
   mass_matrix_St(1,1) = mu2(2,2) + 0.1*sqr(g1)*sqr(vd) + 0.5*sqr(Yu(2
      ,2))*sqr(vu) - 0.1*sqr(g1)*sqr(vu);

   hermitianize(mass_matrix_St);

   return mass_matrix_St;
}

void CLASSNAME::calculate_MSt()
{
   const auto mass_matrix_St(get_mass_matrix_St());
   fs_diagonalize_hermitian(mass_matrix_St, MSt, ZT);

   if (MSt.minCoeff() < 0.) {
      problems.flag_tachyon("St");
   }

   MSt = sqrt(MSt.cwiseAbs());
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_hh() const
{
   Eigen::Matrix<double,2,2> mass_matrix_hh;

   mass_matrix_hh(0,0) = mHd2 + sqr(Mu) + 0.225*sqr(g1)*sqr(vd) +
      0.375*sqr(g2)*sqr(vd) - 0.075*sqr(g1)*sqr(vu) - 0.125*sqr(g2)*sqr(vu);
   mass_matrix_hh(0,1) = -0.5*BMu - 0.5*BMu - 0.15*vd*vu*sqr(g1) -
      0.25*vd*vu*sqr(g2);
   mass_matrix_hh(1,1) = mHu2 + sqr(Mu) - 0.075*sqr(g1)*sqr(vd) -
      0.125*sqr(g2)*sqr(vd) + 0.225*sqr(g1)*sqr(vu) + 0.375*sqr(g2)*sqr(vu);

   symmetrize(mass_matrix_hh);

   return mass_matrix_hh;
}

void CLASSNAME::calculate_Mhh()
{
   const auto mass_matrix_hh(get_mass_matrix_hh());
   fs_diagonalize_hermitian(mass_matrix_hh, Mhh, ZH);

   if (Mhh.minCoeff() < 0.) {
      problems.flag_tachyon("hh");
   }

   Mhh = sqrt(Mhh.cwiseAbs());
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Ah() const
{
   Eigen::Matrix<double,2,2> mass_matrix_Ah;

   mass_matrix_Ah(0,0) = mHd2 + sqr(Mu) + 0.3872983346207417*g1*g2*std::cos
      (ThetaW())*std::sin(ThetaW())*sqr(vd) + 0.075*sqr(g1)*sqr(vd) + 0.125*sqr(g2)*
      sqr(vd) - 0.075*sqr(g1)*sqr(vu) - 0.125*sqr(g2)*sqr(vu) + 0.25*sqr(g2)*
      sqr(vd)*sqr(std::cos(ThetaW())) + 0.15*sqr(g1)*sqr(vd)*sqr(std::sin(ThetaW()));
   mass_matrix_Ah(0,1) = 0.5*BMu + 0.5*BMu - 0.3872983346207417*g1*
      g2*vd*vu*std::cos(ThetaW())*std::sin(ThetaW()) - 0.25*vd*vu*sqr(g2)*sqr(std::cos(ThetaW(
      ))) - 0.15*vd*vu*sqr(g1)*sqr(std::sin(ThetaW()));
   mass_matrix_Ah(1,1) = mHu2 + sqr(Mu) - 0.075*sqr(g1)*sqr(vd) -
      0.125*sqr(g2)*sqr(vd) + 0.3872983346207417*g1*g2*std::cos(ThetaW())*std::sin(ThetaW
      ())*sqr(vu) + 0.075*sqr(g1)*sqr(vu) + 0.125*sqr(g2)*sqr(vu) + 0.25*sqr(g2
      )*sqr(vu)*sqr(std::cos(ThetaW())) + 0.15*sqr(g1)*sqr(vu)*sqr(std::sin(ThetaW()));

   symmetrize(mass_matrix_Ah);

   return mass_matrix_Ah;
}

void CLASSNAME::calculate_MAh()
{
   const auto mass_matrix_Ah(get_mass_matrix_Ah());
   fs_diagonalize_hermitian(mass_matrix_Ah, MAh, ZA);

   if (MAh.minCoeff() < 0.) {
      problems.flag_tachyon("Ah");
   }

   MAh = sqrt(MAh.cwiseAbs());
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Hpm() const
{
   Eigen::Matrix<double,2,2> mass_matrix_Hpm;

   mass_matrix_Hpm(0,0) = mHd2 + sqr(Mu) + 0.075*sqr(g1)*sqr(vd) +
      0.375*sqr(g2)*sqr(vd) - 0.075*sqr(g1)*sqr(vu) + 0.125*sqr(g2)*sqr(vu);
   mass_matrix_Hpm(0,1) = BMu;
   mass_matrix_Hpm(1,1) = mHu2 + sqr(Mu) - 0.075*sqr(g1)*sqr(vd) +
      0.125*sqr(g2)*sqr(vd) + 0.075*sqr(g1)*sqr(vu) + 0.375*sqr(g2)*sqr(vu);

   hermitianize(mass_matrix_Hpm);

   return mass_matrix_Hpm;
}

void CLASSNAME::calculate_MHpm()
{
   const auto mass_matrix_Hpm(get_mass_matrix_Hpm());
   fs_diagonalize_hermitian(mass_matrix_Hpm, MHpm, ZP);

   if (MHpm.minCoeff() < 0.) {
      problems.flag_tachyon("Hpm");
   }

   MHpm = sqrt(MHpm.cwiseAbs());
}

Eigen::Matrix<double,4,4> CLASSNAME::get_mass_matrix_Chi() const
{
   Eigen::Matrix<double,4,4> mass_matrix_Chi;

   mass_matrix_Chi(0,0) = MassB;
   mass_matrix_Chi(0,1) = 0;
   mass_matrix_Chi(0,2) = -0.3872983346207417*g1*vd;
   mass_matrix_Chi(0,3) = 0.3872983346207417*g1*vu;
   mass_matrix_Chi(1,1) = MassWB;
   mass_matrix_Chi(1,2) = 0.5*g2*vd;
   mass_matrix_Chi(1,3) = -0.5*g2*vu;
   mass_matrix_Chi(2,2) = 0;
   mass_matrix_Chi(2,3) = -Mu;
   mass_matrix_Chi(3,3) = 0;

   symmetrize(mass_matrix_Chi);

   return mass_matrix_Chi;
}

void CLASSNAME::calculate_MChi()
{
   const auto mass_matrix_Chi(get_mass_matrix_Chi());
   fs_diagonalize_symmetric(mass_matrix_Chi, MChi, ZN);
}

Eigen::Matrix<double,2,2> CLASSNAME::get_mass_matrix_Cha() const
{
   Eigen::Matrix<double,2,2> mass_matrix_Cha;

   mass_matrix_Cha(0,0) = MassWB;
   mass_matrix_Cha(0,1) = 0.7071067811865475*g2*vu;
   mass_matrix_Cha(1,0) = 0.7071067811865475*g2*vd;
   mass_matrix_Cha(1,1) = Mu;

   return mass_matrix_Cha;
}

void CLASSNAME::calculate_MCha()
{
   const auto mass_matrix_Cha(get_mass_matrix_Cha());
   fs_svd(mass_matrix_Cha, MCha, UM, UP);
}

double CLASSNAME::get_mass_matrix_VWm() const
{
   const double mass_matrix_VWm = 0.25*sqr(g2)*(sqr(vd) + sqr(vu));

   return mass_matrix_VWm;
}

void CLASSNAME::calculate_MVWm()
{
   const auto mass_matrix_VWm = get_mass_matrix_VWm();
   MVWm = std::sqrt(std::abs(mass_matrix_VWm));
}


double CLASSNAME::get_ewsb_eq_hh_1() const
{
   double result = mHd2*vd + vd*sqr(Mu) - 0.5*vu*BMu - 0.5*vu*BMu +
      0.075*std::pow(vd,3)*sqr(g1) + 0.125*std::pow(vd,3)*sqr(g2) - 0.075*vd*sqr(g1)*
      sqr(vu) - 0.125*vd*sqr(g2)*sqr(vu);

   return result;
}

double CLASSNAME::get_ewsb_eq_hh_2() const
{
   double result = mHu2*vu + vu*sqr(Mu) - 0.5*vd*BMu - 0.5*vd*BMu +
      0.075*std::pow(vu,3)*sqr(g1) + 0.125*std::pow(vu,3)*sqr(g2) - 0.075*vu*sqr(g1)*
      sqr(vd) - 0.125*vu*sqr(g2)*sqr(vd);

   return result;
}

double CLASSNAME::ThetaW() const
{
   return std::atan((0.7745966692414834*g1)/g2);
}

double CLASSNAME::v() const
{
   return 2*std::sqrt(sqr(MVWm)/sqr(g2));
}


std::ostream& operator<<(std::ostream& ostr, const MSSMNoFV_onshell_mass_eigenstates& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace gm2calc
