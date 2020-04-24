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

#include "gm2calc/MSSMNoFV_onshell_physical.hpp"

#include "gm2_log.hpp"
#include "gm2_numerics.hpp"

#include <iostream>
#include <limits>

namespace gm2calc {

namespace {

template<int N>
void convert_symmetric_fermion_mixings_to_hk(Eigen::Array<double, N, 1>& m,
                                             Eigen::Matrix<std::complex<double>, N, N>& z)
{
   for (int i = 0; i < N; i++) {
      if (m(i) < 0.) {
         z.row(i) *= std::complex<double>(0.0,1.0);
         m(i) *= -1;
      }
   }
}

template<int N>
void convert_symmetric_fermion_mixings_to_slha(Eigen::Array<double, N, 1>& m,
                                               Eigen::Matrix<std::complex<double>, N, N>& z)
{
   const double eps = std::numeric_limits<double>::epsilon();

   for (int i = 0; i < N; i++) {
      // check if i'th row contains non-zero imaginary parts
      if (!is_zero(z.row(i).imag().cwiseAbs().maxCoeff(), eps)) {
         z.row(i) *= std::complex<double>(0.0,1.0);
         m(i) *= -1;
#ifdef ENABLE_DEBUG
         if (!is_zero(z.row(i).imag().cwiseAbs().maxCoeff(), eps)) {
            WARNING("Row " << i << " of the following fermion mixing matrix"
                    " contains entries which have non-zero real and imaginary"
                    " parts:\nZ = " << z);
         }
#endif
      }
   }
}

} // anonymous namespace

/**
 * Convert masses and mixing matrices to Haber-Kane convention:
 * Fermion masses are always positive and mixing matrices are allowed
 * to be complex.
 */
void MSSMNoFV_onshell_physical::convert_to_hk()
{
   convert_symmetric_fermion_mixings_to_hk(MChi, ZN);
}

/**
 * Convert masses and mixing matrices to SLHA convention: Fermion
 * mixing matrices are always real and fermion masses are allowed to
 * be negative.
 */
void MSSMNoFV_onshell_physical::convert_to_slha()
{
   convert_symmetric_fermion_mixings_to_slha(MChi, ZN);
}

void MSSMNoFV_onshell_physical::print(std::ostream& ostr) const
{
   ostr << "----------------------------------------\n"
           "pole masses:\n"
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
           "pole mass mixing matrices:\n"
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

}

std::ostream& operator<<(std::ostream& ostr, const MSSMNoFV_onshell_physical& phys_pars)
{
   phys_pars.print(ostr);
   return ostr;
}

} // namespace gm2calc
