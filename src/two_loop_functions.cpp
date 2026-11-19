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

#include "two_loop_functions.hpp"

extern "C" {
// tadpole functions
extern void pole2l_ctfss(double&, const double&, const double&, const double&,
                         const double&);
extern void pole2l_ctfsss(double&, const double&, const double&, const double&,
                          const double&);
extern void pole2l_ctfssss(double&, const double&, const double&, const double&,
                           const double&, const double&);
extern void pole2l_ctfssff(double&, const double&, const double&, const double&,
                           const double&, const double&);
extern void pole2l_ctfssfbfb(double&, const double&, const double&,
                             const double&, const double&, const double&);
extern void pole2l_ctffffbs(double&, const double&, const double&,
                            const double&, const double&, const double&);
extern void pole2l_ctfffbfs(double&, const double&, const double&,
                            const double&, const double&, const double&);
extern void pole2l_ctffbfbfbs(double&, const double&, const double&,
                              const double&, const double&, const double&);
extern void pole2l_ctfsv(double&, const double&, const double&);
extern void pole2l_ctffv(double&, const double&, const double&);
extern void pole2l_ctffvmsbar(double&, const double&, const double&);

// self-energy functions
extern void pole2l_cwfssss(double&, const double&, const double&, const double&,
                           const double&, const double&, const double&);
extern void pole2l_cxfsss(double&, const double&, const double&, const double&,
                          const double&, const double&);
extern void pole2l_cyfssss(double&, const double&, const double&, const double&,
                           const double&, const double&, const double&);
extern void pole2l_czfssss(double&, const double&, const double&, const double&,
                           const double&, const double&, const double&);
extern void pole2l_csfsss(double&, const double&, const double&, const double&,
                          const double&, const double&);
extern void pole2l_cufssss(double&, const double&, const double&, const double&,
                           const double&, const double&, const double&);
extern void pole2l_cvfsssss(double&, const double&, const double&,
                            const double&, const double&, const double&,
                            const double&, const double&);
extern void pole2l_cmfsssss(double&, const double&, const double&,
                            const double&, const double&, const double&,
                            const double&, const double&);
extern void pole2l_cwfsssv(double&, const double&, const double&,
                           const double&);
extern void pole2l_cmfssssv(double&, const double&, const double&,
                            const double&, const double&);
extern void pole2l_cwfssff(double&, const double&, const double&, const double&,
                           const double&, const double&, const double&);
extern void pole2l_cwfssfbfb(double&, const double&, const double&,
                             const double&, const double&, const double&,
                             const double&);
extern void pole2l_cmffbfbfbfbs(double&, const double&, const double&,
                                const double&, const double&, const double&,
                                const double&, const double&);
extern void pole2l_cmfffbfbfs(double&, const double&, const double&,
                              const double&, const double&, const double&,
                              const double&, const double&);
extern void pole2l_cmfffbffbs(double&, const double&, const double&,
                              const double&, const double&, const double&,
                              const double&, const double&);
extern void pole2l_cmffffbfbs(double&, const double&, const double&,
                              const double&, const double&, const double&,
                              const double&, const double&);
extern void pole2l_cmfffffs(double&, const double&, const double&,
                            const double&, const double&, const double&,
                            const double&, const double&);
extern void pole2l_cmfsfbsfbfb(double&, const double&, const double&,
                               const double&, const double&, const double&,
                               const double&, const double&);
extern void pole2l_cmfsfsfbf(double&, const double&, const double&,
                             const double&, const double&, const double&,
                             const double&, const double&);
extern void pole2l_cmfsfsffb(double&, const double&, const double&,
                             const double&, const double&, const double&,
                             const double&, const double&);
extern void pole2l_cvfsssfbfb(double&, const double&, const double&,
                              const double&, const double&, const double&,
                              const double&, const double&);
extern void pole2l_cvfsssff(double&, const double&, const double&,
                            const double&, const double&, const double&,
                            const double&, const double&);
extern void pole2l_cvffbfbfbfbs(double&, const double&, const double&,
                                const double&, const double&, const double&,
                                const double&, const double&);
extern void pole2l_cvffbffbfs(double&, const double&, const double&,
                              const double&, const double&, const double&,
                              const double&, const double&);
extern void pole2l_cvffbfffbs(double&, const double&, const double&,
                              const double&, const double&, const double&,
                              const double&, const double&);
extern void pole2l_cvfffbfbfs(double&, const double&, const double&,
                              const double&, const double&, const double&,
                              const double&, const double&);
extern void pole2l_cvffffbfbs(double&, const double&, const double&,
                              const double&, const double&, const double&,
                              const double&, const double&);
extern void pole2l_cvfffffs(double&, const double&, const double&,
                            const double&, const double&, const double&,
                            const double&, const double&);
extern void pole2l_cgfffv(double&, const double&, const double&, const double&,
                          const double&);
extern void pole2l_cgfffvmsbar(double&, const double&, const double&,
                               const double&, const double&);
extern void pole2l_cgffbfbv(double&, const double&, const double&,
                            const double&, const double&);
extern void pole2l_cgffbfbvmsbar(double&, const double&, const double&,
                                 const double&, const double&);
}

namespace flexiblesusy {
namespace two_loop_functions {

double TfSS(double x, double y, double z, double scl2) noexcept
{
   double out = 0.;
   pole2l_ctfss(out, x, y, z, scl2);
   return out;
}

double TfSSS(double x, double y, double z, double scl2) noexcept
{
   double out = 0.;
   pole2l_ctfsss(out, x, y, z, scl2);
   return out;
}

double TfSSSS(double x, double y, double z, double u, double scl2) noexcept
{
   double out = 0.;
   pole2l_ctfssss(out, x, y, z, u, scl2);
   return out;
}

double TfSSFF(double x, double y, double z, double u, double scl2) noexcept
{
   double out = 0.;
   pole2l_ctfssff(out, x, y, z, u, scl2);
   return out;
}

double TfSSFbFb(double x, double y, double z, double u, double scl2) noexcept
{
   double out = 0.;
   pole2l_ctfssfbfb(out, x, y, z, u, scl2);
   return out;
}

double TfFFFbS(double x, double y, double z, double u, double scl2) noexcept
{
   double out = 0.;
   pole2l_ctffffbs(out, x, y, z, u, scl2);
   return out;
}

double TfFFbFS(double x, double y, double z, double u, double scl2) noexcept
{
   double out = 0.;
   pole2l_ctfffbfs(out, x, y, z, u, scl2);
   return out;
}

double TfFbFbFbS(double x, double y, double z, double u, double scl2) noexcept
{
   double out = 0.;
   pole2l_ctffbfbfbs(out, x, y, z, u, scl2);
   return out;
}

double TfSV(double x, double scl2) noexcept
{
   double out = 0.;
   pole2l_ctfsv(out, x, scl2);
   return out;
}

double TfFV(double x, double scl2, Scheme scheme) noexcept
{
   double out = 0.;
   if (scheme == Scheme::MSbar) {
      pole2l_ctffvmsbar(out, x, scl2);
   } else {
      pole2l_ctffv(out, x, scl2);
   }
   return out;
}

// self energy functions

double WfSSSS(double p2, double x, double y, double z, double u, double scl2) noexcept
{
   double out = 0.;
   pole2l_cwfssss(out, p2, x, y, z, u, scl2);
   return out;
}

double XfSSS(double p2, double x, double y, double z, double scl2) noexcept
{
   double out = 0.;
   pole2l_cxfsss(out, p2, x, y, z, scl2);
   return out;
}

double YfSSSS(double p2, double x, double y, double z, double u, double scl2) noexcept
{
   double out = 0.;
   pole2l_cyfssss(out, p2, x, y, z, u, scl2);
   return out;
}

double ZfSSSS(double p2, double x, double y, double z, double u, double scl2) noexcept
{
   double out = 0.;
   pole2l_czfssss(out, p2, x, y, z, u, scl2);
   return out;
}

double SfSSS(double p2, double x, double y, double z, double scl2) noexcept
{
   double out = 0.;
   pole2l_csfsss(out, p2, x, y, z, scl2);
   return out;
}

double UfSSSS(double p2, double x, double y, double z, double u, double scl2) noexcept
{
   double out = 0.;
   pole2l_cufssss(out, p2, x, y, z, u, scl2);
   return out;
}

double VfSSSSS(double p2, double x, double y, double z, double u, double v, double scl2) noexcept
{
   double out = 0.;
   pole2l_cvfsssss(out, p2, x, y, z, u, v, scl2);
   return out;
}

double MfSSSSS(double p2, double x, double y, double z, double u, double v, double scl2) noexcept
{
   double out = 0.;
   pole2l_cmfsssss(out, p2, x, y, z, u, v, scl2);
   return out;
}

double WfSSSV(double p2, double x, double scl2) noexcept
{
   double out = 0.;
   pole2l_cwfsssv(out, p2, x, scl2);
   return out;
}

double MfSSSSV(double p2, double x, double y, double scl2) noexcept
{
   double out = 0.;
   pole2l_cmfssssv(out, p2, x, y, scl2);
   return out;
}

double WfSSFF(double p2, double x, double y, double z, double u, double scl2) noexcept
{
   double out = 0.;
   pole2l_cwfssff(out, p2, x, y, z, u, scl2);
   return out;
}

double WfSSFbFb(double p2, double x, double y, double z, double u, double scl2) noexcept
{
   double out = 0.;
   pole2l_cwfssfbfb(out, p2, x, y, z, u, scl2);
   return out;
}

double MfFbFbFbFbS(double p2, double x, double y, double z, double u, double v, double scl2) noexcept
{
   double out = 0.;
   pole2l_cmffbfbfbfbs(out, p2, x, y, z, u, v, scl2);
   return out;
}

double MfFFbFbFS(double p2, double x, double y, double z, double u, double v, double scl2) noexcept
{
   double out = 0.;
   pole2l_cmfffbfbfs(out, p2, x, y, z, u, v, scl2);
   return out;
}

double MfFFbFFbS(double p2, double x, double y, double z, double u, double v, double scl2) noexcept
{
   double out = 0.;
   pole2l_cmfffbffbs(out, p2, x, y, z, u, v, scl2);
   return out;
}

double MfFFFbFbS(double p2, double x, double y, double z, double u, double v, double scl2) noexcept
{
   double out = 0.;
   pole2l_cmffffbfbs(out, p2, x, y, z, u, v, scl2);
   return out;
}

double MfFFFFS(double p2, double x, double y, double z, double u, double v, double scl2) noexcept
{
   double out = 0.;
   pole2l_cmfffffs(out, p2, x, y, z, u, v, scl2);
   return out;
}

double MfSFbSFbFb(double p2, double x, double y, double z, double u, double v, double scl2) noexcept
{
   double out = 0.;
   pole2l_cmfsfbsfbfb(out, p2, x, y, z, u, v, scl2);
   return out;
}

double MfSFSFbF(double p2, double x, double y, double z, double u, double v, double scl2) noexcept
{
   double out = 0.;
   pole2l_cmfsfsfbf(out, p2, x, y, z, u, v, scl2);
   return out;
}

double MfSFSFFb(double p2, double x, double y, double z, double u, double v, double scl2) noexcept
{
   double out = 0.;
   pole2l_cmfsfsffb(out, p2, x, y, z, u, v, scl2);
   return out;
}

double VfSSSFbFb(double p2, double x, double y, double z, double u, double v, double scl2) noexcept
{
   double out = 0.;
   pole2l_cvfsssfbfb(out, p2, x, y, z, u, v, scl2);
   return out;
}

double VfSSSFF(double p2, double x, double y, double z, double u, double v, double scl2) noexcept
{
   double out = 0.;
   pole2l_cvfsssff(out, p2, x, y, z, u, v, scl2);
   return out;
}

double VfFbFbFbFbS(double p2, double x, double y, double z, double u, double v, double scl2) noexcept
{
   double out = 0.;
   pole2l_cvffbfbfbfbs(out, p2, x, y, z, u, v, scl2);
   return out;
}

double VfFbFFbFS(double p2, double x, double y, double z, double u, double v, double scl2) noexcept
{
   double out = 0.;
   pole2l_cvffbffbfs(out, p2, x, y, z, u, v, scl2);
   return out;
}

double VfFbFFFbS(double p2, double x, double y, double z, double u, double v, double scl2) noexcept
{
   double out = 0.;
   pole2l_cvffbfffbs(out, p2, x, y, z, u, v, scl2);
   return out;
}

double VfFFbFbFS(double p2, double x, double y, double z, double u, double v, double scl2) noexcept
{
   double out = 0.;
   pole2l_cvfffbfbfs(out, p2, x, y, z, u, v, scl2);
   return out;
}

double VfFFFbFbS(double p2, double x, double y, double z, double u, double v, double scl2) noexcept
{
   double out = 0.;
   pole2l_cvffffbfbs(out, p2, x, y, z, u, v, scl2);
   return out;
}

double VfFFFFS(double p2, double x, double y, double z, double u, double v, double scl2) noexcept
{
   double out = 0.;
   pole2l_cvfffffs(out, p2, x, y, z, u, v, scl2);
   return out;
}

double GfFFV(double p2, double x, double y, double scl2, Scheme scheme) noexcept
{
   double out = 0.;
   if (scheme == Scheme::MSbar) {
      pole2l_cgfffvmsbar(out, p2, x, y, scl2);
   } else {
      pole2l_cgfffv(out, p2, x, y, scl2);
   }
   return out;
}

double GfFbFbV(double p2, double x, double y, double scl2, Scheme scheme) noexcept
{
   double out = 0.;
   if (scheme == Scheme::MSbar) {
      pole2l_cgffbfbvmsbar(out, p2, x, y, scl2);
   } else {
      pole2l_cgffbfbv(out, p2, x, y, scl2);
   }
   return out;
}

} // namespace two_loop_functions
} // namespace flexiblesusy
