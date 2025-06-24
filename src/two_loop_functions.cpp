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
extern void __pole2lfunctions_MOD_ctfss(double&, const double&, const double&, const double&, const double&);
extern void __pole2lfunctions_MOD_ctfsss(double&, const double&, const double&, const double&, const double&);
extern void __pole2lfunctions_MOD_ctfssss(double&, const double&, const double&, const double&, const double&, const double&);
extern void __pole2lfunctions_MOD_ctfssff(double&, const double&, const double&, const double&, const double&, const double&);
extern void __pole2lfunctions_MOD_ctfssfbfb(double&, const double&, const double&, const double&, const double&, const double&);
extern void __pole2lfunctions_MOD_ctffffbs(double&, const double&, const double&, const double&, const double&, const double&);
extern void __pole2lfunctions_MOD_ctfffbfs(double&, const double&, const double&, const double&, const double&, const double&);
extern void __pole2lfunctions_MOD_ctffbfbfbs(double&, const double&, const double&, const double&, const double&, const double&);
extern void __pole2lfunctions_MOD_ctfsv(double&, const double&, const double&);
extern void __pole2lfunctions_MOD_ctffv(double&, const double&, const double&);

// self-energy functions
extern void __pole2lfunctions_MOD_cwfssss(double&, const double&, const double&, const double&, const double&, const double&, const double&);
extern void __pole2lfunctions_MOD_cxfsss(double&, const double&, const double&, const double&, const double&, const double&);
extern void __pole2lfunctions_MOD_cyfssss(double&, const double&, const double&, const double&, const double&, const double&, const double&);
extern void __pole2lfunctions_MOD_czfssss(double&, const double&, const double&, const double&, const double&, const double&, const double&);
extern void __pole2lfunctions_MOD_csfsss(double&, const double&, const double&, const double&, const double&, const double&);
extern void __pole2lfunctions_MOD_cufssss(double&, const double&, const double&, const double&, const double&, const double&, const double&);
extern void __pole2lfunctions_MOD_cvfsssss(double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&);
extern void __pole2lfunctions_MOD_cmfsssss(double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&);
extern void __pole2lfunctions_MOD_cwfsssv(double&, const double&, const double&, const double&);
extern void __pole2lfunctions_MOD_cmfssssv(double&, const double&, const double&, const double&, const double&);
extern void __pole2lfunctions_MOD_cwfssff(double&, const double&, const double&, const double&, const double&, const double&, const double&);
extern void __pole2lfunctions_MOD_cwfssfbfb(double&, const double&, const double&, const double&, const double&, const double&, const double&);
extern void __pole2lfunctions_MOD_cmffbfbfbfbs(double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&);
extern void __pole2lfunctions_MOD_cmfffbfbfs(double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&);
extern void __pole2lfunctions_MOD_cmfffbffbs(double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&);
extern void __pole2lfunctions_MOD_cmffffbfbs(double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&);
extern void __pole2lfunctions_MOD_cmfffffs(double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&);
extern void __pole2lfunctions_MOD_cmfsfbsfbfb(double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&);
extern void __pole2lfunctions_MOD_cmfsfsfbf(double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&);
extern void __pole2lfunctions_MOD_cmfsfsffb(double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&);
extern void __pole2lfunctions_MOD_cvfsssfbfb(double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&);
extern void __pole2lfunctions_MOD_cvfsssff(double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&);
extern void __pole2lfunctions_MOD_cvffbfbfbfbs(double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&);
extern void __pole2lfunctions_MOD_cvffbffbfs(double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&);
extern void __pole2lfunctions_MOD_cvffbfffbs(double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&);
extern void __pole2lfunctions_MOD_cvfffbfbfs(double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&);
extern void __pole2lfunctions_MOD_cvffffbfbs(double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&);
extern void __pole2lfunctions_MOD_cvfffffs(double&, const double&, const double&, const double&, const double&, const double&, const double&, const double&);
extern void __pole2lfunctions_MOD_cgfffv(double&, const double&, const double&, const double&, const double&);
extern void __pole2lfunctions_MOD_cgffbfbv(double&, const double&, const double&, const double&, const double&);
}

namespace flexiblesusy {
namespace two_loop_functions {

double TfSS(double x, double y, double z, double scl2) noexcept
{
   double out = 0.;
   __pole2lfunctions_MOD_ctfss(out, x, y, z, scl2);
   return out;
}

double TfSSS(double x, double y, double z, double scl2) noexcept
{
   double out = 0.;
   __pole2lfunctions_MOD_ctfsss(out, x, y, z, scl2);
   return out;
}

double TfSSSS(double x, double y, double z, double u, double scl2) noexcept
{
   double out = 0.;
   __pole2lfunctions_MOD_ctfssss(out, x, y, z, u, scl2);
   return out;
}

double TfSSFF(double x, double y, double z, double u, double scl2) noexcept
{
   double out = 0.;
   __pole2lfunctions_MOD_ctfssff(out, x, y, z, u, scl2);
   return out;
}

double TfSSFbFb(double x, double y, double z, double u, double scl2) noexcept
{
   double out = 0.;
   __pole2lfunctions_MOD_ctfssfbfb(out, x, y, z, u, scl2);
   return out;
}

double TfFFFbS(double x, double y, double z, double u, double scl2) noexcept
{
   double out = 0.;
   __pole2lfunctions_MOD_ctffffbs(out, x, y, z, u, scl2);
   return out;
}

double TfFFbFS(double x, double y, double z, double u, double scl2) noexcept
{
   double out = 0.;
   __pole2lfunctions_MOD_ctfffbfs(out, x, y, z, u, scl2);
   return out;
}

double TfFbFbFbS(double x, double y, double z, double u, double scl2) noexcept
{
   double out = 0.;
   __pole2lfunctions_MOD_ctffbfbfbs(out, x, y, z, u, scl2);
   return out;
}

double TfSV(double x, double scl2) noexcept
{
   double out = 0.;
   __pole2lfunctions_MOD_ctfsv(out, x, scl2);
   return out;
}

double TfFV(double x, double scl2) noexcept
{
   double out = 0.;
   __pole2lfunctions_MOD_ctffv(out, x, scl2);
   return out;
}

// self energy functions

double WfSSSS(double p2, double x, double y, double z, double u, double scl2) noexcept
{
   double out = 0.;
   __pole2lfunctions_MOD_cwfssss(out, p2, x, y, z, u, scl2);
   return out;
}

double XfSSS(double p2, double x, double y, double z, double scl2) noexcept
{
   double out = 0.;
   __pole2lfunctions_MOD_cxfsss(out, p2, x, y, z, scl2);
   return out;
}

double YfSSSS(double p2, double x, double y, double z, double u, double scl2) noexcept
{
   double out = 0.;
   __pole2lfunctions_MOD_cyfssss(out, p2, x, y, z, u, scl2);
   return out;
}

double ZfSSSS(double p2, double x, double y, double z, double u, double scl2) noexcept
{
   double out = 0.;
   __pole2lfunctions_MOD_czfssss(out, p2, x, y, z, u, scl2);
   return out;
}

double SfSSS(double p2, double x, double y, double z, double scl2) noexcept
{
   double out = 0.;
   __pole2lfunctions_MOD_csfsss(out, p2, x, y, z, scl2);
   return out;
}

double UfSSSS(double p2, double x, double y, double z, double u, double scl2) noexcept
{
   double out = 0.;
   __pole2lfunctions_MOD_cufssss(out, p2, x, y, z, u, scl2);
   return out;
}

double VfSSSSS(double p2, double x, double y, double z, double u, double v, double scl2) noexcept
{
   double out = 0.;
   __pole2lfunctions_MOD_cvfsssss(out, p2, x, y, z, u, v, scl2);
   return out;
}

double MfSSSSS(double p2, double x, double y, double z, double u, double v, double scl2) noexcept
{
   double out = 0.;
   __pole2lfunctions_MOD_cmfsssss(out, p2, x, y, z, u, v, scl2);
   return out;
}

double WfSSSV(double p2, double x, double scl2) noexcept
{
   double out = 0.;
   __pole2lfunctions_MOD_cwfsssv(out, p2, x, scl2);
   return out;
}

double MfSSSSV(double p2, double x, double y, double scl2) noexcept
{
   double out = 0.;
   __pole2lfunctions_MOD_cmfssssv(out, p2, x, y, scl2);
   return out;
}

double WfSSFF(double p2, double x, double y, double z, double u, double scl2) noexcept
{
   double out = 0.;
   __pole2lfunctions_MOD_cwfssff(out, p2, x, y, z, u, scl2);
   return out;
}

double WfSSFbFb(double p2, double x, double y, double z, double u, double scl2) noexcept
{
   double out = 0.;
   __pole2lfunctions_MOD_cwfssfbfb(out, p2, x, y, z, u, scl2);
   return out;
}

double MfFbFbFbFbS(double p2, double x, double y, double z, double u, double v, double scl2) noexcept
{
   double out = 0.;
   __pole2lfunctions_MOD_cmffbfbfbfbs(out, p2, x, y, z, u, v, scl2);
   return out;
}

double MfFFbFbFS(double p2, double x, double y, double z, double u, double v, double scl2) noexcept
{
   double out = 0.;
   __pole2lfunctions_MOD_cmfffbfbfs(out, p2, x, y, z, u, v, scl2);
   return out;
}

double MfFFbFFbS(double p2, double x, double y, double z, double u, double v, double scl2) noexcept
{
   double out = 0.;
   __pole2lfunctions_MOD_cmfffbffbs(out, p2, x, y, z, u, v, scl2);
   return out;
}

double MfFFFbFbS(double p2, double x, double y, double z, double u, double v, double scl2) noexcept
{
   double out = 0.;
   __pole2lfunctions_MOD_cmffffbfbs(out, p2, x, y, z, u, v, scl2);
   return out;
}

double MfFFFFS(double p2, double x, double y, double z, double u, double v, double scl2) noexcept
{
   double out = 0.;
   __pole2lfunctions_MOD_cmfffffs(out, p2, x, y, z, u, v, scl2);
   return out;
}

double MfSFbSFbFb(double p2, double x, double y, double z, double u, double v, double scl2) noexcept
{
   double out = 0.;
   __pole2lfunctions_MOD_cmfsfbsfbfb(out, p2, x, y, z, u, v, scl2);
   return out;
}

double MfSFSFbF(double p2, double x, double y, double z, double u, double v, double scl2) noexcept
{
   double out = 0.;
   __pole2lfunctions_MOD_cmfsfsfbf(out, p2, x, y, z, u, v, scl2);
   return out;
}

double MfSFSFFb(double p2, double x, double y, double z, double u, double v, double scl2) noexcept
{
   double out = 0.;
   __pole2lfunctions_MOD_cmfsfsffb(out, p2, x, y, z, u, v, scl2);
   return out;
}

double VfSSSFbFb(double p2, double x, double y, double z, double u, double v, double scl2) noexcept
{
   double out = 0.;
   __pole2lfunctions_MOD_cvfsssfbfb(out, p2, x, y, z, u, v, scl2);
   return out;
}

double VfSSSFF(double p2, double x, double y, double z, double u, double v, double scl2) noexcept
{
   double out = 0.;
   __pole2lfunctions_MOD_cvfsssff(out, p2, x, y, z, u, v, scl2);
   return out;
}

double VfFbFbFbFbS(double p2, double x, double y, double z, double u, double v, double scl2) noexcept
{
   double out = 0.;
   __pole2lfunctions_MOD_cvffbfbfbfbs(out, p2, x, y, z, u, v, scl2);
   return out;
}

double VfFbFFbFS(double p2, double x, double y, double z, double u, double v, double scl2) noexcept
{
   double out = 0.;
   __pole2lfunctions_MOD_cvffbffbfs(out, p2, x, y, z, u, v, scl2);
   return out;
}

double VfFbFFFbS(double p2, double x, double y, double z, double u, double v, double scl2) noexcept
{
   double out = 0.;
   __pole2lfunctions_MOD_cvffbfffbs(out, p2, x, y, z, u, v, scl2);
   return out;
}

double VfFFbFbFS(double p2, double x, double y, double z, double u, double v, double scl2) noexcept
{
   double out = 0.;
   __pole2lfunctions_MOD_cvfffbfbfs(out, p2, x, y, z, u, v, scl2);
   return out;
}

double VfFFFbFbS(double p2, double x, double y, double z, double u, double v, double scl2) noexcept
{
   double out = 0.;
   __pole2lfunctions_MOD_cvffffbfbs(out, p2, x, y, z, u, v, scl2);
   return out;
}

double VfFFFFS(double p2, double x, double y, double z, double u, double v, double scl2) noexcept
{
   double out = 0.;
   __pole2lfunctions_MOD_cvfffffs(out, p2, x, y, z, u, v, scl2);
   return out;
}

double GfFFV(double p2, double x, double y, double scl2) noexcept
{
   double out = 0.;
   __pole2lfunctions_MOD_cgfffv(out, p2, x, y, scl2);
   return out;
}

double GfFbFbV(double p2, double x, double y, double scl2) noexcept
{
   double out = 0.;
   __pole2lfunctions_MOD_cgffbfbv(out, p2, x, y, scl2);
   return out;
}

} // namespace two_loop_functions
} // namespace flexiblesusy
