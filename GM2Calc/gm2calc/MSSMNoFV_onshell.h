/* ====================================================================
 * This file is part of GM2Calc.
 *
 * GM2Calc is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * GM2Calc is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GM2Calc.  If not, see
 * <http://www.gnu.org/licenses/>.
 * ==================================================================== */

#ifndef GM2_MSSMNoFV_ONSHELL_H
#define GM2_MSSMNoFV_ONSHELL_H

/**
 * @file MSSMNoFV_onshell.h
 * @brief contains declarations of C interface functions for the model
 *
 * This file contains the declarations for the C interface functions
 * used to set and retrieve the model parameters and masses.
 */

#ifdef __cplusplus
extern "C" {
#endif

/** handle */
struct MSSMNoFV_onshell;
typedef struct MSSMNoFV_onshell MSSMNoFV_onshell;

/** error codes */
typedef enum {
   gm2calc_NoError = 0,
   gm2calc_InvalidInput,
   gm2calc_PhysicalProblem,
   gm2calc_UnknownError
} gm2calc_error;

/** translate error codes into a string */
const char* gm2calc_error_str(gm2calc_error);

/** allocate new MSSMNoFV model */
MSSMNoFV_onshell* gm2calc_mssmnofv_new(void);

/** delete MSSMNoFV model */
void gm2calc_mssmnofv_free(MSSMNoFV_onshell*);

/** set alpha_em(MZ) */
void gm2calc_mssmnofv_set_alpha_MZ(MSSMNoFV_onshell*, double);

/** set alpha_em(0) in the Thomson limit */
void gm2calc_mssmnofv_set_alpha_thompson(MSSMNoFV_onshell*, double);

/** set soft-breaking trilinear coupling Ae(i,k) */
void gm2calc_mssmnofv_set_Ae(MSSMNoFV_onshell*, unsigned, unsigned, double);

/** set soft-breaking trilinear coupling Au(i,k) */
void gm2calc_mssmnofv_set_Au(MSSMNoFV_onshell*, unsigned, unsigned, double);

/** set soft-breaking trilinear coupling Ad(i,k) */
void gm2calc_mssmnofv_set_Ad(MSSMNoFV_onshell*, unsigned, unsigned, double);

/** set gauge coupling g3 */
void gm2calc_mssmnofv_set_g3(MSSMNoFV_onshell*, double);

/** set bino mass */
void gm2calc_mssmnofv_set_MassB(MSSMNoFV_onshell*, double);

/** set wino mass */
void gm2calc_mssmnofv_set_MassWB(MSSMNoFV_onshell*, double);

/** set gluino mass */
void gm2calc_mssmnofv_set_MassG(MSSMNoFV_onshell*, double);

/** set soft-breaking squared mass parameter mq2(i,k) */
void gm2calc_mssmnofv_set_mq2(MSSMNoFV_onshell*, unsigned, unsigned, double);

/** set soft-breaking squared mass parameter mu2(i,k) */
void gm2calc_mssmnofv_set_mu2(MSSMNoFV_onshell*, unsigned, unsigned, double);

/** set soft-breaking squared mass parameter md2(i,k) */
void gm2calc_mssmnofv_set_md2(MSSMNoFV_onshell*, unsigned, unsigned, double);

/** set soft-breaking squared mass parameter ml2(i,k) */
void gm2calc_mssmnofv_set_ml2(MSSMNoFV_onshell*, unsigned, unsigned, double);

/** set soft-breaking squared mass parameter me2(i,k) */
void gm2calc_mssmnofv_set_me2(MSSMNoFV_onshell*, unsigned, unsigned, double);

/** set soft-breaking squared mass parameter Mu parameter */
void gm2calc_mssmnofv_set_Mu(MSSMNoFV_onshell*, double);

/** set tan(beta) */
void gm2calc_mssmnofv_set_TB(MSSMNoFV_onshell*, double);

/** set renormalization scale */
void gm2calc_mssmnofv_set_scale(MSSMNoFV_onshell*, double);

/** set CP-odd Higgs pole mass */
void gm2calc_mssmnofv_set_MAh_pole(MSSMNoFV_onshell*, double);

/** set Z boson pole mass */
void gm2calc_mssmnofv_set_MZ_pole(MSSMNoFV_onshell*, double);

/** set W boson pole mass */
void gm2calc_mssmnofv_set_MW_pole(MSSMNoFV_onshell*, double);

/** set top-quark pole mass */
void gm2calc_mssmnofv_set_MT_pole(MSSMNoFV_onshell*, double);

/** set MS-bar bottom-quark mass mb at the scale mb */
void gm2calc_mssmnofv_set_MB_running(MSSMNoFV_onshell*, double);

/** set tau-lepton pole mass */
void gm2calc_mssmnofv_set_ML_pole(MSSMNoFV_onshell*, double);

/** set muon pole mass */
void gm2calc_mssmnofv_set_MM_pole(MSSMNoFV_onshell*, double);

/** set smuon pole masses */
void gm2calc_mssmnofv_set_MSm_pole(MSSMNoFV_onshell*, unsigned, double);

/** set muon sneutrino pole masses */
void gm2calc_mssmnofv_set_MSvmL_pole(MSSMNoFV_onshell*, double);

/** set chargino pole masses */
void gm2calc_mssmnofv_set_MCha_pole(MSSMNoFV_onshell*, unsigned, double);

/** set neutralino pole masses */
void gm2calc_mssmnofv_set_MChi_pole(MSSMNoFV_onshell*, unsigned, double);

/** enable/disable verbose output */
void gm2calc_mssmnofv_set_verbose_output(MSSMNoFV_onshell*, int);

/** get Ae(i,k) DR-bar */
double gm2calc_mssmnofv_get_Ae(const MSSMNoFV_onshell*, unsigned, unsigned);

/** get Ad(i,k) */
double gm2calc_mssmnofv_get_Ad(const MSSMNoFV_onshell*, unsigned, unsigned);

/** get Au(i,k) */
double gm2calc_mssmnofv_get_Au(const MSSMNoFV_onshell*, unsigned, unsigned);

/** get electromagnetic gauge coupling at MZ w/o hadronic corrections */
double gm2calc_mssmnofv_get_EL(const MSSMNoFV_onshell*);

/** get electromagnetic gauge coupling in Thomson limit */
double gm2calc_mssmnofv_get_EL0(const MSSMNoFV_onshell*);

/** get Hypercharge gauge coupling (not GUT normalized) */
double gm2calc_mssmnofv_get_gY(const MSSMNoFV_onshell*);

/** get Hypercharge gauge coupling (GUT normalized) */
double gm2calc_mssmnofv_get_g1(const MSSMNoFV_onshell*);

/** get left gauge coupling */
double gm2calc_mssmnofv_get_g2(const MSSMNoFV_onshell*);

/** get strong gauge coupling */
double gm2calc_mssmnofv_get_g3(const MSSMNoFV_onshell*);

/** get tan(beta) */
double gm2calc_mssmnofv_get_TB(const MSSMNoFV_onshell*);

/** get soft-breaking on-shell bino mass parameter */
double gm2calc_mssmnofv_get_MassB(const MSSMNoFV_onshell*);

/** get soft-breaking on-shell wino mass parameter */
double gm2calc_mssmnofv_get_MassWB(const MSSMNoFV_onshell*);

/** get soft-breaking gluino mass parameter */
double gm2calc_mssmnofv_get_MassG(const MSSMNoFV_onshell*);

/** get on-shell superpotential mu parameter */
double gm2calc_mssmnofv_get_Mu(const MSSMNoFV_onshell*);

/** get left-handed up-Squark soft-breaking squared mass */
double gm2calc_mssmnofv_get_mq2(const MSSMNoFV_onshell*, unsigned, unsigned);

/** get right-handed down-Squark soft-breaking squared mass */
double gm2calc_mssmnofv_get_md2(const MSSMNoFV_onshell*, unsigned, unsigned);

/** get right-handed up-Squark soft-breaking squared mass */
double gm2calc_mssmnofv_get_mu2(const MSSMNoFV_onshell*, unsigned, unsigned);

/** get left-handed up-lepton soft-breaking on-shell squared mass */
double gm2calc_mssmnofv_get_ml2(const MSSMNoFV_onshell*, unsigned, unsigned);

/** get right-handed down-lepton soft-breaking on-shell squared mass */
double gm2calc_mssmnofv_get_me2(const MSSMNoFV_onshell*, unsigned, unsigned);

/** get vacuum expectation value */
double gm2calc_mssmnofv_get_vev(const MSSMNoFV_onshell*);

/** get renormalization scale */
double gm2calc_mssmnofv_get_scale(const MSSMNoFV_onshell*);

/** get W boson pole mass */
double gm2calc_mssmnofv_get_MW(const MSSMNoFV_onshell*);

/** get Z boson pole mass */
double gm2calc_mssmnofv_get_MZ(const MSSMNoFV_onshell*);

/** get electron mass */
double gm2calc_mssmnofv_get_ME(const MSSMNoFV_onshell*);

/** get muon pole mass */
double gm2calc_mssmnofv_get_MM(const MSSMNoFV_onshell*);

/** get tau mass */
double gm2calc_mssmnofv_get_ML(const MSSMNoFV_onshell*);

/** get up-quark mass */
double gm2calc_mssmnofv_get_MU(const MSSMNoFV_onshell*);

/** get charm-quark mass */
double gm2calc_mssmnofv_get_MC(const MSSMNoFV_onshell*);

/** get top-quark mass */
double gm2calc_mssmnofv_get_MT(const MSSMNoFV_onshell*);

/** get down-quark mass */
double gm2calc_mssmnofv_get_MD(const MSSMNoFV_onshell*);

/** get strange-quark mass */
double gm2calc_mssmnofv_get_MS(const MSSMNoFV_onshell*);

/** get bottom-quark DR-bar mass mb(MZ) */
double gm2calc_mssmnofv_get_MB(const MSSMNoFV_onshell*);

/** get bottom-quark MS-bar mass mb(mb) */
double gm2calc_mssmnofv_get_MBMB(const MSSMNoFV_onshell*);

/** get CP-odd Higgs mass */
double gm2calc_mssmnofv_get_MAh(const MSSMNoFV_onshell*);

/** get CP-even Higgs masses */
double gm2calc_mssmnofv_get_Mhh(const MSSMNoFV_onshell*, unsigned);

/** get chargino pole masses */
double gm2calc_mssmnofv_get_MCha(const MSSMNoFV_onshell*, unsigned);

/** get chargino pole mixing matrix */
double gm2calc_mssmnofv_get_UM(const MSSMNoFV_onshell*, unsigned, unsigned, double*);

/** get chargino pole mixing matrix */
double gm2calc_mssmnofv_get_UP(const MSSMNoFV_onshell*, unsigned, unsigned, double*);

/** get neutralino pole masses */
double gm2calc_mssmnofv_get_MChi(const MSSMNoFV_onshell*, unsigned);

/** get neutralino pole mixing matrix */
double gm2calc_mssmnofv_get_ZN(const MSSMNoFV_onshell*, unsigned, unsigned, double*);

/** get selectron masses */
double gm2calc_mssmnofv_get_MSe(const MSSMNoFV_onshell*, unsigned);

/** get electron sneutrino mass */
double gm2calc_mssmnofv_get_MSveL(const MSSMNoFV_onshell*);

/** get smuon pole masses */
double gm2calc_mssmnofv_get_MSm(const MSSMNoFV_onshell*, unsigned);

/** get muon sneutrino pole mass */
double gm2calc_mssmnofv_get_MSvmL(const MSSMNoFV_onshell*);

/** get stau masses */
double gm2calc_mssmnofv_get_MStau(const MSSMNoFV_onshell*, unsigned);

/** get tau sneutrino mass */
double gm2calc_mssmnofv_get_MSvtL(const MSSMNoFV_onshell*);

/** get sup masses */
double gm2calc_mssmnofv_get_MSu(const MSSMNoFV_onshell*, unsigned);

/** get sdown masses */
double gm2calc_mssmnofv_get_MSd(const MSSMNoFV_onshell*, unsigned);

/** get scharm masses */
double gm2calc_mssmnofv_get_MSc(const MSSMNoFV_onshell*, unsigned);

/** get sstrange masses */
double gm2calc_mssmnofv_get_MSs(const MSSMNoFV_onshell*, unsigned);

/** get stop masses */
double gm2calc_mssmnofv_get_MSt(const MSSMNoFV_onshell*, unsigned);

/** get sbottom mass */
double gm2calc_mssmnofv_get_MSb(const MSSMNoFV_onshell*, unsigned);

/** get selectron mixing matrix */
double gm2calc_mssmnofv_get_USe(const MSSMNoFV_onshell*, unsigned, unsigned);

/** get smuon pole mixing matrix */
double gm2calc_mssmnofv_get_USm(const MSSMNoFV_onshell*, unsigned, unsigned);

/** get stau mixing matrix */
double gm2calc_mssmnofv_get_UStau(const MSSMNoFV_onshell*, unsigned, unsigned);

/** get sup mixing matrix */
double gm2calc_mssmnofv_get_USu(const MSSMNoFV_onshell*, unsigned, unsigned);

/** get sdown mixing matrix */
double gm2calc_mssmnofv_get_USd(const MSSMNoFV_onshell*, unsigned, unsigned);

/** get scharm mixing matrix */
double gm2calc_mssmnofv_get_USc(const MSSMNoFV_onshell*, unsigned, unsigned);

/** get sstrange mixing matrix */
double gm2calc_mssmnofv_get_USs(const MSSMNoFV_onshell*, unsigned, unsigned);

/** get sstop mixing matrix */
double gm2calc_mssmnofv_get_USt(const MSSMNoFV_onshell*, unsigned, unsigned);

/** get sbottom mixing matrix */
double gm2calc_mssmnofv_get_USb(const MSSMNoFV_onshell*, unsigned, unsigned);

/** get lepton Yukawa coupling */
double gm2calc_mssmnofv_get_Ye(const MSSMNoFV_onshell*, unsigned, unsigned);

/** get down-quark Yukawa coupling */
double gm2calc_mssmnofv_get_Yd(const MSSMNoFV_onshell*, unsigned, unsigned);

/** get up-quark Yukawa coupling */
double gm2calc_mssmnofv_get_Yu(const MSSMNoFV_onshell*, unsigned, unsigned);


/** convert parameters to mixed on-shell/DR-bar scheme */
gm2calc_error gm2calc_mssmnofv_convert_to_onshell(MSSMNoFV_onshell*);

/** convert parameters to mixed on-shell/DR-bar scheme */
gm2calc_error gm2calc_mssmnofv_convert_to_onshell_params(
   MSSMNoFV_onshell*, double precision, unsigned max_iterations);

/** calculate mass spectrum */
gm2calc_error gm2calc_mssmnofv_calculate_masses(MSSMNoFV_onshell*);

/** check for problems */
int gm2calc_mssmnofv_have_problem(MSSMNoFV_onshell*);

/** check for warnings */
int gm2calc_mssmnofv_have_warning(MSSMNoFV_onshell*);

/** get problem descriptions */
void gm2calc_mssmnofv_get_problems(MSSMNoFV_onshell*, char*, unsigned);

/** get warning descriptions */
void gm2calc_mssmnofv_get_warnings(MSSMNoFV_onshell*, char*, unsigned);

/** print model */
void print_mssmnofv(const MSSMNoFV_onshell*);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
