/* ====================================================================
 * This file is part of FlexibleSUSY.
 *
 * FlexibleSUSY is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * FlexibleSUSY is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with FlexibleSUSY.  If not, see
 * <http://www.gnu.org/licenses/>.
 *  ====================================================================
 */

#ifndef _MSSM_TWOLOOPHIGGS_H_
#define _MSSM_TWOLOOPHIGGS_H_

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @file mssm_twoloophiggs.h
 * @brief C interface for two-loop routines from Slavich et al.
 */

/**
 * Two-loop O(a_t a_s) corrections to the CP-even Higgs mass matrix.
 * Routine written by P. Slavich (e-mail: slavich@pd.infn.it).
 * Based on G. Degrassi, P. Slavich and F. Zwirner,
 * Nucl. Phys. B611 (2001) 403 [hep-ph/0105096].
 *
 * Last update:  13/12/2001.
 *
 * I/O PARAMETERS:
 * t = m_top^2, g = m_gluino^2, T1 = m_stop1^2, T2 = m_stop2^2,
 * st = sin(theta_stop), ct = cos(theta_stop), q = Q^2 (ren. scale),
 * mu = Higgs mixing parameter, tanb = tan(beta), vv = v^2,
 * OS = renormalization scheme for 1-loop (0 = DRbar, 1 = On-Shell),
 * Sij = 2-loop corrections to the CP-even Higgs mass matrix elements.
 */
int dszhiggs_(double * t, double * g, double * T1, double * T2,
              double * st, double * ct, double * q, double * mu,
              double * tanb, double *vv, double * gs, int * OS,
              double * S11, double * S22, double * S12);

/**
 * Two-loop O(a_t^2 + at ab + ab^2) corrections to the CP-even Higgs masses
 * Written by P. Slavich (e-mail: slavich@mppmu.mpg.de).
 * Based on A. Dedes, G. Degrassi and P. Slavich, hep-ph/0305127.
 *
 * Last update:  13/05/2003.
 *
 * I/O PARAMETERS:
 * t = m_top^2, b = m_bot^2, A0 = m_A^2, T1 = m_stop1^2, T2 = m_stop2^2,
 * B1 = m_sbot1^2, B2 = m_sbot2^2, st = sin(theta_stop),
 * ct = cos(theta_stop), sb = sin(theta_sbot), cb = cos(theta_sbot),
 * q = Q^2 (ren. scale), mu = Higgs mixing parameter, tanb = tan(beta),
 * vv = v^2,
 * Sij = 2-loop corrections to the CP-even Higgs mass matrix elements,
 *
 * Notice: we assume that the 1-loop part is computed in terms of
 *         running (DRbar) parameters, evaluated at the scale Q. The
 *         parameters in the bottom/sbottom sector should be computed
 *         in terms of the "resummed" bottom Yukawa coupling.
 */
int ddshiggs_(double * t, double * b, double * A0, double * T1,
              double * T2, double * B1, double * B2, double * st,
              double * ct, double * sb, double * cb, double * q,
              double * mu, double * tanb, double * vv,
              double * S11, double * S12, double * S22);

/**
 * Two-loop O(a_tau^2) corrections to the CP-even Higgs mass matrix.
 * Routine written by P. Slavich (e-mail: slavich@mppmu.mpg.de).
 * Based on A. Brignole, G. Degrassi, P. Slavich and F. Zwirner,
 * hep-ph/0112177 with appropriate substitutions for top -> tau.
 *
 * Last update:  19/06/2003.
 *
 * I/O PARAMETERS:
 * t = m_tau^2, A0 = m_A^2, BL = m_snu^2, T1 = m_stau1^2, T2 = m_stau2^2,
 * st = sin(theta_stau), ct = cos(theta_stau), q = Q^2 (ren. scale),
 * mu = Higgs mixing parameter, tb = tan(beta), vv = v^2,
 * OS = renormalization scheme for 1-loop (0 = DRbar, 1 = On-Shell),
 * Sij = 2-loop corrections to the CP-even Higgs mass matrix elements.
 */
int tausqhiggs_(double * t, double * A0, double * BL, double * T1,
                double * T2, double * st, double * ct, double * q,
                double * mu, double * tanb, double * vv, int * OS,
                double * S11, double * S22, double * S12);

/**
 * Two-loop O(a_tau * a_bottom) corrections to the CP-even Higgs mass
 * matrix.  Routine written by P. Slavich (e-mail:
 * slavich@mppmu.mpg.de).  Based on Allanach et.al., JHEP 0409 (2004)
 * 044 [hep-ph/0406166].
 *
 * I/O PARAMETERS:
 * t = m_tau^2, b = m_b^2,
 * T1 = m_stau1^2, T2 = m_stau2^2, B1 = m_sbottom1^2, B2 = m_sbottom2^2,
 * st = sin(theta_stau), ct = cos(theta_stau),
 * sb = sin(theta_sbottom), cb = cos(theta_sbottom),
 * q = Q^2 (ren. scale), mu = Higgs mixing parameter, tb = tan(beta),
 * vv = v^2, Sij = 2-loop corrections to the CP-even Higgs mass matrix
 * elements.
 */
int taubot_(double * t, double * b,
            double * T1, double * T2, double * B1, double * B2,
            double * st, double * ct, double * sb, double * cb,
            double * q, double * mu, double * tanb, double * vv,
            double * S11, double * S22, double * S12);

/**
 * Two-loop O(a_t a_s) corrections to the CP-odd Higgs mass in the
 * DRbar scheme.  Written by P. Slavich (e-mail: slavich@pd.infn.it).
 * Based on G. Degrassi, P. Slavich and F. Zwirner,
 * Nucl. Phys. B611 (2001) 403 [hep-ph/0105096].
 *
 * Last update:  13/12/2001.
 *
 * I/O PARAMETERS:
 * t = m_top^2, g = m_gluino^2, T1 = m_stop1^2, T2 = m_stop2^2,
 * st = sin(theta_stop), ct = cos(theta_stop), q = Q^2 (ren. scale),
 * mu = Higgs mixing parameter, tanb = tan(beta), vv = v^2,
 * gs = strong coupling,
 * dma = 2-loop corrections to the CP-odd Higgs mass.
 */
int dszodd_(double * t, double * g, double * T1,
            double * T2, double * st, double * ct,
            double * q, double * mu, double * tanb,
            double * vv, double * gs, double * dma);

/**
 * Two-loop O(a_t^2 + at ab + ab^2) corrections to the CP-odd Higgs mass
 * Written by P. Slavich (e-mail: slavich@mppmu.mpg.de).
 * Based on A. Dedes, G. Degrassi and P. Slavich, hep-ph/0305127.
 *
 * Last update:  13/05/2003.
 *
 * I/O PARAMETERS:
 * t = m_top^2, b = m_bot^2, A0 = m_A^2, T1 = m_stop1^2, T2 = m_stop2^2,
 * B1 = m_sbot1^2, B2 = m_sbot2^2, st = sin(theta_stop),
 * ct = cos(theta_stop), sb = sin(theta_sbot), cb = cos(theta_sbot),
 * q = Q^2 (ren. scale), mu = Higgs mixing parameter, tanb = tan(beta),
 * vv = v^2,
 * dma = 2-loop corrections to the CP-odd Higgs mass.
 *
 * Notice: we assume that the 1-loop part is computed in terms of
 *         running (DRbar) parameters, evaluated at the scale Q. The
 *         parameters in the bottom/sbottom sector should be computed
 *         in term of the "resummed" bottom Yukawa coupling.
*/
int ddsodd_(double * t, double * b, double * A0, double * T1,
            double * T2, double * B1, double * B2, double * st,
            double * ct, double * sb, double * cb, double * q,
            double * mu, double * tanb, double * vv, double * dma);

/**
 * Two-loop O(a_tau^2) corrections to the CP-odd Higgs mass in the
 * DRbar scheme. Written by P. Slavich (e-mail: slavich@mppmu.mpg.de).
 * Based on A. Brignole, G. Degrassi, P. Slavich and F. Zwirner,
 * hep-ph/0112177 with appropriate substitutions for top -> tau.
 *
 * Last update:  19/06/2003.
 *
 * I/O PARAMETERS:
 * t = m_tau^2, A0 = m_A^2, BL = m_snu^2, T1 = m_stau1^2, T2 = m_stau2^2,
 * st = sin(theta_stau), ct = cos(theta_stau), q = Q^2 (ren. scale),
 * mu = Higgs mixing parameter, tb = tan(beta), vv = v^2,
 * dma = 2-loop corrections to the CP-odd Higgs mass.
 */
int tausqodd_(double * t, double * A0, double * BL, double * T1,
              double * T2, double * st, double * ct, double * q,
              double * mu, double * tanb, double * vv,
              double * dma);

/**
 * Two-loop O(a_tau * a_bottom) corrections to the CP-odd Higgs mass
 * matrix.  Routine written by P. Slavich (e-mail:
 * slavich@mppmu.mpg.de).  Based on Allanach et.al., JHEP 0409 (2004)
 * 044 [hep-ph/0406166].
 *
 * I/O PARAMETERS:
 * t = m_tau^2, b = m_b^2,
 * T1 = m_stau1^2, T2 = m_stau2^2, B1 = m_sbottom1^2, B2 = m_sbottom2^2,
 * st = sin(theta_stau), ct = cos(theta_stau),
 * sb = sin(theta_sbottom), cb = cos(theta_sbottom),
 * q = Q^2 (ren. scale), mu = Higgs mixing parameter, tb = tan(beta),
 * vv = v^2, dma = 2-loop corrections to the CP-odd Higgs mass.
 */
int taubotodd_(double * t, double * b,
            double * T1, double * T2, double * B1, double * B2,
            double * st, double * ct, double * sb, double * cb,
            double * q, double * mu, double * tanb, double * vv,
            double * dma);

/**
 * Two-loop O(a_t a_s) corrections to the Higgs tadpoles.
 * Written by P. Slavich (e-mail: slavich@mppmu.mpg.de).
 * Based on A. Dedes and P. Slavich, hep-ph/0212132.
 *
 * Last update:  12/12/2002.
 *
 * I/O PARAMETERS:
 * t = m_top^2, g = m_gluino^2, T1 = m_stop1^2, T2 = m_stop2^2,
 * st = sin(theta_stop), ct = cos(theta_stop), q = Q^2 (ren. scale),
 * mu = Higgs mixing parameter, tanb = tan(beta), vv = v^2,
 * gs = strong coupling constant
 * Si = 1/vi*dVeff/dvi = 2-loop corrections to the Higgs tadpoles.
 *
 * Notice: we assume that the 1-loop part is computed in terms of
 *         running (DRbar) parameters, evaluated at the scale Q.
 */
int ewsb2loop_(double * t, double * g, double *  T1, double * T2,
               double * st, double * ct, double * q, double * mu,
               double * tanb, double * vv, double * gs,
               double * s1, double * s2);

/**
 * Two-loop O(a_t^2 + at ab + ab^2) corrections to the
 * minimization conditions of the effective potential.
 * Written by P. Slavich (e-mail: slavich@mppmu.mpg.de).
 * Based on A. Dedes, G. Degrassi and P. Slavich, hep-ph/0305127.
 *
 * Last update:  13/05/2003.
 *
 * I/O PARAMETERS:
 * t = m_top^2, b = m_bot^2, A0 = m_A^2, T1 = m_stop1^2, T2 = m_stop2^2,
 * B1 = m_sbot1^2, B2 = m_sbot2^2, st = sin(theta_stop),
 * ct = cos(theta_stop), sb = sin(theta_sbot), cb = cos(theta_sbot),
 * q = Q^2 (ren. scale), mu = Higgs mixing parameter, tanb = tan(beta),
 * vv = v^2,
 * Si = 1/vi*dVeff/dvi = 2-loop corrections to the Higgs tadpoles,
 *
 * Notice: we assume that the 1-loop part is computed in terms of
 *         running (DRbar) parameters, evaluated at the scale Q. The
 *         parameters in the bottom/sbottom sector should be computed
 *         in term of the "resummed" bottom Yukawa coupling.
 */
int ddstad_(double * t, double * b, double * A0, double * T1,
            double * T2, double * B1, double * B2, double * st,
            double * ct, double * sb, double * cb, double * q,
            double * mu, double * tanb, double * vv,
            double * s1, double * s2);

/**
 * Two-loop O(a_t^2) corrections to the Higgs tadpoles.
 * Written by P. Slavich (e-mail: slavich@mppmu.mpg.de).
 * Based on A. Dedes and P. Slavich, hep-ph/0212132
 * with appropriate substitutions for top -> tau.
 *
 * Last update:  19/06/2003.
 *
 * I/O PARAMETERS:
 * t = m_tau^2, A0 = m_A^2, BL = m_snu^2, T1 = m_stau1^2, T2 = m_stau2^2,
 * st = sin(theta_stau), ct = cos(theta_stau), q = Q^2 (ren. scale),
 * mu = Higgs mixing parameter, tb = tan(beta), vv = v^2,
 * Si = 1/vi*dVeff/dvi = 2-loop corrections to the Higgs tadpoles.
 *
 * Notice: we assume that the 1-loop part is computed in terms of
 * running (DRbar) parameters, evaluated at the scale Q.
 */
int tausqtad_(double * t, double * A0, double * BL, double * T1,
			 double * T2, double * st, double * ct, double * q,
			 double * mu, double * tanb, double * vv,
			 double * s1, double * s2);

/**
 * Two-loop O(a_tau * a_bottom) corrections to the Higgs tadpoles.
 * Routine written by P. Slavich (e-mail: slavich@mppmu.mpg.de).
 * Based on Allanach et.al., JHEP 0409 (2004) 044 [hep-ph/0406166].
 *
 * I/O PARAMETERS:
 * t = m_tau^2, b = m_b^2,
 * T1 = m_stau1^2, T2 = m_stau2^2, B1 = m_sbottom1^2, B2 = m_sbottom2^2,
 * st = sin(theta_stau), ct = cos(theta_stau),
 * sb = sin(theta_sbottom), cb = cos(theta_sbottom),
 * q = Q^2 (ren. scale), mu = Higgs mixing parameter, tanbb = tan(beta),
 * vv = v^2,
 * si = 1/vi*dVeff/dvi = 2-loop corrections to the Higgs tadpoles.
 */
int taubottad_(double * t, double * b,
               double * T1, double * T2, double * B1, double * B2,
               double * st, double * ct, double * sb, double * cb,
               double * q, double * mu, double * tanb, double * vv,
               double * s1, double * s2);

#ifdef __cplusplus
}
#endif

#endif
