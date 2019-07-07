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

#ifndef SM_TWO_LOOP_MT_H
#define SM_TWO_LOOP_MT_H

#ifdef TSIL_SIZE_DOUBLE
typedef double TSIL_REAL;
#else
typedef long double TSIL_REAL;
#endif

/**
 * @file sm_twoloop_mt.hpp
 * @brief declaration of 2-loop corrections to mt(MS-bar,SM) from [arXiv:1604.01134]
 * @author Harun Acaroglu
 * @author Alexander Voigt
 *
 * Conventions:
 *
 * g3: SM MS-bar strong gauge coupling
 * yt: SM MS-bar top Yukawa coupling
 * s : squared momentum
 * t : squared SM MS-bar top mass
 * h : squared SM MS-bar higgs mass
 * qq: squared MS-bar renormalization scale
 */

namespace flexiblesusy {
namespace sm_twoloop_mt {

/* ******************** 1-loop ******************** */

/// returns 1-loop O(as) correction to Mt
TSIL_REAL delta_Mt_1loop_as(
   TSIL_REAL g3, TSIL_REAL t, TSIL_REAL qq);

/// returns 1-loop O(at) correction to Mt
TSIL_REAL delta_Mt_1loop_at(
   TSIL_REAL t, TSIL_REAL h, TSIL_REAL yt, TSIL_REAL qq);

/// returns 1-loop O(as) correction to mt
TSIL_REAL delta_mt_1loop_as(
   TSIL_REAL g3, TSIL_REAL t, TSIL_REAL qq);

/// returns scalar part of 1-loop O(at) correction to mt
TSIL_REAL delta_mt_1loop_at_S(
   TSIL_REAL yt, TSIL_REAL t, TSIL_REAL h, TSIL_REAL s, TSIL_REAL qq);

/// returns left part of 1-loop O(at) correction to mt
TSIL_REAL delta_mt_1loop_at_L(
   TSIL_REAL yt, TSIL_REAL t, TSIL_REAL h, TSIL_REAL s, TSIL_REAL qq);

/// returns right part of 1-loop O(at) correction to mt
TSIL_REAL delta_mt_1loop_at_R(
   TSIL_REAL yt, TSIL_REAL t, TSIL_REAL h, TSIL_REAL s, TSIL_REAL qq);

/* ******************** 2-loop ******************** */

/// returns 2-loop O(as^2) correction to Mt, Eq.(2.1) of [1604.01134]
TSIL_REAL delta_Mt_2loop_as_as(
   TSIL_REAL g3, TSIL_REAL t, TSIL_REAL qq);

/// returns 2-loop O(as^2) correction to mt in FlexibleSUSY convention, Eq.(9) of [1710.03760]
TSIL_REAL delta_mt_2loop_as_as_flexiblesusy(
   TSIL_REAL g3, TSIL_REAL t, TSIL_REAL qq);

/// returns scalar part of 2-loop O(as*at) correction to mt in FlexibleSUSY convention
TSIL_REAL delta_mt_2loop_as_at_S_flexiblesusy(
   TSIL_REAL g3, TSIL_REAL yt, TSIL_REAL t, TSIL_REAL h, TSIL_REAL s, TSIL_REAL qq);

/// returns left+right part of 2-loop O(as*at) correction to mt in FlexibleSUSY convention
TSIL_REAL delta_mt_2loop_as_at_LR_flexiblesusy(
   TSIL_REAL g3, TSIL_REAL yt, TSIL_REAL t, TSIL_REAL h, TSIL_REAL s, TSIL_REAL qq);

/// returns scalar part of 2-loop O(at^2) correction to mt in FlexibleSUSY convention
TSIL_REAL delta_mt_2loop_at_at_S_flexiblesusy(
   TSIL_REAL yt, TSIL_REAL t, TSIL_REAL h, TSIL_REAL s, TSIL_REAL qq);

/// returns left+right part of 2-loop O(at^2) correction to mt in FlexibleSUSY convention
TSIL_REAL delta_mt_2loop_at_at_LR_flexiblesusy(
   TSIL_REAL yt, TSIL_REAL t, TSIL_REAL h, TSIL_REAL s, TSIL_REAL qq);

/// returns 2-loop O(as^2) correction to mt in SPheno convention
TSIL_REAL delta_mt_2loop_as_as_spheno(
   TSIL_REAL g3, TSIL_REAL t, TSIL_REAL qq);

/// returns scalar part of 2-loop O(as*at) correction to mt in SPheno convention
TSIL_REAL delta_mt_2loop_as_at_S_spheno(
   TSIL_REAL g3, TSIL_REAL yt, TSIL_REAL t, TSIL_REAL h, TSIL_REAL s, TSIL_REAL qq);

/// returns left+right part of 2-loop O(as*at) correction to mt in SPheno convention
TSIL_REAL delta_mt_2loop_as_at_LR_spheno(
   TSIL_REAL g3, TSIL_REAL yt, TSIL_REAL t, TSIL_REAL h, TSIL_REAL s, TSIL_REAL qq);

/// returns scalar part of 2-loop O(at^2) correction to mt in SPheno convention
TSIL_REAL delta_mt_2loop_at_at_S_spheno(
   TSIL_REAL yt, TSIL_REAL t, TSIL_REAL h, TSIL_REAL s, TSIL_REAL qq);

/// returns left+right part of 2-loop O(at^2) correction to mt in SPheno convention
TSIL_REAL delta_mt_2loop_at_at_LR_spheno(
   TSIL_REAL yt, TSIL_REAL t, TSIL_REAL h, TSIL_REAL s, TSIL_REAL qq);

} // namespace sm_twoloop_mt
} // namespace flexiblesusy

#endif
