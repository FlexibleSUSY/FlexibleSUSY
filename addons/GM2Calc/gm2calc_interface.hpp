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

#ifndef GM2CALC_INTERFACE_H
#define GM2CALC_INTERFACE_H

#include <Eigen/Core>

/**
 * @file gm2calc_interface.hpp
 * @brief contains declarations of GM2Calc interface functions
 */

namespace flexiblesusy {

/**
 * @class GM2Calc_data
 * @brief data to be passed to GM2Calc
 */
struct GM2Calc_data {
   GM2Calc_data();                ///< initializes members to GM2Calc default values
   void initialize();             ///< initializes members to GM2Calc default values

   double scale{0.0};             ///< renormalization scale
   double alpha_em_MZ{0.0};       ///< alpha_em(MZ)
   double alpha_em_0{0.0};        ///< alpha_em(0)
   double alpha_s_MZ{0.0};        ///< alpha_s(MZ) SM MS-bar
   double MZ{0.0};                ///< Z pole mass
   double MW{0.0};                ///< W pole mass
   double mb_mb{0.0};             ///< mb(mb) SM MS-bar
   double MT{0.0};                ///< top quark pole mass
   double MTau{0.0};              ///< tau lepton pole mass
   double MM{0.0};                ///< muon pole mass
   double MA0{0.0};               ///< CP-odd Higgs pole mass
   double MSvm{0.0};              ///< muon sneutrino pole mass
   double TB{0.0};                ///< tan(beta) DR-bar
   double Mu{0.0};                ///< mu parameter (initial guess)
   double M1{0.0};                ///< bino mass parameter (initial guess)
   double M2{0.0};                ///< wino mass parameter (initial guess)
   double M3{0.0};                ///< gluino mass parameter
   Eigen::Array<double,2,1> MSm{Eigen::Array<double,2,1>::Zero()};   ///< smuon pole masses
   Eigen::Array<double,2,1> MCha{Eigen::Array<double,2,1>::Zero()};  ///< chargino pole masses
   Eigen::Array<double,4,1> MChi{Eigen::Array<double,4,1>::Zero()};  ///< neutralino pole masses
   Eigen::Matrix<double,3,3> mq2{Eigen::Matrix<double,3,3>::Zero()}; ///< left-handed squark mass parameters squared
   Eigen::Matrix<double,3,3> mu2{Eigen::Matrix<double,3,3>::Zero()}; ///< right-handed up-type squark mass parameters squared
   Eigen::Matrix<double,3,3> md2{Eigen::Matrix<double,3,3>::Zero()}; ///< right-handed down-type squark mass parameters squared
   Eigen::Matrix<double,3,3> ml2{Eigen::Matrix<double,3,3>::Zero()}; ///< left-handed slepton mass parameters squared
   Eigen::Matrix<double,3,3> me2{Eigen::Matrix<double,3,3>::Zero()}; ///< right-handed down-type slepton mass parameters squared
   Eigen::Matrix<double,3,3> Au{Eigen::Matrix<double,3,3>::Zero()};  ///< up-type squark trilinear coupling
   Eigen::Matrix<double,3,3> Ad{Eigen::Matrix<double,3,3>::Zero()};  ///< down-type squark trilinear coupling
   Eigen::Matrix<double,3,3> Ae{Eigen::Matrix<double,3,3>::Zero()};  ///< down-type slepton trilinear coupling
};

/// calculates amu using GM2Calc
double gm2calc_calculate_amu(const GM2Calc_data&);

/// calculates uncertainty of amu using GM2Calc
double gm2calc_calculate_amu_uncertainty(const GM2Calc_data&);

} // namespace flexiblesusy

#endif
