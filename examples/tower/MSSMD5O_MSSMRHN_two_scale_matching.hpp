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

#ifndef MSSMD5O_MSSMRHN_TWO_SCALE_MATCHING_H
#define MSSMD5O_MSSMRHN_TWO_SCALE_MATCHING_H

#include "MSSMD5O_MSSMRHN_matching.hpp"
#include "MSSMD5O_input_parameters.hpp"
#include "single_scale_matching.hpp"

namespace flexiblesusy {

template<class T>
class MSSMD5O;

template<class T>
class MSSMRHN;

class Two_scale;

class MSSMD5O_MSSMRHN_matching {
public:
    MSSMD5O_MSSMRHN_matching();
    MSSMD5O_MSSMRHN_matching(const MSSMD5O_input_parameters&);
    void match_low_to_high_scale_model();
    void match_high_to_low_scale_model();
    double get_scale() const;
    void set_models(Model *lower, Model *upper);
    double get_initial_scale_guess() const;
    void set_lower_input_parameters(const MSSMD5O_input_parameters&);
    void set_scale(double) {}
    void reset();

private:
    double scale;
    double initial_scale_guess;
    double fixed_scale;
    MSSMD5O<Two_scale> *lower;
    MSSMRHN<Two_scale> *upper;
    MSSMD5O_input_parameters inputPars;

    void make_initial_scale_guess();
    void update_scale();
    void invert_seesaw_formula
    (const Eigen::Matrix3d& WOp, const Eigen::Vector3d& YvDiag,
     Eigen::Matrix3d& Yv, Eigen::Matrix3d& Mv);
};

template<>
class MSSMD5O_MSSMRHN_matching_up<Two_scale> : public Single_scale_matching {
public:
    MSSMD5O_MSSMRHN_matching_up();
    MSSMD5O_MSSMRHN_matching_up(const MSSMD5O_input_parameters&);
    void match();
    double get_scale() const;
    void set_models(Model *lower, Model *upper);
    double get_initial_scale_guess() const;
    void set_lower_input_parameters(const MSSMD5O_input_parameters&);
    void set_scale(double);
    void reset();

private:
    MSSMD5O_MSSMRHN_matching matching;
};

template<>
class MSSMD5O_MSSMRHN_matching_down<Two_scale> : public Single_scale_matching {
public:
    MSSMD5O_MSSMRHN_matching_down();
    MSSMD5O_MSSMRHN_matching_down(const MSSMD5O_input_parameters&);
    void match();
    double get_scale() const;
    void set_models(Model *upper, Model *lower);
    double get_initial_scale_guess() const;
    void set_lower_input_parameters(const MSSMD5O_input_parameters&);
    void set_scale(double);
    void reset();

private:
    MSSMD5O_MSSMRHN_matching matching;
};

}

#endif
