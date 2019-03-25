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

#ifndef TEST_munuSSM_H
#define TEST_munuSSM_H

#include <tuple>

#include "lowe.h"
#include "munuSSM_two_scale_spectrum_generator.hpp"

using namespace flexiblesusy;

munuSSM_slha<munuSSM<Two_scale>>
setup_munuSSM(const munuSSM_input_parameters& input, softsusy::QedQcd& qedqcd)
{

    Spectrum_generator_settings settings;
    settings.set(Spectrum_generator_settings::precision, 1e-4);
    settings.set(Spectrum_generator_settings::max_iterations, 0);
    settings.set(Spectrum_generator_settings::calculate_sm_masses, 1);
    settings.set(Spectrum_generator_settings::calculate_bsm_masses, 1);
    settings.set(Spectrum_generator_settings::ewsb_loop_order, 1);
    settings.set(Spectrum_generator_settings::pole_mass_loop_order, 1);
    settings.set(Spectrum_generator_settings::beta_loop_order, 2);
    settings.set(Spectrum_generator_settings::threshold_corrections_loop_order, 2);
    settings.set(Spectrum_generator_settings::top_pole_qcd_corrections, 1);
    settings.set(Spectrum_generator_settings::threshold_corrections, 123111321);
    settings.set(Spectrum_generator_settings::higgs_2loop_correction_at_as, 1);
    settings.set(Spectrum_generator_settings::higgs_2loop_correction_ab_as, 1);
    settings.set(Spectrum_generator_settings::higgs_2loop_correction_at_at, 1);
    settings.set(Spectrum_generator_settings::higgs_2loop_correction_atau_atau, 1);
    settings.set(Spectrum_generator_settings::higgs_3loop_ren_scheme_atb_as2, 1);
    settings.set(Spectrum_generator_settings::higgs_3loop_correction_at_as2, 1);
    settings.set(Spectrum_generator_settings::higgs_3loop_correction_ab_as2, 1);
    settings.set(Spectrum_generator_settings::higgs_3loop_correction_at2_as, 1);
    settings.set(Spectrum_generator_settings::higgs_3loop_correction_at3, 1);
    settings.set(Spectrum_generator_settings::higgs_4loop_correction_at_as3, 1);
    settings.set(Spectrum_generator_settings::beta_zero_threshold, 1e-16);

    qedqcd.setAlpha(softsusy::ALPHA, 1./1.279340000e+02);
    qedqcd.setAlphaEmInput(1./1.279340000e+02);
    qedqcd.setFermiConstant(1.166378700e-05);
    qedqcd.setAlphaSInput(1.176000000e-01);
    qedqcd.setAlpha(softsusy::ALPHAS, 1.176000000e-01);
    qedqcd.setPoleMZ(9.118760000e+01);
    qedqcd.setMbMb(4.200000000e+00);
    qedqcd.setPoleMt(1.733000000e+02);
    qedqcd.setPoleMtau(1.777000000e+00);
    qedqcd.setNeutrinoPoleMass(3, 0);
    qedqcd.setPoleMW(80.404);
    qedqcd.setPoleMel(5.109989020e-04);
    qedqcd.setNeutrinoPoleMass(1, 0);
    qedqcd.setPoleMmuon(1.056583570e-01);
    qedqcd.setMass(softsusy::mMuon, 1.056583570e-01);
    qedqcd.setNeutrinoPoleMass(2, 0);
    qedqcd.setMd2GeV(4.750000000e-03);
    qedqcd.setMu2GeV(2.400000000e-03);
    qedqcd.setMs2GeV(1.040000000e-01);
    qedqcd.setMcMc(1.270000000e+00);

    munuSSM_spectrum_generator<Two_scale> spectrum_generator;
    spectrum_generator.set_settings(settings);
    spectrum_generator.run(qedqcd, input);

    return std::get<0>(spectrum_generator.get_models_slha());
}

#endif
