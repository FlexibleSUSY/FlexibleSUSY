# ====================================================================
# This file is part of FlexibleSUSY.
#
# FlexibleSUSY is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version.
#
# FlexibleSUSY is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with FlexibleSUSY.  If not, see
# <http://www.gnu.org/licenses/>.
# ====================================================================

"""Model-independent settings and Standard Model inputs.

Every field is optional.  A field left as None is not passed on, so
FlexibleSUSY uses its own default; this package doesn't duplicate the
default values.  Field names follow the C++ enums:

- Settings:       Spectrum_generator_settings  (SLHA block FlexibleSUSY)
- DecaySettings:  FlexibleDecay_settings       (SLHA block FlexibleDecay)
- PhysicalInputs: Physical_input               (SLHA block FlexibleSUSYInput)
- SMInputs:       softsusy::QedQcd             (SLHA block SMINPUTS)
"""

from __future__ import annotations

import enum
from dataclasses import dataclass, field, fields
from typing import Dict, Optional


class LoopLibrary(enum.IntEnum):
    """Loop library, FlexibleSUSY[31].  Decays need COLLIER or LOOPTOOLS."""

    SOFTSUSY = 0
    COLLIER = 1
    LOOPTOOLS = 2
    FFLITE = 3


def _entry(index: int):
    return field(default=None, metadata={"index": index})


class _SLHABlockSettings:
    """Common part: conversion of the set fields to SLHA block entries."""

    _block = ""

    def to_slha_entries(self) -> Dict[int, float]:
        entries = {}
        for f in fields(self):
            value = getattr(self, f.name)
            if value is not None:
                entries[f.metadata["index"]] = float(value)
        return entries


@dataclass
class Settings(_SLHABlockSettings):
    """Spectrum generator settings (SLHA block FlexibleSUSY), see doc/slha_input.rst."""

    _block = "FlexibleSUSY"

    precision: Optional[float] = _entry(0)
    max_iterations: Optional[int] = _entry(1)
    solver: Optional[int] = _entry(2)
    calculate_sm_masses: Optional[bool] = _entry(3)
    pole_mass_loop_order: Optional[int] = _entry(4)
    ewsb_loop_order: Optional[int] = _entry(5)
    beta_loop_order: Optional[int] = _entry(6)
    threshold_corrections_loop_order: Optional[int] = _entry(7)
    higgs_2loop_correction_at_as: Optional[bool] = _entry(8)
    higgs_2loop_correction_ab_as: Optional[bool] = _entry(9)
    higgs_2loop_correction_at_at: Optional[bool] = _entry(10)
    higgs_2loop_correction_atau_atau: Optional[bool] = _entry(11)
    force_output: Optional[bool] = _entry(12)
    top_pole_qcd_corrections: Optional[int] = _entry(13)
    beta_zero_threshold: Optional[float] = _entry(14)
    calculate_observables: Optional[bool] = _entry(15)
    force_positive_masses: Optional[bool] = _entry(16)
    pole_mass_scale: Optional[float] = _entry(17)
    eft_pole_mass_scale: Optional[float] = _entry(18)
    eft_matching_scale: Optional[float] = _entry(19)
    eft_matching_loop_order_up: Optional[int] = _entry(20)
    eft_matching_loop_order_down: Optional[int] = _entry(21)
    eft_higgs_index: Optional[int] = _entry(22)
    calculate_bsm_masses: Optional[bool] = _entry(23)
    threshold_corrections: Optional[int] = _entry(24)
    higgs_3loop_ren_scheme_atb_as2: Optional[int] = _entry(25)
    higgs_3loop_correction_at_as2: Optional[bool] = _entry(26)
    higgs_3loop_correction_ab_as2: Optional[bool] = _entry(27)
    higgs_3loop_correction_at2_as: Optional[bool] = _entry(28)
    higgs_3loop_correction_at3: Optional[bool] = _entry(29)
    higgs_4loop_correction_at_as3: Optional[bool] = _entry(30)
    loop_library: Optional[int] = _entry(31)
    calculate_amm: Optional[float] = _entry(32)


@dataclass
class DecaySettings(_SLHABlockSettings):
    """FlexibleDecay settings (SLHA block FlexibleDecay), see doc/slha_input.rst."""

    _block = "FlexibleDecay"

    calculate_decays: Optional[bool] = _entry(0)
    min_br_to_print: Optional[float] = _entry(1)
    include_higher_order_corrections: Optional[int] = _entry(2)
    use_Thomson_alpha_in_Phigamgam_and_PhigamZ: Optional[bool] = _entry(3)
    offshell_VV_decays: Optional[int] = _entry(4)
    print_effc_block: Optional[bool] = _entry(5)
    calculate_normalized_effc: Optional[bool] = _entry(6)
    call_higgstools: Optional[bool] = _entry(7)
    call_lilith: Optional[bool] = _entry(8)
    use_pole_higgs_mixings: Optional[bool] = _entry(9)


@dataclass
class PhysicalInputs(_SLHABlockSettings):
    """Extra physical inputs (SLHA block FlexibleSUSYInput)."""

    _block = "FlexibleSUSYInput"

    alpha_em_0: Optional[float] = _entry(0)
    mh_pole: Optional[float] = _entry(1)


@dataclass
class SMInputs(_SLHABlockSettings):
    """Standard Model inputs (SLHA block SMINPUTS).

    Note that entry 1 is 1/alpha_em(MZ), as in SLHA.
    """

    _block = "SMINPUTS"

    alpha_em_inv_mz: Optional[float] = _entry(1)
    g_fermi: Optional[float] = _entry(2)
    alpha_s_mz: Optional[float] = _entry(3)
    mz_pole: Optional[float] = _entry(4)
    mb_mb: Optional[float] = _entry(5)
    mt_pole: Optional[float] = _entry(6)
    mtau_pole: Optional[float] = _entry(7)
    mnu3_pole: Optional[float] = _entry(8)
    mw_pole: Optional[float] = _entry(9)
    me_pole: Optional[float] = _entry(11)
    mnu1_pole: Optional[float] = _entry(12)
    mmu_pole: Optional[float] = _entry(13)
    mnu2_pole: Optional[float] = _entry(14)
    md_2gev: Optional[float] = _entry(21)
    mu_2gev: Optional[float] = _entry(22)
    ms_2gev: Optional[float] = _entry(23)
    mc_mc: Optional[float] = _entry(24)
