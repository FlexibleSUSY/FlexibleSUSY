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

"""Results of a calculation."""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from typing import Any, Dict, Iterator, List, Optional, Tuple

from . import arrays, slha
from .model import Metadata, ParameterInfo


@dataclass
class Channel:
    """A decay channel.  final_state holds FlexibleSUSY names such as "Fd(3)"
    (1-based, as in the SLHA output) when available."""

    pdg: Tuple[int, ...]
    width: float
    br: float
    final_state: Tuple[str, ...] = ()


@dataclass
class Decay:
    """Decays of one particle (mass eigenstate) with 0-based index."""

    particle: str
    index: int
    pdg: int
    width: float
    channels: List[Channel] = field(default_factory=list)


class Decays:
    """Decays by (particle name, 0-based index) or by PDG code."""

    def __init__(self, decays: List[Decay], problems: Optional[List[str]] = None,
                 warnings: Optional[List[str]] = None):
        self._list = decays
        self.problems = list(problems or [])
        self.warnings = list(warnings or [])
        self._by_name = {(d.particle, d.index): d for d in decays}
        self._by_pdg = {d.pdg: d for d in decays}

    def __getitem__(self, key: Tuple[str, int]) -> Decay:
        return self._by_name[key]

    def __contains__(self, key: Tuple[str, int]) -> bool:
        return key in self._by_name

    def by_pdg(self, pdg: int) -> Decay:
        return self._by_pdg[pdg]

    def __iter__(self) -> Iterator[Decay]:
        return iter(self._list)

    def __len__(self) -> int:
        return len(self._list)


@dataclass
class RunningParameters:
    """Running (DR-bar/MS-bar) model parameters at one scale."""

    scale: Optional[float]
    values: Dict[str, Any]

    def __getitem__(self, name: str) -> Any:
        return self.values[name]

    def __contains__(self, name: str) -> bool:
        return name in self.values


@dataclass
class SignalResult:
    """Result of HiggsSignals or Lilith: chi^2 of the model and of the SM
    (with Higgs mass mh_ref), and p-values for 1 and 2 degrees of freedom."""

    ndof: int
    chi2: float
    chi2_sm: float
    mh_ref: float
    pval_1d: float
    pval_2d: float


@dataclass
class HiggsBoundsResult:
    """HiggsBounds result for one particle: observed and expected ratio of
    the most sensitive analysis (> 1: excluded)."""

    pdg: int
    obs_ratio: float
    exp_ratio: float
    analysis: str


@dataclass
class HiggsToolsResult:
    higgssignals: Optional[SignalResult]
    higgsbounds: List[HiggsBoundsResult]


_COMPLEX_PART = re.compile(r"^(Re|Im)\((\w+)(?:\(([\d,]+)\))?\)$")


class Observables:
    """Observables by name, as in <M>_observables of the generated code.

    Complex observables, which the C++ code lists as Re(name(i)) and
    Im(name(i)), are combined into complex numbers or arrays.
    """

    def __init__(self, names: List[str], values: List[float], problems: List[str]):
        self.problems = list(problems)
        self.values: Dict[str, Any] = {}
        parts: Dict[str, Dict[Tuple[int, ...], complex]] = {}
        for name, value in zip(names, values):
            m = _COMPLEX_PART.match(name)
            if not m:
                self.values[name] = value
                continue
            base = m.group(2)
            idx = tuple(int(i) for i in m.group(3).split(",")) if m.group(3) else ()
            entry = parts.setdefault(base, {})
            z = entry.get(idx, 0j)
            entry[idx] = complex(value, z.imag) if m.group(1) == "Re" else complex(z.real, value)
        for base, entries in parts.items():
            if list(entries) == [()]:
                self.values[base] = entries[()]
            else:
                shape = tuple(max(i[k] for i in entries) + 1 for k in range(len(next(iter(entries)))))
                self.values[base] = arrays.from_entries(shape, entries, 0j)

    @property
    def ok(self) -> bool:
        return not self.problems

    def __getitem__(self, name: str) -> Any:
        return self.values[name]

    def __contains__(self, name: str) -> bool:
        return name in self.values

    def __iter__(self):
        return iter(self.values)

    def __len__(self) -> int:
        return len(self.values)


@dataclass
class Result:
    """Result of Model.calculate().

    Masses and mixing matrices are those of the SLHA output: masses of
    Majorana fermions can be negative (with real mixing matrices), and
    indices are 0-based.

    observables() and decays() are calculated on demand by the compiled
    backend, from the spectrum kept in memory.  The SLHA backend returns
    decays only if they were calculated in the same run (DecaySettings
    calculate_decays) and has no observables.
    """

    model: str
    problems: List[str]
    warnings: List[str]
    pole_masses: Dict[str, Any]
    mixings: Dict[str, Any]
    parameters: RunningParameters
    slha: str
    blocks: slha.SLHA
    _decays: Optional[Decays] = field(default=None, repr=False)
    _observables: Optional[Observables] = field(default=None, repr=False)
    _higgstools: Dict[Tuple[str, str], HiggsToolsResult] = field(default_factory=dict, repr=False)
    _lilith: Dict[str, Optional[SignalResult]] = field(default_factory=dict, repr=False)
    _backend: Any = field(default=None, repr=False, compare=False)
    _state: Any = field(default=None, repr=False, compare=False)

    @property
    def ok(self) -> bool:
        """True if FlexibleSUSY reported no problems."""
        return not self.problems

    def decays(self) -> Decays:
        """Decays, calculated on first call if the backend supports it."""
        if self._decays is None:
            self._decays = self._backend.calculate_decays(self)
        return self._decays

    def observables(self) -> Observables:
        """Observables, calculated on first call if the backend supports it."""
        if self._observables is None:
            self._observables = self._backend.calculate_observables(self)
        return self._observables

    def higgstools(self, higgsbounds_dataset: Optional[str] = None,
                   higgssignals_dataset: Optional[str] = None) -> HiggsToolsResult:
        """HiggsBounds and HiggsSignals results (needs FlexibleSUSY with HiggsTools).

        The datasets default to those given to Model().  The compiled backend
        calculates the decays first if needed."""
        hb, hs = self._backend.higgstools_datasets(higgsbounds_dataset, higgssignals_dataset)
        if (hb, hs) not in self._higgstools:
            self._higgstools[(hb, hs)] = self._backend.call_higgstools(self, hb, hs)
        return self._higgstools[(hb, hs)]

    def lilith(self, database: Optional[str] = None) -> Optional[SignalResult]:
        """Lilith result (needs FlexibleSUSY with Lilith); None if no Higgs
        boson has a mass between 123 and 128 GeV.  The database defaults to
        the one given to Model()."""
        db = self._backend.lilith_database(database)
        if db not in self._lilith:
            self._lilith[db] = self._backend.call_lilith(self, db)
        return self._lilith[db]

    def release(self) -> None:
        """Free the model kept in memory for on-demand calculations."""
        self._state = None

    def to_dict(self) -> Dict[str, Any]:
        """JSON-serializable dictionary (complex numbers as [re, im])."""
        return {
            "model": self.model,
            "ok": self.ok,
            "problems": self.problems,
            "warnings": self.warnings,
            "pole_masses": {k: arrays.to_list(v) for k, v in self.pole_masses.items()},
            "mixings": {k: arrays.to_list(v) for k, v in self.mixings.items()},
            "parameters": {"scale": self.parameters.scale,
                           **{k: arrays.to_list(v) for k, v in self.parameters.values.items()}},
            "decays": None if self._decays is None else {
                "problems": self._decays.problems,
                "warnings": self._decays.warnings,
                "particles": [{"particle": d.particle, "index": d.index, "pdg": d.pdg,
                               "width": d.width,
                               "channels": [{"pdg": list(c.pdg), "width": c.width, "br": c.br,
                                             "final_state": list(c.final_state)}
                                            for c in d.channels]}
                              for d in self._decays]},
            "observables": None if self._observables is None else {
                "problems": self._observables.problems,
                **{k: arrays.to_list(v) for k, v in self._observables.values.items()}},
        }


# ---------------------------------------------------------------------------
# building a result from SLHA output

def _info_block(data: slha.SLHA, name: str, key: int) -> List[str]:
    block = data.block(name)
    if block is None:
        return []
    return [str(v) for v in block.get_all(key)]


def _pole_masses(data: slha.SLHA, metadata: Metadata) -> Dict[str, Any]:
    mass = data.block("MASS")
    result = {}
    if mass is None:
        return result
    for p in metadata.particles:
        entries = {}
        for i, code in enumerate(p.pdg):
            value = mass.get(abs(code)) if code != 0 else None
            entries[(i,)] = float(value) if isinstance(value, float) else math.nan
        if any(not math.isnan(v) for v in entries.values()):
            result[p.name] = arrays.from_entries((p.multiplicity,), entries, math.nan)
    return result


def _read_parameter(data: slha.SLHA, info: ParameterInfo) -> Any:
    """Value of an output parameter, or None if it isn't in the output."""
    if info.slha is None:
        return None
    block = data.block(info.slha.block)
    im_block = data.block("IM" + info.slha.block)
    if block is None:
        return None

    def get(b, key):
        v = b.get(key) if b is not None else None
        return v if isinstance(v, float) else None

    if info.slha.index is not None:
        re_part = get(block, info.slha.index)
        if re_part is None:
            return None
        im_part = get(im_block, info.slha.index)
        value = re_part if im_part is None else complex(re_part, im_part)
        return value / info.slha_factor
    entries = {}
    for key, value in block:
        if not isinstance(value, float):
            continue
        idx = key if isinstance(key, tuple) else (key,)
        if len(idx) != len(info.shape):
            continue
        im_part = get(im_block, key)
        entries[tuple(i - 1 for i in idx)] = value if im_part is None else complex(value, im_part)
    if not entries:
        return None
    if not info.shape:
        return entries.get((0,))
    return arrays.from_entries(info.shape, entries)


_DECAY_COMMENT = re.compile(r"BR\((.+?)\s*->\s*(.+)\)\s*$")
_PARTICLE_COMMENT = re.compile(r"^(\w+)(?:\((\d+)\))?\s+decays$")


def decays_from_slha(data: slha.SLHA, metadata: Metadata) -> Decays:
    # PDG code -> (particle, 0-based index), from the metadata
    by_pdg = {}
    for p in metadata.particles:
        for i, code in enumerate(p.pdg):
            if code != 0:
                by_pdg.setdefault(code, (p.name, i))
                by_pdg.setdefault(-code, (p.name, i))
    result = []
    for code, d in data.decays.items():
        m = _PARTICLE_COMMENT.match(d.comment)
        if m:
            particle = m.group(1)
            index = int(m.group(2)) - 1 if m.group(2) else 0
        elif code in by_pdg:
            particle, index = by_pdg[code]
        else:
            particle, index = str(code), 0
        channels = []
        for c in d.channels:
            final_state: Tuple[str, ...] = ()
            cm = _DECAY_COMMENT.search(c.comment)
            if cm:
                final_state = tuple(cm.group(2).split())
            channels.append(Channel(pdg=c.pdg, width=c.br * d.width, br=c.br,
                                    final_state=final_state))
        result.append(Decay(particle=particle, index=index, pdg=code, width=d.width,
                            channels=channels))
    return Decays(result, problems=_info_block(data, "DCINFO", 4),
                  warnings=_info_block(data, "DCINFO", 3))


_MH_REF = re.compile(r"mh\s*=\s*([0-9.eE+-]+)")


def signal_result_from_slha(data: slha.SLHA, block_name: str) -> Optional[SignalResult]:
    """Block HIGGSSIGNALS or LILITH of FlexibleSUSY's SLHA output."""
    block = data.block(block_name)
    if block is None or 2 not in block:
        return None
    mh_ref = math.nan
    for (key, _), comment in zip(block.entries, block.comments):
        m = _MH_REF.search(comment)
        if key == 3 and m:
            mh_ref = float(m.group(1))
    return SignalResult(ndof=int(block.get(1, 0)), chi2=float(block[2]),
                        chi2_sm=float(block.get(3, math.nan)), mh_ref=mh_ref,
                        pval_1d=float(block.get(4, math.nan)),
                        pval_2d=float(block.get(5, math.nan)))


def higgsbounds_from_slha(data: slha.SLHA) -> List[HiggsBoundsResult]:
    """Block HIGGSBOUNDS: entries (pdg, 1) observed ratio with the analysis as
    comment, (pdg, 2) expected ratio."""
    block = data.block("HIGGSBOUNDS")
    if block is None:
        return []
    obs, exp, analysis = {}, {}, {}
    for (key, value), comment in zip(block.entries, block.comments):
        if isinstance(key, tuple) and len(key) == 2 and isinstance(value, float):
            pdg, kind = key
            if kind == 1:
                obs[pdg], analysis[pdg] = value, comment
            elif kind == 2:
                exp[pdg] = value
    return [HiggsBoundsResult(pdg=pdg, obs_ratio=obs[pdg], exp_ratio=exp.get(pdg, math.nan),
                              analysis=analysis[pdg]) for pdg in obs]


def from_slha(text: str, metadata: Metadata, backend=None, state=None,
              read_decays: bool = True) -> Result:
    data = slha.read(text)
    mixings = {}
    for info in metadata.mixing_matrices:
        value = _read_parameter(data, info)
        if value is not None:
            mixings[info.name] = value
    values = {}
    scale = None
    for info in metadata.model_parameters:
        value = _read_parameter(data, info)
        if value is not None:
            values[info.name] = value
            block = data.block(info.slha.block) if info.slha else None
            if scale is None and block is not None and block.scale is not None:
                scale = block.scale
    return Result(
        model=metadata.model,
        problems=_info_block(data, "SPINFO", 4),
        warnings=_info_block(data, "SPINFO", 3),
        pole_masses=_pole_masses(data, metadata),
        mixings=mixings,
        parameters=RunningParameters(scale=scale, values=values),
        slha=text,
        blocks=data,
        _decays=decays_from_slha(data, metadata) if read_decays and data.decays else None,
        _backend=backend,
        _state=state,
    )
