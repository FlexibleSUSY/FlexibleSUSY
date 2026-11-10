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

"""Models, their metadata, and parameter points."""

from __future__ import annotations

import dataclasses
import difflib
import json
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

from . import arrays, slha
from .errors import BackendError, InputError
from .settings import DecaySettings, PhysicalInputs, Settings, SMInputs

SUPPORTED_METADATA_FORMAT = 1


# ---------------------------------------------------------------------------
# locating the FlexibleSUSY directory

def _looks_like_root(path: Path) -> bool:
    return (path / "models").is_dir()


def find_root(root: Union[str, os.PathLike, None] = None) -> Path:
    """Directory of FlexibleSUSY (or of an exported model package).

    Order: the given root, the environment variable FLEXIBLESUSY_DIR, then
    the checkout this package is in (<root>/python/flexiblesusy).
    """
    candidates: List[Tuple[str, Path]] = []
    if root is not None:
        candidates.append(("root argument", Path(root)))
    elif os.environ.get("FLEXIBLESUSY_DIR"):
        candidates.append(("FLEXIBLESUSY_DIR", Path(os.environ["FLEXIBLESUSY_DIR"])))
    else:
        candidates.append(("package location", Path(__file__).resolve().parents[2]))
    source, path = candidates[0]
    path = path.expanduser().resolve()
    if not _looks_like_root(path):
        raise BackendError(
            f"{path} (from {source}) is not a FlexibleSUSY directory: it has no "
            "models/ directory.  Pass root=... or set FLEXIBLESUSY_DIR.")
    return path


# ---------------------------------------------------------------------------
# metadata

@dataclass(frozen=True)
class SLHALocation:
    block: str
    index: Optional[int]    # None: the whole block (vectors, matrices)


@dataclass(frozen=True)
class ParameterInfo:
    name: str
    kind: str               # "real" | "complex" | "integer"
    shape: Tuple[int, ...]
    slha: Optional[SLHALocation]
    slha_factor: float = 1.0


@dataclass(frozen=True)
class ParticleInfo:
    name: str
    multiplicity: int
    pdg: Tuple[int, ...]
    sm: bool


def _location(data: Optional[Dict[str, Any]]) -> Optional[SLHALocation]:
    if data is None:
        return None
    return SLHALocation(block=str(data["block"]), index=data.get("index"))


def _parameter(data: Dict[str, Any]) -> ParameterInfo:
    return ParameterInfo(name=data["name"], kind=data["kind"], shape=tuple(data["shape"]),
                         slha=_location(data.get("slha")),
                         slha_factor=float(data.get("slha_factor", 1.0)))


@dataclass(frozen=True)
class Metadata:
    """Content of models/<M>/<M>_metadata.json, written by the meta code."""

    model: str
    flexiblesusy_version: str
    input_parameters: Tuple[ParameterInfo, ...]
    particles: Tuple[ParticleInfo, ...]
    mixing_matrices: Tuple[ParameterInfo, ...]
    model_parameters: Tuple[ParameterInfo, ...]

    @classmethod
    def load(cls, path: Path) -> "Metadata":
        try:
            data = json.loads(path.read_text())
        except FileNotFoundError:
            raise BackendError(
                f"{path} is missing.  It is written by the code generation; "
                "regenerate the model with a FlexibleSUSY version that has the "
                "Python interface.") from None
        version = data.get("format_version")
        if version != SUPPORTED_METADATA_FORMAT:
            raise BackendError(f"{path}: unsupported metadata format {version!r} "
                               f"(this package reads format {SUPPORTED_METADATA_FORMAT})")
        return cls(
            model=data["model"],
            flexiblesusy_version=data["flexiblesusy_version"],
            input_parameters=tuple(_parameter(p) for p in data["input_parameters"]),
            particles=tuple(ParticleInfo(name=p["name"], multiplicity=p["multiplicity"],
                                         pdg=tuple(p["pdg"]), sm=p["sm"])
                            for p in data["particles"]),
            mixing_matrices=tuple(_parameter(p) for p in data["mixing_matrices"]),
            model_parameters=tuple(_parameter(p) for p in data["model_parameters"]),
        )


# ---------------------------------------------------------------------------
# input parameters

def _is_imaginary_part(info: ParameterInfo) -> bool:
    return info.slha is not None and info.slha.block.upper() in ("IMMINPAR", "IMEXTPAR")


class _InputMap:
    """Maps input parameter names to the SLHA entries they are read from."""

    def __init__(self, infos: Sequence[ParameterInfo]):
        self.real: Dict[str, ParameterInfo] = {}
        self.imag: Dict[str, ParameterInfo] = {}
        for info in infos:
            target = self.imag if _is_imaginary_part(info) else self.real
            target[info.name] = info

    def names(self) -> List[str]:
        return list(self.real)

    def settable(self) -> List[str]:
        """Inputs that can be set through SLHA."""
        return [n for n, i in self.real.items() if i.slha is not None]

    def check_name(self, name: str) -> ParameterInfo:
        if name not in self.real:
            close = difflib.get_close_matches(name, self.names(), n=3)
            hint = f"; did you mean {', '.join(close)}?" if close else ""
            raise InputError(f"unknown input parameter {name!r}{hint}")
        return self.real[name]


def _check_value(info: ParameterInfo, value: Any) -> Dict[Tuple[int, ...], Any]:
    try:
        entries = arrays.to_entries(value, info.shape, info.name)
    except ValueError as e:
        raise InputError(str(e)) from None
    for idx, v in entries.items():
        if isinstance(v, complex) and info.kind != "complex":
            if v.imag != 0:
                raise InputError(f"{info.name}: complex value for a {info.kind} parameter")
            entries[idx] = v = v.real
        if isinstance(v, bool):
            entries[idx] = v = int(v)
        if info.kind == "integer" and float(v) != int(v):
            raise InputError(f"{info.name}: expected an integer, got {v!r}")
    return entries


def _add_entries(blocks: Dict[str, Dict[slha.Key, slha.Value]], location: SLHALocation,
                 entries: Dict[Tuple[int, ...], Any], part) -> None:
    block = blocks.setdefault(location.block, {})
    for idx, v in entries.items():
        value = float(part(v))
        if location.index is not None:
            block[location.index] = value
        else:
            key = tuple(i + 1 for i in idx)  # SLHA is 1-based
            block[key[0] if len(key) == 1 else key] = value


# ---------------------------------------------------------------------------
# points

@dataclass(frozen=True)
class Point:
    """A validated parameter point of a model.  Create it with Model.point()."""

    model: str
    inputs: Dict[str, Any]
    settings: Settings = field(default_factory=Settings)
    sm_inputs: SMInputs = field(default_factory=SMInputs)
    physical_inputs: PhysicalInputs = field(default_factory=PhysicalInputs)
    decay_settings: DecaySettings = field(default_factory=DecaySettings)
    extra_blocks: Dict[str, Dict[slha.Key, slha.Value]] = field(default_factory=dict)


class Model:
    """A model generated and built with FlexibleSUSY.

    Model("MSSM") finds models/MSSM in the FlexibleSUSY directory (see
    find_root) and reads its metadata.  backend="slha" runs models/MSSM/run_MSSM.x.
    """

    def __init__(self, name: str, root: Union[str, os.PathLike, None] = None,
                 backend: str = "auto"):
        self.name = name
        self.root = find_root(root)
        self.directory = self.root / "models" / name
        if not self.directory.is_dir():
            available = sorted(p.name for p in (self.root / "models").iterdir() if p.is_dir())
            raise BackendError(f"model {name!r} not found in {self.root / 'models'}; "
                               f"available: {', '.join(available) or 'none'}")
        self.metadata = Metadata.load(self.directory / f"{name}_metadata.json")
        self._inputs = _InputMap(self.metadata.input_parameters)
        self.backend = self._make_backend(backend)

    def _make_backend(self, backend: str):
        from .backends.slha_backend import SLHABackend
        if backend in ("auto", "slha"):
            return SLHABackend(self)
        raise InputError(f"unknown backend {backend!r}; available: 'auto', 'slha'")

    def __repr__(self) -> str:
        return f"Model({self.name!r}, root={str(self.root)!r}, backend={self.backend.name!r})"

    # -- inputs -------------------------------------------------------------

    @property
    def input_parameters(self) -> List[ParameterInfo]:
        """The model's input parameters (real parts for complex inputs)."""
        return list(self._inputs.real.values())

    def point(self, inputs: Dict[str, Any], settings: Optional[Settings] = None,
              sm_inputs: Optional[SMInputs] = None,
              physical_inputs: Optional[PhysicalInputs] = None,
              decay_settings: Optional[DecaySettings] = None,
              extra_blocks: Optional[Dict[str, Dict[slha.Key, slha.Value]]] = None) -> Point:
        """Create a validated point.  All input parameters that can be set
        through SLHA are required; settings and SM inputs default to the
        FlexibleSUSY defaults."""
        checked = {}
        for name, value in inputs.items():
            info = self._inputs.check_name(name)
            if info.slha is None:
                raise InputError(f"input parameter {name!r} is not read from SLHA and can't be "
                                 f"set with backend {self.backend.name!r}")
            _check_value(info, value)
            checked[name] = value
        missing = [n for n in self._inputs.settable() if n not in checked]
        if missing:
            raise InputError(f"missing input parameters: {', '.join(missing)}")
        return Point(model=self.name, inputs=checked,
                     settings=settings or Settings(),
                     sm_inputs=sm_inputs or SMInputs(),
                     physical_inputs=physical_inputs or PhysicalInputs(),
                     decay_settings=decay_settings or DecaySettings(),
                     extra_blocks=dict(extra_blocks or {}))

    def updated(self, point: Point, **changes: Any) -> Point:
        """A copy of point with some input parameters or other fields changed.

        Keyword arguments named like Point fields (settings, sm_inputs, ...)
        replace them; all others are input parameters."""
        point_fields = {f.name for f in dataclasses.fields(Point)} - {"model", "inputs"}
        inputs = dict(point.inputs)
        other = {}
        for key, value in changes.items():
            if key in point_fields:
                other[key] = value
            else:
                inputs[key] = value
        args = {f: getattr(point, f) for f in point_fields}
        args.update(other)
        return self.point(inputs, **args)

    def default_point(self) -> Point:
        """The point of the shipped models/<M>/LesHouches.in.<M>."""
        path = self.directory / f"LesHouches.in.{self.name}"
        if not path.is_file():
            raise BackendError(f"{path} not found")
        return self.point_from_slha(path.read_text())

    def point_from_slha(self, text: str) -> Point:
        """Create a point from an SLHA input file (text)."""
        data = slha.read(text)
        used = set()
        inputs = {}
        for name, info in self._inputs.real.items():
            if info.slha is None:
                continue
            value = self._read_input(data, info)
            imag = self._inputs.imag.get(name)
            if imag is not None and value is not None:
                im = self._read_input(data, imag)
                if im is not None:
                    value = value + 1j * im
            if value is not None:
                inputs[name] = value
            used.add(info.slha.block.upper())
            if imag is not None and imag.slha is not None:
                used.add(imag.slha.block.upper())

        def settings_from(cls):
            block = data.block(cls._block)
            used.add(cls._block.upper())
            if block is None:
                return cls()
            by_index = {f.metadata["index"]: f.name for f in dataclasses.fields(cls)}
            values = {by_index[k]: v for k, v in block if k in by_index}
            return cls(**values)

        settings = settings_from(Settings)
        sm_inputs = settings_from(SMInputs)
        physical_inputs = settings_from(PhysicalInputs)
        decay_settings = settings_from(DecaySettings)
        extra = {name: dict(block.entries) for name, block in data.blocks.items()
                 if name not in used}
        return self.point(inputs, settings=settings, sm_inputs=sm_inputs,
                          physical_inputs=physical_inputs, decay_settings=decay_settings,
                          extra_blocks=extra)

    @staticmethod
    def _read_input(data: slha.SLHA, info: ParameterInfo) -> Any:
        assert info.slha is not None
        block = data.block(info.slha.block)
        if block is None:
            return None
        if info.slha.index is not None:
            return block.get(info.slha.index)
        entries = {}
        for key, value in block:
            idx = key if isinstance(key, tuple) else (key,)
            entries[tuple(i - 1 for i in idx)] = value
        if not info.shape:
            return entries.get((0,))
        return arrays.from_entries(info.shape, entries)

    def to_slha_blocks(self, point: Point) -> Dict[str, Dict[slha.Key, slha.Value]]:
        """SLHA input blocks for a point."""
        if point.model != self.name:
            raise InputError(f"point of model {point.model!r} used with model {self.name!r}")
        blocks: Dict[str, Dict[slha.Key, slha.Value]] = {}
        for obj in (point.settings, point.decay_settings, point.physical_inputs, point.sm_inputs):
            entries = obj.to_slha_entries()
            if entries:
                blocks[obj._block] = entries
        for name, value in point.inputs.items():
            info = self._inputs.real[name]
            entries = _check_value(info, value)
            assert info.slha is not None
            _add_entries(blocks, info.slha, entries, lambda v: complex(v).real)
            if any(complex(v).imag != 0 for v in entries.values()):
                imag = self._inputs.imag.get(name)
                if imag is None or imag.slha is None:
                    raise InputError(f"{name}: complex values are not supported for this "
                                     "input parameter with the SLHA backend")
                _add_entries(blocks, imag.slha, entries, lambda v: complex(v).imag)
        for name, entries in point.extra_blocks.items():
            blocks.setdefault(name, {}).update(entries)
        return blocks

    # -- calculation ----------------------------------------------------------

    def calculate(self, point: Point, raise_on_problems: bool = False, timeout=None):
        """Calculate the spectrum (and decays/observables if enabled) for a point."""
        from .errors import CalculationError
        result = self.backend.calculate(point, timeout=timeout)
        if raise_on_problems and not result.ok:
            raise CalculationError(result.problems, result)
        return result
