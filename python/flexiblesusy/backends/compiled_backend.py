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

"""Backend that calls the model's Python module models/<M>/flexiblesusy_<M>*.so.

The module is built by FlexibleSUSY when configured with
--with-python=<interpreter>.  The spectrum stays in memory (one session per
Result), so observables and decays are calculated on demand from it.
"""

from __future__ import annotations

import dataclasses
import importlib.machinery
import importlib.util
import sys
import warnings
from pathlib import Path
from typing import Optional

from .. import result as result_module
from .. import slha
from ..errors import BackendError, InputError
from ..settings import LoopLibrary


def module_path(model) -> Optional[Path]:
    """Path of the model's Python module for this interpreter, or None."""
    for suffix in importlib.machinery.EXTENSION_SUFFIXES:
        path = model.directory / f"flexiblesusy_{model.name}{suffix}"
        if path.is_file():
            return path
    return None


def _load_module(model):
    path = module_path(model)
    if path is None:
        others = sorted(p.name for p in model.directory.glob(f"flexiblesusy_{model.name}.*"))
        hint = (f"; found {', '.join(others)}, built for another Python version"
                if others else
                f"; configure FlexibleSUSY with --with-python={Path(sys.executable).name} "
                f"and build it")
        raise BackendError(f"the Python module of {model.name} for this interpreter "
                           f"(Python {sys.version_info.major}.{sys.version_info.minor}) "
                           f"is not built{hint}")
    name = f"flexiblesusy_{model.name}"
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, path)
    try:
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    except ImportError as e:
        raise BackendError(f"could not load {path}: {e}") from None
    sys.modules[name] = module
    return module


class CompiledBackend:
    """Calculates with the model's Python module, keeping the model in memory."""

    name = "compiled"

    def __init__(self, model):
        self.model = model
        self.module = _load_module(model)

    @staticmethod
    def _settings_for_session(point):
        """The settings run_<M>.x would use (cf. templates/run.cpp.in), so
        that both backends give the same result: observables need BSM pole
        masses, decays SM and BSM pole masses."""
        settings = point.settings
        changes = {}
        if settings.calculate_observables:
            changes["calculate_bsm_masses"] = True
        if point.decay_settings.calculate_decays:
            changes["calculate_sm_masses"] = True
            changes["calculate_bsm_masses"] = True
        return dataclasses.replace(settings, **changes) if changes else settings

    def _check_loop_library(self, point):
        requested = point.settings.loop_library
        current = self.module.Session.loop_library()
        if requested is not None and requested >= 0 and current >= 0 and int(requested) != current:
            warnings.warn(
                f"loop library {LoopLibrary(int(requested)).name} requested, but this process "
                f"already uses {LoopLibrary(current).name}: the loop library is fixed by the "
                "first calculation in a process", RuntimeWarning, stacklevel=3)

    def calculate(self, point, timeout=None):
        if timeout is not None:
            raise InputError("timeout is not supported by the compiled backend "
                             "(the calculation runs in this process)")
        point = dataclasses.replace(point, settings=self._settings_for_session(point))
        self._check_loop_library(point)
        text = slha.write(self.model.to_slha_blocks(point))
        session = self.module.Session()
        try:
            session.read_slha_input(text)
        except RuntimeError as e:
            raise InputError(f"invalid input: {e}") from None
        try:
            session.calculate_spectrum()
            output = session.slha_output()
        except RuntimeError as e:
            raise BackendError(f"spectrum calculation failed: {e}") from None
        # the output contains (empty) decay blocks even if no decays were
        # calculated: they are read by calculate_decays() instead
        result = result_module.from_slha(output, self.model.metadata, backend=self, state=session,
                                         read_decays=False)
        # as with run_<M>.x: observables and decays requested in the
        # settings are calculated right away
        if point.settings.calculate_observables and result.ok:
            result.observables()
        if point.decay_settings.calculate_decays and result.ok:
            result.decays()
        return result

    def _session(self, result):
        if result._state is None:
            raise BackendError("the model of this result was released (Result.release())")
        return result._state

    def calculate_observables(self, result):
        session = self._session(result)
        try:
            session.calculate_observables()
        except RuntimeError as e:
            raise BackendError(f"calculation of observables failed: {e}") from None
        return result_module.Observables(session.observable_names(), session.observable_values(),
                                         session.observable_problems())

    def calculate_decays(self, result):
        session = self._session(result)
        try:
            session.calculate_decays()
            output = session.slha_output()
        except RuntimeError as e:
            raise BackendError(f"calculation of decays failed: {e}") from None
        decays = result_module.decays_from_slha(slha.read(output), self.model.metadata)
        decays.problems = list(session.decay_problems())
        decays.warnings = list(session.decay_warnings())
        return decays
