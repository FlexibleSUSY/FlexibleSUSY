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

"""Backend that runs the spectrum generator executable models/<M>/run_<M>.x."""

from __future__ import annotations

import os
import subprocess
import tempfile
from pathlib import Path

from .. import result as result_module
from .. import slha
from ..errors import BackendError


class SLHABackend:
    """Writes an SLHA input file, runs run_<M>.x and reads its SLHA output.

    run_<M>.x exits with a non-zero code when the point has problems; it
    still writes the SLHA output with the problems in Block SPINFO.  So the
    exit code alone isn't treated as an error: only a missing or unreadable
    output is.
    """

    name = "slha"

    def __init__(self, model):
        self.model = model
        self.executable = model.directory / f"run_{model.name}.x"

    def _check_executable(self) -> None:
        if not self.executable.is_file():
            raise BackendError(f"{self.executable} not found; build it with "
                               f"'make {self.executable.relative_to(self.model.root)}' "
                               f"in {self.model.root}")
        if not os.access(self.executable, os.X_OK):
            raise BackendError(f"{self.executable} is not executable")

    def input_text(self, point) -> str:
        return slha.write(self.model.to_slha_blocks(point))

    def run(self, input_text: str, timeout=None) -> str:
        """Run the executable on SLHA input text; return the SLHA output text."""
        self._check_executable()
        with tempfile.TemporaryDirectory(prefix="flexiblesusy-") as tmp:
            input_file = Path(tmp) / "input.slha"
            output_file = Path(tmp) / "output.slha"
            input_file.write_text(input_text)
            try:
                proc = subprocess.run(
                    [str(self.executable),
                     f"--slha-input-file={input_file}",
                     f"--slha-output-file={output_file}"],
                    cwd=str(self.model.root), capture_output=True, text=True,
                    timeout=timeout)
            except subprocess.TimeoutExpired as e:
                raise BackendError(f"{self.executable.name} timed out after {timeout} s",
                                   e.stderr if isinstance(e.stderr, str) else None) from None
            if proc.returncode < 0:
                raise BackendError(f"{self.executable.name} was killed by signal "
                                   f"{-proc.returncode}", proc.stderr)
            if not output_file.is_file() or output_file.stat().st_size == 0:
                raise BackendError(f"{self.executable.name} wrote no SLHA output "
                                   f"(exit code {proc.returncode})", proc.stderr)
            text = output_file.read_text()
        if "SPINFO" not in slha.read(text):
            raise BackendError(f"{self.executable.name} wrote no Block SPINFO "
                               f"(exit code {proc.returncode})", proc.stderr)
        return text

    def calculate(self, point, timeout=None):
        output = self.run(self.input_text(point), timeout=timeout)
        return result_module.from_slha(output, self.model.metadata, backend=self, state=point)

    def calculate_decays(self, result):
        """Decays calculated in the same run; this backend can't add them later."""
        point = result._state
        if point is not None and point.decay_settings.calculate_decays:
            data = result.blocks
            return result_module.Decays([], problems=result_module._info_block(data, "DCINFO", 4),
                                        warnings=result_module._info_block(data, "DCINFO", 3))
        raise BackendError(
            "decays were not calculated: with the SLHA backend, enable them before the "
            "calculation (decay_settings=DecaySettings(calculate_decays=True)), or use the "
            "compiled backend, which calculates them on demand")

    def calculate_observables(self, result):
        raise BackendError("observables are only available with the compiled backend "
                           "(configure FlexibleSUSY with --with-python=<interpreter>)")
