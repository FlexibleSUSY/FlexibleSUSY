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

"""Exceptions raised by the FlexibleSUSY Python interface."""

from __future__ import annotations

from typing import List, Optional


class FlexibleSUSYError(Exception):
    """Base class of all exceptions of this package."""


class InputError(FlexibleSUSYError, ValueError):
    """Invalid input: unknown or missing parameter, wrong shape or type."""


class BackendError(FlexibleSUSYError):
    """The backend failed: missing executable or metadata, crash, unreadable output."""

    def __init__(self, message: str, stderr: Optional[str] = None):
        if stderr:
            message = f"{message}\n--- stderr ---\n{stderr.rstrip()}"
        super().__init__(message)
        self.stderr = stderr


class CalculationError(FlexibleSUSYError):
    """The calculation reported problems (only raised with raise_on_problems=True)."""

    def __init__(self, problems: List[str], result=None):
        super().__init__("; ".join(problems) if problems else "calculation failed")
        self.problems = problems
        self.result = result
