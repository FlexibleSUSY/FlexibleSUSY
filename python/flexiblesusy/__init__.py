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

"""Python interface to spectrum generators created with FlexibleSUSY.

Example::

    import flexiblesusy as fs

    mssm = fs.Model("MSSM")
    point = mssm.default_point()                      # models/MSSM/LesHouches.in.MSSM
    point = mssm.updated(point, TanBeta=20)
    res = mssm.calculate(point)
    if res.ok:
        print(res.pole_masses["hh"])                  # 0-based indices
"""

from .errors import BackendError, CalculationError, FlexibleSUSYError, InputError
from .model import Metadata, Model, Point, find_root
from .result import Channel, Decay, Decays, Result, RunningParameters
from .settings import DecaySettings, LoopLibrary, PhysicalInputs, Settings, SMInputs

__all__ = [
    "BackendError", "CalculationError", "FlexibleSUSYError", "InputError",
    "Metadata", "Model", "Point", "find_root",
    "Channel", "Decay", "Decays", "Result", "RunningParameters",
    "DecaySettings", "LoopLibrary", "PhysicalInputs", "Settings", "SMInputs",
]
