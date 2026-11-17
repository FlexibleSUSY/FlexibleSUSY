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

"""Backends that perform the calculation for a Model.

Each backend has a name and a method calculate(point, timeout) that
returns a flexiblesusy.result.Result.
"""

from ..errors import InputError


class Backend:
    """Common part of the backends: default paths for HiggsTools and Lilith."""

    name = ""

    def __init__(self, model):
        self.model = model

    def higgstools_datasets(self, higgsbounds_dataset=None, higgssignals_dataset=None):
        hb = str(higgsbounds_dataset) if higgsbounds_dataset else self.model.higgsbounds_dataset
        hs = str(higgssignals_dataset) if higgssignals_dataset else self.model.higgssignals_dataset
        if not hb or not hs:
            raise InputError("HiggsTools needs the HiggsBounds and HiggsSignals datasets: pass "
                             "them to Model(higgsbounds_dataset=..., higgssignals_dataset=...) "
                             "or to Result.higgstools()")
        return hb, hs

    def lilith_database(self, database=None):
        db = str(database) if database else self.model.lilith_database
        if not db:
            raise InputError("Lilith needs a database: pass it to Model(lilith_database=...) "
                             "or to Result.lilith()")
        return db
