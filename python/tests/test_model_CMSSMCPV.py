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

"""Tests with the generated CMSSMCPV (CP-violating Higgs sector, complex inputs);
skipped if models/CMSSMCPV isn't built."""

import dataclasses
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import flexiblesusy as fs  # noqa: E402
from flexiblesusy import slha  # noqa: E402
from flexiblesusy.backends import compiled_backend  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
MODEL_DIR = ROOT / "models" / "CMSSMCPV"
HAVE_MODEL = (MODEL_DIR / "run_CMSSMCPV.x").is_file() and (MODEL_DIR / "CMSSMCPV_metadata.json").is_file()


def have_module():
    return HAVE_MODEL and compiled_backend.module_path(
        fs.Model("CMSSMCPV", root=ROOT, backend="slha")) is not None


@unittest.skipUnless(have_module(), "CMSSMCPV and its Python module not built for this interpreter")
class CMSSMCPVTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.compiled = fs.Model("CMSSMCPV", root=ROOT, backend="compiled")
        cls.slha_model = fs.Model("CMSSMCPV", root=ROOT, backend="slha")
        p = cls.compiled.default_point()
        cls.point = cls.compiled.updated(
            p, settings=dataclasses.replace(p.settings, loop_library=fs.LoopLibrary.COLLIER),
            decay_settings=fs.DecaySettings(calculate_decays=True, print_effc_block=True))

    def test_inputs(self):
        names = [p.name for p in self.compiled.input_parameters]
        # inputs of their own in IMMINPAR, not imaginary parts of others
        for name in ("Imm12", "SinPhiMu", "ImAzero"):
            self.assertIn(name, names)
        # PhaseMu isn't in the shipped input file: zero, as with run_CMSSMCPV.x
        self.assertEqual(self.point.inputs["PhaseMu"], 0j)

    def test_backends_agree_and_no_effc_block(self):
        res = self.compiled.calculate(self.point)
        ref = self.slha_model.calculate(self.point)
        self.assertTrue(res.ok, res.problems)
        self.assertEqual(res.blocks.block("MASS").entries, ref.blocks.block("MASS").entries)
        # hh(1) is the neutral Goldstone boson (index 0)
        self.assertEqual(res.decays()["hh", 1].width, ref.decays()["hh", 1].width)
        # as run_CMSSMCPV.x: no EFFHIGGSCOUPLINGS for a CP-violating Higgs sector
        self.assertNotIn("EFFHIGGSCOUPLINGS", ref.blocks)
        output = slha.read(self.compiled.backend._session(res).slha_output())
        self.assertNotIn("EFFHIGGSCOUPLINGS", output)


if __name__ == "__main__":
    unittest.main()
