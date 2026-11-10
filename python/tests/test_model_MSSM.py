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

"""Tests with the generated MSSM; skipped if models/MSSM isn't built."""

import math
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import flexiblesusy as fs  # noqa: E402
from flexiblesusy import slha  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
MODEL_DIR = ROOT / "models" / "MSSM"
EXE = MODEL_DIR / "run_MSSM.x"
METADATA = MODEL_DIR / "MSSM_metadata.json"


def run_directly(input_file: Path) -> slha.SLHA:
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "out.slha"
        subprocess.run([str(EXE), f"--slha-input-file={input_file}",
                        f"--slha-output-file={out}"], cwd=str(ROOT), capture_output=True)
        return slha.read(out.read_text())


@unittest.skipUnless(EXE.is_file() and METADATA.is_file(), "MSSM not generated and built")
class MSSMTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = fs.Model("MSSM", root=ROOT)
        cls.point = cls.model.default_point()

    def test_metadata(self):
        names = [p.name for p in self.model.input_parameters]
        self.assertIn("TanBeta", names)
        self.assertIn("mq2Input", names)
        particles = {p.name: p for p in self.model.metadata.particles}
        self.assertEqual(particles["hh"].multiplicity, 2)
        self.assertEqual(particles["hh"].pdg, (25, 35))

    def test_same_as_executable(self):
        """Python round trip of the shipped input gives the executable's spectrum."""
        direct = run_directly(MODEL_DIR / "LesHouches.in.MSSM")
        res = self.model.calculate(self.point)
        self.assertTrue(res.ok, res.problems)
        mass = direct.block("MASS")
        self.assertEqual(sorted(mass.keys()), sorted(res.blocks.block("MASS").keys()))
        for key, value in mass:
            self.assertEqual(res.blocks.block("MASS")[key], value, key)
        self.assertAlmostEqual(res.pole_masses["hh"][0], mass[25])
        self.assertAlmostEqual(res.mixings["ZH"][0][1], direct.block("SCALARMIX")[(1, 2)])

    def test_gauge_couplings_not_gut_normalized(self):
        res = self.model.calculate(self.point)
        gY = res.blocks.block("GAUGE")[1]
        self.assertAlmostEqual(res.parameters["g1"], gY * math.sqrt(5.0 / 3.0))
        self.assertGreater(res.parameters.scale, 0)

    def test_problems(self):
        # one iteration with a very tight precision goal: no convergence
        point = self.model.updated(self.point,
                                   settings=fs.Settings(max_iterations=1, precision=1e-15))
        res = self.model.calculate(point)
        self.assertFalse(res.ok)
        self.assertTrue(res.problems)
        with self.assertRaises(fs.CalculationError):
            self.model.calculate(point, raise_on_problems=True)

    def test_decays(self):
        point = self.model.updated(
            self.point,
            settings=fs.Settings(loop_library=fs.LoopLibrary.COLLIER),
            decay_settings=fs.DecaySettings(calculate_decays=True, min_br_to_print=0))
        res = self.model.calculate(point)
        self.assertTrue(res.ok, res.problems)
        h = res.decays["hh", 0]
        self.assertGreater(h.width, 0)
        self.assertAlmostEqual(sum(c.br for c in h.channels), 1.0, places=6)
        bb = [c for c in h.channels if sorted(c.pdg) == [-5, 5]]
        self.assertEqual(len(bb), 1)


if __name__ == "__main__":
    unittest.main()
