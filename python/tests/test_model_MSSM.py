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

"""Tests with the generated MSSM; skipped if models/MSSM isn't built.

The tests run with the SLHA backend and, if the Python module of the MSSM
is built for this interpreter, with the compiled backend.
"""

import dataclasses
import math
import os
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import flexiblesusy as fs  # noqa: E402
from flexiblesusy import slha  # noqa: E402
from flexiblesusy.backends import compiled_backend  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
MODEL_DIR = ROOT / "models" / "MSSM"
EXE = MODEL_DIR / "run_MSSM.x"
METADATA = MODEL_DIR / "MSSM_metadata.json"
HAVE_MSSM = EXE.is_file() and METADATA.is_file()


def have_module():
    return HAVE_MSSM and compiled_backend.module_path(fs.Model("MSSM", root=ROOT,
                                                              backend="slha")) is not None


def run_directly(input_file: Path) -> slha.SLHA:
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "out.slha"
        subprocess.run([str(EXE), f"--slha-input-file={input_file}",
                        f"--slha-output-file={out}"], cwd=str(ROOT), capture_output=True)
        return slha.read(out.read_text())


def with_decays(model, point):
    return model.updated(
        point,
        settings=dataclasses.replace(point.settings, loop_library=fs.LoopLibrary.COLLIER),
        decay_settings=fs.DecaySettings(calculate_decays=True, min_br_to_print=0))


class MSSMBackendTests:
    """Tests run with each backend (BACKEND set by the subclasses)."""

    BACKEND = ""

    @classmethod
    def setUpClass(cls):
        cls.model = fs.Model("MSSM", root=ROOT, backend=cls.BACKEND)
        cls.point = cls.model.default_point()

    def test_backend(self):
        self.assertEqual(self.model.backend.name, self.BACKEND)

    def test_metadata(self):
        names = [p.name for p in self.model.input_parameters]
        self.assertIn("TanBeta", names)
        self.assertIn("mq2Input", names)
        particles = {p.name: p for p in self.model.metadata.particles}
        self.assertEqual(particles["hh"].multiplicity, 2)
        self.assertEqual(particles["hh"].pdg, (25, 35))

    def test_same_as_executable(self):
        """The spectrum equals the one of run_MSSM.x for the same input."""
        with tempfile.TemporaryDirectory() as tmp:
            input_file = Path(tmp) / "in.slha"
            input_file.write_text(slha.write(self.model.to_slha_blocks(self.point)))
            direct = run_directly(input_file)
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
        self.assertTrue(any("convergence" in p for p in res.problems), res.problems)
        with self.assertRaises(fs.CalculationError):
            self.model.calculate(point, raise_on_problems=True)

    def test_decays_requested(self):
        res = self.model.calculate(with_decays(self.model, self.point))
        self.assertTrue(res.ok, res.problems)
        h = res.decays()["hh", 0]
        self.assertGreater(h.width, 0)
        self.assertAlmostEqual(sum(c.br for c in h.channels), 1.0, places=6)
        bb = [c for c in h.channels if sorted(c.pdg) == [-5, 5]]
        self.assertEqual(len(bb), 1)


@unittest.skipUnless(HAVE_MSSM, "MSSM not generated and built")
class MSSMSLHABackendTest(MSSMBackendTests, unittest.TestCase):
    BACKEND = "slha"

    def test_roundtrip_of_shipped_input(self):
        """Reading, validating and writing the shipped input changes nothing."""
        direct = run_directly(MODEL_DIR / "LesHouches.in.MSSM")
        res = self.model.calculate(self.model.default_point())
        self.assertEqual(res.blocks.block("MASS").entries, direct.block("MASS").entries)


@unittest.skipUnless(have_module(), "Python module of the MSSM not built for this interpreter")
class MSSMCompiledBackendTest(MSSMBackendTests, unittest.TestCase):
    BACKEND = "compiled"

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.slha_model = fs.Model("MSSM", root=ROOT, backend="slha")
        # decays need Collier, and the loop library is fixed per process: set
        # it in the first calculation
        cls.point = cls.model.updated(
            cls.point,
            settings=dataclasses.replace(cls.point.settings, loop_library=fs.LoopLibrary.COLLIER))

    def test_auto_selects_compiled(self):
        self.assertEqual(fs.Model("MSSM", root=ROOT).backend.name, "compiled")

    def test_same_as_slha_backend(self):
        res = self.model.calculate(self.point)
        ref = self.slha_model.calculate(self.point)
        self.assertEqual(res.blocks.block("MASS").entries, ref.blocks.block("MASS").entries)
        self.assertEqual(res.problems, ref.problems)

    def test_observables_on_demand(self):
        res = self.model.calculate(self.point)
        self.assertIsNone(res._observables)
        masses = dict(res.blocks.block("MASS").entries)

        start = time.perf_counter()
        obs = res.observables()
        on_demand = time.perf_counter() - start
        start = time.perf_counter()
        self.model.calculate(self.point)
        spectrum = time.perf_counter() - start
        # the spectrum isn't recalculated: much faster than a spectrum
        self.assertLess(on_demand, 0.5 * spectrum, (on_demand, spectrum))

        self.assertTrue(obs.ok, obs.problems)
        self.assertIn("amm_Fe_1", obs)
        self.assertNotEqual(obs["amm_Fe_1"], 0.0)
        self.assertIs(res.observables(), obs)                    # cached
        self.assertEqual(dict(res.blocks.block("MASS").entries), masses)

    def test_observables_same_as_in_one_run(self):
        later = self.model.calculate(self.point).observables()
        point = self.model.updated(
            self.point, settings=dataclasses.replace(self.point.settings, calculate_observables=True))
        res = self.model.calculate(point)
        self.assertIsNotNone(res._observables)                    # calculated right away
        self.assertEqual(res.observables()["amm_Fe_1"], later["amm_Fe_1"])

    def test_decays_on_demand(self):
        # same channels printed as with with_decays(), but not calculated right
        # away; decays need SM pole masses in the spectrum
        point = self.model.updated(
            self.point,
            settings=dataclasses.replace(self.point.settings, calculate_sm_masses=True),
            decay_settings=fs.DecaySettings(min_br_to_print=0))
        res = self.model.calculate(point)
        self.assertIsNone(res._decays)
        decays = res.decays()
        ref = self.slha_model.calculate(with_decays(self.slha_model, self.point)).decays()
        self.assertEqual(len(decays), len(ref))
        for d in ref:
            mine = decays.by_pdg(d.pdg)
            self.assertEqual(mine.width, d.width, d.particle)
            self.assertEqual(sorted((c.pdg, c.br) for c in mine.channels),
                             sorted((c.pdg, c.br) for c in d.channels))

    def test_decays_need_sm_masses(self):
        point = self.model.updated(
            self.point,
            settings=dataclasses.replace(self.point.settings, calculate_sm_masses=False),
            decay_settings=fs.DecaySettings())
        with self.assertRaisesRegex(fs.BackendError, "SM and BSM pole masses"):
            self.model.calculate(point).decays()

    def test_release(self):
        res = self.model.calculate(self.point)
        res.release()
        with self.assertRaisesRegex(fs.BackendError, "released"):
            res.observables()


HB = os.environ.get("HIGGSBOUNDS_DATASET")
HS = os.environ.get("HIGGSSIGNALS_DATASET")
LILITH = os.environ.get("LILITH_DATABASE")


def point_125(model):
    """MSSM point with mh = 124.8 GeV (Lilith needs a Higgs between 123 and 128 GeV)."""
    import numpy as np
    p = model.default_point()
    ms = 3000
    m2 = ms**2 * np.eye(3)
    return model.updated(
        p, TanBeta=20, mHd2IN=ms**2, mHu2IN=ms**2, mq2Input=m2, mu2Input=m2, md2Input=m2,
        ml2Input=m2, me2Input=m2, Auij=np.zeros((3, 3)),
        MassBInput=1000, MassWBInput=1000, MassGInput=3000,
        settings=dataclasses.replace(p.settings, loop_library=fs.LoopLibrary.COLLIER,
                                     calculate_sm_masses=True),
        decay_settings=fs.DecaySettings(calculate_decays=True, call_higgstools=True,
                                        call_lilith=True))


@unittest.skipUnless(HB and HS and LILITH, "set HIGGSBOUNDS_DATASET, HIGGSSIGNALS_DATASET "
                     "and LILITH_DATABASE to test HiggsTools and Lilith")
@unittest.skipUnless(have_module(), "Python module of the MSSM not built for this interpreter")
class MSSMHiggsToolsLilithTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        paths = dict(higgsbounds_dataset=HB, higgssignals_dataset=HS, lilith_database=LILITH)
        cls.model = fs.Model("MSSM", root=ROOT, backend="compiled", **paths)
        cls.slha_model = fs.Model("MSSM", root=ROOT, backend="slha", **paths)
        if not cls.model.backend.module.Session.has_higgstools():
            raise unittest.SkipTest("FlexibleSUSY configured without HiggsTools")
        cls.point = point_125(cls.model)

    def test_same_as_slha_backend(self):
        # on demand: spectrum first, without decays or HiggsTools
        point = self.model.updated(self.point, decay_settings=fs.DecaySettings())
        res = self.model.calculate(point)
        ht = res.higgstools()
        ref = self.slha_model.calculate(self.point)
        ref_ht = ref.higgstools()
        self.assertAlmostEqual(ht.higgssignals.chi2, ref_ht.higgssignals.chi2, places=4)
        self.assertEqual(ht.higgssignals.ndof, ref_ht.higgssignals.ndof)
        self.assertEqual([b.pdg for b in ht.higgsbounds], [b.pdg for b in ref_ht.higgsbounds])
        self.assertAlmostEqual(ht.higgsbounds[0].obs_ratio, ref_ht.higgsbounds[0].obs_ratio,
                               places=6)
        ll = res.lilith()
        ref_ll = ref.lilith()
        self.assertIsNotNone(ll)
        self.assertEqual(ll.ndof, ref_ll.ndof)
        self.assertAlmostEqual(ll.chi2, ref_ll.chi2, places=4)

    def test_lilith_leaves_the_interpreter_alone(self):
        import __main__
        flags = sys.getdlopenflags()
        names = set(vars(__main__))
        res = self.model.calculate(self.model.updated(self.point, decay_settings=fs.DecaySettings()))
        self.assertIsNotNone(res.lilith())
        self.assertEqual(sys.getdlopenflags(), flags)
        self.assertEqual(set(vars(__main__)), names)

    def test_lilith_without_higgs_near_125(self):
        res = self.model.calculate(dataclasses.replace(self.model.default_point(), settings=dataclasses.replace(
            self.model.default_point().settings, loop_library=fs.LoopLibrary.COLLIER,
            calculate_sm_masses=True)))
        self.assertIsNone(res.lilith())

    def test_bad_lilith_database(self):
        res = self.model.calculate(self.model.updated(self.point, decay_settings=fs.DecaySettings()))
        with self.assertRaisesRegex(fs.BackendError, "is not a file"):
            res.lilith("/nonexistent/latest.list")


if __name__ == "__main__":
    unittest.main()
