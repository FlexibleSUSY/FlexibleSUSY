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

"""Unit tests that don't need a generated model (run: python3 -m unittest)."""

import json
import math
import os
import stat
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import flexiblesusy as fs  # noqa: E402
from flexiblesusy import arrays, slha  # noqa: E402

METADATA = {
    "format_version": 1,
    "model": "Toy",
    "flexiblesusy_version": "2.9.0",
    "input_parameters": [
        {"name": "TanBeta", "kind": "real", "shape": [], "slha": {"block": "MINPAR", "index": 3}},
        {"name": "SignMu", "kind": "integer", "shape": [], "slha": {"block": "MINPAR", "index": 4}},
        {"name": "MuInput", "kind": "complex", "shape": [], "slha": {"block": "EXTPAR", "index": 23}},
        {"name": "MuInput", "kind": "complex", "shape": [], "slha": {"block": "IMEXTPAR", "index": 23}},
        {"name": "mq2Input", "kind": "real", "shape": [2, 2], "slha": {"block": "MSQ2IN", "index": None}},
        {"name": "SinPhi", "kind": "real", "shape": [], "slha": {"block": "IMMINPAR", "index": 7}},
        {"name": "Phase", "kind": "complex", "shape": [], "slha": {"block": "MINPAR", "index": 100}},
        {"name": "Internal", "kind": "real", "shape": [], "slha": None},
    ],
    "particles": [
        {"name": "hh", "multiplicity": 2, "pdg": [25, 35], "sm": False},
        {"name": "Ah", "multiplicity": 2, "pdg": [0, 36], "sm": False},
        {"name": "Fd", "multiplicity": 1, "pdg": [5], "sm": True},
    ],
    "mixing_matrices": [
        {"name": "ZH", "kind": "real", "shape": [2, 2], "slha": {"block": "SCALARMIX", "index": None}},
    ],
    "model_parameters": [
        {"name": "g1", "kind": "real", "shape": [], "slha": {"block": "gauge", "index": 1},
         "slha_factor": 0.7745966692414834},
        {"name": "Yd", "kind": "real", "shape": [2, 2], "slha": {"block": "Yd", "index": None}},
    ],
}

OUTPUT = textwrap.dedent("""\
    Block SPINFO
         1   FlexibleSUSY
         2   2.9.0
         3   a warning with 3 words
         4   no convergence
         5   Toy
    Block MASS
            25     1.25000000E+02   # hh(1)
            35     5.00000000E+02   # hh(2)
            36     5.10000000E+02   # Ah(2)
             5     3.00000000E+00   # Fd(1)
    Block SCALARMIX
      1  1     1.00000000E-01   # ZH(1,1)
      1  2     9.94987437E-01   # ZH(1,2)
      2  1     9.94987437E-01   # ZH(2,1)
      2  2    -1.00000000E-01   # ZH(2,2)
    Block gauge Q= 1.00000000E+03
         1     3.62397102E-01   # g1 * 0.7745966692414834
    Block Yd Q= 1.00000000E+03
      1  1     1.00000000E-03   # Yd(1,1)
      2  2     2.00000000E-02   # Yd(2,2)
    DECAY        25     4.00000000E-03   # hh(1) decays
         7.50000000E-01   2          -5         5  # BR(hh(1) -> barFd(3) Fd(3))
         2.50000000E-01   2          21        21  # BR(hh(1) -> VG VG)
    """)


class ToyModelTest(unittest.TestCase):
    """A fake FlexibleSUSY directory with metadata and a fake executable."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        root = Path(self.tmp.name)
        model_dir = root / "models" / "Toy"
        model_dir.mkdir(parents=True)
        (model_dir / "Toy_metadata.json").write_text(json.dumps(METADATA))
        (model_dir / "LesHouches.in.Toy").write_text(textwrap.dedent("""\
            Block MODSEL
                6   0
            Block FlexibleSUSY
                4   2      # pole mass loop order
               31   0      # loop library
            Block SMINPUTS
                3   1.181e-01
            Block MINPAR
                3   10.0
                4   1
            Block EXTPAR
               23   500.0
            Block IMMINPAR
                7   0.25
            Block IMEXTPAR
               23   5.0
            Block MSQ2IN
              1 1   1.0e6
              1 2   0.0
              2 1   0.0
              2 2   2.0e6
            """))
        # fake executable: copies a fixed output, records its input
        exe = model_dir / "run_Toy.x"
        (model_dir / "output.slha").write_text(OUTPUT)
        exe.write_text(textwrap.dedent(f"""\
            #!/bin/sh
            for a in "$@"; do
              case "$a" in
                --slha-input-file=*) cp "${{a#*=}}" "{model_dir}/last_input.slha" ;;
                --slha-output-file=*) cp "{model_dir}/output.slha" "${{a#*=}}" ;;
              esac
            done
            exit 1
            """))
        exe.chmod(exe.stat().st_mode | stat.S_IXUSR)
        self.root = root
        self.model_dir = model_dir
        self.model = fs.Model("Toy", root=root)

    def tearDown(self):
        self.tmp.cleanup()

    def inputs(self, **changes):
        inputs = {"TanBeta": 10.0, "SignMu": 1, "MuInput": 500.0,
                  "mq2Input": [[1e6, 0.0], [0.0, 2e6]], "SinPhi": 0.5, "Phase": 1.0}
        inputs.update(changes)
        return inputs

    # -- inputs ---------------------------------------------------------------

    def test_missing_input_is_an_error(self):
        inputs = self.inputs()
        del inputs["TanBeta"]
        with self.assertRaisesRegex(fs.InputError, "missing input parameters: TanBeta"):
            self.model.point(inputs)

    def test_unknown_input_suggests_names(self):
        with self.assertRaisesRegex(fs.InputError, "did you mean TanBeta"):
            self.model.point(self.inputs(TanBta=10))

    def test_wrong_shape(self):
        with self.assertRaisesRegex(fs.InputError, r"expected shape \(2, 2\)"):
            self.model.point(self.inputs(mq2Input=[1.0, 2.0]))

    def test_integer_input(self):
        with self.assertRaisesRegex(fs.InputError, "expected an integer"):
            self.model.point(self.inputs(SignMu=0.5))

    def test_input_not_in_slha(self):
        with self.assertRaisesRegex(fs.InputError, "not read from SLHA"):
            self.model.point(self.inputs(Internal=1.0))

    def test_settings_default_to_flexiblesusy(self):
        blocks = self.model.to_slha_blocks(self.model.point(self.inputs()))
        self.assertNotIn("FlexibleSUSY", blocks)
        self.assertNotIn("SMINPUTS", blocks)

    def test_slha_blocks(self):
        point = self.model.point(
            self.inputs(MuInput=500 + 5j),
            settings=fs.Settings(loop_library=fs.LoopLibrary.COLLIER, pole_mass_loop_order=2),
            sm_inputs=fs.SMInputs(alpha_s_mz=0.118))
        blocks = self.model.to_slha_blocks(point)
        self.assertEqual(blocks["FlexibleSUSY"], {31: 1.0, 4: 2.0})
        self.assertEqual(blocks["SMINPUTS"], {3: 0.118})
        self.assertEqual(blocks["MINPAR"], {3: 10.0, 4: 1.0, 100: 1.0})
        self.assertEqual(blocks["EXTPAR"], {23: 500.0})
        self.assertEqual(blocks["IMEXTPAR"], {23: 5.0})
        self.assertEqual(blocks["MSQ2IN"], {(1, 1): 1e6, (1, 2): 0.0, (2, 1): 0.0, (2, 2): 2e6})

    def test_default_point_roundtrip(self):
        point = self.model.default_point()
        self.assertEqual(point.inputs["TanBeta"], 10.0)
        self.assertEqual(point.inputs["MuInput"], 500 + 5j)
        self.assertEqual(point.settings.loop_library, 0.0)
        self.assertEqual(point.sm_inputs.alpha_s_mz, 0.1181)
        self.assertEqual(point.extra_blocks, {"MODSEL": {6: 0.0}})
        blocks = self.model.to_slha_blocks(point)
        self.assertEqual(blocks["MSQ2IN"][(2, 2)], 2e6)
        self.assertEqual(blocks["MODSEL"], {6: 0.0})

    def test_separate_input_in_im_block(self):
        # SinPhi is an input of its own in IMMINPAR, not an imaginary part
        self.assertIn("SinPhi", [p.name for p in self.model.input_parameters])
        blocks = self.model.to_slha_blocks(self.model.point(self.inputs()))
        self.assertEqual(blocks["IMMINPAR"], {7: 0.5})
        self.assertEqual(blocks["MINPAR"][100], 1.0)
        self.assertEqual(self.model.default_point().inputs["SinPhi"], 0.25)

    def test_default_point_missing_input_is_zero(self):
        # Phase is not in LesHouches.in.Toy: zero, as with run_<M>.x
        self.assertEqual(self.model.default_point().inputs["Phase"], 0j)
        inputs = self.inputs()
        del inputs["Phase"]
        with self.assertRaisesRegex(fs.InputError, "missing input parameters: Phase"):
            self.model.point(inputs)

    def test_updated(self):
        point = self.model.updated(self.model.default_point(), TanBeta=20,
                                   settings=fs.Settings(precision=1e-5))
        self.assertEqual(point.inputs["TanBeta"], 20)
        self.assertEqual(point.settings.precision, 1e-5)
        with self.assertRaises(fs.InputError):
            self.model.updated(point, TanBta=20)

    # -- calculation ----------------------------------------------------------

    def test_calculate(self):
        res = self.model.calculate(self.model.point(self.inputs()))
        written = slha.read((self.model_dir / "last_input.slha").read_text())
        self.assertEqual(written.block("MINPAR")[3], 10.0)

        self.assertFalse(res.ok)
        self.assertEqual(res.problems, ["no convergence"])
        self.assertEqual(res.warnings, ["a warning with 3 words"])
        self.assertEqual(list(res.pole_masses["hh"]), [125.0, 500.0])
        ah = list(res.pole_masses["Ah"])
        self.assertTrue(math.isnan(ah[0]))      # Goldstone: PDG code 0
        self.assertEqual(ah[1], 510.0)
        self.assertAlmostEqual(res.mixings["ZH"][1][0], 0.994987437)
        self.assertAlmostEqual(res.parameters["g1"], 0.362397102 / 0.7745966692414834)
        self.assertEqual(res.parameters.scale, 1000.0)
        self.assertEqual(res.parameters["Yd"][1][1], 0.02)
        self.assertEqual(res.parameters["Yd"][0][1], 0.0)
        self.assertEqual(self.model.backend.name, "slha")
        decay = res.decays()["hh", 0]
        self.assertIs(res.decays().by_pdg(25), decay)
        self.assertEqual(decay.width, 4e-3)
        self.assertEqual(decay.channels[0].pdg, (-5, 5))
        self.assertAlmostEqual(decay.channels[0].width, 3e-3)
        self.assertEqual(decay.channels[0].final_state, ("barFd(3)", "Fd(3)"))
        json.dumps(res.to_dict())

    def test_slha_backend_on_demand(self):
        (self.model_dir / "output.slha").write_text(OUTPUT.split("DECAY")[0])
        res = self.model.calculate(self.model.point(self.inputs()))
        with self.assertRaisesRegex(fs.BackendError, "decays were not calculated"):
            res.decays()
        with self.assertRaisesRegex(fs.BackendError, "only available with the compiled backend"):
            res.observables()
        point = self.model.point(self.inputs(),
                                 decay_settings=fs.DecaySettings(calculate_decays=True))
        self.assertEqual(len(self.model.calculate(point).decays()), 0)

    def test_slha_backend_higgstools_not_called(self):
        model = fs.Model("Toy", root=self.root, higgsbounds_dataset="hb",
                         higgssignals_dataset="hs", lilith_database="db")
        res = model.calculate(model.point(self.inputs()))
        with self.assertRaisesRegex(fs.BackendError, "HiggsTools was not called"):
            res.higgstools()
        with self.assertRaisesRegex(fs.InputError, "needs a database"):
            self.model.calculate(self.model.point(self.inputs())).lilith()
        point = model.point(self.inputs(), decay_settings=fs.DecaySettings(call_lilith=True))
        res = model.calculate(point)
        self.assertIsNone(res.lilith())                      # no Block LILITH in the output
        with self.assertRaisesRegex(fs.InputError, "can't call it again"):
            res.lilith("other.list")

    def test_compiled_backend_missing(self):
        with self.assertRaisesRegex(fs.BackendError, "is not built"):
            fs.Model("Toy", root=self.root, backend="compiled")

    def test_raise_on_problems(self):
        with self.assertRaises(fs.CalculationError) as ctx:
            self.model.calculate(self.model.point(self.inputs()), raise_on_problems=True)
        self.assertEqual(ctx.exception.problems, ["no convergence"])

    def test_backend_error_without_output(self):
        (self.model_dir / "output.slha").write_text("")
        with self.assertRaisesRegex(fs.BackendError, "wrote no SLHA output"):
            self.model.calculate(self.model.point(self.inputs()))

    def test_missing_model(self):
        with self.assertRaisesRegex(fs.BackendError, "available: Toy"):
            fs.Model("NoSuchModel", root=self.root)


class SLHATest(unittest.TestCase):
    def test_read_entries(self):
        data = slha.read(OUTPUT)
        self.assertEqual(data.block("spinfo").get_all(3), ["a warning with 3 words"])
        self.assertEqual(data.block("SCALARMIX")[(2, 2)], -0.1)
        self.assertEqual(data.block("GAUGE").scale, 1000.0)
        self.assertEqual(data.decays[25].channels[1].pdg, (21, 21))

    def test_spinfo_text_starting_with_a_number(self):
        data = slha.read("Block SPINFO\n   4   2 iterations failed\n")
        self.assertEqual(data.block("SPINFO")[4], "2 iterations failed")

    def test_write_read_roundtrip(self):
        text = slha.write({"MINPAR": {3: 10.0}, "MSQ2IN": {(1, 2): 1.0 / 3.0}})
        data = slha.read(text)
        self.assertEqual(data.block("MINPAR")[3], 10.0)
        self.assertEqual(data.block("MSQ2IN")[(1, 2)], 1.0 / 3.0)   # full precision

    def test_write_rejects_nan(self):
        with self.assertRaises(ValueError):
            slha.write({"MINPAR": {3: float("nan")}})


HIGGSTOOLS_OUTPUT = textwrap.dedent("""\
    Block HIGGSSIGNALS
         1     1.59000000E+02   # number of degrees of freedom
         2     1.53787720E+02   # chi2
         3     1.52552773E+02   # SM chi2 for mh = 125.090000 GeV
         4     2.66447195E-01   # p-value (1d)
         5     5.39305381E-01   # p-value (2d)
    Block HIGGSBOUNDS
     25  1     6.44167889E-01   # LHC8 [HW,HZ]>[bb] from 1409.6212 (ATLAS 25fb-1, M=(110, 140))
     25  2     1.10349794E+00   # expRatio
    Block LILITH
         1     6.60000000E+01   # number of degrees of freedom
         2     5.88464012E+01   # chi2
         3     5.47819557E+01   # SM chi2 for mh = 125.090000 GeV
         4     4.37950616E-02   # p-value (1d)
         5     1.31043920E-01   # p-value (2d)
    """)


class HiggsToolsSLHATest(unittest.TestCase):
    def test_parse(self):
        from flexiblesusy import result
        data = slha.read(HIGGSTOOLS_OUTPUT)
        hs = result.signal_result_from_slha(data, "HIGGSSIGNALS")
        self.assertEqual((hs.ndof, hs.chi2, hs.mh_ref), (159, 153.78772, 125.09))
        self.assertEqual(hs.pval_2d, 0.539305381)
        hb = result.higgsbounds_from_slha(data)
        self.assertEqual(len(hb), 1)
        self.assertEqual((hb[0].pdg, hb[0].obs_ratio, hb[0].exp_ratio), (25, 0.644167889, 1.10349794))
        self.assertTrue(hb[0].analysis.startswith("LHC8"))
        self.assertEqual(result.signal_result_from_slha(data, "LILITH").ndof, 66)
        self.assertIsNone(result.signal_result_from_slha(slha.read(""), "LILITH"))


class ObservablesTest(unittest.TestCase):
    def test_complex_parts_are_combined(self):
        obs = fs.Observables(
            ["amm_Fe_1", "Re(c(0))", "Im(c(0))", "Re(c(1))", "Im(c(1))", "Re(z)", "Im(z)"],
            [1e-9, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0], [])
        self.assertTrue(obs.ok)
        self.assertEqual(obs["amm_Fe_1"], 1e-9)
        self.assertEqual(list(obs["c"]), [1 + 2j, 3 + 4j])
        self.assertEqual(obs["z"], 5 + 6j)
        self.assertEqual(sorted(obs), ["amm_Fe_1", "c", "z"])


class ArraysTest(unittest.TestCase):
    def test_to_entries(self):
        self.assertEqual(arrays.to_entries(2.0, (), "x"), {(): 2.0})
        self.assertEqual(arrays.to_entries([[1, 2]], (1, 2), "x"), {(0, 0): 1, (0, 1): 2})
        with self.assertRaises(ValueError):
            arrays.to_entries([1, 2, 3], (2,), "x")
        with self.assertRaises(ValueError):
            arrays.to_entries("a", (), "x")

    @unittest.skipUnless(arrays.HAVE_NUMPY, "numpy not installed")
    def test_numpy(self):
        import numpy as np
        self.assertEqual(arrays.to_entries(np.eye(2), (2, 2), "x")[(1, 1)], 1.0)
        a = arrays.from_entries((2,), {(0,): 1.0, (1,): 1j})
        self.assertEqual(a.dtype, complex)
        self.assertEqual(arrays.to_list(a), [[1.0, 0.0], [0.0, 1.0]])


class FindRootTest(unittest.TestCase):
    def test_environment_variable(self):
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / "models").mkdir()
            old = os.environ.get("FLEXIBLESUSY_DIR")
            os.environ["FLEXIBLESUSY_DIR"] = tmp
            try:
                self.assertEqual(fs.find_root(), Path(tmp).resolve())
            finally:
                if old is None:
                    del os.environ["FLEXIBLESUSY_DIR"]
                else:
                    os.environ["FLEXIBLESUSY_DIR"] = old

    def test_not_a_root(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(fs.BackendError, "not a FlexibleSUSY directory"):
                fs.find_root(tmp)


if __name__ == "__main__":
    unittest.main()
