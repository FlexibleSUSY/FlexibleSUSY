# Python interface to FlexibleSUSY

`flexiblesusy` runs spectrum generators created with FlexibleSUSY from
Python 3 (3.8 or newer).  It needs only the standard library; with numpy
installed, masses and matrices are returned as numpy arrays.

## Use from a FlexibleSUSY checkout

Generate and build a model as usual, e.g.

    ./createmodel --name=MSSM
    ./configure --with-models=MSSM
    make

then

    PYTHONPATH=<FlexibleSUSY>/python python3

```python
import flexiblesusy as fs

mssm = fs.Model("MSSM")                       # finds <FlexibleSUSY>/models/MSSM
point = mssm.default_point()                  # models/MSSM/LesHouches.in.MSSM
point = mssm.updated(point, TanBeta=20,
                     settings=fs.Settings(pole_mass_loop_order=2))
res = mssm.calculate(point)

res.ok, res.problems, res.warnings
res.pole_masses["hh"]                         # 0-based: [m_h1, m_h2]
res.mixings["ZH"]
res.parameters.scale, res.parameters["Yu"]    # running parameters
res.to_dict()                                 # JSON-serializable
```

With the compiled backend (see below), observables and decays are calculated
on demand, from the spectrum kept in memory, without calculating it again:

```python
res = mssm.calculate(point)                   # spectrum only
if res.ok and res.pole_masses["hh"][0] > 110:
    obs = res.observables()                   # e.g. obs["amm_Fe_1"], obs.problems
    dec = res.decays()                        # e.g. dec["hh", 0].channels, dec.by_pdg(25)
```

Both are cached in the result.  Decays need Collier or LoopTools
(`Settings(loop_library=fs.LoopLibrary.COLLIER)`); the loop library is fixed
by the first calculation in a Python process.  `res.release()` frees the
model kept in memory.

HiggsBounds/HiggsSignals (HiggsTools) and Lilith work the same way, if
FlexibleSUSY was configured with them (`--enable-higgstools`, `--with-lilith`):

```python
mssm = fs.Model("MSSM", higgsbounds_dataset=".../HBDataSet", higgssignals_dataset=".../HSDataSet",
                lilith_database=".../Lilith/data/latest.list")
res = mssm.calculate(point)
ht = res.higgstools()        # ht.higgssignals.chi2, ht.higgsbounds[0].obs_ratio, ...
ll = res.lilith()            # None if no Higgs boson has a mass between 123 and 128 GeV
```

Both use the normalized effective couplings, so the compiled backend calculates
the decays first if needed.  With the SLHA backend they must be enabled before
the calculation (`DecaySettings(calculate_decays=True, call_higgstools=True,
call_lilith=True)`), and the paths are passed to `run_<M>.x`.

A point can also be created from scratch.  All input parameters of the
model are required; settings and SM inputs that aren't given use the
FlexibleSUSY defaults:

```python
point = mssm.point(
    inputs={"TanBeta": 10, "SignMu": 1, "Qin": 1.9422e16, ...},
    settings=fs.Settings(loop_library=fs.LoopLibrary.COLLIER),
    sm_inputs=fs.SMInputs(alpha_s_mz=0.1181),
    decay_settings=fs.DecaySettings(calculate_decays=True))
```

`mssm.input_parameters` lists the input parameters with their shapes.

## Install into a virtual environment (optional)

    pip install -e <FlexibleSUSY>/python

An editable install keeps pointing at the checkout, so the models are still
found.  Otherwise pass `root=` to `fs.Model` or set `FLEXIBLESUSY_DIR`.

## Errors

- `fs.InputError`: unknown or missing input parameter, wrong shape or type.
- Problems of the calculation (no convergence, tachyons, ...) don't raise:
  `res.ok` is False and `res.problems` lists them.  Use
  `mssm.calculate(point, raise_on_problems=True)` to get `fs.CalculationError`.
- `fs.BackendError`: the model isn't built, its metadata is missing, or the
  spectrum generator crashed.

## Backends

`fs.Model(name, backend=...)`:

- `compiled`: the model's Python module `models/<M>/flexiblesusy_<M>*.so`,
  built when FlexibleSUSY is configured with `--with-python=<interpreter>`
  (needs pybind11 for that interpreter, e.g. `--with-python=python3.14`).
  Calculations run in the Python process, and observables and decays are
  calculated on demand.  Calculations must not run in parallel threads of one
  process (FlexibleSUSY has global state); use processes.
- `slha`: writes an SLHA input file, runs `models/<M>/run_<M>.x` and reads
  its SLHA output.  Decays are only available if enabled before the
  calculation (`DecaySettings(calculate_decays=True)`), and observables not
  at all.
- `auto` (default): `compiled` if the module is built for the running
  interpreter, `slha` otherwise.

Both return the same `Result`: masses, mixings and parameters come from the
SLHA output in both cases (9 significant digits).

The model-specific information (input parameters and their SLHA blocks,
particles, output blocks) comes from `models/<M>/<M>_metadata.json`, which the
code generation writes.

## Tests

    cd <FlexibleSUSY>/python && python3 -m unittest discover -s tests

The MSSM tests are skipped unless `models/MSSM` is generated and built; the
HiggsTools and Lilith tests also need `HIGGSBOUNDS_DATASET`,
`HIGGSSIGNALS_DATASET` and `LILITH_DATABASE` in the environment.
