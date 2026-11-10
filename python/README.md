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
res.decays["hh", 0].channels                  # needs decays enabled, see below
res.to_dict()                                 # JSON-serializable
```

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

- `slha` (current): writes an SLHA input file, runs `models/<M>/run_<M>.x`
  and reads its SLHA output.
- A compiled backend with the same API is planned.

The model-specific information (input parameters and their SLHA blocks,
particles, output blocks) comes from `models/<M>/<M>_metadata.json`, which the
code generation writes.

## Tests

    cd <FlexibleSUSY>/python && python3 -m unittest discover -s tests

The MSSM tests are skipped unless `models/MSSM` is generated and built.
