# FlexibleSUSY: notes for Claude

FlexibleSUSY generates C++ spectrum generators for BSM models from SARAH model
files. The Mathematica meta code reads SARAH's output and fills C++ templates, which
gives a model directory with its own library and executable. FlexibleDecay, the decay
module, is part of it.

## Layout

| Path | Contents |
|---|---|
| `meta/` | Mathematica meta code (`FlexibleSUSY.m` drives it; `Decays.m`, `CXXDiagrams.m`, `Vertices.m`, `SelfEnergies.m`, `Unitarity.m`, ...) |
| `templates/` | C++ templates with `@placeholders@`, filled by the meta code |
| `src/` | Model-independent C++ library (`libflexisusy`), including `src/decays/` and `src/cxx_qft/` |
| `model_files/<M>/` | Shipped FlexibleSUSY model files (`FlexibleSUSY.m.in`) and SLHA inputs (`LesHouches.in.<M>`) |
| `sarah/` | SARAH model files shipped with FlexibleSUSY (others come from SARAH's `Models/`) |
| `Output/` | SARAH's cached output for each model and eigenstate set |
| `model_specific/` | Higher-order corrections from the literature, shared by several models (SM, MSSM, NMSSM, ...) |
| `models/<M>/` | A created model: inputs `FlexibleSUSY.m` and `start.m`, and the **generated** C++ code |
| `test/` | Tests: Boost.Test `.cpp`, Mathematica `.m`, shell `.sh`; registered in `test/module.mk` |
| `utils/` | Reusable scripts for developers and users (SLHA processing, scans, static analysis, ...); new reusable scripts go here |
| `utils/loop_decays/` | FeynArts/FormCalc scripts that produced the generic 1-loop decay diagrams (`make` regenerates `meta/generic_loop_decay_diagram_classes.m` and `src/decays/one_loop_decay_diagrams.*`; it needs FormCalc <= 9.6, because it relies on FormCalc's internal representation of the amplitudes, which changed after 9.6) |
| `doc/` | Documentation (see "Further documentation" below) |

**Never edit the generated C++ files in `models/<M>/` to fix a bug.** They are
overwritten when the code is regenerated. Fix the meta code, the template or `src/`,
then regenerate. Only `models/<M>/FlexibleSUSY.m` (copied from `model_files/` by
`createmodel`) and `start.m` are meant to be edited, to change a model's setup.

## Creating, generating and building a model

```sh
./createmodel --name=THDMII            # copies model_files/THDMII to models/THDMII
./configure --with-models=THDMII,MSSM  # every model to build must be listed
make -j4                               # runs the meta code, then compiles
```

- **Running the meta code:** it runs when `models/<M>/00_DELETE_ME_TO_RERUN_METACODE` is
  missing or older than its inputs (meta code, templates, SARAH files). To rerun it for
  one model: `rm models/<M>/00_DELETE_ME_TO_RERUN_METACODE && make models/<M>/00_DELETE_ME_TO_RERUN_METACODE`.
  Because of this, a generated file that you patched by hand for an experiment is
  overwritten by the next `make` after any change to `meta/` or `templates/` (even a
  comment). Let the code generation run first, then patch.
- **Run times:** code generation takes minutes per model, and longer for large models.
  The slowest steps are usually the vertices, decays, `BrLTo3L` and unitarity.
- **`make -n` isn't a safe dry run:** it can still trigger the meta code, because the
  generated `module.mk` files are make prerequisites. Avoid it in a tree where the
  generated code matters.
- **Progress counter (`[step/total] Building ...`):** `MSG` in `config/Makefile.in` counts
  the steps; the total comes from a dry run (`make -n FS_PROGRESS_COUNT=1 ENABLE_META=no`,
  which doesn't remake dependency files) started when the first counted step runs.
  Rules that can run while make remakes its makefiles (code generation, generated test
  sources) use `MSG_NO_PROGRESS`, as make restarts and resets the counter afterwards.
  The total is only right if `make -n` predicts the real build: a phony prerequisite of a
  file, or a recipe that doesn't update its target (Fortran `.mod` files), makes `make -n`
  count steps that the real build skips. To compare, run
  `make -n FS_PROGRESS_COUNT=1 ENABLE_META=no 'MSG=echo STEP $@'` before `make`.
- **One target generates all models:** for the same reason, `make` of any target
  (even one model's `00_DELETE_ME_TO_RERUN_METACODE`) first generates every configured
  model whose code is missing or out of date. Without `-j` they run one after another;
  with `-j` several code generations can run at once.
- **Models sharing a SARAH model must not be generated in parallel:** e.g. `MSSM` and a
  copy of it with different flags both read and write SARAH's cache in `Output/<SARAH
  model>/`. Generated at the same time (`make -j`), they can leave the cache broken, and
  every later run reuses it ("Found up-to-date vertex file"): errors such as `Drop::drop`
  and `Part::partw` at "Converting self energies", and a model class without some
  self-energies (e.g. `self_energy_VZ_1loop`), which doesn't compile. Generate such models
  with `make -j1` (then compile with `-j`); to recover, move `Output/<SARAH model>` away
  and regenerate. Parallel generations also compete for the licensed subkernels, which can
  make one of them very slow.
- **Reconfiguring:** the first line of `config.status` holds the previous options, but
  without quotes, so an option with spaces (`--with-cxxflags="-O1 -g"`) must be quoted
  again by hand.
- **Sanitizers:** compile and link with them through `CXX`, as not every link rule uses
  `CXXFLAGS`: `make CXX="clang++ -fsanitize=address,undefined -fno-omit-frame-pointer"`,
  with `--with-cxxflags="-std=c++20 -O1 -g -fPIC -D_GLIBCXX_ASSERTIONS"`. The Fortran loop
  libraries aren't instrumented. Known false positives: `test_gsl_vector` allocates
  `SIZE_MAX` elements on purpose (use `ASAN_OPTIONS=allocator_may_return_null=1`), and
  LoopTools' `clooptools.h` addresses its cache relative to a common-block array, which
  UBSan reports as out of bounds (suppress `bounds:clooptools.h`). Tests with tolerances
  near 1e-16 can fail at `-O1` by rounding.
- **Memory:** large generated files such as `models/MRSSM2/decays/MRSSM2_decays.cpp`
  (about 185k lines) need several GB each to compile, so choose `-j` with the available
  memory in mind.
- **Mathematica kernels:** depending on the licence, the number of concurrent kernels
  can be limited. When it's exhausted, code generation or a Mathematica test can fail
  with "Invalid password" or SARAH not loading. The meta code also starts parallel
  subkernels (`FSEnableParallelism`) with `LaunchKernels[]`. Licences limit these
  separately (`$MaxLicenseSubprocesses`), shared by all running kernels. When none is
  free, `LaunchKernels[]` prints a licence error ("No valid password found") for each
  kernel it couldn't start, and `ParallelMap` continues sequentially
  (`ParallelMap::nopar`): the result is correct, only slower.
  A main kernel that gets no licence prints "cannot find a valid password ... Please press
  Enter to quit" and exits with a non-zero status (63). With `make -j`, the Mathematica-generated
  test sources (`test/test_<M>_*.meta`) start together while make remakes its makefiles; a failure
  there is ignored (the `.d` files are `-include`d), and make generates the missing sources
  again after it restarts, so they can appear several times in the log.
- **Generating and compiling separately:** `./configure --disable-compile` only generates
  the code. `./configure --disable-meta` only compiles code that was generated before,
  without Mathematica or SARAH.
- **Stale objects in static libraries (older configurations):** until #769, `configure` set
  `MODULE_MAKE_LIB_CMD` to `ar cruU`, which replaces a member only if the object file is newer than
  it, to the second; an object rebuilt within the same second (e.g. after swapping a generated file
  back, with ccache) stayed outdated in the library. Build trees configured before that still use it
  (see `Makefile`); reconfigure them, or delete the library before rebuilding. If results look
  inconsistent, compare the archived object with `ar x` + `md5sum`.
- **Moving generated code:** `make pack-<M>-src` packs one model's generated code. With
  `--with-install-dir=<dir>`, `make install-src` exports a standalone package without the
  meta code. It only exports files listed in the `module.mk` files (headers and `.inc`
  files too), so list new files there. Files in subdirectories are installed with
  `$(INSTALL_FILES)` (`config/install_files.sh`), which keeps their relative paths. CI runs
  `install-src` without an install directory, where it does nothing, so test an export
  by configuring and building it.
- **Faster development builds:** `--enable-shared-libs` with
  `--with-cxxflags="-std=c++20 -O0 ..."` shortens compile times. `make CXX="ccache clang++"`
  avoids recompiling unchanged files after the meta code rewrites them. With shared
  libraries, executables link them by absolute path, so after moving the tree run
  `make clean-executables && make allexec`.
- **Each model gets** `models/<M>/run_<M>.x`, the command-line spectrum generator, and
  `models/<M>/run_<M>.m` with the LibraryLink library, the Mathematica interface.

### Useful configure options
- **Tests:** `--with-optional-modules=test` enables them (not `--with-addons`).
- **Loop libraries:** `--with-loop-libraries=collier,looptools`, together with
  `--with-collier-incdir/libdir` and `--with-looptools-incdir/libdir`.
- **Compiler:** `--with-cxx=clang++` or `g++`.
- **External tools:** `--enable-higgstools` with `--with-higgstools-incdir/libdir`, and
  `--with-lilith=<dir>`.
- **Boost libraries:** `--with-boost-libdir=...` if `configure` doesn't find
  `boost_unit_test_framework`; otherwise tests compile but don't link.
- **fflite:** it can't be configured together with LoopTools (`configure` refuses). The
  loop-function cross-check (`test_pv_crosschecks.sh`) builds its own fflite, so it
  works in a Collier/LoopTools build.

### Compiler requirements
The code uses C++20 (concepts). It needs g++ >= 10.1, or clang++ >= 10 with
libstdc++ >= 10; on macOS Apple clang with libc++ (CI job `no-meta-macos`, #748). `configure` checks the compiler version and, with a test program, that
the standard library has `<concepts>`, `<span>` and `<numbers>`
(`check_cxx_std_library`). The classic Intel compiler icpc isn't supported (only a
warning); Intel's icpx is detected as clang. Libraries: Boost >= 1.62 except 1.64
(Hana, Boost.Math with C++20; 1.64's odeint doesn't compile), Eigen >= 3.3
(`Eigen::Index`), GSL >= 2.0, all checked by `configure`.

Before using a newer standard library feature, check it against these minimum
versions:
- **Usable:** `std::span`, `std::ssize`, `std::string::starts_with`, `map::contains`,
  `std::erase_if`, `std::remove_cvref_t`, `<numbers>`, concepts and `requires`.
  `std::span` was tested with g++ 10 and clang 10-12 on libstdc++ 10.
- **Not usable yet:** `<ranges>` (fails with clang 10 on libstdc++ 10), and features
  that need g++ >= 11 or later: `std::format` (13), `std::to_chars` for `double`,
  `std::bit_cast`, `std::source_location`.

**Raising the requirements:** when a change (e.g. replacing Boost by a std feature) needs
a newer compiler or standard library, update both places that state the requirement in
the same change: the minimum versions in `configure` (`required_*_compiler_version`, or a
feature test with `try_compile_cpp_program` for things the version number doesn't
show, such as the standard library) and the "Requirements" list in `README.rst`. Also
update this section and mention the new requirement in `CHANGES.rst`.

## Running a spectrum generator

```sh
models/<M>/run_<M>.x --slha-input-file=models/<M>/LesHouches.in.<M> --slha-output-file=out.spc
```

- **Input blocks:** `Block FlexibleSUSY` configures the calculation (precision, loop
  orders, solver, ...); `Block FlexibleSUSYInput` holds extra physical inputs; and
  `Block FlexibleDecay` configures the decays. The entries are documented in
  `doc/slha_input.rst`, those of `FlexibleDecay` in
  `doc/observables/FlexibleDecay.rst`. The EWSB loop order (`FlexibleSUSY[5]`) should
  equal the pole mass loop order (`FlexibleSUSY[4]`).
- **SLHA output follows the SLHA and SLHA2 standards:** these fix which parameters are
  written, in which blocks and in which order. Don't propose changing the standard
  output to suit another program. Quantities outside the standard, such as SARAH's
  `ANGLES 10` (ThetaWp) read by CalcHEP/micrOMEGAs models, belong in the model file's
  `ExtraSLHAOutputBlocks`. Differences in how another program interprets standard
  entries are that program's convention, to be documented in the comparison.
- **Loop library:** `FlexibleSUSY[31]` in the input selects it (0 softsusy, 1 Collier,
  2 LoopTools, 3 fflite). The environment variable `FLEXIBLESUSY_LOOP_LIBRARY` is only
  read when this entry is -1, which is the default when it's missing
  (`src/loop_libraries/loop_library.cpp`). An explicit `31 0`, as in the generated
  inputs, wins over the environment variable.
- **Decays are only calculated with Collier or LoopTools.** With softsusy or fflite they
  are skipped, with only a message on the terminal.
- **Block `FlexibleDecay` (`src/decays/flexibledecay_settings.cpp`):**

  | Entry | Setting |
  |---|---|
  | 0 | calculate decays |
  | 1 | minimum BR to print (set 0 to see all channels) |
  | 2 | higher-order corrections (digits SM\|BSM, default 41) |
  | 7 | HiggsTools |
  | 8 | Lilith |
  | 9 | pole Higgs mixings (default 1) |

## Pole masses and tadpoles

- **Tadpoles enter only through EWSB:** the VEVs vd, vu are fixed and are the minimum of the
  loop-corrected potential, so no tadpole diagrams are added to self-energies. The loop tadpoles
  enter through the `EWSBOutputParameters` (soft Higgs masses, or μ, κ, vS, ... depending on the
  model), solved at `FlexibleSUSY[5]` loops.
- **Which parameters go where:** `get_mass_matrix_X()` in the pole-mass equation uses the
  loop-level EWSB solution. Masses and mixings inside the loops come from
  `calculate_DRbar_masses()`, which re-solves the EWSB at tree level and restores the
  loop-level values afterwards. W, Z, gluino and top start from those stored running masses.
  This is consistent at one loop; other choices differ at two loops.
- **Small VEVs:** if an EWSB output multiplies a small VEV (e.g. the singlet and triplet of the
  MRSSM, vS, vT ~ 1 GeV), the tadpoles enter as t/v and the two-loop differences between codes
  that solve for different parameters (SPheno: mS2, mT2; FS: vS, vT) can be large. Compare at
  the same Lagrangian parameters.

## Tests

- **Building a test:** `make test/test_<name>.x`, then run `./test/test_<name>.x`.
  Running it directly shows the Boost output.
- **Running test sets:** `make alltest` only builds the tests. `make execute-tests` runs
  all of them, and `make execute-compiled-tests`, `make execute-meta-tests` and
  `make execute-shell-tests` run one kind. They write `test/*.xml` and
  `test/*.log`, summarized in `test/test.xml`. `make -B <target>` forces a rebuild, and
  `make clean-test-log` or `make clean-test` clean up.
- **A single Mathematica test:** `make test/test_<name>.m.xml`. It needs the `meta/`
  directory on the Mathematica path, which the rule adds. A meta test loads SARAH,
  FlexibleSUSY, `TestSuite` and the package under test, and ends with
  `PrintTestSummary[]`.
- **Registering a test:** add it to `TEST_SRC`, `TEST_META` or `TEST_SH` in
  `test/module.mk`. Model tests are inside `ifeq ($(WITH_<M>),yes)` blocks, and a model
  test needs a rule `$(DIR)/test_<name>.x: $(LIB<M>)`.
- **FlexibleDecay tests** (`test_<M>_FlexibleDecay`) are only added when
  `FLEXIBLESUSY_LOOP_LIBRARY` is 1 or 2 at **make** time, and need it at run time too:
  `FLEXIBLESUSY_LOOP_LIBRARY=1 make test/test_SM_FlexibleDecay.x`.
- **Tolerances:** many tests compare with stored reference values at tight tolerances
  (down to 1e-16). A change of the last digit can fail them without a real error.
- **CI:** the `with-meta` jobs in `.github/workflows/tests.yml` activate a Wolfram Engine
  first. If they fail at "Activate Wolfram Engine", no test ran.

## The decay module (FlexibleDecay)

FlexibleDecay ([2106.05038]) calculates 2-body decays of scalars, mainly Higgs bosons,
at **leading order**:
- **Tree level or 1-loop:** if a decay exists at tree level, that amplitude is used;
  otherwise the 1-loop amplitude is generated, so loop-induced decays such as
  h -> gamma gamma include all 1-loop BSM contributions.
- **Higher-order corrections:** known SM QED, QCD and electroweak corrections are added
  where possible.
- **Scheme:** the decays use a modified MS-bar scheme in which heavy BSM physics
  decouples, while BSM parameters stay MS-bar or DR-bar.

The code:
- **Generated part (`meta/Decays.m`):** builds, for each decaying particle and final
  state, the tree-level amplitude or, if there is none, the 1-loop amplitude. 1-loop
  diagrams are matched to the generic diagrams in
  `meta/generic_loop_decay_diagram_classes.m`, generated with FeynArts/FormCalc.
  - Topologies are named T1 to T10 as in FeynArts.
  - Each insertion is multiplied by `S_gen/S_class` (FeynArts' `RelativeCF`), as the
    generic expressions contain `1/S_gen`.
- **Fixed part (`src/decays/`):**
  - `decay_amplitudes.cpp`: squared amplitudes per Lorentz structure (`SSS`, `SSV`,
    `SVV`, `SFF`, `FFS`, `FFV`); `SVV` has a Ward identity check. They are spin-averaged
    over the decaying particle. References to check them against: `SquareAmp_*` in the
    SPheno code SARAH writes (`<SPheno output>/LoopDecays/Decay*.f90`; SARAH's overall
    factor and sign of momentum structures can differ), and SPheno's `DecayFunctions.f90`.
  - `one_loop_decay_diagrams.cpp`: generated generic diagrams.
  - `decay.cpp`: decay tables and effective couplings.
  - `experimental_constraints.cpp`: HiggsTools and Lilith.
  - `specializations/`: hand-written channels such as H -> qq, gg, gamma gamma.
- **Couplings and masses:** decays use the `*_mass_eigenstates_decoupling_scheme` class
  (`flag = 1`). It solves EWSB at tree level and, with FlexibleDecay[9] = 1, replaces
  the Higgs mixing matrices by pole mixing matrices.
- **Width form factors:** the width of S -> V V with massless vectors uses
  `form_factor_21` (and `form_factor_eps`). `form_factor_g` only enters through the
  Ward check.
- **PDG codes:** each physical field has a unique PDG code; 0 marks unphysical fields (Goldstone
  bosons, ghosts), which never appear in decay tables. Decay channels are identified by their PDG
  codes (`Decays_list::set_decay`), so (i, j) and (j, i) of identical fields replace each other.
- **Decays of fermions and vectors** are tree level only (`GetMaxLoopOrderForDecay`).
  Vector decays use `Decay_amplitude_{VSS,VFF,VVS,VVV}`; their `square()` matches
  SPheno's `VectorBosonTo*` routines in `DecayFunctions.f90`. They aren't advertised yet,
  but can be enabled for testing (e.g. `VZp`, `Fv` in `FSDecayParticles`): in BLSM,
  Z' -> f f and heavy neutrino N -> W l, Z nu, h nu agree with SARAH's SPheno to a few
  percent at equal masses (Γ_Z'/M_Z' to 1%).
- **Decay scale:** each particle's decays are calculated after `model.run_to(mass)`
  (`run_to_decay_particle_scale` in `templates/decays/decays.hpp.in`). Parameters fixed
  only at the input scale, such as the triangular kinetic-mixing basis (gX1 = 0) of
  U(1)xU(1) models, no longer hold there. Vertices written with mixing angles
  (`ThetaWp`) then break, e.g. Z' -> W W. Check the equivalence theorem
  (Γ(Z' -> WW) ≈ Γ(Z' -> Zh)) when changing Z' decays.
- **Masses of external vectors:** final-state vectors are on-shell, so `SVV`, `VVS` and
  `VVV` use their physical (1-loop) masses on purpose, while the couplings come from the
  model at the decay scale. With heavy decaying particles this breaks the equivalence
  theorem at the level of (m_V,run/m_V,phys)^4: in BLSM with M_Z' = 3.7 TeV,
  Γ(Z' -> WW) is 0.80 of micrOMEGAs's value, while Z' -> f f agree to 1e-4. This is a
  deliberate choice, not a bug; the deviation is of higher order.

## Unitarity constraints

`meta/Unitarity.m` builds the full 2 -> 2 scalar scattering amplitudes with SARAH
(`InitUnitarity`, `GetScatteringDiagrams`: s-, t-, u-channel and contact diagrams) and then takes
s -> infinity, where only the quartic contact term survives. This is a deliberate stepping stone
towards the finite-s calculation of [1805.07306]: changes (also speed-ups) must keep the full
amplitudes and the trilinear vertices, and only touch the s -> infinity step separately.

## Conventions

- **Branches:** development happens on `development`; open pull requests against it.
- **`CHANGES.rst`:** bug fixes cite the commit as `[commit 1234abc]`. These hashes
  change if a branch is rebased. Pull requests are often squash-merged, which also
  changes them. The changelog is reviewed and the hashes fixed when a release is
  prepared, so don't fix them one at a time on `development`.
- **`CHANGES.rst` is for end users:** describe what was wrong from the user's point of
  view (which models or settings, which observables, how large the effect) and what
  changes. The technical cause of a bug belongs in the commit message, not in the
  changelog.
- **Tests with fixes:** add a test with a bug fix where possible, as a Boost.Test
  `.cpp` registered in `test/module.mk`.
- **C++ style:** C++ code follows `.clang-format`
  (`clang-format -style=file -i <file>` from the root). Don't run it on templates:
  it breaks the `@placeholders@`.
  `git clang-format <base>` refuses to touch files with unstaged changes ("Please commit, stage,
  or stash them first"); use `git clang-format --diff <base> -- src test | git apply` instead. It also
  treats `.m` files as Objective-C, so restrict it to C++ paths.
- **clang-tidy:** `.clang-tidy` enables all checks (`*`) minus a list. Run it only on
  changed lines (`git diff -U0 <base> | clang-tidy-diff -p1 ...` with `-std=c++20` and
  the include paths of the build), as existing code has many findings. Newer clang-tidy
  versions add checks that conflict with the code base: unchecked container access
  (`v[i]`, Eigen `m(i,j)`), identifier length (physics names such as `m`, `g2`), math
  parentheses, linkage of Boost.Test cases, and `*-use-ranges` (`<ranges>` isn't usable
  with the minimum compilers). Worth fixing: `bugprone-*` and missing direct includes
  (`misc-include-cleaner`; ignore Boost.Test and Eigen, which are used through umbrella
  headers).
- **Mathematica style:** functions should fail loudly on wrong input, because garbage
  ends up in the generated C++. So:
  - restrict arguments with patterns, such as `i_Integer /; i > 0`
  - report errors with `Utils`AssertOrQuit[...]` and a message
  - add a catch-all with `Utils`MakeUnknownInputDefinition`
  - warn with `Utils`FSFancyWarning`
- **Checking Mathematica code:** `CodeInspector`CodeInspect[File["meta/X.m"]]` finds
  syntax errors and slips such as a missing `;` at the end of a line (which multiplies
  two statements), string juxtaposition (`"a" "b"` is a product, not a concatenation) or
  `x_NumericQ` instead of `x_?NumericQ`. On the existing meta code it reports hundreds of
  low-confidence findings (`PatternRule`, `ImplicitTimesPseudoCall`), so look at
  Fatal/Error findings with `ConfidenceLevel >= 0.9` on the changed lines.

## Pitfalls

- **SARAH:** `conj` isn't mapped over `sum`: `conj[sum[i, 1, 2, a[i]]]` stays
  `sum[i, 1, 2, a[i]]`.
- **Global symmetries in SARAH models:** `Global` lists them (e.g. `{Z[2], RParity}` in the MSSM,
  `{U[1], RSymmetry}` in the MRSSM), with charges per superfield component; `SA`ChargeGlobal[field,
  sym]` gives them for mass eigenstates too. Z[N] charges are phases (`-1`, `Exp[2 Pi I/3]`), and
  `SA`ChargeGlobal[conj[f], sym]` returns the charge of `f` for every Z[N] (only right for Z[2]).
  A model may declare a symmetry that it breaks on purpose: the MRSSM's R-symmetry holds physically,
  but the model file has the R-violating term mu Hu Hd, because SARAH generates the B-term B_mu only
  with it, and mu is set to 0 later (hence `Superpotential::ViolationGlobal` when loading it). Code
  generation can't assume such a parameter is 0, so such a symmetry counts as broken:
  `Vertices`GetExactGlobalSymmetries[]` only returns symmetries that no superpotential or Lagrangian
  term violates.
- **Contexts of model file flags:** the model file is read with `Get` in the `Global``
  context, so an unqualified name like `UseHiggs2LoopMSSM` is set in the first context on
  `$ContextPath` that has the symbol (`FlexibleSUSY`` comes before `SARAH``). Declare every flag
  with its default in the public part of `FlexibleSUSY.m`, and refer to it with the same context
  everywhere in the meta code. Otherwise one explicit `FlexibleSUSY`Flag` (or `SARAH`Flag`)
  elsewhere silently changes which symbol the model file sets.
- **`For` loop bounds that don't evaluate:** `For[l = 1, l <= f[x], l++, ...]` runs zero times if
  `f[x]` stays unevaluated (e.g. a function with no matching definition), without any message.
- **Conditions inside list patterns:** a condition on one element can't use the pattern variable of
  another element: `f[{a_ /; a >= 2, b_ /; b >= a}]` never matches, even `f[{2, 3}]`. Put the
  condition on the whole list: `f[{a_Integer, b_Integer} /; 2 <= a <= b]`. Test such definitions
  with valid input too, not only with input that should be rejected.
- **CKM matrix:** the built-in SM (`model_specific/SM/standard_model.cpp`) has a unit CKM matrix:
  `initialise_from_input(qedqcd)` puts the running masses on the diagonal of the Yukawa matrices and
  ignores the CKM matrix of `qedqcd`. The same holds for every model with `{Yu, Automatic}` (SM, MSSM,
  THDM, ...). Only models whose low-scale input uses it (`{Yu, Sqrt[2] MatMul[upQuarksDRbar, CKM]/vu}`:
  CMSSMCKM, MRSSM2, MRSSM2CKM, S3) have quark flavour mixing. Masses and the Higgs-mass matching don't
  depend on the CKM matrix (it doesn't change the singular values of the Yukawa matrices), but never take
  the CKM matrix of a BSM model from the built-in SM: use `qedqcd.get_complex_ckm()`, as the low-scale
  constraint does (the decays in the decoupling scheme don't yet: PR #776).
- **Quark phases of mass eigenstates:** the phases of the model's quark mass eigenstates come from
  the diagonalization and are arbitrary, so a flavour-changing amplitude (b -> s) can't be normalized
  with the CKM matrix of `qedqcd` (PDG convention): in CMSSMCKM V_tb V_ts^* from `ZUL ZDL^dagger` is
  +0.041, in MRSSM2CKM -0.041 (PDG: -0.042 + 0.0008 i). Rephase (PR 633 adds `ckm_rephasing_phase()`
  in `src/ckm.hpp`) or normalize with the model's CKM.
- **Soft scalar masses aren't all hermitian:** `SARAH`ListSoftBreakingScalarMasses` also contains the
  soft masses that mix two different fields with the same quantum numbers (SARAH's `SA`ListM2ijOFF`,
  e.g. `mlHd2` in the RPV MSSM, the square 3x3 `meEX2` and `muUX2` in U1xMSSM3G). Only the others,
  m^2_{ij} phi_i^* phi_j of one field, are hermitian. SARAH nevertheless sets `Adj[meEX2] = meEX2`.
- **RG running and rounding noise:** the RK stepper (`src/rk.hpp`) controls the relative error of
  each parameter separately (yscal = |y| + |h dy/dx| + 1e-16). A parameter that is analytically zero
  but filled with rounding noise can't be integrated: it runs out of steps, reported as "RG
  integration failed (too many steps) ..., step size limited by <par>" (#786; before, "non-perturbative
  <par>"). Keep such parameters exactly zero rather than nearly zero: the beta functions of hermitian
  soft masses are made exactly hermitian for this reason (#785; before, most CP-violating points failed).
  Don't fix this with an absolute floor relative to the whole parameter vector: it mixes mass
  dimensions, so GUT- or seesaw-scale masses would ruin the precision of small couplings; a floor would
  have to be per parameter matrix.
- **SARAH's vertex lists (`SARAH`VertexList3/4`):** each entry equals `SARAH`Vertex[]`
  **without** dependences, but SARAH keeps a vertex only if it doesn't vanish **with** the
  dependences (`subDependences`, e.g. `g1 -> e Sec[ThetaW]`), so h-γ-γ, χ-χ-γ etc. are missing
  although `SARAH`Vertex[]` without dependences gives them non-zero. A missing vertex is zero
  with dependences (checked for all SM vertices). `Vertices`FSVertex` and
  `Vertices`IsNonZeroVertex` (used by the decays, always with dependences) are built on this.
  `SARAH`Vertex[..., UseDependences -> True]` is slow for large vertices (seconds for 4-squark
  vertices).
- **Eigen comma initialization:** it assigns in order, so `a << a[0] + a[1], a[0];` uses
  the new `a[0]` for the second element.
- **Shell scripts:** prefer `bash` for scripts with variable expansion. Some interactive
  shells, such as zsh, parse `$VAR:path` as a modifier.

## Cross-checks with SARAH and SPheno

FlexibleSUSY and SARAH/SPheno start from the same SARAH model, so they should give
similar results. Every difference in a predicted observable (masses, mixings, decays,
couplings) needs a physics explanation, such as:
- a different perturbative order or set of higher-order corrections
- a different renormalization scheme or scale
- different boundary conditions or input conventions
- a different setup of the calculation, such as the decoupling scheme in FlexibleDecay

A difference without such an explanation is probably a bug, either in FlexibleSUSY or in
SARAH/SPheno. When in doubt about a result or a model definition:
- **Model level:** query SARAH directly (`Start[...]`, `MassMatrix[...]`, `Vertex[...]`,
  `MakeVertexList[EWSB]`). To leave the checkout's `Output/` untouched, copy the model
  to a scratch directory and set `SARAH[InputDirectories]` and `SARAH[OutputDirectory]`
  as in `templates/start.m.in`. `Start` alone doesn't compute vertices: messages such as
  `Vertex::ChargeViolating` only appear with `MakeVertexList` (or `ModelOutput`).
- **Full calculation:** compare with a SPheno spectrum generated by SARAH for the same
  model and inputs.
- **Tree-level decays, exact:** CalcHEP (micrOMEGAs) with SARAH's CalcHEP model reads masses,
  mixings, Yukawas and trilinears from FlexibleSUSY's SLHA output, but computes `el` from
  α(0), M_W from (α, G_F, M_Z), g1, g2, v from those and g3 from α_s(Q). Read g1, g2, g3
  (Block GAUGE) and v (HMIX 3) in its `func1.mdl` instead, run FlexibleSUSY with tree-level
  masses (`FlexibleSUSY[4] = [5] = 0`) and LO decays (`FlexibleDecay[2] = 0`), and switch off
  the decoupling scheme (`flag = 2`) and the running to the decay scale
  (`run_to_decay_particle_scale`) in the generated decays: then all tree-level widths agree to
  CalcHEP's printing precision (MSSM: 282 channels within 7e-7). Call CalcHEP's
  `slhaDecayPrint` directly instead of micrOMEGAs' `main` (relic density, SModelS and cross
  sections take long).
- **CP phases:** a rephasing test checks phase handling without a second code. An R rotation
  (Hu, Hd R-charge 0, matter 1) moves a phase from M_i and A_f onto mu (M_i, A_f -> e^{-i phi},
  mu -> e^{i phi}, B mu and the vevs unchanged); all masses must stay the same. Native SPheno
  isn't invariant, and FeynHiggs and NMSSMCALC give chargino and neutralino masses only at tree
  level. NMSSMCALCEW (complex Z3 NMSSM, OS or DR-bar electroweakino scheme) gives 1-loop
  electroweakino masses, each from the diagonal self-energies in the mass eigenbasis, which is
  independent of the phase convention. SARAH's SPheno uses the conjugated self-energies for Majorana fermions but not for Dirac
  fermions (`mat1a - SigSL - ...`), so its charginos aren't invariant either.
- **SARAH's diagrammatic 2-loop Higgs masses** (Goodsell et al., [1503.03098]) are in the
  gaugeless limit and the Landau gauge. SPheno solves the tree-level tadpole equations
  with the full gauge couplings and only then sets the broken gauge couplings to zero,
  without solving the tadpoles again (`GaugelessLimit` in `SPhenoLoopMasses.m`,
  `SubSolutionsTadpolesLoop /. subZeroGaugeLess` in `SPheno2LPole.m`). With D-terms, the
  Goldstone bosons are then not massless (in the MSSM, m² = -(g_Y² + g2²) v² cos²2β / 8 (g_Y: hypercharge coupling)),
  and the `EffPot` mass routines keep negative m² with their sign. The defaults are
  `GaugelessLimit = True`, `NewGBC = True` (flag 151, the generalised effective potential
  of [1706.05372]), `include1l2lshift = False` (flag 150, consistent tadpole solution) and
  `TwoLoopRegulatorMass = 0` (flag 410).
  The tachyonic Goldstones are intended: [1411.4665] app. A (eq. A15) derives
  m²(G0) = -M_Z² cos²2β / 2, which keeps the loop functions away from m² = 0, and says the
  loop functions must be evaluated correctly for negative m². [1706.05372] sec. 2.1 calls
  this the original "D-term workaround" for the Goldstone boson catastrophe and says the
  imaginary parts were dropped. In the gaugeless MSSM the Goldstones don't couple to the
  Higgs bosons, so there is no catastrophe there in any case ([1411.0675] sec. 2.3).

## Debugging and profiling

- **Crashes:** `gdb --args models/<M>/run_<M>.x --slha-input-file=...`, then `run` and
  `bt`.
- **Which exception caused a problem:** `gdb -batch -ex "catch throw" -ex run -ex bt --args ...`
  shows where it was thrown (e.g. `integrateOdes` for a failed RG running).
- **"RG integration failed" / "non-perturbative <par>" in SPINFO:** the RK integrator (`src/rk.hpp`)
  throws `RGIntegrationError` (derived from `NonPerturbativeRunningError`) when the step size gets too
  small or it runs out of steps (400), and `NonPerturbativeRunningError` itself only for a non-finite
  step. The parameter named in an integration failure is the one that limited the step size, not
  necessarily a large one: a tiny value (1e-15) means rounding noise, not a Landau pole. To see the
  integration, compile a copy of `src/betafunction.cpp` (with copies of `rk.hpp`,
  `basic_rk_integrator.hpp` and a `config.h` with `ENABLE_VERBOSE` in the same directory, so the
  quoted includes find them) and link it before `libflexisusy.a`.
- **Dependency files and missing headers:** `depgen -MI` silently skips headers that don't exist
  yet, so a `.d` file created before the code generation lists no generated headers (nor the `src/`
  headers included through them). Test `.d` files therefore depend on every model's
  `00_DELETE_ME_TO_RERUN_METACODE` and are recreated after each code generation. In older trees
  (or if results look stale: crashes in `free()` in only some tests), delete `test/*.d`.
- **Name lookups aren't bounds-checked:** `<M>_info::parameter_names`/`particle_names` getters index
  plain arrays. `Problems` maps out-of-range parameter indices to -1 (the scale) since #786, but other
  callers must check the index themselves.
- **Profiling:** `valgrind --tool=callgrind --callgrind-out-file=callgrind.log <exe> <args>`.
  Read the result with `callgrind_annotate` or KCachegrind. `--inclusive=yes` shows where
  the time goes, including callees.
- **Where decay time goes:** in the generated amplitude functions
  (`Decay_amplitude_type<...>`, the vertex evaluations and loop functions inlined in
  them) and in Eigen (SVDs and matrix products). The CXXDiagrams index loops
  (`index_range`, `multiindex.hpp`) are inlined at `-O2` and cost about 2 ns per index
  combination, which is negligible.

## Further documentation

| Topic | Where |
|---|---|
| Model file variables (`FSModelName`, `EWSBOutputParameters`, boundary conditions, `FSDecayParticles`, ...) | `doc/model_file.rst` |
| SLHA input entries (except `FlexibleDecay`) | `doc/slha_input.rst` |
| Tadpoles, EWSB, the parameters used inside the loops (`FSSolveEWSBTreeLevelFor`) | `doc/tadpoles.rst` |
| FlexibleDecay: options, LibraryLink use, comparison with NMSSMCALC | `doc/observables/FlexibleDecay.rst` |
| The Mathematica interface | `doc/librarylink.rst` |
| FlexibleEFTHiggs | `doc/FlexibleEFTHiggs.rst` |
| Adding a loop library | `doc/add_loop_library.rst` |
| Reusable meta code packages (loop functions, RGE integrator, SM/MSSM beta functions) | `doc/meta_code.rst` |
| Observables | `doc/observables.rst` |
| Wiki: code formatting, development build profile, debugging, profiling, unit tests, Wolfram Language coding guide, SARAH pitfalls | https://github.com/FlexibleSUSY/FlexibleSUSY/wiki |

Main references:
- FlexibleSUSY [1406.2319], [1710.03760]
- FlexibleEFTHiggs [1609.00371]
- FlexibleDecay [2106.05038]
- NPointFunctions [2402.14630]
- Unitarity constraints [1805.07306]

The full list, with which papers to cite for which feature, is in `README.rst`.
