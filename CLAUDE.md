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
| `utils/loop_decays/` | FeynArts/FormCalc scripts that produced the generic 1-loop decay diagrams |
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
- **Run times:** code generation takes minutes per model, and longer for large models.
  The slowest steps are usually the vertices, decays, `BrLTo3L` and unitarity.
- **`make -n` isn't a safe dry run:** it can still trigger the meta code, because the
  generated `module.mk` files are make prerequisites. Avoid it in a tree where the
  generated code matters.
- **One target generates all models:** for the same reason, `make` of any target
  (even one model's `00_DELETE_ME_TO_RERUN_METACODE`) first generates every configured
  model whose code is missing or out of date. Without `-j` they run one after another;
  with `-j` several code generations can run at once.
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
  subkernels (`FSEnableParallelism`).
- **Generating and compiling separately:** `./configure --disable-compile` only generates
  the code. `./configure --disable-meta` only compiles code that was generated before,
  without Mathematica or SARAH.
- **Moving generated code:** `make pack-<M>-src` packs one model's generated code. With
  `--with-install-dir=<dir>`, `make install-src` exports a standalone package without the
  meta code.
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
libstdc++ >= 10. `configure` checks only the version numbers.

Before using a newer standard library feature, check it against these minimum
versions:
- **Usable:** `std::span`, `std::ssize`, `std::string::starts_with`, `map::contains`,
  `std::erase_if`, `std::remove_cvref_t`, `<numbers>`, concepts and `requires`.
  `std::span` was tested with g++ 10 and clang 10-12 on libstdc++ 10.
- **Not usable yet:** `<ranges>` (fails with clang 10 on libstdc++ 10), and features
  that need g++ >= 11 or later: `std::format` (13), `std::to_chars` for `double`,
  `std::bit_cast`, `std::source_location`.

## Running a spectrum generator

```sh
models/<M>/run_<M>.x --slha-input-file=models/<M>/LesHouches.in.<M> --slha-output-file=out.spc
```

- **Input blocks:** `Block FlexibleSUSY` configures the calculation (precision, loop
  orders, solver, ...); `Block FlexibleSUSYInput` holds extra physical inputs; and
  `Block FlexibleDecay` configures the decays. Every entry is documented in
  `doc/slha_input.rst`. The EWSB loop order (`FlexibleSUSY[5]`) should equal the pole
  mass loop order (`FlexibleSUSY[4]`).
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
    `SVV`, `SFF`, `FFS`, `FFV`); `SVV` has a Ward identity check.
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

## Conventions

- **Branches:** development happens on `development`; open pull requests against it.
- **`CHANGES.rst`:** bug fixes cite the commit as `[commit 1234abc]`. These hashes
  change if a branch is rebased.
- **Tests with fixes:** add a test with a bug fix where possible, as a Boost.Test
  `.cpp` registered in `test/module.mk`.
- **C++ style:** C++ code follows `.clang-format`
  (`clang-format -style=file -i <file>` from the root). Don't run it on templates:
  it breaks the `@placeholders@`.
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

## Debugging and profiling

- **Crashes:** `gdb --args models/<M>/run_<M>.x --slha-input-file=...`, then `run` and
  `bt`.
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
| Every SLHA input entry | `doc/slha_input.rst` |
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
