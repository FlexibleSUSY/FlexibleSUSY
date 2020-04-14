GM2Calc-1.6.0 [March, 23 2020]
==============================

New features
------------

 * Feature: New `cmake/FindGM2Calc.cmake` file.

Fixed bugs
----------

 * Bugfix (commit c06f713): The function `find_right_like_smuon()` has
   been corrected, which finds the most right-like smuon.  This bugfix
   affects only the SLHA input (the GM2Calc-specific input is not
   affected).  On tested parameter points that are affected, the value
   of amu changes by around 0.01%.

 * Bugfix (commit 5e020ee): In case of SLHA input, the neutralino
   mixing matrix is now read from the NMIX block and used to determine
   the bino-like neutralino.  This change fixes convergence problems
   when converting the DR-bar input parameters to the on-shell scheme.

Changes
-------

 * Change: Separate public from private headers.  The public headers
   are located in `include/gm2calc/`, while the private headers are
   located nearby the respective source files.  The public headers
   must be included as

       #include "gm2calc/gm2_version.h"

 * Change: Private header and source files renamed:

    - Model-specific source files are named as
      `MSSMNoFV_onshell*.{h,hpp,cpp}`.

    - Files specific to the calculation of (g-2) are named as
      `gm2_*.{h,hpp,cpp}`.

 * Change (commit f276703): In case the SLHA input contains multiple
   HMIX blocks, use the *last* one to define the renormalization
   scale.  Input parameters are then searched for in blocks with that
   renormalization scale.

 * Change (commit 9907e86): The Mathematica interface now returns
   `Indeterminate` for amu in case there is a problem with the
   calculation.

 * Change (commit aad392c): Update
   [FindMathematica](https://github.com/sakra/FindMathematica) to
   version 3.2.7

 * Change (commit 059744f): Performance improvement by avoiding
   frequent re-calculation of mb(MZ,DR-bar).

 * Change (commit b38eab7): Update of default fine-structure constant
   from PDG (2019) and
   [[arXiv:1802.02995](https://arxiv.org/abs/1802.02995)].

 * Change (commit 67b39a5): New error message when tan(beta) is too
   large, resulting in floating point overflows.

 * Change (commit 9fdc594, e73c63a): Performance improvement of
   reading SLHA input.

 * Change (commits 4ad9d95, 445ea1c, c7a5e9d): More stingent test of
   GM2Calc configuration options.  Bail out if invalid options are
   given.

 * Change (commits ea03df9, 32df289, 6325c07, 438b421): Performance
   improvement of conversion of SLHA input to GM2Calc-specific
   on-shell scheme.

 * Change (commits d54714f, cf8e702, ec60044, da8d948): Make GM2Calc
   compatible to be used as CMake sub-project.

 * Change (commit 77615d7): Verbose, error and warning messages are
   now written to `cerr`.  This allows to separate the physics output
   from the informational messages.

 * Change (commits ee1eab3, 0635231, 576afcb): Avoid redundant error
   messages to stdout/stderr.  This significantly improves the
   performance of parameter scans using the Mathematica interface,
   where frequent flushing of the stdout/stderr buffers is expensive.

 * Change: Many stylistic internal improvements; fixes of clang-tidy
   warnings.

GM2Calc-1.5.2 [Feature, 08 2020]
================================

 * Change: Updated FindMathematica to version 3.2.6, which can detect
   the Wolfram Engine.

GM2Calc-1.5.1 [January, 21 2020]
================================

 * Feature: New Mathematica file `math/ffunctions.m` with analytic
   expressions for the loop functions.

 * Feature: New Mathematica file `math/amu2Lapprox.m` with analytic
   expressions for the 1- and 2-loop fit formulas.

 * Change: Improved performance of 1-loop functions.

GM2Calc-1.5.0 [May, 11 2019]
============================

 * Feature: New `make install` target, so GM2Calc can be installed.
   The installation includes the executable `gm2calc.x`, the GM2Calc
   library, the public headers and a corresponding `gm2calc.pc` file
   for `pkg-config`.

   Example:

       cmake -DCMAKE_INSTALL_PREFIX=${HOME}/.local/ ..
       make
       make install

 * Feature: The required Boost and Eigen libraries can now be
   installed with [Conan](https://conan.io/):

       mkdir -p build
       cd build

       # install dependencies
       conan install ..

       # invoke cmake as usual
       cmake ..

 * Change (commit 577815e): `gm2calc::Error` now inherits from
   `std::runtime_error`.

 * Change (commit 674cc13): If there is more than one entry with the
   same key in an SLHA block, use the last one.

 * Bugfix (commit 2ea83e8): Catch unphysical input parameter where
   `MW >= MZ`.

GM2Calc-1.4.3 [October, 07 2018]
================================

 * Bugfix (commit cd218cd): Properly initialize trilinear A parameters
   in copy constructor of `MSSMNoFV_onshell`.  
   Note: This copy constructor is never used in GM2Calc.

GM2Calc-1.4.2 [July, 14 2018]
=============================

 * Bugfix (commit e784837c): Fix cmake error from `FindDoxygen.cmake`
   when building with cmake < 3.3, by enabling policy CMP0057.
   Thanks to Sho Iwamoto.

GM2Calc-1.4.1 [June, 13 2018]
=============================

 * Change (commit 2114aef): The default output format of the command
   line program `gm2calc.x` has been changed: If no `GM2CalcConfig[0]`
   entry is provided, then the output is

    * written to the `GM2CalcOutput` block for SLHA input
    * written in detailed form to stdout for GM2Calc input

 * Bugfix (commit 304d771): Catch non-numeric SLHA input.
   Thanks to Peter Athron and the GAMBIT collaboration.

GM2Calc-1.4.0 [March, 27 2018]
==============================

 * Feature: Replace GNU make build system by cmake to improve platform
   independence.  See the `README.md` file for build instructions.

 * Bugfix (commit 8d0bac6): Reame `quad()` function to fix a
   compilation error on Windows/Cygwin.

GM2Calc-1.3.3 [July, 19 2017]
=============================

 * Optimization (commit aa7afc0): Avoid redundant calls to complicated
   `lambda_mu_cha()` function.

 * Bugfix (commit cbb4df2): Fix compilation error on Cygwin where
   `M_PI` might not be defined in `<math.h>`.

 * Bugfix (commits 9d4f768, d820a53, b9a6320, 54970a6): Catch floating
   point overflow/underflow from diagonalization of mass matrices
   during DR-bar to on-shell conversion of right-handed smuon mass
   parameter.

 * Bugfix (commit 54b86e8): A small coefficient in the complex dilog
   has been corrected.  However, the complex dilog is not used so far.

GM2Calc-1.3.2 [February, 22 2017]
=================================

 * Bugfix (commit e8af0c9): Catch potential exceptions from C
   interface functions which calculate amu w/o tan(beta) resummation.

 * Bugfix (commit 654276e): Allow user to pass numeric values to
   `GM2CalcAmuSLHAScheme[]`, `GM2CalcAmuGM2CalcScheme[]` and
   `GM2CalcSetSMParameters[]` which are not numbers, but which would
   evaluate to numbers.

GM2Calc-1.3.1 [January, 31 2017]
================================

 * Bugfix (commit ab73c49): Workaround bug in `mcc` 11.0.0.

GM2Calc-1.3.0 [July, 21 2016]
=============================

 * Feature: A Mathematica interface for GM2Calc has been added.  To
   build it, run `make mathlink`.  The compiled MathLink executable
   can then be found in `bin/gm2calc.mx`.  To use it, the MathLink
   executable must be installed in Mathematica via

       Install["bin/gm2calc.mx"];

   Two examples using the Mathematica interface of GM2Calc can be found in

       examples/example-gm2calc.m
       examples/example-slha.m

   These two examples behave exactly like their C/C++ counterparts.

 * Change (commit b3c7357): Abort calculation if tan(beta) is
   undefined or zero.

 * Change (commit 910b1bd): Abort calculation if the mu parameter is
   zero.  For mu = 0 the approximate fermion/sfermion 2-loop
   corrections are ill-defined.

 * Change (commit 16b181b): Abort calculation if the lightest chargino
   mass is zero.  In that case the photonic 2-loop contribution
   originating from loop function `F3C` is ill-defined.

 * Change (commits 6d2a8b3, 2a8b2ed): Create shared library
   `src/libgm2calc.so` (in addition to the static library) for
   convenience.

 * Change (commit c0b5923): Compile all `.cpp` files in `src/` and put
   all generated object files into the libraries, except for
   `src/gm2calc.o`.  In this way, there is no need for GAMBIT to
   modify `src/module.mk` to add further .cpp files for compilation.

 * Bugfix (commit b3d3c0f): Implementation of the x = 0 limits of the
   `F1N[x]`, `F2N[x]`, `F3N[x]`, `F4N[x]` and `F1C[x]` functions.

 * Bugfix (commit 3ec5df4): Catch potential exception during iteration
   to determine the QCD scale, Lambda_QCD.  This exception is thrown,
   for example, if the Z pole mass is set to be larger than 5 TeV.

GM2Calc-1.2.0 [June, 21 2016]
=============================

 * Feature (commit b65d75d): Adding C interface functions to
   retrieving warnings and problems in form of C strings.

 * Feature (commit ffe9b50): Perform the calculation of the
   uncertainty of a_mu(0-loop) and a_mu(1-loop).  The uncertainty of
   a_mu(0-loop) is the magnitude of a_mu(1-loop).  The uncertainty of
   a_mu(1-loop) is the sum of the magnitudes of a_mu(2-loop,best) and
   the uncertainty of a_mu(2-loop,best).

 * Change (commit fbee2e3): Adding calculation of amu uncertainty to
   the C/C++ example programs.

 * Bugfix (commit 0665a09): Correcting the used sbottom masses in
   Delta_b corrections.

 * Bugfix (commit 71b0cf3): Do not allow the calculation of amu for
   negative soft-breaking squared sfermion mass parameters.

 * Bugfix (commit 7d4d4e7): Use better initial guess for the
   root-finding algorithm, which determines the soft-breaking squared
   mass parameter of the right-handed smuon from the mostly
   right-handed smuon pole mass.

 * Bugfix (commit 9eecf61): Catching potential exception during the
   DR-bar to on-shell conversion of the soft-breaking squared mass
   parameter of the right-handed smuon.

GM2Calc-1.1.2 [April, 08 2016]
==============================

 * Bugfix (commit a374bf7): Reformulate Delta_b corrections to avoid
   numerical problems when Mu, M3 or M1 are zero.

 * Bugfix (commit 4177460): Implement limits of `Fa()`, `Fb()` and
   `amuBmuLmuR()` functions when one or more masses go to zero.

GM2Calc-1.1.1 [March, 29 2016]
==============================

 * Bugfix (commit 9ca7890, 33e68a4, c047d46): Implement limit of
   `Iabc()` function when one or more masses go to zero.

GM2Calc-1.1.0 [December, 14 2015]
=================================

 * Change (commit 5cbbc50): The example programs have been moved to a
   separate `examples/` directory.

 * Change (commit 873133f): Adding a C interface.  The new headers

       src/gm2_1loop.h
       src/gm2_2loop.h
       src/gm2_uncertainty.h
       src/MSSMNoFV_onshell.h

   declare C interface functions for GM2Calc routines.  In addition,
   two example C programs

       examples/example-gm2calc_c.c
       examples/example-slha_c.c

   have been added to illustrate the C interface.  These two C example
   programs behave exactly like their C++ counterparts.  They can be
   compiled by running `make examples`.

GM2Calc-1.0.0 [October, 29 2015]
================================

 * Official release 1.0.0.

GM2Calc-0.2.17 [October, 26 2015]
=================================

 * Feature (commit 132851): The uncertainty of amu(2L,TB resummed) can
   be calculated by setting `GM2CalcConfig[5]` to `1`.  Depending on the
   chosen output format (`GM2CalcConfig[0]`) the uncertainty is written

   * as a single number to stdout   in case of minimal output,
   * to the first line              in case of detailed output,
   * to `GM2CalcOutput[1]`          in case of NMSSMTools output,
   * to `GM2CalcOutput[1]`          in case of SPheno output,
   * to `GM2CalcOutput[1]`          in case of GM2Calc output.

 * Feature (commit 02841d): A new SLHA-compliant output format
   *GM2Calc* has been added to avoid interference with SPheno or
   NMSSMTools: If `GM2CalcConfig[0]` is set to `4`, the value of amu
   is written to `GM2CalcOutput[0]`.  If uncertainty estimation has
   been enabled in addition, the uncertainty is written to
   `GM2CalcOutput[1]`.

 * Feature (commits 25eee8d, 148adaf): Adding two C++ interface
   examples for the SLHA-compliant and the GM2Calc-specific input
   format.  The examples can be found in the files
   `src/example-slha.cpp` and `src/example-gm2calc.cpp`.

 * Change (commit 3bf585a): Memorize all tachyonic particles (not only
   one).

 * Change (commit af21077): Print warning if tachyons exist in
   detailed output without tan(beta) resummation.

 * Change (commits 1426497, f558f22): Reduce number of digits after
   the decimal point to 8 digits in the printed result for amu and its
   uncertainty.

 * Bugfix (commit 4b3f634): Check for tachyons if tan(beta)
   resummation is disabled.

GM2Calc-0.2.16 [October, 05 2015]
=================================

 * Change (commit 7a4f62a): The example input files have been renamed
   such that their file name extension reflects their input format.
   `input/slha.in -> input/example.slha`
   `input/gm2calc.in -> input/example.gm2`

 * Bugfix (commits e784353, c194f9c): Fix compilation problem with
   clang++ on Mac due to different STL implementation of
   `std::conj()`.  Thanks to Björn Sarrazin.

GM2Calc-0.2.15 [October, 01 2015]
=================================

 * Bugfix (commit 585614f): Be less strict in checking if scales of
   two SLHA blocks are the same, as SPheno might write the value of
   the scale with different numerical precision to the block header.
   Thanks to Björn Sarrazin.

GM2Calc-0.2.14 [September, 30 2015]
===================================

 * Bugfix (commit 53187a9): Fix compilation error with g++ 4.6.3.

GM2Calc-0.2.13 [September, 30 2015]
===================================

 * Change (commit 9230209): Print neutralino and chargino mixing
   matrices in verbose output.

 * Change (commit 810c56b): Use mb(MZ) in the DR-bar scheme instead of
   mb(mb) in the MS-bar scheme.

GM2Calc-0.2.12 [September, 27 2015]
===================================

 * Change (commit 6205c0b): The strong gauge coupling g3 is now set to
   the non-zero PDG default value of 0.1184.

 * Change (commit 45b914c): The C++ user interface for changing the
   values of vu and vd has been simplified: Now the user only needs to
   set TB.  The VEV v = sqrt(vu^2 + vd^2) is calculated internally
   using the W and Z pole masses.

 * Change (commit 0659005): Regenerate default SLHA input file (CMSSM
   10.1.1, [[arxiv:1109.3859](https://arxiv.org/abs/1109.3859)]) with
   FlexibleSUSY 1.2.2.  Note: slightly updated SM input parameters are
   used compared to SoftSUSY's CMSSM 10.1.1 version
   inOutFiles/lesHouchesInput .

 * Change (commit 5e993dd): 2-loop fermion/sfermion contributions from
   1st and 2nd generation sleptons have been added to Delta_g1 and
   Delta_g2 (Eqs. (6.6a)-(6.6b)
   [[arxiv:1311.1775](https://arxiv.org/abs/1311.1775)]).  Patch:
   Dominik Stöckinger.

GM2Calc-0.2.11 [September, 24 2015]
===================================

 * Change (commit 16a8118): Don't throw an exception if gluino mass,
   M3, is zero, as it contributes at the 2-loop level.

 * Change (commit 53cbee3): Always print model parameters in verbose
   mode.  In addition, in verbose mode the DR-bar to on-shell
   conversion iteration steps are printed.

 * Change (commit 5854532): Replace Fortran implementation of `Li2(z)`
   by C++ implementation.  The C++ implementation is faster and
   eliminates the need of a Fortran compiler.

 * Change (commit d31daa0): Make use of LAPACK library optional.  By
   default the Eigen library is used for diagonalization of mass
   matrices.

 * Change (commit f105385): If SLHA input format has been chosen and
   no `GM2CalcConfig` input block is provided, the default output
   format will be SLHA and the value of amu will be written to the
   SPheno bock `SPhenoLowEnergy[21]`.

 * Change (commit 3f6d5db): If the conversion of the right-handed
   soft-breaking smuon mass parameter from the DR-bar to OS scheme
   fails using a FPI, a root finding algorithm is tried.

 * Bugfix (commit 2f481f5): The sorting of the smuon pole masses got
   lost during the FPI, which lead to non-convergence of the
   right-handed soft-breaking smuon mass parameter with FPI.

 * Bugfix (commits 47c8237, bd15f57, 64cc024): Read DR-bar parameters
   from SLHA blocks with the same scale.  As renormalization scale the
   renormalization scale of the `HMIX` block is chosen.

GM2Calc-0.2.10 [September, 02 2015]
===================================

 * Bugfix (commit d4baa95): fix compilation with g++ 4.7.4

 * Feature (commit 8993e26): New flag `GM2CalcConfig[4]` to
   enable/disable verbose output.

 * Change: Nicer formatting of detailed output.

GM2Calc-0.2.9 [September, 01 2015]
==================================

 * Bugfix (commit 4a5725a): Use default muon mass if no muon mass is
   given in the `SMINPUTS` block.

 * Bugfix (commit c12dfe3): Refine and implement limits for the
   functions `F*C`, `F*N`, `Fa`, `Fb` and `Iabc`.  Patch: Markus Bach.

 * Bugfix (commit d3cc2de): Determine soft-breaking left-handed smuon
   mass parameter from the muon sneutrino pole mass.

 * Change (commit 6fc7ed4): add 2-loop 2L(a) corrections to best
   approximation for amu.  Patch: Markus Bach.

GM2Calc-0.2.8 [August, 06 2015]
===============================

 * Change (commit 2260ac4): If a physical problem has occured, the
   problem description is added to `SPINFO[4]` in case of
   SLHA-compliant output.

 * Feature (commit bfc22cc): Allow user to force output even if
   physical problem has been spotted.

GM2Calc-0.2.7 [August, 06 2015]
===============================

 * Bugfix (commit c514eca): Implementation of the limit x=1 for the
   functions F3C, F4C, F3N, F4N.

GM2Calc-0.2.6 [August, 06 2015]
===============================

 * Bugfix (commit e16898f): Implementation of the limit x=1 for the
   functions `F1C`, `F2C`, `F1N`, `F2N`.  Patch: Markus Bach.

GM2Calc-0.2.5 [August, 04 2015]
===============================

 * Bugfix (commit 8acb0b7): Make loop order selectable in case of
   minimal output.

 * Change: Do not print any output if the spectrum contains a tachyon.

GM2Calc-0.2.4 [August, 04 2015]
===============================

 * Bugfix: Implementation of functions `Fa(x,y)`, `Fb(x,y)` and
   `H2(x,y)` in the limit x=1 and y=1.

 * Bugfix: Implementation of `I(a,b,c)` in the limit of equal
   arguments.
