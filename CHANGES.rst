FlexibleSUSY 2.5.0 [June, 20 2020]
==================================

New features
------------

* Added support for the `COLLIER <https://collier.hepforge.org/>`_
  [`arXiv:1604.06792 <https://arxiv.org/abs/1604.06792>`_]
  Passarino-Veltman loop function library.  The COLLIER support can be
  activated by the following ``configure`` command::

      ./configure --with-loop-libraries=collier [...]

  See ``./configure --help`` and `README.rst <README.rst>`_ for
  further information.

  Note: COLLIER can be automatically installed via::

      conan install . --build=missing

  The loop function library to use at run-time can be selected by
  setting the flag ``FlexibleSUSY[31]`` accordingly in the SLHA input
  record (0 = SOFTSUSY, 1 = COLLIER, 2 = LoopTools, 3 = FFLite).  In
  the Mathematica interface, chose ``loopLibrary -> [...]``
  accordingly.

  Thanks to Uladzimir Khasianevich.

* Added calculation of :math:`$b \to s \gamma$`.  Currently only
  diagrams with scalars and fermions in the loop are supported.  See
  `doc/observables/b_physics.rst <doc/observables/b_physics.rst>`_ for
  further details.

  Thanks to Kien Dang Tran.

* New calculation of the W boson pole mass with decoupling behaviour
  for large BSM masses.

  Thanks to Markus Bach.

Changes
-------

* [commit c5b7e4a8b]: Rename FlexibleSUSY symbol ``Temporary`` to
  ``FSTemporary`` in order to avoid a conflict with an internal
  mathematica symbol.

* The C++ interface for the one-loop integrals has been changed to
  allow switching between the used loop library at run-time. See
  `doc/add_loop_library.rst <doc/add_loop_library.rst>`_ and
  `src/loop_libraries/loop_library_interface.hpp
  <src/loop_libraries/loop_library_interface.hpp>`_ for further
  technical details.

  Thanks to Uladzimir Khasianevich.

* Changed code organization of ``NPointFunctions`` module: improved
  speed of ``C++`` calculations, improved maintainability of the
  metacode.

* Improved performance of ``flexiblesusy-config`` script.

* Improved performance of 1-loop threshold functions from
  [`arXiv:1407.4081 <https://arxiv.org/abs/1407.4081>`_], used in
  HSSUSY.

* ``make all-test`` now returns early and with a non-zero exit code when a
  test fails.  Use ``make -k all-test`` to force running of all tests.

* When installing the dependencies with Conan_, the `Eigen 3`_ library
  from the Conan repository is preferred over the one installed in the
  system directories.

* Improved performance and compile-time of the 2-loop MSSM threshold
  corrections.

* Improved compile time.

* Updated GM2Calc to version 1.6.0.

Fixed bugs
----------

* [commit 6d4310f6a]: Fix linking error with LoopTools on some
  platforms by linking with libquadmath when necessary.

* Fixed numerical instability of SOFTSUSY's B0 function.

* Fixed run-time error on 32-bit ARM platforms

FlexibleSUSY 2.4.2 [April, 10 2020]
===================================

Fixed bugs
----------

* [commit de7091b0d]: Fixed setting of threshold correction flags with
  clang++ 7.0.

* [commit 23f66e54c]: Fixed compilation error of LibraryLink on
  platforms where `mint = long long`.

FlexibleSUSY 2.4.1 [October, 16 2019]
=====================================

New features
------------

* The 4-loop SM-QCD threshold corrections O(αs^4) to the strong
  coupling `[hep-ph/0512060] <https://arxiv.org/abs/hep-ph/0512060>`_
  can be added by setting ::

      UseSMAlphaS4Loop = True

  in the model file.

* New module ``meta/SM/as_4loop_qcd.m`` with 4-loop SM-QCD threshold
  corrections O(αs^4) to the strong coupling `[hep-ph/0512060]
  <https://arxiv.org/abs/hep-ph/0512060>`_.

* New module ``meta/LoopFunctionsZeroMomentum.m`` with
  Passarino-Veltman 1-loop functions for vanishing external momenta.

Fixed bugs
----------

* [commit c06e57497]: The sign of 2- and 3-loop pure QCD threshold
  corrections for αs in the Standard Model has been corrected.  The
  effect is of the order 50 MeV w.r.t. the Higgs pole mass.

* [commit bedc5b83f]: ``./createmodel`` returned an error when the
  ``models`` directory was empty.


FlexibleSUSY 2.4.0 [August, 04 2019]
====================================

New features
------------

* Implementation of the 4-loop O(αs^4) contributions to the running
  MS-bar top mass of the Standard Model from [`1604.01134
  <https://arxiv.org/abs/1604.01134>`_, `1502.01030
  <https://arxiv.org/abs/1502.01030>`_, `1606.06754
  <https://arxiv.org/abs/1606.06754>`_].  The contributions can be
  enabled in SM-like models by setting the flag::

      UseYukawa4LoopQCD = True

  or::

      UseYukawa4LoopQCD = Automatic

  The 4-loop threshold correction is taken into account at run-time if
  both the global threshold correction loop order flag
  (``FlexibleSUSY[7]`` or ``thresholdCorrectionsLoopOrder``) and the
  individual top Yukawa coupling threshold correction flag
  (``FlexibleSUSY[24]`` or ``thresholdCorrections``) are set to a
  value > 3.

  Example (SLHA input file)::

      Block FlexibleSUSY
          7   4          # global threshold corrections loop order flag
         24   124111321  # individual threshold correction loop orders

  Example (Mathematica interface)::

      fsSettings = {
          thresholdCorrectionsLoopOrder -> 4,
          thresholdCorrections -> 124111321,
          ...
      }

* Implementation of 3-loop contributions O(αb,ατ) to the Standard
  Model beta functions from [`1604.00853
  <https://arxiv.org/abs/1604.00853>`_].

* Implementation of the 2-loop O(αt αs + αt^2) contributions to the
  running MS-bar top mass of the Standard Model from [`1604.01134
  <https://arxiv.org/abs/1604.01134>`_].  The contributions can be
  enabled in SM-like models by setting the flag::

      UseSMYukawa2Loop = True

  Note that FlexibleSUSY must be configured with TSIL_ to use these
  corrections, see `README.rst <README.rst>`_.  Furthermore TSIL_ must
  be compiled with ``-fPIC``, which can be achieved by setting in the
  TSIL_ ``Makefile``::

      TSIL_OPT = -O3 -funroll-loops -fPIC

* The libraries required to build FlexibleSUSY can now be installed
  via the Conan_ package manager, see the `README.rst <README.rst>`_
  for more details.

  Example::

      # install Conan (if not already installed)
      pip install conan

      # add remote repository conan-hep (if not already done)
      conan remote add conan-hep https://api.bintray.com/conan/expander/conan-hep

      # install required libraries
      conan install . --build=missing

* The output of ``make`` is now non-verbose by default.  To enable
  verbose ``make`` output run::

      make VERBOSE=1

* New non-SUSY model LeptoSplitMSSM with light 1st and 2nd generation
  sleptons and light charginos and neutrinos.

  Thanks to Fabian Esser.

Changes
-------

* The C++ language version has been increased to C++14.  As a result,
  a C++14-compatible compiler is required to compile FlexibleSUSY.
  This is the case for

  - g++ >= 5.0.0
  - clang++ >= 3.8.1
  - icpc >= 17.0.0

* The support for BLAS/LAPACK as linear algebra libraries has been
  dropped.

Fixed bugs
----------

* [commit 3b417122]: MSSMD5O model is fixed so that the initial guess
  of WOp does not depend on uninitialized vu.

  Thanks to Andrew Miller.

* [commit c47ef34a]: In function ``SLHA_io::read_entry``, if there is
  more than one entry with the same key in an SLHA block, use the last
  one.  Note, that ``SLHA_io::read_entry`` has not been used in
  FlexibleSUSY so far.

* [commit eac58a54]: Correcting 2-loop O(ατ^2) threshold corrections
  to the quartic Higgs coupling in HSSUSY.

  Thanks to Emanuele Bagnaschi.

* [commits 01c9471e, e9814ffc] Fix linking problem due to libpthread
  not linked on some platforms.

* [commit 41e13c3f] Fix compatibility with SARAH 4.14.2.  The issue
  arose due to a name clash regarding the Mathematica function
  ``CreateParameterList[]``.

FlexibleSUSY 2.3.0 [January, 22 2019]
=====================================

New features
------------

* Implementation of the 5-loop beta function of the strong gauge
  coupling in the SM from [`1606.08659
  <https://arxiv.org/abs/1606.08659>`_].  The 5-loop contribution is
  enabled in all FlexibleEFTHiggs models by default and can be enabled
  in other SM-like models by setting the flag ::

      UseSM5LoopRGEs = True

  in the corresponding model file.

* An internal FeynArts_/FormCalc_ interface has been added, which
  allows for loop calculations inside FlexibleSUSY's meta code.  This
  interface is currently optional and FlexibleSUSY can be run without
  a FeynArts_/FormCalc_ installation.

Changes
-------

* The documentation of FlexibleSUSY has been extended and changed to
  the `reStructuredText <http://docutils.sourceforge.net/rst.html>`_
  format for easier access.  The documentation root file is
  `README.rst <README.rst>`_.  It can be read online at `github
  <https://github.com/FlexibleSUSY/FlexibleSUSY/blob/development/README.rst>`_
  or locally using for example `restview
  <https://mg.pov.lt/restview/>`_::

      restview README.rst

* The unused file ``test/SOFTSUSY/nmssm1loop.f`` has been removed.

* The calculation of the vertices with the ``CXXDiagrams`` module has
  been improved and is now significantly faster.

Fixed bugs
----------

* [commit e5473865]: Take non-standard normalization of VEVs into
  account in FlexibleEFTHiggs models.

* [commit 79651844]: Avoid linker-specific ``--start-group`` and
  ``--end-group`` in order to make the tests build on MacOS.

* [commits 6a4a32324, 2cdd71861, 90ca05d70]: Compatibility fixes for
  SARAH 4.14.0.

* [commits 65aeb9dc, 89b4000b, e7c87c6d]: Ensure phase factors have
  unit modulus when converting a CKM matrix to PDG conventions in the
  case that cos(theta13) vanishes, and add missing Majorana phases in
  the definition of the PMNS matrix.

* [commits 05664d66c, b6073b112]: Refining criterion for the selection
  of the degenerate mass limit of the 2-loop SQCD correction to the
  top mass in the MSSM.  This change improves the numerical precision
  and the stability of the correction for large SUSY scales above 10
  TeV.

* [commits 41d704f05, e0b468e3a]: Correcting implementation of
  analytic ``B00`` function in ``meta/LoopFunctions.m`` for vanishing
  momentum.

* [commits 4a8b249e0, ff0ca140b]: The speed of the conversion of the
  SARAH-generated beta functions to FlexibleSUSY format has been
  improved.  This change is significant for complicated BSM models
  with many couplings.

* [commit e88c1c8ab]: Fix linking with ifort compiled LoopTools_.


FlexibleSUSY 2.2.0 [August, 26 2018]
====================================

New features
------------

* The symbols ``upQuarksDRbar``, ``downQuarksDRbar``, ``downLeptonsDRbar``
  and ``neutrinoDRbar`` are now accessible in all the individual
  low-scale constraint settings (not only in the ones that set the
  SM-like Yukawa couplings ``Yu``, ``Yd``, ``Ye`` and ``Yv``).

* 3-loop corrections of O(αt^2 αs^2) from [`1807.03509
  <https://arxiv.org/abs/1807.03509>`_] to the quartic Higgs coupling
  of the SM can be used in HSSUSY.  To use the corrections,
  FlexibleSUSY must be configured with `Himalaya 2.0.0
  <https://github.com/Himalaya-Library/Himalaya>`_::

      HIMALAYA_DIR=/path/to/Himalaya-2.0.0

      ./configure --with-models=[...] --enable-himalaya \
         --with-himalaya-incdir=$HIMALAYA_DIR/source/include \
         --with-himalaya-libdir=$HIMALAYA_DIR/build

      export LD_LIBRARY_PATH="$HIMALAYA_DIR/build:$LD_LIBRARY_PATH"

Changes
-------

* If unspecified, the pole mass and EWSB loop order is set to 4 and
  the threshold correction loop order is set to 3.  In this way all
  available loop corrections are enable by default.

Fixed bugs
----------

* [commit a9860038a]: Properly treat ``Re[p]`` and ``Im[p]`` in the
  parameter list of ``FSFindRoot[]`` and ``FSMinimize[]``, when the
  parameter ``p`` has indices.

* [commit 85f145a72]: Speed-up the generation of the C++ code in
  models with complicated boundary conditions, like HSSUSY for
  example.

* [commit 00b6798a7]: Speed-up the generation of the C++ code with
  Mathematica 11.3.

  Thanks to Wojciech Kotlarski.

* [commit f2b91c358]: Flag error when non-perturbative corrections
  appear in calculation of Weinberg angle.

* [commits 38e5b30, 5c64b0c]: Pick correct neutrino mass eigenstate
  in Delta_VB calculation in models with neutrino mixing.  This
  avoids a division by zero if neutrinos are strongly mixed.

* [commit 8d508521d]: Correcting 4-loop beta function of the Standard
  Model top Yukawa coupling.  Note: A factor yt is missing in the
  expression betaytl4 in the file
  `ttp16_008.m <https://www.ttp.kit.edu/Progdata/ttp16/ttp16-008/ttp16_008.m>`_.
  This factor yt is present in Eq.(3.5) of
  `1604.00853 <https://arxiv.org/abs/1604.00853>`_.

* [commit a2de8a30d]: Stripping leading/trailing whitespace from
  system directory paths used in the configure script.

* [commit 7da3f384d]: Don't extract local const references from a sum
  of expressions because the sum may accidentally be zero.

* [commit 35d4f1952]: Correcting re-scaling factor of the 2-loop
  singlet tadpole in the NMSSM.

  Thanks to Sebastian Pögel.

* [commit a9047718d]: Disable multi-threading when LoopTools is used
  to avoid race conditions.  Note that LoopTools is thread-unsafe,
  because it accepts the renormalization scale via a global variable.

* [commit e31280702]: Correcting 2-loop O(αt αs) threshold
  corrections in the THDM-like EFTs of the MSSM (THDMIIMSSMBC\*,
  HTHDMIIMSSMBC, HGTHDMIIMSSMBC\*).

  Note: The distribution of the sum

      lambda_{345} = lambda_3 + lambda_4 + lambda_5

  (see Eq.(61) of `1508.00576 <https://arxiv.org/abs/1508.00576>`_)
  onto the individual lambda_{3,4,5} is not unique.  In FlexibleSUSY's
  THDM-like models we now chose the same distribution as in MhEFT_ 1.0
  and 1.1.

  Thanks to Jobst Ziebell.

* [commit 8eba91256]: Correcting ``IsMassless[]`` function for ghost
  fields.

  Thanks to Wojciech Kotlarski.


FlexibleSUSY 2.1.0 [March, 05 2018]
===================================

New features
------------

* Allow user to perform replacements on beta functions,
  self-energies/tadpoles and vertices.  The replacement rules are
  specified as::

      FSBetaFunctionRules = {
          {g1 -> 0, g2 -> 0}, (* applied to 1L beta functions *)
          {g1 -> 0, g2 -> 0}, (* applied to 2L beta functions *)
          {g1 -> 0, g2 -> 0}  (* applied to 3L beta functions *)
      };

      FSSelfEnergyRules = {
          (* applied to 1L self-energies/tadpoles *)
          { (Mass|Mass2)[VZ|gZ] -> 0 }
      };

      (* applied to all vertices *)
      FSVertexRules = { g1 -> 0, g2 -> 0 };

* Adding three new input parameters to the HSSUSY model file which
  can be used to estimate the theoretical uncertainty:

  - By setting ``EXTPAR[201]`` to ``0`` or ``1``, the parametrization of the
    2-loop threshold correction to lambda can be switched between
    yt(SM,MS-bar) and yt(MSSM,DR-bar).

  - By setting ``EXTPAR[202]`` to ``0`` or ``1``, the parametrization of the
    2-loop threshold correction to lambda can be switched between
    DR-bar or on-shell stop mass parameters.

  - A non-zero value of ``EXTPAR[203]`` is interpreted as matching
    scale Q_match.  ``EXTPAR[203] = 0`` corresponds to Q_match =
    MSUSY.

  ================== ==================== ============================
   SLHA input field   Mathematica symbol   Description
  ================== ==================== ============================
   ``EXTPAR[201]``    ``DeltaYt``          0 = yt(SM), 1 = yt(MSSM)
   ``EXTPAR[202]``    ``DeltaOS``          0 = OS stops, 1 = DR stops
   ``EXTPAR[203]``    ``Qmatch``           matching scale
  ================== ==================== ============================

* The Mathematica script
  ``model_files/HSSUSY/HSSUSY_uncertainty_estimate.m`` has been added.
  The script defines the ``CalcHSSUSYDMh[]`` function, which performs an
  uncertainty estimate of the predicted SM-like Higgs mass.  The
  three sources of uncertainty defined in [1504.05200] are taken into
  account: SM uncertainty, EFT uncertainty and SUSY uncertainty.  See
  ``?CalcHSSUSYDMh`` for more information.

* The Mathematica script
  ``model_files/MSSMEFTHiggs/MSSMEFTHiggs_uncertainty_estimate.m`` has
  been added.  The script defines the ``CalcMSSMEFTHiggsDMh[]`` function,
  which performs an uncertainty estimate of the predicted SM-like
  Higgs mass.  Two sources of uncertainty defined in [1609.00371] are
  taken into account: SM uncertainty and SUSY uncertainty.  See
  ``?CalcMSSMEFTHiggsDMh`` for more information.  Note, that there is no
  "EFT uncertainty" in MSSMEFTHiggs, because all 1-loop v^n/MS^n
  terms are included.

* The Mathematica script
  ``model_files/NUHMSSMNoFVHimalaya/NUHMSSMNoFVHimalaya_uncertainty_estimate.m``
  has been added.  The script defines the
  ``CalcNUHMSSMNoFVHimalayaDMh[]`` function, which performs an
  uncertainty estimate of the predicted SM-like Higgs mass.  The
  uncertainty is estimated by:
  1) varying the renormalisation scale,
  2) changing the top Yukawa coupling by higher orders (if available) and
  3) changing the stong coupling by higher orders.

* Implementing 2-loop effective potential contributions to the Higgs
  pole mass in the Standard Model of O((αt+αb)^2 + αb αs + ατ^2).
  These corrections are enabled in all FlexibleEFTHiggs models by
  default and can be enabled in other SM-like models by setting the
  flag
  ::

      UseHiggs2LoopSM = True

  in the corresponding model file.  At run-time these corrections
  (and the O(αt αs) contributions) are enabled when the following
  switches are set in the SLHA input file::

      Block FlexibleSUSY
          4   2   # pole mass loop order
          5   2   # EWSB loop order
          8   1   # Higgs 2-loop corrections O(alpha_t alpha_s)
          9   1   # Higgs 2-loop corrections O(alpha_b alpha_s)
         10   1   # Higgs 2-loop corrections O((alpha_t + alpha_b)^2)
         11   1   # Higgs 2-loop corrections O(alpha_tau^2)

  or in the Mathematica interface::

      FSSMOpenHandle[
         fsSettings -> {
           poleMassLoopOrder -> 2,            (* FlexibleSUSY[4] *)
           ewsbLoopOrder -> 2,                (* FlexibleSUSY[5] *)
           higgs2loopCorrectionAtAs -> 1,     (* FlexibleSUSY[8] *)
           higgs2loopCorrectionAbAs -> 1,     (* FlexibleSUSY[9] *)
           higgs2loopCorrectionAtAt -> 1,     (* FlexibleSUSY[10] *)
           higgs2loopCorrectionAtauAtau -> 1, (* FlexibleSUSY[11] *)
         }, ...
      ];

* Implementation of the strong corrections to the 4-loop beta
  functions of the strong gauge, Yukawa and quartic Higgs couplings
  in the SM from [1508.00912, 1604.00853, 1508.02680].  These 4-loop
  beta function contributions are enabled in all FlexibleEFTHiggs
  models by default and can be enabled in other SM-like models by
  setting the flag
  ::

      UseSM4LoopRGEs = True

  in the corresponding model file.

  Many thanks to Alexander Bednyakov for providing the expression for
  the 4-loop beta function for the strong gauge coupling.

* Implementing 4-loop effective potential contributions to the Higgs
  pole mass in the Standard Model from
  `1508.00912 <https://arxiv.org/abs/1508.00912>`_.  The 4-loop
  contributions are enabled in all FlexibleEFTHiggs models by default
  and can be enabled in other SM-like models by setting the flag
  ::

      UseHiggs4LoopSM = True

  in the corresponding model file.

Fixed bugs
----------

* [commit b85e723]: Fix diagonalization of complex symmetric
  matrices when some of the eigenvalues are degenerate.

* [commit d00b4f5]: Extract electric charge from
  particle-anti-particle-photon vertex if the user has not defined
  the electric charge in SARAH's particles.m file.

* [commits 398bb00, d1699bf]: Correcting CKM mixing of up-type
  quarks at the low-energy scale and of the soft-breaking parameters
  at the GUT scale in the CMSSMCKM.

* [commit c95a9b81]: Fix compilation error in 4-scalar
  couplings in models with scalar color triplets due to an unresolved
  SARAH Clebsch-Gordan coefficient ``CG[__][__]``.

* [commit c053fbab]: Adding missing particle multiplicity
  factor w.r.t. unbroken non-abelian non-SM gauge groups in the
  low-energy 1-loop threshold corrections for αem and αs.
  This bugfix affects BSM models with additional unbroken non-abelian
  gauge groups.

* [commit a1fa5e7f]: Fix multiple local definitions of EWSB
  output parameters in models with complex parameters.

  Thanks to Simonas Drauksas.

* [commit c275b60]: Fix linking error on Darwin platforms.

* [commit dd4292a]: Adding the user-defined ``$(CXXFLAGS)`` to
  the command that creates ``config/depgen.x`` in order to avoid
  linking errors on machines where the ``-parallel`` flag is needed
  during linking.

* [commit 996caef]: Workaround intel compiler / Eigen bug in
  ``allFinite()`` function, which may affect FS output.


FlexibleSUSY 2.0.1 [October, 20 2017]
=====================================

New features
------------

* For each FlexibleSUSY ``<model>`` the TeX file
  ``<model>/<model>_references.tex`` is created, which contains
  ``\cite{}`` commands with references to be cited.  Note, that the
  references to be cited are model-specific due to the different
  switches in the FlexibleSUSY model files.

Fixed bugs
----------

* [commit 682de11c]: Include 2-loop gluino contribution in the
  extraction of yt from the top pole mass in the split-MSSM, Eq.(4.7)
  of `1312.5220 <https://arxiv.org/abs/1312.5220>`_, if
  ``UseHiggs3LoopSplit == True``.

  Thanks to Pietro Slavich.

* [commit a783e318]: Distinguish between tree- and loop-level EWSB
  failures, so problem points where only one fails (but not the
  other) get handled properly.

* [commits 88009cac, 5ac9366c]: Now configure script does not hang
  even if Mathematica fails to find a valid license.  The script does
  not quit even if Mathematica does not meet the version requirement,
  unless ``--enable-meta`` or ``--enable-librarylink`` is given.

  Thanks to Anders Kvellestad.


FlexibleSUSY 2.0.0 [October, 10 2017]
=====================================

New features
------------

* The weak mixing angle can now be calculated from the muon decay
  constant at the full 1-loop level (including flavour mixing
  effects) in a wide range of models.  2-loop corrections of the
  order O(αem αs + αt^2) are taken into account, if
  applicable.

  The method to calculate the weak mixing angle can be chosen in the
  model file by setting the variable ``FSWeakMixingAngleInput`` to
  either Automatic, ``FSFermiConstant`` or ``FSMassW``.  If
  ``FSWeakMixingAngleInput == FSFermiConstant``, then the muon decay
  constant will be used to determine the weak mixing angle.  If
  ``FSWeakMixingAngleInput == FSMassW``, then the W mass will be used.
  If ``FSWeakMixingAngleInput == Automatic`` (this is the default),
  then most precise applicable method is chosen automatically.

  Example::

      FSWeakMixingAngleInput = Automatic; (* recommended *)

  The variable ``FSWeakMixingAngleOptions`` has been removed and can no
  longer be used.

* BSM contributions to the anomalous magnetic moment of the muon, aµ,
  at the 1L level in any given BSM model.  Note: Diagrams with non-SM
  vector bosons are not taken into account.

  In order to let FlexibleSUSY calculate aµ, the symbol
  ````FlexibleSUSYObservable``aMuon```` must be written into an SLHA
  output block in the ``ExtraSLHAOutputBlocks`` variable in the
  FlexibleSUSY model file.

  Example::

      ExtraSLHAOutputBlocks = {
         {FlexibleSUSYLowEnergy,
            {{21, FlexibleSUSYObservable``aMuon}}}
      };

  Thanks to Jobst Ziebell.

* BSM contributions to the electric dipole moment of fermions at the
  1L level in any given BSM model.  Note: Diagrams with non-SM vector
  bosons are not taken into account.

  In order to let FlexibleSUSY calculate the EDM of a fermion F, the
  symbol ````FlexibleSUSYObservable``EDM[F]```` must be written into an
  SLHA output block in the ExtraSLHAOutputBlocks variable in the
  FlexibleSUSY model file.  If F is a multiplet, the field index must
  be specified, for example ````FlexibleSUSYObservable``EDM[F[1]]```` for
  the first field in the multiplet.

  Example::

      ExtraSLHAOutputBlocks = {
         {FlexibleSUSYLowEnergy,
            {{23, FlexibleSUSYObservable``EDM[Fe[1]]},
             {24, FlexibleSUSYObservable``EDM[Fe[2]]},
             {25, FlexibleSUSYObservable``EDM[Fe[3]]} } }
      };

  Thanks to Jobst Ziebell.

* New functions, ``FS<model>GetProblems[]``, ``FS<model>GetWarnings[]``
  and ``FS<model>ToSLHA[]``, have been added to FlexibleSUSY's spectrum
  generator Mathematica interface.  The first two functions return
  details about problems / warnings for the given parameter point.
  The third one formats the output according to the SLHA standard.

* 3-loop beta functions (if available) are calculated in parallel if
  multi-threading is enabled.  This leads to a ~25% speed improvement
  in the MSSM when 3-loop RG running is used.

* Support for SLHA-2 input block ``IMEXTPAR``.

* The full 2-loop O(αs^2) corrections to the DR-bar top and bottom
  Yukawa couplings [hep-ph/0210258, hep-ph/0507139, hep-ph/0707.0650]
  can be added by setting
  ::

      UseMSSMYukawa2Loop = True

  in the model file.

  Thanks to Alexander Bednyakov for providing the expressions.

* The full 2-loop O(αs^2 + αt αs + αb αs) corrections to the strong
  coupling [hep-ph/0509048, arXiv:0810.5101, arXiv:1009.5455] can be
  added by setting
  ::

      UseMSSMAlphaS2Loop = True

  in the model file.

  Thanks to Ben Allanach for providing the expressions, which have
  been extracted from SOFTSUSY 4.0.1.

* The 2- and 3-loop SM-QCD threshold corrections O(αs^2 + αs^3)
  to the strong coupling
  `[hep-ph/0004189] <https://arxiv.org/abs/hep-ph/0004189>`_ can be
  added by setting
  ::

      UseSMAlphaS3Loop = True

  in the model file.

* The SQLite database output now contains the MS-bar/DR-bar mass
  spectrum and mixing matrices, in addition to the pole mass
  spectrum.

* The loop orders of the threshold corrections of the SM(5)
  parameters to the BSM model can now be selected individually by
  using the flag ``FlexibleSUSY[24]`` in the SLHA input file or the
  thresholdCorrections variable in the Mathematica interface.  The
  given value consists of 9 digits, each one representing the
  threshold correction loop order of a parameter, as shown in the
  following table.  The default value is ``123111321``, which
  corresponds to the loop orders given in the table.

  ================== =================================== ==============
   digit position n   default value (prefactor of 10^n)   parameter
  ================== =================================== ==============
   0                  1 (1-loop)                          αem
   1                  2 (2-loop)                          sin(theta\_W)
   2                  3 (3-loop)                          αs
   3                  1 (1-loop)                          m\_Z
   4                  1 (1-loop)                          m\_W
   5                  1 (1-loop)                          m\_h
   6                  3 (3-loop)                          m\_t
   7                  2 (2-loop)                          m\_b
   8                  1 (1-loop)                          m\_τ
  ================== =================================== ==============

* An additional boundary value problem solution algorithm, based on
  expanding the soft SUSY breaking or dimensionful parameters in
  terms of semi-analytic solutions to the RGEs, can now be used to
  calculate the spectrum in a model.

  The boundary value solver algorithms to be used in a model can be
  specified by setting the variable ``FSBVPSolvers`` to be a list
  containing all of the desired solvers in the model file.  By
  default, this is set to ``FSBVPSolvers = { TwoScaleSolver }``,
  corresponding to only the two-scale solver being enabled.

  Example: To enable only the semi-analytic solver instead, the
  model file should contain the setting
  ::

      FSBVPSolvers = { SemiAnalyticSolver };

  Currently, the semi-analytic solver can be used in models where
  all of the parameters to be expanded are fixed in the same
  boundary condition, such as the CMSSM or CNMSSM.

* The 3-loop corrections to the Standard Model Higgs mass of the
  order O(αt^3 + αt^2 αs + αt αs^2) of
  `1407.4336 <https://arxiv.org/abs/1407.4336>`_ can be taken into
  account by setting
  ::

      UseHiggs3LoopSM = True;

  in the FlexibleSUSY model file.  In addition, the pole mass loop
  order must be set to a value greater or equal to 3 to switch the
  corrections on (SLHA input: ``FlexibleSUSY[4]``, Mathematica
  interface: poleMassLoopOrder).  To switch on/off the individual
  3-loop contributions, the SLHA input flags ``FlexibleSUSY[26-29]`` or
  the Mathematica symbols
  ::

      higgs3loopCorrectionAtAsAs
      higgs3loopCorrectionAbAsAs
      higgs3loopCorrectionAtAtAs
      higgs3loopCorrectionAtAtAt

  can be used.

* In the MSSM, the 3-loop corrections O(αt αs^2) and O(αb αs^2) to the
  Higgs pole mass from Ref. `1005.5709
  <https://arxiv.org/abs/1005.5709>`_ can be taken into account.  The
  corrections are taken from the Himalaya package `1708.05720
  <https://arxiv.org/abs/1708.05720>`_.  Himalaya can be downloaded
  from https://github.com/jklappert/Himalaya .

  To build Himalaya, run::

      cd $HIMALAY_PATH
      mkdir build
      cd build
      cmake ..
      make

  where ``$HIMALAY_PATH`` is the path to the Himalaya package.

  To enable the 3-loop corrections in a FlexibleSUSY model, set the
  following flag in the FlexibleSUSY model file::

      UseHiggs3LoopMSSM = True;

  In addition, we strongly recommend to set::

      UseHiggs2LoopMSSM = True;
      EffectiveMu = \[Mu]; (* chose sign convention for mu parameter *)
      UseMSSMYukawa2Loop = True;
      UseMSSMAlphaS2Loop = True;
      UseMSSM3LoopRGEs = True;

  There are already three model files with all these corrections
  enabled: MSSMNoFVatMGUTHimalaya, MSSMNoFVHimalaya,
  NUHMSSMNoFVHimalaya.

  To build the FlexibleSUSY spectrum generator with the 3-loop
  corrections from Himalaya, the location of the Himalaya library and
  the Himalaya header files must be passed to the configure script::

      ./configure --with-models=[...] \
         --enable-himalaya \
         --with-himalaya-incdir=$HIMALAY_PATH/source/include \
         --with-himalaya-libdir=$HIMALAY_PATH/build
      make

  To enable the 3-loop corrections at run-time, the following flags
  should be set in the SLHA input::

      Block FlexibleSUSY
          4   3          # pole mass loop order
          5   3          # EWSB loop order
          6   3          # beta-functions loop order
          7   2          # threshold corrections loop order
          8   1          # Higgs 2-loop corrections O(alpha_t alpha_s)
          9   1          # Higgs 2-loop corrections O(alpha_b alpha_s)
         10   1          # Higgs 2-loop corrections O((alpha_t + alpha_b)^2)
         11   1          # Higgs 2-loop corrections O(alpha_tau^2)
         24   123111221  # individual threshold correction loop orders
         25   0          # ren. scheme for Higgs 3L corrections (0 = DR, 1 = MDR)
         26   1          # Higgs 3-loop corrections O(alpha_t alpha_s^2)
         27   1          # Higgs 3-loop corrections O(alpha_b alpha_s^2)

  In FlexibleSUSY's Mathematica interface, the following settings
  should be used::

      fsSettings -> {
          poleMassLoopOrder -> 3,            (* FlexibleSUSY[4] *)
          ewsbLoopOrder -> 3,                (* FlexibleSUSY[5] *)
          betaFunctionLoopOrder -> 3,        (* FlexibleSUSY[6] *)
          thresholdCorrectionsLoopOrder -> 2,(* FlexibleSUSY[7] *)
          higgs2loopCorrectionAtAs -> 1,     (* FlexibleSUSY[8] *)
          higgs2loopCorrectionAbAs -> 1,     (* FlexibleSUSY[9] *)
          higgs2loopCorrectionAtAt -> 1,     (* FlexibleSUSY[10] *)
          higgs2loopCorrectionAtauAtau -> 1, (* FlexibleSUSY[11] *)
          thresholdCorrections -> 123111221, (* FlexibleSUSY[24] *)
          higgs3loopCorrectionRenScheme -> 0,(* FlexibleSUSY[25] *)
          higgs3loopCorrectionAtAsAs -> 1,   (* FlexibleSUSY[26] *)
          higgs3loopCorrectionAbAsAs -> 1,   (* FlexibleSUSY[27] *)
      }

* Adding complete 1-loop O(ατ + αb) and complete 2-loop O((αt +
  αb)^2 + ατ^2) threshold corrections for lambda(MSUSY) to the HSSUSY
  model file from [arXiv:1703.08166].  Many thanks to Pietro Slavich
  and Emanuele Bagnaschi for providing the expressions.

  Note: 5 new flags are introduced to enable/disable the individual
  2-loop corrections.  In the SLHA input the flags which control the
  inclusion of 2-loop corrections are::

      Block EXTPAR                 # Input parameters
        100   2                    # LambdaLoopOrder
        101   1                    # TwoLoopAtAs
        102   1                    # TwoLoopAbAs
        103   1                    # TwoLoopAtAb
        104   1                    # TwoLoopAtauAtau
        105   1                    # TwoLoopAtAt

  In the Mathematica interface the flags which control the inclusion
  of 2-loop corrections are::

      handle = FSHSSUSYOpenHandle[
         fsModelParameters -> {
            ...
            LambdaLoopOrder -> 2,
            TwoLoopAtAs -> 1,
            TwoLoopAbAs -> 1,
            TwoLoopAtAb -> 1,
            TwoLoopAtauAtau -> 1,
            TwoLoopAtAt -> 1
         }
      ]

* Adding the new input parameter DeltaEFT to the HSSUSY spectrum
  generator to allow the user to estimate the EFT uncertainty.  Each
  1-loop term in the threshold correction for lambda(MS) is
  multiplied by the factor (1 + DeltaEFT v^2/MS^2).  Thus, the
  standard calculation of HSSUSY is obtained by DeltaEFT = 0
  (default).  Set ``DeltaEFT = 1`` to obtain an alternative Higgs pole
  mass with a shifted threshold correction to estimate the effect of
  the missing terms of O(v^2/MS^2).

Changes
-------

* The interface for adding constraints and matching conditions to the
  Two_scale_solver class has been simplified.  Constrains and
  matching conditions are now added using the ``add()`` function.  The
  added constraints and matching conditions are imposed in the given
  order.

  Example: To first impose the low-scale, then the high-scale and
  finally the susy-scale constraint call::

      solver.add(&low_scale_constraint, &model);
      solver.add(&high_scale_constraint, &model);
      solver.add(&susy_scale_constraint, &model);

* The readability of the vertices and the self-energies has been
  improved by using the ``SUM()`` and ``IF()`` macros.

* In FlexibleEFTHiggs models, the Standard Model parameters are
  written to the SLHA output in addition to the BSM parameters.

* The return value of ``FS<model>CalculateSpectrum[]`` and
  ``FS<model>CalculateObservables[]`` has been changed.  They have now
  the structure::

      { <model> -> { model parameters ... } }

  In FlexibleEFTHiggs models, the Standard Model parameters are
  returned in addition and the output has the form::

      { <model> -> { model parameters ... },
        StandardModel -> { ... } }

* In multi-threading mode pole masses are calculated using a thread
  pool instead of spawning threads manually.  This avoids over
  subscription on machines with very few CPU cores.

* The performance of FlexibleEFTHiggs has been improved by around 20%
  by using a faster implementation of the B1 loop function in the
  limit of vanishing momentum and non-zero masses.

* [commit b5cad9e]: Automatically chose the maximum number of EWSB
  and pole mass iterations based on the user-defined precision goal.
  This change leads to a performance improvement for some parameter
  points where the EWSB equations cannot be fulfilled.

* The legacy module and most SOFTSUSY components have been removed.

* All gauge couplings written to the GAUGE block in the SLHA output
  are now unnormalized.  Before, only the hypercharge gauge coupling
  was written as unnormalized in the GAUGE block.  Note: Internally,
  FlexibleSUSY uses normalized gauge couplings only.  In order to
  write the normalized gauge couplings to the SLHA output, a separate
  output block should be created.

  Example::

      ExtraSLHAOutputBlocks = {
         {GAUGENORM, (* contains normalized gauge couplings *)
                 {{1, g1},
                  {2, g2},
                  {3, g3},
                  {4, gN} } }
      };

* [commit 9bfd5f1]: tan(beta)-enhanced contributions to the
  down-lepton Yukawa couplings are now resummed.

* [commit 52faaa7]: The symbol SUSYScaleMatching has been renamed to
  MatchingScaleInput, to express that it is imposed at the matching
  scale, which is in general not equal to the SUSYScale.

* Update to GM2Calc_ 1.3.3.

Fixed bugs
----------

* [commits be8f35b, 3432967]: Improve the stability of the algorithm
  which solves the EWSB conditions.  The more stable algorithm leads
  to a faster convergece of the overall iteration in scenarios where
  the EWSB conditions don't have a solution or the solution is hard
  to find.

  Thanks to Peter Athron and Pat Scott.

* [commit 0cb4042]: Matrix products of the form ``A*B`` in the boundary
  conditions are now interpreted as element-wise products (as in
  Mathematica).  For genuine matrix products use ``MatMul[A,B]`` or
  ``A.B``.

* [commits 9738ba1, 118a9a70]: Catch potential numerical
  instabilities during mass matrix diagonalization, which may result
  in eigenvectors with magnitude larger than 1.


FlexibleSUSY-1.7.5 [September, 05 2017]
=======================================

* Bugfix [commit 03e9265]: Correcting wrong vertex between
  chargino-smuon-neutrino and chargino-muon-sneutrino in muon decay.
  Thanks to Markus Bach.

* Bugfix [commit f3f3850]: Correcting coefficient in complex
  dilogarithm.

* Bugfix [commit d8d8c0c]: Make LibraryLink ``Set[]`` functions accept
  matrix-valued parameters.

* Bugfix [commit 072be7e]: Enable 3-loop RGEs in HSSUSY by default.


FlexibleSUSY-1.7.4 [April, 12 2017]
===================================

* Bugfix [commit f434e30]: Rename internal IndexSum symbol which
  conflicts with SARAH version 4.11.0 and higher.

* Bugfix [commit b8d5dcf]: Correcting gauge-dependent term in 2-loop
  beta function of SM vacuum expectation value after a corresponding
  bugfix in SARAH 4.11.0.  This bugfix affects the Higgs mass
  prediction with FlexibleEFTHiggs by around 10 MeV.


FlexibleSUSY-1.7.3 [February, 27 2017]
======================================

* Change [commit 43bb03a]: FlexibleSUSY now aborts the code
  generation if the user tries to fix an unknown parameter in a
  constraint.  (Before this commit, FlexibleSUSY did only print a
  warning.)

* Change [commit cff40dd]: Catch non-numeric user input to the
  LibraryLink interface functions.

* Bugfix [commit 4a5ada7]: Adding missing return statement in
  function ``recalculate_mw_pole()``.  This bug was only present if the
  W pole mass is used as input (not GF).

* Bugfix [commit bd5ee68]: Correctly handle whitespace in directory
  names inside the configure script and search for headers in
  ``$CPATH`` and ``$CPLUS_INCLUDE_PATH`` .

  Thanks to Joshua Ellis.

* Bugfix [commit bc770ae]: Ensure that phase of (complex) mu
  parameter has magnitude 1 in the CMSSMCPV.
  Thanks to Jobst Ziebell.

* Bugfix [commit beb4683]: Accept SLHA output blocks specified as
  strings (not symbols).
  Thanks to Joshua Ellis.

* Bugfix: Implement missing limits of threshold correction functions
  from arXiv:1407.4081.

* Bugfix [commit 581080f]: Catch further NaNs from inside the MSSM 2L
  Higgs mass routines of Pietro Slavich.


FlexibleSUSY-1.7.2 [December, 15 2016]
======================================

* Feature [commit b052e35]: New flag ``FlexibleSUSY[23]`` to disable
  the pole mass calculation of the non-SM particles.  This flag is
  useful in FlexibleEFTHiggs, when the SUSY scale is so high that the
  non-SM particle masses become unreliable or tachyonic: If a non-SM
  pole mass becomes tachyonic (maybe because the loop corrections
  become too large) FlexibleSUSY would flag the given parameter point
  as unphysical.  However, one might still be interested in the value
  of the SM-like Higgs mass, which is valid in FlexibleEFTHiggs even
  for very large SUSY scales.  In such a case ``FlexibleSUSY[23]``
  could be set to 0 to suppress the calculation of the non-SM pole
  masses.

* Feature [commit 998f11e]: Slightly improved speed of the RG
  running.

* Change [commit 189f508]: Speed-up the calculation of the 2L Higgs
  mass corrections in the MSSM and NMSSM, if multi-threading is used,
  by locking the mutex only for the O(αt αt) corrections.

* Change: The limits sin(2 theta) = 0 and m\_stop1 = m\_stop2 have
  been implemented for the 2L O(αt αs) Higgs pole mass corrections in
  the MSSM to avoid numeric instabilities.

* Bugfix [commit 20f169f]: Re-calculate W pole mass in
  FlexibleEFTHiggs.  Before this commit the electroweak gauge
  couplings in FlexibleEFTHiggs are wrong in scenarios with very
  small αem(MZ) (< 1/1000) and/or a small Z pole mass (< 10 GeV).

* Bugfix [commit 38d17ca]: More reliable convergence criterion for
  FlexibleEFTHiggs for large SUSY scales.  Before this commit, only
  the running BSM masses (at the SUSY scale) have been used as
  convergence criterion.  However, they tend to converge very fast,
  compared to the running SM masses at the electroweak scale.  For a
  more reliable convergence criterion, now both the running BSM and
  SM masses are used.

* Bugfix [commits 5e1b6b3, cc5bfae]: Correction of the 2-loop and
  3-loop QCD corrections to the top pole mass in the Standard Model
  in the MS-bar scheme.  Refs. [hep-ph/9803493, hep-ph/9912391,
  hep-ph/9911434] have expressed the relation between the top pole
  mass and the MS-bar mass in terms of Log[Q^2/Mt^2], where Mt is the
  top pole mass.  Before these commits, FlexibleSUSY used the
  expressions from theses references, but wrote result in terms of
  Log[Q^2/mt^2], where mt is the MS-bar mass, while not accounting
  for the difference between Mt and mt in the logarithms.  This
  bugfix affects the Higgs pole mass at the 3-loop level.

* Bugfix [commit cecff4b]: Flag scalar or vector boson gauge singlet
  tachyons.

* Bugfix [commit 4a3fb5b]: Input tan(beta) at the SUSY scale, instead
  of at the matching scale in the FlexibleEFTHiggs model files.  This
  difference matters when the (unphysical) matching scale is varied
  through ``FlexibleSUSY[19]``.

* Bugfix [commit c35dcb2]: Fixed linking problem of the LibraryLink
  on Mac.

* Bugfix [commits a643be5, cc9ebf1]: Avoid function call ambiguities
  when multiple LibraryLink libraries are loaded into Mathematica at
  the same time.

* Bugfix [commit 1f8e135]: Correcting ``FS<model>Set[]`` function for
  models with matrix-valued parameters.

* Bugfix [commit 4097708]: The generated LibraryLink files are now
  added to the model tarball created by ``make pack-<model>-src``.


FlexibleSUSY-1.7.1 [October, 15 2016]
=====================================

* Change [commit b1efa8c]: Updated to GM2Calc 1.3.0.

* Change [commit 05d8e11]: The loop order of the BSM top Yukawa
  coupling at the scale M_SUSY in FlexibleEFTHiggs is now set
  automatically to match the loop order of the matching condition
  from the SM to the BSM model.

  Before this commit, the user had to set ``FlexibleEFTHiggs[13] = 0``
  and ``FlexibleEFTHiggs[20] = 1`` when yt(BSM) should be calculated
  using 1L QCD corrections.  Analogous, the user had to set
  ``FlexibleEFTHiggs[13] = 1`` and ``FlexibleEFTHiggs[20] = 2`` when
  yt(BSM) should be calculated using 2L QCD corrections.  Now,
  ``FlexibleEFTHiggs[13]`` is set automatically to
  ``FlexibleEFTHiggs[20] - 1`` when yt(BSM) is calculated in
  FlexibleEFTHiggs.

* Change [commit b533d67]: Faster calculation of effective vertices h
  -> photon photon and h -> gluon gluon.

* Bugfix [commit 8b04191]: Improve numerical stability of low-scale
  iteration which determines the SM(5) parameters by using a higher
  RG running precision than the precision goal for the convergence.

* Bugfix [commit 44d2f01]: Print SLHA output even if QedQcd class
  throws an exception.


FlexibleSUSY-1.7.0 [September, 19 2016]
=======================================

* Feature: FlexibleSUSY is now able to generate custom spectrum
  generators using the FlexibleEFTHiggs method described in
  [arXiv:1609.00371].  The following FlexibleEFTHiggs example models
  are provided: CMSSMEFTHiggs, MSSMEFTHiggs, MSSMNoFVEFTHiggs,
  NMSSMEFTHiggs, NUHMSSMaltEFTHiggs, MRSSMEFTHiggs, E6SSMEFTHiggs.  A
  documentation of the new model file options to create a custom
  FlexibleEFTHiggs spectrum generator can be found in
  doc/html/FlexibleEFTHiggs.html .

* Feature: FlexibleSUSY now provides a Mathematica interface for the
  generated spectrum generators.  For each model, an example
  Mathematica script

      models/<model>/run_<model>.m

  is generated, which illustrates the usage.  The documentation of
  the Mathematica interface and several examples can be found in
  FlexibleSUSY's HTML documentation.  Please see the section
  "Creating the soucre code documentation" in the README file for a
  description about how to generate the documentation.

* Change: The configure options for creating dynamic libraries and
  statically linked executable have been changed.  By default, static
  FlexibleSUSY libraries and dynamically linked executables are
  created.

  To generate shared FlexibleSUSY libraries, run::

      ./configure --enable-shared-libs ...

  To generate statically linked executables, run::

      ./configure --enable-static ...

  Please refer to the README file for more information.

* Bugfix [commit 39f8d36]: Fix segfault when multi-threading is used
  in statically linked executables.

* Bugfix [commit 3126ac1]: Catch NaNs from inside the MSSM 2L Higgs
  mass routines of Pietro Slavich.

* Bugfix [commit b6db614]: Correcting 2-loop self energy O(αt^2) in
  the Standard Model.  Before, Eq. (20) of
  `1205.6497 <https://arxiv.org/abs/1205.6497>`_ has been used.
  However, this is incorrect, because it includes 2-loop
  contributions from the momentum iteration of the 1-loop self
  energy, which would be double counted, because FlexibleSUSY already
  does a momentum iteration of the 1-loop self energy.  To fix this,
  Eq. (20) of `1504.05200 <https://arxiv.org/abs/1504.05200>`_ has been
  used, which does not include these 2-loop pieces.


FlexibleSUSY-1.6.1 [August, 28 2016]
====================================

* Bugfix [commit db67c81]: Fix compilation with --disable-threads .


FlexibleSUSY-1.6.0 [August, 27 2016]
====================================

* Feature [commit 4e9ef56]: Allow user to access the beta-functions
  of the model parameters on the r.h.s. of the constraints.  BETA[p]
  represents the beta function of the parameter p using the loop
  level given in the SLHA input.  BETA[l,p] represents the l-loop
  beta function of the parameter p.

  Example in the SM::

      HighScaleInput = {
          {\[Lambda], BETA[g1] + BETA[g2] + BETA[1,Yu][3,3]}
      };

* Feature [commit 5e0bca1]: Allow user to add 3-loop QCD corrections
  of `hep-ph/9912391 <https://arxiv.org/abs/hep-ph/9912391>`_ when
  calculating the top pole mass in non-SUSY models.  The 3-loop QCD
  corrections are added if the flag ``FlexibleSUSY[13]`` is set to 2
  and the pole mass loop order, ``FlexibleSUSY[4]``, is set to a value
  > 2.

  * ``FlexibleSUSY[13] = 0`` and ``FlexibleSUSY[4] > 0``: 1L QCD correction
  * ``FlexibleSUSY[13] = 1`` and ``FlexibleSUSY[4] > 1``: 2L QCD correction
  * ``FlexibleSUSY[13] = 2`` and ``FlexibleSUSY[4] > 2``: 3L QCD correction

* Feature [commits 98bc536, e8fd56a]: Speed up of the RG running in
  models with very complicated beta functions.

* Change [commit 728b5ea]: ``make clean`` no longer removes generated
  source files to avoid the need to re-generate them.  To remove the
  generated files use either::

      make clean-<model>-src # deletes generated files for <model>

  or::

      make clean-generated   # deletes all generated files

* Bugfix [commit a5342eb]: Avoid non-portable use of sed in
  createmodel.  This fixes ``make install-src`` on Mac.

* Bugfix [commit 44b31fa]: Fix potential race condition when
  different model classes that make use of the (N)MSSM 2-loop Higgs
  mass routines of P. Slavich call ``calculate_spectrum()`` at the same
  time.

* Bugfix [commit 0d08b99]: Do not try to generate non-squared unit
  matrices for beta function expressions that must be splitted.
  Non-squared unit matrices did appear for non-squared matrix-valued
  parameters, as for example T[hE] in the SE6SSM.

  Thanks to Dylan Harries.


FlexibleSUSY-1.5.1 [July, 12 2016]
==================================

* Bugfix [commit 63f5361]: Fix numerical instability of SOFTSUSY's B1
  function in parameter regions with p << m1,m2 and m1 close to m2.

* Bugfix [commit fc6d509]: Fix makefile bug in the tarball by
  shipping all .m files that appear in the list of dependencies for
  the generated C++ code.


FlexibleSUSY-1.5.0 [June, 29 2016]
==================================

* Feature: Write phases to SLHA output if a SLHA output block is
  defined for them in the SARAH model file.
  Thanks to Dylan Harries.

* Feature: Allow the user to calculate the pole masses at a fixed
  renormalisation scale at run-time, which is different from the one
  set by the SUSYScale model file variable.  The fixed
  renormalisation scale can be given via the ``FlexibleSUSY[17]`` entry
  in the SLHA input.  ``FlexibleSUSY[17]`` is equivalent to
  ``SPhenoInput[33]`` in SPheno.

* Feature: Updated to GM2Calc 1.2.0.

* Bugfix [commit 9a2d576]: Fix compilation error due to ambiguous
  overload of operator<< .
  Thanks to Dylan Harries.

* Bugfix [commits fc748be, 9654a52]: Fix compilation in case Greek
  Symbols appear in ``If[]`` or ``Which[]`` functions in the model file.
  Thanks to Dylan Harries.

* Bugfix: Fix compilation with g++ 4.4.7.
  Thanks to Dylan Harries.

* Bugfix [commit 6f5e38e]: Correcting convergence criterion in the
  iteration which determines the 1st and 2nd generation running
  fermion masses in the SM(5) at the low-energy scale.  After this
  correction, the running 1st and 2nd generation SM(5) fermion masses
  differ from SOFTSUSY by less than 0.5% at the electroweak scale.


FlexibleSUSY-1.4.2 [May, 09 2016]
=================================

* Bugfix: Correcting handling of spaces in configure script if
  ``/bin/sh`` is ``/bin/dash``.


FlexibleSUSY-1.4.1 [May, 09 2016]
=================================

* Feature: Tab-completion for FlexibleSUSY's spectrum generators and
  scripts in the bash.

  Usage::

      . utils/install-bash_completions.bash

* Feature: For each model an example SLHA input file is generated,
  which is located at models/<model>/LesHouches.in.<model>_generated

* Feature [commit 2b95522]: Allow user to provide specific location
  to libpthread using the --with-pthread-libdir= option.

* Change: The algorithm to determine the running fermion masses and
  gauge couplings has been replaced by a more secure one.  The new
  algorithm performs an iteration between 2 GeV and MZ to fix all
  input parameters at their scale.  The new algorithm leads to
  running 1st and 2nd generation quark masses, which differ from
  SOFTSUSY by around 3%.

* Bugfix [commit 59b867d]: Avoid singularity in the limit MSU^2 /
  M3^2 -> MSQ^2 / M3^2 in HSSUSY.

* Bugfix [commit f3864b8]: Catch exception from SOFTSUSY's QedQcd
  class which are triggered when the input value of Mt_pole is chosen
  to be smaller than MZ_pole.

* Bugfix [commit 077c5b9]: Fixing check for SARAH installation with
  Mathematica 10.

* Bugfix [commit e9954d6]: Fixing numerical instability of SOFTSUSY's
  B0 and B22 functions for very heavy spectra and external small
  momenta.

* Bugfix [commits bcb99bc - 8b5d87e]: Fixing compilation error for
  models which don't have input parameters.

* Bugfix [commits 637d099, 8b3a94f, 2e3a972]: Fixing ``make
  install-src`` in case the path to the FlexibleSUSY contains spaces.

* Bugfix [commits ced2072, 8bc8fdd]: Adding support for further
  debian-based multi-architecture linux distributions in the
  configure script.


FlexibleSUSY-1.4.0 [March, 08 2016]
===================================

* Feature: Allow the user to chose the loop order of the RGEs to be
  generated by SARAH.  This is useful in pure low-energy models,
  where no RGE running is needed, or in very complex models, where
  the generation of the RGEs takes a very long time.

  The RGE loop order can be set in the model file using the
  ``FSRGELoopOrder`` variable.

  Example::

      FSRGELoopOrder = 0; (* no RGEs generated *)
      FSRGELoopOrder = 1; (* only 1-loop RGEs generated *)
      FSRGELoopOrder = 2; (* 1- and 2-loop RGEs generated (default) *)

  Note: The RGE loop order can also be specified at run-time in the
  SLHA input block ``FlexibleSUSY[6]``.

* Feature: FlexibleSUSY no longer requires that the weak mixing angle
  and potential Z-Z' mixing angles are provided in terms of
  Lagrangian density parameters (gauge couplings etc.).  Instead,
  FlexibleSUSY makes use of the DependenceSPheno specification given
  in the SARAH model file to calculate these mixing angles
  numerically.  In this way the effect of gauge boson mixings in
  models with extended gauge groups can be taken into account
  automatically.

  Note: If the weak mixing angle is to be fixed at the low-energy
  scale by the running W and Z masses (see ``FSWeakMixingAngleOptions``
  option) in order to determine the electroweak gauge couplings, then
  an expression for it has to be given in either DependenceNum or
  ``FSWeakMixingAngleExpr`` .

  Example for the MRSSM::

      (* determine weak mixing angle from W and Z masses *)
      FSWeakMixingAngleOptions = FSSetOption[
          FSWeakMixingAngleOptions,
          FSWeakMixingAngleInput -> FSMassW
      ];
      (* need to provide expression for weak mixing angle *)
      FSWeakMixingAngleOptions = FSSetOption[
          FSWeakMixingAngleOptions,
          FSWeakMixingAngleExpr  -> ArcSin[Sqrt[1 - (Mass[VWm]^2 - g2^2*vT^2)/Mass[VZ]^2]]
      ];

  Important note: In the SARAH model file a mass ordering of the
  vector bosons is assumed.  For example, the statement
  ::

      DEFINITION[EWSB][GaugeSector] = {
          {{VB,VWB[3],VBp}, {VP,VZ,VZp}, ZZ},
          ...
      };

  assumes MVP < MZ < MZp.  Thus, the user has to make sure that the
  studied parameter region leads to Photon, Z and Z' masses which are
  in agreement with the relation MVP < MZ < MZp.  Otherwise, the
  calculated Z and Z' masses will be incorrect.  If a parameter
  region shall be studied where MVP < MZp < MZ, then the ordering of
  vector bosons in the SARAH model file has to be changed to
  ::

      DEFINITION[EWSB][GaugeSector] = {
          {{VB,VWB[3],VBp}, {VP,VZp,VZ}, ZZ},
          ...
      };

* Feature: By setting the entry ``FlexibleSUSY[16] = 1`` in the SLHA
  input file, the user can force majorana fermion masses to be
  positive.  In this case, the corresponding mixing matrix is not
  purely real and its imaginary part will be written to the output in
  addition.  Note, that setting ``FlexibleSUSY[16] = 1`` is therefore a
  violation of the SLHA standard.

* Feature: FlexibleSUSY calculates the effective 1-loop couplings of
  the CP-even and CP-odd Higgs -> photon + photon and Higgs -> gluon
  + gluon.
  Author: Dylan Harries

  For each model the <model>_effective_couplings class is generated
  and can be used at the C++ level to calculate the effective
  couplings.  In order to write the effective couplings to the SLHA
  output, extra SLHA output blocks have to defined in the
  FlexibleSUSY model file, which contain the symbols
  ::

      FlexibleSUSYObservable``CpHiggsPhotonPhoton
      FlexibleSUSYObservable``CpHiggsGluonGluon
      FlexibleSUSYObservable``CpPseudoScalarPhotonPhoton
      FlexibleSUSYObservable``CpPseudoScalarGluonGluon

  Example:

  Definition of an extra SLHA output block named ``EFFHIGGSCOUPLINGS``,
  containing the effective 1-loop CP-even and CP-odd Higgs -> photon
  + photon and Higgs -> gluon + gluon couplings::

      ExtraSLHAOutputBlocks = {
         {EFFHIGGSCOUPLINGS,
                 {{1, FlexibleSUSYObservable``CpHiggsPhotonPhoton},
                  {2, FlexibleSUSYObservable``CpHiggsGluonGluon},
                  {3, FlexibleSUSYObservable``CpPseudoScalarPhotonPhoton},
                  {4, FlexibleSUSYObservable``CpPseudoScalarGluonGluon} } }
      };

  The calculation of the effective couplings can be disabled (or
  enabled) by setting the flag ``FlexibleSUSY[15]`` to ``0`` (or ``1``) in
  the SLHA input file.

* Feature: Allow the user to temporarily re-define model parameters
  in the boundary conditions, which are restored to their original
  values after the calculations in the boundary condition has been
  finished.

  Example: Temporarily scale the gauge coupling g1 by a factor 1/2
  and set the up-quark Yukawa coupling to zero::

      LowScaleInput = {
         {FSTemporary[g1], g1 / 2},
         {FSTemporary[Yu[1,1]], 0},
         ...
      };

* Feature: The three THDM-like models, which have been used in
  `1512.07761 <https://arxiv.org/abs/1512.07761>`_, are provided.  The
  models implement the 1- and 2-loop threshold corrections of
  `1508.00576 <https://arxiv.org/abs/1508.00576>`_ and
  `hep-ph/9307201 <https://arxiv.org/abs/hep-ph/9307201>`_.  The models
  are named:

  * THDMIIMSSMBC (THDM with boundary condition to the MSSM)
  * HTHDMIIMSSMBC (THDM + Higgsinos with boundary condition to the
       MSSM)
  * HGTHDMIIMSSMBC (THDM + Higgsinos + gauginos with boundary
       condition to the MSSM)

* Feature: In non-SUSY models the 3-loop (Standard Model) QCD
  corrections to the MS-bar Yukawa coupling of the order O(αs^3)
  [hep-ph/9911434, hep-ph/9912391] are added automatically.  They are
  taken into account at run-time if the threshold correction loop
  (``FlexibleSUSY[7]``) order is set to a value > 2 in the SLHA input
  file.

  The generation of 3-loop QCD corrections can be disabled by setting
  in the model file
  ::

      UseYukawa3LoopQCD = False;

* Change [commit f2f913e, 002c904]: When threshold corrections are
  disabled, the charged lepton and top quark pole masses are used to
  determine the corresponding Yukawa couplings.  Before commit
  f2f913e, the running Standard Model masses were used.  This change
  makes it easier to compare the mass spectrum with SPheno when
  threshold corrections are disabled.

* Change [commit 1c7e4a7]: The 2-loop QCD contribution to the top
  Yukawa coupling [hep-ph/0210258 Eq. (60)-(61), hep-ph/9803493
  Eq. (17)] is taken into account only if the threshold correction
  loop order (flag ``FlexibleSUSY[7]``) is set to a value > 1.  Before
  commit 1c7e4a7 the 2-loop QCD contribution was always taken into
  account and could not be disabled.  This change allows the user to
  consistently disable 2-loop contributions.

* Bugfix [commit f7ff872]: Support models which have couplings
  proportional to the epsilon tensor in color space.

* Bugfix [commit 8c1ca39]: Enabling support to use
  LowEnergyConstant[MZ] as scale for the susy-scale contraint.
  LowEnergyConstant[MZ] will be replaced in the C++ code by the
  user-defined SLHA input value of the Z pole mass.

  Example::

      SUSYScale = LowEnergyConstant[MZ];

* Bugfix [commit 0a7934e]: Fix compilation error in models in which a
  multiplet exists, which consists only of Goldstone bosons.

* Bugfix [commit a87042f]: Rename enum entries for matrices to
  prevent compilation errors in models which have mixing matrices
  larger than 10x10.

* Bugfix [commit 61fb1ca]: Fix compilation errors in models which
  don't contain SM-like neutrinos.

* Bugfix [commit 919347d]: Correcting the phase of Dirac fermion
  singlets if their mass is less than zero: Before commit 919347d,
  the phase of Dirac fermion singlets was set to e^(i Pi/2) if their
  mass is less than zero, which is wrong, because in SARAH only one
  Weyl component of the Dirac spinor receives a phase.  After this
  commit, the phase of Dirac fermion singlets is set to e^(i Pi) if
  their mass is less than zero.

* Bugfix [commits 060b492, a6f7741, 306385b]: Implement massless
  limits in C0, D0 and D27 functions.

* Bugfix [commits d62886d]: Ensure that only Standard Model goldstone
  bosons are removed to obtain "heavy" W and Z self-energies.

* Bugfix [commits 60d68af]: Fix compilation error in models where the
  left-handed electron and neutrino mass matrices are of equal size,
  but larger than 3x3.


FlexibleSUSY-1.3.2 [January, 10 2016]
=====================================

* Bugfix [commit d76ca79]: Fix compilation error with Eigen
  3.3-beta1.


FlexibleSUSY-1.3.1 [January, 08 2016]
=====================================

* Bugfix [commit aa8dc76]: Re-enable the output of gauge eigenstate
  masses of 1st and 2nd generation sfermions in the CMSSMNoFV for
  SLHA-1 compatibility.


FlexibleSUSY-1.3.0 [January, 08 2016]
=====================================

* Feature: The output of the spectrum generator can be written into
  an SQLite database using the ``--database-output-file=`` option.  At
  the C++ level, a ``to_database()`` and ``from_database()`` function is
  provided for each model, which write/read a model object (including
  the DR-bar parameters and the pole mass spectrum) to/from a
  database file.

  Example::

      models/CMSSM/run_CMSSM.x \
        --slha-input-file=model_files/CMSSM/LesHouches.in.CMSSM \
        --slha-output-file= --database-output-file=point.db

  Example using the scan script::

      utils/scan-database.sh \
        --spectrum-generator=models/CMSSM/run_CMSSM.x \
        --slha-input-file=model_files/CMSSM/LesHouches.in.CMSSM \
        --scan-range=MINPAR[3]=1~30:21 \
        --database-output-file=scan.db

* Feature: Models can now be matched to the Standard Model at Q =
  MZ_pole, Q = MT_pole or any other dynamically calculated scale, as
  MT_DRbar for example.
  To match at MZ_pole set in the model file: LowScale = LowEnergyConstant[MZ].
  To match at MT_pole set in the model file: LowScale = LowEnergyConstant[MT].
  To match at MT_DRbar set in the model file: LowScale = M[Fu[3]],
  depending on your chosen name for the top quark.

* Feature: 3-loop beta-functions can now be used in the real MSSM.
  To enable the 3-loop MSSM beta-functions, set UseMSSM3LoopRGEs =
  True; in the model file (enabled by default in all real MSSM models
  that are shipped with FlexibleSUSY).  The expressions have been
  obtained from http://www.liv.ac.uk/~dij/betas/allgennb.log and
  include family mixing.

  Note: The 3-loop beta-functions for the vacuum expectation values
  vu and vd are not available so far.  Furthermore, the 3-loop MSSM
  beta-functions miss the "tadpole" contributions corresponding to
  the renormalisation of the Fayet-Iliopoulos D-term, see the note in
  Section 2, page 4 of
  `hep-ph/0308231 <https://arxiv.org/abs/hep-ph/0308231>`_.

* Feature: The anomalous magnetic moment of the muon, (g-2)/2, can be
  calculated in all MSSM models without sfermion flavour violation
  (e.g. the MSSMNoFV and CMSSMNoFV).  The calculation is performed
  with GM2Calc 1.1.0 [arXiv:1510.08071] up to the 2-loop level
  including tan(beta) resummation.

  In order to enable the calculation of (g-2)/2, the symbols
  ::

     FlexibleSUSYObservable``aMuonGM2Calc
     FlexibleSUSYObservable``aMuonGM2CalcUncertainty

  have to be added to ExtraSLHAOutputBlocks variable in the
  FlexibleSUSY model file (they are already added in the MSSMNoFV and
  CMSSMNoFV example models).  In addition, the SLHA input file entry
  ``FlexibleSUSY[15]`` has to be set to 1 to perform the calculation.
  If ``FlexibleSUSY[15]`` is set to 0, (g-2)/2 is not calculated.

* Change [commit d553af8]: No SLHA output is written if the option
  --slha-output-file= is set to the empty string.  To write the SLHA
  output to stdout, set --slha-output-file=- (this is the default).

* Change [commit ac70fec]: In the SM the Higgs pole mass is no longer
  calculated at the scale Qin (= the scale where lambda is input),
  but at the scale Q = M_top.

* Bugfix [commit 1b4fc20]: Correcting W contribution in beta-function
  of α_em in the SM with 5 active quark flavours.
  Imported from SOFTSUSY [commit 0139daa).

* Bugfix [commit d7dbeb6]: Adding neutrino charge, Qv, to list of
  input parameters in the UMSSM.  This fixes a compilation error with
  SARAH 4.6.0.

* Bugfix [commit f1752a7]: Correcting the trilinear couplings and the
  effective mu parameter in the NMSSMRUN SLHA output block in the
  models: NMSSM, NMSSMCPV, NUTNMSSM, SMSSM and NUTSMSSM.

* Bugfix [commit 9ccdb4d]: Workaround a SARAH issue where the list
  SARAH``Masses[EWSB] contains replacement rules of the form ``0 ->
  MassGiven[X]``, instead of ``Mass[X] -> MassGiven[X]``.  Due to this
  issue some massless particles have been missing in FlexibleSUSY
  before commit 9ccdb4d.


FlexibleSUSY-1.2.4 [October, 27 2015]
=====================================

* Change [commit 33af37c]: The spectrum generator, run_<model>.x,
  will no longer overwrite the user-given input parameters of the
  SMINPUTS block.

* Bugfix [commit 9067f3a]: There was an internal programming error in
  the meta code concerning the assignment of tadpole diagrams to the
  Higgs fields, which resulted in a compilation error in the SSM.
  Thanks to John McDowall.

* Bugfix [commit 77d2a86]: Ensure that in the calculation of the pole
  mass of a fermion singlet the prefactor of the self-energies is the
  positive tree-level mass.  Before commit ce1ef83, the prefactor of
  the gluino self-energies in MSSM for example was the soft-breaking
  parameter M3.  If M3 < 0 the gluino pole mass was not calculated
  correctly.
  Thanks to Dylan Harries and Roman Nevzorov.


FlexibleSUSY-1.2.3 [October, 18 2015]
=====================================

* Feature: Adding support for ``If[]`` and ``Which[]`` statements at the
  r.h.s. of contraints.  In addition, the IsClose[a,b,eps] and
  IsCloseRel[a,b,eps] functions have been added to allow for a
  comparison of parameters.

* Feature: New model SplitMSSM, which implements low-energy EFT of
  the MSSM where the sfermions and one Higgs doublet have been
  integrated out.  The model implements the 1- and 2-loop matching
  conditions from `1407.4081 <https://arxiv.org/abs/1407.4081>`_.  The
  Higgs pole mass is calculated at complete 1-loop order plus 2-loop
  contributions O(αt^2) and O(αt αs) from
  `1205.6497 <https://arxiv.org/abs/1205.6497>`_ plus 3-loop
  leading-log contribution from the gluino O(αt αs^2)
  `1312.5220 <https://arxiv.org/abs/1312.5220>`_.

* Feature: New model HSSUSY, which implements a high-scale SUSY
  scenario, where the sfermions, the gauginos, the Higgsinos and one
  Higgs doublet have been integrated out, leaving the Standard Model
  as low-energy EFT.  The model uses the 3-loop Standard Model RGEs
  [1303.4364, 1307.3536] and implements the 1- and 2-loop matching
  conditions to lambda(MSUSY) from
  `1407.4081 <https://arxiv.org/abs/1407.4081>`_.  Furthermore, the
  1-loop matching conditions O(αb) and O(ατ) as well as the 2-loop
  matching condition O(αt^2) from SUSYHD
  `1504.05200 <https://arxiv.org/abs/1504.05200>`_ are implemented.
  The Higgs pole mass is calculated at complete 1-loop order plus
  2-loop contributions O(αt^2) and O(αt αs) from
  `1205.6497 <https://arxiv.org/abs/1205.6497>`_.  The calculation of
  the Higgs pole mass in the HSSUSY model coincides with the one
  obtained with SUSYHD 1.0.2 with a relative deviation of < 0.06%.

* Feature: Allow adding 3-loop gluino contribution to Higgs
  self-energy in split-SUSY models with a physical singlet Higgs.
  The 3-loop gluino contribution is enabled by default in the
  SplitMSSM.

* Change [commit f7cd242]: The ``test`` and ``examples`` modules are no
  longer loaded into the makefile by default.  To load them, run
  ./configure --with-optional-modules="test,examples"

* Change [commit e86d23a]: The FlexibleSUSY test suite is no longer
  shipped with the release tarball.  It can be obtained from the
  official git repository at
  https://github.com/FlexibleSUSY/FlexibleSUSY .

* Change [commit 372bb96]: Use FlexibleSUSY's own dependency file
  generator instead of using the corresponding compiler capabilities.

* Bugfix [commit 20e88db]: Use correct self-energy for 1st and 2nd
  generation charged leptons in \*NoFV models.  Before commit 20e88db,
  the (heavy) tau self-energy was used to convert the running MS-bar
  electron and muon masses to DR-bar masses in \*NoFV models.
  Corresponding test case:
  ``test_CMSSMNoFV_low_scale_constraint::test_delta_Yf()``


FlexibleSUSY-1.2.2 [September, 08 2015]
=======================================

* Feature: The scale at which the EWSB output parameters are fixed
  can now be chosen by the user via the ``FSSolveEWSBFor[{...}]``
  symbol.  By default, the susy-scale is chosen.

* Change [commit 5b9d653]: If ./configure is run without the
  ``--with-models=<models>`` argument, no models will be build.  In
  former FlexibleSUSY versions if the ``--with-models=<models>``
  argument was missing, all models were build.

* Bugfix [commit 5530bf9]: Defining a scale to be a running mass, for
  example SUSYScale = M[hh], resulted in a compilation error.

* Bugfix [commits 2d6c0d2, 87cfe28]: use SLHA input value of the Z
  pole mass as low-energy scale, instead of the hard-coded value MZ =
  91.1876 GeV.

* Bugfix [commit 1ac0aa0]: Use math/physics index convention (index
  starting with 1) in the comments of the extra user-defined SLHA
  output blocks.

* Bugfix [commit 0737c4d]: Properly convert greek symbols in function
  arguments.  Fixes #5.  Thanks to Dylan Harries.

* Bugfix [commit f4eed5d]: Put class Complex into softsusy namespace
  to avoid ambiguities in ``operator*()``.  Fixes #6.  Thanks to Dylan
  Harries.


FlexibleSUSY-1.2.1 [July, 07 2015]
==================================

* Feature: The model name is printed in SPINFO[5] and the SARAH
  version is printed in SPINFO[9].

* Bugfix (fea4d59]: The MODSEL block was not read if SLHA input is
  passed to the spectrum generator via stdin.  Thanks to Peter
  Drechsel.


FlexibleSUSY-1.2.0 [June, 26 2015]
==================================

* Feature: Allow the user to add 3-loop beta-functions in the SM.
  The beta-functions are taken from SUSYHD v1.0.1 (arXiv:1504.05200)
  and `1303.4364 <https://arxiv.org/abs/1303.4364>`_.

* Feature: Allow the user to add 2-loop Higgs self-energy corrections
  O(αt^2 + αt αs) in the SM.  The self-energy corrections were taken
  from `1205.6497 <https://arxiv.org/abs/1205.6497>`_.

* Feature: Allow the user to provide SLHA input via stdin if the SLHA
  input file name is set to - .

  Example::

     cat model_files/CMSSM/LesHouches.in.CMSSM | \
        models/CMSSM/run_CMSSM.x --slha-input-file=-

* Feature: Allow the user create standalone executables that don't
  depend on dynamically linked libraries.  See README for more
  details.

* Bugfix [commit 3843ea7]: Rewrite pole mass tachyon check to fix a
  confusion between goldstone and Higgs bosons in the CP-violating
  MSSM.

* Bugfix [commit e2009f7]: Adding missing declaration of input
  parameters in the generated DependenceNum functions.  This fixes a
  compilation error in the NE6SSM or the UMSSM if ThetaWp is set to
  an expression that involves the charges.

* Bugfix [commits d80c30f, e6c8dda]: Correcting input scale of
  tan(beta) in the lowNMSSM according to SLHA-2 convention.  The
  model file lowNMSSMTanBetaAtMZ has been added, where tan(beta) is
  input at MZ.


FlexibleSUSY-1.1.1 [June, 08 2015]
==================================

* Bugfix [commit e1ea433]: Catch NaNs from Slavich's NMSSM 2-loop
  self-energies.


FlexibleSUSY-1.1.0 [May, 31 2015]
=================================

* Feature: Calculation of DR-bar weak mixing angle from Fermi
  constant and Z pole mass.  The implementation is based on
  expressions from SOFTSUSY and works for the SM, MSSM, NMSSM and
  their variants.  The method for the calculation of the weak mixing
  angle can be selected via the ``FSWeakMixingAngleInput`` variable in
  the FlexibleSUSY model file.

  Example::

      FSWeakMixingAngleInput = FSFermiConstant; (* or FSMassW *)

  Note: To achieve the maximum accuracy available, set the threshold
  corrections loop order to 2 (FlexibleSUSY block entry 7)

* Feature: Support for non-SUSY models, renormalized in the MS-bar
  scheme.

* Feature: 2-loop QCD corrections can be added when calculating the
  top pole mass from the top DR-bar mass.  These 2-loop contributions
  can be enabled/disabled using entries 13 or 4 of the FlexibleSUSY
  block in the SLHA input file.

* Feature: In the shipped FlexibleSUSY model files, the corresponding
  default SARAH model file is specified.  This allows a user to
  create a new model with the simplified command::

      ./createmodel --name=CMSSM

  The default SARAH model file to be used with a given FlexibleSUSY
  model file can be set via ``FSDefaultSARAHModel = <model>``

* Feature: Complex model parameters are now supported.

* Feature: The CKM and PMNS matrix can now be used as low-energy
  inputs.  They are read from the VCKMIN and UPMNSIN input blocks,
  respectively.  Linked to this, the new model file CMSSMCKM was
  added to demonstrate the input of the CKM matrix at low energies.

* Feature: Mark parameter points as invalid, for which the
  calculation of one of the pole masses failed due to
  non-convergence.

* Feature: New (non-templated) intermediate model class
  <model>_mass_eigenstates, which is able to calculate the pole and
  running mass spectrum.  <model>_mass_eigenstates is derived from
  <model>_soft_parameters.  The templated model class
  <model><Two_scale> is now derived from <model>_mass_eigenstates .

* Bugfix [commit 6da2cbd, 8113e32a]: ensure that the MSSM-like CP-odd
  Higgs mass is used in the two-loop Higgs self-energies and
  tadpoles.  Before, there were cases where a Goldstone boson mass or
  a singlet-like pseudoscalar mass was used.

* Bugfix [commit 29a0833]: incorporate tadpole contributions in pole
  masses of singlets

* Bugfix [commit c64a333]: Softsusy's B1 function is now thread-save.
  Before commit c64a333, the τ pole mass was varying due to a race
  condition, if multi-threading is enabled and neither fflite nor
  looptoos is used.

* Bugfix [commit d035544]: Ignore trivial EWSB eqs.  Makes the MRSSM
  work in FlexibleSUSY with SARAH 4.5.x.

* Bugfix [commit d8a1521]: The ``SM()`` preprocessor macro has been
  renamed to ``LowEnergyConstant()`` in order to avoid collisions with
  the copy constructor of the SM model class.

* Bugfix [commit 0c7a7ac]: chop beta-function values smaller than the
  zero-threshold to avoid failures of the RK integrator.  The
  zero-threshold is 1e-11 by default and can be changed via
  ``Beta_function::set_zero_threshold()`` or entry 14 in the SLHA input
  file.

* Bugfix [commit 29a1578]: Ignore goldstone boson "pole mass"
  tachyons.


FlexibleSUSY-1.0.4 [January, 15 2015]
=====================================

* Add new user example program run_cmd_line_<model>.x to run a
  parameter point using command line parameters instead of an SLHA
  input file.

* Allow input parameters in first guesses of scale definitions, for
  example
  SUSYScaleFirstGuess = Sqrt[Sqrt[LHInput[mq2[3,3]] * LHInput[mu2[3,3]]]]

* Adding support for FlexibleSUSY addons.  They are placed inside the
  addons/ directory and can be configued and compiled via
  ``./configure --with-addons=<addon> && make``

* Adding EWSB solvers using a fixed-point iteration (FPIRelative,
  FPIAbsolute, FPITadpole).  FPIRelative is now the first default
  solver used.  Thanks to Dylan Harries.

* Adding new NMSSM model file ``NUTNMSSM`` with non-universal soft
  Higgs masses (EWSB output) and non-universal trilinear couplings
  A_lambda and A_kappa at MX.

* Read user input W boson pole mass form SMINPUTS block entry 9.

* Read user input Z boson pole mass from SMINPUTS block entry 4.

* Automatic check for non-perturbative dimensionless model parameters
  at the high-scale.  The check can be disabled by stetting
  ``FSCheckPerturbativityOfDimensionlessParameters = False`` in the
  model file.  The threshold can be set via the
  ``FSPerturbativityThreshold`` variable.  The default threshold is
  ``N[Sqrt[4 Pi]] = 3.54491``.

* Check for tree-level tachyons at each scale (MZ and M_SUSY)

* Allow to force SLHA output for unphysical points (for example where
  tachyons exist) in FlexibleSUSY block, entry 12.

* Bugfix [commit 6f7d3de]: allow plain model parameters for scale
  definition, for example in the form SUSYScale = vu .

* Bugfix [commit 44baa73]: allow model parameters in first guesses of
  scale definitions, for example
  SUSYScaleFirstGuess = Sqrt[mq2[3,3] * mu2[3,3]]

* Bugfix [commit 77dce8b]: correct momentum guess for the calculation
  of the self-energies with LowPrecision.

* Bugfix [commit fb0b906]: Fix compilation with g++ 4.5.3.

* Support Intel C++ compiler versions 12.1 and 13.x [commits 78d73e7
  and bf5a08e)

* More descriptive error message when an exception is thrown.


FlexibleSUSY-1.0.3 [November, 21 2014]
======================================

* Allow selection of Higgs 2-loop contributions in SLHA input file

* Allow extra user-defined SLHA output blocks

* Allow user-defined matrix- or vector-like SLHA input parameters

* Support low-energy quark flavour violation via CKM matrix

* Bugfix [commit 5f78968]: perform residual color contractions before
  stripping group factors.  Thanks to Philip Diessner and Wojciech
  Kotlarski.

* Bugfix [commit 7160095]: Correcting check for tachyons in pole
  masses of scalar particles, calculated with LowPrecision

* Bugfix [commit a7a33d3]: Implement reading of data from multiple
  SLHA blocks with the same name.  Subsequent block entries will
  overwrite former entries.


FlexibleSUSY-1.0.2 [July, 15 2014]
==================================

* Bugfix [commit 689141da]: Enable non-quadratic superpotential
  coupling matrices.

* Bugfix [commit d0e9cdb]: Correctly set low-energy data (read from
  the SLHA input file) in the low-energy constraint.

* Bugfix [commit 6414e46]: Convert fermion masses and mixing matrices
  to SLHA convention in the SLHA output.

* Install specimen SLHA input files in the model directory when one
  runs the createmodel script.

* Work around fields in ``Cp[]`` carrying an invalid index that cause
  Part::partw when passed to ````SARAH``Vertex[]````.

* Support Cygwin on MS Windows

* New model file for the TMSSM (triplet Higgs model)


FlexibleSUSY-1.0.1 [June, 11 2014]
==================================

* Bugfix [commit 4dc897e]: consts.hpp is not distributed but appears
  in the list of installed headers


FlexibleSUSY-1.0.0 [June, 10 2014]
==================================

* Bugfix [commit 399a1c8]: renaming SoftsusyMSSM and SoftsusyNMSSM
  model classes and files to make ``make all-test`` work on HFS (fixes
  #2).

* Bugfix [commit cfc2562]: correcting MS-bar to DR-bar conversion of
  fermion masses mb and mtau.

* Bugfix [commit ceecc4a]: fixing compilation error with Intel icpc
  14.0, Build 20130728 (and GNU STL 4.6.4 and 4.8.1).

* Bugfix [commit db60205]: fixing linking error of
  test/test_MSSM_NMSSM_linking.x in case LoopTools is used.

* Bugfix [commit 32c3222]: generalizing color summation routine to
  handle single-generation fields and non-fundamental
  representations.

* Bugfix [commit 3fd2699]: Correcting the determination of the number
  of EWSB eqs. in case of CP violating models.

* Bugfix [commit c9cc34f]: Reset fermion phases when ``clear()`` is called.

* Bugfix [commit faa0fb6]: adding boost include directory to
  ``CPPFLAGS`` in the src module.

* Bugfix [commit ac8e38e]: impose EWSB before calculating the
  spectrum.

* Set minimum required SARAH version to 4.0.4, because it implements
  the full two-loop VEV beta functions from arXiv:1310.7629 .

* Add stand-alone examples to illustrate how to use FlexibleSUSY's
  classes and libraries independently of FlexibleSUSY's build system.

* Add tower example to illustrate how to glue multiple models to form
  a stack of effective field theories.

* Add customized-betas example to illustrate how to replace an
  auto-generated C++ component by something of an alternative origin.

* Implement leading two-loop MSSM and NMSSM tadpoles from Slavich
  (used in the EWSB conditions).

* Implement leading two-loop MSSM and NMSSM CP-even and CP-odd Higgs
  self-energy contributions from Slavich.

* Allow to constrain the boundary condition scale via the model file
  variables ``{Low,SUSY,High}ScaleMinimum`` and
  ``{Low,SUSY,High}ScaleMaximum``.

* Allow explicite setting (and disabling) of the Yukawa couplings in
  the constraints.

* Enable/disable multi-threading at the configure level

* lower required g++ version to 4.4.7

* Enable source code export without the meta code via ``make
  install-src``.

* Add FFLite module as a thread-safe alternative to LoopTools

* Create helper function to find the LSP.

* Allow to select beta-function loop order in the SLHA input file.

* Allow disable/enable threshold corrections in the SLHA input file.

* Rename pole mass calculation precision option and set them in the
  model file.


FlexibleSUSY-0.5.3 [January, 21 2014]
=====================================

* Bugfix [commit 44903c]: correcting malformed print out in
  config/list_sarah_model_files.sh in case model files do not exist

* Bugfix [commit 3aae11]: Prevent hard-coding of the running Weinberg
  angle in terms of the gauge couplings

* Bugfix [commit ce4a73]: Generalize calculation of gauge couplings
  at the low-scale (fixes #1)

* Vertices are saved in a file to avoid repeating same calculation.


FlexibleSUSY-0.5.2 [January 14, 2014]
=====================================

* Bugfix [commit 58f8f9]: Convert beta functions which are identical
  zero to the data type of the corresponding parameter.

* Bugfix [commit e5f937]: Correcting check of SARAH patch level
  against minimum required patch level.

* Bugfix [commit e2d43b]: Adapting free phases of fermion fields if
  mass is less than zero.

* Bugfix [commit e777e1]: Converting indices to C convention in
  tree-level EWSB equations.

* Set minimum required SARAH version to 4.0.3, because it includes a
  bug fix in the index structure of the charged Higgs self-energies.

* Allow setting of single matrix/ vector elements in the constraints.

* Model files are now in the directory model_files/
  (instead of templates/)

* The command line arguments of the createmodel script changed.
  Please see ``./createmodel --help`` for more details.

* Add support for the ``LHInput[p]`` command in constraints, which reads
  the parameter ``p`` from the SLHA input file.

* Constrain time used to simplify the beta functions (default: 120
  seconds per beta function).  To change the time constraint, set
  ````FlexibleSUSY``FSSimplifyBetaFunctionsTimeConstraint````.

* Avoid swapping by distributing the calculation of the two-scale
  beta functions among multiple .cpp files.

* Introduce separate meta code stamp (triggers running of the meta
  code) with name ``models/<model-name>/00_DELETE_ME_TO_RERUN_METACODE``


FlexibleSUSY-0.5.1 [November 23, 2013]
======================================

* Handle parameters of type vector in the beta functions.


FlexibleSUSY-0.5 [November 18, 2013]
====================================

* Store particle masses as Eigen::Array and mixing matrices as
  ``Eigen::Matrix``.

.. _Conan: https://conan.io/
.. _Eigen 3: http://eigen.tuxfamily.org
.. _GM2Calc: https://arxiv.org/abs/1510.08071
.. _MhEFT: https://gabrlee.com/code/
.. _FeynArts: http://www.feynarts.de
.. _FormCalc: http://www.feynarts.de/formcalc
.. _LoopTools: http://www.feynarts.de/looptools/
.. _TSIL: https://www.niu.edu/spmartin/tsil/
