Mathematica interface
=====================

The FlexibleSUSY-generated spectrum generators can be called from
within Mathematica using Wolfram's LibraryLink interface.

.. contents:: Table of Contents

Quick start
-----------

The following example calculates the pole mass spectrum and the
observables in the CMSSM for a given parameter point::

    Get["models/CMSSM/CMSSM_librarylink.m"];
    
    (* Create a handle to a model given the input parameters.
       See Options[FSCMSSMOpenHandle] for all default options. *)
    handle = FSCMSSMOpenHandle[
      fsSettings -> { precisionGoal -> 1.*^-4 },
      fsSMParameters -> { Mt -> 173.3 },
      fsModelParameters -> {
          m0 -> 125, m12 -> 500, TanBeta -> 10, SignMu -> 1, Azero -> 0 }
    ];
    
    (* calculate pole mass spectrum *)
    FSCMSSMCalculateSpectrum[handle]
    
    (* calculate further observables *)
    FSCMSSMCalculateObservables[handle]
    
    (* close the model handle *)
    FSCMSSMCloseHandle[handle];

Note: For each model, FlexibleSUSY creates an example Mathematica
script which illustrates the use of the Mathematica interface.  The
generated example can be found in ``models/<model>/run_<model>.m``,
which can be run for example as
::

    math -run "<< \"models/<model>/run_<model>.m\""

Before running the example Mathematica script, the model parameters in
the script should be set to reasonable values.

In the following more advanced example the Higgs pole mass is
calculated in the HSSUSY_ model (an EFT calculation of the SM-like
Higgs mass, assuming that the high-energy completion is the MSSM) as a
function of :math:`X_t / M_S` for :math:`\tan\beta = 5` and different
values of the SUSY scale.  The example also illustrates how
parallelization can be used to exploit the performance of multi-core
systems.
::

    Get["models/HSSUSY/HSSUSY_librarylink.m"];
    
    CalcMh[TB_, Xtt_, MS_] := Module[{handle, spec},
        handle = FSHSSUSYOpenHandle[
            fsSettings -> {
                precisionGoal -> 1.*^-5,
                calculateStandardModelMasses -> 1,
                poleMassLoopOrder -> 4,
                ewsbLoopOrder -> 4,
                betaFunctionLoopOrder -> 4,
                thresholdCorrectionsLoopOrder -> 3,
                thresholdCorrections -> 123111321,
                parameterOutputScale -> 173.34
            },
            fsModelParameters -> {
                TanBeta -> TB,
                MEWSB -> 173.34,
                MSUSY -> MS,
                M1Input -> MS,
                M2Input -> MS,
                M3Input -> MS,
                MuInput -> MS,
                mAInput -> MS,
                AtInput -> (Xtt + 1/TB) * MS,
                AbInput -> 0,
                AtauInput -> 0,
                msq2 -> MS^2 IdentityMatrix[3],
                msu2 -> MS^2 IdentityMatrix[3],
                msd2 -> MS^2 IdentityMatrix[3],
                msl2 -> MS^2 IdentityMatrix[3],
                mse2 -> MS^2 IdentityMatrix[3],
                LambdaLoopOrder -> 2,
                TwoLoopAtAs -> 1,
                TwoLoopAbAs -> 1,
                TwoLoopAtAb -> 1,
                TwoLoopAtauAtau -> 1,
                TwoLoopAtAt -> 1
            }
        ];
        spec = HSSUSY /. FSHSSUSYCalculateSpectrum[handle];
        FSHSSUSYCloseHandle[handle];
        If[spec =!= $Failed, Pole[M[hh]] /. spec, 0]
    ];
    
    LaunchKernels[];
    DistributeDefinitions[CalcMh];
    
    data = {
        ParallelMap[{#, CalcMh[5, #, 1000 ]}&, Range[-3.5, 3.5, 0.1]],
        ParallelMap[{#, CalcMh[5, #, 2000 ]}&, Range[-3.5, 3.5, 0.1]],
        ParallelMap[{#, CalcMh[5, #, 10000]}&, Range[-3.5, 3.5, 0.1]]
    };
    
    plot = ListPlot[data,
                    PlotLegends -> {"MS = 1 TeV", "MS = 2 TeV", "MS = 10 TeV"},
                    Axes -> False, Frame -> True,
                    FrameLabel -> {"Xt / MS", "Mh / GeV"}];
    
    Export["HSSUSY_Mh_Xt.png", plot, ImageSize -> 1000];


Output:

.. image:: images/HSSUSY_Mh_Xt.png
   :align: center

Using the Mathematica interface functions
-----------------------------------------

Building the LibraryLink library
````````````````````````````````

In order to build the LibraryLink library, FlexibleSUSY must be
configured with ``--enable-meta`` (enabled by default).

Example::

    ./configure --with-models=CMSSM
    make

The LibraryLink library can be found in
``models/<model>/<model>_librarylink.so``, where ``<model>`` is the model
name.  In order to use FlexibleSUSY's generated ``<model>`` spectrum
generator at the Mathematica level, the library functions must be
loaded using the ``models/<model>/<model>_librarylink.m`` script.

Example::

    Get["models/CMSSM/CMSSM_librarylink.m"];

FS<model>OpenHandle
```````````````````

First, a handle to the model must be created using the
``FS<model>OpenHandle[]`` function.  The function takes as arguments

- the spectrum generator settings via the ``fsSettings`` variable
- the Standard Model input parameters via the ``fsSMParameters`` variable
- the model input parameters via the ``fsModelParameters`` variable

Example::

    Get["models/CMSSM/CMSSM_librarylink.m"];
    handle = FSCMSSMOpenHandle[
      fsSettings -> { precisionGoal -> 1.*^-4 },
      fsSMParameters -> { Mt -> 173.3 },
      fsModelParameters -> {
          m0 -> 125, m12 -> 500, TanBeta -> 10, SignMu -> 1, Azero -> 0 }
    ];
    FSCMSSMGetSettings[handle]
    FSCMSSMGetSMInputParameters[handle]
    FSCMSSMGetInputParameters[handle]

The ``FS<model>OpenHandle[]`` fixes all settings and input parameters at
once.  Unspecified parameters are set to their default values.  The
default values are stored in the variables ``fsDefaultSettings``,
``fsDefaultSMParameters`` and ``fs<model>DefaultInputParameters``::

    Get["models/CMSSM/CMSSM_librarylink.m"];
    Print[fsDefaultSettings];
    Print[fsDefaultSMParameters];
    Print[fsCMSSMDefaultInputParameters];

The settings associated to a ``handle`` can be listed using the
``FS<model>GetSettings[]`` function.  Please refer to the
`FlexibleSUSY run-time configuration`_ for more information on the
spectrum generator settings.

Example::

    Get["models/CMSSM/CMSSM_librarylink.m"];
    handle = FSCMSSMOpenHandle[
      fsSettings -> { precisionGoal -> 1.*^-5, betaFunctionLoopOrder -> 3 }
    ];
    FSCMSSMGetSettings[handle]

Output::

    { precisionGoal -> 0.00001,
      maxIterations -> 0,
      calculateStandardModelMasses -> 0,
      poleMassLoopOrder -> 2,
      ewsbLoopOrder -> 2,
      betaFunctionLoopOrder -> 3,
      thresholdCorrectionsLoopOrder -> 2,
      higgs2loopCorrectionAtAs -> 1,
      higgs2loopCorrectionAbAs -> 1,
      higgs2loopCorrectionAtAt -> 1,
      higgs2loopCorrectionAtauAtau -> 1,
      forceOutput -> 0,
      top2loopCorrectionsQCD -> 1,
      betaZeroThreshold -> 1.*10^-11,
      forcePositiveMasses -> 0,
      poleMassScale -> 0.,
      parameterOutputScale -> 0. }

The Standard Model input parameters associated to a ``handle`` can be
listed using the ``FS<model>GetSMInputParameters[]`` function.

Example::

    Get["models/CMSSM/CMSSM_librarylink.m"];
    handle = FSCMSSMOpenHandle[
      fsSMParameters -> { Mt -> 173.34 }
    ];
    FSCMSSMGetSMInputParameters[handle]

Output::

    { alphaEmMZ -> 0.00781763, (* alpha_em(MZ) in the SM(5), MS-bar *)
      GF -> 0.000011663787,    (* Fermi constant *)
      alphaSMZ -> 0.1184,      (* alpha_s(MZ) in the SM(5), MS-bar *)
      MZ -> 91.1876,           (* Z pole mass *)
      mbmb -> 4.18,            (* MS-bar bottom mass at Q = mb *)
      Mt -> 173.34,            (* top pole mass *)
      Mtau -> 1.777,           (* tau pole mass *)
      Mv3 -> 0.,               (* 3rd heaviest neutrino mass *)
      MW -> 80.385,            (* W pole mass *)
      Me -> 0.000510999,       (* electron pole mass *)
      Mv1 -> 0.,               (* 1st neutrino mass *)
      Mm -> 0.105658,          (* muon pole masss *)
      Mv2 -> 0.,               (* 2nd neutrino mass *)
      md2GeV -> 0.00475,       (* MS-bar down quark mass at Q = 2 GeV *)
      mu2GeV -> 0.0024,        (* MS-bar up quark mass at Q = 2 GeV *)
      ms2GeV -> 0.104,         (* MS-bar strange quark mass at Q = 2 GeV *)
      mcmc -> 1.27,            (* MS-bar charm quark mass at Q = mc *)
      alphaEm0 -> 0.00729735,  (* alpha_em in the Thompson limit *)
      Mh -> 125.09 }           (* Higgs pole mass *)

The model input parameters associated to a ``handle`` can be listed
using the ``FS<model>GetInputParameters[]`` function.

Example::

    Get["models/CMSSM/CMSSM_librarylink.m"];
    handle = FSCMSSMOpenHandle[
      fsModelParameters -> { m0 -> 125, m12 -> 500, TanBeta -> 10, SignMu -> 1 }
    ];
    FSCMSSMGetInputParameters[handle]

Output::

    { m0 -> 125.,
      m12 -> 500.,
      TanBeta -> 10.,
      SignMu -> 1,
      Azero -> 0. }

FS<model>Set
````````````

Using the ``FS<model>Set[]`` function, the input parameters and settings
associated to a certain handle can be modified.  The ``FS<model>Set[]``
function takes first as argument the handle, and as second argument
the replacement list of new parameters / settings.

Example::

    Get["models/CMSSM/CMSSM_librarylink.m"];
    handle = FSCMSSMOpenHandle[
      fsSettings -> { precisionGoal -> 1.*^-4 },
      fsSMParameters -> { Mt -> 173.3 },
      fsModelParameters -> {
          m0 -> 125, m12 -> 500, TanBeta -> 10, SignMu -> 1, Azero -> 0 }
    ];
    
    FSCMSSMGetInputParameters[handle]
    
    FSCMSSMSet[handle, TanBeta -> 20];
    
    FSCMSSMGetInputParameters[handle]

Output::

    {m0 -> 125., m12 -> 500., TanBeta -> 10., SignMu -> 1, Azero -> 0.}
    
    {m0 -> 125., m12 -> 500., TanBeta -> 20., SignMu -> 1, Azero -> 0.}

FS<model>CalculateSpectrum
``````````````````````````

For each ``<model>``, the ``FS<model>CalculateSpectrum[handle]`` function
solves the boundary value problem and calculates the pole mass
spectrum.  The function takes a model handle as arguments, referring
to the settings and input parameters

The function returns all running model parameters at the parameter
output scale (either the SUSY scale or the scale set via ``fsSettings
-> { parameterOutputScale -> 1000. }``) and the running masses at the
same scale.  The running masses are denoted by ``M[p]`` where ``p`` is the
particle name.  The parameter output scale appears in the returned
list with the symbol ``SCALE``.  The calculated pole masses are denoted
by ``Pole[M[p]]``, respectively.  The mixing matrices which correspond
to the pole masses are denoted by ``Pole[Z]``, where Z is the name of
the mixing matrix.
::

    Get["models/CMSSM/CMSSM_librarylink.m"];
    handle = FSCMSSMOpenHandle[
      fsModelParameters -> { m0 -> 125, m12 -> 500, TanBeta -> 10, SignMu -> 1 }
    ];
    FSCMSSMCalculateSpectrum[handle]

Output::

    {CMSSM ->
       {M[VG] -> 0., M[Glu] -> 1117.18, M[Fv] -> {0., 0., 0.},
        M[Sd] -> {942.251, 977.989, 980.297, 980.3, 1023.94, 1023.94},
        M[Sv] -> {347.371, 348.42, 348.424},
        M[Su] -> {782.7, 983.889, 983.894, 987.561, 1021., 1021.},
        M[Se] -> {219.073, 226.223, 226.248, 356.971, 356.976, 358.335},
        M[hh] -> {88.1593, 732.573}, M[Ah] -> {90.0927, 732.337},
        M[Hpm] -> {78.4808, 736.531},
        M[Chi] -> {207.439, 376.528, 633.944, 647.755},
        M[Cha] -> {376.365, 647.464},
        M[Fe] -> {0.000520523, 0.107628, 1.81042},
        M[Fd] -> {0.00243143, 0.0532355, 2.32379},
        M[Fu] -> {0.00122119, 0.549091, 147.438}, M[VWm] -> 78.4808,
        M[VP] -> 0., M[VZ] -> 90.0927,
        ZD -> {{0., 0., -0.965619, 0., 0., -0.259961}, {0., 0., 0.259961, 0.,
            0., -0.965619}, {0., -0.00456672, 0., 0., -0.99999,
           0.}, {0.000208583, 0., 0., 1., 0., 0.}, {0., -0.99999, 0., 0.,
           0.00456672, 0.}, {1., 0., 0., -0.000208583, 0., 0.}},
        ZV -> {{0., 0., 1.}, {0., 1., 0.}, {1., 0., 0.}},
        ZU -> {{0., 0., 0.430138, 0., 0., 0.902763}, {0., 0.00896415, 0., 0.,
            0.99996, 0.}, {0.000019939, 0., 0., 1., 0., 0.}, {0., 0.,
           0.902763, 0., 0., -0.430138}, {1., 0., 0., -0.000019939, 0.,
           0.}, {0., 0.99996, 0., 0., -0.00896415, 0.}},
        ZE -> {{0., 0., 0.145606, 0., 0., 0.989343}, {0., -0.00903329, 0.,
           0., -0.999959, 0.}, {0.0000436949, 0., 0., 1., 0., 0.}, {1., 0.,
           0., -0.0000436949, 0., 0.}, {0., -0.999959, 0., 0., 0.00903329,
           0.}, {0., 0., 0.989343, 0., 0., -0.145606}},
        ZH -> {{0.105881, 0.994379}, {0.994379, -0.105881}},
        ZA -> {{-0.102825, 0.994699}, {0.994699, 0.102825}},
        ZP -> {{-0.102825, 0.994699}, {0.994699, 0.102825}},
        ZN -> {{-0.995744, 0.018728, -0.0832596, 0.0348113}, {0.0389752,
           0.971833, -0.194009, 0.127995}, {0. - 0.0331609 I,
           0. + 0.0485202 I, 0. + 0.703592 I,
           0. + 0.70817 I}, {0.0766551, -0.229862, -0.678518, 0.69347}},
        UM -> {{0.960661, -0.277725}, {0.277725, 0.960661}},
        UP -> {{0.983012, -0.183543}, {0.183543, 0.983012}},
        ZEL -> {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}},
        ZER -> {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}},
        ZDL -> {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}},
        ZDR -> {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}},
        ZUL -> {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}},
        ZUR -> {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}},
        ZZ -> {{-0.871112, 0.491084}, {-0.491084, -0.871112}},
        Pole[M[VG]] -> 0., Pole[M[Glu]] -> 1151.38,
        Pole[M[Fv]] -> {0., 0., 0.},
        Pole[M[Sd]] -> {970.999, 1012.32, 1015.42, 1015.42, 1059.73,
          1059.73}, Pole[M[Sv]] -> {351.491, 352.69, 352.694},
        Pole[M[Su]] -> {809.283, 1015.61, 1018.71, 1019.46, 1056.91,
          1056.91},
        Pole[M[Se]] -> {222.482, 229.821, 229.847, 361.599, 361.604,
          362.781}, Pole[M[hh]] -> {114.781, 719.259},
        Pole[M[Ah]] -> {88.5742, 718.986},
        Pole[M[Hpm]] -> {77.7605, 723.723},
        Pole[M[Chi]] -> {204.267, 385.936, 636.143, 649.77},
        Pole[M[Cha]] -> {385.949, 650.096}, Pole[M[Fe]] -> {0., 0., 0.},
        Pole[M[Fd]] -> {0., 0., 0.}, Pole[M[Fu]] -> {0., 0., 0.},
        Pole[M[VWm]] -> 80.3924, Pole[M[VP]] -> 0., Pole[M[VZ]] -> 0.,
        Pole[ZD] -> {{0., 0., -0.977566, 0., 0., -0.210631}, {0., 0.,
           0.210631, 0., 0., -0.977566}, {0., -0.0045424, 0., 0., -0.99999,
           0.}, {0.000207472, 0., 0., 1., 0., 0.}, {0., -0.99999, 0., 0.,
           0.0045424, 0.}, {1., 0., 0., -0.000207472, 0., 0.}},
        Pole[ZV] -> {{0., 0., 1.}, {0., 1., 0.}, {1., 0., 0.}},
        Pole[ZU] -> {{0., 0., 0.427999, 0., 0., 0.903779}, {0., 0., 0.903779,
            0., 0., -0.427999}, {0., 0.00911132, 0., 0., 0.999958,
           0.}, {0.0000202664, 0., 0., 1., 0., 0.}, {1., 0.,
           0., -0.0000202664, 0., 0.}, {0., 0.999958, 0., 0., -0.00911132,
           0.}}, Pole[
          ZE] -> {{0., 0., 0.144271, 0., 3.02431*10^-15,
           0.989538}, {0., -0.00895024, 2.08714*10^-14, 0., -0.99996,
           0.}, {0.0000432932, 0., 0., 1., 0., 0.}, {1., 0.,
           0., -0.0000432932, 0., 0.}, {0., -0.99996, -1.86811*10^-16, 0.,
           0.00895024, 0.}, {0., 0., -0.989538, 0., -2.06711*10^-14,
           0.144271}},
        Pole[ZH] -> {{0.106581, 0.994304}, {0.994304, -0.106581}},
        Pole[ZA] -> {{-0.0989827, 0.995089}, {0.995089, 0.0989827}},
        Pole[ZP] -> {{-0.0995943, 0.995028}, {0.995028, 0.0995943}},
        Pole[ZN] -> {{-0.995819, 0.0174686, -0.082821,
           0.0343646}, {0.0380335, 0.970567, -0.197841,
           0.131955}, {0. - 0.0332126 I, 0. + 0.0483916 I, 0. + 0.703447 I,
           0. + 0.70832 I}, {0.0761299, -0.235272, -0.677615, 0.692596}},
        Pole[UM] -> {{0.95912, -0.283001}, {0.283001, 0.95912}},
        Pole[UP] -> {{0.981917, -0.189314}, {0.189314, 0.981917}},
        Pole[ZEL] -> {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}},
        Pole[ZER] -> {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}},
        Pole[ZDL] -> {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}},
        Pole[ZDR] -> {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}},
        Pole[ZUL] -> {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}},
        Pole[ZUR] -> {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}},
        Pole[ZZ] -> {{0., 0.}, {0., 0.}},
        Yd -> {{0.000136987, 0., 0.}, {0., 0.0029993, 0.}, {0., 0.,
           0.130923}},
        Ye -> {{0.0000293264, 0., 0.}, {0., 0.00606377, 0.}, {0., 0.,
           0.102}},
        Yu -> {{7.1123*10^-6, 0., 0.}, {0., 0.00319794, 0.}, {0., 0.,
           0.858685}}, \[Mu] -> 630.611, g1 -> 0.467953, g2 -> 0.642978,
        g3 -> 1.06483, vd -> 25.1013, vu -> 242.823,
        T[Yd] -> {{-0.192259, 0., 0.}, {0., -4.20945, 0.}, {0.,
           0., -171.869}},
        T[Ye] -> {{-0.00878455, 0., 0.}, {0., -1.81633, 0.}, {0.,
           0., -30.3818}},
        T[Yu] -> {{-0.00817412, 0., 0.}, {0., -3.67535, 0.}, {0.,
           0., -764.191}}, B[\[Mu]] -> 54854.6,
        mq2 -> {{1.04513*10^6, 0., 0.}, {0., 1.04512*10^6, 0.}, {0., 0.,
           889135.}},
        ml2 -> {{125372., 0., 0.}, {0., 125369., 0.}, {0., 0., 124639.}},
        mHd2 -> 109915., mHu2 -> -385101.,
        md2 -> {{960350., 0., 0.}, {0., 960345., 0.}, {0., 0., 951180.}},
        mu2 -> {{969326., 0., 0.}, {0., 969321., 0.}, {0., 0., 659257.}},
        me2 -> {{49272.2, 0., 0.}, {0., 49266.9, 0.}, {0., 0., 47778.5}},
        MassB -> 209.358, MassWB -> 388.421, MassG -> 1117.18,
        SCALE -> 879.186}
    }

FS<model>CalculateObservables
`````````````````````````````

For each ``<model>``, the ``FS<model>CalculateObservables[handle]``
function calculates further observables, such as effective Higgs
couplings to two photons or gluons.  See the section on Observables_
for a list of all available observables.

Note: The ``FS<model>CalculateObservables[handle]`` function assumes,
that the pole mass spectrum has been calculated before, using the
``FS<model>CalculateSpectrum[handle]`` function.
::

    Get["models/CMSSM/CMSSM_librarylink.m"];
    handle = FSCMSSMOpenHandle[
      fsModelParameters -> { m0 -> 125, m12 -> 500, TanBeta -> 10, SignMu -> 1 }
    ];
    FSCMSSMCalculateSpectrum[handle]
    FSCMSSMCalculateObservables[handle]

Output::

    {CMSSM ->
       { FlexibleSUSYObservable``CpHiggsPhotonPhoton ->
           {0.0000296409 - 2.1245*10^-7 I, 7.82123*10^-7 + 9.1076*10^-7 I},
         FlexibleSUSYObservable``CpHiggsGluonGluon ->
           {-0.0000670724 - 2.65658*10^-6 I, 2.72135*10^-6 + 4.91993*10^-6 I},
         FlexibleSUSYObservable``CpPseudoScalarPhotonPhoton ->
           1.05105*10^-6 - 8.33068*10^-7 I,
         FlexibleSUSYObservable``CpPseudoScalarGluonGluon ->
           6.71448*10^-6 + 8.41625*10^-7 I }
    }

FS<model>GetProblems and FS<model>GetWarnings
`````````````````````````````````````````````

After the spectrum has been calculated, one should check for problems
or warnings.  They can be obtained for a given handle using the
``FS<model>GetProblems[handle]`` and ``FS<model>GetWarnings[handle]``
functions, respectively.  These functions return the empty list if no
problems / warnings occurred.
::

    Get["models/CMSSM/CMSSM_librarylink.m"];
    handle = FSCMSSMOpenHandle[
      fsModelParameters -> { m0 -> 1000, m12 -> 500, Azero -> -10000, TanBeta -> 2, SignMu -> 1 }
    ];
    FSCMSSMCalculateSpectrum[handle];
    FSCMSSMGetProblems[handle]

Output::

    {CMSSM ->
      { Tachyons -> {M[Sd], M[Su]},
        NoPoleMassConvergence -> {Pole[M[hh]]} }
    }

This list of problems states, that the running up-type and down-type
squarks are tachyonic for this parameter point.  Thus, the spectrum
calculated by FlexibleSUSY for this point cannot be trusted.
Furthermore, the iteration to determine the Higgs pole mass did not
converge.  Thus, the calculated Higgs pole mass cannot be trusted
either for this parameter point.

FS<model>ToSLHA
```````````````

The running parameters, the mass spectrum and/or the observables can
be converted to SLHA format using the ``FS<model>ToSLHA[handle]``
function.  The function returns a string formatted according to
[SLHA1_, SLHA2_].

.. _SLHA1: https://inspirehep.net/record/632863
.. _SLHA2: https://inspirehep.net/record/777216

Example::

    Get["models/CMSSM/CMSSM_librarylink.m"];
    handle = FSCMSSMOpenHandle[
      fsModelParameters -> { m0 -> 1000, m12 -> 500, Azero -> 0, TanBeta -> 10, SignMu -> 1 }
    ];
    FSCMSSMCalculateSpectrum[handle];
    FSCMSSMCalculateObservables[handle];
    Export["spectrum.slha", FSCMSSMToSLHA[handle], "String"];

Output: ``spectrum.slha``
::

    Block SPINFO
         1   FlexibleSUSY
         2   1.7.1
         5   CMSSM
         9   4.9.1
    Block FlexibleSUSY
         0     1.00000000E-04   # precision goal
         1     0.00000000E+00   # max. iterations (0 = automatic)
         2     0.00000000E+00   # algorithm (0 = two_scale)
         3     0.00000000E+00   # calculate SM pole masses
         4     2.00000000E+00   # pole mass loop order
         5     2.00000000E+00   # EWSB loop order
         6     3.00000000E+00   # beta-functions loop order
         7     2.00000000E+00   # threshold corrections loop order
         8     1.00000000E+00   # Higgs 2-loop corrections O(alpha_t alpha_s)
         9     1.00000000E+00   # Higgs 2-loop corrections O(alpha_b alpha_s)
        10     1.00000000E+00   # Higgs 2-loop corrections O((alpha_t + alpha_b)^2)
        11     1.00000000E+00   # Higgs 2-loop corrections O(alpha_tau^2)
        12     0.00000000E+00   # force output
        13     1.00000000E+00   # Top pole mass QCD corrections (0 = 1L, 1 = 2L, 2 = 3L)
        14     1.00000000E-11   # beta-function zero threshold
        15     0.00000000E+00   # calculate observables (a_muon, ...)
        16     0.00000000E+00   # force positive majorana masses
        17     0.00000000E+00   # pole mass renormalization scale (0 = SUSY scale)
        18     0.00000000E+00   # pole mass renormalization scale in the EFT (0 = min(SUSY scale, Mt))
        19     0.00000000E+00   # EFT matching scale (0 = SUSY scale)
        20     2.00000000E+00   # EFT loop order for upwards matching
        21     1.00000000E+00   # EFT loop order for downwards matching
        22     0.00000000E+00   # EFT index of SM-like Higgs in the BSM model
        23     1.00000000E+00   # calculate BSM pole masses
    Block SMINPUTS
         1     1.27916000E+02   # alpha^(-1) SM MSbar(MZ)
         2     1.16637000E-05   # G_Fermi
         3     1.18400000E-01   # alpha_s(MZ) SM MSbar
         4     9.11876000E+01   # MZ(pole)
         5     4.18000000E+00   # mb(mb) SM MSbar
         6     1.73340000E+02   # mtop(pole)
         7     1.77700000E+00   # mtau(pole)
         8     0.00000000E+00   # mnu3(pole)
         9     8.03850000E+01   # MW(pole)
        11     5.10998902E-04   # melectron(pole)
        12     0.00000000E+00   # mnu1(pole)
        13     1.05658372E-01   # mmuon(pole)
        14     0.00000000E+00   # mnu2(pole)
        21     4.75000000E-03   # md
        22     2.40000000E-03   # mu
        23     1.04000000E-01   # ms
        24     1.27000000E+00   # mc
    Block FlexibleSUSYInput
         0     7.29735257E-03   # alpha_em(0)
         1     1.25090000E+02   # mh_pole
    Block MODSEL
         6                  0   # quark/lepton flavour violation
        12     0.00000000E+00   # DRbar parameter output scale (GeV)
    Block MINPAR
         1     1.00000000E+03   # m0
         2     5.00000000E+02   # m12
         3     1.00000000E+01   # TanBeta
         4                  1   # SignMu
         5     0.00000000E+00   # Azero
    Block gauge Q= 1.08941472E+03
         1     3.62448909E-01   # gY
         2     6.41812871E-01   # g2
         3     1.05659964E+00   # g3
    Block Yu Q= 1.08941472E+03
      1  1     7.21657623E-06   # Yu(1,1)
      1  2     0.00000000E+00   # Yu(1,2)
      1  3     0.00000000E+00   # Yu(1,3)
      2  1     0.00000000E+00   # Yu(2,1)
      2  2     3.28521207E-03   # Yu(2,2)
      2  3     0.00000000E+00   # Yu(2,3)
      3  1     0.00000000E+00   # Yu(3,1)
      3  2     0.00000000E+00   # Yu(3,2)
      3  3     8.57291451E-01   # Yu(3,3)
    Block Yd Q= 1.08941472E+03
      1  1     1.38593672E-04   # Yd(1,1)
      1  2     0.00000000E+00   # Yd(1,2)
      1  3     0.00000000E+00   # Yd(1,3)
      2  1     0.00000000E+00   # Yd(2,1)
      2  2     3.03447401E-03   # Yd(2,2)
      2  3     0.00000000E+00   # Yd(2,3)
      3  1     0.00000000E+00   # Yd(3,1)
      3  2     0.00000000E+00   # Yd(3,2)
      3  3     1.31773557E-01   # Yd(3,3)
    Block Ye Q= 1.08941472E+03
      1  1     2.86516322E-05   # Ye(1,1)
      1  2     0.00000000E+00   # Ye(1,2)
      1  3     0.00000000E+00   # Ye(1,3)
      2  1     0.00000000E+00   # Ye(2,1)
      2  2     5.92424993E-03   # Ye(2,2)
      2  3     0.00000000E+00   # Ye(2,3)
      3  1     0.00000000E+00   # Ye(3,1)
      3  2     0.00000000E+00   # Ye(3,2)
      3  3     9.96402101E-02   # Ye(3,3)
    Block Te Q= 1.08941472E+03
      1  1    -8.51105194E-03   # TYe(1,1)
      1  2     0.00000000E+00   # TYe(1,2)
      1  3     0.00000000E+00   # TYe(1,3)
      2  1     0.00000000E+00   # TYe(2,1)
      2  2    -1.75978274E+00   # TYe(2,2)
      2  3     0.00000000E+00   # TYe(2,3)
      3  1     0.00000000E+00   # TYe(3,1)
      3  2     0.00000000E+00   # TYe(3,2)
      3  3    -2.94405366E+01   # TYe(3,3)
    Block Td Q= 1.08941472E+03
      1  1    -1.88977041E-01   # TYd(1,1)
      1  2     0.00000000E+00   # TYd(1,2)
      1  3     0.00000000E+00   # TYd(1,3)
      2  1     0.00000000E+00   # TYd(2,1)
      2  2    -4.13759191E+00   # TYd(2,2)
      2  3     0.00000000E+00   # TYd(2,3)
      3  1     0.00000000E+00   # TYd(3,1)
      3  2     0.00000000E+00   # TYd(3,2)
      3  3    -1.68021836E+02   # TYd(3,3)
    Block Tu Q= 1.08941472E+03
      1  1    -8.04938899E-03   # TYu(1,1)
      1  2     0.00000000E+00   # TYu(1,2)
      1  3     0.00000000E+00   # TYu(1,3)
      2  1     0.00000000E+00   # TYu(2,1)
      2  2    -3.66431893E+00   # TYu(2,2)
      2  3     0.00000000E+00   # TYu(2,3)
      3  1     0.00000000E+00   # TYu(3,1)
      3  2     0.00000000E+00   # TYu(3,2)
      3  3    -7.39719203E+02   # TYu(3,3)
    Block MSQ2 Q= 1.08941472E+03
      1  1     1.92793919E+06   # mq2(1,1)
      1  2     0.00000000E+00   # mq2(1,2)
      1  3     0.00000000E+00   # mq2(1,3)
      2  1     0.00000000E+00   # mq2(2,1)
      2  2     1.92792577E+06   # mq2(2,2)
      2  3     0.00000000E+00   # mq2(2,3)
      3  1     0.00000000E+00   # mq2(3,1)
      3  2     0.00000000E+00   # mq2(3,2)
      3  3     1.46411510E+06   # mq2(3,3)
    Block MSE2 Q= 1.08941472E+03
      1  1     1.02965849E+06   # me2(1,1)
      1  2     0.00000000E+00   # me2(1,2)
      1  3     0.00000000E+00   # me2(1,3)
      2  1     0.00000000E+00   # me2(2,1)
      2  2     1.02959546E+06   # me2(2,2)
      2  3     0.00000000E+00   # me2(2,3)
      3  1     0.00000000E+00   # me2(3,1)
      3  2     0.00000000E+00   # me2(3,2)
      3  3     1.01183185E+06   # me2(3,3)
    Block MSL2 Q= 1.08941472E+03
      1  1     1.09811858E+06   # ml2(1,1)
      1  2     0.00000000E+00   # ml2(1,2)
      1  3     0.00000000E+00   # ml2(1,3)
      2  1     0.00000000E+00   # ml2(2,1)
      2  2     1.09808726E+06   # ml2(2,2)
      2  3     0.00000000E+00   # ml2(2,3)
      3  1     0.00000000E+00   # ml2(3,1)
      3  2     0.00000000E+00   # ml2(3,2)
      3  3     1.08926184E+06   # ml2(3,3)
    Block MSU2 Q= 1.08941472E+03
      1  1     1.85944942E+06   # mu2(1,1)
      1  2     0.00000000E+00   # mu2(1,2)
      1  3     0.00000000E+00   # mu2(1,3)
      2  1     0.00000000E+00   # mu2(2,1)
      2  2     1.85943524E+06   # mu2(2,2)
      2  3     0.00000000E+00   # mu2(2,3)
      3  1     0.00000000E+00   # mu2(3,1)
      3  2     0.00000000E+00   # mu2(3,2)
      3  3     9.41968516E+05   # mu2(3,3)
    Block MSD2 Q= 1.08941472E+03
      1  1     1.85122889E+06   # md2(1,1)
      1  2     0.00000000E+00   # md2(1,2)
      1  3     0.00000000E+00   # md2(1,3)
      2  1     0.00000000E+00   # md2(2,1)
      2  2     1.85121598E+06   # md2(2,2)
      2  3     0.00000000E+00   # md2(2,3)
      3  1     0.00000000E+00   # md2(3,1)
      3  2     0.00000000E+00   # md2(3,2)
      3  3     1.82729635E+06   # md2(3,3)
    Block Phases Q= 1.08941472E+03
         1     1.00000000E+00   # Re(PhaseGlu)
    Block IMPhases Q= 1.08941472E+03
         1     0.00000000E+00   # Im(PhaseGlu)
    Block MASS
       1000021     1.19858229E+03   # Glu
            24     8.03923382E+01   # VWm
       1000024     3.92101428E+02   # Cha(1)
       1000037     6.38758479E+02   # Cha(2)
            25     1.15429842E+02   # hh(1)
            35     1.20308783E+03   # hh(2)
            37     1.20580267E+03   # Hpm(2)
            36     1.20306748E+03   # Ah(2)
       1000012     1.04330120E+03   # Sv(1)
       1000014     1.04765714E+03   # Sv(2)
       1000016     1.04767257E+03   # Sv(3)
       1000022     2.07437446E+02   # Chi(1)
       1000023     3.92113022E+02   # Chi(2)
       1000025    -6.23279987E+02   # Chi(3)
       1000035     6.38634777E+02   # Chi(4)
       1000001     1.23941353E+03   # Sd(1)
       1000003     1.38335815E+03   # Sd(2)
       1000005     1.39207702E+03   # Sd(3)
       2000001     1.39208249E+03   # Sd(4)
       2000003     1.42196160E+03   # Sd(5)
       2000005     1.42196587E+03   # Sd(6)
       1000011     1.00717444E+03   # Se(1)
       1000013     1.01697963E+03   # Se(2)
       1000015     1.01701468E+03   # Se(3)
       2000011     1.04734697E+03   # Se(4)
       2000013     1.05094581E+03   # Se(5)
       2000015     1.05095815E+03   # Se(6)
       1000002     9.94140079E+02   # Su(1)
       1000004     1.26241854E+03   # Su(2)
       1000006     1.39428015E+03   # Su(3)
       2000002     1.39428785E+03   # Su(4)
       2000004     1.41992874E+03   # Su(5)
       2000006     1.41993106E+03   # Su(6)
    Block UMIX
      1  1     9.53292717E-01   # Re(UM(1,1))
      1  2    -3.02048004E-01   # Re(UM(1,2))
      2  1     3.02048004E-01   # Re(UM(2,1))
      2  2     9.53292717E-01   # Re(UM(2,2))
    Block VMIX
      1  1     9.78125635E-01   # Re(UP(1,1))
      1  2    -2.08015004E-01   # Re(UP(1,2))
      2  1     2.08015004E-01   # Re(UP(2,1))
      2  2     9.78125635E-01   # Re(UP(2,2))
    Block PSEUDOSCALARMIX
      1  1    -1.00061419E-01   # ZA(1,1)
      1  2     9.94981262E-01   # ZA(1,2)
      2  1     9.94981262E-01   # ZA(2,1)
      2  2     1.00061419E-01   # ZA(2,2)
    Block DSQMIX
      1  1    -0.00000000E+00   # ZD(1,1)
      1  2    -0.00000000E+00   # ZD(1,2)
      1  3    -9.99033786E-01   # ZD(1,3)
      1  4    -0.00000000E+00   # ZD(1,4)
      1  5    -0.00000000E+00   # ZD(1,5)
      1  6    -4.39487765E-02   # ZD(1,6)
      2  1     0.00000000E+00   # ZD(2,1)
      2  2     0.00000000E+00   # ZD(2,2)
      2  3     4.39487765E-02   # ZD(2,3)
      2  4     0.00000000E+00   # ZD(2,4)
      2  5     0.00000000E+00   # ZD(2,5)
      2  6    -9.99033786E-01   # ZD(2,6)
      3  1     0.00000000E+00   # ZD(3,1)
      3  2     5.00403353E-03   # ZD(3,2)
      3  3     0.00000000E+00   # ZD(3,3)
      3  4     0.00000000E+00   # ZD(3,4)
      3  5     9.99987480E-01   # ZD(3,5)
      3  6     0.00000000E+00   # ZD(3,6)
      4  1     2.28556802E-04   # ZD(4,1)
      4  2     0.00000000E+00   # ZD(4,2)
      4  3     0.00000000E+00   # ZD(4,3)
      4  4     9.99999974E-01   # ZD(4,4)
      4  5     0.00000000E+00   # ZD(4,5)
      4  6     0.00000000E+00   # ZD(4,6)
      5  1     0.00000000E+00   # ZD(5,1)
      5  2     9.99987480E-01   # ZD(5,2)
      5  3     0.00000000E+00   # ZD(5,3)
      5  4     0.00000000E+00   # ZD(5,4)
      5  5    -5.00403353E-03   # ZD(5,5)
      5  6     0.00000000E+00   # ZD(5,6)
      6  1     9.99999974E-01   # ZD(6,1)
      6  2     0.00000000E+00   # ZD(6,2)
      6  3     0.00000000E+00   # ZD(6,3)
      6  4    -2.28556802E-04   # ZD(6,4)
      6  5     0.00000000E+00   # ZD(6,5)
      6  6     0.00000000E+00   # ZD(6,6)
    Block SELMIX
      1  1     0.00000000E+00   # ZE(1,1)
      1  2     0.00000000E+00   # ZE(1,2)
      1  3     1.38267119E-01   # ZE(1,3)
      1  4     0.00000000E+00   # ZE(1,4)
      1  5     0.00000000E+00   # ZE(1,5)
      1  6     9.90394974E-01   # ZE(1,6)
      2  1     0.00000000E+00   # ZE(2,1)
      2  2    -9.59493336E-03   # ZE(2,2)
      2  3     0.00000000E+00   # ZE(2,3)
      2  4     0.00000000E+00   # ZE(2,4)
      2  5    -9.99953968E-01   # ZE(2,5)
      2  6     0.00000000E+00   # ZE(2,6)
      3  1     4.64326774E-05   # ZE(3,1)
      3  2     0.00000000E+00   # ZE(3,2)
      3  3     0.00000000E+00   # ZE(3,3)
      3  4     9.99999999E-01   # ZE(3,4)
      3  5     0.00000000E+00   # ZE(3,5)
      3  6     0.00000000E+00   # ZE(3,6)
      4  1     0.00000000E+00   # ZE(4,1)
      4  2     0.00000000E+00   # ZE(4,2)
      4  3     9.90394974E-01   # ZE(4,3)
      4  4     0.00000000E+00   # ZE(4,4)
      4  5     0.00000000E+00   # ZE(4,5)
      4  6    -1.38267119E-01   # ZE(4,6)
      5  1     0.00000000E+00   # ZE(5,1)
      5  2    -9.99953968E-01   # ZE(5,2)
      5  3     0.00000000E+00   # ZE(5,3)
      5  4     0.00000000E+00   # ZE(5,4)
      5  5     9.59493336E-03   # ZE(5,5)
      5  6     0.00000000E+00   # ZE(5,6)
      6  1     9.99999999E-01   # ZE(6,1)
      6  2     0.00000000E+00   # ZE(6,2)
      6  3     0.00000000E+00   # ZE(6,3)
      6  4    -4.64326774E-05   # ZE(6,4)
      6  5     0.00000000E+00   # ZE(6,5)
      6  6     0.00000000E+00   # ZE(6,6)
    Block SCALARMIX
      1  1     1.04543664E-01   # ZH(1,1)
      1  2     9.94520298E-01   # ZH(1,2)
      2  1     9.94520298E-01   # ZH(2,1)
      2  2    -1.04543664E-01   # ZH(2,2)
    Block NMIX
      1  1     9.95491987E-01   # Re(ZN(1,1))
      1  2    -1.79498849E-02   # Re(ZN(1,2))
      1  3     8.57113475E-02   # Re(ZN(1,3))
      1  4    -3.64289843E-02   # Re(ZN(1,4))
      2  1     4.08078351E-02   # Re(ZN(2,1))
      2  2     9.66067231E-01   # Re(ZN(2,2))
      2  3    -2.10338253E-01   # Re(ZN(2,3))
      2  4     1.44245090E-01   # Re(ZN(2,4))
      3  1    -3.37628172E-02   # Re(ZN(3,1))
      3  2     4.88174736E-02   # Re(ZN(3,2))
      3  3     7.03404869E-01   # Re(ZN(3,3))
      3  4     7.08306796E-01   # Re(ZN(3,4))
      4  1     7.86797123E-02   # Re(ZN(4,1))
      4  2    -2.52999529E-01   # Re(ZN(4,2))
      4  3    -6.73522809E-01   # Re(ZN(4,3))
      4  4     6.90049105E-01   # Re(ZN(4,4))
    Block CHARGEMIX
      1  1    -1.00190943E-01   # ZP(1,1)
      1  2     9.94968228E-01   # ZP(1,2)
      2  1     9.94968228E-01   # ZP(2,1)
      2  2     1.00190943E-01   # ZP(2,2)
    Block USQMIX
      1  1     0.00000000E+00   # ZU(1,1)
      1  2     0.00000000E+00   # ZU(1,2)
      1  3     2.45143979E-01   # ZU(1,3)
      1  4     0.00000000E+00   # ZU(1,4)
      1  5     7.81067721E-14   # ZU(1,5)
      1  6     9.69486683E-01   # ZU(1,6)
      2  1     0.00000000E+00   # ZU(2,1)
      2  2     0.00000000E+00   # ZU(2,2)
      2  3    -9.69486683E-01   # ZU(2,3)
      2  4     0.00000000E+00   # ZU(2,4)
      2  5    -3.08861911E-13   # ZU(2,5)
      2  6     2.45143979E-01   # ZU(2,6)
      3  1     0.00000000E+00   # ZU(3,1)
      3  2    -1.04715235E-02   # ZU(3,2)
      3  3     3.18616538E-13   # ZU(3,3)
      3  4     0.00000000E+00   # ZU(3,4)
      3  5    -9.99945172E-01   # ZU(3,5)
      3  6     0.00000000E+00   # ZU(3,6)
      4  1     2.30067962E-05   # ZU(4,1)
      4  2     0.00000000E+00   # ZU(4,2)
      4  3     0.00000000E+00   # ZU(4,3)
      4  4     1.00000000E+00   # ZU(4,4)
      4  5     0.00000000E+00   # ZU(4,5)
      4  6     0.00000000E+00   # ZU(4,6)
      5  1     0.00000000E+00   # ZU(5,1)
      5  2    -9.99945172E-01   # ZU(5,2)
      5  3    -3.33658350E-15   # ZU(5,3)
      5  4     0.00000000E+00   # ZU(5,4)
      5  5     1.04715235E-02   # ZU(5,5)
      5  6     0.00000000E+00   # ZU(5,6)
      6  1     1.00000000E+00   # ZU(6,1)
      6  2     0.00000000E+00   # ZU(6,2)
      6  3     0.00000000E+00   # ZU(6,3)
      6  4    -2.30067962E-05   # ZU(6,4)
      6  5     0.00000000E+00   # ZU(6,5)
      6  6     0.00000000E+00   # ZU(6,6)
    Block SNUMIX
      1  1     0.00000000E+00   # ZV(1,1)
      1  2     0.00000000E+00   # ZV(1,2)
      1  3     1.00000000E+00   # ZV(1,3)
      2  1     0.00000000E+00   # ZV(2,1)
      2  2     1.00000000E+00   # ZV(2,2)
      2  3     0.00000000E+00   # ZV(2,3)
      3  1     1.00000000E+00   # ZV(3,1)
      3  2     0.00000000E+00   # ZV(3,2)
      3  3     0.00000000E+00   # ZV(3,3)
    Block FlexibleSUSYOutput
         0     2.04206021E+16   # HighScale
         1     1.08941472E+03   # SUSYScale
         2     9.11876000E+01   # LowScale
    Block FlexibleSUSYLowEnergy Q= 1.08941472E+03
        21     2.25853630E-10   # Delta(g-2)_muon/2 FlexibleSUSY 1L
    Block EFFHIGGSCOUPLINGS
           25       22       22     2.99452411E-05   # Abs(effective H-Photon-Photon coupling)
           35       22       22     1.10853075E-06   # Abs(effective H-Photon-Photon coupling)
           25       21       21     6.71211022E-05   # Abs(effective H-Gluon-Gluon coupling)
           35       21       21     2.79047785E-06   # Abs(effective H-Gluon-Gluon coupling)
           36       22       22     1.73035166E-06   # Abs(effective A-Photon-Photon coupling)
           36       21       21     3.59315156E-06   # Abs(effective A-Gluon-Gluon coupling)
    Block ALPHA
              -1.04735039E-01   # ArcSin(Pole(ZH(2,2)))
    Block HMIX Q= 1.08941472E+03
         1     6.15787614E+02   # Mu
         2     9.64402142E+00   # vu/vd
         3     2.43666046E+02   # Sqrt(Sqr(vd) + Sqr(vu))
         4     1.47855510E+06   # Sqr(MAh(2))
       101     1.51682262E+05   # BMu
       102     2.51312780E+01   # vd
       103     2.42366584E+02   # vu
    Block Au Q= 1.08941472E+03
      1  1    -1.11540275E+03   # TYu(1,1)/Yu(1,1)
      2  2    -1.11539799E+03   # TYu(2,2)/Yu(2,2)
      3  3    -8.62856153E+02   # TYu(3,3)/Yu(3,3)
    Block Ad Q= 1.08941472E+03
      1  1    -1.36353297E+03   # TYd(1,1)/Yd(1,1)
      2  2    -1.36352854E+03   # TYd(2,2)/Yd(2,2)
      3  3    -1.27508007E+03   # TYd(3,3)/Yd(3,3)
    Block Ae Q= 1.08941472E+03
      1  1    -2.97052953E+02   # TYe(1,1)/Ye(1,1)
      2  2    -2.97047351E+02   # TYe(2,2)/Ye(2,2)
      3  3    -2.95468431E+02   # TYe(3,3)/Ye(3,3)
    Block MSOFT Q= 1.08941472E+03
         1     2.10560904E+02   # MassB
         2     3.89213480E+02   # MassWB
         3     1.10405452E+03   # MassG
        21     1.05177189E+06   # mHd2
        22    -3.46241804E+05   # mHu2
        31     1.04791153E+03   # SignedAbsSqrt(ml2(1,1))
        32     1.04789659E+03   # SignedAbsSqrt(ml2(2,2))
        33     1.04367708E+03   # SignedAbsSqrt(ml2(3,3))
        34     1.01472089E+03   # SignedAbsSqrt(me2(1,1))
        35     1.01468983E+03   # SignedAbsSqrt(me2(2,2))
        36     1.00589853E+03   # SignedAbsSqrt(me2(3,3))
        41     1.38850250E+03   # SignedAbsSqrt(mq2(1,1))
        42     1.38849767E+03   # SignedAbsSqrt(mq2(2,2))
        43     1.21000624E+03   # SignedAbsSqrt(mq2(3,3))
        44     1.36361630E+03   # SignedAbsSqrt(mu2(1,1))
        45     1.36361110E+03   # SignedAbsSqrt(mu2(2,2))
        46     9.70550625E+02   # SignedAbsSqrt(mu2(3,3))
        47     1.36059873E+03   # SignedAbsSqrt(md2(1,1))
        48     1.36059398E+03   # SignedAbsSqrt(md2(2,2))
        49     1.35177526E+03   # SignedAbsSqrt(md2(3,3))

Redirecting info messages to a file
```````````````````````````````````

When FlexibleSUSY is configured with ``--enable-verbose``, a lot of
additional debug output is written to ``stdout`` and ``stderr`` if
FlexibleSUSY is used at the command line.  When the Mathematica
interface is used, this output is redirected to the notebook and
printed if form of messages of type ``FS<model>::info``, where ``<model>``
is the model name.

By default, no more than three messages of the same type are witten to
the notebook.  In order to write all messages to the notebook, set
::

    Off[General::stop];

The function, which writes the messages is called
``FS<model>Message`` and is defined as
::

    FS<model>Message[s_] := Message[FS<model>::info, s]

where ``s`` is the message string.  If one would like to write the
messages to a file, the function can be re-defined to
::

    FS<model>Message[s_] := WriteString["info.txt", s <> "\n"];

Example::

    Get["models/CMSSM/CMSSM_librarylink.m"];
    
    handle = FSCMSSMOpenHandle[
      fsModelParameters -> { m0 -> 125, m12 -> 500, TanBeta -> 10, SignMu -> 1 }
    ];
    
    (* write all messages to "info.txt" *)
    FSCMSSMMessage[s_] := WriteString["info.txt", s <> "\n"];
    
    FSCMSSMCalculateSpectrum[handle]

.. _HSSUSY: hssusy.rst
.. _Observables: model_file.rst
.. _`FlexibleSUSY run-time configuration`: slha_input.rst
