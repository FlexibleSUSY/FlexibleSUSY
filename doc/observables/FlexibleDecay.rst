.. role:: raw-latex(raw)
    :format: latex

.. role:: mathematica(code)
   :language: mathematica

FlexibleDecay
=============

Creating model with decays
++++++++++++++++++

Whether decays are created for a given model is controled by the viariable

.. code-block:: mathematica
  
  FSCalculateDecays = True;

in the model's ``FlexibleSUSY.m.in`` file.

By defaul this variable is initialized to :mathematica:`False`.
In models distributed with FlexibleSUSY which support decays we do set it explicitly to :mathematica:`True`.

Which decays are included is controled by :mathematica:`FSDecayParticles` variable.
It's possible values are :mathematica:`Automatic` or a list of scalars present in the model.
By default it's set to 

.. code-block:: mathematica

  FSDecayParticles = Automatic;

which means CP-even and -odd neutral Higgs and charged Higgs bosons.
As an example, in SARAH's THDM-II, :mathematica:`Automatic` is equivalen to setting

.. code-block:: mathematica

  FSDecayParticles = {hh, Ah, Hpm};
  
One can also add other scalar.
For example, to add squark decays in the SARAH's MSSM one would write

.. code-block:: mathematica

  FSDecayParticles = {hh, Ah, Hpm, Su, Sd};

Configuration
+++++++++++++

FlexibleDecay **requires** FlexibleSUSY to be configured with a dedicaded loop library (currently this is COLLIER or LoopTools).
See `here`__ for an instruction on how to do it.

__ https://github.com/FlexibleSUSY/FlexibleSUSY/tree/development#support-for-alternative-loop-libraries

Runtime options
+++++++

SLHA input
~~~~~~~~~~

Runtime options are set in ``FlexibleDecay`` block in the SLHA_ input file

.. _SLHA: https://inspirehep.net/record/632863

.. code-block::

  Block FlexibleDecay
     0   1       # calculate decays (0 = no, 1 = yes)
     1   1e-5    # minimum BR to print
     2   4       # include higher order corrections in decays (0 = LO, 1 = NLO, 2 = NNLO, 3 = N^3LO, 4 = N^4LO )
     3   1       # use Thomson alpha(0) instead of alpha(m) in decays to γγ and γZ
     4   2       # off-shell decays into VV pair

The flags are:

0. Turn calculation of decay on/off (default = 1).
#. Minimal branching ratio to print (default = 1e-5).
#. Maximal order of included higher order corrections (default = 4). Note that not all such corrections. 
#. Use α in the Thomson limit instead of a running one in decays to γγ and γZ. This should minimize higher order corrections.

Flag 4 controls treatment of Higgs decay to gauge bosons

0. no off-shell decays
1. on-shell decays above the VV threshold, single off-shell decays above the VV* threshold (V = W, Z), double off-shell below it
2. double off-shell decays also between a VV* and VV thresholds

Flags 2, 3 and 4 influnce only the calculation of Higgs decays.

Finally, the loop library used by decays is controlled by flag 31 in block FlexibleSUSY.
For decays the allowed vlues are 1 and 2.
To compute decays one could set it for example as

.. code-block::

   31    1    # loop library (1 = COLLIER, 2 = LoopTools)

LibraryLink interface
~~~~~~~~~~~~~~~~~

FlexibleDecay can also be used via the LibraryLink interface (see `here`__).
The same options as in the case of SLHA input can be passed as (using CMSSM as an example)

__ https://github.com/FlexibleSUSY/FlexibleSUSY#mathematica-interface

.. code-block:: mathematica

    Get["models/CMSSM/CMSSM_librarylink.m"];

    (* Create a handle to a model given the input parameters.
       See Options[FSCMSSMOpenHandle] for all default options. *)
    handle = FSCMSSMOpenHandle[
      fsSettings -> { precisionGoal -> 1.*^-4 },
      fsSMParameters -> { Mt -> 173.3 },
      fsModelParameters -> {
      m0 -> 125, m12 -> 500, TanBeta -> 10, SignMu -> 1, Azero -> 0 },
      fdSettings -> { minBRtoPrint -> 1*^-5, maxHigherOrderCorrections -> 4, alphaThomson -> 1, offShellVV -> 2 }
    ];

The meaning of symbols in :mathematica:`fdSettings` is the same as options in the `FlexibleDecay` block (as explained in the previous section).
After computing the spectrum via

.. code-block:: mathematica

    FSCMSSMCalculateSpectrum[handle];

Decays can be computed as

.. code-block:: mathematica

    FSCMSSMCalculateDecays[handle];

Example output
++++++++++++++

SLHA
~~~~

.. code-block::

    Block DCINFO
        1   FlexibleSUSY
        2   2.6.0
        5   SM
        9   4.14.3
    DECAY        25     4.01909364E-03   # hh decays
        5.88154048E-01   2          -5         5  # BR(hh -> barFd(3) Fd(3))
        2.04644925E-01   2         -24        24  # BR(hh -> conjVWp VWp)
        8.64458085E-02   2          21        21  # BR(hh -> VG VG)
        6.21678883E-02   2         -15        15  # BR(hh -> barFe(3) Fe(3))
        2.84471939E-02   2          -4         4  # BR(hh -> barFu(2) Fu(2))
        2.59621707E-02   2          23        23  # BR(hh -> VZ VZ)
        2.25173904E-03   2          22        22  # BR(hh -> VP VP)
        1.44211112E-03   2          22        23  # BR(hh -> VP VZ)
        2.63348187E-04   2          -3         3  # BR(hh -> barFd(2) Fd(2))
        2.20054695E-04   2         -13        13  # BR(hh -> barFe(2) Fe(2))
     
The output conforms to the SLHA standard.

LibraryLink
~~~~~~~~

.. code-block:: mathematica

    {
      SM -> {
        Decays[hh] -> { 
          25, 0.00198076, {
            {25, {-15,15}, 0.000157635}, 
            {25, {23,23},  3.16863*10^-7}, 
            {25, {-24,24}, 1.14636*10^-6}, 
            {25, {-3,3},   7.44681*10^-7},
            {25, {22,22},  1.8801*10^-6}, 
            {25, {-13,13}, 5.58985*10^-7}, 
            {25, {-5,5},   0.00164052}, 
            {25, {-4,4},   0.0000812031}, 
            {25, {21,21},  0.0000967487}
          }
        }
      }
    }

At the top of the block we get a PDG id of particle whose with we are computing as well as its total width.
The output for every channel, e.g.

.. code-block:: mathematica

    {25, {-15,15}, 0.000157635}
    
contains PDG identifiers for in and out particles and a partial width in GeV.
