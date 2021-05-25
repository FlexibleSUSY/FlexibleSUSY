.. role:: raw-latex(raw)
    :format: latex

.. role:: mathematica(code)
   :language: mathematica

FlexibleDecay
=============

Prerequisites
+++++++++++++

FlexibleDecay **requires** FlexibleSUSY to be configured with a dedicaded loop library.
See `here`__ for an instruction on how to do it.

__ https://github.com/FlexibleSUSY/FlexibleSUSY/tree/development#support-for-alternative-loop-libraries

Creating model with decays
++++++++++++++++++

Whether decays are created for a given model is controled by viariable

.. code-block:: mathematica
  
  FSCalculateDecays = True;

in model's ``FlexibleSUSY.m.in`` file.

By defaul this variable is initialized to :mathematica:`False`.
In models distributed with FlexibleSUSY which support decays we do set it explicitly to :mathematica:`True`.

Which decays are included is controled by :mathematica:`FSDecayParticles` variable.
By default, it's set to include CP-even and -odd neutral Higgs and charged Higgs.
For example, in SARAH's THDM-II this is equivalen to 

.. code-block:: mathematica

  FSDecayParticles = {hh, Ah, Hpm};
  
Runtime options
+++++++

Runtime options are set in ``FlexibleDecay`` block in the SLHA file

.. code-block::

  Block FlexibleDecay
     0   1       # calculate decays (0 = no, 1 = yes)
     1   1e-5    # minimum BR to print
     2   4       # include higher order corrections in decays (0 = LO, 1 = NLO, 2 = NNLO, 3 = N^3LO, 4 = N^4LO )
     3   1       # use Thomson alpha(0) instead of alpha(m) in decays to γγ and γZ
     4   2       # off-shell decays into VV pair

The options are

0. turn calculation of decay on/off (default = 1)
#. minimal branching ratio to print (default = 1e-5)
#. Maximal order of included higher order corrections (default = 4). Note that not all such corrections. 
#. Use α in the Thomson limit instead of a running one in decays to γγ and γZ. This should minimize higher order corrections.

Flag 4 controls treatment of Higgs decay to gauge bosons

0. no off-shell decays
1. single off-shell decay above VV* threshold (V = W, Z), double offshell below it
2. double off-shell decays also above a VV* threshold

Finally, the loop library used by decays is controlled by flag 31 in block FlexibleSUSY.
For decays the allowed options are 1 and 2.

For example:

.. code-block::

   31   1                    # loop library (1 = COLLIER, 2 = LoopTools)

Example output
++++++++++++++

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
