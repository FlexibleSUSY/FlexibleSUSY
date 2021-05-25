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
#. This item is auto-enumerated

The loop library used by decays is controlled by flag 31 in block FlexibleSUSY

Example output
++++++++++++++
