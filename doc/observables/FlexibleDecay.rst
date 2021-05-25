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

Whether decaysare created for a model is controled by viariable

.. code-block:: mathematica
  
  FSCalculateDecays = True;

in model's ``FlexibleSUSY.m.in`` file.
By defaul this variable is initialized to :mathematica:`False`.
In models distributed with FlexibleSUSY which support decays we set it explicitly to :mathematica:`True`.
Which decays are included is controled by :mathematica:`FSDecayParticles` variable.
By default, it's set to include CP-even and -odd neutral Higgs and charged Higgs.
For example, in SARAH's THDM-II this is equivalen to 

.. code-block:: mathematica

  FSDecayParticles = {hh, Ah, Hpm};
  
Runtime options
+++++++

.. code-block::

  Block FlexibleDecay
     0   1       # calculate decays (0 = no, 1 = yes)
     1   1e-5    # minimum BR to print
     2   4       # include higher order corrections in decays (0 = LO, 1 = NLO, 2 = NNLO, 3 = N^3LO, 4 = N^4LO )
     3   1       # use Thomson alpha(0) instead of alpha(m) in decays to γγ and γZ
     4   2       # off-shell decays into VV pair

The area of a circle is \[ \alpha \]

0. This is the first item
#. This is the second item
#. Enumerators are arabic numbers,
   single letters, or roman numerals
#. List items should be sequentially
   numbered, but need not start at 1
   (although not all formatters will
   honour the first index).
#. This item is auto-enumerated

The loop library used by decays is controlled by flag 31 in block FlexibleSUSY
