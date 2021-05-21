.. role:: raw-latex(raw)
    :format: latex

FlexibleDecay
=============

The code is documented in 

Prerequisites
+++++++++++++

FlexibleDecay 

Models with decays
++++++++++++++++++

.. code-block:: mathematica

  FSCalculateDecays = True;
  FSDecayParticles = {hh, Ah, Hpm};
  
Runtime
+++++++

.. code-block::

  Block FlexibleDecay
     0   1       # calculate decays (0 = no, 1 = yes)
     1   1e-5    # minimum BR to print
     2   4       # include higher order corrections in decays (0 = LO, 1 = NLO, 2 = NNLO, 3 = N^3LO, 4 = N^4LO )
     3   1       # use Thomson alpha(0) instead of alpha(m) in decays to γγ and γZ
     4   2       # off-shell decays into VV pair

The area of a circle is \[ \alpha \]
