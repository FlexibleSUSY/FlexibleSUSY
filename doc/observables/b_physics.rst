===========
b-physics
===========

+++++++++++++++++++++++
:math:`b \to s \gamma`
+++++++++++++++++++++++

To compute :math:`b \to s \gamma` add ``FlexibleSUSYObservable`bsgamma`` into ``ExtraSLHAOutputBlocks`` block in  ``FlexibleSUSY.m.in`` for a model of interest.
As an example, below is the corresponding part from ``MRSSM2CKM``

.. code-block:: mathematica

    ExtraSLHAOutputBlocks = {
        {FlexibleSUSYLowEnergy, {
                {21, FlexibleSUSYObservable`aMuon},
                {23, FlexibleSUSYObservable`EDM[Fe[1]]},
                {24, FlexibleSUSYObservable`EDM[Fe[2]]},
                {25, FlexibleSUSYObservable`EDM[Fe[3]]},
                {26, FlexibleSUSYObservable`BrLToLGamma[Fe[2] -> {Fe[1], VP}]},
                {27, FlexibleSUSYObservable`bsgamma}
            }
        }
    };

By default, calculation of b-physics observables is enabled only for models with a non-trivial CKM matrix.
Currently, these are: CMSSMCKM and MRSSM2CKM (by convention all such model names have suffix ``CKM``).
After a normal ``FlexibleSUSY`` run you'll find in the run directory a file called ``WC_MODELNAME.json``.
This is an input to the flavio_ package by D. Straub.
To use it, first, install ``flavio`` as expained here_.
In the ``utils`` directory of ``FlexibleSUSY`` there's a script ``calc_b_to_s_gamma.py``.
Call it from the directory where the ``.json`` file is located

.. code-block:: bash

    python3 /path/to/calc_b_to_s_gamma.py WC_MRSSM2.json

This will result in an output like

.. code-block:: bash

    BR(B->Xsgamma) =  0.0006411523880676971


.. _flavio: https://flav-io.github.io/

.. _here: https://flav-io.github.io/docs/installation.html