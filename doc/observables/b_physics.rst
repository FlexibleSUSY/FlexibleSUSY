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

+++++++++++++++++++++++++++++++++++++++++++++++
:math:`b \to s \ell \ell` Wilson coefficients
+++++++++++++++++++++++++++++++++++++++++++++++

The Wilson coefficients of :math:`d_i \ell \to d_j \ell` (e.g.
:math:`b \to s \mu\mu`) are calculated at one loop with
NPointFunctions [2402.14630]_ by adding
``FlexibleSUSYObservable`BrDLToDL`` to the block ``FWCOEF`` of
``ExtraSLHAOutputBlocks``, as in ``CMSSMCKM`` and ``MRSSM2CKM``:

.. code-block:: mathematica

    ExtraSLHAOutputBlocks = {
        ...
        {FWCOEF, {
                {1, FlexibleSUSYObservable`BrDLToDL[{Fd[3], Fe[2]} -> {Fd[2], Fe[2]}, {Vectors, Scalars, Boxes}, 1]},
                {2, FlexibleSUSYObservable`BrDLToDL[{Fd[3], Fe[1]} -> {Fd[2], Fe[1]}, {Vectors, Scalars, Boxes}, 1]}
            }
        },
        {IMFWCOEF, { (* the same entries *) }}
    };

calculates :math:`b \to s \mu\mu` and :math:`b \to s ee`. Both are needed
for lepton-universality ratios such as :math:`R_K` and :math:`R_{K^*}`:
coefficients that are missing are taken as zero by flavio.

The first argument gives the incoming and outgoing down-type quark and the
lepton, the second the contributions to calculate and the third the loop
order (0 or 1):

* ``Boxes``: box diagrams, including those with a W boson and a new
  particle (e.g. a charged Higgs boson) in the loop,
* ``Scalars``: penguins with a scalar (e.g. Higgs bosons) between the quark
  and the lepton line,
* ``Vectors``: penguins with a vector boson (Z) between the quark and the
  lepton line.

The photon penguin is always included at one loop (through the form factors
also used for :math:`b \to s \gamma`). In the penguins with a self-energy on
an external quark leg the quark masses are kept, and the Z propagator is
evaluated at :math:`q^2 = 0` as in the vertex penguins, so that the UV
divergences cancel and the Z penguins decouple.

The Wilson coefficients contain only the new physics contributions, as
flavio and the FLHA entries written here expect (flavio adds the Standard
Model part itself; it has no Standard Model contribution to
:math:`C_S, C_P`). Precisely:

* **Removed as Standard Model part:** diagrams in which all fields (in the
  loop and in the tree-level propagator) are Standard Model-like: quarks,
  leptons, neutrinos, gluon, photon, Z, W, the Goldstone bosons and the
  SM-like CP-even Higgs boson, which is taken as the SM Higgs boson (as in
  SARAH/FlavorKit). The SM-like Higgs boson is the one whose pole mass is
  closest to the SM Higgs mass of the input (``FlexibleSUSYInput[1]``,
  default 125.09 GeV; its DR-bar mass if the pole masses aren't
  calculated). These diagrams are removed with the couplings of the BSM
  model, so the deviations of the couplings of the SM-like Higgs boson from
  the Standard Model in them are not included. In the MSSM these are of order
  :math:`M_Z^2/M_A^2` times the Standard Model Higgs penguin, whose
  contribution to :math:`C_S` is about :math:`10^{-4}\,\text{GeV}^{-1}`.
  In models where the SM-like Higgs boson has flavour-changing couplings to
  quarks, the diagrams with these couplings and only Standard Model-like
  fields are missing as well. In the photon penguin and :math:`C_7, C_8`
  (form factors of :math:`b \to s \gamma`) only the quarks, leptons, gauge
  bosons and Goldstone bosons count as Standard Model fields.
* **Included:** all other one-loop diagrams with box, Z-, photon- and
  scalar-penguin topologies, including those with a W boson and a new
  particle in the loop (e.g. the W-charged Higgs boxes and penguins) and
  those with Standard Model particles in the loop and a new particle in the
  tree-level propagator (e.g. the heavy Higgs bosons H, A).
* **Missing:** box diagrams with two vector bosons in the loop (pure
  Standard Model W-W boxes in the MSSM, but e.g. W-W' or Z-Z' boxes in models
  with new vector bosons), and QCD or electroweak corrections beyond one loop
  (e.g. the :math:`\tan\beta`-enhanced corrections to the Higgs penguins).

The new physics contributions to the Wilson coefficients
:math:`C_7, C_8, C_9, C_{10}, C_S, C_P` and their primed counterparts are
written to the blocks ``FWCOEF`` (real parts) and ``IMFWCOEF`` (imaginary
parts) of the Flavour Les Houches Accord (FLHA) [1008.0762]_, e.g. for
:math:`b \to s \mu\mu`:

.. code-block::

    Block FWCOEF Q= 1.66022932E+03
         0305     4422   00   1     ...   # C7 1 loop
         0305     4322   00   1     ...   # C7' 1 loop
         0305     6421   00   1     ...   # C8 1 loop
         0305     6321   00   1     ...   # C8' 1 loop
     03051313     4133   00   1     ...   # C9 Vectors+Scalars+Boxes 1 loop
     03051313     4233   00   1     ...   # C9' Vectors+Scalars+Boxes 1 loop
     03051313     4137   00   1     ...   # C10 Vectors+Scalars+Boxes 1 loop
     03051313     4237   00   1     ...   # C10' Vectors+Scalars+Boxes 1 loop
     03051313     3230   00   1     ...   # CS Vectors+Scalars+Boxes 1 loop
     03051313     3130   00   1     ...   # CS' Vectors+Scalars+Boxes 1 loop
     03051313     3234   00   1     ...   # CP Vectors+Scalars+Boxes 1 loop
     03051313     3134   00   1     ...   # CP' Vectors+Scalars+Boxes 1 loop

The fermion code 0305 stands for :math:`\bar{s} \ldots b` and 1313 for
:math:`\bar{\mu} \ldots \mu`; the operator codes are those of the FLHA
(e.g. 4133 for :math:`(\bar{s} \gamma^\mu P_L b)(\bar{\mu}
\gamma_\mu \mu)`, 3230 for :math:`(\bar{s} P_R b)(\bar{\mu} \mu)`),
and ``00 1`` means leading order, new physics only. The coefficients are
normalized as in the flavio basis (see below), which is also the
normalization of the FLHA and of SuperIso:
:math:`\mathcal{L} = \frac{4 G_F}{\sqrt{2}} V_{tb} V_{ts}^* \sum_i C_i O_i`
with :math:`O_7 = \frac{e}{16\pi^2} m_b (\bar{s} \sigma^{\mu\nu} P_R b)
F_{\mu\nu}`, :math:`O_9 = \frac{e^2}{16\pi^2} (\bar{s} \gamma^\mu P_L
b)(\bar{\ell} \gamma_\mu \ell)`, :math:`O_S = \frac{e^2}{16\pi^2} m_b
(\bar{s} P_R b)(\bar{\ell} \ell)` etc. They are given at the scale Q of the
block in the quark phase convention of the CKM matrix of the input
(``VCKMIN``, PDG convention). SuperIso reads :math:`C_{7,8,9,10}^{(\prime)}`
from these blocks (its program ``flha.x``).

++++++++++++++++++++
Interface to flavio
++++++++++++++++++++

By default, calculation of b-physics observables is enabled only for models with a non-trivial CKM matrix.
Currently, these are: CMSSMCKM and MRSSM2CKM (by convention all such model names have suffix ``CKM``).
After a normal ``FlexibleSUSY`` run you'll find in the run directory a file called ``WC_MODELNAME.json``.
It contains the new physics contributions to the Wilson coefficients in the
flavio basis of the weak effective theory in the WCxf format [1712.05298]_:
:math:`C_7, C_7', C_8, C_8'` from ``bsgamma``, and
:math:`C_9, C_9', C_{10}, C_{10}', C_S, C_S', C_P, C_P'` from ``BrDLToDL``
for :math:`b \to s \ell\ell` and :math:`b \to d \ell\ell`
(e.g. ``C9_bsmumu``).
The file is rewritten in each run (the coefficients of both observables are
collected if they are calculated at the same scale).
This is an input to the flavio_ package by D. Straub.
To use it, first, install ``flavio`` as expained here_.
In the ``utils`` directory of ``FlexibleSUSY`` there's a script
``calc_flavio_observables.py``, which calculates flavour observables with
flavio from this file.
Call it from the directory where the ``.json`` file is located (here for
CMSSMCKM with m0 = 500 GeV, m1/2 = 1 TeV, A0 = 0, tan(beta) = 40, mu > 0)

.. code-block:: bash

    python3 /path/to/calc_flavio_observables.py WC_CMSSMCKM.json

This will result in an output like

.. code-block:: bash

    BR(B->Xsgamma)           = 3.188819e-04   SM: 3.309692e-04   ratio: 0.963479
    BR(Bs->mumu)             = 4.714078e-09   SM: 3.628931e-09   ratio: 1.299026
    <Rmue>(B+->Kll):1.1:6    = 1.000760e+00   SM: 1.000416e+00   ratio: 1.000344
    <Rmue>(B0->K*ll):1.1:6   = 9.967367e-01   SM: 9.964532e-01   ratio: 1.000284

Other flavio observables can be given after the file name; arguments such as
a :math:`q^2` bin are appended with colons, e.g. ``'<Rmue>(B+->Kll):1.1:6'``.
The older script
``calc_b_to_s_gamma.py`` only calculates BR(B->Xsgamma).

.. [1008.0762] Flavour Les Houches Accord, `arxiv:1008.0762 <https://arxiv.org/abs/1008.0762>`_
.. [1712.05298] WCxf, `arxiv:1712.05298 <https://arxiv.org/abs/1712.05298>`_
.. [2402.14630] NPointFunctions, `CPC 302 (2024) 109244 <https://inspirehep.net/literature/2760779>`_ [`arxiv:2402.14630 <https://arxiv.org/abs/2402.14630>`_]

.. _flavio: https://flav-io.github.io/

.. _here: https://flav-io.github.io/docs/installation.html