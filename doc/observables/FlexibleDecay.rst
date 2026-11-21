.. role:: raw-latex(raw)
    :format: latex

.. role:: mathematica(code)
   :language: mathematica

FlexibleDecay
=============

FlexibleDecay [2106.05038] calculates two-body decays of the particles of a
model at leading order. A decay that exists at tree level is calculated at
tree level; for scalars, a decay without a tree-level amplitude is calculated
at one loop with all particles of the model in the loop (e.g.
:math:`h \to \gamma\gamma, gg, Z\gamma`). Known higher-order SM QED, QCD and
electroweak corrections are added where possible. Decays of fermions and
vectors are calculated at tree level only.

The decays are calculated in a modified :math:`\overline{\text{MS}}` scheme in
which heavy BSM particles decouple: the SM-like couplings are fixed by the SM
:math:`\overline{\text{MS}}` parameters, while the BSM parameters stay
:math:`\overline{\text{MS}}` or :math:`\overline{\text{DR}}`. Each particle's
decays are calculated at the scale of its mass.

Creating a model with decays
++++++++++++++++++++++++++++

Whether decays are created for a model is controlled by the variable

.. code-block:: mathematica

  FSCalculateDecays = True;

in the model's ``FlexibleSUSY.m.in`` file.
By default it is :mathematica:`False`.
The models distributed with FlexibleSUSY that support decays set it explicitly to :mathematica:`True`.

The decaying particles are selected with :mathematica:`FSDecayParticles`.
Its possible values are :mathematica:`Automatic`, :mathematica:`All` or a list of particles of the model.
The default is

.. code-block:: mathematica

  FSDecayParticles = Automatic;

which means the CP-even, CP-odd and charged Higgs bosons.
For example, in SARAH's THDM-II, :mathematica:`Automatic` is equivalent to

.. code-block:: mathematica

  FSDecayParticles = {hh, Ah, Hpm};

:mathematica:`All` selects all particles of the model.
Other particles can also be added to the list.
For example, to add squark decays in SARAH's MSSM one would write

.. code-block:: mathematica

  FSDecayParticles = {hh, Ah, Hpm, Su, Sd};

Configuration
+++++++++++++

FlexibleDecay **requires** FlexibleSUSY to be configured with a dedicated loop library (currently COLLIER or LoopTools).
See `here`__ for instructions.
Without one, the decays are switched off when the code is generated.

__ https://github.com/FlexibleSUSY/FlexibleSUSY/tree/development#support-for-alternative-loop-libraries

Higgs-physics constraints from HiggsTools (HiggsBounds and HiggsSignals) and from Lilith are optional.
They need FlexibleSUSY to be configured with them, see
`HiggsTools support`__ and `Lilith support`__ in the README.

__ https://github.com/FlexibleSUSY/FlexibleSUSY/tree/development#higgstools-support
__ https://github.com/FlexibleSUSY/FlexibleSUSY/tree/development#lilith-support

Runtime options
+++++++++++++++

SLHA input
~~~~~~~~~~

The runtime options are set in the ``FlexibleDecay`` block of the SLHA_ input file

.. _SLHA: https://inspirehep.net/record/632863

.. code-block::

  Block FlexibleDecay
     0   1       # calculate decays (0 = no, 1 = yes)
     1   1e-5    # minimum BR to print
     2   41      # include higher order corrections in decays SM|BSM (0 = LO, 1 = NLO, 2 = NNLO, 3 = N^3LO, 4 = N^4LO)
     3   1       # use Thomson alpha(0) instead of alpha(m) in decays to γγ and γZ
     4   2       # off-shell decays into VV pair
     5   1       # print loop-induced Higgs couplings
     6   0       # calculate Higgs effc normalized to SM
     7   0       # call HiggsTools
     8   0       # call Lilith
     9   1       # use pole Higgs mixings in decays

The entries are:

0. Calculate the decays (0 = no, 1 = yes; default 1).
#. Minimum branching ratio to print (default 1e-5).
   Use 0 to print all channels.
#. Maximum order of the included higher-order corrections, as two digits
   (default 41).
   The tens digit sets the order of the SM corrections (0 to 4: LO, NLO, NNLO,
   N\ :sup:`3`\ LO, N\ :sup:`4`\ LO), the units digit that of the BSM
   corrections (0 or 1).
   For example, 41 means SM corrections up to N\ :sup:`4`\ LO and NLO BSM
   corrections, 10 means NLO SM corrections only, and 0 means leading order.
   Not every process has corrections up to every order.
   The BSM corrections are currently the two-loop QCD corrections to the
   colored scalar loops in :math:`h \to \gamma\gamma`.
#. Use :math:`\alpha` in the Thomson limit instead of the running one in decays
   to :math:`\gamma\gamma` and :math:`\gamma Z` (0 = no, 1 = yes; default 1).
   This minimizes the higher-order corrections.
#. Treatment of Higgs decays into gauge bosons :math:`V = W, Z` (default 2):

   0. no off-shell decays
   1. on-shell decays above the :math:`VV` threshold, single off-shell decays
      above the :math:`VV^*` threshold, double off-shell decays below it
   2. double off-shell decays also between the :math:`VV^*` and :math:`VV`
      thresholds

#. Print the loop-induced Higgs couplings to :math:`gg`, :math:`\gamma\gamma`
   and :math:`Z\gamma` in the block ``EFFHIGGSCOUPLINGS`` (0 = no, 1 = yes;
   default 1).
#. Calculate the Higgs couplings normalized to those of an SM Higgs boson of the
   same mass, and print them in the blocks ``NORMALIZEDEFFHIGGSCOUPLINGS`` and
   ``IMNORMALIZEDEFFHIGGSCOUPLINGS`` (0 = no, 1 = yes; default 0).
   They are calculated anyway when HiggsTools or Lilith is called.
#. Call HiggsTools (0 = no, 1 = yes; default 0).
   The results are printed in the blocks ``HIGGSBOUNDS`` and ``HIGGSSIGNALS``.
   The HiggsBounds and HiggsSignals datasets are given on the command line with
   ``--higgsbounds-dataset=<directory>`` and ``--higgssignals-dataset=<directory>``.
#. Call Lilith (0 = no, 1 = yes; default 0).
   The Lilith database is given on the command line with ``--lilith-db=<file>``.
#. Use the Higgs mixing matrices of the pole masses instead of the tree-level ones
   in the decays (0 = no, 1 = yes; default 1).
   With 1, the mixing matrices that diagonalize the loop-corrected Higgs mass
   matrices replace the tree-level ones in the couplings.

Entries 2 to 4 and 9 influence only the calculation of Higgs decays.

The loop library used by the decays is selected by entry 31 of the block
``FlexibleSUSY``.
Decays are only calculated with 1 (COLLIER) or 2 (LoopTools); with the other
loop libraries they are skipped.
For example:

.. code-block::

   31    1    # loop library (1 = COLLIER, 2 = LoopTools)

LibraryLink interface
~~~~~~~~~~~~~~~~~~~~~

FlexibleDecay can also be used via the LibraryLink interface (see `here`__).
The options of the ``FlexibleDecay`` block (except entry 0) are passed in
:mathematica:`fdSettings`, for example for the SM:

__ https://github.com/FlexibleSUSY/FlexibleSUSY#mathematica-interface

.. code-block:: mathematica

    Get["models/SM/SM_librarylink.m"];

    (* Create a handle to a model given the input parameters.
       See Options[FSSMOpenHandle] for all default options. *)
    handle = FSSMOpenHandle[
      fsSettings -> { precisionGoal -> 1.*^-4, loopLibrary -> 1 },
      fsModelParameters -> { LambdaIN -> 0.192, Qin -> 1000, QEWSB -> 173.34 },
      fdSettings -> {
        minBRtoPrint -> 1*^-5,              (* FlexibleDecay[1] *)
        maxHigherOrderCorrections -> 41,    (* FlexibleDecay[2] *)
        alphaThomson -> 1,                  (* FlexibleDecay[3] *)
        offShellVV -> 2,                    (* FlexibleDecay[4] *)
        printEffCBlock -> 1,                (* FlexibleDecay[5] *)
        calcNormalizedEffC -> 0,            (* FlexibleDecay[6] *)
        callHiggsTools -> 0,                (* FlexibleDecay[7] *)
        callLilith -> 0,                    (* FlexibleDecay[8] *)
        usePoleHiggsMixings -> 1            (* FlexibleDecay[9] *)
      }
    ];

After computing the spectrum via

.. code-block:: mathematica

    FSSMCalculateSpectrum[handle];

the decays are computed as

.. code-block:: mathematica

    FSSMCalculateDecays[handle];

The normalized effective couplings and the Higgs-physics constraints are
available as

.. code-block:: mathematica

    FSSMCalculateNormalizedEffectiveCouplings[handle];
    FSSMCallHiggsTools[handle, "/path/to/HSDataSet", "/path/to/HBDataSet"];
    FSSMCallLilith[handle, "/path/to/Lilith/data/latestRun2.list"];

Example output
++++++++++++++

SLHA
~~~~

.. code-block::

    Block DCINFO
         1   FlexibleSUSY
         2   2.9.0
         5   SM
         9   4.15.5
    DECAY        25     4.02018349E-03   # hh decays
         5.88051533E-01   2          -5         5  # BR(hh -> barFd(3) Fd(3))
         2.04761572E-01   2         -24        24  # BR(hh -> conjVWp VWp)
         8.64417393E-02   2          21        21  # BR(hh -> VG VG)
         6.21570890E-02   2         -15        15  # BR(hh -> barFe(3) Fe(3))
         2.84500847E-02   2          -4         4  # BR(hh -> barFu(2) Fu(2))
         2.59561973E-02   2          23        23  # BR(hh -> VZ VZ)
         2.25131718E-03   2          22        22  # BR(hh -> VP VP)
         1.44643725E-03   2          22        23  # BR(hh -> VP VZ)
         2.63302268E-04   2          -3         3  # BR(hh -> barFd(2) Fd(2))
         2.20016469E-04   2         -13        13  # BR(hh -> barFe(2) Fe(2))
    Block EFFHIGGSCOUPLINGS
           25          21          21     6.70326269E-05   # hh-VG-VG
           25          22          22     3.05976327E-05   # hh-VP-VP
           25          22          23     5.44958649E-05   # hh-VP-VZ

The output conforms to the SLHA standard.

LibraryLink
~~~~~~~~~~~

The example above gives

.. code-block:: mathematica

    {
      SM -> {
        hh -> {
          25, 0.00399381, {
            {25, {-3,3},   1.05865*10^-6},
            {25, {-4,4},   0.000114401},
            {25, {22,22},  9.05069*10^-6},
            {25, {-13,13}, 8.84506*10^-7},
            {25, {23,23},  0.000104363},
            {25, {-5,5},   0.00233759},
            {25, {-24,24}, 0.000823168},
            {25, {22,23},  5.81436*10^-6},
            {25, {-15,15}, 0.000249883},
            {25, {21,21},  0.000347594}
          }
        }
      }
    }

For each decaying particle the output contains its PDG identifier and its total
width in GeV, followed by the channels with a branching ratio above
``minBRtoPrint``.
The output for every channel, e.g.

.. code-block:: mathematica

    {25, {-15,15}, 0.000249883}

contains the PDG identifiers of the initial and final particles and the partial width in GeV.

Comparison with NMSSMCALC
+++++++++++++++++++++++++

NMSSMCALC_ [1312.4788] calculates the Higgs masses and decays of the
:math:`\mathbb{Z}_3`-invariant NMSSM. Its decay part, BHDECAY, is an
extension of HDECAY [hep-ph/9704448, 1801.09506]; the corrections it includes
are described in section 5 of [1312.4788]. With ``iewh`` in ``bhdecay.in`` it
also includes the NLO SUSY-EW and SUSY-QCD corrections to the two-body
on-shell decays of NMSSMCALCEW [1907.12060, 2012.14889]; the default
``iewh = 0`` switches them off. When comparing NMSSMCALC with FlexibleDecay
for the same parameter point, note the following differences in setup and in
the higher-order corrections. The comparison refers to the code of
NMSSMCALC 3.00 with BHDECAY_rEWN 2.00 and the default settings of both
programs. Where the code and the publications differ, this is noted.

.. _NMSSMCALC: https://www.itp.kit.edu/~maggie/NMSSMCALC/

Setup
~~~~~

* **Parameters:** FlexibleDecay calculates the decays in a modified
  :math:`\overline{\text{MS}}` scheme in which heavy BSM particles decouple:
  the gauge and Yukawa couplings are fixed by the SM :math:`\overline{\text{MS}}`
  parameters at the decay scale, and the remaining BSM parameters are
  :math:`\overline{\text{MS}}` or :math:`\overline{\text{DR}}` [2106.05038].
  NMSSMCALC uses the parameters of its own spectrum calculation (on-shell or
  :math:`\overline{\text{DR}}`, depending on the parameter).
* **Higgs mixing:** with ``FlexibleDecay[9] = 1`` (the default) FlexibleDecay
  replaces the tree-level Higgs mixing matrices by those that diagonalize the
  loop-corrected mass matrices at the pole masses. NMSSMCALC uses the
  loop-corrected mixing matrix of its mass calculation, evaluated at zero
  external momentum so that it is unitary [1312.4788, sec. 3.2], and, for
  ``iewh >= 1``, the finite wave-function normalization factors
  :math:`\mathbf{Z}^H` [0710.5320] of the external Higgs bosons
  [1907.12060].
* **Masses:** both programs use the pole mass of the decaying Higgs boson.

Corrections in FlexibleDecay but not in NMSSMCALC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* :math:`\phi \to q\bar{q}`: the :math:`\mathcal{O}(\alpha)` QED and the mixed
  :math:`\mathcal{O}(\alpha\alpha_s)` corrections. NMSSMCALC includes QED
  corrections only with ``iewh >= 1`` [1907.12060, sec. 4]; [1312.4788]
  includes no electroweak corrections beyond the SUSY-EW part of
  :math:`\Delta_b`.
* :math:`\phi \to \ell^+\ell^-`: the :math:`\mathcal{O}(\alpha)` QED correction
  (in NMSSMCALC only with ``iewh >= 1``).
* :math:`\phi \to gg`: the NLO and NNLO QCD corrections to the top loop include
  terms of higher order in :math:`\tau = m_\phi^2/(4 m_t^2)` (up to
  :math:`\tau^5` at NLO and :math:`\tau^2` at NNLO for CP-even states); NMSSMCALC
  uses the heavy-top limit.
* :math:`h \to Z\gamma`: the QCD correction to the quark loops (in the
  heavy-quark limit, for :math:`m_h < 0.8\, m_q` and :math:`m_Z < 0.75\, m_h`).
  NMSSMCALC leaves it out as numerically small [1312.4788, sec. 5.1].
* :math:`\phi \to Z\gamma`: all particles of the model run in the loop,
  including charginos, sfermions and (for CP-even states) charged Higgs bosons.
  The NMSSMCALC code (``bhdecay_rEWN.F``, and likewise ``bhdecay.f``) includes
  only the top, bottom and charm quarks and the :math:`\tau` lepton, and, for
  CP-even states, the :math:`W` boson and the charged Higgs boson. This
  differs from [1312.4788, sec. 5.1], which lists sfermion and chargino loops
  for CP-even and chargino loops for CP-odd states. For a singlet-like CP-odd
  state, whose :math:`Z\gamma` width comes mainly from chargino loops, this
  makes a large difference.
* Every loop-induced decay is calculated at one loop with the full particle
  content of the model, including channels that NMSSMCALC doesn't calculate,
  such as :math:`H^\pm \to W^\pm\gamma` and :math:`H^\pm \to W^\pm Z`.

Corrections in NMSSMCALC but not in FlexibleDecay
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* :math:`\phi \to b\bar{b}`: the resummation of the
  :math:`\tan\beta`-enhanced SUSY-QCD and top-induced SUSY-EW corrections
  (:math:`\Delta_b`) in the effective bottom Yukawa couplings
  [hep-ph/9912516, hep-ph/0305101], with the higher-order corrections to them
  [0808.0087, 1001.1935, 1007.0693], plus the remaining
  :math:`\mathcal{O}(\alpha_s)` SUSY-QCD corrections (real NMSSM only)
  [1312.4788, secs. 5.1 and 5.3]. In FlexibleDecay the
  bottom Yukawa coupling is fixed by the SM :math:`\overline{\text{MS}}` bottom
  mass with the tree-level relation, so :math:`\Delta_b` is not included.
* :math:`\phi \to \ell^+\ell^-` and :math:`s\bar{s}`: the corresponding
  resummation (:math:`\Delta_\ell`, from bino and wino loops, and
  :math:`\Delta_s`) in the lepton and strange-quark Yukawa couplings
  [1312.4788, sec. 5.3].
* :math:`\phi \to gg`:

  * the QCD correction factors (derived for the top loop) multiply the whole
    LO amplitude squared, including the bottom, charm and squark loops; in
    FlexibleDecay the corrections are
    added only for the top loop (as :math:`\Gamma^\text{LO}_{\text{top}}
    \delta`), and they are switched off for :math:`\tau > 0.7`
  * an additional QCD correction for squark loops in the heavy-squark limit
    [1312.4788, sec. 5.1]
  * the separation :math:`\phi \to g g^* \to g q\bar{q}` into the
    :math:`q\bar{q}` channels (``NF-GG``)

* :math:`\phi \to Z\gamma`: NMSSMCALC returns zero below the :math:`Z`
  threshold, so this is not a correction but a difference for light states.
* :math:`H^\pm` decays: QCD corrections, e.g. to :math:`H^+ \to t\bar{b}`
  [1312.4788, sec. 5.1].
  FlexibleDecay calculates charged Higgs decays at tree level without
  corrections.
* Decays into squarks: SUSY-QCD corrections (real NMSSM)
  [1312.4788, sec. 5.1]. FlexibleDecay calculates decays
  into BSM particles at tree level.
* Three-body decays below threshold: :math:`\phi \to t t^*`,
  :math:`A_i \to Z^* h_j`, :math:`\phi \to W^{\mp *} H^\pm`,
  :math:`H^+ \to t^* b` and similar (``ON-SHELL = 0``), and off-shell Higgs
  pairs [1312.4788, sec. 5.1]. FlexibleDecay includes
  off-shell decays only for :math:`\phi \to VV^*` and :math:`V^*V^*`
  (:math:`V = W, Z`).
* With ``iewh = 2``: NLO SUSY-EW and SUSY-QCD corrections to the two-body
  on-shell decays of the neutral Higgs bosons into fermions, :math:`VV`,
  :math:`hZ`, Higgs pairs, charginos, neutralinos and squarks [1907.12060],
  and of the charged Higgs bosons into fermions, electroweakinos, sleptons
  and squarks [2012.14889]. The loop-induced decays into :math:`gg`,
  :math:`\gamma\gamma` and :math:`Z\gamma` get no such corrections, as they
  would be of two-loop order [1907.12060].

Corrections in both
~~~~~~~~~~~~~~~~~~~

* :math:`\phi \to q\bar{q}`: massless QCD corrections up to
  :math:`\mathcal{O}(\alpha_s^4)`, the massive :math:`\mathcal{O}(\alpha_s)`
  correction with an interpolation between the running- and pole-mass
  regimes, and the top-induced :math:`\mathcal{O}(\alpha_s^2)` correction
  (for NMSSMCALC see [1312.4788, sec. 5.1]).
* :math:`\phi \to gg`: the QCD corrections to the top loop in the heavy-top
  limit up to N\ :sup:`3`\ LO for CP-even and NNLO for CP-odd states
  (NMSSMCALC: [1312.4788, sec. 5.1], which states N\ :sup:`3`\ LO without
  distinguishing the CP-odd case).
* :math:`\phi \to \gamma\gamma`: the two-loop QCD corrections to the quark
  loops with the full mass dependence [hep-ph/9504378], and to the squark
  loops of CP-even states [hep-ph/0612254]. Both use :math:`\alpha` in the Thomson limit for the photon couplings
  (in FlexibleDecay with ``FlexibleDecay[3] = 1``, the default).
* :math:`\phi \to VV`: double and single off-shell decays below the threshold.
  Neither includes electroweak corrections to them by default; NMSSMCALC adds
  the NLO SUSY-EW corrections with ``iewh = 2`` [1907.12060].

References
~~~~~~~~~~

* [0710.5320] K. Williams, G. Weiglein, Phys. Lett. B 660 (2008) 217,
  `arXiv:0710.5320 <https://arxiv.org/abs/0710.5320>`__
* [0808.0087] D. Noth, M. Spira, Phys. Rev. Lett. 101 (2008) 181801,
  `arXiv:0808.0087 <https://arxiv.org/abs/0808.0087>`__
* [1001.1935] D. Noth, M. Spira, JHEP 1106 (2011) 084,
  `arXiv:1001.1935 <https://arxiv.org/abs/1001.1935>`__
* [1007.0693] L. Mihaila, C. Reisser, JHEP 1008 (2010) 021,
  `arXiv:1007.0693 <https://arxiv.org/abs/1007.0693>`__
* [1312.4788] J. Baglio et al., *NMSSMCALC: A Program Package for the
  Calculation of Loop-Corrected Higgs Boson Masses and Decay Widths in the
  (Complex) NMSSM*, Comput. Phys. Commun. 185 (2014) 3372,
  `arXiv:1312.4788 <https://arxiv.org/abs/1312.4788>`__
* [1801.09506] A. Djouadi, J. Kalinowski, M. Mühlleitner, M. Spira,
  *HDECAY: Twenty++ Years After*,
  `arXiv:1801.09506 <https://arxiv.org/abs/1801.09506>`__
* [1907.12060] J. Baglio, T. N. Dao, M. Mühlleitner, *One-Loop Corrections to
  the Two-Body Decays of the Neutral Higgs Bosons in the Complex NMSSM*,
  Eur. Phys. J. C 80 (2020) 960,
  `arXiv:1907.12060 <https://arxiv.org/abs/1907.12060>`__
* [2012.14889] T. N. Dao, M. Mühlleitner, S. Patel, K. Sakurai, *One-loop
  Corrections to the Two-Body Decays of the Charged Higgs Bosons in the Real
  and Complex NMSSM*, `arXiv:2012.14889 <https://arxiv.org/abs/2012.14889>`__
* [2106.05038] FlexibleDecay, Comput. Phys. Commun. 283 (2023) 108584,
  `arXiv:2106.05038 <https://arxiv.org/abs/2106.05038>`__
* [hep-ph/0305101] J. Guasch, P. Häfliger, M. Spira, Phys. Rev. D 68 (2003)
  115001, `arXiv:hep-ph/0305101 <https://arxiv.org/abs/hep-ph/0305101>`__
* [hep-ph/0612254] M. Mühlleitner, M. Spira, Nucl. Phys. B 790 (2008) 1,
  `arXiv:hep-ph/0612254 <https://arxiv.org/abs/hep-ph/0612254>`__
* [hep-ph/9504378] M. Spira, A. Djouadi, D. Graudenz, P. M. Zerwas,
  Nucl. Phys. B 453 (1995) 17,
  `arXiv:hep-ph/9504378 <https://arxiv.org/abs/hep-ph/9504378>`__
* [hep-ph/9704448] A. Djouadi, J. Kalinowski, M. Spira, *HDECAY*,
  Comput. Phys. Commun. 108 (1998) 56,
  `arXiv:hep-ph/9704448 <https://arxiv.org/abs/hep-ph/9704448>`__
* [hep-ph/9912516] M. Carena, D. Garcia, U. Nierste, C. E. M. Wagner,
  Nucl. Phys. B 577 (2000) 88,
  `arXiv:hep-ph/9912516 <https://arxiv.org/abs/hep-ph/9912516>`__
