.. sectnum::

SLHA input parameters
=====================

.. contents:: Table of Contents

FlexibleSUSY configuration block (FlexibleSUSY)
-----------------------------------------------

**Block name**: ``FlexibleSUSY``

**Default values**::

    Block FlexibleSUSY
        0   1.0e-04   # precision goal
        1   0         # max. iterations (0 = automatic)
        2   0         # algorithm (0 = all, 1 = two_scale, 2 = semi_analytic)
        3   0         # calculate SM pole masses
        4   4         # pole mass loop order
        5   4         # EWSB loop order
        6   4         # beta-functions loop order
        7   4         # threshold corrections loop order
        8   1         # Higgs 2-loop corrections O(alpha_t alpha_s)
        9   1         # Higgs 2-loop corrections O(alpha_b alpha_s)
       10   1         # Higgs 2-loop corrections O((alpha_t + alpha_b)^2)
       11   1         # Higgs 2-loop corrections O(alpha_tau^2)
       12   0         # force output
       13   3         # Top pole mass QCD corrections (0 = 1L, 1 = 2L, 2 = 3L, 3 = 4L)
       14   1.0e-11   # beta-function zero threshold
       15   0         # calculate observables (a_muon, effective couplings)
       16   0         # force positive majorana masses
       17   0         # pole mass scale
       18   0         # pole mass scale in the EFT (0 = min(SUSY scale, Mt))
       19   0         # EFT matching scale (0 = SUSY scale)
       20   2         # EFT loop order for upwards matching (SM -> BSM)
       21   1         # EFT loop order for downwards matching (BSM -> SM)
       22   0         # EFT index of SM-like Higgs in the BSM model (0 = lightest Higgs)
       23   1         # calculate BSM pole masses
       24   124111421 # individual threshold correction loop orders
       25   0         # ren. scheme for Higgs 3L corrections (0 = DR', 1 = MDR', 2 = H3m)
       26   1         # Higgs 3-loop corrections O(alpha_t alpha_s^2)
       27   1         # Higgs 3-loop corrections O(alpha_b alpha_s^2)
       28   1         # Higgs 3-loop corrections O(alpha_t^2 alpha_s)
       29   1         # Higgs 3-loop corrections O(alpha_t^3)
       30   1         # Higgs 4-loop corrections O(alpha_t alpha_s^3)
       31   0         # Loop library to use (0 = softsusy, 1 = COLLIER, 2 = LoopTools, 3 = fflite)

**Description**:

The ``FlexibleSUSY`` block contains fields to configure the spectrum
calculation at run-time.  For example, in the ``FlexibleSUSY`` block the
renormalization group running precsision, the beta function loop order
or the loop order of the pole mass calculation can be selected.

======= ======================================== ===================================================== ======================================
 index   description                              possible values                                       default value
======= ======================================== ===================================================== ======================================
  0      precision goal                           any positive double                                   1.0e-4
  1      max. number of iterations                any positive double                                   0 (= automatic)
  2      BC solver                                0 (all), 1 (two-scale) or 2 (semi-analytic)           0 (= all)
  3      calculate SM pole masses                 0 (no) or 1 (yes)                                     0 (= no)
  4      pole mass loop order                     0, 1, 2, 3, 4                                         4 (= 4-loop)
  5      EWSB loop order                          0, 1, 2, 3, 4                                         4 (= 4-loop)
  6      beta function loop order                 0, 1, 2, 3, 4, 5                                      4 (= 4-loop)
  7      threshold corrections loop order         0, 1, 2, 3, 4                                         4 (= 4-loop)
  8      higgs 2-loop correction O(at as)         0, 1                                                  1 (= enabled)
  9      higgs 2-loop correction O(ab as)         0, 1                                                  1 (= enabled)
 10      higgs 2-loop correction O(at at)         0, 1                                                  1 (= enabled)
 11      higgs 2-loop correction O(atau atau)     0, 1                                                  1 (= enabled)
 12      force output                             0 (no) or 1 (yes)                                     0 (= no)
 13      top quark pole QCD corrections           0 (1L), 1 (2L), 2 (3L), 3 (4L)                        3 (= 4L QCD)
 14      beta function zero threshold             any positive double                                   1.0e-11
 15      calculate observables                    0 (no) or 1 (yes)                                     0 (= no)
 16      force positive Majorana masses           0 (no) or 1 (yes)                                     0 (= no)
 17      pole mass scale                          any positive double                                   0 (= SUSY scale)
 18      EFT pole mass scale                      any positive double                                   0 (= minimum of {Mt, SUSY scale})
 19      EFT matching scale                       any positive double                                   0 (= SUSY scale)
 20      EFT loop order for upwards matching      0, 1, 2                                               2 (= 2-loop)
 21      EFT loop order for downwards matching    0, 1                                                  1 (= 1-loop)
 22      EFT Higgs index                          any integer >= 0                                      0 (= lightest)
 23      calculate pole masses of BSM particles   0 (no) or 1 (yes)                                     1 (= yes)
 24      individual threshold corrections         positive integer                                      124111421
 25      ren. scheme for higgs 3L corrections     0 (DR'), 1 (MDR'), 2 (H3m)                            0 (= DR')
 26      higgs 3-loop correction O(at as^2)       0, 1                                                  1 (= enabled)
 27      higgs 3-loop correction O(ab as^2)       0, 1                                                  1 (= enabled)
 28      higgs 3-loop correction O(at^2 as)       0, 1                                                  1 (= enabled)
 29      higgs 3-loop correction O(at^3)          0, 1                                                  1 (= enabled)
 30      higgs 4-loop correction O(at as^3)       0, 1                                                  1 (= enabled)
 31      Loop library to use                      0 (softsusy), 1 (COLLIER), 2 (LoopTools), 3 (fflite)  0 (= softsusy)
======= ======================================== ===================================================== ======================================

Precision goal (``FlexibleSUSY[0]``)
````````````````````````````````````

FlexibleSUSY solves the given boundary value problem (BVP) by running
all model parameters to each scale and imposing the corresponding
boundary conditions until a convergent solution has been found or the
maximum number of iterations has been reached.  In ``FlexibleSUSY[0]``,
precision goal of the BVP solver can be specified.  The precision goal
determines

- the precision of the numerical solution of the RGEs,

- the precision of the numerical solution of the EWSB equations and

- to test whether the BVP solver has found a convergent solution.


Maximum number of iterations (``FlexibleSUSY[1]``)
``````````````````````````````````````````````````

FlexibleSUSY solves the given boundary value problem (BVP) by running
to each scale and imposing the corresponding boundary conditions until
a convergent solution has been found or the maximum number of
iterations, :math:`N_{\text{max.it.}}`, has been reached.  In
``FlexibleSUSY[1]``, the maximum number of iterations
:math:`N_{\text{max.it.}}` used to solve the BVP can be specified.  If
:math:`N_{\text{max.it.}}` is set to ``0``, the maximum number of
iterations is set to :math:`N_{\text{max.it.}} = -10 \log_{10}(p),`
where :math:`p` is the precision goal specified in
``FlexibleSUSY[0]``.

BVP solver (``FlexibleSUSY[2]``)
````````````````````````````````

Choses the boundary value problem (BVP) solver: 0 = all that are
enabled (starting with the two-scale solver, if present), 1 =
two-scale solver (if present), 2 = semi-analytic solver (if present).

Calculate pole masses of Standard Model particles (``FlexibleSUSY[3]``)
```````````````````````````````````````````````````````````````````````

Calculate pole masses of Standard Model particles: 0 = do not
calculate Standard Model pole masses, 1 = calculate the Standard Model
pole masses.

Pole mass loop order (``FlexibleSUSY[4]``)
``````````````````````````````````````````

Maximum pole mass loop order.  0 = tree-level, 1 = 1-loop, 2 = 2-loop
(if available), 3 = 3-loop (if available).

EWSB loop order (``FlexibleSUSY[5]``)
`````````````````````````````````````

Maximum loop order of the electroweak symmetry breaking (EWSB)
equations.  0 = tree-level, 1 = 1-loop, 2 = 2-loop (if available), 3 =
3-loop (if available).

.. important:: The EWSB loop order should always be set to the same
               value as the pole mass loop order!

beta-function loop order (``FlexibleSUSY[6]``)
``````````````````````````````````````````````

Loop order of the renormalization group running.  0 = no running, 1 =
1-loop running, 2 = 2-loop running, 3 = 3-loop running (if available),
etc.

Threshold correction loop order (``FlexibleSUSY[7]``)
`````````````````````````````````````````````````````

Using the flag ``FlexibleSUSY[7]`` the "global" loop order of the
threshold corrections of the SM to the full BSM model can be selected.
The threshold corrections affect the determination of the running BSM
model parameters :math:`\alpha_{\text{em}}`, :math:`\alpha_s`,
:math:`\sin(\theta_W)`, :math:`y_e`, :math:`y_\mu`, :math:`y_\tau`,
:math:`y_b`, :math:`y_t`, :math:`v` at the low-energy scale
:math:`Q_{\text{low}}` in the :math:`\overline{\text{MS}}` or
:math:`\overline{\text{DR}}` scheme.

.. note:: The individual loop orders of the threshold corrections can
          be specified using ``FlexibleSUSY[24]``.

- :math:`\alpha_{\text{em}}(Q_{\text{low}})`: If the threshold
  correction loop order is set to ``0``,
  :math:`\alpha_{\text{em}}(Q_{\text{low}})` is set to
  :math:`\alpha_{\text{em}}^{\text{SM}(5)}(Q_{\text{low}})` in the
  Standard Model with 5 active quark flavours.  If the threshold
  correction loop order is set to ``1``,
  :math:`\alpha_{\text{em}}(Q_{\text{low}})` is calculated from
  :math:`\alpha_{\text{em}}^{\text{SM}(5)}(Q_{\text{low}})` using the
  full 1-loop threshold correction.

- :math:`\alpha_s(Q_{\text{low}})`: If the threshold correction loop
  order is set to ``0``, :math:`\alpha_s(Q_{\text{low}})` is set to
  :math:`\alpha_s^{\text{SM}(5)}(Q_{\text{low}})` in the Standard
  Model with 5 active quark flavours.  If the threshold correction
  loop order is set to ``1``, :math:`\alpha_s(Q_{\text{low}})` is
  calculated from :math:`\alpha_s^{\text{SM}(5)}(Q_{\text{low}})`
  using the full 1-loop threshold correction.

- :math:`\sin(\theta_W)(Q_{\text{low}})`: If the threshold correction
  loop order is set to ``0``, the weak mixing angle is calculated from
  either (i) :math:`\{G_F,M_Z\}` or (ii) :math:`\{M_W,M_Z\}`
  (depending on the choice of the weak mixing angle calculation in the
  FlexibleSUSY model file, see `FlexibleSUSY model file`_) using the
  corresponding tree-level relation.

  If the threshold correction loop order is set to ``1``, the the weak
  mixing angle is calculated at the 1-loop level, taking into account

  - (i): complete 1-loop corrections to the W and Z self-energies
    :math:`\Pi_{ZZ}^T, \Pi_{ZZ}^T` as well as 1-loop corrections to
    :math:`\Delta r`, which includes vertex and box contributions
    :math:`\delta_{\text{VG}}` from neutralinos, charginos, selectrons
    and smuons.

  - (ii): complete 1-loop corrections to the W and Z self-energies
    :math:`\Pi_{ZZ}^T, \Pi_{ZZ}^T`.

  If the threshold correction loop order is set to ``2``, the weak
  mixing angle is calculated at the 1-loop level, as above, and the
  following 2-loop correction is taken into account:

  - (i): 2-loop corrections to :math:`\Delta r` of the order
    :math:`O(\alpha_{\text{em}} \alpha_s + y_t^4)` from
    [hep-ph:9606211]_ Eqs. (C.5)-(C.6).

- :math:`y_e(Q_{\text{low}})`, :math:`y_\mu(Q_{\text{low}})`,
  :math:`y_\tau(Q_{\text{low}})`: If the threshold correction loop order
  is set to ``0``, the lepton Yukawa couplings
  :math:`y_e(Q_{\text{low}})`, :math:`y_\mu(Q_{\text{low}})`,
  :math:`y_\tau(Q_{\text{low}})` are calculated from the lepton pole
  masses in the Standard Model with 5 active quark flavours using the
  tree-level relation.

  If the threshold correction loop order is set to ``1``,
  :math:`y_e(Q_{\text{low}})`, :math:`y_\mu(Q_{\text{low}})`,
  :math:`y_\tau(Q_{\text{low}})` are calculated at the scale
  :math:`Q_{\text{low}}` at the 1-loop level from the running lepton
  masses in Standard Model with 5 active quark flavours.

- :math:`y_b(Q_{\text{low}})`: If the threshold correction loop order is
  set to ``0``, the bottom Yukawa couplings :math:`y_b(Q_{\text{low}})` is
  calculated from the running bottom mass in the Standard Model with 5
  active quark flavours, :math:`m_b^{(5)}(Q_{\text{low}})`, using the
  tree-level relation.

  If the threshold correction loop order is set to ``1``,
  :math:`y_b(Q_{\text{low}})` is calculated at the scale
  :math:`Q_{\text{low}}` from :math:`m_b^{(5)}(Q_{\text{low}})` taking the
  complete 1-loop correction into account.

- :math:`y_t(Q_{\text{low}})`: If the threshold correction loop order is
  set to ``0``, the running top Yukawa coupling
  :math:`y_t(Q_{\text{low}})` is calculated from the top pole mass,
  :math:`M_t`, using the tree-level relation.

  If the threshold correction loop order is set to ``1``, the running
  :math:`y_t(Q_{\text{low}})` is calculated at the scale
  :math:`Q_{\text{low}}` from :math:`M_t` taking the complete 1-loop
  correction into account.

  .. math::

    m_t(Q) &= M_t +
    \text{Re\;}\Sigma_{t}^{S}(M_t)
    + M_t
    \left[ \text{Re\;}\Sigma_{t}^{L}(M_t) +
      \text{Re\;}\Sigma_{t}^{R}(M_t) + \Delta
      m_t^{(1),\text{QCD}} \right] ,

  where :math:`\Sigma_{t}^{S}(p)`, :math:`\Sigma_{t}^{L}(p)`,
  :math:`\Sigma_{t}^{R}(p)` denote the scalar, left- and right-handed
  parts of the top self-energy without the gluon contribution.  The
  1-loop SM-QCD contribution :math:`m_t^{(1),\text{QCD}}` reads in the
  :math:`\overline{\text{DR}}` scheme

  .. math::

    \Delta m_t^{(1),\text{QCD}} &=
       -\frac{g_3^2}{12 \pi^2} \left[5-3 \log\left(\frac{m_t^2}{Q^2}\right)\right],

  and in the :math:`\overline{\text{MS}}` scheme

  .. math::

    \Delta m_t^{(1),\text{QCD}} &=
       -\frac{g_3^2}{12 \pi^2} \left[4-3 \log\left(\frac{m_t^2}{Q^2}\right)\right].

  If the threshold correction loop order is set to ``2``,
  2-loop SM-QCD corrections are taken into count as

  .. math::

    m_t(Q) &= M_t +
    \text{Re\;}\Sigma_{t}^{S}(M_t)
    + M_t
    \left[ \text{Re\;}\Sigma_{t}^{L}(M_t) +
      \text{Re\;}\Sigma_{t}^{R}(M_t) + \Delta
      m_t^{(1),\text{QCD}} + \Delta m_t^{(2),\text{QCD}} \right] ,

  where :math:`\Delta m_t^{(2),\text{QCD}}` reads in the
  :math:`\overline{\text{DR}}` scheme [hep-ph:0210258]_

  .. math::

    \Delta m_t^{(2),\text{QCD}} &= \left(\Delta
      m_t^{(1),\text{QCD}}\right)^2
    - \frac{g_3^4}{4608 \pi^4} \Bigg[396
    \log^2\left(\frac{m_t^2}{Q^2}\right)-1476
    \log\left(\frac{m_t^2}{Q^2}\right)
    -48 \zeta(3)+2011+16 \pi^2 (1+\log 4)\Bigg] \,,

  and in the :math:`\overline{\text{MS}}` scheme [hep-ph:9803493]_

  .. math::

    \Delta m_t^{(2),\text{QCD}} &= \left(\Delta
      m_t^{(1),\text{QCD}}\right)^2 - \frac{g_3^4}{4608 \pi^4}
    \Bigg[396 \log^2\left(\frac{m_t^2}{Q^2}\right)
    - 2028 \log\left(\frac{m_t^2}{Q^2}\right)
    - 48 \zeta(3) + 2821 + 16 \pi^2 (1+\log 4)\Bigg] \,.

  If the threshold correction loop order is set to ``3`` in *non-SUSY*
  models, the 3-loop SM-QCD corrections from Refs. [hep-ph:9912391]_,
  [hep-ph:9911434]_ are taken into count as

  .. math::

    m_t(Q) &= M_t +
    \text{Re\;}\Sigma_{t}^{S}(M_t)
    + M_t
    \left[ \text{Re\;}\Sigma_{t}^{L}(M_t) +
      \text{Re\;}\Sigma_{t}^{R}(M_t) + \Delta
      m_t^{(1),\text{QCD}} + \Delta m_t^{(2),\text{QCD}} + \Delta m_t^{(3),\text{QCD}} \right] ,

  where :math:`\Delta m_t^{(3),\text{QCD}}` reads in the
  :math:`\overline{\text{MS}}` scheme

  .. math::

     \Delta m_t^{(3),\text{QCD}} =
     -\frac{g_3^6 \left\{2700 \left[-312 \zeta (3)+1645+8 \pi ^2
        (1+\log (4))\right] \log \left(\frac{Q^2}{m^2}\right)+48600 \log
        ^3\left(\frac{Q^2}{m^2}\right)+714420 \log
        ^2\left(\frac{Q^2}{m^2}\right)-15 \left[69120
        \text{Li}_4\left(\frac{1}{2}\right)+116496 \zeta(3)-94800 \zeta
        (5)-531197+2880 \log^4(2)\right] - 4 \pi^2 [129510 \zeta
        (3)-393101+240 \log(2) (697+24 \log(2))] + 10500 \pi
        ^4\right\}}{9953280 \pi^6}

.. note:: The 1-, 2-, and 3-loop QCD corrections can be found in
          Mathematica form in ``meta/TwoLoopQCD.m`` and
          ``meta/ThreeLoopQCD.m``.

2-loop Higgs pole mass contributions (``FlexibleSUSY[8-11]``)
`````````````````````````````````````````````````````````````

Selects (on/off = 1/0) the individual 2-loop Higgs pole mass
contributions (if available).

Force output (``FlexibleSUSY[12]``)
```````````````````````````````````

If set to 1, an output is always printed, even if a problem has
occurred during the calculation.

.. WARNING:: Be careful with this option!  Check the problems and
             warnings that have occurred!

Top pole mass loop order (``FlexibleSUSY[13]``)
```````````````````````````````````````````````

Loop order of contributions to the top pole mass.  0 = full 1-loop, 1
= 2-loop QCD, 2 = 3-loop QCD.

.. note:: The top pole mass is only calculated if ``FlexibleSUSY[3] = 1``.

Beta-function zero threshold (``FlexibleSUSY[14]``)
```````````````````````````````````````````````````

Below this threshold, beta-functions are treated as being exactly
zero.  Setting this threshold to a non-zero value can avoid numerical
problems / non-convergence problems in models with complex parameters.

Calculate observables (``FlexibleSUSY[15]``)
````````````````````````````````````````````

Enable/disable (1/0) the calculation of the observables specified in
the FlexibleSUSY model file.  See the section on observables in
`FlexibleSUSY model file`_ for further details about how to select the
calculation of observables in FlexibleSUSY.

Force positive Majorana masses (``FlexibleSUSY[16]``)
`````````````````````````````````````````````````````

If set to 1, the masses of Majorana fermions will always be positive.
In this case, the corresponding mixing matrices may be complex.

.. WARNING:: Setting ``FlexibleSUSY[6] = 1`` violates the SLHA standard.

Pole mass scale (``FlexibleSUSY[17]``)
``````````````````````````````````````

Using ``FlexibleSUSY[17]``, the renormalization scale at which the
pole mass spectrum is calculated can be overwritten.  By default the
renormalization scale is the SUSY scale (``SUSYScale`` variable in the
model file).  If ``FlexibleSUSY[17]`` is set to ``0``, the value given
by the ``SUSYScale`` variable is used.  If ``FlexibleSUSY[17]`` is set
to a non-zero value, then this value is used as renormalization scale.

EFT pole mass scale (``FlexibleSUSY[18]``)
``````````````````````````````````````````

.. note:: Only used if ``FlexibleEFTHiggs == True``

Using ``FlexibleSUSY[18]``, the renormalization scale at which the
Standard Model pole mass spectrum is calculated in the EFT can be
overwritten.  If unspecified or set to ``0``, the minimum of the top
pole mass and the ``SUSYScale`` is used.

EFT matching scale (``FlexibleSUSY[19]``)
`````````````````````````````````````````

.. note:: Only used if ``FlexibleEFTHiggs == True``

Using ``FlexibleSUSY[19]``, the renormalization scale at which the full
model is matched to the Standard Model can be overwritten.  If
unspecified or set to ``0``, the ``SUSYScale`` is used.

EFT upwards matching loop order (``FlexibleSUSY[20]``)
``````````````````````````````````````````````````````

.. note:: Only used if ``FlexibleEFTHiggs == True``

Using ``FlexibleSUSY[20]``, the loop order for the matching of the
Standard Model to the full BSM model can be selected ("upwards
matching").  If unspecified, the loop order is set to ``2``.

EFT downwards matching loop order (``FlexibleSUSY[21]``)
````````````````````````````````````````````````````````

.. note:: Only used if ``FlexibleEFTHiggs == True``

Using ``FlexibleSUSY[21]``, the loop order for the matching of the BSM
model to the Standard Model can be selected ("downwards matching").
If unspecified, the loop order is set to ``1``.

EFT index of SM-like Higgs (``FlexibleSUSY[22]``)
`````````````````````````````````````````````````

.. note:: Only used if ``FlexibleEFTHiggs == True``

Using ``FlexibleSUSY[22]``, the user can specify which Higgs in the BSM
model should be interpreted to be the SM-like one.  If unspecified,
the index is set to ``0``, i.e. the lightest Higgs eigenstate in the BSM
model is interpreted as the SM-like Higgs.

Calculate pole masses of BSM particles (``FlexibleSUSY[23]``)
`````````````````````````````````````````````````````````````

Enable/disable (1/0) the calculation of the pole masses of
non-Standard Model particles.

Individual threshold corrections (``FlexibleSUSY[24]``)
```````````````````````````````````````````````````````

The entry ``FlexibleSUSY[24]`` can be used for a fine-grained control to
specify the loop orders of the low-energy threshold corrections of the
SM(5) parameters to the parameters of the BSM model.  The given number
is composed of several digits, each one specifying a threshold
correction loop order of a parameter.  The following table shows which
digit is associated with which parameter.

========================== =========================================== ===========================
 digit position :math:`n`   default value (prefactor of :math:`10^n`)   parameter
========================== =========================================== ===========================
 0                          1 (1-loop)                                  :math:`\alpha_{\text{em}}`
 1                          2 (2-loop)                                  :math:`\sin\theta_W`
 2                          4 (4-loop)                                  :math:`\alpha_{s}`
 3                          1 (1-loop)                                  :math:`m_Z`
 4                          1 (1-loop)                                  :math:`m_W`
 5                          1 (1-loop)                                  :math:`m_h`
 6                          4 (4-loop)                                  :math:`m_t`
 7                          2 (2-loop)                                  :math:`m_b`
 8                          1 (1-loop)                                  :math:`m_{\tau}`
========================== =========================================== ===========================

Note, that the threshold correction loop order of a parameter is not
higher than the "global" threshold correction loop order, specified by
``FlexibleSUSY[7]``.

3-loop corrections to the Higgs pole mass (``FlexibleSUSY[25-29]``)
```````````````````````````````````````````````````````````````````

In the MSSM, the 3-loop corrections to the Higgs pole mass of the
order :math:`O(\alpha_t \alpha_s^2 + \alpha_b \alpha_s^2)`
[1005.5709]_ can be taken into account.  To include them, the variable
``UseHiggs3LoopMSSM`` must be set to ``True`` in the model file::

    UseHiggs3LoopMSSM = True;

.. important:: It is strongly recommended to also set ``UseMSSMYukawa2Loop = True;`` and ``UseMSSM3LoopRGEs = True;`` for consistency.

To enable the 3-loop corrections at run-time in general, set both
``FlexibleSUSY[4]`` and ``FlexibleSUSY[5]`` to ``3``.  To enable the
specific :math:`O(\alpha_t \alpha_s^2)` correction at run-time, set the
flag ``FlexibleSUSY[26]`` to ``1``.  To enable the 3-loop correction
:math:`O(\alpha_b \alpha_s^2)` at run-time, set the flag
``FlexibleSUSY[27]`` to ``1``.

The 3-loop corrections from [1005.5709]_ can be calculated in the
:math:`\overline{DR}'`, :math:`\overline{MDR}'` or H3m scheme.  To use
the :math:`\overline{DR}'` scheme, set ``FlexibleSUSY[25]`` to ``0``.
To use the :math:`\overline{MDR}'` scheme, set ``FlexibleSUSY[25]`` to
``1``.  To use the H3m scheme, set ``FlexibleSUSY[25]`` to ``2``.

We recommend to set the following model file options to enable the
3-loop Higgs pole mass corrections in the MSSM::

    UseHiggs2LoopMSSM = True;      (* enable 2-loop corrections *)
    EffectiveMu = \[Mu];           (* sign convention for MSSM mu parameter *)
    UseMSSM3LoopRGEs = True;       (* enable 3-loop RGEs *)
    UseHiggs3LoopMSSM = True;      (* enable 3-loop corrections *)
    UseMSSMYukawa2Loop = True;     (* enable 2-loop SQCD corrections to yt and yb *)
    UseMSSMAlphaS2Loop = True;     (* enable 2-loop SQCD corrections to alpha_s *)

To run FlexibleSUSY with the 3-loop corrections, we recommend the
settings in the SLHA input::

    Block FlexibleSUSY
        4   3                    # pole mass loop order
        5   3                    # EWSB loop order
        6   3                    # beta-functions loop order
        7   2                    # threshold corrections loop order
        8   1                    # Higgs 2-loop corrections O(alpha_t alpha_s)
        9   1                    # Higgs 2-loop corrections O(alpha_b alpha_s)
       10   1                    # Higgs 2-loop corrections O((alpha_t + alpha_b)^2)
       11   1                    # Higgs 2-loop corrections O(alpha_tau^2)
       24   123111221            # individual threshold correction loop orders
       25   0                    # ren. scheme for Higgs 3L corrections (0 = DR', 1 = MDR', 2 = H3m)
       26   1                    # Higgs 3-loop corrections O(alpha_t alpha_s^2)
       27   1                    # Higgs 3-loop corrections O(alpha_b alpha_s^2)
       28   1                    # Higgs 3-loop corrections O(alpha_t^2 alpha_s)
       29   1                    # Higgs 3-loop corrections O(alpha_t^3)
       30   1                    # Higgs 4-loop corrections O(alpha_t alpha_s^3)

In FlexibleSUSY's Mathematica interface, the following settings should
be used::

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

.. note:: In [1708.05720]_ the individual threshold corrections
          (``FlexibleSUSY[24]``, ``thresholdCorrections``) were set to
          ``123111121``, i.e. the 2-loop SQCD threshold corrections to
          :math:`\alpha_s(M_Z)` have not been taken into account for
          clarity, because they would correspond to a partial 4-loop
          contribution to the light CP-even Higgs pole mass.

Loop library to use (``FlexibleSUSY[31]``)
``````````````````````````````````````````

One can configure ``FlexibleSUSY`` with multiple libraries for one
loop Passarino-Veltman integrals via the following command::

   ./configure --with-loop-libraries=<desired libraries>

Where ``<desired libraries>`` is a list of comma separated names for desired libraries.
Currently the following set is available for usage:

======= =========== =========== ============== =================================
 index   name        library     thread-safety   commentary
======= =========== =========== ============== =================================
  0      softsusy    SOFTSUSY_   yes            default value; always enabled
  1      collier     COLLIER_    no             optional; see **specific** below
  2      looptools   LoopTools_  no             optional; see **specific** below
  3      fflite      FFlite      yes            optional; build in; see **specific** below
======= =========== =========== ============== =================================

If the entry ``FlexibleSUSY[31]`` is absent, the loop functions from
SOFTSUSY are used.  If the entry ``FlexibleSUSY[31]`` is set to
``-1``, FlexibleSUSY uses the value from the environment variable
``FLEXIBLESUSY_LOOP_LIBRARY``.

**COLLIER specific**:
To use the COLLIER_ library and header files from a specific directory configure via::

    COLLIER_DIR=/path/to/COLLIER-x.y.z

     ./configure --with-loop-libraries=collier \
                 --with-collier-incdir=$COLLIER_DIR/modules \
                 --with-collier-libdir=$COLLIER_DIR

**LoopTools specific**:
To use the LoopTools_ library and header files from a specific directory configure via::

    LOOPTOOL_DIR=/path/to/looptools/build

    ./configure --with-loop-libraries=looptools \
                --with-looptools-incdir=$LOOPTOOLS_DIR \
                --with-looptools-libdir=$LOOPTOOLS_DIR

As a replacement of ``--with-loop-libraries=<libraries and looptools>`` one can use::

   ./configure --with-loop-libraries=<libraries> --enable-looptools

**fflite specific**:
To use fflite library one can also (as a replacement of ``--with-loop-libraries=<libraries and fflite>``) use::

    ./configure --with-loop-libraries=<libraries> --enable-fflite


Additional physical input parameters (FlexibleSUSYInput)
--------------------------------------------------------

**Block name**: ``FlexibleSUSYInput``

**Default values**::

    Block FlexibleSUSYInput
        0   0.00729735           # alpha_em(0)
        1   125.09               # Mh pole

**Description**:

The ``FlexibleSUSYInput`` block contains fields for additional known
physical input parameters, which are not contained in a SLHA-compliant
``SMINPUTS`` block.

======= ====================================== ============================== ==================
 index   description                            possible values                default value
======= ====================================== ============================== ==================
  0      alpha_em(0) in the Thompson limit      any positive double            1./137.035999074
  1      SM Higgs pole mass                     any positive double            125.09
======= ====================================== ============================== ==================


MODSEL block (MODSEL)
---------------------

**Block name**: ``MODSEL``

**Default values**::

    Block MODSEL
        6    0     # Quark/Lepton flavour violation
       12    0     # running parameter output scale (GeV)

**Description**:

FlexibleSUSYInput supports the following fields of the ``MODSEL``
block, as defined in SLHA-2:

======= ====================================== ========================================= ===========================
 index   description                            possible values                           default value
======= ====================================== ========================================= ===========================
  6      Quark/Lepton flavour violation         0 (no), 1 (quark), 2 (lepton), 3 (both)   0 (= no flavour violation)
 12      Output scale for running parameters    any positive, non-zero double             0 (= SUSYScale)
======= ====================================== ========================================= ===========================


Output blocks
-------------

In FlexibleSUSY the user can define additional SLHA output blocks.
Please refer to the section on output blocks in `FlexibleSUSY model
file`_ section for more information.


References
----------

.. _`FlexibleSUSY model file`: model_file.rst
.. _LoopTools: http://www.feynarts.de/looptools/
.. _COLLIER: https://collier.hepforge.org/
.. _SOFTSUSY: http://softsusy.hepforge.org

.. [1708.05720] `Eur.Phys.J. C77 (2017) no.12, 814 <https://inspirehep.net/record/1617767>`_ [`arxiv:1708.05720 <https://arxiv.org/abs/1708.05720>`_]
.. [1005.5709] `JHEP 1008 (2010) 104 <https://inspirehep.net/record/856612>`_ [`arxiv:1005.5709 <https://arxiv.org/abs/1005.5709>`_]
.. [hep-ph:9606211] `Nucl.Phys. B491 (1997) 3-67 <https://inspirehep.net/record/419242>`_ [`arxiv:hep-ph/9606211 <https://arxiv.org/abs/hep-ph/9606211>`_]
.. [hep-ph:9803493] `Nucl.Phys. B539 (1999) 671-690 <https://inspirehep.net/record/468752>`_ [`arxiv:hep-ph/9803493 <https://arxiv.org/abs/hep-ph/9803493>`_]
.. [hep-ph:9911434] `Nucl.Phys. B573 (2000) 617-651 <https://inspirehep.net/record/510551>`_ [`arxiv:hep-ph/9911434 <https://arxiv.org/abs/hep-ph/9911434>`_]
.. [hep-ph:9912391] `Phys.Lett. B482 (2000) 99-108 <https://inspirehep.net/record/522686>`_ [`arxiv:hep-ph/9912391 <https://arxiv.org/abs/hep-ph/9912391>`_]
.. [hep-ph:0210258] `Eur.Phys.J. C29 (2003) 87-101 <https://inspirehep.net/record/600038>`_ [`arxiv:hep-ph/0210258 <https://arxiv.org/abs/hep-ph/0210258>`_]
