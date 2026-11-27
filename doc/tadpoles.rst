Tadpoles and electroweak symmetry breaking
==========================================

This note describes how FlexibleSUSY treats the tadpoles: how the
electroweak symmetry breaking (EWSB) conditions are imposed, which
parameters enter the loop calculations, and why the choice of these
parameters matters at small VEVs.  The function names refer to the
generated class ``<M>_mass_eigenstates`` (``models/<M>/``) and its EWSB
solver ``<M>_<solver>_ewsb_solver``.

.. contents::
   :local:


The VEVs are fixed, the tadpoles fix other parameters
-----------------------------------------------------

In FlexibleSUSY the VEVs are parameters of the model: they are input
(e.g. :math:`v_d`, :math:`v_u` from :math:`M_Z` and
:math:`\tan\beta`) or are set by a boundary condition.  The EWSB
conditions are then imposed by fixing other parameters, the
``EWSBOutputParameters`` of the model file (e.g. :math:`\mu` and
:math:`B\mu` in the MSSM, or :math:`m_{H_d}^2`, :math:`m_{H_u}^2`,
:math:`m_S^2` in the GNMSSM-like model of the example below), such
that the VEVs are the minimum of the loop-corrected effective
potential.

As a consequence, the VEVs don't receive tadpole corrections, and no
tadpole diagrams are added to the self-energies: the loop tadpoles
enter the results only through the values of the EWSB output
parameters.


The EWSB conditions
-------------------

For each VEV :math:`v_i` there is one EWSB condition (one tadpole
equation)

.. math::

   T_i(P) - t_i^{(1)} - t_i^{(2)} = 0 ,
   \qquad
   T_i(P) = \frac{\partial V^{(0)}}{\partial v_i} ,

where :math:`P` denotes the set of all model parameters,
:math:`V^{(0)}` is the tree-level potential and :math:`t_i^{(n)}` are
the :math:`n`-loop tadpoles.  In the generated code

* :math:`T_i(P)` is ``get_ewsb_eq_<hh>_<i>()``,
* :math:`t_i^{(1)}` is ``tadpole_<hh>_1loop(i)`` (all one-loop
  diagrams of the model, generated from SARAH),
* :math:`t_i^{(2)}` is ``tadpole_<hh>_2loop()``, the two-loop tadpoles
  from the literature (e.g. :math:`O(\alpha_t\alpha_s)` and
  :math:`O(\alpha_t^2)` in the MSSM and NMSSM), consistent with the
  corresponding two-loop Higgs self-energies,
* the left-hand side is ``tadpole_equations()``.

The loop order of the EWSB conditions is set by ``FlexibleSUSY[5]``.
It should be equal to the loop order of the pole masses
(``FlexibleSUSY[4]``), see below.

The EWSB conditions are solved for the ``EWSBOutputParameters``

* at the SUSY scale (``<M>_susy_scale_constraint``) in each iteration
  of the boundary value problem (two-scale solver), and
* once more at the scale where the pole masses are calculated, just
  before the pole masses (``<M>_spectrum_generator``).

The solver first tries a fixed-point iteration: ``ewsb_step()``
evaluates the analytic tree-level solution of the EWSB conditions with
the loop tadpoles taken as constants, the tadpoles are then
recalculated with the new parameters, and so on.  If this doesn't
converge, the GSL root finders (hybrid and Broyden methods) are
applied to ``tadpole_equations()``.  If no solver converges, the
solution that minimizes the tadpole equations is used and the problem
"no EWSB" is flagged.  The tree-level solution, which is also used to
find the initial guess, is found by Mathematica when the code is
generated, or it is given in the model file
(``TreeLevelEWSBSolution``).


Two parameter points: loop-level and tree-level vacuum
------------------------------------------------------

The solution of the EWSB conditions, denoted :math:`P_\text{loop}`,
does **not** fulfil the tree-level EWSB conditions,
:math:`T_i(P_\text{loop}) = t_i \neq 0`.  The tree-level spectrum at
:math:`P_\text{loop}` is therefore not at a tree-level minimum of the
potential: the Goldstone bosons would get masses different from those
of the corresponding gauge bosons in the 't Hooft-Feynman gauge of the
one-loop self-energies, and scalar masses can be tachyonic just
because a loop correction was absorbed into a parameter.

For this reason the running (tree-level) masses and mixing matrices are
calculated at a second parameter point, :math:`P_\text{tree}`, which
fulfils the tree-level EWSB conditions, :math:`T_i(P_\text{tree}) = 0`.
``calculate_DRbar_masses()`` does this:

1. it saves the values of the parameters listed in
   ``FSSolveEWSBTreeLevelFor`` (see below),
2. it solves the tree-level EWSB conditions for these parameters
   (``solve_ewsb_tree_level_custom()``),
3. it calculates all running masses and mixing matrices,
4. it restores the saved parameters, so that the model is at
   :math:`P_\text{loop}` again.

The model therefore holds the parameters :math:`P_\text{loop}`, which
are written to the SLHA output, together with the running masses and
mixings of :math:`P_\text{tree}`.  All loop calculations use the latter:
the self-energies, the loop tadpoles themselves (in each step of the
EWSB iteration), and the threshold corrections at the low scale.  A
tachyonic running mass at :math:`P_\text{tree}` is flagged as a
"running tachyon".


Pole masses
-----------

The pole masses are the solutions of

.. math::

   \det\left[p^2 - M^2_\text{tree}(P_\text{loop})
             + \Sigma(p^2; P_\text{tree})\right] = 0 ,

(``calculate_M<X>_pole()``), where the tree-level mass matrix
``get_mass_matrix_<X>()`` is evaluated with the loop-level EWSB
solution, while the self-energy :math:`\Sigma` is evaluated with the
running masses and couplings of :math:`P_\text{tree}`.  The tadpoles
are contained in the difference between the two points.

**Example: soft Higgs masses.**  Let :math:`X = m_{H_i}^2` be the
parameter that is re-solved in the EWSB condition of :math:`v_i`, where
it appears as :math:`v_i X`.  From :math:`T_i(P_\text{loop}) = t_i` and
:math:`T_i(P_\text{tree}) = 0` one finds

.. math::

   X_\text{tree} = X_\text{loop} - \frac{t_i}{v_i} .

As :math:`m_{H_i}^2` appears only on the diagonal of the CP-even,
CP-odd and charged Higgs mass matrices,

.. math::

   M^2_\text{tree}(P_\text{loop}) = M^2_\text{tree}(P_\text{tree})
   + \operatorname{diag}\left(\frac{t_i}{v_i}\right) .

This is equivalent to starting from the tree-level vacuum and adding
the tadpole contributions :math:`t_i/v_i` to the diagonal of the
self-energies (the default choice was made to agree with SOFTSUSY).

This also shows why the EWSB loop order (``FlexibleSUSY[5]``) should
equal the pole mass loop order (``FlexibleSUSY[4]``): the
:math:`n`-loop tadpoles enter the pole masses only through
:math:`P_\text{loop}`.  If they are missing there, the :math:`n`-loop
self-energies are not accompanied by the corresponding tadpoles.

**General case.**  If the parameter :math:`X` appears as :math:`c_i X`
in the EWSB condition of :math:`v_i`, then
:math:`X_\text{tree} = X_\text{loop} - t_i/c_i` and the running masses in
the loops are shifted by

.. math::

   \delta m^2 = -\frac{\partial M^2_\text{tree}}{\partial X}\,
   \frac{t_i}{c_i}

relative to the masses at :math:`P_\text{loop}`.  This shift is of
one-loop order, so it changes the pole masses only at two-loop order,
beyond the accuracy of a one-loop calculation: every choice of
:math:`X` is correct at one loop.  The size of the two-loop term,
however, depends on the choice.  It is unusually large if
:math:`t_i/c_i` is large, e.g. if :math:`c_i` is a small VEV.

Note that the two points can't be merged by re-solving "the other way
round": :math:`P_\text{loop}` doesn't solve the tree-level conditions,
so there are many tree-level points near it, one for each parameter
that absorbs the tadpole :math:`t_i`, and they differ in their masses.


Which parameters are re-solved at tree level
--------------------------------------------

The parameters re-solved in ``calculate_DRbar_masses()`` are given by
``FSSolveEWSBTreeLevelFor`` in the model file (see the `model file
documentation <model_file.rst>`_).  They are **not** the
``EWSBOutputParameters`` unless the model file says so.  If
``FSSolveEWSBTreeLevelFor`` is empty (the default), they are chosen
automatically (``WriteModelClass`` in ``meta/FlexibleSUSY.m``):

1. In SUSY models, the soft-breaking scalar masses that appear in
   exactly one EWSB condition; in non-SUSY models, the real parameters
   of mass dimension 2 that appear in exactly one EWSB condition.
   Soft masses enter the EWSB conditions linearly, so the tree-level
   solution is unique and always exists, and they change only the
   scalar mass matrices.  The EWSB output parameters are often less
   suited: e.g. :math:`\mu` is solved from :math:`|\mu|^2` with a sign
   choice, may have no real tree-level solution, and changes the
   chargino, neutralino and sfermion masses inside the loops.

A parameter that enters an EWSB condition linearly but none of the
tree-level mass matrices, such as the linear soft-breaking term
:math:`\xi_S` (``L[L1]``) of a gauge singlet, is a special choice:
:math:`\delta m^2 = 0`, so the running masses at :math:`P_\text{tree}`
are exactly the tree-level masses at the loop-level vacuum, and the
two-loop ambiguity of the previous section is absent in this sector.
It has to be chosen in the model file, e.g. ::

    FSSolveEWSBTreeLevelFor = { mHd2, mHu2, L[L1] };

This is not always better than the soft mass.  If the loop corrections
are what makes a direction of the potential stable, the tree-level
mass at the loop-level vacuum is tachyonic.  Re-solving the soft mass
absorbs the tadpole into the curvature and hides this; with the linear
term the running mass is tachyonic and the point is flagged as
"running tachyon".  E.g. the example point of DiracGauginos
(:math:`v_S = -400` GeV) is valid with the default choice
(:math:`m_{H_d}^2`, :math:`m_{H_u}^2`, :math:`m_S^2`, :math:`m_T^2`), but
has a running tachyon with :math:`\xi_S`; its pole masses differ by up
to 1.4 GeV between the two choices.  For this reason the linear term
is not chosen automatically.

Examples:

========== ================================ ========================================
model      ``EWSBOutputParameters``         re-solved for the running masses
========== ================================ ========================================
MSSM       :math:`B\mu`, :math:`\mu`        :math:`m_{H_d}^2`, :math:`m_{H_u}^2`
NMSSM      :math:`\kappa`, :math:`v_S`,     :math:`m_{H_d}^2`, :math:`m_{H_u}^2`,
           :math:`m_S^2`                    :math:`m_S^2`
SMSSM,     :math:`\mu`, :math:`B\mu`,       :math:`m_{H_d}^2`, :math:`m_{H_u}^2`,
NUTSMSSM   :math:`\xi_S` (SMSSM)            :math:`\xi_S` (set in the model files)
MRSSM2     :math:`m_{H_d}^2`,               same (``FSSolveEWSBTreeLevelFor =
           :math:`m_{H_u}^2`, :math:`v_S`,  EWSBOutputParameters``)
           :math:`v_T`
SM         :math:`\mu^2`                    :math:`\mu^2`
========== ================================ ========================================


Small VEVs and linear terms
---------------------------

In the SMSSM, the GNMSSM and other models with linear terms of a gauge
singlet :math:`S`, the singlet mass :math:`m_S^2` enters the EWSB
condition of :math:`v_S` as :math:`v_S m_S^2`.  If it is re-solved
(the automatic choice), the singlet masses inside the loops are
shifted by :math:`t_S/v_S`.  Unlike in the MSSM, where the
VEVs are of the order of the electroweak scale, :math:`v_S` can be
small or negative at physical points of these models.  Moreover, the
tadpole :math:`t_S` doesn't vanish for :math:`v_S \to 0`, as the linear
terms don't leave :math:`S = 0` invariant.

These models are invariant under the shift :math:`S \to S + \delta` of
the singlet superfield (with :math:`\mu \to \mu + \lambda\delta`,
:math:`M_S \to M_S + 2\kappa\delta`, :math:`L_1 \to L_1 + M_S\delta +
\kappa\delta^2`, the corresponding soft-breaking parameters, and
:math:`v_S \to v_S - \sqrt{2}\delta`), which changes :math:`v_S` but
not the physics.  The tadpole :math:`t_S` is invariant, but
:math:`v_S` is not, so with :math:`m_S^2` re-solved the result depended
on the frame: at a GNMSSM point with :math:`v_S = 611` GeV, shifting to
:math:`v_S = 32` GeV moved a singlet-like Higgs mass from 230.7 to
221.6 GeV, and for :math:`v_S \lesssim 10` GeV or :math:`v_S < 0` the
running singlet masses became tachyonic.  With :math:`\xi_S`
re-solved, the masses are the same in every frame (to
:math:`3\cdot 10^{-5}` at that point, down to :math:`|v_S| = 10` GeV,
and for :math:`v_S < 0`).  Therefore the SMSSM and NUTSMSSM model files
set ``FSSolveEWSBTreeLevelFor = { mHd2, mHu2, L[L1] }``, and models with
a linear singlet term in which :math:`v_S` can be small or negative
should do the same.  ``test/test_SMSSM_singlet_shift.cpp`` checks the
invariance.

**Small** :math:`|v_S|` **with** :math:`m_S^2` **as EWSB output
parameter.**  This choice of ``FSSolveEWSBTreeLevelFor`` doesn't help
if :math:`m_S^2` is itself an EWSB output parameter (as in NUTSMSSM),
because then the loop-level EWSB conditions are ill-conditioned for
small :math:`|v_S|`.  :math:`m_S^2` enters the condition of :math:`v_S`
only as :math:`v_S m_S^2`, whose slope :math:`v_S` becomes comparable
to the :math:`m_S^2` dependence of the one-loop tadpole :math:`t_S`
(through the singlet masses in the loops).  The condition

.. math::

   v_S m_S^2 + \ldots - t_S(m_S^2) = 0

can then have two solutions for the same input parameters, and every
solver (the fixed-point iteration and the GSL root finders) finds the
one nearer its starting point.  At :math:`v_S = 0`, :math:`m_S^2` drops
out of the condition altogether.  Both solutions are genuine one-loop
vacua, so no problem is flagged.  Note that the shift
:math:`S \to S + \delta` also shifts the linear soft-breaking term by
:math:`m_S^2 \delta`, so with the linear term fixed, solving for
:math:`m_S^2` is a different problem in each frame.

E.g. at the GNMSSM point above, shifted to :math:`v_S = 3` GeV, the
one-loop EWSB conditions have the solutions
:math:`m_S^2 = -2.02\cdot 10^4~\text{GeV}^2` (the point at
:math:`v_S = 611` GeV) and :math:`m_S^2 = +1.25\cdot 10^5~\text{GeV}^2`
(a CP-even singlet at 442 instead of 231 GeV).  The spectrum generator
starts from the tree-level solution and finds the second one.  It
finds similarly wrong points for :math:`v_S = 4, \ldots, 1` GeV, and
flags running tachyons for :math:`|v_S| \lesssim 1` GeV, while
:math:`v_S \ge 5` GeV and :math:`v_S \le -3` GeV are fine at this point.
To avoid this, choose a frame with a larger :math:`|v_S|` (the physics
is the same), or a linear term (:math:`\xi_S`, whose coefficient in the
condition is constant) as EWSB output parameter, as the SMSSM does.


Decays
------

The decays use the class ``<M>_mass_eigenstates_decoupling_scheme``.
It copies the parameters :math:`P_\text{loop}` (run to the scale of the
decaying particle), replaces the SM-like parameters (gauge couplings,
VEVs and Yukawa couplings) by those obtained from the Standard Model in
the decoupling scheme, and then solves the **tree-level** EWSB
conditions for the same parameters as ``calculate_DRbar_masses()``.  No loop
tadpoles enter there.  The change of the SM-like parameters is
absorbed by the re-solved parameters, i.e. divided by the coefficient
:math:`c_i` of the EWSB condition, which again is large for a small VEV
if a soft mass is re-solved.


Comparing with other codes
--------------------------

Different codes may solve the EWSB conditions for different parameters
(e.g. SPheno for :math:`m_S^2`, :math:`m_T^2` in the MRSSM, where
FlexibleSUSY solves for :math:`v_S`, :math:`v_T`) and use different
parameter points inside the loops.  The results then differ at two-loop
order, by an amount that grows like :math:`t/v` for small VEVs.  Compare
at the same Lagrangian parameters, and expect such differences when a
small VEV multiplies an EWSB output parameter.
