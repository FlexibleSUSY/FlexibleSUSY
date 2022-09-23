.. sectnum::

FlexibleSUSY model file
=======================

.. contents:: Table of Contents

Model information
-----------------

**Symbol**: ``FSModelName``

**Default value**: *unset*

**Description**:

Name of the model class within the generated code.  If ``FSModelName``
is set to the string ``"@CLASSNAME@"``, it will be replaced by the
``createmodel`` script to the name of the FlexibleSUSY model given
during the ``./createmodel --name=<model_name>`` command.

_____________________________________________________________________

**Symbol**: ``FSEigenstates``

**Default value**: ``SARAH`EWSB``

**Description**:

The name of the particle eigenstates in SARAH.  Usually,
````SARAH``EWSB```` corresponds to the mass eigenstates after breaking of
the electroweak symmetry.

_____________________________________________________________________

**Symbol**: ``FSDefaultSARAHModel``

**Default value**: *unset*

**Description**:

Name of the SARAH model to be used.  A sub-model can be specified
after a ``/``.

Example: In the constrained CP-conserving MSSM (``CMSSM``) the SARAH
model ``MSSM`` is used as::

    FSDefaultSARAHModel = MSSM;

Example: In the constrained CP-violating MSSM (``CMSSMCPV``) the SARAH
model ``MSSM`` together with the sub-model ``CPV`` is used::

    FSDefaultSARAHModel = MSSM/CPV;

_____________________________________________________________________

**Symbol**: ``FSBVPSolvers``

**Default value**: ``{ TwoScaleSolver }``

**Description**:

A list of algorithms to use for solving the boundary value problem.
One or both of ``TwoScaleSolver`` or ``SemiAnalyticSolver`` may be
specified in the list.

Input parameters
----------------

**Symbol**: ``MINPAR``

**Default value**: ``{}``

**Description**:

In the ``MINPAR`` variable a list of input parameters for the spectrum
generator can be given, which is read of the ``MINPAR`` block of the
SLHA input file.

``MINPAR`` is supposed to contain a list.  The list elements are
two-component lists, where the first in element is an integer number
representing the index inside the ``MINPAR`` block.  The second element
is the input parameter.  The input parameter must be either a symbol
or a sign of the form ``Sign[p]``, where ``p`` is the name of a model
parameter.

**Example**: In the CMSSM the ``MINPAR`` block has the form::

    MINPAR = {
        {1, m0},
        {2, m12},
        {3, TanBeta},
        {4, Sign[\[Mu]]},
        {5, Azero}
    };

In this case the input parameters can be given in the SLHA input file
as::

    Block MINPAR                 # Input parameters
        1   125                  # m0
        2   500                  # m12
        3   10                   # TanBeta
        4   1                    # SignMu
        5   0                    # Azero

.. note:: Unspecified parameters are assumed to be zero.

_____________________________________________________________________

**Symbol**: ``EXTPAR``

**Default value**: ``{}``

**Description**:

The ``EXTPAR`` variable is a list of input parameters for the spectrum
generator, which is read of the ``EXTPAR`` block of the SLHA input file.
The list assigned to the ``EXTPAR`` variable must have the same form as
the ``MINPAR`` variable.

**Example**: In the NUTNMSSM the ``EXTPAR`` block has the form::

    EXTPAR = {
        {61, LambdaInput},
        {62, KappaInput},
        {63, ALambdaInput},
        {64, AKappaInput},
        {65, MuEff}
    };

In this case the input parameters can be given in the SLHA input file
as::

    Block EXTPAR                 # Input parameters
       61   0.650                # LambdaInput
       62   0.164                # KappaInput
       63   763.8                # ALambdaInput
       64   1268.2               # AKappaInput
       65   265.2                # MuEff

.. note:: Unspecified parameters are assumed to be zero.

_____________________________________________________________________

**Symbol**: ``IMMINPAR``

**Default value**: ``{}``

**Description**:

The ``IMMINPAR`` variable is a list of input parameters for the spectrum
generator, which is read of the ``IMMINPAR`` block of the SLHA input
file.  The list assigned to the ``IMMINPAR`` variable must have the same
form as the ``MINPAR`` variable.

**Example**: In the CP-violating MSSM (``CMSSMCPV``) the ``IMMINPAR`` block
has the form::

    IMMINPAR = {
        {2, Imm12},
        {5, ImAzero}
    };

In this case the input parameters can be given in the SLHA input file
as::

    Block IMMINPAR
        2   10                   # Imm12
        5   10                   # ImAzero

.. note:: Unspecified parameters are assumed to be zero.

_____________________________________________________________________

**Symbol**: ``IMEXTPAR``

**Default value**: ``{}``

**Description**:

The ``IMEXTPAR`` variable is a list of input parameters for the spectrum
generator, which is read of the ``IMEXTPAR`` block of the SLHA input
file.  The list assigned to the ``IMEXTPAR`` variable must have the same
form as the ``MINPAR`` variable.

**Example**: In the CP-violating MSSM (``MSSMCPV``) the ``IMEXTPAR`` block
has the form::

    IMEXTPAR = {
        {1, ImM1Input},
        {2, ImM2Input},
        {3, ImM3Input},
        {23, ImMuInput}
    };

In this case the input parameters can be given in the SLHA input file
as::

    Block IMEXTPAR
        1    100                 # Im(M1(MSUSY))
        2    100                 # Im(M2(MSUSY))
        3    100                 # Im(M3(MSUSY))
       23    100                 # Im(Mu(MSUSY))

.. note:: Unspecified parameters are assumed to be zero.

_____________________________________________________________________

**Symbol**: ``FSAuxiliaryParameterInfo``

**Default value**: ``{}``

**Description**:

In the ``FSAuxiliaryParameterInfo`` variable additional input or extra
parameters can be defined, and extra information provided can be
provided about existing input parameters.  ``FSAuxiliaryParameterInfo``
is expected to be a list, whose element are two-component lists.  The
first element of this list is a symbol representing the parameter.
The second element is a list of properties for that parameter,
specified as replacement rules.  The supported properties are

 - ``InputParameter``: A value of ``True`` or ``False`` indicating if the
   parameter is an input parameter.
 - ``LesHouches``: The name of the SLHA block from which the
   parameter should be read, if it is an input parameter.
 - ``MassDimension``: A number specifying the mass dimension of the
   parameter.
 - ``ParameterDimensions``: A list specifying the vector- or
   matrix-type of the input parameter.  A list of the form ``{N,M}``
   with ``N`` and ``M`` being integer numbers defines a NxM matrix.  A
   list of the form ``{N}``, with ``N`` > 1 defines a vector with ``N``
   rows.  A list of the form ``{1}`` or ``{}`` defines a scalar.

**Example**: In the MSSM the ``FSAuxiliaryParameterInfo`` variable has
the form::

    FSAuxiliaryParameterInfo = {
        {Aeij, { LesHouches -> AeijIN,
                 ParameterDimensions -> {3,3},
                 InputParameter -> True
               } },
        {Adij, { LesHouches -> AdijIN,
                 ParameterDimensions -> {3,3},
                 InputParameter -> True
               } },
        {Auij, { LesHouches -> AuijIN,
                 ParameterDimensions -> {3,3},
                 InputParameter -> True
               } }
    };

Here, three 3x3 matrix-valued parameters are specified: ``Aeij``,
``Adij`` and ``Auij``.  They are defined as input parameters.  These
matrices are read from the blocks ``AeijIN``, ``AdijIN`` and ``AuijIN``,
respectively.

These input parameters can be given in the SLHA input file as::

    Block AeijIN
        1   1   100
        1   2   100
        1   3   100
        2   1   100
        2   2   100
        2   3   100
        3   1   100
        3   2   100
        3   3   100
    Block AdijIN
        1   1   200
        1   2   200
        1   3   200
        2   1   200
        2   2   200
        2   3   200
        3   1   200
        3   2   200
        3   3   200
    Block AuijIN
        1   1   300
        1   2   300
        1   3   300
        2   1   300
        2   2   300
        2   3   300
        3   1   300
        3   2   300
        3   3   300

.. note:: Unspecified parameters are assumed to be zero.

_____________________________________________________________________

**Symbol**: ``RealParameters``

**Default value**: ``{ All }``

**Description**:

``RealParameters`` is a list, which contains the names of all model
parameters, which should be treated as real parameters.  By default,
``RealParameters`` is set to ``{ All }``, meaning that by default all
paramerters are treated to be real.  If ``RealParameters`` is set to the
empty list ``{}``, FlexibleSUSY takes the information which paramerters
are real and which are complex from the SARAH model file.

Example: In the complex Standard Model (``cSM``), the parameters ``mu2``
and ``\[Lambda]`` should be defined to be real::

    RealParameters = { mu2, \[Lambda] };

Note: The gauge couplings and VEVs are always assumed to be real in
SARAH.

Example: In the CP-violating MSSM (``CMSSMCPV``), the ``B[\[Mu]]``
parameter should be defined to be real::

    RealParameters = { B[\[Mu]] };

Boundary conditions
-------------------

In FlexibleSUSY, spectrum generators with maximum 3 boundary
conditions can be generated.  These boundary conditions are named
"high-scale", "susy-scale" and "low-scale" boundary condition and are
described in the following.

However, it is possible to disable the high-scale boundary condition.
In order to do so, set::

    OnlyLowEnergyFlexibleSUSY = True;  (* disable high-scale BC, default: False *)

_____________________________________________________________________

**Symbol**: ``LowScale``

**Default value**: *unset*

**Description**:

The scale of the low-scale boundary condition, at which the model is
matched to the Standard Model.

.. note:: ``LowScale`` is ignored if ``FlexibleEFTHiggs == True``

Example: In the CMSSM the low-energy scale should be set to the Z or
top pole mass.  This choice is achieved by the following expression::

    LowScale = LowEnergyConstant[MZ];

_____________________________________________________________________

**Symbol**: ``LowScaleFirstGuess``

**Default value**: *unset*

**Description**:

First guess of the low-energy scale.

.. note:: ``LowScaleFirstGuess`` is ignored if ``FlexibleEFTHiggs == True``

Example: In the CMSSM the first guess for the low-energy scale should
be set to the Z or top pole mass::

    LowScaleFirstGuess = LowEnergyConstant[MZ];

_____________________________________________________________________

**Symbol**: ``LowScaleInput``

**Default value**: ``{}``

**Description**:

With the ``LowScaleInput`` variable boundary conditions at the
low-energy scale can be specified.  ``LowScaleInput`` is a list.  Please
refer to \ref input_format for details about the list format.

At the low-energy scale, FlexibleSUSY automatically determines the
three gauge couplings from the SLHA input parameters
:math:`\alpha_{em}`, :math:`M_Z` and :math:`G_F` or :math:`M_W`.

.. note:: ``LowScaleInput`` is ignored if ``FlexibleEFTHiggs == True``

Example: In the CMSSM ``LowScaleInput`` is given as follows::

    LowScaleInput = {
       {Yu, Automatic},
       {Yd, Automatic},
       {Ye, Automatic},
       {vd, 2 MZDRbar / Sqrt[GUTNormalization[g1]^2 g1^2 + g2^2] Cos[ArcTan[TanBeta]]},
       {vu, 2 MZDRbar / Sqrt[GUTNormalization[g1]^2 g1^2 + g2^2] Sin[ArcTan[TanBeta]]}
    };

The method to determine the weak mixing angle can be chosen by setting
the variable ``FSWeakMixingAngleInput`` to either ``Automatic``,
``FSFermiConstant`` or ``FSMassW``.  ``FSWeakMixingAngleInput`` is set to
``Automatic`` by default.

====================================== =======================================================
 Value of ``FSWeakMixingAngleInput``    Parameters from which weak mixing angle is determined  
====================================== =======================================================
 ``FSFermiConstant``                    :math:`G_F` and :math:`M_Z`                                
 ``FSMassW``                            :math:`M_W` and :math:`M_Z`                                
 ``Automatic`` (default) (recommended)  chose most precise method automatically                
====================================== =======================================================

Example: Automatically chose most precise method to determine the weak
mixing angle::

    FSWeakMixingAngleInput = Automatic; (* recommended *)

.. note:: If ``FSWeakMixingAngleInput = FSMassW;`` is chosen,
          FlexibleSUSY looks for the definition of the weak mixing
          angle in the symbol ``SARAH`Weinberg``.  If
          ``SARAH`Weinberg`` is not defined, FlexibleSUSY uses the
          expression assigned to ``FSWeakMixingAngleExpr``, which is
          by default set to
          ``ArcSin[Sqrt[1-Mass[SARAH`VectorW]^2/Mass[SARAH`VectorZ]^2]]``.

_____________________________________________________________________

**Symbol**: ``SUSYScale``

**Default value**: *unset*

**Description**:

The scale of the susy-scale boundary condition, which is defined to be
between the low-scale and the high-scale.  This is the scale at which
the electroweak symmetry breaking conditions are imposed by default,
see \ref input_format.

Example: In the CMSSM the SUSY scale should be set to the geometric
average of the two stop masses.  This choice is achieved by the
following expression::

    SUSYScale = Sqrt[Product[M[Su[i]]^(Abs[ZU[i,3]]^2 + Abs[ZU[i,6]]^2), {i,6}]];

_____________________________________________________________________

**Symbol**: ``SUSYScaleFirstGuess``

**Default value**: *unset*

**Description**:

First guess of the SUSY scale.

Example: In the CMSSM a reasonable first guess for the SUSY scale can
be given by the following combination of the mSUGRA parameters::

    SUSYScaleFirstGuess = Sqrt[m0^2 + 4 m12^2];

_____________________________________________________________________

**Symbol**: ``SUSYScaleInput``

**Default value**: ``{}``

**Description**:

With the ``SUSYScaleInput`` variable boundary conditions at the SUSY
scale can be specified.  ``SUSYScaleInput`` is a list.  Please refer to
\ref input_format for details about the list format.

Example: In the NUTNMSSM ``SUSYScaleInput`` is given as follows::

    SUSYScaleInput = {
       {\[Lambda], LambdaInput},
       {\[Kappa], KappaInput},
       {vS, Sqrt[2] MuEff / LambdaInput}
    };

_____________________________________________________________________

**Symbol**: ``HighScale``

**Default value**: *unset*

**Description**:

This is the scale of the high-scale boundary condition.

Example: In the CMSSM the high-energy scale, :math:`M_X`, is given by
the equality of the gauge couplings :math:`g_1(M_X)` and :math:`g_2(M_X)`::

    HighScale = g1 == g2;

_____________________________________________________________________

**Symbol**: ``HighScaleFirstGuess``

**Default value**: *unset*

**Description**:

First guess of the high-energy scale.

Example: In the CMSSM a reasonable initial guess for the high-energy
scale is::

    HighScaleFirstGuess = 2.0 10^16;

_____________________________________________________________________

**Symbol**: ``HighScaleMinimum``

**Default value**: *unset*

**Description**:

Minimum value of the high-energy scale during the iteration.

Example: In the E6SSM the high-energy scale can vary a lot between the
iteration steps.  For this reason, it makes sense to use a minimum
high-energy scale in intermediate steps as::

    HighScaleMinimum = 1.0 10^4;

_____________________________________________________________________

**Symbol**: ``HighScaleMaximum``

**Default value**: *unset*

**Description**:

Maximum value of the high-energy scale during the iteration.

Example: In the E6SSM the high-energy scale can vary a lot between the
iteration steps.  For this reason, it makes sense to use a maximum
high-energy scale in intermediate steps as::

    HighScaleMaximum = 5.0 10^17;

_____________________________________________________________________

**Symbol**: ``HighScaleInput``

**Default value**: ``{}``

**Description**:

With the ``HighScaleInput`` variable boundary conditions at the
high-energy scale can be specified.  ``HighScaleInput`` is a list.
Please refer to \ref input_format for details about the list format.

Example: In the CMSSM ``HighScaleInput`` is set to the mSUGRA boundary
conditions::

    HighScaleInput = {
       {T[Ye], Azero Ye},
       {T[Yd], Azero Yd},
       {T[Yu], Azero Yu},
       {mHd2, m0^2},
       {mHu2, m0^2},
       {mq2, UNITMATRIX[3] m0^2},
       {ml2, UNITMATRIX[3] m0^2},
       {md2, UNITMATRIX[3] m0^2},
       {mu2, UNITMATRIX[3] m0^2},
       {me2, UNITMATRIX[3] m0^2},
       {MassB, m12},
       {MassWB,m12},
       {MassG, m12}
    };

_____________________________________________________________________

**Symbol**: ``InitialGuessAtLowScale``

**Default value**: ``{}``

**Description**:

With the ``InitialGuessAtLowScale`` variable initial values for the
model MS-bar/DR-bar parameters can be given at the low-energy scale
``LowScale``.

.. note:: ``InitialGuessAtLowScale`` is ignored if ``FlexibleEFTHiggs == True``

Example: In the CMSSM ``InitialGuessAtLowScale`` is given as follows::

    InitialGuessAtLowScale = {
       {vd, LowEnergyConstant[vev] Cos[ArcTan[TanBeta]]},
       {vu, LowEnergyConstant[vev] Sin[ArcTan[TanBeta]]},
       {Yu, Automatic},
       {Yd, Automatic},
       {Ye, Automatic}
    };

_____________________________________________________________________

**Symbol**: ``InitialGuessAtSUSYScale``

**Default value**: ``{}``

**Description**:

.. note:: ``InitialGuessAtSUSYScale`` is only used if ``FlexibleEFTHiggs == True``

With the ``InitialGuessAtSUSYScale`` variable initial values for the
model MS-bar/DR-bar parameters can be given at the SUSY scale
``SUSYScale``.

Example: In the MSSMEFTHiggs ``InitialGuessAtSUSYScale`` is given as follows::

    InitialGuessAtSUSYScale = {
        {Yu, Automatic},
        {Yd, Automatic},
        {Ye, Automatic}
        {MassB, Ms},
        {MassWB, Ms},
        {MassG, Ms},
        {mq2, UNITMATRIX[3] Ms^2},
        {mu2, UNITMATRIX[3] Ms^2},
        {md2, UNITMATRIX[3] Ms^2},
        {ml2, UNITMATRIX[3] Ms^2},
        {me2, UNITMATRIX[3] Ms^2},
        {\[Mu], Ms},
        {B[\[Mu]], Sqr[Ms]/(TanBeta + 1/TanBeta)},
        {T[Yu], Ms/TanBeta Yu},
        {T[Yd], Ms TanBeta Yd},
        {T[Ye], Ms TanBeta Ye},
        {T[Yu][3,3], (Ms/TanBeta + Xtt Ms) Yu[3,3]}
    };

_____________________________________________________________________

**Symbol**: ``InitialGuessAtHighScale``

**Default value**: ``{}``

**Description**:

With the ``InitialGuessAtHighScale`` variable initial values for the
model MS-bar/DR-bar parameters can be given at the high-energy scale
``HighScale``.

Example: In the CMSSM ``InitialGuessAtHighScale`` is given as
follows::

    InitialGuessAtHighScale = {
       {\[Mu]   , 1.0},
       {B[\[Mu]], 0.0}
    };

_____________________________________________________________________

**Symbol**: ``EWSBOutputParameters``

**Default value**: ``{}``

**Description**:

In the ``EWSBOutputParameters`` variable the model parameters must be
specified, which are fixed by the electroweak symmetry breaking (EWSB)
conditions, :math:`\partial V_\text{Higgs}/\partial v_i = 0`.  The
length of the ``EWSBOutputParameters`` list must be equal to the number
of EWSB conditions.

Example: In the CMSSM ``EWSBOutputParameters`` is given as follows::

    EWSBOutputParameters = { B[\[Mu]], \[Mu] };

The elements of the ``EWSBOutputParameters`` must be _real_ parameters.
In a model with complex parameters, as in the CMSSMCPV for example,
``EWSBOutputParameters`` is set to be::

    EWSBOutputParameters = { Re[B[\[Mu]]], Im[B[\[Mu]]], \[Mu] };

_____________________________________________________________________

**Symbol**: ``EWSBInitialGuess``

**Default value**: ``{}``

**Description**:

In the ``EWSBInitialGuess`` variable initial guesses for some or all
of the EWSB output parameters can be specified.

Example: In the VCMSSM ``EWSBInitialGuess`` is defined as::

    EWSBInitialGuess = {
       {TanBeta, vu / vd},
       {MuSq, \[Mu]^2}
    };

_____________________________________________________________________

**Symbol**: ``EWSBSubstitutions``

**Default value**: ``{}``

**Description**:

In the ``EWSBSubstitutions`` variable, substitutions for model
parameters in terms of other parameters can be given.
``EWSBSubstitutions`` should be a list of two-component lists, in which
the first element is the parameter to be substituted for, and the
second element is the expression to be substituted in its place.

Example: In the VCMSSM ``EWSBSubstitutions`` is defined as::

    EWSBSubstitutions = {
       {vd, vMSSM Cos[ArcTan[TanBeta]]},
       {vu, vMSSM Sin[ArcTan[TanBeta]]},
       {\[Mu], Sign[\[Mu]] Sqrt[MuSq]}
    };

_____________________________________________________________________

**Symbol**: ``FSSolveEWSBTreeLevelFor``

**Default value**: ``{}``

**Description**:

In the ``FSSolveEWSBTreeLevelFor`` variable the model parameters can be
specified, which are fixed by the tree-level electroweak symmetry
breaking (EWSB) conditions when the running (tree-level) masses are
calculated.  The length of the ``FSSolveEWSBTreeLevelFor`` list must be
either zero (default) or equal to the number of EWSB conditions.  If
``FSSolveEWSBTreeLevelFor`` is the empty list, then the temporary EWSB
output parameters are chosen automatically as follows:

- In SUSY models, by default the soft-breaking squared Higgs mass
  parameters are fixed by the tree-level EWSB equation temporarily
  when the running (tree-level) masses are calculated.

- In non-SUSY models, by default the parameters given in
  ``EWSBOutputParameters`` are fixed by the tree-level EWSB equation
  temporarily when the running (tree-level) masses are calculated.

_____________________________________________________________________

**Symbol**: ``MatchingScaleInput``

**Default value**: ``{}``

**Description**:

.. note:: ``MatchingScaleInput`` is only used if ``FlexibleEFTHiggs == True``

In the ``MatchingScaleInput`` variable, relations between the parameters
of the full model and the Standard Model (the EFT) at the ``SUSYScale``
can be specified.

An important application is the relation between the vacuum
expectation values (VEVs) in a SUSY model and :math:`v` in the Standard
Model: In ``FlexibleEFTHiggs`` the running Yukawa couplings of the full
model are determined from a pole mass matching of the Standard Model
fermions (which need to be present in both models).  For this
determination the running VEVs of the full model must be known and
non-zero.  ``MatchingScaleInput`` allows the user for example to fix the
running VEVs of the full model as a function of the running SM-like
VEV :math:`v` in the full model.

Example: In the MSSM the vacuum expectation values :math:`v_u` and
:math:`v_d` are related to the MSSM SM-like VEV :math:`v = \sqrt{v_u^2 +
v_d^2}` as

.. math::

   v_u &= v \sin\beta , \\
   v_d &= v \cos\beta .

To fix :math:`v_u` and :math:`v_d` in the MSSM in this way,
``MatchingScaleInput`` can be used::

    MatchingScaleInput = {
        {vu, VEV Sin[ArcTan[TanBeta]]},
        {vd, VEV Cos[ArcTan[TanBeta]]}
    };

where ``TanBeta`` is an input parameter.  The symbol ``VEV`` is a
FlexibleSUSY constant which is assigned the value

.. math::

   \text{VEV} = \frac{2 m_Z}{\sqrt{g_Y^2 + g_2^2}} ,

where :math:`m_Z` is the running Z boson mass in the full model,
detetermined by requiring the equality of the Z boson pole masses of
the full model and the Standard Model.  :math:`g_Y` and :math:`g_2`
are the running gauge couplings of :math:`U(1)_Y` and :math:`SU(2)_L`
in the full model, respectively.  These two gauge couplings are
calculated using the 1-loop threshold correction for
:math:`\alpha_{\text{em}}` and the running weak mixing angle,
:math:`\cos\theta_W = m_W / m_Z`.  :math:`m_W` is the running W boson
mass in the full model, detetermined by requiring the equality of the
W boson pole masses of the full model and the Standard Model.

Boundary condition format
`````````````````````````

The variables ``LowScaleInput``, ``SUSYScaleInput`` and ``HighScaleInput``
are lists which specify the boundary conditions for the running model
parameters at the corresponding scale.  The boundary conditions can be
expressed as follows.

Setting a running model parameter to a value or expression
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

A running model parameter can be assigned at the corresponding scale
to a fixed numerical value or a value which is the result of the
evaluation of an expression.  Such an assignment is made by a
two-component list, ``{p, v}``, where the first list element must be the
model parameter (``p`` in this case) and the second list element is a
numerical value or an expression.

Example: An example is the mSUGRA boundary condition in the CMSSM at
the GUT scale::

    HighScaleInput = {
       {T[Ye], Azero*Ye},
       {T[Yd], Azero*Yd},
       {T[Yu], Azero*Yu},
       {mHd2, m0^2},
       {mHu2, m0^2},
       {mq2, UNITMATRIX[3] m0^2},
       {ml2, UNITMATRIX[3] m0^2},
       {md2, UNITMATRIX[3] m0^2},
       {mu2, UNITMATRIX[3] m0^2},
       {me2, UNITMATRIX[3] m0^2},
       {MassB, m12},
       {MassWB,m12},
       {MassG, m12}
    };

The model parameters in the expression in the second list element are
running parameters at the corresponding scale.  I.e. the setting
``{T[Ye], Azero*Ye}`` means :math:`T_{y_e}(Q) := A_0 y_e(Q)`, where
:math:`Q` is the scale.

For matrix- or vector-valued assignments, the following convenient
symbols can be used in the second list element::

    UNITMATRIX[rows]              (* quadratic unit matrix with ``rows' rows *)
    UNITMATRIXCOMPLEX[rows]       (* complex quadratic unit matrix with ``rows' rows *)
    ZEROMATRIX[rows,cols]         (* zero matrix with ``rows' rows and ``cols' columns *)
    ZEROMATRIXCOMPLEX[rows,cols]  (* complex zero matrix with ``rows' rows and ``cols' columns *)
    ZEROVECTOR[rows]              (* zero vector with ``rows' rows *)
    ZEROVECTORCOMPLEX[rows]       (* complex zero vector with ``rows' rows *)

On the r.h.s. of the assignment it is possible to refer to a model
parameter, which is read from an SLHA input block.  These model
parameter input blocks are named after the model parameter output
blocks concatenated with an additionan "IN" (see the SLHA-2 standard,
arXiv:0801.0045, Section 4.1.3).  To refer to such an input model
parameter on the r.h.s. of an assignment one can either add an entry
in ``FSAuxiliaryParameterInfo`` or use the ``LHInput[p]`` symbol,
where ``p`` is the name of the model parameter.

Example::

    SUSYScaleFirstGuess = Sqrt[Sqrt[LHInput[mq2[3,3]] * LHInput[mu2[3,3]]]];
    
    SUSYScaleInput = {
       {mq2, 2 g2^2 LHInput[mq2]}
    };

It is also possible to access the :math:`\beta` functions on the
r.h.s. of an assignment using the ``BETA`` head: ``BETA[p]``
represents the :math:`\beta` function of the parameter ``p`` using the
loop level given in the SLHA input.  ``BETA[l,p]`` represents the
``l``-loop :math:`\beta` function of the parameter ``p``.

Example::

    HighScaleInput = {
        {\[Lambda], BETA[g1] + BETA[g2] + BETA[1,Yu][3,3]}
    };

Temporary parameter re-definitions
''''''''''''''''''''''''''''''''''

Since FlexibleSUSY 1.4.0, the user can perform a temporary parameter
definition to be used in the boundary conditions using the
``FSTemporary[]`` head.

If a parameter ``p`` set in a boundary conditions in the form
``FSTemporary[p,<expr>]``, the following happens: Immediately after the RG
running the value of the parameter is saved locally.  Afterwards, the
parameter is assigned to ``<expr>``.  Now, all further boundary
conditions are imposed and calculations are performed (calculation of
running masses, solution of the EWSB conditions, etc.).  Finally, the
parameter ``p`` is restored to the locally saved value.

Example in ``U1xMSSM3G``: Temporarily rotate the gauge couplings to the
triangular basis::

    g1T  = (g1*gX - g1X*gX1)/Sqrt[gX^2 + gX1^2];
    gXT  = Sqrt[gX^2 + gX1^2];
    g1XT = (g1X*gX + g1*gX1)/Sqrt[gX^2 + gX1^2];
    
    SUSYScaleInput = {
        {FSTemporary[g1], g1T},
        {FSTemporary[gX], gXT},
        {FSTemporary[g1X], g1XT},
        {FSTemporary[gX1], 0},
        {xS, vSInput},
        {x2, Sqrt[4*MZpInput^2 - gX^2*(vu^2 + vd^2)]/(2*gX*Sqrt[1 + TanBetaX^2])},
        {x1, (TanBetaX*Sqrt[4*MZpInput^2 - gX^2*(vu^2 + vd^2)])/(2*gX*Sqrt[1 + TanBetaX^2])},
        {L[lw], 0},
        FSSolveEWSBFor[{mHd2, mHu2, mC12, lw, mS2}]
    };

In this example the gauge couplings, defined in the triangular basis,
are used in every calculation performed at the SUSY scale.  This
includes the calculation of ``x1`` and ``x2`` as well as solving the EWSB
conditions.

Imposing the electroweak symmetry breaking conditions
'''''''''''''''''''''''''''''''''''''''''''''''''''''

The scale, at which the electroweak symmetry breaking (EWSB)
conditions are imposed can be specified by adding
``FSSolveEWSBFor[parameters]`` to the corresponding boundary condition.
The argument ``parameters`` must be the list of model parameters which
are fixed by the electroweak symmetry breaking conditions.

Example: Impose the EWSB conditions at the low-energy scale::

    LowScaleInput = {
       FSSolveEWSBFor[EWSBOutputParameters]
    };

If ``FSSolveEWSBFor[EWSBOutputParameters]`` is not given in any boundary
condition, then it is added to ``SUSYScaleInput``.  This implies, that
by default, the EWSB conditions are imposed at the scale ``SUSYScale``.

Automatic input of unspecified model parameters
'''''''''''''''''''''''''''''''''''''''''''''''

In low-energy models (models where ``OnlyLowEnergyFlexibleSUSY ===
True``) parameters, which are _not_ set in any boundary condition are
automatically input at the ``SUSYScale``.  The values of these
parameters are automatically read from the corresponding SLHA input
blocks.

To disable the automatic input of unspecified parameters, set::

    AutomaticInputAtMSUSY = False;   (* default: True *)

Sub-iterations
''''''''''''''

It is possible to fix model parameters at a scale by performing an
iteration.  Two kinds of iterations are supported:

Minimization
""""""""""""

Model parameters can be fixed by requiring that a function is minimal.
The parameters to be fixed and the function to be minimized must be
specified by the symbol ``FSMinimize[parameters, f]``, where
``parameters`` is the list of parameters to be fixed and ``f`` is the
scalar function to be minimized.

Example::

    SUSYScaleInput = {
       FSMinimize[{vd,vu}, (LowEnergyConstant[MZ] - Pole[M[VZ]])^2 / STANDARDDEVIATION[MZ]^2
                         + (LowEnergyConstant[MH] - Pole[M[hh[1]]])^2 / STANDARDDEVIATION[MH]^2]
    };

Root finding
""""""""""""

Model parameters can be fixed by requiring that a function is zero.
The parameters to be fixed and the function whose zero should be found
must be specified by the symbol ``FSFindRoot[parameters, f]``, where
``parameters`` is the list of parameters to be fixed and ``f`` is the
vector-valued function to be zero.

Example::

    SUSYScaleInput = {
       FSFindRoot[{vd,vu}, {LowEnergyConstant[MZ] - Pole[M[VZ]], LowEnergyConstant[MH] - Pole[M[hh[1]]]}]
    };

Convergence tester
------------------

FlexibleSUSY solves the given boundary value problem (BVP) by running
to each scale and imposing the corresponding boundary conditions until
a convergent solution has been found.

The convergence criterion can be customized using the
``FSConvergenceCheck`` variable.  The default is::

    FSConvergenceCheck = Automatic; (* default *)

If ``FSConvergenceCheck`` is set to ``Automatic``, then the following
convergence criteria are used:

- In SUSY models the BVP solver stops if the maximum number of
  iterations has been reached (``FlexibleSUSY[1]``, see `SLHA input
  file`_ or the maximum relative difference of the DR-bar masses of
  the SUSY particles at the SUSY scale between two successive
  iterations is less than the precision goal (``FlexibleSUSY[0]``, see
  `SLHA input file`_).

- In non-SUSY models the BVP solver stops if the maximum number of
  iterations has been reached (``FlexibleSUSY[1]``, see `SLHA input
  file`_ or the maximum relative difference of all MS-bar masses of
  the model at the SUSY scale between two successive iterations is
  less than the precision goal (``FlexibleSUSY[0]``, see `SLHA input
  file`_).

To create a custom convergence tester, the ``FSConvergenceCheck``
variable must be set to a list containing the running masses and/or
running parameters to be compared between two successive iterations.
The BVP solver stops if the maximum number of iterations has been
reached (``FlexibleSUSY[1]``) or the maximum relative difference of all
running masses and/or parameters given in the ``FSConvergenceCheck``
list at the SUSY scale between two successive iterations is less than
the precision goal (``FlexibleSUSY[0]``).

Example: In the following MSSM example the running masses of all
massive particles as well as the running parameters ``g1, g2, g3, Yu,
Yd[3,3], Ye, B[\[Mu]], \[Mu]`` are tested for convergence.
::

    FSConvergenceCheck = {
        M[hh], M[Ah], M[Hpm],
        M[Su], M[Sd], M[Se],
        M[Chi], M[Cha], M[Glu],
        M[Fu], M[Fd], M[Fe],
        M[VZ], M[VWm],
        g1, g2, g3, Yu, Yd[3,3], Ye, B[\[Mu]], \[Mu]
    };

.. note:: For matrix- or vector-valued parameters every component is
          used in the convergence test, if the matrix/vector indices
          are omitted.


Renormalization group equations (RGEs)
--------------------------------------

The loop order of the RGEs to be used can be selected in the model
file using the ``FSRGELoopOrder`` variable: By setting ``FSRGELoopOrder =
0;`` no RGEs will be generated by SARAH.  By setting ``FSRGELoopOrder =
1;`` only one-loop RGEs will be generated by SARAH.  By setting
``FSRGELoopOrder = 2;`` the two-loop RGEs will be generated by SARAH
(this is the default).

Example::

    FSRGELoopOrder = 2; (* generate two-loop RGEs using SARAH *)

Pole masses
-----------

In order to tune the spectrum generator for speed, the precision of
the pole mass calculation can be selected for each particle.  There
are three different pole mass calculation algorithms available:
``LowPrecision``, ``MediumPrecision`` and ``HighPrecision``.  Please
refer to Section 6.5 of Ref. [1406.2319]_ for details.

By default, the pole masses of all particles are calculated with
``MediumPrecision``, except for the CP-even, CP-odd and charged Higgs
bosons, which are calculated with ``HighPrecision`` in order to include
some momentum-dependent 2-loop corrections.

Example::

    DefaultPoleMassPrecision = MediumPrecision;
    HighPoleMassPrecision    = {hh, Ah, Hpm};
    MediumPoleMassPrecision  = {};
    LowPoleMassPrecision     = {};

Lightest supersymmetric particle (LSP)
--------------------------------------

FlexibleSUSY can generate the helper function ``get_lsp()``, which
returns the mass of the lightest supersymmetric particle (LSP) as well
as the particle type.  The particles which are candidates for being an
LSP must be specified in the ``PotentialLSPParticles`` variable.

Example: In the MSSM the lightest supersymmetric particles might be::

    PotentialLSPParticles = { Chi, Sv, Su, Sd, Se, Cha, Glu };

User-defined replacement rules
------------------------------

User-defined replacement rules can be applied to the beta functions,
self-energies/ tadpoles and vertices.  The rules are specified by the
``FSBetaFunctionRules``, ``FSSelfEnergyRules`` and ``FSVertexRules``
variables, respectively.

Example: Set the gauge couplings ``g1`` and ``g2`` to zero in all 1-loop,
2-loop and 3-loop beta functions::

    FSBetaFunctionRules = {
        {g1 -> 0, g2 -> 0}, (* applied to 1L beta functions *)
        {g1 -> 0, g2 -> 0}, (* applied to 2L beta functions *)
        {g1 -> 0, g2 -> 0}  (* applied to 3L beta functions *)
    };

Example: Set the mass of the Z boson and the corresponding ghost field
to zero in the 1-loop self-energies/ tadpoles::

    FSSelfEnergyRules = {
        { (Mass|Mass2)[VZ|gZ] -> 0 } (* applied to 1L self-energies/tadpoles *)
    };

Example: Set the gauge couplings ``g1`` and ``g2`` to zero in all
vertices::

    FSVertexRules = {
        g1 -> 0,
        g2 -> 0
    };

Observables
-----------

FlexibleSUSY can calculate various observables.  To enable the
calculation of a specific observable, the corresponding symbol must be
added to an extra SLHA output block, see `Output blocks`_ .  In the
following the supported observables are listed.

Effective Higgs-Photon-Photon and Higgs-Gluon-Gluon couplings
`````````````````````````````````````````````````````````````

In the context of [1602.05581]_, FlexibleSUSY has been extended to
calculate the effective couplings of CP-even and CP-odd Higgs bosons
to two photons or two gluons up to NNNLO.  The following table lists
the Mathematica symbols to enable the calculation of these effective
couplings.

================================================================ ==========================================================
 Coupling                                                         Symbol                                                
================================================================ ==========================================================
 CP-even Higgs to two photons, :math:`h\rightarrow\gamma\gamma`   ``FlexibleSUSYObservable`CpHiggsPhotonPhoton``
 CP-odd  Higgs to two photons, :math:`A\rightarrow\gamma\gamma`   ``FlexibleSUSYObservable`CpPseudoScalarPhotonPhoton``
 CP-even Higgs to two gluons,  :math:`h\rightarrow gg`            ``FlexibleSUSYObservable`CpHiggsGluonGluon``
 CP-odd  Higgs to two gluons,  :math:`A\rightarrow gg`            ``FlexibleSUSYObservable`CpPseudoScalarGluonGluon``
================================================================ ==========================================================

Example::

    ExtraSLHAOutputBlocks = {
       {EFFHIGGSCOUPLINGS, NoScale,
               {{1, FlexibleSUSYObservable``CpHiggsPhotonPhoton},
                {2, FlexibleSUSYObservable``CpHiggsGluonGluon},
                {3, FlexibleSUSYObservable``CpPseudoScalarPhotonPhoton},
                {4, FlexibleSUSYObservable``CpPseudoScalarGluonGluon} } }
    };

BSM contributions to the anomalous magnetic moment of the muon
``````````````````````````````````````````````````````````````

Since version 2.0, FlexibleSUSY can calculate the BSM contributions to
the anomalous magnetic moment of the muon, :math:`a_\mu^{\text{BSM}}` at
the 1-loop level, including the leading 2-loop QED logarithmic
corrections.  The following table lists the Mathematica symbols to
enable the calculation of :math:`a_\mu^{\text{BSM}}`.

=================================== =============================================
 Observable                          Symbol
=================================== =============================================
 :math:`a_\mu^{\text{BSM}}`          ``FlexibleSUSYObservable`aMuon``
 :math:`\Delta a_\mu^{\text{BSM}}`   ``FlexibleSUSYObservable`aMuonUncertainty``
=================================== =============================================

:math:`\Delta a_\mu^{\text{BSM}}` is obtained by varying the
renormalization scale by a factor 2.  It therefore represents a \a
lower \a bound of the theoretical uncertainty.

Example::

    ExtraSLHAOutputBlocks = {
       {FlexibleSUSYLowEnergy,
               {{0, FlexibleSUSYObservable``aMuon},
                {1, FlexibleSUSYObservable``aMuonUncertainty} } }
    };

MSSM contributions to the anomalous magnetic moment of the muon with GM2Calc
````````````````````````````````````````````````````````````````````````````

FlexibleSUSY contains an interface to GM2Calc_, which can be used to
calculate the MSSM contributions to the anomalous magnetic moment of
the muon, :math:`a_\mu^{\text{MSSM}}`.  GM2Calc calculates
:math:`a_\mu^{\text{MSSM}}` at the 1-loop level, takes all known
2-loop contributions into account and performs a resummation of
:math:`\tan\beta`-enhanced contributions.

.. note:: GM2Calc version 1.*.* is restricted to CP-conserving MSSM
          without sfermion flavour violation.  Thus, the GM2Calc
          interface in FlexibleSUSY can only be used for MSSM models
          with CP and sfermion flavour conservation.

=================================== ====================================================
 Observable                          Symbol
=================================== ====================================================
 :math:`a_\mu^{\text{MSSM}}`         ``FlexibleSUSYObservable`aMuonGM2Calc``
 :math:`\Delta a_\mu^{\text{MSSM}}`  ``FlexibleSUSYObservable`aMuonGM2CalcUncertainty``
=================================== ====================================================

Example::

    ExtraSLHAOutputBlocks = {
       {FlexibleSUSYLowEnergy,
               {{2, FlexibleSUSYObservable``aMuonGM2Calc},
                {3, FlexibleSUSYObservable``aMuonGM2CalcUncertainty} } }
    };

BSM contributions to the electric dipole moment (EDM)
`````````````````````````````````````````````````````

Since version 2.0 FlexibleSUSY can calculate the BSM contributions to
the electric dipole moments (EDM) of fermions at the 1-loop level in
models with complex parameters.  The following table lists the
Mathematica symbols to enable the calculation of the EDM
:math:`d_f^{\text{BSM}}` of the fermion :math:`f`.

=================================== ====================================================
 Observable                          Symbol
=================================== ====================================================
 :math:`d_f^{\text{BSM}}`            ``FlexibleSUSYObservable`EDM[f]``
=================================== ====================================================

Example: To calculate the EDMs of the electron, muon and tau lepton in
the CP-violating MSSM, add the following to the FlexibleSUSY model file::

    ExtraSLHAOutputBlocks = {
       {FlexibleSUSYLowEnergy,
               {{23, FlexibleSUSYObservable``EDM[Fe[1]]},
                {24, FlexibleSUSYObservable``EDM[Fe[2]]},
                {25, FlexibleSUSYObservable``EDM[Fe[3]]} } }
    };


Output blocks
-------------

The user can define additional SLHA output blocks.  These additional
blocks must be defined in the FlexibleSUSY model file using the
``ExtraSLHAOutputBlocks`` variable.  The ``ExtraSLHAOutputBlocks``
variable is a nested list of the following form::

    ExtraSLHAOutputBlocks = {
       {<blockname>, [<scale>,]
          {{<index>, <expression>},
           {<index>, <expression>},
           {<index>, <expression>}}
       },
       ...
    };

``<blockname>`` is the name of the output block.

Optionally, the renormalization scale can be added to the block head.
``NoScale`` (default) specifies that the block head should have no
scale.  ``CurrentScale`` specifies that the scale written in the block
head should be the current scale of the model.  Otherwise, ``<scale>``
can be numeric value.

The fields inside the block are specified in form of a list of
2-component lists, where the first entry is an integer number
representing the field index.  The second entry is an expression to be
evaluated and whose numeric result is written to the field value.

Example: In the MSSM mode file the following additional output blocks
are defined::

    ExtraSLHAOutputBlocks = {
       {FlexibleSUSYOutput, NoScale,
               {{0, Hold[HighScale]},
                {1, Hold[SUSYScale]},
                {2, Hold[LowScale]} } },
       {FlexibleSUSYLowEnergy,
               {{21, FlexibleSUSYObservable``aMuon} } },
       {EFFHIGGSCOUPLINGS, NoScale,
               {{1, FlexibleSUSYObservable``CpHiggsPhotonPhoton},
                {2, FlexibleSUSYObservable``CpHiggsGluonGluon},
                {3, FlexibleSUSYObservable``CpPseudoScalarPhotonPhoton},
                {4, FlexibleSUSYObservable``CpPseudoScalarGluonGluon} } },
       {ALPHA, NoScale,
               {{ArcSin[Pole[ZH[2,2]]]}}},
       {HMIX , {{1, \[Mu]},
                {2, vu / vd},
                {3, Sqrt[vu^2 + vd^2]},
                {4, M[Ah[2]]^2},
                {101, B[\[Mu]]},
                {102, vd},
                {103, vu} } },
       {Au,    {{1, 1, T[Yu][1,1] / Yu[1,1]},
                {2, 2, T[Yu][2,2] / Yu[2,2]},
                {3, 3, T[Yu][3,3] / Yu[3,3]} } },
       {Ad,    {{1, 1, T[Yd][1,1] / Yd[1,1]},
                {2, 2, T[Yd][2,2] / Yd[2,2]},
                {3, 3, T[Yd][3,3] / Yd[3,3]} } },
       {Ae,    {{1, 1, T[Ye][1,1] / Ye[1,1]},
                {2, 2, T[Ye][2,2] / Ye[2,2]},
                {3, 3, T[Ye][3,3] / Ye[3,3]} } },
       {MSOFT, {{1, MassB},
                {2, MassWB},
                {3, MassG},
                {21, mHd2},
                {22, mHu2},
                {31, SignedAbsSqrt[ml2[1,1]]},
                {32, SignedAbsSqrt[ml2[2,2]]},
                {33, SignedAbsSqrt[ml2[3,3]]},
                {34, SignedAbsSqrt[me2[1,1]]},
                {35, SignedAbsSqrt[me2[2,2]]},
                {36, SignedAbsSqrt[me2[3,3]]},
                {41, SignedAbsSqrt[mq2[1,1]]},
                {42, SignedAbsSqrt[mq2[2,2]]},
                {43, SignedAbsSqrt[mq2[3,3]]},
                {44, SignedAbsSqrt[mu2[1,1]]},
                {45, SignedAbsSqrt[mu2[2,2]]},
                {46, SignedAbsSqrt[mu2[3,3]]},
                {47, SignedAbsSqrt[md2[1,1]]},
                {48, SignedAbsSqrt[md2[2,2]]},
                {49, SignedAbsSqrt[md2[3,3]]} } }
    };

Model-specific switches
-----------------------

Two- and three-loop corrections to pole masses
``````````````````````````````````````````````

MSSM
''''

In the MSSM the dominant two-loop Higgs pole mass corrections
[arxiv:hep-ph/0105096, arxiv:hep-ph/0112177, arxiv:hep-ph/0212132,
arxiv:hep-ph/0206101, arxiv:hep-ph/0305127] can be used by setting in
the model file
::

    UseHiggs2LoopMSSM = True; (* use 2-loop Higgs corrections *)

The known 3-loop Higgs pole mass corrections of the order
:math:`O(\alpha_t\alpha_s^2 + \alpha_b\alpha_s^2)`
[arxiv:hep-ph/0803.0672, arxiv:hep-ph/1005.5709, arxiv:1409.2297,
arxiv:1708.05720] can be used by setting in the model file
::

    UseHiggs3LoopMSSM = True; (* use 3-loop Higgs corrections *)

.. note:: The Himalaya_ library must be linked to FlexibleSUSY in
          order to enable the 3-loop contributions::

              ./configure \
                 --with-models=MSSMNoFVatMGUTHimalaya \
                 --enable-himalaya \
                 --with-himalaya-incdir=${HIMALAYA_DIR}/source/include \
                 --with-himalaya-libdir=${HIMALAYA_DIR}/build

``MSSMNoFVatMGUTHimalaya`` is a pre-defined FlexibleSUSY model which
includes the 3-loop contributions to the light CP-even Higgs mass from
Himalaya.  ``${HIMALAYA_DIR}`` is the path to the Himalaya directory.

To make use of the 2-loop and/or 3-loop corrections the effective
:math:`\mu` parameter must be specified using the ``EffectiveMu``
variable::

    EffectiveMu = \[Mu];

.. note:: When the 3-loop corrections are used, the following switches
          will be set automatically for consistency::

                SARAH`UseHiggs2LoopMSSM = True;
                UseMSSMYukawa2Loop = True; (* use 2-loop SQCD corrections to yt and yb *)
                UseMSSMAlphaS2Loop = True; (* use 2-loop SQCD corrections to alpha_s *)
                UseMSSM3LoopRGEs = True;   (* use 3-loop RGEs *)

NMSSM
'''''

In the NMSSM the dominant two-loop Higgs pole mass corrections from
Ref.  [arXiv:0907.4682] plus the MSSM-like contributions from Refs.
[hep-ph/0105096, hep-ph/0112177, hep-ph/0212132, hep-ph/0206101,
hep-ph/0305127] can be used by setting in the model file::

    UseHiggs2LoopNMSSM = True; (* use 2-loop Higgs corrections *)

In addition, the effective :math:`\mu` parameter must be specified using
the ``EffectiveMu`` variable, Furthermore, the tree-level value of the
effective CP-odd MSSM-like Higgs must be specified in the
``EffectiveMASqr`` variable::

    EffectiveMu = \[Lambda] vS / Sqrt[2];
    EffectiveMASqr = (T[\[Lambda]] vS / Sqrt[2] + 0.5 \[Lambda] \[Kappa] vS^2) (vu^2 + vd^2) / (vu vd);

Standard Model
''''''''''''''

In the Standard Model the two-loop Higgs pole mass corrections of the
order :math:`O(\alpha_t\alpha_s + \alpha_b\alpha_s)` [arxiv:1407.4336],
:math:`O((\alpha_t + \alpha_b)^2)` [arxiv:1205.6497] and
:math:`O(\alpha_\tau^2)` can be used by setting in the model file::

    UseHiggs2LoopSM = True;

The Standard Model the three-loop Higgs pole mass corrections of the
order :math:`O(\alpha_t\alpha_s^2 + \alpha_t^2\alpha_s + \alpha_t^3)`
[arxiv:1407.4336, Eq.(3.2)] can be used by setting in the model file::

    UseHiggs3LoopSM = True;

.. note:: When the 3-loop corrections are used, the following switches
          will be set automatically for consistency::

              UseHiggs2LoopSM = True;
              UseSMYukawa2Loop = True;    (* use 2-loop non-QCD corrections to m_t *)
              UseSMAlphaS3Loop = True;    (* use 2- and 3-loop QCD corrections to alpha_s *)
              UseYukawa3LoopQCD = True;   (* use 2- and 3-loop QCD corrections to m_t *)
              UseSM3LoopRGEs = True;      (* use 3-loop RGEs *)

The Standard Model the 4-loop Higgs pole mass corrections of the order
:math:`O(\alpha_t\alpha_s^3)` [arxiv:1508.00912, Eq.(5.5)] can be used
by setting in the model file::

    UseHiggs4LoopSM = True;

.. note:: When the 4-loop corrections are used, the following switches
          will be set automatically for consistency::

           UseHiggs2LoopSM = True;
           UseHiggs3LoopSM = True;
           UseSMAlphaS4Loop = True;    (* use 2-, 3- and 4-loop QCD corrections to alpha_s *)
           UseYukawa4LoopQCD = True;   (* use 2-, 3- and 4-loop QCD corrections to m_t *)
           UseSM3LoopRGEs = True;      (* use 3-loop RGEs *)
           UseSM4LoopRGEs = True;      (* use 4-loop RGEs *)
           UseSM5LoopRGEs = True;      (* use 4-loop RGEs *)

Split-MSSM
''''''''''

In the split-MSSM (``SplitMSSM``) the two-loop Higgs pole mass
corrections from [arxiv:1312.5220, Eq. (4.8)] of the order
:math:`O(\alpha_t \alpha_s^2)` can be used by setting in the model
file::

    UseHiggs3LoopSplit = True;

Three-loop RGEs for specific models
```````````````````````````````````

Standard Model
''''''''''''''

In the Standard Model the known three-loop RGEs from [arxiv:1303.4364,
arXiv:1307.3536] can be used by setting in the model file::

    UseSM3LoopRGEs = True; (* use three-loop SM RGEs *)

MSSM
''''

In the MSSM the known three-loop RGEs from [hep-ph:0308231]_,
[http://www.liv.ac.uk/~dij/betas/allgennb.log] can be used by setting
in the model file::

    UseMSSM3LoopRGEs = True; (* use three-loop MSSM RGEs *)

Four-loop RGEs for specific models
``````````````````````````````````

Standard Model
''''''''''''''

In the Standard Model the known four-loop RGEs from [arxiv:1508.00912,
arXiv:1604.00853, arxiv:1508.02680] can be used by setting in the
model file::

    UseSM4LoopRGEs = True; (* use four-loop SM RGEs *)

Five-loop RGEs for specific models
``````````````````````````````````

Standard Model
''''''''''''''

In the Standard Model the known five-loop QCD RGE from
[arxiv:1606.08659] can be used by setting in the model file::

    UseSM5LoopRGEs = True; (* use five-loop SM QCD RGE *)

Two-loop threshold corrections
``````````````````````````````

Standard Model
''''''''''''''

The known SM 2-, 3- and 4-loop QCD threshold corrections of order
:math:`O(\alpha_s^2 + \alpha_s^3 + \alpha_s^4)` to the strong coupling
constant are known by [hep-ph/0004189, hep-ph/0512060].  They can be
taken into account by setting in the model file::

    UseSMAlphaS4Loop = True; (* use 2-, 3- and 4-loop threshold for αs *)

The known SM 2-loop threshold corrections of order :math:`O(\alpha_t
\alpha_s + \alpha_t^2)` to the running top mass are known by
[arXiv:1604.01134].  They can be taken into account by setting in the
model file::

    UseSMYukawa2Loop = True; (* use 2-loop thresholds for mt *)

.. note:: These corrections require FlexibleSUSY to be configured with
          TSIL_ support.


MSSM
''''

In the MSSM the known two-loop SQCD relation between the top pole mass
and the DR-bar top mass from
[arxiv:hep-ph/0210258,arxiv:hep-ph/0507139] as well as between the
MS-bar bottom mass in the Standard Model and the DR-bar bottom mass in
the MSSM [arxiv:0707.0650] can be used by setting in the model file::

    UseMSSMYukawa2Loop = True; (* use two-loop threshold for yt and yb *)

The known MSSM two-loop corrections of order :math:`O(\alpha_s^2 +
\alpha_s\alpha_t + \alpha_s\alpha_b)` to the strong coupling
constant are known by [hep-ph/0509048, arXiv:0810.5101,
arXiv:1009.5455]. They can be taken into account by setting in the
model file::

    UseMSSMAlphaS2Loop = True; (* use two-loop threshold for alpha_s *)

Three-loop threshold corrections
````````````````````````````````

Standard Model
''''''''''''''

In non-SUSY models the known 3-loop (Standard Model) QCD corrections
:math:`O(\alpha_s^3)` [arxiv:hep-ph/9911434, arxiv:hep-ph/9912391] can
be used in the determination of the running :math:`\overline{MS}` top
Yukawa coupling :math:`y_t` at the low-energy scale by setting::

    UseYukawa3LoopQCD = Automatic;

or::

    UseYukawa3LoopQCD = True;

Note, that these 3-loop corrections are only applied at run-time if
the threshold correction loop order (block ``FlexibleSUSY[7]``) is set
to a value > 2.

In addition, the 3-loop (Standard Model) QCD corrections
:math:`O(\alpha_s^3)` [arxiv:hep-ph/0004189] to the running
:math:`\overline{MS}` strong coupling :math:`\alpha_s` can be used at
the low-energy scale by setting::

    UseSMAlphaS3Loop = True;

Note, that these 3-loop corrections are only applied at run-time if
the threshold correction loop order (block ``FlexibleSUSY[7]``) is set
to a value > 2.

Four-loop threshold corrections
```````````````````````````````

Standard Model
''''''''''''''

In non-SUSY models the known 4-loop (Standard Model) QCD corrections
:math:`O(\alpha_s^4)` [1604.01134]_, [1502.01030]_, [1606.06754]_ can
be used in the determination of the running :math:`\overline{MS}` top
Yukawa coupling :math:`y_t` at the low-energy scale by setting::

    UseYukawa4LoopQCD = Automatic;

or::

    UseYukawa4LoopQCD = True;


References
----------

.. _GM2Calc: https://arxiv.org/abs/1510.08071
.. _Himalaya: https://github.com/Himalaya-Library/Himalaya
.. _TSIL: https://www.niu.edu/spmartin/tsil/

.. _`SLHA input file`: slha_input.rst

.. [hep-ph:0308231] `Phys.Lett. B579 (2004) 180-188 <https://inspirehep.net/record/626390>`_ [`arxiv:hep-ph/0308231 <https://arxiv.org/abs/hep-ph/0308231>`_]
.. [1406.2319] `CPC 190 (2015) 139-172 <https://inspirehep.net/record/1299998>`_ [`arxiv:1406.2319 <https://arxiv.org/abs/1406.2319>`_]
.. [1602.05581] `Eur.Phys.J. C76 (2016) no.9, 516 <https://inspirehep.net/record/1422208>`_ [`arxiv:1602.05581 <https://arxiv.org/abs/1602.05581>`_]
.. [1604.01134] `Phys.Rev. D93 (2016) no.9, 094017 <https://inspirehep.net/record/1442368>`_ [`arXiv:1604.01134 <https://arxiv.org/abs/1604.01134>`_]
.. [1502.01030] `Phys.Rev.Lett. 114 (2015) 14, 142002 <https://inspirehep.net/literature/1342942>`_ [`arXiv:1502.01030 <https://arxiv.org/abs/1502.01030>`_]
.. [1606.06754] `Phys.Rev.D 94 (2016) 7, 074025 <https://inspirehep.net/literature/1471728>`_ [`arXiv:1606.06754 <https://arxiv.org/abs/1606.06754>`_]
