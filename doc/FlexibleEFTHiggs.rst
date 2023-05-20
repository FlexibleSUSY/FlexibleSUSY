.. sectnum::

FlexibleEFTHiggs
================

.. contents:: Table of Contents

Overview
--------

By default, FlexibleSUSY creates a "full model" spectrum generator,
where the pole mass spectrum is calculated at the 1- or 2-loop level
in the MS-bar/DR-bar scheme in the full model.

FlexibleEFTHiggs [1609.00371]_, [1710.03760]_ is an approach of
FlexibleSUSY, to perform the 1-loop calculation of the lightest Higgs
pole mass of the given BSM model in an effective field theory (EFT),
which is assumed to be the Standard Model.  FlexibleEFTHiggs combines
the features of full model calculation (the inclusion of all
logarithmic and non-logarithmic Higgs mass contributions) with the
ones of an EFT (the resummation of leading and sub-leading logarithms
to all orders).  In FlexibleEFTHiggs, the quartic Higgs coupling of
the Standard Model is determined at the SUSY scale by requiring that
the lightest CP-even Higgs pole mass in the full model is equarl to
the Standard Model Higgs pole mass.

FlexibleEFTHiggs is exact at the 1-loop level.  In particular, all
power-suppressed 1-loop terms of the order :math:`O(v^2/M^2)`, where
:math:`M` is the scale of heavy new non-Standard Model particles, are
correctly taken into account.  In addition, leading and sub-leading
logarithms of the scale :math:`M` are resummed to all orders.  At the
2-loop level, FlexibleEFTHiggs misses only non-logarithmic BSM
contributions.

FlexibleEFTHiggs spectrum generator
-----------------------------------

In order to create a FlexibleEFTHiggs spectrum generator for a given
model, SARAH model files and a FlexibleSUSY model file must be
provided, just as in the case of a "full model" spectrum generator.

In the FlexibleSUSY model file, the ``FlexibleEFTHiggs`` variable must
be set to ``True``::

    FlexibleEFTHiggs = True;

In FlexibleEFTHiggs, the matching of the full model to the Standard
Model is performed at the ``SUSYScale`` (except, if the value of the
matching scale is overwritten by setting the ``FlexibleEFTHiggs[19]``
to a non-zero value).

The low-scale constaint is completely ignored, i.e. the variables
``LowScale``, ``LowScaleFirstGuess``, ``LowScaleInput`` and
``InitialGuessAtLowScale`` have no effect.

The model parameters must be set in the SUSY-scale or high-scale
constraint.  Initial values for the model parameters can be given in
the ``InitialGuessAtSUSYScale`` or ``InitialGuessAtHighScale``
variables.

\note If ``OnlyLowEnergyFlexibleSUSY = True``, then the high-scale
constraint is ignored.  In this case, only the SUSY-scale constraint
is available in FlexibleEFTHiggs.

The Standard Model gauge and Yukawa couplings as well as the SM-like
VEV of the full model are calculated automatically using a full 1-loop
calculation.  Therefore, the :math:`SU(3)_C\times SU(2)_L\times
U(1)_Y` gauge couplings and the up- and down-Quark and lepton Yukawa
couplings don't need to be specified in any of the constraints.

In many models the determination of the running Yukawa couplings
requires the knowledge of the running VEVs.  These VEVs are sometimes
related to the SM-like VEV :math:`v = \sqrt{v_u^2 + v_d^2}`.  For
example in the MSSM the relation reads,

.. math::

   v_u &= v \sin\beta ,
   v_d &= v \cos\beta .

Such a matching condition can be set using the ``MatchingScaleInput``
variable, see the corresponding section in `FlexibleSUSY model file`_.

Example: in the MSSM
::

    MatchingScaleInput = {
        {vu, VEV Sin[ArcTan[TanBeta]]},
        {vd, VEV Cos[ArcTan[TanBeta]]}
    };

The symbol ``VEV`` is a FlexibleSUSY constant which refers to the
running SM-like vacuum expectation value in the full model.  See
`FlexibleSUSY model file`_ for the precise definition.

MSSM example for FlexibleEFTHiggs
---------------------------------

An example for the general MSSM can be found in
``model_files/MSSMEFTHiggs/FlexibleSUSY.m.in``.  Below, we show a
simplified FlexibleEFTHiggs/MSSM spectrum generator, which takes only
three input parameters: The SUSY scale :math:`M_\text{S}`, the stop
mixing parameter :math:`X_t`, and :math:`\tan\beta`::

    FSModelName = "@CLASSNAME@";
    FSEigenstates = SARAH`EWSB;
    FSDefaultSARAHModel = MSSM;
    OnlyLowEnergyFlexibleSUSY = True;
    FlexibleEFTHiggs = True;
    
    MINPAR = {
        {4, Sign[\[Mu]]}
    };
    
    EXTPAR = {
        {0, Ms},      (* SUSY scale *)
        {14, Xtt},    (* Xt / Ms *)
        {25, TanBeta}
    };
    
    EWSBOutputParameters = { mHd2, mHu2 };
    
    SUSYScale = Ms;
    
    SUSYScaleFirstGuess = Ms;
    
    SUSYScaleInput = {
        {MassB, Ms},
        {MassWB, Ms},
        {MassG, Ms},
        {mq2, UNITMATRIX[3] Ms^2},
        {mu2, UNITMATRIX[3] Ms^2},
        {md2, UNITMATRIX[3] Ms^2},
        {ml2, UNITMATRIX[3] Ms^2},
        {me2, UNITMATRIX[3] Ms^2},
        {\[Mu], Ms},
        {B[\[Mu]], Ms^2/(TanBeta + 1/TanBeta)},
        {T[Yu], Ms/TanBeta Yu},
        {T[Yd], Ms TanBeta Yd},
        {T[Ye], Ms TanBeta Ye},
        {T[Yu][3,3], (Ms/TanBeta + Xtt Ms) Yu[3,3]}
    };
    
    InitialGuessAtSUSYScale = SUSYScaleInput;
    
    MatchingScaleInput = {
        {vu, VEV Sin[ArcTan[TanBeta]]},
        {vd, VEV Cos[ArcTan[TanBeta]]}
    };
    
    UseHiggs2LoopMSSM = True;
    EffectiveMu = \[Mu];


References
----------

.. _`FlexibleSUSY model file`: model_file.rst

.. [1609.00371] `JHEP 1701 (2017) 079 <https://inspirehep.net/record/1484857>`_ [`arXiv:1609.00371 <https://arxiv.org/abs/1609.00371>`_]
.. [1710.03760] `CPC 230 (2018) 145-217 <https://inspirehep.net/record/1629978>`_ [`arXiv:1710.03760 <https://arxiv.org/abs/1710.03760>`_]
