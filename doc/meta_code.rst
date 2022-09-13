.. sectnum::

Meta code
=========

.. contents:: Table of Contents

Re-usable meta code packages
----------------------------

FlexibleSUSY contains several packages with general algorithms needed
in Quantum Field Theory calculations and implementations of analytic
expressions from the literature.


Passarino-Veltman loop functions
````````````````````````````````

``meta/LoopFunctions.m`` contains an implementation of
Passarino-Veltman 1-loop functions :math:`A_0`, :math:`B_0`,
:math:`B_1`, :math:`B_{00}`, :math:`B_{11}`, etc.


Passarino-Veltman loop functions for zero external momenta
``````````````````````````````````````````````````````````

``meta/LoopFunctionsZeroMomentum.m`` contains an implementation of the
Passarino-Veltman 1-loop functions :math:`A_0`, :math:`B_0`,
:math:`\frac{\partial B_0}{\partial p^2}`, :math:`C_0`, :math:`D_0`,
:math:`E_0` and :math:`F_0` for zero external momenta.

Example::

    Get["meta/LoopFunctionsZeroMomentum.m"];
    F0[x,x,x,x,x,x,q] //. loopFunctionsZeroMomentum

yields::

    1/(20*x^4)


Renormalization group equation integrator
`````````````````````````````````````````

``meta/RGIntegrator.m`` contains a routine, which perturbatively
integrates a system of coupled renormalization group equations::

    Get["meta/RGIntegrator.m"];
    ?RGIntegrate

Example 1: Integrating a single RGE up to 2-loop order::

    Get["meta/RGIntegrator.m"];
    betas = { {g, h a g^2} };
    RGIntegrate[betas, Q1, Q2, loopOrder -> 2]

yields::

    {g[Q1] -> g[Q2] + a h g[Q2]^2 Log[Q1/Q2] + a^2 h^2 g[Q2]^3 Log[Q1/Q2]^2}

Example 2: Interation of the following coupled system::

    Get["meta/RGIntegrator.m"];
    betas = { {g, h a g^2, h^2 b g^2 l^2},
              {l, h (c l^2 + d g^2), h^2 (e l^4 + f g^2 l^2)} };
    RGIntegrate[betas, Q1, Q2]

yields::

    {g[Q1] -> g[Q2] + a*h*g[Q2]^2*Log[Q1/Q2] +
        h^2*(b*g[Q2]^2*l[Q2]^2*Log[Q1/Q2] + a^2*g[Q2]^3*Log[Q1/Q2]^2),
      l[Q1] -> l[Q2] + h*(d*g[Q2]^2*Log[Q1/Q2] + c*l[Q2]^2*Log[Q1/Q2]) +
        h^2*(f*g[Q2]^2*l[Q2]^2*Log[Q1/Q2] + e*l[Q2]^4*Log[Q1/Q2] +
          a*d*g[Q2]^3*Log[Q1/Q2]^2 + c*d*g[Q2]^2*l[Q2]*Log[Q1/Q2]^2 +
          c^2*l[Q2]^3*Log[Q1/Q2]^2)}


Beta functions
``````````````

Standard Model
''''''''''''''

``meta/ThreeLoopSM.m`` contains a routine which returns the beta
functions of the Standard Model up to the 5-loop level.  The beta
functions are stored in the files ``meta/SM/beta_*.m``, which have
been obtained at 3-loop level from [1303.4364]_, [1504.05200]_ with
some additional 4- and 5-loop corrections from [1508.00912]_,
[1508.02680]_, [1604.00853]_, [1606.08659]_.

.. note:: The loop factors :math:`1/(4\pi)^2` have been omitted.

Example: Extracting the 5-loop QCD beta function for :math:`g_3`::

    b = Get["meta/SM/beta_g3.m"];

    bg3 = { k^1 b[[1]],
            k^2 b[[2]],
            k^3 b[[3]],
            k^4 b[[4]],
            k^5 b[[5]] };

    bg3 /. { g1 -> 0, g2 -> 0, gb -> 0, gt -> 0 } // Chop

Output::

    { -7*g3^3*k,
      -26*g3^5*k^2,
      (65*g3^7*k^3)/2,
      -2472.2837425797156*g3^9*k^4,
      271.4283824198132*g3^11*k^5 }

MSSM
''''

``meta/ThreeLoopMSSM.m`` contains a routine which returns the beta
functions of the MSSM up to the 3-loop level.  The beta functions are
stored in the files ``meta/MSSM/beta_*.m``, which have been obtained
from http://www.liv.ac.uk/~dij/betas/allgennb.log [hep-ph:0308231]_.

Example: Extracting the 3-loop QCD beta function for :math:`g_3` in
the MSSM::

    trace[args__] := Tr[Dot[args]];
    Adj = ConjugateTranspose;
    Yt = Yb = Ye = Array[0&, {3,3}];

    b = Get["meta/MSSM/beta_g3.m"];

    bg3 = { k^1 b[[1]],
            k^2 b[[2]],
            k^3 b[[3]] };

    bg3 /. { g1 -> 0, g2 -> 0 }

Output::

     {-3*g3^3*k, 14*g3^5*k^2, (347*g3^7*k^3)/3}


Loop corrections to masses
``````````````````````````

Standard Model
''''''''''''''

``meta/SM/Mh2_effpot.m`` contains the QCD contributions to the 4-loop
effective Higgs potential in the Standard Model from [1508.00912]_

Example::

    Get["meta/SM/Mh2_effpot.m"];

    k = 1/(4 Pi)^2;
    yt = 0.9;
    g3 = 1.166;
    v = 247.5;
    \[Lambda] = 0.25;
    mu2 = -8.55 10^3;
    mt = yt v / Sqrt[2];
    Q = 173.34;

    Sqrt[DMh2]

Output::

    124.926

------

``meta/ThreeLoopQCD.m`` contains a routine, which returns the ratio of
the :math:`\overline{\text{MS}}` top mass over the top pole mass in
the SM up to 3-loop level in QCD from [hep-ph:9912391]_, Eq. (10).
The expression contains the full renormalization scale dependence,
which has been taken from [hep-ph:9911434]_.

Example::

    Get["meta/ThreeLoopQCD.m"];
    Start["SM"];
    FlexibleSUSY`M[Fu] = mt;

    h = k (4 Pi)^2;

    MfOvermf = GetMTopPoleOverMTopMSbar[{1,h,h^2,h^3}] /. {
        Log[Q^2/mt^2] -> -Lbar[t],
        Log[mt^2/Q^2] -> Lbar[t]
    };

    Mt = mt N[Collect[MfOvermf, {k, g3, Lbar[__]}, Simplify]]

Output::

    mt*(1 +
        g3^2*k*(5.333333333333333 - 4*Lbar[t]) +
        g3^4*k^2*(131.78498721717762 - 80.66666666666667*Lbar[t] + 22*Lbar[t]^2) +
        g3^6*k^3*(4712.740192659316 - 2031.1382275647934*Lbar[t] + 710*Lbar[t]^2 - 132*Lbar[t]^3))

------


``meta/TwoLoopQCD.m`` contains routines, which return the ratio of the
top pole mass over the running top mass up to the 2-loop level in the
:math:`\overline{\text{MS}}` and :math:`\overline{\text{DR}}` schemes
[hep-ph:0210258]_, [hep-ph:9803493]_.

MSSM
''''

``meta/TwoLoopMSSM.m`` contains routines, which return the
analytic 2-loop corrections to the Higgs masses in the CP-conserving
MSSM [hep-ph:0105096]_.


Threshold corrections
`````````````````````

Standard Model
''''''''''''''

``meta/SM/mf_3loop_qcd.m`` contains the 3-loop relation between a quark pole
mass and the corresponding running :math:`\overline{\text{MS}}` mass
from [hep-ph:9912391]_, [hep-ph:9911434]_.

Example::

    Get["meta/SM/mf_3loop_qcd.m"];

    L = Lbar[t];
    NL = 5; (* number of light quark masses *)
    NH = 1; (* number of heavy quark masses *)

    Mt = mt N[Collect[MfOvermf, {k, g3, Lbar[__]}, Simplify]]

Output::

    mt*(1 +
        g3^2*k*(5.333333333333333 - 4*Lbar[t]) +
        g3^4*k^2*(131.78498721717762 - 80.66666666666667*Lbar[t] + 22*Lbar[t]^2) +
        g3^6*k^3*(4712.740192659316 - 2031.1382275647934*Lbar[t] + 710*Lbar[t]^2 - 132*Lbar[t]^3))

------

``meta/SM/mt_4loop_qcd.m`` contains the 4-loop QCD relation between
the top quark pole mass and the corresponding running
:math:`\overline{\text{MS}}` mass from [1604.01134]_, [1502.01030]_,
[1606.06754]_.

Example::

    Get["meta/SM/mt_4loop_qcd.m"];

    L = Lbar[t];

    Mt = mt N[Collect[MtOvermt, {k, g3, Lbar[__]}, Simplify]]

Output::

    mt*(1 +
        g3^2*k*(5.333333333333333 - 4*Lbar[t]) +
        g3^4*k^2*(131.78498721717762 - 80.66666666666667*Lbar[t] + 22*Lbar[t]^2) +
        g3^6*k^3*(4712.740192659316 - 2031.1382275647934*Lbar[t] + 710*Lbar[t]^2 - 132*Lbar[t]^3) +
        g3^8*k^4*(211681.74421123447 - 104673.38261571848*Lbar[t] + 22162.91142653778*Lbar[t]^2 - 5638*Lbar[t]^3 + 825*Lbar[t]^4))

------

``meta/SM/mt_2loop_gaugeless.m`` contains the 2-loop relation between
the top quark pole mass and the corresponding running
:math:`\overline{\text{MS}}` mass from [1604.01134]_ in the gaugeless
limit.  Contributions of :math:`O(\alpha_s^2, \alpha_s\alpha_t,
\alpha_t^2, \alpha_t\lambda^n)` are included.  The relation is
gauge-independent.

------

``meta/SM/as_4loop_qcd.m`` contains the 4-loop QCD relation between
the strong coupling :math:`\alpha_s^{n_f}` with :math:`n_f` flavours
and :math:`\alpha_s^{n_l}` with :math:`n_l = n_f - 1` flavours in the
:math:`\overline{\text{MS}}` scheme [hep-ph:0512060]_.

Example::

    nl = 5;

    L = Log[Q^2/mf[Q]^2];

    alphaS = Get["meta/SM/as_4loop_qcd.m"];


2HDM
''''

``meta/THDM/Thresholds_1L_full.m`` contains the implementation of the
complete analytic 1-loop threshold corrections of the THDM and the
THDM + Higgsinos + gauginos to the MSSM [0901.2065]_.

Example::

    Get["meta/THDM/Thresholds_1L_full.m"];

    tc = (4 Pi)^2 GetTHDMThresholds1L[];

    $Assumptions = { Element[ht, Reals], Element[Mu, Reals], Element[At, Reals] };

    Yu[i_, k_] := DiagonalMatrix[{0,0,ht}][[i,k]];
    Tu[i_, k_] := DiagonalMatrix[{0,0,ht At}][[i,k]];
    Yd[__] := 0;
    Ye[__] := 0;
    Td[__] := 0;
    Te[__] := 0;
    g2 = gY = 0;

    {l1, l2, l3, l4, l5, l6, l7} = Collect[tc, ht, Simplify]

Output::

    { -3*ht^4*Mu^4*D0[msq[3], msq[3], msu[3], msu[3]],
      -3*ht^4*(B0[msq[3], msq[3], Q] + B0[msu[3], msu[3], Q] + At^2*(2*C0[msq[3], msq[3], msu[3]] + 2*C0[msq[3], msu[3], msu[3]] + At^2*D0[msq[3], msq[3], msu[3], msu[3]])),
      -3*ht^4*Mu^2*(C0[msq[3], msu[3], msu[3]] + At^2*D0[msq[3], msq[3], msu[3], msu[3]]),
      -3*ht^4*Mu^2*(C0[msq[3], msq[3], msu[3]] + At^2*D0[msq[3], msq[3], msu[3], msu[3]]),
      -3*At^2*ht^4*Mu^2*D0[msq[3], msq[3], msu[3], msu[3]],
       3*At*ht^4*Mu^3*D0[msq[3], msq[3], msq[3], msq[3]],
       3*At*ht^4*Mu*(C0[msq[3], msq[3], msu[3]] + C0[msq[3], msu[3], msu[3]] + At^2*D0[msq[3], msq[3], msq[3], msq[3]]) }

MSSM
''''

``meta/MSSM/tquark_2loop_strong.m`` contains the analytic
expression for the 2-loop relation :math:`O(\alpha_s^2)` between the top
quark pole mass and the :math:`\overline{\text{DR}}` top mass in the
MSSM [hep-ph:0210258]_, [hep-ph:0507139]_.

----

``meta/MSSM/bquark_2loop_sqcd_decoupling.m`` contains the analytic
expression for the 2-loop relation :math:`O(\alpha_s^2)` between the
:math:`\overline{\text{MS}}` bottom quark mass in the Standard Model
(without the top quark) and the :math:`\overline{\text{DR}}` bottom
mass in the MSSM [0707.0650]_.

----

``meta/MSSM/dmtauas2.m`` contains the analytic expression for the
2-loop relation between the tau lepton pole mass and the
:math:`\overline{\text{DR}}` tau mass in the MSSM.

----

``meta/MSSM/das2.m`` contains the analytic expression for the 2-loop
relation between the :math:`\overline{\text{MS}}` :math:`\alpha_s` in
the Standard Model (without the top quark) and the
:math:`\overline{\text{DR}}` value in the MSSM [hep-ph:0509048]_,
[0810.5101]_, [1009.5455]_.


C++ source code formatting
``````````````````````````

``meta/TextFormatting.m`` contains routines for text formatting of
long expressions in C/C++ form, see ``WrapText[]`` and
``IndentText[]``.

Example: Formatting long expression::

    Get["meta/TextFormatting.m"];

    (* long expression *)
    dmt = Get["meta/MSSM/tquark_2loop_strong.m"];

    maxWidth = 70;
    indent = 3;

    "dmt = " <> WrapText[ToString[dmt, CForm], maxWidth, indent]

Output::

    dmt = (Power(GS,4)*((-11*colorCA*colorCF*MGl*mmst1*s2t)/(-mmgl + mmst1) + (6
       *Power(colorCF,2)*MGl*mmst1*s2t)/(-mmgl + mmst1) - (6*Power(colorCF
       ,2)*MGl*mmst2*s2t)/(-mmgl + mmst1) + (6*Power(colorCF,2)*MGl*mmst1*
       mmst2*s2t)/((-mmgl + mmst1)*(mmst1 - mmst2)) - (6*Power(colorCF,2)*
       MGl*Power(mmst2,2)*s2t)/((-mmgl + mmst1)*(mmst1 - mmst2)) + (11*
       colorCA*colorCF*MGl*mmst2*s2t)/(-mmgl + mmst2) - (6*Power(colorCF,2
       )*MGl*mmst2*s2t)/(-mmgl + mmst2) - (Power(colorCF,2)*MGl*mmst1*
       Power(s2t,3))/(-mmgl + mmst1) + (7*Power(colorCF,2)*MGl*mmst2*Power
       (s2t,3))/(-mmgl + mmst1) - (6*Power(colorCF,2)*MGl*mmst1*mmst2*[...]


References
----------

.. [hep-ph:9803493] `Nucl.Phys. B539 (1999) 671-690 <https://inspirehep.net/record/468752>`_ [`arXiv:hep-ph/9803493 <https://arxiv.org/abs/hep-ph/9803493>`_]
.. [hep-ph:9911434] `Nucl.Phys. B573 (2000) 617-651 <https://inspirehep.net/record/510551>`_ [`arXiv:hep-ph/9911434 <https://arxiv.org/abs/hep-ph/9911434>`_]
.. [hep-ph:9912391] `Phys.Lett. B482 (2000) 99-108 <https://inspirehep.net/record/522686>`_ [`arXiv:hep-ph/9912391 <https://arxiv.org/abs/hep-ph/9912391>`_]
.. [hep-ph:0105096] `Nucl.Phys. B611 (2001) 403-422 <https://inspirehep.net/record/556417>`_ [`arXiv:hep-ph/0105096 <https://arxiv.org/abs/hep-ph/0105096>`_]
.. [hep-ph:0210258] `Eur.Phys.J. C29 (2003) 87-101 <https://inspirehep.net/record/600038>`_ [`arXiv:hep-ph/0210258 <https://arxiv.org/abs/hep-ph/0210258>`_]
.. [hep-ph:0308231] `Phys.Lett. B579 (2004) 180-188 <https://inspirehep.net/record/626390>`_ [`arXiv:hep-ph/0308231 <https://arxiv.org/abs/hep-ph/0308231>`_]
.. [hep-ph:0507139] `Phys.Atom.Nucl. 71 (2008) 343-350 <https://inspirehep.net/record/687205>`_ [`arXiv:hep-ph/0507139 <https://arxiv.org/abs/hep-ph/0507139>`_]
.. [hep-ph:0509048] `Phys.Rev. D72 (2005) 095009 <https://inspirehep.net/record/691479>`_ [`arXiv:hep-ph/0509048 <https://arxiv.org/abs/hep-ph/0509048>`_]
.. [hep-ph:0512060] `Nucl.Phys. B744 (2006) 121-135 <https://inspirehep.net/record/699609>`_ [`arXiv:hep-ph/0512060 <https://arxiv.org/abs/hep-ph/0512060>`_]
.. [0707.0650] `Int.J.Mod.Phys. A22 (2007) 5245-5277 <https://inspirehep.net/record/755029>`_ [`arXiv:0707.0650 <https://arxiv.org/abs/0707.0650>`_]
.. [0810.5101] `JHEP 0902 (2009) 037 <https://inspirehep.net/record/800842>`_ [`arXiv:0810.5101 <https://arxiv.org/abs/0810.5101>`_]
.. [0901.2065] `Phys.Rev. D84 (2011) 034030 <https://inspirehep.net/record/811006>`_ [`arXiv:0901.2065 <https://arxiv.org/abs/0901.2065>`_]
.. [1009.5455] `C10-06-06.1 <https://inspirehep.net/record/871111>`_ [`arXiv:1009.5455 <https://arxiv.org/abs/1009.5455>`_]
.. [1303.4364] `Nucl.Phys. B875 (2013) 552-565 <https://inspirehep.net/record/1224266>`_ [`arXiv:1303.4364 <https://arxiv.org/abs/1303.4364>`_]
.. [1504.05200] `JHEP 1507 (2015) 159 <https://inspirehep.net/record/1362483>`_ [`arXiv:1504.05200 <https://arxiv.org/abs/1504.05200>`_]
.. [1508.00912] `Phys.Rev. D92 (2015) no.5, 054029 <https://inspirehep.net/record/1386688>`_ [`arXiv:1508.00912 <https://arxiv.org/abs/1508.00912>`_]
.. [1508.02680] `Phys.Lett. B762 (2016) 151-156 <https://inspirehep.net/record/1387530>`_ [`arXiv:1508.02680 <https://arxiv.org/abs/1508.02680>`_]
.. [1604.00853] `JHEP 1606 (2016) 175 <https://inspirehep.net/record/1441223>`_ [`arXiv:1604.00853 <https://arxiv.org/abs/1604.00853>`_]
.. [1604.01134] `Phys.Rev. D93 (2016) no.9, 094017 <https://inspirehep.net/record/1442368>`_ [`arXiv:1604.01134 <https://arxiv.org/abs/1604.01134>`_]
.. [1502.01030] `Phys.Rev.Lett. 114 (2015) 14, 142002 <https://inspirehep.net/literature/1342942>`_ [`arXiv:1502.01030 <https://arxiv.org/abs/1502.01030>`_]
.. [1606.06754] `Phys.Rev.D 94 (2016) 7, 074025 <https://inspirehep.net/literature/1471728>`_ [`arXiv:1606.06754 <https://arxiv.org/abs/1606.06754>`_]
.. [1606.08659] `Phys.Rev.Lett. 118 (2017) no.8, 082002 <https://inspirehep.net/record/1472834>`_ [`arXiv:1606.08659 <https://arxiv.org/abs/1606.08659>`_]
