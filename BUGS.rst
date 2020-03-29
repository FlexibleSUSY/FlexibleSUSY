Reporting a bug
===============

If you think you may have found a bug in FlexibleSUSY, which is not
listed below, please file a `bug report`_.


Known problems
==============

Exception due to multi-threading
--------------------------------

If FlexibleSUSY exits with an exception indicating that threads cannot
be opened (which happens for example with g++ 4.8.x), one can try
adding ``-Wl,--no-as-needed`` to the ``CXXFLAGS``, i.e. for example::

    ./configure --with-cxxflags="-Wl,--no-as-needed -std=c++11 -O2"

If this does not solve the problem, one can disable multi-threading
via::

    ./configure --disable-threads


Tensor couplings
----------------

FlexibleSUSY cannot yet deal with tensor couplings (matrix-like
couplings with rank 3, i.e. with 3 or more indices) in the
Lagrangian/Superpotential.  In order to make such models work in
FlexibleSUSY, please split the tensor coupling terms in the
Lagrangian/Superpotential into a sum of terms with a matrix-like
coupling structure of rank 2 or less.


4th generation
--------------

FlexibleSUSY cannot yet include all effects from fermions, which mix
with the Standard Model fermions, in the extraction of the running
Yukawa couplings.  In such models, the Yukawa couplings have to be
fixed by hand in one boundary condition.  FlexibleSUSY provides the
symbols

* ``upQuarksDRbar`` = (mu, ms, mt)
* ``downQuarksDRbar`` = (md, mc,mb)
* ``downLeptonsDRbar`` = (me, mµ, mτ)

to access the running Standard Model fermion masses in such models.
See the munuSSM model file for an example, where the charginos mix
with the Standard Model fermions.  In this model the down-type lepton
Yukawa matrix is fixed in an approximated form at the low-energy scale
given the running SM charged lepton masses.


Compilation errors related to ``TensorStorage``
-----------------------------------------------

The ``Tensor`` module of the Eigen library before version 3.3.5 may
fail to compile with ``clang++`` due to a bug in Eigen.  The bug has
been fixed in Eigen 3.3.5.


Compilation errors related to complex matrix multiplication
-----------------------------------------------------------

When models with complex parameters are built using versions of
Eigen 3.3.0, it may be the case that compilation of the generated
code fails with an error message indicating that different numeric
types have been mixed without using an appropriate cast.  This arises
because of limitations in mixing different scalar types in earlier
versions of Eigen.  In some cases, this problem can be worked around
by reformulating the boundary conditions in the model file.
For example, in a model such as the MSSMCPV, a boundary condition of
the form::

    {T[Ye], (ReAeInput + I ImAeInput) Ye},
    {T[Yd], (ReAdInput + I ImAdInput) Yd},
    {T[Yu], (ReAuInput + I ImAuInput) Yu}

may trigger this compilation error when compiled against versions of
Eigen below 3.3.0, due to the terms in brackets involving addition of
a real and complex matrix.  In cases such as this, the problem can be
avoided by rewriting the boundary condition to avoid an operation
involving mixed matrix types.  For example, in the equivalent form::

    {T[Ye], ReAeInput Ye + I ImAeInput Ye},
    {T[Yd], ReAdInput Yd + I ImAdInput Yd},
    {T[Yu], ReAuInput Yu + I ImAuInput Yu}

all additions involve only complex matrices.

If it is not possible to solve this compilation error by
reformulating the boundary conditions of the model, it is recommended
that the installed Eigen version be upgraded to at least 3.3.0, in
which this problem no longer occurs due to better support for
operations involving mixed numerical types.


Mathematica crashes on exit when multiple LibraryLinks are loaded
-----------------------------------------------------------------

It can happen that Mathematica crashes on exit when multiple
Mathematica interfaces of FlexibleSUSY (LibraryLinks) were loaded.
The crash occurs only when Mathematica is closed, not during run-time.
Depending on the Mathematica version and the used system, the crash
could be caused by one of the following problems::

    *** Error in ``[...]/WolframKernel': corrupted size vs. prev_size
    *** Error in ``[...]/WolframKernel': corrupted double-linked list
    *** Error in ``[...]/WolframKernel': free(): invalid pointer
    *** Error in ``[...]/MathKernel': double free or corruption (!prev)

In some cases the crash could be avoided by disabling multi-threading
(see ``./configure --help``) or by using a different Mathematica
version.

When such a crash occurs during run-time, please file a `bug report`_.


.. _bug report: https://github.com/FlexibleSUSY/FlexibleSUSY/issues
