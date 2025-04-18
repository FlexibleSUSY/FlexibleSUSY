.. sectnum::

============
FlexibleSUSY
============

|release| |commits-since| |static-analysis| |tests|

.. |release| image:: https://img.shields.io/github/v/release/FlexibleSUSY/FlexibleSUSY
.. |commits-since| image:: https://img.shields.io/github/commits-since/FlexibleSUSY/FlexibleSUSY/latest/development
   :alt: GitHub commits since latest release (by SemVer including pre-releases)
.. |static-analysis| image:: https://github.com/FlexibleSUSY/FlexibleSUSY/workflows/static%20analysis/badge.svg?branch=development
   :target: https://github.com/FlexibleSUSY/FlexibleSUSY/actions/workflows/static-analysis.yml
.. |tests| image:: https://github.com/FlexibleSUSY/FlexibleSUSY/actions/workflows/tests.yml/badge.svg?branch=development
   :target: https://github.com/FlexibleSUSY/FlexibleSUSY/actions/workflows/tests.yml

.. image:: doc/images/FS-logo.png
   :align: right

FlexibleSUSY provides Mathematica and C++ code to create spectrum
generators for supersymmetric and non-supersymmetric models.  It is
designed for generating fast and modular C++ code, allowing for easy
modification, extension and reuse.

* Homepage:                https://flexiblesusy.hepforge.org
* Mailing list:            flexiblesusy@projects.hepforge.org
* Source code repository:  https://github.com/FlexibleSUSY
* Bug reports:             https://github.com/FlexibleSUSY/FlexibleSUSY/issues
* References: [1406.2319]_, [1609.00371]_, [1710.03760]_, [2204.05285]_

  If you use **FlexibleSUSY** in your work please cite [1406.2319]_
  and [1710.03760]_.

  If you use the **FlexibleEFTHiggs** approach in your work please
  cite [1609.00371]_ and [1710.03760]_.

  If you use **FlexibleSUSY+Himalaya** or Himalaya_ in your work,
  please cite [1005.5709]_, [1708.05720]_, [1807.03509]_ and
  [1910.03595]_.

  If you use **FlexibleSUSY+GM2Calc** in your work, please cite
  [1510.08071]_ and (if appropriate) [2110.13238]_.

  If you use **FlexibleDecay** in your work, please cite [2106.05038]_ and
  [1211.2099]_. If you compute Higgs-physics constrainst using HiggsTools or
  Lilith please also cite [2502.XXXXX]_ and (as appropriate) [2210.09332]_ or
  [1502.04138]_, [1908.03952]_, [2012.11408]_.

  If you use the W boson pole mass prediction in FlexibleSUSY 2.7.0
  (or later), please cite [2204.05285]_.

  If you use **NPointFunctions** in your work, please cite [2402.14630]_.

  If you use unitarity constraints please cite [XXXX.XXXXX]_ and necessarily
  [1805.07306_].

  FlexibleSUSY depends on SARAH_ and contains components from
  SOFTSUSY_. Therefore, please also cite the following publications
  along with FlexibleSUSY:

  - **SARAH** [0909.2863]_, [1002.0840]_, [1207.0906]_, [1309.7223]_
  - **SOFTSUSY** [hep-ph:0104145]_, [1311.7659]_

  The list of references in BibTeX format can be found in
  `<doc/references.bib>`_.

.. contents:: Table of Contents
   :depth: 2


Quick start
===========

Install required libraries and packages (if not already done)

.. code-block:: shell

    pip install conan
    conan install . --build=missing
    sudo apt-get install libgsl-dev
    ./install-sarah

Build a spectrum generator (here: HSSUSY [1710.03760]_
[1804.09410]_)

.. code-block:: shell

    ./createmodel --name=HSSUSY
    ./configure --with-models=HSSUSY
    make -j4

Run the spectrum generator

.. code-block:: shell

    ./models/HSSUSY/run_HSSUSY.x --slha-input-file=model_files/HSSUSY/LesHouches.in.HSSUSY


Building FlexibleSUSY
=====================

Requirements
------------

* C++17 compiler (g++ >= 7.0.0 or clang++ >= 4.0.0 or icpc >= 18.0.0)
* Fortran compiler (gfortran, ifort)
* `Mathematica/Wolfram Engine`_ (version 11.0 or higher)
* SARAH_ (version 4.11.0 or higher)
* Boost_ (version 1.37.0 or higher)
* `Eigen 3`_ (version 3.1 or higher)
* `GNU scientific library`_

Optional:

* FeynArts_ (version 3.10 or higher)
* FormCalc_ (version 9.8 or higher)
* GM2Calc_ (version 2.1.0 or higher)
* LoopTools_ (version 2.8 or higher)
* HiggsTools_ (version 1.1.4 or higher)
* Lilith_ (requires Python3 development files)
* COLLIER_
* Himalaya_
* TSIL_

Installation of required/optional libraries
-------------------------------------------

The required libraries Boost_, `Eigen 3`_ and `GNU scientific library`_ and the
optionally required Python3 development files can usually be installed via the
package manager of the operating system.
On Debian/Ubuntu one may run for example

.. code-block:: shell

    sudo apt-get install libboost-all-dev libeigen3-dev libgsl-dev python3-dev

The required and optional libraries Boost_, `Eigen 3`_, GM2Calc_,
LoopTools_, Himalaya_ and TSIL_ can be also installed using the Conan_
package manager.  If not already installed, Conan can be installed
with pip

.. code-block:: shell

    pip install conan

To install the libraries required by FlexibleSUSY, run

.. code-block:: shell

    conan install . --build=missing

The `GNU scientific library`_ can currently not be installed via
Conan_.

If the required libraries are installed via Conan or the operating
system's package manager, they will be found automatically by
FlexibleSUSY's ``configure`` script, see below.

Installation of SARAH
---------------------

FlexibleSUSY requires SARAH to be installed and to be loadable with
the ``Needs["SARAH`"]`` command from inside Mathematica.  We recommend
the following setup

.. code-block:: shell

    SARAH_VERSION=4.15.4
    cd ~/.Mathematica/Applications/
    wget https://sarah.hepforge.org/downloads/SARAH-${SARAH_VERSION}.tar.gz
    tar -xf SARAH-${SARAH_VERSION}.tar.gz
    ln -s ${PWD}/SARAH-${SARAH_VERSION}/ SARAH

    cd ~/.Mathematica/Kernel/
    echo "AppendTo[\$Path, \"${HOME}/.Mathematica/Applications/SARAH/\"];" >> init.m

All the above steps can be executed at once with the ``install-sarah``
script

.. code-block:: shell

    ./install-sarah

See ``./install-sarah --help`` for more options.

Installation of FeynArts/FormCalc (optional)
--------------------------------------------

If you want FlexibleSUSY to use FeynArts_ or FormCalc_ you will need
to install these packages first.  Also — as with SARAH — they need to
be loadable with the ``Needs[]`` command from inside Mathematica.  We
recommend using the installation script ``FeynInstall`` provided on
the FeynArts web page. e.g.

.. code-block:: shell

    cd ~/.local
    wget http://www.feynarts.de/FeynInstall
    chmod 755 FeynInstall
    ./FeynInstall

which will install the latest versions of FeynArts, FormCalc and
LoopTools in the ``~/.local/`` directory as well as configure
Mathematica to find these packages.  Note that running the
``FeynInstall`` script might require user intervention.

Building a FlexibleSUSY model
-----------------------------

0. Before you setup a FlexibleSUSY model, you have to provide a SARAH
   model file.  To make it available in FlexibleSUSY, you can put it
   either into FlexibleSUSY's SARAH model directory
   ``FlexibleSUSY/sarah/<model>/`` or directly into SARAH's own model
   directly ``SARAH/Models/<model>/``.  Here ``<model>`` is the name
   of your model (e.g. MSSM, NMSSM, etc.).  Note, that there are
   already plenty of pre-installed model files in FlexibleSUSY's and
   SARAH's model directories that can be used.

1. Create a new or re-initialize an existing FlexibleSUSY model

   .. code-block:: shell

       ./createmodel --name=<model>

   See ``./createmodel --help`` for more details.  Afterwards there will
   be

   * a model directory ``models/<model>/``
   * a makefile module ``models/<model>/module.mk``
   * a Mathematica start script ``models/<model>/start.m``
   * and a FlexibleSUSY model file ``models/<model>/FlexibleSUSY.m``

   To modify the model details (input parameters, boundary conditions,
   etc.), edit the FlexibleSUSY model file
   ``models/<model>/FlexibleSUSY.m``.  For more details see the
   documentation of the `FlexibleSUSY model file`_ and
   `FlexibleEFTHiggs`_.

2. Create the Makefile and register your model(s)

   .. code-block:: shell

       ./configure --with-models=<model>

   Multiple models can be specified, separated by a comma.  See
   ``./configure --help`` for more options.

3. Compile FlexibleSUSY with your model::

       make

   Use ``make -j<N>`` to use ``<N>`` CPU cores.  When ``make`` is
   executed, Mathematica is called, which generates the C++ code for
   the specified models.  All C++ source files are written to the
   directory ``models/<model>/``.  When ``make`` has finished, the
   following spectrum generator(s) are available for each specified
   model:

   * ``models/<model>/run_<model>.x``: command line spectrum generator
   * ``models/<model>/run_<model>.m``: Mathematica interface

Example

.. code-block:: shell

    ./createmodel --name=HSSUSY
    ./configure --with-models=HSSUSY
    make -j4

    ./models/HSSUSY/run_HSSUSY.x --slha-input-file=model_files/HSSUSY/LesHouches.in.HSSUSY


Using FlexibleSUSY
==================

Available models
----------------

FlexibleSUSY ships with many pre-generated models.  The following
table includes an (incomplete) list of models with a detailed
documentation.

======================== ====================================
 Model                    Description
======================== ====================================
 `HSSUSY`_                high-scale MSSM (pure EFT)
 `MSSMEFTHiggs`_          high-scale MSSM (FlexibleEFTHiggs)
 `NUHMSSMNoFVHimalaya`_   fixed-order MSSM
======================== ====================================

.. _`HSSUSY`: doc/models/HSSUSY.rst
.. _`MSSMEFTHiggs`: doc/models/MSSMEFTHiggs.rst
.. _`NUHMSSMNoFVHimalaya`: doc/models/NUHMSSMNoFVHimalaya.rst


Command line
------------

For each model FlexibleSUSY creates an executable
``models/<model>/run_<model>.x`` that can be run from the command
line.  The executable accepts the input in the SLHA format, for
example in form of a file

.. code-block:: shell

    ./models/MSSM/run_MSSM.x \
       --slha-input-file=models/MSSM/LesHouches.in.MSSM \
       --slha-output-file=LesHouches.out.MSSM

or as a stream

.. code-block:: shell

    cat models/MSSM/LesHouches.in.MSSM \
       | ./models/MSSM/run_MSSM.x --slha-input-file=- --slha-output-file=LesHouches.out.MSSM

For a documentation of FlexibleSUSY-specific switches in the SLHA
input see the section on `SLHA input parameters`_.

By default the executable writes the output in SLHA format to stdout.
The output can also be appended to an SQLite database

.. code-block:: shell

    ./models/MSSM/run_MSSM.x \
       --slha-input-file=models/MSSM/LesHouches.in.MSSM \
       --slha-output-file=LesHouches.out.MSSM \
       --database-output-file=points.db

See ``models/<model>/run_<model>.x --help`` for further options.


Mass spectrum and renormalization group running
```````````````````````````````````````````````

The pole mass spectrum and the RG flow can be written to text files
for easy plotting.  In the MSSM for example these text files can be
generated via

.. code-block:: shell

    ./models/MSSM/run_MSSM.x \
       --slha-input-file=model_files/MSSM/LesHouches.in.MSSM \
       --rgflow-output-file=MSSM_rgflow.dat \
       --spectrum-output-file=MSSM_spectrum.dat

The generated files ``MSSM_rgflow.dat`` and ``MSSM_spectrum.dat`` can
be plotted for example with the gnuplot scripts in the model
directory

.. code-block:: shell

    gnuplot -persist -e "filename='MSSM_spectrum.dat'" \
       models/MSSM/MSSM_plot_spectrum.gnuplot

    gnuplot -persist -e "filename='MSSM_rgflow.dat'" \
       models/MSSM/MSSM_plot_rgflow.gnuplot

The gnuplot scripts are just for illustration and currently plot all
running parameters, regardless of their mass dimension, so the
resulting plot is not particularly informative.  However, one may
easily adapt the scripts to plot any chosen subset of the parameters.


Mathematica interface
---------------------

FlexibleSUSY can be called from within Mathematica using Wolfram's
LibraryLink.  By default, FlexibleSUSY creates a LibraryLink library
for each spectrum generator.  The generated library can be found in
``models/<model>/<model>_librarylink.so``, where ``<model>`` is the
model name.

Example

.. code-block:: mathematica

    Get["models/CMSSM/CMSSM_librarylink.m"];

    (* Create a handle to a model given the input parameters.
       See Options[FSCMSSMOpenHandle] for all default options. *)
    handle = FSCMSSMOpenHandle[
      fsSettings -> { precisionGoal -> 1.*^-4 },
      fsSMParameters -> { Mt -> 173.3 },
      fsModelParameters -> {
          m0 -> 125, m12 -> 500, TanBeta -> 10, SignMu -> 1, Azero -> 0 }
    ];

    (* calculate pole mass spectrum *)
    FSCMSSMCalculateSpectrum[handle];

    (* calculate observables *)
    FSCMSSMCalculateObservables[handle];

    (* close the model handle *)
    FSCMSSMCloseHandle[handle];

For each model, FlexibleSUSY creates an example Mathematica script
which illustrates the use of the Mathematica interface.  The generated
example can be found in ``models/<model>/run_<model>.m`` which can be
run for example as

.. code-block:: shell

    math -run "<< \"models/<model>/run_<model>.m\""

Before running it, the model parameters in the script should be set to
reasonable values.  More advanced examples can be found in the
FlexibleSUSY documentation.

Note: In order to compile the library, Mathematica must be installed.
To disable the LibraryLink interface, configure with
``--disable-librarylink``.

Further details and examples can be found in the `LibraryLink
documentation`_.

.. _`LibraryLink documentation`: doc/librarylink.rst

Parameter scans
---------------

FlexibleSUSY contains two shell scripts aiming to help the user
performing parameter scans based on SLHA files.

Tabular output
``````````````

The script ``utils/scan-slha.sh`` performs a scan over an input
parameter.

Examples:

To perform a scan over :math:`\tan\beta(M_Z)` in the CMSSM (given in
the SLHA input file in the ``MINPAR[3]`` field) and print out the the
values of :math:`\tan\beta(M_Z)`, :math:`M_h` (``MASS[25]``) and
:math:`y_t(M_{\text{SUSY}})` (``YU[2,2]``) run

.. code-block:: shell

     utils/scan-slha.sh \
        --spectrum-generator=models/CMSSM/run_CMSSM.x \
        --slha-input-file=model_files/CMSSM/LesHouches.in.CMSSM \
        --scan-range=MINPAR[3]=1~30:10 \
        --output=MINPAR[3],MASS[25],YU[2:2]

Alternatively, the SLHA input can be piped into the script as

.. code-block:: shell

    cat model_files/CMSSM/LesHouches.in.CMSSM \
       | utils/scan-slha.sh \
         --spectrum-generator=models/CMSSM/run_CMSSM.x \
         --scan-range=MINPAR[3]=1~30:10 \
         --output=MINPAR[3],MASS[25],YU[2:2]

The spectrum generator executable is specified using the
``--spectrum-generator=`` option.  The parameter to be scanned over as
well as the scan range and the number of steps must be specified using
the ``--scan-range=`` option.  The syntax is

.. code-block:: shell

    --scan-range=<block>[<field>]=<start>~<stop>:<number_of_steps>

Here ``<block>`` is the SLHA block in which the input parameter is to
be found and ``<field>`` is the block entry corresponding to the
parameter.  ``<start>`` and ``<stop>`` define the scan range and
``<number_of_steps>`` define the number of steps.  By default the step
size is linear.  Alternatively, a logarithmic step size can be chosen
by passing ``--step-size=log`` to the script.  See also
``utils/scan-slha.sh --help``.  The parameters to print to the output
stream must be defined using the ``--output=`` option.  The syntax
is

.. code-block:: shell

    --output=<block>[<fields>]

where ``<block>`` is the SLHA block in which the output parameter is to
be read from and ``<field>`` is the block entry corresponding to the
parameter.  To read a matrix element from a block, use a colon ``:`` to
specify the matrix element indices.  Multiple output parameters can be
specified by a comma.

Database output
```````````````

As an alternative, all parameters calculated during a scan can be
written to a SQLite database using the ``scan-database.sh`` script.

Examples

.. code-block:: shell

    utils/scan-database.sh \
       --spectrum-generator=models/CMSSM/run_CMSSM.x \
       --slha-input-file=model_files/CMSSM/LesHouches.in.CMSSM \
       --scan-range=MINPAR[3]=1~30:10 \
       --database-output-file=scan.db

or

.. code-block:: shell

    cat model_files/CMSSM/LesHouches.in.CMSSM \
       | ./utils/scan-database.sh \
         --spectrum-generator=models/CMSSM/run_CMSSM.x \
         --scan-range=MINPAR[3]=1~30:10 \
         --database-output-file=scan.db

The name of the database file must be set using the
``--database-output-file=`` option.

Convert SPheno to FlexibleSUSY model file
-----------------------------------------

The script ``utils/convert_SPheno_to_FlexibleSUSY.m`` can help to
convert a SPheno model file (``SPheno.m``) to a FlexibleSUSY model
file (``FlexibleSUSY.m.in``).  The conversion is not perfect, because
it is usually not unique.  Therefore one should check the generated
``FlexibleSUSY.m.in`` file.

Example

.. code-block:: shell

    cat << EOF | math -noprompt > FlexibleSUSY.m.in
    sphenoFile = "~/.Mathematica/Applications/SARAH/Models/MSSM/SPheno.m";
    Get["utils/convert_SPheno_to_FlexibleSUSY.m"];
    EOF


Advanced FlexibleSUSY build options
===================================

Generating source code files only (no compilation)
----------------------------------------------------

If you want to only create the C++ source files for your model, but do
not want to compile the code, you can use the ``--disable-compile``
configure option

.. code-block:: shell

    ./configure --with-models=MSSM --disable-compile
    make

Here, configure will not check for installed compilers or libraries.
It will only search for Mathematica and SARAH.  The execution of
``make`` will stop as soon as all C++ source code files are generated.
See below for how to export the generated source code.


Compile only (don't generate source code)
-----------------------------------------

If you want to only compile already created the C++ source files for
your model, you can use the ``--disable-meta`` configure option

.. code-block:: shell

    ./configure --with-models=MSSM --disable-meta
    make

Here, configure will only check for installed compilers or libraries.
It will not check for Mathematica and SARAH.

This option is useful if you want to generate the source code on one
computer and then transfer the generated code to another computer to
compile it.  This option can also be used with the pre-generated
FlexibleSUSY models, which are provided at the FlexibleSUSY home page.

Warning: Please make sure all C++ source files of your model are
available in the model directory ``models/<model>/``.  Otherwise the
compilation will fail.


Exporting the generated source code
-----------------------------------

The complete FlexibleSUSY source code, including the generated C++
code for the specified model(s) (but without the Mathematica meta
code), can be exported to a new directory.  The exported source code
is a complete standalone package, with it's own build system.  To
export the code, one has to set the target directory during
configuration via the ``--with-install-dir=`` option.  For example

.. code-block:: shell

    ./configure --with-models=<models> --with-install-dir=/path/to/export/directory

Afterwards

.. code-block:: shell

    make install-src

must be executed, which will copy the generated C++ source code for
all ``<models>`` to ``/path/to/export/directory``, together with the
non-model specific source code from ``config/``, ``doc/``, ``slhaea/``
and ``src/``.  Afterwards, the standalone package can be build like
this

.. code-block:: shell

    cd /path/to/export/directory
    ./configure
    make

It is also possible to create a "model package", which includes only
the generated source code for a given model, but does not contain the
whole FlexibleSUSY build system.  This is useful when the source code
for a model should be generated on one computer and later transferred
to another one to be compiled.  To create such a "model package" run

.. code-block:: shell

    make pack-<model>-src

where ``<model>`` is the name of the model whose generated source code
shall be packed.  After ``make`` has finished, the package file
``<model>.tar.gz`` can be found in the working directory.


Dynamic libraries
-----------------

If you want to create dynamic model libraries (instead of static
libraries, which is the default) you need to pass the
``--enable-shared-libs`` option to the configure script.  The file
name extension for the shared libraries as well as the command to
build them can be overwritten using the ``--with-shared-lib-ext=``
``--with-shared-lib-cmd=``.  parameters.  For example, when Intel
compilers should be used, replace gcc by icc or icpc

.. code-block:: shell

    ./configure --with-models=CMSSM,NMSSM \
       --enable-shared-libs \
       --with-shared-lib-ext=".so" \
       --with-shared-lib-cmd="gcc -shared -o"

**Important remark:**

The libraries are linked to the executables with *absolute* paths.
This means that, if you for example move the FlexibleSUSY directory to
another location, the executables will no longer find the libraries.
To make the executables find the libraries again, you have to relink
them via

.. code-block:: shell

    make clean-executables
    make allexec


Statically linked executables
-----------------------------

External libraries can be linked statically to the spectrum generator
executables by passing ``--enable-static`` to configure.  This is
useful when the executable should be transferred to another computer,
where some libraries are not available.

Example

.. code-block:: shell

    ./configure --with-models=CMSSM --enable-static

If ``--enable-static`` is used, the following linker flags and
additional libraries will be used

.. code-block:: shell

    LDFLAGS = -static
    LDLIBS  = -ldl

These linker-specific flags and additional libraries can be
overwritten using ``--with-static-ldflags=`` and
``--with-static-ldlibs=``

Example

.. code-block:: shell

    ./configure --with-models=CMSSM \
       --enable-static \
       --with-static-ldflags="-static" \
       --with-static-ldlibs="-lquadmath -ldl"

In case of dynamic linking (``--disable-static``, which is the default),
the options ``--with-shared-ldflags=`` and ``--with-shared-ldlibs=`` must
be used to set ``LDFLAGS`` and ``LDLIBS``.


Support for alternative loop libraries
--------------------------------------

FlexibleSUSY ships with its own implementation of the
Passarino-Veltman 1-loop functions, which have been translated from
SOFTSUSY_.  However, alternative implementations of the 1-loop
functions can be used:

* LoopTools_
* COLLIER_
* FFlite (a thread-safe variant of LoopTools_, shipped with FlexibleSUSY)

The loop function libraries can be enabled by passing
``--with-loop-libraries=`` to the ``configure`` script

.. code-block:: shell

    ./configure --with-loop-libraries=<libraries>

where ``<libraries>`` can be any (or a combination) of ``collier``,
``looptools`` or ``fflite``.

Example

.. code-block:: shell

    ./configure --with-loop-libraries=collier,looptools

When the SLHA input is used, the loop library to use can be selected
by setting the entry of ``FlexibleSUSY[31]`` to ``0`` (= SOFTSUSY),
``1`` ( = COLLIER), ``2`` (= LoopTools) or ``3`` (= FFlite).  See
`SLHA input parameters`_ for details.

Example::

    Block FlexibleSUSY
       31   0    # loop library (0 = SOFTSUSY, 1 = COLLIER, 2 = LoopTools, 3 = FFlite)

When the Mathematica interface is used, the loop library to use can be
selected by setting the value of ``loopLibrary`` appropriately

.. code-block:: mathematica

    FS@ModelName@OpenHandle[
        fsSettings -> {
            loopLibrary -> 0   (* 0 = SOFTSUSY, 1 = COLLIER, 2 = LoopTools, 3 = FFlite *)
        }
    ]

In the following it is described in more detail how to enable these
alternative loop function libraries in FlexibleSUSY.

LoopTools support
`````````````````

It is possible to use LoopTools_ for calculating the loop functions,
instead of using SOFTSUSY's loop functions.  To enable LoopTools,
configure FlexibleSUSY via

.. code-block:: shell

    ./configure --enable-looptools

or

.. code-block:: shell

    ./configure --with-loop-libraries=looptools

If LoopTools has been installed via Conan_, the configure will
automatically find the paths to the LoopTools library.

To use the LoopTools library and header files from a specific
directory, run ``configure`` via

.. code-block:: shell

    LOOPTOOL_DIR=/path/to/looptools/build

    ./configure --enable-looptools \
       --with-looptools-incdir=$LOOPTOOLS_DIR \
       --with-looptools-libdir=$LOOPTOOLS_DIR

Note: LoopTools 2.8 or higher is required.
Also, if FlexibleSUSY is compiled with LibraryLink (default) then LoopTools has to be compiled with ``-fPIC`` option.
This is achieved by setting the ``FFLAGS`` variable during LoopTools configuration as

.. code-block:: shell

    FFLAGS="-O3 -fPIC" ./configure

COLLIER support
```````````````

It is possible to use COLLIER_ for calculating the loop functions,
instead of using SOFTSUSY's loop functions.  To enable COLLIER
configure FlexibleSUSY via

.. code-block:: shell

   ./configure --with-loop-libraries=collier

To use the COLLIER library and header files from a specific
directory configure via

.. code-block:: shell

    COLLIER_DIR=/path/to/COLLIER-x.y.z

    ./configure --with-loop-libraries=collier \
       --with-collier-incdir=$COLLIER_DIR/modules \
       --with-collier-libdir=$COLLIER_DIR

Note: versions since COLLIER-1.2.3 were tested so far.
Also, COLLIER static library should be configured with
``-Dstatic=ON -DCMAKE_POSITION_INDEPENDENT_CODE=ON`` flags.

TSIL support
````````````

Some models of FlexibleSUSY require TSIL_, for example `HSSUSY`_.  When
such models are activated (via ``./configure --with-models=<model>``),
FlexibleSUSY requires TSIL to be available.  If TSIL is installed in a
system directory or installed via Conan_, FlexibleSUSY will find the
TSIL automatically.  To use TSIL from a a non-standard directory,
configure FlexibleSUSY like this

.. code-block:: shell

    $TSIL_DIR=/path/to/tsil

    ./configure --enable-tsil \
       --with-tsil-incdir=$TSIL_DIR \
       --with-tsil-libdir=$TSIL_DIR

Note also that TSIL must be compiled with ``-fPIC``, which can be
achieved by setting in the TSIL ``Makefile``

.. code-block:: shell

    TSIL_OPT = -O3 -funroll-loops -fPIC

Linking with other external libraries
-------------------------------------

HiggsTools support
``````````````````

It is possible to check consistency of a BSM Higgs sector with experimental
data using HiggsTools_. To enable HiggsTools, configure FlexibleSUSY as

.. code-block:: shell

    HIGGSTOOLS_DIR=/path/to/higgstools
    ./configure --with-higgstools-incdir=$HIGGSTOOLS_DIR/include \
       --with-higgstools-libdir=$HIGGSTOOLS_DIR/lib64

Note that depending on the system, HiggsTools library is installed either in
`lib` or `lib64` directory. At runtime user has to provide directories
containing HiggsSignals and HiggsBounds datasets (which are downloaded
separatelly from HiggsTools) as

.. code-block:: shell

    models/<model>/run_<model>.x --higgsbounds-dataset=/path/to/HB/datasetdir \
       --higgssignals-dataset=/path/to/HS/datasetdir

Lilith support
``````````````

It is possible to check consistency of a BSM Higgs sector with experimental
data using Lilith_. To enable Lilith, configure FlexibleSUSY as

.. code-block:: shell

    ./configure --with-lilith=/path/to/Lilith

At runtime user can provide the Lilith database `.list` file

.. code-block:: shell

    models/<model>/run_<model>.x --lilith-db=/path/to/Lilith/database/file

If not provided, the default file will be used.

Creating an addon
-----------------

A FlexibleSUSY addon is a program or library, which uses parts of the
FlexibleSUSY libraries or the generated models or is integrated into
FlexibleSUSY.  An addon can be created via

.. code-block:: shell

    ./createaddon --name=<addon>

where ``<addon>`` is the name of the addon.  The createaddon script
creates the directory ``addons/<addon>/`` and the corresponding makefile
module ``addons/<addon>/module.mk``.  If an addon has been created with
the above script, the user may edit the makefile module
(``addons/<addon>/module.mk``) to add source files in to the three
variables

.. code-block:: shell

    LIB@ADDON@_SRC  # list of source files to be included in library
    EXE@ADDON@_SRC  # list of source files with a main()
    LIB@ADDON@_HDR  # list of header files

Example

.. code-block:: shell

    LIB@ADDON@_SRC := $(DIR)/file1.cpp
    EXE@ADDON@_SRC := $(DIR)/run.cpp
    LIB@ADDON@_HDR := $(DIR)/file1.hpp

To configure and compile the addon run

.. code-block:: shell

    ./configure --with-addons=<addon>
    make

make compiles all source files and creates the addon library
``addons/<addon>/lib<addon>.a`` (including the object file ``file1.o`` in
the above example) and an executable (``addons/<addon>/run.x`` in the
above example).


Creating the source code documentation
--------------------------------------

FlexibleSUSY's source code documentation (including the generated
source code files) can be generated with Doxygen in HTML or man
format.  To generate the HTML documentation please run

.. code-block:: shell

    make doc-html

The generated HTML index file can then be found in
``doc/html/index.html`` and can be viewed with any HTML browser, e.g.

.. code-block:: shell

    firefox doc/html/index.html

To generate the man documentation please run

.. code-block:: shell

    make doc-man

The generated man pages can then be found in ``doc/man/man3/`` and can
be viewed as

.. code-block:: shell

    man doc/man/man3/model_file_options.3


Cleaning
--------

There are several make targets to remove generated files, compiled
object files, libraries or executables::

    make clean      # deletes all .d .o .a .x files

    make distclean  # does `clean` and `clean-generated`
                    # and deletes in addition:
                    # Makefile flexiblesusy-config config.*
                    # config/list_sarah_model_files.sh

    make clean-dep  # deletes all .d files

    make clean-executables # deletes all .x files

    make clean-generated   # deletes generated files

    make clean-lib  # deletes all libraries

    make clean-obj  # deletes all .o files

For each model ``<model>`` or addon there are specific clean targets
to remove model-specific files::

    make clean-<model>     # deletes .d .o .a .x files

    make distclean-<model> # same as `make clean-<model> clean-<model>-src`

    make clean-<model>-dep # deletes .d files

    make clean-<model>-lib # deletes model library

    make clean-<model>-obj # deletes .o files

    make clean-<model>-src # deletes generated files


Package content
===============

In the following all sub-directories within the FlexibleSUSY package
are listed:

* ``addons/`` contains addons for FlexibleSUSY

* ``config/`` contains helper scripts and makefile modules for the
  build system

* ``doc/`` contains the FlexibleSUSY documentation

* ``examples/`` contains examples how to build you own spectrum
  generator based on FlexibleSUSY

* ``fflite/`` contains an alternative implementation of the
  Passarino-Veltman loop functions, based on FF

* ``meta/`` contains the Mathematica meta code which generates the
  spectrum generators.  See the `meta code documentation`_ for more
  details.

* ``model_files/`` contains default model files for some frequently
  used models (SM, SplitMSSM, MSSM, NMSSM, SMSSM, UMSSM, etc.)

* ``model_specific/`` contains model-specific higher order corrections
  for the MSSM, NMSSM, SM and SplitMSSM from the literature

* ``models/`` This is the output directory where the generated C++
  code for the spectrum generators will be stored.

* ``nlohmann/`` contains the nlohmann_ JSON library

* ``Output/`` contains SARAHs model-specific output files

* ``sarah/`` contains SARAH model files shipped with FlexibleSUSY

* ``slhaea/`` contains the slhaea_ SLHA reader library

* ``src/`` contains model-independent FlexibleSUSY C++ source code

* ``templates/`` contains C++ template files for the spectrum generators

* ``test/`` contains the FlexibleSUSY test suite

* ``utils/`` contains some utility scripts to perform scans or extract
  data from SLHA files


Further reading
===============

* `FlexibleSUSY model file`_
* `FlexibleEFTHiggs`_
* `LibraryLink documentation`_
* `meta code documentation`_
* `SLHA input parameters`_
* `Observables`_


References
==========

.. _slhaea: https://github.com/fthomas/slhaea
.. _GM2Calc: https://github.com/GM2Calc/GM2Calc
.. _SARAH: http://sarah.hepforge.org
.. _SOFTSUSY: http://softsusy.hepforge.org
.. _Boost: http://www.boost.org
.. _Conan: https://conan.io/
.. _Eigen 3: http://eigen.tuxfamily.org
.. _FeynArts: http://www.feynarts.de
.. _FormCalc: http://www.feynarts.de/formcalc
.. _GNU scientific library: http://www.gnu.org/software/gsl
.. _HiggsTools: https://gitlab.com/higgsbounds/higgstools
.. _Lilith: https://github.com/sabinekraml/Lilith-2
.. _`Mathematica/Wolfram Engine`: https://www.wolfram.com
.. _LoopTools: http://www.feynarts.de/looptools
.. _COLLIER: https://collier.hepforge.org
.. _Himalaya: https://github.com/Himalaya-Library/Himalaya
.. _TSIL: https://www.niu.edu/spmartin/tsil
.. _nlohmann: https://github.com/nlohmann/json

.. _`FlexibleSUSY model file`: doc/model_file.rst
.. _`FlexibleEFTHiggs`: doc/FlexibleEFTHiggs.rst
.. _`meta code documentation`: doc/meta_code.rst
.. _`SLHA input parameters`: doc/slha_input.rst
.. _`Observables`: doc/observables.rst

.. [hep-ph:0104145] `CPC 143 (2002) 305-331 <https://inspirehep.net/record/555481>`_ [`arxiv:hep-ph/0104145 <http://arxiv.org/abs/hep-ph/0104145>`_]
.. [0909.2863] `CPC 181 (2010) 1077-1086 <https://inspirehep.net/record/831371>`_ [`arxiv:0909.2863 <http://arxiv.org/abs/0909.2863>`_]
.. [1002.0840] `CPC 182 (2011) 808-833 <https://inspirehep.net/record/845241>`_   [`arxiv:1002.0840 <http://arxiv.org/abs/1002.0840>`_]
.. [1005.5709] `JHEP 1008 (2010) 104 <https://inspirehep.net/record/856612>`_  [`arxiv:1005.5709 <https://arxiv.org/abs/1005.5709>`_]
.. [1207.0906] `CPC 184 (2013) 1792-1809 <https://inspirehep.net/record/1121136>`_ [`arxiv:1207.0906 <http://arxiv.org/abs/1207.0906>`_]
.. [1211.2099] `Eur. Phys. J. C73 (2013) no. 2, 2310 <https://inspirehep.net/literature/1201957>`_ [`arxiv:1211.2099 <https://arxiv.org/abs/1211.2099>`_]
.. [1309.7223] `CPC 185 (2014) 1773-1790 <https://inspirehep.net/record/1255845>`_ [`arxiv:1309.7223 <http://arxiv.org/abs/1309.7223>`_]
.. [1311.7659] `CPC 185 (2014) 2322 <https://inspirehep.net/record/1266808>`_  [`arxiv:1311.7659 <http://arxiv.org/abs/1311.7659>`_]
.. [1406.2319] `CPC 190 (2015) 139-172 <https://inspirehep.net/record/1299998>`_ [`arxiv:1406.2319 <https://arxiv.org/abs/1406.2319>`_]
.. [1510.08071] `Eur. Phys. J. C76 (2016) no. 2, 62 <https://inspirehep.net/literature/1401235>`_ [`arXiv:1510.08071 <https://arxiv.org/abs/1510.08071>`_]
.. [1609.00371] `JHEP 1701 (2017) 079 <https://inspirehep.net/record/1484857>`_ [`arxiv:1609.00371 <https://arxiv.org/abs/1609.00371>`_]
.. [1708.05720] `Eur. Phys. J. C77 (2017) no. 12, 814 <https://inspirehep.net/record/1617767>`_ [`arxiv:1708.05720 <https://arxiv.org/abs/1708.05720>`_]
.. [1710.03760] `CPC 230 (2018) 145-217 <https://inspirehep.net/record/1629978>`_ [`arXiv:1710.03760 <https://arxiv.org/abs/1710.03760>`_]
.. [1804.09410] `Eur. Phys. J. C78 (2018) no. 7, 573 <https://inspirehep.net/record/1670032>`_ [`arxiv:1804.09410 <https://arxiv.org/abs/1804.09410>`_]
.. [1805.07306] `Eur. Phys. J. C78 (2018) no. 8, 649 <https://inspirehep.net/literature/1673989>`_ [`arxiv:1805.07306 <https://arxiv.org/abs/1805.07306>`_]
.. [1807.03509] `Eur. Phys. J. C78 (2018) no. 10, 874 <https://inspirehep.net/record/1681658>`_ [`arxiv:1807.03509 <https://arxiv.org/abs/1807.03509>`_]
.. [1910.03595] `Eur. Phys. J. C80 (2020) no. 3, 186 <https://inspirehep.net/record/1758261>`_ [`arxiv:1910.03595 <https://arxiv.org/abs/1910.03595>`_]
.. [2106.05038] `CPC 283 (2023) 108584 <https://inspirehep.net/literature/1867840>`_ [`arxiv:2106.05038 <http://arxiv.org/abs/2106.05038>`_]
.. [2110.13238] `Eur. Phys. J. C82 (2022) no. 3, 229 <https://inspirehep.net/literature/1952331>`_ [`arxiv:2110.13238 <https://arxiv.org/abs/2110.13238>`_]
.. [2204.05285] `Phys. Rev. D 106 (2022) 9, 095023 <https://inspirehep.net/literature/2065408>`_ [`arxiv:2204.05285 <https://arxiv.org/abs/2204.05285>`_]
.. [2402.14630] `CPC 302 (2024) 109244 <https://inspirehep.net/literature/2760779>`_ [`arxiv:2402.14630 <https://arxiv.org/pdf/2402.14630.pdf>`_]
.. [2210.09332] `CPC 291 (2023) 108803 <https://inspirehep.net/literature/2166780>`_ [`arxiv:2210.09332 <https://arxiv.org/pdf/2210.09332.pdf>`_]
.. [1502.04138] `Eur. Phys. J. C75 (2015) no. 9, 440 <https://inspirehep.net/literature/1345039>`_ [`arxiv:1502.04138 <https://arxiv.org/pdf/1502.04138.pdf>`_]
.. [1908.03952] `SciPost Phys. 7 (2019) 4, 052 <https://inspirehep.net/literature/1749361>`_ [`arxiv:1908.03952 <https://arxiv.org/pdf/1908.03952.pdf>`_]
.. [2012.11408] [`arxiv:2012.11408 <https://arxiv.org/pdf/2012.11408.pdf>`_]
.. [2502.XXXXX] `CPC`_ [`arxiv:2502.XXXX <>`_]
