List of potential future FlexibleSUSY features
==============================================


Explicit path to SARAH installation
-----------------------------------

Allow the user to specify full path to the SARAH installation
(suggested by Florian).  This will modify the `Needs[]` command used
in all the scripts (`createmodel`, `configure`, `start.m`).


Decays
------

Comming soon.


Parameter output names
----------------------

Currently we're creating the parametr output names from their
Mathematica symbols.  For example in the MSSM we convert
::

    \[Mu]     ->   "Mu"
    B[\[Mu]]  ->   "BMu"

However, in SARAH the user can chose the output name in the model file
via::

    {{  Description -> "Mu-parameter",
        LaTeX -> "\\mu",
        ...
        OutputName-> Mu }},

    {{  Description -> "Bmu-parameter",
        LaTeX -> "B_{\\mu}",
        ...
        OutputName-> Bmu }},

We should use the user-defined output name, i.e.
::

    \[Mu]     ->   "Mu"
    B[\[Mu]]  ->   "Bmu"


Create a function which provides particle information
-----------------------------------------------------

Jae-hyeon proposed the following: FS could provide a function which
returns a list of (mass ordered) particles and their properies
(R-parity, Hypercharge, ...).  From this list a user can easily
extract for example the LSP (which would be the first particle in the
list with R-charge -1) or the weakly interacting particles etc.


Allow to not run VEVs up to the GUT scale
-----------------------------------------

In the E6SSM for example the VEV running from MZ to MX can become
non-perturbative.  A user interface for disabling the running of some
parameters between some scales might be useful here.  Maybe like this::

    DisableRGRunning = {
        {vu, LowScale, HighScale},
        {vd, LowScale, HighScale}
    };


Spectrum in SLHA convention in Mathematica interface
----------------------------------------------------

The `FS<model>CalculateSpectrum[handle]` function could be extended to
output the masses, mixing matrices and parameters in SLHA-2
convention.  A possible user interface could be::

    FS<model>CalculateSpectrum[handle, Convention -> HK]
    FS<model>CalculateSpectrum[handle, Convention -> SLHA]
