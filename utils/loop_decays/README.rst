This directory contains utility needed to re-generate 1-loop decay amplitudes.
In normal circumstance there should be no reason to use it ever again.
The code uses FeynArts and FormCalc and is not compatible with FormCalc version > 9.6.
It should work with version < 9.6 but has not been tested with them.

# Technical details
Makefile expects that FormCalc can loaded via the
```
<< FormCalc`
```
command.
Please make sure to set it's path in `~/.Mathematica/Kernel/init.m` or apriopriate if you're on a system other than Linux.