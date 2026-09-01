# HolonomicContinuation

Website: [risc.jku.at/sw/holonomiccontinuation](https://risc.jku.at/sw/holonomiccontinuation)

## Overview

HolonomicContinuation is a Mathematica package for analytic continuation of holonomic functions. Given a holonomic function defined by a differential equation and initial values at a specific point, the program can compute a (logarithmically modulated) truncated series expansion of this function at another point. It is specifically designed to efficiently deal with exceptionally large differential equations. Certain functions also have an optional C backend based on FLINT [flintlib.org](https://flintlib.org). 

## Requirements

HolonomicContinuation is developed and tested for Wolfram Mathematica Versions 13-15, but should work on any resent version of Mathematica.
For the optional C backend, FLINT 3.4.0 [flintlib.org](https://flintlib.org) must be installed. 

## Installation

Clone repository and move the directory HolonomicContinuation into one of the directories in $Path where Mathematica can find it. 
To use the C-backends rec_to_val_V1.c or rec_to_val_V2.c navigate to the directory HolonomicContinuation and compile e.g. with

    gcc -o "rec_to_val_V2" "rec_to_val_V2.c" -O3  -march=native  -lflint  -fopenmp -I/path/to/flint/include -L/path/to/flint/lib  -Wl,-rpath,"/path/to/flint/lib"

where the flags -I -L and -Wl,-rpath, give, if necessary, the path to a suitable version of FLINT. Similar for rec_to_val_V1.c. 

## Documentation

See the Mathematica notebook RunningExample.nb and the usage-messages of the individual functions available e.g. via ?GetCoeffSubsFast

## Authors

Abilio De Freitas

Jakob Obrovsky (jakob.obrovsky@risc.jku.at)

Carsten Schneider (carsten.schneider@risc.jku.at) 

## License

HolonomicContinuation is distributed under GPL (GNU General Public License) version 3 or later. See the COPYING file.
