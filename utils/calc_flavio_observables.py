#!/usr/bin/env python3
"""Calculates flavour observables with flavio from the Wilson coefficients
that FlexibleSUSY writes to WC_<model>.json (WCxf format, flavio basis).

Usage:

   calc_flavio_observables.py WC_<model>.json [observable ...]

An observable that depends on arguments, e.g. a q^2 bin, is given as
name:arg1:arg2, e.g. '<Rmue>(B+->Kll):1.1:6'. Without observables,
BR(B->Xsgamma), BR(Bs->mumu) and R_K, R_K* in 1.1 < q^2 < 6 GeV^2 are
calculated. Each observable is printed with its Standard Model prediction
and the ratio of the two. Requires the flavio and wilson packages.

Coefficients that are missing in the file are taken as zero by flavio:
R_K and R_K* need the coefficients of b -> s mu mu and b -> s e e.
"""

import argparse
import sys

DEFAULT_OBSERVABLES = [
    'BR(B->Xsgamma)',
    'BR(Bs->mumu)',
    '<Rmue>(B+->Kll):1.1:6',
    '<Rmue>(B0->K*ll):1.1:6',
]


def parse_observable(text):
    """Splits 'name:arg1:arg2' into the name and the list of arguments"""
    name, *args = text.split(':')
    return name, [float(a) for a in args]


def main():
    """Main program"""
    parser = argparse.ArgumentParser(
        description=__doc__.split('\n\n', maxsplit=1)[0])
    parser.add_argument('json_file', help='WCxf file written by FlexibleSUSY')
    parser.add_argument('observables', nargs='*', default=DEFAULT_OBSERVABLES,
                        help='flavio observables (default: %(default)s)')
    args = parser.parse_args()

    import flavio
    from wilson import Wilson

    try:
        with open(args.json_file, 'r', encoding='utf-8') as json_file:
            coeffs = Wilson.load_wc(json_file)
    except OSError as error:
        sys.exit(f'error: cannot read {args.json_file}: {error}')

    for obs in args.observables:
        name, obs_args = parse_observable(obs)
        np_value = flavio.np_prediction(name, coeffs, *obs_args)
        sm_value = flavio.sm_prediction(name, *obs_args)
        print(f'{obs:24s} = {np_value:.6e}   SM: {sm_value:.6e}'
              f'   ratio: {np_value/sm_value:.6f}')


if __name__ == '__main__':
    main()
