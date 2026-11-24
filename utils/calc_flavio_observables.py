#!/usr/bin/env python3
"""Calculates flavour observables with flavio from the Wilson coefficients
that FlexibleSUSY writes to WC_<model>.json (WCxf format, flavio basis).

Usage:

   calc_flavio_observables.py WC_<model>.json [observable ...]

Without observables, BR(B->Xsgamma) and BR(Bs->mumu) are
calculated. Each observable is printed with its Standard Model prediction
and the ratio of the two. Requires the flavio and wilson packages.
"""

import argparse
import sys

DEFAULT_OBSERVABLES = ['BR(B->Xsgamma)', 'BR(Bs->mumu)']


def main():
    """Main program"""
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('json_file', help='WCxf file written by FlexibleSUSY')
    parser.add_argument('observables', nargs='*', default=DEFAULT_OBSERVABLES,
                        help='flavio observables (default: %(default)s)')
    args = parser.parse_args()

    import flavio
    from wilson import Wilson

    try:
        with open(args.json_file, 'r') as json_file:
            coeffs = Wilson.load_wc(json_file)
    except OSError as error:
        sys.exit(f'error: cannot read {args.json_file}: {error}')

    for obs in args.observables:
        np_value = flavio.np_prediction(obs, coeffs)
        sm_value = flavio.sm_prediction(obs)
        print(f'{obs:20s} = {np_value:.6e}   SM: {sm_value:.6e}'
              f'   ratio: {np_value/sm_value:.6f}')


if __name__ == '__main__':
    main()
