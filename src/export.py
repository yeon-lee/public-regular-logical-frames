#!/usr/bin/env python3
"""Export binary checks and translated canonical pairs, with explicit scope."""
import argparse
import json
from pathlib import Path
from rlf import reconstruct, support, translate

ROOT = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('code', help='Catalogue ID')
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    catalogue = json.loads((ROOT / 'catalog.json').read_text())
    entries = catalogue['entries'] + catalogue.get('comparisons', [])
    match = next((e for e in entries if e['id'] == args.code), None)
    if match is None:
        parser.error('Unknown catalogue ID')
    code = json.loads((ROOT / match['path']).read_text())
    hx, hz, seeds = reconstruct(code)
    P = code['group']['order']
    represented = len(seeds) * P
    complete = represented == code['k']
    if not complete and represented != code['partial_frame']['pairs']:
        raise ValueError('Partial-frame size differs from the catalogue')
    output = {'id': code['id'], 'n': code['n'], 'k': code['k'], 'period': P,
              'complete_canonical_basis': complete,
              'represented_logical_pairs': represented,
              'residual_logical_pairs': code['k'] - represented,
              'frame_scope': 'full' if complete else 'partial regular frame; width does not bound full-code width',
              'coordinates': 'zero-based binary block*P+shift',
              'basis_order': 'seed pair, then shift=0,...,P-1',
              'HX': [support(v) for v in hx], 'HZ': [support(v) for v in hz],
              'X_seeds': [support(x) for x, _ in seeds],
              'Z_seeds': [support(z) for _, z in seeds],
              'X_basis': [support(translate(x, g, P)) for x, _ in seeds for g in range(P)],
              'Z_basis': [support(translate(z, g, P)) for _, z in seeds for g in range(P)]}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + '\n')
    print(f'Exported {code["id"]}: {represented} of {code["k"]} canonical logical pairs to {args.output}')


if __name__ == '__main__':
    main()
