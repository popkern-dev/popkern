"""Emit all 19 regenerated PnS mask CSVs into crates/popkern-core/data/pns-masks/.

Run from the repo root:  ~/venvs/ocp/bin/python tools/reference/pns-tables/emit_all.py

Each table is re-verified against its generator's documented gate before
writing (the star/ngon macro family is self-derived and is NOT compared to
the reference package here — see eop_ngon.py and the data README).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np

import polar
import t_tables
import eop_ngon

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   '..', '..', '..', 'crates', 'popkern-core', 'data', 'pns-masks')

HEADER = ("# Regenerated PnS mask table (popkern, AGPL-3.0). Clean-room derivation\n"
          "# from the published papers; generators live in tools/reference/pns-tables/.\n")


def write(name, table):
    path = os.path.join(OUT, name + '.csv')
    with open(path, 'w') as f:
        f.write(HEADER)
        for row in table:
            f.write(','.join(f'{v:.12g}' for v in row) + '\n')
    print('wrote', name, table.shape)


def main():
    os.makedirs(OUT, exist_ok=True)
    for m in range(3, 9):
        write(f'polarSct{m}', polar.polar_mask(m))
    write('T0', t_tables.t0_mask())
    write('T1', t_tables.t1_mask())
    write('T2', t_tables.t2_mask())
    for n in (3, 5):
        write(f'ngonSct{n}', eop_ngon.ngon_small(n))
        write(f'eopSct{n}', eop_ngon.eop_small(n))
    for n in (6, 7, 8):
        write(f'ngonSct{n}', eop_ngon.ngon_macro(n))
        write(f'eopSct{n}', eop_ngon.eop_macro(n))


if __name__ == '__main__':
    main()
