"""Emit a JSON fixture of the fully-resolved KP15 n=6 bi-6 cap.

For one fixed-seed rational input (C2-consistent shared net + rational
free values) at n = 6, tau = 17/20, this resolves every dependent bi-6
BB-coefficient via the Maple-script construction (see n6_oracle.py) and
writes ALL control points of the deep interior block — SCRIPT indices
i, j in [3, 6] (= paper p[6-i, 6-j], i.e. paper i', j' in [0, 3]) — for
each of the 6 sectors, together with the inputs and free-value choices.

Values are emitted as exact rationals ("num/den") so a future Rust test
can assert bit-for-bit. The G1/G2/sal residuals (all exactly 0) are
recorded as provenance.

Run: <venv>/bin/python gen_rust_fixtures.py   (writes n6_cap_fixture.json)
"""

from __future__ import annotations

import json
import random

import sympy as sp

import n6_oracle as O

SEED = 20260712
N = 6
OUT = "n6_cap_fixture.json"


def rat(x):
    r = sp.nsimplify(x, rational=True)
    r = sp.Rational(r)
    return f"{r.p}/{r.q}"


def main():
    rng = random.Random(SEED)
    c = sp.Rational(1, 2)
    tau = sp.Rational(17, 20)
    free = O.rational_free(N, rng)
    net = O.shared_net(N, rng)
    P, sal, meta = O.build_cap(N, c, tau, net, free, simplify=False)

    # Residual provenance (must be exactly 0).
    g1max = g2max = sp.Integer(0)
    for s in range(N):
        g1, g2 = O.sector_residuals(P, s, N, c, tau)
        g1max = max(g1max, sp.Integer(0) if sp.simplify(g1) == 0 else sp.Integer(1))
        g2max = max(g2max, sp.Integer(0) if sp.simplify(g2) == 0 else sp.Integer(1))

    deep = []
    for s in range(N):
        block = {}
        for i in range(3, 7):
            for j in range(3, 7):
                block[f"{i},{j}"] = rat(P[s][(i, j)])
        deep.append(block)

    free_out = {k: (rat(v) if not isinstance(v, list) else [rat(x) for x in v]) for k, v in free.items()}
    net_out = [[[rat(x) for x in row] for row in blk] for blk in net]

    fixture = {
        "source": "tools/reference/kp15-bi6-cap/n6_oracle.py (Maple script docs/papers/n6.txt)",
        "n": N,
        "c": rat(c),
        "tau": rat(tau),
        "d2": rat(meta["d2"]),
        "seed": SEED,
        "index_convention": "SCRIPT p[i,j] = paper p[6-i,6-j]; [6,6] = extraordinary point",
        "deep_block_indices": "script i,j in [3,6] per sector (paper i',j' in [0,3])",
        "free_set": free_out,
        "input_net_4x4_per_sector": net_out,
        "residuals": {
            "G1_all_zero": bool(g1max == 0),
            "G2_all_zero": bool(g2max == 0),
            "sal_zero": bool(sp.simplify(sal) == 0),
        },
        "deep_block_cp_per_sector": deep,
    }

    with open(OUT, "w") as f:
        json.dump(fixture, f, indent=2)
    print(f"wrote {OUT}")
    print(f"  G1_all_zero={fixture['residuals']['G1_all_zero']} "
          f"G2_all_zero={fixture['residuals']['G2_all_zero']} "
          f"sal_zero={fixture['residuals']['sal_zero']}")
    print(f"  sectors={N} deep CPs/sector={len(deep[0])}")


if __name__ == "__main__":
    main()
