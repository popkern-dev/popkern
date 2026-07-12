"""Cross-check the Rust bi-6 mu / gamma tables against the Maple script.

Transcribes the numeric tables from
  crates/popkern-core/src/geometry/gregory/kp15_bi6_mu.rs      (mu_dot,
      mu_ddot, mu_tdot, mu_tilde, mu_bar, w_scalars)
  crates/popkern-core/src/geometry/gregory/kp15_bi6_reparam.rs (d2, gamma)
and diffs them numerically against the script tables (mu1, mu2, mu3, mut,
mub) and the TBdata row-2 reparameterisation at n = 6, tau = 0.85.

Run: <venv>/bin/python rust_table_xcheck.py
"""

from __future__ import annotations

import sympy as sp

import n6_oracle as O

N = 6
C = sp.Rational(1, 2)
TAU = sp.Rational(17, 20)
D2 = sp.Rational(1, 2) * (3 * TAU**2 + TAU**2 * C - 3 * TAU - 2 * TAU * C + 1) / TAU**2


def approx(x):
    return float(sp.N(x, 20))


# ---- Rust mu tables (transcribed from kp15_bi6_mu.rs) --------------------
def rust_mu_dot(c, tau):
    return [2 * c / 5, 2 - c, c, -2 * c / 5, c / 15]


def rust_mu_ddot(c, tau):
    t2 = tau * tau
    return [
        -4 * c * t2 / 5,
        4 * (1 - c + c * c / 5) + 24 * c * (c - 1) / 5 * tau + 4 * c * c * t2 / 5,
        24 * c / 5 * tau - 4,
        2 * c / 3 * (11 - 3 * c - 6 * c * tau),
        -4 * c,
        2 * c * c,
        -4 * c * c / 5,
        2 * c * c / 15,
    ]


def rust_mu_tdot(c, tau):
    t2 = tau * tau
    t3 = t2 * tau
    return [
        4 * c * t3 / 25,
        4 * c / 25 * (5 - 3 * c * tau + 6 * (1 - c) * t2 - c * t3),
        -24 * c * t2 / 25,
        4 - 2 * c + 6 * c * (c - 2) / 5 * tau + 4 * c * c / 5 * t2,
        12 * c / 5 * tau - 4,
        2 * c / 5 * (5 - 3 * c * tau),
        4 * c / 25 * (3 * c * tau - 5),
        2 * c / 75 * (5 - 3 * c * tau),
        -3 * c / (5 * tau),
        6 * c * (5 - c * tau) / (125 * t2),
        c * (6 * c * tau - 5) / (125 * t3),
    ]


def rust_mu_tilde(c, tau):
    t2 = tau * tau
    t3 = t2 * tau
    return [
        sp.Rational(-11, 20),
        (11 * c - 10) / (8 * c),
        5 / (4 * c),
        sp.Rational(1, 8) - c / 4,
        3 * c / 10 - sp.Rational(1, 20),
        sp.Rational(1, 120) - c / 20,
        3 / (8 * tau),
        (3 * c * tau - 15) / (100 * t2),
        (5 - 6 * c * tau) / (200 * t3),
        -5 / (16 * c),
    ]


def rust_mu_bar(c, tau):
    t2 = tau * tau
    t3 = t2 * tau
    return [
        sp.Integer(1),
        5 * (2 - c) / (2 * c),
        -5 / c,
        sp.Rational(5, 2),
        sp.Integer(-1),
        sp.Rational(1, 6),
        3 / (10 * tau) * (c * tau - 5),
        3 / (125 * t2) * (25 - 20 * c * tau + c * c * t2),
        -1 / (250 * t3) * (25 - 45 * c * tau + 6 * c * c * t2),
        5 / (4 * c),
    ]


def rust_w(c, tau):
    return [-c / (50 * tau), c / 50, c / 15]


def rust_d2(c, tau):
    return (1 - (3 + 2 * c) * tau + (3 + c) * tau * tau) / (2 * tau * tau)


def diff_table(name, script_tab, rust_tab, idx_map):
    """script_tab: dict index->val (script indexing, incl [0] sum entry).
    rust_tab: list (rust storage, [0..]).  idx_map: list of (script_idx,
    rust_idx) pairs to compare."""
    print(f"\n== {name} ==")
    allok = True
    for si, ri in idx_map:
        sv = approx(script_tab[si])
        rv = approx(rust_tab[ri])
        d = abs(sv - rv)
        ok = d < 1e-15
        allok = allok and ok
        flag = "OK " if ok else "MISMATCH"
        print(f"  script[{si}] vs rust[{ri}]: {sv:+.12f}  {rv:+.12f}  d={d:.2e}  {flag}")
    return allok


def main():
    c, tau = C, TAU
    print(f"KP15 Rust-table cross-check at n={N}, c=1/2, tau=17/20, d2={approx(D2):.10f}")

    results = {}

    # mu1 / mu_dot: script[1..5] -> rust[0..4]
    results["mu_dot"] = diff_table(
        "mu_dot (mu1)", O.mu1_tab(c, tau), rust_mu_dot(c, tau),
        [(i, i - 1) for i in range(1, 6)],
    )
    # mu2 / mu_ddot: script[1..8] -> rust[0..7]
    results["mu_ddot"] = diff_table(
        "mu_ddot (mu2)", O.mu2_tab(c, tau), rust_mu_ddot(c, tau),
        [(i, i - 1) for i in range(1, 9)],
    )
    # mu3 / mu_tdot: script[1..11] -> rust[0..10]
    results["mu_tdot"] = diff_table(
        "mu_tdot (mu3)", O.mu3_tab(c, tau), rust_mu_tdot(c, tau),
        [(i, i - 1) for i in range(1, 12)],
    )
    # mut / mu_tilde: script[1..10] -> rust[0..9]
    results["mu_tilde"] = diff_table(
        "mu_tilde (mut)", O.mut_tab(c, tau), rust_mu_tilde(c, tau),
        [(i, i - 1) for i in range(1, 11)],
    )
    # mub / mu_bar: script[1..10] -> rust[0..9]
    results["mu_bar"] = diff_table(
        "mu_bar (mub)", O.mub_tab(c, tau), rust_mu_bar(c, tau),
        [(i, i - 1) for i in range(1, 11)],
    )

    # sum-entry reconstruction: verify rust reconstructs script[0]
    print("\n== sum-entry [0] reconstruction (rust sums stored entries) ==")
    mu1 = O.mu1_tab(c, tau)
    r = rust_mu_dot(c, tau)
    print(f"  mu_dot[0]  script={approx(mu1[0]):+.10f}  rust(2-sum)={approx(2 - sum(r)):+.10f}")
    mu2 = O.mu2_tab(c, tau)
    r = rust_mu_ddot(c, tau)
    print(f"  mu_ddot[0] script={approx(mu2[0]):+.10f}  rust(-sum)={approx(-sum(r)):+.10f}")
    mu3 = O.mu3_tab(c, tau)
    r = rust_mu_tdot(c, tau)
    print(f"  mu_tdot[0] script={approx(mu3[0]):+.10f}  rust(-sum1..8)={approx(-sum(r[0:8])):+.10f}")
    mut = O.mut_tab(c, tau)
    r = rust_mu_tilde(c, tau)
    print(f"  mut[0]     script={approx(mut[0]):+.10f}  rust(1-sum1..6)={approx(1 - sum(r[0:6])):+.10f}")
    mub = O.mub_tab(c, tau)
    r = rust_mu_bar(c, tau)
    print(f"  mub[0]     script={approx(mub[0]):+.10f}  rust(-sum1..6)={approx(-sum(r[0:6])):+.10f}")

    # d2
    print("\n== d2 (reparam.rs vs script eq 13) ==")
    dd = abs(approx(D2) - approx(rust_d2(c, tau)))
    print(f"  script d2={approx(D2):.12f}  rust d2={approx(rust_d2(c,tau)):.12f}  d={dd:.2e}  {'OK' if dd<1e-15 else 'MISMATCH'}")
    results["d2"] = dd < 1e-15

    # w_scalars: NOT present in n6.txt (script stops at S3; eq (11) p41
    # rule with w1/w2/w3 is not in this partial script) -> not cross-checkable.
    print("\n== w_scalars (eq 11) ==")
    print("  NOT in n6.txt (partial script stops at S3). Rust w_scalars =",
          [approx(x) for x in rust_w(c, tau)], "-- unchecked here.")

    # gamma vs TBdata row-2 reparameterisation (best-effort mapping).
    gamma_verdict = gamma_xcheck(c, tau)
    results["gamma"] = gamma_verdict

    print("\n===== SUMMARY =====")
    for k, v in results.items():
        print(f"  {k:10s}: {'MATCH' if v else 'MISMATCH'}")
    print("  ALL:", "MATCH" if all(results.values()) else "SOME MISMATCH")


def gamma_xcheck(c, tau):
    """gamma(k,i,j) [reparam.rs] maps the depth-{0,1,2} bi-3 tensor-border
    b_{ij} into the reparameterised row-2 CP q_{k,2} (script TBdata row 2,
    column k). We verify, symbolically over an independent-b symbol basis,
    that  TBdata_row2[k] == sum_{i,j} gamma(k,i,j) * b[i,j]  for k = 2..6.

    TBdata row 2 also has entries at k=0,1 not covered by the gamma table
    (paper eq (17) is stated for k in 2..6); those are reported separately.
    """
    print("\n== gamma table vs TBdata row 2 (eq 17), symbolic in b_ij ==")
    b = {(i, j): sp.Symbol(f"b_{i}{j}") for i in range(4) for j in range(3)}
    tb = O.TBdata(b, c, sp.Symbol("d2"), tau)  # keep d2 symbolic to match gamma's d2
    d2sym = sp.Symbol("d2")
    row2 = tb[2]  # 7 entries, index = column k (0..6)
    ok = True
    for k in range(7):
        gsum = 0
        for i in range(4):
            for j in range(3):
                g = gamma_rust(c, tau, d2sym, k, i, j)
                gsum += g * b[(i, j)]
        diff = sp.expand(row2[k] - gsum)
        z = sp.simplify(diff) == 0
        covered = 2 <= k <= 6
        tag = "MATCH" if z else ("MISMATCH" if covered else "(k<2: not in eq17 table)")
        if covered and not z:
            ok = False
        print(f"  k={k}: {'covered' if covered else 'uncovered'}  residual==0: {z}  {tag}")
        if covered and not z:
            print("     residual:", diff)
    return ok


def gamma_rust(c, tau, d2, k, i, j):
    """Direct transcription of gamma() from kp15_bi6_reparam.rs."""
    R = sp.Rational
    m = {
        (2, 0, 0): (1 - d2) / 25,
        (2, 0, 1): R(3, 25) + d2 / 50,
        (2, 0, 2): R(1, 25),
        (2, 1, 0): (6 + d2) / 50,
        (2, 2, 0): R(1, 25),
        (2, 1, 1): R(9, 25),
        (2, 1, 2): R(3, 25),
        (2, 2, 1): R(3, 25),
        (2, 2, 2): R(1, 25),
        (3, 0, 0): ((3 + 5 * c - 3 * d2) * tau - 6 * c) / (300 * tau),
        (3, 1, 0): ((9 - 2 * c - 7 * d2) * tau + 2 * c) / (100 * tau),
        (3, 2, 0): (9 + 3 * d2) / 100,
        (3, 3, 0): R(1, 100),
        (3, 0, 1): (18 + 2 * c + 3 * d2) / 600,
        (3, 1, 1): (54 + 9 * d2) / 200,
        (3, 2, 1): R(27, 100),
        (3, 3, 1): R(3, 100),
        (3, 0, 2): R(1, 100),
        (3, 1, 2): R(9, 100),
        (3, 2, 2): R(9, 100),
        (3, 3, 2): R(1, 100),
        (4, 1, 0): ((6 + 6 * c - 5 * d2) * tau - 8 * c) / (150 * tau),
        (4, 2, 0): ((9 - 4 * c - 5 * d2) * tau + 4 * c) / (75 * tau),
        (4, 3, 0): (2 + d2) / 50,
        (4, 1, 1): (18 + 2 * c + 3 * d2) / 150,
        (4, 2, 1): (18 + 3 * d2) / 50,
        (4, 3, 1): R(3, 25),
        (4, 1, 2): R(1, 25),
        (4, 2, 2): R(3, 25),
        (4, 3, 2): R(1, 25),
        (5, 2, 0): ((3 + c - 2 * d2) * tau - 2 * c) / (30 * tau),
        (5, 3, 0): ((3 - 2 * c - d2) * tau + 2 * c) / (30 * tau),
        (5, 2, 1): (18 + 2 * c + 3 * d2) / 60,
        (5, 3, 1): (6 + d2) / 20,
        (5, 2, 2): R(1, 10),
        (5, 3, 2): R(1, 10),
        (6, 3, 0): (6 - 2 * c - 3 * d2) / 30,
        (6, 3, 1): (18 + 2 * c + 3 * d2) / 30,
        (6, 3, 2): R(1, 5),
    }
    return m.get((k, i, j), sp.Integer(0))


if __name__ == "__main__":
    main()
