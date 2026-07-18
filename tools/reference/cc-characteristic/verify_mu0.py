"""PK-134: derive the KP20 eq(23) mu0 table from first principles.

Printed (KP20 p.10): mu0 = 1 + kappa_n/10^5 with
  kappa_n = 185, 1033, 352, -145, -469, -687, -840 for n = 3,5,6,7,8,9,10.

Structure of eq (21)/(22): tau01 = vdot/3, tau02 = 2*vdot/3,
tau03 = mu0*vdot — a cubic row whose first two legs are exactly the
uniform-linear ones, with the whole cubic deviation carried by the
last leg. In monomial form b(t) = t + (mu0 - 1) * t^3 (b(0)=0,
b'(0)=1, b''(0)=0).

Hypotheses tested here, all computed from the characteristic ring
(Stam picking / PR08 cross-checked machinery in this directory):

  H1: mu0 from the ring-1 spoke-edge cubic. Take the corner patch's
      v=0 edge row (B-spline -> Bezier), express the (complex,
      frequency-1) edge cubic in its own parameter, normalize to the
      eq(21) chart and read the cubic coefficient.
  H2: mu0 from the subdominant eigenvalue alone (alpha = -ln lam/ln 2,
      simple rational candidates).
  H3: mu0 from a Taylor-style match of the power-law spoke u^alpha at
      the ring-1 outer boundary.

Run: python3 verify_mu0.py   (needs numpy; kp15venv works)

FINDINGS (2026-07-18, PK-134):
- All simple closed-form candidates FALSIFIED: H2a (power-law Taylor
  coefficient), H2b (lambda rationals), and every ring-1 edge-cubic
  ratio (closest: +0.00193 vs printed +0.00185 at n=3 — wrong by 4%,
  diverges at higher n).
- The paper (p.9, Fig. 17 + eq (21)) defines tau as a TOTAL-DEGREE-3
  triangular jet per sector (10 BB-coefficients: the five uniform
  eq(21) ones + the five mu-parameterized eq(22) ones), constructed
  "G2-connected using rho with the tensor-border of the
  characteristic ring of Catmull-Clark subdivision", G2 across
  sectors, internally C2, diagonally + rotationally symmetric.
  => mu0..mu3 are the NUMERIC SOLUTION of that per-valence G2
  fitting system, printed rounded to 5 decimals (hence the /10^5
  form). Reproducing the table = implementing the tau-sector G2
  solve against the ring tensor border (constraints (8,9) + the
  printed rho) in this oracle. That is the PK-134 continuation.
"""

import numpy as np
from verify_stam import build_A_Abar, char_net, picks

PRINTED = {3: 185e-5, 5: 1033e-5, 6: 352e-5, 7: -145e-5, 8: -469e-5,
           9: -687e-5, 10: -840e-5}

# Uniform cubic B-spline segment -> Bezier (segment between knots).
BSP2BEZ = np.array([
    [1/6, 4/6, 1/6, 0],
    [0,   4/6, 2/6, 0],
    [0,   2/6, 4/6, 0],
    [0,   1/6, 4/6, 1/6],
])


def edge_cubic_bezier(N, which=3):
    """Bezier points of the corner patch's spoke-edge cubic (complex
    frequency-1 representation; global constant cancels in ratios)."""
    lam, C = char_net(N)
    _, Abar = build_A_Abar(N)
    Cbar = Abar @ C
    P = picks(N)
    G = np.array([Cbar[i - 1] for i in P[which]]).reshape(4, 4)
    # v-rows: take the edge row closest to the spoke (see final_xcheck
    # comments; try both first rows and let the caller compare).
    return lam, [BSP2BEZ @ G[r, :] for r in range(2)], [BSP2BEZ @ G[:, c] for c in range(2)]


def mu0_from_cubic(bez):
    """Given 4 (complex) Bezier points of an edge cubic, project onto
    the dominant direction and read mu0 in the eq(21) chart:
    reparameterize so b(0)=0 (subtract), b'(0)=1 (scale), then a
    cubic b(t)=t+q2 t^2+q3 t^3; eq(21) has q2=0 — report both raw q3
    and the q3 after the affine reparameterization that kills q2."""
    b = np.array(bez) - bez[0]
    # dominant direction = end vector
    d = b[3]
    if abs(d) < 1e-14:
        return None
    s = (b / d).real  # scalar profile along the spoke, ends at 1
    # cubic through Bezier scalars s0..s3: monomial coefficients
    # b(t) = 3(s1-s0) t + 3(s2-2s1+s0) t^2 + (s3-3s2+3s1-s0) t^3
    c1 = 3 * (s[1] - s[0])
    c2 = 3 * (s[2] - 2 * s[1] + s[0])
    c3 = s[3] - 3 * s[2] + 3 * s[1] - s[0]
    if abs(c1) < 1e-14:
        return None
    # scale t so linear coeff is 1: t -> t/c1 gives 1, c2/c1^2, c3/c1^3
    q2 = c2 / c1**2
    q3 = c3 / c1**3
    # kill q2 by the quadratic-free reparameterization t = u - q2 u^2 + ...
    # to third order: composed cubic coefficient q3' = q3 - 2 q2^2
    q3_chartfree = q3 - 2 * q2**2
    return q2, q3, q3_chartfree


def main():
    print(f"{'n':>3} {'printed':>9} | {'H2a':>9} {'H2b':>9} | edge-cubic candidates (q3, q3'):")
    for N in sorted(PRINTED):
        lam, C = char_net(N)
        alpha = -np.log(lam) / np.log(0.5)
        # H2 candidates from lambda alone
        h2a = (alpha * (alpha - 1) * (alpha - 2)) / 6.0    # u^alpha cubic Taylor coeff
        h2b = (1 - lam) / (3 * lam) - 1.0 / 3.0
        cands = []
        for which in (1, 2, 3):
            _, rows, cols = edge_cubic_bezier(N, which)
            for bez in rows + cols:
                r = mu0_from_cubic(bez)
                if r is not None:
                    cands.append((which, round(r[1], 5), round(r[2], 5)))
        tgt = PRINTED[N]
        print(f"{N:>3} {tgt:+9.5f} | {h2a:+9.5f} {h2b:+9.5f} | " +
              " ".join(f"[{w}]{a:+.5f}/{b:+.5f}" for w, a, b in cands[:6]))


if __name__ == "__main__":
    main()
