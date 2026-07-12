# KP15 bi-6 G² cap — Maple-script reference oracle (PK-087)

Independent sympy oracle for the Karčiauskas & Peters 2015 bi-6 polar
cap construction — "*Minimal bi-6 G² completion of bicubic spline
surfaces*". It is a faithful, end-to-end port of the paper's companion
Maple script (`docs/papers/n6.txt`, the n = 6 worked example; that
directory is gitignored and never committed — source URL:
<https://www.cise.ufl.edu/research/SurfLab/papers/n6.txt>).

The script symbolically resolves **all** G1/G2 constraints, assigning
every dependent bi-6 BB-coefficient of the cap. That closed-form
resolution was the long-standing blocker (PK-087 / see
`docs/journal/` KP15 notes): the Rust implementation left ~16
deep-interior CPs per sector unpopulated (star-point defaults), which
geometrically guaranteed folded caps. This oracle recovers the full
resolution and cross-checks the Rust tables against it.

## Indexing flip (critical)

The script uses a flipped BB index relative to the paper:

    script  p[i,j]  ==  paper  p[6-i, 6-j]
    script  [6,6]   ==  the extraordinary point (paper [0,0])

**Every CP index in this directory is in SCRIPT convention.** The free
set (`Pfree`) is: `P00`, `TAN[0..1]`, `TWO[0..2]` (central quadratic),
`P220 = p[0][4,4]`, `P210 = p[0][4,5]`, and per-sector
`P23[r] = p[r][4,3]`, `P33[r] = p[r][3,3]`.

## What it verifies

`n6_oracle.py` reproduces the script step for step: `TBdata` (tensor
border from bi-3), `BtoBB`, (S1) central quadratic + eq (6) C² sector
glue, (S2) the eq (7) `2n×2n` circulant (solves `p[r][5,4]`,
`p[r][4,5]`) and the eq (8) `n×n` circulant (solves `p[r][4,4]`;
`lyc[0] = sal` deferred), `d2(τ)` from eq (12,13), the reduced-degree
sector curve eq (5), row-6 glue, and (S3) the `µ̃`/`µ̄` rules for
`p[r][5,3]`, `p[r+1][3,5]`, `p[r+1][3,4]`. It then re-runs the
script's own smoothness verification per sector:

    G1[s] = Fv + Gv - 2c(1-u)² Fu
    G2[s] = -Gvv + Fvv - 4c(1-u)² Fuv + 4c²(1-u)⁴ Fuu
            + c(1-u)³ (A u - 8c - A) Fu + A(1-u)² Fv,
      A = 24c(τ-1)/(1-u+τu),   evaluated at v = 0

plus the deferred `sal` residual. Acceptance: all residuals identically
zero.

## How to run

The venv is **not** in the repo. Create it in scratch:

```bash
python3 -m venv <scratch>/kp15venv
<scratch>/kp15venv/bin/pip install sympy numpy
PY=<scratch>/kp15venv/bin/python

cd tools/reference/kp15-bi6-cap
$PY n6_oracle.py --mode symbolic          # symbolic acceptance (zero net)
$PY n6_oracle.py --mode symbolic --net shared   # + symbolic C²-consistent net
$PY n6_oracle.py --mode numeric --n 6     # exact-rational + float residuals
$PY n6_oracle.py --mode valence           # n = 5..10 generalization probe
$PY n6_oracle.py --mode all               # everything
$PY rust_table_xcheck.py                  # diff Rust µ/γ/d2 tables vs script
$PY gen_rust_fixtures.py                  # write n6_cap_fixture.json
```

Modes:
- **symbolic** — free values kept symbolic, `c = 1/2`, `τ = 17/20`
  exact-rational. Proves `G1 ≡ G2 ≡ sal ≡ 0` as identities. With
  `--net shared` the C²-consistent input net (`a1..a6`, `E`) is also
  symbolic, so the `TBdata` border-reparameterisation path is proven
  symbolically too.
- **numeric** — random rational free values + C²-consistent shared net.
  At n = 6 the whole thing is exact-rational, so residuals are exactly
  0 (bit-for-bit). Float paths report sampled maxima over
  `u ∈ {0.1,…,0.9}`.
- **valence** — float construction at n ∈ {5,6,7,8,9,10} with the
  paper τ table, shared net.

## Findings

**Symbolic acceptance (n = 6):** `G1[s] ≡ 0`, `G2[s] ≡ 0` for all 6
sectors and `sal ≡ 0`, both with zero net and with a fully symbolic
C²-consistent net. `d2 = 103/1156` at τ = 17/20. The free set above is
confirmed exactly.

**Numeric residual maxima (float construction + float sampling):**

| n  | τ    | max\|G1\| | max\|G2\| | \|sal\| | verdict |
|----|------|-----------|-----------|---------|---------|
| 5  | 0.87 | 2.2e-13   | 1.8e-12   | 5.2e-14 | CLOSES  |
| 6  | 0.85 | 1.5e-13   | 2.3e-12   | 2.0e-13 | CLOSES  |
| 7  | 0.84 | 3.1e-14   | 2.9e-13   | 1.1e-14 | CLOSES  |
| 8  | 0.83 | 6.0e-14   | 5.8e-13   | 4.1e-14 | CLOSES  |
| 9  | 0.82 | 4.7e-13   | 4.5e-12   | 3.2e-14 | CLOSES  |
| 10 | 0.81 | 5.7e-14   | 6.3e-12   | 9.5e-15 | CLOSES  |

At n = 6 the exact-rational run gives residuals **identically 0**.

**Valence-generalization verdict: the closed forms GENERALIZE.** Every
weight in the script is a rational function of `c = cos(2π/n)` and `τ`;
running the identical construction at n = 5..10 closes G1/G2 to
~1e-12 (float noise). n = 6 is **not** special. (n = 3 is the paper
§3.5 special case and is skipped.)

**Free-set is valence-dependent (matches Rust `needs_pin`).** The eq (7)
`2n×2n` circulant is rank-deficient (rank 2n-1) **only** at n ∈ {3, 6};
there `p[0][4,5] = P210` is a genuine free parameter. At every other n
the circulant is full rank and `p[0][4,5]` is *determined* (P210
unused). This matches `needs_pin = n == 3 || n == 6` in
`kp15_bi6_circulant.rs`. `p[0][4,4] = P220` is free at all valences.

### Rust table cross-check (`rust_table_xcheck.py`, n = 6, τ = 0.85)

All coefficient vectors in `kp15_bi6_mu.rs` match the script **exactly**
(diff = 0):

- `mu_dot` (µ̇ / mu1), `mu_ddot` (µ̈ / mu2), `mu_tdot` (⃛µ / mu3),
  `mu_tilde` (µ̃ / mut), `mu_bar` (µ̄ / mub) — all entries match,
  including the reconstructed sum-entry `[0]`.
- `d2` in `kp15_bi6_reparam.rs` matches script eq (13).
- `w_scalars` (eq 11) is **not** in `n6.txt` (the partial script stops
  at S3), so it is not cross-checkable here.

**γ table discrepancy — `kp15_bi6_reparam.rs::gamma`, row k = 2.** The
γ table reproduces `TBdata` row-2 exactly for k = 3,4,5,6, but the
k = 2 row disagrees with the script:

| b entry | script γ²ᵢⱼ | Rust gamma(2,i,j) |
|---------|-------------|-------------------|
| b₀₀     | `(1-d2)/25` | `1 - d2/25`  ← wrong (missing `/25` on the `1`) |
| b₀₁     | `3/25 + d2/50` | `0`  ← missing entry `(2,0,1)` |
| b₀₂     | `1/25`      | `0`  ← missing entry `(2,0,2)` |
| b₁₂     | `3/25`      | `0`  ← missing entry `(2,1,2)` |

All other k = 2 entries (b₁₀,b₁₁,b₂₀,b₂₁,b₂₂) match. **Impact: latent.**
`gamma()` has **no call sites** outside its own unit tests — the KP15
polar reparam path is currently unreachable scaffolding (consistent
with the journal note that the KP15 end-to-end polar cap is blocked/
scaffold-only), so this bug does not presently reach production caps.
It must be fixed before that path is wired in.

### Other Rust contradictions found

Beyond the γ k = 2 row, the Rust `µ`/`d2` formulas and the free-vs-
determined classification agree with the script. The (S1), (S2)
circulant, and (S3) *structure* in `kp15_bi6_circulant.rs` /
`kp15_bi6_rule_s3.rs` were read and their coefficient wiring matches
the script's rules; the authoritative check that they compose to
G²-closed caps is the oracle's `G1 ≡ G2 ≡ 0` result, which the Rust
pipeline should be validated against once the polar path is live. The
`beta_d_bb_coeffs` / `beta_e_bb_coeffs` (eq 12) are not line-item
present in `n6.txt` (only `d2` from eq 13 is), so they are validated
only indirectly via the `TBdata`-embedded reparam closing G1/G2.

## Fixture

`gen_rust_fixtures.py` writes `n6_cap_fixture.json` (committed): a
fixed-seed rational instance with all 16 deep-block CPs (script indices
i,j ∈ [3,6]) per sector, the input net, the free-value choices, and the
(zero) residual provenance. Future Rust tests consume this.
