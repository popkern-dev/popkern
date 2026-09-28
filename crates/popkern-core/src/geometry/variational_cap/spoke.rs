//! Exact G2 across the spokes: Loop & Schaefer 2008, eqs. (8)–(10).
//!
//! Between sector `i` and sector `i+1` the shared spoke is
//! `P_i(0,t) = P_{i+1}(t,0)`. With `c = cos(2π/n)`:
//!
//! ```text
//! (8)  P_i(0,t) = P_{i+1}(t,0)
//! (9)  c(1−t)(∂u P_{i+1}(t,0) + ∂v P_i(0,t)) = ∂v P_{i+1}(t,0) + ∂u P_i(0,t)
//! (10) 2c²(1−t)(∂v P_{i+1} − ∂u P_i)
//!        + 2c(1−t)(1 + c(1−t))(∂uv P_{i+1} − ∂uv P_i)
//!        = (1 + c(1−t))(∂vv P_{i+1} − ∂uu P_i)
//! ```
//!
//! (8) is an identity between control rows (same knot vector); (9) and
//! (10) are piecewise-polynomial identities in `t` of degree at most
//! `degree + 1` per knot span, imposed at `degree + 3` Gauss points per
//! span, which makes them exact. They are the same for every
//! coordinate. Verified symbolically for n = 3..8 (a cubic surface
//! pulled back through φ_n satisfies them to 1e-14), and by
//! `tests::pulled_back_polynomials_satisfy_the_spoke_rows`.

use super::space::{Functional, SectorSpace, gauss_legendre};

/// Every hard interior row, per coordinate.
pub(super) fn spoke_rows(space: &SectorSpace, c: f64) -> Vec<Functional> {
    let n = space.n;
    let m = space.m;
    let mut rows = Vec::new();
    let nodes = gauss_legendre(space.degree + 3);
    for i in 0..n {
        let next = (i + 1) % n;
        // (8): the shared control row.
        for k in 0..m {
            rows.push(vec![
                (space.index(i, 0, k), 1.0),
                (space.index(next, k, 0), -1.0),
            ]);
        }
        for &(a, b) in &space.spans() {
            for &(x, _) in &nodes {
                let t = a + (b - a) * x;
                let s = 1.0 - t;
                let bt = space.basis(t);
                let b0 = space.basis(0.0);
                // P_i(0,t): u = 0, v = t.   P_{i+1}(t,0): u = t, v = 0.
                let pi = |du, dv| space.functional_from(i, &b0, &bt, du, dv);
                let pn = |du, dv| space.functional_from(next, &bt, &b0, du, dv);
                // (9)
                let mut r9 = Vec::new();
                push(&mut r9, &pn(1, 0), c * s);
                push(&mut r9, &pi(0, 1), c * s);
                push(&mut r9, &pn(0, 1), -1.0);
                push(&mut r9, &pi(1, 0), -1.0);
                rows.push(normalise(r9));
                // (10)
                let k1 = 2.0 * c * c * s;
                let k2 = 2.0 * c * s * (1.0 + c * s);
                let k3 = 1.0 + c * s;
                let mut r10 = Vec::new();
                push(&mut r10, &pn(0, 1), k1);
                push(&mut r10, &pi(1, 0), -k1);
                push(&mut r10, &pn(1, 1), k2);
                push(&mut r10, &pi(1, 1), -k2);
                push(&mut r10, &pn(0, 2), -k3);
                push(&mut r10, &pi(2, 0), k3);
                rows.push(normalise(r10));
            }
        }
    }
    rows
}

fn push(row: &mut Functional, f: &Functional, scale: f64) {
    if scale != 0.0 {
        row.extend(f.iter().map(|&(i, w)| (i, w * scale)));
    }
}

/// Merge duplicate indices and scale the row to unit max weight.
fn normalise(mut row: Functional) -> Functional {
    row.sort_by_key(|&(i, _)| i);
    let mut merged: Functional = Vec::with_capacity(row.len());
    for (i, w) in row {
        match merged.last_mut() {
            Some(last) if last.0 == i => last.1 += w,
            _ => merged.push((i, w)),
        }
    }
    let scale = merged.iter().map(|&(_, w)| w.abs()).fold(0.0, f64::max);
    if scale > 0.0 {
        for e in &mut merged {
            e.1 /= scale;
        }
    }
    merged.retain(|&(_, w)| w != 0.0);
    merged
}
