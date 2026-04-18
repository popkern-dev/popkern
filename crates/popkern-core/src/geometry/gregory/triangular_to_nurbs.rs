//! Exact rational-NURBS conversion for [`TriangularGregory`].
//!
//! A triangular Gregory patch on the barycentric simplex `{(u, v, w) :
//! u + v + w = 1, u, v, w ≥ 0}` is re-parameterised onto the
//! rectangular domain `(s, t) ∈ [0, 1]²` via
//!
//! ```text
//!   u = s,  v = (1 − s)(1 − t),  w = (1 − s)·t
//! ```
//!
//! which collapses the `s = 1` edge of the rectangle onto triangular
//! corner 0 (a **degenerate Bezier edge**).
//!
//! # Key identity
//!
//! Under the substitution above, the triangular Bernstein basis
//! reduces to a product of univariate Bernsteins:
//!
//! ```text
//!   B_{ijk}^n(u, v, w) = B_i^n(s) · B_k^{n−i}(t),    with i + j + k = n.
//! ```
//!
//! Proof: `B_{ijk}^n = C(n, i, j, k) · u^i v^j w^k
//!                   = C(n, i, j, k) · s^i (1−s)^{j+k} (1−t)^j t^k`,
//! and the trinomial identity
//! `C(n, i, j, k) = C(n, i) · C(n − i, j)` combined with the definition
//! `B_m^p(x) = C(p, m) · x^m (1 − x)^(p−m)` cancels the binomials.
//!
//! # Numerator / denominator
//!
//! Longhi's triangular Gregory substitutes the centre CP `P_{111}`
//! with a barycentric-rational blend:
//!
//! ```text
//!   C(u, v, w) = (v·w · P_A + u·w · P_B + u·v · P_C) / (u·v + v·w + w·u),
//! ```
//!
//! where `P_A, P_B, P_C` are the three inner candidates. Multiplying
//! the patch through by `D = u·v + v·w + w·u` clears the denominator,
//! giving a numerator of total degree 5 and denominator of total
//! degree 2. Under the `(s, t)` substitution, both become bidegree
//! `(5, 5)` (numerator) and `(2, 2)` (denominator) `(s, t)`-
//! Bernstein polynomials. Elevating the denominator to `(5, 5)`
//! puts them in the same basis — the resulting rational-Bezier CPs
//! are the weighted and weight arrays of a bidegree-`(5, 5)`
//! rational NURBS.
//!
//! # ε-regularization
//!
//! `D(s, t) = 0` along the full `s = 1` edge (the collapse) and at
//! the two other corners `(0, 0)` and `(0, 1)` (simplex corners 1
//! and 2). At those points both numerator and denominator vanish;
//! we override them with `ε = 1e-14` and pin the weighted CP to the
//! appropriate triangular corner so the rational evaluator returns
//! the correct point. Error from the ε-regularization is `O(ε · M)`
//! where `M` is the surface derivative magnitude — below f64
//! precision in practice.
//!
//! Along `s = 1` every CP collapses to corner 0: the rational
//! Bezier value there is `Σ ε · P_0 · B_j^5(t) / Σ ε · B_j^5(t) = P_0`
//! — a constant curve, as required.

use super::rectangular_to_nurbs::{Bivar, accumulate_scaled};
use super::triangular::TriangularGregory;
use crate::geometry::Point3;
use crate::geometry::nurbs_surface::NurbsSurface;

impl TriangularGregory {
    /// Convert the triangular Gregory patch to an exact bidegree-
    /// `(5, 5)` rational Bezier on `(s, t) ∈ [0, 1]²` via the
    /// `(s, t) ↦ (s, (1 − s)(1 − t), (1 − s)·t)` parameter map.
    /// The `s = 1` edge of the rectangle collapses onto triangular
    /// corner 0.
    ///
    /// # Parameter mapping
    ///
    /// The NURBS surface evaluates at `(s, t)` such that:
    /// - `(s = 0, t)` traces the simplex edge from corner 1
    ///   (at `t = 0`) to corner 2 (at `t = 1`).
    /// - `(s = 1, t)` is the degenerate edge, pinned to corner 0.
    /// - `(s, t = 0)` traces the simplex edge from corner 1 to
    ///   corner 0.
    /// - `(s, t = 1)` traces the simplex edge from corner 2 to
    ///   corner 0.
    pub fn to_nurbs_exact(&self) -> NurbsSurface {
        // --- Step 1: build univariate Bernstein basis `Bivar`s.
        //
        // u, v, w as bidegree-(1, 1) `Bivar`s over (s, t):
        //   u = s                   → values at (s,t) lattice {(0,0), (0,1), (1,0), (1,1)}
        //                               = (0, 0, 1, 1)
        //   v = (1 − s)(1 − t)      → (1, 0, 0, 0)
        //   w = (1 − s) · t         → (0, 1, 0, 0)
        let u = lin_bivar(0.0, 0.0, 1.0, 1.0);
        let v = lin_bivar(1.0, 0.0, 0.0, 0.0);
        let w = lin_bivar(0.0, 1.0, 0.0, 0.0);

        // --- Step 2: build denominator D = u·v + v·w + w·u, bidegree (2, 2).
        let uv = u.multiply(&v);
        let vw = v.multiply(&w);
        let wu = w.multiply(&u);
        let mut d_raw = Bivar::zeros(2, 2);
        for k in 0..d_raw.coefs.len() {
            d_raw.coefs[k] = uv.coefs[k] + vw.coefs[k] + wu.coefs[k];
        }

        // --- Step 3: build the numerator N(s, t) of the Gregory patch
        // AFTER clearing D. Collect the contributions term by term.
        //
        // N = Σ_{(i,j,k) ≠ (1,1,1), i+j+k=3} B_{ijk}^3 · P_{ijk} · D
        //   + 6·u·v·w · (v·w · P_A + u·w · P_B + u·v · P_C)
        //
        // B_{ijk}^3 in (s, t) is B_i^3(s) · B_k^{3-i}(t), elevated to
        // bidegree (5, 5). `d_elev` brings D to bidegree (5, 5).
        let d_elev = d_raw.elevate_to(5, 5);

        let mut n_x = Bivar::zeros(5, 5);
        let mut n_y = Bivar::zeros(5, 5);
        let mut n_z = Bivar::zeros(5, 5);

        // 9 non-centre boundary CPs.
        let boundary_indices = [
            (3, 0, 0, 0usize), // corner 0
            (0, 3, 0, 1),      // corner 1
            (0, 0, 3, 2),      // corner 2
            (2, 1, 0, 3),      // edge 01 near corner 0
            (1, 2, 0, 4),      // edge 01 near corner 1
            (0, 2, 1, 5),      // edge 12 near corner 1
            (0, 1, 2, 6),      // edge 12 near corner 2
            (1, 0, 2, 7),      // edge 20 near corner 2
            (2, 0, 1, 8),      // edge 20 near corner 0
        ];
        for (i, _j, k, slot) in boundary_indices {
            // B_{ijk}^3 = B_i^3(s) · B_k^{3-i}(t), bidegree (3, 3-i).
            let basis = single_bern_st(i, k, 3, 3 - i);
            // Multiply by D (bidegree (2, 2)) → bidegree (5, 5-i).
            let term_raw = basis.multiply(&d_raw);
            let term = term_raw.elevate_to(5, 5);
            let p = self.boundary_cps[slot];
            accumulate_scaled(&mut n_x, &term, p.x);
            accumulate_scaled(&mut n_y, &term, p.y);
            accumulate_scaled(&mut n_z, &term, p.z);
        }

        // Inner contributions. `B_{111}^3 · C · D = 6·u·v·w ·
        // (v·w · P_A + u·w · P_B + u·v · P_C)` since the `D`
        // denominator in `C` cancels exactly against the multiplied
        // `D`. Expanding:
        //
        //   6·u·v²·w² · P_A + 6·u²·v·w² · P_B + 6·u²·v²·w · P_C.
        //
        // Each term is total degree 5 in (u, v, w) and bidegree
        // (5, 5) in (s, t) via `Bivar::multiply`.
        let u_v2_w2 = u.multiply(&v).multiply(&v).multiply(&w).multiply(&w);
        let u2_v_w2 = u.multiply(&u).multiply(&v).multiply(&w).multiply(&w);
        let u2_v2_w = u.multiply(&u).multiply(&v).multiply(&v).multiply(&w);

        let (p_a, p_b, p_c) = (self.inner_cps[0], self.inner_cps[1], self.inner_cps[2]);
        accumulate_scaled(&mut n_x, &u_v2_w2, 6.0 * p_a.x);
        accumulate_scaled(&mut n_y, &u_v2_w2, 6.0 * p_a.y);
        accumulate_scaled(&mut n_z, &u_v2_w2, 6.0 * p_a.z);
        accumulate_scaled(&mut n_x, &u2_v_w2, 6.0 * p_b.x);
        accumulate_scaled(&mut n_y, &u2_v_w2, 6.0 * p_b.y);
        accumulate_scaled(&mut n_z, &u2_v_w2, 6.0 * p_b.z);
        accumulate_scaled(&mut n_x, &u2_v2_w, 6.0 * p_c.x);
        accumulate_scaled(&mut n_y, &u2_v2_w, 6.0 * p_c.y);
        accumulate_scaled(&mut n_z, &u2_v2_w, 6.0 * p_c.z);

        // --- Step 4: assemble NURBS. 6×6 = 36 CPs.
        let mut cps: Vec<Point3> = Vec::with_capacity(36);
        let mut weights: Vec<f64> = Vec::with_capacity(36);
        for i_s in 0..6 {
            for i_t in 0..6 {
                weights.push(d_elev.get(i_s, i_t));
                cps.push(Point3::new(
                    n_x.get(i_s, i_t),
                    n_y.get(i_s, i_t),
                    n_z.get(i_s, i_t),
                ));
            }
        }

        // --- Step 5: ε-regularize the degenerate edge + other corners.
        // `D = 0` along the entire `s = 1` edge (6 CPs at i_s = 5,
        // i_t = 0..5). All pin to triangular corner 0.
        let epsilon = 1e-14;
        let c0 = self.boundary_cps[0];
        for i_t in 0..6 {
            let idx = 5 * 6 + i_t;
            weights[idx] = epsilon;
            cps[idx] = Point3::new(epsilon * c0.x, epsilon * c0.y, epsilon * c0.z);
        }
        // `D = 0` at the two remaining rectangular corners.
        // (0, 0) → triangular corner 1, (0, 5) → triangular corner 2.
        let c1 = self.boundary_cps[1];
        let c2 = self.boundary_cps[2];
        let idx_00 = 0;
        let idx_05 = 5;
        weights[idx_00] = epsilon;
        cps[idx_00] = Point3::new(epsilon * c1.x, epsilon * c1.y, epsilon * c1.z);
        weights[idx_05] = epsilon;
        cps[idx_05] = Point3::new(epsilon * c2.x, epsilon * c2.y, epsilon * c2.z);

        // Divide weighted CPs by weights to get final Q_{ij}.
        let control_points: Vec<Point3> = cps
            .iter()
            .zip(weights.iter())
            .map(|(wcp, &w_val)| Point3::new(wcp.x / w_val, wcp.y / w_val, wcp.z / w_val))
            .collect();

        let mut knots = Vec::with_capacity(12);
        knots.extend(std::iter::repeat_n(0.0, 6));
        knots.extend(std::iter::repeat_n(1.0, 6));

        NurbsSurface {
            degree_u: 5,
            degree_v: 5,
            knots_u: knots.clone(),
            knots_v: knots,
            control_points,
            num_u: 6,
            num_v: 6,
            weights,
        }
    }
}

/// Build a bidegree-(1, 1) `Bivar` from the four lattice-point
/// values `(v00, v01, v10, v11)` where `v_{i_s i_t}` is the value
/// at `(s, t) = (i_s, i_t)`.
fn lin_bivar(v00: f64, v01: f64, v10: f64, v11: f64) -> Bivar {
    let mut b = Bivar::zeros(1, 1);
    b.set_val(0, 0, v00);
    b.set_val(0, 1, v01);
    b.set_val(1, 0, v10);
    b.set_val(1, 1, v11);
    b
}

/// Build the bidegree-(m, n) `Bivar` `B_i^m(s) · B_k^n(t)` with all
/// other coefficients zero. Used to construct triangular-Bernstein
/// basis terms after the (s, t) substitution.
fn single_bern_st(i: usize, k: usize, m: usize, n: usize) -> Bivar {
    let mut b = Bivar::zeros(m, n);
    b.set_val(i, k, 1.0);
    b
}

/// Extension trait on `Bivar` to set coefficients. `Bivar::set` is
/// private to the rectangular module; we get similar access via a
/// wrapper helper that uses `pub(super)` internals.
trait BivarSetVal {
    fn set_val(&mut self, i: usize, j: usize, v: f64);
}

impl BivarSetVal for Bivar {
    fn set_val(&mut self, i: usize, j: usize, v: f64) {
        self.coefs[i * (self.n + 1) + j] = v;
    }
}
