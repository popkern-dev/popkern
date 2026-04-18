//! Exact rational-NURBS conversion for `RectangularGregory`.
//!
//! Converts a cubic rectangular rational Gregory patch to a bidegree-(7,7)
//! rational Bezier patch (expressed as a clamped `NurbsSurface`). The
//! conversion is exact by polynomial identity away from the four corners,
//! and exact *at* the corners via ε-regularization (the natural polynomial
//! denominator `D = D_0 · D_1 · D_2 · D_3` vanishes at corners, so the
//! corner weights would be zero — we force them to a small positive ε and
//! set the corresponding control point to the Gregory corner position so
//! the evaluation reproduces the analytic corner exactly).
//!
//! # Algorithm
//!
//! The rectangular Gregory patch evaluator is
//!
//! ```text
//!   S(u,v) = Σ_{(i,j) ∉ I} B_i^3(u) B_j^3(v) · P_{ij}
//!          + Σ_{(i,j) ∈ I} B_i^3(u) B_j^3(v) · C_{ij}(u,v)
//! ```
//!
//! with `I = {(1,1), (1,2), (2,1), (2,2)}` and rational inner control
//! points
//!
//! ```text
//!   C_{1,1}(u,v) = (v·P_{0,u} + u·P_{0,v}) / (u + v)          (= N_0 / D_0)
//!   C_{1,2}(u,v) = ((1-v)·P_{1,u} + u·P_{1,v}) / (u + 1 - v)  (= N_1 / D_1)
//!   C_{2,2}(u,v) = ((1-v)·P_{2,u} + (1-u)·P_{2,v}) / (2-u-v)  (= N_2 / D_2)
//!   C_{2,1}(u,v) = (v·P_{3,u} + (1-u)·P_{3,v}) / (1-u+v)      (= N_3 / D_3)
//! ```
//!
//! Multiplying through by `D = D_0·D_1·D_2·D_3` clears all four
//! denominators, giving a bidegree-(7,7) numerator `N(u,v)` and a
//! bidegree-(4,4) denominator `D(u,v)`. Elevating `D` to bidegree (7,7)
//! puts both on the same basis, and the Bernstein coefficients give the
//! weighted control points `w_{ij} Q_{ij}` and weights `w_{ij}` of a
//! rational Bezier patch that equals `S` pointwise.
//!
//! # References
//!
//! * Gregory, J. A. (1974). "Smooth Interpolation Without Twist Constraints."
//! * Chiyokura, H., & Kimura, F. (1983). "Design of Solids with Free-form Surfaces."
//! * Piegl, L., & Tiller, W. (1997). "The NURBS Book", §9.

use super::rectangular::RectangularGregory;
use crate::geometry::Point3;
use crate::geometry::nurbs_surface::NurbsSurface;

/// Bidegree-(m, n) bivariate Bernstein polynomial with scalar coefficients,
/// indexed `coefs[i * (n + 1) + j]` for term `B_i^m(u) B_j^n(v)`.
#[derive(Clone)]
pub(super) struct Bivar {
    pub m: usize,
    pub n: usize,
    pub coefs: Vec<f64>,
}

impl Bivar {
    pub(super) fn zeros(m: usize, n: usize) -> Self {
        Self {
            m,
            n,
            coefs: vec![0.0; (m + 1) * (n + 1)],
        }
    }

    #[inline]
    pub(super) fn get(&self, i: usize, j: usize) -> f64 {
        self.coefs[i * (self.n + 1) + j]
    }

    #[inline]
    fn set(&mut self, i: usize, j: usize, v: f64) {
        self.coefs[i * (self.n + 1) + j] = v;
    }

    #[inline]
    fn add(&mut self, i: usize, j: usize, v: f64) {
        self.coefs[i * (self.n + 1) + j] += v;
    }

    /// Bernstein-basis product.
    ///
    /// Uses the identity `B_i^m(t) · B_j^n(t) = [C(m,i) C(n,j) / C(m+n,i+j)] ·
    /// B_{i+j}^{m+n}(t)` in each variable; bivariate product factors.
    pub(super) fn multiply(&self, other: &Self) -> Self {
        let m_out = self.m + other.m;
        let n_out = self.n + other.n;
        let mut out = Bivar::zeros(m_out, n_out);
        for i in 0..=self.m {
            for j in 0..=self.n {
                let a = self.get(i, j);
                if a == 0.0 {
                    continue;
                }
                for ip in 0..=other.m {
                    for jp in 0..=other.n {
                        let b = other.get(ip, jp);
                        if b == 0.0 {
                            continue;
                        }
                        let k = i + ip;
                        let l = j + jp;
                        let coef = (binom(self.m, i) * binom(other.m, ip) / binom(m_out, k))
                            * (binom(self.n, j) * binom(other.n, jp) / binom(n_out, l));
                        out.add(k, l, coef * a * b);
                    }
                }
            }
        }
        out
    }

    /// Raise the bidegree of `self` to `(m_target, n_target)` via Bernstein
    /// degree elevation in each direction separately.
    pub(super) fn elevate_to(&self, m_target: usize, n_target: usize) -> Self {
        assert!(m_target >= self.m && n_target >= self.n);
        let mut cur = self.clone();
        while cur.m < m_target {
            cur = elevate_u(&cur);
        }
        while cur.n < n_target {
            cur = elevate_v(&cur);
        }
        cur
    }
}

/// Binomial coefficient `C(n, k)` for small `n`, returned as `f64` so it
/// composes into ratios directly. Exact for `n ≤ 34` in f64.
fn binom(n: usize, k: usize) -> f64 {
    if k > n {
        return 0.0;
    }
    let mut acc = 1.0;
    let k = k.min(n - k);
    for i in 0..k {
        acc = acc * ((n - i) as f64) / ((i + 1) as f64);
    }
    acc
}

/// Elevate u-degree by 1: `a_i B_i^m = (i/(m+1)) a_{i-1} B_i^{m+1} +
/// ((m+1-i)/(m+1)) a_i B_i^{m+1}`, summed.
fn elevate_u(src: &Bivar) -> Bivar {
    let m1 = src.m + 1;
    let mut out = Bivar::zeros(m1, src.n);
    for i in 0..=m1 {
        for j in 0..=src.n {
            let left = if i > 0 {
                (i as f64 / m1 as f64) * src.get(i - 1, j)
            } else {
                0.0
            };
            let right = if i <= src.m {
                ((m1 - i) as f64 / m1 as f64) * src.get(i, j)
            } else {
                0.0
            };
            out.set(i, j, left + right);
        }
    }
    out
}

/// Elevate v-degree by 1.
fn elevate_v(src: &Bivar) -> Bivar {
    let n1 = src.n + 1;
    let mut out = Bivar::zeros(src.m, n1);
    for i in 0..=src.m {
        for j in 0..=n1 {
            let left = if j > 0 {
                (j as f64 / n1 as f64) * src.get(i, j - 1)
            } else {
                0.0
            };
            let right = if j <= src.n {
                ((n1 - j) as f64 / n1 as f64) * src.get(i, j)
            } else {
                0.0
            };
            out.set(i, j, left + right);
        }
    }
    out
}

/// The four corner denominators `D_0..D_3` of the Gregory rational blend,
/// each as a bidegree-(1,1) Bernstein polynomial on `[0,1]²`.
///
/// Recall that a bidegree-(1,1) polynomial `p(u,v)` has Bernstein
/// coefficients equal to its vertex values: `c_{ij} = p(i, j)` for
/// `i, j ∈ {0, 1}`.
pub(super) fn denom_bernstein() -> [Bivar; 4] {
    // D_0 = u + v      → vertex values (0, 1, 1, 2)
    // D_1 = u + 1 - v  → vertex values (1, 0, 2, 1)  (remember j=0 is v=0)
    // D_2 = 2 - u - v  → vertex values (2, 1, 1, 0)
    // D_3 = 1 - u + v  → vertex values (1, 2, 0, 1)
    let mk = |v00: f64, v01: f64, v10: f64, v11: f64| {
        let mut b = Bivar::zeros(1, 1);
        b.set(0, 0, v00);
        b.set(0, 1, v01);
        b.set(1, 0, v10);
        b.set(1, 1, v11);
        b
    };
    [
        mk(0.0, 1.0, 1.0, 2.0),
        mk(1.0, 0.0, 2.0, 1.0),
        mk(2.0, 1.0, 1.0, 0.0),
        mk(1.0, 2.0, 0.0, 1.0),
    ]
}

/// Build a bidegree-(3,3) Bernstein polynomial from a single `B_i^3(u) B_j^3(v)`
/// term with unit coefficient (all other coefficients zero).
fn bernstein_basis_33(i: usize, j: usize) -> Bivar {
    let mut b = Bivar::zeros(3, 3);
    b.set(i, j, 1.0);
    b
}

/// Build a bidegree-(5,5) Bernstein polynomial from a single `B_i^5(u) B_j^5(v)`
/// term with unit coefficient.
pub(super) fn bernstein_basis_55(i: usize, j: usize) -> Bivar {
    let mut b = Bivar::zeros(5, 5);
    b.set(i, j, 1.0);
    b
}

/// Linear-coefficient helper: build a bidegree-(1,1) Bernstein polynomial
/// from its four vertex values.
pub(super) fn linear_bernstein(v00: f64, v01: f64, v10: f64, v11: f64) -> Bivar {
    let mut b = Bivar::zeros(1, 1);
    b.set(0, 0, v00);
    b.set(0, 1, v01);
    b.set(1, 0, v10);
    b.set(1, 1, v11);
    b
}

impl RectangularGregory {
    /// Convert the Gregory patch to an exact bidegree-(7,7) rational Bezier
    /// surface, represented as a clamped [`NurbsSurface`] on `[0,1]²`.
    ///
    /// The conversion is polynomial-identity exact on the open interior.
    /// At the four corners, where the natural denominator `D` vanishes,
    /// weights are ε-regularized and control points are pinned to the
    /// analytic Gregory corner positions so that evaluation at the
    /// corners returns the correct point. The resulting surface agrees
    /// with [`RectangularGregory::point_at`] to floating-point precision.
    pub fn to_nurbs_exact(&self) -> NurbsSurface {
        let denoms = denom_bernstein();

        // D(u,v) = D_0 · D_1 · D_2 · D_3 — bidegree (4,4).
        let d_full = denoms[0]
            .multiply(&denoms[1])
            .multiply(&denoms[2])
            .multiply(&denoms[3]);

        // For each (i,j) index in the 3×3 inner block we build the
        // bidegree-(7,7) numerator contribution. Start with all zero for
        // each of the three coordinate numerators.
        let mut n_x = Bivar::zeros(7, 7);
        let mut n_y = Bivar::zeros(7, 7);
        let mut n_z = Bivar::zeros(7, 7);

        // --- Non-inner contribution: B_i^3 B_j^3 · P_{ij} · D.
        for i in 0..4 {
            for j in 0..4 {
                if matches!((i, j), (1, 1) | (1, 2) | (2, 1) | (2, 2)) {
                    continue;
                }
                let p = self.boundary_cps[i][j];
                // (B_i^3 B_j^3) × D  is bidegree (3+4, 3+4) = (7,7).
                let term = bernstein_basis_33(i, j).multiply(&d_full);
                accumulate_scaled(&mut n_x, &term, p.x);
                accumulate_scaled(&mut n_y, &term, p.y);
                accumulate_scaled(&mut n_z, &term, p.z);
            }
        }

        // --- Inner contributions: B_{i_k}^3 B_{j_k}^3 · (w_v P_{k,u} + w_u P_{k,v}) · (D/D_k).
        // D/D_k = product of the other three denominators (bidegree (3,3)).
        let d_over_dk: [Bivar; 4] = [
            denoms[1].multiply(&denoms[2]).multiply(&denoms[3]),
            denoms[0].multiply(&denoms[2]).multiply(&denoms[3]),
            denoms[0].multiply(&denoms[1]).multiply(&denoms[3]),
            denoms[0].multiply(&denoms[1]).multiply(&denoms[2]),
        ];

        // Weight-pair coefficient tables per corner, bidegree (1,1).
        //
        // Corner 0 inner (1,1): w_v = v, w_u = u.
        //   numerator_coord = v·P_0u + u·P_0v
        //   v in (1,1)-Bernstein has coefs [0,1,0,1] (for j=0→0, j=1→1 at each i)
        //   u in (1,1)-Bernstein has coefs [0,0,1,1]
        //
        // Corner 1 inner (1,2): w_v = 1-v, w_u = u. (from (1,2)=>(1,u,1-v) where wu=u, wv=1-v)
        //   numerator = (1-v)·P_1u + u·P_1v
        //
        // Corner 2 inner (2,2): w_v = 1-v, w_u = 1-u.
        //   numerator = (1-v)·P_2u + (1-u)·P_2v
        //
        // Corner 3 inner (2,1): w_v = v, w_u = 1-u.
        //   numerator = v·P_3u + (1-u)·P_3v
        let corner_data = [
            // (inner index i,j, k, w_v_bernstein, w_u_bernstein)
            (
                1usize,
                1usize,
                0usize,
                linear_bernstein(0.0, 1.0, 0.0, 1.0), // v
                linear_bernstein(0.0, 0.0, 1.0, 1.0), // u
            ),
            (
                1,
                2,
                1,
                linear_bernstein(1.0, 0.0, 1.0, 0.0), // 1 - v
                linear_bernstein(0.0, 0.0, 1.0, 1.0), // u
            ),
            (
                2,
                2,
                2,
                linear_bernstein(1.0, 0.0, 1.0, 0.0), // 1 - v
                linear_bernstein(1.0, 1.0, 0.0, 0.0), // 1 - u
            ),
            (
                2,
                1,
                3,
                linear_bernstein(0.0, 1.0, 0.0, 1.0), // v
                linear_bernstein(1.0, 1.0, 0.0, 0.0), // 1 - u
            ),
        ];

        for (i, j, k, wv, wu) in &corner_data {
            let (p_u, p_v) = self.inner_cps[*k];
            // base = B_i^3 B_j^3 · (D/D_k)  — bidegree (3+3, 3+3) = (6,6).
            let base = bernstein_basis_33(*i, *j).multiply(&d_over_dk[*k]);
            // w_v · base — bidegree (6+1, 6+1) = (7,7).
            let term_u = wv.multiply(&base);
            // w_u · base — bidegree (7,7).
            let term_v = wu.multiply(&base);
            accumulate_scaled(&mut n_x, &term_u, p_u.x);
            accumulate_scaled(&mut n_y, &term_u, p_u.y);
            accumulate_scaled(&mut n_z, &term_u, p_u.z);
            accumulate_scaled(&mut n_x, &term_v, p_v.x);
            accumulate_scaled(&mut n_y, &term_v, p_v.y);
            accumulate_scaled(&mut n_z, &term_v, p_v.z);
        }

        // Elevate D from (4,4) to (7,7) so the weight polynomial is in the
        // same Bernstein basis as N.
        let d_elev = d_full.elevate_to(7, 7);

        // Assemble the 8×8 weighted control points and weights.
        let mut cps: Vec<Point3> = Vec::with_capacity(64);
        let mut weights: Vec<f64> = Vec::with_capacity(64);
        for i_u in 0..8 {
            for i_v in 0..8 {
                let w = d_elev.get(i_u, i_v);
                let wx = n_x.get(i_u, i_v);
                let wy = n_y.get(i_u, i_v);
                let wz = n_z.get(i_u, i_v);
                cps.push(Point3::new(wx, wy, wz)); // store weighted CP for now
                weights.push(w);
            }
        }

        // ε-regularize the four corner grid entries. D vanishes exactly at
        // (u,v) = (0,0), (0,1), (1,0), (1,1) — those map to grid index
        // (0,0), (0,7), (7,0), (7,7). We set the weight to EPS and the
        // weighted CP to EPS · P_corner so that Q = P_corner.
        let epsilon = 1e-14;
        let corners = [
            (0usize, 0usize, self.boundary_cps[0][0]),
            (0, 7, self.boundary_cps[0][3]),
            (7, 0, self.boundary_cps[3][0]),
            (7, 7, self.boundary_cps[3][3]),
        ];
        for (i_u, i_v, p) in corners {
            let idx = i_u * 8 + i_v;
            weights[idx] = epsilon;
            cps[idx] = Point3::new(epsilon * p.x, epsilon * p.y, epsilon * p.z);
        }

        // Divide weighted CPs by weights to get final Q_{ij}.
        let control_points: Vec<Point3> = cps
            .iter()
            .zip(weights.iter())
            .map(|(wcp, &w)| Point3::new(wcp.x / w, wcp.y / w, wcp.z / w))
            .collect();

        // Clamped degree-7 knot vectors in each direction (pure Bezier).
        let mut knots = Vec::with_capacity(16);
        knots.extend(std::iter::repeat_n(0.0, 8));
        knots.extend(std::iter::repeat_n(1.0, 8));

        NurbsSurface {
            degree_u: 7,
            degree_v: 7,
            knots_u: knots.clone(),
            knots_v: knots,
            control_points,
            num_u: 8,
            num_v: 8,
            weights,
        }
    }
}

/// `out += scale · src`. Bidegrees must match.
pub(super) fn accumulate_scaled(out: &mut Bivar, src: &Bivar, scale: f64) {
    assert_eq!(out.m, src.m);
    assert_eq!(out.n, src.n);
    for k in 0..out.coefs.len() {
        out.coefs[k] += scale * src.coefs[k];
    }
}
