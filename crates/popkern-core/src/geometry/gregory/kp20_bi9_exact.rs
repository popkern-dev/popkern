//! The exact composition `h = g ∘ τ` of one KP20 cap piece and its
//! derivatives through second order, in the piece's own `(u, v)`
//! parameters (PK-599).
//!
//! `τ` is the planar bi-6 piece (rotated into its sector); the sector
//! guide `g` is a total-degree-6 polynomial in its own chart. The chart
//! map `L` (global plane → guide `(a, b)`, spoke `k` along `a`, spoke
//! `k+1` along `b`) and the sector rotation `R` are linear, so with
//! `T = L R τ`:
//!
//! ```text
//! h_u  = g_a T¹_u + g_b T²_u
//! h_uv = g_aa T¹_u T¹_v + g_ab (T¹_u T²_v + T²_u T¹_v) + g_bb T²_u T²_v
//!        + g_a T¹_uv + g_b T²_uv
//! ```
//!
//! and likewise for the other orders — closed form, no differencing.

// Coefficient nets are addressed by (i, j) index pairs throughout.
#![allow(clippy::needless_range_loop)]

use super::kp20_bi6_cap::Kp20Bi6Patch;
use super::kp20_bi6_guide_sample::Kp20Tri6GuideSector;
use crate::geometry::{Point3, Vec3};
use popkern_math::PortableMath;

/// Value and derivatives through order two at one parameter point.
#[derive(Debug, Clone, Copy)]
pub(super) struct Jet2 {
    pub p: Vec3,
    pub u: Vec3,
    pub v: Vec3,
    pub uu: Vec3,
    pub vv: Vec3,
}

/// One piece's exact composition `h = g ∘ τ`.
pub(super) struct ExactPiece<'a> {
    guide: &'a Kp20Tri6GuideSector,
    /// Planar τ coefficients `(x, y)`.
    tau: [[[f64; 2]; 7]; 7],
    /// `L R` as a 2×2 matrix (row-major).
    map: [[f64; 2]; 2],
}

fn bernstein(n: usize, i: usize, t: f64) -> f64 {
    if i > n {
        return 0.0;
    }
    let c = (0..i).fold(1.0, |c, k| c * (n - k) as f64 / (k + 1) as f64);
    c * t.powi(i as i32) * (1.0 - t).powi((n - i) as i32)
}

/// Degree-6 Bernstein values and first/second derivatives at `t`.
fn basis6(t: f64) -> [[f64; 7]; 3] {
    let mut out = [[0.0; 7]; 3];
    for i in 0..7 {
        out[0][i] = bernstein(6, i, t);
        let b5 = |k: isize| {
            if k < 0 {
                0.0
            } else {
                bernstein(5, k as usize, t)
            }
        };
        let b4 = |k: isize| {
            if k < 0 {
                0.0
            } else {
                bernstein(4, k as usize, t)
            }
        };
        let k = i as isize;
        out[1][i] = 6.0 * (b5(k - 1) - b5(k));
        out[2][i] = 30.0 * (b4(k - 2) - 2.0 * b4(k - 1) + b4(k));
    }
    out
}

impl<'a> ExactPiece<'a> {
    /// The composition for piece `tau_piece` of sector `sector` of an
    /// `n`-sided cap.
    pub fn new(
        n: usize,
        sector: usize,
        guide: &'a Kp20Tri6GuideSector,
        tau_piece: &Kp20Bi6Patch,
    ) -> Self {
        let alpha = core::f64::consts::TAU / n as f64;
        let (ca, sa) = (
            (sector as f64 * alpha).pcos(),
            (sector as f64 * alpha).psin(),
        );
        let (a0, a1) = (sector as f64 * alpha, (sector + 1) as f64 * alpha);
        let (s0, s1) = ([a0.pcos(), a0.psin()], [a1.pcos(), a1.psin()]);
        let det = s1[0] * s0[1] - s1[1] * s0[0];
        // L: [a; b] = [[-s1y, s1x], [s0y, -s0x]] / det · [x; y].
        let l = [[-s1[1] / det, s1[0] / det], [s0[1] / det, -s0[0] / det]];
        let r = [[ca, -sa], [sa, ca]];
        let map = [
            [
                l[0][0] * r[0][0] + l[0][1] * r[1][0],
                l[0][0] * r[0][1] + l[0][1] * r[1][1],
            ],
            [
                l[1][0] * r[0][0] + l[1][1] * r[1][0],
                l[1][0] * r[0][1] + l[1][1] * r[1][1],
            ],
        ];
        let tau = core::array::from_fn(|i| {
            core::array::from_fn(|j| [tau_piece.cps[i][j].x, tau_piece.cps[i][j].y])
        });
        Self { guide, tau, map }
    }

    /// `h` and its derivatives through order two at `(u, v)`.
    pub fn jet(&self, u: f64, v: f64) -> Jet2 {
        let (bu, bv) = (basis6(u), basis6(v));
        // τ partials [value, u, v, uu, uv, vv] × (x, y).
        let orders = [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2)];
        let mut t = [[0.0f64; 2]; 6];
        for (k, &(du, dv)) in orders.iter().enumerate() {
            for i in 0..7 {
                for j in 0..7 {
                    let w = bu[du][i] * bv[dv][j];
                    t[k][0] += w * self.tau[i][j][0];
                    t[k][1] += w * self.tau[i][j][1];
                }
            }
        }
        let m = |p: [f64; 2]| {
            [
                self.map[0][0] * p[0] + self.map[0][1] * p[1],
                self.map[1][0] * p[0] + self.map[1][1] * p[1],
            ]
        };
        let tt: [[f64; 2]; 6] = core::array::from_fn(|k| m(t[k]));
        let g = self.guide.jet3(tt[0][0], tt[0][1]).h;
        let gv = |i: usize, j: usize| Vec3::new(g[i][j].x, g[i][j].y, g[i][j].z);
        let (ga, gb, gaa, gab, gbb) = (gv(1, 0), gv(0, 1), gv(2, 0), gv(1, 1), gv(0, 2));
        let first = |d: [f64; 2]| ga * d[0] + gb * d[1];
        let second = |d: [f64; 2], e: [f64; 2], de: [f64; 2]| {
            gaa * (d[0] * e[0])
                + gab * (d[0] * e[1] + d[1] * e[0])
                + gbb * (d[1] * e[1])
                + first(de)
        };
        let (tu, tv) = (tt[1], tt[2]);
        Jet2 {
            p: gv(0, 0),
            u: first(tu),
            v: first(tv),
            uu: second(tu, tu, tt[3]),
            vv: second(tv, tv, tt[5]),
        }
    }

    /// The exact point `h(u, v)`.
    pub fn point(&self, u: f64, v: f64) -> Point3 {
        let p = self.jet(u, v).p;
        Point3::new(p.x, p.y, p.z)
    }
}
