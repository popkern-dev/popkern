//! The rim approximant: each side half as a B-spline control row in
//! the shared knot vector.
//!
//! Side `s` runs `C(t0)` (corner `c_{s-1}`) → `C(tm)` (`m_s`) → `C(t1)`
//! (corner `c_s`) with `tm` the parameter midpoint. Its trailing half is
//! sector `s`'s row `i = m−1` (`S_s(1, v)`, `v: m_s → c_s`); its leading
//! half is sector `s−1`'s row `j = m−1` (`S_{s−1}(u, 1)`, `u: m_s →
//! c_{s−1}`). Both rows interpolate `m_s` and the side's first and
//! second derivatives there (scaled by the half length `h`, opposite
//! first-derivative signs), which is exactly what the spoke conditions
//! (9)–(10) demand at `t = 1`; they interpolate their corner, and fit
//! the side in least squares between. A polynomial side of degree
//! `<= degree` per half (a strip's section spline) is reproduced
//! exactly.

use faer::Mat;

use super::linalg::lstsq;
use super::space::{SectorSpace, gauss_legendre};
use super::{VariationalCapError, VariationalCapSide};
use crate::geometry::curve::CurveEval;
use crate::geometry::{Point3, Vec3};

/// One half row: `m` control points from `m_s` to the corner.
pub(super) type HalfRow = Vec<Point3>;

/// Both half rows of one side.
#[derive(Debug, Clone)]
pub(super) struct SideRows {
    /// Leading half (sector s−1, row j = m−1), `m_s` first.
    pub leading: HalfRow,
    /// Trailing half (sector s, row i = m−1), `m_s` first.
    pub trailing: HalfRow,
    /// Side-curve parameters `(t0, tm, t1)`.
    pub params: (f64, f64, f64),
}

/// The map from a half row's parameter `x ∈ [0, 1]` (measured from
/// `m_s`) to the side curve: `t = tm ± h (x ± a x²)`, `h` the mean half
/// length and `a = (H_t − H_l) / (H_t + H_l)`. Both halves share the
/// 2-jet at `m_s` (first derivative negated), so the spoke conditions
/// hold exactly for an off-centre split; `a = 0` is the affine split.
#[derive(Debug, Clone, Copy)]
pub(super) struct HalfMap {
    pub tm: f64,
    pub h: f64,
    pub a: f64,
}

impl HalfMap {
    /// The split `(t0, tm, t1)`.
    pub fn new(t0: f64, tm: f64, t1: f64) -> Self {
        let (ht, hl) = (t1 - tm, tm - t0);
        HalfMap {
            tm,
            h: 0.5 * (ht + hl),
            a: (ht - hl) / (ht + hl),
        }
    }

    /// Side parameter of `x` on the trailing (`true`) or leading half.
    pub fn t(&self, trailing: bool, x: f64) -> f64 {
        if trailing {
            self.tm + self.h * (x + self.a * x * x)
        } else {
            self.tm - self.h * (x - self.a * x * x)
        }
    }
}

/// Fit both half rows of `side` in `space`.
pub(super) fn side_rows(
    space: &SectorSpace,
    side: &VariationalCapSide,
    index: usize,
    tm: f64,
) -> Result<SideRows, VariationalCapError> {
    let (t0, t1) = side.curve.parameter_range();
    let map = HalfMap::new(t0, tm, t1);
    let h = map.h;
    let (pm, c1, c2) = side.curve_jet(tm);
    let second = c2 * (h * h) + c1 * (2.0 * map.a * h);
    let trailing = half_row(
        space,
        |x| side.curve_at(map.t(true, x)),
        (pm, c1 * h, second),
        index,
    )?;
    let leading = half_row(
        space,
        |x| side.curve_at(map.t(false, x)),
        (pm, c1 * -h, second),
        index,
    )?;
    Ok(SideRows {
        leading,
        trailing,
        params: (t0, tm, t1),
    })
}

impl VariationalCapSide {
    /// `(C, C', C'')` at `t`.
    pub(super) fn curve_jet(&self, t: f64) -> (Point3, Vec3, Vec3) {
        (
            self.curve.point_at(t),
            self.curve.derivative_at(t, 1),
            self.curve.derivative_at(t, 2),
        )
    }

    fn curve_at(&self, t: f64) -> Point3 {
        self.curve.point_at(t)
    }
}

/// One half row for `r(x) = C(tm ± x·h)`, `x ∈ [0, 1]`, with the
/// exact 2-jet `(r(0), r'(0), r''(0))`. Both halves of a side receive
/// the same jet (first derivative negated), bit for bit, so the spoke
/// conditions at `m_s` are met exactly.
fn half_row(
    space: &SectorSpace,
    r: impl Fn(f64) -> Point3,
    jet: (Point3, Vec3, Vec3),
    index: usize,
) -> Result<HalfRow, VariationalCapError> {
    let d = space.degree;
    let m = space.m;
    let k = &space.knots;
    let (p0, d1, d2) = jet;
    let q0 = d1;
    let p1 = p0 + q0 * (k[d + 1] / d as f64);
    let q1 = q0 + d2 * (k[d + 1] / (d as f64 - 1.0));
    let p2 = p1 + q1 * (k[d + 2] / d as f64);
    let end = r(1.0);
    // Free control points 3..m−2 by least squares over Gauss samples.
    let free = m.saturating_sub(4);
    let mut row = vec![Point3::new(0.0, 0.0, 0.0); m];
    row[0] = p0;
    row[1] = p1;
    row[2] = p2;
    row[m - 1] = end;
    if free == 0 {
        return Ok(row);
    }
    let nodes = gauss_legendre(d + 3);
    let samples: Vec<f64> = space
        .spans()
        .iter()
        .flat_map(|&(a, b)| nodes.iter().map(move |&(x, _)| a + (b - a) * x))
        .collect();
    let mut a = Mat::<f64>::zeros(samples.len(), free);
    let mut rhs = Mat::<f64>::zeros(samples.len(), 3);
    for (s, &x) in samples.iter().enumerate() {
        let basis = space.basis(x);
        let mut known = Vec3::new(0.0, 0.0, 0.0);
        let target = r(x);
        for (ri, f) in basis.funcs.iter().enumerate() {
            let g = basis.first + ri;
            match g {
                0..=2 => known = known + to_vec(row[g]) * f[0],
                g if g == m - 1 => known = known + to_vec(row[g]) * f[0],
                g => a.write(s, g - 3, f[0]),
            }
        }
        let t = to_vec(target) - known;
        rhs.write(s, 0, t.x);
        rhs.write(s, 1, t.y);
        rhs.write(s, 2, t.z);
    }
    let sol = lstsq(a, rhs, 1e-12).ok_or(VariationalCapError::Singular {
        stage: "rim approximant",
    })?;
    for (c, slot) in row.iter_mut().enumerate().take(m - 1).skip(3) {
        let c = c - 3;
        *slot = Point3::new(sol.read(c, 0), sol.read(c, 1), sol.read(c, 2));
    }
    if row
        .iter()
        .any(|p| !(p.x.is_finite() && p.y.is_finite() && p.z.is_finite()))
    {
        return Err(VariationalCapError::DegenerateSide {
            side: index,
            length: 0.0,
        });
    }
    Ok(row)
}

fn to_vec(p: Point3) -> Vec3 {
    Vec3::new(p.x, p.y, p.z)
}
