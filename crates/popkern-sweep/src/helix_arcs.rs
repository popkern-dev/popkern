//! The constant-pitch helix path (PK-675).
//!
//! A helix is not a rational curve, so its NURBS is an approximation:
//! the question is where the error sits and how small it is.
//! - **What it was.** Rational quadratic quarter arcs with linearly
//!   placed z control points. They were exact circles in x and y, but a
//!   rational arc's parameter does not advance linearly in angle, so
//!   z(θ) oscillated by up to 2.64e-3 at pitch 0.5 (0.5 % of the pitch,
//!   measured ±2.4e-3 on the M3 nut's tool wall). Every thread built on
//!   it was off its helicoid by that much.
//! - **What it is.** Each arc keeps the exact circle and is degree-
//!   elevated to [`HELIX_DEGREE`], which is still the same circle in
//!   x and y. That frees control points for z, which are fitted by
//!   least squares to the helix's linear z(θ) with both ends pinned.
//!   Every arc is the same rotated shape, so joints meet with matching
//!   tangents. With [`ARCS_PER_TURN`] = 4 the largest z error is
//!   3.2e-7 at pitch 0.5, about 10⁴ times smaller, and it scales with
//!   the pitch.
//! - **The swept wall.** The path is only half of it: the sweep places
//!   its frames through an arc-length table, which must resolve the
//!   speed variation inside each arc (see
//!   `frame_budget::arc_length_samples`). Together, the thread tool
//!   walls sit within 6e-6 of their helicoids axially (M3..M10,
//!   `popkern-fastener` `helicoid_fidelity_tests`).
//! - **Why four arcs per turn.** Eight halve the z error again
//!   (2.5e-9) and the radial skinning residue, but the frame budget
//!   gives every knot span six frames, so they double a thread tool
//!   wall's stations: bolt and nut builds ran 1.4 to 4.4 times
//!   slower (M10 bolt 13.8 s -> 60.6 s) for an axial gain far below
//!   any thread tolerance.

use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_math::PortableMath;
use popkern_math::linalg::{DenseMat, HouseholderQr};

use crate::error::SweepError;
use crate::types::HelixSpec;

/// Degree of each helix arc.
pub(crate) const HELIX_DEGREE: usize = 5;
/// Arcs per turn.
pub(crate) const ARCS_PER_TURN: usize = 4;
/// Samples per arc for the z fit.
const FIT_SAMPLES: usize = 64;

/// The Bernstein basis of degree `p` at `t`.
fn bernstein(p: usize, t: f64) -> Vec<f64> {
    let mut b = vec![0.0; p + 1];
    b[0] = 1.0;
    for j in 1..=p {
        let mut saved = 0.0;
        for bk in b.iter_mut().take(j) {
            let temp = *bk;
            *bk = saved + (1.0 - t) * temp;
            saved = t * temp;
        }
        b[j] = saved;
    }
    b
}

/// One template arc from angle 0 to `2·half` (`half` signed: negative
/// winds clockwise) on the unit circle, as homogeneous control points
/// `(w x, w y, w)` of degree [`HELIX_DEGREE`], and the homogeneous z
/// coefficients `w z` of a rise from 0 to 1 that tracks the angle.
fn template_arc(half: f64) -> Option<(Vec<[f64; 3]>, Vec<f64>)> {
    let w = half.pcos();
    // The rational quadratic arc, homogeneous.
    let mut hom: Vec<[f64; 3]> = vec![
        [1.0, 0.0, 1.0],
        [w, w * half.ptan(), w],
        [(2.0 * half).pcos(), (2.0 * half).psin(), 1.0],
    ];
    // Degree elevation (exact: the same rational curve).
    for degree in 2..HELIX_DEGREE {
        let mut next = Vec::with_capacity(degree + 2);
        next.push(hom[0]);
        for i in 1..=degree {
            let a = i as f64 / (degree + 1) as f64;
            next.push([
                a * hom[i - 1][0] + (1.0 - a) * hom[i][0],
                a * hom[i - 1][1] + (1.0 - a) * hom[i][1],
                a * hom[i - 1][2] + (1.0 - a) * hom[i][2],
            ]);
        }
        next.push(hom[degree]);
        hom = next;
    }
    let p = HELIX_DEGREE;
    // z(t)·W(t) = Σ Z_i B_i(t) with Z_0 = 0, Z_p = w_p (z = 1 at the
    // end), the interior Z fitted to the angle's share of the rise.
    let free = p - 1;
    let mut rows = DenseMat::zeros(FIT_SAMPLES, free);
    let mut rhs = DenseMat::zeros(FIT_SAMPLES, 1);
    for k in 0..FIT_SAMPLES {
        let t = (k as f64 + 0.5) / FIT_SAMPLES as f64;
        let b = bernstein(p, t);
        let (mut x, mut y, mut big_w) = (0.0, 0.0, 0.0);
        for (bi, h) in b.iter().zip(&hom) {
            x += bi * h[0];
            y += bi * h[1];
            big_w += bi * h[2];
        }
        let theta = (y / big_w).patan2(x / big_w);
        let target = theta / (2.0 * half) * big_w - b[p] * hom[p][2];
        for (j, bj) in b[1..p].iter().enumerate() {
            rows.set(k, j, *bj);
        }
        rhs.set(k, 0, target);
    }
    let qr = HouseholderQr::factor_left(rows, free);
    if (0..free).any(|i| qr.r_diag(i).abs() < 1e-12) {
        return None;
    }
    let solved = qr.solve(rhs);
    let mut z = vec![0.0; p + 1];
    for (j, zj) in z[1..p].iter_mut().enumerate() {
        *zj = solved.get(j, 0);
    }
    z[p] = hom[p][2];
    Some((hom, z))
}

/// The constant-pitch, constant-radius helix of `spec` as degree
/// [`HELIX_DEGREE`] arcs, [`ARCS_PER_TURN`] per turn.
pub(crate) fn build_constant_helix(spec: &HelixSpec, pitch: f64) -> Result<NurbsCurve, SweepError> {
    let axis = spec
        .axis_dir
        .normalized()
        .ok_or_else(|| SweepError::InvalidPath("helix axis is zero-length".into()))?;
    let radius = spec.start_radius;
    let (u_dir, v_dir) = super::helix::orthonormal_pair(&axis);
    let arcs = (spec.turns * ARCS_PER_TURN as f64).ceil() as usize;
    if arcs == 0 {
        return Err(SweepError::InvalidPath("helix has no turns".into()));
    }
    // A left-hand helix winds clockwise about the axis while still
    // advancing along it.
    let winding = if spec.left_handed { -1.0 } else { 1.0 };
    let total_angle = winding * spec.turns * 2.0 * std::f64::consts::PI;
    let total_height = spec.turns * pitch;
    let sweep = total_angle / arcs as f64;
    let rise = total_height / arcs as f64;
    let (template, z_template) = template_arc(0.5 * sweep)
        .ok_or_else(|| SweepError::InvalidPath("helix arc fit is singular".into()))?;
    let p = HELIX_DEGREE;
    let mut control_points = Vec::with_capacity(arcs * p + 1);
    let mut weights = Vec::with_capacity(arcs * p + 1);
    for a in 0..arcs {
        let (s, c) = (spec.start_angle + a as f64 * sweep).psin_cos();
        let base = a as f64 * rise;
        let first = if a == 0 { 0 } else { 1 };
        for i in first..=p {
            let [hx, hy, w] = template[i];
            let (x, y) = (hx / w, hy / w);
            let (xr, yr) = (c * x - s * y, s * x + c * y);
            let z = base + rise * z_template[i] / w;
            let pos: Point3 =
                spec.axis_origin + u_dir * (radius * xr) + v_dir * (radius * yr) + axis * z;
            control_points.push(pos);
            weights.push(w);
        }
    }
    let mut knots = vec![0.0; p + 1];
    for a in 1..arcs {
        let k = a as f64 / arcs as f64;
        knots.extend(std::iter::repeat_n(k, p));
    }
    knots.extend(std::iter::repeat_n(1.0, p + 1));
    Ok(NurbsCurve {
        degree: p,
        knots,
        control_points,
        weights,
    })
}

#[cfg(test)]
#[path = "helix_arcs_tests.rs"]
mod tests;
