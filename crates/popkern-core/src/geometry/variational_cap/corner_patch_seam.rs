//! The seam between a corner patch and the main cap (PK-619): the cut
//! curve K and the ribbon of tangent planes both pieces meet along it.

use super::VariationalCapSide;
use crate::geometry::curve::CurveEval;
use crate::geometry::fitting::interpolate_curve;
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::nurbs_surface::NurbsSurface;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::uv_projection::project_point_to_uv;
use crate::geometry::{Point3, Vec3};

/// Samples along K that the seam curve and the ribbon interpolate.
const SEAM_SAMPLES: usize = 24;

/// The seam over thin side `t`: K from p_B (fraction `tail` of the next
/// side) to p_A (fraction 1 − `tail` of the previous side), and the
/// ruled ribbon `K(s) + w·(N(s) × K̂(s))` whose natural normal along K is
/// the seam normal N. N follows the thin side's oriented neighbour
/// normal backwards along T, corrected linearly to the neighbours'
/// normals at p_B and p_A, so each end of the ribbon agrees with the
/// side it meets.
pub(super) fn corner_seam(
    sides: &[VariationalCapSide],
    t: usize,
    tail: f64,
    bulge: f64,
) -> Result<(NurbsCurve, SurfaceGeometry), String> {
    let n = sides.len();
    let (a, thin, b) = (&sides[(t + n - 1) % n], &sides[t], &sides[(t + 1) % n]);
    let pa = at(&a.curve, 1.0 - tail);
    let pb = at(&b.curve, tail);
    let na = oriented_normal(a, &pa).ok_or("the previous side's neighbour has no normal at p_A")?;
    let nb = oriented_normal(b, &pb).ok_or("the next side's neighbour has no normal at p_B")?;
    // Sides run counter-clockwise about the outward normal: N × t
    // points into the opening.
    let into_a = na
        .cross(&unit_tangent(&a.curve, 1.0 - tail))
        .normalized()
        .ok_or("the previous side is tangent to its neighbour's normal")?;
    let into_b = nb
        .cross(&unit_tangent(&b.curve, tail))
        .normalized()
        .ok_or("the next side is tangent to its neighbour's normal")?;
    let handle = bulge * pa.distance_to(&pb);
    let ctrl = [pb, pb + into_b * handle, pa + into_a * handle, pa];
    let start_normal = oriented_normal(thin, &at(&thin.curve, 1.0))
        .ok_or("the thin side's neighbour has no normal at its end")?;
    let end_normal = oriented_normal(thin, &at(&thin.curve, 0.0))
        .ok_or("the thin side's neighbour has no normal at its start")?;
    let (fix_b, fix_a) = (nb - start_normal, na - end_normal);
    let params: Vec<f64> = (0..=SEAM_SAMPLES)
        .map(|i| i as f64 / SEAM_SAMPLES as f64)
        .collect();
    let mut points = Vec::with_capacity(params.len());
    let mut across = Vec::with_capacity(params.len());
    for &s in &params {
        let (p, d) = bezier(&ctrl, s);
        let k = d.normalized().ok_or("the seam curve is singular")?;
        let base = oriented_normal(thin, &at(&thin.curve, 1.0 - s))
            .ok_or_else(|| format!("the thin side's neighbour has no normal at {s:.3}"))?;
        let guess = base + fix_b * (1.0 - s) + fix_a * s;
        let normal = (guess - k * guess.dot(&k))
            .normalized()
            .ok_or_else(|| format!("the seam normal is tangent to K at {s:.3}"))?;
        let w = normal.cross(&k);
        points.push(p);
        across.push(Point3::new(w.x, w.y, w.z));
    }
    let seam = interpolate_curve(&points, &params, 3).map_err(|e| format!("seam curve: {e}"))?;
    let width = interpolate_curve(&across, &params, 3).map_err(|e| format!("ribbon: {e}"))?;
    // Same parameters and degree: the two interpolants share their knots,
    // so the ruled net is exact.
    let h = pa.distance_to(&pb);
    let m = seam.control_points.len();
    let mut net = Vec::with_capacity(2 * m);
    for (k, w) in seam.control_points.iter().zip(&width.control_points) {
        let w = Vec3::new(w.x, w.y, w.z);
        net.push(*k + w * -h);
        net.push(*k + w * h);
    }
    let ribbon = NurbsSurface {
        degree_u: 3,
        degree_v: 1,
        knots_u: seam.knots.clone(),
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points: net,
        num_u: m,
        num_v: 2,
        weights: vec![1.0; 2 * m],
    };
    Ok((seam, SurfaceGeometry::Nurbs(ribbon)))
}

/// `curve` restricted to the fractions `[lo, hi]` of its parameter range.
pub(super) fn restrict(curve: &NurbsCurve, lo: f64, hi: f64) -> Result<NurbsCurve, String> {
    let (a, b) = curve.parameter_range();
    let mut piece = curve.clone();
    if hi < 1.0 {
        piece = piece
            .split_at(a + (b - a) * hi)
            .map_err(|e| format!("split at {hi:.3}: {e}"))?
            .0;
    }
    if lo > 0.0 {
        piece = piece
            .split_at(a + (b - a) * lo)
            .map_err(|e| format!("split at {lo:.3}: {e}"))?
            .1;
    }
    Ok(piece)
}

/// The same curve traversed backwards over the same parameter range.
pub(super) fn reversed(curve: &NurbsCurve) -> NurbsCurve {
    let (a, b) = curve.parameter_range();
    NurbsCurve {
        degree: curve.degree,
        knots: curve.knots.iter().rev().map(|k| a + b - k).collect(),
        control_points: curve.control_points.iter().rev().copied().collect(),
        weights: curve.weights.iter().rev().copied().collect(),
    }
}

/// Polyline length of `curve` over 64 chords.
pub(super) fn length(curve: &NurbsCurve) -> f64 {
    const CHORDS: usize = 64;
    (0..CHORDS)
        .map(|k| {
            at(curve, k as f64 / CHORDS as f64)
                .distance_to(&at(curve, (k + 1) as f64 / CHORDS as f64))
        })
        .sum()
}

/// The point at fraction `f` of `curve`'s parameter range.
fn at(curve: &NurbsCurve, f: f64) -> Point3 {
    let (a, b) = curve.parameter_range();
    curve.point_at(a + (b - a) * f)
}

/// The unit tangent at fraction `f`.
fn unit_tangent(curve: &NurbsCurve, f: f64) -> Vec3 {
    let (a, b) = curve.parameter_range();
    curve
        .derivative_at(a + (b - a) * f, 1)
        .normalized()
        .unwrap_or(Vec3::new(0.0, 0.0, 0.0))
}

/// The side's neighbour normal at the foot of `p`, oriented outward.
fn oriented_normal(side: &VariationalCapSide, p: &Point3) -> Option<Vec3> {
    let (u, v) = project_point_to_uv(&side.neighbour, p)?;
    let normal = side.neighbour.normal_at(u, v);
    Some(if side.neighbour_reversed {
        normal * -1.0
    } else {
        normal
    })
}

/// The cubic Bézier `ctrl` and its derivative at `s`.
fn bezier(ctrl: &[Point3; 4], s: f64) -> (Point3, Vec3) {
    let m = 1.0 - s;
    let w = [m * m * m, 3.0 * m * m * s, 3.0 * m * s * s, s * s * s];
    let o = ctrl[0];
    let mut p = Vec3::new(0.0, 0.0, 0.0);
    for (c, wi) in ctrl.iter().zip(w) {
        p = p + (*c - o) * wi;
    }
    let d = (ctrl[1] - ctrl[0]) * (3.0 * m * m)
        + (ctrl[2] - ctrl[1]) * (6.0 * m * s)
        + (ctrl[3] - ctrl[2]) * (3.0 * s * s);
    (o + p, d)
}
