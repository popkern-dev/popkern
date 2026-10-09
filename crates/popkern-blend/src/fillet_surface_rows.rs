//! Exact cross-section row primitives for the fillet surface sweep:
//! the rational-quadratic arc row, the sampled arc fallback, and the
//! clamped-uniform knot builder, the collapsed row and the strip row
//! points. Split out of `fillet_surface.rs` (500-line cap).

use popkern_core::geometry::{Point3, Vec3};
use popkern_math::PortableMath;

use super::FilletSurfaceResult;
use crate::center_curve::CenterCurveResult;

/// One EXACT rational-quadratic arc row from `start` to `end` about
/// `center`: control points [start, tangent-intersection, end] with
/// mid weight cos(θ/2). Returns `None` when the section is outside
/// the well-conditioned single-segment range (θ > ~115°, radii
/// mismatched, or degenerate).
pub(crate) type ArcRow = ([Point3; 3], f64);

pub(crate) fn rational_arc_row(center: &Point3, start: &Point3, end: &Point3) -> Option<ArcRow> {
    let v1 = *start - *center;
    let v2 = *end - *center;
    let (r1, r2) = (v1.length(), v2.length());
    if r1 < 1e-12 || r2 < 1e-12 {
        return None;
    }
    if (r1 - r2).abs() > 1e-6 * r1.max(r2) {
        return None;
    }
    let u1 = v1 * (1.0 / r1);
    let u2 = v2 * (1.0 / r2);
    let cos_theta = u1.dot(&u2).clamp(-1.0, 1.0);
    let theta = cos_theta.pacos();
    if !(1e-9..=2.0).contains(&theta) {
        return None;
    }
    let half = 0.5 * theta;
    let w_mid = half.pcos();
    let mid_dir = (u1 + u2).normalized()?;
    let r = 0.5 * (r1 + r2);
    let mid = *center + mid_dir * (r / w_mid);
    Some(([*start, mid, *end], w_mid))
}

/// Sample a circular arc from `start` to `end` centered at `center`.
///
/// The arc lies in the plane defined by the three points. If `start` and
/// `end` are equidistant from `center` at approximately `radius`, the arc
/// sweeps the shorter angle between them.
pub(crate) fn sample_arc(
    center: &Point3,
    start: &Point3,
    end: &Point3,
    _radius: f64,
    n: usize,
) -> Vec<Point3> {
    let v1 = *start - *center;
    let v2 = *end - *center;
    let r1 = v1.length();
    let r2 = v2.length();

    // Handle degenerate cases.
    if r1 < 1e-15 || r2 < 1e-15 || n < 2 {
        return vec![*start; n];
    }

    // Normalize to unit vectors.
    let u = Vec3::new(v1.x / r1, v1.y / r1, v1.z / r1);
    let w_raw = v2 - u * v2.dot(&u);
    let w_len = w_raw.length();

    if w_len < 1e-15 {
        // start and end are nearly coincident or antipodal — just linearly interpolate.
        let mut pts = Vec::with_capacity(n);
        for i in 0..n {
            let t = i as f64 / (n - 1).max(1) as f64;
            pts.push(Point3::new(
                start.x + (end.x - start.x) * t,
                start.y + (end.y - start.y) * t,
                start.z + (end.z - start.z) * t,
            ));
        }
        return pts;
    }

    let v = Vec3::new(w_raw.x / w_len, w_raw.y / w_len, w_raw.z / w_len);

    // Compute sweep angle.
    let cos_theta = (v1.dot(&v2) / (r1 * r2)).clamp(-1.0, 1.0);
    let theta = cos_theta.pacos();

    // Use average radius for the arc.
    let r = (r1 + r2) * 0.5;

    let mut pts = Vec::with_capacity(n);
    for i in 0..n {
        let t = i as f64 / (n - 1).max(1) as f64;
        let angle = theta * t;
        let x = r * angle.pcos();
        let y = r * angle.psin();
        pts.push(Point3::new(
            center.x + u.x * x + v.x * y,
            center.y + u.y * x + v.y * y,
            center.z + u.z * x + v.z * y,
        ));
    }

    pts
}

/// Build clamped uniform knot vector for `n` control points at given degree.
pub(crate) fn build_clamped_uniform_knots(n: usize, degree: usize) -> Vec<f64> {
    let n_interior = n.saturating_sub(degree);
    let mut knots = Vec::with_capacity(n + degree + 1);
    knots.extend(std::iter::repeat_n(0.0, degree + 1));
    for i in 1..n_interior {
        knots.push(i as f64 / n_interior as f64);
    }
    knots.extend(std::iter::repeat_n(1.0, degree + 1));
    knots
}

/// Push the exact degenerate arc of a COLLAPSED section — both contacts
/// one point, a strip's vanishing end at a tangent contact (PK-764):
/// one point, unit weights. Declining it threw the whole strip onto the
/// sampled skin, whose rows sat 6.5e-3 off the exact junction section.
pub(crate) fn push_collapsed_row(
    point: &Point3,
    rational_rows: &mut Option<Vec<ArcRow>>,
    bi_rows: &mut Option<Vec<([Point3; 5], f64)>>,
    sampled_points: &mut Vec<Point3>,
    n_arc: usize,
) {
    if let Some(rows) = rational_rows.as_mut() {
        rows.push(([*point; 3], 1.0));
    }
    if let Some(rows) = bi_rows.as_mut() {
        rows.push(([*point; 5], 1.0));
    }
    sampled_points.extend(std::iter::repeat_n(*point, n_arc));
}

/// The points the strip's boundary rows pass through, which the trims
/// must follow. The G1 sweeps start their sections at the clipped
/// centre-curve footprints. The G2 sweep (and its exact revolution)
/// starts them at the walker's contact solutions, a different
/// construction: footprint trims left a G2 pentagonal-pyramid base
/// strip's rail 1.25e-2 off (PK-531). The G2 sweep ran exactly when G2
/// was requested, the walker supplied contacts, and no downgrade was
/// recorded.
pub(crate) fn strip_row_points<'a>(
    center_result: &'a CenterCurveResult,
    fillet_result: &FilletSurfaceResult,
    wants_g2: bool,
) -> (
    std::borrow::Cow<'a, [Point3]>,
    std::borrow::Cow<'a, [Point3]>,
) {
    use std::borrow::Cow;
    match center_result
        .contact_solutions
        .as_deref()
        .filter(|s| wants_g2 && s.len() >= 2 && fillet_result.continuity_downgrade.is_none())
    {
        Some(solutions) => (
            Cow::Owned(solutions.iter().map(|c| c.contact0).collect()),
            Cow::Owned(solutions.iter().map(|c| c.contact1).collect()),
        ),
        None => (
            Cow::Borrowed(center_result.footprint_1.as_slice()),
            Cow::Borrowed(center_result.footprint_2.as_slice()),
        ),
    }
}
