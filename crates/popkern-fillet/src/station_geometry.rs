//! Exact fillet station geometry (PK-177).
//!
//! Shared solvers for EXACT station rows: the alternating-offset
//! centre solve (equidistant `r` from both supports, ~1e-10) and
//! the two-segment piecewise-rational arc row for sweeps beyond the
//! single-segment rational range. The chart extension and the loft
//! builder both consume these so every station — interior,
//! extension, wide-arc — is exact station geometry rather than
//! walker-accuracy samples.

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::SurfaceGeometry;

/// Re-solve one station's support contacts exactly: alternating
/// offset projection to the centre equidistant `radius` from both
/// supports (seeded between the walker footprints), then the exact
/// foot points. `None` when the solve does not converge or wanders
/// more than half a radius from the seeds.
pub(crate) fn refine_station_contacts(
    s1: &SurfaceGeometry,
    s2: &SurfaceGeometry,
    fp1: &Point3,
    fp2: &Point3,
    radius: f64,
) -> Option<(Point3, Point3)> {
    use popkern_core::geometry::SurfaceEval;
    let project = |sg: &SurfaceGeometry, p: Point3| -> Option<Point3> {
        let (u, v) = popkern_core::geometry::uv_projection::project_point_to_uv(sg, &p)?;
        Some(sg.point_at(u, v))
    };
    let seed = Point3::new(
        0.5 * (fp1.x + fp2.x),
        0.5 * (fp1.y + fp2.y),
        0.5 * (fp1.z + fp2.z),
    );
    let mut centre = seed;
    let mut converged = false;
    for _ in 0..200 {
        let q1 = project(s1, centre)?;
        let d1 = centre.distance_to(&q1);
        if d1 > 1e-12 {
            centre = q1 + (centre - q1) * (radius / d1);
        }
        let q2 = project(s2, centre)?;
        let d2 = centre.distance_to(&q2);
        if d2 > 1e-12 {
            centre = q2 + (centre - q2) * (radius / d2);
        }
        let r1 = (centre.distance_to(&project(s1, centre)?) - radius).abs();
        let r2 = (centre.distance_to(&project(s2, centre)?) - radius).abs();
        if r1.max(r2) < 1e-10 {
            converged = true;
            break;
        }
    }
    if !converged {
        return None;
    }
    let c1 = project(s1, centre)?;
    let c2 = project(s2, centre)?;
    // Containment: an exact solve should land within footprint
    // accuracy of the walker's contacts; a slide along the support
    // means a different station solved.
    if c1.distance_to(fp1) > radius * 0.5 || c2.distance_to(fp2) > radius * 0.5 {
        return None;
    }
    Some((c1, c2))
}

/// One EXACT two-segment piecewise-rational arc row from `start` to
/// `end` about `center`: control points [start, M1, mid, M2, end]
/// with segment-mid weights cos(θ/4) and a doubled interior knot at
/// 0.5. Exact for sweeps the single-segment row declines
/// (~115°..170°); `None` outside its own well-conditioned range.
pub(crate) fn bi_rational_arc_row(
    center: &Point3,
    start: &Point3,
    end: &Point3,
) -> Option<([Point3; 5], f64)> {
    let v1 = *start - *center;
    let v2 = *end - *center;
    let r1 = v1.length();
    let r2 = v2.length();
    if r1 < 1e-12 || r2 < 1e-12 {
        return None;
    }
    // Row is exact for the MEAN circle; a footprint radius mismatch
    // of ε leaves an ε/2 residual. 1e-3 relative admits walker
    // footprints (~1e-5 class) at a sub-precision residual — still
    // an order under the sampled skin's chordal error.
    if (r1 - r2).abs() > 1e-3 * r1.max(r2) {
        return None;
    }
    let u1 = v1 * (1.0 / r1);
    let u2 = v2 * (1.0 / r2);
    let cos_theta = u1.dot(&u2).clamp(-1.0, 1.0);
    let theta = cos_theta.acos();
    if !(1e-9..=3.0).contains(&theta) {
        return None;
    }
    let r = 0.5 * (r1 + r2);
    let bis = (u1 + u2).normalized()?;
    let mid = *center + bis * r;
    let w = (theta * 0.25).cos();
    let m1 = *center + (u1 + bis).normalized()? * (r / w);
    let m2 = *center + (u2 + bis).normalized()? * (r / w);
    Some(([*start, m1, mid, m2, *end], w))
}
