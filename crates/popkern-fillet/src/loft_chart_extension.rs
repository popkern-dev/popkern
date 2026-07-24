//! PK-141 chart overlap: extend a rational-loft strip chart by one
//! EXACT extrapolated station row at each end, so the chart covers
//! the corner mutual-trim chords beyond the chain stations. Rows
//! are re-projected onto the supports and the station centre is
//! re-solved (equidistant `radius` by alternating offset
//! projection), so the extension is exact station geometry rather
//! than extrapolation error; trims and anchors are untouched
//! (walk-based window extension was falsified three ways).

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{NurbsSurface, Point3};

use crate::fillet_surface::build_clamped_uniform_knots;

/// Extend a rational-loft strip chart by one exact row at each end:
/// linearly extrapolate the boundary station, re-project the rail
/// control points onto their supports, re-solve the station centre
/// (equidistant `r` by alternating offset projection) and rebuild
/// the exact rational arc row.
pub(crate) fn extend_loft_rows(
    loft: &mut NurbsSurface,
    surface1: &SurfaceGeometry,
    surface2: &SurfaceGeometry,
    radius: f64,
) {
    use popkern_core::geometry::SurfaceEval;
    let n = loft.num_u;
    if n < 3 || loft.num_v != 3 || loft.control_points.len() != n * 3 {
        return;
    }
    let row =
        |cps: &[Point3], i: usize| -> [Point3; 3] { [cps[i * 3], cps[i * 3 + 1], cps[i * 3 + 2]] };
    let project = |sg: &SurfaceGeometry, p: Point3| -> Option<Point3> {
        let (u, v) = popkern_core::geometry::uv_projection::project_point_to_uv(sg, &p)?;
        Some(sg.point_at(u, v))
    };
    let solve_row = |c0_seed: Point3, c2_seed: Point3| -> Option<([Point3; 3], f64)> {
        let c0 = project(surface1, c0_seed)?;
        let c2 = project(surface2, c2_seed)?;
        // Station centre: equidistant r from both supports, seeded
        // from the rail midpoint pushed outward.
        let mid_seed = Point3::new(
            0.5 * (c0.x + c2.x),
            0.5 * (c0.y + c2.y),
            0.5 * (c0.z + c2.z),
        );
        let mut centre = mid_seed;
        for _ in 0..200 {
            let q1 = project(surface1, centre)?;
            let d1 = centre.distance_to(&q1);
            if d1 > 1e-12 {
                centre = q1 + (centre - q1) * (radius / d1);
            }
            let q2 = project(surface2, centre)?;
            let d2 = centre.distance_to(&q2);
            if d2 > 1e-12 {
                centre = q2 + (centre - q2) * (radius / d2);
            }
            let r1 = (centre.distance_to(&project(surface1, centre)?) - radius).abs();
            let r2 = (centre.distance_to(&project(surface2, centre)?) - radius).abs();
            if r1.max(r2) < 1e-10 {
                break;
            }
        }
        // Exact contacts for the solved centre.
        let f0 = project(surface1, centre)?;
        let f2 = project(surface2, centre)?;
        let na = (f0 - centre).normalized()?;
        let nb = (f2 - centre).normalized()?;
        let theta = na.dot(&nb).clamp(-1.0, 1.0).acos();
        if !(1e-9..=std::f64::consts::PI * 0.9).contains(&theta) {
            return None;
        }
        let w = (theta * 0.5).cos();
        let bis = (na + nb).normalized()?;
        let mid = centre + bis * (radius / w);
        Some(([f0, mid, f2], w))
    };
    let first = row(&loft.control_points, 0);
    let second = row(&loft.control_points, 1);
    let last = row(&loft.control_points, n - 1);
    let prev = row(&loft.control_points, n - 2);
    let step = |a: Point3, b: Point3| -> Point3 {
        Point3::new(2.0 * a.x - b.x, 2.0 * a.y - b.y, 2.0 * a.z - b.z)
    };
    let head = solve_row(step(first[0], second[0]), step(first[2], second[2]));
    let tail = solve_row(step(last[0], prev[0]), step(last[2], prev[2]));
    let (Some((head_row, head_w)), Some((tail_row, tail_w))) = (head, tail) else {
        return;
    };
    let mut cps = Vec::with_capacity((n + 2) * 3);
    let mut ws = Vec::with_capacity((n + 2) * 3);
    cps.extend_from_slice(&head_row);
    ws.extend_from_slice(&[1.0, head_w, 1.0]);
    cps.extend_from_slice(&loft.control_points);
    ws.extend_from_slice(&loft.weights);
    cps.extend_from_slice(&tail_row);
    ws.extend_from_slice(&[1.0, tail_w, 1.0]);
    loft.control_points = cps;
    loft.weights = ws;
    loft.num_u = n + 2;
    loft.knots_u = build_clamped_uniform_knots(n + 2, loft.degree_u);
}
