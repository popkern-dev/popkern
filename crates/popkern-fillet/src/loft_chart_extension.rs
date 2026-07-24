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
    let step = |a: Point3, b: Point3| -> Point3 {
        Point3::new(2.0 * a.x - b.x, 2.0 * a.y - b.y, 2.0 * a.z - b.z)
    };
    // ONE exact extrapolated row per end. Deeper extension was
    // measured and REGRESSED: 3-4 rows per end moved the boss
    // corner-chord gauge from 3.5e-3 to 3.7e-2 (the wider chart
    // shifts later refit/region choices), and the residual it
    // targeted is not a U-coverage gap — the corner vertex sits off
    // the AS-BUILT twin strip laterally (see corner_chord_settle).
    let mut head_rows: Vec<([Point3; 3], f64)> = Vec::new();
    let mut tail_rows: Vec<([Point3; 3], f64)> = Vec::new();
    for (rows, edge_index, inner_index) in [
        (&mut head_rows, 0usize, 1usize),
        (&mut tail_rows, n - 1, n - 2),
    ] {
        let outer = row(&loft.control_points, edge_index);
        let inner = row(&loft.control_points, inner_index);
        if let Some((new_row, w)) = solve_row(step(outer[0], inner[0]), step(outer[2], inner[2])) {
            rows.push((new_row, w));
        }
    }
    if head_rows.is_empty() && tail_rows.is_empty() {
        return;
    }
    let added = head_rows.len() + tail_rows.len();
    let mut cps = Vec::with_capacity((n + added) * 3);
    let mut ws = Vec::with_capacity((n + added) * 3);
    for (r, w) in head_rows.iter().rev() {
        cps.extend_from_slice(r);
        ws.extend_from_slice(&[1.0, *w, 1.0]);
    }
    cps.extend_from_slice(&loft.control_points);
    ws.extend_from_slice(&loft.weights);
    for (r, w) in &tail_rows {
        cps.extend_from_slice(r);
        ws.extend_from_slice(&[1.0, *w, 1.0]);
    }
    loft.control_points = cps;
    loft.weights = ws;
    loft.num_u = n + added;
    loft.knots_u = build_clamped_uniform_knots(n + added, loft.degree_u);
}
