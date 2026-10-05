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
use popkern_math::PortableMath;

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
    let m = loft.num_v;
    if n < 3 || m < 3 || loft.control_points.len() != n * m {
        return;
    }
    // Two supported layouts: the exact rational 3-row loft (mid CP
    // weighted cos(θ/2)) and the sampled-skin fallback (m unit-
    // weight arc samples per station).
    let nonunit = loft.weights.iter().any(|w| (w - 1.0).abs() > 1e-12);
    let rational3 = m == 3 && nonunit;
    let bi5 = m == 5 && nonunit;
    let row = |cps: &[Point3], i: usize| -> [Point3; 3] {
        [cps[i * m], cps[i * m + m / 2], cps[i * m + m - 1]]
    };
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
        let theta = na.dot(&nb).clamp(-1.0, 1.0).pacos();
        if !(1e-9..=std::f64::consts::PI * 0.9).contains(&theta) {
            return None;
        }
        let w = (theta * 0.5).pcos();
        let bis = (na + nb).normalized()?;
        let mid = centre + bis * (radius / w);
        Some(([f0, mid, f2], w))
    };
    let step = |a: Point3, b: Point3| -> Point3 {
        Point3::new(2.0 * a.x - b.x, 2.0 * a.y - b.y, 2.0 * a.z - b.z)
    };
    // Extend each end until the chart reaches PAST the corner-ball
    // station: cumulative added rail arc >= 1.2 r per end (measured:
    // one row ≈ 0.5 r left the corner contact W OUTSIDE both charts,
    // clamped at 2.4e-2/3.7e-2; the earlier one-sided multi-row
    // experiment regressed because the sampled-skin twin could not
    // extend at all — both layouts extend now). Cap 4 rows per end.
    let target = radius * 1.2;
    let expand_row = |contacts: &([Point3; 3], f64)| -> (Vec<Point3>, Vec<f64>) {
        let ([f0, mid, f2], w) = contacts;
        if rational3 {
            (vec![*f0, *mid, *f2], vec![1.0, *w, 1.0])
        } else if bi5 {
            // Exact two-segment rational row (see
            // bi_rational_arc_row): recover the station centre from
            // the solved contacts, rebuild the piecewise mids.
            let Some(d) = ((*f0 - *mid) + (*f2 - *mid)).normalized() else {
                return (Vec::new(), Vec::new());
            };
            let centre = *mid + d * (radius / *w);
            let (Some(u1), Some(u2)) = ((*f0 - centre).normalized(), (*f2 - centre).normalized())
            else {
                return (Vec::new(), Vec::new());
            };
            let theta = u1.dot(&u2).clamp(-1.0, 1.0).pacos();
            let Some(bis) = (u1 + u2).normalized() else {
                return (Vec::new(), Vec::new());
            };
            let arc_mid = centre + bis * radius;
            let wq = (theta * 0.25).pcos();
            let (Some(m1d), Some(m2d)) = ((u1 + bis).normalized(), (u2 + bis).normalized()) else {
                return (Vec::new(), Vec::new());
            };
            let m1 = centre + m1d * (radius / wq);
            let m2 = centre + m2d * (radius / wq);
            (vec![*f0, m1, arc_mid, m2, *f2], vec![1.0, wq, 1.0, wq, 1.0])
        } else {
            // Sampled-skin layout: re-sample the exact station arc
            // at the row's arc density. The station centre is
            // recoverable from the exact contacts and mid.
            let centre_dir = ((*f0 - *mid) + (*f2 - *mid)).normalized();
            let centre = match centre_dir {
                Some(d) => {
                    // mid = centre + bis * (r/w) with bis unit: walk
                    // back from mid toward the chord's centre side.
                    let bis_len = radius / *w;
                    *mid + d * (bis_len - radius).max(0.0) + d * radius
                }
                None => return (Vec::new(), Vec::new()),
            };
            let pts = crate::fillet_surface::sample_arc(&centre, f0, f2, radius, m);
            if pts.len() != m {
                return (Vec::new(), Vec::new());
            }
            (pts, vec![1.0; m])
        }
    };
    let mut head_rows: Vec<(Vec<Point3>, Vec<f64>)> = Vec::new();
    let mut tail_rows: Vec<(Vec<Point3>, Vec<f64>)> = Vec::new();
    for (rows, edge_index, inner_index) in [
        (&mut head_rows, 0usize, 1usize),
        (&mut tail_rows, n - 1, n - 2),
    ] {
        let mut outer = row(&loft.control_points, edge_index);
        let mut inner = row(&loft.control_points, inner_index);
        let mut added = 0.0;
        for _ in 0..4 {
            let Some(contacts) = solve_row(step(outer[0], inner[0]), step(outer[2], inner[2]))
            else {
                break;
            };
            let (pts, ws) = expand_row(&contacts);
            if pts.is_empty() {
                break;
            }
            added += contacts.0[0]
                .distance_to(&outer[0])
                .max(contacts.0[2].distance_to(&outer[2]));
            inner = outer;
            outer = contacts.0;
            rows.push((pts, ws));
            if added >= target {
                break;
            }
        }
    }
    if head_rows.is_empty() && tail_rows.is_empty() {
        return;
    }
    let added = head_rows.len() + tail_rows.len();
    let mut cps = Vec::with_capacity((n + added) * m);
    let mut ws = Vec::with_capacity((n + added) * m);
    for (r, w) in head_rows.iter().rev() {
        cps.extend_from_slice(r);
        ws.extend_from_slice(w);
    }
    cps.extend_from_slice(&loft.control_points);
    ws.extend_from_slice(&loft.weights);
    for (r, w) in &tail_rows {
        cps.extend_from_slice(r);
        ws.extend_from_slice(w);
    }
    loft.control_points = cps;
    loft.weights = ws;
    loft.num_u = n + added;
    loft.knots_u = build_clamped_uniform_knots(n + added, loft.degree_u);
}
