//! Per-edge offset-curve solve dispatch: analytic plane×plane line,
//! seam ladder, and the curved-arm ladder (exact arms first, general
//! numeric SSI last). Split from `edge_recompute.rs` (PK-329) to honour
//! the 500-line structural cap.

use crate::error::ShellError;
use popkern_core::geometry::Vec3;
use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle};

/// Compute the edge curve from the intersection of two offset surfaces.
///
/// For plane-plane intersections, computes the line analytically.
/// For other surface types, uses the general surface-surface intersection.
pub(crate) fn compute_edge_curve(
    inner: &BRep,
    edge_handle: EdgeHandle,
    surf1: &SurfaceGeometry,
    surf2: &SurfaceGeometry,
    tolerance: &Tolerance,
) -> Result<Option<CurveOnEdge>, ShellError> {
    // Get edge endpoint positions for trimming.
    let edge = &inner.edges[edge_handle];
    let he = &inner.half_edges[edge.half_edge];
    let twin = &inner.half_edges[he.twin];
    let start_pos = inner.vertices[he.origin].position;
    let end_pos = inner.vertices[twin.origin].position;

    match (surf1, surf2) {
        // Same-surface edges (seams, free rims of open sheets) go to
        // the seam ladder even for planes — the plane×plane arm below
        // would read an identical pair as "parallel, no intersection"
        // and leave the clone stale.
        _ if surf1 == surf2 => {
            let Some(old) = inner.edges[edge_handle].curve.clone() else {
                return Ok(None);
            };
            if let Some(c) = crate::curved_recompute::straight_ruling_seam_curve(
                &old, surf1, start_pos, end_pos, tolerance,
            ) {
                return Ok(Some(c));
            }
            let mapped = crate::curved_recompute_similarity::parallel_map_edge_curve(
                &old, surf1, surf2, start_pos, end_pos, tolerance,
            );
            if mapped.is_none() && std::env::var_os("PK_RECOMP_TRACE").is_some() {
                eprintln!(
                    "[recomp-seam-decline] {edge_handle:?} start=({:.3},{:.3},{:.3}) \
                     end=({:.3},{:.3},{:.3})",
                    start_pos.x, start_pos.y, start_pos.z, end_pos.x, end_pos.y, end_pos.z
                );
            }
            Ok(mapped)
        }
        (
            SurfaceGeometry::Plane {
                origin: o1,
                normal: n1,
                ..
            },
            SurfaceGeometry::Plane {
                origin: o2,
                normal: n2,
                ..
            },
        ) => {
            // Plane-plane intersection: line direction is n1 x n2.
            let direction = n1.cross(n2);
            let dir_len = direction.length();
            if dir_len < tolerance.solver() {
                // Parallel planes — no intersection line (shouldn't happen
                // for adjacent faces of a valid solid).
                return Ok(None);
            }
            let direction = direction * (1.0 / dir_len);

            // Find a point on the intersection line.
            let d1 = n1.dot(&Vec3::new(o1.x, o1.y, o1.z));
            let d2 = n2.dot(&Vec3::new(o2.x, o2.y, o2.z));
            let point_on_line = solve_point_on_plane_plane(*n1, d1, *n2, d2, direction)?;

            // Trim to edge vertex endpoints by projecting onto the line.
            let t_start = (start_pos - point_on_line).dot(&direction);
            let t_end = (end_pos - point_on_line).dot(&direction);

            Ok(Some(CurveOnEdge {
                geometry: CurveGeometry::Line {
                    origin: point_on_line,
                    direction,
                },
                t_start,
                t_end,
            }))
        }
        _ => {
            // Curved-arm ladder (PK-322 exact arms + PK-324 general
            // arms), cheapest/most-exact first. Every arm verifies its
            // candidate by sampled residuals and recomputed-vertex
            // endpoints, and declines to the next; when the whole
            // ladder declines the clone curve stays stale and the
            // serve gate declines the body.
            let Some(old) = inner.edges[edge_handle].curve.clone() else {
                return Ok(None);
            };
            if let Some(c) = crate::curved_recompute::coaxial_circle_edge_curve(
                &old, surf1, surf2, start_pos, tolerance,
            ) {
                return Ok(Some(c));
            }
            // G1-tangential pairs and remaining seams: both offset
            // surfaces displace the edge identically along the shared
            // normal (projection-agreement witness).
            if let Some(c) = crate::curved_recompute_similarity::parallel_map_edge_curve(
                &old, surf1, surf2, start_pos, end_pos, tolerance,
            ) {
                return Ok(Some(c));
            }
            // Transversal coaxial revolution pairs (bore×countersink
            // junction circles): 1D profile-plane solve.
            if let Some(c) = crate::curved_recompute_similarity::coaxial_profile_edge_curve(
                &old, surf1, surf2, start_pos, end_pos, tolerance,
            ) {
                return Ok(Some(c));
            }
            // General numeric SSI (design §5 Stage 2): pairs with no
            // analytic form, e.g. fitted NURBS offsets against planes.
            if let Some(c) = crate::curved_recompute_similarity::general_ssi_edge_curve(
                &old, surf1, surf2, start_pos, end_pos, tolerance,
            ) {
                return Ok(Some(c));
            }
            Ok(None)
        }
    }
}

/// Solve for a point on the intersection of two planes.
///
/// Given `n1 . p = d1` and `n2 . p = d2`, fixes the coordinate with the
/// largest component of `direction` to zero and solves the 2x2 system.
fn solve_point_on_plane_plane(
    n1: Vec3,
    d1: f64,
    n2: Vec3,
    d2: f64,
    direction: Vec3,
) -> Result<popkern_core::geometry::Point3, ShellError> {
    let ax = direction.x.abs();
    let ay = direction.y.abs();
    let az = direction.z.abs();

    let (x, y, z) = if ax >= ay && ax >= az {
        let (y, z) = solve_2x2(n1.y, n1.z, d1, n2.y, n2.z, d2)?;
        (0.0, y, z)
    } else if ay >= ax && ay >= az {
        let (x, z) = solve_2x2(n1.x, n1.z, d1, n2.x, n2.z, d2)?;
        (x, 0.0, z)
    } else {
        let (x, y) = solve_2x2(n1.x, n1.y, d1, n2.x, n2.y, d2)?;
        (x, y, 0.0)
    };

    Ok(popkern_core::geometry::Point3::new(x, y, z))
}

/// Solve a 2x2 linear system: `a*x + b*y = c1`, `d*x + e*y = c2`.
pub(crate) fn solve_2x2(
    a: f64,
    b: f64,
    c1: f64,
    d: f64,
    e: f64,
    c2: f64,
) -> Result<(f64, f64), ShellError> {
    let det = a * e - b * d;
    if det.abs() < f64::EPSILON * 1e6 {
        return Err(ShellError::Topology(
            "degenerate 2x2 system in plane-plane intersection".into(),
        ));
    }
    let x = (c1 * e - b * c2) / det;
    let y = (a * c2 - c1 * d) / det;
    Ok((x, y))
}
