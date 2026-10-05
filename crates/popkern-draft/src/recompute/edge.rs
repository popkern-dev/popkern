//! Edge curve recomputation for drafted faces.
//!
//! For each edge between a drafted and an undrafted face (or between two
//! drafted faces), computes the new edge curve from the intersection of
//! the two adjacent surfaces. Pairs without a recompute path keep the
//! original curve and record [`DraftWarning::EdgeCurveStale`] — the
//! explicit fallback boundary (measured: a drafted cylinder→cone kept
//! its r-of-cylinder rim circle while the cone flared past it, and the
//! inconsistent body meshed with 638 cracked seam edges).

use popkern_math::PortableMath;
use std::collections::HashMap;

use crate::error::{DraftError, DraftWarning};
use popkern_core::geometry::Vec3;
use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

/// Recompute edge curves affected by drafted face surfaces.
///
/// For each edge adjacent to at least one drafted face, computes the new
/// curve from the intersection of the two adjacent surfaces. Returns the
/// edges whose curves were recomputed (input for the vertex snap pass).
///
/// # Errors
///
/// Returns [`DraftError::Topology`] if adjacency lookups fail.
pub fn recompute_draft_edges(
    brep: &mut BRep,
    drafted: &HashMap<FaceHandle, SurfaceGeometry>,
    moved_vertices: &[popkern_core::topology::VertexHandle],
    tolerance: &Tolerance,
    warnings: &mut Vec<DraftWarning>,
) -> Result<Vec<EdgeHandle>, DraftError> {
    let mut recomputed = Vec::new();
    if drafted.is_empty() {
        return Ok(recomputed);
    }

    let edge_handles: Vec<EdgeHandle> = brep.edges.keys().collect();

    for eh in edge_handles {
        let edge = &brep.edges[eh];
        let he_handle = edge.half_edge;

        if !brep.half_edges.contains_key(he_handle) {
            continue;
        }

        let he = &brep.half_edges[he_handle];
        let twin_handle = he.twin;

        if !brep.half_edges.contains_key(twin_handle) {
            continue;
        }

        let twin = &brep.half_edges[twin_handle];

        let loop1 = he.loop_;
        let loop2 = twin.loop_;

        if !brep.loops.contains_key(loop1) || !brep.loops.contains_key(loop2) {
            continue;
        }

        let face1 = brep.loops[loop1].face;
        let face2 = brep.loops[loop2].face;

        // An edge needs its curve rebuilt when a drafted face flanks
        // it — OR when either endpoint vertex moved. The second-ring
        // case (both faces undrafted, a corner vertex shared with the
        // drafted face relocated) previously kept its STALE curve and
        // shipped a validate-red ShellGap the size of the draft
        // displacement (measured: tan(angle) on every plane-neutral
        // cube draft; invisible to open-edge counting, which is what
        // the historical tests checked).
        let f1_drafted = drafted.contains_key(&face1);
        let f2_drafted = drafted.contains_key(&face2);
        let endpoint_moved = {
            let origin = he.origin;
            let twin_origin = twin.origin;
            moved_vertices.contains(&origin) || moved_vertices.contains(&twin_origin)
        };
        if !f1_drafted && !f2_drafted && !endpoint_moved {
            continue;
        }

        // Get surfaces (drafted if available, otherwise original).
        let surf1 = get_surface(brep, face1, drafted);
        let surf2 = get_surface(brep, face2, drafted);

        let (surf1, surf2) = match (surf1, surf2) {
            (Some(s1), Some(s2)) => (s1, s2),
            _ => continue,
        };

        // SEAM edge (both half-edges on one face): on a drafted
        // cone the seam is the apex ruling at the seam azimuth —
        // azimuth is preserved because the draft moves points
        // within their meridian plane.
        if face1 == face2 {
            if let SurfaceGeometry::Cone {
                apex,
                axis,
                half_angle,
            } = &surf1
            {
                if let Some(curve) = cone_seam_ruling(brep, eh, apex, axis, *half_angle) {
                    brep.edges[eh].curve = Some(curve);
                    recomputed.push(eh);
                } else {
                    warnings.push(DraftWarning::EdgeCurveStale { edge: eh });
                }
            }
            continue;
        }

        // Compute new edge curve.
        match compute_edge_curve(brep, eh, &surf1, &surf2, tolerance)? {
            EdgeCurveOutcome::Recomputed(curve_on_edge) => {
                brep.edges[eh].curve = Some(curve_on_edge);
                recomputed.push(eh);
            }
            EdgeCurveOutcome::Degenerate => {}
            EdgeCurveOutcome::UnsupportedPair => {
                warnings.push(DraftWarning::EdgeCurveStale { edge: eh });
            }
        }
    }

    Ok(recomputed)
}

/// Outcome of a single edge-curve recompute attempt.
enum EdgeCurveOutcome {
    /// A new curve was derived from the adjacent surfaces.
    Recomputed(CurveOnEdge),
    /// The pair is degenerate (e.g. parallel planes); the old curve
    /// is geometrically unaffected and kept without a warning.
    Degenerate,
    /// No recompute path exists for this surface pair — the old
    /// curve is kept and MAY be stale (warned by the caller).
    UnsupportedPair,
}

/// The apex ruling of a cone at the seam azimuth of the edge's
/// current curve. Returns `None` when the edge has no curve to take
/// the azimuth from.
fn cone_seam_ruling(
    brep: &BRep,
    edge_handle: EdgeHandle,
    apex: &popkern_core::geometry::Point3,
    axis: &Vec3,
    half_angle: f64,
) -> Option<CurveOnEdge> {
    let edge = &brep.edges[edge_handle];
    let old = edge.curve.as_ref()?;
    // Azimuth and axial station of the old seam's midpoint. The
    // azimuth and station survive the draft (points move within
    // their meridian); only the radius changes.
    let (p1, p2) = popkern_core::geometry::analytic_curves::circle_basis(*axis);
    let sample = old.geometry.point_at(0.5 * (old.t_start + old.t_end));
    let d = sample - *apex;
    let phi = d.dot(&p2).patan2(d.dot(&p1));
    let station = d.dot(axis);
    if station.abs() <= 0.0 {
        return None; // sample at the apex — azimuth undefined
    }
    // Point of the NEW cone at that station and azimuth; the ruling
    // is apex → that point. Deriving it from the surface makes the
    // direction correct for either apex-side convention (draft's
    // cone need not point its axis toward growing radius).
    let radial = p1 * phi.pcos() + p2 * phi.psin();
    let on_cone = *apex + *axis * station + radial * (station.abs() * half_angle.ptan());
    let direction = (on_cone - *apex).normalized()?;
    // Preserve the old parameter direction.
    let old_dir = old.geometry.point_at(old.t_end) - old.geometry.point_at(old.t_start);
    let direction = if old_dir.dot(&direction) < 0.0 {
        direction * -1.0
    } else {
        direction
    };
    let he = &brep.half_edges[edge.half_edge];
    let twin = &brep.half_edges[he.twin];
    let start = brep.vertices[he.origin].position;
    let end = brep.vertices[twin.origin].position;
    Some(CurveOnEdge {
        geometry: CurveGeometry::Line {
            origin: *apex,
            direction,
        },
        t_start: (start - *apex).dot(&direction),
        t_end: (end - *apex).dot(&direction),
    })
}

/// Get the surface for a face, preferring the drafted version.
fn get_surface(
    brep: &BRep,
    face: FaceHandle,
    drafted: &HashMap<FaceHandle, SurfaceGeometry>,
) -> Option<SurfaceGeometry> {
    if let Some(surf) = drafted.get(&face) {
        return Some(surf.clone());
    }
    brep.faces
        .get(face)
        .and_then(|f| f.surface.as_ref())
        .map(|sof| sof.geometry.clone())
}

/// Compute the edge curve from the intersection of two surfaces.
fn compute_edge_curve(
    brep: &BRep,
    edge_handle: EdgeHandle,
    surf1: &SurfaceGeometry,
    surf2: &SurfaceGeometry,
    tolerance: &Tolerance,
) -> Result<EdgeCurveOutcome, DraftError> {
    let edge = &brep.edges[edge_handle];
    let he = &brep.half_edges[edge.half_edge];
    let twin = &brep.half_edges[he.twin];
    let start_pos = brep.vertices[he.origin].position;
    let end_pos = brep.vertices[twin.origin].position;

    match (surf1, surf2) {
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
            let direction = n1.cross(n2);
            let dir_len = direction.length();
            if dir_len < tolerance.solver() {
                return Ok(EdgeCurveOutcome::Degenerate);
            }
            let direction = direction * (1.0 / dir_len);

            let d1 = n1.dot(&Vec3::new(o1.x, o1.y, o1.z));
            let d2 = n2.dot(&Vec3::new(o2.x, o2.y, o2.z));
            let point_on_line = solve_point_on_plane_plane(*n1, d1, *n2, d2, direction)?;

            let t_start = (start_pos - point_on_line).dot(&direction);
            let t_end = (end_pos - point_on_line).dot(&direction);

            Ok(EdgeCurveOutcome::Recomputed(CurveOnEdge {
                geometry: CurveGeometry::Line {
                    origin: point_on_line,
                    direction,
                },
                t_start,
                t_end,
            }))
        }
        (
            SurfaceGeometry::Cone {
                apex,
                axis,
                half_angle,
            },
            SurfaceGeometry::Plane {
                origin: po,
                normal: pn,
                ..
            },
        )
        | (
            SurfaceGeometry::Plane {
                origin: po,
                normal: pn,
                ..
            },
            SurfaceGeometry::Cone {
                apex,
                axis,
                half_angle,
            },
        ) => {
            // Perpendicular case only: plane normal parallel to the
            // cone axis → the intersection is the exact rim circle
            // at that axial station. Oblique cone×plane (a conic) has
            // no recompute path yet.
            let n = match pn.normalized() {
                Some(n) => n,
                None => return Ok(EdgeCurveOutcome::Degenerate),
            };
            if n.cross(axis).length() > 1.0e-9 {
                return Ok(EdgeCurveOutcome::UnsupportedPair);
            }
            let denom = axis.dot(&n);
            if denom.abs() < tolerance.solver() {
                return Ok(EdgeCurveOutcome::Degenerate);
            }
            let t_axis = (*po - *apex).dot(&n) / denom;
            let radius = t_axis.abs() * half_angle.ptan();
            if radius < tolerance.linear() {
                return Ok(EdgeCurveOutcome::Degenerate);
            }
            let center = *apex + *axis * t_axis;
            // Preserve the old curve's winding: the Circle evaluator
            // runs counter-clockwise about its normal, so pick the
            // normal that keeps the edge's parameter direction.
            let old_ccw = edge
                .curve
                .as_ref()
                .map(|old| {
                    let tmid = 0.5 * (old.t_start + old.t_end);
                    let p = old.geometry.point_at(tmid);
                    let tangent = old.geometry.derivative_at(tmid, 1);
                    axis.cross(&(p - center)).dot(&tangent) >= 0.0
                })
                .unwrap_or(true);
            let normal = if old_ccw { *axis } else { *axis * -1.0 };
            // Anchor the angular range at the origin vertex's
            // azimuth (azimuth survives the draft; radius does not).
            let (p1, p2) = popkern_core::geometry::analytic_curves::circle_basis(normal);
            let d = start_pos - center;
            let a0 = d.dot(&p2).patan2(d.dot(&p1));
            let closed = he.origin == twin.origin;
            let a1 = if closed {
                a0 + std::f64::consts::TAU
            } else {
                let d_end = end_pos - center;
                let mut a = d_end.dot(&p2).patan2(d_end.dot(&p1));
                while a <= a0 {
                    a += std::f64::consts::TAU;
                }
                a
            };
            Ok(EdgeCurveOutcome::Recomputed(CurveOnEdge {
                geometry: CurveGeometry::Circle {
                    center,
                    normal,
                    radius,
                },
                t_start: a0,
                t_end: a1,
            }))
        }
        _ => Ok(EdgeCurveOutcome::UnsupportedPair),
    }
}

/// Find a point on the intersection of two planes.
fn solve_point_on_plane_plane(
    n1: Vec3,
    d1: f64,
    n2: Vec3,
    d2: f64,
    direction: Vec3,
) -> Result<popkern_core::geometry::Point3, DraftError> {
    let ax = direction.x.abs();
    let ay = direction.y.abs();
    let az = direction.z.abs();

    let (x, y, z) = if ax >= ay && ax >= az {
        let det = n1.y * n2.z - n1.z * n2.y;
        if det.abs() < 1e-10 {
            return Err(DraftError::Topology(
                "degenerate plane-plane intersection".into(),
            ));
        }
        let y = (d1 * n2.z - n1.z * d2) / det;
        let z = (n1.y * d2 - d1 * n2.y) / det;
        (0.0, y, z)
    } else if ay >= ax && ay >= az {
        let det = n1.x * n2.z - n1.z * n2.x;
        if det.abs() < 1e-10 {
            return Err(DraftError::Topology(
                "degenerate plane-plane intersection".into(),
            ));
        }
        let x = (d1 * n2.z - n1.z * d2) / det;
        let z = (n1.x * d2 - d1 * n2.x) / det;
        (x, 0.0, z)
    } else {
        let det = n1.x * n2.y - n1.y * n2.x;
        if det.abs() < 1e-10 {
            return Err(DraftError::Topology(
                "degenerate plane-plane intersection".into(),
            ));
        }
        let x = (d1 * n2.y - n1.y * d2) / det;
        let y = (n1.x * d2 - d1 * n2.x) / det;
        (x, y, 0.0)
    };

    Ok(popkern_core::geometry::Point3::new(x, y, z))
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::Point3;
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    fn tol() -> Tolerance {
        Tolerance::default()
    }

    #[test]
    fn recompute_empty_drafted_is_noop() {
        let mut brep = BRep::unit_cube(tol()).unwrap();
        let drafted = HashMap::new();
        recompute_draft_edges(&mut brep, &drafted, &[], &tol(), &mut Vec::new()).unwrap();
    }

    #[test]
    fn recompute_drafted_cube_face_edges() {
        let mut brep = BRep::unit_cube(tol()).unwrap();
        let face = brep.faces.keys().next().unwrap();

        let tilted = SurfaceGeometry::Plane {
            origin: Point3::new(0.5, 0.5, 0.0),
            normal: Vec3::new(0.0, 0.087, 0.996),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 0.996, -0.087),
        };

        let mut drafted = HashMap::new();
        drafted.insert(face, tilted);

        let result = recompute_draft_edges(&mut brep, &drafted, &[], &tol(), &mut Vec::new());
        assert!(result.is_ok());
    }

    #[test]
    fn solve_point_on_plane_plane_basic() {
        let n1 = Vec3::new(0.0, 0.0, 1.0);
        let n2 = Vec3::new(0.0, 1.0, 0.0);
        let direction = n1.cross(&n2);
        let dir_norm = direction.normalized().unwrap();

        let result = solve_point_on_plane_plane(n1, 0.5, n2, 0.3, dir_norm).unwrap();
        // n1.p = 0.5 => z = 0.5
        // n2.p = 0.3 => y = 0.3
        assert!((result.z - 0.5).abs() < 1e-10);
        assert!((result.y - 0.3).abs() < 1e-10);
    }
}
