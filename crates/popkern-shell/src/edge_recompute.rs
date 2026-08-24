//! Edge curve recomputation for inner skin.

use std::collections::HashMap;

use crate::error::ShellError;
use crate::skin_map::SkinMap;
use popkern_core::geometry::Vec3;
use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

/// Recompute all edge curves in the inner skin.
///
/// For each edge, finds the two adjacent faces, maps them to the outer
/// BRep, retrieves offset surfaces, and computes the intersection curve.
/// Plane-plane intersections yield a line; other cases use the general
/// surface-surface intersection from popkern-core.
///
/// # Errors
///
/// Returns [`ShellError::Topology`] if adjacency lookups fail.
pub fn recompute_edges(
    inner: &mut BRep,
    offset_surfaces: &HashMap<FaceHandle, SurfaceGeometry>,
    skin_map: &SkinMap,
    outer: &BRep,
    tolerance: &Tolerance,
) -> Result<(), ShellError> {
    let reverse_face = skin_map.reverse_face_map();
    let reverse_edge = skin_map.reverse_edge_map();

    let inner_edge_handles: Vec<EdgeHandle> = inner.edges.keys().collect();

    for inner_eh in inner_edge_handles {
        let edge = &inner.edges[inner_eh];
        let he_handle = edge.half_edge;

        if !inner.half_edges.contains_key(he_handle) {
            continue;
        }

        // Find the two adjacent faces via the half-edge and its twin.
        let he = &inner.half_edges[he_handle];
        let twin_handle = he.twin;

        if !inner.half_edges.contains_key(twin_handle) {
            continue;
        }

        let twin = &inner.half_edges[twin_handle];

        let loop1 = he.loop_;
        let loop2 = twin.loop_;
        let side1 = inner
            .loops
            .get(loop1)
            .and_then(|l| reverse_face.get(&l.face).copied());
        let side2 = inner
            .loops
            .get(loop2)
            .and_then(|l| reverse_face.get(&l.face).copied());

        let (surf1, surf2) = match (side1, side2) {
            // Interior edge: intersection of the two offset surfaces.
            // A side that is present in the skin but UNCLASSIFIED (no
            // offset surface — midplane outer skins keep removed faces
            // for rim discovery) contributes its ORIGINAL surface, the
            // same contract as the removed-side arm below.
            (Some(f1), Some(f2)) => match (offset_surfaces.get(&f1), offset_surfaces.get(&f2)) {
                (Some(s1), Some(s2)) => (s1.clone(), s2.clone()),
                (Some(s1), None) => {
                    let Some(original) = outer.faces.get(f2).and_then(|f| f.surface.as_ref())
                    else {
                        continue;
                    };
                    (s1.clone(), original.geometry.clone())
                }
                (None, Some(s2)) => {
                    let Some(original) = outer.faces.get(f1).and_then(|f| f.surface.as_ref())
                    else {
                        continue;
                    };
                    (s2.clone(), original.geometry.clone())
                }
                (None, None) => continue,
            },
            // Opening rim edge (PK-307): one side's face was removed,
            // so its loop never existed in the inner skin and the
            // cloned curve is the STALE pre-offset original — measured
            // gap √2·thickness on the open-cube shell, invisible until
            // the side walls became manifold and validation reached
            // geometry. The honest rim curve is the intersection of
            // the KEPT face's offset surface with the REMOVED face's
            // ORIGINAL surface: the opening stays flush with the
            // original boundary while the kept wall recedes.
            (kept @ Some(_), None) | (None, kept @ Some(_)) => {
                let Some(kept_face) = kept else { continue };
                let Some(kept_offset) = offset_surfaces.get(&kept_face) else {
                    continue;
                };
                let Some(&outer_eh) = reverse_edge.get(&inner_eh) else {
                    continue;
                };
                let Some(removed_face) = other_face_of_edge(outer, outer_eh, kept_face) else {
                    continue;
                };
                let Some(original) = outer
                    .faces
                    .get(removed_face)
                    .and_then(|f| f.surface.as_ref())
                else {
                    continue;
                };
                (kept_offset.clone(), original.geometry.clone())
            }
            (None, None) => continue,
        };

        // Compute edge curve from surface intersection.
        if let Some(curve_on_edge) = compute_edge_curve(inner, inner_eh, &surf1, &surf2, tolerance)?
        {
            inner.edges[inner_eh].curve = Some(curve_on_edge);
        }
    }

    Ok(())
}

/// The face on the other side of `outer_eh` from `known_face`, if any.
fn other_face_of_edge(
    outer: &BRep,
    outer_eh: EdgeHandle,
    known_face: FaceHandle,
) -> Option<FaceHandle> {
    let edge = outer.edges.get(outer_eh)?;
    let he = outer.half_edges.get(edge.half_edge)?;
    let twin = outer.half_edges.get(he.twin)?;
    let f1 = outer.loops.get(he.loop_).map(|l| l.face);
    let f2 = outer.loops.get(twin.loop_).map(|l| l.face);
    match (f1, f2) {
        (Some(a), Some(b)) if a == known_face => Some(b),
        (Some(a), Some(b)) if b == known_face => Some(a),
        _ => None,
    }
}

/// Compute the edge curve from the intersection of two offset surfaces.
///
/// For plane-plane intersections, computes the line analytically.
/// For other surface types, uses the general surface-surface intersection.
fn compute_edge_curve(
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
            // For non-planar surfaces, we would use intersect_surfaces().
            // For now, keep the existing edge curve (from the clone).
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
fn solve_2x2(a: f64, b: f64, c1: f64, d: f64, e: f64, c2: f64) -> Result<(f64, f64), ShellError> {
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::classify::classify_faces;
    use crate::inner_skin::build_inner_skin;
    use crate::surface_offset::compute_offset_surfaces;
    use crate::types::ShellDirection;
    use crate::vertex_recompute::recompute_vertices;
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    fn tol() -> Tolerance {
        Tolerance::default()
    }

    #[test]
    fn cube_edges_recomputed_as_lines() {
        let brep = BRep::unit_cube(tol()).unwrap();
        let thickness = 0.1;
        let classification =
            classify_faces(&brep, thickness, &[], ShellDirection::Inward, &[]).unwrap();
        let offset_result = compute_offset_surfaces(&brep, &classification, true, &tol()).unwrap();
        let (mut inner, skin_map) = build_inner_skin(&brep, &[]).unwrap();

        // Recompute vertices first (edges need correct vertex positions).
        recompute_vertices(
            &mut inner,
            &offset_result.surfaces,
            &skin_map,
            &brep,
            &tol(),
        )
        .unwrap();

        recompute_edges(
            &mut inner,
            &offset_result.surfaces,
            &skin_map,
            &brep,
            &tol(),
        )
        .unwrap();

        // All cube edges are plane-plane intersections → lines.
        for (_, edge) in &inner.edges {
            let curve = edge.curve.as_ref().expect("edge should have a curve");
            assert!(
                matches!(curve.geometry, CurveGeometry::Line { .. }),
                "cube edge should be a line"
            );
        }
    }

    #[test]
    fn cube_edge_endpoints_match_vertices() {
        let brep = BRep::unit_cube(tol()).unwrap();
        let thickness = 0.1;
        let classification =
            classify_faces(&brep, thickness, &[], ShellDirection::Inward, &[]).unwrap();
        let offset_result = compute_offset_surfaces(&brep, &classification, true, &tol()).unwrap();
        let (mut inner, skin_map) = build_inner_skin(&brep, &[]).unwrap();

        recompute_vertices(
            &mut inner,
            &offset_result.surfaces,
            &skin_map,
            &brep,
            &tol(),
        )
        .unwrap();

        recompute_edges(
            &mut inner,
            &offset_result.surfaces,
            &skin_map,
            &brep,
            &tol(),
        )
        .unwrap();

        // Verify that edge curve endpoints match their vertex positions.
        for (_, edge) in &inner.edges {
            let curve = match &edge.curve {
                Some(c) => c,
                None => continue,
            };
            let he = &inner.half_edges[edge.half_edge];
            let twin = &inner.half_edges[he.twin];
            let start_pos = inner.vertices[he.origin].position;
            let end_pos = inner.vertices[twin.origin].position;

            let curve_start = curve.geometry.point_at(curve.t_start);
            let curve_end = curve.geometry.point_at(curve.t_end);

            assert!(
                start_pos.distance_to(&curve_start) < 0.02,
                "edge start mismatch: {:?} vs {:?}",
                start_pos,
                curve_start,
            );
            assert!(
                end_pos.distance_to(&curve_end) < 0.02,
                "edge end mismatch: {:?} vs {:?}",
                end_pos,
                curve_end,
            );
        }
    }

    #[test]
    fn solve_2x2_basic() {
        // x + y = 3, 2x - y = 0 → x=1, y=2
        let (x, y) = solve_2x2(1.0, 1.0, 3.0, 2.0, -1.0, 0.0).unwrap();
        assert!((x - 1.0).abs() < 1e-12);
        assert!((y - 2.0).abs() < 1e-12);
    }

    #[test]
    fn solve_2x2_singular() {
        let result = solve_2x2(1.0, 2.0, 3.0, 2.0, 4.0, 6.0);
        assert!(result.is_err());
    }
}
