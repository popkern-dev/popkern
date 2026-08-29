//! Edge curve recomputation for inner skin.

use std::collections::HashMap;

use crate::edge_curve_dispatch::compute_edge_curve;
use crate::error::ShellError;
use crate::skin_map::SkinMap;
use popkern_core::geometry::curve::CurveOnEdge;
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
    recompute_edges_with_skips(
        inner,
        offset_surfaces,
        skin_map,
        outer,
        tolerance,
        &std::collections::BTreeSet::new(),
    )
    .map(|_| ())
}

/// [`recompute_edges`] with gap-join integration (PK-325): edges in
/// `skip_edges` (arc-join gap edges — their clone curves anchor the
/// canal construction) are left untouched, and the returned list names
/// every interior edge whose whole recompute ladder declined silently
/// (the clone curve stayed stale) so the caller can convert gap-edge
/// declines into the typed `JoinFailed`.
pub fn recompute_edges_with_skips(
    inner: &mut BRep,
    offset_surfaces: &HashMap<FaceHandle, SurfaceGeometry>,
    skin_map: &SkinMap,
    outer: &BRep,
    tolerance: &Tolerance,
    skip_edges: &std::collections::BTreeSet<EdgeHandle>,
) -> Result<Vec<EdgeHandle>, ShellError> {
    let reverse_face = skin_map.reverse_face_map();
    let reverse_edge = skin_map.reverse_edge_map();

    let inner_edge_handles: Vec<EdgeHandle> = inner.edges.keys().collect();

    // Per-edge verdict from the parallel compute phase. `FreeRim`
    // declines are silent (the clone stays, no gap-join follow-up);
    // `Interior` declines feed the caller's `undone` list.
    enum EdgeOutcome {
        Skip,
        FreeRim(Option<CurveOnEdge>),
        Interior(Option<CurveOnEdge>),
    }

    // Phase 1 (parallel, PK-329): each edge's recompute reads the skin
    // immutably — its own clone curve, endpoint vertices (finalized in
    // the prior stage), and the offset surfaces — and never another
    // edge's freshly computed curve, so the per-edge solves dispatch in
    // parallel with an ordered collect and the sequential apply below
    // reproduces the single-loop behavior exactly (first failing edge's
    // error included). The per-edge timing trace forces serial dispatch
    // so thread interleaving cannot garble it.
    let trace = popkern_core::trace::trace_enabled("recompute");
    let compute_one = |inner_eh: EdgeHandle| -> Result<EdgeOutcome, ShellError> {
        if skip_edges.contains(&inner_eh) {
            return Ok(EdgeOutcome::Skip);
        }
        let edge = &inner.edges[inner_eh];
        let he_handle = edge.half_edge;

        if !inner.half_edges.contains_key(he_handle) {
            return Ok(EdgeOutcome::Skip);
        }

        // Find the two adjacent faces via the half-edge and its twin.
        let he = &inner.half_edges[he_handle];
        let twin_handle = he.twin;

        if !inner.half_edges.contains_key(twin_handle) {
            return Ok(EdgeOutcome::Skip);
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
                        return Ok(EdgeOutcome::Skip);
                    };
                    (s1.clone(), original.geometry.clone())
                }
                (None, Some(s2)) => {
                    let Some(original) = outer.faces.get(f1).and_then(|f| f.surface.as_ref())
                    else {
                        return Ok(EdgeOutcome::Skip);
                    };
                    (s2.clone(), original.geometry.clone())
                }
                (None, None) => return Ok(EdgeOutcome::Skip),
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
                let Some(kept_face) = kept else {
                    return Ok(EdgeOutcome::Skip);
                };
                let Some(kept_offset) = offset_surfaces.get(&kept_face) else {
                    return Ok(EdgeOutcome::Skip);
                };
                let Some(&outer_eh) = reverse_edge.get(&inner_eh) else {
                    return Ok(EdgeOutcome::Skip);
                };
                let Some(removed_face) = other_face_of_edge(outer, outer_eh, kept_face) else {
                    // FREE rim (open-sheet thicken, PK-324): no face was
                    // ever on the other side. The rim edge lies ON the
                    // kept face alone, so its offset image is the
                    // normal-offset map within that face — dispatch as
                    // a single-surface (seam-like) edge; the
                    // parallel-map arm serves it exactly.
                    let (s1, s2) = (kept_offset.clone(), kept_offset.clone());
                    return Ok(EdgeOutcome::FreeRim(compute_edge_curve(
                        inner, inner_eh, &s1, &s2, tolerance,
                    )?));
                };
                let Some(original) = outer
                    .faces
                    .get(removed_face)
                    .and_then(|f| f.surface.as_ref())
                else {
                    return Ok(EdgeOutcome::Skip);
                };
                (kept_offset.clone(), original.geometry.clone())
            }
            (None, None) => return Ok(EdgeOutcome::Skip),
        };

        // Compute edge curve from surface intersection.
        let trace_started = trace.then(std::time::Instant::now);
        let computed = compute_edge_curve(inner, inner_eh, &surf1, &surf2, tolerance)?;
        if let Some(t0) = trace_started {
            eprintln!(
                "[recomp-edge] {inner_eh:?} {}x{} {:.1} ms served={}",
                surface_kind(&surf1),
                surface_kind(&surf2),
                t0.elapsed().as_secs_f64() * 1e3,
                computed.is_some()
            );
        }
        Ok(EdgeOutcome::Interior(computed))
    };

    let outcomes = crate::parallel_dispatch::dispatch_ordered(&inner_edge_handles, trace, |&eh| {
        compute_one(eh)
    });

    // Phase 2 (sequential): apply curve writes and collect declines in
    // edge order — identical to the original single loop.
    let mut undone: Vec<EdgeHandle> = Vec::new();
    for (inner_eh, outcome) in inner_edge_handles.iter().zip(outcomes) {
        match outcome? {
            EdgeOutcome::Skip => {}
            EdgeOutcome::FreeRim(computed) => {
                if let Some(curve_on_edge) = computed {
                    inner.edges[*inner_eh].curve = Some(curve_on_edge);
                }
            }
            EdgeOutcome::Interior(computed) => {
                if let Some(curve_on_edge) = computed {
                    inner.edges[*inner_eh].curve = Some(curve_on_edge);
                } else {
                    undone.push(*inner_eh);
                }
            }
        }
    }

    Ok(undone)
}

/// Surface-class label for the the `recompute` trace topic per-edge timing dump.
fn surface_kind(s: &SurfaceGeometry) -> &'static str {
    match s {
        SurfaceGeometry::Plane { .. } => "plane",
        SurfaceGeometry::Cylinder { .. } => "cylinder",
        SurfaceGeometry::Cone { .. } => "cone",
        SurfaceGeometry::Sphere { .. } => "sphere",
        SurfaceGeometry::Torus { .. } => "torus",
        SurfaceGeometry::Nurbs(_) => "nurbs",
        SurfaceGeometry::Gregory(_) => "gregory",
        SurfaceGeometry::PolyhedralSpline(_) => "polyhedral_spline",
    }
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::classify::classify_faces;
    use crate::edge_curve_dispatch::solve_2x2;
    use crate::inner_skin::build_inner_skin;
    use crate::surface_offset::compute_offset_surfaces;
    use crate::types::ShellDirection;
    use crate::vertex_recompute::recompute_vertices;
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::curve::CurveGeometry;
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
