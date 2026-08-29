//! §8.3 end-to-end integration tests for cap insertion on open-shell
//! topology (`BRep::open_box`). Split out of `cap_insert_tests.rs`
//! to stay under the 500-line structural cap.
//!
//! Unit-level tests for `endpoint_is_closed`, `insert_cap_face`, and
//! the closed-body cap-warning paths live in the sibling file.

use super::*;
use crate::types::{FilletOptions, FilletRadius, FilletSpec};
use crate::{CapKind, fillet_edges};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Return all interior edges of `open_box` with at least one
/// endpoint on the boundary, sorted by handle for determinism.
fn boundary_adjacent_interior_edges(ob: &BRep) -> Vec<popkern_core::topology::EdgeHandle> {
    let mut edges: Vec<_> = ob
        .edges
        .iter()
        .filter(|(eh, _)| {
            if ob.boundary_edges.contains(eh) {
                return false;
            }
            let he = ob.edges[*eh].half_edge;
            let twin = ob.half_edges[he].twin;
            if twin == popkern_core::topology::HalfEdgeHandle::default() {
                return false;
            }
            let v_start = ob.half_edges[he].origin;
            let v_end = ob.half_edges[twin].origin;
            let touches_boundary = |v: popkern_core::topology::VertexHandle| -> bool {
                ob.half_edges
                    .values()
                    .any(|h| h.origin == v && ob.boundary_edges.contains(&h.edge))
            };
            touches_boundary(v_start) || touches_boundary(v_end)
        })
        .map(|(h, _)| h)
        .collect();
    edges.sort();
    edges
}

fn pick_boundary_adjacent_interior_edge(ob: &BRep) -> popkern_core::topology::EdgeHandle {
    boundary_adjacent_interior_edges(ob)
        .into_iter()
        .next()
        .expect("open_box has a boundary-adjacent interior edge")
}

/// Fillet a single interior edge of `open_box` with a given
/// `CapKind`, asserting the output validates end-to-end.
fn fillet_open_box_with_cap(kind: CapKind) -> crate::BlendResult {
    let ob = BRep::open_box(tol()).unwrap();
    let edge = pick_boundary_adjacent_interior_edge(&ob);
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];
    let opts = FilletOptions {
        cap_kind: kind,
        ..FilletOptions::default()
    };
    let faces_before = ob.faces.len();
    let result = fillet_edges(&ob, &specs, &opts, tol(), None)
        .unwrap_or_else(|e| panic!("fillet on open_box (kind={kind:?}) failed: {e:?}"));
    assert!(
        result.brep.faces.len() > faces_before,
        "{kind:?} cap: face count should grow (before={}, after={})",
        faces_before,
        result.brep.faces.len()
    );
    result
        .brep
        .validate()
        .unwrap_or_else(|e| panic!("open_box fillet with {kind:?} cap should validate: {e:?}"));
    result
}

/// `CapKind::Flat` on open_box: cap insertion fires at the free
/// boundary endpoint, output validates, and at least one
/// CapApplied warning carries `topology_inserted = true`.
#[test]
fn fillet_open_box_interior_edge_with_flat_cap() {
    let result = fillet_open_box_with_cap(CapKind::Flat);
    let any_inserted = result.edge_warnings.iter().any(|(_, w)| {
        matches!(
            w,
            BlendWarning::CapApplied {
                topology_inserted: true,
                kind: CapKind::Flat,
                ..
            }
        )
    });
    assert!(
        any_inserted,
        "§8.3 end-to-end: expected at least one Flat cap to be inserted on open_box"
    );
}

/// Open cap on open_box: no insertion warnings, output validates.
/// Confirms that the boundary-edge split path produces valid
/// topology before cap insertion enters the picture.
#[test]
fn fillet_open_box_interior_edge_with_open_cap_validates() {
    let result = fillet_open_box_with_cap(CapKind::Open);
    let cap_warns = result
        .edge_warnings
        .iter()
        .filter(|(_, w)| matches!(w, BlendWarning::CapApplied { .. }))
        .count();
    assert_eq!(
        cap_warns, 0,
        "Open cap must not emit CapApplied warnings on open_box"
    );
}

/// `CapKind::Rounded` on open_box: the inserted cap face carries
/// a `Sphere` surface (matching §8.3's G1-with-blend spec).
#[test]
fn fillet_open_box_interior_edge_with_rounded_cap() {
    let result = fillet_open_box_with_cap(CapKind::Rounded);
    let inserted = result.edge_warnings.iter().any(|(_, w)| {
        matches!(
            w,
            BlendWarning::CapApplied {
                topology_inserted: true,
                kind: CapKind::Rounded,
                ..
            }
        )
    });
    assert!(inserted, "Rounded cap should be inserted on open_box");
    use popkern_core::geometry::surface::SurfaceGeometry;
    let sphere_count = result
        .brep
        .faces
        .values()
        .filter(|f| {
            f.surface
                .as_ref()
                .is_some_and(|sof| matches!(sof.geometry, SurfaceGeometry::Sphere { .. }))
        })
        .count();
    assert!(
        sphere_count >= 1,
        "at least one face should carry a Sphere surface (cap or vertex blend)"
    );
}

/// `CapKind::TangentMatch` on open_box: the inserted cap face
/// carries a `Nurbs` surface (bilinear with a spine-tangent bulge
/// so the (0.5, 0.5) normal stays non-degenerate).
#[test]
fn fillet_open_box_interior_edge_with_tangent_match_cap() {
    let result = fillet_open_box_with_cap(CapKind::TangentMatch);
    let inserted = result.edge_warnings.iter().any(|(_, w)| {
        matches!(
            w,
            BlendWarning::CapApplied {
                topology_inserted: true,
                kind: CapKind::TangentMatch,
                ..
            }
        )
    });
    assert!(inserted, "TangentMatch cap should be inserted on open_box");
}

// ---------------------------------------------------------------------------
// Multi-edge open-shell integration
// ---------------------------------------------------------------------------

/// Fillet the first N boundary-adjacent interior edges of open_box
/// in one batch call and assert the output validates.
fn fillet_open_box_multi(n_edges: usize, kind: CapKind) -> crate::BlendResult {
    let ob = BRep::open_box(tol()).unwrap();
    let edges = boundary_adjacent_interior_edges(&ob);
    assert!(
        edges.len() >= n_edges,
        "open_box has {} such edges",
        edges.len()
    );
    let specs: Vec<_> = edges
        .into_iter()
        .take(n_edges)
        .map(|edge| FilletSpec {
            edge,
            radius: FilletRadius::Constant(0.1),
        })
        .collect();
    let opts = FilletOptions {
        cap_kind: kind,
        ..FilletOptions::default()
    };
    let result = fillet_edges(&ob, &specs, &opts, tol(), None)
        .unwrap_or_else(|e| panic!("multi-edge fillet ({n_edges} edges) failed: {e:?}"));
    result
        .brep
        .validate()
        .unwrap_or_else(|e| panic!("multi-edge ({n_edges}, {kind:?}) validate: {e:?}"));
    result
}

fn count_inserted(result: &crate::BlendResult) -> usize {
    result
        .edge_warnings
        .iter()
        .filter(|(_, w)| {
            matches!(
                w,
                BlendWarning::CapApplied {
                    topology_inserted: true,
                    ..
                }
            )
        })
        .count()
}

/// Two boundary-adjacent interior edges filleted together: at
/// least two caps insert, output validates.
#[test]
fn fillet_open_box_two_interior_edges_with_flat_caps() {
    let result = fillet_open_box_multi(2, CapKind::Flat);
    let inserted = count_inserted(&result);
    assert!(
        inserted >= 2,
        "expected ≥2 inserted caps across 2 edges, got {inserted}"
    );
}

/// Same two-edge configuration with Rounded caps — sphere faces
/// appear in the output.
#[test]
fn fillet_open_box_two_interior_edges_with_rounded_caps() {
    let result = fillet_open_box_multi(2, CapKind::Rounded);
    use popkern_core::geometry::surface::SurfaceGeometry;
    let sphere_faces = result
        .brep
        .faces
        .values()
        .filter(|f| {
            f.surface
                .as_ref()
                .is_some_and(|sof| matches!(sof.geometry, SurfaceGeometry::Sphere { .. }))
        })
        .count();
    assert!(sphere_faces >= 2, "got {sphere_faces} sphere faces");
}

/// All boundary-adjacent interior edges filleted together — the
/// aggressive open-shell batch case where every boundary vertex
/// is touched by two filleted edges.
#[test]
fn fillet_open_box_all_boundary_edges() {
    let ob = BRep::open_box(tol()).unwrap();
    let n = boundary_adjacent_interior_edges(&ob).len();
    let result = fillet_open_box_multi(n, CapKind::Flat);
    let inserted = count_inserted(&result);
    assert!(
        inserted >= 1,
        "expected ≥1 inserted caps across {n} edges, got {inserted}"
    );
}
