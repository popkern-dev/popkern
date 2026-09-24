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

/// The edges filleted by the open-box fixtures: the first `n`
/// boundary-adjacent interior edges, at r = 0.1.
fn open_box_specs(ob: &BRep, n: usize) -> Vec<FilletSpec> {
    let edges = boundary_adjacent_interior_edges(ob);
    assert!(edges.len() >= n, "open_box has {} such edges", edges.len());
    edges
        .into_iter()
        .take(n)
        .map(|edge| FilletSpec {
            edge,
            radius: FilletRadius::Constant(0.1),
        })
        .collect()
}

fn fillet_open_box(n_edges: usize, kind: CapKind) -> Result<crate::BlendResult, crate::BlendError> {
    let ob = BRep::open_box(tol()).unwrap();
    let opts = FilletOptions {
        cap_kind: kind,
        ..FilletOptions::default()
    };
    fillet_edges(&ob, &open_box_specs(&ob, n_edges), &opts, tol(), None)
}

/// Open cap on open_box: the strip runs out at the sheet boundary and
/// its free end settles onto the strip's section by the rim plane
/// (PK-530). The rim edges split at the fillet's anchors keep their
/// curves. The whole body holds check 9.
#[test]
fn fillet_open_box_interior_edge_with_open_cap_validates() {
    let result = fillet_open_box(1, CapKind::Open).expect("open-cap fillet serves");
    let cap_warns = result
        .edge_warnings
        .iter()
        .filter(|(_, w)| matches!(w, BlendWarning::CapApplied { .. }))
        .count();
    assert_eq!(
        cap_warns, 0,
        "Open cap must not emit CapApplied warnings on open_box"
    );
    let errors = result.brep.validate_all();
    assert!(errors.is_empty(), "validate_all: {errors:?}");
    assert!(
        result.brep.edges.values().all(|e| e.curve.is_some()),
        "every edge carries a curve"
    );
}

/// Inserted caps on an open sheet (Flat, Rounded, TangentMatch) do
/// not yet author their strip|cap edge. `insert_cap_face` splits the
/// loop with an edge that has no curve, and check 9 had silently
/// skipped it while `validate()` passed. The territory gate (PK-530)
/// now declines them typed, with the input handed back. PK-532 builds
/// the exact cap boundaries and flips these pins back to serving.
#[test]
fn inserted_open_sheet_caps_decline_until_their_boundaries_are_exact() {
    let ob = BRep::open_box(tol()).unwrap();
    let all = boundary_adjacent_interior_edges(&ob).len();
    for (kind, n) in [
        (CapKind::Flat, 1),
        (CapKind::Rounded, 1),
        (CapKind::TangentMatch, 1),
        (CapKind::Flat, 2),
        (CapKind::Rounded, 2),
        (CapKind::Flat, all),
    ] {
        let Err(crate::BlendError::PartialFailure {
            partial_result,
            failures,
            ..
        }) = fillet_open_box(n, kind)
        else {
            panic!("{kind:?} x{n}: expected the typed territory decline");
        };
        assert_eq!(
            partial_result.faces.len(),
            ob.faces.len(),
            "input handed back"
        );
        assert!(
            failures.iter().any(|(_, f)| matches!(
                f,
                crate::BlendFailure::InconsistentGeometry { reason, .. }
                    if reason.contains("without a curve")
            )),
            "{kind:?} x{n}: {failures:?}"
        );
    }
}
