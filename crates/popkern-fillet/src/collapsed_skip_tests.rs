//! Tests for the §13.1(5) collapsed-edge soft-skip pre-pass.

use super::*;
use crate::types::FilletRadius;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default()
}

// ---------------------------------------------------------------------------
// negative path: healthy cube
// ---------------------------------------------------------------------------

#[test]
fn skip_collapsed_edges_is_noop_on_healthy_cube() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let specs: Vec<FilletSpec> = brep
        .edges
        .keys()
        .map(|e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(0.1),
        })
        .collect();

    let (kept, warnings) = skip_collapsed_edges(&brep, &specs, &tol());
    assert_eq!(
        kept.len(),
        specs.len(),
        "healthy cube should not drop any specs"
    );
    assert!(warnings.is_empty(), "healthy cube should emit no warnings");
    for (orig, new) in specs.iter().zip(kept.iter()) {
        assert_eq!(orig.edge, new.edge, "spec order should be preserved");
    }
}

// ---------------------------------------------------------------------------
// positive path: synthetic collapsed edge
// ---------------------------------------------------------------------------

#[test]
fn skip_collapsed_edges_drops_collapsed_spec_and_emits_warning() {
    // Build a cube, collapse one edge by moving its start vertex
    // to coincide with its end vertex, then run the pre-pass on
    // the full spec list. The collapsed edge should be dropped
    // and a CollapsedEdgeSkipped warning emitted, while the
    // other 11 specs should pass through unchanged.
    let mut brep = BRep::unit_cube(tol()).unwrap();
    let target_edge = brep.edges.keys().next().unwrap();
    let he = brep.edges[target_edge].half_edge;
    let twin = brep.half_edges[he].twin;
    let v_start = brep.half_edges[he].origin;
    let v_end = brep.half_edges[twin].origin;
    brep.vertices[v_start].position = brep.vertices[v_end].position;

    let specs: Vec<FilletSpec> = brep
        .edges
        .keys()
        .map(|e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(0.1),
        })
        .collect();
    let original_count = specs.len();

    let (kept, warnings) = skip_collapsed_edges(&brep, &specs, &tol());

    assert_eq!(
        kept.len(),
        original_count - 1,
        "expected one spec dropped from the batch"
    );
    assert_eq!(warnings.len(), 1, "expected exactly one warning");
    match &warnings[0] {
        BlendWarning::CollapsedEdgeSkipped { edge, length } => {
            assert_eq!(*edge, target_edge);
            assert!(
                *length < tol().linear(),
                "reported length {length} should be below tolerance"
            );
        }
        other => panic!("expected CollapsedEdgeSkipped, got {other:?}"),
    }
    // The kept specs should NOT contain the collapsed edge.
    assert!(
        !kept.iter().any(|s| s.edge == target_edge),
        "collapsed edge should not appear in kept spec list"
    );
}

// ---------------------------------------------------------------------------
// edge case: closed edge (v_start == v_end) is NOT collapsed
// ---------------------------------------------------------------------------

#[test]
fn skip_collapsed_edges_preserves_closed_intentional_edges() {
    // A closed edge (v_start == v_end) is the intentional
    // closed-rim case used by cylinder seams and closed-fillet
    // annuli. It must NOT be dropped by the collapsed-edge
    // pre-pass even though its endpoint positions are identical.
    //
    // Synthesise this by overriding the edge's half-edge chain
    // so both origins point at the same vertex — simpler than
    // constructing a real closed-rim BRep and sufficient to
    // exercise the short-circuit branch.
    //
    // NB: this leaves the cube's half-edge chain structurally
    // invalid, but the pre-pass only reads positions, not the
    // loop chain, so the invariant-breaking is local to the test.
    let mut brep = BRep::unit_cube(tol()).unwrap();
    let target_edge = brep.edges.keys().next().unwrap();
    let he = brep.edges[target_edge].half_edge;
    let twin = brep.half_edges[he].twin;
    let v_start = brep.half_edges[he].origin;
    // Override the twin's origin to match the first half-edge's
    // origin — simulates a closed edge.
    brep.half_edges[twin].origin = v_start;

    let specs = vec![FilletSpec {
        edge: target_edge,
        radius: FilletRadius::Constant(0.1),
    }];

    let (kept, warnings) = skip_collapsed_edges(&brep, &specs, &tol());
    assert_eq!(
        kept.len(),
        1,
        "closed edge should NOT be dropped as collapsed"
    );
    assert_eq!(warnings.len(), 0, "closed edge should emit no warning");
}
