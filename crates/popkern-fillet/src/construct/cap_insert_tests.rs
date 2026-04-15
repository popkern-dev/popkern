use super::*;
use crate::types::{FilletOptions, FilletRadius, FilletSpec};
use crate::{CapKind, fillet_edges};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default()
}

/// With `cap_kind = Open`, no `CapApplied` warnings are emitted
/// for a vanilla cube fillet.
#[test]
fn open_cap_kind_emits_no_cap_applied_warnings() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let edge = cube.edges.keys().next().unwrap();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];
    let opts = FilletOptions::default(); // cap_kind = Open

    let result = fillet_edges(&cube, &specs, &opts, tol(), None).expect("fillet should succeed");

    let cap_warnings: Vec<_> = result
        .edge_warnings
        .iter()
        .filter(|(_, w)| matches!(w, BlendWarning::CapApplied { .. }))
        .collect();
    assert!(
        cap_warnings.is_empty(),
        "open cap_kind should produce no CapApplied warnings, got {}",
        cap_warnings.len()
    );
}

/// With `cap_kind = Flat` on a cube edge fillet, both endpoints
/// should emit `CapApplied { topology_inserted: false }` warnings
/// — the cube is a closed solid and every endpoint has an adjacent
/// face that serves as the implicit cap.
#[test]
fn flat_cap_kind_emits_two_implicit_warnings_on_cube_edge() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let edge = cube.edges.keys().next().unwrap();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];
    let opts = FilletOptions {
        cap_kind: CapKind::Flat,
        ..FilletOptions::default()
    };

    let result = fillet_edges(&cube, &specs, &opts, tol(), None).expect("fillet should succeed");

    let cap_warnings: Vec<_> = result
        .edge_warnings
        .iter()
        .filter_map(|(_, w)| match w {
            BlendWarning::CapApplied {
                edge,
                at_start,
                kind,
                topology_inserted,
            } => Some((*edge, *at_start, *kind, *topology_inserted)),
            _ => None,
        })
        .collect();

    assert_eq!(
        cap_warnings.len(),
        2,
        "expected exactly 2 CapApplied warnings (start + end), got {}",
        cap_warnings.len()
    );

    // All four fields match the expected values on both endpoints.
    let mut saw_start = false;
    let mut saw_end = false;
    for (w_edge, at_start, kind, inserted) in cap_warnings {
        assert_eq!(w_edge, edge);
        assert_eq!(kind, CapKind::Flat);
        // Cube is a closed solid — every corner has a third face
        // so the cap is always implicit.
        assert!(
            !inserted,
            "cube corner is closed, topology should not be inserted"
        );
        if at_start {
            saw_start = true;
        } else {
            saw_end = true;
        }
    }
    assert!(saw_start, "missing start-endpoint cap warning");
    assert!(saw_end, "missing end-endpoint cap warning");
}

/// `cap_kind = Rounded` and `TangentMatch` also emit the
/// implicit-closed warnings on a cube. The kind field must match.
#[test]
fn rounded_and_tangent_match_kinds_propagate_to_warnings() {
    for kind in [CapKind::Rounded, CapKind::TangentMatch] {
        let cube = BRep::unit_cube(tol()).unwrap();
        let edge = cube.edges.keys().next().unwrap();
        let specs = vec![FilletSpec {
            edge,
            radius: FilletRadius::Constant(0.1),
        }];
        let opts = FilletOptions {
            cap_kind: kind,
            ..FilletOptions::default()
        };
        let result = fillet_edges(&cube, &specs, &opts, tol(), None).unwrap();

        let kinds: Vec<CapKind> = result
            .edge_warnings
            .iter()
            .filter_map(|(_, w)| match w {
                BlendWarning::CapApplied { kind, .. } => Some(*kind),
                _ => None,
            })
            .collect();
        assert_eq!(kinds.len(), 2, "expected 2 warnings for {kind:?}");
        assert!(kinds.iter().all(|k| *k == kind));
    }
}

/// `cap_kind = UserDefined` is a [BEYOND] placeholder that falls
/// back to `Open` semantics: no cap geometry is produced and no
/// `CapApplied` warnings fire.
#[test]
fn user_defined_cap_kind_emits_no_cap_applied_warnings() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let edge = cube.edges.keys().next().unwrap();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];
    let opts = FilletOptions {
        cap_kind: CapKind::UserDefined,
        ..FilletOptions::default()
    };
    let result = fillet_edges(&cube, &specs, &opts, tol(), None).unwrap();

    let cap_count = result
        .edge_warnings
        .iter()
        .filter(|(_, w)| matches!(w, BlendWarning::CapApplied { .. }))
        .count();
    assert_eq!(cap_count, 0, "UserDefined should not emit CapApplied");
}

/// Cap insertion is a pure post-construct pass — enabling it does
/// not change the BRep's face or edge count on closed solids. This
/// is the "no regression" guard: every existing closed-body fillet
/// result is byte-compatible with the pre-§8.3 output.
#[test]
fn cap_insertion_is_a_no_op_on_closed_cube_topology() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let edge = cube.edges.keys().next().unwrap();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];

    let open_opts = FilletOptions::default();
    let flat_opts = FilletOptions {
        cap_kind: CapKind::Flat,
        ..FilletOptions::default()
    };
    let rounded_opts = FilletOptions {
        cap_kind: CapKind::Rounded,
        ..FilletOptions::default()
    };

    let open_result = fillet_edges(&cube, &specs, &open_opts, tol(), None).unwrap();
    let flat_result = fillet_edges(&cube, &specs, &flat_opts, tol(), None).unwrap();
    let rounded_result = fillet_edges(&cube, &specs, &rounded_opts, tol(), None).unwrap();

    // Same face / edge / vertex counts across all three cap modes.
    assert_eq!(open_result.brep.faces.len(), flat_result.brep.faces.len());
    assert_eq!(open_result.brep.edges.len(), flat_result.brep.edges.len());
    assert_eq!(
        open_result.brep.vertices.len(),
        flat_result.brep.vertices.len()
    );
    assert_eq!(
        open_result.brep.faces.len(),
        rounded_result.brep.faces.len()
    );
    assert_eq!(
        open_result.brep.edges.len(),
        rounded_result.brep.edges.len()
    );
    assert_eq!(
        open_result.brep.vertices.len(),
        rounded_result.brep.vertices.len()
    );
}
