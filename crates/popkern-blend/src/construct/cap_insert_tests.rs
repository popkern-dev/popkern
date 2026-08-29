use super::*;
use crate::types::{FilletOptions, FilletRadius, FilletSpec};
use crate::{CapKind, fillet_edges};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default()
}

// ---------------------------------------------------------------------------
// Free-boundary detection
// ---------------------------------------------------------------------------

/// Cube corner has a 3rd face, so cap is implicit.
#[test]
fn endpoint_is_closed_returns_true_on_cube_corner() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    let mut faces = cube.faces.keys();
    let f1 = faces.next().unwrap();
    let f2 = faces.next().unwrap();
    assert!(endpoint_is_closed(&cube, vertex, f1, f2));
}

/// Sheet-quad corner has only 2 faces (front + backside). Passing
/// both as supports leaves no "other" — returns `false`.
#[test]
fn endpoint_is_closed_returns_false_on_sheet_quad_corner() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let vertex = sheet.vertices.keys().next().unwrap();
    let mut faces = sheet.faces.keys();
    let f1 = faces.next().unwrap();
    let f2 = faces.next().unwrap();
    assert!(!endpoint_is_closed(&sheet, vertex, f1, f2));
}

/// Sheet-quad vertex with a bogus second support: the real second
/// face counts as "other" — returns `true`.
#[test]
fn endpoint_is_closed_detects_third_face_on_sheet_quad_with_one_support() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let vertex = sheet.vertices.keys().next().unwrap();
    let f1 = sheet.faces.keys().next().unwrap();
    let fake = popkern_core::topology::FaceHandle::default();
    assert!(endpoint_is_closed(&sheet, vertex, f1, fake));
}

// ---------------------------------------------------------------------------
// Sheet-body fillet — current-state spike
// ---------------------------------------------------------------------------

/// Filleting an edge of a planar sheet body is blocked by the
/// current fillet analyzer (the front and backside faces of a
/// sheet are coincident planes with opposite normals, which the
/// convexity classifier treats as degenerate). This test pins the
/// behaviour: running `fillet_edges` on a `sheet_quad` edge either
/// fails cleanly at analysis or returns a result with no new faces,
/// and critically **does not panic**.
///
/// When §5.7 2D fillets / sheet-body fillets are implemented,
/// update this test to assert success plus the presence of the
/// cap face — that's where cap topology insertion actually gets
/// exercised end-to-end.
#[test]
fn sheet_quad_fillet_is_currently_blocked_at_analysis() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let edge = sheet.edges.keys().next().unwrap();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];
    let opts = FilletOptions {
        cap_kind: CapKind::Flat,
        ..FilletOptions::default()
    };
    // Any outcome except a panic is acceptable — this test is a
    // regression guard against the pipeline crashing on unusual
    // topology, not a correctness claim about sheet filleting.
    let _ = fillet_edges(&sheet, &specs, &opts, tol(), None);
}

// ---------------------------------------------------------------------------
// Cap-warning shape on closed solids (cube)
// ---------------------------------------------------------------------------

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
    let opts = FilletOptions::default();
    let result = fillet_edges(&cube, &specs, &opts, tol(), None).unwrap();
    let count = result
        .edge_warnings
        .iter()
        .filter(|(_, w)| matches!(w, BlendWarning::CapApplied { .. }))
        .count();
    assert_eq!(count, 0);
}

/// Cube edge with `Flat` cap: both endpoints emit implicit
/// (topology_inserted: false) warnings.
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
    let result = fillet_edges(&cube, &specs, &opts, tol(), None).unwrap();
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
    assert_eq!(cap_warnings.len(), 2);
    let mut saw_start = false;
    let mut saw_end = false;
    for (w_edge, at_start, kind, inserted) in cap_warnings {
        assert_eq!(w_edge, edge);
        assert_eq!(kind, CapKind::Flat);
        assert!(!inserted);
        if at_start {
            saw_start = true;
        } else {
            saw_end = true;
        }
    }
    assert!(saw_start && saw_end);
}

/// `Rounded` and `TangentMatch` on a cube emit the same shape of
/// warnings — kind propagates through.
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

/// `UserDefined` is a [BEYOND] placeholder — falls back to Open
/// semantics, no cap warnings.
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
    let count = result
        .edge_warnings
        .iter()
        .filter(|(_, w)| matches!(w, BlendWarning::CapApplied { .. }))
        .count();
    assert_eq!(count, 0);
}

// ---------------------------------------------------------------------------
// Direct topology surgery unit tests
// ---------------------------------------------------------------------------

/// `insert_cap_face` on a sheet_quad creates one new face + one
/// new diagonal edge; the new face carries the cap's plane.
#[test]
fn insert_cap_face_creates_new_face_on_sheet_quad() {
    use popkern_core::geometry::surface::SurfaceGeometry;
    let mut sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let v0 = sheet
        .vertices
        .iter()
        .find(|(_, v)| v.position.x == 0.0 && v.position.y == 0.0)
        .map(|(h, _)| h)
        .unwrap();
    let cap = crate::cap::CapGeometry {
        kind: CapKind::Flat,
        surface: Some(SurfaceGeometry::Plane {
            origin: popkern_core::geometry::Point3::new(0.0, 0.0, 0.0),
            normal: popkern_core::geometry::Vec3::new(0.0, 0.0, 1.0),
            u_axis: popkern_core::geometry::Vec3::new(1.0, 0.0, 0.0),
            v_axis: popkern_core::geometry::Vec3::new(0.0, 1.0, 0.0),
        }),
        boundary_points: vec![
            popkern_core::geometry::Point3::new(1.0, 0.0, 0.0),
            popkern_core::geometry::Point3::new(0.0, 1.0, 0.0),
        ],
    };
    let faces_before = sheet.faces.len();
    let edges_before = sheet.edges.len();
    assert!(super::insert_cap_face(&mut sheet, v0, &cap));
    assert_eq!(sheet.faces.len(), faces_before + 1);
    assert_eq!(sheet.edges.len(), edges_before + 1);
    assert!(sheet.faces.values().any(|f| {
        f.surface
            .as_ref()
            .is_some_and(|sof| matches!(sof.geometry, SurfaceGeometry::Plane { .. }))
    }));
}

/// `insert_cap_face` returns false when no surrounding face has
/// surface geometry.
#[test]
fn insert_cap_face_returns_false_when_all_faces_surfaceless() {
    let mut sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    for fh in sheet.faces.keys().collect::<Vec<_>>() {
        sheet.faces[fh].surface = None;
    }
    let v0 = sheet.vertices.keys().next().unwrap();
    let cap = crate::cap::CapGeometry {
        kind: CapKind::Flat,
        surface: Some(popkern_core::geometry::surface::SurfaceGeometry::Plane {
            origin: popkern_core::geometry::Point3::new(0.0, 0.0, 0.0),
            normal: popkern_core::geometry::Vec3::new(0.0, 0.0, 1.0),
            u_axis: popkern_core::geometry::Vec3::new(1.0, 0.0, 0.0),
            v_axis: popkern_core::geometry::Vec3::new(0.0, 1.0, 0.0),
        }),
        boundary_points: vec![],
    };
    let faces_before = sheet.faces.len();
    assert!(!super::insert_cap_face(&mut sheet, v0, &cap));
    assert_eq!(sheet.faces.len(), faces_before);
}

/// Closed-solid cube fillet: cap insertion is a no-op — face /
/// edge / vertex counts match across all cap modes.
#[test]
fn cap_insertion_is_a_no_op_on_closed_cube_topology() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let edge = cube.edges.keys().next().unwrap();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];
    let with_kind = |k: CapKind| {
        let opts = FilletOptions {
            cap_kind: k,
            ..FilletOptions::default()
        };
        fillet_edges(&cube, &specs, &opts, tol(), None).unwrap()
    };
    let open_r = with_kind(CapKind::Open);
    for k in [CapKind::Flat, CapKind::Rounded] {
        let r = with_kind(k);
        assert_eq!(open_r.brep.faces.len(), r.brep.faces.len());
        assert_eq!(open_r.brep.edges.len(), r.brep.edges.len());
        assert_eq!(open_r.brep.vertices.len(), r.brep.vertices.len());
    }
}
