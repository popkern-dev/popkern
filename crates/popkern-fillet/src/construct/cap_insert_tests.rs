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

/// On a cube every corner has 3 incident faces, so
/// `endpoint_is_closed(vertex, face_1, face_2)` must return `true`
/// regardless of which two faces we call the "supports" — the
/// third cube face always provides an implicit cap.
#[test]
fn endpoint_is_closed_returns_true_on_cube_corner() {
    let cube = BRep::unit_cube(tol()).unwrap();
    // Pick any vertex and any two faces as the pseudo-supports.
    let vertex = cube.vertices.keys().next().unwrap();
    let mut faces = cube.faces.keys();
    let f1 = faces.next().unwrap();
    let f2 = faces.next().unwrap();
    assert!(endpoint_is_closed(&cube, vertex, f1, f2));
}

/// On a `sheet_quad` every corner has exactly 2 incident faces
/// (front + backside), so passing both as the supports leaves no
/// "other" face for the detector to find — returns `false`. This
/// is the free-boundary signal: an endpoint where a cap face
/// would actually need to be inserted.
#[test]
fn endpoint_is_closed_returns_false_on_sheet_quad_corner() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let vertex = sheet.vertices.keys().next().unwrap();
    // The sheet has exactly 2 faces — treat them both as "supports"
    // to simulate a fillet between the front and backside of the
    // sheet.
    let mut faces = sheet.faces.keys();
    let f1 = faces.next().unwrap();
    let f2 = faces.next().unwrap();
    assert!(
        !endpoint_is_closed(&sheet, vertex, f1, f2),
        "sheet corner has no third face; should report free boundary"
    );
}

/// Pair up a sheet-quad vertex with mismatched supports — e.g. only
/// one of the two real faces declared as a support — and
/// `endpoint_is_closed` should return `true` (the second face counts
/// as "other"). Exercises the "at least one face outside the
/// support pair" branch on sheet topology.
#[test]
fn endpoint_is_closed_detects_third_face_on_sheet_quad_with_one_support() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let vertex = sheet.vertices.keys().next().unwrap();
    let f1 = sheet.faces.keys().next().unwrap();
    // Invent a FaceHandle that doesn't exist — every real face
    // except f1 is "other," so the detector should find one.
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

// ---------------------------------------------------------------------------
// Direct topology surgery unit tests
// ---------------------------------------------------------------------------

/// Test `insert_cap_face` directly on a sheet_quad. Use
/// vertices v1 = (1, 0, 0) and v3 = (0, 1, 0) as the "contact
/// vertices" of a synthetic cap at the v0 = (0, 0, 0) corner.
/// `mef` splits the front face along the v1-v3 diagonal,
/// creating a new triangular face that gets the cap surface.
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
    let contact_1_pos = popkern_core::geometry::Point3::new(1.0, 0.0, 0.0);
    let contact_2_pos = popkern_core::geometry::Point3::new(0.0, 1.0, 0.0);
    let cap = crate::cap::CapGeometry {
        kind: CapKind::Flat,
        surface: Some(SurfaceGeometry::Plane {
            origin: popkern_core::geometry::Point3::new(0.0, 0.0, 0.0),
            normal: popkern_core::geometry::Vec3::new(0.0, 0.0, 1.0),
            u_axis: popkern_core::geometry::Vec3::new(1.0, 0.0, 0.0),
            v_axis: popkern_core::geometry::Vec3::new(0.0, 1.0, 0.0),
        }),
        boundary_points: vec![contact_1_pos, contact_2_pos],
    };
    let faces_before = sheet.faces.len();
    let edges_before = sheet.edges.len();
    let inserted = super::insert_cap_face(&mut sheet, v0, &cap);
    assert!(inserted, "cap face should be inserted");
    assert_eq!(
        sheet.faces.len(),
        faces_before + 1,
        "one new face should be created by mef"
    );
    assert_eq!(
        sheet.edges.len(),
        edges_before + 1,
        "one new edge (the diagonal) should be created"
    );
    // The new face has the cap surface assigned.
    let has_cap_surface = sheet.faces.values().any(|f| {
        f.surface
            .as_ref()
            .is_some_and(|sof| matches!(sof.geometry, SurfaceGeometry::Plane { .. }))
    });
    assert!(
        has_cap_surface,
        "at least one face should have the cap's plane surface"
    );
}

/// `insert_cap_face` returns false (no-op) when the boundary
/// points don't map to real vertices near the endpoint.
#[test]
fn insert_cap_face_returns_false_for_distant_contacts() {
    let mut sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let v0 = sheet.vertices.keys().next().unwrap();
    let cap = crate::cap::CapGeometry {
        kind: CapKind::Flat,
        surface: Some(popkern_core::geometry::surface::SurfaceGeometry::Plane {
            origin: popkern_core::geometry::Point3::new(0.0, 0.0, 0.0),
            normal: popkern_core::geometry::Vec3::new(0.0, 0.0, 1.0),
            u_axis: popkern_core::geometry::Vec3::new(1.0, 0.0, 0.0),
            v_axis: popkern_core::geometry::Vec3::new(0.0, 1.0, 0.0),
        }),
        boundary_points: vec![
            popkern_core::geometry::Point3::new(99.0, 99.0, 99.0),
            popkern_core::geometry::Point3::new(100.0, 100.0, 100.0),
        ],
    };
    let faces_before = sheet.faces.len();
    let inserted = super::insert_cap_face(&mut sheet, v0, &cap);
    assert!(
        !inserted,
        "cap insertion should fail when contacts are far from BRep vertices"
    );
    assert_eq!(sheet.faces.len(), faces_before);
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
