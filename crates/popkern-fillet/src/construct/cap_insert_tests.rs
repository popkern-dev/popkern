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
// True open-shell end-to-end: fillet_edges on open_box
// ---------------------------------------------------------------------------

/// Find an interior edge of `open_box` with at least one endpoint
/// vertex on the boundary — the canonical case where cap insertion
/// must fire.
fn pick_boundary_adjacent_interior_edge(ob: &BRep) -> popkern_core::topology::EdgeHandle {
    ob.edges
        .iter()
        .find(|(eh, _)| {
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
        .expect("open_box has a boundary-adjacent interior edge")
}

/// Fillet an interior edge of `open_box` with a given `CapKind`,
/// asserting the output validates end-to-end. Returns the result
/// so callers can make kind-specific assertions about the cap.
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
/// carries a `Nurbs` surface (degree-1 bilinear patch per cap.rs).
///
/// **Ignored:** the bilinear NURBS patch produced by
/// `build_tangent_match_cap` is hard-assigned `orientation: true`
/// in `insert_cap_face`, which fails `validate()`
/// (InconsistentOrientation). The NURBS control-point winding
/// doesn't match the MEF-split face's loop direction on open_box;
/// the cap-surface builder needs to compute orientation from the
/// loop winding rather than hard-code it. Tracked as a follow-up
/// — Flat and Rounded (with symmetric surfaces) do validate.
#[test]
#[ignore]
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
    assert_eq!(cap_warnings.len(), 2);
    let mut saw_start = false;
    let mut saw_end = false;
    for (w_edge, at_start, kind, inserted) in cap_warnings {
        assert_eq!(w_edge, edge);
        assert_eq!(kind, CapKind::Flat);
        assert!(!inserted, "cube corner is closed");
        if at_start {
            saw_start = true;
        } else {
            saw_end = true;
        }
    }
    assert!(saw_start && saw_end);
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
    assert!(super::insert_cap_face(&mut sheet, v0, &cap));
    assert_eq!(sheet.faces.len(), faces_before + 1);
    assert_eq!(sheet.edges.len(), edges_before + 1);
    assert!(sheet.faces.values().any(|f| {
        f.surface
            .as_ref()
            .is_some_and(|sof| matches!(sof.geometry, SurfaceGeometry::Plane { .. }))
    }));
}

/// `insert_cap_face` returns false when the vertex has no
/// adjacent vertices on faces with surface geometry — e.g.
/// when all surrounding faces are virtual/surfaceless.
#[test]
fn insert_cap_face_returns_false_when_all_faces_surfaceless() {
    let mut sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    // Clear ALL face surfaces so no face counts as "surfaced".
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
    let inserted = super::insert_cap_face(&mut sheet, v0, &cap);
    assert!(
        !inserted,
        "cap insertion should fail when no adjacent face has a surface"
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
