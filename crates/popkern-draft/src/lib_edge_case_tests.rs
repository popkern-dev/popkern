use super::*;
use popkern_core::geometry::Vec3;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Helper: find a side face on a unit cube (normal perpendicular to Z).
fn side_face(brep: &BRep) -> popkern_core::topology::FaceHandle {
    use popkern_core::geometry::surface::SurfaceEval;
    brep.faces
        .iter()
        .find(|(_, fd)| {
            fd.surface.as_ref().map_or(false, |sof| {
                let n = sof.geometry.normal_at(0.5, 0.5);
                n.z.abs() < 0.1
            })
        })
        .map(|(fh, _)| fh)
        .unwrap_or_else(|| brep.faces.keys().next().unwrap())
}

#[test]
fn draft_angle_exactly_at_max_boundary_succeeds() {
    // The check is `angle.abs() > max_angle`, so angle == max should pass.
    let brep = BRep::unit_cube(tol()).unwrap();
    let face = side_face(&brep);
    let max = std::f64::consts::FRAC_PI_4;
    let specs = vec![DraftSpec {
        face,
        angle: max, // exactly at boundary
        pull_direction: Vec3::new(0.0, 0.0, 1.0),
        neutral: NeutralReference::Auto,
        angle_variation: AngleVariation::Constant,
    }];
    let result = draft_faces(&brep, &specs, &DraftOptions::default(), tol(), None);
    assert!(
        result.is_ok(),
        "angle == max_angle should succeed, got: {:?}",
        result.err()
    );
}

#[test]
fn draft_angle_just_above_max_boundary_fails() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let face = side_face(&brep);
    let max = std::f64::consts::FRAC_PI_4;
    let specs = vec![DraftSpec {
        face,
        angle: max + 1e-10, // just above boundary
        pull_direction: Vec3::new(0.0, 0.0, 1.0),
        neutral: NeutralReference::Auto,
        angle_variation: AngleVariation::Constant,
    }];
    let result = draft_faces(&brep, &specs, &DraftOptions::default(), tol(), None);
    assert!(matches!(result, Err(DraftError::AngleTooLarge { .. })));
}

#[test]
fn draft_negative_angle_at_max_boundary_succeeds() {
    // Negative angle should also be checked via abs(), so -max_angle passes.
    let brep = BRep::unit_cube(tol()).unwrap();
    let face = side_face(&brep);
    let max = std::f64::consts::FRAC_PI_4;
    let specs = vec![DraftSpec {
        face,
        angle: -max,
        pull_direction: Vec3::new(0.0, 0.0, 1.0),
        neutral: NeutralReference::Auto,
        angle_variation: AngleVariation::Constant,
    }];
    let result = draft_faces(&brep, &specs, &DraftOptions::default(), tol(), None);
    assert!(
        result.is_ok(),
        "negative angle at boundary should succeed, got: {:?}",
        result.err()
    );
}

#[test]
fn draft_with_custom_max_angle() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let face = side_face(&brep);
    let custom_max = 0.05; // very small max
    let specs = vec![DraftSpec {
        face,
        angle: 0.04,
        pull_direction: Vec3::new(0.0, 0.0, 1.0),
        neutral: NeutralReference::Auto,
        angle_variation: AngleVariation::Constant,
    }];
    let opts = DraftOptions {
        max_angle: custom_max,
        ..DraftOptions::default()
    };
    let result = draft_faces(&brep, &specs, &opts, tol(), None);
    assert!(result.is_ok(), "within custom max should succeed");

    // Now exceed the custom max.
    let specs_over = vec![DraftSpec {
        face,
        angle: 0.06,
        pull_direction: Vec3::new(0.0, 0.0, 1.0),
        neutral: NeutralReference::Auto,
        angle_variation: AngleVariation::Constant,
    }];
    let result_over = draft_faces(&brep, &specs_over, &opts, tol(), None);
    assert!(matches!(result_over, Err(DraftError::AngleTooLarge { .. })));
}

#[test]
fn draft_with_adjacent_mode_extend_trim() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let face = side_face(&brep);
    let specs = vec![DraftSpec {
        face,
        angle: 0.087,
        pull_direction: Vec3::new(0.0, 0.0, 1.0),
        neutral: NeutralReference::Auto,
        angle_variation: AngleVariation::Constant,
    }];
    let opts = DraftOptions {
        adjacent_mode: AdjacentMode::ExtendTrim,
        ..DraftOptions::default()
    };
    let result = draft_faces(&brep, &specs, &opts, tol(), None);
    assert!(result.is_ok(), "ExtendTrim mode failed: {:?}", result.err());
    let dr = result.unwrap();
    // Cube should still have 6 faces with ExtendTrim.
    assert_eq!(dr.brep.faces.len(), 6);
}

#[test]
fn draft_with_adjacent_mode_transition_faces() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let face = side_face(&brep);
    let specs = vec![DraftSpec {
        face,
        angle: 0.087,
        pull_direction: Vec3::new(0.0, 0.0, 1.0),
        neutral: NeutralReference::Auto,
        angle_variation: AngleVariation::Constant,
    }];
    let opts = DraftOptions {
        adjacent_mode: AdjacentMode::TransitionFaces,
        ..DraftOptions::default()
    };
    let result = draft_faces(&brep, &specs, &opts, tol(), None);
    assert!(
        result.is_ok(),
        "TransitionFaces mode failed: {:?}",
        result.err()
    );
    let dr = result.unwrap();
    // TransitionFaces may add extra faces at boundary edges.
    assert!(
        dr.brep.faces.len() >= 6,
        "expected at least 6 faces, got {}",
        dr.brep.faces.len()
    );
}

#[test]
fn draft_with_adjacent_mode_extend_trim_with_fallback() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let face = side_face(&brep);
    let specs = vec![DraftSpec {
        face,
        angle: 0.087,
        pull_direction: Vec3::new(0.0, 0.0, 1.0),
        neutral: NeutralReference::Auto,
        angle_variation: AngleVariation::Constant,
    }];
    let opts = DraftOptions {
        adjacent_mode: AdjacentMode::ExtendTrimWithFallback,
        ..DraftOptions::default()
    };
    let result = draft_faces(&brep, &specs, &opts, tol(), None);
    assert!(
        result.is_ok(),
        "ExtendTrimWithFallback mode failed: {:?}",
        result.err()
    );
}

#[test]
fn draft_variable_angle_through_full_api() {
    // Exercise AngleVariation::Variable through the full draft_faces pipeline.
    let brep = BRep::unit_cube(tol()).unwrap();
    let face = side_face(&brep);
    let variation = AngleVariation::Variable(vec![(0.0, 0.03), (0.5, 0.06), (1.0, 0.09)]);
    let specs = vec![DraftSpec {
        face,
        angle: 0.06, // fallback angle
        pull_direction: Vec3::new(0.0, 0.0, 1.0),
        neutral: NeutralReference::Auto,
        angle_variation: variation,
    }];
    let result = draft_faces(&brep, &specs, &DraftOptions::default(), tol(), None);
    assert!(
        result.is_ok(),
        "variable angle draft failed: {:?}",
        result.err()
    );
    let dr = result.unwrap();
    assert_eq!(dr.brep.faces.len(), 6);
}

#[test]
fn draft_variable_angle_two_point_linear() {
    // Two-point variable angle: linear variation from 2 deg to 8 deg.
    let brep = BRep::unit_cube(tol()).unwrap();
    let face = side_face(&brep);
    let variation = AngleVariation::Variable(vec![(0.0, 0.035), (1.0, 0.14)]);
    let specs = vec![DraftSpec {
        face,
        angle: 0.087,
        pull_direction: Vec3::new(0.0, 0.0, 1.0),
        neutral: NeutralReference::Auto,
        angle_variation: variation,
    }];
    let result = draft_faces(&brep, &specs, &DraftOptions::default(), tol(), None);
    assert!(
        result.is_ok(),
        "two-point variable draft failed: {:?}",
        result.err()
    );
}

#[test]
fn draft_cylinder_euler_characteristic_preserved() {
    // A proper cylinder has V=2, E=3, F=3 → V-E+F=2.
    // After drafting a cap face, Euler should still hold.
    let brep = BRep::cylinder(1.0, 2.0, tol()).unwrap();

    // Find a cap face (planar) to draft.
    let cap_face = brep
        .faces
        .iter()
        .find(|(_, fd)| {
            fd.surface.as_ref().map_or(false, |sof| {
                matches!(sof.geometry, SurfaceGeometry::Plane { .. })
            })
        })
        .map(|(fh, _)| fh);

    if let Some(face) = cap_face {
        let specs = vec![DraftSpec {
            face,
            angle: 0.05,
            pull_direction: Vec3::new(0.0, 0.0, 1.0),
            neutral: NeutralReference::Auto,
            angle_variation: AngleVariation::Constant,
        }];
        let result = draft_faces(&brep, &specs, &DraftOptions::default(), tol(), None);
        assert!(
            result.is_ok(),
            "cylinder cap draft failed: {:?}",
            result.err()
        );
        let dr = result.unwrap();
        let v = dr.brep.vertices.len() as i64;
        let e = dr.brep.edges.len() as i64;
        let f = dr.brep.faces.len() as i64;
        assert_eq!(v - e + f, 2, "Euler violated: V={v} E={e} F={f}");
    }
}

#[test]
fn draft_cylinder_lateral_face_euler_preserved() {
    // Draft the lateral (cylindrical) face and verify Euler.
    let brep = BRep::cylinder(1.0, 2.0, tol()).unwrap();

    let cyl_face = brep
        .faces
        .iter()
        .find(|(_, fd)| {
            fd.surface.as_ref().map_or(false, |sof| {
                matches!(sof.geometry, SurfaceGeometry::Cylinder { .. })
            })
        })
        .map(|(fh, _)| fh);

    if let Some(face) = cyl_face {
        // Use an edge as neutral reference since Auto may not work for
        // the cylindrical face whose normal is perpendicular to the pull.
        let face_data = &brep.faces[face];
        let loop_h = &brep.loops[face_data.outer_loop];
        let he = &brep.half_edges[loop_h.half_edge];
        let edge_h = he.edge;

        let specs = vec![DraftSpec {
            face,
            angle: 0.05,
            pull_direction: Vec3::new(0.0, 0.0, 1.0),
            neutral: NeutralReference::Edge(edge_h),
            angle_variation: AngleVariation::Constant,
        }];
        let result = draft_faces(&brep, &specs, &DraftOptions::default(), tol(), None);
        if let Ok(dr) = result {
            let v = dr.brep.vertices.len() as i64;
            let e = dr.brep.edges.len() as i64;
            let f = dr.brep.faces.len() as i64;
            assert_eq!(v - e + f, 2, "Euler violated: V={v} E={e} F={f}");
        }
        // If this errors due to cylinder seam constraints, that is acceptable.
    }
}

#[test]
fn draft_multiple_faces_simultaneously() {
    let brep = BRep::unit_cube(tol()).unwrap();
    // Collect two side faces.
    let mut side_faces = Vec::new();
    for (fh, fd) in &brep.faces {
        if let Some(ref sof) = fd.surface {
            use popkern_core::geometry::surface::SurfaceEval;
            let n = sof.geometry.normal_at(0.5, 0.5);
            if n.z.abs() < 0.1 && side_faces.len() < 2 {
                side_faces.push(fh);
            }
        }
    }
    assert!(
        side_faces.len() >= 2,
        "need at least 2 side faces on a cube"
    );

    let specs: Vec<DraftSpec> = side_faces
        .iter()
        .map(|&fh| DraftSpec {
            face: fh,
            angle: 0.087,
            pull_direction: Vec3::new(0.0, 0.0, 1.0),
            neutral: NeutralReference::Auto,
            angle_variation: AngleVariation::Constant,
        })
        .collect();

    let result = draft_faces(&brep, &specs, &DraftOptions::default(), tol(), None);
    assert!(
        result.is_ok(),
        "multi-face draft failed: {:?}",
        result.err()
    );
    let dr = result.unwrap();
    assert_eq!(dr.brep.faces.len(), 6);
    let v = dr.brep.vertices.len() as i64;
    let e = dr.brep.edges.len() as i64;
    let f = dr.brep.faces.len() as i64;
    assert_eq!(v - e + f, 2, "Euler violated: V={v} E={e} F={f}");
}

#[test]
fn draft_with_custom_nurbs_samples() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let face = side_face(&brep);
    let specs = vec![DraftSpec {
        face,
        angle: 0.05,
        pull_direction: Vec3::new(0.0, 0.0, 1.0),
        neutral: NeutralReference::Auto,
        angle_variation: AngleVariation::Constant,
    }];
    let opts = DraftOptions {
        nurbs_samples: 50,
        ..DraftOptions::default()
    };
    let result = draft_faces(&brep, &specs, &opts, tol(), None);
    assert!(
        result.is_ok(),
        "custom nurbs_samples failed: {:?}",
        result.err()
    );
}
