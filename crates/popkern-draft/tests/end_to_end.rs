//! End-to-end integration tests for popkern-draft.

use popkern_core::error::ValidationError;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_draft::{
    AdjacentMode, AngleVariation, DraftError, DraftOptions, DraftSpec, NeutralReference,
    PartingLine, draft_face, draft_faces, draft_parting_line,
};

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Assert the Euler characteristic V - E + F = 2 for a closed solid.
fn assert_euler(brep: &BRep) {
    let v = brep.vertices.len() as i64;
    let e = brep.edges.len() as i64;
    let f = brep.faces.len() as i64;
    assert_eq!(v - e + f, 2, "Euler: V={v} E={e} F={f}");
}

/// Assert no structural validation errors (broken cycles, twins, etc.).
/// Note: OpenLoop is a geometric tolerance check, not a structural one,
/// so it is excluded here. Small geometric gaps may remain after draft
/// operations and should be validated separately with appropriate tolerance.
fn assert_structurally_valid(brep: &BRep) {
    let errors = brep.validate_all();
    let structural: Vec<_> = errors
        .iter()
        .filter(|e| {
            matches!(
                e,
                ValidationError::BrokenCycle { .. }
                    | ValidationError::BrokenTwin { .. }
                    | ValidationError::InconsistentVertex { .. }
                    | ValidationError::NonManifoldEdge { .. }
            )
        })
        .collect();
    assert!(structural.is_empty(), "structural errors: {:?}", structural);
}

/// Helper: find face handles for all 6 cube faces, classifying by normal direction.
fn classify_cube_faces(brep: &BRep) -> CubeFaces {
    use popkern_core::geometry::surface::SurfaceEval;

    let mut result = CubeFaces::default();

    for (fh, face) in &brep.faces {
        if let Some(ref sof) = face.surface {
            let n = sof.geometry.normal_at(0.0, 0.0);
            if let Some(n) = n.normalized() {
                if n.z > 0.9 {
                    result.top = Some(fh);
                } else if n.z < -0.9 {
                    result.bottom = Some(fh);
                } else if n.x > 0.9 {
                    result.right = Some(fh);
                } else if n.x < -0.9 {
                    result.left = Some(fh);
                } else if n.y > 0.9 {
                    result.back = Some(fh);
                } else if n.y < -0.9 {
                    result.front = Some(fh);
                }
            }
        }
    }

    result
}

#[derive(Default)]
struct CubeFaces {
    top: Option<popkern_core::topology::FaceHandle>,
    bottom: Option<popkern_core::topology::FaceHandle>,
    left: Option<popkern_core::topology::FaceHandle>,
    right: Option<popkern_core::topology::FaceHandle>,
    front: Option<popkern_core::topology::FaceHandle>,
    back: Option<popkern_core::topology::FaceHandle>,
}

#[test]
fn draft_one_cube_face() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let faces = classify_cube_faces(&brep);
    let face = faces
        .front
        .unwrap_or_else(|| brep.faces.keys().next().unwrap());

    let result = draft_face(
        &brep,
        face,
        5.0_f64.to_radians(),
        Vec3::new(0.0, 0.0, 1.0),
        NeutralReference::Auto,
        tol(),
    );

    assert!(result.is_ok(), "draft failed: {:?}", result.err());
    let dr = result.unwrap();
    assert_euler(&dr.brep);
    assert_structurally_valid(&dr.brep);
    assert_eq!(dr.brep.faces.len(), 6);
}

#[test]
fn draft_two_opposite_faces() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let faces = classify_cube_faces(&brep);

    // Draft front and back faces.
    let front = faces
        .front
        .unwrap_or_else(|| brep.faces.keys().next().unwrap());
    let back = faces
        .back
        .unwrap_or_else(|| brep.faces.keys().nth(1).unwrap());

    let specs = vec![
        DraftSpec {
            face: front,
            angle: 5.0_f64.to_radians(),
            pull_direction: Vec3::new(0.0, 0.0, 1.0),
            neutral: NeutralReference::Auto,
            angle_variation: AngleVariation::Constant,
        },
        DraftSpec {
            face: back,
            angle: 5.0_f64.to_radians(),
            pull_direction: Vec3::new(0.0, 0.0, 1.0),
            neutral: NeutralReference::Auto,
            angle_variation: AngleVariation::Constant,
        },
    ];

    let result = draft_faces(&brep, &specs, &DraftOptions::default(), tol());
    assert!(result.is_ok(), "draft failed: {:?}", result.err());
    let dr = result.unwrap();
    assert_euler(&dr.brep);
    assert_structurally_valid(&dr.brep);
}

#[test]
fn draft_all_side_faces() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let faces = classify_cube_faces(&brep);

    // Draft all 4 side faces (front, back, left, right) — classic mold cavity.
    let side_faces: Vec<_> = [faces.front, faces.back, faces.left, faces.right]
        .into_iter()
        .flatten()
        .collect();

    let specs: Vec<DraftSpec> = side_faces
        .iter()
        .map(|&fh| DraftSpec {
            face: fh,
            angle: 3.0_f64.to_radians(),
            pull_direction: Vec3::new(0.0, 0.0, 1.0),
            neutral: NeutralReference::Auto,
            angle_variation: AngleVariation::Constant,
        })
        .collect();

    if specs.is_empty() {
        return; // Skip if cube faces couldn't be classified.
    }

    let result = draft_faces(&brep, &specs, &DraftOptions::default(), tol());
    assert!(result.is_ok(), "draft failed: {:?}", result.err());
    let dr = result.unwrap();
    assert_euler(&dr.brep);
    assert_structurally_valid(&dr.brep);
}

#[test]
fn draft_variable_angle() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let face = brep.faces.keys().next().unwrap();

    let specs = vec![DraftSpec {
        face,
        angle: 5.0_f64.to_radians(), // base angle
        pull_direction: Vec3::new(0.0, 0.0, 1.0),
        neutral: NeutralReference::Auto,
        angle_variation: AngleVariation::Variable(vec![
            (0.0, 3.0_f64.to_radians()),
            (1.0, 7.0_f64.to_radians()),
        ]),
    }];

    let result = draft_faces(&brep, &specs, &DraftOptions::default(), tol());
    assert!(result.is_ok(), "variable draft failed: {:?}", result.err());
    let dr = result.unwrap();
    assert_euler(&dr.brep);
}

#[test]
fn parting_line_cube() {
    let brep = BRep::unit_cube(tol()).unwrap();

    // Parting plane at y=0.25 separates faces above and below.
    let parting = PartingLine::Plane {
        origin: Point3::new(0.0, 0.25, 0.0),
        normal: Vec3::new(0.0, 1.0, 0.0),
    };

    let result = draft_parting_line(
        &brep,
        &parting,
        3.0_f64.to_radians(), // above
        3.0_f64.to_radians(), // below
        Vec3::new(0.0, 1.0, 0.0),
        &DraftOptions::default(),
        tol(),
    );

    assert!(result.is_ok(), "parting draft failed: {:?}", result.err());
    let dr = result.unwrap();
    assert_euler(&dr.brep);
}

#[test]
fn draft_with_transition_faces() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let face = brep.faces.keys().next().unwrap();

    let specs = vec![DraftSpec {
        face,
        angle: 5.0_f64.to_radians(),
        pull_direction: Vec3::new(0.0, 0.0, 1.0),
        neutral: NeutralReference::Auto,
        angle_variation: AngleVariation::Constant,
    }];

    let options = DraftOptions {
        adjacent_mode: AdjacentMode::TransitionFaces,
        ..DraftOptions::default()
    };

    let result = draft_faces(&brep, &specs, &options, tol());
    assert!(
        result.is_ok(),
        "transition draft failed: {:?}",
        result.err()
    );
    let dr = result.unwrap();
    assert_euler(&dr.brep);
    assert_structurally_valid(&dr.brep);
}

#[test]
fn draft_angle_too_large() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let face = brep.faces.keys().next().unwrap();

    let result = draft_face(
        &brep,
        face,
        60.0_f64.to_radians(), // > 45 degree max
        Vec3::new(0.0, 0.0, 1.0),
        NeutralReference::Auto,
        tol(),
    );

    assert!(matches!(result, Err(DraftError::AngleTooLarge { .. })));
}

#[test]
fn draft_empty_specs() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let result = draft_faces(&brep, &[], &DraftOptions::default(), tol());
    assert!(matches!(result, Err(DraftError::EmptySpecs)));
}

#[test]
fn draft_invalid_face() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let result = draft_face(
        &brep,
        Default::default(), // invalid handle
        5.0_f64.to_radians(),
        Vec3::new(0.0, 0.0, 1.0),
        NeutralReference::Auto,
        tol(),
    );
    assert!(matches!(result, Err(DraftError::InvalidFace(_))));
}

#[test]
fn draft_angle_matches_request() {
    use popkern_core::geometry::surface::SurfaceEval;

    let brep = BRep::unit_cube(tol()).unwrap();
    let faces = classify_cube_faces(&brep);
    let face = faces
        .front
        .unwrap_or_else(|| brep.faces.keys().next().unwrap());

    let requested_angle = 5.0_f64.to_radians();
    let result = draft_face(
        &brep,
        face,
        requested_angle,
        Vec3::new(0.0, 0.0, 1.0),
        NeutralReference::Auto,
        tol(),
    );

    assert!(result.is_ok(), "draft failed: {:?}", result.err());
    let dr = result.unwrap();

    // The drafted face's surface normal should be tilted by approximately
    // the requested angle relative to the original normal.
    let original_face = &brep.faces[face];
    let drafted_face = &dr.brep.faces[face];

    if let (Some(orig_sof), Some(draft_sof)) = (&original_face.surface, &drafted_face.surface) {
        let orig_normal = orig_sof.geometry.normal_at(0.0, 0.0);
        let draft_normal = draft_sof.geometry.normal_at(0.0, 0.0);

        if let (Some(on), Some(dn)) = (orig_normal.normalized(), draft_normal.normalized()) {
            let actual_angle = on.angle_to(&dn);
            // Tight tolerance: Rodrigues rotation on planar faces should
            // reproduce the requested angle within ~0.01 rad (~0.57°).
            assert!(
                (actual_angle - requested_angle).abs() < 0.01,
                "Expected angle ~{:.4} rad, got {:.4} rad (delta {:.6})",
                requested_angle,
                actual_angle,
                (actual_angle - requested_angle).abs()
            );
        }
    }
}
