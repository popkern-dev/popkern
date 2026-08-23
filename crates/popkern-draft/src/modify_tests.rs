use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::modify::modify_draft_angle;
use crate::types::{AngleVariation, DraftOptions, DraftSpec, NeutralReference};

fn tol() -> Tolerance {
    Tolerance::default()
}

/// A cube side face whose plane normal is +X.
fn side_face(brep: &BRep) -> popkern_core::topology::FaceHandle {
    brep.faces
        .iter()
        .find(|(_, fd)| {
            fd.surface.as_ref().is_some_and(|sof| {
                matches!(
                    &sof.geometry,
                    popkern_core::geometry::surface::SurfaceGeometry::Plane { normal, .. }
                        if (normal.x - 1.0).abs() < 1e-9
                )
            })
        })
        .map(|(f, _)| f)
        .expect("cube has a +X side face")
}

fn bottom_neutral() -> NeutralReference {
    NeutralReference::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
    }
}

/// The additive-composition claim, measured: draft 5° then modify to
/// 10° equals a fresh 10° draft (vertex positions within tolerance).
#[test]
fn angle_change_matches_fresh_draft() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let face = side_face(&cube);
    let pull = Vec3::new(0.0, 0.0, 1.0);
    let a1 = 5.0_f64.to_radians();
    let a2 = 10.0_f64.to_radians();

    let spec = |angle: f64| DraftSpec {
        face,
        angle,
        pull_direction: pull,
        neutral: bottom_neutral(),
        angle_variation: AngleVariation::Constant,
    };
    let drafted = crate::draft_faces(&cube, &[spec(a1)], &DraftOptions::default(), tol(), None)
        .expect("initial draft serves")
        .brep;
    let modified = modify_draft_angle(&drafted, face, a1, a2, pull, bottom_neutral(), tol())
        .expect("in-place angle change serves")
        .brep;
    modified.validate().expect("modified body validates");

    let fresh = crate::draft_faces(&cube, &[spec(a2)], &DraftOptions::default(), tol(), None)
        .expect("fresh draft serves")
        .brep;

    let corners = |b: &BRep| -> Vec<[f64; 3]> {
        let mut pts: Vec<[f64; 3]> = b
            .vertices
            .values()
            .map(|v| [v.position.x, v.position.y, v.position.z])
            .collect();
        pts.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        pts
    };
    let got = corners(&modified);
    let want = corners(&fresh);
    assert_eq!(got.len(), want.len());
    for (g, w) in got.iter().zip(want.iter()) {
        for k in 0..3 {
            assert!(
                (g[k] - w[k]).abs() < 1.0e-9,
                "composed-draft vertex {g:?} deviates from fresh-draft vertex {w:?}"
            );
        }
    }
    // Topology untouched.
    assert_eq!(modified.faces.len(), drafted.faces.len());
    assert_eq!(modified.edges.len(), drafted.edges.len());
    assert_eq!(modified.vertices.len(), drafted.vertices.len());
}

/// Shrinking back restores the original within tolerance.
#[test]
fn angle_round_trip_restores_geometry() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let face = side_face(&cube);
    let pull = Vec3::new(0.0, 0.0, 1.0);
    let a1 = 3.0_f64.to_radians();
    let a2 = 8.0_f64.to_radians();
    let drafted = crate::draft_faces(
        &cube,
        &[DraftSpec {
            face,
            angle: a1,
            pull_direction: pull,
            neutral: bottom_neutral(),
            angle_variation: AngleVariation::Constant,
        }],
        &DraftOptions::default(),
        tol(),
        None,
    )
    .expect("initial draft")
    .brep;
    let up = modify_draft_angle(&drafted, face, a1, a2, pull, bottom_neutral(), tol())
        .expect("grow")
        .brep;
    let back = modify_draft_angle(&up, face, a2, a1, pull, bottom_neutral(), tol())
        .expect("shrink")
        .brep;
    for (v, vd) in &drafted.vertices {
        let p = back.vertices.get(v).expect("same vertex set").position;
        assert!(
            vd.position.distance_to(&p) < 1.0e-9,
            "round-trip vertex drifted: {:?} vs {p:?}",
            vd.position
        );
    }
}

/// Auto neutral declines with the named reason.
#[test]
fn auto_neutral_declines_named() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let face = side_face(&cube);
    let err = modify_draft_angle(
        &cube,
        face,
        0.1,
        0.2,
        Vec3::new(0.0, 0.0, 1.0),
        NeutralReference::Auto,
        tol(),
    )
    .expect_err("Auto declines");
    assert!(err.to_string().contains("anchored"), "{err}");
}

/// Non-planar face declines.
#[test]
fn non_planar_face_declines() {
    let cyl = BRep::cylinder(0.5, 2.0, tol()).unwrap();
    let wall = cyl
        .faces
        .iter()
        .find(|(_, fd)| {
            fd.surface.as_ref().is_some_and(|sof| {
                matches!(
                    sof.geometry,
                    popkern_core::geometry::surface::SurfaceGeometry::Cylinder { .. }
                )
            })
        })
        .map(|(f, _)| f)
        .expect("cylinder wall");
    let err = modify_draft_angle(
        &cyl,
        wall,
        0.1,
        0.2,
        Vec3::new(0.0, 0.0, 1.0),
        NeutralReference::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
        },
        tol(),
    )
    .expect_err("non-planar declines");
    assert!(err.to_string().contains("planar"), "{err}");
}
