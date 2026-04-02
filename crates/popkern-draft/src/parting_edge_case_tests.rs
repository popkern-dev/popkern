use super::*;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default()
}

// ---------------------------------------------------------------------------
// Draft angle edge cases
// ---------------------------------------------------------------------------

#[test]
fn very_small_draft_angle_succeeds() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let parting = PartingLine::Plane {
        origin: Point3::new(0.5, 0.5, 0.3),
        normal: Vec3::new(0.0, 0.0, 1.0),
    };
    let dr = draft_parting_line(
        &brep,
        &parting,
        0.001,
        0.001,
        Vec3::new(0.0, 0.0, 1.0),
        &DraftOptions::default(),
        tol(),
    )
    .unwrap();
    assert_eq!(dr.brep.faces.len(), 6);
}

#[test]
fn angle_near_max_limit_succeeds() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let parting = PartingLine::Plane {
        origin: Point3::new(0.5, 0.5, 0.3),
        normal: Vec3::new(0.0, 0.0, 1.0),
    };
    let near_max = std::f64::consts::FRAC_PI_4 - 0.01;
    assert!(
        draft_parting_line(
            &brep,
            &parting,
            near_max,
            near_max,
            Vec3::new(0.0, 0.0, 1.0),
            &DraftOptions::default(),
            tol()
        )
        .is_ok()
    );
}

#[test]
fn angle_exceeding_max_returns_angle_too_large() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let parting = PartingLine::Plane {
        origin: Point3::new(0.5, 0.5, 0.3),
        normal: Vec3::new(0.0, 0.0, 1.0),
    };
    let over_max = std::f64::consts::FRAC_PI_4 + 0.01;
    let err = draft_parting_line(
        &brep,
        &parting,
        over_max,
        over_max,
        Vec3::new(0.0, 0.0, 1.0),
        &DraftOptions::default(),
        tol(),
    )
    .unwrap_err();
    assert!(err.to_string().contains("angle too large"), "got: {err}");
}

#[test]
fn custom_max_angle_allows_larger_draft() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let parting = PartingLine::Plane {
        origin: Point3::new(0.5, 0.5, 0.3),
        normal: Vec3::new(0.0, 0.0, 1.0),
    };
    let opts = DraftOptions {
        max_angle: std::f64::consts::FRAC_PI_2,
        ..DraftOptions::default()
    };
    let large = std::f64::consts::FRAC_PI_4 + 0.1;
    assert!(
        draft_parting_line(
            &brep,
            &parting,
            large,
            large,
            Vec3::new(0.0, 0.0, 1.0),
            &opts,
            tol()
        )
        .is_ok()
    );
}

#[test]
fn negative_draft_angle_succeeds() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let parting = PartingLine::Plane {
        origin: Point3::new(0.5, 0.5, 0.3),
        normal: Vec3::new(0.0, 0.0, 1.0),
    };
    let dr = draft_parting_line(
        &brep,
        &parting,
        -0.05,
        -0.05,
        Vec3::new(0.0, 0.0, 1.0),
        &DraftOptions::default(),
        tol(),
    )
    .unwrap();
    assert_eq!(dr.brep.faces.len(), 6);
}

// ---------------------------------------------------------------------------
// Different input variants and modes
// ---------------------------------------------------------------------------

#[test]
fn parting_plane_unnormalized_normal_succeeds() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let parting = PartingLine::Plane {
        origin: Point3::new(0.5, 0.5, 0.3),
        normal: Vec3::new(0.0, 0.0, 100.0),
    };
    assert!(
        draft_parting_line(
            &brep,
            &parting,
            0.05,
            0.05,
            Vec3::new(0.0, 0.0, 1.0),
            &DraftOptions::default(),
            tol()
        )
        .is_ok()
    );
}

#[test]
fn parting_line_edges_multiple() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let mut it = brep.edges.keys();
    let e1 = it.next().unwrap();
    let e2 = it.next().unwrap();
    let _ = draft_parting_line(
        &brep,
        &PartingLine::Edges(vec![e1, e2]),
        0.05,
        0.05,
        Vec3::new(0.0, 0.0, 1.0),
        &DraftOptions::default(),
        tol(),
    );
}

#[test]
fn face_centroids_are_distinct() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let cs: Vec<_> = brep
        .faces
        .keys()
        .map(|fh| compute_face_centroid(&brep, fh))
        .collect();
    for (i, a) in cs.iter().enumerate() {
        for b in cs.iter().skip(i + 1) {
            assert!(a.distance_to(b) > 0.01, "centroids too close: {a:?} {b:?}");
        }
    }
}

#[test]
fn parting_line_no_warnings_for_cube() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let parting = PartingLine::Plane {
        origin: Point3::new(0.5, 0.5, 0.3),
        normal: Vec3::new(0.0, 0.0, 1.0),
    };
    let dr = draft_parting_line(
        &brep,
        &parting,
        0.05,
        0.05,
        Vec3::new(0.0, 0.0, 1.0),
        &DraftOptions::default(),
        tol(),
    )
    .unwrap();
    assert!(
        dr.warnings.is_empty(),
        "unexpected warnings: {:?}",
        dr.warnings
    );
}

#[test]
fn asymmetric_angles_euler_valid() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let parting = PartingLine::Plane {
        origin: Point3::new(0.25, 0.5, 0.5),
        normal: Vec3::new(1.0, 0.0, 0.0),
    };
    let dr = draft_parting_line(
        &brep,
        &parting,
        1.0_f64.to_radians(),
        10.0_f64.to_radians(),
        Vec3::new(1.0, 0.0, 0.0),
        &DraftOptions::default(),
        tol(),
    )
    .unwrap();
    let v = dr.brep.vertices.len() as i64;
    let e = dr.brep.edges.len() as i64;
    let f = dr.brep.faces.len() as i64;
    assert_eq!(v - e + f, 2, "Euler: V={v} E={e} F={f}");
}

#[test]
fn extend_trim_mode_on_parting_line() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let parting = PartingLine::Plane {
        origin: Point3::new(0.5, 0.5, 0.3),
        normal: Vec3::new(0.0, 0.0, 1.0),
    };
    let opts = DraftOptions {
        adjacent_mode: crate::types::AdjacentMode::ExtendTrim,
        ..DraftOptions::default()
    };
    assert!(
        draft_parting_line(
            &brep,
            &parting,
            0.05,
            0.05,
            Vec3::new(0.0, 0.0, 1.0),
            &opts,
            tol()
        )
        .is_ok()
    );
}

#[test]
fn transition_faces_mode_on_parting_line() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let parting = PartingLine::Plane {
        origin: Point3::new(0.5, 0.5, 0.3),
        normal: Vec3::new(0.0, 0.0, 1.0),
    };
    let opts = DraftOptions {
        adjacent_mode: crate::types::AdjacentMode::TransitionFaces,
        ..DraftOptions::default()
    };
    assert!(
        draft_parting_line(
            &brep,
            &parting,
            0.05,
            0.05,
            Vec3::new(0.0, 0.0, 1.0),
            &opts,
            tol()
        )
        .is_ok()
    );
}

#[test]
fn parting_line_on_cylinder() {
    let brep = BRep::cylinder(1.0, 2.0, tol()).unwrap();
    let parting = PartingLine::Plane {
        origin: Point3::new(0.0, 0.0, 1.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
    };
    // May succeed or fail depending on face classification, but must not panic.
    let _ = draft_parting_line(
        &brep,
        &parting,
        0.05,
        0.05,
        Vec3::new(0.0, 0.0, 1.0),
        &DraftOptions::default(),
        tol(),
    );
}
