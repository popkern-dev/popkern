use super::*;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Find a face whose surface normal aligns with `axis` (dot > 0.9).
fn find_face_by_normal(brep: &BRep, axis: Vec3) -> Option<popkern_core::topology::FaceHandle> {
    brep.faces
        .iter()
        .find(|(_, fd)| {
            fd.surface.as_ref().map_or(false, |sof| {
                sof.geometry.normal_at(0.5, 0.5).dot(&axis).abs() > 0.9
            })
        })
        .map(|(fh, _)| fh)
}

// ---------------------------------------------------------------------------
// Existing tests (migrated from inline module)
// ---------------------------------------------------------------------------

#[test]
fn parting_plane_through_cube_classifies_faces() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let parting = PartingLine::Plane {
        origin: Point3::new(0.0, 0.25, 0.0),
        normal: Vec3::new(0.0, 1.0, 0.0),
    };
    let dr = draft_parting_line(
        &brep,
        &parting,
        0.05,
        0.05,
        Vec3::new(0.0, 1.0, 0.0),
        &DraftOptions::default(),
        tol(),
    )
    .unwrap();
    assert_eq!(dr.brep.faces.len(), 6);
}

#[test]
fn parting_line_zero_normal_returns_error() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let parting = PartingLine::Plane {
        origin: Point3::new(0.0, 0.0, 0.5),
        normal: Vec3::zero(),
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
        .is_err()
    );
}

#[test]
fn parting_line_empty_edges_returns_error() {
    let brep = BRep::unit_cube(tol()).unwrap();
    assert!(
        draft_parting_line(
            &brep,
            &PartingLine::Edges(vec![]),
            0.05,
            0.05,
            Vec3::new(0.0, 0.0, 1.0),
            &DraftOptions::default(),
            tol()
        )
        .is_err()
    );
}

#[test]
fn face_centroid_on_cube() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let face = brep.faces.keys().next().unwrap();
    let c = compute_face_centroid(&brep, face);
    assert!(c.x >= -0.1 && c.x <= 1.1);
    assert!(c.y >= -0.1 && c.y <= 1.1);
    assert!(c.z >= -0.1 && c.z <= 1.1);
}

#[test]
fn compute_perpendicular_is_orthogonal() {
    for v in [
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        Vec3::new(1.0, 1.0, 1.0),
    ] {
        let perp = compute_perpendicular(&v);
        assert!(perp.dot(&v).abs() < 1e-10, "dot={}", perp.dot(&v));
        assert!((perp.length() - 1.0).abs() < 1e-10, "len={}", perp.length());
    }
}

#[test]
fn parting_line_edges_variant_classifies_faces() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let edge_h = brep.edges.keys().next().unwrap();
    let _ = draft_parting_line(
        &brep,
        &PartingLine::Edges(vec![edge_h]),
        0.05,
        0.05,
        Vec3::new(0.0, 0.0, 1.0),
        &DraftOptions::default(),
        tol(),
    );
}

#[test]
fn parting_line_edges_invalid_handle_returns_error() {
    let brep = BRep::unit_cube(tol()).unwrap();
    assert!(
        draft_parting_line(
            &brep,
            &PartingLine::Edges(vec![Default::default()]),
            0.05,
            0.05,
            Vec3::new(0.0, 0.0, 1.0),
            &DraftOptions::default(),
            tol()
        )
        .is_err()
    );
}

#[test]
fn faces_parallel_to_parting_plane_are_skipped() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let parting = PartingLine::Plane {
        origin: Point3::new(0.5, 0.5, 0.3),
        normal: Vec3::new(0.0, 0.0, 1.0),
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
fn all_faces_on_parting_plane_returns_error() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let big_tol = Tolerance {
        linear: 10.0,
        angular: 1e-10,
        parametric: 1e-8,
    };
    let parting = PartingLine::Plane {
        origin: Point3::new(0.5, 0.5, 0.5),
        normal: Vec3::new(1.0, 0.0, 0.0),
    };
    assert!(
        draft_parting_line(
            &brep,
            &parting,
            0.05,
            0.05,
            Vec3::new(1.0, 0.0, 0.0),
            &DraftOptions::default(),
            big_tol
        )
        .is_err()
    );
}

#[test]
fn zero_above_and_below_angles_returns_error() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let parting = PartingLine::Plane {
        origin: Point3::new(0.5, 0.5, 0.25),
        normal: Vec3::new(0.0, 0.0, 1.0),
    };
    assert!(
        draft_parting_line(
            &brep,
            &parting,
            0.0,
            0.0,
            Vec3::new(0.0, 0.0, 1.0),
            &DraftOptions::default(),
            tol()
        )
        .is_err()
    );
}

#[test]
fn parting_curve_resolves_without_panic() {
    use popkern_core::geometry::NurbsCurve;
    let brep = BRep::unit_cube(tol()).unwrap();
    let curve = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.5, 0.0), Point3::new(1.0, 0.5, 0.0)],
        weights: vec![1.0, 1.0],
    };
    let _ = draft_parting_line(
        &brep,
        &PartingLine::Curve(curve),
        0.05,
        0.05,
        Vec3::new(0.0, 0.0, 1.0),
        &DraftOptions::default(),
        tol(),
    );
}

// ---------------------------------------------------------------------------
// Gap 1: draft_parting_line geometric verification
// ---------------------------------------------------------------------------

#[test]
fn parting_plane_at_midpoint_skips_parallel_top_bottom() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let parting = PartingLine::Plane {
        origin: Point3::new(0.5, 0.5, 0.4),
        normal: Vec3::new(0.0, 0.0, 1.0),
    };
    let dr = draft_parting_line(
        &brep,
        &parting,
        0.05,
        0.0,
        Vec3::new(0.0, 0.0, 1.0),
        &DraftOptions::default(),
        tol(),
    )
    .unwrap();
    assert_eq!(dr.brep.faces.len(), 6);
    assert!(find_face_by_normal(&dr.brep, Vec3::new(0.0, 0.0, 1.0)).is_some());
    assert!(find_face_by_normal(&dr.brep, Vec3::new(0.0, 0.0, -1.0)).is_some());
}

#[test]
fn parting_line_preserves_euler_characteristic() {
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
    let v = dr.brep.vertices.len() as i64;
    let e = dr.brep.edges.len() as i64;
    let f = dr.brep.faces.len() as i64;
    assert_eq!(v - e + f, 2, "Euler: V={v} E={e} F={f}");
}

#[test]
fn parting_line_above_below_different_angles() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let parting = PartingLine::Plane {
        origin: Point3::new(0.25, 0.5, 0.5),
        normal: Vec3::new(1.0, 0.0, 0.0),
    };
    let dr = draft_parting_line(
        &brep,
        &parting,
        0.08,
        0.04,
        Vec3::new(1.0, 0.0, 0.0),
        &DraftOptions::default(),
        tol(),
    )
    .unwrap();
    assert_eq!(dr.brep.faces.len(), 6);
}

#[test]
fn parting_line_only_above_nonzero() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let parting = PartingLine::Plane {
        origin: Point3::new(0.5, 0.5, 0.3),
        normal: Vec3::new(0.0, 0.0, 1.0),
    };
    assert!(
        draft_parting_line(
            &brep,
            &parting,
            0.06,
            0.0,
            Vec3::new(0.0, 0.0, 1.0),
            &DraftOptions::default(),
            tol()
        )
        .is_ok()
    );
}

#[test]
fn parting_line_only_below_nonzero() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let parting = PartingLine::Plane {
        origin: Point3::new(0.5, 0.5, 0.7),
        normal: Vec3::new(0.0, 0.0, 1.0),
    };
    assert!(
        draft_parting_line(
            &brep,
            &parting,
            0.0,
            0.06,
            Vec3::new(0.0, 0.0, 1.0),
            &DraftOptions::default(),
            tol()
        )
        .is_ok()
    );
}

#[test]
fn parting_plane_all_parallel_returns_error() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let wide_tol = Tolerance {
        linear: 1e-6,
        angular: 1.5,
        parametric: 1e-8,
    };
    let parting = PartingLine::Plane {
        origin: Point3::new(0.5, 0.5, 0.5),
        normal: Vec3::new(0.0, 0.0, 1.0),
    };
    assert!(
        draft_parting_line(
            &brep,
            &parting,
            0.05,
            0.05,
            Vec3::new(0.0, 0.0, 1.0),
            &DraftOptions::default(),
            wide_tol
        )
        .is_err()
    );
}
