//! Tests for [`super::DirectBoundaryInput`] and [`super::TangentPlane`].

use super::*;

fn demo_spoke(c: Point3, v: Point3, edge_id: &str) -> SpokeBoundaryData {
    // Bidegree-7 line from C to V as 8 equally-spaced CPs.
    let mut cps = [Point3::new(0.0, 0.0, 0.0); 8];
    for (j, cp) in cps.iter_mut().enumerate() {
        let t = (j as f64) / 7.0;
        *cp = Point3::new(
            c.x + t * (v.x - c.x),
            c.y + t * (v.y - c.y),
            c.z + t * (v.z - c.z),
        );
    }
    SpokeBoundaryData {
        boundary_cps: cps,
        xtan_a_at_c: Vec3::new(0.0, 1.0, 0.0),
        xtan_b_at_c: Vec3::new(0.0, -1.0, 0.0),
        xtan_a_at_v: Vec3::new(0.0, 1.0, 0.0),
        xtan_b_at_v: Vec3::new(0.0, -1.0, 0.0),
        radius_at_c: 0.1,
        radius_at_v: 0.1,
        edge_id: edge_id.into(),
    }
}

fn demo_tps() -> TangentPlane {
    TangentPlane::from_normal(Vec3::new(0.0, 0.0, 1.0)).unwrap()
}

#[test]
fn validate_rejects_fewer_than_three_spokes() {
    let input = DirectBoundaryInput {
        central_vertex: Point3::new(0.0, 0.0, 0.0),
        spokes: vec![demo_spoke(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            "e0",
        )],
        outer_corners: vec![OuterCornerData {
            position: Point3::new(1.0, 0.0, 0.0),
            tangent_plane: demo_tps(),
        }],
        central_vertex_tps: vec![demo_tps()],
        seed_patches: None,
    };
    assert!(input.validate().is_err());
}

#[test]
fn validate_rejects_mismatched_outer_corner_count() {
    let c = Point3::new(0.0, 0.0, 0.0);
    let input = DirectBoundaryInput {
        central_vertex: c,
        spokes: vec![
            demo_spoke(c, Point3::new(1.0, 0.0, 0.0), "e0"),
            demo_spoke(c, Point3::new(0.0, 1.0, 0.0), "e1"),
            demo_spoke(c, Point3::new(-1.0, -1.0, 0.0), "e2"),
        ],
        outer_corners: vec![OuterCornerData {
            position: Point3::new(1.0, 0.0, 0.0),
            tangent_plane: demo_tps(),
        }],
        central_vertex_tps: vec![demo_tps(); 3],
        seed_patches: None,
    };
    assert!(input.validate().is_err());
}

#[test]
fn validate_accepts_well_formed_three_sided() {
    let c = Point3::new(0.0, 0.0, 0.0);
    let vs = [
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(-1.0, -1.0, 0.0),
    ];
    let input = DirectBoundaryInput {
        central_vertex: c,
        spokes: (0..3)
            .map(|i| demo_spoke(c, vs[i], &format!("e{i}")))
            .collect(),
        outer_corners: vs
            .iter()
            .map(|&p| OuterCornerData {
                position: p,
                tangent_plane: demo_tps(),
            })
            .collect(),
        central_vertex_tps: vec![demo_tps(); 3],
        seed_patches: None,
    };
    assert!(input.validate().is_ok());
    assert_eq!(input.n(), 3);
}

#[test]
fn tangent_plane_from_normal_is_orthonormal() {
    let n = Vec3::new(1.0, 2.0, 3.0);
    let tps = TangentPlane::from_normal(n).unwrap();
    // u and v unit.
    assert!((tps.u.length() - 1.0).abs() < 1e-12);
    assert!((tps.v.length() - 1.0).abs() < 1e-12);
    // u ⊥ v.
    let dot = tps.u.x * tps.v.x + tps.u.y * tps.v.y + tps.u.z * tps.v.z;
    assert!(dot.abs() < 1e-12, "dot = {dot}");
    // normal || n.
    let n_unit = n.normalized().unwrap();
    let computed = tps.normal();
    let dot_n = n_unit.x * computed.x + n_unit.y * computed.y + n_unit.z * computed.z;
    assert!(dot_n.abs() > 0.9999, "normal not aligned");
}

#[test]
fn tangent_plane_from_basis_gram_schmidt() {
    let u = Vec3::new(2.0, 0.0, 0.0); // not unit
    let v = Vec3::new(1.0, 1.0, 0.0); // not orthogonal to u
    let tps = TangentPlane::from_basis(u, v).unwrap();
    // u should be normalised to (1, 0, 0).
    assert!((tps.u.x - 1.0).abs() < 1e-12);
    assert!(tps.u.y.abs() < 1e-12);
    // v should be orthogonalised to (0, 1, 0).
    assert!(tps.v.x.abs() < 1e-12);
    assert!((tps.v.y - 1.0).abs() < 1e-12);
}

#[test]
fn tangent_plane_from_basis_rejects_parallel() {
    let u = Vec3::new(1.0, 0.0, 0.0);
    let v = Vec3::new(2.0, 0.0, 0.0); // parallel
    assert!(TangentPlane::from_basis(u, v).is_none());
}

#[test]
fn tangent_plane_contains_in_plane_vector() {
    let tps = TangentPlane::from_normal(Vec3::new(0.0, 0.0, 1.0)).unwrap();
    assert!(tps.contains(Vec3::new(1.0, 2.0, 0.0), 1e-9));
    assert!(!tps.contains(Vec3::new(0.0, 0.0, 1.0), 1e-9));
}
