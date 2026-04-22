//! Tests for the KP13 bi-5 macro-patch production dispatch wrapper.

use super::try_kp_bi5_macro_cap;
use popkern_core::geometry::Point3;
use popkern_core::geometry::Vec3;
use popkern_core::geometry::gregory::{
    OuterCornerData, PolarBoundaryInput, SpokeBoundaryData, TangentPlane,
};

fn pt(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

fn tp_xy() -> TangentPlane {
    TangentPlane::from_basis(Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0)).unwrap()
}

fn tp_radial(theta: f64) -> TangentPlane {
    let u = Vec3::new(-theta.sin(), theta.cos(), 0.0);
    let v = Vec3::new(theta.cos(), theta.sin(), 0.0);
    TangentPlane::from_basis(u, v).unwrap()
}

/// Regular N-gon polar input at origin. Mirrors the bi-6 polar
/// test fixture for consistency.
fn regular_polar(n: usize) -> PolarBoundaryInput {
    let center = pt(0.0, 0.0, 0.0);
    let mut spokes = Vec::with_capacity(n);
    let mut outer_corners = Vec::with_capacity(n);
    let mut face_tangent_planes = Vec::with_capacity(n);
    let mut per_spoke_radii = Vec::with_capacity(n);
    for k in 0..n {
        let theta = (2.0 * std::f64::consts::PI) * (k as f64) / (n as f64);
        let v_outer = pt(theta.cos(), theta.sin(), 0.0);
        let mut boundary_cps = [pt(0.0, 0.0, 0.0); 8];
        for (i, cp) in boundary_cps.iter_mut().enumerate() {
            let t = i as f64 / 7.0;
            *cp = pt(t * v_outer.x, t * v_outer.y, 0.0);
        }
        spokes.push(SpokeBoundaryData {
            boundary_cps,
            xtan_a_at_c: Vec3::new(0.0, 1.0, 0.0),
            xtan_b_at_c: Vec3::new(0.0, -1.0, 0.0),
            xtan_a_at_v: Vec3::new(0.0, 1.0, 0.0),
            xtan_b_at_v: Vec3::new(0.0, -1.0, 0.0),
            radius_at_c: 0.1,
            radius_at_v: 0.1,
            edge_id: format!("edge_{k}"),
        });
        outer_corners.push(OuterCornerData {
            position: v_outer,
            tangent_plane: tp_xy(),
        });
        face_tangent_planes.push(tp_radial(theta));
        per_spoke_radii.push(0.1);
    }
    PolarBoundaryInput {
        central_vertex: center,
        central_tangent_plane: tp_xy(),
        spokes,
        outer_corners,
        face_tangent_planes,
        per_spoke_radii,
    }
}

#[test]
fn pentagon_produces_4n_patches() {
    let polar = regular_polar(5);
    let patches = try_kp_bi5_macro_cap(&polar).expect("pentagon KP13 bi-5 succeeds");
    assert_eq!(patches.len(), 4 * 5);
    for patch in &patches {
        for row in patch {
            for cp in row {
                assert!(cp.x.is_finite() && cp.y.is_finite() && cp.z.is_finite());
            }
        }
    }
}

#[test]
fn hexagon_produces_4n_patches() {
    let polar = regular_polar(6);
    let patches = try_kp_bi5_macro_cap(&polar).expect("hexagon KP13 bi-5 succeeds");
    assert_eq!(patches.len(), 4 * 6);
}

#[test]
fn multiple_valences_succeed() {
    for n in [5_usize, 6, 7, 8, 10] {
        let polar = regular_polar(n);
        let patches = try_kp_bi5_macro_cap(&polar).unwrap_or_else(|e| panic!("n={n}: {e}"));
        assert_eq!(patches.len(), 4 * n);
    }
}

#[test]
fn low_valence_n3_n4_supported() {
    // KP09a § 6 pseudo-inverse path flows through end-to-end.
    for n in [3_usize, 4] {
        let polar = regular_polar(n);
        let patches = try_kp_bi5_macro_cap(&polar).unwrap_or_else(|e| panic!("n={n}: {e}"));
        assert_eq!(patches.len(), 4 * n);
    }
}

#[test]
fn deterministic() {
    let polar = regular_polar(5);
    let a = try_kp_bi5_macro_cap(&polar).unwrap();
    let b = try_kp_bi5_macro_cap(&polar).unwrap();
    assert_eq!(a, b);
}

#[test]
fn invalid_polar_rejected() {
    // Build a degenerate polar with mismatched vector lengths.
    let mut polar = regular_polar(5);
    polar.spokes.pop();
    assert!(try_kp_bi5_macro_cap(&polar).is_err());
}
