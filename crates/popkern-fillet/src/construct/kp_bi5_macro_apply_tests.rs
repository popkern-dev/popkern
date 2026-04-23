//! Tests for the KP13 bi-5 macro-patch production dispatch wrapper.

use super::{
    tensor_border_kp_polar_to_kp13, try_kp_bi5_macro_cap, try_kp_bi5_macro_cap_as_gregory,
};
use popkern_core::geometry::Point3;
use popkern_core::geometry::Vec3;
use popkern_core::geometry::gregory::{
    OuterCornerData, PolarBoundaryInput, SpokeBoundaryData, TangentPlane, TensorBorder,
    TensorBorderBi3, TensorBorderBi5,
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

#[test]
fn gregory_wrapper_preserves_count_and_cps() {
    // The Gregory wrapper must return exactly one wrapper per bi-9
    // grid, in the same order, with identical CP contents.
    let polar = regular_polar(6);
    let grids = try_kp_bi5_macro_cap(&polar).unwrap();
    let wrapped = try_kp_bi5_macro_cap_as_gregory(&polar).unwrap();
    assert_eq!(wrapped.len(), grids.len());
    assert_eq!(wrapped.len(), 4 * 6);
    for (g, w) in grids.iter().zip(wrapped.iter()) {
        assert_eq!(&w.grid, g);
    }
}

/// Evaluate a bi-3 Bézier patch at `(u, v)` ∈ `[0, 1]²`.
fn eval_bi3(cps: &[[Point3; 4]; 4], u: f64, v: f64) -> Point3 {
    let bu: [f64; 4] = {
        let s = 1.0 - u;
        [s * s * s, 3.0 * u * s * s, 3.0 * u * u * s, u * u * u]
    };
    let bv: [f64; 4] = {
        let s = 1.0 - v;
        [s * s * s, 3.0 * v * s * s, 3.0 * v * v * s, v * v * v]
    };
    let mut out = Point3::new(0.0, 0.0, 0.0);
    for i in 0..4 {
        for j in 0..4 {
            let w = bu[i] * bv[j];
            out.x += w * cps[i][j].x;
            out.y += w * cps[i][j].y;
            out.z += w * cps[i][j].z;
        }
    }
    out
}

/// Evaluate a bi-2 Bézier column at `v` ∈ `[0, 1]`, given `[C_0, C_1, C_2]`.
fn eval_bi2_col(c: [Point3; 3], v: f64) -> Point3 {
    let s = 1.0 - v;
    Point3::new(
        s * s * c[0].x + 2.0 * v * s * c[1].x + v * v * c[2].x,
        s * s * c[0].y + 2.0 * v * s * c[1].y + v * v * c[2].y,
        s * s * c[0].z + 2.0 * v * s * c[1].z + v * v * c[2].z,
    )
}

#[test]
fn adapter_degree_elevation_preserves_polynomial() {
    // Construct a polar bi-3 tensor-border with each column a known
    // bi-2-in-v polynomial; verify that the adapted KP13 bi-3 patch
    // evaluates to the v-reversed polynomial at any (u, v).
    // Position (C_0) is polar v=0 (outer mouth); second-cross (C_2)
    // is polar v=1 (deepest inward). KP13 v_kp13 = 1 − v_polar.
    let mut sectors = Vec::with_capacity(1);
    // Four columns of distinct polynomials.
    let columns: [[Point3; 3]; 4] = [
        [pt(0.0, 0.0, 0.0), pt(1.0, 0.0, 0.5), pt(2.0, 1.0, 1.0)],
        [pt(0.5, 0.0, 0.0), pt(1.5, 0.5, 0.25), pt(2.5, 1.5, 0.75)],
        [pt(1.0, 0.0, 0.0), pt(2.0, 0.75, 0.0), pt(3.0, 2.0, 0.5)],
        [pt(1.5, 0.0, 0.0), pt(2.5, 1.0, -0.25), pt(3.5, 2.5, 0.25)],
    ];
    let mut rows: [[Point3; 4]; 3] = [[pt(0.0, 0.0, 0.0); 4]; 3];
    for (i, col) in columns.iter().enumerate() {
        rows[0][i] = col[0]; // C_0
        rows[1][i] = col[1]; // C_1
        rows[2][i] = col[2]; // C_2
    }
    sectors.push(rows);
    let polar_tb = TensorBorder::Bi3(TensorBorderBi3 { sectors });
    let kp13 = tensor_border_kp_polar_to_kp13(&polar_tb).unwrap();
    let patch = &kp13.patches[0].cps;

    // Sample at (u, v_kp13) ∈ grid and verify equality with the
    // polar column's value at v_polar = 1 − v_kp13.
    for iu in 0..=4 {
        let u = iu as f64 / 4.0;
        for iv in 0..=4 {
            let v_kp13 = iv as f64 / 4.0;
            let v_polar = 1.0 - v_kp13;
            let got = eval_bi3(patch, u, v_kp13);
            // u-interp across the 4 polar columns (linear-in-CP bi-3
            // boundary — polar extraction is bi-3 in u as well).
            let bu: [f64; 4] = {
                let s = 1.0 - u;
                [s * s * s, 3.0 * u * s * s, 3.0 * u * u * s, u * u * u]
            };
            let mut expect = Point3::new(0.0, 0.0, 0.0);
            for i in 0..4 {
                let pp = eval_bi2_col(columns[i], v_polar);
                expect.x += bu[i] * pp.x;
                expect.y += bu[i] * pp.y;
                expect.z += bu[i] * pp.z;
            }
            let eps = 1e-12;
            assert!(
                (got.x - expect.x).abs() < eps
                    && (got.y - expect.y).abs() < eps
                    && (got.z - expect.z).abs() < eps,
                "(u={u}, v_kp13={v_kp13}): got {got:?} expect {expect:?}"
            );
        }
    }
}

#[test]
fn adapter_rejects_bi5_input() {
    let sectors = vec![[[pt(0.0, 0.0, 0.0); 6]; 3]];
    let polar_tb = TensorBorder::Bi5(TensorBorderBi5 { sectors });
    assert!(tensor_border_kp_polar_to_kp13(&polar_tb).is_none());
}

#[test]
fn gregory_wrapper_low_valence() {
    // Low-valence path (n ∈ {3, 4}) also flows through the wrapper.
    for n in [3_usize, 4, 5, 7, 8] {
        let polar = regular_polar(n);
        let wrapped =
            try_kp_bi5_macro_cap_as_gregory(&polar).unwrap_or_else(|e| panic!("n={n}: {e}"));
        assert_eq!(wrapped.len(), 4 * n);
    }
}
