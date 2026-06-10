use super::{NormalCurvatureError, normal_curvature_along};
use crate::geometry::gregory::{GregoryPatch, RectangularGregory, RectangularGregoryBiNonic};
use crate::geometry::{Point3, Vec3};

const DEGENERACY_TOL: f64 = 1e-14;

/// Bi-9 patch reproducing the plane `z = 0` exactly via linear
/// precision of the Bernstein basis.
fn plane_binonic() -> GregoryPatch {
    let mut grid = [[Point3::new(0.0, 0.0, 0.0); 10]; 10];
    for (i, row) in grid.iter_mut().enumerate() {
        for (j, cp) in row.iter_mut().enumerate() {
            *cp = Point3::new(i as f64 / 9.0, j as f64 / 9.0, 0.0);
        }
    }
    GregoryPatch::RectangularBiNonic(Box::new(RectangularGregoryBiNonic { grid }))
}

/// Bi-9 patch reproducing the parabolic cylinder `z = u²`
/// exactly: the degree-9 Bezier coefficients of `t²` are
/// `b_i = i(i − 1) / (9 · 8)`.
fn parabolic_cylinder_binonic() -> GregoryPatch {
    let mut grid = [[Point3::new(0.0, 0.0, 0.0); 10]; 10];
    for (i, row) in grid.iter_mut().enumerate() {
        for (j, cp) in row.iter_mut().enumerate() {
            let z = (i as f64) * (i as f64 - 1.0) / 72.0;
            *cp = Point3::new(i as f64 / 9.0, j as f64 / 9.0, z);
        }
    }
    GregoryPatch::RectangularBiNonic(Box::new(RectangularGregoryBiNonic { grid }))
}

#[test]
fn plane_has_zero_normal_curvature() {
    let patch = plane_binonic();
    for &(u, v) in &[(0.25, 0.25), (0.5, 0.5), (0.8, 0.3)] {
        let curvature =
            normal_curvature_along(&patch, u, v, Vec3::new(1.0, 1.0, 0.0), DEGENERACY_TOL)
                .expect("plane sample must be measurable");
        assert!(
            curvature.abs() < 1e-12,
            "plane curvature should vanish, got {curvature}"
        );
    }
}

#[test]
fn parabolic_cylinder_matches_analytic_curvature() {
    let patch = parabolic_cylinder_binonic();
    // Along the bent direction at u = 0.5: the planar curve
    // z = x² has κ = z'' / (1 + z'²)^{3/2} = 2 / 2^{3/2} = 1/√2.
    // The patch normal (∂u × ∂v) points toward −x, +z, so the
    // surface bends away from the normal: κ_n = −1/√2... sign
    // check: n = (−1, 0, 1)/√2, S_uu = (0, 0, 2) ⇒ L = √2 > 0,
    // so κ_n = +1/√2 with this normal.
    let along_u =
        normal_curvature_along(&patch, 0.5, 0.5, Vec3::new(1.0, 0.0, 1.0), DEGENERACY_TOL)
            .expect("cylinder sample must be measurable");
    let expected = std::f64::consts::FRAC_1_SQRT_2;
    assert!(
        (along_u - expected).abs() < 1e-10,
        "expected κ ≈ {expected}, got {along_u}"
    );
    // Along the ruling direction the cylinder is flat.
    let along_v =
        normal_curvature_along(&patch, 0.5, 0.5, Vec3::new(0.0, 1.0, 0.0), DEGENERACY_TOL)
            .expect("ruling sample must be measurable");
    assert!(
        along_v.abs() < 1e-12,
        "ruling curvature should vanish, got {along_v}"
    );
}

#[test]
fn direction_magnitude_is_irrelevant() {
    let patch = parabolic_cylinder_binonic();
    let unit = normal_curvature_along(&patch, 0.5, 0.5, Vec3::new(1.0, 0.0, 1.0), DEGENERACY_TOL)
        .expect("unit-scale sample must be measurable");
    let scaled = normal_curvature_along(
        &patch,
        0.5,
        0.5,
        Vec3::new(250.0, 0.0, 250.0),
        DEGENERACY_TOL,
    )
    .expect("scaled sample must be measurable");
    assert!(
        (unit - scaled).abs() < 1e-12,
        "direction scaling changed κ: {unit} vs {scaled}"
    );
}

#[test]
fn normal_component_only_direction_is_degenerate() {
    let patch = plane_binonic();
    let result = normal_curvature_along(&patch, 0.5, 0.5, Vec3::new(0.0, 0.0, 1.0), DEGENERACY_TOL);
    assert_eq!(result, Err(NormalCurvatureError::DegenerateDirection));
}

#[test]
fn zero_patch_has_degenerate_normal() {
    let patch = GregoryPatch::RectangularBiNonic(Box::new(RectangularGregoryBiNonic::zeros()));
    let result = normal_curvature_along(&patch, 0.5, 0.5, Vec3::new(1.0, 0.0, 0.0), DEGENERACY_TOL);
    assert_eq!(result, Err(NormalCurvatureError::DegenerateNormal));
}

#[test]
fn near_parallel_partials_report_degenerate_first_form() {
    // su = (1, 0, 0), sv = (1, 1e-8, 0): the normal still
    // normalizes, but E·G − F² ≈ 1e-16 falls under the guard.
    let mut grid = [[Point3::new(0.0, 0.0, 0.0); 10]; 10];
    for (i, row) in grid.iter_mut().enumerate() {
        for (j, cp) in row.iter_mut().enumerate() {
            *cp = Point3::new(
                i as f64 / 9.0 + j as f64 / 9.0,
                j as f64 / 9.0 * 1.0e-8,
                0.0,
            );
        }
    }
    let patch = GregoryPatch::RectangularBiNonic(Box::new(RectangularGregoryBiNonic { grid }));
    let result = normal_curvature_along(&patch, 0.5, 0.5, Vec3::new(1.0, 0.0, 0.0), DEGENERACY_TOL);
    assert_eq!(
        result,
        Err(NormalCurvatureError::DegenerateFirstFundamentalForm)
    );
}

#[test]
fn unsupported_variant_reports_typed_error() {
    let origin = Point3::new(0.0, 0.0, 0.0);
    let patch = GregoryPatch::Rectangular(Box::new(RectangularGregory {
        boundary_cps: [[origin; 4]; 4],
        inner_cps: [(origin, origin); 4],
    }));
    let result = normal_curvature_along(&patch, 0.5, 0.5, Vec3::new(1.0, 0.0, 0.0), DEGENERACY_TOL);
    assert_eq!(
        result,
        Err(NormalCurvatureError::UnsupportedVariant {
            variant: "Rectangular"
        })
    );
}
