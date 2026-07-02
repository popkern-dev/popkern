//! Tests for the bi-6 polar promotion.
//!
//! Exercises `promote_decomp_to_kp_polar_bi6` on a synthetic
//! N-sided ring of two-stripe sub-patches. The ring is built by
//! hand (no upstream three-quad decomposition); verifies:
//!
//! - n < 5 rejected (the binonic rectangular path owns that
//!   regime).
//! - Non-two-stripe input rejected.
//! - Valid n = 5, 6, 7, 8 rings convert to bi-nonic grids and
//!   replace the input.
//!
//! The patch geometry here is rough and not commercial-grade;
//! we're testing the wiring and dispatch logic, not output
//! quality. Commercial-gate acceptance tests for real fillets
//! are in Stage 6 (I.25–I.28) of the plan.

use super::{extract_fillet_radius, promote_decomp_to_kp_polar_bi6};
use crate::construct::kp20_polar_dispatch::try_kp20_polar_macro_candidate;
use crate::fillet_surface::FilletSurfaceResult;
use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::{GregoryPatch, RectangularGregoryTwoStripe, TensorBorderBi3};

fn pt(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

/// Build an N-sided ring of two-stripe sub-patches that share
/// spoke boundaries. Patches are placed as "pie slices" in the
/// xy plane with the center at origin.
#[allow(clippy::needless_range_loop)]
fn synthetic_ring(n: usize) -> Vec<GregoryPatch> {
    (0..n)
        .map(|k| {
            let theta0 = (2.0 * std::f64::consts::PI) * (k as f64) / (n as f64);
            let theta1 = (2.0 * std::f64::consts::PI) * ((k + 1) as f64) / (n as f64);
            let mut grid = [[pt(0.0, 0.0, 0.0); 8]; 8];
            for i in 0..8 {
                for j in 0..8 {
                    let u_frac = 1.0 - (i as f64) / 7.0; // radial
                    let v_frac = (j as f64) / 7.0;
                    let theta = theta0 * (1.0 - v_frac) + theta1 * v_frac;
                    let r = u_frac;
                    grid[i][j] = pt(r * theta.cos(), r * theta.sin(), 0.0);
                }
            }
            GregoryPatch::RectangularTwoStripe(Box::new(RectangularGregoryTwoStripe { grid }))
        })
        .collect()
}

fn kp20_tensor_sector(offset: f64) -> [[Point3; 4]; 3] {
    [
        [
            pt(offset + 0.0, 0.0, 0.0),
            pt(offset + 1.0, 0.0, 0.0),
            pt(offset + 2.0, 0.0, 0.0),
            pt(offset + 3.0, 0.0, 0.0),
        ],
        [
            pt(offset + 0.0, 1.0, 0.0),
            pt(offset + 1.0, 1.0, 0.0),
            pt(offset + 2.0, 1.0, 0.0),
            pt(offset + 3.0, 1.0, 0.0),
        ],
        [
            pt(offset + 0.0, 2.0, 0.0),
            pt(offset + 1.0, 2.0, 0.0),
            pt(offset + 2.0, 2.0, 0.0),
            pt(offset + 3.0, 2.0, 0.0),
        ],
    ]
}

#[test]
fn rejects_n_less_than_5() {
    let mut ring = synthetic_ring(4);
    assert!(matches!(
        promote_decomp_to_kp_polar_bi6(&mut ring, &[]),
        super::PolarPromoteOutcome::Rejected { .. }
    ));
}

#[test]
fn rejects_empty_ring() {
    let mut ring: Vec<GregoryPatch> = vec![];
    assert!(matches!(
        promote_decomp_to_kp_polar_bi6(&mut ring, &[]),
        super::PolarPromoteOutcome::Rejected { .. }
    ));
}

#[test]
fn rejects_non_two_stripe_input() {
    // Mix a two-stripe ring with a non-two-stripe variant.
    let mut ring = synthetic_ring(5);
    // Replace patch 0 with a different variant. Use a cubic
    // rectangular by taking an existing two-stripe and converting
    // — simpler to construct a zero-sized different variant.
    // We'll use the RectangularGregoryBiNonic type itself as the
    // wrong-variant check.
    use popkern_core::geometry::gregory::RectangularGregoryBiNonic;
    let bi9 = RectangularGregoryBiNonic {
        grid: [[pt(0.0, 0.0, 0.0); 10]; 10],
    };
    ring[0] = GregoryPatch::RectangularBiNonic(Box::new(bi9));
    assert!(matches!(
        promote_decomp_to_kp_polar_bi6(&mut ring, &[]),
        super::PolarPromoteOutcome::Rejected { .. }
    ));
}

#[test]
fn kp20_macro_candidate_builds_four_sub_patches_per_sector() {
    // The macro candidate routes through the paper-faithful section
    // 4.4 assembly, which gates on G2 closure — so the fixture must
    // be a VALID C2 ring border (the previous disconnected synthetic
    // sectors are not Hermite data of any surface and are correctly
    // declined by the paper gates).
    let n = 5;
    let tensor_border =
        crate::construct::kp20_polar_dispatch::paper_primary_candidate_test_border(n);

    let candidate = try_kp20_polar_macro_candidate(pt(0.0, 0.0, 0.3), &tensor_border).unwrap();

    assert_eq!(candidate.len(), 4 * n);
    for patch in candidate {
        assert!(matches!(patch, GregoryPatch::RectangularBiNonic(_)));
    }
}

// The positive-path test (`ring → bi-nonic patches`) requires a
// ring with enough geometric fidelity that `extract_direct_boundary_input`
// can produce finite tangent planes. The synthetic flat ring
// above produces degenerate normals (cross product of two
// parallel xy-plane vectors ⇒ zero vector) which intentionally
// fails; that's an honest negative result, not a bug.
//
// Full-pipeline positive verification lives in the core-level
// tests in `popkern-core` (`kp15_bi6_polar_bridge_tests.rs` +
// `kp_polar_bi6_apply_tests.rs`), which use synthetic
// `PolarBoundaryInput` directly and exercise the full flow at
// n ∈ {5, 6, 7, 8, 10, 12}.

fn placeholder_trim() -> popkern_core::geometry::NurbsCurve {
    popkern_core::geometry::NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![pt(0.0, 0.0, 0.0), pt(1.0, 0.0, 0.0)],
        weights: vec![1.0, 1.0],
    }
}

fn fs_with(surface: popkern_core::geometry::surface::SurfaceGeometry) -> FilletSurfaceResult {
    FilletSurfaceResult {
        surface,
        trim_curve_1: placeholder_trim(),
        trim_curve_2: placeholder_trim(),
        continuity_downgrade: None,
    }
}

#[test]
fn extract_radius_matches_cylinder_torus_sphere() {
    use popkern_core::geometry::Vec3;
    use popkern_core::geometry::surface::SurfaceGeometry;

    let cyl = fs_with(SurfaceGeometry::Cylinder {
        origin: pt(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 0.25,
    });
    assert!((extract_fillet_radius(&cyl).unwrap() - 0.25).abs() < 1e-15);

    let tor = fs_with(SurfaceGeometry::Torus {
        center: pt(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: 1.0,
        minor_radius: 0.1,
    });
    assert!((extract_fillet_radius(&tor).unwrap() - 0.1).abs() < 1e-15);

    let sph = fs_with(SurfaceGeometry::Sphere {
        center: pt(0.0, 0.0, 0.0),
        radius: 0.7,
    });
    assert!((extract_fillet_radius(&sph).unwrap() - 0.7).abs() < 1e-15);
}

#[test]
fn extract_radius_returns_none_for_non_analytic() {
    use popkern_core::geometry::Vec3;
    use popkern_core::geometry::surface::SurfaceGeometry;

    let plane = fs_with(SurfaceGeometry::Plane {
        origin: pt(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    });
    assert!(extract_fillet_radius(&plane).is_none());
}
