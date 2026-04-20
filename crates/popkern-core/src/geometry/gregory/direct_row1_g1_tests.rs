//! Tests for [`super::solve_single_spoke_row1_g1`].
//!
//! The decisive test — **"non-G1 setup becomes G1 to machine
//! precision after one solve"** — verifies dense 100-sample
//! min-cos ≥ 0.9999 after calling the direct formula, which
//! is the claim the whole Phase G infrastructure exists to
//! justify.

#![allow(clippy::needless_range_loop)]

use super::*;
use crate::geometry::Point3;
use crate::geometry::gregory::RectangularGregoryTwoStripe;
use crate::geometry::gregory::direct_blend_coeffs::compute_ring_blend_coeffs;
use crate::geometry::gregory::direct_input::{
    DirectBoundaryInput, OuterCornerData, SpokeBoundaryData, TangentPlane,
};
use crate::geometry::gregory::direct_options::DirectConstructionOptions;

/// Flat unit-square two-stripe patch: grid[i][j] = origin +
/// (i/7)·u_dir + (j/7)·v_dir on z=0.
fn flat_patch(
    origin: (f64, f64),
    u_dir: (f64, f64),
    v_dir: (f64, f64),
) -> RectangularGregoryTwoStripe {
    let mut grid = [[Point3::new(0.0, 0.0, 0.0); 8]; 8];
    for i in 0..8 {
        for j in 0..8 {
            let s = (i as f64) / 7.0;
            let t = (j as f64) / 7.0;
            grid[i][j] = Point3::new(
                origin.0 + s * u_dir.0 + t * v_dir.0,
                origin.1 + s * u_dir.1 + t * v_dir.1,
                0.0,
            );
        }
    }
    RectangularGregoryTwoStripe { grid }
}

/// Package a single-spoke test setup: two adjacent flat
/// patches sharing x=0 from (0, 0, 0) to (0, 1, 0). Patch A
/// occupies [-1, 0] × [0, 1]; patch B [0, 1] × [0, 1]. B's
/// UZero boundary (u=0) is the shared spoke (not VZero as in
/// three-quad — this simplifies the test; the math is identical).
fn flat_spoke_setup(
    b_perturb: f64,
) -> (
    RectangularGregoryTwoStripe,
    RectangularGregoryTwoStripe,
    SingleSpokeRow1Input,
) {
    let a = flat_patch((-1.0, 0.0), (1.0, 0.0), (0.0, 1.0));
    let mut b = flat_patch((0.0, 0.0), (1.0, 0.0), (0.0, 1.0));
    // Break G1 by lifting B's row-1-near-UZero free CPs in z.
    for j in 1..=6 {
        b.grid[1][j].z += b_perturb;
    }
    // Package spoke data. For this test we use A.UOne ↔ B.UZero
    // NOT reversed (simpler than three-quad's VZero reversed
    // case). Spoke boundary_cps = A's grid[7][0..7] = B's grid[0][0..7].
    let mut boundary_cps = [Point3::new(0.0, 0.0, 0.0); 8];
    boundary_cps.copy_from_slice(&a.grid[7]);
    let spoke = SpokeBoundaryData {
        boundary_cps,
        // Pre-perturb xtans (the "G1-valid" state used for blend-coeff computation).
        // On flat coplanar patches: xtan_A = (1, 0, 0), xtan_B = (1, 0, 0).
        xtan_a_at_c: Vec3::new(1.0, 0.0, 0.0),
        xtan_b_at_c: Vec3::new(1.0, 0.0, 0.0),
        xtan_a_at_v: Vec3::new(1.0, 0.0, 0.0),
        xtan_b_at_v: Vec3::new(1.0, 0.0, 0.0),
        radius_at_c: 0.1,
        radius_at_v: 0.1,
        edge_id: "test-spoke".into(),
    };

    // Row-1 CPs in spoke-t order. For A's UOne, row-1 = grid[6][j] for j=0..7 — already in spoke-t order.
    // For B's UZero (NOT reversed), row-1 = grid[1][j] for j=0..7 — also in spoke-t order.
    let mut patch_a_row1 = [Point3::new(0.0, 0.0, 0.0); 8];
    let mut patch_b_row1 = [Point3::new(0.0, 0.0, 0.0); 8];
    patch_a_row1.copy_from_slice(&a.grid[6]);
    patch_b_row1.copy_from_slice(&b.grid[1]);

    // Compute blend coefficients by constructing a synthetic
    // DirectBoundaryInput with the same spoke (the blend-coeff
    // solver only uses spoke + xtans + options).
    let c = Point3::new(0.0, 0.0, 0.0);
    let v = Point3::new(0.0, 1.0, 0.0);
    let tps_xy = TangentPlane::from_normal(Vec3::new(0.0, 0.0, 1.0)).unwrap();
    let input_for_blend = DirectBoundaryInput {
        central_vertex: c,
        spokes: vec![spoke.clone(), spoke.clone(), spoke.clone()],
        outer_corners: vec![
            OuterCornerData {
                position: v,
                tangent_plane: tps_xy,
            };
            3
        ],
        central_vertex_tps: vec![tps_xy; 3],
        seed_patches: None,
    };
    let options = DirectConstructionOptions::with_characteristic_length(1.0);
    let coeffs = compute_ring_blend_coeffs(&input_for_blend, &options).expect("blend coeffs")[0];

    let input = SingleSpokeRow1Input {
        spoke,
        blend_coeffs: coeffs,
        patch_a_row1,
        patch_b_row1,
    };
    (a, b, input)
}

/// Apply the solve output back to the two patches' grids.
/// For the simplified A.UOne ↔ B.UZero (not reversed) case:
///   A.grid[6][j] ← output.patch_a_row1_free[j-1] for j=1..=6
///   B.grid[1][j] ← output.patch_b_row1_free[j-1] for j=1..=6.
fn apply_output(
    a: &mut RectangularGregoryTwoStripe,
    b: &mut RectangularGregoryTwoStripe,
    output: &SingleSpokeRow1Output,
) {
    for (idx, j) in (1..=6).enumerate() {
        a.grid[6][j] = output.patch_a_row1_free[idx];
        b.grid[1][j] = output.patch_b_row1_free[idx];
    }
}

/// Dense 100-sample G1 min-cos along the shared spoke.
fn dense_min_cos(
    a: &RectangularGregoryTwoStripe,
    b: &RectangularGregoryTwoStripe,
    samples: usize,
) -> f64 {
    let mut min_cos = 1.0_f64;
    for k in 1..=samples {
        let t = (k as f64) / ((samples + 1) as f64);
        let na = a.partial_u(1.0, t).cross(&a.partial_v(1.0, t)).normalized();
        let nb = b.partial_u(0.0, t).cross(&b.partial_v(0.0, t)).normalized();
        let (Some(na), Some(nb)) = (na, nb) else {
            return 0.0;
        };
        let cos = (na.x * nb.x + na.y * nb.y + na.z * nb.z).abs();
        if cos < min_cos {
            min_cos = cos;
        }
    }
    min_cos
}

#[test]
fn flat_coplanar_setup_is_already_g1_and_solve_is_a_no_op() {
    let (a, b, input) = flat_spoke_setup(0.0);
    let options = DirectConstructionOptions::with_characteristic_length(1.0);
    let output = solve_single_spoke_row1_g1(&input, &options).expect("solved");
    // Residual ≈ 0 on flat coplanar (already G1).
    assert!(
        output.residual_max < 1e-9,
        "flat residual = {}",
        output.residual_max
    );
    // Applied output shouldn't visibly move CPs.
    let mut a2 = a.clone();
    let mut b2 = b.clone();
    apply_output(&mut a2, &mut b2, &output);
    let min_cos = dense_min_cos(&a2, &b2, 100);
    assert!(
        (min_cos - 1.0).abs() < 1e-9,
        "post-solve min_cos = {min_cos}"
    );
}

#[test]
fn non_g1_setup_recovers_exact_g1_after_single_solve() {
    // Perturb B's row-1 CPs in z → non-G1 initial state.
    let (a, b, input) = flat_spoke_setup(0.2);

    // Before the solve: non-G1.
    let before = dense_min_cos(&a, &b, 50);
    assert!(
        before < 0.99,
        "test setup should be non-G1 before solve (min_cos = {before})"
    );

    let options = DirectConstructionOptions::with_characteristic_length(1.0);
    let output = solve_single_spoke_row1_g1(&input, &options).expect("solved");
    assert!(
        output.residual_max < 1e-9,
        "residual = {} after solve",
        output.residual_max
    );

    let mut a2 = a.clone();
    let mut b2 = b.clone();
    apply_output(&mut a2, &mut b2, &output);

    // The decisive assertion: dense 100-sample G1 min-cos ≈ 1.0
    // to machine precision — this is the claim the whole Phase
    // G infrastructure exists to establish.
    let min_cos = dense_min_cos(&a2, &b2, 100);
    assert!(
        (min_cos - 1.0).abs() < 1e-9,
        "post-solve dense min_cos = {min_cos}, want 1.0 to machine precision"
    );
}

#[test]
fn solved_cps_stay_finite_on_asymmetric_perturbation() {
    // Larger perturbation — make sure the solver doesn't diverge.
    let (_, _, input) = flat_spoke_setup(0.5);
    let options = DirectConstructionOptions::with_characteristic_length(1.0);
    let output = solve_single_spoke_row1_g1(&input, &options).expect("solved");
    for cp in output
        .patch_a_row1_free
        .iter()
        .chain(output.patch_b_row1_free.iter())
    {
        assert!(cp.x.is_finite() && cp.y.is_finite() && cp.z.is_finite());
    }
}

#[test]
fn residual_is_bounded_even_on_larger_perturbation() {
    let (_, _, input) = flat_spoke_setup(1.0);
    let options = DirectConstructionOptions::with_characteristic_length(1.0);
    let output = solve_single_spoke_row1_g1(&input, &options).expect("solved");
    // For flat coplanar patches (with ∑ tangent = const) the
    // system is CONSISTENT regardless of perturbation magnitude.
    // The residual should still vanish.
    assert!(
        output.residual_max < 1e-8,
        "large perturbation residual = {}",
        output.residual_max
    );
}
