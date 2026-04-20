//! Tests for the single-spoke polynomial-identity G1 solver.
//!
//! The decisive tests are the **dense-sampling G1 checks**:
//! after the solve, we evaluate the surface normals on either
//! side of the shared spoke at 100 parameter values and assert
//! that every pair is aligned to machine precision. That's the
//! structural claim the solver makes — not merely at samples the
//! solver saw during optimisation, but at *every* `t ∈ [0, 1]`
//! because the underlying polynomial identity zeroes every
//! Bernstein coefficient of the LHS.

use super::*;
use crate::geometry::{Point3, Vec3};

/// A flat unit square on the `z = 0` plane. Trivially G1-smooth
/// with any neighbour that shares its boundary.
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

/// Unit normal at `(u, v)` on the patch via the cross of the
/// partial derivatives.
fn unit_normal(p: &RectangularGregoryTwoStripe, u: f64, v: f64) -> Option<Vec3> {
    p.partial_u(u, v).cross(&p.partial_v(u, v)).normalized()
}

/// Min |n_A · n_B| over a dense set of parameters on the spoke.
/// The test's ultimate witness — should be ≈ 1.0 after the solve.
fn spoke_min_cos(
    patch_a: &RectangularGregoryTwoStripe,
    side_a: Side,
    patch_b: &RectangularGregoryTwoStripe,
    side_b: Side,
    b_reversed: bool,
    n_samples: usize,
) -> f64 {
    let mut min_cos = 1.0_f64;
    for k in 1..=n_samples {
        let t = (k as f64) / ((n_samples + 1) as f64);
        let (ua, va) = uv_for_side(side_a, t);
        let t_b = if b_reversed { 1.0 - t } else { t };
        let (ub, vb) = uv_for_side(side_b, t_b);
        // Strict: if a normal is ill-defined (degenerate patch
        // with a zero partial), return 0.0. This guards against
        // a false "G1" report when the solver collapses row-1
        // CPs onto the boundary.
        let Some(na) = unit_normal(patch_a, ua, va) else {
            return 0.0;
        };
        let Some(nb) = unit_normal(patch_b, ub, vb) else {
            return 0.0;
        };
        let cos = (na.x * nb.x + na.y * nb.y + na.z * nb.z).abs();
        if cos < min_cos {
            min_cos = cos;
        }
    }
    min_cos
}

fn uv_for_side(side: Side, t: f64) -> (f64, f64) {
    match side {
        Side::UZero => (0.0, t),
        Side::UOne => (1.0, t),
        Side::VZero => (t, 0.0),
        Side::VOne => (t, 1.0),
    }
}

#[test]
fn xtan_bernstein_cps_for_flat_patch() {
    // On the flat unit square with u_dir = (1, 0), v_dir = (0, 1),
    // xtan at UZero (= ∂/∂u at u = 0) points in +x everywhere.
    // Bernstein coefficient magnitude: 7·(P[1][j] − P[0][j]) =
    // 7·(1/7, 0, 0) = (1, 0, 0).
    let p = flat_patch((0.0, 0.0), (1.0, 0.0), (0.0, 1.0));
    let coeffs = xtan_bernstein_cps(&p, Side::UZero);
    for c in coeffs.iter() {
        assert!((c.x - 1.0).abs() < 1e-12, "x = {}", c.x);
        assert!(c.y.abs() < 1e-12);
        assert!(c.z.abs() < 1e-12);
    }
}

#[test]
fn boundary_tangent_poly_for_linear_boundary() {
    // A linear boundary has constant tangent. Verify that
    // every Bernstein coefficient of T(t) equals the step size.
    let p = flat_patch((0.0, 0.0), (1.0, 0.0), (0.0, 1.0));
    let bcps = boundary_bernstein_cps(&p, Side::UZero);
    let t_poly = boundary_tangent_poly(&bcps);
    assert_eq!(t_poly.degree(), 6);
    // Boundary along v at u=0: points are (0, j/7, 0). C'(v) = (0, 1, 0).
    // Every Bernstein coefficient of T = (0, 1, 0).
    for j in 0..=t_poly.degree() {
        let c = t_poly.coeff(j);
        assert!(c.x.abs() < 1e-12);
        assert!((c.y - 1.0).abs() < 1e-12);
        assert!(c.z.abs() < 1e-12);
    }
}

#[test]
fn flat_adjacent_patches_solve_to_zero_residual() {
    // Two coplanar flat patches sharing the line x = 0 in the
    // z = 0 plane. Patch A occupies x ∈ [−1, 0], patch B
    // occupies x ∈ [0, 1]. Both are already G1-compatible, so
    // the solve should return residual ≈ 0 with β ≈ 0.
    //
    //   A: u_dir = (1, 0), v_dir = (0, 1), origin (−1, 0)
    //      → spoke is A's UOne (u=1, x=0)
    //   B: u_dir = (1, 0), v_dir = (0, 1), origin (0, 0)
    //      → spoke is B's UZero (u=0, x=0)
    //   b_reversed = false (both param v → y increasing).
    let a = flat_patch((-1.0, 0.0), (1.0, 0.0), (0.0, 1.0));
    let b = flat_patch((0.0, 0.0), (1.0, 0.0), (0.0, 1.0));
    let result =
        solve_single_spoke_poly_g1(&a, Side::UOne, &b, Side::UZero, false).expect("solved");
    assert!(
        result.residual_max < 1e-10,
        "flat residual should be zero, got {}",
        result.residual_max
    );
}

#[test]
fn non_g1_patches_become_g1_after_solve() {
    // A curved non-G1 setup: patch A on x ≤ 0 is flat, patch B
    // on x ≥ 0 has its row 1 lifted in z. That breaks G1 along
    // x = 0. The polynomial-identity solve adjusts BOTH
    // patches' row-1 CPs + β so that G1 holds everywhere.
    let a = flat_patch((-1.0, 0.0), (1.0, 0.0), (0.0, 1.0));
    let mut b = flat_patch((0.0, 0.0), (1.0, 0.0), (0.0, 1.0));
    // Break G1: shift B's row i = 1 up by 0.1 in z. This is B's
    // near-spoke row for side UZero.
    for j in 0..8 {
        b.grid[1][j].z += 0.1;
    }

    // Before: spoke min-cos < 1.
    let before = spoke_min_cos(&a, Side::UOne, &b, Side::UZero, false, 20);
    assert!(
        before < 0.999,
        "test setup should be non-G1 before solve (min_cos = {before})"
    );

    let result =
        solve_single_spoke_poly_g1(&a, Side::UOne, &b, Side::UZero, false).expect("solved");
    assert!(
        result.residual_max < 1e-9,
        "residual should vanish (got {})",
        result.residual_max
    );

    let mut a2 = a.clone();
    let mut b2 = b.clone();
    apply_single_spoke_solve(&mut a2, &mut b2, Side::UOne, Side::UZero, false, &result);

    // Dense 100-sample min-cos: should be 1.0 to machine precision.
    let min_cos = spoke_min_cos(&a2, Side::UOne, &b2, Side::UZero, false, 100);
    assert!(
        (min_cos - 1.0).abs() < 1e-9,
        "post-solve dense min_cos = {min_cos}, want 1.0 to machine precision"
    );
}

#[test]
fn solver_handles_reversed_orientation() {
    // Same setup but pretend B's spoke parameter runs in the
    // opposite direction. Even with `b_reversed = true`, for
    // patches that are genuinely coherent (and reversing doesn't
    // apply geometrically), the solve still has to yield a
    // consistent result — the identity is linear and the
    // coefficient reversal is a no-op when the geometry is
    // already symmetric.
    let a = flat_patch((-1.0, 0.0), (1.0, 0.0), (0.0, 1.0));
    let b = flat_patch((0.0, 0.0), (1.0, 0.0), (0.0, 1.0));
    // This isn't a physically-meaningful reversed spoke, but it
    // exercises the coefficient-reversal code path.
    let result = solve_single_spoke_poly_g1(&a, Side::UOne, &b, Side::UZero, true).expect("solved");
    // The system is still consistent (the coefficient reversal
    // is symmetric under the flat coplanar geometry), so residual
    // ≈ 0 is still reachable — we just verify the solver doesn't
    // panic or return None.
    assert!(result.residual_max.is_finite());
}

#[test]
fn g1_holds_with_non_trivial_boundary_curvature() {
    // A curved shared boundary (not a straight line). The
    // boundary tangent T(t) is no longer constant, so β·T is a
    // non-trivial vector polynomial — this exercises the
    // Bernstein-product code path.
    let mut a = flat_patch((-1.0, 0.0), (1.0, 0.0), (0.0, 1.0));
    let mut b = flat_patch((0.0, 0.0), (1.0, 0.0), (0.0, 1.0));
    // Bend the shared spoke: shift boundary CPs in z along a
    // sinusoidal profile. The boundary is the j-th column of
    // each patch's A-UOne / B-UZero sides, which coincide on the
    // same x = 0 line — must be set identically in both patches.
    for j in 0..8 {
        let tj = (j as f64) / 7.0;
        let z_shift = 0.2 * (std::f64::consts::PI * tj).sin();
        a.grid[7][j].z = z_shift;
        b.grid[0][j].z = z_shift;
    }
    // Also perturb row 1 on B to break G1.
    for j in 0..8 {
        b.grid[1][j].z += 0.15;
    }

    let before = spoke_min_cos(&a, Side::UOne, &b, Side::UZero, false, 50);
    assert!(
        before < 0.99,
        "setup should be non-G1 before solve (min_cos = {before})"
    );

    let result =
        solve_single_spoke_poly_g1(&a, Side::UOne, &b, Side::UZero, false).expect("solved");
    assert!(
        result.residual_max < 1e-9,
        "residual should vanish (got {})",
        result.residual_max
    );

    let mut a2 = a.clone();
    let mut b2 = b.clone();
    apply_single_spoke_solve(&mut a2, &mut b2, Side::UOne, Side::UZero, false, &result);

    let min_cos = spoke_min_cos(&a2, Side::UOne, &b2, Side::UZero, false, 100);
    assert!(
        (min_cos - 1.0).abs() < 1e-9,
        "curved-spoke post-solve min_cos = {min_cos}, want 1.0"
    );
}

#[test]
fn boundary_cps_unchanged_by_solve() {
    // The solver must NOT move the boundary CPs — only row-1
    // near-spoke CPs plus the (internal) β coefficients change.
    let a = flat_patch((-1.0, 0.0), (1.0, 0.0), (0.0, 1.0));
    let mut b = flat_patch((0.0, 0.0), (1.0, 0.0), (0.0, 1.0));
    for j in 0..8 {
        b.grid[1][j].z += 0.1;
    }

    let a_boundary_before: Vec<Point3> = (0..8).map(|j| a.grid[7][j]).collect();
    let b_boundary_before: Vec<Point3> = (0..8).map(|j| b.grid[0][j]).collect();

    let result =
        solve_single_spoke_poly_g1(&a, Side::UOne, &b, Side::UZero, false).expect("solved");
    let mut a2 = a.clone();
    let mut b2 = b.clone();
    apply_single_spoke_solve(&mut a2, &mut b2, Side::UOne, Side::UZero, false, &result);

    for j in 0..8 {
        let a_after = a2.grid[7][j];
        let b_after = b2.grid[0][j];
        assert!(
            (a_after.x - a_boundary_before[j].x).abs() < 1e-12
                && (a_after.y - a_boundary_before[j].y).abs() < 1e-12
                && (a_after.z - a_boundary_before[j].z).abs() < 1e-12,
            "A boundary at j={j} moved"
        );
        assert!(
            (b_after.x - b_boundary_before[j].x).abs() < 1e-12
                && (b_after.y - b_boundary_before[j].y).abs() < 1e-12
                && (b_after.z - b_boundary_before[j].z).abs() < 1e-12,
            "B boundary at j={j} moved"
        );
    }
}
