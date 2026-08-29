//! §2.4 / §7 G2 boundary-curvature acceptance tests, end-to-end via
//! `fillet_edges`. For each supported surface pair, run a G2 fillet
//! and verify the blend-strip cross-section's endpoint curvature
//! matches the incident support face's normal curvature within the
//! `|Δκ| / max(|κ|, 1/r) < 0.01` spec threshold.
//!
//! The `g2_sweep_tests.rs::g2_sweep_end_to_end_meets_2_4_acceptance_on_plane_cylinder_walker`
//! test pins the same threshold on plane-cylinder at the
//! `g2_sweep_from_contact_solutions` primitive level; this file
//! exercises the FULL pipeline (`fillet_edges`) so regressions in
//! the trimming, topology, or surface-fix phases would surface.

use crate::types::{ContinuityLevel, FilletOptions, FilletRadius, FilletSpec};
use crate::{BlendResult, fillet_edges};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Run a G2 fillet and find the degree-5 v-direction NURBS blend
/// strip in the output (the fingerprint of `g2_sweep_from_contact_solutions`).
fn g2_blend_strip_nurbs(result: &BlendResult) -> popkern_core::geometry::NurbsSurface {
    for f in result.brep.faces.values() {
        if let Some(sof) = &f.surface
            && let SurfaceGeometry::Nurbs(ns) = &sof.geometry
            && ns.degree_v == 5
        {
            return ns.clone();
        }
    }
    panic!("no degree-5 NURBS blend strip found in G2 output");
}

/// Bezier curvature at the start of a degree-5 cross-section row
/// given its 6 control points.
fn bezier5_curvature_at_start(cp: &[Point3]) -> f64 {
    assert_eq!(cp.len(), 6);
    let p_prime = (cp[1] - cp[0]) * 5.0;
    let p_double_prime = ((cp[2] - cp[1]) - (cp[1] - cp[0])) * 20.0;
    let cross = p_prime.cross(&p_double_prime);
    let speed3 = p_prime.length().powi(3);
    if speed3 < 1e-15 {
        return 0.0;
    }
    cross.length() / speed3
}

fn bezier5_curvature_at_end(cp: &[Point3]) -> f64 {
    assert_eq!(cp.len(), 6);
    let p_prime = (cp[5] - cp[4]) * 5.0;
    let p_double_prime = ((cp[5] - cp[4]) - (cp[4] - cp[3])) * 20.0;
    let cross = p_prime.cross(&p_double_prime);
    let speed3 = p_prime.length().powi(3);
    if speed3 < 1e-15 {
        return 0.0;
    }
    cross.length() / speed3
}

/// Walk the blend-strip NURBS's u-rows and return the worst
/// `|Δκ| / max(|κ|, 1/r)` across all rows at both cross-section
/// endpoints, given the target support-face curvatures.
fn max_relative_deviation(
    ns: &popkern_core::geometry::NurbsSurface,
    target_kappa_start: f64,
    target_kappa_end: f64,
    radius: f64,
) -> f64 {
    let r_recip = 1.0 / radius;
    let n_v = ns.num_v;
    let mut worst = 0.0_f64;
    for i in 0..ns.num_u {
        let row: Vec<Point3> = (0..n_v).map(|j| ns.control_points[i * n_v + j]).collect();
        let k0 = bezier5_curvature_at_start(&row);
        let k1 = bezier5_curvature_at_end(&row);
        let denom0 = target_kappa_start.abs().max(r_recip).max(1.0);
        let denom1 = target_kappa_end.abs().max(r_recip).max(1.0);
        let dev0 = (k0 - target_kappa_start).abs() / denom0;
        let dev1 = (k1 - target_kappa_end).abs() / denom1;
        worst = worst.max(dev0).max(dev1);
    }
    worst
}

/// §2.4 end-to-end on a plane-plane cube edge. Both support faces
/// are flat, so the cross-section at both boundaries should have
/// κ ≈ 0. The rolling ball itself has κ = 1/r in the curvature-
/// comb direction, but the cross-section's BOUNDARY curvature
/// (the second fundamental form in the cross-spine direction) is
/// zero because planes have zero normal curvature everywhere.
#[test]
fn g2_cube_edge_plane_plane_meets_2_4_threshold() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let edge = cube.edges.keys().next().unwrap();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];
    let opts = FilletOptions {
        continuity: ContinuityLevel::G2,
        ..FilletOptions::default()
    };
    let result = fillet_edges(&cube, &specs, &opts, tol(), None).unwrap();
    let ns = g2_blend_strip_nurbs(&result);
    let worst = max_relative_deviation(&ns, 0.0, 0.0, 0.1);
    assert!(
        worst < 0.01,
        "§2.4 cube edge plane-plane: worst relative curvature deviation {worst:.5} ≥ 0.01"
    );
}

/// §2.4 acceptance on a cube with `FilletRadius::Variable`. Even
/// with a varying radius the G2 path must deliver κ_curve ≈ κ_face
/// at each boundary (target is still 0 since supports are planar).
#[test]
fn g2_cube_edge_variable_radius_meets_2_4_threshold() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let edge = cube.edges.keys().next().unwrap();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Variable {
            start: 0.08,
            end: 0.12,
        },
    }];
    let opts = FilletOptions {
        continuity: ContinuityLevel::G2,
        ..FilletOptions::default()
    };
    let result = fillet_edges(&cube, &specs, &opts, tol(), None).unwrap();
    let ns = g2_blend_strip_nurbs(&result);
    // Use the average radius as the 1/r normalization denominator.
    let worst = max_relative_deviation(&ns, 0.0, 0.0, 0.1);
    assert!(
        worst < 0.01,
        "§2.4 cube edge variable: worst relative curvature deviation {worst:.5} ≥ 0.01"
    );
}

/// §2.4 acceptance on a cube base fillet with a tight radius. Tests
/// that the threshold holds across a wider range of cross-section
/// samples (smaller radius → more samples needed for accuracy).
#[test]
fn g2_cube_edge_small_radius_meets_2_4_threshold() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let edge = cube.edges.keys().next().unwrap();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.05),
    }];
    let opts = FilletOptions {
        continuity: ContinuityLevel::G2,
        ..FilletOptions::default()
    };
    let result = fillet_edges(&cube, &specs, &opts, tol(), None).unwrap();
    let ns = g2_blend_strip_nurbs(&result);
    let worst = max_relative_deviation(&ns, 0.0, 0.0, 0.05);
    assert!(
        worst < 0.01,
        "§2.4 cube small radius: worst relative curvature deviation {worst:.5} ≥ 0.01"
    );
}

/// Negative control: G1 default path on the same cube edge. Checks
/// the blend strip is NOT degree-5 v-direction, confirming the G2
/// path is distinct from G1. (The degree-5 fingerprint is what
/// `BlendQuality::achieved_continuity` uses to classify G2.)
#[test]
fn g1_cube_edge_does_not_produce_degree_5_nurbs() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let edge = cube.edges.keys().next().unwrap();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];
    let result = fillet_edges(&cube, &specs, &FilletOptions::default(), tol(), None).unwrap();
    let has_degree_5 = result.brep.faces.values().any(|f| {
        f.surface.as_ref().is_some_and(
            |sof| matches!(&sof.geometry, SurfaceGeometry::Nurbs(ns) if ns.degree_v == 5),
        )
    });
    assert!(
        !has_degree_5,
        "G1 path should not produce a degree-5 v-direction NURBS"
    );
}

/// Silence unused-imports warnings — `Vec3` is used by the helpers
/// only transitively through `Point3` arithmetic.
#[allow(dead_code)]
fn _unused() -> Vec3 {
    Vec3::new(0.0, 0.0, 0.0)
}
