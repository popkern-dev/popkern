//! PK-225: N=3 setback vertex blends are production for corners the
//! sphere does not own.
//!
//! Three root fixes landed the proven KP20 setback fan on trihedral
//! cube corners (measured on the uniform conic corner, r=0.15
//! rho=0.5):
//! - crossings graze-merge: a setback arc that TERMINATES on a rail
//!   approaches it near-tangentially and the polyline test read the
//!   tangent wiggle as extra crossings (same-edge grazes at 2.1e-6
//!   next to exact endpoint contacts at 1e-17);
//! - rim-refit row pick by endpoint-corner identity: midpoint
//!   proximity was ambiguous on wide 120-degree pieces and moved a
//!   spoke-mid vertex 4.31e-3 onto the WRONG grid row (fan seams
//!   g0 2.3e-2 / g1 1.0 — the whole historical N=3 decline class);
//! - tangent-coherent corner-pair selection: on a trihedral corner
//!   every rail pair is coplanar and TWO pairs intersect exactly; the
//!   creased on-edge pick fed the assembly a 3.6e-2 corner band
//!   mismatch (min regularity 8.7e-5, three folded near-pole pieces),
//!   while the shared-support contact corner closes the assembly at
//!   machine zero (g2 2.7e-12).
//!
//! Honest gates added with the flip: a fan FOLD postcondition (the
//! fast structural validate() skips self-intersection) and a no-new-
//! face-intersection commit gate (flat conic sections below the
//! circular-equivalent rho at small radii genuinely interpenetrate
//! near the shared contact corner once the strips extend to the
//! setback arcs).
//!
//! PK-226 added the fourth root fix: a 3-D margin gate in
//! `arc_boundary_crossings` (the legacy carve's off-surface axis
//! chords projected across the strip crown and the detector snapped
//! the arc midpoint onto the raw cube edge), which flips the
//! variable-radius class from decline to serve for both profiles.

#[path = "kp_gate_helpers.rs"]
mod helpers;

use helpers::{incident_edges, measure_cap_outer_boundary};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle};
use popkern_fillet::fillet_edges;
use popkern_fillet::types::{BlendProfile, FilletOptions, FilletRadius, FilletSpec};

fn cube_corner(
    profile: BlendProfile,
    radius_fn: impl Fn(usize) -> FilletRadius,
) -> popkern_fillet::BlendResult {
    let tol = Tolerance::default();
    let brep = BRep::unit_cube(tol).unwrap();
    let vertex = brep.vertices.keys().next().unwrap();
    let edges: Vec<EdgeHandle> = incident_edges(&brep, vertex);
    let specs: Vec<FilletSpec> = edges
        .iter()
        .enumerate()
        .map(|(k, &e)| FilletSpec {
            edge: e,
            radius: radius_fn(k),
        })
        .collect();
    let options = FilletOptions {
        profile,
        ..FilletOptions::default()
    };
    fillet_edges(&brep, &specs, &options, tol, None).expect("cube corner fillet succeeds")
}

/// The uniform conic cube corner serves the KP20 setback fan at the
/// commercial cap/support bars with a fully valid body.
#[test]
fn uniform_conic_cube_corner_serves_setback_fan_commercially() {
    let result = cube_corner(BlendProfile::Conic { rho: 0.5 }, |_| {
        FilletRadius::Constant(0.15)
    });
    assert!(
        !result.corner_patches.is_empty() && result.corner_patches.iter().all(|cp| cp.setback_cap),
        "uniform conic corner must serve setback caps, got {:?}",
        result
            .corner_patches
            .iter()
            .map(|cp| (cp.method, cp.setback_cap))
            .collect::<Vec<_>>()
    );
    let errors = result.brep.validate_all();
    assert!(
        errors.is_empty(),
        "setback serve must validate clean, got {errors:?}"
    );
    let cap_faces: Vec<_> = result.corner_patches.iter().map(|cp| cp.face).collect();
    let metrics = measure_cap_outer_boundary(&result.brep, &cap_faces, 100);
    assert!(
        metrics.worst_g0_gap <= 5.0e-3,
        "cap_outer_g0 {:.3e} exceeds the commercial gate 5e-3",
        metrics.worst_g0_gap
    );
    assert!(
        metrics.worst_g1_coplanarity <= 5.0e-2,
        "cap_outer_g1 {:.3e} exceeds the commercial gate 5e-2",
        metrics.worst_g1_coplanarity
    );
}

/// Flat conic sections (rho below the circular equivalent) at small
/// radii interpenetrate at the corner once the strips extend to the
/// setback arcs — the commit must decline (no-new-intersection gate)
/// and the serve stays structurally clean on the legacy path.
#[test]
fn flat_rho_small_radius_corner_declines_setback_and_stays_clean() {
    let result = cube_corner(BlendProfile::Conic { rho: 0.35 }, |_| {
        FilletRadius::Constant(0.1)
    });
    assert!(
        result.corner_patches.iter().all(|cp| !cp.setback_cap),
        "flat-rho small-radius corner must decline the setback commit"
    );
    use popkern_core::error::ValidationError;
    for error in result.brep.validate_all() {
        assert!(
            matches!(error, ValidationError::CurveSurfaceDeviation { .. }),
            "declined corner must keep the legacy serve's class, got {error:?}"
        );
    }
}

/// Sphere-owned corners (uniform radii, analytic strips) keep the
/// exact rolling-ball serve — the N=3 setback flip must not touch
/// them.
#[test]
fn circular_uniform_corner_keeps_sphere_serve() {
    let result = cube_corner(BlendProfile::Circular, |_| FilletRadius::Constant(0.15));
    assert!(
        result
            .corner_patches
            .iter()
            .any(|cp| cp.method == popkern_fillet::types::CornerMethod::RollingBallSphere),
        "uniform circular corner must keep the rolling-ball sphere"
    );
    // PK-233: the serving sphere must be the TRUE rolling-ball
    // sphere (radius = strip radius), not a boundary-ring fit. The
    // measured impostor — the arc+bridge hexagon fits a sphere of
    // radius 1.48r at 1.4e-4 — passed the G0 gate and shipped a rim
    // edge 6.2e-2 off the cap before the tangency certification.
    let sphere_radius = result
        .corner_patches
        .iter()
        .find_map(|cp| {
            match result
                .brep
                .faces
                .get(cp.face)
                .and_then(|f| f.surface.as_ref())
                .map(|s| &s.geometry)
            {
                Some(popkern_core::geometry::surface::SurfaceGeometry::Sphere {
                    radius, ..
                }) => Some(*radius),
                _ => None,
            }
        })
        .expect("rolling-ball corner must serve a Sphere face");
    assert!(
        (sphere_radius - 0.15).abs() < 1e-6,
        "serving sphere radius {sphere_radius} must equal the strip radius 0.15"
    );
    assert!(result.brep.validate_all().is_empty());
}

/// Variable-radius corners serve the setback fan at the commercial
/// bars (PK-226). The historical decline was a false UV crossing:
/// closest-point projection folded the legacy carve's off-surface
/// axis chords onto the strip crown, and the detector snapped the
/// arc MIDPOINT onto the raw cube edge (margin 4.9e-2 accepted
/// silently). The 3-D margin gate in `arc_boundary_crossings`
/// rejects the artifact and the existing endpoint completion
/// restores the true t=1 contact.
#[test]
fn variable_radius_corner_serves_setback_fan_commercially() {
    for profile in [BlendProfile::Conic { rho: 0.5 }, BlendProfile::Circular] {
        let result = cube_corner(profile, |k| FilletRadius::Variable {
            start: 0.10 + 0.02 * k as f64,
            end: 0.16 - 0.02 * k as f64,
        });
        assert!(
            !result.corner_patches.is_empty()
                && result.corner_patches.iter().all(|cp| cp.setback_cap),
            "variable {profile:?} corner must serve setback caps, got {:?}",
            result
                .corner_patches
                .iter()
                .map(|cp| (cp.method, cp.setback_cap))
                .collect::<Vec<_>>()
        );
        let errors = result.brep.validate_all();
        assert!(
            errors.is_empty(),
            "variable {profile:?} serve must validate clean, got {errors:?}"
        );
        let cap_faces: Vec<_> = result.corner_patches.iter().map(|cp| cp.face).collect();
        let metrics = measure_cap_outer_boundary(&result.brep, &cap_faces, 100);
        assert!(
            metrics.worst_g0_gap <= 5.0e-3,
            "variable {profile:?} cap_outer_g0 {:.3e} exceeds the commercial gate 5e-3",
            metrics.worst_g0_gap
        );
        assert!(
            metrics.worst_g1_coplanarity <= 5.0e-2,
            "variable {profile:?} cap_outer_g1 {:.3e} exceeds the commercial gate 5e-2",
            metrics.worst_g1_coplanarity
        );
    }
}
