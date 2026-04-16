//! Tests for `blend_diagnostics` — the one-call UI helpers.

use super::*;
use crate::fillet_edges;
use crate::types::{ContinuityLevel, FilletOptions, FilletRadius, FilletSpec};
use popkern_core::geometry::Vec3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Locate the primary blend strip face in a fillet result — the
/// degree-5 v-direction NURBS when G2 was achieved, or the analytic
/// torus/cylinder on the G1 path.
fn find_blend_face(
    brep: &BRep,
    pre_faces: &std::collections::HashSet<popkern_core::topology::FaceHandle>,
) -> popkern_core::topology::FaceHandle {
    use popkern_core::geometry::surface::SurfaceGeometry;
    for (fh, f) in brep.faces.iter() {
        if pre_faces.contains(&fh) {
            continue;
        }
        if let Some(sof) = &f.surface {
            match &sof.geometry {
                SurfaceGeometry::Nurbs(_)
                | SurfaceGeometry::Torus { .. }
                | SurfaceGeometry::Cylinder { .. } => return fh,
                _ => continue,
            }
        }
    }
    panic!("no blend strip face found");
}

/// `blend_curvature_comb_u` on a cube-edge G1 fillet: the output is
/// an analytic Cylinder, for which the `u` parameter is angular
/// (across-spine). The comb along u-midline should read κ ≈ 1/r.
#[test]
fn comb_u_on_g1_cube_fillet_reads_ball_curvature() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let edge = cube.edges.keys().next().unwrap();
    let pre: std::collections::HashSet<_> = cube.faces.keys().collect();
    let radius = 0.15;
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(radius),
    }];
    let result = fillet_edges(&cube, &specs, &FilletOptions::default(), tol(), None).unwrap();
    let blend = find_blend_face(&result.brep, &pre);
    let comb = blend_curvature_comb_u(&result.brep, blend, 12).expect("comb should build");
    assert_eq!(comb.samples.len(), 12);
    let expected = 1.0 / radius;
    // Interior samples only (loose 10% tolerance for the 3-point
    // stencil on an analytic cylinder).
    for s in &comb.samples[3..9] {
        let rel_err = (s.curvature - expected).abs() / expected;
        assert!(
            rel_err < 0.1,
            "cylinder cross-section κ should be ≈ 1/r = {expected}, got {} (rel_err {rel_err:.3})",
            s.curvature
        );
    }
}

/// `blend_curvature_comb_v` on a cube-edge G2 fillet: the output is
/// a degree-5 v-direction NURBS sweep, where `v` IS the cross-
/// section direction. The comb should show non-zero curvature
/// (the cross-section bulges) across the width.
#[test]
fn comb_v_on_g2_cube_fillet_is_populated() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let edge = cube.edges.keys().next().unwrap();
    let pre: std::collections::HashSet<_> = cube.faces.keys().collect();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];
    let opts = FilletOptions {
        continuity: ContinuityLevel::G2,
        ..FilletOptions::default()
    };
    let result = fillet_edges(&cube, &specs, &opts, tol(), None).unwrap();
    let blend = find_blend_face(&result.brep, &pre);
    let comb = blend_curvature_comb_v(&result.brep, blend, 10).expect("comb builds");
    let max = comb.max_curvature().unwrap();
    assert!(
        max > 0.0,
        "G2 cross-section comb must have positive peak curvature, got {max}"
    );
}

/// `blend_highlight_kink_max` on a cube edge produces a finite
/// number (or None if no samples match). A G2 fillet should produce
/// a smaller kink angle than the same-geometry G1 fillet; strict
/// comparison is hard without a controlled reference so we just
/// assert the call returns a finite value.
#[test]
fn highlight_kink_max_returns_finite_on_cube_fillet() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let edge = cube.edges.keys().next().unwrap();
    let pre: std::collections::HashSet<_> = cube.faces.keys().collect();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];
    let result = fillet_edges(&cube, &specs, &FilletOptions::default(), tol(), None).unwrap();
    let blend = find_blend_face(&result.brep, &pre);
    let kink = blend_highlight_kink_max(
        &result.brep,
        blend,
        Vec3::new(0.0, 0.0, 1.0),
        Vec3::new(0.0, 0.0, 1.0),
        15,
        0.02,
    );
    // Either the analysis matched some pairs (Some with a finite
    // angle) or it didn't (None). Both are acceptable — we're
    // verifying the one-call wrapper plumbs correctly, not the
    // isophote algorithm's coverage.
    if let Some(angle) = kink {
        assert!(angle.is_finite(), "kink angle should be finite");
        assert!(angle >= 0.0, "kink angle should be non-negative");
    }
}

/// `blend_curvature_comb_v` returns `None` for a face with no
/// surface geometry.
#[test]
fn comb_v_returns_none_for_surfaceless_face() {
    let mut cube = BRep::unit_cube(tol()).unwrap();
    let face = cube.faces.keys().next().unwrap();
    cube.faces[face].surface = None;
    assert!(blend_curvature_comb_v(&cube, face, 10).is_none());
}
