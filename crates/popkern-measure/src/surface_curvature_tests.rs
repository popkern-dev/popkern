//! Analytic-fixture tests for the curvature field: sphere, cylinder,
//! and plane have closed-form K and H everywhere.

use super::*;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_mesh::TessellationOptions;

fn field_stats(brep: &BRep, kind: CurvatureKind) -> CurvatureField {
    let mesh = popkern_mesh::tessellate(brep, &TessellationOptions::default()).expect("mesh");
    curvature_field(brep, &mesh, kind).expect("field")
}

#[test]
fn sphere_curvatures_are_uniform() {
    let tol = Tolerance::default();
    let r = 0.75;
    let sphere = BRep::sphere(r, tol).expect("sphere");
    let k = field_stats(&sphere, CurvatureKind::Gaussian);
    let expected_k = 1.0 / (r * r);
    // Evaluated vertices must read K = 1/r² within analytic slack;
    // polar-degenerate chart rows may decline.
    assert!(
        (k.min - expected_k).abs() < 1e-3 * expected_k
            && (k.max - expected_k).abs() < 1e-3 * expected_k,
        "sphere K range [{:.6}, {:.6}] vs {:.6}",
        k.min,
        k.max,
        expected_k
    );
    let h = field_stats(&sphere, CurvatureKind::Mean);
    let expected_h = 1.0 / r;
    assert!(
        (h.min.abs() - expected_h).abs() < 1e-3 * expected_h
            && (h.max.abs() - expected_h).abs() < 1e-3 * expected_h,
        "sphere |H| range [{:.6}, {:.6}] vs {:.6}",
        h.min,
        h.max,
        expected_h
    );
}

#[test]
fn cube_faces_are_flat() {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).expect("cube");
    for kind in [
        CurvatureKind::Gaussian,
        CurvatureKind::Mean,
        CurvatureKind::MaxPrincipalAbs,
    ] {
        let f = field_stats(&cube, kind);
        assert!(
            f.min.abs() < 1e-9 && f.max.abs() < 1e-9,
            "flat faces must read zero {kind:?}: [{:.3e}, {:.3e}]",
            f.min,
            f.max
        );
        assert_eq!(f.unevaluated, 0, "planar charts never decline");
    }
}

#[test]
fn cylinder_principal_curvature_is_reciprocal_radius() {
    let tol = Tolerance::default();
    let r = 0.4;
    let cyl = BRep::cylinder(r, 1.2, tol).expect("cylinder");
    let f = field_stats(&cyl, CurvatureKind::MaxPrincipalAbs);
    let expected = 1.0 / r;
    // The caps are flat (0); the wall reads 1/r.
    assert!(
        (f.max - expected).abs() < 1e-3 * expected,
        "cylinder wall max |κ| {:.6} vs {:.6}",
        f.max,
        expected
    );
    let k = field_stats(&cyl, CurvatureKind::Gaussian);
    assert!(
        k.max.abs() < 1e-6 && k.min.abs() < 1e-6,
        "cylinder K must be ~0 everywhere: [{:.3e}, {:.3e}]",
        k.min,
        k.max
    );
}
