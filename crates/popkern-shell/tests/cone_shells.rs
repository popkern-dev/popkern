//! Shells of the cone primitive (PK-881).

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_shell::*;

fn tol() -> Tolerance {
    Tolerance::default()
}

#[test]
fn outward_cone_shell_serves_clean_with_honest_pcurves() {
    // PK-881: the cone primitive wired its base rim backwards for both
    // faces; the outward shell then served a misrun pcurve and a
    // canonical-backward rim (OpenCASCADE: 0.766 vs popkern 0.783; now
    // equal to 1e-10).
    use popkern_core::geometry::surface::SurfaceGeometry;
    let brep = BRep::cone(1.0, 2.0, tol()).unwrap();
    let base_cap = brep
        .faces
        .iter()
        .find(|(_, f)| {
            matches!(
                f.surface.as_ref().map(|s| &s.geometry),
                Some(SurfaceGeometry::Plane { .. })
            )
        })
        .map(|(h, _)| h)
        .unwrap();
    let result = shell(
        &brep,
        0.1,
        ShellDirection::Outward,
        &[base_cap],
        &[],
        &ShellOptions::default(),
        tol(),
        None,
    )
    .unwrap();
    let errors = result.brep.validate_all();
    assert!(errors.is_empty(), "{errors:?}");
}
