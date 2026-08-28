//! PK-340: opening-rim gap corners serve on the classic shell path.
//! An outward shell with a removed face has gap edges reaching the
//! opening rim; the boundary two-sector peel closes each rim corner
//! with a connector authored as the strip × removed-ORIGINAL-surface
//! transversal (the flush contract), and the opening walls grow
//! planar corner sectors between the quads.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_shell::{ShellDirection, ShellOptions, shell};

fn tol() -> Tolerance {
    Tolerance::default()
}

#[test]
fn outward_shell_with_opening_serves_rounded_rim_corners() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let face = cube.faces.keys().next().unwrap();
    let result = shell(
        &cube,
        0.1,
        ShellDirection::Outward,
        &[face],
        &[],
        &ShellOptions::default(),
        tol(),
        None,
    )
    .expect("outward shell with opening serves");

    assert!(
        !result
            .warnings
            .iter()
            .any(|w| matches!(w, popkern_shell::ShellWarning::JoinStyleDowngraded { .. })),
        "opening-rim corners must serve without downgrade: {:?}",
        result.warnings
    );
    let errors = result.brep.validate_all();
    assert!(errors.is_empty(), "served invalid: {errors:?}");

    // Closed: no loopless half-edges survive the wall build.
    let free = result
        .brep
        .half_edges
        .iter()
        .filter(|(_, he)| !result.brep.loops.contains_key(he.loop_))
        .count();
    assert_eq!(free, 0, "free half-edges remain after opening walls");

    // The moved (void) skin is rounded: 8 gap strips (4 verticals +
    // 4 far-face edges), 4 all-gap sphere corners at the far face,
    // 4 rim arc corners. Faces: outer body 5 (opening removed) +
    // moved skin 17 + wall quads 4 + corner sectors 4.
    assert_eq!(result.brep.faces.len(), 30, "faces");
    let spheres = result
        .brep
        .faces
        .values()
        .filter(|f| {
            matches!(
                f.surface.as_ref().map(|s| &s.geometry),
                Some(popkern_core::geometry::surface::SurfaceGeometry::Sphere { .. })
            )
        })
        .count();
    assert_eq!(spheres, 4, "four sphere corner windows on the moved skin");
}

#[test]
fn midplane_shell_with_opening_serves_rounded_rim_corners() {
    // The midplane outer skin moves outward with the removed face
    // RETAINED (rim discovery): its wedge acts as the boundary cut and
    // the connector lands in the removed face's original surface.
    let cube = BRep::unit_cube(tol()).unwrap();
    let face = cube.faces.keys().next().unwrap();
    let result = shell(
        &cube,
        0.2,
        ShellDirection::Midplane,
        &[face],
        &[],
        &ShellOptions::default(),
        tol(),
        None,
    )
    .expect("midplane shell with opening serves");
    assert!(
        !result
            .warnings
            .iter()
            .any(|w| matches!(w, popkern_shell::ShellWarning::JoinStyleDowngraded { .. })),
        "midplane opening-rim corners must serve without downgrade: {:?}",
        result.warnings
    );
    assert!(result.brep.validate_all().is_empty());
    assert_eq!(result.brep.faces.len(), 30, "rounded outer + inner + walls");
}
