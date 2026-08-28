//! PK-346 gate: a removed face pierced by a bore carries the bore rim
//! on an INNER loop; that rim is an opening boundary and gets its own
//! wall ring — the outer-loop-only walk left its half-edges loopless
//! (structural BrokenCycle) and the whole opening declined.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_shell::{ShellDirection, ShellOptions, shell};

/// A 1×1×0.5 block with a through-bore: the top face carries the bore
/// rim as an inner loop.
fn bored_block(tol: Tolerance) -> BRep {
    use popkern_core::geometry::transform::Transform3;
    let block = BRep::prism(
        &[
            popkern_core::geometry::Point3::new(0.0, 0.0, 0.0),
            popkern_core::geometry::Point3::new(1.0, 0.0, 0.0),
            popkern_core::geometry::Point3::new(1.0, 1.0, 0.0),
            popkern_core::geometry::Point3::new(0.0, 1.0, 0.0),
        ],
        popkern_core::geometry::Vec3::new(0.0, 0.0, 1.0),
        0.5,
        tol,
    )
    .expect("block");
    let bore = BRep::cylinder(0.2, 0.7, tol)
        .expect("bore")
        .transform(&Transform3::translation(0.5, 0.5, -0.1))
        .expect("place bore");
    let a = popkern_part::Part::from_brep(block, tol, "block").expect("part");
    let b = popkern_part::Part::from_brep(bore, tol, "bore").expect("part");
    a.subtract(&b).expect("subtract").brep().clone()
}

#[test]
fn opening_with_inner_rim_serves_walls_on_both_loops() {
    let tol = Tolerance::default();
    let brep = bored_block(tol);
    let top = brep
        .faces
        .iter()
        .find(|(fh, face)| {
            face.surface.as_ref().is_some_and(|s| {
                matches!(
                    &s.geometry,
                    SurfaceGeometry::Plane { normal, origin, .. }
                        if normal.z.abs() > 0.9 && origin.z > 0.4
                )
            }) && !brep.face_outer_loop_positions(*fh).is_empty()
                && !face.inner_loops.is_empty()
        })
        .map(|(fh, _)| fh)
        .expect("bored top face");
    let result = shell(
        &brep,
        0.06,
        ShellDirection::Inward,
        &[top],
        &[],
        &ShellOptions::default(),
        tol,
        None,
    )
    .expect("bore-pierced opening must serve");
    assert!(
        result.brep.validate_all().is_empty(),
        "opening with inner rim serves validate_all-clean"
    );
    // Both rims walled: the serve gains wall faces beyond the plain
    // outer+inner skins.
    assert!(
        result.brep.faces.len() > (brep.faces.len() - 1) * 2,
        "wall faces present: {}",
        result.brep.faces.len()
    );
}

/// Visual verification (geometry-affecting change rule): the
/// bore-pierced top opening of the corpus bracket, walls on both rims.
#[test]
fn render_pk346_bore_pierced_opening() {
    use glam::Vec3 as GVec3;
    use popkern_examples::offset_corpus::bracket_bore_countersink;
    use popkern_viz::camera::Camera;

    let tol = Tolerance::default();
    let f = bracket_bore_countersink(tol).expect("fixture");
    let top = f
        .brep
        .faces
        .iter()
        .find(|(fh, face)| {
            face.surface.as_ref().is_some_and(|s| {
                matches!(
                    &s.geometry,
                    SurfaceGeometry::Plane { normal, origin, .. }
                        if normal.z.abs() > 0.9 && origin.z > 0.4
                )
            }) && !f.brep.face_outer_loop_positions(*fh).is_empty()
        })
        .map(|(fh, _)| fh)
        .expect("top face");
    let cup = shell(
        &f.brep,
        0.06,
        ShellDirection::Inward,
        &[top],
        &[],
        &ShellOptions::default(),
        tol,
        None,
    )
    .expect("serves");
    let out_dir = std::path::Path::new("../../examples/output");
    std::fs::create_dir_all(out_dir).expect("output dir");
    let camera = Camera {
        eye: GVec3::new(2.6, 2.2, 2.4),
        target: GVec3::new(1.0, 0.5, 0.2),
        up: GVec3::new(0.0, 0.0, 1.0),
        fov_y: 45.0_f32.to_radians(),
        aspect: 800.0 / 600.0,
        near: 0.01,
        far: 100.0,
    };
    let path = out_dir.join("bore_pierced_opening.png");
    popkern_viz::headless::render_to_png(&cup.brep, &camera, 800, 600, &path, true)
        .expect("render");
}
