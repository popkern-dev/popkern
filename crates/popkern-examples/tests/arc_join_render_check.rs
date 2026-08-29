//! PK-325 visual verification renders (geometry-affecting change
//! rule): arc-join results to PNG for eyeballing. Doubles as a
//! render-path smoke test; artifacts land in the target tmpdir (path printed).

use glam::Vec3 as GVec3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_examples::offset_corpus::corpus;
use popkern_shell::{ShellOptions, offset::offset_solid, thicken::thicken};
use popkern_viz::camera::Camera;

fn tol() -> Tolerance {
    Tolerance::default()
}

fn render(brep: &BRep, name: &str, eye: (f32, f32, f32), target: (f32, f32, f32)) {
    let out_dir = std::path::Path::new(env!("CARGO_TARGET_TMPDIR"));
    std::fs::create_dir_all(out_dir).expect("output dir");
    let camera = Camera {
        eye: GVec3::new(eye.0, eye.1, eye.2),
        target: GVec3::new(target.0, target.1, target.2),
        up: GVec3::new(0.0, 0.0, 1.0),
        fov_y: 45.0_f32.to_radians(),
        aspect: 800.0 / 600.0,
        near: 0.01,
        far: 100.0,
    };
    let path = out_dir.join(format!("{name}.png"));
    popkern_viz::headless::render_to_png(brep, &camera, 800, 600, &path, true).expect("render");
    println!("[render] {}", path.display());
}

#[test]
fn render_pk325_arc_join_bodies() {
    // Rounded cube: 12 cylinder strips + 8 sphere octants.
    let cube = BRep::unit_cube(tol()).unwrap();
    let rounded = offset_solid(&cube, 0.15, &ShellOptions::default(), tol()).unwrap();
    render(
        &rounded.brep,
        "rounded_cube",
        (2.6, 2.0, 1.8),
        (0.5, 0.5, 0.5),
    );

    // Cylinder thicken: torus quarter-round rims.
    let cyl = BRep::cylinder(1.0, 2.0, tol()).unwrap();
    let tc = thicken(&cyl, 0.3, &ShellOptions::default(), tol()).unwrap();
    render(
        &tc.brep,
        "cylinder_arc_thicken",
        (3.4, 2.6, 2.6),
        (0.0, 0.0, 1.0),
    );

    // Cone offset: torus base rim + (recorded) sharp apex.
    let cone = BRep::cone(1.0, 2.0, tol()).unwrap();
    let oc = offset_solid(&cone, 0.15, &ShellOptions::default(), tol()).unwrap();
    render(
        &oc.brep,
        "cone_arc_offset",
        (3.2, 2.4, 2.2),
        (0.0, 0.0, 0.8),
    );

    // The acceptance body: wavy plate thicken (cylinder strips along
    // its straight boundary edges, sphere octant corners).
    let fixtures = corpus(tol()).unwrap();
    let wavy = fixtures.iter().find(|f| f.name == "wavy_plate").unwrap();
    let tw = thicken(&wavy.brep, 0.06, &ShellOptions::default(), tol()).unwrap();
    render(
        &tw.brep,
        "wavy_plate_arc_thicken",
        (1.9, 1.5, 1.3),
        (0.5, 0.5, 0.25),
    );
}
