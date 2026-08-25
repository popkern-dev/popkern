//! PK-324 visual verification renders (geometry-affecting change
//! rule): serves the newly-serving bodies to PNG for eyeballing.
//! Doubles as a render-path smoke test; artifacts land in examples/output.

use glam::Vec3 as GVec3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_examples::offset_corpus::{bracket_bore_countersink, filleted_cube};
use popkern_shell::{ShellDirection, ShellOptions, shell, thicken::thicken};
use popkern_viz::camera::Camera;

fn tol() -> Tolerance {
    Tolerance::default()
}

fn render(brep: &BRep, name: &str, eye: (f32, f32, f32)) {
    let out_dir = std::path::Path::new("../../examples/output");
    std::fs::create_dir_all(out_dir).expect("output dir");
    let camera = Camera {
        eye: GVec3::new(eye.0, eye.1, eye.2),
        target: GVec3::new(0.0, 0.0, 0.0),
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
fn render_pk324_serving_bodies() {
    // Hollow shelled sphere and torus (validate_all-clean since PK-324).
    let sphere = BRep::sphere(0.5, tol()).unwrap();
    let shelled_sphere = shell(
        &sphere,
        0.06,
        ShellDirection::Inward,
        &[],
        &[],
        &ShellOptions::default(),
        tol(),
        None,
    )
    .unwrap();
    render(
        &shelled_sphere.brep,
        "pk324_shelled_sphere",
        (1.4, 1.0, 0.9),
    );

    let torus = BRep::torus(0.6, 0.25, tol()).unwrap();
    let shelled_torus = shell(
        &torus,
        0.06,
        ShellDirection::Inward,
        &[],
        &[],
        &ShellOptions::default(),
        tol(),
        None,
    )
    .unwrap();
    render(&shelled_torus.brep, "pk324_shelled_torus", (1.8, 1.4, 1.2));

    let fc = filleted_cube(tol()).unwrap().brep;
    let thick_fc = thicken(&fc, 0.05, &ShellOptions::default(), tol()).unwrap();
    render(&thick_fc.brep, "pk324_thick_filleted_cube", (2.2, 1.7, 1.5));

    let bracket = bracket_bore_countersink(tol()).unwrap().brep;
    let thick_bracket = thicken(&bracket, 0.02, &ShellOptions::default(), tol()).unwrap();
    render(&thick_bracket.brep, "pk324_thick_bracket", (3.2, 2.4, 2.0));

    let open = BRep::open_box(tol()).unwrap();
    let thick_open = thicken(&open, 0.06, &ShellOptions::default(), tol()).unwrap();
    render(&thick_open.brep, "pk324_thick_open_box", (2.0, 1.5, -1.8));
}
