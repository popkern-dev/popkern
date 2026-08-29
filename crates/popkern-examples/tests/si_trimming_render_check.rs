//! PK-327 visual verification renders (geometry-affecting change
//! rule): the globally trimmed grow of the grooved cylinder next to
//! its sub-reach sibling. Artifacts land in the target tmpdir (path printed).

use glam::Vec3 as GVec3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_examples::offset_corpus::{grooved_cylinder, torus_lens};
use popkern_shell::{ShellOptions, offset::offset_solid};
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
fn render_pk327_trimmed_grow_bodies() {
    let fixture = grooved_cylinder(tol()).unwrap().brep;
    render(
        &fixture,
        "grooved_cylinder_source",
        (2.2, 1.8, 1.6),
        (0.0, 0.0, 0.5),
    );
    // Sub-reach grow: the groove survives, narrowed, with mouth tori.
    let sub = offset_solid(&fixture, 0.03, &ShellOptions::default(), tol()).unwrap();
    render(
        &sub.brep,
        "grooved_cylinder_grow_sub_reach",
        (2.2, 1.8, 1.6),
        (0.0, 0.0, 0.5),
    );
    // Past-reach grow: the groove closes into a shallow V-crease where
    // the two mouth tori weld along their exact crossing circle.
    let over = offset_solid(&fixture, 0.08, &ShellOptions::default(), tol()).unwrap();
    render(
        &over.brep,
        "grooved_cylinder_grow_trimmed",
        (2.2, 1.8, 1.6),
        (0.0, 0.0, 0.5),
    );
}

#[test]
fn render_pk327_torus_lens_apple() {
    let fixture = torus_lens(tol()).unwrap().brep;
    render(
        &fixture,
        "torus_lens_source",
        (2.6, 2.2, 1.8),
        (0.0, 0.0, 0.0),
    );
    let apple = offset_solid(&fixture, 0.2, &ShellOptions::default(), tol()).unwrap();
    render(
        &apple.brep,
        "torus_lens_apple",
        (2.6, 2.2, 1.8),
        (0.0, 0.0, 0.0),
    );
    // Dimple view from above the pinch.
    render(
        &apple.brep,
        "torus_lens_apple_dimple",
        (0.9, 0.7, 2.6),
        (0.0, 0.0, 0.1),
    );
}
