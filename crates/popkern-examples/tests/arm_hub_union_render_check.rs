//! PK-439 visual verification render (geometry-affecting change
//! rule): the flush-cap carrier-arm ∪ round-hub union — the spider
//! carrier class the planetary gearbox routes around — to PNG for
//! eyeballing. Doubles as a render-path smoke test; the artifact
//! lands in the target tmpdir (path printed).

use glam::Vec3 as GVec3;
use popkern_part::prelude::*;
use popkern_viz::camera::Camera;

#[test]
fn render_pk439_flush_arm_hub_union() {
    let hub = Profile::on(Plane::xy_at(0.0))
        .circle(0.24)
        .extrude(0.44)
        .expect("hub");
    let t0: f64 = (0.24f64 * 0.24 - 0.09 * 0.09).sqrt();
    let arm = Profile::on(Plane::xy_at(0.0))
        .polygon(&[
            (t0 - 0.05, -0.09),
            (0.9, -0.09),
            (0.9, 0.09),
            (t0 - 0.05, 0.09),
        ])
        .extrude(0.44)
        .expect("arm");
    let out = hub.union(&arm).expect("flush arm ∪ hub union serves");

    let out_dir = std::path::Path::new(env!("CARGO_TARGET_TMPDIR"));
    std::fs::create_dir_all(out_dir).expect("output dir");
    let camera = Camera {
        eye: GVec3::new(1.4, -1.1, 1.0),
        target: GVec3::new(0.3, 0.0, 0.22),
        up: GVec3::new(0.0, 0.0, 1.0),
        fov_y: 45.0_f32.to_radians(),
        aspect: 800.0 / 600.0,
        near: 0.01,
        far: 100.0,
    };
    let path = out_dir.join("pk439_flush_arm_hub_union.png");
    popkern_viz::headless::render_to_png(out.brep(), &camera, 800, 600, &path, true)
        .expect("render");
    println!("[render] {}", path.display());
}
