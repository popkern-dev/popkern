//! PK-475 visual verification render (geometry-affecting change
//! rule): M8x1 external thread cut onto an existing r = 4 shank —
//! the flush-crest-band class whose seam-ladder loops gaped one
//! pitch while the groove tool was swept inside-out. Renders the
//! served body to PNG for eyeballing; the artifact lands in the
//! target tmpdir (path printed).

use glam::Vec3 as GVec3;
use popkern_core::geometry::{Point3, Vec3};
use popkern_part::fastener::ThreadSpec;
use popkern_part::prelude::*;
use popkern_viz::camera::Camera;

#[test]
fn render_pk475_external_thread_on_shank() {
    let spec = ThreadSpec::parse("M8x1").expect("M8 fine");
    let shank = Profile::on(Plane::XY)
        .circle(4.0)
        .extrude(10.0)
        .expect("shank");
    let threaded = shank
        .external_thread(
            &spec,
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
            6.0,
            10.0,
        )
        .expect("external thread serves");
    assert_eq!(
        threaded.open_edge_count().expect("mesh"),
        0,
        "watertight before rendering"
    );

    let out_dir = std::path::Path::new(env!("CARGO_TARGET_TMPDIR"));
    std::fs::create_dir_all(out_dir).expect("output dir");
    let camera = Camera {
        eye: GVec3::new(16.0, -14.0, 14.0),
        target: GVec3::new(0.0, 0.0, 6.0),
        up: GVec3::new(0.0, 0.0, 1.0),
        fov_y: 45.0_f32.to_radians(),
        aspect: 800.0 / 600.0,
        near: 0.1,
        far: 200.0,
    };
    let path = out_dir.join("pk475_external_thread_shank.png");
    popkern_viz::headless::render_to_png(threaded.brep(), &camera, 800, 600, &path, true)
        .expect("render");
    println!("[render] {}", path.display());
}
