//! PK-480 visual verification: render the angled tapped hole (30°
//! off vertical, M6, blind, entry chamfer) that the fixture
//! `angled_tapped_hole_serves_watertight` gates numerically.

use popkern_core::geometry::{Point3, Vec3};
use popkern_part::Part;
use popkern_part::fastener::{HoleDepth, TapOptions, ThreadSpec};
use popkern_part::prelude::*;

fn main() {
    let spec = ThreadSpec::parse("M6").expect("M6");
    let plate: Part = Profile::on(Plane::XY)
        .rect(20.0, 20.0)
        .extrude(8.0)
        .expect("block");
    let tapped = plate
        .tapped_hole(
            &spec,
            Point3::new(1.0, -0.5, 8.0),
            Vec3::new(0.35, 0.2, -1.0),
            HoleDepth::Blind(4.5),
            &TapOptions {
                clearance: 0.15,
                entry_chamfer: true,
            },
        )
        .expect("angled tapped hole");
    println!(
        "open edges: {:?}",
        tapped.open_edge_count().expect("mesh oracle")
    );
    let brep = tapped.brep();
    let out = std::path::Path::new("examples/output/pk480_angled_tap.png");
    std::fs::create_dir_all("examples/output").expect("output dir");
    let camera = popkern_viz::camera::Camera {
        eye: glam::Vec3::new(16.0, -22.0, 24.0),
        target: glam::Vec3::new(1.0, -0.5, 6.0),
        up: glam::Vec3::new(0.0, 0.0, 1.0),
        fov_y: 40.0_f32.to_radians(),
        aspect: 1000.0 / 800.0,
        near: 0.01,
        far: 200.0,
    };
    popkern_viz::headless::render_to_png(brep, &camera, 1000, 800, out, true).expect("render");
    // A second view straight down the hole axis for the chamfer rim.
    let camera2 = popkern_viz::camera::Camera {
        eye: glam::Vec3::new(1.0 + 0.35 * 14.0, -0.5 + 0.2 * 14.0, 8.0 + 14.0),
        target: glam::Vec3::new(1.0, -0.5, 6.0),
        up: glam::Vec3::new(0.0, 1.0, 0.0),
        fov_y: 40.0_f32.to_radians(),
        aspect: 1.0,
        near: 0.01,
        far: 200.0,
    };
    let out2 = std::path::Path::new("examples/output/pk480_angled_tap_axis.png");
    popkern_viz::headless::render_to_png(brep, &camera2, 900, 900, out2, true).expect("render");
    println!("wrote {} and {}", out.display(), out2.display());
}
