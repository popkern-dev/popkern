//! PK-326 visual verification renders (geometry-affecting change
//! rule): mixed-corner arc-join results to PNG for eyeballing.
//! Artifacts land in examples/output.

use glam::Vec3 as GVec3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_examples::offset_corpus::{l_bracket, open_sheet, thin_rib};
use popkern_shell::{ShellOptions, offset::offset_solid, thicken::thicken};
use popkern_viz::camera::Camera;

fn tol() -> Tolerance {
    Tolerance::default()
}

fn render(brep: &BRep, name: &str, eye: (f32, f32, f32), target: (f32, f32, f32)) {
    let out_dir = std::path::Path::new("../../examples/output");
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
fn render_pk326_mixed_corner_bodies() {
    // L-bracket outward: saddle corners at the reflex-edge endpoints —
    // the two top-edge strips trim each other along an ellipse.
    let lb = l_bracket(tol()).unwrap().brep;
    let lo = offset_solid(&lb, 0.1, &ShellOptions::default(), tol()).unwrap();
    render(
        &lo.brep,
        "pk326_l_bracket_arc_offset",
        (3.4, 3.0, 2.4),
        (0.7, 0.6, 0.3),
    );

    // L-bracket inward: the reflex edge erodes to a rounded concave
    // strip dead-ending into the caps along circle connectors.
    let li = offset_solid(&lb, -0.1, &ShellOptions::default(), tol()).unwrap();
    render(
        &li.brep,
        "pk326_l_bracket_arc_erosion",
        (3.4, 3.0, 2.4),
        (0.7, 0.6, 0.3),
    );

    // Thin rib outward: rib-base T-corners — vertical quarter-cylinder
    // strips dead-end into the raised plate-top plane along quarter
    // circles.
    let tr = thin_rib(tol()).unwrap().brep;
    let to = offset_solid(&tr, 0.1, &ShellOptions::default(), tol()).unwrap();
    render(
        &to.brep,
        "pk326_thin_rib_arc_offset",
        (4.2, -2.2, 2.6),
        (1.0, 0.5, 0.3),
    );

    // Dilated cone: sphere apex cap tangent to the offset cone along
    // an exact circle, torus base rim.
    let cone = BRep::cone(1.0, 1.5, tol()).unwrap();
    let oc = offset_solid(&cone, 0.15, &ShellOptions::default(), tol()).unwrap();
    render(
        &oc.brep,
        "pk326_cone_apex_cap",
        (2.8, 2.2, 2.6),
        (0.0, 0.0, 0.8),
    );

    // Open-sheet thicken: rounded outer skin with rim arc corners,
    // flat rim walls.
    let sheet = open_sheet(tol()).unwrap().brep;
    let ts = thicken(&sheet, 0.08, &ShellOptions::default(), tol()).unwrap();
    render(
        &ts.brep,
        "pk326_open_sheet_thicken",
        (2.4, 1.9, 2.2),
        (0.5, 0.5, 0.5),
    );
}
