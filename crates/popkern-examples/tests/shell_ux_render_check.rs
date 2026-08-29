//! PK-328 visual verification renders (geometry-affecting change
//! rule): the two serve classes whose construction changed with the
//! per-face material-sense resolver — the sphere-capped rod's inward
//! hollow (previously built its wall skin OUTWARD) and the bracket
//! hollow (previously declined on phantom trim regions). Both render
//! as HALF bodies cut by dropping faces on one side of the camera so
//! the void skin is visible. Artifacts land in the target tmpdir (path printed).

use glam::Vec3 as GVec3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_examples::offset_corpus::{bracket_bore_countersink, sphere_capped_rod};
use popkern_shell::{ShellDirection, ShellOptions, shell};
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
fn render_pk328_sense_fixed_hollows() {
    let rod = sphere_capped_rod(tol()).unwrap().brep;
    let hollow = shell(
        &rod,
        0.1,
        ShellDirection::Inward,
        &[],
        &[],
        &ShellOptions::default(),
        tol(),
        None,
    )
    .expect("rod hollow serves");
    assert!(hollow.brep.validate_all().is_empty());
    render(&hollow.brep, "rod_hollow", (1.4, 1.1, 1.8), (0.0, 0.0, 0.7));

    let bracket = bracket_bore_countersink(tol()).unwrap().brep;
    let hollow = shell(
        &bracket,
        0.1,
        ShellDirection::Inward,
        &[],
        &[],
        &ShellOptions::default(),
        tol(),
        None,
    )
    .expect("bracket hollow serves");
    assert!(hollow.brep.validate_all().is_empty());
    render(
        &hollow.brep,
        "bracket_hollow",
        (3.2, 2.4, 2.0),
        (1.0, 0.5, 0.25),
    );

    // OPEN-SIDE shell (remove the hole-free x=0 wall) so the hollow
    // interior and the grown bore skin are directly visible.
    let side = bracket
        .faces
        .iter()
        .find(|(fh, f)| {
            f.surface.as_ref().is_some_and(|s| {
                matches!(
                    &s.geometry,
                    popkern_core::geometry::surface::SurfaceGeometry::Plane { normal, .. }
                        if normal.x.abs() > 0.9
                ) && f.inner_loops.is_empty()
                    && bracket
                        .face_outer_loop_positions(*fh)
                        .iter()
                        .all(|p| p.x < 0.5)
            })
        })
        .map(|(fh, _)| fh)
        .expect("bracket x=0 side face");
    let cup = shell(
        &bracket,
        0.06,
        ShellDirection::Inward,
        &[side],
        &[],
        &ShellOptions::default(),
        tol(),
        None,
    )
    .expect("bracket open cup serves");
    render(
        &cup.brep,
        "bracket_open_cup",
        (-1.8, 1.6, 1.6),
        (0.7, 0.5, 0.25),
    );
}
