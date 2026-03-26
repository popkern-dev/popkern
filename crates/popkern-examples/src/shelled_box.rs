//! Shelled box: a unit cube with one face removed, hollowed out.
//!
//! Demonstrates the `shell` operation from popkern-shell.

mod helpers;

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_shell::{ShellDirection, ShellOptions, shell};

fn main() {
    println!("=== Shelled Box ===");
    let tol = Tolerance::default();

    let cube = BRep::unit_cube(tol).expect("cube");
    println!(
        "  Cube: {} faces, {} edges",
        cube.faces.len(),
        cube.edges.len()
    );

    // Remove the first face (which face we get depends on iteration order,
    // but any face works for the demo).
    let face_to_remove = cube.faces.keys().next().expect("cube has faces");
    println!("  Removing one face and shelling with thickness 0.1");

    match shell(
        &cube,
        0.1,
        ShellDirection::Inward,
        &[face_to_remove],
        &[],
        &ShellOptions::default(),
        tol,
    ) {
        Ok(result) => {
            let brep = &result.brep;
            println!(
                "  Shelled: {} faces, {} edges, {} vertices",
                brep.faces.len(),
                brep.edges.len(),
                brep.vertices.len()
            );
            if !result.warnings.is_empty() {
                println!("  Warnings: {:?}", result.warnings);
            }
            helpers::export(brep, "Shelled Box", "shelled_box.step");
            helpers::try_render(brep, "shelled_box.png");
            helpers::render_multi_angle(brep, "shelled_box");
        }
        Err(e) => {
            eprintln!("  Shell operation failed: {e}");
            eprintln!("  Exporting plain cube instead.");
            helpers::export(&cube, "Shelled Box", "shelled_box.step");
            helpers::try_render(&cube, "shelled_box.png");
            helpers::render_multi_angle(&cube, "shelled_box");
        }
    }

    println!("=== Done ===");
}
