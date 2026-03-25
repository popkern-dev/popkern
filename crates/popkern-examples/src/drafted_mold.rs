//! Drafted mold: a unit cube with one face drafted for mold withdrawal.
//!
//! Demonstrates the `draft_face` operation from popkern-draft.

mod helpers;

use popkern_core::geometry::Vec3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_draft::{NeutralReference, draft_face};

fn main() {
    println!("=== Drafted Mold ===");
    let tol = Tolerance::default();

    let cube = BRep::unit_cube(tol).expect("cube");
    println!(
        "  Cube: {} faces, {} edges",
        cube.faces.len(),
        cube.edges.len()
    );

    // Draft a face with a 5-degree taper (0.087 radians).
    let face = cube.faces.keys().next().expect("cube has faces");
    let pull_dir = Vec3::new(0.0, 0.0, 1.0);
    let angle: f64 = 0.087; // ~5 degrees

    println!("  Drafting one face at {:.1} degrees", angle.to_degrees());

    match draft_face(&cube, face, angle, pull_dir, NeutralReference::Auto, tol) {
        Ok(result) => {
            let brep = &result.brep;
            println!(
                "  Drafted: {} faces, {} edges, {} vertices",
                brep.faces.len(),
                brep.edges.len(),
                brep.vertices.len()
            );
            if !result.warnings.is_empty() {
                println!("  Warnings: {:?}", result.warnings);
            }
            helpers::export(brep, "Drafted Mold", "drafted_mold.step");
            helpers::try_render(brep, "drafted_mold.png");
        }
        Err(e) => {
            eprintln!("  Draft operation failed: {e}");
            eprintln!("  Exporting plain cube instead.");
            helpers::export(&cube, "Drafted Mold", "drafted_mold.step");
            helpers::try_render(&cube, "drafted_mold.png");
        }
    }

    println!("=== Done ===");
}
