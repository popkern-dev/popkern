//! Boolean showcase: cube with cylindrical holes and sphere cutout.
//!
//! Demonstrates union, subtraction, and intersection of primitive solids.

mod helpers;

use popkern_bool::{BoolOperation, boolean, boolean_ref};
use popkern_core::geometry::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn main() {
    println!("=== Boolean Showcase ===");
    let tol = Tolerance::default();

    // 1. Build a unit cube at origin [0,1]^3.
    let cube = BRep::unit_cube(tol).expect("cube");
    println!(
        "  Cube: {} faces, {} edges",
        cube.faces.len(),
        cube.edges.len()
    );

    // 2. Create a cylinder and position it to punch a hole through the top face.
    //    Cylinder at origin, radius 0.2, height 2.0 (oversized to fully penetrate).
    let hole1 = BRep::cylinder(0.2, 2.0, tol).expect("cylinder");
    //    Translate so the cylinder is centered at (0.5, 0.5, -0.5) -- punching through Z.
    let hole1 = hole1
        .transform(&Transform3::translation(0.5, 0.5, -0.5))
        .expect("translate hole1");
    println!("  Hole cylinder: {} faces", hole1.faces.len());

    // 3. Subtract the cylinder from the cube.
    let with_hole = match boolean(BoolOperation::Subtraction, cube, hole1, tol) {
        Ok(brep) => {
            println!(
                "  After hole subtraction: {} faces, {} edges",
                brep.faces.len(),
                brep.edges.len()
            );
            brep
        }
        Err(e) => {
            eprintln!("  Hole subtraction failed: {e}");
            eprintln!("  Falling back to plain cube.");
            BRep::unit_cube(tol).expect("cube")
        }
    };

    // 4. Create a sphere cutout -- centered at corner (1, 1, 1), radius 0.4.
    let sphere = BRep::sphere(0.4, tol).expect("sphere");
    let sphere = sphere
        .transform(&Transform3::translation(1.0, 1.0, 1.0))
        .expect("translate sphere");
    println!("  Sphere: {} faces", sphere.faces.len());

    // 5. Subtract the sphere (use boolean_ref to keep with_hole for fallback).
    let final_model = match boolean_ref(BoolOperation::Subtraction, &with_hole, &sphere, tol) {
        Ok(brep) => {
            println!(
                "  After sphere cutout: {} faces, {} edges",
                brep.faces.len(),
                brep.edges.len()
            );
            brep
        }
        Err(e) => {
            eprintln!("  Sphere subtraction failed: {e}");
            eprintln!("  Using cube-with-hole result instead.");
            with_hole
        }
    };

    // Export STEP and render PNG.
    helpers::export(&final_model, "Boolean Showcase", "boolean_showcase.step");
    helpers::try_render(&final_model, "boolean_showcase.png");

    println!("=== Done ===");
}
