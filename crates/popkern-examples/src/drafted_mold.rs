//! Drafted mold: a unit cube with one face drafted for mold withdrawal.
//!
//! Demonstrates the `draft_face` operation from popkern-draft.

mod helpers;

use popkern_core::geometry::Vec3;
use popkern_core::geometry::surface::SurfaceEval;
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

    // Draft a side face (one whose normal is NOT parallel to the pull direction)
    // so that the neutral plane intersects the face cleanly.
    let pull_dir = Vec3::new(0.0, 0.0, 1.0);
    let face = cube
        .faces
        .iter()
        .find(|(_, fd)| {
            fd.surface.as_ref().is_some_and(|sof| {
                let n = sof.geometry.normal_at(0.5, 0.5);
                // Pick a face whose normal is roughly perpendicular to pull_dir
                // (i.e., a side face, not top or bottom).
                n.dot(&pull_dir).abs() < 0.5
            })
        })
        .map(|(fh, _)| fh)
        .expect("cube should have a side face");
    let angle: f64 = 0.087; // ~5 degrees

    println!("  Drafting one face at {:.1} degrees", angle.to_degrees());

    match draft_face(
        &cube,
        face,
        angle,
        pull_dir,
        NeutralReference::Auto,
        tol,
        None,
    ) {
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
            // Verify no vertex positions are NaN before exporting.
            let nan_count = brep
                .vertices
                .values()
                .filter(|v| {
                    !v.position.x.is_finite()
                        || !v.position.y.is_finite()
                        || !v.position.z.is_finite()
                })
                .count();
            if nan_count > 0 {
                eprintln!("  WARNING: {nan_count} vertices have non-finite positions");
            }
            helpers::export(brep, "Drafted Mold", "drafted_mold.step");
            helpers::try_render(brep, "drafted_mold.png");
            helpers::render_multi_angle(brep, "drafted_mold");
        }
        Err(e) => {
            eprintln!("  Draft operation failed: {e}");
            eprintln!("  Exporting plain cube instead.");
            helpers::export(&cube, "Drafted Mold", "drafted_mold.step");
            helpers::try_render(&cube, "drafted_mold.png");
            helpers::render_multi_angle(&cube, "drafted_mold");
        }
    }

    println!("=== Done ===");
}
