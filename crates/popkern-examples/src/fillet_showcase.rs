//! Fillet showcase: unit cube with edges filleted.
//!
//! Demonstrates the `fillet_edges` operation from popkern-fillet.

mod helpers;

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};

fn main() {
    println!("=== Fillet Showcase ===");
    let tol = Tolerance::default();

    let cube = BRep::unit_cube(tol).expect("cube");
    println!(
        "  Cube: {} faces, {} edges",
        cube.faces.len(),
        cube.edges.len()
    );

    // Collect all 12 edges of the cube.
    let edge_handles: Vec<_> = cube.edges.keys().collect();
    println!("  Filleting {} edges with radius 0.1", edge_handles.len());

    // Build fillet specs for all edges with a small radius.
    let specs: Vec<FilletSpec> = edge_handles
        .iter()
        .map(|&edge| FilletSpec {
            edge,
            radius: FilletRadius::Constant(0.05),
        })
        .collect();

    let options = FilletOptions::default();

    match fillet_edges(&cube, &specs, &options, tol, None) {
        Ok(result) => {
            let brep = &result.brep;
            println!(
                "  Filleted: {} faces, {} edges, {} vertices",
                brep.faces.len(),
                brep.edges.len(),
                brep.vertices.len()
            );
            if !result.edge_warnings.is_empty() {
                println!("  Warnings: {}", result.edge_warnings.len());
                for (edge, warning) in &result.edge_warnings {
                    println!("    edge {:?}: {:?}", edge, warning);
                }
            }
            helpers::export(brep, "Fillet Showcase", "fillet_showcase.step");
            helpers::try_render(brep, "fillet_showcase.png");
            helpers::render_multi_angle(brep, "fillet_showcase");
        }
        Err(e) => {
            eprintln!("  Full fillet failed: {e}");
            eprintln!("  Trying single-edge fillet as fallback...");

            // Try filleting just one edge as a minimal demo.
            let single_spec = vec![FilletSpec {
                edge: edge_handles[0],
                radius: FilletRadius::Constant(0.1),
            }];
            match fillet_edges(&cube, &single_spec, &options, tol, None) {
                Ok(result) => {
                    let brep = &result.brep;
                    println!(
                        "  Single fillet: {} faces, {} edges",
                        brep.faces.len(),
                        brep.edges.len()
                    );
                    helpers::export(brep, "Fillet Showcase", "fillet_showcase.step");
                    helpers::try_render(brep, "fillet_showcase.png");
                    helpers::render_multi_angle(brep, "fillet_showcase");
                }
                Err(e2) => {
                    eprintln!("  Single fillet also failed: {e2}");
                    eprintln!("  Exporting plain cube instead.");
                    helpers::export(&cube, "Fillet Showcase", "fillet_showcase.step");
                    helpers::try_render(&cube, "fillet_showcase.png");
                    helpers::render_multi_angle(&cube, "fillet_showcase");
                }
            }
        }
    }

    println!("=== Done ===");
}
