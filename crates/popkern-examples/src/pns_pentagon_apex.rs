//! PnS N=5 showcase: pentagonal pyramid with the apex edges
//! filleted. Gated by `POPKERN_PNS=1` to force the PnS dispatch
//! rung; without it the fillet falls back through the existing
//! ladder. Produces STEP + multi-angle PNGs under
//! `examples/output/pns_pentagon_apex*`.

mod helpers;

use std::path::{Path, PathBuf};

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, BrepRepairDiagnosticExport};
use popkern_examples::repair_artifacts::write_whole_brep_repair_knowledge_artifact;
use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};

fn main() {
    println!("=== PnS Pentagonal Pyramid Apex ===");
    let tol = Tolerance::default();

    let pyramid = BRep::n_gonal_pyramid(5, 1.0, 1.2, tol).expect("pentagonal pyramid");
    println!(
        "  Pyramid: {} faces, {} edges, {} vertices",
        pyramid.faces.len(),
        pyramid.edges.len(),
        pyramid.vertices.len()
    );

    // Apex is the vertex at (0, 0, 1.2) — collect the 5 edges
    // incident on it.
    let apex = pyramid
        .vertices
        .iter()
        .find(|(_, v)| v.position.x.abs() < 1e-9 && v.position.y.abs() < 1e-9 && v.position.z > 0.0)
        .map(|(h, _)| h)
        .expect("apex vertex");

    let apex_edges: Vec<_> = pyramid
        .edges
        .iter()
        .filter_map(|(eh, _)| {
            let he = pyramid.edges[eh].half_edge;
            let origin = pyramid.half_edges[he].origin;
            let twin_origin = pyramid.half_edges[pyramid.half_edges[he].twin].origin;
            if origin == apex || twin_origin == apex {
                Some(eh)
            } else {
                None
            }
        })
        .collect();
    println!("  Apex has {} incident edges", apex_edges.len());

    let specs: Vec<FilletSpec> = apex_edges
        .iter()
        .map(|&edge| FilletSpec {
            edge,
            radius: FilletRadius::Constant(0.1),
        })
        .collect();
    let options = FilletOptions::default();

    match fillet_edges(&pyramid, &specs, &options, tol, None) {
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
            helpers::export(brep, "PnS Pentagon Apex", "pns_pentagon_apex.step");
            write_apex_repair_knowledge(brep).expect("write apex repair knowledge artifact");
            helpers::try_render(brep, "pns_pentagon_apex.png");
            helpers::render_multi_angle(brep, "pns_pentagon_apex");
        }
        Err(e) => {
            eprintln!("  Fillet failed: {e}");
            helpers::export(&pyramid, "PnS Pentagon Apex", "pns_pentagon_apex.step");
            write_apex_repair_knowledge(&pyramid).expect("write apex repair knowledge artifact");
            helpers::try_render(&pyramid, "pns_pentagon_apex.png");
            helpers::render_multi_angle(&pyramid, "pns_pentagon_apex");
        }
    }

    println!("=== Done ===");
}

fn write_apex_repair_knowledge(brep: &BRep) -> Result<(), String> {
    let path = output_path("pns_pentagon_apex_repair_knowledge.json");
    write_whole_brep_repair_knowledge_artifact(
        &path,
        "pns_pentagon_apex",
        brep,
        9,
        BrepRepairDiagnosticExport::FullJson,
    )
}

fn output_path(filename: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .join("examples")
        .join("output")
        .join(filename)
}
