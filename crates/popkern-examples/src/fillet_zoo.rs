//! Fillet zoo: showcase fillets at every PnS-supported valence
//! (N ∈ {3, 5, 6, 7, 8, 10, 12}) on n-gonal pyramids, plus a
//! filleted hexagonal prism for a non-pyramid sample.
//!
//! Each model gets STEP + multi-angle PNG output under
//! `examples/output/fillet_zoo/`. Useful as a demo, and as a
//! regression catch — anything that crashes here lands a real bug.

mod helpers;

use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, VertexHandle};
use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};

fn main() {
    println!("=== Fillet Zoo ===");
    let tol = Tolerance::default();

    let valences: &[usize] = &[3, 4, 5, 6, 7, 8, 10, 12];
    for &n in valences {
        showcase_pyramid(n, tol);
    }

    showcase_hex_prism(tol);
    showcase_oct_prism(tol);

    println!("=== Done ===");
}

fn showcase_pyramid(n: usize, tol: Tolerance) {
    let name = format!("pyramid_n{n:02}");
    println!("\n--- {} ---", name);

    let pyramid = match BRep::n_gonal_pyramid(n, 1.0, 1.2, tol) {
        Ok(p) => p,
        Err(e) => {
            eprintln!("  build failed: {e}");
            return;
        }
    };
    println!(
        "  BRep: {} faces, {} edges, {} vertices",
        pyramid.faces.len(),
        pyramid.edges.len(),
        pyramid.vertices.len()
    );

    let apex = pyramid
        .vertices
        .iter()
        .find(|(_, v)| v.position.x.abs() < 1e-9 && v.position.y.abs() < 1e-9 && v.position.z > 0.0)
        .map(|(h, _)| h)
        .expect("apex vertex");

    let apex_edges = collect_incident_edges(&pyramid, apex);
    println!("  Apex valence: {}", apex_edges.len());

    // Mixed radii ensure the rolling-ball sphere fit doesn't snap
    // and the dispatch lands on the Gregory / PnS rung. Walks
    // through 0.07–0.13 across the apex spokes so each pyramid is
    // a real test of the high-valence cap path.
    let specs: Vec<FilletSpec> = apex_edges
        .iter()
        .enumerate()
        .map(|(i, &edge)| FilletSpec {
            edge,
            radius: FilletRadius::Constant(
                0.07 + 0.06 * (i as f64) / ((apex_edges.len() - 1).max(1) as f64),
            ),
        })
        .collect();
    let options = FilletOptions::default();

    match fillet_edges(&pyramid, &specs, &options, tol, None) {
        Ok(result) => {
            let brep = &result.brep;
            println!(
                "  Filleted: {} faces, {} edges, {} vertices, {} corner patches",
                brep.faces.len(),
                brep.edges.len(),
                brep.vertices.len(),
                result.corner_patches.len()
            );
            let methods: std::collections::BTreeMap<String, usize> = result
                .corner_patches
                .iter()
                .fold(std::collections::BTreeMap::new(), |mut m, p| {
                    *m.entry(format!("{:?}", p.method)).or_insert(0) += 1;
                    m
                });
            for (method, count) in &methods {
                println!("    method {method}: {count} patches");
            }
            if !result.edge_warnings.is_empty() {
                let kinds: std::collections::BTreeMap<_, usize> = result.edge_warnings.iter().fold(
                    std::collections::BTreeMap::new(),
                    |mut m, (_, w)| {
                        *m.entry(format!("{:?}", std::mem::discriminant(w)))
                            .or_insert(0) += 1;
                        m
                    },
                );
                println!("  Warnings: {} ({:?})", result.edge_warnings.len(), kinds);
            }
            let step_path = format!("fillet_zoo/{name}.step");
            let png_path = format!("fillet_zoo/{name}.png");
            helpers::export(brep, &name, &step_path);
            helpers::try_render(brep, &png_path);
            helpers::render_multi_angle(brep, &format!("fillet_zoo/{name}"));
        }
        Err(e) => eprintln!("  fillet_edges failed: {e}"),
    }
}

fn showcase_hex_prism(tol: Tolerance) {
    let name = "hex_prism";
    println!("\n--- {name} ---");

    let r = 1.0;
    let h = 0.6;
    let profile: Vec<Point3> = (0..6)
        .map(|k| {
            let theta = 2.0 * std::f64::consts::PI * (k as f64) / 6.0;
            Point3::new(r * theta.cos(), r * theta.sin(), 0.0)
        })
        .collect();
    let prism = match BRep::prism(&profile, Vec3::new(0.0, 0.0, 1.0), h, tol) {
        Ok(p) => p,
        Err(e) => {
            eprintln!("  build failed: {e}");
            return;
        }
    };
    println!(
        "  BRep: {} faces, {} edges, {} vertices",
        prism.faces.len(),
        prism.edges.len(),
        prism.vertices.len()
    );

    let edges: Vec<EdgeHandle> = prism.edges.keys().collect();
    let specs: Vec<FilletSpec> = edges
        .iter()
        .map(|&edge| FilletSpec {
            edge,
            radius: FilletRadius::Constant(0.04),
        })
        .collect();
    let options = FilletOptions::default();

    match fillet_edges(&prism, &specs, &options, tol, None) {
        Ok(result) => {
            let brep = &result.brep;
            println!(
                "  Filleted: {} faces, {} edges, {} corner patches",
                brep.faces.len(),
                brep.edges.len(),
                result.corner_patches.len()
            );
            if !result.edge_warnings.is_empty() {
                println!("  Warnings: {}", result.edge_warnings.len());
            }
            helpers::export(brep, name, &format!("fillet_zoo/{name}.step"));
            helpers::try_render(brep, &format!("fillet_zoo/{name}.png"));
            helpers::render_multi_angle(brep, &format!("fillet_zoo/{name}"));
        }
        Err(e) => eprintln!("  fillet_edges failed: {e}"),
    }
}

fn showcase_oct_prism(tol: Tolerance) {
    let name = "oct_prism";
    println!("\n--- {name} ---");

    // Octagonal prism — purely convex edges, every base/top
    // vertex is N=3 trihedral, mid-plane edges are 2-edge meets.
    // Showcases the rolling-ball-sphere + rational-Gregory rungs
    // on a smooth wide-angle profile.
    let r = 1.0;
    let h = 0.5;
    let profile: Vec<Point3> = (0..8)
        .map(|k| {
            let theta = 2.0 * std::f64::consts::PI * (k as f64) / 8.0;
            Point3::new(r * theta.cos(), r * theta.sin(), 0.0)
        })
        .collect();
    let prism = match BRep::prism(&profile, Vec3::new(0.0, 0.0, 1.0), h, tol) {
        Ok(p) => p,
        Err(e) => {
            eprintln!("  build failed: {e}");
            return;
        }
    };
    println!(
        "  BRep: {} faces, {} edges, {} vertices",
        prism.faces.len(),
        prism.edges.len(),
        prism.vertices.len()
    );

    let edges: Vec<EdgeHandle> = prism.edges.keys().collect();
    let specs: Vec<FilletSpec> = edges
        .iter()
        .map(|&edge| FilletSpec {
            edge,
            radius: FilletRadius::Constant(0.05),
        })
        .collect();
    let options = FilletOptions::default();

    match fillet_edges(&prism, &specs, &options, tol, None) {
        Ok(result) => {
            let brep = &result.brep;
            println!(
                "  Filleted: {} faces, {} edges, {} corner patches",
                brep.faces.len(),
                brep.edges.len(),
                result.corner_patches.len()
            );
            if !result.edge_warnings.is_empty() {
                println!("  Warnings: {}", result.edge_warnings.len());
            }
            helpers::export(brep, name, &format!("fillet_zoo/{name}.step"));
            helpers::try_render(brep, &format!("fillet_zoo/{name}.png"));
            helpers::render_multi_angle(brep, &format!("fillet_zoo/{name}"));
        }
        Err(e) => eprintln!("  fillet_edges failed: {e}"),
    }
}

fn collect_incident_edges(brep: &BRep, vertex: VertexHandle) -> Vec<EdgeHandle> {
    brep.edges
        .iter()
        .filter_map(|(eh, _)| {
            let he = brep.edges[eh].half_edge;
            let origin = brep.half_edges[he].origin;
            let twin_origin = brep.half_edges[brep.half_edges[he].twin].origin;
            (origin == vertex || twin_origin == vertex).then_some(eh)
        })
        .collect()
}
