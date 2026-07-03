//! PK-094 stage-2: tessellation export of the setback vertex blend
//! (commit mode) — the honest end-to-end replacing the showcase's
//! from-parts writer.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{ContinuityLevel, FilletOptions, FilletRadius, FilletSpec, fillet_edges};

#[test]
#[ignore = "PK-094 stage-2: enable with the setback commit arm"]
fn setback_blend_tessellates_closed() {
    popkern_fillet::construct_test_hooks::enable_setback_topology_commit();
    let tol = Tolerance::default();
    let pyramid = BRep::n_gonal_pyramid(5, 1.0, 1.2, tol.clone()).expect("pyramid");
    let apex = pyramid
        .vertices
        .iter()
        .max_by(|a, b| {
            a.1.position
                .z
                .partial_cmp(&b.1.position.z)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|(h, _)| h)
        .expect("apex");
    let apex_edges: Vec<_> = pyramid
        .edges
        .iter()
        .filter_map(|(eh, _)| {
            let he = pyramid.edges[eh].half_edge;
            let origin = pyramid.half_edges[he].origin;
            let twin_origin = pyramid.half_edges[pyramid.half_edges[he].twin].origin;
            (origin == apex || twin_origin == apex).then_some(eh)
        })
        .collect();
    let specs: Vec<FilletSpec> = apex_edges
        .iter()
        .map(|&edge| FilletSpec {
            edge,
            radius: FilletRadius::Constant(0.16),
        })
        .collect();
    let options = FilletOptions {
        continuity: ContinuityLevel::G2,
        ..Default::default()
    };
    let result = fillet_edges(&pyramid, &specs, &options, tol, None).expect("pipeline");
    let brep = &result.brep;
    brep.validate().expect("blend must validate");

    let mesh =
        crate::tessellate(brep, &crate::TessellationOptions::default()).expect("tessellation");
    eprintln!(
        "tessellated: {} vertices, {} triangles",
        mesh.vertices.len(),
        mesh.triangles.len()
    );
    assert!(mesh.triangles.len() > 100, "implausibly small mesh");

    // PK-095 closure gate: the whole-object tessellation must be
    // watertight — after position welding, every edge is used by
    // exactly two triangles.
    let mut canon: std::collections::HashMap<(i64, i64, i64), usize> =
        std::collections::HashMap::new();
    let mut remap = Vec::with_capacity(mesh.vertices.len());
    for v in &mesh.vertices {
        let k = (
            (v.position.x * 1.0e7).round() as i64,
            (v.position.y * 1.0e7).round() as i64,
            (v.position.z * 1.0e7).round() as i64,
        );
        let next = canon.len();
        remap.push(*canon.entry(k).or_insert(next));
    }
    let mut counts: std::collections::HashMap<(usize, usize), usize> =
        std::collections::HashMap::new();
    for t in &mesh.triangles {
        let w = [remap[t[0]], remap[t[1]], remap[t[2]]];
        if w[0] == w[1] || w[1] == w[2] || w[0] == w[2] {
            continue;
        }
        for (a, b) in [(w[0], w[1]), (w[1], w[2]), (w[2], w[0])] {
            *counts.entry((a.min(b), a.max(b))).or_insert(0) += 1;
        }
    }
    let open: usize = counts.values().filter(|&&c| c == 1).count();
    assert_eq!(open, 0, "tessellation must be watertight (open edges)");
}
