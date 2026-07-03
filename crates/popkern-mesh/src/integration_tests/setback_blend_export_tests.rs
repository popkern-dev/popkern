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
}
