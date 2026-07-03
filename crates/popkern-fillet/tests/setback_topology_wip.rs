//! PK-094 stage-2 WIP integration probe: drive the gated setback
//! boundary surgery through the production pipeline on the
//! pentagonal-pyramid fixture and report how far it gets.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{ContinuityLevel, FilletOptions, FilletRadius, FilletSpec, fillet_edges};

#[test]
#[ignore = "PK-094 stage-2 WIP: enable while developing the setback topology arm"]
fn setback_surgery_reaches_ngon_checkpoint_on_pentagon() {
    popkern_fillet::construct_test_hooks::enable_setback_topology_wip();
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
    let result = fillet_edges(&pyramid, &specs, &options, tol, None);
    assert!(result.is_ok(), "pipeline must stay green: {result:?}");
    let checkpoint = popkern_fillet::construct_test_hooks::take_setback_topology_checkpoint()
        .expect("WIP arm should have recorded a checkpoint");
    eprintln!("checkpoint: {checkpoint}");
    assert!(
        checkpoint.contains("setback surgery WIP checkpoint"),
        "surgery did not reach the n-gon checkpoint: {checkpoint}"
    );
    assert!(
        checkpoint.contains("arcs=5 corners=5"),
        "unexpected topology shape: {checkpoint}"
    );
    assert!(
        checkpoint.contains("blend_loop_edges=Some(15)"),
        "blend loop must be the pure setback n-gon (3 segments x 5 arcs): {checkpoint}"
    );
    assert!(
        checkpoint.contains("validate=ok"),
        "surgery result failed validation: {checkpoint}"
    );
}
