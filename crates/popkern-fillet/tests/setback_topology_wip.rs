//! PK-094 stage-2 WIP integration probe: drive the gated setback
//! boundary surgery through the production pipeline on the
//! pentagonal-pyramid fixture and report how far it gets.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{ContinuityLevel, FilletOptions, FilletRadius, FilletSpec, fillet_edges};

#[test]
#[ignore = "PK-094 stage-2: commit-mode end-to-end (enable with the WIP arm)"]
fn setback_commit_produces_valid_g2_vertex_blend() {
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
    brep.validate().expect("blended BRep must validate");
    assert_eq!(
        result.corner_patches.len(),
        20,
        "expected the 4n KP20 piece faces as corner patches"
    );
    for cp in &result.corner_patches {
        assert_eq!(cp.achieved_continuity, ContinuityLevel::G2);
    }
    let gregory_faces = brep
        .faces
        .iter()
        .filter(|(_, f)| {
            matches!(
                f.surface.as_ref().map(|s| &s.geometry),
                Some(popkern_core::geometry::SurfaceGeometry::Gregory(_))
            )
        })
        .count();
    assert_eq!(gregory_faces, 20, "expected 20 Gregory cap faces");
    eprintln!(
        "commit-mode blend: faces={} edges={} vertices={}",
        brep.faces.len(),
        brep.edges.len(),
        brep.vertices.len()
    );
}
