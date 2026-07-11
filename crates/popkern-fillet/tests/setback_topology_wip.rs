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
    let pyramid = BRep::n_gonal_pyramid(5, 1.0, 1.2, tol).expect("pyramid");
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
    if let Some(reason) = popkern_fillet::construct_test_hooks::take_setback_topology_checkpoint() {
        eprintln!("decline: {reason}");
    }
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
    // Base-corner runoff edges must be genuinely curved (the true
    // strip∩base intersection, not the anchor chord — PK-097/098).
    {
        use popkern_core::geometry::curve::CurveEval;
        let mut checked = 0;
        for (_, edge) in brep.edges.iter() {
            let Some(curve) = edge.curve.as_ref() else {
                continue;
            };
            let ceval: &dyn CurveEval = &curve.geometry;
            let a = ceval.point_at(curve.t_start);
            let b = ceval.point_at(curve.t_end);
            let span = a.distance_to(&b);
            if a.z.abs() > 1e-6 || b.z.abs() > 1e-6 || !(0.05..=0.5).contains(&span) {
                continue;
            }
            let mid = ceval.point_at(0.5 * (curve.t_start + curve.t_end));
            let chord_mid = popkern_core::geometry::Point3::new(
                0.5 * (a.x + b.x),
                0.5 * (a.y + b.y),
                0.5 * (a.z + b.z),
            );
            assert!(
                mid.distance_to(&chord_mid) > 0.05 * span,
                "base runoff edge must be curved (sagitta {:.4} over span {:.4})",
                mid.distance_to(&chord_mid),
                span
            );
            checked += 1;
        }
        assert_eq!(checked, 5, "one curved runoff edge per base corner");
    }
    eprintln!(
        "commit-mode blend: faces={} edges={} vertices={}",
        brep.faces.len(),
        brep.edges.len(),
        brep.vertices.len()
    );
}
