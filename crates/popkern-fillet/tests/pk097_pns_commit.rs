//! PK-097 Phase B gate: commit-mode G1 blend — the setback surgery
//! with the PnS fan on the pentagon.
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

#[test]
#[ignore = "PK-097 Phase B: commit-mode G1 vertex blend end-to-end"]
fn pns_commit_produces_valid_g1_vertex_blend() {
    // SAFETY: test-local env; routes the G1 apex through the PnS
    // bi-4 comparison family (the production default is KP20).
    unsafe { std::env::set_var("POPKERN_G1_PNS_CAP", "1") };
    popkern_fillet::construct_test_hooks::enable_setback_topology_commit();
    let tol = Tolerance::default();
    let pyramid = BRep::n_gonal_pyramid(5, 1.0, 1.2, tol.clone()).expect("pyramid");
    let apex = pyramid
        .vertices
        .iter()
        .max_by(|a, b| a.1.position.z.partial_cmp(&b.1.position.z).unwrap())
        .map(|(h, _)| h)
        .expect("apex");
    let apex_edges: Vec<_> = pyramid
        .edges
        .iter()
        .filter_map(|(eh, _)| {
            let he = pyramid.edges[eh].half_edge;
            let o = pyramid.half_edges[he].origin;
            let t = pyramid.half_edges[pyramid.half_edges[he].twin].origin;
            (o == apex || t == apex).then_some(eh)
        })
        .collect();
    let specs: Vec<popkern_fillet::FilletSpec> = apex_edges
        .iter()
        .map(|&edge| popkern_fillet::FilletSpec {
            edge,
            radius: popkern_fillet::FilletRadius::Constant(0.16),
        })
        .collect();
    let options = popkern_fillet::FilletOptions {
        continuity: popkern_fillet::ContinuityLevel::G1,
        ..Default::default()
    };
    let result =
        popkern_fillet::fillet_edges(&pyramid, &specs, &options, tol, None).expect("pipeline");
    if let Some(reason) = popkern_fillet::construct_test_hooks::take_setback_topology_checkpoint() {
        eprintln!("decline: {reason}");
    }
    let brep = &result.brep;
    brep.validate().expect("G1 blend must validate");
    eprintln!(
        "G1 commit blend: faces={} edges={} vertices={}",
        brep.faces.len(),
        brep.edges.len(),
        brep.vertices.len()
    );
    let mut kinds = std::collections::BTreeMap::new();
    for (_, f) in brep.faces.iter() {
        let k = match f.surface.as_ref().map(|s| &s.geometry) {
            Some(popkern_core::geometry::SurfaceGeometry::Plane { .. }) => "Plane",
            Some(popkern_core::geometry::SurfaceGeometry::Nurbs(_)) => "Nurbs",
            Some(popkern_core::geometry::SurfaceGeometry::Cylinder { .. }) => "Cylinder",
            Some(popkern_core::geometry::SurfaceGeometry::Torus { .. }) => "Torus",
            Some(_) => "Other",
            None => "None",
        };
        *kinds.entry(k).or_insert(0usize) += 1;
    }
    eprintln!("face kinds: {kinds:?}");
    // Spoke audit: every cap↔cap edge's curve mid vs both faces.
    {
        use popkern_core::geometry::curve::CurveEval;
        use popkern_core::geometry::surface::SurfaceEval;
        use popkern_core::geometry::uv_projection::project_point_to_uv;
        for (eh, edge) in brep.edges.iter() {
            let Some(curve) = edge.curve.as_ref() else {
                continue;
            };
            let he = edge.half_edge;
            let f1 = brep.loops[brep.half_edges[he].loop_].face;
            let f2 = brep.loops[brep.half_edges[brep.half_edges[he].twin].loop_].face;
            let is_nurbs = |fh| {
                matches!(
                    brep.faces[fh].surface.as_ref().map(|s| &s.geometry),
                    Some(popkern_core::geometry::SurfaceGeometry::Nurbs(_))
                )
            };
            if !(is_nurbs(f1) && is_nurbs(f2)) {
                continue;
            }
            let ceval: &dyn CurveEval = &curve.geometry;
            let mid = ceval.point_at(0.5 * (curve.t_start + curve.t_end));
            let d = |fh: popkern_core::topology::FaceHandle| {
                let sfc = &brep.faces[fh].surface.as_ref().unwrap().geometry;
                match project_point_to_uv(sfc, &mid) {
                    Some((u, v)) => sfc.point_at(u, v).distance_to(&mid),
                    None => f64::NAN,
                }
            };
            eprintln!(
                "spoke {eh:?}: mid ({:.3},{:.3},{:.3}) dA {:.2e} dB {:.2e}",
                mid.x,
                mid.y,
                mid.z,
                d(f1),
                d(f2)
            );
        }
    }
    assert_eq!(
        result.corner_patches.len(),
        10,
        "expected the 2n PnS sector faces as corner patches"
    );
    for cp in &result.corner_patches {
        assert_eq!(cp.achieved_continuity, popkern_fillet::ContinuityLevel::G1);
    }
}
