//! PK-097 Phase C gate: G1 fillets with the KP20 (watertight-class)
//! setback cap.
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

#[test]
#[ignore = "PK-097: G1 production route (KP20 cap family)"]
fn g1_kp20_commit_produces_valid_blend() {
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
    brep.validate().expect("G1/KP20 blend must validate");
    let mut kinds = std::collections::BTreeMap::new();
    for (_, f) in brep.faces.iter() {
        let k = match f.surface.as_ref().map(|s| &s.geometry) {
            Some(popkern_core::geometry::SurfaceGeometry::Plane { .. }) => "Plane",
            Some(popkern_core::geometry::SurfaceGeometry::Nurbs(_)) => "Nurbs",
            Some(popkern_core::geometry::SurfaceGeometry::Cylinder { .. }) => "Cylinder",
            Some(popkern_core::geometry::SurfaceGeometry::Gregory(_)) => "Gregory",
            Some(_) => "Other",
            None => "None",
        };
        *kinds.entry(k).or_insert(0usize) += 1;
    }
    eprintln!(
        "G1/KP20 blend: faces={} edges={} vertices={} kinds={kinds:?}",
        brep.faces.len(),
        brep.edges.len(),
        brep.vertices.len()
    );
    assert_eq!(result.corner_patches.len(), 20, "the 4n KP20 piece faces");
    assert_eq!(*kinds.get("Gregory").unwrap_or(&0), 20);
    assert_eq!(*kinds.get("Cylinder").unwrap_or(&0), 5, "5 strips");
    assert_eq!(
        *kinds.get("Plane").unwrap_or(&0),
        6,
        "base + 5 laterals (PK-097: the legacy sliver-absorb pass used \
         to consume 3 of them on setback-built blends)"
    );
}
