//! PK-097 KP20-commit gate: G1 fillets with the KP20 (watertight-class)
//! setback cap.
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

#[test]
fn g1_kp20_commit_produces_valid_blend() {
    popkern_fillet::construct_test_hooks::enable_setback_topology_commit();
    let tol = Tolerance::default();
    let pyramid = BRep::n_gonal_pyramid(5, 1.0, 1.2, tol).expect("pyramid");
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
    // Base-corner runoff edges must be the true cylinder∩base
    // ellipse (bulging toward the corner), not a straight chord —
    // the axial runoff march degenerated at symmetric corners and
    // produced visibly flat chamfers (PK-097).
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
            if a.z.abs() > 1e-6 || b.z.abs() > 1e-6 || a.distance_to(&b) < 0.05 {
                continue;
            }
            // Base-plane edges: original pentagon edges are LONG
            // (side ~1.18); the corner runoff spans ~0.156.
            if a.distance_to(&b) > 0.5 {
                continue;
            }
            let mid = ceval.point_at(0.5 * (curve.t_start + curve.t_end));
            let chord_mid = popkern_core::geometry::Point3::new(
                0.5 * (a.x + b.x),
                0.5 * (a.y + b.y),
                0.5 * (a.z + b.z),
            );
            let sagitta = mid.distance_to(&chord_mid);
            assert!(
                sagitta > 0.05 * a.distance_to(&b),
                "base runoff edge must be curved (sagitta {:.4} over span {:.4})",
                sagitta,
                a.distance_to(&b)
            );
            checked += 1;
        }
        assert_eq!(checked, 5, "one runoff edge per base corner");
    }
    assert_eq!(
        *kinds.get("Plane").unwrap_or(&0),
        6,
        "base + 5 laterals (PK-097: the earlier sliver-absorb pass used \
         to consume 3 of them on setback-built blends)"
    );
}
