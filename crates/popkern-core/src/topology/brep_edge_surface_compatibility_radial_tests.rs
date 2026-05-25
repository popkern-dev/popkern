use super::*;

#[test]
fn detects_radial_face_use_edge_surface_incompatibility() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let edge = first_repairable_edge(&brep);
    let base = brep.edges[edge].half_edge;
    let twin = brep.half_edges[base].twin;
    let radial_use = brep.half_edges.insert(HalfEdge {
        origin: brep.half_edges[base].origin,
        twin,
        next: brep.half_edges[base].next,
        prev: brep.half_edges[base].prev,
        edge,
        loop_: brep.half_edges[base].loop_,
        pcurve: brep.half_edges[base].pcurve.clone(),
    });
    brep.nonmanifold_radial.insert(base, twin);
    brep.nonmanifold_radial.insert(twin, radial_use);
    brep.nonmanifold_radial.insert(radial_use, base);
    corrupt_edge_curve(&mut brep, edge, Vec3::new(0.11, 0.07, 0.13));

    let report = analyze_brep_edge_surface_compatibility(&brep, options());

    assert!(
        report.inspected_half_edges >= 3,
        "radial face uses must be included in edge/surface validation: {report:#?}"
    );
    assert!(
        report
            .findings
            .iter()
            .any(|finding| finding.half_edge == radial_use),
        "radial face-use incompatibility must not be hidden behind manifold twin traversal: {report:#?}"
    );
    let pcurve_report = brep
        .measure_edge_pcurve_deviation(edge)
        .expect("corrupted radial edge remains measurable");
    assert!(
        pcurve_report
            .half_edges
            .iter()
            .any(|finding| finding.half_edge == radial_use),
        "radial face-use pcurves must be included in pcurve fidelity measurement: {pcurve_report:#?}"
    );
}
