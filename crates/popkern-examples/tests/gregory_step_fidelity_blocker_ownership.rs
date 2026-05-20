use popkern_core::topology::BRep;

pub(crate) fn blocker_ownership_summary(brep: &BRep, limit: usize) -> Vec<String> {
    let mut reports = brep
        .edges
        .keys()
        .filter_map(|edge| brep.measure_edge_pcurve_deviation(edge))
        .filter(|report| report.max_deviation > 1e-12)
        .collect::<Vec<_>>();
    reports.sort_by(|a, b| b.max_deviation.total_cmp(&a.max_deviation));
    reports
        .into_iter()
        .take(limit)
        .map(|edge| {
            let role = brep
                .edges
                .get(edge.edge)
                .map(|edge| format!("{:?}", edge.role))
                .unwrap_or_else(|| "unknown".to_string());
            let faces = edge
                .half_edges
                .iter()
                .map(|half_edge| {
                    let face_role = brep
                        .faces
                        .get(half_edge.face)
                        .map(|face| format!("{:?}", face.role))
                        .unwrap_or_else(|| "unknown".to_string());
                    format!(
                        "{:?} on {:?} role={} deviation={:.6}",
                        half_edge.half_edge, half_edge.face, face_role, half_edge.max_deviation
                    )
                })
                .collect::<Vec<_>>()
                .join(", ");
            format!(
                "{:?} role={} deviation={:.6}: {}",
                edge.edge, role, edge.max_deviation, faces
            )
        })
        .collect()
}
