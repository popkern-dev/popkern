use super::*;

/// Chamfer the specified edges of a BRep solid. See [`BlendError`] for failure modes.
pub fn chamfer_edges(
    brep: &BRep,
    specs: &[ChamferSpec],
    options: &ChamferOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<BlendResult, BlendError> {
    // §10.1 — reject invalid parameters at the API boundary.
    validate::validate_chamfer_specs(brep, specs)?;

    // Capture pre-operation faces for naming.
    let pre_faces: std::collections::HashSet<FaceHandle> = brep.faces.keys().collect();

    // 1. Analyze edges.
    let edges: Vec<EdgeHandle> = specs.iter().map(|s| s.edge).collect();
    let analyses = analysis::analyze_edges(brep, &edges, options.propagate_tangent)
        .map_err(BlendError::NoValidEdges)?;

    // 2. Compute geometry for each edge.
    let mut edge_geometries = Vec::new();
    let mut failures = Vec::new();

    for (i, spec) in specs.iter().enumerate() {
        if i >= analyses.len() {
            continue;
        }
        match compute_chamfer_geometry(&analyses[i], &spec.kind, brep, &tolerance) {
            Ok(geom) => edge_geometries.push(geom),
            Err(failure) => failures.push((spec.edge, failure)),
        }
    }

    if edge_geometries.is_empty() && !failures.is_empty() {
        return Err(BlendError::PartialFailure {
            failed: failures.len(),
            total: specs.len(),
            partial_result: Box::new(brep.clone()),
            failures,
        });
    }

    // 3. Apply topology. Chamfers always use mitred vertex mode + no rational-Gregory.
    let mut apply_result = construct::apply_blend_geometry(
        brep,
        &edge_geometries,
        &tolerance,
        &VertexMode::Mitred,
        false,
        false,
        false,
        false,
        false,
        false,
        false,
        false, // chamfers don't go through the DirectApex apex builder
    )
    .map_err(|edge_failures| BlendError::PartialFailure {
        failed: edge_failures.len(),
        total: specs.len(),
        partial_result: Box::new(brep.clone()),
        failures: edge_failures,
    })?;

    // Assign stable names when tracking is active.
    if let Some(ctx) = tracking {
        assign_fillet_names(&mut apply_result.brep, ctx, &pre_faces);
    }

    Ok(BlendResult {
        brep: apply_result.brep,
        edge_warnings: apply_result.warnings,
        corner_patches: apply_result.corner_patches,
        repair_report: apply_result.repair_report,
    })
}
