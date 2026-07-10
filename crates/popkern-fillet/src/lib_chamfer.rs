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

    if std::env::var_os("PK_FILLET_TRACE").is_some() {
        eprintln!("  [chamfer] pre-apply faces={}", brep.faces.len());
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
        crate::types::ContinuityLevel::G1,
    )
    .map_err(|edge_failures| BlendError::PartialFailure {
        failed: edge_failures.len(),
        total: specs.len(),
        partial_result: Box::new(brep.clone()),
        failures: edge_failures,
    })?;

    if std::env::var_os("PK_FILLET_TRACE").is_some() {
        eprintln!(
            "  [chamfer] post-apply faces={}",
            apply_result.brep.faces.len()
        );
    }
    // Post-apply pcurve maintenance — the same ladder the fillet
    // path runs. The surgery shortens the edges adjacent to each
    // chamfered edge but leaves their pcurves covering the OLD span;
    // without the recompute a retrimmed face tessellates to its
    // pre-chamfer extent and overlaps the bevel strips (measured: a
    // tilted boolean-cut top face meshed clear across both setbacks
    // — 75 open mesh edges on a 9-face body; flat charts survived
    // the staleness by accident, which is why the bracket suite
    // never caught it).
    // Order matters: stale EDGE WINDOWS first (or the pcurve
    // recompute faithfully reproduces the stale extent from the
    // curve samples), then null stale pcurves, then refill —
    // `recompute_all_pcurves` is FILL-ONLY (existing pcurves are
    // trusted).
    apply_result.brep.fix_stale_edge_windows();
    apply_result.brep.fix_stale_pcurves();
    apply_result.brep.recompute_all_pcurves_lenient();
    let _ = construct::clamp_bounded_pcurves_post_phase5(&mut apply_result.brep);

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
