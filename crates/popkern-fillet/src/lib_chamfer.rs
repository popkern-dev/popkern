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

    // 3a. Closed rim chains (PK-135): when the chamfered edges form
    // complete outer loops of planar cap faces, the per-edge strip
    // pipeline cannot serve them (each edge would retrim the shared
    // cap independently and the strips never meet at the corners).
    // The chain surgery builds the mitred inset by construction. It
    // only runs when EVERY requested edge belongs to a chain — a
    // mixed request falls through whole to the per-edge path so the
    // operation stays all-or-nothing either way.
    let chains = construct::rim_chain_chamfer::detect_closed_rim_chains(brep, &edge_geometries);
    let chain_edge_count: usize = chains.iter().map(|c| c.geometry_indices().count()).sum();
    if !chains.is_empty() && chain_edge_count == edge_geometries.len() && failures.is_empty() {
        let mut result = brep.clone();
        let mut chain_error: Option<(FaceHandle, String)> = None;
        for chain in &chains {
            if let Err(reason) = construct::rim_chain_chamfer::apply_rim_chain_chamfer(
                &mut result,
                chain,
                &edge_geometries,
                &tolerance,
            ) {
                chain_error = Some((chain.cap_face(), reason));
                break;
            }
        }
        match chain_error {
            None => {
                result.fix_stale_pcurves();
                result.recompute_all_pcurves_lenient();
                let structural = result.validate();
                if structural.is_ok() {
                    if let Some(ctx) = tracking {
                        assign_fillet_names(&mut result, ctx, &pre_faces);
                    }
                    return Ok(BlendResult {
                        brep: result,
                        edge_warnings: Vec::new(),
                        corner_patches: Vec::new(),
                        repair_report: Default::default(),
                    });
                }
                // Post-surgery validation failed: transactional decline
                // with the validation evidence (never ship the clone).
                return Err(BlendError::PartialFailure {
                    failed: specs.len(),
                    total: specs.len(),
                    partial_result: Box::new(brep.clone()),
                    failures: specs
                        .iter()
                        .map(|s| {
                            (
                                s.edge,
                                BlendFailure::InvalidTopology {
                                    reason: format!(
                                        "rim chain chamfer failed post-validation: {}",
                                        structural
                                            .as_ref()
                                            .err()
                                            .map(|e| e.to_string())
                                            .unwrap_or_default()
                                    ),
                                    location: None,
                                },
                            )
                        })
                        .collect(),
                });
            }
            Some((cap, reason)) => {
                // Declined chain construction is a hard, transactional
                // error: the legacy per-edge path is measured-broken
                // for closed chains, so falling through would ship
                // corrupt topology.
                return Err(BlendError::PartialFailure {
                    failed: specs.len(),
                    total: specs.len(),
                    partial_result: Box::new(brep.clone()),
                    failures: specs
                        .iter()
                        .map(|s| {
                            (
                                s.edge,
                                BlendFailure::InvalidTopology {
                                    reason: format!(
                                        "rim chain chamfer around {cap:?} declined: {reason}"
                                    ),
                                    location: None,
                                },
                            )
                        })
                        .collect(),
                });
            }
        }
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
    // Terminal chord trim (the planar counterpart of PK-148): settle
    // chain-end corners so the band's end stubs carry ONE geometric
    // authority — without it a chamfered cube round-trips through
    // STEP with forked face-local end edges (χ=6 vs 2).
    apply_result.repair_report.terminal_arc_trim =
        construct::terminal_arc_trim::settle_terminal_arc_stubs(
            brep,
            &mut apply_result.brep,
            &tolerance,
            construct::terminal_arc_trim::PlanarStubRecognition::SurfaceScoped,
        );
    // Second reconcile pass AFTER the terminal settle: the settle
    // snaps anchor vertices onto their true triple points, which is
    // exactly what the first pass declined on (AnchorOffSurface on
    // every stale gap chord of an interrupted rim). The stale-only
    // gate keeps already-exact curves untouched, so the second pass
    // only serves the edges the settle just made reconcilable.
    construct::edge_reconcile::reconcile_edge_curves(&mut apply_result.brep);

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
