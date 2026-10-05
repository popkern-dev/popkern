use super::*;

/// Chamfer the specified edges of a BRep solid. See [`BlendError`] for failure modes.
pub fn chamfer_edges(
    brep: &BRep,
    specs: &[ChamferSpec],
    options: &ChamferOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<BlendResult, BlendError> {
    chamfer_edges_internal(brep, specs, options, tolerance, tracking, None)
}

/// Cancellable variant of [`chamfer_edges`] — returns
/// [`BlendError::Cancelled`] within one second of the token firing
/// (§17: clean termination < 1 s, model unchanged), via the shared
/// detach-on-cancel layer (PK-276; same contract as
/// [`crate::fillet_edges_cancellable`], including the
/// [`TrackingContext`] handoff: on cancellation the caller holds a
/// fresh context for the same operation and nothing is recorded).
pub fn chamfer_edges_cancellable(
    brep: &BRep,
    specs: &[ChamferSpec],
    options: &ChamferOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
    cancel: &crate::cancellation::CancellationToken,
) -> Result<BlendResult, BlendError> {
    if cancel.is_cancelled() {
        return Err(BlendError::Cancelled);
    }
    let mut tracking = tracking;
    let worker_tracking: Option<TrackingContext> = tracking.as_deref_mut().map(|ctx| {
        let operation = ctx.operation;
        let scheme = ctx.naming_scheme;
        std::mem::replace(
            ctx,
            TrackingContext::new(operation).with_naming_scheme(scheme),
        )
    });
    let worker_brep = brep.clone();
    let worker_specs = specs.to_vec();
    let worker_options = options.clone();
    let worker_token = cancel.clone();
    match popkern_core::cancellation::run_detach_on_cancel(cancel, move || {
        let mut worker_tracking = worker_tracking;
        let result = chamfer_edges_internal(
            &worker_brep,
            &worker_specs,
            &worker_options,
            tolerance,
            worker_tracking.as_mut(),
            Some(&worker_token),
        );
        (result, worker_tracking)
    }) {
        popkern_core::cancellation::DetachOutcome::Completed((result, recorded)) => {
            if let (Some(caller_ctx), Some(recorded)) = (tracking.take(), recorded) {
                *caller_ctx = recorded;
            }
            result
        }
        popkern_core::cancellation::DetachOutcome::Detached => Err(BlendError::Cancelled),
    }
}

/// Shared chamfer pipeline: cooperative polls between stages and the
/// per-edge loop; `cancel` also rides into the apply phase's per-stage
/// and per-corner polls.
fn chamfer_edges_internal(
    brep: &BRep,
    specs: &[ChamferSpec],
    options: &ChamferOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
    cancel: Option<&cancellation::CancellationToken>,
) -> Result<BlendResult, BlendError> {
    macro_rules! check_cancel {
        () => {
            if let Some(tok) = cancel
                && tok.is_cancelled()
            {
                return Err(BlendError::Cancelled);
            }
        };
    }
    check_cancel!();
    // §10.1 — reject invalid parameters at the API boundary.
    validate::validate_chamfer_specs(brep, specs)?;

    // Capture pre-operation faces, edges and vertices for naming.
    let pre = crate::naming::PreOperationSnapshot::capture(brep);
    let pre_faces: &std::collections::HashSet<FaceHandle> = &pre.faces;

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
        // Per-edge cancellation point (§17).
        check_cancel!();
        match compute_chamfer_geometry(&analyses[i], &spec.kind, brep, &tolerance) {
            Ok(geom) => edge_geometries.push(geom),
            Err(failure) => failures.push((spec.edge, failure)),
        }
    }

    if edge_geometries.is_empty() && !failures.is_empty() {
        return Err(BlendError::partial_failure(
            specs.len(),
            brep.clone(),
            failures,
        ));
    }

    // 3a. Closed rim chains (PK-135): when the chamfered edges form
    // complete outer loops of planar cap faces, the per-edge strip
    // pipeline cannot serve them (each edge would retrim the shared
    // cap independently and the strips never meet at the corners).
    // The chain surgery builds the mitred inset by construction. It
    // only runs when EVERY requested edge belongs to a chain — a
    // mixed request falls through whole to the per-edge path so the
    // operation stays all-or-nothing either way.
    check_cancel!();
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
                // PK-244 serving-boundary flag reconcile (see
                // fillet_pipeline): flags are derived state.
                let _ = result.reconcile_orientation_flags();
                // Check 9 (curve-surface consistency) is part of this
                // serve gate: a chamfer that leaves a stale edge piece
                // off its re-trimmed neighbour passes structural
                // validation and ships an invalid body the NEXT
                // boolean declines on (measured: chamfering a rounded
                // plate's corner-fillet top arcs left the vertical
                // tangent edge's top segment unconsumed, 7.5e-3 off
                // the chamfer surface; the body served here and the
                // following subtract failed its result validation).
                let structural = result
                    .validate()
                    .and_then(|()| result.validate_curve_surface());
                if structural.is_ok() {
                    if let Some(ctx) = tracking {
                        assign_fillet_names(&mut result, ctx, pre_faces);
                        crate::naming::assign_output_edge_vertex_names(
                            &mut result,
                            ctx,
                            &pre,
                            &tolerance,
                        );
                    }
                    // PK-478 stage 2 (see fillet_pipeline exit).
                    result.repair_canonical_forward();
                    return Ok(BlendResult {
                        brep: result,
                        edge_warnings: Vec::new(),
                        corner_patches: Vec::new(),
                        repair_report: Default::default(),
                    });
                }
                // Post-surgery validation failed: transactional decline
                // with the validation evidence (never ship the clone).
                return Err(BlendError::partial_failure(
                    specs.len(),
                    brep.clone(),
                    specs
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
                ));
            }
            Some((cap, reason)) => {
                // Declined chain construction is a hard, transactional
                // error: the pre-batch per-edge path is measured-broken
                // for closed chains, so falling through would ship
                // corrupt topology.
                return Err(BlendError::partial_failure(
                    specs.len(),
                    brep.clone(),
                    specs
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
                ));
            }
        }
    }

    if popkern_core::trace::trace_enabled("fillet") {
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
        &std::collections::BTreeSet::new(),
        crate::types::CornerCapPreference::default(), // mitred corners have no cap
        &crate::types::CornerSetback::default(),
        cancel,
        None, // chamfers report no per-corner progress (mitred corners are cheap)
    )
    .map_err(|apply_error| apply_error.into_blend_error(specs.len(), brep))?;
    // Trihedral mitre corners (PK-378): where all three edges at a
    // valence-3 vertex are chamfered with planar bevels, replace the
    // per-edge cluster with the exact mitre apex. Transactional.
    let mitre = construct::trihedral_mitre_corner::settle_trihedral_mitre_corners(
        brep,
        &mut apply_result.brep,
        &edge_geometries,
        &tolerance,
    )
    .map_err(|failure| BlendError::partial_failure(specs.len(), brep.clone(), vec![failure]))?;
    apply_result.repair_report.trihedral_mitre_corners = mitre.corners_settled;

    if popkern_core::trace::trace_enabled("fillet") {
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
    let _ = construct::clamp_bounded_pcurves(&mut apply_result.brep);
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
        assign_fillet_names(&mut apply_result.brep, ctx, pre_faces);
        crate::naming::assign_output_edge_vertex_names(
            &mut apply_result.brep,
            ctx,
            &pre,
            &tolerance,
        );
    }
    // PK-244 serving-boundary flag reconcile: flags are derived state.
    apply_result.repair_report.orientation_flag_reconciliations +=
        apply_result.brep.reconcile_orientation_flags();

    // Check 9 (curve-surface consistency) arms this serve gate too: a
    // chamfer terminating AT a tangent junction (a corner fillet's
    // top arc, whose ends meet the straight rim edges G1) can leave
    // the neighbour's stale edge piece off the re-trimmed surface —
    // structural validation passes and the invalid body surfaces as
    // the NEXT boolean's decline (measured: rounded plate, corner-arc
    // chamfer 0.03, vertical tangent edge's top segment 7.5e-3 off).
    // Declining here is transactional; construction support for the
    // tangent-termination class is the filed follow-up.
    if let Err(deviation) = apply_result.brep.validate_curve_surface() {
        return Err(BlendError::partial_failure(
            specs.len(),
            brep.clone(),
            specs
                .iter()
                .map(|s| {
                    (
                        s.edge,
                        BlendFailure::InvalidTopology {
                            reason: format!(
                                "chamfer result failed curve-surface validation: {deviation}"
                            ),
                            location: None,
                        },
                    )
                })
                .collect(),
        ));
    }

    let mut served = apply_result.brep;
    // PK-478 stage 2 (see fillet_pipeline exit).
    served.repair_canonical_forward();
    Ok(BlendResult {
        brep: served,
        edge_warnings: apply_result.warnings,
        corner_patches: apply_result.corner_patches,
        repair_report: apply_result.repair_report,
    })
}
