use std::collections::HashMap;
use std::time::{Duration, Instant};

use super::complex::PolyhedralSplineSurfaceComplex;
use super::complex_weld_scope::{
    PolyhedralBoundaryWeldScope, build_boundary_weld_scope, build_post_repair_weld_quality_scope,
};
use super::coverage::PolyhedralCoverageOptions;
use crate::geometry::Point3;
use crate::topology::{
    BRep, BrepCurveTrimReconstructionOptions, BrepDegenerateEdgeKind, BrepDegenerateEdgeOptions,
    BrepEdgeCurveCompletionDeclineReason, BrepEdgeCurveCompletionOptions,
    BrepEdgeCurveCompletionOutcome, BrepEdgeCurveCompletionReport,
    BrepEdgeSurfaceCompatibilityOptions, BrepExportQualityOptions, BrepRepairActionArtifact,
    BrepRepairPipelineOptions, BrepScopedExportQualityReport, BrepSewingClassificationOptions,
    BrepSewingEdgePair, BrepSewingOptions, BrepSewingProvenance, BrepSnapSewingOptions,
    BrepTopologyOutputPolicy, BrepTrimRepairOptions, EdgeHandle, FaceHandle, HalfEdgeHandle, Shell,
    ShellHandle, Solid, analyze_brep_export_quality_with_pcurve_fidelity,
    analyze_scoped_brep_export_quality, build_brep_repair_action_artifact, edge_segment,
    repair_brep_boundary_edges_transactionally,
};

#[derive(Debug, Clone)]
pub struct PolyhedralWeldedSurfaceComplex {
    pub brep: BRep,
    pub report: PolyhedralBoundaryWeldReport,
    pub repair_actions: BrepRepairActionArtifact,
}

#[derive(Debug, Clone, PartialEq)]
pub struct PolyhedralBoundaryWeldReport {
    pub candidate_pairs: usize,
    pub welded_pairs: usize,
    pub declined_pairs: usize,
    pub remaining_boundary_edges: usize,
    pub max_weld_gap: f64,
    pub trim_repairs: usize,
    pub edge_curve_repairs: usize,
    pub edge_curve_declines: usize,
    pub edge_curve_decline_reasons: Vec<(BrepEdgeCurveCompletionDeclineReason, usize)>,
    pub degenerate_edges: usize,
    pub collapsed_surface_trims: usize,
    pub zero_length_vertex_edges: usize,
    pub ambiguous_degenerate_edges: usize,
    pub topology_errors: usize,
    pub topology_warnings: usize,
    pub export_ready: bool,
    pub expected_open_sheet_boundaries: bool,
    pub export_blockers: usize,
    pub missing_edge_curves: usize,
    pub trim_blockers: usize,
    /// Pre-repair topology scope touched by accepted weld candidates.
    pub scope: PolyhedralBoundaryWeldScope,
    /// Edge-local export-quality diagnostics over the pre-repair weld scope.
    pub pre_repair_scoped_quality: BrepScopedExportQualityReport,
    /// Edge-local export-quality diagnostics over repaired topology touched by repair stages.
    pub post_repair_scoped_quality: BrepScopedExportQualityReport,
    pub timing: PolyhedralBoundaryWeldTiming,
}

/// Per-stage wall-clock timing for PnS boundary welding.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct PolyhedralBoundaryWeldTiming {
    pub assemble_ms: u128,
    pub candidate_match_ms: u128,
    pub repair_pipeline_ms: u128,
    pub export_quality_ms: u128,
    pub repair_artifact_ms: u128,
    pub total_ms: u128,
}

#[derive(Debug, Clone, Copy)]
struct BoundarySegment {
    component: usize,
    edge: EdgeHandle,
    start: Point3,
    end: Point3,
}

struct AssembledSurfaceComplex {
    brep: BRep,
    edge_component: HashMap<EdgeHandle, usize>,
}

pub fn weld_clean_surface_complex_boundaries(
    complex: &PolyhedralSplineSurfaceComplex,
    options: PolyhedralCoverageOptions,
) -> Result<PolyhedralWeldedSurfaceComplex, String> {
    let started_at = Instant::now();
    let assembled = assemble_components(complex)?;
    let assembled_at = Instant::now();
    let brep = assembled.brep;
    let boundaries = brep_boundary_segments(&brep, &assembled.edge_component);
    let mut used = vec![false; boundaries.len()];
    let mut candidates = Vec::new();

    for a in 0..boundaries.len() {
        if used[a] {
            continue;
        }
        let mut best: Option<(usize, f64)> = None;
        for b in (a + 1)..boundaries.len() {
            if used[b] || boundaries[a].component == boundaries[b].component {
                continue;
            }
            let endpoint_gap = reversed_endpoint_gap(boundaries[a], boundaries[b]);
            if endpoint_gap > options.c0_tolerance {
                continue;
            }
            if best.is_none_or(|(_, current)| endpoint_gap < current) {
                best = Some((b, endpoint_gap));
            }
        }
        if let Some((b, gap)) = best {
            used[a] = true;
            used[b] = true;
            let _ = gap;
            candidates.push(BrepSewingEdgePair {
                first: boundaries[a].edge,
                second: boundaries[b].edge,
                provenance: BrepSewingProvenance {
                    first: format!("pns component {}", boundaries[a].component),
                    second: format!("pns component {}", boundaries[b].component),
                },
            });
        }
    }
    let scope = build_boundary_weld_scope(&brep, &candidates, &assembled.edge_component);
    let pre_repair_scoped_quality = analyze_scoped_brep_export_quality(
        &brep,
        &scope.export_quality_scope(),
        export_quality_options(options),
    );
    let candidates_at = Instant::now();

    let result =
        repair_brep_boundary_edges_transactionally(&brep, &candidates, pipeline_options(options))
            .map_err(|failure| {
            format!(
                "post-weld repair failed at {:?}: {:?}",
                failure.stage, failure.reason
            )
        })?;
    let repaired_at = Instant::now();
    let snap_report = result.report.snap_sewing.as_ref();
    let clean_report = snap_report.and_then(|report| report.clean_report.as_ref());
    let trim_repairs = result
        .report
        .trim_repair
        .as_ref()
        .map_or(0, |report| report.repaired_half_edges);
    let edge_curve_repairs = result
        .report
        .edge_curve_completion
        .as_ref()
        .map_or(0, |report| report.completed_edges);
    let edge_curve_declines = result
        .report
        .edge_curve_completion
        .as_ref()
        .map_or(0, |report| report.declined_edges);
    let edge_curve_decline_reasons = result
        .report
        .edge_curve_completion
        .as_ref()
        .map_or_else(Vec::new, edge_curve_decline_reasons);
    let degenerate_report = result.report.degenerate_edges.as_ref();
    let degenerate_edges = degenerate_report.map_or(0, |report| report.degenerate_edges);
    let collapsed_surface_trims =
        degenerate_report.map_or(0, |report| report.collapsed_surface_trims);
    let zero_length_vertex_edges =
        degenerate_report.map_or(0, |report| report.zero_length_vertex_edges);
    let ambiguous_degenerate_edges = degenerate_report.map_or(0, |report| {
        report
            .findings
            .iter()
            .filter(|finding| finding.kind == BrepDegenerateEdgeKind::AmbiguousCollapsedTrim)
            .count()
    });
    let topology_errors = result.report.loop_shell.as_ref().map_or(0, |report| {
        report
            .findings
            .iter()
            .filter(|finding| {
                finding.severity == crate::topology::BrepTopologyDiagnosticSeverity::Error
            })
            .count()
    });
    let topology_warnings = result.report.loop_shell.as_ref().map_or(0, |report| {
        report
            .findings
            .iter()
            .filter(|finding| {
                finding.severity == crate::topology::BrepTopologyDiagnosticSeverity::Warning
            })
            .count()
    });
    let post_repair_scope =
        build_post_repair_weld_quality_scope(&result.brep, &candidates, &result.report);
    let post_repair_scoped_quality = analyze_scoped_brep_export_quality(
        &result.brep,
        &post_repair_scope,
        export_quality_options(options),
    );
    let export_quality = analyze_brep_export_quality_with_pcurve_fidelity(
        &result.brep,
        export_quality_options(options),
        (
            post_repair_scoped_quality.pcurve_blockers,
            post_repair_scoped_quality.max_pcurve_deviation,
        ),
    );
    let quality_at = Instant::now();
    let export_blockers = export_quality.topology_blockers
        + export_quality.trim_blockers
        + export_quality.edge_curve_blockers
        + export_quality.pcurve_blockers
        + usize::from(!export_quality.validation_passed);
    let repair_actions = build_brep_repair_action_artifact(&result.report.actions);
    let actions_at = Instant::now();
    let report = PolyhedralBoundaryWeldReport {
        candidate_pairs: result.report.candidate_pairs,
        welded_pairs: clean_report.map_or(0, |report| report.accepted_pairs),
        declined_pairs: snap_report.map_or(0, |report| {
            report.declined_pairs.len()
                + report
                    .clean_report
                    .as_ref()
                    .map_or(0, |clean| clean.declined_pairs.len())
        }),
        remaining_boundary_edges: result.brep.boundary_edges.len(),
        max_weld_gap: clean_report.map_or(0.0, |report| {
            report.max_endpoint_gap.max(report.max_surface_gap)
        }),
        trim_repairs,
        edge_curve_repairs,
        edge_curve_declines,
        edge_curve_decline_reasons,
        degenerate_edges,
        collapsed_surface_trims,
        zero_length_vertex_edges,
        ambiguous_degenerate_edges,
        topology_errors,
        topology_warnings,
        export_ready: export_quality.is_ready(),
        expected_open_sheet_boundaries: export_quality.boundary_edges > 0
            && export_quality.topology_blockers == 0,
        export_blockers,
        missing_edge_curves: export_quality.missing_edge_curves,
        trim_blockers: export_quality.trim_blockers,
        scope,
        pre_repair_scoped_quality,
        post_repair_scoped_quality,
        timing: PolyhedralBoundaryWeldTiming {
            assemble_ms: elapsed_ms(assembled_at.duration_since(started_at)),
            candidate_match_ms: elapsed_ms(candidates_at.duration_since(assembled_at)),
            repair_pipeline_ms: elapsed_ms(repaired_at.duration_since(candidates_at)),
            export_quality_ms: elapsed_ms(quality_at.duration_since(repaired_at)),
            repair_artifact_ms: elapsed_ms(actions_at.duration_since(quality_at)),
            total_ms: elapsed_ms(actions_at.duration_since(started_at)),
        },
    };
    Ok(PolyhedralWeldedSurfaceComplex {
        brep: result.brep,
        report,
        repair_actions,
    })
}

fn elapsed_ms(duration: Duration) -> u128 {
    duration.as_millis()
}

fn pipeline_options(options: PolyhedralCoverageOptions) -> BrepRepairPipelineOptions {
    let clean_sewing = BrepSewingOptions {
        tolerance: options.c0_tolerance,
        samples: options.samples,
        require_surface_agreement: true,
        validate_output: true,
        allow_nonmanifold_edges: false,
    };
    BrepRepairPipelineOptions {
        classification: BrepSewingClassificationOptions {
            exact_tolerance: options.c0_tolerance,
            snap_tolerance: options.c0_tolerance,
            samples: options.samples,
            require_surface_agreement: true,
            allow_nonmanifold_edges: false,
        },
        snap_sewing: BrepSnapSewingOptions {
            max_vertex_motion: 0.0,
            clean_sewing,
        },
        trim_repair: trim_options(options),
        edge_curve_completion: BrepEdgeCurveCompletionOptions {
            tolerance: options.c0_tolerance,
            samples: options.samples,
            allow_vertex_segments_without_surface: true,
            require_all_completed: false,
            validate_output: true,
        },
        curve_trim_reconstruction: BrepCurveTrimReconstructionOptions {
            tolerance: options.c0_tolerance,
            samples: options.samples,
            replace_existing: false,
            replace_stale_existing: true,
            expand_bounded_surface_domains: true,
            require_all_reconstructed: false,
            validate_output: true,
        },
        edge_surface_compatibility: BrepEdgeSurfaceCompatibilityOptions {
            tolerance: options.c0_tolerance,
            max_anchor_motion: 0.0,
            samples: options.samples,
            require_all_compatible: false,
            validate_output: true,
        },
        degenerate_edges: BrepDegenerateEdgeOptions {
            tolerance: options.c0_tolerance,
            samples: options.samples,
        },
        require_all_pairs_repaired: true,
        complete_edge_curves: true,
        reconstruct_curve_trims: true,
        repair_edge_surfaces: false,
        repair_trims: true,
        diagnose_loop_shell: true,
        topology_output: BrepTopologyOutputPolicy::OpenSheet,
    }
}

fn trim_options(options: PolyhedralCoverageOptions) -> BrepTrimRepairOptions {
    BrepTrimRepairOptions {
        endpoint_tolerance: options.c0_tolerance,
        sampled_tolerance: options.c0_tolerance,
        samples: options.samples,
        validate_output: true,
    }
}

fn export_quality_options(options: PolyhedralCoverageOptions) -> BrepExportQualityOptions {
    BrepExportQualityOptions {
        topology_output: BrepTopologyOutputPolicy::OpenSheet,
        trim: trim_options(options),
        require_edge_curves: true,
        allow_collapsed_trim_edge_curves: true,
        max_pcurve_deviation: Some(options.c0_tolerance),
        max_edge_surface_deviation: Some(options.c0_tolerance),
        validate_brep: true,
    }
}

fn edge_curve_decline_reasons(
    report: &BrepEdgeCurveCompletionReport,
) -> Vec<(BrepEdgeCurveCompletionDeclineReason, usize)> {
    let mut counts = Vec::<(BrepEdgeCurveCompletionDeclineReason, usize)>::new();
    for decision in &report.decisions {
        let BrepEdgeCurveCompletionOutcome::Declined(reason) = decision.outcome else {
            continue;
        };
        if let Some((_, count)) = counts.iter_mut().find(|(kind, _)| *kind == reason) {
            *count += 1;
        } else {
            counts.push((reason, 1));
        }
    }
    counts
}

fn assemble_components(
    complex: &PolyhedralSplineSurfaceComplex,
) -> Result<AssembledSurfaceComplex, String> {
    let mut out = BRep::new(crate::tolerance::Tolerance::default());
    let mut edge_component = HashMap::new();
    let solid = out.solids.insert(Solid {
        shell: ShellHandle::default(),
    });
    let shell = out.shells.insert(Shell {
        face: FaceHandle::default(),
        solid,
    });
    out.solids[solid].shell = shell;
    for (component_index, component) in complex.components.iter().enumerate() {
        copy_component_brep(
            &mut out,
            shell,
            component_index,
            &component.brep,
            &mut edge_component,
        )?;
    }
    Ok(AssembledSurfaceComplex {
        brep: out,
        edge_component,
    })
}

fn copy_component_brep(
    out: &mut BRep,
    shell: ShellHandle,
    component_index: usize,
    src: &BRep,
    edge_component: &mut HashMap<EdgeHandle, usize>,
) -> Result<(), String> {
    let mut vertices = HashMap::new();
    let mut edges = HashMap::new();
    let mut half_edges = HashMap::new();
    let mut loops = HashMap::new();
    let mut faces = HashMap::new();

    for (old, vertex) in &src.vertices {
        vertices.insert(old, out.vertices.insert(vertex.clone()));
    }
    for (old, edge) in &src.edges {
        let new = out.edges.insert(edge.clone());
        edges.insert(old, new);
        edge_component.insert(new, component_index);
        if src.boundary_edges.contains(&old) {
            out.boundary_edges.insert(new);
        }
    }
    for (old, half_edge) in &src.half_edges {
        half_edges.insert(old, out.half_edges.insert(half_edge.clone()));
    }
    for (old, loop_) in &src.loops {
        loops.insert(old, out.loops.insert(loop_.clone()));
    }
    for (old, face) in &src.faces {
        let mut new_face = face.clone();
        new_face.shell = shell;
        let new = out.faces.insert(new_face);
        if out.shells[shell].face == FaceHandle::default() {
            out.shells[shell].face = new;
        }
        faces.insert(old, new);
    }

    for (old, edge) in &src.edges {
        out.edges[edges[&old]].half_edge = half_edges[&edge.half_edge];
    }
    for (old, half_edge) in &src.half_edges {
        let new = half_edges[&old];
        out.half_edges[new].origin = vertices[&half_edge.origin];
        out.half_edges[new].twin = remap_half_edge(half_edge.twin, &half_edges);
        out.half_edges[new].next = half_edges[&half_edge.next];
        out.half_edges[new].prev = half_edges[&half_edge.prev];
        out.half_edges[new].edge = edges[&half_edge.edge];
        out.half_edges[new].loop_ = loops[&half_edge.loop_];
    }
    for (old, loop_) in &src.loops {
        let new = loops[&old];
        out.loops[new].half_edge = half_edges[&loop_.half_edge];
        out.loops[new].face = faces[&loop_.face];
    }
    for (old, face) in &src.faces {
        let new = faces[&old];
        out.faces[new].outer_loop = loops[&face.outer_loop];
        out.faces[new].inner_loops = face.inner_loops.iter().map(|loop_| loops[loop_]).collect();
    }
    Ok(())
}

fn remap_half_edge(
    old: HalfEdgeHandle,
    half_edges: &HashMap<HalfEdgeHandle, HalfEdgeHandle>,
) -> HalfEdgeHandle {
    half_edges.get(&old).copied().unwrap_or_default()
}

fn brep_boundary_segments(
    brep: &BRep,
    edge_component: &HashMap<EdgeHandle, usize>,
) -> Vec<BoundarySegment> {
    brep.boundary_edges
        .iter()
        .filter_map(|&edge| {
            edge_segment(brep, edge).map(|(start, end)| BoundarySegment {
                component: edge_component.get(&edge).copied().unwrap_or_default(),
                edge,
                start,
                end,
            })
        })
        .collect()
}

fn reversed_endpoint_gap(a: BoundarySegment, b: BoundarySegment) -> f64 {
    a.start.distance_to(&b.end).max(a.end.distance_to(&b.start))
}
