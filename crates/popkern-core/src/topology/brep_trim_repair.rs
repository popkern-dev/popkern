//! Trim and pcurve diagnostics plus transactional repair for BReps.

use std::collections::HashSet;

use crate::geometry::SurfaceEval;
use crate::topology::{BRep, FaceHandle, HalfEdgeHandle};

/// Options controlling trim diagnostics and repair.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BrepTrimRepairOptions {
    /// Maximum allowed 3D mismatch between pcurve/surface and trim vertices.
    pub endpoint_tolerance: f64,
    /// Maximum allowed sampled mismatch between edge curves and pcurve/surface.
    pub sampled_tolerance: f64,
    /// Number of samples per half-edge for sampled trim checks.
    pub samples: usize,
    /// Whether repaired output must pass full BRep validation.
    pub validate_output: bool,
}

/// Aggregate trim diagnostic report.
#[derive(Debug, Clone, PartialEq)]
pub struct BrepTrimDiagnosticReport {
    /// Number of half-edges inspected on surfaced faces.
    pub inspected_half_edges: usize,
    /// Number of half-edges with missing pcurves.
    pub missing_pcurves: usize,
    /// Number of half-edges with stale pcurve endpoint mappings.
    pub stale_endpoint_pcurves: usize,
    /// Number of half-edges with sampled edge/surface mismatch.
    pub sampled_mismatches: usize,
    /// Maximum endpoint mismatch observed.
    pub max_endpoint_deviation: f64,
    /// Maximum sampled edge/surface mismatch observed.
    pub max_sampled_deviation: f64,
    /// Per-half-edge findings.
    pub findings: Vec<BrepTrimFinding>,
}

impl BrepTrimDiagnosticReport {
    /// True when no repairable trim findings were produced.
    pub fn is_clean(&self) -> bool {
        self.findings.is_empty()
    }
}

/// One trim diagnostic finding.
#[derive(Debug, Clone, PartialEq)]
pub struct BrepTrimFinding {
    /// Face carrying the trim.
    pub face: FaceHandle,
    /// Half-edge carrying the pcurve.
    pub half_edge: HalfEdgeHandle,
    /// Finding kind.
    pub kind: BrepTrimFindingKind,
    /// Maximum deviation attached to this finding.
    pub deviation: f64,
}

/// Machine-readable trim finding kind.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BrepTrimFindingKind {
    /// A surfaced face half-edge has no pcurve.
    MissingPcurve,
    /// Pcurve endpoints do not map to trim vertices on the face surface.
    StaleEndpoint,
    /// Interior samples disagree between edge curve and pcurve/surface.
    SampledMismatch,
}

/// Successful transactional trim repair result.
#[derive(Debug, Clone)]
pub struct BrepTrimRepairResult {
    /// Repaired clone of the input BRep.
    pub brep: BRep,
    /// Repair report.
    pub report: BrepTrimRepairReport,
}

/// Aggregate trim repair report.
#[derive(Debug, Clone, PartialEq)]
pub struct BrepTrimRepairReport {
    /// Diagnostics before repair.
    pub before: BrepTrimDiagnosticReport,
    /// Diagnostics after repair.
    pub after: BrepTrimDiagnosticReport,
    /// Number of faces whose pcurves were recomputed.
    pub repaired_faces: usize,
    /// Number of half-edges cleared and recomputed.
    pub repaired_half_edges: usize,
}

/// Transactional trim repair failure.
#[derive(Debug, Clone, PartialEq)]
pub struct BrepTrimRepairFailure {
    /// Failure reason.
    pub reason: BrepTrimRepairFailureReason,
    /// Report accumulated before failure.
    pub report: BrepTrimRepairReport,
}

/// Reason trim repair failed.
#[derive(Debug, Clone, PartialEq)]
pub enum BrepTrimRepairFailureReason {
    /// Recomputed trims still exceed tolerance.
    UnresolvedTrimFindings,
    /// Repaired clone failed validation.
    ValidationFailed(String),
}

/// Analyze pcurves on surfaced faces without mutating the BRep.
pub fn analyze_brep_trims(brep: &BRep, options: BrepTrimRepairOptions) -> BrepTrimDiagnosticReport {
    let mut report = BrepTrimDiagnosticReport {
        inspected_half_edges: 0,
        missing_pcurves: 0,
        stale_endpoint_pcurves: 0,
        sampled_mismatches: 0,
        max_endpoint_deviation: 0.0,
        max_sampled_deviation: 0.0,
        findings: Vec::new(),
    };

    for (face_handle, face) in &brep.faces {
        if face.surface.is_none() {
            continue;
        }
        for half_edge in face_half_edges(brep, face_handle) {
            report.inspected_half_edges += 1;
            inspect_half_edge_trim(brep, face_handle, half_edge, options, &mut report);
        }
    }
    report
}

/// Repair missing or stale pcurves on a cloned BRep.
pub fn repair_brep_trims_transactionally(
    brep: &BRep,
    options: BrepTrimRepairOptions,
) -> Result<BrepTrimRepairResult, Box<BrepTrimRepairFailure>> {
    let before = analyze_brep_trims(brep, options);
    let mut repaired = brep.clone();
    let mut affected_faces = HashSet::new();
    let mut affected_half_edges = HashSet::new();

    for finding in &before.findings {
        affected_faces.insert(finding.face);
        affected_half_edges.insert(finding.half_edge);
        if let Some(half_edge) = repaired.half_edges.get_mut(finding.half_edge) {
            half_edge.pcurve = None;
        }
    }

    for face in &affected_faces {
        repaired.recompute_face_pcurves(*face).map_err(|err| {
            trim_failure(
                BrepTrimRepairFailureReason::ValidationFailed(err.to_string()),
                before.clone(),
                &repaired,
                options,
                affected_faces.len(),
                affected_half_edges.len(),
            )
        })?;
    }
    repaired.fix_stale_pcurves();

    let after = analyze_brep_trims(&repaired, options);
    let report = BrepTrimRepairReport {
        before,
        after,
        repaired_faces: affected_faces.len(),
        repaired_half_edges: affected_half_edges.len(),
    };
    if !report.after.is_clean() {
        return Err(Box::new(BrepTrimRepairFailure {
            reason: BrepTrimRepairFailureReason::UnresolvedTrimFindings,
            report,
        }));
    }
    if options.validate_output
        && let Err(err) = repaired.validate()
    {
        return Err(Box::new(BrepTrimRepairFailure {
            reason: BrepTrimRepairFailureReason::ValidationFailed(err.to_string()),
            report,
        }));
    }
    Ok(BrepTrimRepairResult {
        brep: repaired,
        report,
    })
}

fn inspect_half_edge_trim(
    brep: &BRep,
    face: FaceHandle,
    half_edge: HalfEdgeHandle,
    options: BrepTrimRepairOptions,
    report: &mut BrepTrimDiagnosticReport,
) {
    let Some(he) = brep.half_edges.get(half_edge) else {
        return;
    };
    if he.pcurve.is_none() {
        report.missing_pcurves += 1;
        report.findings.push(BrepTrimFinding {
            face,
            half_edge,
            kind: BrepTrimFindingKind::MissingPcurve,
            deviation: f64::INFINITY,
        });
        return;
    }
    let endpoint_deviation =
        trim_endpoint_deviation(brep, face, half_edge).unwrap_or(f64::INFINITY);
    report.max_endpoint_deviation = report.max_endpoint_deviation.max(endpoint_deviation);
    if endpoint_deviation > options.endpoint_tolerance {
        report.stale_endpoint_pcurves += 1;
        report.findings.push(BrepTrimFinding {
            face,
            half_edge,
            kind: BrepTrimFindingKind::StaleEndpoint,
            deviation: endpoint_deviation,
        });
    }
    let sampled_deviation = trim_sampled_deviation(brep, half_edge, options.samples).unwrap_or(0.0);
    report.max_sampled_deviation = report.max_sampled_deviation.max(sampled_deviation);
    if sampled_deviation > options.sampled_tolerance {
        report.sampled_mismatches += 1;
        report.findings.push(BrepTrimFinding {
            face,
            half_edge,
            kind: BrepTrimFindingKind::SampledMismatch,
            deviation: sampled_deviation,
        });
    }
}

fn face_half_edges(brep: &BRep, face: FaceHandle) -> Vec<HalfEdgeHandle> {
    let Some(face_data) = brep.faces.get(face) else {
        return Vec::new();
    };
    let mut loops = vec![face_data.outer_loop];
    loops.extend(face_data.inner_loops.iter().copied());
    let mut out = Vec::new();
    for loop_handle in loops {
        let Some(loop_data) = brep.loops.get(loop_handle) else {
            continue;
        };
        let start = loop_data.half_edge;
        if !brep.half_edges.contains_key(start) {
            continue;
        }
        let mut cursor = start;
        let bound = brep.half_edges.len().max(1);
        for _ in 0..=bound {
            out.push(cursor);
            cursor = brep.half_edges[cursor].next;
            if cursor == start || !brep.half_edges.contains_key(cursor) {
                break;
            }
        }
    }
    out
}

fn trim_endpoint_deviation(
    brep: &BRep,
    face: FaceHandle,
    half_edge: HalfEdgeHandle,
) -> Option<f64> {
    let he = brep.half_edges.get(half_edge)?;
    let pcurve = he.pcurve.as_ref()?;
    let surface = &brep.faces.get(face)?.surface.as_ref()?.geometry;
    let origin = brep.vertices.get(he.origin)?.position;
    let dest = brep
        .vertices
        .get(brep.half_edges.get(he.next)?.origin)?
        .position;
    let (t0, t1) = pcurve.parameter_range();
    let start_uv = pcurve.point_at(t0);
    let end_uv = pcurve.point_at(t1);
    let start = surface.point_at(start_uv.x, start_uv.y);
    let end = surface.point_at(end_uv.x, end_uv.y);
    Some(start.distance_to(&origin).max(end.distance_to(&dest)))
}

fn trim_sampled_deviation(brep: &BRep, half_edge: HalfEdgeHandle, _samples: usize) -> Option<f64> {
    brep.measure_half_edge_pcurve_deviation(half_edge)
        .map(|report| report.max_deviation)
}

fn trim_failure(
    reason: BrepTrimRepairFailureReason,
    before: BrepTrimDiagnosticReport,
    brep: &BRep,
    options: BrepTrimRepairOptions,
    repaired_faces: usize,
    repaired_half_edges: usize,
) -> Box<BrepTrimRepairFailure> {
    Box::new(BrepTrimRepairFailure {
        reason,
        report: BrepTrimRepairReport {
            before,
            after: analyze_brep_trims(brep, options),
            repaired_faces,
            repaired_half_edges,
        },
    })
}

#[cfg(test)]
#[path = "brep_trim_repair_tests.rs"]
mod tests;
