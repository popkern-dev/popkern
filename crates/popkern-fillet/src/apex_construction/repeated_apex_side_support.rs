//! Read-only side-support region planner for repeated-apex trim rewrites.

use popkern_core::topology::{BRep, FaceHandle};

use super::error::ApexBuildError;
use super::repeated_apex_region_candidates::CandidateSetSnapshot;
use super::repeated_apex_region_graph::RegionGraph;
use super::repeated_apex_side_support_bridge::FaceFamilyBridgePrimitiveDerivation;
use super::repeated_apex_side_support_execution::SideSupportRegionExecutionReport;
use super::repeated_apex_side_support_facts::SideSupportRegionFacts;
use super::repeated_apex_side_support_merge::SupportRegionMergePrimitiveDerivation;
use super::repeated_apex_side_support_outside::OutsideOutsideSupportPrimitiveDerivation;
use super::repeated_apex_side_support_plan::SideSupportRegionPlanSelection;
use super::repeated_apex_side_support_plan::SideSupportRegionRefusalReason;
use super::repeated_apex_side_support_result::SideSupportRegionContinuationResult;
use super::repeated_apex_topology_knowledge::TopologyKnowledgeGraph;
use super::repeated_apex_trim_subcomplex::TrimRewriteSubcomplexPlan;

#[derive(Clone, Debug)]
pub(super) struct SideSupportRegionDiagnosticSnapshot {
    knowledge: TopologyKnowledgeGraph,
    facts: SideSupportRegionFacts,
    region_graph: RegionGraph,
    candidates: CandidateSetSnapshot,
    selection: SideSupportRegionPlanSelection,
    execution: SideSupportRegionExecutionReport,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) enum SideSupportRegionContinuationOutcome {
    Applied(SideSupportRegionContinuationResult),
    NotReady,
}

#[derive(Clone, Debug, PartialEq, Eq)]
enum SideSupportRegionContinuationFailure {
    StaleKnowledgeGraph,
    PlannerRefused { reason: String },
    BridgeTopologyOperator { reason: String },
    MergeTopologyAndGeometry { reason: String },
}

impl SideSupportRegionDiagnosticSnapshot {
    pub(super) fn build_read_only(brep: &BRep, plan: &TrimRewriteSubcomplexPlan) -> Self {
        let knowledge = TopologyKnowledgeGraph::build_local(brep, plan);
        let facts = SideSupportRegionFacts::classify(&knowledge);
        let region_graph = RegionGraph::build(&facts);
        let candidates = CandidateSetSnapshot::enumerate(&region_graph);
        let selection = SideSupportRegionPlanSelection::select(&candidates);
        let execution =
            SideSupportRegionExecutionReport::preflight(&selection, &knowledge, brep, true);
        Self {
            knowledge,
            facts,
            region_graph,
            candidates,
            selection,
            execution,
        }
    }

    pub(super) fn validate_fresh(&self, brep: &BRep) -> bool {
        self.knowledge.validate_fresh(brep)
    }

    fn apply_executable_continuation(
        &self,
        brep: &mut BRep,
    ) -> Result<SideSupportRegionContinuationOutcome, SideSupportRegionContinuationFailure> {
        if !self.validate_fresh(brep) {
            return Err(SideSupportRegionContinuationFailure::StaleKnowledgeGraph);
        }
        if let SideSupportRegionPlanSelection::Refused(refusal) = &self.selection
            && self.facts.trim_endpoint_pairs().next().is_some()
        {
            return Err(SideSupportRegionContinuationFailure::PlannerRefused {
                reason: format!("stage={:?} reason={:?}", refusal.stage(), refusal.reason()),
            });
        }
        if let SideSupportRegionPlanSelection::Blocked(_, refusal) = &self.selection
            && refusal.reason()
                == SideSupportRegionRefusalReason::MissingValidatedOutsideOutsideSupportContinuationPrimitive
        {
            let outside_primitive =
                OutsideOutsideSupportPrimitiveDerivation::build(&self.selection, &self.knowledge, brep);
            let topology_report = outside_primitive.topology_operator_report(brep);
            if outside_primitive.is_complete() && topology_report.is_ready() {
                let rollback = brep.clone();
                match outside_primitive.apply_topology_operator(brep) {
                    Ok(result) => {
                        // The core operator authors TOPOLOGY only; the
                        // bridge faces (planar continuations of their
                        // supports, bounded by tangency generators)
                        // must receive verified surface geometry and
                        // pcurves or every downstream cap/support seam
                        // measurement declines the whole apex
                        // (measured: cube_extreme — 12 pcurve-less
                        // seams, native Gregory declined, the legacy
                        // fallback committed g0 = 14.2).
                        // Best-effort: bridge faces that continue a
                        // twin-adjacent surface get it (verified) plus
                        // pcurves; the rest stay for the result-level
                        // support-cell completion, which owns their
                        // geometry (the pipeline that kept cube_extreme
                        // green before the strict candidate acceptance
                        // landed). A refusal here was measured to fail
                        // the whole continuation EARLIER than the old
                        // pipeline ever did.
                        if let Err(reason) =
                            assign_bridge_face_geometry(brep, result.bridge_faces())
                            && crate::direct_apex_trace::enabled()
                        {
                            eprintln!(
                                "[direct-apex-trace] bridge geometry authoring left \
                                 residue for downstream completion: {reason}"
                            );
                        }
                        return Ok(SideSupportRegionContinuationOutcome::Applied(result));
                    }
                    Err(err) => {
                        *brep = rollback;
                        return Err(
                            SideSupportRegionContinuationFailure::BridgeTopologyOperator {
                                reason: err.to_string(),
                            },
                        );
                    }
                }
            }
            return Err(SideSupportRegionContinuationFailure::PlannerRefused {
                reason: format!("stage={:?} reason={:?}", refusal.stage(), refusal.reason()),
            });
        }
        let bridge_primitive = FaceFamilyBridgePrimitiveDerivation::build(&self.selection, brep);
        let topology_report = bridge_primitive.topology_operator_report(brep);
        if bridge_primitive.is_complete() && topology_report.is_ready() {
            let rollback = brep.clone();
            match bridge_primitive.apply_topology_operator(brep) {
                Ok(result) => {
                    return Ok(SideSupportRegionContinuationOutcome::Applied(result));
                }
                Err(err) => {
                    *brep = rollback;
                    return Err(
                        SideSupportRegionContinuationFailure::BridgeTopologyOperator {
                            reason: err,
                        },
                    );
                }
            }
        }
        let merge_primitive = SupportRegionMergePrimitiveDerivation::build(&self.selection, brep);
        if !merge_primitive.exact_geometry_ready() {
            return Ok(SideSupportRegionContinuationOutcome::NotReady);
        }
        let rollback = brep.clone();
        if let Err(err) = merge_primitive.apply_topology_and_geometry(brep) {
            *brep = rollback;
            return Err(
                SideSupportRegionContinuationFailure::MergeTopologyAndGeometry { reason: err },
            );
        }
        Ok(SideSupportRegionContinuationOutcome::Applied(
            SideSupportRegionContinuationResult::not_applied(),
        ))
    }

    pub(super) fn diagnostic_label(&self) -> String {
        format!(
            "{}; {}; {}; {}; {}; {}",
            self.knowledge.diagnostic_label(),
            self.facts.diagnostic_label(),
            self.region_graph.diagnostic_label(),
            self.candidates.diagnostic_label(),
            self.selection.diagnostic_label(),
            self.execution.diagnostic_label()
        )
    }
}

pub(super) fn apply_selected_side_support_continuation(
    brep: &mut BRep,
    subcomplex_plan: &TrimRewriteSubcomplexPlan,
) -> Result<SideSupportRegionContinuationOutcome, ApexBuildError> {
    let diagnostics = SideSupportRegionDiagnosticSnapshot::build_read_only(brep, subcomplex_plan);
    diagnostics
        .apply_executable_continuation(brep)
        .map_err(|failure| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition target-cylinder trim rewrite could not apply \
                 executable side-support continuation: {}. Side-support diagnostics: {}",
                failure.diagnostic_label(),
                diagnostics.diagnostic_label()
            ),
        })
}

/// Author geometry for freshly-bridged faces: each bridge face is a
/// planar continuation of one of its twin-adjacent surfaced faces
/// (bounded by the mediating cylinder's tangency generators, so the
/// support surface contains it exactly). Assign the best-fitting
/// adjacent surface — VERIFIED against every loop vertex and edge-curve
/// midpoint — then author the face's missing pcurves and fill missing
/// pcurves on its twin-adjacent faces (the new seam edges' other
/// sides). A bridge face that fits no adjacent surface is a refusal:
/// the caller rolls the candidate back and records the reason.
fn assign_bridge_face_geometry(brep: &mut BRep, bridge_faces: &[FaceHandle]) -> Result<(), String> {
    use popkern_core::geometry::surface::SurfaceEval;

    let fit_gate = brep.tolerance.linear().max(1.0e-9) * 10.0;
    let mut refresh: std::collections::BTreeSet<FaceHandle> = std::collections::BTreeSet::new();
    for &bridge in bridge_faces {
        let Some(face) = brep.faces.get(bridge) else {
            return Err(format!("bridge face {bridge:?} is stale"));
        };
        if face.surface.is_some() {
            refresh.insert(bridge);
            continue;
        }
        // Loop samples: vertices plus each edge curve's midpoint.
        let Some(loop_) = brep.loops.get(face.outer_loop) else {
            return Err(format!("bridge face {bridge:?} has no outer loop"));
        };
        let start = loop_.half_edge;
        let mut samples = Vec::new();
        let mut neighbors: Vec<FaceHandle> = Vec::new();
        let mut he = start;
        for _ in 0..10_000 {
            let Some(data) = brep.half_edges.get(he) else {
                break;
            };
            samples.push(brep.vertices[data.origin].position);
            if let Some(edge) = brep.edges.get(data.edge)
                && let Some(curve) = edge.curve.as_ref()
            {
                use popkern_core::geometry::curve::CurveEval;
                samples.push(curve.geometry.point_at(0.5 * (curve.t_start + curve.t_end)));
            }
            if let Some(twin) = brep.half_edges.get(data.twin)
                && let Some(twin_loop) = brep.loops.get(twin.loop_)
                && twin_loop.face != bridge
                && !neighbors.contains(&twin_loop.face)
            {
                neighbors.push(twin_loop.face);
            }
            he = data.next;
            if he == start {
                break;
            }
        }
        let mut best: Option<(f64, FaceHandle)> = None;
        for neighbor in &neighbors {
            let Some(sof) = brep.faces.get(*neighbor).and_then(|f| f.surface.as_ref()) else {
                continue;
            };
            let mut worst = 0.0f64;
            for p in &samples {
                let d =
                    popkern_core::geometry::uv_projection::project_point_to_uv(&sof.geometry, p)
                        .map(|(u, v)| sof.geometry.point_at(u, v).distance_to(p))
                        .unwrap_or(f64::INFINITY);
                worst = worst.max(d);
                if worst > fit_gate {
                    break;
                }
            }
            if worst <= fit_gate && best.is_none_or(|(bd, _)| worst < bd) {
                best = Some((worst, *neighbor));
            }
        }
        let Some((_, source)) = best else {
            if crate::direct_apex_trace::enabled() {
                for neighbor in &neighbors {
                    let Some(sof) = brep.faces.get(*neighbor).and_then(|f| f.surface.as_ref())
                    else {
                        continue;
                    };
                    let mut worst = 0.0f64;
                    for p in &samples {
                        let d = popkern_core::geometry::uv_projection::project_point_to_uv(
                            &sof.geometry,
                            p,
                        )
                        .map(|(u, v)| sof.geometry.point_at(u, v).distance_to(p))
                        .unwrap_or(f64::INFINITY);
                        worst = worst.max(d);
                    }
                    eprintln!(
                        "[direct-apex-trace] bridge {bridge:?} vs neighbor {neighbor:?}: worst={worst:.4e}"
                    );
                }
                // Planarity of the sample set: best-fit plane residual
                // via the centroid + Newell normal.
                let n_pts = samples.len().max(1) as f64;
                let mut c = popkern_core::geometry::Vec3::new(0.0, 0.0, 0.0);
                for p in &samples {
                    c = c + popkern_core::geometry::Vec3::new(p.x, p.y, p.z);
                }
                let c = c * (1.0 / n_pts);
                let mut normal = popkern_core::geometry::Vec3::new(0.0, 0.0, 0.0);
                for i in 0..samples.len() {
                    let a = samples[i];
                    let b = samples[(i + 1) % samples.len()];
                    let a = popkern_core::geometry::Vec3::new(a.x - c.x, a.y - c.y, a.z - c.z);
                    let b = popkern_core::geometry::Vec3::new(b.x - c.x, b.y - c.y, b.z - c.z);
                    normal = normal + a.cross(&b);
                }
                let flat = normal.normalized().map(|n| {
                    samples
                        .iter()
                        .map(|p| {
                            (popkern_core::geometry::Vec3::new(p.x - c.x, p.y - c.y, p.z - c.z)
                                .dot(&n))
                            .abs()
                        })
                        .fold(0.0, f64::max)
                });
                eprintln!(
                    "[direct-apex-trace] bridge {bridge:?} planarity worst={flat:?} samples={}",
                    samples.len()
                );
            }
            return Err(format!(
                "bridge face {bridge:?} fits no twin-adjacent surface within {fit_gate:.2e} \
                 over {} samples ({} neighbors tried)",
                samples.len(),
                neighbors.len()
            ));
        };
        let source_surface = brep
            .faces
            .get(source)
            .and_then(|f| f.surface.clone())
            .ok_or_else(|| format!("source face {source:?} lost its surface"))?;
        brep.faces[bridge].surface = Some(source_surface);
        refresh.insert(bridge);
        for neighbor in neighbors {
            refresh.insert(neighbor);
        }
    }
    for face in refresh {
        brep.recompute_face_pcurves(face)
            .map_err(|err| format!("pcurve authoring on {face:?} failed: {err}"))?;
    }
    Ok(())
}

impl SideSupportRegionContinuationOutcome {
    pub(super) fn applied_diagnostic(
        &self,
        cap_face: FaceHandle,
    ) -> Option<crate::types::ApexConstructionDiagnostic> {
        let SideSupportRegionContinuationOutcome::Applied(result) = self else {
            return None;
        };
        if !result.was_applied() {
            return None;
        }
        Some(
            crate::types::ApexConstructionDiagnostic::SideSupportContinuation(
                crate::types::SideSupportContinuationDiagnostic {
                    stage: "side_support_continuation".to_string(),
                    code: format!(
                        "side_support_continuation.{}.applied",
                        result.primitive_code()
                    ),
                    cap_face,
                    primitive: result.primitive_code().to_string(),
                    applied: true,
                    bridge_face_count: result.bridge_face_count(),
                    bridge_loop_count: result.bridge_loop_count(),
                },
            ),
        )
    }

    pub(super) fn diagnostic_label(&self) -> String {
        match self {
            SideSupportRegionContinuationOutcome::Applied(result) => {
                if result.was_applied() {
                    format!("applied-with-support-payload {}", result.diagnostic_label())
                } else {
                    "applied".to_string()
                }
            }
            SideSupportRegionContinuationOutcome::NotReady => "not-ready".to_string(),
        }
    }

    pub(super) fn bridge_loops(&self) -> &[popkern_core::topology::LoopHandle] {
        match self {
            SideSupportRegionContinuationOutcome::Applied(result) => result.bridge_loops(),
            SideSupportRegionContinuationOutcome::NotReady => &[],
        }
    }
}

impl SideSupportRegionContinuationFailure {
    pub(super) fn diagnostic_label(&self) -> String {
        match self {
            SideSupportRegionContinuationFailure::StaleKnowledgeGraph => {
                "stale-knowledge-graph".to_string()
            }
            SideSupportRegionContinuationFailure::PlannerRefused { reason } => {
                format!("planner-refused {reason}")
            }
            SideSupportRegionContinuationFailure::BridgeTopologyOperator { reason } => {
                format!("bridge-topology-operator-failed reason={reason}")
            }
            SideSupportRegionContinuationFailure::MergeTopologyAndGeometry { reason } => {
                format!("merge-topology-and-geometry-failed reason={reason}")
            }
        }
    }
}

#[cfg(test)]
#[path = "repeated_apex_side_support_tests.rs"]
mod tests;
