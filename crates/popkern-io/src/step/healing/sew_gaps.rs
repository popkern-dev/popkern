//! Sew-gaps healing pass: merges coincident vertices and reconciles edge curves.
//!
//! Adapted from the boolean-merge stitching code. After STEP import, models may
//! have duplicate vertices at the same position and edge curves whose endpoints
//! don't match vertex positions. This pass repairs both defects.

use std::collections::HashMap;

use popkern_core::{
    geometry::curve::CurveEval,
    tolerance::Tolerance,
    topology::{BRep, EdgeHandle, VertexHandle},
};

use crate::error::IoError;
use crate::step::diagnostics::{DiagnosticEntry, DiagnosticLog, Phase};

use super::HealingPass;

/// Healing pass that merges coincident vertices and reconciles edge curves.
///
/// Vertices within geometric tolerance are unified, then edge curve endpoints
/// are adjusted to match the canonical vertex positions.
pub struct SewGapsPass;

impl HealingPass for SewGapsPass {
    fn name(&self) -> &str {
        "sew_gaps"
    }

    fn applies_to(&self, entry: &DiagnosticEntry) -> bool {
        // Applies to topology validation issues (gaps, duplicate vertices).
        entry.phase == Phase::Validate || entry.phase == Phase::MapTopology
    }

    fn attempt_repair(
        &self,
        brep: &mut BRep,
        _entry: &DiagnosticEntry,
        log: &mut DiagnosticLog,
    ) -> Result<bool, IoError> {
        let tol = brep.tolerance;
        let merged = merge_coincident_vertices(brep, &tol);
        let reconciled = reconcile_edge_curves(brep);

        if merged > 0 || reconciled > 0 {
            log.warn(
                Phase::Heal,
                None,
                format!(
                    "sew_gaps: merged {merged} coincident vertices, reconciled {reconciled} edge curves"
                ),
            );
            Ok(true)
        } else {
            Ok(false)
        }
    }
}

/// Merge coincident vertices in a BRep.
///
/// Finds vertex pairs within tolerance, keeps the lower-index vertex, and
/// rewires all half-edge origin references. Returns the number of vertices
/// removed. Only vertices of one piece (joined through its edges) merge:
/// a gap lies within a shell, while coincident vertices of different
/// pieces are where they touch, one copy per piece (a self-touching
/// shell imported as pieces); welding those made the pieces one body
/// again by an origin rewire no glue checks (PK-980, PK-1115).
pub fn merge_coincident_vertices(brep: &mut BRep, tolerance: &Tolerance) -> usize {
    let piece = vertex_pieces(brep);
    let verts: Vec<(VertexHandle, _)> = brep
        .vertices
        .iter()
        .map(|(vh, v)| (vh, v.position))
        .collect();

    let mut remap: HashMap<VertexHandle, VertexHandle> = HashMap::new();
    let n = verts.len();
    for i in 0..n {
        let (vh_i, pos_i) = verts[i];
        if remap.contains_key(&vh_i) {
            continue;
        }
        for &(vh_j, pos_j) in &verts[(i + 1)..n] {
            if remap.contains_key(&vh_j) {
                continue;
            }
            // A vertex no fin uses (an orphan) lies in no piece.
            let one_piece = match (piece.get(&vh_i), piece.get(&vh_j)) {
                (Some(a), Some(b)) => a == b,
                _ => true,
            };
            if one_piece && tolerance.points_coincident(&pos_i, &pos_j) {
                remap.insert(vh_j, vh_i);
            }
        }
    }

    if remap.is_empty() {
        return 0;
    }

    let count = remap.len();

    // Rewire half-edge origins.
    for (_, he) in &mut brep.half_edges {
        if let Some(&canonical) = remap.get(&he.origin) {
            he.origin = canonical;
        }
    }

    // Fix vertex outgoing pointers.
    for heh in brep.half_edges.keys().collect::<Vec<_>>() {
        let origin = brep.half_edges[heh].origin;
        brep.vertices[origin].outgoing = heh;
    }

    // Remove merged vertices.
    for old_vh in remap.keys() {
        brep.vertices.remove(*old_vh);
    }

    count
}

/// Each fin-using vertex's piece: vertices joined through an edge (a
/// fin's origin and its successor's) are one.
fn vertex_pieces(brep: &BRep) -> HashMap<VertexHandle, VertexHandle> {
    let mut parent: HashMap<VertexHandle, VertexHandle> = brep
        .half_edges
        .values()
        .map(|h| (h.origin, h.origin))
        .collect();
    fn root(parent: &mut HashMap<VertexHandle, VertexHandle>, mut v: VertexHandle) -> VertexHandle {
        while let Some(&p) = parent.get(&v) {
            if p == v {
                break;
            }
            let up = parent.get(&p).copied().unwrap_or(p);
            parent.insert(v, up);
            v = p;
        }
        v
    }
    for he in brep.half_edges.values() {
        let Some(next) = brep.half_edges.get(he.next) else {
            continue;
        };
        let (a, b) = (root(&mut parent, he.origin), root(&mut parent, next.origin));
        if a != b {
            parent.insert(a, b);
        }
    }
    let keys: Vec<VertexHandle> = parent.keys().copied().collect();
    keys.into_iter()
        .map(|v| (v, root(&mut parent, v)))
        .collect()
}

/// Reconcile edge curve endpoints with vertex positions.
///
/// For each edge with a curve, checks whether the curve's start and end
/// points match the vertex positions. If not, attempts to swap or refit
/// the curve parameters. Returns the number of edges reconciled.
pub fn reconcile_edge_curves(brep: &mut BRep) -> usize {
    let tolerance = brep.tolerance;
    let edge_handles: Vec<EdgeHandle> = brep.edges.keys().collect();
    let mut reconciled = 0usize;

    for eh in edge_handles {
        let curve = match &brep.edges[eh].curve {
            Some(c) => c.clone(),
            None => continue,
        };

        let he_h = brep.edges[eh].half_edge;
        let he = match brep.half_edges.get(he_h) {
            Some(h) => h,
            None => continue,
        };

        let v_start_pos = match brep.vertices.get(he.origin) {
            Some(v) => v.position,
            None => continue,
        };

        let is_self_twin = he.mate() == Some(he_h);

        let twin = match he.mate().and_then(|m| brep.half_edges.get(m)) {
            Some(t) => t,
            None => continue,
        };
        let v_end_pos = match brep.vertices.get(twin.origin) {
            Some(v) => v.position,
            None => continue,
        };

        let c_start = curve.geometry.point_at(curve.t_start);
        let c_end = curve.geometry.point_at(curve.t_end);

        let gap_start = v_start_pos.distance_to(&c_start);
        let gap_end = v_end_pos.distance_to(&c_end);

        if gap_start > tolerance.linear() || gap_end > tolerance.linear() {
            if is_self_twin {
                // Self-twin edges: remove curve so validation uses vertex fallback.
                brep.edges[eh].curve = None;
                reconciled += 1;
            } else {
                // Try swapping t_start/t_end.
                let alt_gap_start = v_start_pos.distance_to(&c_end);
                let alt_gap_end = v_end_pos.distance_to(&c_start);
                if alt_gap_start <= tolerance.linear() && alt_gap_end <= tolerance.linear() {
                    let mut new_curve = curve.clone();
                    std::mem::swap(&mut new_curve.t_start, &mut new_curve.t_end);
                    brep.edges[eh].curve = Some(new_curve);
                    reconciled += 1;
                } else if let Some(new_curve) =
                    refit_line_curve(&curve, v_start_pos, v_end_pos, &tolerance)
                {
                    brep.edges[eh].curve = Some(new_curve);
                    reconciled += 1;
                }
            }
        }
    }

    reconciled
}

/// Try to refit a line curve's parameters so it evaluates to the given
/// start and end positions. Returns the adjusted curve on success.
fn refit_line_curve(
    curve: &popkern_core::geometry::CurveOnEdge,
    start_pos: popkern_core::geometry::Point3,
    end_pos: popkern_core::geometry::Point3,
    tolerance: &Tolerance,
) -> Option<popkern_core::geometry::CurveOnEdge> {
    use popkern_core::geometry::curve::CurveGeometry;

    match &curve.geometry {
        CurveGeometry::Line { origin, direction } => {
            let dir_len_sq = direction.dot(direction);
            if dir_len_sq < 1e-30 {
                return None;
            }
            let diff_start = start_pos - *origin;
            let t_start = diff_start.dot(direction) / dir_len_sq;
            let diff_end = end_pos - *origin;
            let t_end = diff_end.dot(direction) / dir_len_sq;

            let check_start = curve.geometry.point_at(t_start);
            let check_end = curve.geometry.point_at(t_end);
            if start_pos.distance_to(&check_start) <= tolerance.linear()
                && end_pos.distance_to(&check_end) <= tolerance.linear()
            {
                let mut new_curve = curve.clone();
                new_curve.t_start = t_start;
                new_curve.t_end = t_end;
                Some(new_curve)
            } else {
                None
            }
        }
        _ => None,
    }
}

#[cfg(test)]
#[path = "sew_gaps_tests.rs"]
mod sew_gaps_tests;
