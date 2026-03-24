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
/// removed.
pub fn merge_coincident_vertices(brep: &mut BRep, tolerance: &Tolerance) -> usize {
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
            if tolerance.points_coincident(&pos_i, &pos_j) {
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

        let is_self_twin = he.twin == he_h;

        let twin = match brep.half_edges.get(he.twin) {
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

        if gap_start > tolerance.linear || gap_end > tolerance.linear {
            if is_self_twin {
                // Self-twin edges: remove curve so validation uses vertex fallback.
                brep.edges[eh].curve = None;
                reconciled += 1;
            } else {
                // Try swapping t_start/t_end.
                let alt_gap_start = v_start_pos.distance_to(&c_end);
                let alt_gap_end = v_end_pos.distance_to(&c_start);
                if alt_gap_start <= tolerance.linear && alt_gap_end <= tolerance.linear {
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
            if start_pos.distance_to(&check_start) <= tolerance.linear
                && end_pos.distance_to(&check_end) <= tolerance.linear
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
mod tests {
    use super::*;
    use popkern_core::{
        geometry::{Point3, Vec3},
        topology::{BRep, Vertex},
    };

    fn tol() -> Tolerance {
        Tolerance::default()
    }

    #[test]
    fn merge_removes_coincident_vertices() {
        let mut brep = BRep::new(tol());
        let _v0 = brep.vertices.insert(Vertex {
            position: Point3::new(0.0, 0.0, 0.0),
            outgoing: Default::default(),
        });
        let _v1 = brep.vertices.insert(Vertex {
            position: Point3::new(0.0, 0.0, 0.0),
            outgoing: Default::default(),
        });
        let _v2 = brep.vertices.insert(Vertex {
            position: Point3::new(1.0, 0.0, 0.0),
            outgoing: Default::default(),
        });

        assert_eq!(brep.vertices.len(), 3);
        let removed = merge_coincident_vertices(&mut brep, &tol());
        assert_eq!(removed, 1);
        assert_eq!(brep.vertices.len(), 2);
    }

    #[test]
    fn merge_noop_when_no_duplicates() {
        let mut brep = BRep::new(tol());
        let _v0 = brep.vertices.insert(Vertex {
            position: Point3::new(0.0, 0.0, 0.0),
            outgoing: Default::default(),
        });
        let _v1 = brep.vertices.insert(Vertex {
            position: Point3::new(1.0, 0.0, 0.0),
            outgoing: Default::default(),
        });

        let removed = merge_coincident_vertices(&mut brep, &tol());
        assert_eq!(removed, 0);
        assert_eq!(brep.vertices.len(), 2);
    }

    #[test]
    fn sew_gaps_pass_name() {
        let pass = SewGapsPass;
        assert_eq!(pass.name(), "sew_gaps");
    }

    #[test]
    fn sew_gaps_applies_to_validate_phase() {
        let pass = SewGapsPass;
        let entry = DiagnosticEntry {
            severity: crate::step::diagnostics::Severity::Warning,
            phase: Phase::Validate,
            entity_id: None,
            message: "test".to_string(),
        };
        assert!(pass.applies_to(&entry));
    }

    #[test]
    fn sew_gaps_does_not_apply_to_parse_phase() {
        let pass = SewGapsPass;
        let entry = DiagnosticEntry {
            severity: crate::step::diagnostics::Severity::Warning,
            phase: Phase::Parse,
            entity_id: None,
            message: "test".to_string(),
        };
        assert!(!pass.applies_to(&entry));
    }

    #[test]
    fn reconcile_edge_curves_empty_brep() {
        let mut brep = BRep::new(tol());
        let count = reconcile_edge_curves(&mut brep);
        assert_eq!(count, 0);
    }

    #[test]
    fn sew_gaps_attempt_repair_no_changes() {
        let pass = SewGapsPass;
        let mut brep = BRep::new(tol());
        let _ = brep.vertices.insert(Vertex {
            position: Point3::new(0.0, 0.0, 0.0),
            outgoing: Default::default(),
        });
        let entry = DiagnosticEntry {
            severity: crate::step::diagnostics::Severity::Warning,
            phase: Phase::Validate,
            entity_id: None,
            message: "gap".to_string(),
        };
        let mut log = DiagnosticLog::new();
        let result = pass.attempt_repair(&mut brep, &entry, &mut log);
        assert_eq!(result.ok(), Some(false));
    }

    #[test]
    fn refit_line_curve_degenerate_returns_none() {
        use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
        let curve = CurveOnEdge {
            geometry: CurveGeometry::Line {
                origin: Point3::new(0.0, 0.0, 0.0),
                direction: Vec3::new(0.0, 0.0, 0.0),
            },
            t_start: 0.0,
            t_end: 1.0,
        };
        let result = refit_line_curve(
            &curve,
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            &tol(),
        );
        assert!(result.is_none());
    }

    #[test]
    fn refit_line_curve_non_line_returns_none() {
        use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
        let curve = CurveOnEdge {
            geometry: CurveGeometry::Circle {
                center: Point3::new(0.0, 0.0, 0.0),
                normal: Vec3::new(0.0, 0.0, 1.0),
                radius: 1.0,
            },
            t_start: 0.0,
            t_end: 1.0,
        };
        let result = refit_line_curve(
            &curve,
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            &tol(),
        );
        assert!(result.is_none());
    }

    #[test]
    fn refit_line_curve_successful() {
        use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
        // Line from (0,0,0) in direction (1,0,0)
        let curve = CurveOnEdge {
            geometry: CurveGeometry::Line {
                origin: Point3::new(0.0, 0.0, 0.0),
                direction: Vec3::new(1.0, 0.0, 0.0),
            },
            t_start: 0.0,
            t_end: 1.0,
        };
        // Refit to start=(2,0,0) end=(5,0,0) — these lie on the line
        let result = refit_line_curve(
            &curve,
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(5.0, 0.0, 0.0),
            &tol(),
        );
        assert!(result.is_some());
        let new_curve = result.unwrap();
        assert!((new_curve.t_start - 2.0).abs() < 1e-10);
        assert!((new_curve.t_end - 5.0).abs() < 1e-10);
    }

    #[test]
    fn sew_gaps_attempt_repair_with_merges() {
        let pass = SewGapsPass;
        let mut brep = BRep::new(tol());
        // Insert two coincident vertices
        let _ = brep.vertices.insert(Vertex {
            position: Point3::new(0.0, 0.0, 0.0),
            outgoing: Default::default(),
        });
        let _ = brep.vertices.insert(Vertex {
            position: Point3::new(0.0, 0.0, 0.0),
            outgoing: Default::default(),
        });
        let entry = DiagnosticEntry {
            severity: crate::step::diagnostics::Severity::Warning,
            phase: Phase::Validate,
            entity_id: None,
            message: "gap".to_string(),
        };
        let mut log = DiagnosticLog::new();
        let result = pass.attempt_repair(&mut brep, &entry, &mut log);
        assert_eq!(result.ok(), Some(true));
        assert_eq!(brep.vertices.len(), 1);
    }

    #[test]
    fn sew_gaps_applies_to_map_topology_phase() {
        let pass = SewGapsPass;
        let entry = DiagnosticEntry {
            severity: crate::step::diagnostics::Severity::Warning,
            phase: Phase::MapTopology,
            entity_id: None,
            message: "test".to_string(),
        };
        assert!(pass.applies_to(&entry));
    }

    #[test]
    fn merge_rewires_half_edge_origins() {
        use popkern_core::topology::{Edge, HalfEdge, HalfEdgeHandle};
        let mut brep = BRep::new(tol());
        let v0 = brep.vertices.insert(Vertex {
            position: Point3::new(0.0, 0.0, 0.0),
            outgoing: Default::default(),
        });
        // v1 is coincident with v0
        let v1 = brep.vertices.insert(Vertex {
            position: Point3::new(0.0, 0.0, 0.0),
            outgoing: Default::default(),
        });
        let v2 = brep.vertices.insert(Vertex {
            position: Point3::new(1.0, 0.0, 0.0),
            outgoing: Default::default(),
        });

        let edge_handle = brep.edges.insert(Edge {
            half_edge: HalfEdgeHandle::default(),
            curve: None,
        });
        // half-edge from v1 (will be remapped to v0)
        let he = brep.half_edges.insert(HalfEdge {
            origin: v1,
            twin: HalfEdgeHandle::default(),
            next: HalfEdgeHandle::default(),
            prev: HalfEdgeHandle::default(),
            edge: edge_handle,
            loop_: Default::default(),
            pcurve: None,
        });
        brep.half_edges[he].twin = he; // self-twin for simplicity
        brep.edges[edge_handle].half_edge = he;
        brep.vertices[v1].outgoing = he;

        let removed = merge_coincident_vertices(&mut brep, &tol());
        assert_eq!(removed, 1);
        // The half-edge origin should now point to v0
        assert_eq!(brep.half_edges[he].origin, v0);
        assert_eq!(brep.vertices.len(), 2);
    }
}
