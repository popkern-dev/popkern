//! Watertight repair pass for boolean results.
//!
//! Detects and closes small gaps between adjacent faces along shared edges.
//! Uses vertex position averaging for gaps smaller than the repair tolerance.
//! Reports gaps that are too large to repair.

use popkern_core::topology::BRep;

/// Result of a repair pass.
#[derive(Debug, Clone, Default)]
pub struct RepairReport {
    /// Number of vertices adjusted.
    pub vertices_adjusted: usize,
    /// Maximum gap that was closed.
    pub max_gap_closed: f64,
    /// Gaps that were too large to repair.
    pub unrepairable_gaps: usize,
}

impl RepairReport {
    /// Returns true if all gaps were repaired (or none existed).
    pub fn fully_repaired(&self) -> bool {
        self.unrepairable_gaps == 0
    }
}

/// Attempt to close small gaps in a BRep by adjusting vertex positions.
///
/// For each edge, evaluates the edge curve at its endpoints and checks
/// whether the vertex positions match. If the deviation is within
/// `repair_tolerance`, the vertex is moved to the curve endpoint.
///
/// Returns a report of what was repaired.
pub fn repair_vertex_positions(brep: &mut BRep, repair_tolerance: f64) -> RepairReport {
    use popkern_core::geometry::curve::CurveEval;

    let mut report = RepairReport::default();

    let edge_handles: Vec<_> = brep.edges.keys().collect();

    for eh in edge_handles {
        let edge = &brep.edges[eh];
        let curve = match &edge.curve {
            Some(c) => c,
            None => continue,
        };

        let he = match brep.half_edges.get(edge.half_edge) {
            Some(h) => h,
            None => continue,
        };
        let twin = match brep.half_edges.get(he.twin) {
            Some(t) => t,
            None => continue,
        };

        let vh_start = he.origin;
        let vh_end = twin.origin;

        // Evaluate curve at start and end parameters
        let curve_start = curve.geometry.point_at(curve.t_start);
        let curve_end = curve.geometry.point_at(curve.t_end);

        // Check start vertex
        let v_start_pos = brep.vertices[vh_start].position;
        let start_dev = v_start_pos.distance_to(&curve_start);
        if start_dev > 1e-15 && start_dev < repair_tolerance {
            brep.vertices[vh_start].position = curve_start;
            report.vertices_adjusted += 1;
            report.max_gap_closed = report.max_gap_closed.max(start_dev);
        } else if start_dev >= repair_tolerance {
            report.unrepairable_gaps += 1;
        }

        // Check end vertex
        let v_end_pos = brep.vertices[vh_end].position;
        let end_dev = v_end_pos.distance_to(&curve_end);
        if end_dev > 1e-15 && end_dev < repair_tolerance {
            brep.vertices[vh_end].position = curve_end;
            report.vertices_adjusted += 1;
            report.max_gap_closed = report.max_gap_closed.max(end_dev);
        } else if end_dev >= repair_tolerance {
            report.unrepairable_gaps += 1;
        }
    }

    report
}
