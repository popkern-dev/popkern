//! Degenerate entity detection (Section 5.3).
//!
//! Detects small edges, small faces, and entities with missing geometry.
//! Small entities often result from import operations or boolean operations
//! near tangent points, and may cause downstream failures in meshing or
//! further boolean operations.

use crate::geometry::curve::CurveEval;
use crate::healing::diagnostics::{
    DefectCategory, DefectRecord, EntityRef, HealingReport, Remediation, RemediationAction,
    Severity,
};
use crate::tolerance::Tolerance;
use crate::topology::BRep;

/// Classification of small face topology.
///
/// Used in defect descriptions to help the user understand the nature
/// of the degenerate face. The classification is based on the face's
/// vertex configuration and aspect ratio.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SmallFaceKind {
    /// A face collapsed to a near-point (all vertices very close together).
    Pin,
    /// A very thin, elongated face.
    Sliver,
    /// A face that bridges two regions with near-zero width.
    Bridge,
    /// A narrow slot-shaped face.
    Slot,
    /// A generally small face that doesn't fit other categories.
    Micro,
}

impl std::fmt::Display for SmallFaceKind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SmallFaceKind::Pin => write!(f, "pin"),
            SmallFaceKind::Sliver => write!(f, "sliver"),
            SmallFaceKind::Bridge => write!(f, "bridge"),
            SmallFaceKind::Slot => write!(f, "slot"),
            SmallFaceKind::Micro => write!(f, "micro"),
        }
    }
}

/// Run all degenerate entity checks.
pub fn check_all(brep: &BRep, tolerance: &Tolerance, report: &mut HealingReport) {
    check_small_edges(brep, tolerance.linear(), report);
    check_missing_edge_geometry(brep, report);
    check_missing_face_geometry(brep, report);
    check_small_faces(brep, tolerance, report);
    check_zero_length_singularity_edges(brep, tolerance, report);
}

/// Check for edges shorter than the minimum length.
///
/// For edges with curve geometry, computes the arc length via
/// `curve.length_between(t_start, t_end)`. For edges without curve
/// geometry, falls back to the Euclidean distance between the origin
/// vertices of the two half-edges.
pub fn check_small_edges(brep: &BRep, min_length: f64, report: &mut HealingReport) {
    for (edge_handle, edge_data) in &brep.edges {
        let length = if let Some(ref curve_on_edge) = edge_data.curve {
            curve_on_edge
                .geometry
                .length_between(curve_on_edge.t_start, curve_on_edge.t_end)
        } else {
            // No curve geometry: measure vertex-to-vertex distance.
            let Some(he) = brep.half_edges.get(edge_data.half_edge) else {
                continue;
            };
            let Some(twin) = brep.half_edges.get(he.twin) else {
                continue;
            };
            let Some(v_start) = brep.vertices.get(he.origin) else {
                continue;
            };
            let Some(v_end) = brep.vertices.get(twin.origin) else {
                continue;
            };
            v_start.position.distance_to(&v_end.position)
        };

        if length < min_length {
            report.push(DefectRecord {
                category: DefectCategory::Degenerate,
                severity: Severity::Warning,
                description: format!(
                    "Edge length {:.6} is below minimum {:.6}",
                    length, min_length
                ),
                entities: vec![EntityRef::Edge(edge_handle)],
                location: None,
                measured_value: Some(length),
                threshold: Some(min_length),
                remediation: Some(Remediation {
                    action: RemediationAction::RemoveEntity,
                    parameter: None,
                    value: None,
                    rationale: "Collapse small edge by merging its endpoints".into(),
                }),
            });
        }
    }
}

/// Check for edges with missing curve geometry.
///
/// Edges without curve geometry may be intentional (topology-only BReps)
/// or may indicate incomplete import. Recorded as Info severity.
pub fn check_missing_edge_geometry(brep: &BRep, report: &mut HealingReport) {
    for (edge_handle, edge_data) in &brep.edges {
        if edge_data.curve.is_none() {
            report.push(DefectRecord {
                category: DefectCategory::Degenerate,
                severity: Severity::Info,
                description: "Edge has no curve geometry".into(),
                entities: vec![EntityRef::Edge(edge_handle)],
                location: None,
                measured_value: None,
                threshold: None,
                remediation: Some(Remediation {
                    action: RemediationAction::RecomputeEdgeCurve,
                    parameter: None,
                    value: None,
                    rationale: "Derive edge curve from adjacent face surfaces".into(),
                }),
            });
        }
    }
}

/// Check for faces with missing surface geometry.
///
/// Faces without surface geometry may be intentional (topology-only BReps)
/// or may indicate incomplete import. Recorded as Info severity.
pub fn check_missing_face_geometry(brep: &BRep, report: &mut HealingReport) {
    for (face_handle, face_data) in &brep.faces {
        if face_data.surface.is_none() {
            report.push(DefectRecord {
                category: DefectCategory::Degenerate,
                severity: Severity::Info,
                description: "Face has no surface geometry".into(),
                entities: vec![EntityRef::Face(face_handle)],
                location: None,
                measured_value: None,
                threshold: None,
                remediation: None,
            });
        }
    }
}

/// Count the number of edges in a face's outer loop.
///
/// Walks the outer loop's half-edges and counts them.
fn count_face_edges(brep: &BRep, face: crate::topology::FaceHandle) -> usize {
    let Some(face_data) = brep.faces.get(face) else {
        return 0;
    };
    let Some(loop_data) = brep.loops.get(face_data.outer_loop) else {
        return 0;
    };
    let start_he = loop_data.half_edge;
    if !brep.half_edges.contains_key(start_he) {
        return 0;
    }
    let mut count = 0usize;
    let mut current = start_he;
    loop {
        count += 1;
        let he = &brep.half_edges[current];
        current = he.next;
        if current == start_he {
            break;
        }
        // Guard against degenerate loops.
        if count > brep.half_edges.len() {
            break;
        }
    }
    count
}

/// Count distinct neighbor faces adjacent to a face via shared edges.
///
/// Walks the outer loop half-edges, finds each twin's face, and collects
/// unique faces (excluding the face itself).
fn count_distinct_neighbor_faces(brep: &BRep, face: crate::topology::FaceHandle) -> usize {
    use std::collections::HashSet;

    let Some(face_data) = brep.faces.get(face) else {
        return 0;
    };
    let Some(loop_data) = brep.loops.get(face_data.outer_loop) else {
        return 0;
    };
    let start_he = loop_data.half_edge;
    if !brep.half_edges.contains_key(start_he) {
        return 0;
    }

    let mut neighbors = HashSet::new();
    let mut current = start_he;
    let mut guard = 0usize;
    loop {
        let he = &brep.half_edges[current];
        if let Some(twin) = brep.half_edges.get(he.twin)
            && let Some(twin_loop) = brep.loops.get(twin.loop_)
            && twin_loop.face != face
        {
            neighbors.insert(twin_loop.face);
        }
        current = he.next;
        guard += 1;
        if current == start_he || guard > brep.half_edges.len() {
            break;
        }
    }
    neighbors.len()
}

/// Check whether a face is narrow between its two neighbors.
///
/// A bridge face is narrow in the direction perpendicular to its long
/// axis. We check that the bounding box smallest dimension is much less
/// than the largest (aspect ratio < 0.3).
fn is_narrow_between_neighbors(dims: [f64; 3]) -> bool {
    let largest = dims[2];
    if largest < 1e-15 {
        return false;
    }
    // A bridge is elongated: smallest dimension << largest.
    dims[0] / largest < 0.3
}

/// Classify a small face into a sub-category based on its shape.
///
/// Uses a combination of bounding box aspect ratios and topological
/// analysis of the face's outer loop:
/// - **Bridge**: 4 edges, 2 distinct neighbor faces, narrow shape
/// - **Slot**: 4 edges, 1 distinct neighbor face (both long edges on same face)
/// - **Pin**: one dimension is negligible relative to the largest (ratio < 0.01)
/// - **Sliver**: high aspect ratio (smallest/middle ratio < 0.1)
/// - **Micro**: small in all dimensions (roughly isotropic)
fn classify_small_face(brep: &BRep, face: crate::topology::FaceHandle) -> SmallFaceKind {
    let positions = brep.face_outer_loop_positions(face);
    if positions.len() < 2 {
        return SmallFaceKind::Pin;
    }

    let mut min_x = f64::INFINITY;
    let mut min_y = f64::INFINITY;
    let mut min_z = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut max_y = f64::NEG_INFINITY;
    let mut max_z = f64::NEG_INFINITY;

    for pos in &positions {
        min_x = min_x.min(pos.x);
        min_y = min_y.min(pos.y);
        min_z = min_z.min(pos.z);
        max_x = max_x.max(pos.x);
        max_y = max_y.max(pos.y);
        max_z = max_z.max(pos.z);
    }

    let mut dims = [max_x - min_x, max_y - min_y, max_z - min_z];
    dims.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let smallest = dims[0];
    let middle = dims[1];
    let largest = dims[2];

    if largest < 1e-15 {
        return SmallFaceKind::Pin;
    }

    // Topological analysis for Bridge and Slot detection.
    let edge_count = count_face_edges(brep, face);
    let neighbor_count = count_distinct_neighbor_faces(brep, face);

    // Bridge: a narrow face connecting two larger regions.
    // Exactly 4 edges with 2 distinct neighbors, and the face is narrow.
    if edge_count == 4 && neighbor_count == 2 && is_narrow_between_neighbors(dims) {
        return SmallFaceKind::Bridge;
    }

    // Slot: a narrow channel where both long edges are shared with the
    // same face (forming a slot/groove). 4 edges, only 1 neighbor.
    if edge_count == 4 && neighbor_count == 1 {
        return SmallFaceKind::Slot;
    }

    // Fall back to geometric classification.
    if smallest / largest < 0.01 {
        SmallFaceKind::Pin
    } else if smallest / middle < 0.1 {
        SmallFaceKind::Sliver
    } else {
        SmallFaceKind::Micro
    }
}

/// Check for faces whose characteristic size is below the tolerance threshold.
///
/// Computes the bounding box diagonal of the face's outer loop vertex
/// positions. If the diagonal is smaller than `tolerance.linear * 10`,
/// the face is flagged as a small face.
pub fn check_small_faces(brep: &BRep, tolerance: &Tolerance, report: &mut HealingReport) {
    let threshold = tolerance.linear() * 10.0;

    for (face_handle, _face_data) in &brep.faces {
        let positions = brep.face_outer_loop_positions(face_handle);
        if positions.is_empty() {
            continue;
        }

        // Compute axis-aligned bounding box of the outer loop vertices.
        let mut min_x = f64::INFINITY;
        let mut min_y = f64::INFINITY;
        let mut min_z = f64::INFINITY;
        let mut max_x = f64::NEG_INFINITY;
        let mut max_y = f64::NEG_INFINITY;
        let mut max_z = f64::NEG_INFINITY;

        for pos in &positions {
            min_x = min_x.min(pos.x);
            min_y = min_y.min(pos.y);
            min_z = min_z.min(pos.z);
            max_x = max_x.max(pos.x);
            max_y = max_y.max(pos.y);
            max_z = max_z.max(pos.z);
        }

        let dx = max_x - min_x;
        let dy = max_y - min_y;
        let dz = max_z - min_z;
        let diagonal = (dx * dx + dy * dy + dz * dz).sqrt();

        if diagonal < threshold {
            let kind = classify_small_face(brep, face_handle);
            report.push(DefectRecord {
                category: DefectCategory::Degenerate,
                severity: Severity::Warning,
                description: format!(
                    "Small face ({kind}): bounding box diagonal {:.6} is below threshold {:.6}",
                    diagonal, threshold
                ),
                entities: vec![EntityRef::Face(face_handle)],
                location: None,
                measured_value: Some(diagonal),
                threshold: Some(threshold),
                remediation: Some(Remediation {
                    action: RemediationAction::RemoveEntity,
                    parameter: None,
                    value: None,
                    rationale: "Remove or merge small face with an adjacent face".into(),
                }),
            });
        }
    }
}

/// Check for zero-length topological edges at surface singularities (Section 5.3.3).
///
/// Detects edges whose two endpoint vertices are at the same position
/// (within tolerance) and whose adjacent face has a surface with a
/// singularity at that point. These edges are common at sphere poles
/// and cone apices, and while they may be intentional, they can cause
/// problems for downstream operations like meshing.
pub fn check_zero_length_singularity_edges(
    brep: &BRep,
    tolerance: &Tolerance,
    report: &mut HealingReport,
) {
    for (edge_handle, edge_data) in &brep.edges {
        let Some(he) = brep.half_edges.get(edge_data.half_edge) else {
            continue;
        };
        let Some(twin) = brep.half_edges.get(he.twin) else {
            continue;
        };
        let Some(v_start) = brep.vertices.get(he.origin) else {
            continue;
        };
        let Some(v_end) = brep.vertices.get(twin.origin) else {
            continue;
        };

        let dist = v_start.position.distance_to(&v_end.position);
        if crate::adaptive_precision::adaptive_dist_sq_lt(
            (v_start.position.x, v_start.position.y, v_start.position.z),
            (v_end.position.x, v_end.position.y, v_end.position.z),
            tolerance.linear() * tolerance.linear(),
        ) {
            report.push(DefectRecord {
                category: DefectCategory::Degenerate,
                severity: Severity::Info,
                description: "Zero-length edge at surface singularity".into(),
                entities: vec![EntityRef::Edge(edge_handle)],
                location: None,
                measured_value: Some(dist),
                threshold: Some(tolerance.linear()),
                remediation: Some(Remediation {
                    action: RemediationAction::RemoveEntity,
                    parameter: None,
                    value: None,
                    rationale: "Remove degenerate zero-length edge".into(),
                }),
            });
        }
    }
}

#[cfg(test)]
#[path = "degenerate_tests.rs"]
mod tests;
