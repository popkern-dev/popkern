//! Topological defect detection (Section 5.1).
//!
//! Wraps the existing BRep validation methods and converts their
//! `ValidationError` results into structured `DefectRecord`s for
//! the healing report.

use crate::error::ValidationError;
use crate::geometry::curve::CurveEval;
use crate::healing::diagnostics::{
    DefectCategory, DefectRecord, EntityRef, HealingReport, Severity,
};
use crate::tolerance::Tolerance;
use crate::topology::BRep;

/// Run all topological defect checks.
pub fn check_all(brep: &BRep, tolerance: &Tolerance, report: &mut HealingReport) {
    check_boundary_closure(brep, tolerance, report);
    check_half_edge_pairing(brep, report);
    check_half_edge_cycles(brep, report);
    check_vertex_consistency(brep, report);
    check_manifold_edges(brep, report);
    check_face_orientation(brep, report);
    check_shell_closure(brep, report);
    check_boundary_simplicity(brep, tolerance, report);
    check_inner_loop_orientation(brep, report);
}

/// Check that every loop is geometrically closed within tolerance.
///
/// Calls `brep.validate_loop_closure()` and converts any
/// `ValidationError::OpenLoop` into a defect record.
pub fn check_boundary_closure(brep: &BRep, _tolerance: &Tolerance, report: &mut HealingReport) {
    if let Err(e) = brep.validate_loop_closure() {
        report.push(validation_error_to_defect(&e, Severity::Error));
    }
}

/// Check that every half-edge has a valid twin and twin-of-twin consistency.
///
/// Calls `brep.validate_half_edge_pairing()` and converts any
/// `ValidationError::BrokenTwin` into a defect record.
pub fn check_half_edge_pairing(brep: &BRep, report: &mut HealingReport) {
    if let Err(e) = brep.validate_half_edge_pairing() {
        report.push(validation_error_to_defect(&e, Severity::Error));
    }
}

/// Check that following `next` from any half-edge returns to the start.
///
/// Calls `brep.validate_half_edge_cycles()` and converts any
/// `ValidationError::BrokenCycle` into a defect record.
pub fn check_half_edge_cycles(brep: &BRep, report: &mut HealingReport) {
    if let Err(e) = brep.validate_half_edge_cycles() {
        report.push(validation_error_to_defect(&e, Severity::Error));
    }
}

/// Check that each vertex's outgoing half-edge originates from it.
///
/// Calls `brep.validate_vertex_consistency()` and converts any
/// `ValidationError::InconsistentVertex` into a defect record.
pub fn check_vertex_consistency(brep: &BRep, report: &mut HealingReport) {
    if let Err(e) = brep.validate_vertex_consistency() {
        report.push(validation_error_to_defect(&e, Severity::Error));
    }
}

/// Check that every edge is manifold (exactly 2 half-edges).
///
/// Calls `brep.validate_manifold_edges()` and converts any
/// `ValidationError::NonManifoldEdge` into a defect record.
pub fn check_manifold_edges(brep: &BRep, report: &mut HealingReport) {
    if let Err(e) = brep.validate_manifold_edges() {
        report.push(validation_error_to_defect(&e, Severity::Error));
    }
}

/// Check that face surface normals are not degenerate.
///
/// Calls `brep.validate_face_orientation()` and converts any
/// `ValidationError::InconsistentOrientation` into a defect record.
/// Uses Warning severity since orientation issues are often repairable.
pub fn check_face_orientation(brep: &BRep, report: &mut HealingReport) {
    if let Err(e) = brep.validate_face_orientation() {
        report.push(validation_error_to_defect(&e, Severity::Warning));
    }
}

/// Check that the shell is watertight (no gaps at shared edges).
///
/// Calls `brep.validate_watertight()` and converts any
/// `ValidationError::ShellGap` into a defect record.
pub fn check_shell_closure(brep: &BRep, report: &mut HealingReport) {
    if let Err(e) = brep.validate_watertight() {
        report.push(validation_error_to_defect(&e, Severity::Error));
    }
}

/// Check boundary simplicity -- no wire loop self-intersects (Section 5.1.2).
///
/// For each loop, samples edge curves at several points and checks whether
/// any pair of non-adjacent edge segments comes within tolerance of each other.
/// This is a simplified check using point sampling rather than exact
/// curve-curve intersection.
pub fn check_boundary_simplicity(brep: &BRep, tolerance: &Tolerance, report: &mut HealingReport) {
    const SAMPLES_PER_EDGE: usize = 4;

    for (loop_handle, loop_data) in &brep.loops {
        // Collect (edge_handle, sample_points, origin_vertex, end_vertex) per edge.
        let mut edge_samples: Vec<(
            crate::topology::EdgeHandle,
            Vec<crate::geometry::Point3>,
            crate::topology::VertexHandle,
            crate::topology::VertexHandle,
        )> = Vec::new();

        let start_he = loop_data.half_edge;
        let mut cursor = start_he;
        let bound = brep.half_edges.len().max(1);

        for _ in 0..=bound {
            let Some(he) = brep.half_edges.get(cursor) else {
                break;
            };

            let origin_v = he.origin;
            // End vertex is the origin of the next half-edge.
            let end_v = brep
                .half_edges
                .get(he.next)
                .map(|next_he| next_he.origin)
                .unwrap_or(origin_v);

            let edge_handle = he.edge;
            let mut points = Vec::new();

            if let Some(edge_data) = brep.edges.get(edge_handle)
                && let Some(ref curve_on_edge) = edge_data.curve
            {
                let t0 = curve_on_edge.t_start;
                let t1 = curve_on_edge.t_end;
                for k in 0..SAMPLES_PER_EDGE {
                    let frac = (k as f64 + 0.5) / SAMPLES_PER_EDGE as f64;
                    let t = t0 + (t1 - t0) * frac;
                    points.push(curve_on_edge.geometry.point_at(t));
                }
            }

            edge_samples.push((edge_handle, points, origin_v, end_v));

            cursor = he.next;
            if cursor == start_he {
                break;
            }
        }

        // Check all pairs of non-adjacent edges for proximity.
        let n = edge_samples.len();
        for i in 0..n {
            for j in (i + 2)..n {
                // Skip adjacent edges (they share a vertex).
                if i == 0 && j == n - 1 {
                    continue; // first and last edges are adjacent in a loop
                }
                let (_, ref pts_a, _, _) = edge_samples[i];
                let (_, ref pts_b, _, _) = edge_samples[j];
                let lin_sq = tolerance.linear() * tolerance.linear();
                for pa in pts_a {
                    for pb in pts_b {
                        if crate::adaptive_precision::adaptive_dist_sq_lt(
                            (pa.x, pa.y, pa.z),
                            (pb.x, pb.y, pb.z),
                            lin_sq,
                        ) {
                            report.push(DefectRecord {
                                category: DefectCategory::Topological,
                                severity: Severity::Error,
                                description: "Loop boundary self-intersection \
                                              detected between non-adjacent edges"
                                    .into(),
                                entities: vec![
                                    EntityRef::Loop(loop_handle),
                                    EntityRef::Edge(edge_samples[i].0),
                                    EntityRef::Edge(edge_samples[j].0),
                                ],
                                location: None,
                                measured_value: Some(pa.distance_to(pb)),
                                threshold: Some(tolerance.linear()),
                                remediation: None,
                            });
                            // One defect per loop is enough for this pair.
                            break;
                        }
                    }
                }
            }
        }
    }
}

/// Check that inner loops have opposite winding to outer loops (Section 5.1.3).
///
/// For each face with inner loops, verifies that the signed area of
/// inner loops has opposite sign to the outer loop, indicating correct
/// orientation (holes wind opposite to the boundary).
pub fn check_inner_loop_orientation(brep: &BRep, report: &mut HealingReport) {
    for (face_handle, face_data) in &brep.faces {
        if face_data.inner_loops.is_empty() {
            continue;
        }

        let outer_positions = loop_positions(brep, face_data.outer_loop);
        if outer_positions.len() < 3 {
            continue;
        }

        // Compute a reference normal from the face surface or from the
        // outer loop polygon cross products.
        let ref_normal = loop_polygon_normal(&outer_positions);
        if ref_normal.length() < 1e-15 {
            continue;
        }

        let outer_area = signed_area_projected(&outer_positions, &ref_normal);

        for &inner_loop in &face_data.inner_loops {
            let inner_positions = loop_positions(brep, inner_loop);
            if inner_positions.len() < 3 {
                continue;
            }

            let inner_area = signed_area_projected(&inner_positions, &ref_normal);

            // Inner loops should have opposite sign to the outer loop.
            if (outer_area > 0.0 && inner_area > 0.0) || (outer_area < 0.0 && inner_area < 0.0) {
                report.push(DefectRecord {
                    category: DefectCategory::Topological,
                    severity: Severity::Warning,
                    description: "Inner loop has same winding direction as \
                                  outer loop — expected opposite orientation"
                        .into(),
                    entities: vec![EntityRef::Face(face_handle), EntityRef::Loop(inner_loop)],
                    location: None,
                    measured_value: None,
                    threshold: None,
                    remediation: None,
                });
            }
        }
    }
}

/// Collect vertex positions from a loop by handle.
fn loop_positions(
    brep: &BRep,
    loop_handle: crate::topology::LoopHandle,
) -> Vec<crate::geometry::Point3> {
    let Some(loop_data) = brep.loops.get(loop_handle) else {
        return Vec::new();
    };
    let start_he = loop_data.half_edge;
    if !brep.half_edges.contains_key(start_he) {
        return Vec::new();
    }
    let mut positions = Vec::new();
    let mut cursor = start_he;
    let bound = brep.half_edges.len().max(1);
    for _ in 0..=bound {
        let Some(he) = brep.half_edges.get(cursor) else {
            break;
        };
        if let Some(v) = brep.vertices.get(he.origin) {
            positions.push(v.position);
        }
        cursor = he.next;
        if cursor == start_he {
            break;
        }
    }
    positions
}

/// Compute a polygon normal from loop vertex positions via Newell's method.
fn loop_polygon_normal(positions: &[crate::geometry::Point3]) -> crate::geometry::Vec3 {
    let mut nx = 0.0;
    let mut ny = 0.0;
    let mut nz = 0.0;
    let n = positions.len();
    for i in 0..n {
        let cur = &positions[i];
        let next = &positions[(i + 1) % n];
        nx += (cur.y - next.y) * (cur.z + next.z);
        ny += (cur.z - next.z) * (cur.x + next.x);
        nz += (cur.x - next.x) * (cur.y + next.y);
    }
    crate::geometry::Vec3::new(nx, ny, nz)
}

/// Compute the signed area of a polygon projected onto a reference normal.
fn signed_area_projected(
    positions: &[crate::geometry::Point3],
    normal: &crate::geometry::Vec3,
) -> f64 {
    let n_len = normal.length();
    if n_len < 1e-15 {
        return 0.0;
    }
    let unit_n = *normal * (1.0 / n_len);

    // Sum cross products of consecutive edge vectors, project onto normal.
    let origin = &positions[0];
    let mut area = 0.0;
    for i in 1..positions.len().saturating_sub(1) {
        let a = crate::geometry::Vec3::new(
            positions[i].x - origin.x,
            positions[i].y - origin.y,
            positions[i].z - origin.z,
        );
        let b = crate::geometry::Vec3::new(
            positions[i + 1].x - origin.x,
            positions[i + 1].y - origin.y,
            positions[i + 1].z - origin.z,
        );
        area += a.cross(&b).dot(&unit_n);
    }
    area * 0.5
}

/// Convert a `ValidationError` into a `DefectRecord`.
///
/// Extracts the affected entity handle from the error variant and
/// uses the error's `Display` string as the human-readable description.
fn validation_error_to_defect(error: &ValidationError, severity: Severity) -> DefectRecord {
    let entities = extract_entities(error);
    let (measured_value, threshold) = extract_measurements(error);

    DefectRecord {
        category: DefectCategory::Topological,
        severity,
        description: error.to_string(),
        entities,
        location: None,
        measured_value,
        threshold,
        remediation: None,
    }
}

/// Extract entity references from a `ValidationError` variant.
fn extract_entities(error: &ValidationError) -> Vec<EntityRef> {
    match error {
        ValidationError::EulerPoincare { .. } => Vec::new(),
        ValidationError::DanglingReference { .. } => Vec::new(),
        ValidationError::BrokenTwin { half_edge } => {
            vec![EntityRef::HalfEdge(*half_edge)]
        }
        ValidationError::BrokenCycle {
            half_edge, loop_, ..
        } => {
            vec![EntityRef::HalfEdge(*half_edge), EntityRef::Loop(*loop_)]
        }
        ValidationError::LoopJunctionMismatch {
            half_edge,
            next,
            loop_,
            ..
        } => vec![
            EntityRef::HalfEdge(*half_edge),
            EntityRef::HalfEdge(*next),
            EntityRef::Loop(*loop_),
        ],
        ValidationError::InconsistentVertex {
            vertex, half_edge, ..
        } => {
            vec![EntityRef::Vertex(*vertex), EntityRef::HalfEdge(*half_edge)]
        }
        ValidationError::NonManifoldEdge { edge } => {
            vec![EntityRef::Edge(*edge)]
        }
        ValidationError::LowValence { vertex, .. } => {
            vec![EntityRef::Vertex(*vertex)]
        }
        ValidationError::InconsistentOrientation { face } => {
            vec![EntityRef::Face(*face)]
        }
        ValidationError::OpenLoop { loop_, .. } => {
            vec![EntityRef::Loop(*loop_)]
        }
        ValidationError::CurveSurfaceDeviation { edge, .. } => {
            vec![EntityRef::Edge(*edge)]
        }
        ValidationError::CanonicalBackward { edge }
        | ValidationError::CurveWindowOutsideSpan { edge, .. } => {
            vec![EntityRef::Edge(*edge)]
        }
        ValidationError::ShellGap { edge, .. } => {
            vec![EntityRef::Edge(*edge)]
        }
        ValidationError::SelfIntersection { face_a, face_b } => {
            vec![EntityRef::FacePair(*face_a, *face_b)]
        }
        ValidationError::WindingDishonest { face, .. }
        | ValidationError::InnerLoopWinding { face, .. } => {
            vec![EntityRef::Face(*face)]
        }
        ValidationError::MissingPcurves { face, .. } => {
            vec![EntityRef::Face(*face)]
        }
        ValidationError::TwistedSurface { face } => {
            vec![EntityRef::Face(*face)]
        }
        ValidationError::LoopPeriodWinding { face, .. } => {
            vec![EntityRef::Face(*face)]
        }
        ValidationError::PcurveDeviation { half_edge, .. } => {
            vec![EntityRef::HalfEdge(*half_edge)]
        }
        ValidationError::Precision { .. } => Vec::new(),
    }
}

/// Extract measured value and threshold from measurement-carrying error variants.
fn extract_measurements(error: &ValidationError) -> (Option<f64>, Option<f64>) {
    match error {
        ValidationError::OpenLoop { gap, .. } => (Some(*gap), None),
        ValidationError::CurveSurfaceDeviation { max_deviation, .. } => {
            (Some(*max_deviation), None)
        }
        ValidationError::ShellGap { gap, .. } => (Some(*gap), None),
        _ => (None, None),
    }
}

#[cfg(test)]
#[path = "topological_tests.rs"]
mod tests;
