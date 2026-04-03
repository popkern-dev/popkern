//! Snap-to-existing logic for the imprinting phase.
//!
//! When an intersection curve endpoint lies within the snapping tolerance of
//! an existing vertex or edge, the endpoint is merged to the existing entity
//! rather than creating a new vertex. This prevents micro-edges and is the
//! foundation of the watertightness guarantee.
//!
//! Snapping can cascade: merging vertex A to B may bring another vertex C
//! into snapping range. The [`snap_points`] function uses a worklist to
//! handle cascading with a bounded iteration count.

use popkern_core::geometry::Point3;

/// Result of a snap decision.
#[derive(Debug, Clone)]
pub struct SnapDecision {
    /// Index of the point that was snapped.
    pub source_index: usize,
    /// Index of the target it was snapped to.
    pub target_index: usize,
    /// Distance that was snapped (before merge).
    pub distance: f64,
}

/// Log of all snap decisions made during an imprinting pass.
#[derive(Debug, Clone, Default)]
pub struct SnapLog {
    /// All snap decisions in order.
    pub decisions: Vec<SnapDecision>,
}

impl SnapLog {
    /// Number of snaps performed.
    pub fn count(&self) -> usize {
        self.decisions.len()
    }
}

/// Snap a set of points to each other using a worklist approach.
///
/// Any two points within `snap_tolerance` are merged: the second point
/// is moved to the first point's position. Cascading is handled by
/// re-checking after each merge, up to `max_iterations`.
///
/// Returns the deduplicated point set and the snap log.
///
/// # Arguments
///
/// * `points` — mutable slice of points to snap
/// * `snap_tolerance` — maximum distance for snapping
/// * `max_iterations` — bound on cascading iterations (typically 10)
pub fn snap_points(points: &mut [Point3], snap_tolerance: f64, max_iterations: usize) -> SnapLog {
    let mut log = SnapLog::default();
    let snap_tol_sq = snap_tolerance * snap_tolerance;

    for _iteration in 0..max_iterations {
        let mut any_snapped = false;

        for i in 0..points.len() {
            for j in (i + 1)..points.len() {
                let dx = points[i].x - points[j].x;
                let dy = points[i].y - points[j].y;
                let dz = points[i].z - points[j].z;
                let dist_sq = dx * dx + dy * dy + dz * dz;

                if dist_sq < snap_tol_sq && dist_sq > 0.0 {
                    let dist = dist_sq.sqrt();
                    // Snap j to i (keep i's position)
                    points[j] = points[i];
                    log.decisions.push(SnapDecision {
                        source_index: j,
                        target_index: i,
                        distance: dist,
                    });
                    any_snapped = true;
                }
            }
        }

        if !any_snapped {
            break;
        }
    }

    log
}

/// Snap a single point to the nearest existing vertex if within tolerance.
///
/// Returns `Some(index)` if the point was snapped to vertex at that index,
/// or `None` if no vertex is close enough.
pub fn snap_to_nearest(
    point: &Point3,
    existing_vertices: &[Point3],
    snap_tolerance: f64,
) -> Option<usize> {
    let snap_tol_sq = snap_tolerance * snap_tolerance;
    let mut best_idx = None;
    let mut best_dist_sq = snap_tol_sq;

    for (i, v) in existing_vertices.iter().enumerate() {
        let dx = point.x - v.x;
        let dy = point.y - v.y;
        let dz = point.z - v.z;
        let dist_sq = dx * dx + dy * dy + dz * dz;

        if dist_sq < best_dist_sq {
            best_dist_sq = dist_sq;
            best_idx = Some(i);
        }
    }

    best_idx
}

/// Check if a proposed edge (between two points) would be a micro-edge.
///
/// Returns `true` if the edge length is below `micro_edge_tolerance`.
pub fn is_micro_edge(p1: &Point3, p2: &Point3, micro_edge_tolerance: f64) -> bool {
    p1.distance_to(p2) < micro_edge_tolerance
}
