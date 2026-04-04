//! Edge validity checks per spec §7.5.
//!
//! Before using any edge in the pipeline, validate:
//! - Parameter range: t_start != t_end
//! - Geometric length: arc length > micro_edge threshold
//! - Endpoint consistency: curve endpoints match vertex positions

use popkern_core::geometry::curve::CurveEval;
use popkern_core::topology::{BRep, EdgeHandle};

/// Result of edge validation.
#[derive(Debug, Clone)]
pub struct EdgeValidityReport {
    /// Edges with collapsed parameter range.
    pub collapsed_range: Vec<EdgeHandle>,
    /// Edges shorter than the threshold.
    pub too_short: Vec<EdgeHandle>,
    /// Edges where endpoints don't match vertices.
    pub endpoint_mismatch: Vec<(EdgeHandle, f64)>,
    /// Total edges checked.
    pub checked: usize,
}

impl EdgeValidityReport {
    /// True if all edges are valid.
    pub fn is_valid(&self) -> bool {
        self.collapsed_range.is_empty()
            && self.too_short.is_empty()
            && self.endpoint_mismatch.is_empty()
    }
}

/// Validate all edges in a BRep.
pub fn validate_edges(brep: &BRep, micro_tol: f64, endpoint_tol: f64) -> EdgeValidityReport {
    let mut report = EdgeValidityReport {
        collapsed_range: Vec::new(),
        too_short: Vec::new(),
        endpoint_mismatch: Vec::new(),
        checked: 0,
    };

    for (eh, edge) in &brep.edges {
        report.checked += 1;

        let curve = match &edge.curve {
            Some(c) => c,
            None => continue,
        };

        // 1. Parameter range: t_start != t_end
        if (curve.t_start - curve.t_end).abs() < 1e-15 {
            report.collapsed_range.push(eh);
            continue;
        }

        // 2. Geometric length
        let he = match brep.half_edges.get(edge.half_edge) {
            Some(h) => h,
            None => continue,
        };
        let twin = match brep.half_edges.get(he.twin) {
            Some(t) => t,
            None => continue,
        };

        let p_start = brep.vertices[he.origin].position;
        let p_end = brep.vertices[twin.origin].position;
        let chord_len = p_start.distance_to(&p_end);

        if chord_len < micro_tol {
            report.too_short.push(eh);
        }

        // 3. Endpoint consistency
        let curve_start = curve.geometry.point_at(curve.t_start);
        let curve_end = curve.geometry.point_at(curve.t_end);

        let start_dev = curve_start.distance_to(&p_start);
        let end_dev = curve_end.distance_to(&p_end);
        let max_dev = start_dev.max(end_dev);

        if max_dev > endpoint_tol {
            report.endpoint_mismatch.push((eh, max_dev));
        }
    }

    report
}

/// Validate intersection curves before use.
///
/// Checks each curve in a set of intersection records for:
/// - Minimum segment length
/// - Curve doesn't oscillate wildly (max deviation from chord)
pub fn validate_intersection_curves(
    curves: &[popkern_core::geometry::curve_on_surface::CurveOnSurface],
    min_length: f64,
) -> Vec<usize> {
    let mut invalid = Vec::new();

    for (i, cos) in curves.iter().enumerate() {
        let (t0, t1) = cos.curve_3d.parameter_range();
        let p0 = cos.curve_3d.point_at(t0);
        let p1 = cos.curve_3d.point_at(t1);
        let chord = p0.distance_to(&p1);

        if chord < min_length {
            invalid.push(i);
            continue;
        }

        // Check for oscillation: sample midpoint, verify it's near the chord
        let t_mid = (t0 + t1) * 0.5;
        let p_mid = cos.curve_3d.point_at(t_mid);
        let chord_mid = popkern_core::geometry::Point3::new(
            (p0.x + p1.x) * 0.5,
            (p0.y + p1.y) * 0.5,
            (p0.z + p1.z) * 0.5,
        );
        let deviation = p_mid.distance_to(&chord_mid);

        // If the midpoint deviates more than 10x the chord length, it's wild
        if deviation > chord * 10.0 {
            invalid.push(i);
        }
    }

    invalid
}
