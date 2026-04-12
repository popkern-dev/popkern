//! Quality scoring report for boolean results.
//!
//! Implements §11.4 of the boolean requirements: aggregates the
//! per-aspect detector outputs (sliver counts, micro-edge counts,
//! watertightness gaps, tolerance maxima) into a single
//! [`QualityReport`] and produces a normalized 0–100 [`QualityReport::score`]
//! against configurable thresholds.
//!
//! The score is intended for regression detection and CI dashboards,
//! not as a hard pass/fail gate. A "perfect" boolean result on
//! well-conditioned input scores 100; each detected issue subtracts a
//! configurable penalty.
//!
//! Quality reporting is read-only — it never mutates the BRep. Callers
//! that want repair should drive that from the report's individual
//! sub-counts (e.g. "if `micro_edges > 0`, run the collapse pass").

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::clip::face_vertex_positions;
use crate::imprint::micro_edge::find_micro_edges;
use crate::imprint::shared_geometry::validate_watertightness;
use crate::imprint::sliver::{SliverThresholds, analyze_face_polygon};

/// Aggregated quality metrics for a boolean result (§11.4).
///
/// All counts are non-negative; all "min" fields default to `f64::INFINITY`
/// when no sample is available, and "max" fields default to `0.0`.
#[derive(Debug, Clone)]
pub struct QualityReport {
    /// Number of faces classified as slivers by `analyze_face_polygon`.
    pub sliver_count: usize,
    /// Number of edges shorter than `tolerance.linear()`.
    pub micro_edge_count: usize,
    /// Number of watertightness gaps reported by `validate_watertightness`.
    pub watertightness_gap_count: usize,
    /// Maximum watertightness gap measurement (0 when watertight).
    pub max_watertightness_gap: f64,
    /// Maximum entity tolerance assigned anywhere in the result. The
    /// requirement (§11.3) is that this should not exceed the configured
    /// boolean tolerance; this field exposes the achieved value.
    pub max_entity_tolerance: f64,
    /// Length of the shortest edge in the result. `f64::INFINITY` when
    /// the result has no edges.
    pub min_edge_length: f64,
    /// Approximate area of the smallest face. `f64::INFINITY` when the
    /// result has no faces.
    pub min_face_area: f64,
    /// Worst aspect ratio (smallest value across all faces).
    /// 1.0 = perfectly square; close to 0 = highly elongated.
    pub worst_aspect_ratio: f64,
    /// Total face count.
    pub face_count: usize,
    /// Total edge count.
    pub edge_count: usize,
    /// Total vertex count.
    pub vertex_count: usize,
    /// Composite quality score in 0–100. 100 = clean result.
    /// Computed by [`QualityReport::compute_score`] from the configured
    /// thresholds; included on the report for convenience.
    pub score: u8,
}

/// Configurable thresholds and penalty weights for [`QualityReport::score`].
///
/// Defaults match the §11.4 acceptance criteria: any sliver, any
/// micro-edge, or any watertightness gap is treated as a quality issue
/// worth deducting points for. Tune the per-issue penalty when
/// integrating with a CI dashboard.
#[derive(Debug, Clone, Copy)]
pub struct QualityThresholds {
    /// Subtract this many score points per sliver face.
    pub sliver_penalty: u8,
    /// Subtract this many score points per micro-edge.
    pub micro_edge_penalty: u8,
    /// Subtract this many score points per watertightness gap.
    pub gap_penalty: u8,
    /// Subtract this many score points if the maximum entity
    /// tolerance exceeds the configured boolean tolerance.
    pub tolerance_violation_penalty: u8,
}

impl Default for QualityThresholds {
    fn default() -> Self {
        Self {
            sliver_penalty: 5,
            micro_edge_penalty: 3,
            gap_penalty: 10,
            tolerance_violation_penalty: 15,
        }
    }
}

/// Compute a [`QualityReport`] for a boolean result BRep.
///
/// All sub-detectors are read-only; the BRep is unchanged on return.
///
/// # Arguments
///
/// * `result` — the boolean result to score.
/// * `tolerance` — the boolean tolerance the result was produced under.
/// * `thresholds` — score weighting; pass `QualityThresholds::default()`
///   for the standard penalties.
pub fn quality_report(
    result: &BRep,
    tolerance: &Tolerance,
    thresholds: QualityThresholds,
) -> QualityReport {
    let sliver_thresholds = SliverThresholds::from_tolerance(tolerance.linear() * 5.0);

    let mut sliver_count = 0;
    let mut min_face_area = f64::INFINITY;
    let mut worst_aspect_ratio = 1.0_f64;
    for fh in result.faces.keys() {
        let positions = face_vertex_positions(result, fh);
        if positions.len() < 3 {
            continue;
        }
        let report = analyze_face_polygon(&positions, &sliver_thresholds);
        if report.is_sliver {
            sliver_count += 1;
        }
        if report.area < min_face_area {
            min_face_area = report.area;
        }
        if report.aspect_ratio < worst_aspect_ratio {
            worst_aspect_ratio = report.aspect_ratio;
        }
    }

    let micro_edges = find_micro_edges(result, tolerance.linear());
    let micro_edge_count = micro_edges.len();

    // Minimum edge length: walk every edge once, evaluate the
    // curve geometry's parameter range and take the chord length
    // between endpoints. We treat edges with no curve as having
    // length 0 (they shouldn't exist in a valid result, but the
    // detector must not panic on them).
    use popkern_core::geometry::curve::CurveEval;
    let mut min_edge_length = f64::INFINITY;
    for (_, edge) in &result.edges {
        let Some(curve_on_edge) = &edge.curve else {
            min_edge_length = 0.0;
            continue;
        };
        let p0 = curve_on_edge.geometry.point_at(curve_on_edge.t_start);
        let p1 = curve_on_edge.geometry.point_at(curve_on_edge.t_end);
        let len = p0.distance_to(&p1);
        if len < min_edge_length {
            min_edge_length = len;
        }
    }

    let wt = validate_watertightness(result, tolerance.linear() * 10.0, 3);
    let watertightness_gap_count = wt.gaps.len();
    let max_watertightness_gap = wt.max_gap();

    // Maximum entity tolerance: scan vertices, edges, faces and
    // pick the largest explicit tolerance value. Entities with
    // `None` inherit from `BRep::tolerance` and are excluded.
    let mut max_entity_tolerance = 0.0_f64;
    for (_, v) in &result.vertices {
        if let Some(t) = v.tolerance
            && t > max_entity_tolerance
        {
            max_entity_tolerance = t;
        }
    }
    for (_, e) in &result.edges {
        if let Some(t) = e.tolerance
            && t > max_entity_tolerance
        {
            max_entity_tolerance = t;
        }
    }
    for (_, f) in &result.faces {
        if let Some(t) = f.tolerance
            && t > max_entity_tolerance
        {
            max_entity_tolerance = t;
        }
    }

    let mut report = QualityReport {
        sliver_count,
        micro_edge_count,
        watertightness_gap_count,
        max_watertightness_gap,
        max_entity_tolerance,
        min_edge_length,
        min_face_area,
        worst_aspect_ratio,
        face_count: result.faces.len(),
        edge_count: result.edges.len(),
        vertex_count: result.vertices.len(),
        score: 100,
    };
    report.score = report.compute_score(tolerance, thresholds);
    report
}

impl QualityReport {
    /// Compute a 0–100 score from the report and configured thresholds.
    ///
    /// Subtracts a per-issue penalty for each detected defect. The
    /// score saturates at 0 (it never wraps around).
    pub fn compute_score(&self, tolerance: &Tolerance, thresholds: QualityThresholds) -> u8 {
        let mut deductions: u32 = 0;
        deductions += self.sliver_count as u32 * thresholds.sliver_penalty as u32;
        deductions += self.micro_edge_count as u32 * thresholds.micro_edge_penalty as u32;
        deductions += self.watertightness_gap_count as u32 * thresholds.gap_penalty as u32;
        if self.max_entity_tolerance > tolerance.linear() {
            deductions += thresholds.tolerance_violation_penalty as u32;
        }
        if deductions >= 100 {
            0
        } else {
            (100 - deductions) as u8
        }
    }

    /// Returns true if the result is "clean" — no slivers, no micro-edges,
    /// no watertightness gaps. Callers can short-circuit on this before
    /// inspecting the more granular fields.
    pub fn is_clean(&self) -> bool {
        self.sliver_count == 0 && self.micro_edge_count == 0 && self.watertightness_gap_count == 0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::Transform3;

    fn tol() -> Tolerance {
        Tolerance::default()
    }

    #[test]
    fn clean_cube_scores_high() {
        let cube = BRep::unit_cube(tol()).unwrap();
        let report = quality_report(&cube, &tol(), QualityThresholds::default());
        assert!(report.is_clean(), "unit cube should be clean: {report:?}");
        assert_eq!(report.score, 100, "perfect cube should score 100");
        assert_eq!(report.face_count, 6);
        assert_eq!(report.edge_count, 12);
        assert_eq!(report.vertex_count, 8);
    }

    #[test]
    fn boolean_result_quality_scoring() {
        let a = BRep::unit_cube(tol()).unwrap();
        let b = BRep::unit_cube(tol())
            .unwrap()
            .transform(&Transform3::translation(0.5, 0.3, 0.2))
            .unwrap();
        if let Ok(result) = crate::union(a, b, tol(), None) {
            let report = quality_report(&result, &tol(), QualityThresholds::default());
            // Score is in valid range
            assert!(report.score <= 100);
            assert!(report.face_count > 0);
            // Boolean results may have minor issues; score should reflect them.
        }
    }

    #[test]
    fn compute_score_saturates_at_zero() {
        let report = QualityReport {
            sliver_count: 100,
            micro_edge_count: 0,
            watertightness_gap_count: 0,
            max_watertightness_gap: 0.0,
            max_entity_tolerance: 0.0,
            min_edge_length: 1.0,
            min_face_area: 1.0,
            worst_aspect_ratio: 1.0,
            face_count: 100,
            edge_count: 0,
            vertex_count: 0,
            score: 0,
        };
        // 100 slivers × 5 penalty = 500 deductions, saturates at 0
        let score = report.compute_score(&tol(), QualityThresholds::default());
        assert_eq!(score, 0);
    }

    #[test]
    fn tolerance_violation_penalty_applied() {
        let mut cube = BRep::unit_cube(tol()).unwrap();
        // Force an entity tolerance much higher than the boolean tolerance
        let vh = cube.vertices.keys().next().unwrap();
        cube.vertices[vh].tolerance = Some(tol().linear() * 100.0);
        let report = quality_report(&cube, &tol(), QualityThresholds::default());
        // Score should be deducted for the tolerance violation.
        let default_thresholds = QualityThresholds::default();
        assert!(
            report.score < 100,
            "tolerance violation should reduce score from 100, got {}",
            report.score
        );
        assert_eq!(
            report.score as i32,
            100_i32 - default_thresholds.tolerance_violation_penalty as i32
        );
    }
}
