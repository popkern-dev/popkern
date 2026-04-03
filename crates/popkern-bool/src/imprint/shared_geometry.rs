//! Shared geometry enforcement for watertight boolean output.
//!
//! The single most important implementation discipline for watertight output:
//! both sides of every intersection edge SHALL reference the exact same
//! geometric curve object — the same data, the same floating-point values.
//!
//! This module provides validation functions that check whether the shared
//! geometry principle is maintained in a BRep result, and reports any
//! violations.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::topology::BRep;

/// A detected gap between two incident faces at a shared edge.
#[derive(Debug, Clone)]
pub struct WatertightnessGap {
    /// The edge where the gap was detected.
    pub edge: popkern_core::topology::EdgeHandle,
    /// Maximum deviation between the two faces at the edge.
    pub max_deviation: f64,
    /// The 3D point where the maximum deviation occurs.
    pub worst_point: Point3,
}

/// Result of a watertightness validation.
#[derive(Debug, Clone, Default)]
pub struct WatertightnessReport {
    /// All detected gaps.
    pub gaps: Vec<WatertightnessGap>,
    /// Number of edges checked.
    pub edges_checked: usize,
}

impl WatertightnessReport {
    /// Returns true if no gaps were found.
    pub fn is_watertight(&self) -> bool {
        self.gaps.is_empty()
    }

    /// Maximum gap found across all edges.
    pub fn max_gap(&self) -> f64 {
        self.gaps
            .iter()
            .map(|g| g.max_deviation)
            .fold(0.0_f64, f64::max)
    }
}

/// Validate watertightness of a BRep by checking edge-surface consistency.
///
/// For each edge that has curve geometry and whose incident faces have
/// surface geometry, samples the edge at several points and evaluates
/// both incident face surfaces. The deviation between the two evaluations
/// is the gap.
///
/// # Arguments
///
/// * `brep` — the BRep to validate
/// * `tolerance` — gaps below this are acceptable
/// * `n_samples` — number of sample points per edge (default: 5)
pub fn validate_watertightness(
    brep: &BRep,
    tolerance: f64,
    n_samples: usize,
) -> WatertightnessReport {
    let mut report = WatertightnessReport::default();
    let n_samples = n_samples.max(2);

    for (eh, edge) in &brep.edges {
        let curve_on_edge = match &edge.curve {
            Some(c) => c,
            None => continue,
        };

        // Get the two incident faces via half-edges
        let he = match brep.half_edges.get(edge.half_edge) {
            Some(h) => h,
            None => continue,
        };
        let twin = match brep.half_edges.get(he.twin) {
            Some(t) => t,
            None => continue,
        };

        let loop_a = match brep.loops.get(he.loop_) {
            Some(l) => l,
            None => continue,
        };
        let loop_b = match brep.loops.get(twin.loop_) {
            Some(l) => l,
            None => continue,
        };
        let face_a = loop_a.face;
        let face_b = loop_b.face;

        let surf_a = match brep.faces.get(face_a).and_then(|f| f.surface.as_ref()) {
            Some(s) => &s.geometry,
            None => continue,
        };
        let surf_b = match brep.faces.get(face_b).and_then(|f| f.surface.as_ref()) {
            Some(s) => &s.geometry,
            None => continue,
        };

        report.edges_checked += 1;

        let t_start = curve_on_edge.t_start;
        let t_end = curve_on_edge.t_end;
        let mut max_dev = 0.0_f64;
        let mut worst_pt = Point3::new(0.0, 0.0, 0.0);

        for i in 0..n_samples {
            let t = t_start + (t_end - t_start) * (i as f64) / ((n_samples - 1) as f64);
            let edge_pt = curve_on_edge.geometry.point_at(t);

            // Project edge point onto both surfaces
            let proj_a = surf_a.closest_point(&edge_pt);
            let proj_b = surf_b.closest_point(&edge_pt);

            if let (Some((ua, va)), Some((ub, vb))) = (proj_a, proj_b) {
                let pt_a = surf_a.point_at(ua, va);
                let pt_b = surf_b.point_at(ub, vb);
                let dev = pt_a.distance_to(&pt_b);
                if dev > max_dev {
                    max_dev = dev;
                    worst_pt = edge_pt;
                }
            }
        }

        if max_dev > tolerance {
            report.gaps.push(WatertightnessGap {
                edge: eh,
                max_deviation: max_dev,
                worst_point: worst_pt,
            });
        }
    }

    report
}

/// Check if a single edge has consistent geometry with its incident faces.
///
/// Returns the maximum deviation found along the edge, or `None` if the
/// edge or faces lack geometry.
pub fn edge_face_deviation(
    brep: &BRep,
    edge: popkern_core::topology::EdgeHandle,
    n_samples: usize,
) -> Option<f64> {
    let edge_data = brep.edges.get(edge)?;
    let curve = edge_data.curve.as_ref()?;

    let he = brep.half_edges.get(edge_data.half_edge)?;
    let twin = brep.half_edges.get(he.twin)?;

    let face_a = brep.loops.get(he.loop_)?.face;
    let face_b = brep.loops.get(twin.loop_)?.face;

    let surf_a = brep.faces.get(face_a)?.surface.as_ref()?;
    let surf_b = brep.faces.get(face_b)?.surface.as_ref()?;

    let t_start = curve.t_start;
    let t_end = curve.t_end;
    let mut max_dev = 0.0_f64;

    for i in 0..n_samples.max(2) {
        let t = t_start + (t_end - t_start) * (i as f64) / ((n_samples.max(2) - 1) as f64);
        let edge_pt = curve.geometry.point_at(t);

        let dist_a = surf_a
            .geometry
            .closest_point(&edge_pt)
            .map(|(u, v)| surf_a.geometry.point_at(u, v).distance_to(&edge_pt));
        let dist_b = surf_b
            .geometry
            .closest_point(&edge_pt)
            .map(|(u, v)| surf_b.geometry.point_at(u, v).distance_to(&edge_pt));

        if let (Some(da), Some(db)) = (dist_a, dist_b) {
            max_dev = max_dev.max(da).max(db);
        }
    }

    Some(max_dev)
}
