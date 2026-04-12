//! Continuity enforcement — adjust geometry to achieve target continuity.
//!
//! Supports G0 enforcement (close positional gaps), G1 enforcement
//! (tangent-plane continuity via NURBS control-point adjustment),
//! G2 enforcement (curvature continuity via second-row control-point
//! adjustment), and dispatches to [`super::g3_enforcement`] for G3
//! (curvature-rate continuity via third-row control-point adjustment).

use crate::error::CoreError;
use crate::geometry::Point3;
use crate::geometry::curve::CurveEval;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::healing::curvature::surface_curvature::surface_curvature_at;
use crate::healing::quality::ContinuityLevel;
use crate::tolerance::Tolerance;
use crate::topology::{BRep, EdgeHandle, FaceHandle};

/// Result of a continuity enforcement operation.
#[derive(Debug, Clone)]
pub struct EnforcementResult {
    /// Edge that was enforced.
    pub edge: EdgeHandle,
    /// Maximum displacement introduced.
    pub max_displacement: f64,
    /// Continuity level achieved after enforcement.
    pub achieved_level: ContinuityLevel,
}

/// Enforce G0 continuity at a shared edge by adjusting vertex positions.
///
/// Measures the gap between the edge curve endpoints and their topological
/// vertices. If the gap exceeds tolerance, moves vertices toward the curve
/// endpoints. Reports maximum displacement. Returns error if displacement
/// would exceed `max_displacement`.
///
/// # Errors
///
/// Returns [`CoreError::InvalidHandle`] if the edge is missing from the BRep.
/// Returns [`CoreError::Healing`] if the required displacement exceeds
/// `max_displacement`.
pub fn enforce_g0(
    brep: &mut BRep,
    edge: EdgeHandle,
    tolerance: &Tolerance,
    max_displacement: f64,
) -> Result<EnforcementResult, CoreError> {
    let edge_data = brep
        .edges
        .get(edge)
        .ok_or_else(|| CoreError::InvalidHandle(format!("edge {edge:?}")))?;

    // If no curve geometry, nothing to enforce — already trivially G0.
    let curve_on_edge = match &edge_data.curve {
        Some(c) => c.clone(),
        None => {
            return Ok(EnforcementResult {
                edge,
                max_displacement: 0.0,
                achieved_level: ContinuityLevel::G0,
            });
        }
    };

    // Get start vertex via the half-edge.
    let he_handle = edge_data.half_edge;
    let he = brep
        .half_edges
        .get(he_handle)
        .ok_or_else(|| CoreError::InvalidHandle("half-edge".to_string()))?;
    let start_vh = he.origin;
    let twin_handle = he.twin;

    // Get end vertex via the twin half-edge.
    let twin = brep
        .half_edges
        .get(twin_handle)
        .ok_or_else(|| CoreError::InvalidHandle("twin half-edge".to_string()))?;
    let end_vh = twin.origin;

    // Read current vertex positions.
    let start_pos = brep
        .vertices
        .get(start_vh)
        .ok_or_else(|| CoreError::InvalidHandle("start vertex".to_string()))?
        .position;
    let end_pos = brep
        .vertices
        .get(end_vh)
        .ok_or_else(|| CoreError::InvalidHandle("end vertex".to_string()))?
        .position;

    // Evaluate curve endpoints.
    let curve_start = curve_on_edge.geometry.point_at(curve_on_edge.t_start);
    let curve_end = curve_on_edge.geometry.point_at(curve_on_edge.t_end);

    let gap_start = start_pos.distance_to(&curve_start);
    let gap_end = end_pos.distance_to(&curve_end);

    if gap_start < tolerance.linear() && gap_end < tolerance.linear() {
        return Ok(EnforcementResult {
            edge,
            max_displacement: 0.0,
            achieved_level: ContinuityLevel::G0,
        });
    }

    let mut actual_max_displacement = 0.0_f64;

    let start_adjustment = if gap_start >= tolerance.linear() {
        let mid = midpoint(&start_pos, &curve_start);
        let disp = start_pos.distance_to(&mid);
        actual_max_displacement = actual_max_displacement.max(disp);
        Some((start_vh, mid, disp))
    } else {
        None
    };

    let end_adjustment = if gap_end >= tolerance.linear() {
        let mid = midpoint(&end_pos, &curve_end);
        let disp = end_pos.distance_to(&mid);
        actual_max_displacement = actual_max_displacement.max(disp);
        Some((end_vh, mid, disp))
    } else {
        None
    };

    if actual_max_displacement > max_displacement {
        return Err(CoreError::Healing(
            crate::error::HealingError::DisplacementExceeded {
                edge,
                displacement: actual_max_displacement,
                limit: max_displacement,
            },
        ));
    }

    if let Some((vh, new_pos, _)) = start_adjustment
        && let Some(v) = brep.vertices.get_mut(vh)
    {
        v.position = new_pos;
    }
    if let Some((vh, new_pos, _)) = end_adjustment
        && let Some(v) = brep.vertices.get_mut(vh)
    {
        v.position = new_pos;
    }

    Ok(EnforcementResult {
        edge,
        max_displacement: actual_max_displacement,
        achieved_level: ContinuityLevel::G0,
    })
}

/// Enforce G0 continuity across all edges in a BRep.
///
/// Iterates every edge and attempts [`enforce_g0`]. Edges that fail
/// (e.g. displacement exceeds limit or missing topology) are silently
/// skipped. Returns results only for edges where enforcement succeeded
/// and had a non-zero displacement.
pub fn enforce_g0_all(
    brep: &mut BRep,
    tolerance: &Tolerance,
    max_displacement: f64,
) -> Vec<EnforcementResult> {
    let edge_handles: Vec<EdgeHandle> = brep.edges.keys().collect();
    let mut results = Vec::new();
    for eh in edge_handles {
        if let Ok(result) = enforce_g0(brep, eh, tolerance, max_displacement)
            && result.max_displacement > 0.0
        {
            results.push(result);
        }
    }
    results
}

/// Enforce G1 (tangent-plane) continuity at a shared edge.
///
/// For each sample point along the shared boundary, evaluates normals on
/// both adjacent surfaces. If the angle between normals exceeds tolerance
/// and at least one surface is NURBS, adjusts the first row of control
/// points perpendicular to the boundary to align the normals. For purely
/// analytical surfaces, reports that G1 cannot be enforced without
/// changing surface type (returns with `achieved_level: G0`).
pub fn enforce_g1(
    brep: &mut BRep,
    edge: EdgeHandle,
    tolerance: &Tolerance,
    max_displacement: f64,
) -> Result<EnforcementResult, CoreError> {
    let (f1, f2, sa, sb) = get_adjacent_surfaces(brep, edge)?;

    let is_nurbs = |sg: &SurfaceGeometry| matches!(sg, SurfaceGeometry::Nurbs(_));
    if !is_nurbs(&sa.geometry) && !is_nurbs(&sb.geometry) {
        return Ok(EnforcementResult {
            edge,
            max_displacement: 0.0,
            achieved_level: ContinuityLevel::G0,
        });
    }

    let samples = sample_edge_uvs(brep, edge, 8);
    let max_angle = measure_normal_angle(&sa.geometry, &sb.geometry, &samples);

    if max_angle <= tolerance.angular() {
        return Ok(EnforcementResult {
            edge,
            max_displacement: 0.0,
            achieved_level: ContinuityLevel::G1,
        });
    }

    let disp = super::nurbs_adjust::adjust_for_g1(brep, f1, f2, &samples, max_displacement);

    // Re-evaluate after adjustment.
    let (_, _, sa2, sb2) = get_adjacent_surfaces(brep, edge)?;
    let post = measure_normal_angle(&sa2.geometry, &sb2.geometry, &samples);
    let achieved = if post <= tolerance.angular() {
        ContinuityLevel::G1
    } else {
        ContinuityLevel::G0
    };

    Ok(EnforcementResult {
        edge,
        max_displacement: disp,
        achieved_level: achieved,
    })
}

/// Enforce G2 (curvature) continuity at a shared edge.
///
/// Requires both adjacent surfaces to be NURBS with degree >= 3.
/// Adjusts the second row of control points perpendicular to the
/// boundary to match curvature across the junction. First ensures G1
/// is achieved; if G1 cannot be reached, returns at that level.
pub fn enforce_g2(
    brep: &mut BRep,
    edge: EdgeHandle,
    tolerance: &Tolerance,
    max_displacement: f64,
) -> Result<EnforcementResult, CoreError> {
    let g1 = enforce_g1(brep, edge, tolerance, max_displacement)?;
    if g1.achieved_level < ContinuityLevel::G1 {
        return Ok(g1);
    }
    let budget = max_displacement - g1.max_displacement;
    if budget <= 0.0 {
        return Ok(g1);
    }

    let (f1, f2, sa, sb) = get_adjacent_surfaces(brep, edge)?;
    let need_g2 = |sg: &SurfaceGeometry| matches!(sg, SurfaceGeometry::Nurbs(n) if n.degree_u >= 3 && n.degree_v >= 3);
    if !need_g2(&sa.geometry) || !need_g2(&sb.geometry) {
        return Ok(g1);
    }

    let samples = sample_edge_uvs(brep, edge, 8);
    let curv_tol = 1e-3_f64;
    let diff = measure_curvature_diff(&sa.geometry, &sb.geometry, &samples);
    if diff <= curv_tol {
        return Ok(EnforcementResult {
            edge,
            max_displacement: g1.max_displacement,
            achieved_level: ContinuityLevel::G2,
        });
    }

    let g2_disp = super::nurbs_adjust::adjust_for_g2(brep, f1, f2, &samples, budget);

    let (_, _, sa2, sb2) = get_adjacent_surfaces(brep, edge)?;
    let post = measure_curvature_diff(&sa2.geometry, &sb2.geometry, &samples);
    let achieved = if post <= curv_tol {
        ContinuityLevel::G2
    } else {
        ContinuityLevel::G1
    };

    Ok(EnforcementResult {
        edge,
        max_displacement: g1.max_displacement + g2_disp,
        achieved_level: achieved,
    })
}

/// Enforce continuity at all edges up to the specified level.
///
/// Iterates every edge and calls the appropriate enforcement function
/// based on `target_level`. Returns a result for every edge where
/// enforcement succeeded.
pub fn enforce_continuity_all(
    brep: &mut BRep,
    target_level: ContinuityLevel,
    tolerance: &Tolerance,
    max_displacement: f64,
) -> Vec<EnforcementResult> {
    let handles: Vec<EdgeHandle> = brep.edges.keys().collect();
    let mut results = Vec::new();
    for eh in handles {
        let r = match target_level {
            ContinuityLevel::G0 => enforce_g0(brep, eh, tolerance, max_displacement),
            ContinuityLevel::G1 => enforce_g1(brep, eh, tolerance, max_displacement),
            ContinuityLevel::G2 => enforce_g2(brep, eh, tolerance, max_displacement),
            ContinuityLevel::G3 => {
                super::g3_enforcement::enforce_g3(brep, eh, tolerance, max_displacement)
            }
        };
        if let Ok(result) = r {
            results.push(result);
        }
    }
    results
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// Retrieve the two adjacent face surfaces across an edge.
pub(super) fn get_adjacent_surfaces(
    brep: &BRep,
    edge: EdgeHandle,
) -> Result<
    (
        FaceHandle,
        FaceHandle,
        crate::geometry::SurfaceOnFace,
        crate::geometry::SurfaceOnFace,
    ),
    CoreError,
> {
    let ed = brep
        .edges
        .get(edge)
        .ok_or_else(|| CoreError::InvalidHandle(format!("edge {edge:?}")))?;
    let he1 = brep
        .half_edges
        .get(ed.half_edge)
        .ok_or_else(|| CoreError::InvalidHandle("half-edge".into()))?;
    let he2 = brep
        .half_edges
        .get(he1.twin)
        .ok_or_else(|| CoreError::InvalidHandle("twin half-edge".into()))?;
    let l1 = brep
        .loops
        .get(he1.loop_)
        .ok_or_else(|| CoreError::InvalidHandle("loop".into()))?;
    let l2 = brep
        .loops
        .get(he2.loop_)
        .ok_or_else(|| CoreError::InvalidHandle("loop".into()))?;
    let f1 = brep
        .faces
        .get(l1.face)
        .ok_or_else(|| CoreError::InvalidHandle("face".into()))?;
    let f2 = brep
        .faces
        .get(l2.face)
        .ok_or_else(|| CoreError::InvalidHandle("face".into()))?;
    let s1 = f1
        .surface
        .clone()
        .ok_or_else(|| CoreError::PreconditionFailed("face missing surface".into()))?;
    let s2 = f2
        .surface
        .clone()
        .ok_or_else(|| CoreError::PreconditionFailed("face missing surface".into()))?;
    Ok((l1.face, l2.face, s1, s2))
}

/// Sample UV coordinates on both adjacent surfaces along an edge.
pub(super) fn sample_edge_uvs(
    brep: &BRep,
    edge: EdgeHandle,
    count: usize,
) -> Vec<(f64, f64, f64, f64)> {
    let Some(ed) = brep.edges.get(edge) else {
        return Vec::new();
    };
    let Some(he1) = brep.half_edges.get(ed.half_edge) else {
        return Vec::new();
    };
    let Some(he2) = brep.half_edges.get(he1.twin) else {
        return Vec::new();
    };
    let Some(l1) = brep.loops.get(he1.loop_) else {
        return Vec::new();
    };
    let Some(l2) = brep.loops.get(he2.loop_) else {
        return Vec::new();
    };
    let Some(f1) = brep.faces.get(l1.face) else {
        return Vec::new();
    };
    let Some(f2) = brep.faces.get(l2.face) else {
        return Vec::new();
    };

    (0..count)
        .map(|i| {
            let t = (i as f64 + 0.5) / count as f64;
            let (ua, va) = uv_at_fraction(he1, f1, t);
            let (ub, vb) = uv_at_fraction(he2, f2, t);
            (ua, va, ub, vb)
        })
        .collect()
}

/// Get UV at a fraction t in [0,1] along a half-edge's pcurve.
fn uv_at_fraction(
    he: &crate::topology::HalfEdge,
    face: &crate::topology::Face,
    t_frac: f64,
) -> (f64, f64) {
    if let Some(ref pc) = he.pcurve {
        let (t0, t1) = pc.parameter_range();
        let pt = pc.point_at(t0 + (t1 - t0) * t_frac);
        return (pt.x, pt.y);
    }
    if let Some(ref s) = face.surface {
        let ((u0, u1), (v0, v1)) = s.geometry.parameter_range();
        let u = if u0.is_finite() && u1.is_finite() {
            (u0 + u1) * 0.5
        } else {
            0.0
        };
        let v = if v0.is_finite() && v1.is_finite() {
            (v0 + v1) * 0.5
        } else {
            0.0
        };
        return (u, v);
    }
    (0.0, 0.0)
}

/// Measure maximum normal angle (radians) between two surfaces at samples.
fn measure_normal_angle(
    sa: &SurfaceGeometry,
    sb: &SurfaceGeometry,
    samples: &[(f64, f64, f64, f64)],
) -> f64 {
    let mut worst = 0.0_f64;
    for &(ua, va, ub, vb) in samples {
        let na = sa.normal_at(ua, va);
        let nb = sb.normal_at(ub, vb);
        let la = na.length();
        let lb = nb.length();
        if la < 1e-15 || lb < 1e-15 {
            continue;
        }
        let cos = (na.dot(&nb) / (la * lb)).clamp(-1.0, 1.0);
        worst = worst.max(cos.acos());
    }
    worst
}

/// Measure maximum Gaussian curvature difference between surfaces at samples.
fn measure_curvature_diff(
    sa: &SurfaceGeometry,
    sb: &SurfaceGeometry,
    samples: &[(f64, f64, f64, f64)],
) -> f64 {
    let mut worst = 0.0_f64;
    for &(ua, va, ub, vb) in samples {
        if let (Some(a), Some(b)) = (
            surface_curvature_at(sa, ua, va),
            surface_curvature_at(sb, ub, vb),
        ) {
            worst = worst.max((a.gaussian - b.gaussian).abs());
        }
    }
    worst
}

/// Compute the midpoint of two 3D points.
fn midpoint(a: &Point3, b: &Point3) -> Point3 {
    Point3::new((a.x + b.x) * 0.5, (a.y + b.y) * 0.5, (a.z + b.z) * 0.5)
}

#[cfg(test)]
#[path = "enforcement_tests.rs"]
mod tests;
