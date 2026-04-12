//! G3 (curvature-rate) continuity enforcement.
//!
//! Adjusts the third row of NURBS control points perpendicular to
//! the boundary to match the derivative of curvature across the
//! junction. Requires both surfaces to be NURBS with degree >= 5.

use crate::error::CoreError;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::healing::curvature::surface_curvature::surface_curvature_at;
use crate::healing::quality::ContinuityLevel;
use crate::tolerance::Tolerance;
use crate::topology::{BRep, EdgeHandle};

use super::enforcement::{self, EnforcementResult};

/// Enforce G3 continuity at a shared edge.
///
/// First ensures G2 is achieved, then adjusts the third row of control
/// points to match curvature derivatives. Requires both surfaces to be
/// NURBS with degree >= 5 (sufficient freedom for third-order matching).
/// If degree is insufficient, returns with G2 achieved and a graceful
/// degradation message.
///
/// **Acceptance (from spec):** curvature-derivative deviation < 1% of
/// maximum curvature magnitude; graceful degradation if representation
/// lacks sufficient freedom.
pub fn enforce_g3(
    brep: &mut BRep,
    edge: EdgeHandle,
    tolerance: &Tolerance,
    max_displacement: f64,
) -> Result<EnforcementResult, CoreError> {
    // Step 1: Ensure G2 is achieved first.
    let g2 = enforcement::enforce_g2(brep, edge, tolerance, max_displacement)?;
    if g2.achieved_level < ContinuityLevel::G2 {
        // Graceful degradation: cannot reach G3 without G2.
        return Ok(g2);
    }

    // Step 2: Compute remaining displacement budget.
    let budget = max_displacement - g2.max_displacement;
    if budget <= 0.0 {
        return Ok(g2);
    }

    // Step 3: Check both surfaces are NURBS with degree >= 5.
    let (f1, f2, sa, sb) = enforcement::get_adjacent_surfaces(brep, edge)?;
    let has_g3_freedom = |sg: &SurfaceGeometry| matches!(sg, SurfaceGeometry::Nurbs(n) if n.degree_u >= 5 && n.degree_v >= 5);
    if !has_g3_freedom(&sa.geometry) || !has_g3_freedom(&sb.geometry) {
        // Graceful degradation: insufficient degree for G3.
        return Ok(g2);
    }

    // Step 4: Sample UV points along the edge.
    let samples = enforcement::sample_edge_uvs(brep, edge, 8);
    if samples.is_empty() {
        return Ok(g2);
    }

    // Step 5: Measure curvature-rate difference before adjustment.
    let (rate_diff, max_curv_mag) =
        measure_curvature_rate_diff(&sa.geometry, &sb.geometry, &samples);
    let threshold = curvature_rate_threshold(max_curv_mag);

    if rate_diff <= threshold {
        // Already G3 continuous.
        return Ok(EnforcementResult {
            edge,
            max_displacement: g2.max_displacement,
            achieved_level: ContinuityLevel::G3,
        });
    }

    // Step 6: Adjust third row of control points.
    let g3_disp = super::nurbs_adjust::adjust_for_g3(brep, f1, f2, &samples, budget);

    // Step 7: Re-evaluate curvature rate after adjustment.
    let (_, _, sa2, sb2) = enforcement::get_adjacent_surfaces(brep, edge)?;
    let (post_diff, post_mag) = measure_curvature_rate_diff(&sa2.geometry, &sb2.geometry, &samples);
    let post_threshold = curvature_rate_threshold(post_mag);

    let achieved = if post_diff <= post_threshold {
        ContinuityLevel::G3
    } else {
        ContinuityLevel::G2
    };

    Ok(EnforcementResult {
        edge,
        max_displacement: g2.max_displacement + g3_disp,
        achieved_level: achieved,
    })
}

/// Compute the curvature-rate acceptance threshold.
///
/// The spec requires curvature-derivative deviation < 1% of
/// the maximum curvature magnitude. Uses a minimum floor to
/// avoid zero-threshold on flat surfaces.
fn curvature_rate_threshold(max_curvature_magnitude: f64) -> f64 {
    let relative = 0.01 * max_curvature_magnitude;
    // Floor at 1e-6 to handle near-flat surfaces.
    relative.max(1e-6)
}

/// Measure maximum curvature-rate difference between two surfaces.
///
/// For consecutive sample pairs along the edge, computes the Gaussian
/// curvature at each point on both surfaces and estimates dk/ds via
/// finite difference. Returns `(max_rate_diff, max_curvature_magnitude)`.
fn measure_curvature_rate_diff(
    sa: &SurfaceGeometry,
    sb: &SurfaceGeometry,
    samples: &[(f64, f64, f64, f64)],
) -> (f64, f64) {
    if samples.len() < 2 {
        return (0.0, 0.0);
    }

    let mut max_rate_diff = 0.0_f64;
    let mut max_curv_mag = 0.0_f64;

    // Collect curvatures at each sample on both surfaces.
    let curvatures_a: Vec<Option<f64>> = samples
        .iter()
        .map(|&(ua, va, _ub, _vb)| surface_curvature_at(sa, ua, va).map(|c| c.gaussian))
        .collect();
    let curvatures_b: Vec<Option<f64>> = samples
        .iter()
        .map(|&(_ua, _va, ub, vb)| surface_curvature_at(sb, ub, vb).map(|c| c.gaussian))
        .collect();

    // Track max curvature magnitude.
    for k in curvatures_a.iter().chain(curvatures_b.iter()).flatten() {
        max_curv_mag = max_curv_mag.max(k.abs());
    }

    // Estimate curvature rate via central differences between consecutive samples.
    for i in 0..samples.len() - 1 {
        let (ka_cur, ka_next) = match (curvatures_a[i], curvatures_a[i + 1]) {
            (Some(a), Some(b)) => (a, b),
            _ => continue,
        };
        let (kb_cur, kb_next) = match (curvatures_b[i], curvatures_b[i + 1]) {
            (Some(a), Some(b)) => (a, b),
            _ => continue,
        };

        // Approximate arc-length step from the 3D positions.
        let ds = estimate_arc_step(sa, samples[i], samples[i + 1]);
        if ds < 1e-15 {
            continue;
        }

        let rate_a = (ka_next - ka_cur) / ds;
        let rate_b = (kb_next - kb_cur) / ds;
        max_rate_diff = max_rate_diff.max((rate_a - rate_b).abs());
    }

    (max_rate_diff, max_curv_mag)
}

/// Estimate the arc-length step between two sample points using surface A.
fn estimate_arc_step(
    sa: &SurfaceGeometry,
    s1: (f64, f64, f64, f64),
    s2: (f64, f64, f64, f64),
) -> f64 {
    let p1 = sa.point_at(s1.0, s1.1);
    let p2 = sa.point_at(s2.0, s2.1);
    p1.distance_to(&p2)
}
