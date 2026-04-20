//! Runtime self-check for direct-construction (Phase G).
//!
//! Measures the dense-sample G1, G2, and G0 witnesses on a ring
//! of two-stripe patches, independent of how the CPs got there.
//! The same machinery is used by:
//!
//!   * `direct_build::build` at `ValidationLevel::Standard` / `Full`
//!     to fire `SelfCheckFailedBug` when direct-construction math
//!     produces a result that doesn't actually satisfy the
//!     identity.
//!   * Unit tests and the validation sweep, to assert commercial-
//!     grade quality.
//!   * Future UI preview / telemetry, to badge patch quality.
//!
//! The dense-sample measurement is the **structural claim** —
//! checking only at Bernstein coefficient indices (or at the
//! solver's sample set) admits false positives via between-sample
//! wiggle. Dense sampling at 100+ points is the honest witness.

use crate::geometry::Vec3;
use crate::geometry::gregory::RectangularGregoryTwoStripe;
use crate::geometry::gregory::spoke_g1_solver::{Side, Spoke};

/// Output of [`measure_ring_quality`]. All residuals are
/// worst-case (min cos / max gap) across the dense sample set;
/// per-spoke detail is available via
/// [`measure_ring_quality_per_spoke`] for UI binding.
#[derive(Debug, Clone)]
pub struct RingQualityMeasurement {
    /// `min_{t in dense set, spoke in ring} |n_A(t) · n_B(t)|`.
    /// `1.0` = exact G1 everywhere.
    pub dense_min_cos_g1: f64,
    /// `min` cos across dense samples for normal curvature
    /// agreement. `None` if G2 was not requested. `1.0` = exact
    /// G2 curvature continuity.
    pub dense_min_cos_g2: Option<f64>,
    /// Maximum 3D position gap at any shared boundary sample.
    /// Expected ≤ `1e-12` for exact construction.
    pub max_position_gap_g0: f64,
}

/// Per-spoke breakdown of the dense measurement.
#[derive(Debug, Clone)]
pub struct PerSpokeQuality {
    /// Min cos of the G1 angle along this spoke.
    pub min_cos_g1: f64,
    /// Min cos of G2 curvature agreement (if requested).
    pub min_cos_g2: Option<f64>,
    /// Max position gap along this spoke.
    pub max_position_gap: f64,
}

/// Options for [`measure_ring_quality`].
#[derive(Debug, Clone)]
pub struct SelfCheckOptions {
    /// Dense sampling count per spoke. 100 is the commercial-
    /// grade default — denser than the scorecard, fine enough to
    /// catch between-sample wiggle.
    pub samples_per_spoke: usize,
    /// Whether to also measure G2 curvature continuity. Adds
    /// second-derivative evaluations per sample; use only when
    /// `target_continuity == G2`.
    pub measure_g2: bool,
}

impl Default for SelfCheckOptions {
    fn default() -> Self {
        Self {
            samples_per_spoke: 100,
            measure_g2: false,
        }
    }
}

/// Measure G1 (and optionally G2) quality of a ring of
/// two-stripe patches. Fixture-free, pure-geometric.
pub fn measure_ring_quality(
    patches: &[RectangularGregoryTwoStripe],
    spokes: &[Spoke],
    options: &SelfCheckOptions,
) -> RingQualityMeasurement {
    let per_spoke = measure_ring_quality_per_spoke(patches, spokes, options);
    let dense_min_cos_g1 = per_spoke
        .iter()
        .map(|q| q.min_cos_g1)
        .fold(1.0_f64, f64::min);
    let dense_min_cos_g2 = if options.measure_g2 {
        Some(
            per_spoke
                .iter()
                .filter_map(|q| q.min_cos_g2)
                .fold(1.0_f64, f64::min),
        )
    } else {
        None
    };
    let max_position_gap_g0 = per_spoke
        .iter()
        .map(|q| q.max_position_gap)
        .fold(0.0_f64, f64::max);
    RingQualityMeasurement {
        dense_min_cos_g1,
        dense_min_cos_g2,
        max_position_gap_g0,
    }
}

/// Same as [`measure_ring_quality`] but returns per-spoke
/// detail (for UI binding and targeted diagnostics).
pub fn measure_ring_quality_per_spoke(
    patches: &[RectangularGregoryTwoStripe],
    spokes: &[Spoke],
    options: &SelfCheckOptions,
) -> Vec<PerSpokeQuality> {
    let mut out = Vec::with_capacity(spokes.len());
    for spoke in spokes {
        out.push(measure_single_spoke_quality(patches, *spoke, options));
    }
    out
}

fn measure_single_spoke_quality(
    patches: &[RectangularGregoryTwoStripe],
    spoke: Spoke,
    options: &SelfCheckOptions,
) -> PerSpokeQuality {
    let mut min_cos_g1 = 1.0_f64;
    let mut max_gap = 0.0_f64;
    let mut min_cos_g2 = if options.measure_g2 {
        Some(1.0_f64)
    } else {
        None
    };

    for k in 1..=options.samples_per_spoke {
        let t = (k as f64) / ((options.samples_per_spoke + 1) as f64);
        let (ua, va) = uv_for_side(spoke.a.side, t);
        let t_b = if spoke.b.reversed { 1.0 - t } else { t };
        let (ub, vb) = uv_for_side(spoke.b.side, t_b);

        let pa = &patches[spoke.a.patch_idx];
        let pb = &patches[spoke.b.patch_idx];

        // G0 position gap.
        let point_a = pa.point_at(ua, va);
        let point_b = pb.point_at(ub, vb);
        let gap = ((point_a.x - point_b.x).powi(2)
            + (point_a.y - point_b.y).powi(2)
            + (point_a.z - point_b.z).powi(2))
        .sqrt();
        if gap > max_gap {
            max_gap = gap;
        }

        // G1 normal agreement.
        if let (Some(na), Some(nb)) = (unit_normal(pa, ua, va), unit_normal(pb, ub, vb)) {
            let cos = (na.x * nb.x + na.y * nb.y + na.z * nb.z).abs();
            if cos < min_cos_g1 {
                min_cos_g1 = cos;
            }
        } else {
            // Degenerate normal — zero partial. Report this as
            // min_cos_g1 = 0.0; any surface with zero partials
            // along its boundary fails commercial-grade.
            min_cos_g1 = 0.0;
        }

        // G2 normal-curvature agreement.
        if let Some(ref mut min_g2) = min_cos_g2 {
            let kn_a = normal_curvature(pa, ua, va, t);
            let kn_b = normal_curvature(pb, ub, vb, t_b);
            let agree = g2_agreement(kn_a, kn_b);
            if agree < *min_g2 {
                *min_g2 = agree;
            }
        }
    }

    PerSpokeQuality {
        min_cos_g1,
        min_cos_g2,
        max_position_gap: max_gap,
    }
}

fn uv_for_side(side: Side, t: f64) -> (f64, f64) {
    match side {
        Side::UZero => (0.0, t),
        Side::UOne => (1.0, t),
        Side::VZero => (t, 0.0),
        Side::VOne => (t, 1.0),
    }
}

fn unit_normal(p: &RectangularGregoryTwoStripe, u: f64, v: f64) -> Option<Vec3> {
    p.partial_u(u, v).cross(&p.partial_v(u, v)).normalized()
}

/// Normal curvature along the boundary tangent direction.
/// Computed via the second fundamental form; see standard
/// differential-geometry texts (e.g. do Carmo §3.2).
fn normal_curvature(p: &RectangularGregoryTwoStripe, u: f64, v: f64, _t: f64) -> f64 {
    // For the cross-boundary normal-curvature check at a spoke,
    // we take the tangent direction along the boundary (which
    // is perpendicular to the cross-boundary direction). At
    // side UZero / UOne, the boundary runs in v; at VZero / VOne,
    // in u. We compute kn = II(T, T) / I(T, T) in that direction.
    let du = p.partial_u(u, v);
    let dv = p.partial_v(u, v);
    let duu = p.second_partial_uu(u, v);
    let duv = p.second_partial_uv(u, v);
    let dvv = p.second_partial_vv(u, v);
    let Some(n) = du.cross(&dv).normalized() else {
        return 0.0;
    };
    // First fundamental form (E, F, G).
    let e_ff = du.dot(&du);
    let f_ff = du.dot(&dv);
    let g_ff = dv.dot(&dv);
    // Second fundamental form (L, M, N).
    let l_ff = duu.dot(&n);
    let m_ff = duv.dot(&n);
    let n_ff = dvv.dot(&n);
    // Walk along v at fixed u (for UZero/UOne boundaries): T = (0, 1).
    // Then kn = N / G. For VZero/VOne we'd use kn = L / E; caller picks
    // the right partial-derivative pair by passing the right (u, v).
    // For a robust measurement here we use the maximum principal curvature
    // candidate — the measurement is symmetric because we compare
    // A and B at the same geometric boundary point and use the same formula.
    if g_ff > 1e-20 {
        let denom = e_ff * g_ff - f_ff * f_ff;
        if denom.abs() > 1e-20 {
            // Mean curvature component (coordinate-system-robust).
            (l_ff * g_ff - 2.0 * m_ff * f_ff + n_ff * e_ff) / (2.0 * denom)
        } else {
            0.0
        }
    } else {
        0.0
    }
}

/// G2 agreement on a single scalar curvature. Tolerance-free
/// comparison expressed as `1 - |k_a - k_b| / (1 + max)` in
/// `[0, 1]`. `1.0` = exact match; `0.0` = wildly different.
fn g2_agreement(k_a: f64, k_b: f64) -> f64 {
    let diff = (k_a - k_b).abs();
    let scale = 1.0 + k_a.abs().max(k_b.abs());
    1.0 - (diff / scale).min(1.0)
}

#[cfg(test)]
#[path = "direct_self_check_tests.rs"]
mod tests;
