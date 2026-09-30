//! Interior fairness of a setback corner cap (PK-612).
//!
//! The seam gates read a cap's rims and internal seams only. The KP20
//! fan's star pieces can pinch or ripple well inside the cap while every
//! seam passes. On the census's KP20 serves the sampled interior
//! curvature reached 4.5e5 / r_min (hexagon h1.2 variable, a near-fold
//! at piece regularity 2.4e-3), 51 (pyramid4 h0.8 wide) and 7.5 on an
//! equal-radius octagon h0.4, each a visible zebra ripple or dimple,
//! while the variational fit on the same corners measured 1.3–2.6.
//!
//! The measure is the largest principal curvature over interior samples
//! of every cap face, times the corner's smallest fillet radius: how far
//! the cap bends past its sharpest strip.

use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::healing::curvature::principal_curvatures;
use popkern_core::topology::{BRep, FaceHandle};

/// The largest interior curvature a served cap may reach, in units of
/// `1 / r_min`. Fair caps measure 1.1–2.6 (the fit on every census
/// corner it serves cleanly) and at most 3.85 (the pinned KP20 caps:
/// flagship 1.84, blunt 1.70, cube 1.76, pentagon h0.8 variable 3.75).
/// Every KP20 cap over 4 that was rendered showed a visible ripple or
/// dimple (octagon h0.4 at 7.5, hexagon h0.4 variable at 18, pyramid4
/// h0.8 wide at 51).
pub(crate) const CAP_FAIRNESS_BOUND: f64 = 4.0;

/// Interior samples per face direction (the endpoints, where a Gregory
/// corner blend is 0/0, are excluded).
const SAMPLES: usize = 32;

/// Smallest interior regularity (min/mean of `|S_u × S_v|` per face) a
/// served KP20 cap may have: the variational fit's fold floor
/// (`variational_cap::measure::SERVED_MIN_REGULARITY`, E6). On the
/// census only the hexagon h1.2 variable fan sits under it (2.4e-3, its
/// 4.5e5 pinch); every other KP20 serve reads 3.05e-2 or more.
pub(crate) const CAP_MIN_REGULARITY: f64 = 2.0e-2;

/// A cap's interior readings.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct CapInterior {
    /// Largest principal curvature times `r_min`.
    pub overshoot: f64,
    /// Smallest per-face min/mean of `|S_u × S_v|`.
    pub regularity: f64,
}

impl CapInterior {
    /// Whether the cap nearly folds somewhere inside.
    pub(crate) fn near_fold(&self) -> bool {
        self.regularity < CAP_MIN_REGULARITY
    }

    /// Whether the cap bends past the fairness bound or nearly folds.
    pub(crate) fn unfair(&self) -> bool {
        self.near_fold() || self.overshoot > CAP_FAIRNESS_BOUND
    }
}

/// The interior readings over `SAMPLES - 1` squared interior samples of
/// every face in `faces`. A face without a finite surface is skipped; a
/// sample at a degenerate point contributes no curvature.
pub(crate) fn cap_interior(brep: &BRep, faces: &[FaceHandle], r_min: f64) -> CapInterior {
    let mut worst = 0.0_f64;
    let mut regularity = f64::INFINITY;
    for &face in faces {
        let Some(surface) = brep.faces.get(face).and_then(|f| f.surface.as_ref()) else {
            continue;
        };
        let geometry: &dyn SurfaceEval = &surface.geometry;
        let ((u0, u1), (v0, v1)) = geometry.parameter_range();
        if ![u0, u1, v0, v1].iter().all(|x| x.is_finite()) {
            continue;
        }
        let (mut jac_min, mut jac_sum) = (f64::INFINITY, 0.0_f64);
        for i in 1..SAMPLES {
            for j in 1..SAMPLES {
                let u = u0 + (u1 - u0) * i as f64 / SAMPLES as f64;
                let v = v0 + (v1 - v0) * j as f64 / SAMPLES as f64;
                if let Some((k1, k2)) = principal_curvatures(geometry, u, v) {
                    worst = worst.max(k1.abs().max(k2.abs()));
                }
                let jac = geometry
                    .partial_u(u, v)
                    .cross(&geometry.partial_v(u, v))
                    .length();
                jac_min = jac_min.min(jac);
                jac_sum += jac;
            }
        }
        let mean = jac_sum / ((SAMPLES - 1) * (SAMPLES - 1)) as f64;
        if mean > 0.0 {
            regularity = regularity.min(jac_min / mean);
        }
    }
    CapInterior {
        overshoot: worst * r_min,
        regularity,
    }
}
