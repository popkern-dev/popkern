//! Generic trim-domain classification for bounded surface pcurves.
//!
//! This module answers a narrow question before topology repair or export trusts
//! a pcurve: does the pcurve live inside the owning surface parameter domain,
//! and if so is it a canonical boundary subspan or an interior trim?

use serde::{Deserialize, Serialize};

use crate::geometry::edge_triple::BoundarySide;
use crate::geometry::pcurve::Pcurve;
use crate::geometry::{Point3, SurfaceEval, SurfaceGeometry};

/// Accepted trim-domain class for a pcurve on a bounded surface.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum TrimDomainCandidateKind {
    /// The pcurve lies on one rectangular parameter-domain boundary.
    BoundarySubspan {
        /// Canonical rectangular boundary containing the trim subspan.
        side: BoundarySide,
    },
    /// The pcurve is inside the bounded parameter domain but not on a canonical
    /// boundary. It is valid trim geometry, not boundary ownership evidence.
    Interior,
}

/// Validated trim-domain candidate.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct TrimDomainCandidate {
    /// Accepted domain class.
    pub kind: TrimDomainCandidateKind,
    /// Largest sampled distance outside the finite UV parameter box.
    pub max_domain_overshoot: f64,
    /// Closest rectangular boundary by sampled UV distance.
    pub best_boundary: Option<BoundarySide>,
    /// Sampled deviation from [`Self::best_boundary`].
    pub best_boundary_deviation: f64,
    /// Largest 3-D gap between evaluated pcurve endpoints and topology anchors.
    pub endpoint_gap: f64,
    /// Number of pcurve samples used for classification.
    pub samples: usize,
}

/// Machine-readable reason a trim-domain candidate was declined.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum TrimDomainDeclineReason {
    /// The owning surface does not expose a finite rectangular parameter box.
    NonFiniteDomain,
    /// The pcurve parameter interval is invalid or near-zero.
    DegeneratePcurveRange,
    /// At least one sampled UV point lies outside the finite parameter box.
    OutsideDomain,
    /// Pcurve endpoints do not evaluate to the supplied topology anchors.
    EndpointMismatch,
}

/// Structured trim-domain decline.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct TrimDomainDecline {
    /// Machine-readable decline reason.
    pub reason: TrimDomainDeclineReason,
    /// Largest sampled distance outside the finite UV parameter box.
    pub max_domain_overshoot: f64,
    /// Closest rectangular boundary by sampled UV distance.
    pub best_boundary: Option<BoundarySide>,
    /// Sampled deviation from [`Self::best_boundary`].
    pub best_boundary_deviation: f64,
    /// Largest 3-D gap between evaluated pcurve endpoints and topology anchors.
    pub endpoint_gap: f64,
    /// Number of pcurve samples used for classification.
    pub samples: usize,
}

/// Classify and validate a pcurve against a bounded surface parameter domain.
///
/// The tolerance is local to the caller. Boundary ownership is accepted only
/// when the sampled pcurve lies on a canonical side. Interior trims are accepted
/// only when every sample is inside the finite parameter box and the pcurve
/// endpoints evaluate back to the supplied topology anchors.
pub fn classify_trim_domain_candidate(
    surface: &SurfaceGeometry,
    pcurve: &Pcurve,
    anchors: Option<(Point3, Point3)>,
    tolerance: f64,
    samples: usize,
) -> Result<TrimDomainCandidate, TrimDomainDecline> {
    let samples = samples.max(2);
    let tolerance = tolerance.max(0.0);
    let (t0, t1) = pcurve.parameter_range();
    if !t0.is_finite() || !t1.is_finite() || (t1 - t0).abs() <= f64::EPSILON {
        return Err(empty_decline(
            TrimDomainDeclineReason::DegeneratePcurveRange,
            samples,
        ));
    }
    let Some(domain) = finite_domain(surface) else {
        return Err(empty_decline(
            TrimDomainDeclineReason::NonFiniteDomain,
            samples,
        ));
    };
    let metrics = sample_trim_domain(surface, pcurve, anchors, domain, samples);
    if metrics.max_domain_overshoot > tolerance {
        return Err(metrics.decline(TrimDomainDeclineReason::OutsideDomain, samples));
    }
    if metrics.endpoint_gap > tolerance {
        return Err(metrics.decline(TrimDomainDeclineReason::EndpointMismatch, samples));
    }
    let kind = if metrics.best_boundary_deviation <= tolerance {
        let Some(side) = metrics.best_boundary else {
            return Err(metrics.decline(TrimDomainDeclineReason::OutsideDomain, samples));
        };
        TrimDomainCandidateKind::BoundarySubspan { side }
    } else {
        TrimDomainCandidateKind::Interior
    };
    Ok(TrimDomainCandidate {
        kind,
        max_domain_overshoot: metrics.max_domain_overshoot,
        best_boundary: metrics.best_boundary,
        best_boundary_deviation: metrics.best_boundary_deviation,
        endpoint_gap: metrics.endpoint_gap,
        samples,
    })
}

#[derive(Copy, Clone)]
struct Domain {
    u_min: f64,
    u_max: f64,
    v_min: f64,
    v_max: f64,
}

#[derive(Copy, Clone)]
struct TrimMetrics {
    max_domain_overshoot: f64,
    best_boundary: Option<BoundarySide>,
    best_boundary_deviation: f64,
    endpoint_gap: f64,
}

impl TrimMetrics {
    fn decline(self, reason: TrimDomainDeclineReason, samples: usize) -> TrimDomainDecline {
        TrimDomainDecline {
            reason,
            max_domain_overshoot: self.max_domain_overshoot,
            best_boundary: self.best_boundary,
            best_boundary_deviation: self.best_boundary_deviation,
            endpoint_gap: self.endpoint_gap,
            samples,
        }
    }
}

fn finite_domain(surface: &SurfaceGeometry) -> Option<Domain> {
    let ((u_min, u_max), (v_min, v_max)) = surface.parameter_range();
    [u_min, u_max, v_min, v_max]
        .into_iter()
        .all(f64::is_finite)
        .then_some(Domain {
            u_min,
            u_max,
            v_min,
            v_max,
        })
}

fn sample_trim_domain(
    surface: &SurfaceGeometry,
    pcurve: &Pcurve,
    anchors: Option<(Point3, Point3)>,
    domain: Domain,
    samples: usize,
) -> TrimMetrics {
    let (t0, t1) = pcurve.parameter_range();
    let mut overshoot = 0.0_f64;
    let mut u_min_dev = 0.0_f64;
    let mut u_max_dev = 0.0_f64;
    let mut v_min_dev = 0.0_f64;
    let mut v_max_dev = 0.0_f64;
    for index in 0..samples {
        let s = index as f64 / (samples - 1) as f64;
        let uv = pcurve.point_at(t0 + (t1 - t0) * s);
        overshoot = overshoot.max(domain_overshoot(uv, domain));
        u_min_dev = u_min_dev.max((uv.x - domain.u_min).abs());
        u_max_dev = u_max_dev.max((uv.x - domain.u_max).abs());
        v_min_dev = v_min_dev.max((uv.y - domain.v_min).abs());
        v_max_dev = v_max_dev.max((uv.y - domain.v_max).abs());
    }
    let (best_boundary, best_boundary_deviation) = [
        (BoundarySide::UMin, u_min_dev),
        (BoundarySide::UMax, u_max_dev),
        (BoundarySide::VMin, v_min_dev),
        (BoundarySide::VMax, v_max_dev),
    ]
    .into_iter()
    .min_by(|a, b| a.1.total_cmp(&b.1))
    .map(|(side, deviation)| (Some(side), deviation))
    .unwrap_or((None, f64::INFINITY));
    TrimMetrics {
        max_domain_overshoot: overshoot,
        best_boundary,
        best_boundary_deviation,
        endpoint_gap: endpoint_gap(surface, pcurve, anchors),
    }
}

fn domain_overshoot(uv: Point3, domain: Domain) -> f64 {
    (domain.u_min - uv.x)
        .max(uv.x - domain.u_max)
        .max(domain.v_min - uv.y)
        .max(uv.y - domain.v_max)
        .max(0.0)
}

fn endpoint_gap(
    surface: &SurfaceGeometry,
    pcurve: &Pcurve,
    anchors: Option<(Point3, Point3)>,
) -> f64 {
    let Some((start, end)) = anchors else {
        return 0.0;
    };
    let start_uv = pcurve.point_at(pcurve.t_start);
    let end_uv = pcurve.point_at(pcurve.t_end);
    let surface_start = surface.point_at(start_uv.x, start_uv.y);
    let surface_end = surface.point_at(end_uv.x, end_uv.y);
    start
        .distance_to(&surface_start)
        .max(end.distance_to(&surface_end))
}

fn empty_decline(reason: TrimDomainDeclineReason, samples: usize) -> TrimDomainDecline {
    TrimDomainDecline {
        reason,
        max_domain_overshoot: f64::INFINITY,
        best_boundary: None,
        best_boundary_deviation: f64::INFINITY,
        endpoint_gap: f64::INFINITY,
        samples,
    }
}

#[cfg(test)]
#[path = "trim_domain_tests.rs"]
mod tests;
