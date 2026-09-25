//! Certified seam continuity: rigorous upper bounds on the position
//! gap, normal-angle sine and second-fundamental-form difference over
//! a whole shared edge (PK-499 design §Algorithms).
//!
//! The seam fraction `[0, 1]` is covered by s-intervals. On each
//! interval `T` with midpoint `m`, every quantity component `q` —
//! a truncated Taylor series with coefficients `c_k` — is enclosed by
//! the intersection of the natural enclosure `c_0(T)` and the Taylor
//! forms `Σ_{k<j} c_k(m)(T − m)^k + c_j(T)(T − m)^j` for every order
//! `j` (order 1 is the mean-value form), except where a side's box
//! straddles a knot or closure, where a quantity may jump and only
//! the natural enclosure is sound. Refinement is branch and bound on
//! the supremum (see [`certify_seam`]).

use super::CertifyDecline;
use super::seam_interval::SideB;
use super::seam_projection::ProjectedSide;
use super::seam_refine::certify_with;
use super::seam_state::SeamSideSource;

/// Refinement targets: an interval is settled once each upper bound
/// exceeds the best certified witness `w` on the seam by at most
/// `target + relative · w`.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct CertifyTargets {
    /// Position gap (length units).
    pub gap: f64,
    /// Normal-angle sine (dimensionless).
    pub normal_sine: f64,
    /// Second-fundamental-form difference (1/length).
    pub curvature: f64,
    /// Slack relative to the best witness (dimensionless), so a
    /// seam with a real residual is bounded to a fraction of it
    /// instead of to the absolute floor.
    pub relative: f64,
}

/// Options for [`certify_seam`]. Every tolerance is explicit.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct SeamCertifyOptions {
    /// Refinement targets.
    pub targets: CertifyTargets,
    /// Largest certified gap under which the two sides are accepted
    /// as fraction-co-parameterized.
    pub correspondence_tolerance: f64,
    /// Most s-intervals evaluated per seam.
    pub budget: usize,
    /// Narrowest s-interval that is still bisected.
    pub min_width: f64,
    /// Narrowest s-interval bisected around a degenerate chart (a
    /// Jacobian enclosure reaching zero): a pole cannot be certified
    /// however fine the split, so its neighbourhood is reported as
    /// uncertified at this resolution.
    pub singular_min_width: f64,
    /// Stop at the first structural failure (a side that cannot be
    /// enclosed, a degenerate chart at the singular floor): such a
    /// certificate can never cover the seam, so a consumer that only
    /// uses covering bounds (the blend evidence) skips the remaining
    /// work. Default `false` keeps the partial bounds as diagnostics.
    pub stop_at_structural_failure: bool,
    /// Decision gates for (gap, sine, curvature). When set, refinement
    /// serves the decision against each gate rather than tight bounds:
    /// a quantity whose certified lower bound exceeds its gate
    /// (certified failure) stops refining, and an interval whose upper
    /// bound is within the gate (certified pass) is settled for it.
    /// `None` refines every bound to its target.
    pub decision_gates: Option<[f64; 3]>,
}

impl SeamCertifyOptions {
    /// Defaults for a model of characteristic length `scale` and
    /// point-coincidence (linear) tolerance `linear`: sides within
    /// `linear` of each other are one seam (correspondence), the gap
    /// is refined to `linear/10`, the sine to `1e-6`, curvature to
    /// `1e-3/scale`, relative slack `1e-3`, 4096 intervals.
    pub fn for_model(scale: f64, linear: f64) -> Self {
        Self {
            targets: CertifyTargets {
                gap: 0.1 * linear,
                normal_sine: 1e-6,
                curvature: 1e-3 / scale,
                relative: 1e-3,
            },
            correspondence_tolerance: linear,
            budget: 4096,
            min_width: 1e-12,
            singular_min_width: 1e-7,
            stop_at_structural_failure: false,
            decision_gates: None,
        }
    }
}

/// How the two sides were put in correspondence.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum SeamCorrespondence {
    /// Equal seam fractions are the same point (certified gap).
    Fraction,
    /// Side b is compared at the certified foot point of side a's
    /// seam point on side b's seam curve (Krawczyk; see
    /// [`super::seam_projection`]).
    Projected,
    /// Side b (closed) is compared at the phase-shifted fraction
    /// `σ = φ ± s (mod 1)`; like fraction correspondence, certified
    /// by the gap bound (PK-542).
    PhaseShifted,
}

/// Why a certificate is not a full certification.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum InconclusiveReason {
    /// A side could not be enclosed (index 0 = side a).
    Side(usize, CertifyDecline),
    /// A side's Jacobian enclosure reaches zero (a pole or collapsed
    /// chart).
    Jacobian,
    /// The sides are not fraction-co-parameterized within the
    /// correspondence tolerance.
    Correspondence,
    /// The interval budget ran out before every target was met.
    Budget,
    /// An interval reached the minimum width without meeting its
    /// target.
    Resolution,
    /// The edge is not a two-face seam with curve, pcurves and
    /// surfaces (reported by the topology wrapper).
    NotASeam,
}

/// Outcome status of a certificate.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub enum CertificateStatus {
    /// Every bound holds over the whole seam and met its target.
    Certified,
    /// Something could not be certified; `at` is the first failing
    /// s-interval. The bounds cover every interval that was
    /// enclosed, and `uncertified_length` says how much was not.
    Inconclusive {
        /// Why.
        reason: InconclusiveReason,
        /// Where (seam fractions).
        at: (f64, f64),
    },
}

/// Certified bounds for one seam.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct SeamCertificate {
    /// Upper bound of `|P_a − P_b|`.
    pub g0_upper: f64,
    /// Upper bound of the normal-angle sine.
    pub g1_upper: f64,
    /// Upper bound of `max |II_a − II_b|` over the shared frame's
    /// components.
    pub g2_upper: f64,
    /// Certified lower bounds of the sup of (gap, sine, curvature):
    /// values attained at some interval midpoint. A lower bound above
    /// a gate is a certified failure.
    pub lower: [f64; 3],
    /// Lower bound of each side's `|S_u × S_v|` along the seam.
    pub jacobian_lower: [f64; 2],
    /// Seam fractions of the intervals attaining `g0`, `g1`, `g2`.
    pub worst_at: [(f64, f64); 3],
    /// Correspondence used.
    pub correspondence: SeamCorrespondence,
    /// s-intervals evaluated.
    pub intervals: usize,
    /// Total seam fraction not covered by the bounds.
    pub uncertified_length: f64,
    /// Status.
    pub status: CertificateStatus,
}

impl SeamCertificate {
    /// A certificate that bounds nothing: every upper bound infinite,
    /// status `Inconclusive(NotASeam)` until a caller sets its reason.
    pub fn unbounded() -> Self {
        Self::inconclusive(InconclusiveReason::NotASeam)
    }

    pub(super) fn inconclusive(reason: InconclusiveReason) -> Self {
        Self {
            g0_upper: f64::INFINITY,
            g1_upper: f64::INFINITY,
            g2_upper: f64::INFINITY,
            jacobian_lower: [0.0, 0.0],
            lower: [0.0; 3],
            worst_at: [(0.0, 1.0); 3],
            correspondence: SeamCorrespondence::Fraction,
            intervals: 0,
            uncertified_length: 1.0,
            status: CertificateStatus::Inconclusive {
                reason,
                at: (0.0, 1.0),
            },
        }
    }

    /// True only for a full certification.
    pub fn is_certified(&self) -> bool {
        self.status == CertificateStatus::Certified
    }

    /// True when the upper bounds hold over the whole seam: fully
    /// certified, or stopped by budget or resolution with nothing left
    /// uncertified (the bounds are rigorous, only looser than the
    /// refinement targets). A structural failure (a side that cannot be
    /// enclosed, a degenerate chart, failed correspondence, not a seam)
    /// never covers the seam.
    pub fn bounds_cover_seam(&self) -> bool {
        match self.status {
            CertificateStatus::Certified => true,
            CertificateStatus::Inconclusive { reason, .. } => {
                matches!(
                    reason,
                    InconclusiveReason::Budget | InconclusiveReason::Resolution
                ) && self.uncertified_length == 0.0
                    && self.g0_upper.is_finite()
            }
        }
    }
}

/// Certify the seam between two prepared sides under fraction
/// correspondence.
///
/// Branch and bound on the supremum: every enclosed interval carries
/// upper bounds, and its midpoint carries a certified lower witness.
/// An interval is settled once each upper bound is within
/// `target + relative·witness` of the best witness found anywhere on
/// the seam; the interval overshooting most is bisected first. So an
/// exact seam refines until its bounds reach the absolute targets,
/// and a seam with a real residual refines only near its maximum.
pub fn certify_seam(
    a: &mut SeamSideSource,
    b: &mut SeamSideSource,
    options: &SeamCertifyOptions,
) -> SeamCertificate {
    // Fraction correspondence is a short probe (an eighth of the
    // budget, at least 64 intervals): exact co-parameterized seams
    // certify in ≤ 60 (PK-499 E1, PK-542), while sides co-parameterized only
    // up to a small real mismatch would spend the budget resolving a
    // parameterization difference that is not a geometric gap —
    // projection measures the geometric gap and settles quickly
    // (measured: boss-plate ring rails spent 250 of 256 evidence
    // intervals in fraction mode before projection certified them).
    let fraction_options = SeamCertifyOptions {
        budget: (options.budget / 8).max(64).min(options.budget),
        ..*options
    };
    let fraction = certify_with(a, &mut SideB::Fraction(b), &fraction_options);
    // Certified first, then bounds that cover the seam, then the
    // tighter gap (the better correspondence: a probe stopped by its
    // short budget otherwise shadows a refined alternative); a
    // certified fraction mismatch leaves only the other modes.
    let rank = |c: &SeamCertificate| (c.is_certified(), c.bounds_cover_seam());
    let loose_gap = fraction.g0_upper > options.correspondence_tolerance;
    let failed_correspondence = matches!(
        fraction.status,
        CertificateStatus::Inconclusive {
            reason: InconclusiveReason::Correspondence,
            ..
        }
    );
    if fraction.is_certified() {
        return fraction;
    }
    let probe_ran_out = matches!(
        fraction.status,
        CertificateStatus::Inconclusive {
            reason: InconclusiveReason::Budget,
            ..
        }
    );
    // The sides correspond at equal fractions when the gap bound is
    // already within tolerance, or when the probe merely ran out of
    // budget while every midpoint witness shows the sides coinciding
    // (a certified gap witness under a tenth of the tolerance): finish
    // in fraction mode. A genuine mismatch (a phase shift, a
    // rational-arc parameter) has witnesses above tolerance.
    let witnesses_coincide =
        probe_ran_out && fraction.lower[0] <= 0.1 * options.correspondence_tolerance;
    if !failed_correspondence && (!loose_gap || witnesses_coincide) {
        return if probe_ran_out && options.budget > fraction_options.budget {
            certify_with(a, &mut SideB::Fraction(b), options)
        } else {
            fraction
        };
    }
    let Ok(mut projected_side) = ProjectedSide::new(b.clone(), options.correspondence_tolerance)
    else {
        return fraction;
    };
    // A closed rail parameterized by angle on both sides differs only
    // by a phase (and direction): certify that affine correspondence
    // directly, at fraction-mode cost, before full projection
    // (measured: the corpus boss rings needed 339 projected intervals
    // for the G0 decision).
    let mut spent = fraction.intervals;
    let mut best = fraction;
    if let Some((phase, dir)) = projected_side.phase_estimate(a) {
        let shifted_options = SeamCertifyOptions {
            budget: options.budget.saturating_sub(spent).max(16) / 2,
            ..*options
        };
        let mut shifted = certify_with(
            a,
            &mut SideB::Shifted {
                side: &mut projected_side,
                phase,
                dir,
            },
            &shifted_options,
        );
        shifted.correspondence = SeamCorrespondence::PhaseShifted;
        spent += shifted.intervals;
        shifted.intervals = spent;
        if shifted.is_certified() {
            return shifted;
        }
        if rank(&shifted) >= rank(&best) && shifted.g0_upper < best.g0_upper {
            best = shifted;
        }
    }
    let fraction = best;
    let projected_options = SeamCertifyOptions {
        budget: options.budget.saturating_sub(spent).max(16),
        ..*options
    };
    let mut projected = certify_with(
        a,
        &mut SideB::Projected(&mut projected_side),
        &projected_options,
    );
    projected.correspondence = SeamCorrespondence::Projected;
    projected.intervals += spent;
    let better = rank(&projected) > rank(&fraction)
        || (rank(&projected) == rank(&fraction) && projected.g0_upper < fraction.g0_upper);
    if failed_correspondence || better {
        projected
    } else {
        fraction
    }
}
