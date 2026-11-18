//! The atlas's cap rows (design §6a, step 1): around each chart
//! singularity, the parameter row bounding the cap a pole projection
//! charts, where the surface normal has turned 45° from the singular
//! point's. Between the distortion radius (speed ratio ten, ~5.7° on a
//! sphere) and the cap's 60° bound, with room either side to move the
//! row off a near-tangency with a loop.

use popkern_math::PortableMath;

use crate::geometry::Vec3;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};

use super::descriptor::{Axis, ChartDescriptor, SingularKind, SingularRow};
use super::regular::{RegularChart, pole_projection};

/// The normal turn at which a cap row is placed.
pub const CAP_TURN: f64 = std::f64::consts::FRAC_PI_4;

/// The largest normal turn a cap may reach (its projection's Jacobian at
/// least cos 60°, so chord height in the cap chart is within two of
/// space).
pub const CAP_LIMIT: f64 = std::f64::consts::FRAC_PI_3;

/// Longitudes at which a row's normals are read.
const ROW_SAMPLES: usize = 16;

/// A cap row: the chart singularity it bounds and the row parameter.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CapRow {
    /// The singular row the cap is about.
    pub singular: SingularRow,
    /// The cut row's value of the singular row's coordinate (`v` for a
    /// row `v = at`, `u` for a row `u = at`).
    pub at: f64,
    /// Whether the whole surface lies within the cap (its normal never
    /// turns `CAP_TURN` across the domain): no band beyond the row.
    pub whole: bool,
}

/// Why no cap row could be placed.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum CapRowError {
    /// The singular row is an apex (no pole projection charts it).
    #[error("singular row at {at} is a geometric singularity")]
    Geometric {
        /// The row parameter.
        at: f64,
    },
    /// The row's normal could not be read.
    #[error("singular row at {at}: no normal")]
    NoNormal {
        /// The row parameter.
        at: f64,
    },
    /// The domain across the row is unbounded.
    #[error("singular row at {at}: unbounded domain across it")]
    Unbounded {
        /// The row parameter.
        at: f64,
    },
}

/// The worst normal turn from `n0` along the row `coordinate = t`.
fn row_turn(d: &ChartDescriptor, s: &SurfaceGeometry, axis: Axis, t: f64, n0: Vec3) -> f64 {
    let ((u0, u1), (v0, v1)) = d.domain;
    let (lo, hi) = match axis {
        Axis::V => (u0, u1),
        Axis::U => (v0, v1),
    };
    let (lo, hi) = if lo.is_finite() && hi.is_finite() {
        (lo, hi)
    } else {
        (0.0, std::f64::consts::TAU)
    };
    (0..ROW_SAMPLES)
        .map(|k| {
            let a = lo + (hi - lo) * (k as f64 + 0.5) / ROW_SAMPLES as f64;
            let (u, v) = match axis {
                Axis::V => (a, t),
                Axis::U => (t, a),
            };
            let (su, sv) = s.first_partials(u, v);
            su.cross(&sv)
                .normalized()
                .map_or(std::f64::consts::PI, |n| {
                    n.dot(&n0).clamp(-1.0, 1.0).pacos()
                })
        })
        .fold(0.0, f64::max)
}

/// The cap row about `row` (design §6a step 1): the row parameter, moving
/// away from the singular row across the domain, at which the worst
/// normal turn reaches [`CAP_TURN`]; bisected to a millionth of the span.
///
/// # Errors
///
/// [`CapRowError`]: an apex, an unreadable normal, an unbounded span.
pub fn cap_row(
    d: &ChartDescriptor,
    s: &SurfaceGeometry,
    row: &SingularRow,
) -> Result<CapRow, CapRowError> {
    cap_row_at(d, s, row, CAP_TURN)
}

/// The turns a consumer may place a cap row at, preferred first: the cap
/// row may move off a near-tangency with a loop (a loop running along the
/// row flips its crossings as its samples refine), within the distortion
/// radius's speed ratio and the cap's 60° bound.
pub const CAP_TURNS: [f64; 5] = [
    std::f64::consts::FRAC_PI_4,
    40.0 * std::f64::consts::PI / 180.0,
    50.0 * std::f64::consts::PI / 180.0,
    35.0 * std::f64::consts::PI / 180.0,
    55.0 * std::f64::consts::PI / 180.0,
];

/// [`cap_row`] at the normal turn `turn` (one of [`CAP_TURNS`]).
///
/// # Errors
///
/// [`CapRowError`]: an apex, an unreadable normal, an unbounded span.
pub fn cap_row_at(
    d: &ChartDescriptor,
    s: &SurfaceGeometry,
    row: &SingularRow,
    turn_at: f64,
) -> Result<CapRow, CapRowError> {
    if row.kind != SingularKind::Chart {
        return Err(CapRowError::Geometric { at: row.at });
    }
    let Some(RegularChart::PoleProjection { e1, e2, .. }) = pole_projection(d, s, row) else {
        return Err(CapRowError::NoNormal { at: row.at });
    };
    let n0 = e1.cross(&e2);
    let ((u0, u1), (v0, v1)) = d.domain;
    let (lo, hi) = match row.axis {
        Axis::V => (v0, v1),
        Axis::U => (u0, u1),
    };
    if !(lo.is_finite() && hi.is_finite()) {
        return Err(CapRowError::Unbounded { at: row.at });
    }
    // Away from the singular row: toward the domain's far end.
    let far = if (row.at - lo).abs() <= (row.at - hi).abs() {
        hi
    } else {
        lo
    };
    let turn = |t: f64| row_turn(d, s, row.axis, t, n0);
    // The first of 64 steps across the span past the turn brackets it.
    let steps = 64;
    let mut inside = row.at;
    let mut outside = None;
    for k in 1..=steps {
        let t = row.at + (far - row.at) * k as f64 / steps as f64;
        if turn(t) >= turn_at {
            outside = Some(t);
            break;
        }
        inside = t;
    }
    let Some(mut outside) = outside else {
        return Ok(CapRow {
            singular: *row,
            at: far,
            whole: true,
        });
    };
    let span = (far - row.at).abs();
    while (outside - inside).abs() > 1e-6 * span {
        let mid = 0.5 * (inside + outside);
        if turn(mid) >= turn_at {
            outside = mid;
        } else {
            inside = mid;
        }
    }
    Ok(CapRow {
        singular: *row,
        at: 0.5 * (inside + outside),
        whole: false,
    })
}

/// The regular chart a cap is meshed in: the projection about its
/// singular row (`None` when the row's normal cannot be read).
pub fn cap_chart(d: &ChartDescriptor, s: &SurfaceGeometry, cap: &CapRow) -> Option<RegularChart> {
    pole_projection(d, s, &cap.singular)
}

/// How far `to` lies from `from` walking along a row in the direction
/// `walk` (`±1`), whole periods removed when the row closes (`period`):
/// in `[0, period)` then, else the plain signed distance in that sense.
pub fn along_offset(period: Option<f64>, from: f64, to: f64, walk: f64) -> f64 {
    let d = (to - from) * walk;
    match period {
        Some(p) => popkern_math::periodic::modulo(d, p),
        None => d,
    }
}

/// `a`'s representative on a closed row's canonical branch
/// `[start, start + period)` (the row's own parameter range), or `a` when
/// the row does not close.
pub fn canonical_along(period: Option<f64>, start: f64, a: f64) -> f64 {
    match period {
        Some(p) => start + popkern_math::periodic::modulo(a - start, p),
        None => a,
    }
}
