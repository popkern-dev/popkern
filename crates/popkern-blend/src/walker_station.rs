//! Station certificate for the spine walker's first sample (PK-419).
//!
//! The contact solver closes the rolling-ball system with the spine
//! plane — the plane through the spine point normal to the spine
//! tangent. That plane is a station only up to the spine's own
//! symmetry: on a closed circular rim the antipode lies in the same
//! plane, and a Newton solve seeded from the far side of the circle
//! converges to the antipodal ball with a clean residual. Every later
//! sample warm-starts from the first, so the whole walk — and the
//! surface of revolution built on it — then originates half a turn
//! from where the caller asked (measured on the flagship boss plate:
//! the closed-fillet surgery seamed the ring at the rail vertex while
//! the strip's `u = 0` sat at the antipode, 70 open mesh edges).
//!
//! The certificate is spine-only: the ball at station `t` is nearer
//! the spine point at `t` than to the spine point at any station a
//! good fraction of the span away. "Nearer" reads the closest of the
//! ball centre and its two contact points, so the test is
//! insensitive to the dihedral angle, and the mismatch must be gross
//! (a factor of two) so a legitimately tight spine never trips it.

use popkern_core::geometry::Point3;

use crate::contact::ContactSolution;
use crate::walker::{WalkRequest, WalkerError};

/// Stations sampled along the spine.
const STATIONS: usize = 32;
/// Stations closer than this fraction of the span are the same
/// neighbourhood as `t`.
const NEIGHBOURHOOD: f64 = 1.0 / 8.0;
/// A foreign station must be this much nearer to flag a mismatch.
const MARGIN: f64 = 0.5;

/// `Ok` when `solution` is the ball at spine parameter `t`;
/// [`WalkerError::StationMismatch`] naming the station it belongs to
/// otherwise.
pub(crate) fn certify_station(
    request: &WalkRequest,
    t: f64,
    solution: &ContactSolution,
) -> Result<(), WalkerError> {
    let (t0, t1) = request.t_range;
    let span = t1 - t0;
    if span == 0.0 || !span.is_finite() {
        return Ok(());
    }
    let reach = |p: Point3| -> f64 {
        solution
            .center
            .distance_to(&p)
            .min(solution.contact0.distance_to(&p))
            .min(solution.contact1.distance_to(&p))
    };
    let own = reach(request.spine.point_at(t));
    let ends = (request.spine.point_at(t0), request.spine.point_at(t1));
    let scale = ends
        .0
        .distance_to(&ends.1)
        .max(own)
        .max(solution.center.distance_to(&solution.contact0));
    let closed = ends.0.distance_to(&ends.1) <= 1.0e-9 * scale.max(1.0e-12);
    let mut nearest: Option<(f64, f64)> = None;
    for k in 0..=STATIONS {
        let tk = t0 + span * k as f64 / STATIONS as f64;
        let mut dt = ((tk - t) / span).abs();
        if closed {
            dt = dt.min(1.0 - dt);
        }
        if dt < NEIGHBOURHOOD {
            continue;
        }
        let d = reach(request.spine.point_at(tk));
        if d < MARGIN * own && nearest.is_none_or(|(_, best)| d < best) {
            nearest = Some((tk, d));
        }
    }
    match nearest {
        Some((nearest_t, _)) => Err(WalkerError::StationMismatch {
            t_spine: t,
            nearest_t,
            center: solution.center,
        }),
        None => Ok(()),
    }
}

#[cfg(test)]
#[path = "walker_station_tests.rs"]
mod tests;
