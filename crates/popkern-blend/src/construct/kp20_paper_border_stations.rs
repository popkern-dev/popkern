//! Strip stations for the bridged setback border (PK-611).
//!
//! On the corners-at-crossings border a strip's arc runs between the
//! rail crossings its two neighbours fix. A thin strip between wide
//! neighbours then gets corners far apart along its axis, and its arc
//! takes up that displacement mid-arc, across its own trimmed boundary
//! (the "setback arcs overrun a narrow strip" decline).
//!
//! The bridged border instead sets each strip back to one station of
//! its own, a pure cross-section. That station lies past both of the
//! strip's rail crossings, so the rails diverge there onto clean
//! support. A bridge on each support face then joins consecutive
//! strips' feet (see `junctions`), and the stations balance each
//! bridge's two legs (see [`balanced_stations`]).

use popkern_core::geometry::curve::{CurveEval, closest_point_on_curve};
use popkern_core::geometry::{NurbsCurve, Point3};

use super::junctions::shared_rails;
use crate::fillet_surface::FilletSurfaceResult;

/// Each bridge leg's minimum length past its junction's crossing, as
/// a fraction of the junction's setback scale (the farther of its
/// crossing's arclengths along the two rails). Measured on the four
/// PK-594 overrun apexes (pentagon h0.4 variable and wide, h0.8 wide,
/// hexagon h0.4 wide), each served by the tangent fit in 4–8 s:
///
/// - 0.25 and 0.5 serve all four, but leave 4 and 1 of the hexagon's
///   48 gated seams uncertified;
/// - every value in 0.75..=2.0 certifies G1 on every seam of all four.
///
/// 1.0 sits inside that window: each leg is as long as the crossing
/// distance it clears.
const BRIDGE_LEG_FRACTION: f64 = 1.0;

/// Polyline samples along a rail from its vertex station.
const RAIL_SAMPLES: usize = 512;

/// Each strip's feet at its station, in `order`: `[on the rail shared
/// with the previous strip, on the rail shared with the next]`.
/// `corners[k]` is the crossing of strip `order[k]`'s and strip
/// `order[k + 1]`'s rails.
pub(super) fn strip_station_feet(
    fillet_surfaces: &[FilletSurfaceResult],
    order: &[usize],
    corners: &[(f64, Point3)],
    vertex: Point3,
) -> Result<Vec<[Point3; 2]>, String> {
    let n = order.len();
    let rails = |a: usize, b: usize| {
        shared_rails(&fillet_surfaces[order[a]], &fillet_surfaces[order[b]])
            .ok_or_else(|| "adjacent strips have no coplanar rail pair".to_string())
    };
    let walks: Vec<[RailWalk<'_>; 2]> = (0..n)
        .map(|k| {
            let (prev, next) = ((k + n - 1) % n, (k + 1) % n);
            Ok([
                RailWalk::new(rails(prev, k)?.1, vertex),
                RailWalk::new(rails(k, next)?.0, vertex),
            ])
        })
        .collect::<Result<_, String>>()?;
    let behind = || "a rail crossing lies behind the strip's vertex station".to_string();
    // Junction k's crossing along strip k's rail (`along[k].0`) and
    // strip k+1's (`along[k].1`).
    let along: Vec<(f64, f64)> = (0..n)
        .map(|k| {
            let next = (k + 1) % n;
            Ok((
                walks[k][1].arclength_to(corners[k].1).ok_or_else(behind)?,
                walks[next][0]
                    .arclength_to(corners[k].1)
                    .ok_or_else(behind)?,
            ))
        })
        .collect::<Result<_, String>>()?;
    let station = balanced_stations(&along);
    (0..n)
        .map(|k| {
            let foot = |rail: &RailWalk<'_>| {
                rail.point_at_arclength(station[k]).ok_or_else(|| {
                    format!(
                        "strip {} rail ends before its station {:.4}",
                        order[k], station[k]
                    )
                })
            };
            Ok([foot(&walks[k][0])?, foot(&walks[k][1])?])
        })
        .collect()
}

/// One station per strip (arclength from each rail's vertex station)
/// whose bridge legs balance. Junction k's legs are
/// `station[k] − along[k].0` and `station[k+1] − along[k].1`; their
/// difference sums to a fixed total around the ring, so the stations
/// that minimise the squared differences give every junction the same
/// difference, the mean offset. Equal legs keep both feet convex cap
/// corners: each chord
/// makes the triangle angle `(180° − α) / 2` with its rails, for a face
/// apex angle α. A uniform shift then gives every leg at least
/// [`BRIDGE_LEG_FRACTION`] of its junction's setback scale, the
/// farther crossing.
fn balanced_stations(along: &[(f64, f64)]) -> Vec<f64> {
    let n = along.len();
    let offset = along.iter().map(|(a, b)| a - b).sum::<f64>() / n as f64;
    let mut station = vec![0.0_f64; n];
    for k in 0..n - 1 {
        station[k + 1] = station[k] - (along[k].0 - along[k].1) + offset;
    }
    let shift = (0..n)
        .map(|k| {
            let (a, b) = along[k];
            let legs = (station[k] - a).min(station[(k + 1) % n] - b);
            BRIDGE_LEG_FRACTION * a.max(b) - legs
        })
        .fold(f64::NEG_INFINITY, f64::max);
    station.iter().map(|s| s + shift).collect()
}

/// A rail walked from its vertex station into the strip. Trims run the
/// swept edge's length plus overshoot, so the strip lies toward the
/// domain midpoint.
struct RailWalk<'a> {
    rail: &'a NurbsCurve,
    t_vertex: f64,
    t_end: f64,
}

impl<'a> RailWalk<'a> {
    fn new(rail: &'a NurbsCurve, vertex: Point3) -> Self {
        let (t0, t1) = rail.parameter_range();
        let t_vertex = closest_point_on_curve(rail, &vertex).unwrap_or(t0);
        let t_end = if 0.5 * (t0 + t1) >= t_vertex {
            t0.max(t1)
        } else {
            t0.min(t1)
        };
        RailWalk {
            rail,
            t_vertex,
            t_end,
        }
    }

    /// Arclength from the vertex station to the rail's foot of `p`, or
    /// `None` when the foot is not inside the strip.
    fn arclength_to(&self, p: Point3) -> Option<f64> {
        let t = closest_point_on_curve(self.rail, &p)?;
        if (t - self.t_vertex) * (self.t_end - self.t_vertex) <= 0.0 {
            return None;
        }
        let mut acc = 0.0;
        let mut prev = self.rail.point_at(self.t_vertex);
        for i in 1..=RAIL_SAMPLES {
            let s = self.t_vertex + (t - self.t_vertex) * i as f64 / RAIL_SAMPLES as f64;
            let q = self.rail.point_at(s);
            acc += prev.distance_to(&q);
            prev = q;
        }
        Some(acc)
    }

    /// The rail point at `arclength` from the vertex station, or `None`
    /// when the rail ends first.
    fn point_at_arclength(&self, arclength: f64) -> Option<Point3> {
        let mut acc = 0.0;
        let mut t_prev = self.t_vertex;
        let mut prev = self.rail.point_at(t_prev);
        for i in 1..=RAIL_SAMPLES {
            let t = self.t_vertex + (self.t_end - self.t_vertex) * i as f64 / RAIL_SAMPLES as f64;
            let q = self.rail.point_at(t);
            let step = prev.distance_to(&q);
            if acc + step >= arclength {
                let f = if step > 0.0 {
                    (arclength - acc) / step
                } else {
                    0.0
                };
                // On the rail, not its chord.
                return Some(self.rail.point_at(t_prev + (t - t_prev) * f));
            }
            acc += step;
            t_prev = t;
            prev = q;
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A thin strip (index 1) between wide neighbours: its crossings
    /// lie far along its rails, theirs close to the vertex.
    const ALONG: [(f64, f64); 3] = [(0.04, 0.2), (0.22, 0.05), (0.06, 0.07)];

    #[test]
    fn balanced_stations_give_every_junction_the_same_leg_difference() {
        let station = balanced_stations(&ALONG);
        let offset = ALONG.iter().map(|(a, b)| a - b).sum::<f64>() / 3.0;
        for k in 0..3 {
            let (a, b) = ALONG[k];
            let (leg_a, leg_b) = (station[k] - a, station[(k + 1) % 3] - b);
            assert!(
                ((leg_a - leg_b) + offset).abs() < 1e-12,
                "junction {k}: legs {leg_a} {leg_b}"
            );
        }
    }

    #[test]
    fn balanced_stations_clear_every_crossing_by_the_leg_fraction() {
        let station = balanced_stations(&ALONG);
        let tight = (0..3)
            .map(|k| {
                let (a, b) = ALONG[k];
                let legs = (station[k] - a).min(station[(k + 1) % 3] - b);
                legs - BRIDGE_LEG_FRACTION * a.max(b)
            })
            .fold(f64::INFINITY, f64::min);
        assert!(
            tight.abs() < 1e-12,
            "the tightest junction sits on the bound: {tight}"
        );
    }
}
