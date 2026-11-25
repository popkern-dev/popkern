//! A run's direction read where it leaves its node's zone (tangency
//! design §4.8).
//!
//! The vertex table grows a vertex's zone to the reach its members lie
//! within: at a tangent crossing the records' own ends, a bending
//! length apart, are one vertex, and each record is re-anchored onto it
//! with a short hook from its own end. Read inside the zone, a run's
//! direction is the hook's: two records of a sphere crossing a fillet at
//! their touch, re-anchored 1.1e-3 rad along the chart, both departed
//! the vertex the wrong way round it, the wedge matching paired them
//! with the wrong neighbours, and three regions of the fillet walked as
//! one loop (generated corpus SL-V 614). Past twice the zone the run is
//! the curve's own.

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};

/// A chart point.
type Uv = (f64, f64);

/// The chart point where `uv` (from its first sample, which lies at
/// `node`) first leaves twice `zone` in space, or `None` when the whole
/// run stays inside. Bisected on the chart segment that crosses.
pub(super) fn zone_exit(
    surface: &SurfaceGeometry,
    uv: &[Uv],
    node: Point3,
    zone: f64,
) -> Option<Uv> {
    let radius = 2.0 * zone;
    let out = |q: Uv| surface.point_at(q.0, q.1).distance_to(&node) > radius;
    let k = (1..uv.len()).find(|&k| out(uv[k]))?;
    let (mut lo, mut hi) = (uv[k - 1], uv[k]);
    for _ in 0..40 {
        let mid = (0.5 * (lo.0 + hi.0), 0.5 * (lo.1 + hi.1));
        if out(mid) {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    Some(hi)
}

/// A refined sample this far off the zone exit's bearing (cosine) reads
/// a hook, not the chord error the sample corrects (≤ ~0.06 rad).
const HOOK_COS: f64 = 0.95;

/// Replace the run's endpoint-refined samples with their zone exits
/// where they read a hook (`zones` at the start and end node): inside the
/// zone and off the exit's bearing past [`HOOK_COS`]; returns which ends
/// were replaced. Every other refined sample keeps its precision: every
/// ordinary vertex grows a zone of a few tolerances, and reading all its
/// runs off the coarse run misordered the near-tangent wedges of a torus
/// less a slab and of the tolerance envelope's torus x box.
pub(super) fn read_past_zones(
    surface: &SurfaceGeometry,
    uv: &[Uv],
    ends: (Point3, Point3),
    zones: (Option<f64>, Option<f64>),
    refined: (&mut Option<Uv>, &mut Option<Uv>),
) -> (bool, bool) {
    let hook_exit = |run: &[Uv], q: Option<Uv>, node: Point3, z: Option<f64>| -> Option<Uv> {
        let (z, q, o) = (z?, q?, *run.first()?);
        if surface.point_at(q.0, q.1).distance_to(&node) >= 2.0 * z {
            return None;
        }
        let e = zone_exit(surface, run, node, z)?;
        let (a, b) = ((q.0 - o.0, q.1 - o.1), (e.0 - o.0, e.1 - o.1));
        let la = (a.0 * a.0 + a.1 * a.1).sqrt();
        let lb = (b.0 * b.0 + b.1 * b.1).sqrt();
        (la > 0.0 && lb > 0.0 && (a.0 * b.0 + a.1 * b.1) / (la * lb) < HOOK_COS).then_some(e)
    };
    let mut read = (false, false);
    if let Some(e) = hook_exit(uv, *refined.0, ends.0, zones.0) {
        *refined.0 = Some(e);
        read.0 = true;
    }
    let back: Vec<Uv> = uv.iter().rev().copied().collect();
    if let Some(e) = hook_exit(&back, *refined.1, ends.1, zones.1) {
        *refined.1 = Some(e);
        read.1 = true;
    }
    read
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::Vec3;

    fn xy_plane() -> SurfaceGeometry {
        SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        }
    }

    /// A run re-anchored onto its node hooks +u for 1e-3 before running
    /// -u: its first segment departs +u, its zone exit -u.
    #[test]
    fn a_hooked_run_departs_where_it_leaves_its_zone() {
        let run = [(0.0, 0.0), (1e-3, 1e-4), (-4e-3, -1e-3), (-2e-2, -5e-3)];
        let node = Point3::new(0.0, 0.0, 0.0);
        assert!(run[1].0 > 0.0, "the first segment reads the hook");
        let exit = zone_exit(&xy_plane(), &run, node, 1e-3).expect("leaves");
        assert!(exit.0 < 0.0, "the exit {exit:?} reads the run's course");
        let d = (exit.0 * exit.0 + exit.1 * exit.1).sqrt();
        assert!((d - 2e-3).abs() < 1e-9, "exit at twice the zone, {d:e}");
        // A run that never leaves keeps its refined samples.
        assert_eq!(zone_exit(&xy_plane(), &run[..2], node, 1e-3), None);
    }

    /// A refined sample on the hook is replaced; one inside the zone on
    /// the run's course keeps its precision.
    #[test]
    fn only_a_hooked_refined_sample_is_read_past_the_zone() {
        let run = [(0.0, 0.0), (1e-3, 1e-4), (-4e-3, -1e-3), (-2e-2, -5e-3)];
        let far = Point3::new(-2e-2, -5e-3, 0.0);
        let ends = (Point3::new(0.0, 0.0, 0.0), far);
        let (mut hooked, mut last) = (Some((5e-4, 5e-5)), None);
        let read = read_past_zones(
            &xy_plane(),
            &run,
            ends,
            (Some(1e-3), None),
            (&mut hooked, &mut last),
        );
        assert_eq!(read, (true, false));
        assert!(hooked.is_some_and(|q| q.0 < 0.0), "replaced: {hooked:?}");
        let mut on_course = Some((-1e-3, -2.5e-4));
        let read = read_past_zones(
            &xy_plane(),
            &run,
            ends,
            (Some(1e-3), None),
            (&mut on_course, &mut last),
        );
        assert_eq!(read, (false, false));
        assert_eq!(on_course, Some((-1e-3, -2.5e-4)));
    }
}
