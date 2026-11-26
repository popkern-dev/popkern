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
        space_cos(surface, o, a, b)
            .is_some_and(|c| c < HOOK_COS)
            .then_some(e)
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

/// The cosine between the chords from `at` by chart steps `a` and `b`,
/// measured in space between the points they reach (not through the
/// first partials, which a pole collapses: a run 1e-5 beside a ball's
/// pole read its refined sample off its chord and lost its wedges);
/// `None` where either chord is a point.
fn space_cos(surface: &SurfaceGeometry, at: Uv, a: Uv, b: Uv) -> Option<f64> {
    let o = surface.point_at(at.0, at.1);
    let x = surface.point_at(at.0 + a.0, at.1 + a.1) - o;
    let y = surface.point_at(at.0 + b.0, at.1 + b.1) - o;
    let l = x.length() * y.length();
    (l > 0.0).then(|| x.dot(&y) / l)
}

/// A refined sample this far off its run's first chord in space
/// (cosine) is an artifact, never the chord error it corrects (≤ ~0.06
/// rad): the chart test in `ChartEmbedding::start_dir` alone is blind
/// in a stretched chart.
const REFINED_COS: f64 = 0.95;

/// Drop each endpoint-refined sample (not a zone exit) whose direction
/// strays off its run's first chord in space past [`REFINED_COS`], so
/// the run departs along its chord. A chart stretched along one
/// parameter folds directions toward it, and the chart's own cosine
/// passes what space does not (measured, YL-V 1189: on a 1 mm cylinder,
/// u an angle and v a length, the four records leaving a tangent
/// crossing hooked inside its touch band 1.6e-7 from the node; each
/// refined sample read 50 degrees off its chord in space, mirrored
/// across u, a cosine of 1 - 1e-7 in the chart, and the wedge matching
/// paired every arm with the wrong neighbour).
pub(super) fn drop_refined_off_course(
    surface: &SurfaceGeometry,
    uv: &[Uv],
    zone_ends: (bool, bool),
    refined: (&mut Option<Uv>, &mut Option<Uv>),
) {
    let off_course = |mut run: Box<dyn Iterator<Item = &Uv> + '_>, q: Uv| -> bool {
        let Some(&o) = run.next() else {
            return false;
        };
        let Some(&k) = run.find(|k| (k.0 - o.0).powi(2) + (k.1 - o.1).powi(2) > 1e-26) else {
            return false;
        };
        let (a, b) = ((q.0 - o.0, q.1 - o.1), (k.0 - o.0, k.1 - o.1));
        space_cos(surface, o, a, b).is_some_and(|c| c < REFINED_COS)
    };
    if !zone_ends.0
        && refined
            .0
            .is_some_and(|q| off_course(Box::new(uv.iter()), q))
    {
        *refined.0 = None;
    }
    if !zone_ends.1
        && refined
            .1
            .is_some_and(|q| off_course(Box::new(uv.iter().rev()), q))
    {
        *refined.1 = None;
    }
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

    /// On a 1 mm cylinder (u an angle, v a length) a refined sample
    /// mirrored across u from its run's chord is dropped, though the
    /// chart's own cosine passes it (the negative control); one along
    /// the chord is kept. The run (YL-V 1189) leaves at -27 degrees in
    /// space, the refined sample at +27.
    #[test]
    fn a_refined_sample_off_its_chord_in_space_is_dropped_in_a_stretched_chart() {
        let cylinder = SurfaceGeometry::Cylinder {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 1e-3,
        };
        let run = [(2.8, 2.3e-3), (2.9, 2.25e-3), (3.0, 2.2e-3)];
        let mirrored = (2.8 + 1e-4, 2.3e-3 + 5e-8);
        let chart_cos = {
            let (a, b): ((f64, f64), (f64, f64)) = ((1e-4, 5e-8), (0.1, -5e-5));
            (a.0 * b.0 + a.1 * b.1) / (a.0.hypot(a.1) * b.0.hypot(b.1))
        };
        assert!(chart_cos > REFINED_COS, "the chart passes it: {chart_cos}");
        let (mut first, mut last) = (Some(mirrored), Some((3.0 - 1e-4, 2.2e-3 + 5e-8)));
        drop_refined_off_course(&cylinder, &run, (false, false), (&mut first, &mut last));
        assert_eq!(first, None, "the mirrored sample is dropped");
        assert!(last.is_some(), "a sample along the chord is kept: {last:?}");
        // A zone exit is deliberate and kept whatever its bearing.
        let mut exit = Some(mirrored);
        drop_refined_off_course(&cylinder, &run, (true, false), (&mut exit, &mut None));
        assert_eq!(exit, Some(mirrored));
    }
}
