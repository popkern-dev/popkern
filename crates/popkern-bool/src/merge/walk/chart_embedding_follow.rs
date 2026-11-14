//! Continuation of a projection walk through a chart's fast sweeps
//! (PK-1184).
//!
//! A projection walk continues each sample from its neighbour and
//! unwraps it by the nearest period image, which is faithful only while
//! consecutive samples lie under half a period apart. Near a sphere's
//! pole a curve sweeps longitude fast: a circle passing 1e-4 from the
//! pole and round it turns a whole period within a stretch as short as
//! that distance, and stepped station to station its run lost the turn
//! (a 17-sample run spanning 2.99 rad where the circle spans 2 pi), the
//! seam's borders found no window and the walk kept the ball's cap
//! inside out. A step turning a periodic coordinate past a quarter
//! period is halved instead, down to the curve's parameter resolution.

use popkern_core::geometry::Point3;

/// A projection continued from a seed point of the chart.
pub(super) type Project<'a> = dyn Fn(&Point3, (f64, f64)) -> Option<(f64, f64)> + 'a;

/// The most halvings one step may take.
const MAX_STEPS: usize = 4096;

/// The chart point of the curve at `tb`, continued from `seed` at `ta`
/// through steps that each turn a periodic coordinate (periods
/// `periods`) by at most a quarter period; `None` where a projection
/// fails.
pub(super) fn follow(
    point: &dyn Fn(f64) -> Point3,
    project: &Project<'_>,
    periods: (Option<f64>, Option<f64>),
    (ta, tb): (f64, f64),
    seed: (f64, f64),
) -> Option<(f64, f64)> {
    let turns = |a: (f64, f64), b: (f64, f64)| {
        periods.0.is_some_and(|p| (b.0 - a.0).abs() > 0.25 * p)
            || periods.1.is_some_and(|p| (b.1 - a.1).abs() > 0.25 * p)
    };
    let floor = 1e-12 * (tb - ta).abs().max(1e-300);
    let (mut uv, mut t) = (seed, ta);
    let mut targets = vec![tb];
    let mut steps = 0;
    while let Some(&target) = targets.last() {
        let next = project(&point(target), uv)?;
        steps += 1;
        if turns(uv, next) && (target - t).abs() > floor && steps < MAX_STEPS {
            targets.push(0.5 * (t + target));
            continue;
        }
        (uv, t) = (next, target);
        targets.pop();
    }
    Some(uv)
}
