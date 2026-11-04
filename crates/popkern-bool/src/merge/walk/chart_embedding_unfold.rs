//! Period-fold unfolding for pcurve-sourced chart embedding runs
//! (split from `chart_embedding.rs` under the 500-line file budget).

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};

/// PERIOD-FOLD UNFOLDING (PK-491): stored pcurves are normalized
/// into one canonical period branch, so a record crossing the
/// chart's period lattice mid-span carries a whole-period jump
/// between consecutive samples — chart-space discontinuity the
/// midpoint image check cannot see (a period shift maps to the
/// same 3D point). An embedding is a CONTINUOUS uv run by
/// contract; unfold each sample to the branch nearest its
/// predecessor, the same minimal-image rule `projection_run`
/// already enforces (measured: the LH two-start's tip-chamfer
/// runout record folded at u = 0, the anchored interior read
/// 7.75 > 2π, ChartOverwrap declined the cone to the 3D-Newell
/// fallback, and its ring parity stranded 24 record sides).
/// Samples are 1/16 of a record apart — well under half a period
/// — so the minimal image is faithful (the wrap detector's
/// established density argument). The δ-samples sit next to the
/// run's endpoints and snap onto the endpoint's branch so
/// directions stay coherent.
pub(super) fn unfold_pcurve_run(
    surface: &SurfaceGeometry,
    uv: &mut [(f64, f64)],
    refined_first: &mut Option<(f64, f64)>,
    refined_last: &mut Option<(f64, f64)>,
) {
    let (u_analytic, v_period) = surface.parameter_periodicity();
    let (u_chart, v_chart) = crate::clip::classification_chart_periods(surface);
    let u_period = u_analytic.or(u_chart);
    let v_period = v_period.or(v_chart);
    if u_period.is_none() && v_period.is_none() {
        return;
    }
    for i in 1..uv.len() {
        if let Some(p) = u_period {
            let du = uv[i].0 - uv[i - 1].0;
            uv[i].0 -= (du / p).round() * p;
        }
        if let Some(p) = v_period {
            let dv = uv[i].1 - uv[i - 1].1;
            uv[i].1 -= (dv / p).round() * p;
        }
    }
    if let Some(r) = refined_first
        && let Some(&(u0, v0)) = uv.first()
    {
        if let Some(p) = u_period {
            r.0 -= ((r.0 - u0) / p).round() * p;
        }
        if let Some(p) = v_period {
            r.1 -= ((r.1 - v0) / p).round() * p;
        }
    }
    if let Some(r) = refined_last
        && let Some(&(u1, v1)) = uv.last()
    {
        if let Some(p) = u_period {
            r.0 -= ((r.0 - u1) / p).round() * p;
        }
        if let Some(p) = v_period {
            r.1 -= ((r.1 - v1) / p).round() * p;
        }
    }
}

/// Bisection levels [`run_parameters`] may add to one uniform step.
const MAX_BISECTIONS: usize = 8;

/// Parameters of `pc` over `(s0, s1)` to sample a run at: `n` uniform
/// steps, each bisected while the step between its samples (to the
/// nearest period image) exceeds a quarter period. The unfold above
/// reads each step as the minimal image, which is faithful only for
/// steps under half a period; near a pole the longitude sweeps half a
/// turn over a short stretch of the curve (PK-1055: a loop 0.02 from
/// a sphere's pole swept more than pi within one sixteenth of it, the
/// unfold wrapped it the short way through the seam, the anchored
/// interior straddled the seam, and the face fell to the 3D walk,
/// which joined the two regions the loop separates).
pub(super) fn run_parameters(
    pc: &popkern_core::geometry::pcurve::Pcurve,
    surface: &SurfaceGeometry,
    (s0, s1): (f64, f64),
    n: usize,
) -> Vec<f64> {
    let (u_analytic, v_analytic) = surface.parameter_periodicity();
    let (u_chart, v_chart) = crate::clip::classification_chart_periods(surface);
    let periods = [u_analytic.or(u_chart), v_analytic.or(v_chart)];
    let wide = |a: f64, b: f64| -> bool {
        let (pa, pb) = (pc.point_at(a), pc.point_at(b));
        [(pb.x - pa.x, periods[0]), (pb.y - pa.y, periods[1])]
            .iter()
            .any(|&(d, p)| p.is_some_and(|p| (d - (d / p).round() * p).abs() > 0.25 * p))
    };
    let mut params = vec![s0];
    for i in 1..=n {
        let a = s0 + (s1 - s0) * (i - 1) as f64 / n as f64;
        let b = s0 + (s1 - s0) * i as f64 / n as f64;
        // Depth-first bisection, emitted in order.
        let mut stack = vec![(a, b, 0)];
        while let Some((lo, hi, depth)) = stack.pop() {
            if depth < MAX_BISECTIONS && wide(lo, hi) {
                let mid = 0.5 * (lo + hi);
                stack.push((mid, hi, depth + 1));
                stack.push((lo, mid, depth + 1));
            } else {
                params.push(hi);
            }
        }
    }
    params
}
