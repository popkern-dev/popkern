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
