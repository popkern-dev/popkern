//! Sample-sequence fitting and closed-loop welding for marched SSI
//! traces (split from `surface_surface_fitting.rs` under the
//! 500-line file budget).

use crate::error::CoreError;
use crate::geometry::Point3;
use crate::geometry::curve_on_surface::CurveOnSurface;
use crate::geometry::fitting::fit_curve;
use crate::geometry::intersection::marching::MarchSample;

/// Fit NURBS curves (3D + UV1 + UV2) to a sequence of march samples.
///
/// The fit is a pure function of the samples. The MARCHING caller
/// welds its closed loops with [`weld_closed_loop`] before fitting
/// (its tip gap is a marching artifact); analytic arms (spiric
/// sections, cross-drill branches) sample their loops exactly and
/// must NOT be welded — the perturbed tail parameterization trips
/// their sub-tolerance fitted-curve gates (measured: the spiric
/// dispatch pin declined to marching and emitted fragment chains).
pub(crate) fn fit_march_samples(
    samples: &[MarchSample],
    is_closed: bool,
) -> Result<CurveOnSurface, CoreError> {
    if samples.len() < 2 {
        return Err(CoreError::PreconditionFailed(
            "Too few samples to fit".into(),
        ));
    }

    // Chord-length parameterization
    let mut params = vec![0.0];
    let mut total_length = 0.0;
    for i in 1..samples.len() {
        total_length += samples[i].point.distance_to(&samples[i - 1].point);
        params.push(total_length);
    }
    if total_length < 1e-30 {
        return Err(CoreError::PreconditionFailed(
            "Degenerate curve (zero length)".into(),
        ));
    }
    // Normalize to [0, 1]
    for p in &mut params {
        *p /= total_length;
    }

    let degree = 3.min(samples.len() - 1);
    let m = samples.len();

    // 3D curve, residual-adaptive: a fixed 1-in-5 control ratio (and
    // a hard 200-point cap) under-resolved tight curvature — a blade
    // girdle rounding a 6e-3-radius edge fitted 3.7e-4 off the true
    // intersection and the imprinted edge failed the curve-surface
    // band (PK-440). Double the control count until the worst sample
    // residual is far inside any serving band (1e-7 of the curve's
    // own scale), up to full interpolation.
    let pts_3d: Vec<Point3> = samples.iter().map(|s| s.point).collect();
    let target = 1.0e-7 * total_length.max(1.0);
    let mut n_cp = (m / 5).clamp(degree + 1, m);
    let curve_3d = loop {
        use crate::geometry::curve::CurveEval as _;
        let c = fit_curve(&pts_3d, &params, degree, n_cp)?;
        let worst = pts_3d
            .iter()
            .zip(&params)
            .map(|(p, &t)| p.distance_to(&c.point_at(t)))
            .fold(0.0_f64, f64::max);
        if worst <= target || n_cp >= m {
            break c;
        }
        // Past half the sample count, least squares leaves knot spans
        // holding barely one sample each: the clamped end spans are
        // then nearly unconstrained and the curve RINGS between its
        // first samples while every sample residual stays under the
        // target (measured, PK-673: the M3 nut's crest × flank
        // records at 204 control points for 258 exact samples bowed
        // 1.0e-4 off the flank inside their first interval — ten
        // times the curve-surface band — against 4.1e-7 for the
        // interpolant). The next rung is full interpolation, whose
        // knot averaging satisfies Schoenberg–Whitney.
        n_cp = if n_cp * 2 > m / 2 { m } else { n_cp * 2 };
    };

    // UV1 curve (using Point3(u1, v1, 0.0) convention)
    let pts_uv1: Vec<Point3> = samples
        .iter()
        .map(|s| Point3::new(s.u1, s.v1, 0.0))
        .collect();
    let curve_uv1 = fit_curve(&pts_uv1, &params, degree, n_cp)?;

    // UV2 curve
    let pts_uv2: Vec<Point3> = samples
        .iter()
        .map(|s| Point3::new(s.u2, s.v2, 0.0))
        .collect();
    let curve_uv2 = fit_curve(&pts_uv2, &params, degree, n_cp)?;

    Ok(CurveOnSurface {
        curve_3d,
        curve_uv1,
        curve_uv2,
        is_closed,
    })
}

/// Weld a marched closed loop: drop trailing samples within half an
/// average step of the start, then append the closing sample — the
/// exact start point, with the start's UVs when the loop also closes
/// in UV space (per surface, judged against the UV step scale) and
/// the tail's UVs otherwise (seam-winding loop: 3D closes, UV does
/// not).
pub(crate) fn weld_closed_loop(samples: &[MarchSample]) -> Vec<MarchSample> {
    let mut out = samples.to_vec();
    let first = out[0].clone();
    let avg_step = {
        let total: f64 = out
            .windows(2)
            .map(|w| w[0].point.distance_to(&w[1].point))
            .sum();
        (total / (out.len() - 1) as f64).max(1e-30)
    };
    while out.len() > 4 {
        let last = &out[out.len() - 1];
        if last.point.distance_to(&first.point) < avg_step * 0.5 {
            out.pop();
        } else {
            break;
        }
    }
    let tail = out[out.len() - 1].clone();
    let uv_gap = |a: (f64, f64), b: (f64, f64)| ((a.0 - b.0).powi(2) + (a.1 - b.1).powi(2)).sqrt();
    let uv_step = |pick: fn(&MarchSample) -> (f64, f64)| {
        let total: f64 = out
            .windows(2)
            .map(|w| uv_gap(pick(&w[0]), pick(&w[1])))
            .sum();
        (total / (out.len() - 1) as f64).max(1e-30)
    };
    let uv1_closes =
        uv_gap((tail.u1, tail.v1), (first.u1, first.v1)) < uv_step(|s| (s.u1, s.v1)) * 4.0;
    let uv2_closes =
        uv_gap((tail.u2, tail.v2), (first.u2, first.v2)) < uv_step(|s| (s.u2, s.v2)) * 4.0;
    out.push(MarchSample {
        point: first.point,
        u1: if uv1_closes { first.u1 } else { tail.u1 },
        v1: if uv1_closes { first.v1 } else { tail.v1 },
        u2: if uv2_closes { first.u2 } else { tail.u2 },
        v2: if uv2_closes { first.v2 } else { tail.v2 },
    });
    out
}

#[cfg(test)]
#[path = "surface_surface_fit_samples_tests.rs"]
mod tests;
