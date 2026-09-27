//! Trim pcurves for the setback vertex blend's sweep and cap-internal
//! half-edges (PK-576): a degree-3 interpolant of the edge curve's
//! point-inverted samples, with the degree-1 polyline through the same
//! samples as the measured fallback.
//!
//! Why: the polyline is exact where the trim is straight in the
//! strip's chart, which is every constant-radius setback arc (an
//! iso-line). On a VARIABLE-radius strip the arc bends in uv, and the
//! chord sagitta over 32 segments put the strip-side trim 1.1e-6 off
//! the edge curve while the strip surface itself met the cap to 3.5e-8
//! (saddle_setback variable radius; cubic: 3.5e-8). Measured: seam G0
//! only. The cap-strip G1 there did not move (1.75e-5), because
//! it lives in the cap's own rim tangent, not the trim.
//!
//! Samples are point-inverted by CONTINUATION (each seeded from the
//! previous sample's uv), with the global projection kept only when it
//! is measurably closer. Global projection per sample jumped sheets on
//! the step's G2 strips (PK-583): near the saddle, the NURBS strip's
//! extension maps a second uv region onto the same rail, and one sample
//! landed there, so the trim ran 0.225 past the edge mid-span while
//! every sample was on the surface and both pcurve ends were exact.
//!
//! Fallback boundary: the cubic is kept only when its worst midpoint
//! distance to the edge curve does not exceed the polyline's. A
//! projection that jumped basins or a genuine kink in the samples
//! (where a cubic overshoots) keeps the polyline, which never
//! overshoots between samples.

use popkern_core::geometry::fitting::interpolate_curve;
use popkern_core::geometry::pcurve::{Pcurve, PcurveGeometry};
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::uv_projection::{project_point_to_uv, project_point_to_uv_near};
use popkern_core::geometry::{NurbsCurve, Point3};

/// The trim pcurve through `uv` (chart points, `z = 0`) sampled at
/// uniform fractions `i / (N − 1)` of the edge curve `curve_at`.
pub(super) fn sweep_trim_pcurve(
    uv: Vec<Point3>,
    surface: &SurfaceGeometry,
    curve_at: impl Fn(f64) -> Point3,
) -> Pcurve {
    let n = uv.len();
    let params: Vec<f64> = (0..n).map(|i| i as f64 / (n - 1).max(1) as f64).collect();
    let polyline = polyline_pcurve(uv.clone(), &params);
    let Ok(cubic) = interpolate_curve(&uv, &params, 3) else {
        return polyline;
    };
    let cubic = Pcurve {
        geometry: PcurveGeometry::Nurbs(cubic),
        t_start: 0.0,
        t_end: 1.0,
    };
    let miss = |pc: &Pcurve| -> f64 {
        (0..n.saturating_sub(1))
            .map(|i| {
                let f = (params[i] + params[i + 1]) * 0.5;
                let q = pc.point_at(f);
                surface.point_at(q.x, q.y).distance_to(&curve_at(f))
            })
            .fold(0.0, f64::max)
    };
    if miss(&cubic) <= miss(&polyline) {
        cubic
    } else {
        polyline
    }
}

fn polyline_pcurve(cps: Vec<Point3>, params: &[f64]) -> Pcurve {
    let mut knots = Vec::with_capacity(params.len() + 2);
    knots.push(0.0);
    knots.extend_from_slice(params);
    knots.push(1.0);
    let weights = vec![1.0; cps.len()];
    Pcurve {
        geometry: PcurveGeometry::Nurbs(NurbsCurve {
            degree: 1,
            knots,
            control_points: cps,
            weights,
        }),
        t_start: 0.0,
        t_end: 1.0,
    }
}

/// Chart images of `points` along a trim, by continuation: the first
/// sample projects globally, each later one from its predecessor's uv.
/// A seeded result stands unless the global one is closer by more than
/// `CONTINUATION_SLACK` (the seeded solve stalled). `None` when a
/// sample has no preimage.
pub(super) fn project_trim_samples(
    surface: &SurfaceGeometry,
    points: &[Point3],
) -> Option<Vec<Point3>> {
    const CONTINUATION_SLACK: f64 = 1.0e-9;
    let miss = |uv: (f64, f64), p: &Point3| surface.point_at(uv.0, uv.1).distance_to(p);
    let mut out: Vec<Point3> = Vec::with_capacity(points.len());
    for p in points {
        let global = project_point_to_uv(surface, p);
        let seeded = out
            .last()
            .and_then(|prev| project_point_to_uv_near(surface, p, prev.x, prev.y));
        let uv = match (seeded, global) {
            (Some(s), Some(g)) if miss(g, p) + CONTINUATION_SLACK < miss(s, p) => g,
            (Some(s), _) => s,
            (None, g) => g?,
        };
        out.push(Point3::new(uv.0, uv.1, 0.0));
    }
    Some(out)
}
