//! Rail pcurves on a NURBS blend strip for an exact-arc rail
//! (PK-417). A closed strip is a ring in `u`: the generic point
//! projection seeds from a fixed guess and converges to the nearest
//! local minimum, which for a rail whose arc starts a quarter turn
//! from the strip's `u = 0` lands 0.5–0.7 off the rail (measured on
//! the G2 boss plate once the boolean rim reached the fillet as an
//! `Arc`). The rail IS one of the strip's four boundary iso-curves, so
//! the pcurve is found by searching that iso-curve exhaustively: for
//! every arc sample, a dense scan of each boundary followed by a
//! golden-section refinement; the boundary with the smallest worst
//! deviation wins, provided it meets the gate.

use popkern_core::geometry::curve::{CurveEval, CurveGeometry};
use popkern_core::geometry::pcurve::{Pcurve, PcurveGeometry};
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{NurbsCurve, Point3};

const SAMPLES: usize = 65;
const SCAN: usize = 256;

/// One boundary iso-curve of a surface's parameter rectangle.
#[derive(Clone, Copy)]
struct Boundary {
    /// `true`: `u` varies, `v` fixed; `false`: `v` varies, `u` fixed.
    along_u: bool,
    fixed: f64,
    lo: f64,
    hi: f64,
}

impl Boundary {
    fn uv(&self, s: f64) -> (f64, f64) {
        if self.along_u {
            (s, self.fixed)
        } else {
            (self.fixed, s)
        }
    }
}

/// Nearest parameter on `boundary` to `target`, by dense scan then
/// golden-section refinement; returns `(s, distance)`.
fn nearest_on_boundary(
    surface: &SurfaceGeometry,
    boundary: &Boundary,
    target: &Point3,
) -> (f64, f64) {
    let distance = |s: f64| {
        let (u, v) = boundary.uv(s);
        surface.point_at(u, v).distance_to(target)
    };
    let mut best = (boundary.lo, f64::INFINITY);
    for i in 0..=SCAN {
        let s = boundary.lo + (boundary.hi - boundary.lo) * i as f64 / SCAN as f64;
        let d = distance(s);
        if d < best.1 {
            best = (s, d);
        }
    }
    let step = (boundary.hi - boundary.lo) / SCAN as f64;
    let (mut a, mut b) = (
        (best.0 - step).max(boundary.lo),
        (best.0 + step).min(boundary.hi),
    );
    let phi = 0.618_033_988_749_895;
    for _ in 0..48 {
        let c = b - phi * (b - a);
        let d = a + phi * (b - a);
        if distance(c) < distance(d) {
            b = d;
        } else {
            a = c;
        }
    }
    let s = 0.5 * (a + b);
    (s, distance(s))
}

/// The co-parameterised pcurve of `arc` over `range` on the strip
/// boundary it lies on, or `None` when no boundary carries it within
/// `gate`.
pub(crate) fn rail_pcurve_on_strip(
    surface: &SurfaceGeometry,
    arc: &CurveGeometry,
    range: (f64, f64),
    gate: f64,
) -> Option<Pcurve> {
    let ((u0, u1), (v0, v1)) = surface.parameter_range();
    if !(u0.is_finite() && u1.is_finite() && v0.is_finite() && v1.is_finite()) {
        return None;
    }
    let boundaries = [
        Boundary {
            along_u: true,
            fixed: v0,
            lo: u0,
            hi: u1,
        },
        Boundary {
            along_u: true,
            fixed: v1,
            lo: u0,
            hi: u1,
        },
        Boundary {
            along_u: false,
            fixed: u0,
            lo: v0,
            hi: v1,
        },
        Boundary {
            along_u: false,
            fixed: u1,
            lo: v0,
            hi: v1,
        },
    ];
    let params: Vec<f64> = (0..SAMPLES)
        .map(|i| range.0 + (range.1 - range.0) * i as f64 / (SAMPLES - 1) as f64)
        .collect();
    let targets: Vec<Point3> = params.iter().map(|&t| arc.point_at(t)).collect();
    let (period_u, period_v) = surface.periodicity();
    let (chart_u, chart_v) = surface.parameter_periodicity();
    let mut best: Option<(f64, Vec<(f64, f64)>)> = None;
    for boundary in &boundaries {
        let mut worst: f64 = 0.0;
        let mut uvs: Vec<(f64, f64)> = Vec::with_capacity(SAMPLES);
        // A closed rail on a periodic strip crosses the strip's seam:
        // the principal parameter wraps there, so consecutive samples
        // are unwrapped by the period — the pcurve runs continuously
        // past the chart edge (the periodic-pcurve convention every
        // seam-crossing rim uses), never jumping back. A NURBS ring
        // that reports no period but whose boundary closes on itself
        // is unwrapped by its parameter span.
        let declared = if boundary.along_u {
            period_u.or(chart_u)
        } else {
            period_v.or(chart_v)
        };
        let closes = {
            let (u_a, v_a) = boundary.uv(boundary.lo);
            let (u_b, v_b) = boundary.uv(boundary.hi);
            surface
                .point_at(u_a, v_a)
                .distance_to(&surface.point_at(u_b, v_b))
                <= gate
        };
        let period = declared.or_else(|| closes.then_some(boundary.hi - boundary.lo));
        for target in &targets {
            let (mut s, d) = nearest_on_boundary(surface, boundary, target);
            worst = worst.max(d);
            if let (Some(period), Some(prev)) = (period, uvs.last()) {
                let prev_s = if boundary.along_u { prev.0 } else { prev.1 };
                s -= period * ((s - prev_s) / period).round();
            }
            uvs.push(boundary.uv(s));
        }
        if popkern_core::trace::trace_enabled("fillet") {
            let show = |i: usize| format!("({:.3},{:.3})", uvs[i].0, uvs[i].1);
            eprintln!(
                "  [closed-fillet] rail on strip boundary along_u={} fixed={:.3} lo={:.3} hi={:.3} period={period:?} closes={closes}: worst {worst:.3e} (gate {gate:.1e}) uv {} {} {} {}",
                boundary.along_u,
                boundary.fixed,
                boundary.lo,
                boundary.hi,
                show(0),
                show(SAMPLES / 3),
                show(2 * SAMPLES / 3),
                show(SAMPLES - 1)
            );
        }
        if best.as_ref().is_none_or(|(w, _)| worst < *w) {
            best = Some((worst, uvs));
        }
    }
    let (worst, uvs) = best?;
    if worst > gate {
        return None;
    }
    let mut knots = Vec::with_capacity(SAMPLES + 2);
    knots.push(params[0]);
    knots.extend(params.iter().copied());
    knots.push(params[SAMPLES - 1]);
    Some(Pcurve {
        geometry: PcurveGeometry::Nurbs(NurbsCurve {
            degree: 1,
            knots,
            control_points: uvs
                .into_iter()
                .map(|(u, v)| Point3::new(u, v, 0.0))
                .collect(),
            weights: vec![1.0; SAMPLES],
        }),
        t_start: range.0,
        t_end: range.1,
    })
}
