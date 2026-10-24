//! Exact area of a trimmed cylinder, cone or sphere face (PK-935).
//!
//! On these surfaces the area element depends on v alone, `g(v) du dv`,
//! so by Green's theorem the area of a chart region is the boundary
//! integral `−∮ G(v) du` with `G' = g`. Three facts make it exact on
//! any trim, folded seams and pole visits included:
//! - Shifting a pcurve by a period leaves `∫ G(v) du` unchanged, so
//!   each half-edge's pcurve integrates as written, on whatever branch.
//! - Between one pcurve's end and the next one's start the loop runs
//!   along an iso-row in the chart: nothing at a regular vertex (the
//!   wrapped gap is zero), the pole row at a cone apex or sphere pole.
//!   The row contributes `−G(v)·Δu`.
//! - `G` is taken to vanish at a pole, so a cone's apex row costs
//!   nothing. A sphere's other pole row costs `G(π)·P = 4πr²`, the
//!   whole sphere, per period it runs: the region's area is exact
//!   modulo the sphere, and a face's area lies in (0, 4πr²].
//!
//! Each edge integrates by Gauss–Legendre on each knot span (on
//! several panels for analytic curves) along whichever of its pcurve
//! and its 3D curve sits nearer the true trim, the intersection with
//! the neighbouring face's surface: neither is always the better (a
//! plane-cut sphere's hole pcurve is a polyline whose chords sag 1e-7
//! where its circle is exact; a cone wall's notch pcurves are exact where its arcs sit
//! 1e-10 off). The area is exact to that trim's accuracy. Measured before: a notched cone wall and a 7/8
//! sphere face read their measuring mesh, 1.4e-4 relative off.

use crate::planar_area::traverses_forward;
use popkern_core::geometry::Point3;
use popkern_core::geometry::analytic_surfaces::analytic_signed_distance;
use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::uv_projection::project_point_to_uv_analytic;
use popkern_core::topology::{BRep, FaceHandle};
use popkern_math::PortableMath;

/// Gauss–Legendre nodes and weights on [−1, 1], 16 points.
const GL16: [(f64, f64); 8] = [
    (0.095_012_509_837_637_44, 0.189_450_610_455_068_5),
    (0.281_603_550_779_258_9, 0.182_603_415_044_923_6),
    (0.458_016_777_657_227_4, 0.169_156_519_395_002_5),
    (0.617_876_244_402_643_8, 0.149_595_988_816_576_7),
    (0.755_404_408_355_003, 0.124_628_971_255_533_9),
    (0.865_631_202_387_831_8, 0.095_158_511_682_492_78),
    (0.944_575_023_073_232_6, 0.062_253_523_938_647_89),
    (0.989_400_934_991_649_9, 0.027_152_459_411_754_09),
];
/// Panels per analytic (knotless) pcurve.
const ANALYTIC_PANELS: usize = 16;

/// The antiderivative `G` of a surface's area element, vanishing at the
/// chart's pole.
#[derive(Clone, Copy)]
enum Antiderivative {
    /// `r·v`.
    Cylinder(f64),
    /// `k·v²/2`, `k = tan α / cos α`.
    Cone(f64),
    /// `r²·(1 − cos v)`.
    Sphere(f64),
}

impl Antiderivative {
    fn of(surface: &SurfaceGeometry) -> Option<Self> {
        match *surface {
            SurfaceGeometry::Cylinder { radius, .. } => Some(Self::Cylinder(radius)),
            SurfaceGeometry::Cone { half_angle, .. } => {
                Some(Self::Cone(half_angle.ptan() / half_angle.pcos()))
            }
            SurfaceGeometry::Sphere { radius, .. } => Some(Self::Sphere(radius)),
            _ => None,
        }
    }

    fn at(self, v: f64) -> f64 {
        match self {
            Self::Cylinder(r) => r * v,
            Self::Cone(k) => 0.5 * k * v * v,
            Self::Sphere(r) => r * r * (1.0 - v.pcos()),
        }
    }

    /// The area modulus: the closed surface's area when its chart has a
    /// second pole, whose row costs that much per period.
    fn modulus(self) -> Option<f64> {
        match self {
            Self::Sphere(r) => Some(4.0 * std::f64::consts::PI * r * r),
            _ => None,
        }
    }
}

/// `du/dt` along `pc`. `Pcurve::tangent_at` is the derivative for the
/// analytic variants but the UNIT tangent for NURBS, so NURBS pcurves
/// take their first derivative directly.
fn du_dt(pc: &Pcurve, t: f64) -> f64 {
    match pc.as_nurbs() {
        Some(n) => n.derivative_at(t, 1).x,
        None => pc.tangent_at(t).x,
    }
}

/// `∫ f` over `[t0, t1]` (either order), by Gauss–Legendre on each
/// span between `knots` inside the interval, or on uniform panels when
/// there are none.
fn integrate(t0: f64, t1: f64, knots: Option<&[f64]>, f: impl Fn(f64) -> f64) -> f64 {
    let (lo, hi) = (t0.min(t1), t0.max(t1));
    let mut breaks: Vec<f64> = match knots {
        Some(k) => k.iter().copied().filter(|&k| k > lo && k < hi).collect(),
        None => (1..ANALYTIC_PANELS)
            .map(|i| lo + (hi - lo) * i as f64 / ANALYTIC_PANELS as f64)
            .collect(),
    };
    breaks.push(lo);
    breaks.push(hi);
    breaks.sort_by(f64::total_cmp);
    breaks.dedup();
    let sum: f64 = breaks
        .windows(2)
        .map(|w| {
            let (mid, half) = (0.5 * (w[0] + w[1]), 0.5 * (w[1] - w[0]));
            GL16.iter()
                .flat_map(|&(x, wt)| [(mid - half * x, wt), (mid + half * x, wt)])
                .map(|(t, wt)| wt * half * f(t))
                .sum::<f64>()
        })
        .sum();
    if t1 < t0 { -sum } else { sum }
}

/// `∫ G(v) du` along `pc`.
fn pcurve_integral(pc: &Pcurve, g: &dyn Fn(f64) -> f64) -> f64 {
    let (t0, t1) = pc.parameter_range();
    let knots = pc.as_nurbs().map(|n| n.knots.as_slice());
    integrate(t0, t1, knots, |t| g(pc.point_at(t).y) * du_dt(pc, t))
}

/// `∫ G(v) du` along the edge's 3D curve, traversed as `forward` says,
/// and its u travel: `v` from the surface's inverse of each point and
/// `du/ds = C′·S_u / |S_u|²` (the charts are orthogonal; the rate does
/// not depend on the u branch). `None` where the inverse fails or a
/// point sits on the pole.
fn edge_integral(
    surface: &SurfaceGeometry,
    edge: &CurveOnEdge,
    forward: bool,
    g: &dyn Fn(f64) -> f64,
) -> Option<(f64, f64)> {
    let (t0, t1) = if forward {
        (edge.t_start, edge.t_end)
    } else {
        (edge.t_end, edge.t_start)
    };
    let knots = match &edge.geometry {
        CurveGeometry::Nurbs(n) => Some(n.knots.as_slice()),
        _ => None,
    };
    let failed = std::cell::Cell::new(false);
    let rates = |s: f64| -> (f64, f64) {
        let p = edge.geometry.point_at(s);
        let Some((u, v)) = project_point_to_uv_analytic(surface, &p, None) else {
            failed.set(true);
            return (0.0, 0.0);
        };
        let su = surface.partial_u(u, v);
        let len2 = su.dot(&su);
        if len2 <= 0.0 {
            failed.set(true);
            return (0.0, 0.0);
        }
        (v, edge.geometry.derivative_at(s, 1).dot(&su) / len2)
    };
    let value = integrate(t0, t1, knots, |s| {
        let (v, du) = rates(s);
        g(v) * du
    });
    let travel = integrate(t0, t1, knots, |s| rates(s).1);
    (!failed.get()).then_some((value, travel))
}

/// Samples per representation when scoring a trim against the surfaces:
/// a prime count at half-step offsets, so they fall between the
/// vertices of a uniform polyline pcurve (a plane-cut sphere's hole
/// pcurve is a 1024-chord polyline; samples on its vertices scored it
/// 4e-13 where its chords sag 1e-7).
const MISFIT_SAMPLES: usize = 17;
/// Largest disagreement (radians) between an edge's u travel on its 3D
/// curve and on its pcurve for the curve to integrate the edge: a
/// pcurve over a pole runs the pole row inside the edge, which the 3D
/// curve cannot see (half a period).
const TRAVEL_AGREEMENT: f64 = 1e-6;

/// How far `pc`'s image and `edge` sit off the true trim, the
/// intersection of the face's surface `own` with the neighbour's
/// `other`: the pcurve's image lies on `own`, so it scores its distance
/// to `other`; the 3D curve scores its distance to either. `None` when
/// `other` has no implicit form.
fn misfits(
    pc: &Pcurve,
    edge: &CurveOnEdge,
    own: &SurfaceGeometry,
    other: &SurfaceGeometry,
) -> Option<(f64, f64)> {
    let off = |s: &SurfaceGeometry, p: &Point3| analytic_signed_distance(s, p).map(f64::abs);
    let (t0, t1) = pc.parameter_range();
    let mut on_pcurve: f64 = 0.0;
    let mut on_curve: f64 = 0.0;
    for k in 0..MISFIT_SAMPLES {
        let f = (k as f64 + 0.5) / MISFIT_SAMPLES as f64;
        let q = pc.point_at(t0 + (t1 - t0) * f);
        on_pcurve = on_pcurve.max(off(other, &own.point_at(q.x, q.y))?);
        let p = edge
            .geometry
            .point_at(edge.t_start + (edge.t_end - edge.t_start) * f);
        on_curve = on_curve.max(off(other, &p)?).max(off(own, &p)?);
    }
    Some((on_pcurve, on_curve))
}

/// The face's area by its boundary integral, or `None` when its surface
/// has no v-only area element (plane, torus, NURBS) or a half-edge has
/// no pcurve.
pub(crate) fn boundary_integral_area(brep: &BRep, face: FaceHandle) -> Option<f64> {
    let f = brep.faces.get(face)?;
    let surface = f.surface.as_ref()?;
    let antiderivative = Antiderivative::of(&surface.geometry)?;
    let g = |v: f64| antiderivative.at(v);
    let period = surface.geometry.parameter_periodicity().0?;
    let wrap = |d: f64| d - period * (d / period).round();
    let mut total = 0.0;
    for l in std::iter::once(f.outer_loop).chain(f.inner_loops.iter().copied()) {
        let start = brep.loops.get(l)?.half_edge;
        let mut uses = Vec::new();
        let mut h = start;
        for _ in 0..=brep.half_edges.len() {
            let he = brep.half_edges.get(h)?;
            let pc = brep.face_use_pcurve(face, h).or(he.pcurve.as_ref())?;
            uses.push((h, he, pc));
            h = he.next;
            if h == start {
                break;
            }
        }
        for (k, &(h, he, pc)) in uses.iter().enumerate() {
            // The pcurve, unless the 3D curve sits nearer the true trim.
            let on_curve = brep.edges.get(he.edge).and_then(|e| {
                let curve = e.curve.as_ref()?;
                let other = brep
                    .half_edges
                    .get(he.twin)
                    .and_then(|t| brep.loops.get(t.loop_))
                    .and_then(|tl| brep.faces.get(tl.face))
                    .and_then(|tf| tf.surface.as_ref())?;
                let (pcurve_off, curve_off) =
                    misfits(pc, curve, &surface.geometry, &other.geometry)?;
                if curve_off >= pcurve_off {
                    return None;
                }
                // The sense is the vertices', not the edge's owner: a
                // drilled sphere's split hole arcs ran their curves
                // backwards from their owning half-edge (PK-937).
                let from = brep.vertices.get(he.origin)?.position;
                let to = brep
                    .vertices
                    .get(brep.half_edges.get(he.next)?.origin)?
                    .position;
                let forward = traverses_forward(brep, h, curve, from, to);
                let (value, travel) = edge_integral(&surface.geometry, curve, forward, &g)?;
                let (t0, t1) = pc.parameter_range();
                let pcurve_travel = pc.point_at(t1).x - pc.point_at(t0).x;
                ((travel - pcurve_travel).abs() < TRAVEL_AGREEMENT).then_some(value)
            });
            total -= on_curve.unwrap_or_else(|| pcurve_integral(pc, &g));
            // The iso-row from this pcurve's end to the next one's start.
            let next = uses[(k + 1) % uses.len()].2;
            let end = pc.point_at(pc.parameter_range().1);
            let resume = next.point_at(next.parameter_range().0);
            total -= g(end.y) * wrap(resume.x - end.x);
        }
    }
    // The kernel flag: loops run counter-clockwise in the chart when
    // set, so the integral is the area's own sign.
    let signed = if surface.orientation { total } else { -total };
    let area = match antiderivative.modulus() {
        Some(m) => {
            let a = signed.rem_euclid(m);
            if a < 1e-12 * m { m } else { a }
        }
        None => signed.abs(),
    };
    area.is_finite().then_some(area)
}

#[cfg(test)]
#[path = "face_green_tests.rs"]
mod tests;
