//! Chart trim gate for the classify ray: cone (PK-669) and sphere
//! (PK-877) hits are tested against the face's own uv boundary.
//!
//! The boundary-derived cone gates (`axial_span_gate`) bound a hit by
//! the face's AXIAL extent only. That is the whole trim for a band
//! between two circles (a tip or mouth chamfer), but not for a cone
//! whose boundary varies with azimuth. A hex corner chamfer's cone
//! is trimmed by the six flats along hyperbolic arcs, and a ray
//! leaving the hex through a flat just above the chamfer edge then
//! also "hits" the cone's extension outside the hex, inside its
//! axial span (measured on the M6 nut: flat hit at (5.00, 2.42,
//! 0.47) and cone phantom at (5.01, 2.42, 0.47) — three crossings,
//! and a bore-void point read Inside).
//!
//! A sphere face had no trim test at all: its solver window is the
//! whole chart and the cap gate (`sphere_cap_gate`) abstains on a
//! hemisphere (its boundary directions are not all on one side). The
//! face's AABB prefilter hid that on an upright dome, whose rays from
//! the wall below the equator never enter the cap's box; turned 37°,
//! a ray from inside the dome's wall crossed the sphere's lower,
//! untrimmed half next to its real exit through the wall, the parity
//! read two, and a drill's wall fragment classified Outside (PK-877).
//!
//! The face's own uv boundary is the trim truth. The hit is tested
//! by exact even-odd parity over the period lifts of the boundary
//! loops (`periodic_chart_parity`). The gate abstains — keeps the
//! tuned behavior — whenever that truth is unavailable or untrusted:
//! a boundary half-edge without a pcurve, a pcurve whose image
//! leaves its edge's 3D curve, or (spheres) a loop that winds around
//! the chart without closing through a pole. A cone's chart has one
//! end (the apex, v = 0), so such a loop bounds the region on the
//! apex side, which is what the +v parity ray assumes. A sphere's chart
//! has a pole at each end, and its loops close THROUGH them: a step
//! along a pole row (a seam reaching the pole and the next leaving it
//! a period over) is a real chart segment, not a branch jump — a full
//! sphere's seam spine bounds the whole chart that way. Counted
//! literally, a v = π row adds one +v crossing per period and a v = 0
//! row none; a loop that still winds leaves which pole the region
//! holds undecided.
//!
//! A hit within the EDGE BAND of the face's 3D boundary (either side)
//! is a graze of the edge it shares with a neighbour: the ray crosses
//! the solid once there, but the two faces' trim tests cannot both be
//! exact so close to the edge — a plane neighbour's outer-loop gate
//! deliberately KEEPS hits inside its own sampling dead-band (0.2 %
//! of its uv extent), so the neighbour's extension hit survives and
//! the crossing counts twice (measured on the M6 nut: flat hit at
//! (1.57, −4.87, 5.00), 0.01 above its trim, and the cone hit 0.007
//! inside the hex). Such a ray is inconclusive; another direction
//! decides. The band mirrors the plane gate's: 0.2 % of the face's
//! boundary extent, floored at `EDGE_BAND_FLOOR_PER_LINEAR` linear
//! tolerances.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{Point3, pcurve::Pcurve};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle};

use super::gates::HitGate;

/// Polygon stations per boundary pcurve.
const STATIONS: usize = 64;
/// Edge band as a fraction of the face's boundary-sample extent (the
/// plane outer-loop gate's dead-band fraction).
const EDGE_BAND_FRACTION: f64 = 0.002;
/// Edge-band floor, in linear tolerances.
const EDGE_BAND_FLOOR_PER_LINEAR: f64 = 100.0;
/// Edge-band ceiling as a fraction of the boundary-sample extent.
const EDGE_BAND_CEILING_FRACTION: f64 = 0.05;
/// Image-honesty bound, in linear tolerances, between a boundary
/// pcurve's midpoint image and its edge's 3D curve.
const IMAGE_HONESTY_PER_LINEAR: f64 = 1.0e3;

/// Accept the hit when `(u, v)` lies on the cone or sphere face's trim,
/// skip it when it lies clearly off, poison the ray at the trim's
/// edge; [`HitGate::Accept`] for other faces and whenever the face's
/// uv boundary cannot be trusted. FINAL bodies only (caller-guarded).
pub(super) fn chart_trim_gate(
    brep: &BRep,
    face_handle: FaceHandle,
    surface: &SurfaceGeometry,
    uv: (f64, f64),
    hit: &Point3,
    tolerance: &Tolerance,
) -> HitGate {
    if !matches!(
        surface,
        SurfaceGeometry::Sphere { .. } | SurfaceGeometry::Cone { .. }
    ) {
        return HitGate::Accept;
    }
    chart_trim_verdict(brep, face_handle, surface, uv, hit, tolerance).unwrap_or(HitGate::Accept)
}

/// The trim verdict for a hit on a TORUS face whose chart cell cannot
/// be resolved (PK-431): a revolved groove's torus is full in u and
/// bounded in v by two circles and a seam, `face_chart_cell` declines
/// it, and the cell test poisoned every ray that crossed it — points
/// beside a bearing race's groove classified OnBoundary on all six
/// directions. Its loops do not wind in v, so the u-periodic parity
/// decides. `None` when the trim cannot be trusted (the caller keeps
/// poisoning the ray).
fn torus_trim_verdict(
    brep: &BRep,
    face_handle: FaceHandle,
    surface: &SurfaceGeometry,
    uv: (f64, f64),
    hit: &Point3,
    tolerance: &Tolerance,
) -> Option<HitGate> {
    matches!(surface, SurfaceGeometry::Torus { .. })
        .then(|| chart_trim_verdict(brep, face_handle, surface, uv, hit, tolerance))
        .flatten()
}

/// The trim gate for a hit on a torus face: fold-containment against
/// the face's chart cell when it resolves, else the chart-trim parity
/// ([`torus_trim_verdict`]); a hit neither can judge poisons the ray.
pub(super) fn torus_cell_gate(
    brep: &BRep,
    face_handle: FaceHandle,
    surface: &SurfaceGeometry,
    uv: (f64, f64),
    hit: &Point3,
    tolerance: &Tolerance,
) -> HitGate {
    match crate::clip::torus_hit_in_face_cell(brep, face_handle, uv.0, uv.1, tolerance) {
        Some(true) => HitGate::Accept,
        Some(false) => HitGate::SkipHit,
        None => torus_trim_verdict(brep, face_handle, surface, uv, hit, tolerance)
            .unwrap_or(HitGate::PoisonRay),
    }
}

/// The chart-trim verdict for a sphere, cone or torus hit (module
/// docs); `None` when the face's uv boundary cannot be trusted or does
/// not decide (an unclosed pole winding, a torus loop winding in v).
fn chart_trim_verdict(
    brep: &BRep,
    face_handle: FaceHandle,
    surface: &SurfaceGeometry,
    (u, v): (f64, f64),
    hit: &Point3,
    tolerance: &Tolerance,
) -> Option<HitGate> {
    let sphere = matches!(surface, SurfaceGeometry::Sphere { .. });
    let torus = matches!(surface, SurfaceGeometry::Torus { .. });
    let period = surface.parameter_periodicity().0?;
    let outer = crate::clip::face_uv_boundary(brep, face_handle)?;
    let inner = crate::clip::face_uv_inner_boundaries(brep, face_handle);
    let honesty = tolerance.linear() * IMAGE_HONESTY_PER_LINEAR;
    let mut rings: Vec<Vec<(f64, f64)>> = Vec::with_capacity(1 + inner.len());
    for contour in std::iter::once(&outer).chain(inner.iter()) {
        let mut ring = Vec::with_capacity(contour.len() * STATIONS);
        for (he, pcurve) in contour {
            if !image_honest(brep, *he, pcurve, surface, honesty) {
                return None;
            }
            let (t0, t1) = pcurve.parameter_range();
            // Endpoints included: a seam reaching a pole then puts a
            // sample ON the pole row (junction duplicates have zero
            // length and the parity skips them).
            for i in 0..=STATIONS {
                let q = pcurve.point_at(t0 + (t1 - t0) * i as f64 / STATIONS as f64);
                ring.push((q.x, q.y));
            }
        }
        rings.push(ring);
    }
    if rings.iter().map(Vec::len).sum::<usize>() < 3 {
        return None;
    }
    // Sphere rings close through the poles: a step along a pole row
    // is a real full-period chart segment (a full sphere's seam spine
    // is two meridians joined there), not a branch jump. At v = 0 it
    // lies below every +v parity ray; at v = π it crosses each once.
    let mut pole_crossings = 0usize;
    if sphere {
        for ring in &rings {
            pole_crossings += pole_row_closure(ring, period)?;
        }
    }
    // A torus is periodic in v too: the parity runs on the rings' own v
    // branch, which needs loops that do not wind in v.
    let mut v = v;
    if torus {
        let v_period = surface.parameter_periodicity().1?;
        if rings.iter().any(|r| winds(r.iter().map(|q| q.1), v_period)) {
            return None;
        }
        let n = rings.iter().map(Vec::len).sum::<usize>() as f64;
        let mean = rings.iter().flatten().map(|q| q.1).sum::<f64>() / n;
        v -= v_period * ((v - mean) / v_period).round();
    }
    let loops: Vec<&[(f64, f64)]> = rings.iter().map(Vec::as_slice).collect();
    let on_trim = crate::periodic_chart_parity::contains_on_u_periodic_chart(&loops, u, v, period)
        != (pole_crossings % 2 == 1);
    let boundary = brep.face_boundary_sample_points(face_handle, STATIONS);
    let (lo, hi) = boundary.iter().fold(
        (
            Point3::new(f64::MAX, f64::MAX, f64::MAX),
            Point3::new(f64::MIN, f64::MIN, f64::MIN),
        ),
        |(lo, hi), p| {
            (
                Point3::new(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z)),
                Point3::new(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z)),
            )
        },
    );
    let extent = (hi.x - lo.x) + (hi.y - lo.y) + (hi.z - lo.z);
    // Ceiling (PK-683): at a coarse tolerance the floor alone exceeds
    // the face (1e-2 x 100 = 1.0 on a radius-0.4 cone), every hit read
    // as an edge graze, and points inside the cone could not be
    // classified at all — the vertex-fan winding fallback had been
    // inventing their verdicts. A band past a few percent of the face
    // no longer separates edge grazes from interior crossings.
    let band = (EDGE_BAND_FRACTION * extent)
        .max(tolerance.linear() * EDGE_BAND_FLOOR_PER_LINEAR)
        .min(EDGE_BAND_CEILING_FRACTION * extent);
    let near_edge = boundary
        .windows(2)
        .any(|w| segment_distance(hit, &w[0], &w[1]) < band);
    if near_edge {
        if popkern_core::trace::trace_enabled("ray-classify") {
            eprintln!(
                "  [ray-band-chart] {face_handle:?} uv=({u:.4},{v:.4}) within {band:.2e} of the trim edge (on trim {on_trim}) — ray inconclusive"
            );
        }
        return Some(HitGate::PoisonRay);
    }
    if on_trim {
        return Some(HitGate::Accept);
    }
    if popkern_core::trace::trace_enabled("ray-classify") {
        eprintln!("  [ray-rej-chart] {face_handle:?} uv=({u:.4},{v:.4}) off the trim");
    }
    Some(HitGate::SkipHit)
}

/// Whether the closed sequence `values` winds round `period` (net
/// short-step travel of half a period or more).
fn winds(values: impl Iterator<Item = f64>, period: f64) -> bool {
    let vs: Vec<f64> = values.collect();
    let n = vs.len();
    let net: f64 = (0..n)
        .map(|i| {
            let d = vs[(i + 1) % n] - vs[i];
            d - period * (d / period).round()
        })
        .sum();
    net.abs() > 0.5 * period
}

/// Distance in v, in radians, within which a ring sample is on a pole.
const POLE_TOUCH: f64 = 1.0e-6;

/// The full periods a sphere ring's v = π pole-row steps traverse
/// (each crosses a +v parity ray once per period), or `None` when the
/// ring still winds around the chart with its pole-row steps taken
/// literally: it then separates the poles without closing through
/// either, and which one the region holds is undecided.
fn pole_row_closure(ring: &[(f64, f64)], period: f64) -> Option<usize> {
    let pi = std::f64::consts::PI;
    let row = |v: f64| {
        if v < POLE_TOUCH {
            Some(false)
        } else if v > pi - POLE_TOUCH {
            Some(true)
        } else {
            None
        }
    };
    let n = ring.len();
    let (mut net, mut high) = (0.0_f64, 0usize);
    for i in 0..n {
        let (a, b) = (ring[i], ring[(i + 1) % n]);
        let du = b.0 - a.0;
        match (row(a.1), row(b.1)) {
            (Some(x), Some(y)) if x == y => {
                net += du;
                if x {
                    high += (du.abs() / period).round() as usize;
                }
            }
            _ => net += du - period * (du / period).round(),
        }
    }
    (net.abs() <= 0.5 * period).then_some(high)
}

/// Whether `pcurve`'s midpoint image lies on its edge's 3D curve.
fn image_honest(
    brep: &BRep,
    he: HalfEdgeHandle,
    pcurve: &Pcurve,
    surface: &SurfaceGeometry,
    bound: f64,
) -> bool {
    let Some(curve) = brep
        .half_edges
        .get(he)
        .and_then(|h| brep.edges.get(h.edge))
        .and_then(|e| e.curve.as_ref())
    else {
        return false;
    };
    let (t0, t1) = pcurve.parameter_range();
    let uv = pcurve.point_at(0.5 * (t0 + t1));
    let image = surface.point_at(uv.x, uv.y);
    let eval: &dyn CurveEval = &curve.geometry;
    let Some(t) = eval.closest_point(&image) else {
        return false;
    };
    let (lo, hi) = (
        curve.t_start.min(curve.t_end),
        curve.t_start.max(curve.t_end),
    );
    eval.point_at(t.clamp(lo, hi)).distance_to(&image) <= bound
}

/// Distance from `p` to the segment `a`–`b`.
fn segment_distance(p: &Point3, a: &Point3, b: &Point3) -> f64 {
    let ab = *b - *a;
    let len2 = ab.dot(&ab);
    let s = if len2 > 0.0 {
        ((*p - *a).dot(&ab) / len2).clamp(0.0, 1.0)
    } else {
        0.0
    };
    (*a + ab * s).distance_to(p)
}
