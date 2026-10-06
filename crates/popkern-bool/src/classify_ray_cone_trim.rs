//! CONE trim gate for the classify ray (PK-669).
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
//! The face's own uv boundary is the trim truth. The hit is tested
//! by exact even-odd parity over the period lifts of the boundary
//! loops (`periodic_chart_parity`). The gate abstains — keeps the
//! tuned behavior — whenever that truth is unavailable or untrusted:
//! a boundary half-edge without a pcurve, or a pcurve whose image
//! leaves its edge's 3D curve.
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

/// Accept the hit when `(u, v)` lies on the cone face's trim, skip it
/// when it lies clearly off, poison the ray at the trim's edge;
/// [`HitGate::Accept`] for non-cone faces and whenever the face's uv
/// boundary cannot be trusted. FINAL bodies only (caller-guarded).
pub(super) fn cone_trim_gate(
    brep: &BRep,
    face_handle: FaceHandle,
    surface: &SurfaceGeometry,
    (u, v): (f64, f64),
    hit: &Point3,
    tolerance: &Tolerance,
) -> HitGate {
    if !matches!(surface, SurfaceGeometry::Cone { .. }) {
        return HitGate::Accept;
    }
    let Some(period) = surface.parameter_periodicity().0 else {
        return HitGate::Accept;
    };
    let Some(outer) = crate::clip::face_uv_boundary(brep, face_handle) else {
        return HitGate::Accept;
    };
    let inner = crate::clip::face_uv_inner_boundaries(brep, face_handle);
    let honesty = tolerance.linear() * IMAGE_HONESTY_PER_LINEAR;
    let mut rings: Vec<Vec<(f64, f64)>> = Vec::with_capacity(1 + inner.len());
    for contour in std::iter::once(&outer).chain(inner.iter()) {
        let mut ring = Vec::with_capacity(contour.len() * STATIONS);
        for (he, pcurve) in contour {
            if !image_honest(brep, *he, pcurve, surface, honesty) {
                return HitGate::Accept;
            }
            let (t0, t1) = pcurve.parameter_range();
            for i in 0..STATIONS {
                let q = pcurve.point_at(t0 + (t1 - t0) * i as f64 / STATIONS as f64);
                ring.push((q.x, q.y));
            }
        }
        rings.push(ring);
    }
    if rings.iter().map(Vec::len).sum::<usize>() < 3 {
        return HitGate::Accept;
    }
    let loops: Vec<&[(f64, f64)]> = rings.iter().map(Vec::as_slice).collect();
    let on_trim = crate::periodic_chart_parity::contains_on_u_periodic_chart(&loops, u, v, period);
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
                "  [ray-band-cone] {face_handle:?} uv=({u:.4},{v:.4}) within {band:.2e} of the trim edge (on trim {on_trim}) — ray inconclusive"
            );
        }
        return HitGate::PoisonRay;
    }
    if on_trim {
        return HitGate::Accept;
    }
    if popkern_core::trace::trace_enabled("ray-classify") {
        eprintln!("  [ray-rej-cone] {face_handle:?} uv=({u:.4},{v:.4}) off the trim");
    }
    HitGate::SkipHit
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
