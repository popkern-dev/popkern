//! UV-loop containment for ray classification: pcurve-sampled
//! even-odd polygon tests for outer and inner face loops.
//!
//! Split from `classify.rs` to keep it within the structural line
//! ratchet.

use popkern_core::topology::{BRep, FaceHandle};
use popkern_math::PortableMath;

/// Polygon stations per CURVED pcurve in [`uv_point_in_loop`]: a full
/// circle at 64 stations has sagitta r·(1 − cos(π/64)) ≈ 0.0012·r,
/// the same density the outer-loop gate uses.
const CURVED_LOOP_STATIONS: usize = 64;

/// Sample a loop's pcurves into a UV polygon and even-odd test
/// `(u, v)` against it. Returns `None` when the loop has no usable
/// pcurve polygon (fewer than 3 samples) — the caller must treat
/// that as "no verdict", not "outside".
///
/// Linear pcurves are exact at their stations; CURVED pcurves are
/// sampled at [`CURVED_LOOP_STATIONS`] so the polygon's sagitta stays
/// a small fraction of the rim (PK-669: at 8 stations a bore circle
/// of radius 2.61 sampled to an octagon whose apothem is 2.41, so
/// ray hits up to 0.2 inside the M6 nut's bore read as hits ON the
/// holed end face, every ray rising through the bore void flipped
/// parity, and a tool fragment 0.01 below the nut classified Inside).
pub(crate) fn uv_point_in_loop(
    brep: &BRep,
    lh: popkern_core::topology::LoopHandle,
    u: f64,
    v: f64,
) -> Option<bool> {
    let loop_data = brep.loops.get(lh)?;
    let start = loop_data.half_edge;
    if !brep.half_edges.contains_key(start) {
        return None;
    }
    let mut ring: Vec<(f64, f64)> = Vec::new();
    let mut he = start;
    let cap = brep.half_edges.len() + 1;
    for _ in 0..cap {
        let d = &brep.half_edges[he];
        if let Some(pc) = d.pcurve.as_ref() {
            let (t0, t1) = pc.parameter_range();
            let linear_content = match &pc.geometry {
                popkern_core::geometry::pcurve::PcurveGeometry::Linear { .. } => true,
                popkern_core::geometry::pcurve::PcurveGeometry::Nurbs(c) => c.degree == 1,
                _ => false,
            };
            let stations = if linear_content {
                8
            } else {
                CURVED_LOOP_STATIONS
            };
            for i in 0..stations {
                let t = t0 + (t1 - t0) * i as f64 / stations as f64;
                let p = pc.point_at(t);
                ring.push((p.x, p.y));
            }
        }
        he = d.next;
        if he == start {
            break;
        }
    }
    if ring.len() < 3 {
        return None;
    }
    let mut crossings = 0;
    let n = ring.len();
    for i in 0..n {
        let (x0, y0) = ring[i];
        let (x1, y1) = ring[(i + 1) % n];
        if (y0 > v) != (y1 > v) {
            let xc = x0 + (v - y0) / (y1 - y0) * (x1 - x0);
            if xc > u {
                crossings += 1;
            }
        }
    }
    Some(crossings % 2 == 1)
}

/// True when the surface hit at `(u, v)` lands inside the face's
/// OUTER pcurve polygon. The verdict is trusted when the outer ring
/// is COMPLETE (every half-edge carries a pcurve). Linear pcurves
/// contribute their exact endpoints; CURVED pcurves are sampled
/// densely (64 stations — sagitta ≲ 0.2% of the span), and hits
/// within a dead-band of the polygon boundary abstain (keep the
/// hit) so sampling error can never reject a rim-grazing crossing
/// (measured: an 8-station circle misjudged a 2.7 mm band at r=35
/// and broke star-socket caps). A hit clearly outside — like the
/// 0.67-off phantom on an ellipse-fronted top face that flipped a
/// tool fragment to Inside — is rejected. Partial rings
/// mid-pipeline still abstain entirely.
///
/// Near the rim the decision is EXACT, not abstained: the curved
/// polygon's error is its measured chord sagitta, so outside twice that
/// (plus linear tolerance) polygon parity is decisive, and inside it the
/// hit is projected onto the true pcurve over the nearest segment's
/// window and its side read against the ring's winding. Only a hit
/// within linear tolerance of the rim abstains. The previous blanket
/// band of 0.2% of the face's span (0.036 on a 9-wide cap) accepted
/// rays crossing the cap plane up to 0.01 OUTSIDE a rod's rim: a point
/// 1e-3 below the rim read Inside the rod, a coplanar cap's outward
/// normal calibrated inverted, and a flush union dropped both copies of
/// the overlap (PK-750).
pub(crate) fn hit_inside_outer_loop(brep: &BRep, face: FaceHandle, u: f64, v: f64) -> bool {
    brep.faces
        .get(face)
        .and_then(|f| point_in_loop_exact(brep, f.outer_loop, u, v))
        .unwrap_or(true)
}

/// The exact-near-the-rim containment of [`hit_inside_outer_loop`] for
/// any loop `lh`: `None` when the ring is incomplete (a half-edge
/// without a pcurve) or has fewer than three segments — no verdict. A
/// point within linear tolerance of the rim reads inside.
pub(crate) fn point_in_loop_exact(
    brep: &BRep,
    lh: popkern_core::topology::LoopHandle,
    u: f64,
    v: f64,
) -> Option<bool> {
    let start = brep.loops.get(lh)?.half_edge;
    if !brep.half_edges.contains_key(start) {
        return None;
    }
    // Each half-edge's pcurve as polygon segments, each remembering the
    // pcurve window it chords (curved pcurves only).
    let mut segs: Vec<Seg<'_>> = Vec::new();
    let mut sagitta: f64 = 0.0;
    let mut he = start;
    let cap = brep.half_edges.len() + 1;
    for _ in 0..cap {
        let d = &brep.half_edges[he];
        let pc = d.pcurve.as_ref()?; // incomplete ring: no verdict
        let linear_content = match &pc.geometry {
            popkern_core::geometry::pcurve::PcurveGeometry::Linear { .. } => true,
            popkern_core::geometry::pcurve::PcurveGeometry::Nurbs(c) => c.degree == 1,
            _ => false,
        };
        let (t0, t1) = pc.parameter_range();
        if linear_content {
            let p0 = pc.point_at(t0);
            let p1 = pc.point_at(t1);
            segs.push(Seg {
                a: (p0.x, p0.y),
                b: (p1.x, p1.y),
                curve: None,
            });
        } else {
            // Curved pcurve: densify, measuring each chord's sagitta.
            let n = CURVED_LOOP_STATIONS;
            let mut prev = pc.point_at(t0);
            for i in 1..=n {
                let ta = t0 + (t1 - t0) * (i - 1) as f64 / n as f64;
                let tb = t0 + (t1 - t0) * i as f64 / n as f64;
                let q = pc.point_at(tb);
                let m = pc.point_at(0.5 * (ta + tb));
                let chord_mid = ((prev.x + q.x) * 0.5, (prev.y + q.y) * 0.5);
                sagitta = sagitta.max((m.x - chord_mid.0).phypot(m.y - chord_mid.1));
                segs.push(Seg {
                    a: (prev.x, prev.y),
                    b: (q.x, q.y),
                    curve: Some((pc, ta, tb)),
                });
                prev = q;
            }
        }
        he = d.next;
        if he == start {
            break;
        }
    }
    if segs.len() < 3 {
        return None;
    }
    // Chain the segments head to tail (a stored pcurve may run against
    // the traversal); `ring[k] -> ring[k + 1]` is segment k.
    let d2 = |a: (f64, f64), b: (f64, f64)| (a.0 - b.0).powi(2) + (a.1 - b.1).powi(2);
    let mut ring: Vec<(f64, f64)> = Vec::with_capacity(segs.len());
    let mut chained: Vec<&Seg<'_>> = Vec::with_capacity(segs.len());
    let mut tail = {
        let (s0, s1) = (&segs[0], &segs[1]);
        let fwd = d2(s0.b, s1.a).min(d2(s0.b, s1.b));
        let rev = d2(s0.a, s1.a).min(d2(s0.a, s1.b));
        ring.push(if fwd <= rev { s0.a } else { s0.b });
        chained.push(s0);
        if fwd <= rev { s0.b } else { s0.a }
    };
    for seg in segs.iter().skip(1) {
        if d2(tail, seg.a) <= d2(tail, seg.b) {
            ring.push(seg.a);
            tail = seg.b;
        } else {
            ring.push(seg.b);
            tail = seg.a;
        }
        chained.push(seg);
    }
    let n = ring.len();
    let crossings = (0..n)
        .filter(|&i| {
            let ((x0, y0), (x1, y1)) = (ring[i], ring[(i + 1) % n]);
            (y0 > v) != (y1 > v) && x0 + (v - y0) / (y1 - y0) * (x1 - x0) > u
        })
        .count();
    let polygon_inside = crossings % 2 == 1;
    if sagitta <= 0.0 {
        return Some(polygon_inside); // linear ring: exact
    }
    // Nearest segment; outside the sampling band the polygon decides.
    let tol = brep.tolerance.linear();
    let band = 2.0 * sagitta + tol;
    let (mut near, mut near_d) = (0, f64::MAX);
    for i in 0..n {
        let ((x0, y0), (x1, y1)) = (ring[i], ring[(i + 1) % n]);
        let (dx, dy) = (x1 - x0, y1 - y0);
        let l2 = (dx * dx + dy * dy).max(1e-300);
        let t = (((u - x0) * dx + (v - y0) * dy) / l2).clamp(0.0, 1.0);
        let d = (u - x0 - dx * t).phypot(v - y0 - dy * t);
        if d < near_d {
            (near, near_d) = (i, d);
        }
    }
    if near_d >= band {
        return Some(polygon_inside);
    }
    let Some((pc, ta, tb)) = chained[near].curve else {
        return Some(polygon_inside); // nearest is an exact linear piece
    };
    // Exact foot on the pcurve over the segment's window widened by one
    // segment each side (ternary search on a locally unimodal distance).
    let w = (tb - ta).abs();
    let (mut lo, mut hi) = (ta.min(tb) - w, ta.max(tb) + w);
    let (r0, r1) = pc.parameter_range();
    (lo, hi) = (lo.max(r0.min(r1)), hi.min(r0.max(r1)));
    let dist = |t: f64| {
        let p = pc.point_at(t);
        (p.x - u).phypot(p.y - v)
    };
    for _ in 0..60 {
        let (m1, m2) = (lo + (hi - lo) / 3.0, hi - (hi - lo) / 3.0);
        if dist(m1) <= dist(m2) {
            hi = m2;
        } else {
            lo = m1;
        }
    }
    let tf = 0.5 * (lo + hi);
    if dist(tf) <= tol {
        return Some(true); // on the rim within tolerance: abstain (keep)
    }
    // Side of the hit against the ring's traversal: the segment's
    // traversal direction (ring order) and the ring's winding.
    let foot = pc.point_at(tf);
    let (a, b) = (ring[near], ring[(near + 1) % n]);
    let along = (b.0 - a.0, b.1 - a.1);
    let left = along.0 * (v - foot.y) - along.1 * (u - foot.x) > 0.0;
    let area2: f64 = (0..n)
        .map(|i| {
            let ((x0, y0), (x1, y1)) = (ring[i], ring[(i + 1) % n]);
            x0 * y1 - x1 * y0
        })
        .sum();
    Some(left == (area2 > 0.0))
}

/// Exact REGION containment on a non-periodic chart: inside the outer
/// loop and outside every inner loop, each decided by
/// [`point_in_loop_exact`]. `None` when any ring gives no verdict.
pub(crate) fn point_in_face_region_exact(
    brep: &BRep,
    face: FaceHandle,
    u: f64,
    v: f64,
) -> Option<bool> {
    let f = brep.faces.get(face)?;
    if !point_in_loop_exact(brep, f.outer_loop, u, v)? {
        return Some(false);
    }
    for &lh in &f.inner_loops {
        if point_in_loop_exact(brep, lh, u, v)? {
            return Some(false);
        }
    }
    Some(true)
}

/// One polygon segment of a loop ring and, for a curved pcurve, the
/// pcurve window it chords.
struct Seg<'a> {
    a: (f64, f64),
    b: (f64, f64),
    curve: Option<(&'a popkern_core::geometry::pcurve::Pcurve, f64, f64)>,
}

/// True when the surface hit at `(u, v)` lands inside one of the
/// face's INNER-loop polygons (pcurve-sampled, even-odd).
pub(crate) fn hit_inside_inner_loop(brep: &BRep, face: FaceHandle, u: f64, v: f64) -> bool {
    let Some(face_data) = brep.faces.get(face) else {
        return false;
    };
    face_data
        .inner_loops
        .iter()
        .any(|&lh| uv_point_in_loop(brep, lh, u, v) == Some(true))
}

// ---------------------------------------------------------------------------
// 4. Point-in-solid with multi-ray voting
// ---------------------------------------------------------------------------

/// True when the hit at `(u, v)` lands strictly inside one of a
/// SINGLY-periodic face's holes (inner loops). Each hole contour is
/// unwrapped to continuity (a compact hole straddling the u seam sits on
/// both branches) and the query is moved to the contour's branch; a
/// contour still wider than half the period is torn and abstains, as do
/// hits within 20 × linear of a hole rim. Measured (PK-758): a web's
/// notch through a ring wall is a hole of the wall; ray hits inside it
/// counted as wall crossings and a sample inside the web read Outside.
pub(crate) fn hit_inside_periodic_hole(
    brep: &BRep,
    face: FaceHandle,
    period: f64,
    u: f64,
    v: f64,
) -> bool {
    let band = brep.tolerance.linear() * 20.0;
    crate::clip::face_uv_inner_boundaries(brep, face)
        .iter()
        .any(|contour| {
            let mut ring: Vec<(f64, f64)> = Vec::new();
            for (_he, pc) in contour {
                let (t0, t1) = pc.parameter_range();
                for i in 0..CURVED_LOOP_STATIONS {
                    let q = pc.point_at(t0 + (t1 - t0) * i as f64 / CURVED_LOOP_STATIONS as f64);
                    let x = match ring.last() {
                        Some(&(pu, _)) => q.x - ((q.x - pu) / period).round() * period,
                        None => q.x,
                    };
                    ring.push((x, q.y));
                }
            }
            let n = ring.len();
            if n < 3 {
                return false;
            }
            let (lo, hi) = ring
                .iter()
                .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), q| {
                    (lo.min(q.0), hi.max(q.0))
                });
            if hi - lo > 0.5 * period {
                return false;
            }
            let centre = 0.5 * (lo + hi);
            let u = u + ((centre - u) / period).round() * period;
            let mut inside = false;
            let mut near = f64::INFINITY;
            for i in 0..n {
                let ((x0, y0), (x1, y1)) = (ring[i], ring[(i + 1) % n]);
                if (y0 > v) != (y1 > v) && x0 + (v - y0) / (y1 - y0) * (x1 - x0) > u {
                    inside = !inside;
                }
                let (dx, dy) = (x1 - x0, y1 - y0);
                let l2 = (dx * dx + dy * dy).max(1e-300);
                let t = (((u - x0) * dx + (v - y0) * dy) / l2).clamp(0.0, 1.0);
                near = near.min((u - x0 - dx * t).phypot(v - y0 - dy * t));
            }
            inside && near > band
        })
}

/// [`hit_inside_periodic_hole`] on a singly-periodic chart, false on
/// any other surface.
pub(crate) fn hit_in_periodic_hole(
    brep: &BRep,
    face: FaceHandle,
    surface: &popkern_core::geometry::surface::SurfaceGeometry,
    (u, v): (f64, f64),
) -> bool {
    match crate::clip::classification_chart_periods(surface) {
        (Some(p), None) => hit_inside_periodic_hole(brep, face, p, u, v),
        _ => false,
    }
}
