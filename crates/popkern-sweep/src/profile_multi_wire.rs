//! Multi-wire and face-loop profile analysis (split from
//! `profile.rs` under the 500-line budget): outer/inner wire roles,
//! per-role winding normalization, curve-faithful loop adoption
//! (PK-477, PK-479).

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, LoopHandle};

use super::{
    AnalyzedProfile, PROFILE_SAMPLES, analyze_wire, check_closed, compute_centroid,
    detect_and_split, reversed_curve, signed_area_2d, translate_curve,
};
use crate::error::SweepError;

/// Walk a loop's half-edges into an ordered, traversal-oriented span
/// chain of NURBS curves.
pub(super) fn loop_to_wire_chain(
    brep: &BRep,
    loop_handle: LoopHandle,
) -> Result<Vec<NurbsCurve>, SweepError> {
    let loop_data = brep
        .loops
        .get(loop_handle)
        .ok_or_else(|| SweepError::InvalidProfile("loop handle not found in BRep".into()))?;
    let start_he = loop_data.half_edge;
    let mut spans = Vec::new();
    let mut current = start_he;
    loop {
        let he = brep.half_edges.get(current).ok_or_else(|| {
            SweepError::InvalidProfile("half-edge handle not found in BRep".into())
        })?;
        let twin = brep
            .half_edges
            .get(he.twin)
            .ok_or_else(|| SweepError::InvalidProfile("twin handle not found in BRep".into()))?;
        let origin = brep
            .vertices
            .get(he.origin)
            .ok_or_else(|| SweepError::InvalidProfile("vertex handle not found in BRep".into()))?
            .position;
        let end = brep
            .vertices
            .get(twin.origin)
            .ok_or_else(|| SweepError::InvalidProfile("vertex handle not found in BRep".into()))?
            .position;
        let edge = brep
            .edges
            .get(he.edge)
            .ok_or_else(|| SweepError::InvalidProfile("edge handle not found in BRep".into()))?;
        let span = match &edge.curve {
            Some(c) => {
                // Normalize the adopted span onto the [0,1] knot
                // domain the sweep machinery assumes everywhere
                // (surface corners, cap curves, seam pcurves are all
                // evaluated at 0 and 1). An adopted rim curve carries
                // its SOURCE domain — an arc's angles — and u = 0.0
                // then lands mid-wire: the ring vertices sat 90° off
                // every cap curve start and both caps gaped by the
                // quarter chord (measured: gap 1.414 on an r = 1
                // drilled rim).
                let nurbs = normalized_domain(c.geometry.to_nurbs_in_range((c.t_start, c.t_end)));
                let (t0, t1) = nurbs.parameter_range();
                let s = nurbs.point_at(t0);
                let e = nurbs.point_at(t1);
                let forward = if origin.distance_to(&end) > brep.tolerance.linear() {
                    // Open span: orient by endpoint matching.
                    s.distance_to(&origin) + e.distance_to(&end)
                        <= s.distance_to(&end) + e.distance_to(&origin)
                } else {
                    // Closed span (a full circle rim): endpoints
                    // cannot orient it; the canonical half-edge runs
                    // the stored curve forward (PK-303 convention).
                    edge.half_edge == current
                };
                if forward {
                    nurbs
                } else {
                    reversed_curve(&nurbs)
                }
            }
            None => NurbsCurve {
                degree: 1,
                knots: vec![0.0, 0.0, 1.0, 1.0],
                control_points: vec![origin, end],
                weights: vec![1.0, 1.0],
            },
        };
        spans.push(span);
        current = he.next;
        if current == start_he {
            break;
        }
    }
    if spans.is_empty() {
        return Err(SweepError::InvalidProfile("loop has no half-edges".into()));
    }
    Ok(spans)
}

/// Assemble an [`AnalyzedProfile`] from per-loop span chains: the
/// largest-area chain becomes the outer wire (ordered first), every
/// chain must close, windings are normalized per role (outer CCW,
/// inner CW), and inner chains must lie inside the outer — the same
/// contract as [`analyze_multi_wire`], fed exact spans instead of
/// vertex polylines.
pub(super) fn analyze_wire_chains(
    chains: Vec<Vec<NurbsCurve>>,
    tolerance: &Tolerance,
) -> Result<AnalyzedProfile, SweepError> {
    let chain_points = |chain: &[NurbsCurve]| -> Vec<Point3> {
        let per_span = (PROFILE_SAMPLES / chain.len().max(1)).max(4);
        let mut pts = Vec::with_capacity(per_span * chain.len());
        for span in chain {
            let (t0, t1) = span.parameter_range();
            for i in 0..per_span {
                let t = t0 + (i as f64 / per_span as f64) * (t1 - t0);
                pts.push(span.point_at(t));
            }
        }
        pts
    };
    let chain_area = |chain: &[NurbsCurve]| -> f64 {
        let pts = chain_points(chain);
        let n = pts.len();
        let mut area = 0.0;
        for i in 0..n {
            let (a, b) = (pts[i], pts[(i + 1) % n]);
            area += a.x * b.y - b.x * a.y;
        }
        0.5 * area
    };
    let chain_closed = |chain: &[NurbsCurve]| -> bool {
        let first = &chain[0];
        let last = &chain[chain.len() - 1];
        let (t0, _) = first.parameter_range();
        let (_, t1) = last.parameter_range();
        first.point_at(t0).distance_to(&last.point_at(t1)) < tolerance.linear() * 10.0
    };

    let mut outer_idx = 0;
    let mut best = 0.0_f64;
    for (i, chain) in chains.iter().enumerate() {
        let a = chain_area(chain).abs();
        if a > best {
            best = a;
            outer_idx = i;
        }
    }
    let mut order: Vec<usize> = (0..chains.len()).collect();
    order.retain(|&i| i != outer_idx);
    order.insert(0, outer_idx);

    // Centroid of the outer chain (uniform samples).
    let outer_pts = chain_points(&chains[outer_idx]);
    let inv = 1.0 / outer_pts.len().max(1) as f64;
    let centroid = Point3::new(
        outer_pts.iter().map(|p| p.x).sum::<f64>() * inv,
        outer_pts.iter().map(|p| p.y).sum::<f64>() * inv,
        outer_pts.iter().map(|p| p.z).sum::<f64>() * inv,
    );
    let outer_poly: Vec<(f64, f64)> = outer_pts
        .iter()
        .map(|p| (p.x - centroid.x, p.y - centroid.y))
        .collect();

    let mut spans = Vec::new();
    let mut wires = Vec::with_capacity(chains.len());
    for (rank, &ci) in order.iter().enumerate() {
        let chain: Vec<NurbsCurve> = chains[ci]
            .iter()
            .map(|c| translate_curve(c, centroid))
            .collect();
        if !chain_closed(&chain) {
            return Err(SweepError::InvalidProfile(format!(
                "face profile: boundary loop {ci} does not close"
            )));
        }
        let is_outer = rank == 0;
        if !is_outer {
            let pts = chain_points(&chain);
            let inv = 1.0 / pts.len().max(1) as f64;
            let c = (
                pts.iter().map(|p| p.x).sum::<f64>() * inv,
                pts.iter().map(|p| p.y).sum::<f64>() * inv,
            );
            if !point_in_polygon(&outer_poly, c) {
                return Err(SweepError::InvalidProfile(format!(
                    "face profile: inner loop {ci} lies outside the outer boundary"
                )));
            }
        }
        let area = chain_area(&chain);
        let wants_ccw = is_outer;
        let chain: Vec<NurbsCurve> = if (area < 0.0) == wants_ccw {
            chain.iter().rev().map(reversed_curve).collect()
        } else {
            chain
        };
        wires.push(chain.len());
        spans.extend(chain);
    }

    Ok(AnalyzedProfile {
        spans,
        wires,
        is_closed: true,
        centroid,
    })
}

/// Analyze multiple wires: the largest (by area) is the OUTER wire
/// and is ordered first; the rest are inner (hole) wires.
///
/// Contract (PK-477): every wire must be closed — a mixed open/closed
/// multi-wire has no solid or sheet interpretation and declines as
/// `InvalidProfile`. Winding is normalized per role (outer CCW, inner
/// CW in the profile plane) so `du×dv` faces out of the material on
/// every side face, hole walls included. Each inner wire's centroid
/// must lie inside the outer wire's polygon — a disjoint "inner" wire
/// would sweep a second, cross-wired solid and declines instead.
pub(super) fn analyze_multi_wire(
    curves: &[NurbsCurve],
    tolerance: &Tolerance,
) -> Result<AnalyzedProfile, SweepError> {
    if curves.is_empty() {
        return Err(SweepError::InvalidProfile("empty multi-wire".into()));
    }
    if curves.len() == 1 {
        return analyze_wire(&curves[0], tolerance);
    }

    // The outer wire is the one with the largest absolute area.
    let mut outer_idx = 0;
    let mut best_area = 0.0_f64;
    for (i, curve) in curves.iter().enumerate() {
        let area = signed_area_2d(curve).abs();
        if area > best_area {
            best_area = area;
            outer_idx = i;
        }
    }

    let centroid = compute_centroid(&curves[outer_idx]);
    let mut order: Vec<usize> = (0..curves.len()).collect();
    order.retain(|&i| i != outer_idx);
    order.insert(0, outer_idx);

    // Outer polygon (centered frame) for the containment check.
    let outer_centered = translate_curve(&curves[outer_idx], centroid);
    let outer_poly: Vec<(f64, f64)> = {
        let (t0, t1) = outer_centered.parameter_range();
        (0..PROFILE_SAMPLES)
            .map(|i| {
                let t = t0 + (i as f64 / PROFILE_SAMPLES as f64) * (t1 - t0);
                let p = outer_centered.point_at(t);
                (p.x, p.y)
            })
            .collect()
    };

    let mut all_spans = Vec::new();
    let mut wires = Vec::with_capacity(curves.len());
    for (rank, &i) in order.iter().enumerate() {
        let centered = translate_curve(&curves[i], centroid);
        if !check_closed(&centered, tolerance) {
            return Err(SweepError::InvalidProfile(format!(
                "multi-wire profile: wire {i} is not closed (every wire \
                 of a multi-wire profile must be a closed loop)"
            )));
        }
        let is_outer = rank == 0;
        if !is_outer {
            let c = compute_centroid(&centered);
            if !point_in_polygon(&outer_poly, (c.x, c.y)) {
                return Err(SweepError::InvalidProfile(format!(
                    "multi-wire profile: wire {i} lies outside the outer \
                     wire (inner wires must be holes of the outer)"
                )));
            }
        }
        let area = signed_area_2d(&centered);
        let wants_ccw = is_outer;
        let centered = if (area < 0.0) == wants_ccw {
            reversed_curve(&centered)
        } else {
            centered
        };
        let spans = detect_and_split(&centered)?;
        wires.push(spans.len());
        all_spans.extend(spans.into_iter().map(normalized_domain));
    }

    Ok(AnalyzedProfile {
        spans: all_spans,
        wires,
        is_closed: true,
        centroid,
    })
}

/// Affinely rescale a curve's knot vector onto [0, 1]; geometry is
/// untouched.
pub(super) fn normalized_domain(mut curve: NurbsCurve) -> NurbsCurve {
    let (t0, t1) = curve.parameter_range();
    let span = t1 - t0;
    if span > 0.0 && (t0 != 0.0 || t1 != 1.0) {
        for k in &mut curve.knots {
            *k = (*k - t0) / span;
        }
    }
    curve
}

/// Even-odd point-in-polygon over a sampled outer wire.
pub(super) fn point_in_polygon(poly: &[(f64, f64)], p: (f64, f64)) -> bool {
    let mut inside = false;
    let n = poly.len();
    for i in 0..n {
        let (x0, y0) = poly[i];
        let (x1, y1) = poly[(i + 1) % n];
        if (y0 > p.1) != (y1 > p.1) {
            let x = x0 + (p.1 - y0) / (y1 - y0) * (x1 - x0);
            if p.0 < x {
                inside = !inside;
            }
        }
    }
    inside
}
