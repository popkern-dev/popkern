//! Local helpers for split.rs — kept in a sibling file to fit the
//! 500-line per-file structural cap. Used only by `split.rs` and
//! `split_boundary.rs` siblings; not part of the public API.

use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::pcurve::{Pcurve, PcurveGeometry};

/// Walk `loop_start_he`'s loop and return the first half-edge (other
/// than `skip`) whose origin vertex coincides with `target` under the
/// tolerance — the second-pass endpoint matching of `split_face`,
/// shared by both endpoint orders.
pub(crate) fn find_loop_vertex_coincident(
    brep: &crate::topology::BRep,
    loop_start_he: crate::topology::HalfEdgeHandle,
    skip: Option<crate::topology::HalfEdgeHandle>,
    target: &crate::geometry::Point3,
    tolerance: &crate::tolerance::Tolerance,
) -> Option<crate::topology::HalfEdgeHandle> {
    let max_iters = brep.half_edges.len() + 1;
    let mut cursor = loop_start_he;
    let mut iters = 0;
    loop {
        let he_data = &brep.half_edges[cursor];
        let vertex_pos = brep.vertices[he_data.origin].position;
        if Some(cursor) != skip && tolerance.points_coincident(&vertex_pos, target) {
            return Some(cursor);
        }
        cursor = he_data.next;
        iters += 1;
        if cursor == loop_start_he || iters >= max_iters {
            return None;
        }
    }
}

/// Split a pcurve at parameter `t` (in [0,1] of the edge).
///
/// Returns `(left, right)` where `left` covers the first portion and
/// `right` covers the second. If the pcurve is `None` or splitting
/// fails, returns `(None, None)`. Handles both `Linear` and `Nurbs`
/// pcurve variants.
/// Split a pcurve at the parameter whose IMAGE on `surface` is nearest
/// `split_pos`, falling back to the parameter fraction `t` when the
/// pcurve cannot be evaluated against a surface or its image never
/// comes near the point.
///
/// A fraction split trusts that the pcurve's parameter range agrees
/// with the edge curve's. When an earlier surgery shortened the edge
/// without trimming its pcurve (measured: a strip rail's side-plane
/// pcurve still spanned the full plate height while the edge ended at
/// the rim fillet's rail), the fraction lands the cut 25 % past the
/// split point and the left piece's image overshoots the edge — every
/// distance-to-curve consumer then reads a gap of up to 0.7 r (PK-391).
/// The image is the truth; the parameter is derived from it.
pub(crate) fn split_pcurve_at_point(
    pcurve: &Option<Pcurve>,
    t: f64,
    surface: Option<&crate::geometry::surface::SurfaceGeometry>,
    split_pos: &crate::geometry::Point3,
    tolerance: f64,
) -> (Option<Pcurve>, Option<Pcurve>) {
    use crate::geometry::surface::SurfaceEval;
    let (Some(pc), Some(surface)) = (pcurve, surface) else {
        return split_pcurve_at(pcurve, t);
    };
    if !matches!(pc.geometry, PcurveGeometry::Nurbs(_)) {
        return split_pcurve_at(pcurve, t);
    }
    let (p_start, p_end) = (pc.t_start, pc.t_end);
    let span = p_end - p_start;
    if span.abs() <= 0.0 {
        return split_pcurve_at(pcurve, t);
    }
    let image_distance = |f: f64| {
        let uv = pc.point_at(p_start + f * span);
        surface.point_at(uv.x, uv.y).distance_to(split_pos)
    };
    const COARSE: usize = 64;
    let mut best = (t, image_distance(t));
    for i in 0..=COARSE {
        let f = i as f64 / COARSE as f64;
        let d = image_distance(f);
        if d < best.1 {
            best = (f, d);
        }
    }
    let window = 1.0 / COARSE as f64;
    let (mut lo, mut hi) = (
        (best.0 - window).clamp(0.0, 1.0),
        (best.0 + window).clamp(0.0, 1.0),
    );
    for _ in 0..40 {
        let m1 = lo + (hi - lo) * 0.382;
        let m2 = lo + (hi - lo) * 0.618;
        if image_distance(m1) < image_distance(m2) {
            hi = m2;
        } else {
            lo = m1;
        }
    }
    let refined = 0.5 * (lo + hi);
    let refined_distance = image_distance(refined);
    // Take the image-honest parameter only when the image actually
    // reaches the split point and it is an interior cut; otherwise the
    // pcurve is not a faithful image of this edge at all and the
    // fraction is the least surprising choice.
    if refined_distance <= tolerance * 10.0 && refined > 1e-9 && refined < 1.0 - 1e-9 {
        split_pcurve_at(pcurve, refined)
    } else {
        split_pcurve_at(pcurve, t)
    }
}

pub(crate) fn split_pcurve_at(pcurve: &Option<Pcurve>, t: f64) -> (Option<Pcurve>, Option<Pcurve>) {
    match pcurve {
        Some(pc) => match pc.split_at_fraction(t) {
            Some((left, right)) => (Some(left), Some(right)),
            None => (None, None),
        },
        None => (None, None),
    }
}

/// Reverse a NURBS curve (swap direction of traversal).
///
/// The reversed curve evaluates the same geometry but in the opposite
/// direction: `reversed.point_at(t_start)` equals
/// `original.point_at(t_end)`.
pub(crate) fn reverse_nurbs_curve(curve: &NurbsCurve) -> NurbsCurve {
    use crate::geometry::curve::CurveEval;
    let (t_start, t_end) = curve.parameter_range();
    let mut reversed_knots: Vec<f64> = curve
        .knots
        .iter()
        .rev()
        .map(|&k| t_start + t_end - k)
        .collect();
    reversed_knots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let mut reversed_cp = curve.control_points.clone();
    reversed_cp.reverse();

    let mut reversed_weights = curve.weights.clone();
    reversed_weights.reverse();

    NurbsCurve {
        degree: curve.degree,
        knots: reversed_knots,
        control_points: reversed_cp,
        weights: reversed_weights,
    }
}

/// The chart-branch occurrences of a cut's endpoints. A vertex a loop
/// visits more than once (a periodic face's seam vertex, reached once
/// up the seam and once down its twin) matches the cut's 3D endpoint at
/// every visit; the cut belongs at the one on its pcurve's branch
/// (PK-798: a fin's root rail ending on the housing seam at u = 2π was
/// cut at the u = 0 visit, and the split traded the housing for the
/// territory under the fin). Each visit is measured by the nearest
/// end of its two incident half-edges' pcurves to the trim's uv end.
/// A vertex visited once, or a visit without pcurve evidence, keeps the
/// first match.
pub(crate) fn chart_branch_occurrences(
    brep: &crate::topology::BRep,
    he_entry: crate::topology::HalfEdgeHandle,
    he_exit: crate::topology::HalfEdgeHandle,
    trim_uv: &NurbsCurve,
) -> (
    crate::topology::HalfEdgeHandle,
    crate::topology::HalfEdgeHandle,
) {
    use crate::geometry::curve::CurveEval;
    let (t0, t1) = trim_uv.parameter_range();
    let pick = |found: crate::topology::HalfEdgeHandle, uv: crate::geometry::Point3| {
        let vertex = brep.half_edges[found].origin;
        let visits: Vec<_> = loop_half_edges(brep, found)
            .into_iter()
            .filter(|&he| brep.half_edges[he].origin == vertex)
            .collect();
        if visits.len() < 2 {
            return found;
        }
        let miss = |he: crate::topology::HalfEdgeHandle| -> Option<f64> {
            let prev = brep.half_edges[he].prev;
            let mut best: Option<f64> = None;
            for pc in [&brep.half_edges[he].pcurve, &brep.half_edges[prev].pcurve] {
                let pc = pc.as_ref()?;
                for t in [pc.t_start, pc.t_end] {
                    let p = pc.point_at(t);
                    let d = ((p.x - uv.x).powi(2) + (p.y - uv.y).powi(2)).sqrt();
                    best = Some(best.map_or(d, |b: f64| b.min(d)));
                }
            }
            best
        };
        let measured: Option<Vec<(crate::topology::HalfEdgeHandle, f64)>> =
            visits.iter().map(|&he| miss(he).map(|d| (he, d))).collect();
        measured
            .and_then(|m| m.into_iter().min_by(|a, b| a.1.total_cmp(&b.1)))
            .map_or(found, |(he, _)| he)
    };
    (
        pick(he_entry, trim_uv.point_at(t0)),
        pick(he_exit, trim_uv.point_at(t1)),
    )
}

/// The half-edges of `start`'s loop, in order from `start`.
fn loop_half_edges(
    brep: &crate::topology::BRep,
    start: crate::topology::HalfEdgeHandle,
) -> Vec<crate::topology::HalfEdgeHandle> {
    let mut out = vec![start];
    let mut cursor = brep.half_edges[start].next;
    while cursor != start && out.len() <= brep.half_edges.len() {
        out.push(cursor);
        cursor = brep.half_edges[cursor].next;
    }
    out
}
