//! Exact area of a planar face from its boundary curves.
//!
//! By Stokes, a closed planar loop's vector area is `½ ∮ p × p′ dt`.
//! The face's vector area is the outer loop's, less each hole's
//! magnitude along the outer loop's normal (holes need not wind
//! opposite to the outer loop). Integrating the edge
//! CURVES (not the polygon through the vertices) is what makes a
//! circular hole bounded by one closed edge count (PK-771).

use popkern_core::geometry::{CurveEval, CurveGeometry, CurveOnEdge};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, LoopHandle};

/// 8-point Gauss–Legendre nodes on [-1, 1].
const GL_X: [f64; 8] = [
    -0.960_289_856_497_536_2,
    -0.796_666_477_413_626_7,
    -0.525_532_409_916_329,
    -0.183_434_642_495_649_8,
    0.183_434_642_495_649_8,
    0.525_532_409_916_329,
    0.796_666_477_413_626_7,
    0.960_289_856_497_536_2,
];
/// Matching Gauss–Legendre weights.
const GL_W: [f64; 8] = [
    0.101_228_536_290_376_26,
    0.222_381_034_453_374_47,
    0.313_706_645_877_887_3,
    0.362_683_783_378_362,
    0.362_683_783_378_362,
    0.313_706_645_877_887_3,
    0.222_381_034_453_374_47,
    0.101_228_536_290_376_26,
];

/// Largest angular step one quadrature panel spans on a conic. An
/// eighth of a turn keeps the 8-point rule's error below 1e-15
/// relative on the trigonometric integrand.
const CONIC_PANEL: f64 = std::f64::consts::FRAC_PI_4;

/// Quadrature panels per NURBS knot span. A rational quadratic
/// quarter circle reads 9e-12 relative on one panel and below 1e-15
/// on four.
const NURBS_SPAN_PANELS: usize = 4;

/// Vector area of a planar face: the outer loop's `½ ∮ (p − o) × p′ dt`
/// less each inner loop's magnitude along it; its length is the net
/// area (outer region minus holes) and it points along the outer
/// loop's winding normal. `o`
/// is any point; closed loops make the sum independent of it, and
/// the face's first vertex keeps the products well conditioned.
/// `None` when a loop is not a well-formed ring.
pub(crate) fn planar_face_vector_area(brep: &BRep, face: FaceHandle) -> Option<Vec3> {
    let face_data = brep.faces.get(face)?;
    let origin = brep
        .half_edges
        .get(brep.loops.get(face_data.outer_loop)?.half_edge)
        .and_then(|he| brep.vertices.get(he.origin))?
        .position;
    // Holes subtract by magnitude along the outer loop's normal, not by
    // their winding: an inner loop wound like its outer one (a tube
    // cap's bore rim, measured on the arm's forearm tube) otherwise
    // ADDED its disc and the tube read 3.9x its volume.
    let outer = loop_moment(brep, face_data.outer_loop, origin)?;
    let Some(normal) = outer.normalized() else {
        return Some(outer * 0.5);
    };
    let mut net = outer.dot(&normal);
    for &inner in &face_data.inner_loops {
        net -= loop_moment(brep, inner, origin)?.dot(&normal).abs();
    }
    Some(normal * (0.5 * net))
}

/// A planar face's divergence-theorem volume term `⅓ ∬ p · n dA`:
/// every point of the plane has the same `p · n`, so the term is
/// `⅓ (o · A)` for any vertex `o` and the face's vector area `A`
/// (signed by the loops' winding). Zero for a face whose loops are
/// not well-formed rings, as the polygon fan gave.
pub(crate) fn planar_face_divergence(brep: &BRep, face: FaceHandle) -> f64 {
    let Some(area) = planar_face_vector_area(brep, face) else {
        return 0.0;
    };
    let Some(o) = brep
        .faces
        .get(face)
        .and_then(|f| brep.loops.get(f.outer_loop))
        .and_then(|l| brep.half_edges.get(l.half_edge))
        .and_then(|he| brep.vertices.get(he.origin))
    else {
        return 0.0;
    };
    let p = o.position;
    (p.x * area.x + p.y * area.y + p.z * area.z) / 3.0
}

/// `∮ (p − o) × p′ dt` around one loop.
fn loop_moment(brep: &BRep, loop_handle: LoopHandle, o: Point3) -> Option<Vec3> {
    let start = brep.loops.get(loop_handle)?.half_edge;
    let mut current = start;
    let mut sum = Vec3::zero();
    let mut guard = 0usize;
    loop {
        sum = sum + half_edge_moment(brep, current, o)?;
        current = brep.half_edges.get(current)?.next;
        guard += 1;
        if current == start {
            return Some(sum);
        }
        if guard > brep.half_edges.len() {
            return None;
        }
    }
}

/// One half-edge's contribution, traversed origin → destination.
fn half_edge_moment(brep: &BRep, half_edge: HalfEdgeHandle, o: Point3) -> Option<Vec3> {
    let he = brep.half_edges.get(half_edge)?;
    let from = brep.vertices.get(he.origin)?.position;
    let to = brep
        .vertices
        .get(brep.half_edges.get(he.next)?.origin)?
        .position;
    let edge = brep.edges.get(he.edge)?;
    let chord = (from - o).cross(&(to - o));
    let Some(curve) = edge.curve.as_ref() else {
        // A curveless edge is read as the straight chord between its
        // vertices, the same reading the vertex polygon gave.
        return Some(chord);
    };
    if matches!(curve.geometry, CurveGeometry::Line { .. }) {
        return Some(chord);
    }
    let forward = traverses_forward(brep, half_edge, curve, from, to);
    let (t0, t1) = if forward {
        (curve.t_start, curve.t_end)
    } else {
        (curve.t_end, curve.t_start)
    };
    Some(curve_moment(&curve.geometry, t0, t1, o))
}

/// Whether the half-edge runs its edge's curve in increasing `t`.
/// Open edges answer by their endpoints. A closed edge's endpoints
/// coincide, so it answers by the PK-303 convention: the edge's
/// canonical half-edge runs the stored curve forward.
pub(crate) fn traverses_forward(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    curve: &CurveOnEdge,
    from: Point3,
    to: Point3,
) -> bool {
    let p0 = curve.geometry.point_at(curve.t_start);
    let p1 = curve.geometry.point_at(curve.t_end);
    let span = p0.distance_to(&p1);
    let scale = from.distance_to(&to).max(span);
    if span > 1e-9 * scale.max(1.0) {
        return from.distance_to(&p0) + to.distance_to(&p1)
            <= from.distance_to(&p1) + to.distance_to(&p0);
    }
    brep.half_edges
        .get(half_edge)
        .and_then(|he| brep.edges.get(he.edge))
        .is_some_and(|edge| edge.half_edge == half_edge)
}

/// Composite Gauss–Legendre `∫ (p − o) × p′ dt` over `[t0, t1]`
/// (either order), in panels within each NURBS knot span or of at
/// most an eighth turn of a conic.
fn curve_moment(geometry: &CurveGeometry, t0: f64, t1: f64, o: Point3) -> Vec3 {
    let (lo, hi) = (t0.min(t1), t0.max(t1));
    let mut breaks = vec![lo];
    match geometry {
        CurveGeometry::Nurbs(nurbs) => {
            // A rational span's integrand is not a polynomial, so each
            // knot span takes NURBS_SPAN_PANELS panels.
            let mut knots = vec![lo];
            knots.extend(nurbs.knots.iter().copied().filter(|&k| k > lo && k < hi));
            knots.push(hi);
            knots.dedup();
            for w in knots.windows(2) {
                for i in 1..=NURBS_SPAN_PANELS {
                    breaks.push(w[0] + (w[1] - w[0]) * i as f64 / NURBS_SPAN_PANELS as f64);
                }
            }
            breaks.pop();
        }
        _ => {
            let panels = ((hi - lo) / CONIC_PANEL).ceil().clamp(1.0, 4096.0) as usize;
            for i in 1..panels {
                breaks.push(lo + (hi - lo) * i as f64 / panels as f64);
            }
        }
    }
    breaks.push(hi);
    let mut total = Vec3::zero();
    for w in breaks.windows(2) {
        let (a, b) = (w[0], w[1]);
        let half = 0.5 * (b - a);
        let mid = 0.5 * (a + b);
        for (x, wt) in GL_X.iter().zip(GL_W) {
            let t = mid + half * x;
            let p = geometry.point_at(t) - o;
            let d = geometry.derivative_at(t, 1);
            total = total + p.cross(&d) * (wt * half);
        }
    }
    if t1 < t0 { total * -1.0 } else { total }
}
