//! The boundary half of the chordal certificate (see
//! `chordal_certificate`): each boundary segment against the face edge
//! (or the two-edge corner chain) it lies on.

use std::collections::HashMap;

use popkern_core::geometry::certify::chord::{curve_piece_bound, curve_point_distance};
use popkern_core::geometry::curve::{CurveEval, CurveGeometry, closest_point_in_window};
use popkern_core::geometry::{Point3, PortableMath};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

use super::ChordalInconclusive;
use crate::mesh::TriMesh;

/// The largest bound over the face's boundary segments (triangle
/// edges used once within the face).
pub(super) fn boundary_bound(
    brep: &BRep,
    mesh: &TriMesh,
    face: FaceHandle,
    tris: &[usize],
    tolerance: f64,
) -> Result<f64, ChordalInconclusive> {
    // Edges keyed by welded position, not index: a pole fan's apex and
    // a seam's two sides are duplicated vertices, and their edges are
    // interior although each index pair is used once.
    let mut weld: HashMap<[u64; 3], usize> = HashMap::new();
    let mut id = |i: usize| -> usize {
        let p = mesh.vertices[i].position;
        let key = [p.x, p.y, p.z].map(|x| (x + 0.0).to_bits());
        let next = weld.len();
        *weld.entry(key).or_insert(next)
    };
    let mut uses: HashMap<(usize, usize), (usize, (usize, usize))> = HashMap::new();
    for &t in tris {
        let tri = mesh.triangles[t];
        for k in 0..3 {
            let (a, b) = (tri[k], tri[(k + 1) % 3]);
            let (wa, wb) = (id(a), id(b));
            uses.entry((wa.min(wb), wa.max(wb)))
                .or_insert((0, (a, b)))
                .0 += 1;
        }
    }
    let edges = face_edges(brep, face);
    let mut worst = 0.0_f64;
    let mut segments: Vec<(usize, usize)> = uses
        .into_values()
        .filter(|&(n, _)| n == 1)
        .map(|(_, e)| e)
        .collect();
    segments.sort_unstable();
    // Each boundary vertex ends two segments: its feet on every edge
    // are found once.
    let mut feet: HashMap<usize, Vec<Option<(f64, f64)>>> = HashMap::new();
    let mut feet_of = |i: usize| -> Vec<Option<(f64, f64)>> {
        let p = mesh.vertices[i].position;
        feet.entry(i)
            .or_insert_with(|| edges.iter().map(|&eh| edge_foot(brep, eh, p)).collect())
            .clone()
    };
    for (a, b) in segments {
        let (pa, pb) = (mesh.vertices[a].position, mesh.vertices[b].position);
        if pa.distance_to(&pb) <= f64::EPSILON * (1.0 + magnitude(pa)) {
            continue;
        }
        let ends = (feet_of(a), feet_of(b));
        worst = worst.max(segment_bound(brep, &edges, (pa, pb), &ends, tolerance)?);
    }
    Ok(worst)
}

/// The foot `(t, distance)` of `p` on edge `eh`'s curve within the
/// edge's own window (a curve's declared range can lie on another
/// period: an arc's angles vs a lifted edge window); `None` for an
/// edge without a curve or a finite window.
fn edge_foot(brep: &BRep, eh: EdgeHandle, p: Point3) -> Option<(f64, f64)> {
    let curve = brep.edges.get(eh)?.curve.as_ref()?;
    let window = (curve.t_start, curve.t_end);
    if !(window.0.is_finite() && window.1.is_finite()) {
        return None;
    }
    let t = window_foot(&curve.geometry, p, window);
    Some((t, p.distance_to(&curve.geometry.point_at(t))))
}

/// A foot of `p` on the curve within `window`: closed-form on a line
/// or circle (the angle lifted into the window), sampled and polished
/// otherwise. Any foot is valid for the bound, which adds the ends'
/// distances; the nearest keeps it tight.
pub(super) fn window_foot(g: &CurveGeometry, p: Point3, window: (f64, f64)) -> f64 {
    let (lo, hi) = (window.0.min(window.1), window.0.max(window.1));
    match g {
        CurveGeometry::Line { origin, direction } if direction.dot(direction) > 0.0 => {
            ((p - *origin).dot(direction) / direction.dot(direction)).clamp(lo, hi)
        }
        CurveGeometry::Circle { center, .. } | CurveGeometry::Arc { center, .. } => {
            // The evaluator's own basis, so the angle is its parameter.
            let bu = g.point_at(0.0) - *center;
            let bv = g.point_at(std::f64::consts::FRAC_PI_2) - *center;
            let d = p - *center;
            let t = lo + (d.dot(&bv).patan2(d.dot(&bu)) - lo).rem_euclid(std::f64::consts::TAU);
            if t <= hi {
                t
            } else if g.point_at(lo).distance_to(&p) <= g.point_at(hi).distance_to(&p) {
                lo
            } else {
                hi
            }
        }
        _ => closest_point_in_window(g, &p, lo, hi),
    }
}

fn magnitude(p: Point3) -> f64 {
    p.x.abs().max(p.y.abs()).max(p.z.abs())
}

/// Every edge on the face's loops.
fn face_edges(brep: &BRep, face: FaceHandle) -> Vec<EdgeHandle> {
    let Some(f) = brep.faces.get(face) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for &lh in std::iter::once(&f.outer_loop).chain(&f.inner_loops) {
        let Some(start) = brep.loops.get(lh).map(|l| l.half_edge) else {
            continue;
        };
        let mut he = start;
        for _ in 0..100_000 {
            let Some(h) = brep.half_edges.get(he) else {
                break;
            };
            if !out.contains(&h.edge) {
                out.push(h.edge);
            }
            he = h.next;
            if he == start {
                break;
            }
        }
    }
    out
}

/// One edge's projection of a segment's ends.
struct Projection<'a> {
    curve: &'a CurveGeometry,
    window: (f64, f64),
    ta: f64,
    tb: f64,
    da: f64,
    db: f64,
}

/// The bound of segment `pa pb` against the face's boundary: the edge
/// whose curve both ends are nearest, or, for a segment cutting a
/// corner (its ends on two edges meeting at a vertex), the two-edge
/// chain through that vertex, whichever bound is smaller (both are
/// pieces of the face boundary, so both are valid).
type Feet = Vec<Option<(f64, f64)>>;

fn segment_bound(
    brep: &BRep,
    edges: &[EdgeHandle],
    (pa, pb): (Point3, Point3),
    (feet_a, feet_b): &(Feet, Feet),
    tolerance: f64,
) -> Result<f64, ChordalInconclusive> {
    let mut found: Vec<Projection<'_>> = Vec::new();
    for (k, &eh) in edges.iter().enumerate() {
        let (Some(curve), Some((ta, da)), Some((tb, db))) = (
            brep.edges.get(eh).and_then(|e| e.curve.as_ref()),
            feet_a.get(k).copied().flatten(),
            feet_b.get(k).copied().flatten(),
        ) else {
            continue;
        };
        found.push(Projection {
            curve: &curve.geometry,
            window: (curve.t_start, curve.t_end),
            ta,
            tb,
            da,
            db,
        });
    }
    // The nearest edge, however near: the bound adds the ends'
    // distance from its curve, so a mesh boundary off its edge reads
    // as a proven excess, not as an unmatched segment.
    let Some(single) = found
        .iter()
        .min_by(|x, y| x.da.max(x.db).total_cmp(&y.da.max(y.db)))
    else {
        return Err(ChordalInconclusive::SegmentOnNoEdge { at: pa });
    };
    let mut bound = single_edge_bound(single, pa, pb, tolerance)?;
    if single.da.max(single.db) > tolerance {
        bound = bound.min(corner_chain_bound(&found, pa, pb, tolerance)?);
    }
    Ok(bound)
}

fn single_edge_bound(
    e: &Projection<'_>,
    pa: Point3,
    pb: Point3,
    tolerance: f64,
) -> Result<f64, ChordalInconclusive> {
    let decline = ChordalInconclusive::Decline;
    let (curve, ta, tb) = (e.curve, e.ta, e.tb);
    let ends = curve_point_distance(curve, ta, pa)
        .map_err(decline)?
        .max(curve_point_distance(curve, tb, pb).map_err(decline)?);
    let piece = |a: f64, b: f64| curve_piece_bound(curve, a, b, pa, pb, tolerance).map_err(decline);
    // The piece between the two parameters, or on a closed curve the
    // other way round (split at the seam for a NURBS). Both are pieces
    // of the edge's curve, so either bound is valid; the one whose
    // middle lies nearer the segment's is the piece the segment stands
    // for (bounding both costs a deep bisection of the long way round).
    let mid = Point3::new(
        0.5 * (pa.x + pb.x),
        0.5 * (pa.y + pb.y),
        0.5 * (pa.z + pb.z),
    );
    let near = |t: f64| curve.point_at(t).distance_to(&mid);
    let body = match curve_period(curve) {
        Some((start, period)) => {
            let (lo, hi) = (ta.min(tb), ta.max(tb));
            let around_mid = 0.5 * (hi + lo + period);
            if near(around_mid) < near(0.5 * (lo + hi)) {
                match curve {
                    CurveGeometry::Nurbs(_) => piece(start, lo)?.max(piece(hi, start + period)?),
                    _ => piece(hi, lo + period)?,
                }
            } else {
                piece(ta, tb)?
            }
        }
        None => piece(ta, tb)?,
    };
    Ok(body + ends)
}

/// The segment against two edges meeting at a vertex: `pa` on one,
/// `pb` on the other, each piece running from its end to the shared
/// vertex. Infinite when no such pair exists.
fn corner_chain_bound(
    found: &[Projection<'_>],
    pa: Point3,
    pb: Point3,
    tolerance: f64,
) -> Result<f64, ChordalInconclusive> {
    let decline = ChordalInconclusive::Decline;
    let mut best = f64::INFINITY;
    for (i, ea) in found.iter().enumerate() {
        for (j, eb) in found.iter().enumerate() {
            if i == j || ea.da > tolerance || eb.db > tolerance {
                continue;
            }
            // A window end of each edge at the same point: the vertex.
            for &va in &[ea.window.0, ea.window.1] {
                for &vb in &[eb.window.0, eb.window.1] {
                    let (qa, qb) = (ea.curve.point_at(va), eb.curve.point_at(vb));
                    if qa.distance_to(&qb) > tolerance {
                        continue;
                    }
                    let one = curve_piece_bound(ea.curve, ea.ta, va, pa, pb, tolerance)
                        .map_err(decline)?;
                    let two = curve_piece_bound(eb.curve, vb, eb.tb, pa, pb, tolerance)
                        .map_err(decline)?;
                    let ends = curve_point_distance(ea.curve, ea.ta, pa)
                        .map_err(decline)?
                        .max(curve_point_distance(eb.curve, eb.tb, pb).map_err(decline)?);
                    // The chain is continuous up to the vertex gap.
                    let gap = qa.distance_to(&qb);
                    best = best.min(one.max(two) + ends + gap);
                }
            }
        }
    }
    Ok(best)
}

/// `(start, period)` of a closed curve's parameter: 2π for circles,
/// the knot domain for a NURBS whose ends meet.
fn curve_period(curve: &CurveGeometry) -> Option<(f64, f64)> {
    match curve {
        // Angle parameters repeat every turn.
        CurveGeometry::Circle { .. }
        | CurveGeometry::Arc { .. }
        | CurveGeometry::Ellipse { .. } => Some((0.0, std::f64::consts::TAU)),
        CurveGeometry::Nurbs(n) => {
            let (lo, hi) = (n.knots[n.degree], n.knots[n.control_points.len()]);
            let closed = n.point_at(lo).distance_to(&n.point_at(hi))
                <= 1e-9 * (1.0 + magnitude(n.point_at(lo)));
            closed.then_some((lo, hi - lo))
        }
        _ => None,
    }
}
