//! Hole → region nesting for the 3D-Newell fragment election.
//!
//! The flat assignment hung every hole on the primary fragment; with
//! CONCENTRIC rings (countersunk bore: the stale top rim sits inside
//! the new countersink circle) the annulus piece lost its hole and
//! the primary kept a stale one, leaving a single-use edge that made
//! the cellular assembly non-manifold — the silent graph-merge
//! fallback then shipped an invalid body. A hole moves to a separate
//! region only when that region strictly contains it (positive,
//! curve-sampled proof), winds WITH the outer, and shares no graph
//! edge with it — the twin-sharing exclusion that keeps the tuned
//! paired-ring tool semantics (same rule as assemble_chart_roles).

use std::collections::HashSet;

use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use super::super::graph::IntersectionGraph;
use super::walk_geometry::{self, signed_area};

/// Reassign the election's holes to their smallest strictly
/// containing separate region. `inner_loops_out` keeps the holes that
/// stay on the primary fragment; the returned vector is indexed like
/// `separate_outers` and carries each region's holes.
pub(super) fn assign_holes_to_regions(
    graph: &IntersectionGraph,
    outer_newell: &Vec3,
    inner_loops_out: &mut Vec<Vec<usize>>,
    separate_outers: &[Vec<usize>],
) -> Vec<Vec<Vec<usize>>> {
    let mut separate_inners: Vec<Vec<Vec<usize>>> = vec![Vec::new(); separate_outers.len()];
    let mut kept: Vec<Vec<usize>> = Vec::new();
    for hole in inner_loops_out.drain(..) {
        let hole_twins: HashSet<usize> = hole
            .iter()
            .map(|&he_idx| graph.half_edges[he_idx].twin)
            .collect();
        let mut best: Option<(f64, usize)> = None;
        for (si, sep) in separate_outers.iter().enumerate() {
            if sep.iter().any(|h| hole_twins.contains(h)) {
                continue;
            }
            let sep_newell = walk_geometry::loop_newell(graph, sep);
            if sep_newell.dot(outer_newell) <= 0.0 {
                continue;
            }
            // Positive-proof containment (curve-sampled, so 1-node
            // closed rims vote with real points).
            if !loop_strictly_inside_outer(graph, &hole, sep, &sep_newell) {
                continue;
            }
            let sep_area = signed_area(graph, sep).abs();
            if best.is_none_or(|(ba, _)| sep_area < ba) {
                best = Some((sep_area, si));
            }
        }
        match best {
            Some((_, si)) => separate_inners[si].push(hole),
            None => kept.push(hole),
        }
    }
    *inner_loops_out = kept;
    separate_inners
}

/// True when `lp` is UNAMBIGUOUSLY inside the candidate outer ring:
/// every sampled point inside and at least three of them. Closed
/// single-half-edge rims are sampled along their CURVE (they carry
/// one node, and the node-quorum of [`loop_outside_outer`] can never
/// clear on them — a countersunk bore's stale rim hole failed the
/// outside test against a faraway notch piece and was reassigned to
/// it, deleting the bore from the face). Used by the nesting pass:
/// holes move fragments only on positive proof.
fn loop_strictly_inside_outer(
    graph: &IntersectionGraph,
    lp: &[usize],
    outer: &[usize],
    outer_newell: &Vec3,
) -> bool {
    let n = *outer_newell;
    let (ax, ay) = if n.x.abs() >= n.y.abs() && n.x.abs() >= n.z.abs() {
        (1, 2)
    } else if n.y.abs() >= n.z.abs() {
        (0, 2)
    } else {
        (0, 1)
    };
    let coord = |p: &popkern_core::geometry::Point3, a: usize| match a {
        0 => p.x,
        1 => p.y,
        _ => p.z,
    };
    let sample_loop = |hes: &[usize]| -> Vec<(f64, f64)> {
        let mut pts = Vec::new();
        for &he_idx in hes {
            let he = &graph.half_edges[he_idx];
            let own = he.curve.is_some();
            let curve = he
                .curve
                .as_ref()
                .or_else(|| graph.half_edges[he.twin].curve.as_ref());
            let curved = curve.filter(|c| {
                !matches!(
                    &c.geometry,
                    popkern_core::geometry::curve::CurveGeometry::Line { .. }
                )
            });
            if let Some(c) = curved {
                use popkern_core::geometry::curve::CurveEval;
                let (t0, t1) = (c.t_start, c.t_end);
                for i in 0..8 {
                    let f = i as f64 / 8.0;
                    let t = if own {
                        t0 + (t1 - t0) * f
                    } else {
                        t1 - (t1 - t0) * f
                    };
                    let q = c.geometry.point_at(t);
                    pts.push((coord(&q, ax), coord(&q, ay)));
                }
            } else {
                let p = graph.nodes[he.origin].position;
                pts.push((coord(&p, ax), coord(&p, ay)));
            }
        }
        pts
    };
    let ring = sample_loop(outer);
    if ring.len() < 3 {
        return false;
    }
    let (mut lo, mut hi) = ((f64::MAX, f64::MAX), (f64::MIN, f64::MIN));
    for &(x, y) in &ring {
        lo = (lo.0.min(x), lo.1.min(y));
        hi = (hi.0.max(x), hi.1.max(y));
    }
    let band = 1.0e-3 * ((hi.0 - lo.0) + (hi.1 - lo.1)).max(1e-9);
    let mut inside_n = 0usize;
    for (px, py) in sample_loop(lp) {
        let mut min_d = f64::MAX;
        for i in 0..ring.len() {
            let (x0, y0) = ring[i];
            let (x1, y1) = ring[(i + 1) % ring.len()];
            let (dx, dy) = (x1 - x0, y1 - y0);
            let l2 = (dx * dx + dy * dy).max(1e-300);
            let t = (((px - x0) * dx + (py - y0) * dy) / l2).clamp(0.0, 1.0);
            let d = ((px - x0 - dx * t).powi(2) + (py - y0 - dy * t).powi(2)).sqrt();
            min_d = min_d.min(d);
        }
        if min_d < band {
            continue;
        }
        let mut crossings = 0;
        for i in 0..ring.len() {
            let (x0, y0) = ring[i];
            let (x1, y1) = ring[(i + 1) % ring.len()];
            if (y0 > py) != (y1 > py) {
                let xc = x0 + (py - y0) / (y1 - y0) * (x1 - x0);
                if xc > px {
                    crossings += 1;
                }
            }
        }
        if crossings % 2 == 1 {
            inside_n += 1;
        } else {
            return false; // any outside vote = not contained
        }
    }
    inside_n >= 3
}

/// TURNAROUND PROMOTION (PK-343): a "hole" that DOUBLES BACK through
/// an edge and its twin consecutively (the apex turnaround signature —
/// a cone lateral's seam halves are twins inside one loop) and whose
/// curve-sampled 3D box is strictly disjoint from the primary outer's
/// cannot be a hole of anything: it is a REGION of its own (the
/// cone-top loop above a through-cube union was filed as a hole of the
/// disjoint bottom frustum and the apex piece vanished from the union
/// wholesale). Ordinary holes — including the geometrically-disjoint
/// kept holes the tuned election relies on downstream (measured: a
/// holder-chain slot ring) — are untouched.
pub(super) fn promote_disjoint_turnaround_holes(
    graph: &IntersectionGraph,
    inner_loops_out: &mut Vec<Vec<usize>>,
    outer_loop: &[usize],
) -> Vec<Vec<usize>> {
    let has_turnaround = |lp: &[usize]| -> bool {
        lp.len() >= 2
            && (0..lp.len()).any(|i| {
                let a = lp[i];
                let b = lp[(i + 1) % lp.len()];
                graph.half_edges[a].twin == b
            })
    };
    let mut promoted: Vec<Vec<usize>> = Vec::new();
    inner_loops_out.retain(|hole| {
        if has_turnaround(hole) && loops_bbox_disjoint(graph, hole, outer_loop) {
            promoted.push(hole.clone());
            false
        } else {
            true
        }
    });
    promoted
}

/// Both operands and their face AABBs, computed once per boolean, so
/// the exterior-disjoint promotion's per-loop membership test does not
/// rebuild them for every face of every fragment.
pub(super) struct OppositeAabbs<'a> {
    brep_a: &'a BRep,
    brep_b: &'a BRep,
    a: Vec<(FaceHandle, Aabb)>,
    b: Vec<(FaceHandle, Aabb)>,
}

impl<'a> OppositeAabbs<'a> {
    /// Precompute both operands' face AABBs.
    pub(super) fn new(brep_a: &'a BRep, brep_b: &'a BRep) -> Self {
        Self {
            brep_a,
            brep_b,
            a: crate::intersect::compute_face_aabbs(brep_a),
            b: crate::intersect::compute_face_aabbs(brep_b),
        }
    }

    /// The solid and AABBs OPPOSITE the given face's operand — the
    /// membership target for a face-piece of `operand` (an A-face is
    /// tested against B, a B-face against A).
    fn opposite(&self, operand: super::super::graph::Operand) -> (&BRep, &[(FaceHandle, Aabb)]) {
        match operand {
            super::super::graph::Operand::A => (self.brep_b, &self.b),
            super::super::graph::Operand::B => (self.brep_a, &self.a),
        }
    }
}

/// EXTERIOR-DISJOINT PROMOTION (PK-438): a "hole" whose curve-sampled
/// 3D box is strictly disjoint from the primary outer AND whose
/// interior sample lies strictly OUTSIDE the opposite operand cannot be
/// a hole of this fragment — it is a spatially separate trim component
/// of the SAME surface. A tool wall that protrudes below the workpiece
/// carries an in-body wall piece (kept) and a below-body protrusion
/// piece (culled) on one face; the election glued the protrusion in as
/// an inner loop, so a single classify sample (taken in the in-body
/// piece) dragged the outside-body protrusion into the result, where it
/// triple-claimed the shared flush rim and stranded it (6 open edges).
///
/// The opposite operand for an A-face is B and vice versa. The
/// Outside-vs-opposite gate is what separates this class from the
/// disjoint holes the tuned nesting KEEPS (PK-343's holder-chain slot
/// ring, which classifies Inside/OnBoundary against the opposite
/// operand and stays attached). The membership sample is the loop's
/// curve-sampled centroid; only a decisive `Outside` promotes — an
/// `Inside`, `OnBoundary`, or classification error keeps the loop
/// attached (a curved protrusion whose centroid leaves its own surface
/// fails safe as a kept hole rather than a wrong promotion). A single
/// self-loop rim (a full closed circle) is never promoted: it is a
/// disk-shaped hole boundary by construction, not a stray patch. The
/// turnaround class (PK-343) is promoted separately and unconditionally.
pub(super) fn promote_exterior_disjoint_holes(
    graph: &IntersectionGraph,
    inner_loops_out: &mut Vec<Vec<usize>>,
    outer_loop: &[usize],
    operand: super::super::graph::Operand,
    opposite: &OppositeAabbs,
    tolerance: &Tolerance,
) -> Vec<Vec<usize>> {
    let (opp_brep, opp_aabbs) = opposite.opposite(operand);
    use popkern_core::geometry::curve::CurveEval as _;
    let loop_centroid = |lp: &[usize]| -> Option<Point3> {
        let mut sum = Point3::new(0.0, 0.0, 0.0);
        let mut n = 0.0_f64;
        for &h in lp {
            let he = &graph.half_edges[h];
            if let Some(c) = he
                .curve
                .as_ref()
                .or_else(|| graph.half_edges[he.twin].curve.as_ref())
            {
                for i in 0..8 {
                    let t = c.t_start + (c.t_end - c.t_start) * (i as f64) / 8.0;
                    let p = c.geometry.point_at(t);
                    sum.x += p.x;
                    sum.y += p.y;
                    sum.z += p.z;
                    n += 1.0;
                }
            } else {
                let p = graph.nodes[he.origin].position;
                sum.x += p.x;
                sum.y += p.y;
                sum.z += p.z;
                n += 1.0;
            }
        }
        (n > 0.0).then(|| Point3::new(sum.x / n, sum.y / n, sum.z / n))
    };
    // The pathology is a fragment whose outer piece is KEPT (its wall is
    // inside the opposite operand — a tool wall inside the workpiece) but
    // whose disjoint inner piece is DISCARDED (outside the opposite
    // operand — the same wall protruding past the body). The two pieces
    // sit in DIFFERENT membership regions, so they cannot share a
    // fragment. A workpiece face's genuine through-hole is the opposite
    // configuration — its outer face and its hole are BOTH outside the
    // (thin) tool — so both vote the same way and the loop stays put.
    // Requiring the outer to be strictly Inside is also what keeps the
    // centroid membership honest: a loop whose centroid leaves its own
    // surface (a tapered or non-convex face) fails safe as a kept hole.
    let outer_inside = loop_centroid(outer_loop).is_some_and(|c| {
        matches!(
            crate::classify::point_in_solid(&c, opp_brep, opp_aabbs, tolerance),
            Ok(crate::error::PointClassification::Inside)
        )
    });
    let mut promoted: Vec<Vec<usize>> = Vec::new();
    inner_loops_out.retain(|hole| {
        // A single self-loop half-edge is a full closed rim (a drilled
        // bore's exit circle): by construction the boundary of a
        // disk-shaped face hole, never a stray trimmed-away patch. Its
        // curve-sampled box can read disjoint from an outer on another
        // z-band and its centroid (the circle centre) lands ON a face
        // plane, giving a misleading Outside vote — so exclude it. The
        // below-body protrusion pathology is always a multi-edge trimmed
        // loop (rim arcs + trim edges).
        if hole.len() < 2 {
            return true;
        }
        // Only a fragment whose outer wall is kept (inside the opposite
        // operand) can be carrying a discarded protrusion.
        if !outer_inside {
            return true;
        }
        if !loops_bbox_disjoint(graph, hole, outer_loop) {
            return true;
        }
        let Some(centroid) = loop_centroid(hole) else {
            return true;
        };
        let pc = crate::classify::point_in_solid(&centroid, opp_brep, opp_aabbs, tolerance);
        let outside = matches!(pc, Ok(crate::error::PointClassification::Outside));
        if outside {
            promoted.push(hole.clone());
            false
        } else {
            true
        }
    });
    promoted
}

/// Strict 3D bounding-box separation of two loops — the
/// projection-free arm of the turnaround promotion (PK-343). Boxes
/// are CURVE-SAMPLED, never node-only: a one-anchor closed rim's node
/// box is a point, and a point-box reads "disjoint" from everything
/// (measured: node-only boxes evicted a holder chain's slot rings and
/// χ flipped 0 → 2).
fn loops_bbox_disjoint(graph: &IntersectionGraph, a: &[usize], b: &[usize]) -> bool {
    use popkern_core::geometry::curve::CurveEval as _;
    type P3 = popkern_core::geometry::Point3;
    let bbox = |hes: &[usize]| -> Option<(P3, P3)> {
        let mut lo = P3::new(f64::MAX, f64::MAX, f64::MAX);
        let mut hi = P3::new(f64::MIN, f64::MIN, f64::MIN);
        let mut any = false;
        let mut take = |p: P3| {
            lo.x = lo.x.min(p.x);
            lo.y = lo.y.min(p.y);
            lo.z = lo.z.min(p.z);
            hi.x = hi.x.max(p.x);
            hi.y = hi.y.max(p.y);
            hi.z = hi.z.max(p.z);
            any = true;
        };
        for &h in hes {
            for n in [
                graph.half_edges[h].origin,
                graph.half_edges[graph.half_edges[h].twin].origin,
            ] {
                take(graph.nodes[n].position);
            }
            let he = &graph.half_edges[h];
            if let Some(c) = he
                .curve
                .as_ref()
                .or_else(|| graph.half_edges[he.twin].curve.as_ref())
            {
                for i in 1..8 {
                    let t = c.t_start + (c.t_end - c.t_start) * (i as f64) / 8.0;
                    take(c.geometry.point_at(t));
                }
            }
        }
        any.then_some((lo, hi))
    };
    let (Some((alo, ahi)), Some((blo, bhi))) = (bbox(a), bbox(b)) else {
        return false;
    };
    // STRICTLY separated on some axis (touching boxes do not count —
    // a hole ring can share its outer's box surface).
    alo.x > bhi.x
        || blo.x > ahi.x
        || alo.y > bhi.y
        || blo.y > ahi.y
        || alo.z > bhi.z
        || blo.z > ahi.z
}
