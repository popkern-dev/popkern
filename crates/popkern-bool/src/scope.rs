//! Local boolean scope, stage 1 (docs/architecture/local-boolean-scope.md,
//! PK-813): the active region of operand `a` and the compare check that
//! measures whether serving it alone would be sound.
//!
//! Active = faces of `a` whose padded box meets `b`'s padded body box
//! (no ring: see `active_scope`); far = the rest, in `a`'s face order; frontier = edges between an active and
//! a far face. Stage 1 runs nothing scoped: under
//! `POPKERN_BOOL_SCOPE=compare` the unscoped pipeline's result is checked
//! against the premises a scoped run would rely on (design §3):
//!
//! 1. no intersection record lies on a far face;
//! 2. every far face comes through unchanged (same surface kind, the same
//!    loop vertex positions bit for bit) for a union or a subtraction,
//!    and is absent for an intersection;
//! 3. every frontier edge survives unsplit (an edge between the same two
//!    end positions).
//!
//! A violated premise fails the boolean in compare mode with the first
//! violation named, so a suite run under the variable lists every case
//! the scope would get wrong.

use std::collections::{BTreeSet, HashSet};

use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, LoopHandle, VertexHandle};
use slotmap::SecondaryMap;

use crate::error::BoolOperation;
use crate::intersect::IntersectionResult;

/// Operand `a`'s active region against a tool (see the module docs).
#[derive(Debug, Clone, Default)]
pub(crate) struct ActiveScope {
    /// Faces whose padded box meets the tool's.
    pub(crate) active: BTreeSet<FaceHandle>,
    /// The rest, in `a`'s face order.
    pub(crate) far: Vec<FaceHandle>,
    /// Edges between an active and a far face.
    pub(crate) frontier: BTreeSet<EdgeHandle>,
    /// Inner loops of active PLANAR faces that the tool cannot reach
    /// (stage 2b): carried like far faces, never entering the graph.
    pub(crate) far_holes:
        std::collections::BTreeMap<FaceHandle, Vec<popkern_core::topology::LoopHandle>>,
}

/// Whether compare mode is on (`POPKERN_BOOL_SCOPE=compare`).
pub(crate) fn compare_enabled() -> bool {
    std::env::var("POPKERN_BOOL_SCOPE").is_ok_and(|v| v == "compare")
}

fn padded(
    b: &Aabb,
    pad: f64,
) -> (
    popkern_core::geometry::Point3,
    popkern_core::geometry::Point3,
) {
    use popkern_core::geometry::Point3;
    (
        Point3::new(b.min.x - pad, b.min.y - pad, b.min.z - pad),
        Point3::new(b.max.x + pad, b.max.y + pad, b.max.z + pad),
    )
}

/// The active region of `a` against `b` (face boxes from the
/// intersection stage), with boxes padded by `pad` (the linear
/// tolerance plus the coincidence band).
pub(crate) fn active_scope(
    a: &BRep,
    aabbs_a: &[(FaceHandle, Aabb)],
    aabbs_b: &[(FaceHandle, Aabb)],
    pad: f64,
) -> ActiveScope {
    use popkern_core::geometry::Point3;
    let mut lo = Point3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
    let mut hi = Point3::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
    for (_, bx) in aabbs_b {
        let (l, h) = padded(bx, pad);
        lo = Point3::new(lo.x.min(l.x), lo.y.min(l.y), lo.z.min(l.z));
        hi = Point3::new(hi.x.max(h.x), hi.y.max(h.y), hi.z.max(h.z));
    }
    // Dense maps on `a`'s handles throughout: the scope is computed for
    // every boolean, over the whole body (PK-814: 0.44 ms of SipHash maps
    // at the slot disc's 64th cut).
    let mut boxes: SecondaryMap<FaceHandle, &Aabb> = SecondaryMap::with_capacity(aabbs_a.len());
    for (f, b) in aabbs_a {
        boxes.insert(*f, b);
    }
    // A face with no box (no surface) stays active: nothing proves it far.
    let active: BTreeSet<FaceHandle> = a
        .faces
        .keys()
        .filter(|&f| {
            boxes.get(f).is_none_or(|bx| {
                let (l, h) = padded(bx, pad);
                l.x <= hi.x
                    && h.x >= lo.x
                    && l.y <= hi.y
                    && h.y >= lo.y
                    && l.z <= hi.z
                    && h.z >= lo.z
            })
        })
        .collect();
    // No ring is needed (a refinement of design §3 rule 3, measured on
    // the 64-slot disc: a ring through shared vertices pulled every slot
    // wall in through the disc's top and bottom faces, far = 0). Every
    // record lies on the tool's surfaces, so inside its box; a face whose
    // padded box misses the tool's carries none, and neither do its
    // edges and vertices, which lie in that box. A face sharing an edge
    // whose box meets the tool's has a box meeting it too: it is a seed.
    let far: Vec<FaceHandle> = a.faces.keys().filter(|f| !active.contains(f)).collect();
    let mut far_set: SecondaryMap<FaceHandle, ()> = SecondaryMap::with_capacity(far.len());
    for &f in &far {
        far_set.insert(f, ());
    }
    let mut frontier = BTreeSet::new();
    for (eh, e) in &a.edges {
        let face_of = |he| {
            a.half_edges
                .get(he)
                .and_then(|h| a.loops.get(h.loop_))
                .map(|l| l.face)
        };
        let twin = a.half_edges.get(e.half_edge).map(|h| h.twin);
        let (f0, f1) = (face_of(e.half_edge), twin.and_then(face_of));
        if let (Some(f0), Some(f1)) = (f0, f1)
            && far_set.contains_key(f0) != far_set.contains_key(f1)
        {
            frontier.insert(eh);
        }
    }
    // The walk's peel pad: 100 linear tolerances, ten scope pads.
    let far_holes = far_holes(a, &active, (lo, hi), 10.0 * pad);
    ActiveScope {
        active,
        far,
        frontier,
        far_holes,
    }
}

/// The inner loops of active planar faces whose conservative curve box
/// misses the tool's box padded by `pad`, sharing no vertex with the
/// face's other loops (the walk's peel premises, `identity_peel`): on a
/// plane a loop enclosed by another lies in that loop's box, and every
/// record lies in the tool's box, so no record touches or encloses such
/// a hole.
fn far_holes(
    a: &BRep,
    active: &BTreeSet<FaceHandle>,
    (lo, hi): (
        popkern_core::geometry::Point3,
        popkern_core::geometry::Point3,
    ),
    pad: f64,
) -> std::collections::BTreeMap<FaceHandle, Vec<popkern_core::topology::LoopHandle>> {
    use popkern_core::geometry::surface::SurfaceGeometry;
    let mut out = std::collections::BTreeMap::new();
    for &fh in active {
        let Some(face) = a.faces.get(fh) else {
            continue;
        };
        if face.inner_loops.is_empty()
            || !matches!(
                face.surface.as_ref().map(|s| &s.geometry),
                Some(SurfaceGeometry::Plane { .. })
            )
        {
            continue;
        }
        // Loops each vertex bounds: a hole sharing a vertex with
        // another loop has one counted twice (one pass, not per hole).
        let loops: Vec<_> = std::iter::once(face.outer_loop)
            .chain(face.inner_loops.iter().copied())
            .collect();
        // Each vertex's loop count, and the last loop that counted it.
        let mut bounds: SecondaryMap<VertexHandle, (usize, LoopHandle)> = SecondaryMap::new();
        for &lh in &loops {
            for h in loop_half_edges(a, lh) {
                let Some(v) = a.half_edges.get(h).map(|x| x.origin) else {
                    continue;
                };
                match bounds.get_mut(v) {
                    Some((n, last)) if *last != lh => (*n, *last) = (*n + 1, lh),
                    Some(_) => {}
                    None => {
                        bounds.insert(v, (1, lh));
                    }
                }
            }
        }
        let mut holes = Vec::new();
        let mut points = Vec::new();
        for &lh in &face.inner_loops {
            let shares = loop_half_edges(a, lh).iter().any(|&h| {
                a.half_edges
                    .get(h)
                    .is_some_and(|x| bounds.get(x.origin).is_some_and(|&(n, _)| n > 1))
            });
            let clear = loop_box(a, lh, &mut points).is_some_and(|(l, h)| {
                l.x - pad > hi.x
                    || h.x + pad < lo.x
                    || l.y - pad > hi.y
                    || h.y + pad < lo.y
                    || l.z - pad > hi.z
                    || h.z + pad < lo.z
            });
            if clear && !shares {
                holes.push(lh);
            }
        }
        if !holes.is_empty() {
            out.insert(fh, holes);
        }
    }
    out
}

fn loop_half_edges(
    a: &BRep,
    lh: popkern_core::topology::LoopHandle,
) -> Vec<popkern_core::topology::HalfEdgeHandle> {
    let mut out = Vec::new();
    let Some(start) = a.loops.get(lh).map(|l| l.half_edge) else {
        return out;
    };
    let mut he = start;
    for _ in 0..=a.half_edges.len() {
        let Some(h) = a.half_edges.get(he) else { break };
        out.push(he);
        he = h.next;
        if he == start {
            break;
        }
    }
    out
}

/// A loop's conservative 3D box: vertices plus each edge curve's
/// conservative extent (as the walk's peel uses).
fn loop_box(
    a: &BRep,
    lh: popkern_core::topology::LoopHandle,
    points: &mut Vec<popkern_core::geometry::Point3>,
) -> Option<(
    popkern_core::geometry::Point3,
    popkern_core::geometry::Point3,
)> {
    points.clear();
    for h in loop_half_edges(a, lh) {
        let he = a.half_edges.get(h)?;
        points.push(a.vertices.get(he.origin)?.position);
        if let Some(c) = a.edges.get(he.edge).and_then(|e| e.curve.as_ref()) {
            crate::face_box::conservative_curve_extent(c, points);
        }
    }
    let b = Aabb::from_points(points)?;
    Some((b.min, b.max))
}

/// Every half-edge of the face's loops.
fn face_half_edges(brep: &BRep, face: FaceHandle) -> Vec<popkern_core::topology::HalfEdgeHandle> {
    let mut out = Vec::new();
    let Some(f) = brep.faces.get(face) else {
        return out;
    };
    for lh in std::iter::once(f.outer_loop).chain(f.inner_loops.iter().copied()) {
        let Some(start) = brep.loops.get(lh).map(|l| l.half_edge) else {
            continue;
        };
        let mut he = start;
        for _ in 0..=brep.half_edges.len() {
            let Some(h) = brep.half_edges.get(he) else {
                break;
            };
            out.push(he);
            he = h.next;
            if he == start {
                break;
            }
        }
    }
    out
}

/// A face's identity for the carry check: surface kind and its loop
/// vertex positions (bits), sorted.
fn face_key(brep: &BRep, face: FaceHandle) -> (String, Vec<[u64; 3]>) {
    let kind = brep
        .faces
        .get(face)
        .and_then(|f| f.surface.as_ref())
        .map(|s| format!("{:?}", std::mem::discriminant(&s.geometry)))
        .unwrap_or_default();
    let mut points: Vec<[u64; 3]> = face_half_edges(brep, face)
        .into_iter()
        .filter_map(|he| brep.half_edges.get(he))
        .filter_map(|h| brep.vertices.get(h.origin))
        .map(|v| {
            [
                v.position.x.to_bits(),
                v.position.y.to_bits(),
                v.position.z.to_bits(),
            ]
        })
        .collect();
    points.sort_unstable();
    (kind, points)
}

/// The first premise of the scope the unscoped `result` violates, or
/// `None` (see the module docs).
pub(crate) fn scope_violation(
    op: BoolOperation,
    a: &BRep,
    intersection: &IntersectionResult,
    scope: &ActiveScope,
    result: &BRep,
) -> Option<String> {
    let far: HashSet<FaceHandle> = scope.far.iter().copied().collect();
    if let Some(r) = intersection
        .records
        .iter()
        .find(|r| far.contains(&r.face_a))
    {
        return Some(format!("a record lies on far face {:?}", r.face_a));
    }
    let result_keys: HashSet<(String, Vec<[u64; 3]>)> =
        result.faces.keys().map(|f| face_key(result, f)).collect();
    for &f in &scope.far {
        let carried = result_keys.contains(&face_key(a, f));
        let expected = !matches!(op, BoolOperation::Intersection);
        if carried != expected {
            return Some(format!(
                "far face {f:?} {} the result ({op:?}){}",
                if carried {
                    "is in"
                } else {
                    "is not carried unchanged into"
                },
                if carried {
                    String::new()
                } else {
                    carry_miss(a, f, result)
                }
            ));
        }
    }
    let ends = |brep: &BRep, e: &popkern_core::topology::Edge| {
        let h = brep.half_edges.get(e.half_edge)?;
        let p = brep.vertices.get(h.origin)?.position;
        let q = brep
            .vertices
            .get(brep.half_edges.get(h.twin)?.origin)?
            .position;
        let (p, q) = (
            [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()],
            [q.x.to_bits(), q.y.to_bits(), q.z.to_bits()],
        );
        Some(if p <= q { (p, q) } else { (q, p) })
    };
    if !matches!(op, BoolOperation::Intersection) {
        let result_edges: HashSet<_> = result
            .edges
            .values()
            .filter_map(|e| ends(result, e))
            .collect();
        for &eh in &scope.frontier {
            if let Some(key) = a.edges.get(eh).and_then(|e| ends(a, e))
                && !result_edges.contains(&key)
            {
                return Some(format!("frontier edge {eh:?} did not survive unsplit"));
            }
        }
    }
    None
}

/// Why a far face did not come through: the nearest result face of its
/// surface kind by vertex positions — its vertex count against the far
/// face's and the largest distance from a far-face vertex to the
/// nearest of its vertices.
fn carry_miss(a: &BRep, face: FaceHandle, result: &BRep) -> String {
    let (kind, _) = face_key(a, face);
    let points = |brep: &BRep, f: FaceHandle| -> Vec<popkern_core::geometry::Point3> {
        face_half_edges(brep, f)
            .into_iter()
            .filter_map(|he| brep.half_edges.get(he))
            .filter_map(|h| brep.vertices.get(h.origin))
            .map(|v| v.position)
            .collect()
    };
    let mine = points(a, face);
    let mut best: Option<(f64, usize, FaceHandle)> = None;
    for rf in result.faces.keys() {
        if face_key(result, rf).0 != kind {
            continue;
        }
        let theirs = points(result, rf);
        let worst = mine
            .iter()
            .map(|p| {
                theirs
                    .iter()
                    .map(|q| p.distance_to(q))
                    .fold(f64::INFINITY, f64::min)
            })
            .fold(0.0, f64::max);
        if best.is_none_or(|b| worst < b.0) {
            best = Some((worst, theirs.len(), rf));
        }
    }
    match best {
        Some((d, n, rf)) => format!(
            ": nearest {kind} result face {rf:?} has {n} loop vertices (far face {}), worst vertex miss {d:.3e}",
            mine.len()
        ),
        None => format!(": no {kind} face in the result"),
    }
}

/// Far faces below which a scoped run is not worth its carry.
const MIN_FAR_FACES: usize = 8;

/// The scope a boolean SERVES with (stage 2, PK-814), or `None` to run
/// unscoped: only a union or a subtraction keeps every far face; compare
/// mode runs unscoped (it checks the scope's premises on that result);
/// `POPKERN_BOOL_SCOPE=off` is the kill switch; and a scope with few far
/// faces gains nothing. The pad is the compare harness's (stage 1 found
/// no premise violated with it across the boolean, part, rhai,
/// examples, blend, fastener and hinge suites).
pub(crate) fn serving_scope(
    a: &BRep,
    intersection: &IntersectionResult,
    op: BoolOperation,
    tolerance: &popkern_core::tolerance::Tolerance,
) -> Option<ActiveScope> {
    if !matches!(op, BoolOperation::Union | BoolOperation::Subtraction)
        || compare_enabled()
        || std::env::var("POPKERN_BOOL_SCOPE").is_ok_and(|v| v == "off")
    {
        return None;
    }
    let scope = active_scope(
        a,
        &intersection.aabbs_a,
        &intersection.aabbs_b,
        10.0 * tolerance.linear(),
    );
    if popkern_core::trace::trace_enabled("boolean") {
        eprintln!(
            "  [scope] serving: active {} far {} far holes {}",
            scope.active.len(),
            scope.far.len(),
            scope.far_holes.values().map(Vec::len).sum::<usize>()
        );
    }
    (scope.far.len() >= MIN_FAR_FACES).then_some(scope)
}

/// Compare mode (`POPKERN_BOOL_SCOPE=compare`): the first premise the
/// cellular assembly `assembled` (before validation, pcurve recompute
/// and the serve postprocess, which stage 2 leaves whole-body) violates,
/// or `None`. Off, or with no far face, it costs nothing.
pub(crate) fn compare_violation(
    a: &BRep,
    intersection: &IntersectionResult,
    op: BoolOperation,
    tolerance: &popkern_core::tolerance::Tolerance,
    assembled: &BRep,
) -> Option<String> {
    if !compare_enabled() {
        return None;
    }
    let scope = active_scope(
        a,
        &intersection.aabbs_a,
        &intersection.aabbs_b,
        10.0 * tolerance.linear(),
    );
    if popkern_core::trace::trace_enabled("boolean") {
        eprintln!(
            "  [scope] active {} far {} frontier {}",
            scope.active.len(),
            scope.far.len(),
            scope.frontier.len()
        );
    }
    if scope.far.is_empty() {
        return None;
    }
    scope_violation(op, a, intersection, &scope, assembled)
}

#[cfg(test)]
#[path = "scope_tests.rs"]
mod tests;
