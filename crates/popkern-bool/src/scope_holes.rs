//! The far holes of the active scope (stage 2b, PK-814): inner loops of
//! active planar faces the scoped boolean carries like far faces.

use std::collections::BTreeSet;

use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::topology::{BRep, FaceHandle, LoopHandle, VertexHandle};
use slotmap::SecondaryMap;

/// The inner loops of active planar faces whose conservative curve box
/// misses the tool's box padded by `pad`, sharing no vertex with the
/// face's other loops (the walk's peel premises, `identity_peel`): on a
/// plane a loop enclosed by another lies in that loop's box, and every
/// record lies in the tool's box, so no record touches or encloses such
/// a hole.
pub(super) fn far_holes(
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
        // One buffer for every loop walk below (a disc face's 64 loops,
        // walked per cut).
        let mut hes = Vec::new();
        for &lh in &loops {
            loop_half_edges_into(a, lh, &mut hes);
            for &h in &hes {
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
            loop_half_edges_into(a, lh, &mut hes);
            let shares = hes.iter().any(|&h| {
                a.half_edges
                    .get(h)
                    .is_some_and(|x| bounds.get(x.origin).is_some_and(|&(n, _)| n > 1))
            });
            let clear = !shares
                && loop_box(a, &hes, &mut points).is_some_and(|(l, h)| {
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

/// A loop's half-edges in loop order, into `out` (cleared first).
fn loop_half_edges_into(
    a: &BRep,
    lh: popkern_core::topology::LoopHandle,
    out: &mut Vec<popkern_core::topology::HalfEdgeHandle>,
) {
    out.clear();
    let Some(start) = a.loops.get(lh).map(|l| l.half_edge) else {
        return;
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
}

/// The conservative 3D box of a loop's half-edges `hes`: vertices plus
/// each edge curve's conservative extent (as the walk's peel uses).
fn loop_box(
    a: &BRep,
    hes: &[popkern_core::topology::HalfEdgeHandle],
    points: &mut Vec<popkern_core::geometry::Point3>,
) -> Option<(
    popkern_core::geometry::Point3,
    popkern_core::geometry::Point3,
)> {
    points.clear();
    for &h in hes {
        let he = a.half_edges.get(h)?;
        points.push(a.vertices.get(he.origin)?.position);
        if let Some(c) = a.edges.get(he.edge).and_then(|e| e.curve.as_ref()) {
            crate::face_box::conservative_curve_extent(c, points);
        }
    }
    let b = Aabb::from_points(points)?;
    Some((b.min, b.max))
}
