//! Seam bookkeeping for the closed-fillet surgery: seam-edge lookup,
//! seam-consistent contact anchors, and seam-curve rebuild. Split out
//! of `closed_fillet.rs` (500-line ratchet).

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, EdgeHandle, HalfEdgeHandle, LoopHandle, VertexHandle};

/// Seam-consistent anchor on a contact circle (PK-285): the point on
/// the exact contact circle radially nearest `v_old_pos` — NOT the
/// trim curve's own start point. For a closed rim the trim NURBS
/// inherits its parameterization origin from the offset-surface
/// intersection marcher's seed, which is an arbitrary, ULP-sensitive
/// artifact (measured: a 7e-15 extrude-height difference moved the
/// anchor from the seam angle to the diametrically opposite side, and
/// the rebuilt seam line crossed the cylinder as a chord — check-9
/// caught it at deviation = radius). Anchoring at the old rim-seam
/// vertex keeps the shortened seam edge a generator of its surface.
pub(super) fn seam_consistent_anchor(
    circle: &Option<(Point3, Vec3, f64)>,
    trim: &popkern_core::geometry::NurbsCurve,
    v_old_pos: Point3,
) -> Point3 {
    if let Some((center, normal, radius)) = circle {
        let d = v_old_pos - *center;
        let radial = d - *normal * d.dot(normal);
        if let Some(dir) = radial.normalized() {
            return *center + dir * *radius;
        }
    }
    // No exact circle (non-circular closed trim): closest sampled
    // point on the trim curve to the old vertex — same
    // seam-consistency contract, sampled.
    let (t0, t1) = trim.parameter_range();
    let mut best = trim.point_at(t0);
    let mut best_dist = f64::INFINITY;
    for i in 0..=64 {
        let t = t0 + (t1 - t0) * (i as f64) / 64.0;
        let p = trim.point_at(t);
        let dist = p.distance_to(&v_old_pos);
        if dist < best_dist {
            best_dist = dist;
            best = p;
        }
    }
    best
}

/// Find a half-edge in `loop_handle` other than `skip` whose
/// `origin` is `vertex` — i.e. an edge incident to the rim vertex
/// from a non-rim direction. For the cylinder case this returns
/// the seam edge half-edge.
pub(super) fn find_seam_edge_at_vertex(
    brep: &BRep,
    loop_handle: LoopHandle,
    skip: HalfEdgeHandle,
    vertex: VertexHandle,
) -> Option<EdgeHandle> {
    let start = brep.loops.get(loop_handle)?.half_edge;
    let mut cursor = start;
    let mut iters = 0;
    loop {
        iters += 1;
        if iters > 10_000 {
            return None;
        }
        if cursor != skip {
            let h = brep.half_edges.get(cursor)?;
            if h.origin == vertex {
                return Some(h.edge);
            }
        }
        let next = brep.half_edges.get(cursor)?.next;
        cursor = next;
        if cursor == start {
            break;
        }
    }
    None
}

/// Replace a seam edge's curve geometry with a straight line that
/// preserves the canonical "forward" direction of the edge — the
/// curve starts at `edge.half_edge.origin`'s current position and
/// ends at the twin's `origin` current position. Used to update
/// the seam after one of its endpoint vertices moved to the new
/// contact-circle height.
pub(super) fn rebuild_seam_curve(brep: &mut BRep, edge: EdgeHandle) {
    use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};

    let Some(edata) = brep.edges.get(edge) else {
        return;
    };
    let he = edata.half_edge;
    let Some(he_data) = brep.half_edges.get(he).cloned() else {
        return;
    };
    let Some(twin_data) = brep.half_edges.get(he_data.twin).cloned() else {
        return;
    };

    let v_start_pos = brep.vertices[he_data.origin].position;
    let v_end_pos = brep.vertices[twin_data.origin].position;
    let direction = Vec3::new(
        v_end_pos.x - v_start_pos.x,
        v_end_pos.y - v_start_pos.y,
        v_end_pos.z - v_start_pos.z,
    );
    if direction.length() < 1e-12 {
        return;
    }

    brep.edges[edge].curve = Some(CurveOnEdge {
        geometry: CurveGeometry::Line {
            origin: v_start_pos,
            direction,
        },
        t_start: 0.0,
        t_end: 1.0,
    });
}
