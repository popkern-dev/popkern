//! Material-side audit for contested shared records (PK-494 layer-3
//! diagnostic) — split from `assemble.rs` under the frozen-ratchet
//! budget. Gated behind the `uv-scramble` topic plus
//! `POPKERN_CLAIM_TRACE_Z`.

use popkern_core::topology::BRep;

use crate::merge::fragment_classify::FragmentSelection;
use crate::merge::graph::{IntersectionGraph, Operand};
use crate::merge::walk::Fragment;

/// For each kept fragment traversing a Shared record near the query
/// height, print the OUTPUT-orientation data: traversal tangent T
/// (post-flip), outward normal N (surface normal · orientation,
/// negated on flip), and the material side M = N × T — M must point
/// into the fragment's own region for a consistent orientation.
pub(super) fn trace_material_sides(
    graph: &IntersectionGraph,
    frag: &Fragment,
    sel: &FragmentSelection,
    brep_a: &BRep,
    brep_b: &BRep,
) {
    // PK-494 layer-3 diagnostic (uv-scramble topic +
    // POPKERN_CLAIM_TRACE_Z): the material-side audit at contested
    // shared records. For each kept fragment traversing a Shared
    // record near the query height, print the OUTPUT-orientation
    // data: traversal tangent T (post-flip), outward normal N
    // (surface normal · orientation, negated on flip), and the
    // material side M = N × T — M must point into the fragment's
    // own region for a consistent orientation; the walk whose M
    // points away is the one whose convention is inverted.
    if popkern_core::trace::trace_enabled("uv-scramble")
        && let Some(zq) = std::env::var("POPKERN_CLAIM_TRACE_Z")
            .ok()
            .and_then(|v| v.parse::<f64>().ok())
    {
        use popkern_core::geometry::curve::CurveEval;
        use popkern_core::geometry::surface::SurfaceEval;
        for &raw_ghe in frag.outer_loop.iter() {
            let ghe = if sel.flip {
                graph.half_edges[raw_ghe].twin
            } else {
                raw_ghe
            };
            if !matches!(
                graph.half_edges[ghe].source,
                crate::merge::graph::EdgeSource::Shared { .. }
            ) {
                continue;
            }
            let twin = graph.half_edges[ghe].twin;
            let (na, nb) = (graph.half_edges[ghe].origin, graph.half_edges[twin].origin);
            let (pa, pb) = (graph.nodes[na].position, graph.nodes[nb].position);
            let mid = popkern_core::geometry::Point3::new(
                0.5 * (pa.x + pb.x),
                0.5 * (pa.y + pb.y),
                0.5 * (pa.z + pb.z),
            );
            if (mid.z - zq).abs() > 0.35 && (pa.z - zq).abs() > 0.2 && (pb.z - zq).abs() > 0.2 {
                continue;
            }
            // Traversal tangent from the he's curve (or chord).
            let t_vec = match graph.half_edges[ghe]
                .curve
                .as_ref()
                .or_else(|| graph.half_edges[twin].curve.as_ref())
            {
                Some(c) => {
                    let m = 0.5 * (c.t_start + c.t_end);
                    let d = 1e-4 * (c.t_end - c.t_start).abs().max(1e-9);
                    let (p0, p1) = (c.geometry.point_at(m - d), c.geometry.point_at(m + d));
                    let mut v =
                        popkern_core::geometry::Vec3::new(p1.x - p0.x, p1.y - p0.y, p1.z - p0.z);
                    // Curve param direction vs he direction: align to
                    // the he's chord.
                    let chord =
                        popkern_core::geometry::Vec3::new(pb.x - pa.x, pb.y - pa.y, pb.z - pa.z);
                    if v.dot(&chord) < 0.0 {
                        v = v * -1.0;
                    }
                    v
                }
                None => popkern_core::geometry::Vec3::new(pb.x - pa.x, pb.y - pa.y, pb.z - pa.z),
            };
            let src_brep = match frag.face.operand {
                Operand::A => brep_a,
                Operand::B => brep_b,
            };
            let Some(sof) = src_brep
                .faces
                .get(frag.face.face)
                .and_then(|f| f.surface.as_ref())
            else {
                continue;
            };
            let Some((u, v)) = sof.geometry.closest_point(&mid) else {
                continue;
            };
            let mut n = sof.geometry.normal_at(u, v);
            if !sof.orientation {
                n = n * -1.0;
            }
            if sel.flip {
                n = n * -1.0;
            }
            let m_side = n.cross(&t_vec);
            eprintln!(
                "  [mat-side] frag {:?} flip={} ghe {ghe} mid ({:.3},{:.3},{:.3}) T ({:.3},{:.3},{:.3}) N ({:.3},{:.3},{:.3}) M ({:.3},{:.3},{:.3})",
                frag.face,
                sel.flip,
                mid.x,
                mid.y,
                mid.z,
                t_vec.x,
                t_vec.y,
                t_vec.z,
                n.x,
                n.y,
                n.z,
                m_side.x,
                m_side.y,
                m_side.z
            );
        }
    }
}
