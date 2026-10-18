//! `tangent-chain` trace output of the closed tangent-chain
//! construction (split from `tangent_chain_geometry` for the line cap).

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::topology::{BRep, FaceHandle};

/// `tangent-chain` trace: every edge of the strips' and their
/// neighbours' loops whose curve ends stray from its vertices.
pub(super) fn trace_edge_ends(brep: &BRep, strips: &[FaceHandle]) {
    let mut faces: Vec<FaceHandle> = strips.to_vec();
    for &f in strips {
        let lp = brep.faces[f].outer_loop;
        let start = brep.loops[lp].half_edge;
        let mut he = start;
        loop {
            let tw = brep.half_edges[he].twin;
            let g = brep.loops[brep.half_edges[tw].loop_].face;
            if !faces.contains(&g) {
                faces.push(g);
            }
            he = brep.half_edges[he].next;
            if he == start {
                break;
            }
        }
    }
    for f in faces {
        let face = &brep.faces[f];
        for lp in std::iter::once(face.outer_loop).chain(face.inner_loops.iter().copied()) {
            let start = brep.loops[lp].half_edge;
            let mut he = start;
            loop {
                let h = &brep.half_edges[he];
                let e = &brep.edges[h.edge];
                if let Some(c) = e.curve.as_ref() {
                    let eh = &brep.half_edges[e.half_edge];
                    let a = brep.vertices[eh.origin].position;
                    let b = brep.vertices[brep.half_edges[eh.twin].origin].position;
                    let da = c.geometry.point_at(c.t_start).distance_to(&a);
                    let db = c.geometry.point_at(c.t_end).distance_to(&b);
                    if da.max(db) > 1e-7 {
                        eprintln!(
                            "[tangent-chain] edge {:?} on {f:?} ends off its vertices by {da:.3e}/{db:.3e} ({:?})",
                            h.edge,
                            std::mem::discriminant(&c.geometry)
                        );
                    }
                }
                he = h.next;
                if he == start {
                    break;
                }
            }
        }
    }
}

/// `tangent-chain` trace: how far each strip's two end sections lie
/// from its surface (strip `i` runs from junction `i − 1` to `i`).
pub(super) fn trace_sections(surfaces: &[&SurfaceGeometry], contacts: &[(Point3, Point3)]) {
    let k = surfaces.len();
    for (i, surface) in surfaces.iter().enumerate() {
        // An open chain of one link has no junction: nothing to trace.
        let (Some(&before), Some(&after)) = (contacts.get((i + k - 1) % k), contacts.get(i)) else {
            continue;
        };
        let off = |p: Point3| {
            surface
                .closest_point(&p)
                .map(|(u, v)| surface.point_at(u, v).distance_to(&p))
                .unwrap_or(f64::INFINITY)
        };
        eprintln!(
            "[tangent-chain] strip {i} {} sections off the strip: start {:.3e}/{:.3e} end {:.3e}/{:.3e}",
            surface.kind_name(),
            off(before.0),
            off(before.1),
            off(after.0),
            off(after.1)
        );
    }
}
