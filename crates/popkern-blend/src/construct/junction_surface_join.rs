//! Stale-vertex settling across smooth-through junctions (PK-177).
//!
//! At a smooth-through junction the chain continues over ONE
//! envelope; vertices left at stale original-boundary positions
//! measure as join-chord bulges (2.07e-2 on the boss cross-seam).
//! Settling snaps the junction vertex and its near neighbours onto
//! their incident-surface intersections.

use popkern_core::geometry::Vec3;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::BRep;

pub(crate) struct SmoothJunctionContext {
    pub vertex: popkern_core::topology::VertexHandle,
    pub ta: Vec3,
    pub tb: Vec3,
    pub trace: bool,
}

/// Smooth-through junction tail: settle the junction vertex onto
/// its incident-surface intersection (corner snap fires only for
/// true corners; a stale vertex measured a 2.07e-2 join-chord
/// bulge), then settle stale two-surface neighbours the same way.
pub(crate) fn handle_smooth_junction(
    brep: &mut BRep,
    ctx: SmoothJunctionContext,
    tolerance: &popkern_core::tolerance::Tolerance,
) -> usize {
    let SmoothJunctionContext {
        vertex,
        ta,
        tb,
        trace,
    } = ctx;
    if ta.dot(&tb).abs() < 0.99 {
        return 0;
    }
    let mut absorbed = 0usize;
    if let Some(n) =
        super::corner_chord_settle::settle_vertex_on_incident_surfaces(brep, vertex, tolerance)
        && n > 0
    {
        absorbed += n;
        if trace {
            eprintln!("[junction] {vertex:?} smooth junction vertex settled");
        }
    }
    // Neighbor stale-vertex settle (mirrors the corner pass):
    // two-surface vertices within a few loop hops of the junction
    // can sit at stale original-boundary positions — measured: the
    // cross-seam chord's far endpoint at the original wall position,
    // 2.07e-2 off BOTH mega strips, at the boss top junction.
    let mut neighbors: Vec<popkern_core::topology::VertexHandle> = Vec::new();
    for (_, he) in brep.half_edges.iter() {
        if he.origin != vertex {
            continue;
        }
        let mut cur = he.next;
        for _ in 0..3 {
            let Some(d) = brep.half_edges.get(cur) else {
                break;
            };
            if d.origin != vertex && !neighbors.contains(&d.origin) {
                neighbors.push(d.origin);
            }
            cur = d.next;
        }
    }
    for v in neighbors {
        let mut sigs: Vec<String> = Vec::new();
        let mut nurbs = 0usize;
        let mut cap_like = false;
        for (_, he) in brep.half_edges.iter() {
            if he.origin != v {
                continue;
            }
            // Both sides of every incident edge: a cap on the TWIN
            // side evades an origin-only scan, and moving a
            // cap-boundary vertex broke the cube-uniform cap G1
            // gate to 1.0 (measured).
            for l in [
                brep.loops.get(he.loop_),
                brep.half_edges
                    .get(he.twin)
                    .and_then(|tw| brep.loops.get(tw.loop_)),
            ]
            .into_iter()
            .flatten()
            {
                let Some(sof) = brep.faces.get(l.face).and_then(|fd| fd.surface.as_ref()) else {
                    continue;
                };
                let sig = format!("{:?}", sof.geometry);
                if !sigs.contains(&sig) {
                    match &sof.geometry {
                        SurfaceGeometry::Nurbs(_) => nurbs += 1,
                        SurfaceGeometry::Sphere { .. } | SurfaceGeometry::Gregory(_) => {
                            cap_like = true;
                        }
                        _ => {}
                    }
                    sigs.push(sig);
                }
            }
        }
        // 2 = stale boundary vertex; 3 = the junction TRIPLE with
        // the exact signature 2 strips + 1 analytic support
        // (measured: a 2-surface settle put the shared rail
        // endpoint on ONE strip's intersection, 1.9e-2 off the
        // other; the unconstrained 3-surface settle regressed the
        // cube-uniform cap G1 gate). Cap-adjacent vertices and
        // 4+-surface hole corners never settle here.
        let triple = sigs.len() == 3 && nurbs == 2 && !cap_like;
        if !(sigs.len() == 2 && !cap_like || triple) {
            continue;
        }
        if let Some(n) =
            super::corner_chord_settle::settle_vertex_on_incident_surfaces(brep, v, tolerance)
            && n > 0
        {
            absorbed += n;
            if trace {
                eprintln!("[junction] {vertex:?}: stale neighbor {v:?} settled");
            }
        }
    }
    absorbed
}
