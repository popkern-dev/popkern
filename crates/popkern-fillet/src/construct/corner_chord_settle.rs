//! PK-141: settle the strip x strip corner chord's FAR endpoint.
//!
//! After a corner vertex snaps to the rails' wall crossing, the
//! chord's other endpoint is still the stale planner position;
//! alternating-project it onto ALL its incident surfaces so both
//! chord ends sit on the envelope intersection (transactional —
//! snap + reconcile accepted only when validation passes and the
//! residual strictly improves).

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

/// Settle the far endpoint of the strip x strip corner chord at a
/// snapped corner vertex onto all of ITS incident surfaces
/// (transactional; returns Some(moved_count) when the corner has a
/// chord, None otherwise).
pub(crate) fn settle_corner_chord_far_end(
    brep: &mut BRep,
    corner_vertex: VertexHandle,
    tolerance: &Tolerance,
) -> Option<usize> {
    use popkern_core::geometry::SurfaceEval;
    let chord = brep.edges.iter().find_map(|(_, ed)| {
        let he = brep.half_edges.get(ed.half_edge)?;
        let tw = brep.half_edges.get(he.twin)?;
        let f1 = brep.loops.get(he.loop_)?.face;
        let f2 = brep.loops.get(tw.loop_)?.face;
        let is_strip = |f: FaceHandle| {
            matches!(
                brep.faces
                    .get(f)
                    .and_then(|d| d.surface.as_ref())
                    .map(|s| &s.geometry),
                Some(popkern_core::geometry::surface::SurfaceGeometry::Nurbs(_))
            )
        };
        if f1 == f2 || !is_strip(f1) || !is_strip(f2) {
            return None;
        }
        if he.origin == corner_vertex {
            Some(tw.origin)
        } else if tw.origin == corner_vertex {
            Some(he.origin)
        } else {
            None
        }
    });
    let far = chord?;
    let pos = brep.vertices.get(far)?.position;
    // Incident surfaces of the far endpoint.
    let mut surfaces: Vec<popkern_core::geometry::surface::SurfaceGeometry> = Vec::new();
    for (_, he) in brep.half_edges.iter() {
        if he.origin != far {
            continue;
        }
        if let Some(l) = brep.loops.get(he.loop_)
            && let Some(fd) = brep.faces.get(l.face)
            && let Some(sof) = fd.surface.as_ref()
        {
            let sig = format!("{:?}", sof.geometry);
            if !surfaces.iter().any(|g| format!("{g:?}") == sig) {
                surfaces.push(sof.geometry.clone());
            }
        }
    }
    if surfaces.len() < 2 {
        return Some(0);
    }
    let mut p = pos;
    let residual = |q: &popkern_core::geometry::Point3| -> f64 {
        surfaces
            .iter()
            .filter_map(|g| {
                popkern_core::geometry::uv_projection::project_point_to_uv(g, q)
                    .map(|(u, v)| q.distance_to(&g.point_at(u, v)))
            })
            .fold(0.0_f64, f64::max)
    };
    let before = residual(&pos);
    if before < tolerance.linear().max(1e-9) * 10.0 {
        return Some(0);
    }
    for _ in 0..500 {
        for g in &surfaces {
            if let Some((u, v)) = popkern_core::geometry::uv_projection::project_point_to_uv(g, &p)
            {
                p = g.point_at(u, v);
            }
        }
        if residual(&p) < tolerance.linear().max(1e-9) {
            break;
        }
    }
    let after = residual(&p);
    if after >= before || p.distance_to(&pos) > before * 6.0 {
        return Some(0);
    }
    let mut cand = brep.clone();
    super::edge_reconcile::snap_vertex_and_adjacent_chords(
        &mut cand,
        far,
        p,
        popkern_core::topology::EdgeHandle::default(),
    );
    super::edge_reconcile::reconcile_edge_curves(&mut cand);
    if cand.validate().is_ok() {
        *brep = cand;
        Some(1)
    } else {
        Some(0)
    }
}
