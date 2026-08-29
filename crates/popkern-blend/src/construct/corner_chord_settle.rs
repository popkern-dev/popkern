//! PK-141: settle the strip x strip corner chord's FAR endpoint.
//!
//! After a corner vertex snaps to the rails' wall crossing, the
//! chord's other endpoint is still the stale planner position;
//! alternating-project it onto ALL its incident surfaces so both
//! chord ends sit on the envelope intersection (transactional —
//! snap + reconcile accepted only when validation passes and the
//! residual strictly improves).
//!
//! FALSIFIED mechanisms for the remaining 3.5e-3 mid-span residual
//! (both measured on the boss fixture, both reverted): (a) an
//! arc-length intersection march refit — the residual is ENDPOINT-
//! dominated (the corner vertex sits 3.7e-2 off the as-built
//! region-merged twin strip; every marched curve shares those
//! endpoints, before == after); (b) settling the corner vertex onto
//! its actual incident surfaces — with 4+ as-built surfaces no
//! zero-residual common point exists, gauge regressed to 1.28e-2.
//! The remaining mechanism is UPSTREAM: strip surface accuracy near
//! corners (the region-merged refit deviates from the envelope).

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

/// Locate the strip x strip chord at `corner_vertex`: the edge, its
/// far-end vertex, and the two strip surface geometries.
pub(super) fn corner_strip_chord(
    brep: &BRep,
    corner_vertex: VertexHandle,
) -> Option<(EdgeHandle, VertexHandle, SurfaceGeometry, SurfaceGeometry)> {
    brep.edges.iter().find_map(|(e, ed)| {
        let he = brep.half_edges.get(ed.half_edge)?;
        let tw = brep.half_edges.get(he.twin)?;
        let f1 = brep.loops.get(he.loop_)?.face;
        let f2 = brep.loops.get(tw.loop_)?.face;
        let strip_geometry = |f: FaceHandle| -> Option<SurfaceGeometry> {
            match brep.faces.get(f).and_then(|d| d.surface.as_ref()) {
                Some(s) if matches!(s.geometry, SurfaceGeometry::Nurbs(_)) => {
                    Some(s.geometry.clone())
                }
                _ => None,
            }
        };
        if f1 == f2 {
            return None;
        }
        let (g1, g2) = (strip_geometry(f1)?, strip_geometry(f2)?);
        let far = if he.origin == corner_vertex {
            tw.origin
        } else if tw.origin == corner_vertex {
            he.origin
        } else {
            return None;
        };
        Some((e, far, g1, g2))
    })
}

/// Settle the far endpoint of the strip x strip corner chord at a
/// snapped corner vertex onto all of ITS incident surfaces
/// (transactional; returns Some(moved_count) when the corner has a
/// chord, None otherwise).
pub(crate) fn settle_corner_chord_far_end(
    brep: &mut BRep,
    corner_vertex: VertexHandle,
    tolerance: &Tolerance,
) -> Option<usize> {
    let (_, far, _, _) = corner_strip_chord(brep, corner_vertex)?;
    settle_vertex_on_incident_surfaces(brep, far, tolerance)
}

/// Settle a vertex onto the mutual intersection of ALL its incident
/// face surfaces by alternating projection (transactional — snap +
/// reconcile accepted only when validation passes and the residual
/// strictly improves). NOT applied to the corner vertex itself:
/// with 4+ as-built surfaces no zero-residual common point exists
/// (measured: settling the boss corner vertex regressed the corner
/// chord gauge 3.5e-3 -> 1.28e-2); the analytic corner-ball snap
/// stays the corner target.
pub(crate) fn settle_vertex_on_incident_surfaces(
    brep: &mut BRep,
    far: VertexHandle,
    tolerance: &Tolerance,
) -> Option<usize> {
    use popkern_core::geometry::SurfaceEval;
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
