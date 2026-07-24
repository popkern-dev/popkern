//! PK-177: carve the hole-corner wall boundary onto its honest
//! trims.
//!
//! The served topology closes each hole corner with ONE strip x
//! strip chord that also bounds the WALL — geometrically impossible
//! to satisfy both (the strips' mutual intersection bows off the
//! wall by ~0.75 r mid-span, measured 3.75e-2 on the boss). The
//! honest closure: TWO wall trims (wall ∩ strip_a, wall ∩ strip_b)
//! meeting at the corner-ball wall contact W, with the chord kept
//! purely between the strips. The trim geometry comes from the
//! solved corner cap (`CornerCapGeometry`): the rolling ball's wall
//! contacts along each strip's extension path.
//!
//! Transactional per corner: candidate clone, validate()-gated.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, VertexHandle};

use super::corner_cap_geometry::CornerCapGeometry;

/// Loop context of `vertex` on a wall-surfaced face: the face, and
/// the loop half-edges in order starting at the side whose origin is
/// `vertex`.
fn wall_loop_context(
    brep: &BRep,
    wall: &SurfaceGeometry,
    vertex: VertexHandle,
) -> Option<(FaceHandle, Vec<HalfEdgeHandle>)> {
    let sig = format!("{wall:?}");
    for (h, he) in brep.half_edges.iter() {
        if he.origin != vertex {
            continue;
        }
        let Some(l) = brep.loops.get(he.loop_) else {
            continue;
        };
        let Some(fd) = brep.faces.get(l.face) else {
            continue;
        };
        let Some(sof) = fd.surface.as_ref() else {
            continue;
        };
        if format!("{:?}", sof.geometry) != sig {
            continue;
        }
        let mut sides = vec![h];
        let mut cur = he.next;
        for _ in 0..brep.half_edges.len() {
            if cur == h {
                break;
            }
            sides.push(cur);
            cur = brep.half_edges.get(cur)?.next;
        }
        return Some((l.face, sides));
    }
    None
}

/// Settle the stale wall-boundary vertices around a snapped corner
/// onto their honest wall ∩ strip intersections.
///
/// MEASURED structure (boss, all four corners probed): the -y
/// corners already carry the honest two-trim shape (both rail
/// anchors exist as loop vertices at d = 0); at the +y corners the
/// wall loop still runs along STALE original-edge segments (the
/// pre-fillet boss edge at its original offset, twinning the
/// region-merged strip) instead of the contact rails — the same
/// root as the PK141_TRUST_STRIPS 0.7 mis-side. Every stale vertex
/// sits ON the wall but off its strip twin; the incident-surface
/// settle (wall + strip, a transversal pair) moves it onto the
/// contact rail and the reconcile refits the adjacent curves.
/// Transactional per vertex (inside the settle).
pub(crate) fn carve_corner_wall_trims(
    brep: &mut BRep,
    cap: &CornerCapGeometry,
    corner_vertex: VertexHandle,
    tolerance: &Tolerance,
    trace: bool,
) -> Option<bool> {
    let (_, sides) = wall_loop_context(brep, &cap.shared_surface, corner_vertex)?;
    // Vertices within a few sides of the corner on the wall loop,
    // both directions (the stale runs are short: 1-2 segments).
    let mut candidates: Vec<VertexHandle> = Vec::new();
    let take = 4usize.min(sides.len());
    for h in sides.iter().take(take).chain(sides.iter().rev().take(take)) {
        if let Some(he) = brep.half_edges.get(*h)
            && he.origin != corner_vertex
            && !candidates.contains(&he.origin)
        {
            candidates.push(he.origin);
        }
    }
    let mut moved = 0usize;
    for v in candidates {
        // Only the STALE two-surface class settles: a vertex with 3+
        // distinct incident surfaces is a corner-class vertex — the
        // over-determined settle was falsified there (no common
        // zero-residual point; measured regression 3.5e-3 ->
        // 1.28e-2 when corner vertices re-settled through this
        // pass).
        let mut sigs: Vec<String> = Vec::new();
        for (_, he) in brep.half_edges.iter() {
            if he.origin != v {
                continue;
            }
            if let Some(l) = brep.loops.get(he.loop_)
                && let Some(fd) = brep.faces.get(l.face)
                && let Some(sof) = fd.surface.as_ref()
            {
                let sig = format!("{:?}", sof.geometry);
                if !sigs.contains(&sig) {
                    sigs.push(sig);
                }
            }
        }
        if sigs.len() != 2 {
            continue;
        }
        if let Some(n) =
            super::corner_chord_settle::settle_vertex_on_incident_surfaces(brep, v, tolerance)
            && n > 0
        {
            moved += n;
            if trace {
                eprintln!("[corner-carve] {corner_vertex:?}: stale wall vertex {v:?} settled");
            }
        }
    }
    if trace
        && let Some((chord, far, g1, g2)) =
            super::corner_chord_settle::corner_strip_chord(brep, corner_vertex)
    {
        // Chart-validity probe at the corner: where does each strip
        // chart place W, and how far off is its image?
        use popkern_core::geometry::SurfaceEval;
        let w = brep.vertices.get(corner_vertex)?.position;
        let far_pos = brep.vertices.get(far)?.position;
        for (label, g) in [("g1", &g1), ("g2", &g2)] {
            let (nu, nv) = match g {
                SurfaceGeometry::Nurbs(n) => (n.num_u, n.num_v),
                _ => (0, 0),
            };
            let img = popkern_core::geometry::uv_projection::project_point_to_uv(g, &w)
                .map(|(u, v)| {
                    let q = g.point_at(u, v);
                    format!("uv=({u:.4},{v:.4}) d={:.3e}", w.distance_to(&q))
                })
                .unwrap_or_else(|| "proj FAILED".into());
            eprintln!("[corner-carve]   chord {chord:?} {label} ({nu}x{nv}): W {img}");
        }
        eprintln!(
            "[corner-carve]   chord far {far:?} at ({:.4},{:.4},{:.4})",
            far_pos.x, far_pos.y, far_pos.z
        );
    }
    Some(moved > 0)
}
