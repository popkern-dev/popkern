//! Fan ends of FITTED strips (G2, variable-radius NURBS): crossings
//! and arcs solved on the strip surface itself (PK-530).
//!
//! No analytic section exists. Each crossing `q_j` is the root of the
//! strip's signed distance along its fan edge's line, bracketed from
//! the corner and bisected. The root is accepted only where the point
//! genuinely touches the strip, so a projection clamped at a patch
//! border cannot pass. Each arc is the strip ∩ fan-plane intersection
//! marched between its chain points (`surface_pair_march`). A fitted
//! strip meets its rails only within its fit accuracy, so the arcs
//! inherit the band recorded at the feet and declare it.

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{BRep, HalfEdgeHandle};

use crate::types::TerminalArcTrimDeclineReason as Reason;

use super::super::blend_seam_corner_repair::neighbour_chart_point;
use super::super::surface_pair_march::{
    MarchedIntersection, intersection_curve, intersection_pcurve, march_intersection,
};
use super::super::terminal_arc_author::StagedStubArc;
use super::super::terminal_arc_trim::CornerCandidate;
use super::Plane;

fn plane_geometry(plane: &Plane) -> SurfaceGeometry {
    SurfaceGeometry::Plane {
        origin: plane.0,
        u_axis: plane.1,
        v_axis: plane.2,
        normal: plane.3,
    }
}

/// The band recorded at the fan end: the strip face's and the feet's
/// incident edges' recorded tolerances, at least `pos_tol`.
pub(super) fn recorded_band(brep: &BRep, c: &CornerCandidate, pos_tol: f64) -> f64 {
    let feet = [
        brep.half_edges.get(c.he_in).map(|h| h.origin),
        brep.half_edges
            .get(c.he_out)
            .and_then(|h| brep.half_edges.get(h.twin))
            .map(|tw| tw.origin),
    ];
    let face = brep.faces.get(c.strip_face).and_then(|f| f.tolerance);
    let edges = brep
        .half_edges
        .values()
        .filter(|he| feet.contains(&Some(he.origin)))
        .filter_map(|he| brep.edges.get(he.edge)?.tolerance);
    face.into_iter().chain(edges).fold(pos_tol, f64::max)
}

/// March every arc of the chain: arc `j` on face surface `j`, from
/// `chain[j]` to `chain[j + 1]`.
pub(super) fn march_arcs(
    geom: &SurfaceGeometry,
    surfaces: &[SurfaceGeometry],
    chain: &[Point3],
    pos_tol: f64,
    recorded: f64,
) -> Result<Vec<MarchedIntersection>, Reason> {
    surfaces
        .iter()
        .enumerate()
        .map(|(j, surface)| {
            march_intersection(geom, surface, chain[j], chain[j + 1], pos_tol, recorded)
                .map_err(|_| Reason::ArcBranchUnresolved)
        })
        .collect()
}

/// Stage one marched arc on the strip half-edge `he` (running along
/// the chain) and its cap twin. The first arc seeds its strip chart
/// from the in-foot's neighbouring pcurve.
pub(super) fn stage_marched_arc(
    brep: &BRep,
    he: HalfEdgeHandle,
    geom: &SurfaceGeometry,
    cap_surface: &SurfaceGeometry,
    arc: &MarchedIntersection,
    first: bool,
) -> Result<StagedStubArc, Reason> {
    let h = brep.half_edges.get(he).ok_or(Reason::FanEndRejected)?;
    let (edge, twin) = (h.edge, h.twin);
    let strip_face = brep.loops.get(h.loop_).ok_or(Reason::FanEndRejected)?.face;
    let primary = brep
        .edges
        .get(edge)
        .ok_or(Reason::FanEndRejected)?
        .half_edge;
    let curve = intersection_curve(arc, primary != he).map_err(|_| Reason::FanEndRejected)?;
    let seed = if first {
        neighbour_chart_point(brep, strip_face, h.prev, true)
    } else {
        None
    };
    let strip_pcurve =
        intersection_pcurve(geom, arc, false, seed).map_err(|_| Reason::FanEndRejected)?;
    let cap_pcurve =
        intersection_pcurve(cap_surface, arc, true, None).map_err(|_| Reason::FanEndRejected)?;
    Ok(StagedStubArc {
        edge,
        curve,
        strip_half_edge: he,
        strip_pcurve,
        cap_half_edge: twin,
        cap_pcurve,
    })
}

/// A fitted strip ending on ONE oblique planar cap: the terminal iso
/// does not lie in the cap plane, so the swept circle fit declined.
/// The end is the strip ∩ cap intersection from foot to foot. The
/// corner settles onto it at the marched point nearest the stale
/// corner (the circular rule), and each stub becomes its marched arc.
pub(super) fn settle_single_cap(
    brep: &mut BRep,
    c: &CornerCandidate,
    geom: &SurfaceGeometry,
    tolerance: &popkern_core::tolerance::Tolerance,
) -> Result<(), Reason> {
    let pos_tol = (tolerance.linear() * 10.0).max(1e-9);
    let face_of = |b: &BRep, he: HalfEdgeHandle| -> Option<popkern_core::topology::FaceHandle> {
        Some(b.loops.get(b.half_edges.get(he)?.loop_)?.face)
    };
    let twin_in = brep
        .half_edges
        .get(c.he_in)
        .ok_or(Reason::FanEndRejected)?
        .twin;
    let twin_out = brep
        .half_edges
        .get(c.he_out)
        .ok_or(Reason::FanEndRejected)?
        .twin;
    let cap = face_of(brep, twin_in).ok_or(Reason::NoCommonPlanarCap)?;
    if face_of(brep, twin_out) != Some(cap) {
        return Err(Reason::NoCommonPlanarCap);
    }
    let plane = super::plane_of(brep, cap).ok_or(Reason::NoCommonPlanarCap)?;
    let position = |v| {
        brep.vertices
            .get(v)
            .map(|x: &popkern_core::topology::Vertex| x.position)
    };
    let foot_a = position(
        brep.half_edges
            .get(c.he_in)
            .ok_or(Reason::FanEndRejected)?
            .origin,
    )
    .ok_or(Reason::FanEndRejected)?;
    let foot_b = position(
        brep.half_edges
            .get(twin_out)
            .ok_or(Reason::FanEndRejected)?
            .origin,
    )
    .ok_or(Reason::FanEndRejected)?;
    let stale = position(c.corner).ok_or(Reason::FanEndRejected)?;
    let recorded = recorded_band(brep, c, pos_tol);
    let section = march_intersection(
        geom,
        &plane_geometry(&plane),
        foot_a,
        foot_b,
        pos_tol,
        recorded,
    )
    .map_err(|_| Reason::ArcBranchUnresolved)?;
    // Interior sample nearest the stale corner.
    let n = section.samples.len();
    let settled = section.samples[1..n - 1]
        .iter()
        .copied()
        .min_by(|a, b| a.distance_to(&stale).total_cmp(&b.distance_to(&stale)))
        .ok_or(Reason::ArcBranchUnresolved)?;
    let cap_surface = plane_geometry(&plane);
    let arcs = march_arcs(
        geom,
        &[cap_surface.clone(), cap_surface.clone()],
        &[foot_a, settled, foot_b],
        pos_tol,
        recorded,
    )?;
    if let Some(v) = brep.vertices.get_mut(c.corner) {
        v.position = settled;
    }
    let staged = [
        (
            stage_marched_arc(brep, c.he_in, geom, &cap_surface, &arcs[0], true)?,
            arcs[0].band,
        ),
        (
            stage_marched_arc(brep, c.he_out, geom, &cap_surface, &arcs[1], false)?,
            arcs[1].band,
        ),
    ];
    for (stub, band) in staged {
        super::commit_stub(brep, c.strip_face, stub, Some(band), pos_tol)?;
    }
    brep.validate().map_err(|_| Reason::FanEndRejected)?;
    Ok(())
}
