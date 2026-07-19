//! Terminal CHORD trim for planar blend strips (the chamfer
//! counterpart of the PK-148 terminal arc trim).
//!
//! A chamfer band that terminates at a plain body vertex keeps the
//! two support-edge *stubs* running through the original corner —
//! which sits off the chamfer plane by `d/√2` for a 90° equal
//! chamfer — while the band's own trim describes the straight
//! terminal chord ON the plane. Two geometric authorities on the
//! same edges: STEP export cannot co-parameterize the band-side
//! pcurves with the stub curves, forks face-local edge copies, and
//! the imported body fails Euler-Poincaré (measured: chamfered cube
//! round-trip χ=6 vs 2, four duplicate edges + eight coincident
//! vertices).
//!
//! **Repair, by construction:** the terminal cross-section of a
//! planar strip is the straight chord between the two on-plane
//! anchors. Reposition the corner vertex onto the chord (nearest
//! point to the stale corner), re-author both stub 3-D curves as
//! exact line segments, and write matching `Linear` pcurves on both
//! the strip side and the planar cap side. Same trust model,
//! transactional staging, and decline taxonomy as the arc trim.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle};

use crate::types::TerminalArcTrimDeclineReason;

use super::terminal_arc_author::{StagedStubArc, compute_stub_chord};
use super::terminal_arc_trim::{CornerCandidate, incident_edge_count};

/// A plane frame extracted once per repair (unit axes).
pub(super) struct StripPlane {
    pub(super) origin: Point3,
    pub(super) normal: Vec3,
    pub(super) u_axis: Vec3,
    pub(super) v_axis: Vec3,
}

fn plane_of(brep: &BRep, face: FaceHandle) -> Option<StripPlane> {
    match brep
        .faces
        .get(face)
        .and_then(|f| f.surface.as_ref())
        .map(|s| &s.geometry)
    {
        Some(SurfaceGeometry::Plane {
            origin,
            normal,
            u_axis,
            v_axis,
        }) => Some(StripPlane {
            origin: *origin,
            normal: *normal,
            u_axis: *u_axis,
            v_axis: *v_axis,
        }),
        _ => None,
    }
}

pub(super) fn repair_corner_planar(
    brep: &mut BRep,
    c: &CornerCandidate,
    tolerance: &Tolerance,
) -> Result<(), TerminalArcTrimDeclineReason> {
    let strip_plane =
        plane_of(brep, c.strip_face).ok_or(TerminalArcTrimDeclineReason::NonCylinderStrip)?;

    // Both twins must bound one shared planar cap face.
    let twin_face = |he: HalfEdgeHandle| -> Option<FaceHandle> {
        let tw = brep.half_edges.get(brep.half_edges.get(he)?.twin)?;
        let l = brep.loops.get(tw.loop_)?;
        Some(l.face)
    };
    let cap_face = twin_face(c.he_in).ok_or(TerminalArcTrimDeclineReason::NoCommonPlanarCap)?;
    if twin_face(c.he_out) != Some(cap_face) {
        return Err(TerminalArcTrimDeclineReason::NoCommonPlanarCap);
    }
    let cap_plane =
        plane_of(brep, cap_face).ok_or(TerminalArcTrimDeclineReason::NoCommonPlanarCap)?;

    // The corner vertex must be used by exactly the two stub edges.
    if incident_edge_count(brep, c.corner) != 2 {
        return Err(TerminalArcTrimDeclineReason::VertexNotIsolated);
    }

    // Trusted geometric inputs: the two anchor vertex positions.
    let anchor_of = |he: HalfEdgeHandle, at_origin: bool| -> Option<Point3> {
        let data = brep.half_edges.get(he)?;
        let vertex = if at_origin {
            data.origin
        } else {
            brep.half_edges.get(data.next)?.origin
        };
        Some(brep.vertices.get(vertex)?.position)
    };
    let anchor_in =
        anchor_of(c.he_in, true).ok_or(TerminalArcTrimDeclineReason::ArcBranchUnresolved)?;
    let anchor_out =
        anchor_of(c.he_out, false).ok_or(TerminalArcTrimDeclineReason::ArcBranchUnresolved)?;

    // Anchors must sit ON both the strip plane and the cap plane —
    // the chord between them is then their intersection line.
    let pos_tol = (tolerance.linear() * 10.0).max(1e-9);
    let trace = std::env::var_os("PK_TERMINAL_ARC_TRACE").is_some();
    for anchor in [anchor_in, anchor_out] {
        let off_strip = (anchor - strip_plane.origin).dot(&strip_plane.normal).abs();
        let off_cap = (anchor - cap_plane.origin).dot(&cap_plane.normal).abs();
        if off_strip > pos_tol || off_cap > pos_tol {
            if trace {
                eprintln!(
                    "[tat]   chord anchor {anchor:?} off strip {off_strip:.3e} / cap {off_cap:.3e}"
                );
            }
            return Err(TerminalArcTrimDeclineReason::AnchorsOffSurface);
        }
    }
    let chord = anchor_out - anchor_in;
    let chord_len = chord.length();
    if chord_len <= pos_tol {
        return Err(TerminalArcTrimDeclineReason::ArcBranchUnresolved);
    }

    // Settled corner: the chord point nearest the stale corner (its
    // projection, clamped into the span).
    let old_corner = brep
        .vertices
        .get(c.corner)
        .ok_or(TerminalArcTrimDeclineReason::ArcBranchUnresolved)?
        .position;
    let t = ((old_corner - anchor_in).dot(&chord) / (chord_len * chord_len)).clamp(0.0, 1.0);
    let new_corner = anchor_in + chord * t;
    // Sanity: the settled corner stays near the stale one (the stale
    // offset is d/√2 for a 90° equal chamfer).
    if (new_corner - old_corner).length() > chord_len {
        return Err(TerminalArcTrimDeclineReason::ArcBranchUnresolved);
    }

    // Compute every authored piece first; commit only when the whole
    // corner resolves (transactional).
    let mut staged: Vec<StagedStubArc> = Vec::with_capacity(2);
    for (he_strip, start_pos, end_pos) in [
        (c.he_in, anchor_in, new_corner),
        (c.he_out, new_corner, anchor_out),
    ] {
        staged.push(compute_stub_chord(
            brep,
            he_strip,
            &strip_plane,
            &cap_plane,
            start_pos,
            end_pos,
            tolerance,
        )?);
    }
    if let Some(v) = brep.vertices.get_mut(c.corner) {
        v.position = new_corner;
    }
    for stub in staged {
        if let Some(e) = brep.edges.get_mut(stub.edge) {
            e.curve = Some(stub.curve);
        }
        if let Some(he) = brep.half_edges.get_mut(stub.strip_half_edge) {
            he.pcurve = Some(stub.strip_pcurve.clone());
        }
        brep.set_face_use_pcurve(c.strip_face, stub.strip_half_edge, stub.strip_pcurve);
        if let Some(tw) = brep.half_edges.get_mut(stub.cap_half_edge) {
            tw.pcurve = Some(stub.cap_pcurve.clone());
        }
        brep.set_face_use_pcurve(cap_face, stub.cap_half_edge, stub.cap_pcurve);
    }
    Ok(())
}
