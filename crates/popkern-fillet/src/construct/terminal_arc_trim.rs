//! PK-148 terminal arc trim at open fillet chain ends.
//!
//! **Problem.** When a filleted edge terminates at a plain body
//! vertex (no meeting blend), Phase A splits the two support edges
//! at the tangent stations and wires the leftover *stubs* — straight
//! segments running through the original corner — as the strip's
//! terminal boundary. The strip's own trim describes the true
//! terminal cross-section arc, so the same two edges carry two
//! disagreeing geometric authorities: a line chain through the
//! corner (3-D curve + cap side) and an arc pair (strip side). The
//! corner sits `(√2−1)·r` off the arc; tessellation cracks along
//! the seam and STEP export forks face-local edges (measured: 8–12
//! OCC free closed edges on a plain single-edge cube fillet).
//!
//! **Repair, by construction.** The terminal cross-section of a
//! cylindrical strip at an axial station is an exact circular arc.
//! For each recognized corner: reposition the corner vertex onto
//! that arc, re-author both stub edges' 3-D curves as exact circle
//! arcs (`CurveGeometry::Ellipse` with explicit basis so no frame
//! convention is assumed), re-author the strip-side stub pcurves as
//! exact iso-v `Linear` spans, and write matching cap-side pcurves
//! (`Elliptical` on the planar end face). This *is* the end-face
//! trim: the cap loop now rounds the corner, and strip, cap, and
//! 3-D curves agree exactly. Zero topology change.
//!
//! **Trust model.** Pcurves are used only to *recognize* a stale
//! corner (two consecutive iso-v spans meeting at an off-surface
//! vertex). All repair *positions* come from trusted data — the two
//! anchor vertex positions (verified on-surface) and the cylinder
//! itself — because stale charts lie: an in-place radius modify
//! swaps the strip surface for a differently-parameterized cylinder
//! and the old pcurve stations map far off the body (measured:
//! station v=2.0 in the new chart landed a vertex at x=2 on a unit
//! cube).
//!
//! **Invariants preserved:** topology handles and counts; support
//! trims (stubs are shared exclusively by strip and cap — verified
//! before repair). **Fallback:** any structural doubt records a
//! [`TerminalArcTrimDecline`] and leaves the stale geometry as-is
//! (pre-PK-148 behaviour); the repair is transactional. **Tolerance:**
//! the caller's blend tolerance (`linear` for point checks) — local,
//! explicit, no global state.

use std::collections::HashSet;

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, VertexHandle};

use crate::types::{TerminalArcTrimDecline, TerminalArcTrimDeclineReason, TerminalArcTrimReport};

/// A recognized stale chain-end corner on one strip face.
pub(super) struct CornerCandidate {
    pub(super) strip_face: FaceHandle,
    /// Strip-loop half-edge arriving at the corner.
    pub(super) he_in: HalfEdgeHandle,
    /// Strip-loop half-edge leaving the corner (`he_in.next`).
    pub(super) he_out: HalfEdgeHandle,
    pub(super) corner: VertexHandle,
}

/// Settle stale chain-end corners on every face the blend operation
/// created. Returns the structured repair report.
pub(crate) fn settle_terminal_arc_stubs(
    source: &BRep,
    result: &mut BRep,
    tolerance: &Tolerance,
) -> TerminalArcTrimReport {
    let pre_faces: HashSet<FaceHandle> = source.faces.keys().collect();
    let created: Vec<FaceHandle> = result
        .faces
        .keys()
        .filter(|f| !pre_faces.contains(f))
        .collect();
    settle_terminal_arc_faces(result, &created, tolerance)
}

/// Settle stale chain-end corners on an explicit set of candidate
/// strip faces. Used by the fresh-blend pipeline (scoped to created
/// faces) and by in-place radius modification (scoped to the
/// modified strip, whose mid-arc corners must land on the NEW
/// terminal arc).
pub(crate) fn settle_terminal_arc_faces(
    result: &mut BRep,
    faces: &[FaceHandle],
    tolerance: &Tolerance,
) -> TerminalArcTrimReport {
    let mut report = TerminalArcTrimReport::default();
    for &face in faces {
        for candidate in recognize_corners(result, face, tolerance) {
            let outcome = repair_corner(result, &candidate, tolerance);
            if std::env::var_os("PK_TERMINAL_ARC_TRACE").is_some() {
                eprintln!(
                    "[tat] corner {:?} on {:?} (he {:?},{:?}) -> {:?}",
                    candidate.corner,
                    candidate.strip_face,
                    candidate.he_in,
                    candidate.he_out,
                    outcome
                );
            }
            match outcome {
                Ok(()) => report.repaired_corners += 1,
                Err(reason) => report.declines.push(TerminalArcTrimDecline {
                    strip_face: candidate.strip_face,
                    vertex: candidate.corner,
                    reason,
                }),
            }
        }
    }
    report
}

/// Walk one face's outer loop and recognize stale chain-end corners:
/// two consecutive half-edges whose strip pcurves sit on one iso-v
/// station while the shared vertex lies off the strip surface.
fn recognize_corners(brep: &BRep, face: FaceHandle, tolerance: &Tolerance) -> Vec<CornerCandidate> {
    let mut found = Vec::new();
    let Some(f) = brep.faces.get(face) else {
        return found;
    };
    let Some(surface) = f.surface.as_ref() else {
        return found;
    };
    let Some(l) = brep.loops.get(f.outer_loop) else {
        return found;
    };
    let start = l.half_edge;
    let mut he_in = start;
    // Bounded walk: loops are cycles; cap at half-edge count.
    for _ in 0..brep.half_edges.len() {
        let Some(a) = brep.half_edges.get(he_in) else {
            break;
        };
        let he_out = a.next;
        let Some(b) = brep.half_edges.get(he_out) else {
            break;
        };
        if let Some(candidate) = recognize_corner_pair(
            brep,
            face,
            &surface.geometry,
            he_in,
            he_out,
            b.origin,
            tolerance,
        ) {
            found.push(candidate);
        }
        he_in = he_out;
        if he_in == start {
            break;
        }
    }
    found
}

/// UV endpoints of a half-edge's strip-side pcurve, oriented along
/// the half-edge (origin → destination).
pub(super) fn strip_pcurve_span(
    brep: &BRep,
    face: FaceHandle,
    he: HalfEdgeHandle,
) -> Option<((f64, f64), (f64, f64))> {
    let pc = brep
        .face_use_pcurve(face, he)
        .or(brep.half_edges.get(he)?.pcurve.as_ref())?;
    let (t0, t1) = pc.parameter_range();
    let a = pc.point_at(t0);
    let b = pc.point_at(t1);
    Some(((a.x, a.y), (b.x, b.y)))
}

#[allow(clippy::too_many_arguments)]
fn recognize_corner_pair(
    brep: &BRep,
    face: FaceHandle,
    surface: &SurfaceGeometry,
    he_in: HalfEdgeHandle,
    he_out: HalfEdgeHandle,
    corner: VertexHandle,
    tolerance: &Tolerance,
) -> Option<CornerCandidate> {
    let (a_start, a_end) = strip_pcurve_span(brep, face, he_in)?;
    let (b_start, b_end) = strip_pcurve_span(brep, face, he_out)?;
    // Both stubs on one iso-v station, meeting at the same u. The
    // chart values may be stale — they gate recognition only; all
    // repair positions are computed from trusted geometry.
    // Recognition-only gate: swept-strip stub pcurves projected from
    // a stale corner drift off the exact iso by the corner's chart
    // shadow (measured 1.4e-3 in u on the r=0.10 end of a variable
    // sweep), so the station tolerance is loose; every ACCEPT is
    // still guarded downstream by off-surface, isolation, common-cap,
    // circularity, and anchor checks.
    let eps_uv = 5e-3;
    // Cylinder strips author their chain-end stubs on one iso-v
    // station (u meets); planar chamfer bands author them on one
    // iso-u station (v meets). Accept either — the chart values are
    // recognition-only either way.
    let station_v = a_end.1;
    let iso_v = (a_start.1 - station_v).abs() < eps_uv
        && (b_start.1 - station_v).abs() < eps_uv
        && (b_end.1 - station_v).abs() < eps_uv
        && (a_end.0 - b_start.0).abs() < eps_uv;
    let station_u = a_end.0;
    let iso_u = (a_start.0 - station_u).abs() < eps_uv
        && (b_start.0 - station_u).abs() < eps_uv
        && (b_end.0 - station_u).abs() < eps_uv
        && (a_end.1 - b_start.1).abs() < eps_uv;
    if !(iso_v || iso_u) {
        return None;
    }
    // The vertex must be OFF the strip surface — that is the stale
    // signature. Coherent corners (already settled) are left alone,
    // which keeps the pass idempotent. Measured against the true
    // support geometry, not the possibly-stale chart.
    let corner_pos = brep.vertices.get(corner)?.position;
    let off_surface = match surface {
        SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        } => {
            let d = corner_pos - *origin;
            let radial = d - *axis * d.dot(axis);
            (radial.length() - radius).abs() > tolerance.linear() * 10.0
        }
        SurfaceGeometry::Plane { origin, normal, .. } => {
            (corner_pos - *origin).dot(normal).abs() > tolerance.linear() * 10.0
        }
        SurfaceGeometry::Nurbs(_) => {
            // Swept (variable-radius) strips: the chart + surface are
            // the construction authority; the stale signature is the
            // corner vertex sitting off the surface point at the
            // stubs' shared chart corner.
            use popkern_core::geometry::SurfaceEval;
            let shared = a_end;
            let on_surface = surface.point_at(shared.0, shared.1);
            corner_pos.distance_to(&on_surface) > tolerance.linear() * 10.0
        }
        _ => return None,
    };
    if !off_surface {
        return None;
    }
    Some(CornerCandidate {
        strip_face: face,
        he_in,
        he_out,
        corner,
    })
}

/// Count edges incident to a vertex (via each edge's half-edge pair).
pub(super) fn incident_edge_count(brep: &BRep, vertex: VertexHandle) -> usize {
    brep.edges
        .iter()
        .filter(|(_, e)| {
            brep.half_edges.get(e.half_edge).is_some_and(|he| {
                he.origin == vertex
                    || brep
                        .half_edges
                        .get(he.twin)
                        .is_some_and(|tw| tw.origin == vertex)
            })
        })
        .count()
}

/// Cylinder data extracted once per repair.
pub(super) struct StripCylinder {
    pub(super) origin: Point3,
    pub(super) axis: Vec3,
    pub(super) radius: f64,
}

fn repair_corner(
    brep: &mut BRep,
    c: &CornerCandidate,
    tolerance: &Tolerance,
) -> Result<(), TerminalArcTrimDeclineReason> {
    // Cylindrical strips settle to the terminal ARC; planar strips
    // (chamfer bands) settle to the terminal CHORD.
    let cyl = {
        let geom = brep
            .faces
            .get(c.strip_face)
            .and_then(|f| f.surface.as_ref())
            .map(|s| &s.geometry);
        match geom {
            Some(SurfaceGeometry::Cylinder {
                origin,
                axis,
                radius,
            }) => StripCylinder {
                origin: *origin,
                axis: *axis,
                radius: *radius,
            },
            Some(SurfaceGeometry::Plane { .. }) => {
                return super::terminal_chord_trim::repair_corner_planar(brep, c, tolerance);
            }
            Some(SurfaceGeometry::Nurbs(_)) => {
                return super::terminal_arc_fit_trim::repair_corner_swept(brep, c, tolerance);
            }
            _ => return Err(TerminalArcTrimDeclineReason::NonCylinderStrip),
        }
    };
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
    let plane = match brep
        .faces
        .get(cap_face)
        .and_then(|f| f.surface.as_ref())
        .map(|s| &s.geometry)
    {
        Some(SurfaceGeometry::Plane {
            origin,
            u_axis,
            v_axis,
            ..
        }) => (*origin, *u_axis, *v_axis),
        _ => return Err(TerminalArcTrimDeclineReason::NoCommonPlanarCap),
    };
    // The corner vertex must be used by exactly the two stub edges.
    if incident_edge_count(brep, c.corner) != 2 {
        return Err(TerminalArcTrimDeclineReason::VertexNotIsolated);
    }

    // Trusted geometric inputs: the two anchor vertex positions.
    let anchor_in = {
        let he = brep
            .half_edges
            .get(c.he_in)
            .ok_or(TerminalArcTrimDeclineReason::ArcBranchUnresolved)?;
        brep.vertices
            .get(he.origin)
            .ok_or(TerminalArcTrimDeclineReason::ArcBranchUnresolved)?
            .position
    };
    let anchor_out = {
        let he = brep
            .half_edges
            .get(c.he_out)
            .ok_or(TerminalArcTrimDeclineReason::ArcBranchUnresolved)?;
        let next = brep
            .half_edges
            .get(he.next)
            .ok_or(TerminalArcTrimDeclineReason::ArcBranchUnresolved)?;
        brep.vertices
            .get(next.origin)
            .ok_or(TerminalArcTrimDeclineReason::ArcBranchUnresolved)?
            .position
    };
    // Anchors must sit ON the cylinder at a common axial station.
    let station_of = |p: Point3| -> f64 { (p - cyl.origin).dot(&cyl.axis) };
    let radial_of = |p: Point3| -> Vec3 {
        let d = p - cyl.origin;
        d - cyl.axis * d.dot(&cyl.axis)
    };
    let pos_tol = (tolerance.linear() * 10.0).max(1e-9);
    let trace = std::env::var_os("PK_TERMINAL_ARC_TRACE").is_some();
    for anchor in [anchor_in, anchor_out] {
        if (radial_of(anchor).length() - cyl.radius).abs() > pos_tol {
            if trace {
                eprintln!(
                    "[tat]   anchor {anchor:?} radial {} vs r {} (cyl origin {:?} axis {:?})",
                    radial_of(anchor).length(),
                    cyl.radius,
                    cyl.origin,
                    cyl.axis
                );
            }
            return Err(TerminalArcTrimDeclineReason::AnchorsOffSurface);
        }
    }
    let (sa, sb) = (station_of(anchor_in), station_of(anchor_out));
    if (sa - sb).abs() > pos_tol {
        if trace {
            eprintln!("[tat]   station mismatch {sa} vs {sb} (a {anchor_in:?} b {anchor_out:?})");
        }
        return Err(TerminalArcTrimDeclineReason::AnchorsOffSurface);
    }
    let station = 0.5 * (sa + sb);
    let center = cyl.origin + cyl.axis * station;

    // Settled corner: the arc point radially nearest the stale
    // corner (the tangent-line intersection). This is exact for the
    // convex chain-end configuration and never trusts the chart.
    let old_corner = brep
        .vertices
        .get(c.corner)
        .ok_or(TerminalArcTrimDeclineReason::ArcBranchUnresolved)?
        .position;
    let corner_radial = radial_of(old_corner);
    if corner_radial.length() < cyl.radius * 1e-6 {
        return Err(TerminalArcTrimDeclineReason::ArcBranchUnresolved);
    }
    let new_corner = center
        + corner_radial
            .normalized()
            .ok_or(TerminalArcTrimDeclineReason::ArcBranchUnresolved)?
            * cyl.radius;
    // Sanity: the settled corner stays near the stale one (the
    // stale offset is (√2−1)·r for a 90° corner; the antipodal
    // mis-pick would be ≈2.4·r away).
    if (new_corner - old_corner).length() > cyl.radius {
        return Err(TerminalArcTrimDeclineReason::ArcBranchUnresolved);
    }

    // Compute every authored piece first; commit only when the whole
    // corner resolves (transactional — a decline must not leave a
    // half-repaired body).
    let mut staged: Vec<super::terminal_arc_author::StagedStubArc> = Vec::with_capacity(2);
    for (he_strip, start_pos, end_pos) in [
        (c.he_in, anchor_in, new_corner),
        (c.he_out, new_corner, anchor_out),
    ] {
        staged.push(super::terminal_arc_author::compute_stub_arc(
            brep, he_strip, &cyl, center, start_pos, end_pos, plane, tolerance,
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
