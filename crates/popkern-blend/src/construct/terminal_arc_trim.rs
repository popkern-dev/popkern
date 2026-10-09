//! PK-148 terminal arc trim at open fillet chain ends.
//!
//! **Problem.** When a filleted edge terminates at a plain body
//! vertex (no meeting blend), the anchor-split stage splits the two support edges
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

/// How a PLANAR created face's stubs are recognized. Chamfer bands
/// are planar STRIPS: an off-plane loop corner is a defect by
/// definition, including skewed stubs (bevel meeting its end face
/// obliquely on a yawed+tilted boolean-cut rest plane — PK-168), so
/// recognition is surface-scoped. Fillet pipelines also CREATE
/// planar faces (supports, bridges, cap sectors) whose corners
/// legitimately fail chart gates mid-construction; loosening there
/// mangled the V-free serving (measured: cube_mixed
/// pcurve_pair_count 6 -> 0), so the fillet entry stays chart-gated.
#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum PlanarStubRecognition {
    /// Chart-direction gates apply to planar stubs (fillet entry).
    ChartGated,
    /// Planar strips recognize on the off-plane signature alone
    /// (chamfer entry).
    SurfaceScoped,
}

/// Settle stale chain-end corners on every face the blend operation
/// created. Returns the structured repair report.
pub(crate) fn settle_terminal_arc_stubs(
    source: &BRep,
    result: &mut BRep,
    tolerance: &Tolerance,
    planar_mode: PlanarStubRecognition,
) -> TerminalArcTrimReport {
    let pre_faces: HashSet<FaceHandle> = source.faces.keys().collect();
    let created: Vec<FaceHandle> = result
        .faces
        .keys()
        .filter(|f| !pre_faces.contains(f))
        .collect();
    settle_terminal_faces_with_mode(result, &created, tolerance, planar_mode)
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
    settle_terminal_faces_with_mode(result, faces, tolerance, PlanarStubRecognition::ChartGated)
}

use super::terminal_corner_recognition::recognize_corners;

fn settle_terminal_faces_with_mode(
    result: &mut BRep,
    faces: &[FaceHandle],
    tolerance: &Tolerance,
    planar_mode: PlanarStubRecognition,
) -> TerminalArcTrimReport {
    let mut report = TerminalArcTrimReport::default();
    for &face in faces {
        let candidates = recognize_corners(result, face, tolerance, planar_mode);
        if popkern_core::trace::trace_enabled("terminal-arc") {
            eprintln!(
                "[tat] face {face:?}: {} corner candidates",
                candidates.len()
            );
        }
        for candidate in candidates {
            let outcome = repair_corner(result, &candidate, tolerance);
            if popkern_core::trace::trace_enabled("terminal-arc") {
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
#[derive(Clone)]
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
            Some(SurfaceGeometry::Cone {
                apex,
                axis,
                half_angle,
            }) => {
                // Circular chamfers on interrupted rims produce CONE
                // strips; their terminal boundary against a planar
                // cap is the cone∩plane conic (PK-171: the notched
                // flange's stale rim-end corner shipped a 1.0 lie
                // because Cone declined here as NonCylinderStrip).
                let strip = super::terminal_cone_trim::StripCone {
                    apex: *apex,
                    axis: *axis,
                    half_angle: *half_angle,
                };
                return super::terminal_cone_trim::repair_corner_conic(brep, c, &strip, tolerance);
            }
            Some(SurfaceGeometry::Nurbs(_)) => {
                return super::terminal_arc_fit_trim::repair_corner_swept(brep, c, tolerance);
            }
            // A torus strip on a meridian end plane: its local section
            // cylinder (PK-769).
            Some(SurfaceGeometry::Torus { .. }) => {
                super::terminal_torus_section::meridian_section_cylinder(brep, c, tolerance)?
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
    let cap_face = match (twin_face(c.he_in), twin_face(c.he_out)) {
        // Both stubs on the sheet boundary: a free end (PK-530).
        (None, None) => {
            return super::terminal_arc_oblique::repair_free_end(brep, c, &cyl, tolerance);
        }
        (Some(face), _) => face,
        (None, Some(_)) => return Err(TerminalArcTrimDeclineReason::NoCommonPlanarCap),
    };
    if twin_face(c.he_out) != Some(cap_face) {
        // The strip ends into a fan of faces (a slant edge running into
        // a pyramid apex): one exact arc per fan face (PK-530).
        return super::terminal_fan_end::repair_fan_end(
            brep,
            c,
            super::terminal_fan_end::FanStrip::Cylinder(&cyl),
            tolerance,
        );
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
            normal,
        }) => (*origin, *u_axis, *v_axis, *normal),
        _ => return Err(TerminalArcTrimDeclineReason::NoCommonPlanarCap),
    };
    let (plane, plane_normal) = ((plane.0, plane.1, plane.2), plane.3);
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
    let trace = popkern_core::trace::trace_enabled("terminal-arc");
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
            // One off-surface anchor with the other exact is the
            // signature of an OBSTRUCTED end (an adjacent face rises
            // inside the fillet's terminal footprint and the
            // perp-walk clamped the anchor short) — attempt the
            // topology surgery before declining.
            return super::terminal_obstructed_end::try_repair_obstructed_end(
                brep, c, &cyl, tolerance,
            );
        }
    }
    let old_corner = brep
        .vertices
        .get(c.corner)
        .ok_or(TerminalArcTrimDeclineReason::ArcBranchUnresolved)?
        .position;
    let (sa, sb) = (station_of(anchor_in), station_of(anchor_out));
    // A cap not perpendicular to the axis cuts an ellipse, even when
    // symmetry puts both feet at one station (a pyramid slant edge
    // running out on the base: equal stations, cap pcurve 1.3e-2 off
    // under the circular assumption).
    let oblique_cap = plane_normal
        .normalized()
        .zip(cyl.axis.normalized())
        .is_none_or(|(n, a)| n.cross(&a).length() > 1e-9);
    let (new_corner, staged) = if oblique_cap || (sa - sb).abs() > pos_tol {
        // Feet at different stations or an oblique cap: the section
        // is an ellipse (PK-530).
        super::terminal_arc_oblique::stage_oblique_end(
            brep,
            c,
            &cyl,
            [anchor_in, old_corner, anchor_out],
            (plane, plane_normal),
            tolerance,
        )?
    } else {
        let station = 0.5 * (sa + sb);
        let center = cyl.origin + cyl.axis * station;

        // Settled corner: the arc point radially nearest the stale
        // corner (the tangent-line intersection). This is exact for
        // the convex chain-end configuration and never trusts the
        // chart.
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

        // Compute every authored piece first; commit only when the
        // whole corner resolves (transactional — a decline must not
        // leave a half-repaired body).
        let mut staged: Vec<super::terminal_arc_author::StagedStubArc> = Vec::with_capacity(2);
        for (he_strip, start_pos, end_pos) in [
            (c.he_in, anchor_in, new_corner),
            (c.he_out, new_corner, anchor_out),
        ] {
            staged.push(super::terminal_arc_author::compute_stub_arc(
                brep, he_strip, &cyl, center, start_pos, end_pos, plane, tolerance,
            )?);
        }
        (new_corner, staged)
    };
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
