//! Obstructed strip-end repair: fillet terminal arc interrupted by
//! an adjacent face at a planar end face.
//!
//! **Problem.** A straight-edge fillet terminating on a planar end
//! face normally settles to the exact cylinder ∩ plane terminal arc
//! between its two tangency feet ([`super::terminal_arc_trim`]).
//! When another face rises out of the end plane INSIDE the fillet's
//! footprint (an L-bracket wall overlapping the base within one
//! radius of the filleted corner), one tangency foot has no boundary
//! edge to land on: the support face wraps around the obstruction and
//! the true foot is INTERIOR to it. the anchor-split stage's perp-walk clamps the
//! anchor to the short end-face rim instead (measured on the corpus
//! `offset_bracket`: rail authored as a slant chord 0.47 off the
//! cylinder, terminal settle declining on `AnchorsOffSurface`, and
//! the stale corner edges serving at exactly `(√2−1)·r` off the strip
//! — the check-9 catch that kept the part off the corpus baseline).
//!
//! **Repair, by construction.** The correct local topology is fully
//! determined by the strip cylinder, the end plane, and the
//! obstructing face:
//!
//! 1. the TRUE foot `A` is the cylinder's tangency contact with the
//!    wrapped support plane at the terminal station (foot of the
//!    perpendicular from the arc center — exact);
//! 2. the terminal arc `foot → A` crosses the end face's boundary
//!    edge shared with the obstructing face at a crossing `C`
//!    (found by monotone bisection on the arc, then verified on the
//!    edge segment);
//! 3. the region of the end face between the arc, the obstructing
//!    face, and the support face becomes a NEW end-plane face (the
//!    obstruction's "ceiling": same plane, opposite orientation —
//!    material continues above it);
//! 4. the stale corner vertex becomes the crossing `C`; the clamped
//!    anchor vertex becomes the true foot `A`.
//!
//! Euler sequence (χ-neutral overall): `split_edge_at_point` on the
//! crossed boundary edge (V+1, E+1), `mef` corner→split-vertex
//! (E+1, F+1, L+1) carving the ceiling out of the end face, then
//! `kve` contracting the zero-length mef edge (V−1, E−1) so the
//! crossing is a single vertex.
//!
//! **Scope.** One obstruction crossing per corner (the single-wall
//! class); multi-crossing configurations decline with a recorded
//! reason and keep the stale-but-valid geometry. **Transactional:**
//! the whole surgery runs on a candidate clone and commits only when
//! `validate()` holds. **Tolerance:** the caller's blend tolerance;
//! the tangency gate reuses the terminal-arc anchor gate
//! (`linear × 10`).

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, VertexHandle};

use super::terminal_arc_trim::{CornerCandidate, StripCylinder, incident_edge_count};
use super::terminal_obstructed_crossing::find_arc_crossing;
use crate::types::TerminalArcTrimDeclineReason;

/// Attempt the obstructed-end surgery for a corner whose plain
/// settle failed on an off-surface anchor. Returns `Ok(())` when the
/// repair commits; any decline keeps `brep` untouched.
pub(super) fn try_repair_obstructed_end(
    brep: &mut BRep,
    c: &CornerCandidate,
    cyl: &StripCylinder,
    tolerance: &Tolerance,
) -> Result<(), TerminalArcTrimDeclineReason> {
    let trace = popkern_core::trace::trace_enabled("terminal-arc");
    let pos_tol = (tolerance.linear() * 10.0).max(1e-9);

    let radial_of = |p: Point3| -> Vec3 {
        let d = p - cyl.origin;
        d - cyl.axis * d.dot(&cyl.axis)
    };
    let on_cyl = |p: Point3| (radial_of(p).length() - cyl.radius).abs() <= pos_tol;

    // Stub endpoint vertices away from the corner.
    let v_a_in = brep
        .half_edges
        .get(c.he_in)
        .map(|he| he.origin)
        .ok_or(TerminalArcTrimDeclineReason::ArcBranchUnresolved)?;
    let v_a_out = brep
        .half_edges
        .get(c.he_out)
        .and_then(|he| brep.half_edges.get(he.next))
        .map(|he| he.origin)
        .ok_or(TerminalArcTrimDeclineReason::ArcBranchUnresolved)?;
    let pos_in = brep.vertices[v_a_in].position;
    let pos_out = brep.vertices[v_a_out].position;

    // Exactly one anchor off the cylinder = the obstructed side.
    let (foot_v, foot_pos, off_v, rail_he) = match (on_cyl(pos_in), on_cyl(pos_out)) {
        (true, false) => {
            // Rail leaves the strip loop after the off anchor.
            let rail = brep
                .half_edges
                .get(c.he_out)
                .map(|he| he.next)
                .ok_or(TerminalArcTrimDeclineReason::ArcBranchUnresolved)?;
            (v_a_in, pos_in, v_a_out, rail)
        }
        (false, true) => {
            let rail = brep
                .half_edges
                .get(c.he_in)
                .map(|he| he.prev)
                .ok_or(TerminalArcTrimDeclineReason::ArcBranchUnresolved)?;
            (v_a_out, pos_out, v_a_in, rail)
        }
        _ => return Err(TerminalArcTrimDeclineReason::AnchorsOffSurface),
    };
    let _ = foot_v;

    // Terminal station and arc center from the trusted on-surface foot.
    let station = (foot_pos - cyl.origin).dot(&cyl.axis);
    let center = cyl.origin + cyl.axis * station;

    // Support face across the rail: its plane carries the true foot.
    let support_face =
        twin_face(brep, rail_he).ok_or(TerminalArcTrimDeclineReason::ObstructedEnd)?;
    let (s_origin, s_normal) =
        plane_of(brep, support_face).ok_or(TerminalArcTrimDeclineReason::ObstructedEnd)?;
    let n = s_normal
        .normalized()
        .ok_or(TerminalArcTrimDeclineReason::ObstructedEnd)?;
    // Tangent support of a straight-edge fillet: plane parallel to
    // the axis, one radius from it.
    if n.dot(&cyl.axis).abs() > 1e-6 {
        return Err(TerminalArcTrimDeclineReason::ObstructedEnd);
    }
    let signed = (center - s_origin).dot(&n);
    if (signed.abs() - cyl.radius).abs() > pos_tol {
        if trace {
            eprintln!(
                "[tat-obs] support-plane distance {} vs r {} — not tangent",
                signed.abs(),
                cyl.radius
            );
        }
        return Err(TerminalArcTrimDeclineReason::ObstructedEnd);
    }
    let foot_a = center - n * signed;

    // Both stubs must border one shared planar cap (the end face).
    let cap_face =
        twin_face(brep, c.he_in).ok_or(TerminalArcTrimDeclineReason::NoCommonPlanarCap)?;
    if twin_face(brep, c.he_out) != Some(cap_face) {
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
    // The off anchor joins exactly stub + rail + one cap-boundary
    // edge (the step that will bound the ceiling).
    if incident_edge_count(brep, off_v) != 3 {
        return Err(TerminalArcTrimDeclineReason::ObstructedEnd);
    }

    // Arc basis foot → A (minor sweep).
    let u_basis = foot_pos - center;
    let mut v_basis = cyl.axis.cross(&u_basis);
    let d = foot_a - center;
    let mut sweep = d.dot(&v_basis).atan2(d.dot(&u_basis));
    if sweep < 0.0 {
        v_basis = -v_basis;
        sweep = -sweep;
    }
    if !(1e-9..=std::f64::consts::PI - 1e-9).contains(&sweep) {
        return Err(TerminalArcTrimDeclineReason::ObstructedEnd);
    }
    let u_hat = match u_basis.normalized() {
        Some(u) => u,
        None => return Err(TerminalArcTrimDeclineReason::ObstructedEnd),
    };
    let v_hat = match v_basis.normalized() {
        Some(v) => v,
        None => return Err(TerminalArcTrimDeclineReason::ObstructedEnd),
    };
    let arc_at = |theta: f64| -> Point3 {
        center + (u_hat * theta.cos() + v_hat * theta.sin()) * cyl.radius
    };

    // Find the single crossing of the open arc with the cap
    // boundary (straight edges not incident to the stub chain).
    let crossing = find_arc_crossing(brep, cap_face, c, off_v, &arc_at, sweep, pos_tol)?;
    let (edge_x, cross_pt) = crossing;

    // ---- Mutation on a candidate clone; validate() gates commit ----
    let mut cand = brep.clone();

    let split = cand
        .split_edge_at_point(edge_x, cross_pt, tolerance)
        .map_err(|e| {
            if trace {
                eprintln!("[tat-obs] boundary split failed: {e:?}");
            }
            TerminalArcTrimDeclineReason::ObstructedEnd
        })?;
    let vc = split.vertex;
    if let Some(v) = cand.vertices.get_mut(vc) {
        v.position = cross_pt;
    }

    // True positions: stale corner → crossing, clamped anchor →
    // foot. `snap_vertex_and_adjacent_chords` moves the vertex and
    // rebuilds every incident chord Line from the new positions (the
    // clamp offset can be a large move — the small-move robust snap
    // is not enough). Stub curves are re-authored as exact arcs just
    // below, overwriting the interim chords.
    let stub_in_edge = cand
        .half_edges
        .get(c.he_in)
        .map(|he| he.edge)
        .ok_or(TerminalArcTrimDeclineReason::ArcBranchUnresolved)?;
    let stub_out_edge = cand
        .half_edges
        .get(c.he_out)
        .map(|he| he.edge)
        .ok_or(TerminalArcTrimDeclineReason::ArcBranchUnresolved)?;
    crate::construct::edge_reconcile::snap_vertex_and_adjacent_chords(
        &mut cand,
        c.corner,
        cross_pt,
        stub_in_edge,
    );
    crate::construct::edge_reconcile::snap_vertex_and_adjacent_chords(
        &mut cand,
        off_v,
        foot_a,
        if on_cyl(pos_in) {
            stub_out_edge
        } else {
            stub_in_edge
        },
    );

    // Author both stub edges as exact terminal-arc segments (strip
    // iso pcurve + cap-plane elliptical pcurve). Positions run in
    // half-edge direction: he_in arrives at the corner (now the
    // crossing), he_out leaves it.
    let stubs = if on_cyl(pos_in) {
        [(c.he_in, foot_pos, cross_pt), (c.he_out, cross_pt, foot_a)]
    } else {
        [(c.he_in, foot_a, cross_pt), (c.he_out, cross_pt, pos_out)]
    };
    let mut staged = Vec::with_capacity(2);
    for (he_strip, sp, ep) in stubs {
        staged.push(super::terminal_arc_author::compute_stub_arc(
            &cand, he_strip, cyl, center, sp, ep, plane, tolerance,
        )?);
    }
    for stub in &staged {
        if let Some(e) = cand.edges.get_mut(stub.edge) {
            e.curve = Some(stub.curve.clone());
        }
        if let Some(he) = cand.half_edges.get_mut(stub.strip_half_edge) {
            he.pcurve = Some(stub.strip_pcurve.clone());
        }
        cand.set_face_use_pcurve(
            c.strip_face,
            stub.strip_half_edge,
            stub.strip_pcurve.clone(),
        );
        if let Some(tw) = cand.half_edges.get_mut(stub.cap_half_edge) {
            tw.pcurve = Some(stub.cap_pcurve.clone());
        }
    }

    // Carve the ceiling out of the cap: mef corner → vc (both at the
    // crossing position), then contract the zero-length edge.
    let h_c = cap_loop_he_with_origin(&cand, cap_face, c.corner)
        .ok_or(TerminalArcTrimDeclineReason::ObstructedEnd)?;
    let h_vc = cap_loop_he_with_origin(&cand, cap_face, vc)
        .ok_or(TerminalArcTrimDeclineReason::ObstructedEnd)?;
    let mef = cand.mef(h_c, h_vc).map_err(|e| {
        if trace {
            eprintln!("[tat-obs] mef failed: {e:?}");
        }
        TerminalArcTrimDeclineReason::ObstructedEnd
    })?;
    let ceiling = if loop_contains_vertex(&cand, cand.faces[mef.face].outer_loop, off_v) {
        mef.face
    } else {
        cap_face
    };
    // Ceiling: same end plane, opposite side of the material.
    let cap_surface = cand
        .faces
        .get(cap_face)
        .and_then(|f| f.surface.clone())
        .ok_or(TerminalArcTrimDeclineReason::NoCommonPlanarCap)?;
    let cap_orientation = cap_surface.orientation;
    if let Some(f) = cand.faces.get_mut(ceiling) {
        f.surface = Some(popkern_core::geometry::surface::SurfaceOnFace {
            geometry: cap_surface.geometry.clone(),
            orientation: if ceiling == cap_face {
                cap_orientation
            } else {
                !cap_orientation
            },
        });
    }
    if ceiling == cap_face {
        // The parent kept the ceiling side: the new face is the
        // remaining end face and inherits the parent's orientation;
        // flip the parent instead.
        if let Some(f) = cand.faces.get_mut(cap_face)
            && let Some(s) = f.surface.as_mut()
        {
            s.orientation = !cap_orientation;
        }
        if let Some(f) = cand.faces.get_mut(mef.face) {
            f.surface = Some(popkern_core::geometry::surface::SurfaceOnFace {
                geometry: cap_surface.geometry,
                orientation: cap_orientation,
            });
        }
    }
    cand.kve(mef.half_edge).map_err(|e| {
        if trace {
            eprintln!("[tat-obs] kve contraction failed: {e:?}");
        }
        TerminalArcTrimDeclineReason::ObstructedEnd
    })?;

    // The stub twin on the obstructed side now bounds the ceiling —
    // re-register its elliptical pcurve there (same plane chart).
    for stub in &staged {
        let twf = twin_face(&cand, stub.strip_half_edge);
        if let Some(f) = twf {
            cand.set_face_use_pcurve(f, stub.cap_half_edge, stub.cap_pcurve.clone());
        }
    }
    // Refresh projected pcurves on every touched face (analytic
    // pcurves written above survive the recompute guards).
    let obstructing_face = twin_face(&cand, cand.edges[edge_x].half_edge);
    for f in [
        Some(cap_face),
        Some(mef.face),
        Some(support_face),
        obstructing_face,
    ]
    .into_iter()
    .flatten()
    {
        let _ = cand.recompute_face_pcurves(f);
    }

    if let Err(e) = cand.validate() {
        if trace {
            eprintln!("[tat-obs] candidate failed validation: {e:?}");
        }
        return Err(TerminalArcTrimDeclineReason::ObstructedEnd);
    }
    if trace {
        eprintln!("[tat-obs] committed: foot {foot_a:?} crossing {cross_pt:?} ceiling {ceiling:?}");
    }
    *brep = cand;
    Ok(())
}

fn twin_face(brep: &BRep, he: HalfEdgeHandle) -> Option<FaceHandle> {
    let tw = brep.half_edges.get(brep.half_edges.get(he)?.twin)?;
    Some(brep.loops.get(tw.loop_)?.face)
}

fn plane_of(brep: &BRep, face: FaceHandle) -> Option<(Point3, Vec3)> {
    match brep
        .faces
        .get(face)
        .and_then(|f| f.surface.as_ref())
        .map(|s| &s.geometry)
    {
        Some(SurfaceGeometry::Plane { origin, normal, .. }) => Some((*origin, *normal)),
        _ => None,
    }
}

fn loop_contains_vertex(
    brep: &BRep,
    loop_h: popkern_core::topology::LoopHandle,
    v: VertexHandle,
) -> bool {
    let Some(l) = brep.loops.get(loop_h) else {
        return false;
    };
    let start = l.half_edge;
    let mut cur = start;
    for _ in 0..brep.half_edges.len() + 1 {
        let Some(he) = brep.half_edges.get(cur) else {
            return false;
        };
        if he.origin == v {
            return true;
        }
        cur = he.next;
        if cur == start {
            break;
        }
    }
    false
}

/// Find the half-edge in any loop of `face` whose origin is `v`.
fn cap_loop_he_with_origin(
    brep: &BRep,
    face: FaceHandle,
    v: VertexHandle,
) -> Option<HalfEdgeHandle> {
    let f = brep.faces.get(face)?;
    let loops = std::iter::once(f.outer_loop).chain(f.inner_loops.iter().copied());
    for lh in loops {
        let l = brep.loops.get(lh)?;
        let start = l.half_edge;
        let mut cur = start;
        for _ in 0..brep.half_edges.len() + 1 {
            let he = brep.half_edges.get(cur)?;
            if he.origin == v {
                return Some(cur);
            }
            cur = he.next;
            if cur == start {
                break;
            }
        }
    }
    None
}
