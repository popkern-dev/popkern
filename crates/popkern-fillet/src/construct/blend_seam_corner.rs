//! PK-152 blend×blend corner-seam settlement.
//!
//! **Problem.** Two fillets on adjacent edges meeting at a shared
//! body corner (third edge sharp) are constructed independently, so
//! near the corner each strip runs past the other and the support
//! plane keeps a phantom flap face: a small quad through the original
//! corner whose boundary edges sit `(√2−1)·r` off both strips. The
//! terminal settlement family correctly declines these corners
//! (`VertexNotIsolated`), the body validates, and the defect ships as
//! welded-mesh opens and off-surface corner edges.
//!
//! **Repair, by construction.** The strips' true mutual boundary is
//! their surface×surface intersection. For equal-radius cylinder
//! strips with intersecting axes this is a planar ellipse in a
//! bisector plane of the axes, and both rail-crossing anchor
//! vertices of the corner already lie on it. The settlement:
//!
//! 1. recognizes the corner tangle — a stale edge shared by the two
//!    strips whose far endpoint is off both surfaces, flanked by two
//!    stub chains that enclose one planar flap face between the
//!    trusted anchors;
//! 2. retires the flap and chains with Euler operators only (`kef`
//!    merges the flap into one strip, `kev` dissolves the spur
//!    chain, `kve` contracts the far chain onto the far anchor);
//! 3. re-authors the surviving shared edge as the exact intersection
//!    ellipse arc with exact `Sinusoidal` pcurves in both strip
//!    charts (see `blend_seam_corner_geometry`).
//!
//! **Invariants:** Euler-neutral surgery (χ unchanged), strips keep
//! their face identity, the seam edge records both strips as its
//! faces, support rails and caps are untouched. **Fallback:** every
//! decline restores the pre-corner snapshot and records a
//! [`BlendSeamCornerDecline`] — the stale tangle is never made
//! worse. **Tolerance:** the caller's blend tolerance, local and
//! explicit.

use std::collections::HashSet;

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

use crate::types::{BlendSeamCornerDecline, BlendSeamCornerDeclineReason, BlendSeamCornerReport};

use super::blend_seam_corner_geometry as geometry;
use super::terminal_arc_trim::StripCylinder;

/// Longest stub chain the recognizer will walk on one strip.
const MAX_CHAIN: usize = 6;

/// One recognized pair-corner tangle.
struct PairCornerCandidate {
    /// Strip whose seam half-edge points AT the corner; receives the
    /// flap merge (next-walk side).
    strip_n: FaceHandle,
    /// Strip whose seam half-edge leaves the corner (prev-walk side).
    strip_p: FaceHandle,
    /// The stale strip-shared edge; survives as the seam.
    seam_edge: EdgeHandle,
    /// Seam half-edge in `strip_n`'s loop (destination = corner).
    he_n: HalfEdgeHandle,
    /// Seam half-edge in `strip_p`'s loop (origin = corner).
    he_p: HalfEdgeHandle,
    /// The stale corner vertex (dies).
    corner: VertexHandle,
    /// Seam anchor already shared by the strips (`he_n`'s origin).
    anchor_near: VertexHandle,
    /// Far seam anchor where both stub chains meet (survives).
    anchor_far: VertexHandle,
    /// Stub half-edges in `strip_n`'s loop, corner → far anchor.
    chain_n: Vec<HalfEdgeHandle>,
    /// Stub half-edges in `strip_p`'s loop, each oriented toward the
    /// corner, ordered corner-side first.
    chain_p: Vec<HalfEdgeHandle>,
    /// The phantom flap face enclosed by the chains (dies).
    flap: FaceHandle,
}

/// Settle stale pair corners between blend strips the operation
/// created. Returns the structured repair report.
pub(crate) fn settle_blend_seam_corners(
    source: &BRep,
    result: &mut BRep,
    tolerance: &Tolerance,
) -> BlendSeamCornerReport {
    let pre: HashSet<FaceHandle> = source.faces.keys().collect();
    let created: HashSet<FaceHandle> = result.faces.keys().filter(|f| !pre.contains(f)).collect();
    let mut report = BlendSeamCornerReport::default();
    let trace = std::env::var_os("PK_BLEND_SEAM_TRACE").is_some();
    for candidate in recognize_pair_corners(result, &created, tolerance) {
        // Corners are structurally disjoint, but stay defensive
        // across repairs: every handle must still be alive.
        if result.edges.get(candidate.seam_edge).is_none()
            || result.faces.get(candidate.flap).is_none()
        {
            continue;
        }
        let snapshot = result.clone();
        let outcome = repair_pair_corner(result, &candidate, tolerance);
        if trace {
            eprintln!(
                "[bsc] corner {:?} strips {:?}/{:?} seam {:?} -> {:?}",
                candidate.corner,
                candidate.strip_n,
                candidate.strip_p,
                candidate.seam_edge,
                outcome
            );
        }
        match outcome {
            Ok(()) => report.seams_built += 1,
            Err(reason) => {
                *result = snapshot;
                report.declines.push(BlendSeamCornerDecline {
                    strips: (candidate.strip_n, candidate.strip_p),
                    corner: candidate.corner,
                    reason,
                });
            }
        }
    }
    report
}

fn strip_cylinder(brep: &BRep, face: FaceHandle) -> Option<StripCylinder> {
    match brep.faces.get(face)?.surface.as_ref().map(|s| &s.geometry) {
        Some(SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        }) => Some(StripCylinder {
            origin: *origin,
            axis: *axis,
            radius: *radius,
        }),
        _ => None,
    }
}

fn off_cylinder(p: Point3, cyl: &StripCylinder) -> f64 {
    let d = p - cyl.origin;
    let radial: Vec3 = d - cyl.axis * d.dot(&cyl.axis);
    (radial.length() - cyl.radius).abs()
}

/// Recognize every pair-corner tangle among the created strip faces.
fn recognize_pair_corners(
    brep: &BRep,
    created: &HashSet<FaceHandle>,
    tolerance: &Tolerance,
) -> Vec<PairCornerCandidate> {
    let stale_tol = (tolerance.linear() * 10.0).max(1e-9);
    let mut found = Vec::new();
    for (edge, e) in brep.edges.iter() {
        let Some(he_a) = brep.half_edges.get(e.half_edge) else {
            continue;
        };
        let Some(he_b) = brep.half_edges.get(he_a.twin) else {
            continue;
        };
        let (Some(la), Some(lb)) = (brep.loops.get(he_a.loop_), brep.loops.get(he_b.loop_)) else {
            continue;
        };
        if la.face == lb.face || !created.contains(&la.face) || !created.contains(&lb.face) {
            continue;
        }
        let (Some(cyl_a), Some(cyl_b)) =
            (strip_cylinder(brep, la.face), strip_cylinder(brep, lb.face))
        else {
            continue;
        };
        // The stale corner endpoint is off BOTH strips; the near
        // anchor is on at least one (on both when radii match; on
        // its own strip only across a radius mismatch, which the
        // repair then declines as the structured stage-2 reason).
        let off_all = |v: VertexHandle| -> Option<bool> {
            let p = brep.vertices.get(v)?.position;
            Some(off_cylinder(p, &cyl_a) > stale_tol && off_cylinder(p, &cyl_b) > stale_tol)
        };
        let va = he_a.origin;
        let vb = he_b.origin;
        let (Some(va_off), Some(vb_off)) = (off_all(va), off_all(vb)) else {
            continue;
        };
        let (he_n, he_p, corner, anchor_near) = match (va_off, vb_off) {
            (false, true) => (e.half_edge, he_a.twin, vb, va),
            (true, false) => (he_a.twin, e.half_edge, va, vb),
            _ => continue,
        };
        let Some(candidate) = walk_corner_tangle(
            brep,
            edge,
            he_n,
            he_p,
            corner,
            anchor_near,
            &cyl_a,
            &cyl_b,
            stale_tol,
        ) else {
            continue;
        };
        found.push(candidate);
    }
    found
}

/// Walk both stub chains from the corner and verify they enclose one
/// planar flap face between the trusted anchors.
#[allow(clippy::too_many_arguments)]
fn walk_corner_tangle(
    brep: &BRep,
    seam_edge: EdgeHandle,
    he_n: HalfEdgeHandle,
    he_p: HalfEdgeHandle,
    corner: VertexHandle,
    anchor_near: VertexHandle,
    cyl_a: &StripCylinder,
    cyl_b: &StripCylinder,
    stale_tol: f64,
) -> Option<PairCornerCandidate> {
    let on_both = |v: VertexHandle| -> Option<bool> {
        let p = brep.vertices.get(v)?.position;
        Some(off_cylinder(p, cyl_a) <= stale_tol && off_cylinder(p, cyl_b) <= stale_tol)
    };
    let twin_face = |he: HalfEdgeHandle| -> Option<FaceHandle> {
        let tw = brep.half_edges.get(brep.half_edges.get(he)?.twin)?;
        Some(brep.loops.get(tw.loop_)?.face)
    };
    // Next-walk on the strip whose seam half-edge arrives at the
    // corner: half-edges oriented away from the corner.
    let mut chain_n = Vec::new();
    let mut far_n = None;
    let mut he = brep.half_edges.get(he_n)?.next;
    for _ in 0..MAX_CHAIN {
        let h = brep.half_edges.get(he)?;
        let dest = brep.half_edges.get(h.twin)?.origin;
        chain_n.push(he);
        if on_both(dest)? {
            far_n = Some(dest);
            break;
        }
        he = h.next;
    }
    // Prev-walk on the other strip: half-edges oriented toward the
    // corner, corner-side first.
    let mut chain_p = Vec::new();
    let mut far_p = None;
    let mut he = brep.half_edges.get(he_p)?.prev;
    for _ in 0..MAX_CHAIN {
        let h = brep.half_edges.get(he)?;
        chain_p.push(he);
        if on_both(h.origin)? {
            far_p = Some(h.origin);
            break;
        }
        he = h.prev;
    }
    let anchor_far = far_n?;
    if far_p? != anchor_far || anchor_far == anchor_near || chain_n.is_empty() {
        return None;
    }
    // Every chain half-edge must border one common flap face.
    let strip_n = brep.loops.get(brep.half_edges.get(he_n)?.loop_)?.face;
    let strip_p = brep.loops.get(brep.half_edges.get(he_p)?.loop_)?.face;
    let flap = twin_face(*chain_n.first()?)?;
    if flap == strip_n || flap == strip_p {
        return None;
    }
    for he in chain_n.iter().chain(chain_p.iter()) {
        if twin_face(*he) != Some(flap) {
            return None;
        }
    }
    // The flap must be exactly the region the chains enclose: planar,
    // no inner loops, boundary = the chain edges.
    let f = brep.faces.get(flap)?;
    if !f.inner_loops.is_empty()
        || !matches!(
            f.surface.as_ref().map(|s| &s.geometry),
            Some(SurfaceGeometry::Plane { .. })
        )
    {
        return None;
    }
    let chain_edges: HashSet<EdgeHandle> = chain_n
        .iter()
        .chain(chain_p.iter())
        .filter_map(|he| Some(brep.half_edges.get(*he)?.edge))
        .collect();
    let mut flap_edges = HashSet::new();
    let start = brep.loops.get(f.outer_loop)?.half_edge;
    let mut he = start;
    for _ in 0..brep.half_edges.len() {
        let h = brep.half_edges.get(he)?;
        flap_edges.insert(h.edge);
        he = h.next;
        if he == start {
            break;
        }
    }
    if flap_edges != chain_edges || flap_edges.len() != chain_n.len() + chain_p.len() {
        return None;
    }
    Some(PairCornerCandidate {
        strip_n,
        strip_p,
        seam_edge,
        he_n,
        he_p,
        corner,
        anchor_near,
        anchor_far,
        chain_n,
        chain_p,
        flap,
    })
}

fn repair_pair_corner(
    brep: &mut BRep,
    c: &PairCornerCandidate,
    tolerance: &Tolerance,
) -> Result<(), BlendSeamCornerDeclineReason> {
    use BlendSeamCornerDeclineReason as Reason;
    let cyl_n = strip_cylinder(brep, c.strip_n).ok_or(Reason::EulerRejected)?;
    let cyl_p = strip_cylinder(brep, c.strip_p).ok_or(Reason::EulerRejected)?;
    let pos_of = |v: VertexHandle| -> Option<Point3> { Some(brep.vertices.get(v)?.position) };
    let corner_pos = pos_of(c.corner).ok_or(Reason::EulerRejected)?;
    let a_pos = pos_of(c.anchor_near).ok_or(Reason::EulerRejected)?;
    let b_pos = pos_of(c.anchor_far).ok_or(Reason::EulerRejected)?;

    // Author the seam geometry FIRST — every analytic decline exits
    // before any topology is touched.
    let frame = geometry::solve_equal_radius_seam_frame(&cyl_n, &cyl_p, a_pos, b_pos, tolerance)?;
    let primary = brep
        .edges
        .get(c.seam_edge)
        .ok_or(Reason::EulerRejected)?
        .half_edge;
    // After surgery `he_n` runs anchor_near → anchor_far and `he_p`
    // the reverse; the stored curve follows the primary half-edge.
    let (p_start, p_end) = if primary == c.he_n {
        (a_pos, b_pos)
    } else {
        (b_pos, a_pos)
    };
    let seam_curve = geometry::author_seam_curve(&frame, p_start, p_end, corner_pos, tolerance)?;

    // Euler surgery. kef keeps the given half-edge's face, so the
    // flap (the twin side) dies into strip_n; the remaining chain_n
    // stubs become spurs dissolved tip-first; chain_p contracts onto
    // the far anchor, whose position kve preserves (origin side).
    let last_n = *c.chain_n.last().ok_or(Reason::EulerRejected)?;
    brep.kef(last_n).map_err(|_| Reason::EulerRejected)?;
    for he in c.chain_n[..c.chain_n.len() - 1].iter().rev() {
        brep.kev(*he).map_err(|_| Reason::EulerRejected)?;
    }
    for he in &c.chain_p {
        brep.kve(*he).map_err(|_| Reason::EulerRejected)?;
    }
    // The seam edge must now run between the two anchors.
    let he_n_data = brep.half_edges.get(c.he_n).ok_or(Reason::EulerRejected)?;
    let he_p_data = brep
        .half_edges
        .get(he_n_data.twin)
        .ok_or(Reason::EulerRejected)?;
    if he_n_data.origin != c.anchor_near || he_p_data.origin != c.anchor_far {
        return Err(Reason::EulerRejected);
    }

    // Strip-side pcurves, anchored to the loop neighbours' charts.
    let staged: Vec<(
        FaceHandle,
        HalfEdgeHandle,
        popkern_core::geometry::pcurve::Pcurve,
    )> = [
        (c.strip_n, c.he_n, a_pos, b_pos),
        (c.strip_p, c.he_p, b_pos, a_pos),
    ]
    .into_iter()
    .map(|(face, he, origin_pos, dest_pos)| {
        let surface = brep
            .faces
            .get(face)
            .and_then(|f| f.surface.as_ref())
            .map(|s| s.geometry.clone())
            .ok_or(Reason::EulerRejected)?;
        let h = brep.half_edges.get(he).ok_or(Reason::EulerRejected)?;
        let origin_anchor = geometry::ChartAnchor {
            uv: neighbour_chart_point(brep, face, h.prev, true),
        };
        let dest_anchor = geometry::ChartAnchor {
            uv: neighbour_chart_point(brep, face, h.next, false),
        };
        let pcurve = geometry::author_side_pcurve(
            &surface,
            &frame,
            origin_pos,
            dest_pos,
            &origin_anchor,
            &dest_anchor,
            tolerance,
        )?;
        Ok((face, he, pcurve))
    })
    .collect::<Result<_, Reason>>()?;

    // Commit.
    if let Some(e) = brep.edges.get_mut(c.seam_edge) {
        e.curve = Some(seam_curve);
    }
    for (face, he, pcurve) in staged {
        if let Some(h) = brep.half_edges.get_mut(he) {
            h.pcurve = Some(pcurve.clone());
        }
        brep.set_face_use_pcurve(face, he, pcurve);
    }

    // Postconditions: the seam must sit on both strips and the body
    // must still validate.
    let pos_tol = (tolerance.linear() * 10.0).max(1e-9);
    if brep.measure_edge_deviation(c.seam_edge) > pos_tol {
        return Err(Reason::ValidationFailed);
    }
    brep.validate().map_err(|_| Reason::ValidationFailed)?;
    Ok(())
}

/// Authoritative chart position of a seam endpoint on one face: the
/// adjacent half-edge's pcurve endpoint (its end when the neighbour
/// arrives at the vertex, its start when it leaves).
fn neighbour_chart_point(
    brep: &BRep,
    face: FaceHandle,
    neighbour: HalfEdgeHandle,
    at_end: bool,
) -> Option<(f64, f64)> {
    let pc = brep.face_use_pcurve(face, neighbour).or(brep
        .half_edges
        .get(neighbour)?
        .pcurve
        .as_ref())?;
    let (t0, t1) = pc.parameter_range();
    let p = pc.point_at(if at_end { t1 } else { t0 });
    Some((p.x, p.y))
}
