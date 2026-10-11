//! Reflex pair corners (PK-740): two concave fillets on adjacent root
//! edges of a boss, meeting where the SUPPORT's corner is reflex and
//! the boss's convex corner edge stays sharp (a rib's base loop on a
//! plate, a fin's root loop on a housing).
//!
//! **Residue.** The anchor-split stage lands both strips' support-side
//! rail ends at the rails' mitre point `M` (one spur from the corner
//! `V`, shared), and both wall-side anchors at `T` on the boss's sharp
//! corner edge. After the consumed-merge stage the two strips share a
//! two-edge chain `M → V → T` through the original corner, which lies
//! `(√2−1)·r`-class off both strips — `V` has valence 2.
//!
//! **Repair, by construction.** The strips' true mutual boundary is
//! their intersection from `M` to `T`. For equal-radius cylinder strips
//! it is the exact bisector-plane ellipse (`blend_seam_corner_geometry`,
//! the PK-152 family). The surgery is one `kve` contracting `V → T`
//! into `T` (V−1, E−1: χ-neutral), leaving the spur edge as the single
//! shared seam `M–T`, which is re-authored with exact sinusoid pcurves
//! in both strip charts.
//!
//! **Invariants:** strips keep their handles and surfaces; support and
//! wall faces are untouched. **Fallback:** every decline restores the
//! corner's snapshot and is recorded, typed; the corner-resolution gate
//! then declines the fillet with that reason. **Tolerance:** the
//! caller's blend tolerance (`linear × 10` position gate, PK-152).

use std::collections::HashSet;

use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

use crate::types::{BlendSeamCornerDecline, BlendSeamCornerDeclineReason as Reason};

use popkern_core::geometry::{CurveOnEdge, Point3};

use super::blend_seam_corner::strip_cylinder;
use super::blend_seam_corner_geometry as geometry;
use super::blend_seam_corner_repair::{commit_edge_authoring, neighbour_chart_point};
use super::blend_seam_corner_swept::surface_distance;
use super::surface_pair_march as marched;

#[path = "blend_seam_corner_reflex_wall_end.rs"]
mod wall_end;

/// One recognized reflex pair corner.
struct ReflexChain {
    corner: VertexHandle,
    /// The shared spur edge `M–V`; survives as the seam `M–T`.
    seam: EdgeHandle,
    /// The shared wall stub `V–T`'s half-edge from `T` to `V`.
    stub_from_t: HalfEdgeHandle,
    strips: (FaceHandle, FaceHandle),
    /// The chain's far end when it sits on one strip's wall rail but
    /// off the other strip: an asymmetric mitre on a curved support
    /// (PK-742), moved to the rail's true end before the seam is built.
    wall_end: Option<wall_end::WallEnd>,
}

/// Settle every reflex pair corner among the strips the operation
/// created. Returns `(seams built, declines)`.
pub(super) fn settle_reflex_pair_corners(
    result: &mut BRep,
    created: &HashSet<FaceHandle>,
    tolerance: &Tolerance,
) -> (usize, Vec<BlendSeamCornerDecline>) {
    let mut built = 0;
    let mut declines = Vec::new();
    let trace = popkern_core::trace::trace_enabled("blend-seam");
    for chain in recognize(result, created, tolerance) {
        let snapshot = result.clone();
        let outcome = repair(result, &chain, tolerance);
        if trace {
            eprintln!(
                "[bsc-reflex] corner {:?} strips {:?} -> {outcome:?}",
                chain.corner, chain.strips
            );
        }
        match outcome {
            Ok(()) => built += 1,
            Err(reason) => {
                *result = snapshot;
                declines.push(BlendSeamCornerDecline {
                    strips: chain.strips,
                    corner: chain.corner,
                    reason,
                });
            }
        }
    }
    (built, declines)
}

fn faces_of(brep: &BRep, edge: EdgeHandle) -> Option<(FaceHandle, FaceHandle)> {
    let he = brep.half_edges.get(brep.edges.get(edge)?.half_edge)?;
    let twin = brep.half_edges.get(he.twin)?;
    Some((
        brep.loops.get(he.loop_)?.face,
        brep.loops.get(twin.loop_)?.face,
    ))
}

fn recognize(
    brep: &BRep,
    created: &HashSet<FaceHandle>,
    tolerance: &Tolerance,
) -> Vec<ReflexChain> {
    let pos_tol = (tolerance.linear() * 10.0).max(1e-9);
    let mut out = Vec::new();
    for (v, vd) in brep.vertices.iter() {
        let outgoing: Vec<HalfEdgeHandle> = brep
            .half_edges
            .iter()
            .filter(|(_, h)| h.origin == v)
            .map(|(h, _)| h)
            .collect();
        if outgoing.len() != 2 {
            continue;
        }
        let edges: Vec<EdgeHandle> = outgoing.iter().map(|&h| brep.half_edges[h].edge).collect();
        let (Some(fa), Some(fb)) = (faces_of(brep, edges[0]), faces_of(brep, edges[1])) else {
            continue;
        };
        let same_pair = (fa == fb || (fa.0 == fb.1 && fa.1 == fb.0)) && fa.0 != fa.1;
        if !same_pair || !created.contains(&fa.0) || !created.contains(&fa.1) {
            continue;
        }
        let geom = |f: FaceHandle| {
            brep.faces
                .get(f)
                .and_then(|d| d.surface.as_ref())
                .map(|s| s.geometry.clone())
        };
        let (Some(ga), Some(gb)) = (geom(fa.0), geom(fa.1)) else {
            continue;
        };
        let on_both = |p: Point3| {
            surface_distance(&ga, &p, None) <= pos_tol && surface_distance(&gb, &p, None) <= pos_tol
        };
        let far = |h: HalfEdgeHandle| {
            let twin = brep.half_edges[h].twin;
            brep.half_edges[twin].origin
        };
        let (Some(mut x0), Some(mut x1)) = (
            brep.vertices.get(far(outgoing[0])).map(|d| d.position),
            brep.vertices.get(far(outgoing[1])).map(|d| d.position),
        ) else {
            continue;
        };
        // One end may sit on a strip's wall rail but off the other
        // strip (PK-742): its true position is where the rail meets it.
        let mut asymmetric = None;
        for (i, x) in [&mut x0, &mut x1].into_iter().enumerate() {
            if !on_both(*x)
                && let Some(end) = wall_end::find(brep, far(outgoing[i]), fa, pos_tol)
            {
                *x = end.target;
                asymmetric = Some(end);
                break;
            }
        }
        if !on_both(x0) || !on_both(x1) || x0.distance_to(&x1) <= pos_tol {
            continue;
        }
        // The stale corner is off both strips — or a tail pass already
        // slid it onto one chain end (a degenerate stub; measured on a
        // fin root loop: the corner sat on the mitre point).
        let p = vd.position;
        if on_both(p) && p.distance_to(&x0) > pos_tol && p.distance_to(&x1) > pos_tol {
            continue;
        }
        // Either edge may survive; contract the second (T side).
        out.push(ReflexChain {
            corner: v,
            seam: edges[0],
            stub_from_t: brep.half_edges[outgoing[1]].twin,
            strips: fa,
            wall_end: asymmetric,
        });
    }
    out
}

fn repair(brep: &mut BRep, c: &ReflexChain, tolerance: &Tolerance) -> Result<(), Reason> {
    let pos_tol = (tolerance.linear() * 10.0).max(1e-9);
    if let Some(end) = &c.wall_end {
        brep.vertices
            .get_mut(end.vertex)
            .ok_or(Reason::EulerRejected)?
            .position = end.target;
    }
    let corner_pos = brep
        .vertices
        .get(c.corner)
        .map(|v| v.position)
        .ok_or(Reason::EulerRejected)?;
    // Contract V → T into T: the spur edge becomes the seam M–T.
    brep.kve(c.stub_from_t).map_err(|_| Reason::EulerRejected)?;
    let primary = brep
        .edges
        .get(c.seam)
        .ok_or(Reason::EulerRejected)?
        .half_edge;
    let twin = brep
        .half_edges
        .get(primary)
        .ok_or(Reason::EulerRejected)?
        .twin;
    let pos = |h: HalfEdgeHandle| -> Result<_, Reason> {
        let o = brep.half_edges.get(h).ok_or(Reason::EulerRejected)?.origin;
        Ok(brep.vertices.get(o).ok_or(Reason::EulerRejected)?.position)
    };
    let (p_start, p_end) = (pos(primary)?, pos(twin)?);
    let sides: Vec<(HalfEdgeHandle, FaceHandle)> = [primary, twin]
        .into_iter()
        .map(|he| {
            let h = brep.half_edges.get(he).ok_or(Reason::EulerRejected)?;
            Ok((
                he,
                brep.loops.get(h.loop_).ok_or(Reason::EulerRejected)?.face,
            ))
        })
        .collect::<Result<_, Reason>>()?;
    let equal_cylinders = match (
        strip_cylinder(brep, c.strips.0),
        strip_cylinder(brep, c.strips.1),
    ) {
        (Some(a), Some(b)) => (a.radius - b.radius).abs() <= pos_tol.max(a.radius * 1e-9),
        _ => false,
    };
    let (seam_curve, staged) = if equal_cylinders {
        analytic_seam(brep, c, &sides, p_start, p_end, corner_pos, tolerance)?
    } else {
        marched_seam(brep, &sides, p_start, p_end, pos_tol)?
    };
    commit_edge_authoring(brep, c.seam, seam_curve, staged);
    if brep.measure_edge_deviation(c.seam) > pos_tol {
        return Err(Reason::ValidationFailed);
    }
    if let Some(end) = &c.wall_end {
        wall_end::reauthor_incident(brep, end.vertex, c.seam, pos_tol)?;
    }
    // Both seam ends: a support rail ending at the mitre point may span
    // past the root edge's trim (PK-742).
    for he in [primary, twin] {
        let v = brep.half_edges.get(he).ok_or(Reason::EulerRejected)?.origin;
        wall_end::reauthor_contact_rails(brep, v, c.strips, pos_tol)?;
    }
    brep.validate().map_err(|_| Reason::ValidationFailed)?;
    Ok(())
}

type Staged = Vec<(FaceHandle, HalfEdgeHandle, Pcurve)>;

/// Equal-radius cylinder strips: the exact bisector-plane ellipse.
fn analytic_seam(
    brep: &BRep,
    c: &ReflexChain,
    sides: &[(HalfEdgeHandle, FaceHandle)],
    p_start: Point3,
    p_end: Point3,
    corner_pos: Point3,
    tolerance: &Tolerance,
) -> Result<(CurveOnEdge, Staged), Reason> {
    let cyl_a = strip_cylinder(brep, c.strips.0).ok_or(Reason::EulerRejected)?;
    let cyl_b = strip_cylinder(brep, c.strips.1).ok_or(Reason::EulerRejected)?;
    let frame = geometry::solve_equal_radius_seam_frame(&cyl_a, &cyl_b, p_start, p_end, tolerance)?;
    let seam_curve = geometry::author_seam_curve(&frame, p_start, p_end, corner_pos, tolerance)?;
    let mut staged = Vec::new();
    for (i, &(he, face)) in sides.iter().enumerate() {
        let h = brep.half_edges.get(he).ok_or(Reason::EulerRejected)?;
        let surface = face_geom(brep, face)?;
        let (origin_pos, dest_pos) = if i == 0 {
            (p_start, p_end)
        } else {
            (p_end, p_start)
        };
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
        staged.push((face, he, pcurve));
    }
    Ok((seam_curve, staged))
}

/// Any other strip pair (a cylinder root meeting a torus root on a
/// curved housing): the marched strip×strip intersection.
fn marched_seam(
    brep: &BRep,
    sides: &[(HalfEdgeHandle, FaceHandle)],
    p_start: Point3,
    p_end: Point3,
    pos_tol: f64,
) -> Result<(CurveOnEdge, Staged), Reason> {
    let ga = face_geom(brep, sides[0].1)?;
    let gb = face_geom(brep, sides[1].1)?;
    let seam = marched::march_intersection(&ga, &gb, p_start, p_end, pos_tol, pos_tol)?;
    let curve = marched::intersection_curve(&seam, false)?;
    let mut staged = Vec::new();
    for (i, &(he, face)) in sides.iter().enumerate() {
        let h = brep.half_edges.get(he).ok_or(Reason::EulerRejected)?;
        let g = if i == 0 { &ga } else { &gb };
        // The march's `from` is p_start: the primary side's origin,
        // the twin side's destination.
        let seed = if i == 0 {
            neighbour_chart_point(brep, face, h.prev, true)
        } else {
            neighbour_chart_point(brep, face, h.next, false)
        };
        let pcurve = marched::intersection_pcurve(g, &seam, i == 1, seed)?;
        staged.push((face, he, pcurve));
    }
    Ok((curve, staged))
}

fn face_geom(
    brep: &BRep,
    face: FaceHandle,
) -> Result<popkern_core::geometry::surface::SurfaceGeometry, Reason> {
    brep.faces
        .get(face)
        .and_then(|f| f.surface.as_ref())
        .map(|s| s.geometry.clone())
        .ok_or(Reason::EulerRejected)
}
