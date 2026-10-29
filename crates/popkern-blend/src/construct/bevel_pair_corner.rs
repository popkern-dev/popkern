//! Bevel-pair corners for planar chamfers (PK-994).
//!
//! **The class.** A valence-3 vertex `V` where two edges are chamfered
//! with planar bevels `b1`, `b2` and the third edge stays sharp. The
//! bevels share the support `A`; `b1` also cuts `S1`, `b2` also cuts
//! `S2`, and the sharp edge is `e = S1 ∩ S2`. Every support is planar.
//! The exact result is the half-space configuration of the five
//! planes: the bevels mitre along `b1 ∩ b2` from `P = A ∩ b1 ∩ b2`,
//! and the sharp edge is cut back to whichever bevel reaches deeper
//! along it, `c1 = e ∩ b1` or `c2 = e ∩ b2`.
//!
//! - **Concurrent** (`c1 = c2`, e.g. every orthogonal corner, or equal
//!   face angles at `V` in `S1` and `S2`): one 4-valent vertex `Q`, the
//!   mitre edge `P–Q`. The anchor-split stage already lands this
//!   exactly (the two rails share `P`, both rails meet `e` at `Q`).
//! - **Non-concurrent** (a hexagonal prism's top edge with its
//!   vertical, slanted ramp sides): `c1 ≠ c2`. When `b1` reaches deeper
//!   the corner is `c1 = S1 ∩ S2 ∩ b1`, `R = b1 ∩ b2 ∩ S2`, edges
//!   `c1–R` (`b1 | S2`) and `R–P` (`b1 | b2`); symmetrically with `b2`
//!   deeper (`c2`, `R = b1 ∩ b2 ∩ S1`, connector `S1 | b2`). The
//!   anchor stage leaves the connector's inner end at the SHALLOWER
//!   rail crossing — off `b1` by the depth difference (measured: a
//!   hexagon top edge with its vertical, 9.5e-3 at d = 0.1).
//!
//! **Settlement, by construction.** From `P` (the shared junction
//! vertex) the stage walks the `b1 | b2` edge to its other end `X` and
//! checks the support `X` touches: the anchor stage leaves `X` at the
//! shallower bevel's rail crossing, on that bevel's own support, which
//! is the support the exact configuration puts `R` on. `X` moves to
//! `R` and its three edges are re-authored as exact segments. A
//! residue with `X` on the other support (the connector between the
//! wrong face pair) is not produced by the anchor stage; it declines
//! with a reason naming the corner instead of being guessed at.
//!
//! **Invariants:** at most one `kve` (the stale corner, when an earlier
//! settle left it on the mitre chain); otherwise a vertex moves and its
//! three edges are re-authored. Bevels and supports keep their handles. **Fallback:** a cluster that does not read as above
//! (a face missing, a vertex of the wrong valence) is left to the
//! serve gate, which declines with the structured check-9 reason; an
//! operator failure restores the corner's snapshot and declines the
//! operation with the reason. **Tolerance:** the caller's blend
//! tolerance; positions gate at `linear × 10`.

use super::EdgeGeometry;
use super::mitre_apex_corner::{Plane, destination, outgoing};
use super::mitre_apex_scaffold::three_plane_point;
use crate::error::BlendFailure;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{CurveOnEdge, Point3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};
use std::collections::BTreeMap;

/// Outcome of the stage.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct BevelPairReport {
    /// Non-concurrent corners rebuilt exactly.
    pub corners_settled: usize,
}

/// One candidate corner of the input, its planes resolved.
struct PairCorner {
    vertex: VertexHandle,
    /// One of the two chamfered edges (the failure's subject).
    edge: EdgeHandle,
    position: Point3,
    a: Plane,
    s1: Plane,
    s2: Plane,
    b1: Plane,
    b2: Plane,
}

fn plane_of(surface: &SurfaceGeometry) -> Option<Plane> {
    match surface {
        SurfaceGeometry::Plane { origin, normal, .. } => Some((*origin, normal.normalized()?)),
        _ => None,
    }
}

fn face_plane(brep: &BRep, face: FaceHandle) -> Option<Plane> {
    plane_of(&brep.faces.get(face)?.surface.as_ref()?.geometry)
}

/// Settle every non-concurrent bevel-pair corner of `input` in `out`.
pub(crate) fn settle_bevel_pair_corners(
    input: &BRep,
    out: &mut BRep,
    edge_geometries: &[EdgeGeometry],
    tolerance: &Tolerance,
) -> Result<BevelPairReport, (EdgeHandle, BlendFailure)> {
    let trace = popkern_core::trace::trace_enabled("bevel-pair");
    let gate = (tolerance.linear() * 10.0).max(1e-9);
    let mut report = BevelPairReport::default();
    for corner in candidate_corners(input, edge_geometries) {
        let c1 = three_plane_point(&corner.s1, &corner.s2, &corner.b1);
        let c2 = three_plane_point(&corner.s1, &corner.s2, &corner.b2);
        let (Some(c1), Some(c2)) = (c1, c2) else {
            continue;
        };
        if c1.distance_to(&c2) <= gate {
            continue; // concurrent: the anchor stage lands it exactly
        }
        let Some(p) = three_plane_point(&corner.a, &corner.b1, &corner.b2) else {
            continue;
        };
        // The deeper crossing is further from V along the sharp edge.
        let b1_deeper = c1.distance_to(&corner.position) >= c2.distance_to(&corner.position);
        // The support the connector's inner end R sits on: the other
        // bevel's own support.
        let r_support = if b1_deeper { corner.s2 } else { corner.s1 };
        let Some(r) = three_plane_point(&corner.b1, &corner.b2, &r_support) else {
            continue;
        };
        let snapshot = out.clone();
        match settle_one(out, &corner, p, r, &r_support, gate) {
            Ok(Ok(())) => {
                report.corners_settled += 1;
                if trace {
                    eprintln!(
                        "[bevel-pair] corner {:?}: R ({:.6},{:.6},{:.6}), b1 deeper {b1_deeper}",
                        corner.vertex, r.x, r.y, r.z
                    );
                }
            }
            Ok(Err(why)) => {
                *out = snapshot;
                if trace {
                    eprintln!(
                        "[bevel-pair] corner {:?}: cluster not recognised: {why}",
                        corner.vertex
                    );
                }
            }
            Err(reason) => {
                *out = snapshot;
                return Err((
                    corner.edge,
                    BlendFailure::InvalidTopology {
                        reason: format!("bevel pair corner at {:?}: {reason}", corner.vertex),
                        location: Some(corner.position),
                    },
                ));
            }
        }
    }
    Ok(report)
}

/// Valence-3 vertices of `input` with exactly two incident edges
/// chamfered by planar bevels sharing a support, all supports planar.
fn candidate_corners(input: &BRep, edge_geometries: &[EdgeGeometry]) -> Vec<PairCorner> {
    let by_edge: BTreeMap<EdgeHandle, usize> = edge_geometries
        .iter()
        .enumerate()
        .map(|(i, g)| (g.analysis.edge, i))
        .collect();
    let mut out = Vec::new();
    for (vh, vd) in input.vertices.iter() {
        let mut edges: Vec<EdgeHandle> = outgoing(input, vh)
            .into_iter()
            .filter_map(|he| input.half_edges.get(he).map(|h| h.edge))
            .collect();
        edges.sort();
        edges.dedup();
        if edges.len() != 3 {
            continue;
        }
        let blended: Vec<usize> = edges
            .iter()
            .filter_map(|e| by_edge.get(e).copied())
            .collect();
        if blended.len() != 2 {
            continue;
        }
        let (g1, g2) = (&edge_geometries[blended[0]], &edge_geometries[blended[1]]);
        let (Some(b1), Some(b2)) = (plane_of(&g1.surface), plane_of(&g2.surface)) else {
            continue;
        };
        let f1 = [g1.analysis.face_1, g1.analysis.face_2];
        let f2 = [g2.analysis.face_1, g2.analysis.face_2];
        let Some(&a) = f1.iter().find(|f| f2.contains(f)) else {
            continue;
        };
        let s1 = if f1[0] == a { f1[1] } else { f1[0] };
        let s2 = if f2[0] == a { f2[1] } else { f2[0] };
        if s1 == s2 {
            continue;
        }
        let (Some(pa), Some(p1), Some(p2)) = (
            face_plane(input, a),
            face_plane(input, s1),
            face_plane(input, s2),
        ) else {
            continue;
        };
        out.push(PairCorner {
            vertex: vh,
            edge: g1.analysis.edge,
            position: vd.position,
            a: pa,
            s1: p1,
            s2: p2,
            b1,
            b2,
        });
    }
    out
}

/// The face of `out` incident to `vertex` lying on `plane`.
fn incident_face_on(
    out: &BRep,
    vertex: VertexHandle,
    plane: &Plane,
    gate: f64,
) -> Option<FaceHandle> {
    outgoing(out, vertex).into_iter().find_map(|he| {
        let face = out.loops.get(out.half_edges.get(he)?.loop_)?.face;
        let fp = face_plane(out, face)?;
        let parallel = fp.1.dot(&plane.1).abs() >= 1.0 - 1e-9;
        (parallel && (fp.0 - plane.0).dot(&plane.1).abs() <= gate).then_some(face)
    })
}

fn faces_of(out: &BRep, he: HalfEdgeHandle) -> Option<(FaceHandle, FaceHandle)> {
    let h = out.half_edges.get(he)?;
    let f1 = out.loops.get(h.loop_)?.face;
    let f2 = out.loops.get(out.half_edges.get(h.twin)?.loop_)?.face;
    Some((f1, f2))
}

fn incident_faces(out: &BRep, vertex: VertexHandle) -> Vec<FaceHandle> {
    let mut faces: Vec<FaceHandle> = outgoing(out, vertex)
        .into_iter()
        .filter_map(|he| out.loops.get(out.half_edges.get(he)?.loop_).map(|l| l.face))
        .collect();
    faces.sort();
    faces.dedup();
    faces
}

/// Rebuild one corner: `Ok(Err(why))` when the cluster is not the
/// expected one (left to the serve gate; `why` goes to the trace).
fn settle_one(
    out: &mut BRep,
    corner: &PairCorner,
    p: Point3,
    r: Point3,
    r_support: &Plane,
    gate: f64,
) -> Result<Result<(), &'static str>, String> {
    // P: the shared junction vertex of the two rails on A.
    let Some(pv) = out
        .vertices
        .iter()
        .find(|(_, v)| v.position.distance_to(&p) <= gate)
        .map(|(h, _)| h)
    else {
        return Ok(Err("no vertex at the mitre point P"));
    };
    let (Some(fb1), Some(fb2)) = (
        incident_face_on(out, pv, &corner.b1, gate),
        incident_face_on(out, pv, &corner.b2, gate),
    ) else {
        return Ok(Err("P does not touch both bevel faces"));
    };
    // X: the other end of the b1 | b2 mitre edge.
    let Some(mitre) = outgoing(out, pv).into_iter().find(|&he| {
        faces_of(out, he).is_some_and(|(f, g)| (f, g) == (fb1, fb2) || (f, g) == (fb2, fb1))
    }) else {
        return Ok(Err("no b1 | b2 mitre edge at P"));
    };
    let mut x = destination(out, mitre).ok_or("mitre edge without destination")?;
    // The stale original corner may still sit on the mitre chain as a
    // valence-2 vertex between the bevels (P → V → X): the reflex-pair
    // settle that usually contracts it declines when it cannot bracket
    // its rail root (measured: a hexagon corner at 1/100 scale). One
    // `kve` merges it into X, so this stage does not depend on that.
    let x_out = outgoing(out, x);
    let between_bevels = |he: HalfEdgeHandle| {
        faces_of(out, he).is_some_and(|(f, g)| (f, g) == (fb1, fb2) || (f, g) == (fb2, fb1))
    };
    if x_out.len() == 2 && x_out.iter().all(|&he| between_bevels(he)) {
        let onward = x_out
            .iter()
            .copied()
            .find(|&he| destination(out, he) != Some(pv))
            .ok_or("stale corner has no onward edge")?;
        let far = destination(out, onward).ok_or("onward edge without destination")?;
        let far_to_corner = out.half_edges.get(onward).ok_or("stale onward edge")?.twin;
        out.kve(far_to_corner)
            .map_err(|e| format!("kve stale corner: {e}"))?;
        x = far;
    }
    let x_faces = incident_faces(out, x);
    if popkern_core::trace::trace_enabled("bevel-pair") {
        let q = out.vertices[x].position;
        let nbrs: Vec<String> = outgoing(out, x)
            .into_iter()
            .filter_map(|he| destination(out, he))
            .map(|v| {
                let n = out.vertices[v].position;
                format!("({:.7},{:.7},{:.7})", n.x, n.y, n.z)
            })
            .collect();
        eprintln!(
            "[bevel-pair] X ({:.7},{:.7},{:.7}) faces {} neighbours {}",
            q.x,
            q.y,
            q.z,
            x_faces.len(),
            nbrs.join(" ")
        );
    }
    if x_faces.len() != 3 || !x_faces.contains(&fb1) || !x_faces.contains(&fb2) {
        return Ok(Err(
            "the mitre edge's far end is not a bevel-bevel-support vertex",
        ));
    }
    let Some(&x_support) = x_faces.iter().find(|f| **f != fb1 && **f != fb2) else {
        return Ok(Err("the mitre edge's far end touches no support"));
    };
    // Support faces are matched by plane: the per-edge apply does not
    // keep their handles (measured: a triangular prism's side face
    // 3v1 in the input, 8v1 after apply).
    let on_r_support = face_plane(out, x_support).is_some_and(|fp| {
        fp.1.dot(&r_support.1).abs() >= 1.0 - 1e-9
            && (fp.0 - r_support.0).dot(&r_support.1).abs() <= gate
    });
    if !on_r_support {
        // The anchor stage lands the connector's inner end at the
        // SHALLOWER bevel's rail crossing, on that bevel's own support,
        // which is where R belongs. A residue with it on the other
        // support has the connector between the wrong face pair:
        // decline, naming it, rather than ship or guess.
        return Err(format!(
            "connector end {x:?} lies on the support the exact corner does not put it on \
             (residue topology differs from the non-concurrent mitre)"
        ));
    }
    if let Some(v) = out.vertices.get_mut(x) {
        v.position = r;
    }
    for he in outgoing(out, x) {
        author_segment(out, he, gate / 10.0)?;
    }
    Ok(Ok(()))
}

/// Author the edge of `he` as the straight segment between its
/// endpoints (both faces planar, so the edge is their line).
fn author_segment(out: &mut BRep, he: HalfEdgeHandle, tol: f64) -> Result<(), String> {
    let h = out.half_edges.get(he).ok_or("stale half-edge")?;
    let edge = h.edge;
    let primary = out.edges.get(edge).ok_or("stale edge")?.half_edge;
    let o = out.half_edges.get(primary).ok_or("stale primary")?.origin;
    let d = destination(out, primary).ok_or("edge without destination")?;
    let (po, pd) = (out.vertices[o].position, out.vertices[d].position);
    let segment = CurveOnEdge::line_segment(po, pd, tol).ok_or("degenerate segment")?;
    if let Some(e) = out.edges.get_mut(edge) {
        e.curve = Some(segment);
    }
    Ok(())
}
