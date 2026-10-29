//! Chamfers at mixed corners: fills first (PK-996).
//!
//! **The class.** A vertex where a CONCAVE chamfer (a fill: it adds
//! material) and a CONVEX chamfer (a cut) both end, such as an
//! L-bracket's re-entrant vertical chamfered with a top edge running
//! into it. There is no order-independent corner: bounding the cut by
//! the fill's plane, bounding the fill by the cut's plane, and bounding
//! each by the other on its own side of the shared wall are three valid
//! solids (measured on the L at d = 0.1: 1.4981667, 1.4968333 and
//! 1.4975000, the last leaving a d × d square of the shared wall
//! between the bevels).
//!
//! **Semantics: fills first.** A request with a mixed corner runs as
//! one transactional operation in two stages: every concave chamfer of
//! the request, then every other chamfer on that result. The cut ends
//! against the fill, whose face reaches the top: what chamfering the
//! concave edges and then the convex edges in separate operations
//! gives, now in one feature. Requests without a mixed corner are
//! unchanged (one pass).
//!
//! **Remapping.** The fill shortens the convex edges that end at it.
//! Each convex spec is carried to the intermediate body by geometry:
//! the edges whose ends and midpoint lie on the original edge's curve,
//! within its span (two edges of a valid solid never overlap, so these
//! are exactly its surviving pieces). Two-distance and distance-angle
//! reference faces are carried by surface: the face at the remapped
//! edge lying on the original face's surface.
//!
//! **Naming.** Both stages run untracked and the result is named once
//! against the ORIGINAL input, so the document sees one chamfer
//! feature with one Generated/Modified role sequence.
//!
//! **Fallback.** Either stage declining hands the input back with the
//! failures on the requested edges, restated at the mixed corner as
//! [`BlendFailure::VertexBlendFailed`] naming it; a convex edge the
//! fill swallowed whole declines the same way.
//! **Tolerance:** the caller's, curve containment at `linear × 10`.

use crate::analysis::{Convexity, analyze_edges};
use crate::cancellation::CancellationToken;
use crate::error::{BlendError, BlendFailure};
use crate::types::{BlendResult, ChamferKind, ChamferSpec};
use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::history::TrackingContext;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};
use std::collections::BTreeMap;

/// The reason prefix the remediation mapping recognises.
pub(crate) const MIXED_CORNER: &str = "mixed chamfer corner";

/// Vertices where requested chamfers of both convexities end, with the
/// requested edges at each.
pub(crate) fn mixed_corners(
    brep: &BRep,
    specs: &[ChamferSpec],
) -> Vec<(VertexHandle, Vec<EdgeHandle>)> {
    let edges: Vec<EdgeHandle> = specs.iter().map(|s| s.edge).collect();
    let Ok(analyses) = analyze_edges(brep, &edges, false) else {
        return Vec::new();
    };
    let mut at: BTreeMap<VertexHandle, (Vec<EdgeHandle>, bool, bool)> = BTreeMap::new();
    for a in &analyses {
        for v in [a.vertex_start, a.vertex_end] {
            let entry = at.entry(v).or_default();
            entry.0.push(a.edge);
            entry.1 |= a.convexity == Convexity::Convex;
            entry.2 |= a.convexity == Convexity::Concave;
        }
    }
    at.into_iter()
        .filter(|(_, (_, convex, concave))| *convex && *concave)
        .map(|(v, (edges, _, _))| (v, edges))
        .collect()
}

/// Run `specs` fills first (see the module docs). `pipeline` is the
/// single-pass chamfer (untracked; this names the result).
pub(crate) fn chamfer_fills_first(
    brep: &BRep,
    specs: &[ChamferSpec],
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
    cancel: Option<&CancellationToken>,
    mixed: &[(VertexHandle, Vec<EdgeHandle>)],
    pipeline: impl Fn(
        &BRep,
        &[ChamferSpec],
        Option<&CancellationToken>,
    ) -> Result<BlendResult, BlendError>,
) -> Result<BlendResult, BlendError> {
    let total = specs.len();
    let edges: Vec<EdgeHandle> = specs.iter().map(|s| s.edge).collect();
    let concave: Vec<EdgeHandle> = analyze_edges(brep, &edges, false)
        .map_err(BlendError::NoValidEdges)?
        .into_iter()
        .filter(|a| a.convexity == Convexity::Concave)
        .map(|a| a.edge)
        .collect();
    let (fills, cuts): (Vec<ChamferSpec>, Vec<ChamferSpec>) = specs
        .iter()
        .cloned()
        .partition(|s| concave.contains(&s.edge));
    let decline = |failures: Vec<(EdgeHandle, BlendFailure)>| {
        let failures = failures
            .into_iter()
            .map(
                |(edge, failure)| match mixed.iter().find(|(_, es)| es.contains(&edge)) {
                    Some((vertex, _)) => (edge, restate(brep, *vertex, &failure)),
                    None => (edge, failure),
                },
            )
            .collect();
        BlendError::partial_failure(total, brep.clone(), failures)
    };
    let failures_of = |error: BlendError| match error {
        BlendError::PartialFailure { failures, .. } => Ok(failures),
        other => Err(other),
    };

    // Stage 1: the fills.
    let filled = match pipeline(brep, &fills, cancel) {
        Ok(r) => r,
        Err(e) => return Err(decline(failures_of(e)?)),
    };

    // Stage 2: the cuts, carried onto the filled body.
    let gate = (tolerance.linear() * 10.0).max(1e-12);
    let mut carried: Vec<ChamferSpec> = Vec::new();
    let mut origin: BTreeMap<EdgeHandle, EdgeHandle> = BTreeMap::new();
    for spec in &cuts {
        let pieces = remap_edge(brep, spec.edge, &filled.brep, gate);
        if pieces.is_empty() {
            return Err(decline(vec![(
                spec.edge,
                BlendFailure::InconsistentGeometry {
                    reason: "the concave chamfer consumed this edge".into(),
                    location: edge_midpoint(brep, spec.edge),
                },
            )]));
        }
        for piece in pieces {
            let kind =
                remap_kind(brep, &spec.kind, &filled.brep, piece, gate).ok_or_else(|| {
                    decline(vec![(
                        spec.edge,
                        BlendFailure::InconsistentGeometry {
                            reason: "the chamfer's reference face did not survive the concave \
                                     chamfer"
                                .into(),
                            location: edge_midpoint(brep, spec.edge),
                        },
                    )])
                })?;
            origin.insert(piece, spec.edge);
            carried.push(ChamferSpec { edge: piece, kind });
        }
    }
    let cut = match pipeline(&filled.brep, &carried, cancel) {
        Ok(r) => r,
        Err(e) => {
            let failures = failures_of(e)?
                .into_iter()
                .map(|(edge, f)| (origin.get(&edge).copied().unwrap_or(edge), f))
                .collect();
            return Err(decline(failures));
        }
    };

    let mut result = cut;
    if let Some(ctx) = tracking {
        let pre = crate::naming::PreOperationSnapshot::capture(brep);
        crate::naming::assign_fillet_names(&mut result.brep, ctx, &pre.faces);
        crate::naming::assign_output_edge_vertex_names(&mut result.brep, ctx, &pre, &tolerance);
    }
    for (edge, _) in &mut result.edge_warnings {
        if let Some(o) = origin.get(edge) {
            *edge = *o;
        }
    }
    let mut warnings = filled.edge_warnings;
    warnings.append(&mut result.edge_warnings);
    result.edge_warnings = warnings;
    let mut patches = filled.corner_patches;
    patches.append(&mut result.corner_patches);
    result.corner_patches = patches;
    result.repair_report.mitre_apex_corners += filled.repair_report.mitre_apex_corners;
    result.repair_report.bevel_pair_corners += filled.repair_report.bevel_pair_corners;
    result.repair_report.truncated_apex_corners += filled.repair_report.truncated_apex_corners;
    result.repair_report.bevel_growth_flips += filled.repair_report.bevel_growth_flips;
    result.repair_report.bevel_consumed_faces += filled.repair_report.bevel_consumed_faces;
    Ok(result)
}

/// Restate a failure at a mixed corner, keeping the original reason.
fn restate(brep: &BRep, vertex: VertexHandle, failure: &BlendFailure) -> BlendFailure {
    BlendFailure::VertexBlendFailed {
        vertex,
        reason: format!(
            "{MIXED_CORNER} (fills-first construction declined: {})",
            failure_reason(failure)
        ),
        location: brep.vertices.get(vertex).map(|v| v.position),
    }
}

fn failure_reason(failure: &BlendFailure) -> String {
    match failure {
        BlendFailure::InvalidTopology { reason, .. }
        | BlendFailure::InconsistentGeometry { reason, .. }
        | BlendFailure::VertexBlendFailed { reason, .. }
        | BlendFailure::IntersectionFailed { reason, .. }
        | BlendFailure::OffsetFailed { reason, .. } => reason.clone(),
        other => format!("{other:?}"),
    }
}

fn endpoints(brep: &BRep, edge: EdgeHandle) -> Option<(Point3, Point3)> {
    let he = brep.half_edges.get(brep.edges.get(edge)?.half_edge)?;
    let twin = brep.half_edges.get(he.twin)?;
    Some((
        brep.vertices.get(he.origin)?.position,
        brep.vertices.get(twin.origin)?.position,
    ))
}

fn edge_midpoint(brep: &BRep, edge: EdgeHandle) -> Option<Point3> {
    let e = brep.edges.get(edge)?;
    match e.curve.as_ref() {
        Some(c) => Some(c.geometry.point_at(0.5 * (c.t_start + c.t_end))),
        None => endpoints(brep, edge).map(|(a, b)| a + (b - a) * 0.5),
    }
}

/// Distance from `p` to the original edge (its curve between its
/// ends, or its chord): coarse samples, then ternary refinement.
fn distance_to_edge(brep: &BRep, edge: EdgeHandle, p: &Point3) -> Option<f64> {
    let e = brep.edges.get(edge)?;
    let (a, b) = endpoints(brep, edge)?;
    let at = |s: f64| match e.curve.as_ref() {
        Some(c) => c.geometry.point_at(c.t_start + (c.t_end - c.t_start) * s),
        None => a + (b - a) * s,
    };
    let (mut best_s, mut best) = (0.0, f64::INFINITY);
    for k in 0..=128 {
        let s = f64::from(k) / 128.0;
        let d = at(s).distance_to(p);
        if d < best {
            (best_s, best) = (s, d);
        }
    }
    let (mut lo, mut hi) = (
        (best_s - 1.0 / 128.0).max(0.0),
        (best_s + 1.0 / 128.0).min(1.0),
    );
    for _ in 0..60 {
        let (m1, m2) = (lo + (hi - lo) / 3.0, hi - (hi - lo) / 3.0);
        if at(m1).distance_to(p) < at(m2).distance_to(p) {
            hi = m2;
        } else {
            lo = m1;
        }
    }
    Some(at(0.5 * (lo + hi)).distance_to(p))
}

/// The edges of `staged` that are surviving pieces of `edge` of
/// `original`.
fn remap_edge(original: &BRep, edge: EdgeHandle, staged: &BRep, gate: f64) -> Vec<EdgeHandle> {
    staged
        .edges
        .keys()
        .filter(|&e| {
            let (Some((p, q)), Some(m)) = (endpoints(staged, e), edge_midpoint(staged, e)) else {
                return false;
            };
            p.distance_to(&q) > gate
                && [p, q, m]
                    .iter()
                    .all(|x| distance_to_edge(original, edge, x).is_some_and(|d| d <= gate))
        })
        .collect()
}

fn face_surface(brep: &BRep, face: FaceHandle) -> Option<&SurfaceGeometry> {
    Some(&brep.faces.get(face)?.surface.as_ref()?.geometry)
}

fn surface_distance(surface: &SurfaceGeometry, p: &Point3) -> Option<f64> {
    let (u, v) = surface.closest_point(p)?;
    Some(surface.point_at(u, v).distance_to(p))
}

/// The two faces of `edge`.
fn edge_faces(brep: &BRep, edge: EdgeHandle) -> Option<[FaceHandle; 2]> {
    let he = brep.half_edges.get(brep.edges.get(edge)?.half_edge)?;
    let twin = brep.half_edges.get(he.twin)?;
    Some([
        brep.loops.get(he.loop_)?.face,
        brep.loops.get(twin.loop_)?.face,
    ])
}

/// Carry a chamfer kind's reference face to the face at `piece` lying
/// on the same surface.
fn remap_kind(
    original: &BRep,
    kind: &ChamferKind,
    staged: &BRep,
    piece: EdgeHandle,
    gate: f64,
) -> Option<ChamferKind> {
    let carry = |face: FaceHandle| -> Option<FaceHandle> {
        let reference = face_surface(original, face)?;
        edge_faces(staged, piece)?.into_iter().find(|&f| {
            staged
                .face_outer_loop_positions(f)
                .iter()
                .all(|x| surface_distance(reference, x).is_some_and(|d| d <= gate))
        })
    };
    Some(match kind {
        ChamferKind::Equal(d) => ChamferKind::Equal(*d),
        ChamferKind::TwoDistance {
            face_1,
            dist_1,
            dist_2,
        } => ChamferKind::TwoDistance {
            face_1: carry(*face_1)?,
            dist_1: *dist_1,
            dist_2: *dist_2,
        },
        ChamferKind::DistanceAngle {
            face,
            distance,
            angle,
        } => ChamferKind::DistanceAngle {
            face: carry(*face)?,
            distance: *distance,
            angle: *angle,
        },
    })
}
