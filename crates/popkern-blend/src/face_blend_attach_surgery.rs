//! Trim-surgery stage of the TrimBoth face-face blend attach
//! (PK-142 stage 2): wall and end-face splits, doomed-region flood,
//! kef-merge + spur dissolution, and strip re-surfacing. See
//! `face_blend_attach` for the full transactional contract.

use std::collections::HashSet;

use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry, SurfaceOnFace};
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle};

use crate::contact::ContactSolution;
use crate::error::BlendError;
use crate::face_blend::FaceBlendSpec;
use crate::face_blend_attach::AttachOutcome;
use crate::face_blend_attach_ends::{EndOutcome, split_end_face};
use crate::face_blend_attach_geometry::{
    arc_departure, decline, ensure_boundary_vertex, face_surface, loop_half_edges, pick_doomed,
    project_uv, uniform_polyline_nurbs, uv_polyline_nurbs,
};
use crate::face_blend_sheet::{StripSkin, strip_skin};

struct WallSplit {
    kept: FaceHandle,
    doomed: FaceHandle,
    edge: EdgeHandle,
}

pub(crate) fn surgery(
    brep: &BRep,
    spec: &FaceBlendSpec,
    chain: &[ContactSolution],
    bulge: Option<Vec3>,
    tolerance: Tolerance,
) -> Result<AttachOutcome, BlendError> {
    let mut work = brep.clone();
    let skin = strip_skin(chain, spec.radius, bulge);

    let wall_a = split_wall(
        &mut work,
        spec.face_a,
        &skin.side_a,
        chain,
        true,
        bulge,
        &tolerance,
    )?;
    let wall_b = split_wall(
        &mut work,
        spec.face_b,
        &skin.side_b,
        chain,
        false,
        bulge,
        &tolerance,
    )?;
    let wall_faces = [wall_a.kept, wall_a.doomed, wall_b.kept, wall_b.doomed];

    let first = chain
        .first()
        .ok_or_else(|| decline("end split", "clipped chain is empty".into()))?;
    let last = chain
        .last()
        .ok_or_else(|| decline("end split", "clipped chain is empty".into()))?;
    let end_start = split_end_face(&mut work, &skin.start_arc, first, &wall_faces, &tolerance)?;
    let end_end = split_end_face(&mut work, &skin.end_arc, last, &wall_faces, &tolerance)?;

    let mut blocked = vec![wall_a.edge, wall_b.edge];
    let mut barriers: Vec<FaceHandle> = Vec::new();
    let mut members = vec![("wall B doomed part", wall_b.doomed)];
    let mut kept_faces = vec![wall_a.kept, wall_b.kept];
    let mut through_ends: Vec<&crate::face_blend_attach_through::ThroughEnd> = Vec::new();
    for (label, outcome) in [
        ("start-end doomed part", &end_start),
        ("far-end doomed part", &end_end),
    ] {
        match outcome {
            EndOutcome::Closed(split) => {
                blocked.push(split.edge);
                members.push((label, split.doomed));
                kept_faces.push(split.kept);
            }
            EndOutcome::Through(t) => {
                barriers.push(t.end_face);
                kept_faces.push(t.end_face);
                through_ends.push(t);
            }
        }
    }
    let region = flood_doomed(&work, wall_a.doomed, &blocked, &barriers)?;
    for (label, member) in members {
        if !region.contains(&member) {
            return Err(decline(
                "doomed-region flood",
                format!("{label} {member:?} is not connected to the gap region"),
            ));
        }
    }
    for kept in kept_faces {
        if region.contains(&kept) {
            return Err(decline(
                "doomed-region flood",
                format!("flood escaped past the trim curves into kept face {kept:?}"),
            ));
        }
    }
    if region.len() >= work.faces.len() {
        return Err(decline(
            "doomed-region flood",
            format!(
                "flood covered {} of {} faces — the trim curves do not enclose the gap",
                region.len(),
                work.faces.len()
            ),
        ));
    }

    let mut strip_face = merge_region(&mut work, region, &blocked)?;
    // PK-163: extend the end caps of through ends down to the arcs
    // before the strip resurface (the mef arcs complete the strip's
    // four-border loop).
    for t in &through_ends {
        strip_face = crate::face_blend_attach_through::extend_through_end(
            &mut work, strip_face, t, &tolerance,
        )?;
    }
    resurface_strip(&mut work, strip_face, &skin, chain, &tolerance)?;

    work.validate()
        .map_err(|e| decline("post-surgery validation", format!("{e:?}")))?;
    Ok(AttachOutcome {
        body: work,
        strip_face,
        chain: chain.to_vec(),
    })
}

#[allow(clippy::too_many_arguments)]
fn split_wall(
    work: &mut BRep,
    face: FaceHandle,
    side: &[Point3],
    chain: &[ContactSolution],
    is_side_a: bool,
    bulge: Option<Vec3>,
    tolerance: &Tolerance,
) -> Result<WallSplit, BlendError> {
    let first = *side
        .first()
        .ok_or_else(|| decline("wall split", "contact polyline is empty".into()))?;
    let last = *side
        .last()
        .ok_or_else(|| decline("wall split", "contact polyline is empty".into()))?;
    ensure_boundary_vertex(work, face, first, tolerance)?;
    ensure_boundary_vertex(work, face, last, tolerance)?;

    let surface = face_surface(work, face)?.clone();
    let uv_pts = project_uv(&surface, side)?;
    let uv_curve = uv_polyline_nurbs(&uv_pts);
    let curve3 = uniform_polyline_nurbs(side);
    let split = work
        .split_face(face, &uv_curve, &curve3, tolerance)
        .map_err(|e| decline("wall split", format!("{e}")))?;
    let edge = *split
        .new_edges
        .first()
        .ok_or_else(|| decline("wall split", "split_face produced no trim edge".into()))?;
    reconcile_edge_curve(work, edge, &curve3, tolerance);

    let mid = &chain[chain.len() / 2];
    let (anchor, uv, w3) = if is_side_a {
        (mid.contact0, mid.uv0, arc_departure(mid, true, bulge))
    } else {
        (mid.contact1, mid.uv1, arc_departure(mid, false, bulge))
    };
    let w3 = w3.ok_or_else(|| {
        decline(
            "wall trim side",
            "degenerate section arc at mid-chain".into(),
        )
    })?;
    let n = surface.normal_at(uv.0, uv.1);
    let w = (w3 - n * w3.dot(&n)).normalized().ok_or_else(|| {
        decline(
            "wall trim side",
            "arc departure is normal to the wall — witness degenerate".into(),
        )
    })?;
    let (doomed, kept) = pick_doomed(
        work,
        split.face_a,
        split.face_b,
        anchor,
        w,
        "wall trim side",
    )?;
    Ok(WallSplit { kept, doomed, edge })
}

/// `split_face` may canonicalize a straight trim polyline into an
/// analytic line; keep that recognition when its parameterization
/// agrees with the authored polyline, otherwise restore the polyline
/// as the single geometric authority.
pub(crate) fn reconcile_edge_curve(
    work: &mut BRep,
    edge: EdgeHandle,
    poly: &NurbsCurve,
    tolerance: &Tolerance,
) {
    let Some(curve) = work.edges[edge].curve.as_ref() else {
        return;
    };
    let (p0, p1) = poly.parameter_range();
    let mut mismatch = 0.0_f64;
    for f in [0.25, 0.5, 0.75] {
        let tc = curve.t_start + (curve.t_end - curve.t_start) * f;
        let tp = p0 + (p1 - p0) * f;
        mismatch = mismatch.max(curve.geometry.point_at(tc).distance_to(&poly.point_at(tp)));
    }
    if mismatch > tolerance.linear().max(1.0e-9) * 10.0 {
        work.edges[edge].curve = Some(CurveOnEdge {
            geometry: CurveGeometry::Nurbs(poly.clone()),
            t_start: p0,
            t_end: p1,
        });
    }
}

fn flood_doomed(
    work: &BRep,
    seed: FaceHandle,
    blocked: &[EdgeHandle],
    barriers: &[FaceHandle],
) -> Result<HashSet<FaceHandle>, BlendError> {
    let mut region = HashSet::new();
    region.insert(seed);
    let mut queue = vec![seed];
    while let Some(f) = queue.pop() {
        if !work.faces[f].inner_loops.is_empty() {
            return Err(decline(
                "doomed-region flood",
                format!("gap-region face {f:?} has inner loops — unsupported in stage 2"),
            ));
        }
        for he in loop_half_edges(work, work.faces[f].outer_loop) {
            let edge = work.half_edges[he].edge;
            if blocked.contains(&edge) {
                continue;
            }
            let twin = work.half_edges[he].twin;
            let neighbor = work.loops[work.half_edges[twin].loop_].face;
            if barriers.contains(&neighbor) {
                continue;
            }
            if region.insert(neighbor) {
                queue.push(neighbor);
            }
        }
    }
    Ok(region)
}

fn merge_region(
    work: &mut BRep,
    mut region: HashSet<FaceHandle>,
    blocked: &[EdgeHandle],
) -> Result<FaceHandle, BlendError> {
    for &f in &region {
        work.faces[f].surface = None;
    }
    loop {
        let mut candidate: Option<HalfEdgeHandle> = None;
        'outer: for &f in &region {
            for he in loop_half_edges(work, work.faces[f].outer_loop) {
                let edge = work.half_edges[he].edge;
                if blocked.contains(&edge) {
                    continue;
                }
                let twin = work.half_edges[he].twin;
                let other = work.loops[work.half_edges[twin].loop_].face;
                if other != f && region.contains(&other) {
                    candidate = Some(he);
                    break 'outer;
                }
            }
        }
        let Some(he) = candidate else {
            break;
        };
        let merged = work
            .merge_faces_across_edge(he)
            .map_err(|e| decline("region merge", format!("{e}")))?;
        region.remove(&merged.removed_face);
    }
    if region.len() != 1 {
        return Err(decline(
            "region merge",
            format!("merging left {} disconnected faces", region.len()),
        ));
    }
    let face = *region
        .iter()
        .next()
        .ok_or_else(|| decline("region merge", "merged region vanished".into()))?;

    // Dissolve leftover interior spur trees tip-first.
    loop {
        let mut tip: Option<HalfEdgeHandle> = None;
        for he in loop_half_edges(work, work.faces[face].outer_loop) {
            if work.half_edges[he].next == work.half_edges[he].twin {
                tip = Some(he);
                break;
            }
        }
        let Some(he) = tip else {
            break;
        };
        let tip_vertex = work.half_edges[work.half_edges[he].twin].origin;
        work.dissolve_same_loop_bridge(he)
            .map_err(|e| decline("spur dissolve", format!("{e}")))?;
        if work.half_edges.iter().all(|(_, d)| d.origin != tip_vertex) {
            work.vertices.remove(tip_vertex);
        }
    }
    Ok(face)
}

fn resurface_strip(
    work: &mut BRep,
    face: FaceHandle,
    skin: &StripSkin,
    chain: &[ContactSolution],
    tolerance: &Tolerance,
) -> Result<(), BlendError> {
    let hes = loop_half_edges(work, work.faces[face].outer_loop);
    if hes.len() != 4 {
        return Err(decline(
            "strip seam",
            format!(
                "merged gap boundary has {} half-edges, expected 4",
                hes.len()
            ),
        ));
    }

    let surface = SurfaceGeometry::Nurbs(skin.surface.clone());
    let q = surface.point_at(0.5, 0.5);
    let c = chain[chain.len() / 2].center;
    let n = surface.normal_at(0.5, 0.5);
    let orientation = n.dot(&(c - q)) >= 0.0;
    work.faces[face].surface = Some(SurfaceOnFace {
        geometry: surface,
        orientation,
    });

    // Chart borders by endpoint matching (mirrors the sheet builder;
    // loop order after the merges is an implementation detail).
    type Border<'a> = (&'a [Point3], ((f64, f64), (f64, f64)));
    let sides: [Border; 4] = [
        (&skin.side_a, ((0.0, 0.0), (1.0, 0.0))),
        (&skin.end_arc, ((1.0, 0.0), (1.0, 1.0))),
        (&skin.side_b, ((0.0, 1.0), (1.0, 1.0))),
        (&skin.start_arc, ((0.0, 0.0), (0.0, 1.0))),
    ];
    let tol_match = tolerance.linear().max(1.0e-9) * 100.0;
    for &he in &hes {
        let he_data = &work.half_edges[he];
        let origin = work.vertices[he_data.origin].position;
        let dest = work.vertices[work.half_edges[he_data.twin].origin].position;
        let border = sides.iter().find_map(|(pts, border)| {
            let first = *pts.first()?;
            let last = *pts.last()?;
            if first.distance_to(&origin) < tol_match && last.distance_to(&dest) < tol_match {
                Some((border.0, border.1))
            } else if first.distance_to(&dest) < tol_match && last.distance_to(&origin) < tol_match
            {
                Some((border.1, border.0))
            } else {
                None
            }
        });
        let (b0, b1) = border.ok_or_else(|| {
            decline(
                "strip seam",
                format!("no chart border matches the boundary edge at {origin:?}"),
            )
        })?;
        work.half_edges[he].pcurve = Some(Pcurve::linear(b0, b1, 0.0, 1.0));
    }
    Ok(())
}
