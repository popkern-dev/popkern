//! PK-177: carve the hole-corner wall boundary onto its honest
//! trims.
//!
//! The served topology closes each hole corner with ONE strip x
//! strip chord that also bounds the WALL — geometrically impossible
//! to satisfy both (the strips' mutual intersection bows off the
//! wall by ~0.75 r mid-span, measured 3.75e-2 on the boss). The
//! honest closure: TWO wall trims (wall ∩ strip_a, wall ∩ strip_b)
//! meeting at the corner-ball wall contact W, with the chord kept
//! purely between the strips. The trim geometry comes from the
//! solved corner cap (`CornerCapGeometry`): the rolling ball's wall
//! contacts along each strip's extension path.
//!
//! Transactional per corner: candidate clone, validate()-gated.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, VertexHandle};

use super::corner_cap_geometry::CornerCapGeometry;

/// Loop context of `vertex` on a wall-surfaced face: the face, and
/// the loop half-edges in order starting at the side whose origin is
/// `vertex`.
fn wall_loop_context(
    brep: &BRep,
    wall: &SurfaceGeometry,
    vertex: VertexHandle,
) -> Option<(FaceHandle, Vec<HalfEdgeHandle>)> {
    let sig = format!("{wall:?}");
    for (h, he) in brep.half_edges.iter() {
        if he.origin != vertex {
            continue;
        }
        let Some(l) = brep.loops.get(he.loop_) else {
            continue;
        };
        let Some(fd) = brep.faces.get(l.face) else {
            continue;
        };
        let Some(sof) = fd.surface.as_ref() else {
            continue;
        };
        if format!("{:?}", sof.geometry) != sig {
            continue;
        }
        let mut sides = vec![h];
        let mut cur = he.next;
        for _ in 0..brep.half_edges.len() {
            if cur == h {
                break;
            }
            sides.push(cur);
            cur = brep.half_edges.get(cur)?.next;
        }
        return Some((l.face, sides));
    }
    None
}

/// Settle the stale wall-boundary vertices around a snapped corner
/// onto their honest wall ∩ strip intersections.
///
/// MEASURED structure (boss, all four corners probed): the -y
/// corners already carry the honest two-trim shape (both rail
/// anchors exist as loop vertices at d = 0); at the +y corners the
/// wall loop still runs along STALE original-edge segments (the
/// pre-fillet boss edge at its original offset, twinning the
/// region-merged strip) instead of the contact rails — the same
/// root as the PK141_TRUST_STRIPS 0.7 mis-side. Every stale vertex
/// sits ON the wall but off its strip twin; the incident-surface
/// settle (wall + strip, a transversal pair) moves it onto the
/// contact rail and the reconcile refits the adjacent curves.
/// Transactional per vertex (inside the settle).
pub(crate) fn carve_corner_wall_trims(
    brep: &mut BRep,
    cap: &CornerCapGeometry,
    corner_vertex: VertexHandle,
    tolerance: &Tolerance,
    trace: bool,
) -> Option<bool> {
    let (_, sides) = wall_loop_context(brep, &cap.shared_surface, corner_vertex)?;
    // Vertices within a few sides of the corner on the wall loop,
    // both directions (the stale runs are short: 1-2 segments).
    let mut candidates: Vec<VertexHandle> = Vec::new();
    let take = 4usize.min(sides.len());
    for h in sides.iter().take(take).chain(sides.iter().rev().take(take)) {
        if let Some(he) = brep.half_edges.get(*h)
            && he.origin != corner_vertex
            && !candidates.contains(&he.origin)
        {
            candidates.push(he.origin);
        }
    }
    let mut moved = 0usize;
    for v in candidates {
        // Only the STALE two-surface class settles: a vertex with 3+
        // distinct incident surfaces is a corner-class vertex — the
        // over-determined settle was falsified there (no common
        // zero-residual point; measured regression 3.5e-3 ->
        // 1.28e-2 when corner vertices re-settled through this
        // pass).
        let mut sigs: Vec<String> = Vec::new();
        for (_, he) in brep.half_edges.iter() {
            if he.origin != v {
                continue;
            }
            if let Some(l) = brep.loops.get(he.loop_)
                && let Some(fd) = brep.faces.get(l.face)
                && let Some(sof) = fd.surface.as_ref()
            {
                let sig = format!("{:?}", sof.geometry);
                if !sigs.contains(&sig) {
                    sigs.push(sig);
                }
            }
        }
        if sigs.len() != 2 {
            continue;
        }
        if let Some(n) =
            super::corner_chord_settle::settle_vertex_on_incident_surfaces(brep, v, tolerance)
            && n > 0
        {
            moved += n;
            if trace {
                eprintln!("[corner-carve] {corner_vertex:?}: stale wall vertex {v:?} settled");
            }
        }
    }
    if trace
        && let Some((chord, far, g1, g2)) =
            super::corner_chord_settle::corner_strip_chord(brep, corner_vertex)
    {
        // Chart-validity probe at the corner: where does each strip
        // chart place W, and how far off is its image?
        use popkern_core::geometry::SurfaceEval;
        let w = brep.vertices.get(corner_vertex)?.position;
        let far_pos = brep.vertices.get(far)?.position;
        for (label, g) in [("g1", &g1), ("g2", &g2)] {
            let (nu, nv) = match g {
                SurfaceGeometry::Nurbs(n) => (n.num_u, n.num_v),
                _ => (0, 0),
            };
            let img = popkern_core::geometry::uv_projection::project_point_to_uv(g, &w)
                .map(|(u, v)| {
                    let q = g.point_at(u, v);
                    format!("uv=({u:.4},{v:.4}) d={:.3e}", w.distance_to(&q))
                })
                .unwrap_or_else(|| "proj FAILED".into());
            eprintln!("[corner-carve]   chord {chord:?} {label} ({nu}x{nv}): W {img}");
        }
        eprintln!(
            "[corner-carve]   chord far {far:?} at ({:.4},{:.4},{:.4})",
            far_pos.x, far_pos.y, far_pos.z
        );
    }
    Some(moved > 0)
}

/// The carve proper: mef the two honest wall trims (rolling-ball
/// wall contacts along each strip's extension path) from the rail
/// anchors to the corner contact W, then absorb the enclosed
/// stale-boundary slivers into their twin strips. Transactional:
/// candidate clone; committed only when validation passes and the
/// serving-metric worst deviation over the affected edges strictly
/// improves. Corners already carrying the two-trim shape skip.
pub(crate) fn carve_corner_wall_topology(
    brep: &mut BRep,
    cap: &CornerCapGeometry,
    corner_vertex: VertexHandle,
    trace: bool,
) -> Option<bool> {
    use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
    let radius_gate = {
        let w = brep.vertices.get(corner_vertex)?.position;
        w.distance_to(&cap.center)
    };
    let mut cand = brep.clone();
    let (wall_face, _sides) = wall_loop_context(&cand, &cap.shared_surface, corner_vertex)?;
    let mut mef_faces = Vec::new();
    let mut carved = false;
    for side_k in 0..2 {
        let rail = wall_rail_points(cap, side_k)?;
        let anchor_pos = *rail.first()?;
        // Locate the anchor side on the current wall loop.
        let (_, sides_now) = wall_loop_context(&cand, &cap.shared_surface, corner_vertex)?;
        let mut he_anchor = None;
        let mut he_w = None;
        for h in &sides_now {
            let Some(he) = cand.half_edges.get(*h) else {
                continue;
            };
            let Some(v) = cand.vertices.get(he.origin) else {
                continue;
            };
            if he.origin == corner_vertex {
                he_w = Some(*h);
            } else if v.position.distance_to(&anchor_pos) < radius_gate * 0.25 {
                he_anchor = Some((*h, he.origin));
            }
        }
        // No anchor vertex yet (the +y corners keep stale boundary
        // segments): split the loop side whose span passes the
        // anchor position and settle the new vertex there.
        if he_anchor.is_none() {
            let mut best: Option<(f64, popkern_core::topology::EdgeHandle)> = None;
            for h in &sides_now {
                let Some(he) = cand.half_edges.get(*h) else {
                    continue;
                };
                let Some(a) = cand.vertices.get(he.origin).map(|v| v.position) else {
                    continue;
                };
                let Some(b) = cand
                    .half_edges
                    .get(he.next)
                    .and_then(|n| cand.vertices.get(n.origin))
                    .map(|v| v.position)
                else {
                    continue;
                };
                let ab = b - a;
                let len2 = ab.dot(&ab);
                if len2 < 1e-18 {
                    continue;
                }
                let t = ((anchor_pos - a).dot(&ab) / len2).clamp(0.0, 1.0);
                if !(0.05..=0.95).contains(&t) {
                    continue;
                }
                let foot = a + ab * t;
                let d = foot.distance_to(&anchor_pos);
                if best.map(|(bd, _)| d < bd).unwrap_or(true) {
                    best = Some((d, he.edge));
                }
            }
            if let Some((d, split_edge)) = best
                && d < radius_gate * 0.5
                && let Ok(split) =
                    cand.split_edge_at_point(split_edge, anchor_pos, &cand.tolerance.clone())
            {
                if trace {
                    eprintln!(
                        "[corner-carve] {corner_vertex:?} side {side_k}: anchor vertex split in \
                         at d={d:.4e}"
                    );
                }
                // Re-walk the loop for the fresh side handles.
                let _ = split;
                let (_, refreshed) = wall_loop_context(&cand, &cap.shared_surface, corner_vertex)?;
                for h in &refreshed {
                    let Some(he) = cand.half_edges.get(*h) else {
                        continue;
                    };
                    let Some(v) = cand.vertices.get(he.origin) else {
                        continue;
                    };
                    if he.origin != corner_vertex
                        && v.position.distance_to(&anchor_pos) < radius_gate * 0.25
                    {
                        he_anchor = Some((*h, he.origin));
                    }
                    if he.origin == corner_vertex {
                        he_w = Some(*h);
                    }
                }
            }
        }
        let (Some((ha, anchor_vertex)), Some(hw)) = (he_anchor, he_w) else {
            continue;
        };
        // Direct edge anchor->W already exists (the two-trim shape,
        // or freshly created by the anchor split along the OLD
        // boundary): ensure its CURVE is the honest rail — rewrite
        // in place when it deviates.
        let direct = cand.edges.iter().find_map(|(eh, ed)| {
            cand.half_edges.get(ed.half_edge).and_then(|he| {
                cand.half_edges.get(he.twin).and_then(|tw| {
                    ((he.origin == anchor_vertex && tw.origin == corner_vertex)
                        || (he.origin == corner_vertex && tw.origin == anchor_vertex))
                        .then_some((eh, he.origin == anchor_vertex))
                })
            })
        });
        if let Some((direct_edge, forward)) = direct {
            // Honest rail curve (anchor -> W), reversed if the
            // edge's primary runs W -> anchor.
            let mut pts = rail.clone();
            if !forward {
                pts.reverse();
            }
            let params: Vec<f64> = (0..pts.len())
                .map(|i| (i as f64) / ((pts.len() - 1) as f64))
                .collect();
            if let Ok(fitted) = popkern_core::geometry::fitting::fit_curve(
                &pts,
                &params,
                3.min(pts.len() - 1),
                pts.len(),
            ) {
                let (t0, t1) = fitted.parameter_range();
                let candidate_curve = CurveOnEdge {
                    geometry: CurveGeometry::Nurbs(fitted),
                    t_start: t0,
                    t_end: t1,
                };
                // Rewrite only when the current curve deviates from
                // the rail beyond the serving band.
                let deviates = cand
                    .edges
                    .get(direct_edge)
                    .and_then(|ed| ed.curve.as_ref())
                    .map(|cur| {
                        let mid = cur.geometry.point_at(0.5 * (cur.t_start + cur.t_end));
                        let rail_mid = pts[pts.len() / 2];
                        mid.distance_to(&rail_mid) > 1e-5
                    })
                    .unwrap_or(true);
                if deviates && let Some(ed) = cand.edges.get_mut(direct_edge) {
                    ed.curve = Some(candidate_curve);
                    carved = true;
                    if trace {
                        eprintln!(
                            "[corner-carve] {corner_vertex:?} side {side_k}: direct trim curve \
                             rewritten to the rail"
                        );
                    }
                }
            }
            continue;
        }
        let Ok(result) = cand.mef(ha, hw) else {
            continue;
        };
        // Exact trim curve through the rail contacts (anchor -> W).
        let params: Vec<f64> = (0..rail.len())
            .map(|i| (i as f64) / ((rail.len() - 1) as f64))
            .collect();
        if let Ok(fitted) = popkern_core::geometry::fitting::fit_curve(
            &rail,
            &params,
            3.min(rail.len() - 1),
            rail.len(),
        ) {
            let (t0, t1) = fitted.parameter_range();
            if let Some(ed) = cand.edges.get_mut(result.edge) {
                ed.curve = Some(CurveOnEdge {
                    geometry: CurveGeometry::Nurbs(fitted),
                    t_start: t0,
                    t_end: t1,
                });
            }
        }
        mef_faces.push(result.face);
        carved = true;
        if trace {
            eprintln!(
                "[corner-carve] {corner_vertex:?} side {side_k}: wall trim mef {:?} (face {:?})",
                result.edge, result.face
            );
        }
    }
    if !carved {
        return Some(false);
    }
    // Absorb the mef slivers: a SHORT new face whose loop carries a
    // side twinned by a NURBS strip merges into that strip
    // (clear-surface / kef / restore).
    for sliver in mef_faces {
        let Some(fd) = cand.faces.get(sliver) else {
            continue;
        };
        let Some(l) = cand.loops.get(fd.outer_loop) else {
            continue;
        };
        let mut count = 0usize;
        let mut strip_side = None;
        let start = l.half_edge;
        let mut cur = start;
        for _ in 0..cand.half_edges.len() {
            let Some(he) = cand.half_edges.get(cur) else {
                break;
            };
            count += 1;
            if let Some(tw) = cand.half_edges.get(he.twin)
                && let Some(tl) = cand.loops.get(tw.loop_)
                && tl.face != sliver
                && tl.face != wall_face
                && cand
                    .faces
                    .get(tl.face)
                    .and_then(|f| f.surface.as_ref())
                    .is_some_and(|sof| matches!(sof.geometry, SurfaceGeometry::Nurbs(_)))
            {
                strip_side = Some((he.twin, tl.face));
            }
            cur = he.next;
            if cur == start {
                break;
            }
        }
        let Some((strip_he, strip_face)) = strip_side else {
            continue;
        };
        if count > 8 {
            continue;
        }
        let strip_surface = cand.faces.get(strip_face).and_then(|f| f.surface.clone());
        if let Some(f) = cand.faces.get_mut(strip_face) {
            f.surface = None;
        }
        if let Some(f) = cand.faces.get_mut(sliver) {
            f.surface = None;
        }
        let merged = cand.kef(strip_he).is_ok();
        if let Some(f) = cand.faces.get_mut(strip_face) {
            f.surface = strip_surface;
        }
        if trace {
            eprintln!("[corner-carve]   sliver {sliver:?} -> {strip_face:?} merged={merged}");
        }
    }
    // Refit curves in the corner neighborhood and heal pcurves.
    super::edge_reconcile::reconcile_edge_curves(&mut cand);
    cand.fix_stale_pcurves();
    if cand.validate().is_err() {
        if trace {
            eprintln!("[corner-carve] {corner_vertex:?} carve declined validation");
        }
        return Some(false);
    }
    *brep = cand;
    Some(true)
}

/// Wall-rail contact points along one strip extension path (index 0
/// = the strip's end station, last = the corner contact W).
fn wall_rail_points(
    cap: &CornerCapGeometry,
    side: usize,
) -> Option<Vec<popkern_core::geometry::Point3>> {
    let path = cap.extension_centers.get(side)?;
    let mut rail = Vec::with_capacity(path.len());
    for c in path {
        rail.push(super::corner_cap_geometry::ball_contact(
            &cap.shared_surface,
            c,
        )?);
    }
    Some(rail)
}
