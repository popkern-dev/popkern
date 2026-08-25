//! PK-174 — recordless corner-group unfillet.
//!
//! A trihedral corner filleted on all three edges serves three
//! strips meeting one corner cap. Undoing it without an operation
//! record is a GROUP surgery — the strips' ends at the cap have no
//! planar end face, and the three restored sharp edges meet at a
//! vertex that exists nowhere in the blended body: the
//! triple-intersection of the three support planes.
//!
//! Surgery (inverse of the corner construction, transactional on a
//! clone):
//! 1. Discover the group from any member strip: cap = the
//!    Sphere/Nurbs/Gregory blend neighbour; strips = the cap's
//!    Cylinder/Torus blend neighbours; supports = the strips'
//!    support planes (exactly three distinct planes, else a
//!    structured decline names the class).
//! 2. Merge the whole group into one region face with cross-face
//!    kefs. The group adjacency is a 6-cycle (support-strip ring)
//!    plus the cap, so one strip contact edge survives as a
//!    same-face slit — after the corner weld and curve rebuild that
//!    slit IS one of the three restored edges, re-legitimized when
//!    the two mef splits put its sides in different loops.
//! 3. Weld the corner cluster (cap ring + strip corner arcs + old
//!    corner remnants) into one vertex at the triple intersection;
//!    weld each strip's far cluster onto its restored far endpoint.
//! 4. Two mefs split the region into the three planar sectors;
//!    surfaces are assigned by measured plane deviation, restored
//!    edges carry Line curves, pcurves heal through the standard
//!    passes.

use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

use crate::error::BlendError;
use crate::unfillet_corner_group_geometry::{
    cylinder_support_pair, edge_between, face_diameter, intersect_three_planes, invalid,
    is_cap_surface, is_strip_surface, loop_vertices_near, needs_mef, pick_mef_pair,
    plane_deviation, surface_kind_of, weld_cluster,
};
use crate::unfillet_helpers::{
    PlaneData, collect_end_face_planes, extract_plane, intersect_line_plane, intersect_planes,
};
use crate::unfillet_surgery::{
    rebuild_curves_touching, remove_degenerate_self_loop_edges, remove_slit_loop_edges,
};

/// Outcome of a corner-group restore: the restored body plus the
/// restored sharp edges with the radii their strips carried (the
/// recordless MODIFY leg re-fillets the group from this mapping).
pub(crate) struct CornerGroupRestore {
    pub brep: BRep,
    /// (restored edge, the strip's old radius, strip was the seed).
    pub edges: Vec<(EdgeHandle, f64, bool)>,
}

/// Recordless corner-group restore. `seed` is any strip of the
/// group. Returns the restored body or a structured decline naming
/// what the group discovery or surgery could not establish.
pub(crate) fn unfillet_corner_group(brep: &BRep, seed: FaceHandle) -> Result<BRep, BlendError> {
    restore_corner_group(brep, seed).map(|r| r.brep)
}

/// Full restore with the edge/radius mapping (see
/// [`CornerGroupRestore`]).
pub(crate) fn restore_corner_group(
    brep: &BRep,
    seed: FaceHandle,
) -> Result<CornerGroupRestore, BlendError> {
    // --- 1. Group discovery -------------------------------------------------
    let caps: Vec<FaceHandle> = crate::identify::blend_neighbours(brep, seed)
        .into_iter()
        .filter(|&f| surface_kind_of(brep, f).is_some_and(is_cap_surface))
        .collect();
    let Some(&cap) = caps.first() else {
        return Err(invalid(
            "blend_face",
            format!("{seed:?}"),
            "corner-group restore needs a cap neighbour, found none".to_string(),
        ));
    };
    // Strips adjacent to the cap directly, or through corner
    // FRAGMENTS — the mixed-radius (setback/Gregory) corner class
    // interposes small plane fragments and surfaceless residuals
    // between the cap and the strips. Traversal is bounded to two
    // fragment hops and to faces no larger than a few cap diameters
    // so it can never walk out through a support plane.
    let cap_diameter = face_diameter(brep, cap).max(1e-9);
    let mut strips: Vec<FaceHandle> = vec![seed];
    let mut visited: Vec<FaceHandle> = caps.clone();
    let mut frontier: Vec<FaceHandle> = caps.clone();
    for _ in 0..2 {
        let mut next: Vec<FaceHandle> = Vec::new();
        for &f in &frontier {
            for n in crate::identify::blend_neighbours(brep, f) {
                if visited.contains(&n) || strips.contains(&n) {
                    continue;
                }
                match surface_kind_of(brep, n) {
                    Some(g) if is_strip_surface(g) => strips.push(n),
                    // Sibling cap sectors (a Gregory corner is a fan
                    // of sector faces) are always traversable.
                    Some(g) if is_cap_surface(g) => {
                        visited.push(n);
                        next.push(n);
                    }
                    Some(SurfaceGeometry::Plane { .. }) | None
                        if face_diameter(brep, n) <= cap_diameter * 3.0 =>
                    {
                        visited.push(n);
                        next.push(n);
                    }
                    _ => {}
                }
            }
        }
        frontier = next;
    }
    strips.sort();
    strips.dedup();
    if strips.len() != 3 {
        return Err(invalid(
            "blend_face",
            format!("{seed:?}"),
            format!(
                "corner-group restore covers the trihedral class (3 strips at one cap), found {}",
                strips.len()
            ),
        ));
    }
    // Support pair per strip, and the distinct support set. The
    // supports are derived from the strip's own cylinder geometry
    // (plane neighbours whose normal is perpendicular to the fillet
    // axis at axis distance ~= radius) — the tangency-sample pair
    // from identify_blend picks the cap or a far face on capped
    // corner strips (measured on the trihedral cube fixture).
    let mut strip_supports: Vec<(FaceHandle, FaceHandle)> = Vec::new();
    for &s in &strips {
        let pair = cylinder_support_pair(brep, s).ok_or_else(|| {
            invalid(
                "blend_face",
                format!("{s:?}"),
                "could not determine a strip's tangent support plane pair".to_string(),
            )
        })?;
        strip_supports.push(pair);
    }
    let mut supports: Vec<FaceHandle> = strip_supports.iter().flat_map(|&(a, b)| [a, b]).collect();
    supports.sort();
    supports.dedup();
    if supports.len() != 3 {
        return Err(invalid(
            "blend_face",
            format!("{seed:?}"),
            format!(
                "trihedral corner-group restore needs 3 distinct support planes, found {}",
                supports.len()
            ),
        ));
    }
    let planes: Vec<PlaneData> = supports
        .iter()
        .map(|&f| extract_plane(brep, f))
        .collect::<Result<_, _>>()?;

    // --- 2. Corner vertex = triple plane intersection ----------------------
    let corner_pos = intersect_three_planes(&planes).ok_or_else(|| {
        invalid(
            "supports",
            format!("{supports:?}"),
            "support planes are near-parallel; no triple-intersection corner".to_string(),
        )
    })?;

    // Restored far endpoint + support line per strip.
    let mut restored: Vec<(Point3, Vec3)> = Vec::new(); // (far endpoint, line dir)
    let mut strip_radii: Vec<f64> = Vec::new();
    let mut max_radius = 0.0_f64;
    for (i, &s) in strips.iter().enumerate() {
        let mut this_radius = 0.0_f64;
        if let Some(SurfaceGeometry::Cylinder { radius, .. }) = surface_kind_of(brep, s) {
            max_radius = max_radius.max(*radius);
            this_radius = *radius;
        }
        strip_radii.push(this_radius);
        let (a, b) = strip_supports[i];
        let pa = extract_plane(brep, a)?;
        let pb = extract_plane(brep, b)?;
        let (line_origin, line_dir) = intersect_planes(&pa, &pb)?;
        let end_planes = collect_end_face_planes(brep, s, a, b)?;
        let far = end_planes
            .iter()
            .filter_map(|p| intersect_line_plane(&line_origin, &line_dir, p))
            .max_by(|p, q| {
                corner_pos
                    .distance_to(p)
                    .total_cmp(&corner_pos.distance_to(q))
            })
            .ok_or_else(|| {
                invalid(
                    "blend_face",
                    format!("{s:?}"),
                    "strip has no planar far end to bound the restored edge".to_string(),
                )
            })?;
        restored.push((far, line_dir));
    }

    // --- 3. Transactional surgery ------------------------------------------
    let mut out = brep.clone();
    let mut group: Vec<FaceHandle> = Vec::new();
    group.extend(&supports);
    group.extend(&strips);
    group.push(cap);
    // Corner fragments and residuals living entirely inside the
    // corner ball belong to the group (the mixed-radius class).
    let corner_ball = (max_radius * 3.0).max(1e-6);
    for (fh, _) in brep.faces.iter() {
        if group.contains(&fh) {
            continue;
        }
        let positions = brep.face_outer_loop_positions(fh);
        if !positions.is_empty()
            && positions
                .iter()
                .all(|p| p.distance_to(&corner_pos) <= corner_ball)
        {
            group.push(fh);
        }
    }

    // Merge the group into one face: repeatedly kef an edge whose two
    // sides are distinct faces of the group. kef kills the twin side's
    // face, so both sides stay in the shrinking group set.
    // Surfaceless residual faces (the mixed-radius/setback class
    // ships ghost-ridden 2-3 side residuals) refuse standard kef;
    // remove them with the consolidation's tolerant removal before
    // merging the surfaced group.
    let residuals: Vec<FaceHandle> = group
        .iter()
        .copied()
        .filter(|&f| out.faces.get(f).is_some_and(|d| d.surface.is_none()))
        .collect();
    for f in residuals {
        crate::apex_construction::v_free_shell_consolidation::remove_face(&mut out, f);
        group.retain(|&g| g != f);
    }
    // Ghost sweep: the residual class carries half-edges whose loop_
    // points at the removed loops without being reachable in the
    // ring walk (the measured "ghost twin records"); remove_face's
    // ring walk cannot see them. Anything referencing a dead loop
    // goes, with its edge when both sides are gone.
    let ghosts: Vec<HalfEdgeHandle> = out
        .half_edges
        .iter()
        .filter(|(_, he)| !out.loops.contains_key(he.loop_))
        .map(|(h, _)| h)
        .collect();
    for h in ghosts {
        let Some(data) = out.remove_half_edge_reanchored(h) else {
            continue;
        };
        if out.half_edges.contains_key(data.twin) {
            out.boundary_edges.insert(data.edge);
            if let Some(ed) = out.edges.get_mut(data.edge)
                && ed.half_edge == h
            {
                ed.half_edge = data.twin;
            }
        } else if out.edges.remove(data.edge).is_some() {
            out.boundary_edges.remove(&data.edge);
        }
    }
    let mut merges = 0usize;
    loop {
        // Candidate cross-face edges inside the group. kef kills the
        // TWIN side's face, so orient each candidate to kill the
        // surfaceless residual first (the mixed-radius class ships
        // ghost-ridden 2-side residual faces; killing them into a
        // surfaced neighbour is the only kef direction that holds).
        let mut best: Option<(HalfEdgeHandle, bool)> = None;
        for (_, ed) in out.edges.iter() {
            let Some(he) = out.half_edges.get(ed.half_edge) else {
                continue;
            };
            let Some(tw) = out.half_edges.get(he.twin) else {
                continue;
            };
            let (Some(la), Some(lb)) = (out.loops.get(he.loop_), out.loops.get(tw.loop_)) else {
                continue;
            };
            let (fa, fb) = (la.face, lb.face);
            if fa == fb || !group.contains(&fa) || !group.contains(&fb) {
                continue;
            }
            let fa_res = out.faces.get(fa).is_some_and(|d| d.surface.is_none());
            let fb_res = out.faces.get(fb).is_some_and(|d| d.surface.is_none());
            // he kills fb; he.twin kills fa. Prefer killing a
            // residual side.
            let (cand, kills_residual) = if fb_res {
                (ed.half_edge, true)
            } else if fa_res {
                (he.twin, true)
            } else {
                (ed.half_edge, false)
            };
            if kills_residual {
                best = Some((cand, true));
                break;
            }
            if best.is_none() {
                best = Some((cand, false));
            }
        }
        let Some((he, _)) = best else { break };
        let killed = out
            .half_edges
            .get(out.half_edges[he].twin)
            .and_then(|tw| out.loops.get(tw.loop_).map(|l| l.face));
        out.kef(he).map_err(BlendError::Topology)?;
        if let Some(k) = killed {
            group.retain(|&f| f != k);
        }
        merges += 1;
        if merges > 64 {
            return Err(invalid(
                "group",
                format!("{group:?}"),
                "corner-group merge did not converge".to_string(),
            ));
        }
    }
    let [merged] = group.as_slice() else {
        return Err(invalid(
            "group",
            format!("{group:?}"),
            "corner-group did not merge to a single region".to_string(),
        ));
    };
    let merged = *merged;

    // Weld the corner cluster into one vertex at the triple point.
    let corner_gate = (max_radius * 2.5).max(1e-6);
    let corner_cluster = loop_vertices_near(&out, merged, &corner_pos, corner_gate);
    if corner_cluster.len() < 2 {
        return Err(invalid(
            "corner",
            format!("{corner_pos:?}"),
            format!(
                "corner cluster has {} vertices; the cap ring did not collapse",
                corner_cluster.len()
            ),
        ));
    }
    let corner_v = weld_cluster(&mut out, &corner_cluster, corner_pos)?;

    // Weld each strip's far cluster onto its restored endpoint.
    let mut far_vs: Vec<VertexHandle> = Vec::new();
    let mut welded = vec![corner_v];
    for &(far, _) in &restored {
        let far_gate = (max_radius * 1.5).max(1e-6);
        let cluster = loop_vertices_near(&out, merged, &far, far_gate);
        if cluster.is_empty() {
            return Err(invalid(
                "far_endpoint",
                format!("{far:?}"),
                "no far cluster vertices near the restored endpoint".to_string(),
            ));
        }
        let v = weld_cluster(&mut out, &cluster, far)?;
        far_vs.push(v);
        welded.push(v);
    }

    remove_degenerate_self_loop_edges(&mut out);
    rebuild_curves_touching(&mut out, &welded);

    // --- 4. Two mefs restore the sector split ------------------------------
    // The third restored edge is the surviving slit; the mefs put its
    // two sides into different sector loops.
    let mut mef_edges: Vec<(EdgeHandle, usize)> = Vec::new();
    let mut sector_faces: Vec<FaceHandle> = vec![merged];
    for (i, &far_v) in far_vs.iter().enumerate() {
        if !needs_mef(&out, corner_v, far_v) {
            continue;
        }
        let (he_c, he_f) = pick_mef_pair(&out, corner_v, far_v).ok_or_else(|| {
            invalid(
                "mef",
                format!("{corner_v:?}->{far_v:?}"),
                "no same-loop half-edge pair for the restored-edge mef".to_string(),
            )
        })?;
        let r = out.mef(he_c, he_f).map_err(BlendError::Topology)?;
        mef_edges.push((r.edge, i));
        sector_faces.push(r.face);
    }
    // Residual zero-area slit loops (redundant group contact edges
    // the mefs did not re-legitimize) must go BEFORE the restored-
    // edge mapping: `edge_between` would otherwise pick a slit
    // coincident with a real restored edge.
    remove_slit_loop_edges(&mut out);
    for &(edge, i) in &mef_edges {
        let (far, _) = restored[i];
        let dir = far - corner_pos;
        if dir.length() > 1e-12 {
            out.edges[edge].curve = Some(CurveOnEdge {
                geometry: CurveGeometry::Line {
                    origin: corner_pos,
                    direction: dir,
                },
                t_start: 0.0,
                t_end: 1.0,
            });
        }
    }

    // Assign each sector face (the merged parent + the mef children,
    // which inherit the parent's plane) the support plane it
    // geometrically measures against.
    let source_planes: Vec<SurfaceOnFace> = supports
        .iter()
        .filter_map(|&f| brep.faces.get(f).and_then(|d| d.surface.clone()))
        .collect();
    for f in sector_faces {
        let positions = out.face_outer_loop_positions(f);
        if positions.is_empty() {
            continue;
        }
        let best = source_planes.iter().min_by(|p, q| {
            plane_deviation(&p.geometry, &positions)
                .total_cmp(&plane_deviation(&q.geometry, &positions))
        });
        if let Some(sof) = best {
            out.faces[f].surface = Some(sof.clone());
        }
    }

    out.recompute_all_pcurves_lenient();
    out.fix_stale_pcurves();
    // Transactional postcondition: the restore must serve a fully
    // valid body. The setback/Gregory corner class (spanning sector
    // faces, far-corner residuals) reaches here with an incomplete
    // group and fails — that class's commercial path is the
    // operation record (PK-143), and this decline routes there.
    let errors = out.validate_all();
    if !errors.is_empty() {
        return Err(invalid(
            "restored_body",
            format!("{} validation errors", errors.len()),
            format!(
                "corner-group restore postcondition failed (first: {:?}); the setback-class                  corner is served by the operation record path",
                errors.first()
            ),
        ));
    }
    // Restored-edge mapping: mef-created edges are tracked; the
    // surviving-slit sector edge is found between the corner vertex
    // and its far vertex.
    let mut edges: Vec<(EdgeHandle, f64, bool)> = Vec::new();
    for (i, &far_v) in far_vs.iter().enumerate() {
        let handle = mef_edges
            .iter()
            .find(|(_, j)| *j == i)
            .map(|(e, _)| *e)
            .or_else(|| edge_between(&out, corner_v, far_v));
        let Some(e) = handle else {
            return Err(invalid(
                "restored_edges",
                format!("{corner_v:?}->{far_v:?}"),
                "restored sharp edge not found after surgery".to_string(),
            ));
        };
        edges.push((e, strip_radii[i], strips[i] == seed));
    }
    Ok(CornerGroupRestore { brep: out, edges })
}
