//! Setback vertex-blend topology surgery (PK-094 stage 2).
//!
//! Rebuilds the apex blend region's TOPOLOGY around the classical
//! corner-setback boundary: each fillet face is trimmed at its
//! setback arc, each lateral face keeps its wedge only up to the
//! trim-intersection corner c, and the surrendered apex region is
//! merged into a single n-gon blend face bounded by the n arcs.
//! This is the boundary the KP20 setback cap is built for (stage 1),
//! so downstream cap acceptance measures against the right edges.
//!
//! Work-in-progress gate: the arm only runs when
//! [`SETBACK_TOPOLOGY_WIP`] is set (integration tests); the cap face
//! fan (splitting the blend n-gon into the 4n KP20 piece faces) is
//! the recorded follow-up, so production keeps the legacy path until
//! that lands.

use std::sync::atomic::{AtomicBool, Ordering};

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

use super::kp20_paper_border_extract::{Kp20SetbackArc, kp20_paper_border_from_fillet_setback};
use crate::fillet_surface::FilletSurfaceResult;

/// WIP gate for the setback topology arm (see module docs).
pub(crate) static SETBACK_TOPOLOGY_WIP: AtomicBool = AtomicBool::new(false);

/// WIP checkpoint report channel for the integration probe test.
pub(crate) static SETBACK_TOPOLOGY_LAST: std::sync::Mutex<Option<String>> =
    std::sync::Mutex::new(None);

/// Result of the setback boundary surgery.
#[derive(Debug)]
pub(crate) struct SetbackTopologyOutcome {
    /// The merged blend face bounded by the n setback arcs.
    pub blend_face: FaceHandle,
    /// The n setback-arc edges, fillet order.
    pub arc_edges: Vec<EdgeHandle>,
    /// The n corner vertices (trim intersections), fillet order
    /// (corner k = between fillet k and fillet k+1).
    pub corner_vertices: Vec<VertexHandle>,
}

/// Attempt the setback boundary surgery at `vertex`. On any error the
/// BRep is restored from a snapshot and the error is returned for the
/// caller's structured diagnostics; the legacy path then proceeds.
pub(crate) fn try_setback_boundary_surgery(
    brep: &mut BRep,
    vertex: VertexHandle,
    fillet_surfaces: &[FilletSurfaceResult],
    tolerance: &Tolerance,
) -> Result<SetbackTopologyOutcome, String> {
    let snapshot = brep.clone();
    match surgery_inner(brep, vertex, fillet_surfaces, tolerance) {
        Ok(outcome) => Ok(outcome),
        Err(e) => {
            *brep = snapshot;
            Err(e)
        }
    }
}

fn surgery_inner(
    brep: &mut BRep,
    vertex: VertexHandle,
    fillet_surfaces: &[FilletSurfaceResult],
    tolerance: &Tolerance,
) -> Result<SetbackTopologyOutcome, String> {
    let n = fillet_surfaces.len();
    if n < 3 {
        return Err(format!("setback surgery needs >= 3 fillets, got {n}"));
    }
    let vertex_position = brep
        .vertices
        .get(vertex)
        .map(|v| v.position)
        .ok_or_else(|| "apex vertex not in BRep".to_string())?;

    // Setback arcs from the extractor (the same geometry the stage-1
    // candidate assembles the cap on).
    let (_border, extract) =
        kp20_paper_border_from_fillet_setback(fillet_surfaces, vertex_position)?;
    if extract.setback_arcs.len() != n {
        return Err(format!(
            "extractor returned {} arcs for {} fillets",
            extract.setback_arcs.len(),
            n
        ));
    }
    let mut arcs: Vec<&Kp20SetbackArc> = Vec::with_capacity(n);
    for k in 0..n {
        arcs.push(
            extract
                .setback_arcs
                .iter()
                .find(|a| a.fillet_index == k)
                .ok_or_else(|| format!("no setback arc for fillet {k}"))?,
        );
    }

    // WIP probe: where does the arc midpoint live relative to the
    // current face partition?
    if std::env::var_os("PK_SETBACK_PROBE").is_some() {
        for k in 0..n.min(2) {
            let mid = eval_cubic_mid(&arcs[k].cubic);
            let mut nearest = (f64::INFINITY, String::new());
            for (fh, face) in brep.faces.iter() {
                let Some(surface) = face.surface.as_ref() else {
                    continue;
                };
                let kind = match &surface.geometry {
                    popkern_core::geometry::SurfaceGeometry::Plane { .. } => "Plane",
                    popkern_core::geometry::SurfaceGeometry::Nurbs(_) => "Nurbs",
                    popkern_core::geometry::SurfaceGeometry::Cylinder { .. } => "Cylinder",
                    _ => "Other",
                };
                if let Some((u, v)) = popkern_core::geometry::uv_projection::project_point_to_uv(
                    &surface.geometry,
                    &mid,
                ) {
                    let d = (&surface.geometry as &dyn SurfaceEval)
                        .point_at(u, v)
                        .distance_to(&mid);
                    if d < nearest.0 {
                        nearest = (d, format!("{fh:?} {kind}"));
                    }
                }
            }
            eprintln!(
                "  [setback-probe] arc {k} midpoint nearest surface: {} at {:.3e}",
                nearest.1, nearest.0
            );
            // And the fillet face's u-extent: distance from arc mid to the fillet face's loop edges.
            let ff = find_fillet_face(brep, &fillet_surfaces[k])?;
            let mut best_edge = f64::INFINITY;
            for he in face_loop_half_edges(brep, ff)? {
                let edge = brep.half_edges[he].edge;
                if let Some(c) = brep.edges[edge].curve.as_ref() {
                    best_edge = best_edge.min(curve_distance(&c.geometry, c.t_start, c.t_end, mid));
                }
            }
            eprintln!(
                "  [setback-probe] arc {k} midpoint distance to its fillet face boundary: {best_edge:.3e}"
            );
            eprintln!("  [setback-probe] fillet {k} face = {ff:?}");
            // Corner territory: nearest surfaced face to the corner point.
            let corner = arcs[k].corner_next;
            let mut ranked: Vec<(f64, String)> = Vec::new();
            for (fh, face) in brep.faces.iter() {
                let Some(surface) = face.surface.as_ref() else {
                    continue;
                };
                let kind = match &surface.geometry {
                    popkern_core::geometry::SurfaceGeometry::Plane { .. } => "Plane",
                    popkern_core::geometry::SurfaceGeometry::Nurbs(_) => "Nurbs",
                    popkern_core::geometry::SurfaceGeometry::Cylinder { .. } => "Cylinder",
                    _ => "Other",
                };
                if let Some((u, v)) = popkern_core::geometry::uv_projection::project_point_to_uv(
                    &surface.geometry,
                    &corner,
                ) {
                    let d = (&surface.geometry as &dyn SurfaceEval)
                        .point_at(u, v)
                        .distance_to(&corner);
                    ranked.push((d, format!("{fh:?} {kind}")));
                }
            }
            ranked.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
            for (d, label) in ranked.iter().take(3) {
                eprintln!("  [setback-probe] corner {k} near: {label} d={d:.3e}");
            }
        }
    }

    // ---- discovery: fillet faces --------------------------------------
    let fillet_faces: Vec<FaceHandle> = (0..n)
        .map(|k| find_fillet_face(brep, &fillet_surfaces[k]))
        .collect::<Result<_, _>>()?;

    // ---- surgery A/B/C: three-segment arc insertion ---------------------
    // Measured reality (see probes): each setback arc's MIDDLE is
    // interior to its fillet face, while its two END segments exit
    // through the fillet face's apex-end boundary and run across the
    // lateral PLANE faces to the corner points (which lie exactly on
    // those planes). Per fillet: split the boundary at the two arc
    // crossings, insert the middle sub-arc across the fillet face
    // (mef); per corner: cut the lateral face with the two adjacent
    // end sub-arcs via a spur (mev to c) + close (mef).
    let mut arc_middle_edges = Vec::with_capacity(n);
    let mut apex_slivers = Vec::with_capacity(n);
    let mut crossing_vertices = vec![[None::<VertexHandle>; 2]; n]; // [X1 at t1, X2 at t2]
    let mut end_segments = vec![[None::<[Point3; 4]>; 2]; n]; // [c_prev→X1, X2→c_next]
    for k in 0..n {
        let crossings = arc_boundary_crossings(brep, fillet_faces[k], &arcs[k].cubic)?;
        let [(t1, e1, p1), (t2, e2, p2)] = crossings;
        let v1 = brep
            .split_edge_at_point(e1, p1, tolerance)
            .map_err(|e| format!("split crossing 1 of fillet {k}: {e}"))?
            .vertex;
        let v2 = brep
            .split_edge_at_point(e2, p2, tolerance)
            .map_err(|e| format!("split crossing 2 of fillet {k}: {e}"))?
            .vertex;
        let (mut seg_pre, mut seg_mid, mut seg_post) = split_cubic_three(&arcs[k].cubic, t1, t2);
        // Snap segment endpoints to the topology vertices (the split
        // point is the edge's own foot, ~4e-4 from the raw arc).
        let p1v = brep.vertices[v1].position;
        let p2v = brep.vertices[v2].position;
        seg_pre[3] = p1v;
        seg_mid[0] = p1v;
        seg_mid[3] = p2v;
        seg_post[0] = p2v;
        let he1 = half_edge_at_vertex_on_face(brep, fillet_faces[k], v1)?;
        let he2 = half_edge_at_vertex_on_face(brep, fillet_faces[k], v2)?;
        let mef = brep
            .mef(he1, he2)
            .map_err(|e| format!("mef arc middle on fillet {k}: {e}"))?;
        if let Some(edge) = brep.edges.get_mut(mef.edge) {
            edge.curve = Some(popkern_core::geometry::CurveOnEdge {
                geometry: popkern_core::geometry::CurveGeometry::Nurbs(bezier_cubic_curve(
                    &seg_mid,
                )),
                t_start: 0.0,
                t_end: 1.0,
            });
        }
        let d_old = face_centroid(brep, fillet_faces[k])?.distance_to(&vertex_position);
        let d_new = face_centroid(brep, mef.face)?.distance_to(&vertex_position);
        let sliver = if d_new < d_old {
            mef.face
        } else {
            fillet_faces[k]
        };
        arc_middle_edges.push(mef.edge);
        apex_slivers.push(sliver);
        crossing_vertices[k] = [Some(v1), Some(v2)];
        end_segments[k] = [Some(seg_pre), Some(seg_post)];
    }

    // Per corner c_k: path X2_k → c_k → X1_{k+1} across the shared
    // lateral face.
    let mut corner_vertices = Vec::with_capacity(n);
    let mut wedge_slivers: Vec<FaceHandle> = Vec::new();
    for k in 0..n {
        let k_next = (k + 1) % n;
        let x2 = crossing_vertices[k][1].ok_or("missing X2")?;
        let x1n = crossing_vertices[k_next][0].ok_or("missing X1'")?;
        let corner = arcs[k].corner_next;
        let seg_a = end_segments[k][1].ok_or("missing post segment")?;
        let seg_b = end_segments[k_next][0].ok_or("missing pre segment")?;
        // The lateral face: common face of X2 and X1' that is not a
        // fillet face.
        let face = common_face_excluding(brep, x2, x1n, &fillet_faces)
            .ok_or_else(|| format!("no shared lateral face at corner {k}"))?;
        // Spur to the corner point.
        let he_x2 = half_edge_at_vertex_on_face(brep, face, x2)?;
        let mev = brep
            .mev(he_x2, corner)
            .map_err(|e| format!("mev corner spur at {k}: {e}"))?;
        if let Some(edge) = brep.edges.get_mut(mev.edge) {
            edge.curve = Some(popkern_core::geometry::CurveOnEdge {
                geometry: popkern_core::geometry::CurveGeometry::Nurbs(bezier_cubic_curve(&seg_a)),
                t_start: 0.0,
                t_end: 1.0,
            });
        }
        // Close across the face to X1'.
        let he_c = half_edge_at_vertex_on_face(brep, face, mev.vertex)?;
        let he_x1 = half_edge_at_vertex_on_face(brep, face, x1n)?;
        let mef = brep
            .mef(he_c, he_x1)
            .map_err(|e| format!("mef corner close at {k}: {e}"))?;
        if let Some(edge) = brep.edges.get_mut(mef.edge) {
            edge.curve = Some(popkern_core::geometry::CurveOnEdge {
                geometry: popkern_core::geometry::CurveGeometry::Nurbs(bezier_cubic_curve(&seg_b)),
                t_start: 0.0,
                t_end: 1.0,
            });
        }
        let d_old = face_centroid(brep, face)?.distance_to(&vertex_position);
        let d_new = face_centroid(brep, mef.face)?.distance_to(&vertex_position);
        wedge_slivers.push(if d_new < d_old { mef.face } else { face });
        corner_vertices.push(mev.vertex);
    }
    let arc_edges = arc_middle_edges;

    // ---- surgery D: merge the surrendered region ------------------------
    let mut surrendered: Vec<FaceHandle> = Vec::new();
    surrendered.extend(apex_slivers.iter().copied());
    surrendered.extend(wedge_slivers.iter().copied());
    let blend_face =
        merge_faces(brep, vertex, &surrendered).map_err(|e| format!("blend region merge: {e}"))?;

    Ok(SetbackTopologyOutcome {
        blend_face,
        arc_edges,
        corner_vertices,
    })
}

/// Find the face carrying this fillet's sweep surface (sample match).
fn find_fillet_face(brep: &BRep, fillet: &FilletSurfaceResult) -> Result<FaceHandle, String> {
    let probes = [(0.5, 0.5), (0.25, 0.75), (0.75, 0.25)];
    let targets: Vec<Point3> = probes
        .iter()
        .map(|&(u, v)| fillet.surface.point_at(u, v))
        .collect();
    for (fh, face) in brep.faces.iter() {
        let Some(surface) = face.surface.as_ref() else {
            continue;
        };
        let matches = probes.iter().zip(&targets).all(|(&(u, v), target)| {
            (&surface.geometry as &dyn SurfaceEval)
                .point_at(u, v)
                .distance_to(target)
                < 1.0e-9
        });
        if matches {
            return Ok(fh);
        }
    }
    Err("no face carries the fillet sweep surface".to_string())
}

/// Merge the surrendered faces (apex slivers, wedges, and everything
/// still touching the apex vertex) into a single blend face by
/// kef-ing shared edges.
fn merge_faces(
    brep: &mut BRep,
    vertex: VertexHandle,
    seed_faces: &[FaceHandle],
) -> Result<FaceHandle, String> {
    use std::collections::BTreeSet;
    let mut region: BTreeSet<FaceHandle> = seed_faces.iter().copied().collect();
    // Everything still incident to the apex vertex is surrendered.
    if brep.vertices.get(vertex).is_some() {
        for (fh, face) in brep.faces.iter() {
            let mut touches = false;
            for he in face_loop_half_edges(brep, fh).unwrap_or_default() {
                if brep.half_edges[he].origin == vertex {
                    touches = true;
                    break;
                }
            }
            let _ = face;
            if touches {
                region.insert(fh);
            }
        }
    }
    if region.is_empty() {
        return Err("no faces to merge".to_string());
    }
    // Iteratively kef edges interior to the region.
    loop {
        let mut merged_any = false;
        let edges: Vec<EdgeHandle> = brep.edges.iter().map(|(eh, _)| eh).collect();
        for eh in edges {
            let Some(edge) = brep.edges.get(eh) else {
                continue;
            };
            let he = edge.half_edge;
            let Some(he_data) = brep.half_edges.get(he) else {
                continue;
            };
            let twin = he_data.twin;
            let f1 = brep.loops[he_data.loop_].face;
            let f2 = brep.loops[brep.half_edges[twin].loop_].face;
            if f1 != f2 && region.contains(&f1) && region.contains(&f2) {
                match brep.kef(he) {
                    Ok(kef) => {
                        region.remove(&kef.removed_face);
                        let survivor = brep.loops[kef.surviving_loop].face;
                        region.insert(survivor);
                        merged_any = true;
                        break;
                    }
                    Err(_) => continue,
                }
            }
        }
        if !merged_any {
            break;
        }
    }
    // Dissolve interior remnants of the old apex structure: edges
    // whose two half-edges sit in the SAME loop (bridges the blend
    // loop bounces through twice). Iterate so chains and the old
    // anchor ring unwind from their exposed ends.
    loop {
        let mut dissolved = false;
        let candidates: Vec<HalfEdgeHandle> = brep
            .half_edges
            .iter()
            .filter(|(_, data)| {
                brep.half_edges
                    .get(data.twin)
                    .map(|t| t.loop_ == data.loop_)
                    .unwrap_or(false)
            })
            .map(|(h, _)| h)
            .collect();
        if std::env::var_os("PK_SETBACK_PROBE").is_some() {
            eprintln!(
                "  [merge] same-loop bridge half-edges: {}",
                candidates.len()
            );
        }
        for he in candidates {
            if brep.half_edges.get(he).is_none() {
                continue;
            }
            match brep.dissolve_same_loop_bridge(he) {
                Ok(()) => dissolved = true,
                Err(e) => {
                    if std::env::var_os("PK_SETBACK_PROBE").is_some() {
                        eprintln!("  [merge] dissolve declined: {e}");
                    }
                }
            }
        }
        if !dissolved {
            break;
        }
    }
    super::blend_insert::kev_v_bridge(brep, vertex).map_err(|e| format!("kev_v_bridge: {e}"))?;
    super::blend_insert::remove_orphaned_vertex(brep, vertex);
    // Dissolved bridges strand their interior vertices; sweep all
    // orphans (no incident half-edges).
    let orphans: Vec<VertexHandle> = brep
        .vertices
        .iter()
        .map(|(vh, _)| vh)
        .filter(|&vh| !brep.half_edges.iter().any(|(_, h)| h.origin == vh))
        .collect();
    for vh in orphans {
        super::blend_insert::remove_orphaned_vertex(brep, vh);
    }

    let mut it = region.iter();
    let blend = *it.next().ok_or("merge produced no face")?;
    if it.next().is_some() {
        return Err(format!(
            "blend region did not merge to one face ({} remain)",
            region.len()
        ));
    }
    Ok(blend)
}

// ---- helpers --------------------------------------------------------------

fn face_loop_half_edges(brep: &BRep, face: FaceHandle) -> Result<Vec<HalfEdgeHandle>, String> {
    let f = brep.faces.get(face).ok_or("face handle invalid")?;
    let start = brep
        .loops
        .get(f.outer_loop)
        .ok_or("outer loop invalid")?
        .half_edge;
    let mut out = Vec::new();
    let mut he = start;
    for _ in 0..10_000 {
        out.push(he);
        he = brep.half_edges[he].next;
        if he == start {
            return Ok(out);
        }
    }
    Err("unterminated face loop".to_string())
}

fn half_edge_at_vertex_on_face(
    brep: &BRep,
    face: FaceHandle,
    vertex: VertexHandle,
) -> Result<HalfEdgeHandle, String> {
    for he in face_loop_half_edges(brep, face)? {
        if brep.half_edges[he].origin == vertex {
            return Ok(he);
        }
    }
    Err("vertex not on face loop".to_string())
}

fn common_face(brep: &BRep, a: VertexHandle, b: VertexHandle) -> Option<FaceHandle> {
    let mut faces_a = Vec::new();
    for (fh, _) in brep.faces.iter() {
        if let Ok(hes) = face_loop_half_edges(brep, fh) {
            let has_a = hes.iter().any(|&he| brep.half_edges[he].origin == a);
            let has_b = hes.iter().any(|&he| brep.half_edges[he].origin == b);
            if has_a && has_b {
                faces_a.push(fh);
            }
            let _ = &has_a;
        }
    }
    faces_a.first().copied()
}

fn face_centroid(brep: &BRep, face: FaceHandle) -> Result<Point3, String> {
    let hes = face_loop_half_edges(brep, face)?;
    let mut acc = [0.0_f64; 3];
    let mut count = 0.0;
    for he in hes {
        let p = brep.vertices[brep.half_edges[he].origin].position;
        acc[0] += p.x;
        acc[1] += p.y;
        acc[2] += p.z;
        count += 1.0;
    }
    Ok(Point3::new(acc[0] / count, acc[1] / count, acc[2] / count))
}

fn curve_distance(
    curve: &popkern_core::geometry::CurveGeometry,
    t_start: f64,
    t_end: f64,
    point: Point3,
) -> f64 {
    // Use the edge's own parameter window: curve.parameter_range()
    // is unbounded for analytic lines (the known gotcha).
    let eval: &dyn CurveEval = curve;
    let (t0, t1) = if t_start.is_finite() && t_end.is_finite() {
        (t_start, t_end)
    } else {
        eval.parameter_range()
    };
    if !(t0.is_finite() && t1.is_finite()) {
        return f64::INFINITY;
    }
    let mut best = f64::INFINITY;
    for i in 0..=64 {
        let t = t0 + (t1 - t0) * i as f64 / 64.0;
        best = best.min(eval.point_at(t).distance_to(&point));
    }
    best
}

/// Find where the arc crosses the face's boundary: for each boundary
/// edge, closest approach between the arc and the edge curve; the two
/// sub-tolerance approaches (sorted by arc parameter) are the
/// crossings.
fn arc_boundary_crossings(
    brep: &BRep,
    face: FaceHandle,
    cubic: &[Point3; 4],
) -> Result<[(f64, EdgeHandle, Point3); 2], String> {
    let arc = bezier_cubic_curve(cubic);
    let arc_eval: &dyn CurveEval = &popkern_core::geometry::CurveGeometry::Nurbs(arc);
    let mut hits: Vec<(f64, EdgeHandle, Point3, f64)> = Vec::new();
    for he in face_loop_half_edges(brep, face)? {
        let edge = brep.half_edges[he].edge;
        let Some(curve) = brep.edges[edge].curve.as_ref() else {
            continue;
        };
        let edge_eval: &dyn CurveEval = &curve.geometry;
        let (e0, e1) = if curve.t_start.is_finite() && curve.t_end.is_finite() {
            (curve.t_start, curve.t_end)
        } else {
            continue;
        };
        let mut best = (f64::INFINITY, 0.0_f64, Point3::new(0.0, 0.0, 0.0));
        for i in 0..=128 {
            let t = i as f64 / 128.0;
            let pa = arc_eval.point_at(t);
            for j in 0..=256 {
                let te = e0 + (e1 - e0) * j as f64 / 256.0;
                let pe = edge_eval.point_at(te);
                let d = pe.distance_to(&pa);
                if d < best.0 {
                    // Split at the EDGE's own point (the splitter
                    // requires exact point inversion on the edge).
                    best = (d, t, pe);
                }
            }
        }
        if best.0 < 1.0e-3 {
            hits.push((best.1, edge, best.2, best.0));
        }
    }
    // Dedup near-identical arc parameters (shared vertices between
    // consecutive edges register on both).
    hits.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    hits.dedup_by(|a, b| (a.0 - b.0).abs() < 1.0e-3);
    if hits.len() != 2 {
        return Err(format!(
            "expected 2 arc/boundary crossings, found {} ({:?})",
            hits.len(),
            hits.iter().map(|h| (h.0, h.3)).collect::<Vec<_>>()
        ));
    }
    Ok([
        (hits[0].0, hits[0].1, hits[0].2),
        (hits[1].0, hits[1].1, hits[1].2),
    ])
}

/// De Casteljau three-way split of a cubic at t1 < t2, returning the
/// (pre, mid, post) sub-cubics.
fn split_cubic_three(
    cps: &[Point3; 4],
    t1: f64,
    t2: f64,
) -> ([Point3; 4], [Point3; 4], [Point3; 4]) {
    fn split(cps: &[Point3; 4], t: f64) -> ([Point3; 4], [Point3; 4]) {
        let lerp = |a: Point3, b: Point3| {
            Point3::new(
                a.x + (b.x - a.x) * t,
                a.y + (b.y - a.y) * t,
                a.z + (b.z - a.z) * t,
            )
        };
        let q0 = lerp(cps[0], cps[1]);
        let q1 = lerp(cps[1], cps[2]);
        let q2 = lerp(cps[2], cps[3]);
        let r0 = lerp(q0, q1);
        let r1 = lerp(q1, q2);
        let sm = lerp(r0, r1);
        ([cps[0], q0, r0, sm], [sm, r1, q2, cps[3]])
    }
    let (pre, rest) = split(cps, t1);
    let t2r = if t1 < 1.0 {
        (t2 - t1) / (1.0 - t1)
    } else {
        0.0
    };
    let (mid, post) = split(&rest, t2r);
    (pre, mid, post)
}

/// Common face of two vertices, excluding the given fillet faces.
fn common_face_excluding(
    brep: &BRep,
    a: VertexHandle,
    b: VertexHandle,
    exclude: &[FaceHandle],
) -> Option<FaceHandle> {
    for (fh, _) in brep.faces.iter() {
        if exclude.contains(&fh) {
            continue;
        }
        if let Ok(hes) = face_loop_half_edges(brep, fh) {
            let has_a = hes.iter().any(|&he| brep.half_edges[he].origin == a);
            let has_b = hes.iter().any(|&he| brep.half_edges[he].origin == b);
            if has_a && has_b {
                return Some(fh);
            }
        }
    }
    None
}

fn eval_cubic_mid(cps: &[Point3; 4]) -> Point3 {
    Point3::new(
        0.125 * (cps[0].x + 3.0 * cps[1].x + 3.0 * cps[2].x + cps[3].x),
        0.125 * (cps[0].y + 3.0 * cps[1].y + 3.0 * cps[2].y + cps[3].y),
        0.125 * (cps[0].z + 3.0 * cps[1].z + 3.0 * cps[2].z + cps[3].z),
    )
}

fn chord_distance(a: Point3, b: Point3, p: Point3) -> f64 {
    let ab = b - a;
    let len2 = ab.dot(&ab);
    if len2 < 1.0e-30 {
        return a.distance_to(&p);
    }
    let t = ((p - a).dot(&ab) / len2).clamp(0.0, 1.0);
    Point3::new(a.x + ab.x * t, a.y + ab.y * t, a.z + ab.z * t).distance_to(&p)
}

fn bezier_cubic_curve(cps: &[Point3; 4]) -> NurbsCurve {
    NurbsCurve {
        degree: 3,
        knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
        control_points: cps.to_vec(),
        weights: vec![1.0; 4],
    }
}

/// Loop length of a face's outer loop (WIP checkpoint reporting).
pub(crate) fn blend_loop_len(brep: &BRep, face: FaceHandle) -> Option<usize> {
    face_loop_half_edges(brep, face).ok().map(|v| v.len())
}

/// WIP checkpoint: describe the blend loop's edges.
pub(crate) fn blend_loop_description(brep: &BRep, face: FaceHandle) -> String {
    let Ok(hes) = face_loop_half_edges(brep, face) else {
        return "?".to_string();
    };
    let mut out = String::new();
    for he in hes {
        let v = brep.half_edges[he].origin;
        let p = brep.vertices[v].position;
        let edge = brep.half_edges[he].edge;
        let has_curve = brep.edges[edge].curve.is_some();
        let vn = brep.half_edges[brep.half_edges[he].next].origin;
        let pn = brep.vertices[vn].position;
        let len = p.distance_to(&pn);
        out.push_str(&format!(
            "\n    ({:.3},{:.3},{:.3}) len={:.4} curve={}",
            p.x, p.y, p.z, len, has_curve
        ));
    }
    out
}

/// Returns whether the WIP arm is enabled (integration tests only).
pub(crate) fn setback_topology_enabled() -> bool {
    SETBACK_TOPOLOGY_WIP.load(Ordering::Relaxed)
}
