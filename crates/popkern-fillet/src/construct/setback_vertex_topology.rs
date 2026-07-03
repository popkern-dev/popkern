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

    // ---- discovery: fillet faces --------------------------------------
    let fillet_faces: Vec<FaceHandle> = (0..n)
        .map(|k| find_fillet_face(brep, &fillet_surfaces[k]))
        .collect::<Result<_, _>>()?;

    // ---- surgery A: split trim edges at the corners --------------------
    // Corner c_k joins fillet k (its corner_next) and fillet k+1 (its
    // corner_prev). Each corner lies on two trim edges — one on each
    // adjacent fillet face's loop. Splitting is idempotent per edge.
    // Returns per (fillet, side) the corner vertex handle.
    let mut corner_vertex_at = vec![[None::<VertexHandle>; 2]; n]; // [prev, next]
    for k in 0..n {
        let corner = arcs[k].corner_next;
        let k_next = (k + 1) % n;
        let v_this = split_face_boundary_at(brep, fillet_faces[k], corner, tolerance)
            .map_err(|e| format!("split fillet {k} trim at corner {k}: {e}"))?;
        let v_next = split_face_boundary_at(brep, fillet_faces[k_next], corner, tolerance)
            .map_err(|e| format!("split fillet {k_next} trim at corner {k}: {e}"))?;
        corner_vertex_at[k][1] = Some(v_this);
        corner_vertex_at[k_next][0] = Some(v_next);
    }

    // ---- surgery B: arc edge across each fillet face -------------------
    let mut arc_edges = Vec::with_capacity(n);
    let mut apex_slivers = Vec::with_capacity(n);
    for k in 0..n {
        let v_prev = corner_vertex_at[k][0].ok_or("missing prev corner vertex")?;
        let v_next = corner_vertex_at[k][1].ok_or("missing next corner vertex")?;
        let (arc_edge, sliver) = split_face_with_arc(
            brep,
            fillet_faces[k],
            v_prev,
            v_next,
            &arcs[k].cubic,
            vertex_position,
        )
        .map_err(|e| format!("arc split of fillet face {k}: {e}"))?;
        arc_edges.push(arc_edge);
        apex_slivers.push(sliver);
    }

    // ---- surgery C: pinch lateral faces at the corners ------------------
    // The two corner vertices at c_k (one per adjacent trim edge) both
    // lie on the shared lateral face's loop. A chord between them
    // splits off the face's wedge-above-c; contracting the chord
    // (kev) merges the two vertices into the single corner vertex.
    let mut corner_vertices = Vec::with_capacity(n);
    let mut wedge_slivers = Vec::with_capacity(n);
    for k in 0..n {
        let w1 = corner_vertex_at[k][1].ok_or("missing corner vertex")?;
        let w2 = corner_vertex_at[(k + 1) % n][0].ok_or("missing corner vertex")?;
        let (corner_vertex, wedge) = pinch_lateral_corner(brep, w1, w2, vertex_position)
            .map_err(|e| format!("lateral pinch at corner {k}: {e}"))?;
        corner_vertices.push(corner_vertex);
        wedge_slivers.push(wedge);
    }

    // ---- surgery D: merge the surrendered region ------------------------
    let mut surrendered: Vec<FaceHandle> = Vec::new();
    surrendered.extend(apex_slivers.iter().copied());
    surrendered.extend(wedge_slivers.iter().flatten().copied());
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

/// Split the boundary edge of `face` nearest to `point` at that
/// point. Returns the new (or coincident existing) vertex.
fn split_face_boundary_at(
    brep: &mut BRep,
    face: FaceHandle,
    point: Point3,
    tolerance: &Tolerance,
) -> Result<VertexHandle, String> {
    // Existing vertex within tolerance?
    for he in face_loop_half_edges(brep, face)? {
        let v = brep.half_edges[he].origin;
        if brep.vertices[v].position.distance_to(&point) < 1.0e-9 {
            return Ok(v);
        }
    }
    // Nearest boundary edge: by curve when present, else by the
    // chord between its endpoint vertices.
    let mut best: Option<(f64, EdgeHandle, bool)> = None;
    let mut curved = 0_usize;
    let mut total = 0_usize;
    for he in face_loop_half_edges(brep, face)? {
        let edge = brep.half_edges[he].edge;
        total += 1;
        let d = match brep.edges[edge].curve.as_ref() {
            Some(curve) => {
                curved += 1;
                curve_distance(&curve.geometry, curve.t_start, curve.t_end, point)
            }
            None => {
                let a = brep.vertices[brep.half_edges[he].origin].position;
                let b = brep.vertices[brep.half_edges[brep.half_edges[he].next].origin].position;
                chord_distance(a, b, point)
            }
        };
        let has_curve = brep.edges[edge].curve.is_some();
        if best.map(|(bd, _, _)| d < bd).unwrap_or(true) {
            best = Some((d, edge, has_curve));
        }
    }
    let (d, edge, has_curve) = best.ok_or("face loop has no edges")?;
    if d > 5.0e-3 {
        return Err(format!(
            "corner point is {d:.3e} from the nearest boundary edge \
             (loop edges: {total}, with curves: {curved})"
        ));
    }
    if !has_curve {
        return Err(format!(
            "nearest boundary edge at {d:.3e} has no curve geometry \
             (loop edges: {total}, with curves: {curved})"
        ));
    }
    let result = brep
        .split_edge_at_point(edge, point, tolerance)
        .map_err(|e| format!("split_edge_at_point: {e}"))?;
    Ok(result.vertex)
}

/// Insert the setback arc as an edge across `face` between the two
/// corner vertices, splitting it into the kept fillet face and the
/// apex-side sliver. Returns (arc edge, sliver face).
fn split_face_with_arc(
    brep: &mut BRep,
    face: FaceHandle,
    v_prev: VertexHandle,
    v_next: VertexHandle,
    cubic: &[Point3; 4],
    apex: Point3,
) -> Result<(EdgeHandle, FaceHandle), String> {
    let he1 = half_edge_at_vertex_on_face(brep, face, v_prev)?;
    let he2 = half_edge_at_vertex_on_face(brep, face, v_next)?;
    let result = brep
        .mef(he1, he2)
        .map_err(|e| format!("mef arc edge: {e}"))?;
    // Assign the arc curve (oriented v_prev → v_next).
    let curve = bezier_cubic_curve(cubic);
    if let Some(edge) = brep.edges.get_mut(result.edge) {
        edge.curve = Some(popkern_core::geometry::CurveOnEdge {
            geometry: popkern_core::geometry::CurveGeometry::Nurbs(curve),
            t_start: 0.0,
            t_end: 1.0,
        });
    }
    // Which side is the apex sliver? The face whose loop centroid is
    // nearer the apex.
    let old_face = face;
    let new_face = result.face;
    let d_old = face_centroid(brep, old_face)?.distance_to(&apex);
    let d_new = face_centroid(brep, new_face)?.distance_to(&apex);
    let (sliver, _kept) = if d_new < d_old {
        (new_face, old_face)
    } else {
        (old_face, new_face)
    };
    Ok((result.edge, sliver))
}

/// Split the shared lateral face with a chord w1→w2 (both at the same
/// corner position), then contract the chord so the corner becomes a
/// single vertex. Returns (corner vertex, wedge sliver face if one
/// was split off).
fn pinch_lateral_corner(
    brep: &mut BRep,
    w1: VertexHandle,
    w2: VertexHandle,
    apex: Point3,
) -> Result<(VertexHandle, Option<FaceHandle>), String> {
    if w1 == w2 {
        return Ok((w1, None));
    }
    // The face whose loop contains BOTH vertices (the lateral face).
    let face = common_face(brep, w1, w2).ok_or("no common face for corner pinch")?;
    let he1 = half_edge_at_vertex_on_face(brep, face, w1)?;
    let he2 = half_edge_at_vertex_on_face(brep, face, w2)?;
    let result = brep
        .mef(he1, he2)
        .map_err(|e| format!("mef corner chord: {e}"))?;
    let d_old = face_centroid(brep, face)?.distance_to(&apex);
    let d_new = face_centroid(brep, result.face)?.distance_to(&apex);
    let wedge = if d_new < d_old { result.face } else { face };
    // Contract the chord: kev kills the chord edge and one endpoint.
    let chord_he = brep.edges[result.edge].half_edge;
    brep.kev(chord_he)
        .map_err(|e| format!("kev corner chord: {e}"))?;
    // The surviving vertex is whichever of w1/w2 still exists.
    let corner = if brep.vertices.get(w1).is_some() {
        w1
    } else {
        w2
    };
    Ok((corner, Some(wedge)))
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
    // Remove any leftover interior bridge edges (same face both
    // sides) via kev from their dangling ends, then the apex vertex.
    super::blend_insert::kev_v_bridge(brep, vertex).map_err(|e| format!("kev_v_bridge: {e}"))?;
    super::blend_insert::remove_orphaned_vertex(brep, vertex);

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

/// Returns whether the WIP arm is enabled (integration tests only).
pub(crate) fn setback_topology_enabled() -> bool {
    SETBACK_TOPOLOGY_WIP.load(Ordering::Relaxed)
}
