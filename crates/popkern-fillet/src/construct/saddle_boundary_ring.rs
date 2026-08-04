//! Boundary-driven blend-ring carve for saddle (`OppositeSide`)
//! corners (PK-234).
//!
//! The mitred corner construction leaves the saddle's blend ring on
//! the SHARP model edges: rim chords claim fillet cylinders they sit
//! up to `√2·r − r` away from, while the three-quad decomposition's
//! Gregory patches interpolate the legacy [`BlendBoundary`] 1.3e-1
//! away (measured, notched-cube saddle). The corner-point setback
//! surgery cannot serve this class either — at the reflex junction
//! between the two convex strips the rims never meet, so its
//! 3-corner ring model has no valid corner (measured: the
//! tangent-line crossing sits outside both strips' physical
//! extents).
//!
//! This carve builds the ring the boundary actually describes:
//!
//! * per `FilletEnd` segment — split the strip face's boundary at
//!   the arc's two junction endpoints (they land on the strip's
//!   tangency rails by construction) and `mef` the exact arc curve
//!   across the strip;
//! * per `FaceEdge` segment — `mef` the bridge curve across its
//!   support face between the two junction vertices the arcs
//!   created;
//! * merge the surrendered region (apex slivers + wedges inside the
//!   ring) into a single `2n`-gon blend face via the setback
//!   surgery's region merge.
//!
//! The blend face then flows through the UNCHANGED legacy three-quad
//! dispatch: the poke sees a `2n` ring, the seam re-anchor measures
//! already-anchored (the ring IS the boundary), and the decomposition
//! bonds onto edges that lie on their faces.
//!
//! Transactional: any step error restores the entry snapshot and the
//! caller falls through to the mitred legacy construction with the
//! decline recorded (today's behavior, no regression).

use super::*;

/// Junction acceptance: a boundary junction must land on a strip
/// boundary edge within this band before we split there. Junctions
/// are constructed ON the tangency rails, and the BRep rail curves
/// are the same carve output the boundary read, so genuine landings
/// measure ~1e-9; a 1e-4 gate rejects a boundary/topology mismatch
/// class without admitting mis-splits.
const JUNCTION_ON_EDGE_BAND: f64 = 1.0e-4;

/// Ring postcondition: every ring edge curve sample must sit on both
/// adjacent surfaces within this band (the arcs are exact chart isos
/// and the bridges live on their support faces, so genuine rings
/// measure ≤ 1e-9; check 9's default band is 1e-5).
const RING_ON_SURFACE_BAND: f64 = 1.0e-5;

/// Carve the boundary ring at a saddle vertex and merge the interior
/// into one `2n`-gon blend face. Returns the blend face; on any
/// decline the BRep is exactly the entry state.
pub(crate) fn try_saddle_boundary_ring(
    brep: &mut BRep,
    vertex: VertexHandle,
    boundary: &crate::blend::boundary::BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
    tolerance: &Tolerance,
) -> Result<FaceHandle, String> {
    use crate::blend::boundary::BoundarySegment;

    let n = fillet_surfaces.len();
    if boundary.segments.len() != 2 * n {
        return Err(format!(
            "boundary has {} segments, expected {}",
            boundary.segments.len(),
            2 * n
        ));
    }
    let vertex_position = brep
        .vertices
        .get(vertex)
        .map(|v| v.position)
        .ok_or_else(|| "saddle vertex not in BRep".to_string())?;

    // Segment pattern: alternating FilletEnd / FaceEdge starting at a
    // FilletEnd (the boundary builder's contract).
    for (i, seg) in boundary.segments.iter().enumerate() {
        let ok = if i % 2 == 0 {
            matches!(seg, BoundarySegment::FilletEnd { .. })
        } else {
            matches!(seg, BoundarySegment::FaceEdge { .. })
        };
        if !ok {
            return Err(format!("segment {i} breaks the arc/bridge alternation"));
        }
    }

    let snapshot = brep.clone();
    match carve_ring(
        brep,
        vertex,
        vertex_position,
        boundary,
        fillet_surfaces,
        tolerance,
    ) {
        Ok(face) => Ok(face),
        Err(reason) => {
            *brep = snapshot;
            Err(reason)
        }
    }
}

fn carve_ring(
    brep: &mut BRep,
    vertex: VertexHandle,
    vertex_position: Point3,
    boundary: &crate::blend::boundary::BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
    tolerance: &Tolerance,
) -> Result<FaceHandle, String> {
    use crate::blend::boundary::BoundarySegment;

    let n = fillet_surfaces.len();
    let trace = std::env::var_os("PK_CARVE_TRACE").is_some();

    // ---- arcs: one exact strip-end curve per fillet -----------------
    let mut arc_vertices: Vec<[VertexHandle; 2]> = Vec::with_capacity(n);
    let mut ring_edges: std::collections::BTreeSet<EdgeHandle> = Default::default();
    let mut slivers: Vec<FaceHandle> = Vec::new();
    for k in 0..n {
        let BoundarySegment::FilletEnd {
            fillet_index,
            curve,
        } = &boundary.segments[2 * k]
        else {
            return Err("alternation check missed a non-arc even segment".into());
        };
        let fillet = fillet_surfaces
            .get(*fillet_index)
            .ok_or_else(|| format!("fillet index {fillet_index} out of range"))?;
        let (t0, t1) = curve.parameter_range();
        let p_start = curve.point_at(t0);
        let p_end = curve.point_at(t1);
        // Bind the strip face by the arc's interior samples (the
        // binding cubic is geometric — interior points only).
        let pseudo_cubic = [
            p_start,
            curve.point_at(t0 + (t1 - t0) / 3.0),
            curve.point_at(t0 + 2.0 * (t1 - t0) / 3.0),
            p_end,
        ];
        let strip_face =
            setback_arc_binding::find_fillet_face(brep, fillet, &pseudo_cubic, vertex_position)?;
        let v_start = junction_vertex(brep, strip_face, p_start, tolerance)
            .map_err(|e| format!("arc {k} start junction: {e}"))?;
        let v_end = junction_vertex(brep, strip_face, p_end, tolerance)
            .map_err(|e| format!("arc {k} end junction: {e}"))?;
        if v_start == v_end {
            return Err(format!("arc {k} junctions collapse to one vertex"));
        }
        let he_start = half_edge_at_vertex_on_face(brep, strip_face, v_start)?;
        let he_end = half_edge_at_vertex_on_face(brep, strip_face, v_end)?;
        let mef = brep
            .mef(he_start, he_end)
            .map_err(|e| format!("mef arc {k} across strip: {e}"))?;
        write_ring_curve(brep, mef.edge, curve, v_start, v_end);
        let d_old = face_centroid(brep, strip_face)?.distance_to(&vertex_position);
        let d_new = face_centroid(brep, mef.face)?.distance_to(&vertex_position);
        slivers.push(if d_new < d_old { mef.face } else { strip_face });
        ring_edges.insert(mef.edge);
        arc_vertices.push([v_start, v_end]);
        if trace {
            eprintln!(
                "  [saddle-ring] arc {k} strip {strip_face:?} \
                 ({:.4},{:.4},{:.4})->({:.4},{:.4},{:.4})",
                p_start.x, p_start.y, p_start.z, p_end.x, p_end.y, p_end.z
            );
        }
    }

    // ---- bridges: one support-face curve per junction ---------------
    // Bridge 2k+1 connects arc k's END to arc (k+1)'s START (the
    // boundary builder emits the ring in cyclic order).
    for k in 0..n {
        let BoundarySegment::FaceEdge { face, curve } = &boundary.segments[2 * k + 1] else {
            return Err("alternation check missed a non-bridge odd segment".into());
        };
        let v_a = arc_vertices[k][1];
        let v_b = arc_vertices[(k + 1) % n][0];
        if v_a == v_b {
            // Degenerate bridge: the ring closes through the shared
            // junction. The blend face then has fewer than 2n edges
            // and the poke contract breaks — decline until a fixture
            // exercises this class.
            return Err(format!("bridge {k} is degenerate (shared junction vertex)"));
        }
        // The support face recorded by the boundary must still carry
        // both junction vertices on its loop (the arc mefs only split
        // rail edges, which bound this face).
        let support = *face;
        let (t0, t1) = curve.parameter_range();
        let b_start = curve.point_at(t0);
        let b_end = curve.point_at(t1);
        let pa = brep.vertices[v_a].position;
        // Bridge orientation in the boundary may run either way.
        let (v_from, v_to) = if b_start.distance_to(&pa) <= b_end.distance_to(&pa) {
            (v_a, v_b)
        } else {
            (v_b, v_a)
        };
        let d_from = brep.vertices[v_from].position.distance_to(&b_start);
        let d_to = brep.vertices[v_to].position.distance_to(&b_end);
        if d_from > JUNCTION_ON_EDGE_BAND || d_to > JUNCTION_ON_EDGE_BAND {
            return Err(format!(
                "bridge {k} endpoints miss the junction vertices \
                 (start {d_from:.3e}, end {d_to:.3e})"
            ));
        }
        let he_from = half_edge_at_vertex_on_face(brep, support, v_from)
            .map_err(|e| format!("bridge {k} start not on support {support:?}: {e}"))?;
        let he_to = half_edge_at_vertex_on_face(brep, support, v_to)
            .map_err(|e| format!("bridge {k} end not on support {support:?}: {e}"))?;
        let mef = brep
            .mef(he_from, he_to)
            .map_err(|e| format!("mef bridge {k} across support: {e}"))?;
        write_ring_curve(brep, mef.edge, curve, v_from, v_to);
        let d_old = face_centroid(brep, support)?.distance_to(&vertex_position);
        let d_new = face_centroid(brep, mef.face)?.distance_to(&vertex_position);
        slivers.push(if d_new < d_old { mef.face } else { support });
        ring_edges.insert(mef.edge);
        if trace {
            eprintln!(
                "  [saddle-ring] bridge {k} support {support:?} \
                 ({:.4},{:.4},{:.4})->({:.4},{:.4},{:.4})",
                b_start.x, b_start.y, b_start.z, b_end.x, b_end.y, b_end.z
            );
        }
    }

    // ---- merge the surrendered interior ------------------------------
    let blend_face = merge_faces(brep, vertex, &slivers, &ring_edges)
        .map_err(|e| format!("saddle ring region merge: {e}"))?;

    // ---- postconditions ----------------------------------------------
    let ring_len = face_loop_half_edges(brep, blend_face)?.len();
    if ring_len != 2 * n {
        return Err(format!(
            "blend face ring has {ring_len} edges, expected {}",
            2 * n
        ));
    }
    let worst = worst_ring_on_surface_deviation(brep, &ring_edges)?;
    if worst > RING_ON_SURFACE_BAND {
        return Err(format!(
            "ring edge sits {worst:.3e} off an adjacent surface (band {RING_ON_SURFACE_BAND:.0e})"
        ));
    }
    if trace {
        eprintln!(
            "  [saddle-ring] blend face {blend_face:?} ring {ring_len} edges, \
             worst on-surface deviation {worst:.3e}"
        );
    }
    let _ = tolerance;
    Ok(blend_face)
}

/// Find or create the junction vertex where `point` lands on `face`'s
/// outer-loop boundary: reuse a loop vertex within the junction band,
/// else split the containing edge at the point.
fn junction_vertex(
    brep: &mut BRep,
    face: FaceHandle,
    point: Point3,
    tolerance: &Tolerance,
) -> Result<VertexHandle, String> {
    // Reuse: nearest existing loop vertex.
    let mut nearest_vertex: Option<(f64, VertexHandle)> = None;
    for he in face_loop_half_edges(brep, face)? {
        let v = brep.half_edges[he].origin;
        if let Some(data) = brep.vertices.get(v) {
            let d = data.position.distance_to(&point);
            if nearest_vertex.map(|(bd, _)| d < bd).unwrap_or(true) {
                nearest_vertex = Some((d, v));
            }
        }
    }
    if let Some((d, v)) = nearest_vertex
        && d <= JUNCTION_ON_EDGE_BAND
    {
        return Ok(v);
    }
    // Split: the loop edge whose curve passes nearest the point.
    // Coarse sample, then windowed parametric refinement — the
    // coarse minimum alone overstates the distance by up to half a
    // sample step (measured 1.56e-3 on a 0.3-long rail whose true
    // point-to-line distance is ~0).
    let mut best: Option<(f64, f64, EdgeHandle)> = None;
    for he in face_loop_half_edges(brep, face)? {
        let edge = brep.half_edges[he].edge;
        let Some(curve) = brep.edges[edge].curve.as_ref() else {
            continue;
        };
        let (e0, e1) = (curve.t_start, curve.t_end);
        if !(e0.is_finite() && e1.is_finite()) {
            continue;
        }
        let eval: &dyn CurveEval = &curve.geometry;
        let mut d_min = f64::INFINITY;
        let mut t_min = e0;
        for i in 0..=64 {
            let t = e0 + (e1 - e0) * i as f64 / 64.0;
            let d = eval.point_at(t).distance_to(&point);
            if d < d_min {
                d_min = d;
                t_min = t;
            }
        }
        let mut window = (e1 - e0).abs() / 64.0;
        for _ in 0..24 {
            let (lo, hi) = (
                (t_min - window).max(e0.min(e1)),
                (t_min + window).min(e0.max(e1)),
            );
            for i in 0..=16 {
                let t = lo + (hi - lo) * i as f64 / 16.0;
                let d = eval.point_at(t).distance_to(&point);
                if d < d_min {
                    d_min = d;
                    t_min = t;
                }
            }
            window *= 0.4;
        }
        if best.map(|(bd, _, _)| d_min < bd).unwrap_or(true) {
            best = Some((d_min, t_min, edge));
        }
    }
    let (d, t_foot, edge) = best.ok_or_else(|| "no curved boundary edge to split".to_string())?;
    if d > JUNCTION_ON_EDGE_BAND {
        if std::env::var_os("PK_CARVE_TRACE").is_some() {
            for he in face_loop_half_edges(brep, face)? {
                let hd = &brep.half_edges[he];
                let p = brep.vertices[hd.origin].position;
                let other = brep.loops[brep.half_edges[hd.twin].loop_].face;
                let mid = brep.edges[hd.edge].curve.as_ref().map(|c| {
                    let eval: &dyn CurveEval = &c.geometry;
                    eval.point_at(0.5 * (c.t_start + c.t_end))
                });
                eprintln!(
                    "    [saddle-ring] {face:?} loop v ({:.4},{:.4},{:.4}) e {:?} mid {:?} | other {other:?}",
                    p.x,
                    p.y,
                    p.z,
                    hd.edge,
                    mid.map(|m| (
                        (m.x * 1e4).round() / 1e4,
                        (m.y * 1e4).round() / 1e4,
                        (m.z * 1e4).round() / 1e4
                    ))
                );
            }
        }
        return Err(format!(
            "junction ({:.4},{:.4},{:.4}) is {d:.3e} off face {face:?}'s boundary \
             (band {JUNCTION_ON_EDGE_BAND:.0e})",
            point.x, point.y, point.z
        ));
    }
    // Split at the FOOT on the edge curve (the junction target may
    // sit a hair off the carved edge; the vertex must be on it).
    let foot = {
        let curve = brep.edges[edge]
            .curve
            .as_ref()
            .ok_or_else(|| "split edge lost its curve".to_string())?;
        let eval: &dyn CurveEval = &curve.geometry;
        eval.point_at(t_foot)
    };
    brep.split_edge_at_point(edge, foot, tolerance)
        .map(|r| r.vertex)
        .map_err(|e| format!("split at junction: {e}"))
}

/// Write a boundary segment curve onto a freshly-mef'd ring edge,
/// oriented `v_from → v_to` (the mef argument order) with endpoint
/// control points snapped to the topology vertices.
fn write_ring_curve(
    brep: &mut BRep,
    edge: EdgeHandle,
    curve: &NurbsCurve,
    v_from: VertexHandle,
    v_to: VertexHandle,
) {
    let (t0, t1) = curve.parameter_range();
    let p_from = brep.vertices[v_from].position;
    let mut oriented = curve.clone();
    if oriented.point_at(t1).distance_to(&p_from) < oriented.point_at(t0).distance_to(&p_from) {
        oriented = reverse_nurbs(&oriented);
    }
    if let (Some(first), Some(_)) = (
        oriented.control_points.first().copied(),
        oriented.control_points.last().copied(),
    ) {
        let _ = first;
        let len = oriented.control_points.len();
        oriented.control_points[0] = p_from;
        oriented.control_points[len - 1] = brep.vertices[v_to].position;
    }
    let (t0, t1) = oriented.parameter_range();
    if let Some(edge_data) = brep.edges.get_mut(edge) {
        edge_data.curve = Some(popkern_core::geometry::CurveOnEdge {
            geometry: popkern_core::geometry::CurveGeometry::Nurbs(oriented),
            t_start: t0,
            t_end: t1,
        });
    }
}

/// Reverse a clamped NURBS curve (control points, weights, mirrored
/// knots) so `point_at(t0)` and `point_at(t1)` swap.
fn reverse_nurbs(curve: &NurbsCurve) -> NurbsCurve {
    let (lo, hi) = curve.parameter_range();
    let mut control_points = curve.control_points.clone();
    control_points.reverse();
    let mut weights = curve.weights.clone();
    weights.reverse();
    let mut knots: Vec<f64> = curve.knots.iter().map(|k| lo + hi - k).collect();
    knots.reverse();
    NurbsCurve {
        degree: curve.degree,
        control_points,
        weights,
        knots,
    }
}

/// Worst sampled distance from any ring edge curve to either adjacent
/// face surface (9 samples per edge, cold projection — the adjacent
/// surfaces are analytic strips and support planes).
fn worst_ring_on_surface_deviation(
    brep: &BRep,
    ring_edges: &std::collections::BTreeSet<EdgeHandle>,
) -> Result<f64, String> {
    let mut worst = 0.0_f64;
    for &eh in ring_edges {
        let edge = brep
            .edges
            .get(eh)
            .ok_or_else(|| format!("ring edge {eh:?} vanished during merge"))?;
        let curve = edge
            .curve
            .as_ref()
            .ok_or_else(|| format!("ring edge {eh:?} lost its curve during merge"))?;
        let he = edge.half_edge;
        let he_data = brep
            .half_edges
            .get(he)
            .ok_or_else(|| format!("ring edge {eh:?} has no half-edge"))?;
        for h in [he, he_data.twin] {
            let Some(hd) = brep.half_edges.get(h) else {
                continue;
            };
            let face = brep.loops[hd.loop_].face;
            let Some(surface) = brep.faces[face].surface.as_ref() else {
                continue;
            };
            let eval: &dyn CurveEval = &curve.geometry;
            for i in 0..=8 {
                let t = curve.t_start + (curve.t_end - curve.t_start) * i as f64 / 8.0;
                let p = eval.point_at(t);
                let d = popkern_core::geometry::uv_projection::project_point_to_uv(
                    &surface.geometry,
                    &p,
                )
                .map(|(u, v)| {
                    (&surface.geometry as &dyn SurfaceEval)
                        .point_at(u, v)
                        .distance_to(&p)
                })
                .unwrap_or(f64::INFINITY);
                worst = worst.max(d);
            }
        }
    }
    Ok(worst)
}
