//! Exact-curve authoring and measurement helpers for the torus
//! strip end-section caps (split from `torus_end_section_caps.rs`
//! for the 500-line structural cap; see that module's docs for the
//! class and the repair contract).

use super::*;
use popkern_math::PortableMath;

/// Re-anchor every edge at `vertex` whose curve endpoints no longer
/// match its vertex positions (the kve repair contract). Straight
/// curves are rebuilt exactly; a stale non-line curve declines.
pub(super) fn rebuild_stale_lines_at(brep: &mut BRep, vertex: VertexHandle) -> Result<(), String> {
    let stale: Vec<EdgeHandle> = brep
        .edges
        .iter()
        .filter_map(|(eh, e)| {
            let hd = brep.half_edges.get(e.half_edge)?;
            let td = brep.half_edges.get(hd.twin)?;
            if hd.origin != vertex && td.origin != vertex {
                return None;
            }
            let c = e.curve.as_ref()?;
            let p0 = brep.vertices.get(hd.origin)?.position;
            let p1 = brep.vertices.get(td.origin)?.position;
            let c0 = c.geometry.point_at(c.t_start);
            let c1 = c.geometry.point_at(c.t_end);
            let hit = |p: Point3| c0.distance_to(&p).min(c1.distance_to(&p)) <= 1.0e-9;
            if hit(p0) && hit(p1) {
                return None;
            }
            Some(eh)
        })
        .collect();
    for eh in stale {
        let is_line = matches!(
            brep.edges[eh].curve.as_ref().map(|c| &c.geometry),
            Some(CurveGeometry::Line { .. })
        );
        if !is_line {
            return Err(format!(
                "stale non-line curve on {eh:?} at the merged corner"
            ));
        }
        rebuild_line(brep, eh)?;
    }
    Ok(())
}

/// Rebuild an edge's curve as the exact segment between its current
/// vertex positions.
pub(super) fn rebuild_line(brep: &mut BRep, edge: EdgeHandle) -> Result<(), String> {
    let he = brep
        .edges
        .get(edge)
        .map(|e| e.half_edge)
        .ok_or("edge vanished")?;
    let hd = brep.half_edges.get(he).ok_or("half-edge vanished")?;
    let p0 = brep.vertices[hd.origin].position;
    let p1 = brep.vertices[brep.half_edges[hd.twin].origin].position;
    brep.edges[edge].curve = CurveOnEdge::line_segment(p0, p1, 1.0e-12);
    Ok(())
}

/// Split `edge` at the closest point to `target` (measured along its
/// chord from `anchor`'s side), land the new vertex exactly on
/// `target`, and return the new vertex.
pub(super) fn split_at(
    brep: &mut BRep,
    edge: EdgeHandle,
    target: Point3,
    anchor: VertexHandle,
) -> Result<VertexHandle, String> {
    let he = brep
        .edges
        .get(edge)
        .map(|e| e.half_edge)
        .ok_or("edge vanished")?;
    let hd = brep.half_edges.get(he).ok_or("half-edge vanished")?;
    let origin = hd.origin;
    let p0 = brep.vertices[origin].position;
    let p1 = brep.vertices[brep.half_edges[hd.twin].origin].position;
    let along = p1 - p0;
    let len2 = along.length_squared();
    if len2 < 1.0e-24 {
        return Err("degenerate chord".into());
    }
    let t = ((target - p0).dot(&along) / len2).clamp(1.0e-6, 1.0 - 1.0e-6);
    let split = brep
        .split_edge(edge, t)
        .map_err(|e| format!("chord split failed: {e}"))?;
    if let Some(v) = brep.vertices.get_mut(split.vertex) {
        v.position = target;
    }
    let _ = anchor;
    // Both pieces get exact segments on their landed endpoints.
    rebuild_line(brep, split.edge_a)?;
    rebuild_line(brep, split.edge_b)?;
    Ok(split.vertex)
}

/// The edge whose two endpoints are exactly `v1`/`v2`, if any.
pub(super) fn edge_between(brep: &BRep, v1: VertexHandle, v2: VertexHandle) -> Option<EdgeHandle> {
    for (eh, e) in &brep.edges {
        let Some(hd) = brep.half_edges.get(e.half_edge) else {
            continue;
        };
        let Some(td) = brep.half_edges.get(hd.twin) else {
            continue;
        };
        if (hd.origin == v1 && td.origin == v2) || (hd.origin == v2 && td.origin == v1) {
            return Some(eh);
        }
    }
    None
}

/// Author `edge` as the exact circle arc `center/normal/…` between
/// its endpoints, sweeping through the side nearest `mid_hint`.
pub(super) fn author_circle_arc(
    brep: &mut BRep,
    edge: EdgeHandle,
    center: Point3,
    normal: Vec3,
    p_from: Point3,
    p_to: Point3,
    mid_hint: Option<Point3>,
) -> Result<(), String> {
    let u_basis = p_from - center;
    let radius = u_basis.length();
    if radius < 1.0e-12 {
        return Err("degenerate arc radius".into());
    }
    let mut v_basis = normal.cross(&u_basis);
    let d = p_to - center;
    let mut sweep = d.dot(&v_basis).patan2(d.dot(&u_basis));
    if sweep < 0.0 {
        v_basis = -v_basis;
        sweep = -sweep;
    }
    // Two arcs join the endpoints; pick the one through the hint.
    if let Some(hint) = mid_hint {
        let arc_mid =
            |u: Vec3, v: Vec3, s: f64| center + u * (s * 0.5).pcos() + v * (s * 0.5).psin();
        let mid = arc_mid(u_basis, v_basis, sweep);
        let other_sweep = std::f64::consts::TAU - sweep;
        let other_mid = arc_mid(u_basis, v_basis * -1.0, other_sweep);
        if hint.distance_to(&other_mid) < hint.distance_to(&mid) {
            v_basis = -v_basis;
            sweep = other_sweep;
        }
    }
    let arc = CurveOnEdge {
        geometry: CurveGeometry::Ellipse {
            center,
            u_basis,
            v_basis,
            start_angle: 0.0,
            end_angle: sweep,
        },
        t_start: 0.0,
        t_end: sweep,
    };
    let hit_from = arc.geometry.point_at(arc.t_start).distance_to(&p_from);
    let hit_to = arc.geometry.point_at(arc.t_end).distance_to(&p_to);
    if hit_from > EXACT_BAND || hit_to > EXACT_BAND {
        return Err(format!(
            "authored arc misses its endpoints (from {hit_from:.3e}, to {hit_to:.3e})"
        ));
    }
    brep.edges.get_mut(edge).ok_or("arc edge vanished")?.curve = Some(arc);
    Ok(())
}

/// The single-edge on-torus trim joining the two run-flank vertices
/// directly (the far contact trim), if it exists as one edge.
pub(super) fn on_torus_trim_edge(
    brep: &BRep,
    _strip: FaceHandle,
    entries: &[(HalfEdgeHandle, EdgeHandle, VertexHandle, Point3, bool)],
    idx_prev: usize,
    idx_after: usize,
) -> Option<EdgeHandle> {
    edge_between(brep, entries[idx_prev].2, entries[idx_after].2)
}

/// Midpoint sample of an edge's current curve (arc-side hint).
pub(super) fn far_trim_mid(brep: &BRep, edge: EdgeHandle) -> Option<Point3> {
    let c = brep.edges.get(edge)?.curve.as_ref()?;
    Some(c.geometry.point_at((c.t_start + c.t_end) * 0.5))
}

/// Whether the shell's flag convention is outward: witnessed by the
/// support face across the strip boundary nearest `probe` — its
/// flag-folded natural normal is compared with the known
/// outward-of-material direction `outward` at that point.
pub(super) fn shell_outward_witness(
    brep: &BRep,
    strip: FaceHandle,
    probe: &Point3,
    outward: &Vec3,
) -> bool {
    let mut best: Option<(f64, bool)> = None;
    // Walk the strip's own loop and test twin faces.
    let Some(f) = brep.faces.get(strip) else {
        return true;
    };
    let Some(l) = brep.loops.get(f.outer_loop) else {
        return true;
    };
    let start = l.half_edge;
    let mut cur = start;
    for _ in 0..brep.half_edges.len() + 1 {
        let Some(hd) = brep.half_edges.get(cur) else {
            break;
        };
        let Some(td) = brep.half_edges.get(hd.twin) else {
            break;
        };
        if let Some(other) = brep.loops.get(td.loop_).map(|ol| ol.face)
            && let Some(sof) = brep.faces.get(other).and_then(|fd| fd.surface.as_ref())
            && let Some((u, v)) = (&sof.geometry as &dyn SurfaceEval).closest_point(probe)
        {
            let on = (&sof.geometry as &dyn SurfaceEval)
                .point_at(u, v)
                .distance_to(probe);
            if let Some(nat) = (&sof.geometry as &dyn SurfaceEval)
                .normal_at(u, v)
                .normalized()
            {
                let folded = nat * if sof.orientation { 1.0 } else { -1.0 };
                let agree = folded.dot(outward) > 0.0;
                if best.is_none_or(|(d, _)| on < d) {
                    best = Some((on, agree));
                }
            }
        }
        cur = hd.next;
        if cur == start {
            break;
        }
    }
    best.map(|(_, agree)| agree).unwrap_or(true)
}

/// One end's mef + section arc + cap plane.
#[allow(clippy::too_many_arguments)]
pub(super) fn author_section_cap(
    brep: &mut BRep,
    strip: FaceHandle,
    corner: VertexHandle,
    spine: Point3,
    span_dir: Vec3,
    vt_pos: Point3,
    vb_pos: Point3,
    minor: f64,
    shell_outward: bool,
) -> Result<(), String> {
    // Locate the corner's half-edge on the strip loop; its loop
    // neighbours must be the two exact end points.
    let loop_h = brep
        .faces
        .get(strip)
        .map(|f| f.outer_loop)
        .ok_or("strip vanished")?;
    let start = brep
        .loops
        .get(loop_h)
        .map(|l| l.half_edge)
        .ok_or("strip loop vanished")?;
    let mut he_v = None;
    let mut cur = start;
    for _ in 0..brep.half_edges.len() + 1 {
        let hd = brep.half_edges.get(cur).ok_or("loop walk broke")?;
        if hd.origin == corner {
            he_v = Some(cur);
            break;
        }
        cur = hd.next;
        if cur == start {
            break;
        }
    }
    let he_v = he_v.ok_or("corner not on the strip loop")?;
    let he1 = brep.half_edges[he_v].prev;
    let he2 = brep.half_edges[he_v].next;
    let o1 = brep.half_edges[he1].origin;
    let o2 = brep.half_edges[he2].origin;
    let (q1, q2) = (brep.vertices[o1].position, brep.vertices[o2].position);
    let near = |p: Point3, q: Point3| p.distance_to(&q) < 1.0e-7;
    if !((near(q1, vt_pos) && near(q2, vb_pos)) || (near(q1, vb_pos) && near(q2, vt_pos))) {
        return Err("corner loop neighbours are not the section-arc endpoints".into());
    }
    let mef = brep
        .mef(he1, he2)
        .map_err(|e| format!("section-cap mef failed: {e}"))?;
    author_circle_arc(brep, mef.edge, spine, span_dir, q1, q2, None)?;
    // Verify the minor arc (the fillet profile spans < π).
    if let Some(c) = brep.edges.get(mef.edge).and_then(|e| e.curve.as_ref())
        && (c.t_end - c.t_start) > std::f64::consts::PI - 1.0e-9
    {
        return Err("section arc swept the major side".into());
    }
    let _ = minor;
    // Cap plane: outward-of-material points INTO the removed span
    // (the stopped-fillet end wall); the stored flag then encodes
    // the shell's own convention, witnessed from the neighbours.
    let corner_pos = brep.vertices[corner].position;
    let u_axis = (vb_pos - spine)
        .normalized()
        .ok_or("degenerate cap frame")?;
    let v_axis = span_dir.cross(&u_axis);
    if let Some(f) = brep.faces.get_mut(mef.face) {
        f.surface = Some(SurfaceOnFace {
            geometry: SurfaceGeometry::Plane {
                origin: corner_pos,
                normal: span_dir,
                u_axis,
                v_axis,
            },
            orientation: shell_outward,
        });
    }
    let _ = brep.recompute_face_pcurves(mef.face);
    Ok(())
}

/// Boundary edges of `strip` and (one level out) of its caps.
pub(super) fn touched_edges(brep: &BRep, strip: FaceHandle) -> Vec<EdgeHandle> {
    let mut out = Vec::new();
    let Some(f) = brep.faces.get(strip) else {
        return out;
    };
    let Some(l) = brep.loops.get(f.outer_loop) else {
        return out;
    };
    let start = l.half_edge;
    let mut cur = start;
    for _ in 0..brep.half_edges.len() + 1 {
        let Some(hd) = brep.half_edges.get(cur) else {
            break;
        };
        out.push(hd.edge);
        cur = hd.next;
        if cur == start {
            break;
        }
    }
    out
}

/// Max deviation of an edge's curve from both adjacent face
/// surfaces, sampled at 9 interior stations.
pub(super) fn edge_max_deviation(brep: &BRep, edge: EdgeHandle) -> f64 {
    let Some(e) = brep.edges.get(edge) else {
        return 0.0;
    };
    let Some(c) = e.curve.as_ref() else {
        return 0.0;
    };
    let Some(hd) = brep.half_edges.get(e.half_edge) else {
        return 0.0;
    };
    let Some(td) = brep.half_edges.get(hd.twin) else {
        return 0.0;
    };
    let faces = [
        brep.loops.get(hd.loop_).map(|l| l.face),
        brep.loops.get(td.loop_).map(|l| l.face),
    ];
    let mut worst = 0.0_f64;
    for face in faces.into_iter().flatten() {
        let Some(sof) = brep.faces.get(face).and_then(|f| f.surface.as_ref()) else {
            continue;
        };
        let eval = &sof.geometry as &dyn SurfaceEval;
        for k in 1..=9 {
            let t = c.t_start + (c.t_end - c.t_start) * (k as f64) / 10.0;
            let p = c.geometry.point_at(t);
            if let Some((u, v)) = eval.closest_point(&p) {
                worst = worst.max(eval.point_at(u, v).distance_to(&p));
            }
        }
    }
    worst
}
