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

/// Commit override for tests (three-state: unset → production
/// default ON, overridable by the POPKERN_SETBACK_VERTEX_BLEND=0
/// kill-switch).
pub(crate) static SETBACK_TOPOLOGY_COMMIT: AtomicBool = AtomicBool::new(false);
pub(crate) static SETBACK_TOPOLOGY_COMMIT_SET: AtomicBool = AtomicBool::new(false);

/// Whether the setback vertex blend replaces the legacy corner
/// machinery for eligible apexes. Production default: ON — the
/// surgery declines (with a recorded reason) and falls through to
/// the legacy path for anything out of its class, and the cap gates
/// must close at relative 1e-6 before any topology commits.
/// Kill-switch: POPKERN_SETBACK_VERTEX_BLEND=0.
pub(crate) fn setback_topology_commit_enabled() -> bool {
    if SETBACK_TOPOLOGY_COMMIT_SET.load(Ordering::Relaxed) {
        return SETBACK_TOPOLOGY_COMMIT.load(Ordering::Relaxed);
    }
    std::env::var("POPKERN_SETBACK_VERTEX_BLEND")
        .map(|v| v != "0")
        .unwrap_or(true)
}

/// WIP checkpoint report channel for the integration probe test.
pub(crate) static SETBACK_TOPOLOGY_LAST: std::sync::Mutex<Option<String>> =
    std::sync::Mutex::new(None);

/// Result of the setback boundary surgery.
#[derive(Debug)]
#[allow(dead_code)] // arc_edges/corner_vertices feed the upcoming acceptance sweep
pub(crate) struct SetbackTopologyOutcome {
    /// The merged blend face bounded by the n setback arcs.
    pub blend_face: FaceHandle,
    /// The n setback-arc edges, fillet order.
    pub arc_edges: Vec<EdgeHandle>,
    /// The n corner vertices (trim intersections), fillet order
    /// (corner k = between fillet k and fillet k+1).
    pub corner_vertices: Vec<VertexHandle>,
    /// The 4n KP20 piece faces subdividing the blend region.
    pub cap_faces: Vec<FaceHandle>,
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

    // ---- cap face fan: subdivide the blend n-gon into the 4n KP20
    // piece faces --------------------------------------------------------
    let cap_faces = build_cap_face_fan(
        brep,
        blend_face,
        &arc_edges,
        fillet_surfaces,
        vertex_position,
    )?;

    // ---- pcurves for the new topology -----------------------------------
    // Plane faces: rebuild loop pcurves from the 3D curves (the
    // corner path edges live on the lateral planes). Cap pieces and
    // sweep faces: the new edges follow iso/straight uv paths, so
    // endpoint-projected Linear pcurves are exact.
    author_setback_pcurves(brep, &cap_faces)?;

    Ok(SetbackTopologyOutcome {
        blend_face,
        arc_edges,
        corner_vertices,
        cap_faces,
    })
}

/// Sector chart (measured in the c0-exact campaign): o at (0,0), the
/// two bounding spokes along the chart's i=0 row and j=0 column,
/// m_s at (0,6), m_{s+1} at (6,0), the rim corner c at (6,6); pieces
/// P1/P4/P2/P3 are the chart quarters. Rim quarter-points come from
/// the CAP's own piece corners (the eq (14) band reparameterization
/// means they are NOT at naive arc parameters).
fn build_cap_face_fan(
    brep: &mut BRep,
    blend_face: FaceHandle,
    arc_edges: &[EdgeHandle],
    fillet_surfaces: &[FilletSurfaceResult],
    vertex_position: Point3,
) -> Result<Vec<FaceHandle>, String> {
    use popkern_core::geometry::gregory::Kp20PieceId;
    let n = arc_edges.len();
    // Assemble the cap (the same call the stage-1 candidate makes).
    let (border, _extract) =
        kp20_paper_border_from_fillet_setback(fillet_surfaces, vertex_position)?;
    let centroid = {
        let mut acc = [0.0_f64; 3];
        let mut count = 0.0;
        for sector in &border.sectors {
            for piece in sector {
                for cp in &piece.rows[0] {
                    acc[0] += cp.x;
                    acc[1] += cp.y;
                    acc[2] += cp.z;
                    count += 1.0;
                }
            }
        }
        Point3::new(acc[0] / count, acc[1] / count, acc[2] / count)
    };
    let star = Point3::new(
        centroid.x + (vertex_position.x - centroid.x) * 0.25,
        centroid.y + (vertex_position.y - centroid.y) * 0.25,
        centroid.z + (vertex_position.z - centroid.z) * 0.25,
    );
    let report = popkern_core::geometry::gregory::assemble_kp20_paper_cap_conditioned(
        &border,
        star,
        popkern_core::geometry::gregory::Kp20StripEnforcement::HardG2,
        0.0,
    )
    .map_err(|e| format!("cap assembly: {e}"))?;
    // Acceptance gate: the paper residuals must close at the relative
    // tolerance (the same criterion as the stage-1 candidate) before
    // any topology is committed.
    {
        let scale = border
            .sectors
            .iter()
            .flat_map(|sector| sector.iter())
            .flat_map(|piece| piece.rows.iter().flatten())
            .map(|cp| cp.distance_to(&star))
            .fold(0.0_f64, f64::max);
        let gate = 1.0e-6 * scale.max(f64::MIN_POSITIVE);
        if report.max_g1_residual > gate
            || report.max_g2_residual > gate
            || report.max_c0_residual > gate
            || report.max_border_residual > gate
            || report.min_regularity <= 0.0
        {
            return Err(format!(
                "cap gates open: g1={:.3e} g2={:.3e} c0={:.3e} border={:.3e} min_reg={:.3e} gate={gate:.3e}",
                report.max_g1_residual,
                report.max_g2_residual,
                report.max_c0_residual,
                report.max_border_residual,
                report.min_regularity,
            ));
        }
    }
    let sectors = &report.sectors;
    if sectors.len() != n {
        return Err(format!("cap has {} sectors for {n} arcs", sectors.len()));
    }

    const P1: Kp20PieceId = Kp20PieceId::NearOriginLeft;
    const P2: Kp20PieceId = Kp20PieceId::NearOriginRight;
    const P4: Kp20PieceId = Kp20PieceId::NearMidpointLeft;

    // Key chart positions per sector.
    let o_point = sectors[0].piece(P1).cps[0][0];
    let sector_m = |s: usize| sectors[s].piece(P4).cps[0][6]; // chart (0,6) = m_s
    let sector_center = |s: usize| sectors[s].piece(P1).cps[6][6]; // chart (3,3)
    let sector_spoke_mid_row = |s: usize| sectors[s].piece(P1).cps[0][6]; // chart (0,3)
    let rim_q_row = |s: usize| sectors[s].piece(P4).cps[6][6]; // chart (3,6): rim mid of [m_s→c]
    let rim_q_col = |s: usize| sectors[s].piece(P2).cps[6][6]; // chart (6,3): rim mid of [c→m_{s+1}]

    // 1. Split each arc's middle edge at the cap's rim subdivision
    //    points: m (arc midpoint = chart (0,6)) and the two rim
    //    quarter corners. Sector s's rim = [m_s→c_s] (second half of
    //    arc s) + [c_s→m_{s+1}] (first half of arc s+1), so arc k
    //    carries: q_col of sector k−1 (its first half) and m_k + q_row
    //    of sector k (its second half).
    if std::env::var_os("PK_SETBACK_PROBE").is_some() {
        let dist_to_loop = |pt: Point3| -> f64 {
            let mut best = f64::INFINITY;
            if let Ok(hes) = face_loop_half_edges(brep, blend_face) {
                for he in hes {
                    let edge = brep.half_edges[he].edge;
                    if let Some(c) = brep.edges[edge].curve.as_ref() {
                        best = best.min(curve_distance(&c.geometry, c.t_start, c.t_end, pt));
                    }
                }
            }
            best
        };
        for s in 0..1 {
            use popkern_core::geometry::gregory::Kp20PieceId as Pid;
            let sec = &sectors[s];
            for (label, pt) in [
                ("P1(0,0) o?", sec.piece(Pid::NearOriginLeft).cps[0][0]),
                ("P1(0,6)", sec.piece(Pid::NearOriginLeft).cps[0][6]),
                ("P1(6,0)", sec.piece(Pid::NearOriginLeft).cps[6][0]),
                ("P1(6,6)", sec.piece(Pid::NearOriginLeft).cps[6][6]),
                ("P4(0,6)", sec.piece(Pid::NearMidpointLeft).cps[0][6]),
                ("P4(6,6)", sec.piece(Pid::NearMidpointLeft).cps[6][6]),
                ("P4(6,0)", sec.piece(Pid::NearMidpointLeft).cps[6][0]),
                ("P2(6,6)", sec.piece(Pid::NearOriginRight).cps[6][6]),
                ("P2(6,0)", sec.piece(Pid::NearOriginRight).cps[6][0]),
                ("P2(0,6)", sec.piece(Pid::NearOriginRight).cps[0][6]),
                ("P3(0,0)", sec.piece(Pid::NearMidpointRight).cps[0][0]),
                ("P3(6,6)", sec.piece(Pid::NearMidpointRight).cps[6][6]),
                ("P3(0,6)", sec.piece(Pid::NearMidpointRight).cps[0][6]),
                ("P3(6,0)", sec.piece(Pid::NearMidpointRight).cps[6][0]),
            ] {
                eprintln!(
                    "  [chart] s{s} {label}: ({:+.3},{:+.3},{:+.3}) loop-dist {:.2e}",
                    pt.x,
                    pt.y,
                    pt.z,
                    dist_to_loop(pt)
                );
            }
        }
    }
    // Split at every sector's rim subdivision positions directly —
    // the extractor's sector order is angular and generally differs
    // from the fillet index order, so no arc↔sector pairing is
    // assumed; the splitter locates the right boundary edge globally.
    let mut m_vertices = Vec::with_capacity(n);
    for s in 0..n {
        let m_vertex = split_blend_boundary_at(brep, blend_face, sector_m(s))?;
        split_blend_boundary_at(brep, blend_face, rim_q_row(s))?;
        split_blend_boundary_at(brep, blend_face, rim_q_col(s))?;
        m_vertices.push(m_vertex);
    }

    // 2. Spokes: chains m_s → sm_s → o. The first spoke creates o via
    //    mev; later spokes close onto it with mef.
    let mut o_vertex: Option<VertexHandle> = None;
    let mut spoke_mid_vertices = Vec::with_capacity(n);
    for s in 0..n {
        let m_v = m_vertices[s];
        let he_m = half_edge_at_vertex_on_face_containing(brep, m_v, &o_point)?;
        let mev1 = brep
            .mev(he_m, sector_spoke_mid_row(s))
            .map_err(|e| format!("mev spoke mid {s}: {e}"))?;
        // Curve: chart (0,6)→(0,3) = P4 row i=0 (m → spoke mid).
        if let Some(edge) = brep.edges.get_mut(mev1.edge) {
            edge.curve = Some(bezier6_curve_between(
                piece_row(sectors[s].piece(P4), 0),
                brep.vertices[m_v].position,
                sector_spoke_mid_row(s),
            ));
        }
        spoke_mid_vertices.push(mev1.vertex);
        match o_vertex {
            None => {
                let he_sm = half_edge_at_vertex_on_face_containing(brep, mev1.vertex, &o_point)?;
                let mev2 = brep
                    .mev(he_sm, o_point)
                    .map_err(|e| format!("mev o from spoke {s}: {e}"))?;
                if let Some(edge) = brep.edges.get_mut(mev2.edge) {
                    edge.curve = Some(bezier6_curve_between(
                        piece_row(sectors[s].piece(P1), 0),
                        brep.vertices[mev1.vertex].position,
                        o_point,
                    ));
                }
                o_vertex = Some(mev2.vertex);
            }
            Some(o_v) => {
                let (he_sm, he_o) = half_edge_pair_in_common_loop(brep, mev1.vertex, o_v)?;
                let mef = brep
                    .mef(he_sm, he_o)
                    .map_err(|e| format!("mef spoke {s} to o: {e}"))?;
                if let Some(edge) = brep.edges.get_mut(mef.edge) {
                    edge.curve = Some(bezier6_curve_between(
                        piece_row(sectors[s].piece(P1), 0),
                        brep.vertices[mev1.vertex].position,
                        o_point,
                    ));
                }
            }
        }
    }
    let o_vertex = o_vertex.ok_or("no o vertex")?;

    // 3. Per sector: the interior cross (center vertex + 4 arms).
    for s in 0..n {
        let center = sector_center(s);
        let q_row_v = find_vertex_at(brep, rim_q_row(s))?;
        let q_col_v = find_vertex_at(brep, rim_q_col(s))?;
        let sm_this = spoke_mid_vertices[s];
        let sm_next = spoke_mid_vertices[(s + 1) % n];
        // Spur from the rim quarter (3,6) into the sector face.
        let he_q = half_edge_at_vertex_on_face_containing(brep, q_row_v, &center)?;
        let mev_c = brep
            .mev(he_q, center)
            .map_err(|e| format!("mev sector {s} center: {e}"))?;
        let center_v = mev_c.vertex;
        // Curve: chart (3,6)→(3,3) = P4 row i=6 (q_row → center).
        if let Some(edge) = brep.edges.get_mut(mev_c.edge) {
            edge.curve = Some(bezier6_curve_between(
                piece_row(sectors[s].piece(P4), 6),
                brep.vertices[q_row_v].position,
                center,
            ));
        }
        // Arm curves: center→sm_this = chart col j=3 i∈[0,3] = P1 col 6;
        // center→sm_next = chart row i=3 j∈[0,3] = P1 row 6;
        // center→q_col = chart (3,3)→(6,3) = P2 col j=6.
        let arm_curves = [
            piece_col(sectors[s].piece(P1), 6),
            piece_row(sectors[s].piece(P1), 6),
            piece_col(sectors[s].piece(P2), 6),
        ];
        for ((label, target), cps) in [
            ("sm_this", sm_this),
            ("sm_next", sm_next),
            ("q_col", q_col_v),
        ]
        .into_iter()
        .zip(arm_curves)
        {
            let (he_c, he_t) = half_edge_pair_in_common_loop(brep, center_v, target)
                .map_err(|e| format!("sector {s} arm {label}: {e}"))?;
            let mef = brep
                .mef(he_c, he_t)
                .map_err(|e| format!("mef sector {s} arm {label}: {e}"))?;
            if let Some(edge) = brep.edges.get_mut(mef.edge) {
                edge.curve = Some(bezier6_curve_between(
                    cps,
                    center,
                    brep.vertices[target].position,
                ));
            }
        }
    }
    let _ = o_vertex;

    // 4. Assign piece surfaces by centroid matching.
    let patches = super::kp20_polar_fold_summary::gregory_patches_from_macro_patches(sectors);
    let mut cap_faces = Vec::new();
    let mut assigned = vec![false; patches.len()];
    let face_handles: Vec<FaceHandle> = brep.faces.iter().map(|(fh, _)| fh).collect();
    for fh in face_handles {
        let Ok(c) = face_centroid(brep, fh) else {
            continue;
        };
        // A cap face's centroid is near one piece's mid-point.
        let mut best = (f64::INFINITY, usize::MAX);
        for (idx, patch) in patches.iter().enumerate() {
            if assigned[idx] {
                continue;
            }
            let popkern_core::geometry::gregory::GregoryPatch::RectangularBiNonic(bp) = patch
            else {
                continue;
            };
            let mid = bp.point_at(0.5, 0.5);
            let d = mid.distance_to(&c);
            if d < best.0 {
                best = (d, idx);
            }
        }
        if best.1 != usize::MAX && best.0 < 0.05 {
            let idx = best.1;
            if let Some(face) = brep.faces.get_mut(fh) {
                face.surface = Some(popkern_core::geometry::SurfaceOnFace {
                    geometry: popkern_core::geometry::SurfaceGeometry::Gregory(
                        patches[idx].clone(),
                    ),
                    orientation: true,
                });
                assigned[idx] = true;
                cap_faces.push(fh);
            }
        }
    }
    if cap_faces.len() != 4 * n {
        return Err(format!(
            "assigned {} cap faces, expected {}",
            cap_faces.len(),
            4 * n
        ));
    }
    Ok(cap_faces)
}

/// Author pcurves for the setback vertex-blend topology (see the
/// call site). Plane faces are reprojected wholesale; cap piece and
/// sweep half-edges missing a pcurve get endpoint-projected Linear
/// pcurves (exact for iso-line piece boundaries and the straight-uv
/// setback arcs).
fn author_setback_pcurves(brep: &mut BRep, cap_faces: &[FaceHandle]) -> Result<(), String> {
    use popkern_core::geometry::pcurve::{Pcurve, PcurveGeometry};
    use popkern_core::geometry::uv_projection::project_point_to_uv;
    // 1. Plane faces.
    let plane_faces: Vec<FaceHandle> = brep
        .faces
        .iter()
        .filter(|(_, f)| {
            matches!(
                f.surface.as_ref().map(|s| &s.geometry),
                Some(popkern_core::geometry::SurfaceGeometry::Plane { .. })
            )
        })
        .map(|(fh, _)| fh)
        .collect();
    for fh in plane_faces {
        super::topology_repair::topology_repair_geometry::reproject_loop_pcurves_into_plane(
            brep, fh,
        );
    }
    // 2. Cap pieces and any other curved faces with missing pcurves.
    let curved_faces: Vec<FaceHandle> = brep
        .faces
        .iter()
        .filter(|(fh, f)| {
            cap_faces.contains(fh)
                || matches!(
                    f.surface.as_ref().map(|s| &s.geometry),
                    Some(popkern_core::geometry::SurfaceGeometry::Nurbs(_))
                        | Some(popkern_core::geometry::SurfaceGeometry::Gregory(_))
                )
        })
        .map(|(fh, _)| fh)
        .collect();
    for fh in curved_faces {
        let surface = match brep.faces.get(fh).and_then(|f| f.surface.as_ref()) {
            Some(s) => s.geometry.clone(),
            None => continue,
        };
        let hes = face_loop_half_edges(brep, fh)?;
        for he in hes {
            if brep.half_edges[he].pcurve.is_some() {
                continue;
            }
            let a = brep.vertices[brep.half_edges[he].origin].position;
            let b = brep.vertices[brep.half_edges[brep.half_edges[he].next].origin].position;
            let (Some(ua), Some(ub)) = (
                project_point_to_uv(&surface, &a),
                project_point_to_uv(&surface, &b),
            ) else {
                continue;
            };
            if let Some(he_mut) = brep.half_edges.get_mut(he) {
                he_mut.pcurve = Some(Pcurve {
                    geometry: PcurveGeometry::Linear { start: ua, end: ub },
                    t_start: 0.0,
                    t_end: 1.0,
                });
            }
        }
    }
    Ok(())
}

/// Degree-6 Bezier from 7 control points, oriented start → end (the
/// CPs are reversed when needed so the curve runs from `from` to
/// `to`).
fn bezier6_curve_between(
    cps: [Point3; 7],
    from: Point3,
    to: Point3,
) -> popkern_core::geometry::CurveOnEdge {
    let mut pts = cps.to_vec();
    if pts[0].distance_to(&from) > pts[6].distance_to(&from) {
        pts.reverse();
    }
    // Snap endpoints to the topology vertices (rim vertices sit at
    // nearest-parameter positions ~1e-4 from the cap's ideal points;
    // unsnapped endpoints reopen shell gaps).
    pts[0] = from;
    pts[6] = to;
    popkern_core::geometry::CurveOnEdge {
        geometry: popkern_core::geometry::CurveGeometry::Nurbs(NurbsCurve {
            degree: 6,
            knots: vec![
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            ],
            control_points: pts,
            weights: vec![1.0; 7],
        }),
        t_start: 0.0,
        t_end: 1.0,
    }
}

fn piece_row(piece: &popkern_core::geometry::gregory::Kp20Bi6Patch, i: usize) -> [Point3; 7] {
    core::array::from_fn(|j| piece.cps[i][j])
}

fn piece_col(piece: &popkern_core::geometry::gregory::Kp20Bi6Patch, j: usize) -> [Point3; 7] {
    core::array::from_fn(|i| piece.cps[i][j])
}

/// Split the blend face's boundary at an exact point (which must lie
/// on one of its edge curves). Returns the (new or existing) vertex.
fn split_blend_boundary_at(
    brep: &mut BRep,
    face: FaceHandle,
    point: Point3,
) -> Result<VertexHandle, String> {
    for he in face_loop_half_edges(brep, face)? {
        let v = brep.half_edges[he].origin;
        if brep.vertices[v].position.distance_to(&point) < 1.0e-4 {
            return Ok(v);
        }
    }
    let mut best: Option<(f64, EdgeHandle)> = None;
    for he in face_loop_half_edges(brep, face)? {
        let edge = brep.half_edges[he].edge;
        let Some(curve) = brep.edges[edge].curve.as_ref() else {
            continue;
        };
        let d = curve_distance(&curve.geometry, curve.t_start, curve.t_end, point);
        if best.map(|(bd, _)| d < bd).unwrap_or(true) {
            best = Some((d, edge));
        }
    }
    let (d, edge) = best.ok_or("no curved boundary edge")?;
    if d > 1.0e-3 {
        return Err(format!("rim split point {d:.3e} off the boundary"));
    }
    // The generic-NURBS splitter lacks point inversion; compute the
    // parameter by projection (dense scan + ternary refine) and use
    // the parametric splitter.
    let (t_split, t_dist) = {
        let curve = brep.edges[edge].curve.as_ref().ok_or("edge lost curve")?;
        let eval: &dyn CurveEval = &curve.geometry;
        let (e0, e1) = (curve.t_start, curve.t_end);
        let mut best_t = (f64::INFINITY, e0);
        for i in 0..=512 {
            let t = e0 + (e1 - e0) * i as f64 / 512.0;
            let dd = eval.point_at(t).distance_to(&point);
            if dd < best_t.0 {
                best_t = (dd, t);
            }
        }
        let mut lo = (best_t.1 - (e1 - e0) / 512.0).max(e0.min(e1));
        let mut hi = (best_t.1 + (e1 - e0) / 512.0).min(e0.max(e1));
        for _ in 0..60 {
            let t1 = lo + (hi - lo) / 3.0;
            let t2 = hi - (hi - lo) / 3.0;
            if eval.point_at(t1).distance_to(&point) < eval.point_at(t2).distance_to(&point) {
                hi = t2;
            } else {
                lo = t1;
            }
        }
        let t = 0.5 * (lo + hi);
        (t, eval.point_at(t).distance_to(&point))
    };
    // The cap's rim points differ from the (endpoint-snapped) arc
    // sub-cubics by the band-reparameterization noise (~1e-4); the
    // split vertex lands on the topology's own curve at the nearest
    // parameter, which is the consistent choice.
    if t_dist > 5.0e-3 {
        return Err(format!("rim split projection residual {t_dist:.3e}"));
    }
    // split_edge takes the normalized fraction along the edge's
    // parameter window, not the raw curve parameter.
    let fraction = {
        let curve = brep.edges[edge].curve.as_ref().ok_or("edge lost curve")?;
        (t_split - curve.t_start) / (curve.t_end - curve.t_start)
    };
    let result = brep
        .split_edge(edge, fraction)
        .map_err(|e| format!("rim split: {e}"))?;
    if std::env::var_os("PK_SETBACK_PROBE").is_some() {
        let vp = brep.vertices[result.vertex].position;
        eprintln!(
            "  [rimsplit] target ({:+.4},{:+.4},{:+.4}) edge-dist {d:.2e} t={t_split:.4} t_dist {t_dist:.2e} vertex-off {:.2e}",
            point.x,
            point.y,
            point.z,
            vp.distance_to(&point)
        );
    }
    Ok(result.vertex)
}

/// The half-edge at `vertex` whose loop's face contains `interior`
/// (the face whose centroid is nearest to it).
fn half_edge_at_vertex_on_face_containing(
    brep: &BRep,
    vertex: VertexHandle,
    interior: &Point3,
) -> Result<HalfEdgeHandle, String> {
    let mut best: Option<(f64, HalfEdgeHandle)> = None;
    for (h, data) in brep.half_edges.iter() {
        if data.origin != vertex {
            continue;
        }
        let face = brep.loops[data.loop_].face;
        let Ok(c) = face_centroid(brep, face) else {
            continue;
        };
        let d = c.distance_to(interior);
        if best.map(|(bd, _)| d < bd).unwrap_or(true) {
            best = Some((d, h));
        }
    }
    best.map(|(_, h)| h)
        .ok_or("vertex has no half-edges".into())
}

/// Half-edges at `va` and `vb` that share one loop (for mef).
fn half_edge_pair_in_common_loop(
    brep: &BRep,
    va: VertexHandle,
    vb: VertexHandle,
) -> Result<(HalfEdgeHandle, HalfEdgeHandle), String> {
    for (ha, da) in brep.half_edges.iter() {
        if da.origin != va {
            continue;
        }
        let loop_ = da.loop_;
        let mut he = ha;
        for _ in 0..10_000 {
            he = brep.half_edges[he].next;
            if brep.half_edges[he].origin == vb && brep.half_edges[he].loop_ == loop_ {
                return Ok((ha, he));
            }
            if he == ha {
                break;
            }
        }
    }
    Err("no common loop for mef pair".into())
}

fn find_vertex_at(brep: &BRep, point: Point3) -> Result<VertexHandle, String> {
    // Rim vertices sit on the topology's own curves, within the band
    // reparameterization noise (~1e-4) of the cap's ideal positions.
    let mut best: Option<(f64, VertexHandle)> = None;
    for (vh, v) in brep.vertices.iter() {
        let d = v.position.distance_to(&point);
        if best.map(|(bd, _)| d < bd).unwrap_or(true) {
            best = Some((d, vh));
        }
    }
    match best {
        Some((d, vh)) if d < 1.0e-3 => Ok(vh),
        Some((d, _)) => Err(format!(
            "no vertex near ({:.4},{:.4},{:.4}); nearest {d:.3e}",
            point.x, point.y, point.z
        )),
        None => Err("empty vertex set".into()),
    }
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

fn bezier_cubic_curve(cps: &[Point3; 4]) -> NurbsCurve {
    NurbsCurve {
        degree: 3,
        knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
        control_points: cps.to_vec(),
        weights: vec![1.0; 4],
    }
}
