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
/// Which cap construction fills the blend n-gon.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum CapFamily {
    /// KP20 §4.4 paper cap: 4n bi-9 Gregory pieces, G2.
    Kp20G2,
    /// PnS bi-4 cap on the setback border: 2n Bezier patches, G1.
    PnsG1,
}

pub(crate) fn try_setback_boundary_surgery(
    brep: &mut BRep,
    vertex: VertexHandle,
    fillet_surfaces: &[FilletSurfaceResult],
    tolerance: &Tolerance,
    family: CapFamily,
) -> Result<SetbackTopologyOutcome, String> {
    let snapshot = brep.clone();
    match surgery_inner(brep, vertex, fillet_surfaces, tolerance, family) {
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
    family: CapFamily,
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
        // A crossing can land exactly on an existing vertex — with
        // trim-crossing anchors (PK-095/PK-097) the pre-carve already
        // split the boundary at the trim crossings, so reuse rather
        // than split.
        let split_or_reuse = |brep: &mut BRep, edge, point: Point3, which: &str| {
            if let Ok(v) = find_vertex_at(brep, point) {
                return Ok(v);
            }
            brep.split_edge_at_point(edge, point, tolerance)
                .map(|r| r.vertex)
                .map_err(|e| format!("split {which} of fillet {k}: {e}"))
        };
        let v1 = split_or_reuse(brep, e1, p1, "crossing 1")?;
        let v2 = split_or_reuse(brep, e2, p2, "crossing 2")?;
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

    // ---- surgery E: absorb perpendicular-split remnant bands -----------
    // The legacy pre-carve split each fillet face near the apex; the
    // arc mef leaves a narrow same-surface band between the arc and
    // the old split line. The blended fillet face must extend to the
    // arc (design: [arc, trim, base-end, trim]), so scaffold edges
    // between same-surface fillet pieces are dissolved.
    absorb_same_surface_bands(brep, &arc_edges, blend_face)?;

    // ---- cap face fan (family-specific) ---------------------------------
    let cap_faces = match family {
        CapFamily::Kp20G2 => {
            let faces = build_cap_face_fan(
                brep,
                blend_face,
                &arc_edges,
                fillet_surfaces,
                vertex_position,
            )?;
            // Rim refit: exact cap-row edge curves (PK-096). The PnS
            // family does not need this — its border equals the
            // topology arcs by construction.
            refit_cap_rim(brep, &faces)?;
            faces
        }
        CapFamily::PnsG1 => build_pns_cap_fan(
            brep,
            blend_face,
            fillet_surfaces,
            vertex_position,
            &corner_vertices,
            &arc_edges,
        )?,
    };

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
    // Keep every split's vertex handle — re-finding them later by
    // position fails for the G1 arc class (~2e-3 extraction noise vs
    // the G2 class's ~1e-4; PK-097 Phase C).
    let mut m_vertices = Vec::with_capacity(n);
    let mut q_row_vertices = Vec::with_capacity(n);
    let mut q_col_vertices = Vec::with_capacity(n);
    for s in 0..n {
        let m_vertex = split_blend_boundary_at(brep, blend_face, sector_m(s))?;
        q_row_vertices.push(split_blend_boundary_at(brep, blend_face, rim_q_row(s))?);
        q_col_vertices.push(split_blend_boundary_at(brep, blend_face, rim_q_col(s))?);
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
        let q_row_v = q_row_vertices[s];
        let q_col_v = q_col_vertices[s];
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
    let mut patches = super::kp20_polar_fold_summary::gregory_patches_from_macro_patches(sectors);
    // Sweep reconciliation (PK-096): the border's cubic row0 can only
    // follow the fillet surface to ~4e-4 over an arc half-span, and
    // the cap interpolates it exactly — a genuine C0 shell gap that
    // STEP consumers refuse to auto-sew. The bi-9 boundary has the
    // degrees of freedom the cubic lacks: project each rim-side
    // boundary row onto its fillet surface and translate the two
    // derivative rows in parallel (preserving the G1/G2 structure).
    let lateral_planes: Vec<popkern_core::geometry::SurfaceGeometry> = brep
        .faces
        .iter()
        .filter_map(|(_, f)| {
            let sof = f.surface.as_ref()?;
            matches!(
                sof.geometry,
                popkern_core::geometry::SurfaceGeometry::Plane { .. }
            )
            .then(|| sof.geometry.clone())
        })
        .collect();
    reconcile_cap_rim_to_sweeps(&mut patches, fillet_surfaces, &lateral_planes);
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
        let is_sweep = matches!(surface, popkern_core::geometry::SurfaceGeometry::Nurbs(_));
        let hes = face_loop_half_edges(brep, fh)?;
        for he in hes {
            // Sweep (fillet strip) faces: the recorded trim curves lie
            // INSIDE the sweep's parameter box (measured 0.03-0.08 off
            // the v = 0/1 iso borders), so iso-assuming pcurves are
            // wrong data — replace them with point-inversion
            // projections of the 3D curves. Cap pieces keep exact
            // endpoint-Linear pcurves for missing entries only (their
            // boundaries ARE iso lines).
            if is_sweep {
                let window = brep.half_edges[he].edge;
                let Some(curve) = brep.edges.get(window).and_then(|e| e.curve.clone()) else {
                    continue;
                };
                if !(curve.t_start.is_finite() && curve.t_end.is_finite()) {
                    continue;
                }
                let origin = brep.vertices[brep.half_edges[he].origin].position;
                use popkern_core::geometry::curve::CurveEval;
                let ceval: &dyn CurveEval = &curve.geometry;
                let forward = ceval.point_at(curve.t_start).distance_to(&origin)
                    <= ceval.point_at(curve.t_end).distance_to(&origin);
                const N: usize = 33;
                let mut cps = Vec::with_capacity(N);
                let mut ok = true;
                for i in 0..N {
                    let f = i as f64 / (N - 1) as f64;
                    let f = if forward { f } else { 1.0 - f };
                    let t = curve.t_start + (curve.t_end - curve.t_start) * f;
                    match project_point_to_uv(&surface, &ceval.point_at(t)) {
                        Some((u, v)) => cps.push(popkern_core::geometry::Point3::new(u, v, 0.0)),
                        None => {
                            ok = false;
                            break;
                        }
                    }
                }
                if !ok {
                    continue;
                }
                let mut knots = Vec::with_capacity(N + 2);
                knots.push(0.0);
                for i in 0..N {
                    knots.push(i as f64 / (N - 1) as f64);
                }
                knots.push(1.0);
                if let Some(he_mut) = brep.half_edges.get_mut(he) {
                    he_mut.pcurve = Some(Pcurve {
                        geometry: PcurveGeometry::Nurbs(NurbsCurve {
                            degree: 1,
                            knots,
                            control_points: cps,
                            weights: vec![1.0; N],
                        }),
                        t_start: 0.0,
                        t_end: 1.0,
                    });
                }
                continue;
            }
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

/// Move rim vertices onto the cap boundary rows and assign each rim
/// edge (cap piece ↔ non-cap face) the exact degree-9 Bezier segment
/// of the piece's boundary row between its (re-projected) endpoints.
/// Non-rim curves that end at a moved vertex have their nearest
/// endpoint control point translated to match.
fn refit_cap_rim(brep: &mut BRep, cap_faces: &[FaceHandle]) -> Result<(), String> {
    use popkern_core::geometry::gregory::GregoryPatch;
    // Collect (rim half-edge on cap side, row, ..) jobs first.
    struct RimJob {
        edge: EdgeHandle,
        /// Cap-side half-edge (the loop side on the cap face).
        cap_he: popkern_core::topology::HalfEdgeHandle,
        /// Which grid side of the piece (0: u=0, 1: u=9-row, 2: v=0
        /// col, 3: v=9 col) — determines the pcurve's iso line.
        side: usize,
        row: [Point3; 10],
        va: VertexHandle,
        vb: VertexHandle,
    }
    let mut jobs: Vec<RimJob> = Vec::new();
    for &fh in cap_faces {
        let Some(popkern_core::geometry::SurfaceGeometry::Gregory(
            GregoryPatch::RectangularBiNonic(patch),
        )) = brep
            .faces
            .get(fh)
            .and_then(|f| f.surface.as_ref())
            .map(|s| &s.geometry)
        else {
            continue;
        };
        let grid = patch.grid;
        let rows: [[Point3; 10]; 4] = [
            core::array::from_fn(|j| grid[0][j]),
            core::array::from_fn(|j| grid[9][j]),
            core::array::from_fn(|i| grid[i][0]),
            core::array::from_fn(|i| grid[i][9]),
        ];
        for he in face_loop_half_edges(brep, fh)? {
            let data = &brep.half_edges[he];
            let twin_face = brep.loops[brep.half_edges[data.twin].loop_].face;
            let twin_is_cap = cap_faces.contains(&twin_face);
            if twin_is_cap && fh > twin_face {
                // Interior cap↔cap edges are refit once (from the
                // lower-handle side); both grids agree on the shared
                // row after sweep reconciliation.
                continue;
            }
            let va = data.origin;
            let vb = brep.half_edges[data.next].origin;
            let pa = brep.vertices[va].position;
            let pb = brep.vertices[vb].position;
            // Pick the row by the edge's MIDPOINT (endpoints sit at
            // patch corners where two rows tie; the midpoint is
            // interior to exactly one side).
            let mid = match brep.edges[data.edge].curve.as_ref() {
                Some(c) => {
                    let eval: &dyn CurveEval = &c.geometry;
                    eval.point_at(0.5 * (c.t_start + c.t_end))
                }
                None => Point3::new(
                    0.5 * (pa.x + pb.x),
                    0.5 * (pa.y + pb.y),
                    0.5 * (pa.z + pb.z),
                ),
            };
            let mut best: Option<(f64, usize)> = None;
            for (ri, row) in rows.iter().enumerate() {
                let d = bezier9_nearest(row, mid).1;
                if best.map(|(bd, _)| d < bd).unwrap_or(true) {
                    best = Some((d, ri));
                }
            }
            let (_, ri) = best.ok_or("no rim row")?;
            jobs.push(RimJob {
                edge: data.edge,
                cap_he: he,
                side: ri,
                row: rows[ri],
                va,
                vb,
            });
        }
    }

    // Move rim vertices onto their rows (average when a vertex is
    // claimed by multiple rows — e.g. c between two sectors; the rows
    // agree at corners to extraction accuracy).
    let mut targets: std::collections::HashMap<VertexHandle, Vec<Point3>> =
        std::collections::HashMap::new();
    for job in &jobs {
        for v in [job.va, job.vb] {
            let p = brep.vertices[v].position;
            let (t, _) = bezier9_nearest(&job.row, p);
            targets
                .entry(v)
                .or_default()
                .push(bezier9_eval(&job.row, t));
        }
    }
    // A vertex may move ONLY if every incident non-rim curve can
    // follow it: full-window NURBS whose end control point sits at
    // the vertex. Split trims carry WINDOWED curves (their endpoint
    // is curve(t_start), not a control point) — translating a CP
    // would deform the wrong span, and skipping it breaks loop
    // closure (measured 4.7e-4 gaps), so those vertices stay put and
    // the rim segments snap to them instead.
    let rim_edge_set: std::collections::HashSet<EdgeHandle> = jobs.iter().map(|j| j.edge).collect();
    let movable = |v: VertexHandle| -> bool {
        let old = brep.vertices[v].position;
        for (eh, edge) in brep.edges.iter() {
            let he = edge.half_edge;
            let va = brep.half_edges[he].origin;
            let vb = brep.half_edges[brep.half_edges[he].twin].origin;
            if va != v && vb != v {
                continue;
            }
            if rim_edge_set.contains(&eh) {
                continue;
            }
            let Some(curve) = edge.curve.as_ref() else {
                return false;
            };
            let popkern_core::geometry::CurveGeometry::Nurbs(n) = &curve.geometry else {
                return false;
            };
            let cp_end_at_v = n
                .control_points
                .first()
                .map(|p| p.distance_to(&old) < 1.0e-6)
                .unwrap_or(false)
                || n.control_points
                    .last()
                    .map(|p| p.distance_to(&old) < 1.0e-6)
                    .unwrap_or(false);
            if !cp_end_at_v {
                return false;
            }
        }
        true
    };
    let movable_set: std::collections::HashSet<VertexHandle> =
        targets.keys().copied().filter(|&v| movable(v)).collect();
    let mut moves: std::collections::HashMap<VertexHandle, (Point3, Point3)> =
        std::collections::HashMap::new();
    for (v, list) in targets {
        if !movable_set.contains(&v) {
            continue;
        }
        let old = brep.vertices[v].position;
        let mut acc = [0.0_f64; 3];
        for p in &list {
            acc[0] += p.x;
            acc[1] += p.y;
            acc[2] += p.z;
        }
        let n = list.len() as f64;
        let new = Point3::new(acc[0] / n, acc[1] / n, acc[2] / n);
        if old.distance_to(&new) > 1.0e-12 {
            brep.vertices[v].position = new;
            moves.insert(v, (old, new));
        }
    }

    // Assign exact row segments to the rim edges.
    let rim_edges: std::collections::HashSet<EdgeHandle> = jobs.iter().map(|j| j.edge).collect();
    for job in &jobs {
        let pa = brep.vertices[job.va].position;
        let pb = brep.vertices[job.vb].position;
        let (ta, _) = bezier9_nearest(&job.row, pa);
        let (tb, _) = bezier9_nearest(&job.row, pb);
        // The row's parameterization may run against the half-edge:
        // extract the ordered span and reverse the control points.
        let mut seg = if ta <= tb {
            bezier9_extract(&job.row, ta, tb)
        } else {
            let mut rev = bezier9_extract(&job.row, tb, ta);
            rev.reverse();
            rev
        };
        // Endpoints exactly at the vertices, with the deltas ramped
        // linearly across all control points (a degree-1 correction
        // in Bernstein form) so the segment stays smooth.
        let da = pa - bezier9_eval(&seg, 0.0);
        let db = pb - bezier9_eval(&seg, 1.0);
        for (i, cp) in seg.iter_mut().enumerate() {
            let f = i as f64 / 9.0;
            *cp = Point3::new(
                cp.x + da.x * (1.0 - f) + db.x * f,
                cp.y + da.y * (1.0 - f) + db.y * f,
                cp.z + da.z * (1.0 - f) + db.z * f,
            );
        }
        let mut knots = Vec::with_capacity(20);
        knots.extend(std::iter::repeat_n(0.0, 10));
        knots.extend(std::iter::repeat_n(1.0, 10));
        if let Some(edge) = brep.edges.get_mut(job.edge) {
            edge.curve = Some(popkern_core::geometry::CurveOnEdge {
                geometry: popkern_core::geometry::CurveGeometry::Nurbs(NurbsCurve {
                    degree: 9,
                    knots,
                    control_points: seg.to_vec(),
                    weights: vec![1.0; 10],
                }),
                t_start: 0.0,
                t_end: 1.0,
            });
        }
        // Cap-side pcurve: the row segment IS the iso line of the
        // piece domain between ta and tb, and de Casteljau extraction
        // parameterizes it linearly — so an exact, parameter-
        // synchronized pcurve is the straight domain segment. (The
        // STEP exporter demands trace synchronization; asynchronous
        // pcurves force face-local edges and an OPEN_SHELL export.)
        // grid[i][j] evaluates with B_i(u)·B_j(v): sides 0/1 are u
        // iso lines parameterized by v, sides 2/3 v iso lines by u.
        let (uv_a, uv_b) = match job.side {
            0 => ((0.0, ta), (0.0, tb)),
            1 => ((1.0, ta), (1.0, tb)),
            2 => ((ta, 0.0), (tb, 0.0)),
            _ => ((ta, 1.0), (tb, 1.0)),
        };
        if let Some(he_mut) = brep.half_edges.get_mut(job.cap_he) {
            he_mut.pcurve = Some(popkern_core::geometry::pcurve::Pcurve {
                geometry: popkern_core::geometry::pcurve::PcurveGeometry::Linear {
                    start: uv_a,
                    end: uv_b,
                },
                t_start: 0.0,
                t_end: 1.0,
            });
        }
    }

    // Repair endpoint coherence for NON-rim curves at moved vertices:
    // translate whichever curve end is nearest the old position.
    let edge_keys: Vec<EdgeHandle> = brep.edges.keys().collect();
    for eh in edge_keys {
        if rim_edges.contains(&eh) {
            continue;
        }
        let he = brep.edges[eh].half_edge;
        let va = brep.half_edges[he].origin;
        let vb = brep.half_edges[brep.half_edges[he].twin].origin;
        for v in [va, vb] {
            let Some(&(old, new)) = moves.get(&v) else {
                continue;
            };
            let Some(curve) = brep.edges.get_mut(eh).and_then(|e| e.curve.as_mut()) else {
                continue;
            };
            if let popkern_core::geometry::CurveGeometry::Nurbs(n) = &mut curve.geometry {
                let front = n.control_points.first().map(|p| p.distance_to(&old));
                let back = n.control_points.last().map(|p| p.distance_to(&old));
                match (front, back) {
                    (Some(f), Some(b)) if f <= b && f < 1.0e-6 => {
                        *n.control_points.first_mut().ok_or("empty curve")? = new;
                    }
                    (Some(_), Some(b)) if b < 1.0e-6 => {
                        *n.control_points.last_mut().ok_or("empty curve")? = new;
                    }
                    _ => {}
                }
            }
        }
    }
    Ok(())
}

fn bezier9_eval(cps: &[Point3; 10], t: f64) -> Point3 {
    let mut work = *cps;
    for level in (1..10).rev() {
        for i in 0..level {
            work[i] = Point3::new(
                work[i].x + (work[i + 1].x - work[i].x) * t,
                work[i].y + (work[i + 1].y - work[i].y) * t,
                work[i].z + (work[i + 1].z - work[i].z) * t,
            );
        }
    }
    work[0]
}

/// Nearest parameter on a degree-9 Bezier (dense scan + ternary).
fn bezier9_nearest(cps: &[Point3; 10], p: Point3) -> (f64, f64) {
    let mut best = (f64::INFINITY, 0.0_f64);
    for i in 0..=128 {
        let t = i as f64 / 128.0;
        let d = bezier9_eval(cps, t).distance_to(&p);
        if d < best.0 {
            best = (d, t);
        }
    }
    let (mut lo, mut hi) = (
        (best.1 - 1.0 / 128.0).max(0.0),
        (best.1 + 1.0 / 128.0).min(1.0),
    );
    for _ in 0..60 {
        let m1 = lo + (hi - lo) / 3.0;
        let m2 = hi - (hi - lo) / 3.0;
        if bezier9_eval(cps, m1).distance_to(&p) < bezier9_eval(cps, m2).distance_to(&p) {
            hi = m2;
        } else {
            lo = m1;
        }
    }
    let t = 0.5 * (lo + hi);
    (t, bezier9_eval(cps, t).distance_to(&p))
}

/// De Casteljau extraction of the sub-segment [ta, tb].
fn bezier9_extract(cps: &[Point3; 10], ta: f64, tb: f64) -> [Point3; 10] {
    // Split at ta, keep right; split that at (tb-ta)/(1-ta), keep left.
    let right = bezier9_split_right(cps, ta);
    let f = if (1.0 - ta).abs() < 1e-15 {
        0.0
    } else {
        (tb - ta) / (1.0 - ta)
    };
    bezier9_split_left(&right, f)
}

fn bezier9_split_right(cps: &[Point3; 10], t: f64) -> [Point3; 10] {
    // Right-segment control points are the terminal point of each de
    // Casteljau level: [C(t), b_1^8, ..., b_8^1, b_9].
    let mut work = *cps;
    let mut right = [cps[9]; 10];
    for level in (1..10).rev() {
        for i in 0..level {
            work[i] = Point3::new(
                work[i].x + (work[i + 1].x - work[i].x) * t,
                work[i].y + (work[i + 1].y - work[i].y) * t,
                work[i].z + (work[i + 1].z - work[i].z) * t,
            );
        }
        right[level - 1] = work[level - 1];
    }
    right[0] = work[0];
    right
}

fn bezier9_split_left(cps: &[Point3; 10], t: f64) -> [Point3; 10] {
    let mut work = *cps;
    let mut left = [cps[0]; 10];
    for (level, slot) in (1..10).rev().zip(1..10) {
        for i in 0..level {
            work[i] = Point3::new(
                work[i].x + (work[i + 1].x - work[i].x) * t,
                work[i].y + (work[i + 1].y - work[i].y) * t,
                work[i].z + (work[i + 1].z - work[i].z) * t,
            );
        }
        left[slot] = work[0];
    }
    left
}

/// Project the rim-side boundary rows of the bi-9 cap pieces onto the
/// fillet sweeps (see call site). A grid side is a rim side when every
/// sampled point lies within `RIM_NEAR` of some fillet surface; the
/// projection deltas are applied to the boundary row and, in parallel,
/// to the two rows behind it. Deltas at shared corners are computed
/// from identical input points, so adjacent pieces stay sewn.
fn reconcile_cap_rim_to_sweeps(
    patches: &mut [popkern_core::geometry::gregory::GregoryPatch],
    fillet_surfaces: &[crate::fillet_surface::FilletSurfaceResult],
    lateral_planes: &[popkern_core::geometry::SurfaceGeometry],
) {
    use popkern_core::geometry::gregory::GregoryPatch;
    use popkern_core::geometry::uv_projection::project_point_to_uv;
    const RIM_NEAR: f64 = 2.0e-3;

    // The rim's neighbors are the sweeps along the arcs AND the
    // lateral planes along the corner segments (past the trim
    // crossings X the border leaves the fillet surface). Candidates
    // are enumerated so callers can track WHICH surface a sample
    // used: min-distance alone jumps between the two planes meeting
    // at a corner c, and interpolating the resulting discontinuous
    // field rings at the millimetre level.
    let candidates: Vec<popkern_core::geometry::SurfaceGeometry> = fillet_surfaces
        .iter()
        .map(|f| f.surface.clone())
        .chain(lateral_planes.iter().cloned())
        .collect();
    let project_on = |idx: usize, p: Point3| -> Option<(f64, Point3)> {
        let sfc = &candidates[idx];
        let (u, v) = project_point_to_uv(sfc, &p)?;
        let q = sfc.point_at(u, v);
        Some((q.distance_to(&p), q))
    };
    let project_best = |p: Point3| -> Option<(usize, f64, Point3)> {
        let mut best: Option<(usize, f64, Point3)> = None;
        for idx in 0..candidates.len() {
            if let Some((d, q)) = project_on(idx, p)
                && best.map(|(_, bd, _)| d < bd).unwrap_or(true)
            {
                best = Some((idx, d, q));
            }
        }
        best.filter(|(_, d, _)| *d < RIM_NEAR)
    };
    // Continuation projection: keep the previous sample's surface
    // while it stays near; switch only when it drifts (the genuine
    // arc→plane handoff at X happens where both surfaces agree).
    let project_seq = |pts: &[Point3]| -> Option<Vec<Point3>> {
        let mut out = Vec::with_capacity(pts.len());
        let mut current: Option<usize> = None;
        for &p in pts {
            let pick = match current {
                Some(idx) => match project_on(idx, p) {
                    Some((d, q)) if d < RIM_NEAR => Some((idx, d, q)),
                    _ => project_best(p),
                },
                None => project_best(p),
            };
            let (idx, _, q) = pick?;
            current = Some(idx);
            out.push(q);
        }
        Some(out)
    };

    for patch in patches.iter_mut() {
        let GregoryPatch::RectangularBiNonic(bp) = patch else {
            continue;
        };
        // Four sides: (fixed index, along index, is_row, fixed at 0 or 9).
        // Represent as closures over (k) returning grid coordinates for
        // depth d in 0..3.
        for side in 0..4 {
            // Gather the boundary row for this side.
            let coord = |d: usize, k: usize| -> (usize, usize) {
                match side {
                    0 => (d, k),     // row i=0 boundary, depth into +i
                    1 => (9 - d, k), // row i=9 boundary, depth into -i
                    2 => (k, d),     // col j=0 boundary, depth into +j
                    _ => (k, 9 - d), // col j=9 boundary, depth into -j
                }
            };
            // Rim test on CURVE SAMPLES: control points legitimately
            // sit 1-2e-3 off their own curve, so a CP-based gate
            // wrongly rejects true rim sides.
            let side_row: [Point3; 10] = core::array::from_fn(|k| {
                let (i, j) = coord(0, k);
                bp.grid[i][j]
            });
            let gate_pts: Vec<Point3> = (0..=16)
                .map(|i| bezier9_eval(&side_row, i as f64 / 16.0))
                .collect();
            if project_seq(&gate_pts).is_none() {
                continue;
            }
            // Land the boundary CURVE on the sweep: projecting control
            // points directly leaves the curve ~3e-4 off (CPs are not
            // curve points). Sample the deviation along the curve at
            // 10 parameters and solve the degree-9 Bernstein
            // interpolation for the exact CP deltas; iterate (the
            // projection field is mildly nonlinear).
            for _pass in 0..3 {
                let row: [Point3; 10] = core::array::from_fn(|k| {
                    let (i, j) = coord(0, k);
                    bp.grid[i][j]
                });
                let ts: [f64; 10] = core::array::from_fn(|k| k as f64 / 9.0);
                let node_pts: Vec<Point3> = ts.iter().map(|&t| bezier9_eval(&row, t)).collect();
                let Some(projected) = project_seq(&node_pts) else {
                    break;
                };
                let mut rhs = [[0.0_f64; 3]; 10];
                for k in 0..10 {
                    let (p, q) = (node_pts[k], projected[k]);
                    rhs[k] = [q.x - p.x, q.y - p.y, q.z - p.z];
                }
                // Bernstein Vandermonde B[k][i] = B_i^9(t_k).
                let mut mat = [[0.0_f64; 10]; 10];
                for (k, &t) in ts.iter().enumerate() {
                    mat[k] = bernstein9_all(t);
                }
                let Some(deltas) = solve_10(&mat, &rhs) else {
                    break;
                };
                for k in 0..10 {
                    for d in 0..3 {
                        let (i, j) = coord(d, k);
                        let g = bp.grid[i][j];
                        bp.grid[i][j] =
                            Point3::new(g.x + deltas[k][0], g.y + deltas[k][1], g.z + deltas[k][2]);
                    }
                }
            }
        }
    }
}

fn bernstein9_all(t: f64) -> [f64; 10] {
    const C: [f64; 10] = [1.0, 9.0, 36.0, 84.0, 126.0, 126.0, 84.0, 36.0, 9.0, 1.0];
    let mut out = [0.0_f64; 10];
    let u = 1.0 - t;
    let mut tp = 1.0;
    for (i, slot) in out.iter_mut().enumerate() {
        *slot = C[i] * tp * u.powi(9 - i as i32);
        tp *= t;
    }
    out
}

/// Solve the 10×10 linear system `mat · x = rhs` (3 right-hand sides)
/// by Gaussian elimination with partial pivoting.
fn solve_10(mat: &[[f64; 10]; 10], rhs: &[[f64; 3]; 10]) -> Option<[[f64; 3]; 10]> {
    let mut a = *mat;
    let mut b = *rhs;
    for col in 0..10 {
        let mut piv = col;
        for r in (col + 1)..10 {
            if a[r][col].abs() > a[piv][col].abs() {
                piv = r;
            }
        }
        if a[piv][col].abs() < 1e-14 {
            return None;
        }
        a.swap(col, piv);
        b.swap(col, piv);
        for r in (col + 1)..10 {
            let f = a[r][col] / a[col][col];
            for c in col..10 {
                a[r][c] -= f * a[col][c];
            }
            for c in 0..3 {
                b[r][c] -= f * b[col][c];
            }
        }
    }
    let mut x = [[0.0_f64; 3]; 10];
    for row in (0..10).rev() {
        for c in 0..3 {
            let mut acc = b[row][c];
            for col in (row + 1)..10 {
                acc -= a[row][col] * x[col][c];
            }
            x[row][c] = acc / a[row][row];
        }
    }
    Some(x)
}

/// Subdivide the blend n-gon into the 2n PnS sector faces (PK-097
/// Phase B): one quad per half-arc sector, all meeting at the cap's
/// pole. The candidate's border equals the topology arcs exactly
/// (Phase A), so no rim refit is needed; spoke curves are the
/// patches' shared radial rows.

/// Split a specific edge at the parameter nearest `point` (dense scan
/// + ternary refine; `split_edge` takes the normalized fraction of
/// the window). Reuses an existing vertex when the projection lands
/// on an endpoint.
fn split_edge_at_projected_point(
    brep: &mut BRep,
    edge: EdgeHandle,
    point: Point3,
) -> Result<VertexHandle, String> {
    let curve = brep
        .edges
        .get(edge)
        .and_then(|e| e.curve.as_ref())
        .ok_or("edge has no curve")?
        .clone();
    let eval: &dyn CurveEval = &curve.geometry;
    let (e0, e1) = (curve.t_start, curve.t_end);
    let mut best = (f64::INFINITY, e0);
    for i in 0..=512 {
        let t = e0 + (e1 - e0) * i as f64 / 512.0;
        let d = eval.point_at(t).distance_to(&point);
        if d < best.0 {
            best = (d, t);
        }
    }
    let mut lo = (best.1 - (e1 - e0) / 512.0).max(e0.min(e1));
    let mut hi = (best.1 + (e1 - e0) / 512.0).min(e0.max(e1));
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
    let fraction = (t - e0) / (e1 - e0);
    if !(0.001..0.999).contains(&fraction) {
        // Landed on an endpoint — reuse the nearer vertex.
        let p = eval.point_at(t);
        return nearest_vertex_within(brep, p, 5.0e-3).ok_or_else(|| {
            format!(
                "projection landed on an endpoint with no vertex nearby (fraction {fraction:.4})"
            )
        });
    }
    brep.split_edge(edge, fraction)
        .map(|r| r.vertex)
        .map_err(|e| format!("split_edge: {e}"))
}

/// Nearest live vertex within `gate` of `p`.
fn nearest_vertex_within(brep: &BRep, p: Point3, gate: f64) -> Option<VertexHandle> {
    brep.vertices
        .iter()
        .map(|(vh, vd)| (vd.position.distance_to(&p), vh))
        .min_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal))
        .filter(|(d, _)| *d < gate)
        .map(|(_, vh)| vh)
}

fn build_pns_cap_fan(
    brep: &mut BRep,
    blend_face: FaceHandle,
    fillet_surfaces: &[FilletSurfaceResult],
    vertex_position: Point3,
    corner_vertices: &[VertexHandle],
    arc_edges: &[EdgeHandle],
) -> Result<Vec<FaceHandle>, String> {
    let candidate =
        super::pns_setback_candidate::try_pns_setback_candidate(fillet_surfaces, vertex_position)?;
    // Acceptance before topology (the same discipline as the KP20
    // fan's gate): border exact, seams exact, G1 within the
    // production coplanarity gate.
    let r = &candidate.report;
    if r.g0_border_after > 1.0e-9 || r.seam_g0 > 1.0e-9 || r.g1_angle_after > 5.0e-2 {
        return Err(format!(
            "PnS candidate gates open: g0={:.3e} seam={:.3e} g1={:.3e}",
            r.g0_border_after, r.seam_g0, r.g1_angle_after
        ));
    }
    let sectors = candidate.half_cubics.len();

    // 1. Split the rim at the arc midpoints m (the junction between
    // each arc's two half-cubic sectors). The target EDGE is a known
    // topological fact (the surgery's arc middles, indexed by
    // fillet) — a global nearest-edge search occasionally split the
    // wrong arc's neighbor.
    let mut junction_vertices: Vec<VertexHandle> = Vec::new();
    for k in (0..sectors).step_by(2) {
        let m_pos = candidate.half_cubics[k].0[3];
        let fillet_idx = candidate.half_cubics[k].1;
        let edge = *arc_edges
            .get(fillet_idx)
            .ok_or_else(|| format!("no arc edge for fillet {fillet_idx}"))?;
        let v = split_edge_at_projected_point(brep, edge, m_pos)?;
        junction_vertices.push(v);
    }
    // Corner junctions are TOPOLOGICAL facts (the surgery's corner
    // vertices) — match each even sector's start to the nearest of
    // them; odd sectors start at the m splits made above.
    let mut all_junctions: Vec<(Point3, VertexHandle)> = Vec::new();
    for k in 0..sectors {
        let start = candidate.half_cubics[k].0[0];
        let v = if k % 2 == 0 {
            // Handles can go stale through the merge/dissolve passes;
            // consider only live ones and fall back to a positional
            // search at the family tolerance.
            let live_nearest = corner_vertices
                .iter()
                .filter_map(|&vh| {
                    brep.vertices
                        .get(vh)
                        .map(|vd| (vd.position.distance_to(&start), vh))
                })
                .min_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
            match live_nearest {
                Some((d, vh)) if d < 5.0e-3 => vh,
                _ => nearest_vertex_within(brep, start, 5.0e-3)
                    .ok_or_else(|| format!("no corner vertex near sector {k} start"))?,
            }
        } else {
            junction_vertices[k / 2]
        };
        let pos = brep.vertices[v].position;
        all_junctions.push((pos, v));
    }

    // 2. Pole spur from junction 0, then spokes from every other
    // junction. Spoke curve for junction k = patch k's radial row
    // i=0, running border→pole (grid [0][4] → [0][0]).
    // The assembly's requested pole seeds the fit, but the patches'
    // ACTUAL shared corner is where they meet — use it for the pole
    // vertex and the spoke snaps (measured 3e-2 apart on the G1
    // pentagon; snapping to the seed bent every spoke off the caps).
    let pole = {
        let seed = candidate.pole;
        let p0 = &candidate.patches[0];
        let corners = [p0[0][0], p0[0][4], p0[4][0], p0[4][4]];
        corners
            .into_iter()
            .min_by(|a, b| {
                a.distance_to(&seed)
                    .partial_cmp(&b.distance_to(&seed))
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .unwrap_or(seed)
    };
    let spoke_curve = |patch: &[[Point3; 5]; 5],
                       neighbor: &popkern_core::geometry::SurfaceGeometry,
                       from: Point3,
                       to: Point3|
     -> popkern_core::geometry::CurveOnEdge {
        // The PnS patches alternate grid orientation sector to
        // sector (their radial seams pair row i=0 against col j=0),
        // so the spoke side is selected empirically: the grid side
        // whose endpoints best match (junction, pole).
        let sides: [[Point3; 5]; 4] = [
            core::array::from_fn(|j| patch[0][j]),
            core::array::from_fn(|j| patch[4][j]),
            core::array::from_fn(|i| patch[i][0]),
            core::array::from_fn(|i| patch[i][4]),
        ];
        // The defining property of the seam side: its interior lies
        // ON the neighboring patch (radial seams are surface-exact,
        // 1e-15 measured). Score = endpoint match + midpoint's
        // distance to the neighbor surface.
        let mut best: (f64, [Point3; 5]) = (f64::INFINITY, sides[0]);
        for side in &sides {
            for rev in [false, true] {
                let (a, b) = if rev {
                    (side[4], side[0])
                } else {
                    (side[0], side[4])
                };
                let mid = {
                    let mut w = *side;
                    for level in (1..5).rev() {
                        for i in 0..level {
                            w[i] = Point3::new(
                                w[i].x + (w[i + 1].x - w[i].x) * 0.5,
                                w[i].y + (w[i + 1].y - w[i].y) * 0.5,
                                w[i].z + (w[i + 1].z - w[i].z) * 0.5,
                            );
                        }
                    }
                    w[0]
                };
                let on_neighbor = match popkern_core::geometry::uv_projection::project_point_to_uv(
                    neighbor, &mid,
                ) {
                    Some((u, v)) => neighbor.point_at(u, v).distance_to(&mid),
                    None => f64::INFINITY,
                };
                let d = a.distance_to(&from).max(b.distance_to(&to)) + on_neighbor;
                if d < best.0 {
                    let mut arr = *side;
                    if rev {
                        arr.reverse();
                    }
                    best = (d, arr);
                }
            }
        }
        let mut cps: Vec<Point3> = best.1.to_vec();
        cps[0] = from;
        cps[4] = to;
        popkern_core::geometry::CurveOnEdge {
            geometry: popkern_core::geometry::CurveGeometry::Nurbs(NurbsCurve {
                degree: 4,
                knots: vec![0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0],
                control_points: cps,
                weights: vec![1.0; 5],
            }),
            t_start: 0.0,
            t_end: 1.0,
        }
    };

    let (first_pos, first_v) = all_junctions[0];
    let he_first = half_edge_at_vertex_on_face_containing(brep, first_v, &pole)?;
    let mev = brep
        .mev(he_first, pole)
        .map_err(|e| format!("PnS pole mev: {e}"))?;
    let pole_vertex = mev.vertex;
    let neighbor_of = |k: usize| &candidate.surfaces[(k + sectors - 1) % sectors];
    if let Some(edge) = brep.edges.get_mut(mev.edge) {
        edge.curve = Some(spoke_curve(
            &candidate.patches[0],
            neighbor_of(0),
            first_pos,
            pole,
        ));
    }
    for (k, &(pos, v)) in all_junctions.iter().enumerate().skip(1) {
        let (he_v, he_pole) = half_edge_pair_in_common_loop(brep, v, pole_vertex)?;
        let mef = brep
            .mef(he_v, he_pole)
            .map_err(|e| format!("PnS spoke mef {k}: {e}"))?;
        if let Some(edge) = brep.edges.get_mut(mef.edge) {
            edge.curve = Some(spoke_curve(
                &candidate.patches[k],
                neighbor_of(k),
                pos,
                pole,
            ));
        }
    }

    // 3. Assign sector surfaces. Candidate faces are exactly those
    // whose loop touches the pole vertex; match faces↔patches by
    // globally sorted centroid distances (per-face greedy steals
    // neighbors' patches on strongly curved quads).
    let mut fan_faces: Vec<FaceHandle> = Vec::new();
    for (fh, _) in brep.faces.iter() {
        if fh == blend_face || face_loop_half_edges(brep, fh).is_err() {
            // blend_face handle may be one of the subdivided faces —
            // include it via the loop test below instead.
        }
        let Ok(hes) = face_loop_half_edges(brep, fh) else {
            continue;
        };
        if hes
            .iter()
            .any(|&he| brep.half_edges[he].origin == pole_vertex)
        {
            fan_faces.push(fh);
        }
    }
    if fan_faces.len() != sectors {
        return Err(format!(
            "PnS fan produced {} pole-adjacent faces, expected {sectors}",
            fan_faces.len()
        ));
    }
    // Face↔patch pairing must use INTERIOR evidence: the sector
    // vertices (junctions, pole, X) are shared by several patches,
    // so vertex- or centroid-fit is underdetermined (measured: every
    // spoke curve landed on one adjacent face and 5e-2 off the
    // other). Each fan face's RIM edge midpoint lies on exactly one
    // half-arc — the patch index directly.
    let fan_set: std::collections::HashSet<FaceHandle> = fan_faces.iter().copied().collect();
    let mut cap_faces = Vec::new();
    let mut patch_taken = vec![false; sectors];
    for &fh in &fan_faces {
        // Rim edges: boundary half-edges whose twin face is outside
        // the fan. Use the BEST (mid, half-arc) pair over all of the
        // face's rim edges — a face's first rim edge in loop order
        // can be a short corner sliver.
        // Use the LONGEST rim edge: the short corner segs carry the
        // crossing-foot snap noise (up to 5e-3) and mis-key near the
        // shared corners; the arc portions are long and central.
        let mut longest: Option<(f64, Point3)> = None;
        for he in face_loop_half_edges(brep, fh)? {
            let data = &brep.half_edges[he];
            let twin_face = brep.loops[brep.half_edges[data.twin].loop_].face;
            if fan_set.contains(&twin_face) {
                continue;
            }
            let Some(curve) = brep.edges[data.edge].curve.as_ref() else {
                continue;
            };
            let ceval: &dyn CurveEval = &curve.geometry;
            let a = ceval.point_at(curve.t_start);
            let b = ceval.point_at(curve.t_end);
            let len = a.distance_to(&b);
            let rim_mid = ceval.point_at(0.5 * (curve.t_start + curve.t_end));
            if longest.map(|(l, _)| len > l).unwrap_or(true) {
                longest = Some((len, rim_mid));
            }
        }
        let mut best = (f64::INFINITY, usize::MAX);
        if let Some((_, rim_mid)) = longest {
            for (pi, (cubic, _)) in candidate.half_cubics.iter().enumerate() {
                let mut d = f64::INFINITY;
                for i in 0..=32 {
                    let t = i as f64 / 32.0;
                    let mut w = *cubic;
                    for level in (1..4).rev() {
                        for j in 0..level {
                            w[j] = Point3::new(
                                w[j].x + (w[j + 1].x - w[j].x) * t,
                                w[j].y + (w[j + 1].y - w[j].y) * t,
                                w[j].z + (w[j + 1].z - w[j].z) * t,
                            );
                        }
                    }
                    d = d.min(w[0].distance_to(&rim_mid));
                }
                if d < best.0 {
                    best = (d, pi);
                }
            }
        }
        let pi = best.1;
        if pi == usize::MAX || best.0 > 5.0e-3 {
            return Err(format!(
                "fan face {fh:?}: no half-arc within 5e-3 of its rim (best {:.3e})",
                best.0
            ));
        }
        if patch_taken[pi] {
            return Err(format!(
                "fan face {fh:?}: half-arc {pi} claimed twice — rim keying ambiguous"
            ));
        }
        patch_taken[pi] = true;
        if let Some(face) = brep.faces.get_mut(fh) {
            face.surface = Some(popkern_core::geometry::SurfaceOnFace {
                geometry: candidate.surfaces[pi].clone(),
                orientation: true,
            });
            cap_faces.push(fh);
        }
    }
    if cap_faces.len() != sectors {
        return Err(format!(
            "PnS fan assigned {} faces, expected {sectors}",
            cap_faces.len()
        ));
    }
    Ok(cap_faces)
}

/// Dissolve scaffold edges between fillet-face pieces that carry the
/// SAME sweep surface (perpendicular-split remnants). An edge
/// qualifies when its two faces are distinct, neither is the blend
/// face, and both faces' Nurbs surfaces share the same control-point
/// grid corners (same underlying sweep).
fn absorb_same_surface_bands(
    brep: &mut BRep,
    arc_edges: &[EdgeHandle],
    blend_face: FaceHandle,
) -> Result<(), String> {
    let _ = arc_edges;
    loop {
        let mut target: Option<popkern_core::topology::HalfEdgeHandle> = None;
        'scan: for (he, data) in brep.half_edges.iter() {
            let f1 = brep.loops[data.loop_].face;
            let twin = data.twin;
            let f2 = brep.loops[brep.half_edges[twin].loop_].face;
            if f1 == f2 || f1 == blend_face || f2 == blend_face {
                continue;
            }
            let (Some(s1), Some(s2)) = (
                brep.faces.get(f1).and_then(|f| f.surface.as_ref()),
                brep.faces.get(f2).and_then(|f| f.surface.as_ref()),
            ) else {
                continue;
            };
            let (
                popkern_core::geometry::SurfaceGeometry::Nurbs(n1),
                popkern_core::geometry::SurfaceGeometry::Nurbs(n2),
            ) = (&s1.geometry, &s2.geometry)
            else {
                continue;
            };
            if n1.control_points.len() != n2.control_points.len() {
                continue;
            }
            let corners = |n: &popkern_core::geometry::NurbsSurface| {
                let (nu, nv) = (n.num_u, n.num_v);
                [
                    n.control_points[0],
                    n.control_points[nu - 1],
                    n.control_points[(nv - 1) * nu],
                    n.control_points[nv * nu - 1],
                ]
            };
            let (c1, c2) = (corners(n1), corners(n2));
            if c1
                .iter()
                .zip(c2.iter())
                .all(|(a, b)| a.distance_to(b) < 1.0e-12)
            {
                target = Some(he);
                break 'scan;
            }
        }
        let Some(he) = target else { break };
        brep.kef(he)
            .map_err(|e| format!("absorb same-surface band: {e}"))?;
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
    // Gate matches the projection-residual gate below: the cap rim
    // and the topology boundary disagree by band-reparameterization
    // noise plus the trim-crossing anchor shift (PK-095), ~1e-3; the
    // split lands on the topology's own curve either way.
    if d > 5.0e-3 {
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
        // Acceptance: the split lands at the EDGE's own foot point,
        // so the arc-to-edge distance is a recorded quality margin,
        // not a correctness condition. The G1 cylinder strips carry
        // their setback corners a few 1e-3 from the trim extent
        // (PK-097) — 5e-3 covers both families.
        if best.0 < 5.0e-3 {
            hits.push((best.1, edge, best.2, best.0));
        } else if std::env::var_os("PK_CARVE_TRACE").is_some() {
            eprintln!(
                "    [crossings] edge {edge:?}: best approach {:.3e} at arc t {:.3}",
                best.0, best.1
            );
        }
    }
    if std::env::var_os("PK_CARVE_TRACE").is_some() {
        let a0 = arc_eval.point_at(0.0);
        let a1 = arc_eval.point_at(1.0);
        eprintln!(
            "    [crossings] arc ({:.3},{:.3},{:.3})->({:.3},{:.3},{:.3}): {} hits",
            a0.x,
            a0.y,
            a0.z,
            a1.x,
            a1.y,
            a1.z,
            hits.len()
        );
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
