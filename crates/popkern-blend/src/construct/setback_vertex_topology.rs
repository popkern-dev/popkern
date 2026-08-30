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
//! the recorded follow-up, so production keeps the vertex-blend path until
//! that lands.

use std::sync::atomic::{AtomicBool, Ordering};

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

use super::kp20_paper_border_extract::{Kp20SetbackArc, kp20_paper_border_from_fillet_setback};
use crate::fillet_surface::FilletSurfaceResult;

/// Commit override for tests (unset → production default ON).
pub(crate) static SETBACK_TOPOLOGY_COMMIT: AtomicBool = AtomicBool::new(false);
pub(crate) static SETBACK_TOPOLOGY_COMMIT_SET: AtomicBool = AtomicBool::new(false);

#[path = "setback_arc_binding.rs"]
mod setback_arc_binding;
/// Whether the setback vertex blend replaces the vertex-blend corner
/// machinery for eligible apexes. Production default: ON — the
/// surgery declines (with a recorded reason) and falls through to
/// the vertex-blend path for anything out of its class, and the cap gates
/// must close at relative 1e-6 before any topology commits.
#[path = "setback_rim_refit.rs"]
mod setback_rim_refit;
use setback_rim_refit::refit_cap_rim;
#[path = "setback_arc_crossings.rs"]
mod setback_arc_crossings;
#[path = "setback_bezier9.rs"]
mod setback_bezier9;
#[path = "setback_boundary_tolerance.rs"]
mod setback_boundary_tolerance;
#[path = "setback_fan_seam_refresh.rs"]
mod setback_fan_seam_refresh;
#[path = "setback_rim_jet.rs"]
mod setback_rim_jet;
pub(crate) use setback_bezier9::{bezier9_eval, bezier9_extract, bezier9_nearest};
#[path = "boundary_ring.rs"]
mod boundary_ring;
#[path = "boundary_ring_geom.rs"]
mod boundary_ring_geom;
#[path = "boundary_ring_recut.rs"]
mod boundary_ring_recut;
#[path = "setback_rim_reconcile.rs"]
mod setback_rim_reconcile;
#[path = "setback_shared_net_unify.rs"]
mod setback_shared_net_unify;
pub(crate) use boundary_ring::try_boundary_ring_carve;
use setback_arc_crossings::arc_boundary_crossings;

/// Env-gated (the fan-quality trace topic) pipeline-stage observability for
/// committed cap fans: worst live-seam g0/g1 of the corner-patch
/// faces. This tracer localized the pentagon5_blunt post-pipeline
/// pcurve clobber to `corner_patch_support_completion` (PK-087).
pub(crate) fn fan_quality_stage_probe(
    label: &str,
    brep: &BRep,
    corner_patches: &[crate::types::CornerPatchInfo],
) {
    if !popkern_core::trace::trace_enabled("fan-quality") {
        return;
    }
    let faces: Vec<FaceHandle> = corner_patches.iter().map(|cp| cp.face).collect();
    let (g0, g1, g2) = setback_fan_seam_refresh::worst_fan_seam_quality(brep, &faces);
    eprintln!("  [fan-quality] {label}: g0={g0:.3e} g1={g1:.3e} g2={g2:.3e}");
}

pub(crate) fn setback_topology_commit_enabled() -> bool {
    if SETBACK_TOPOLOGY_COMMIT_SET.load(Ordering::Relaxed) {
        return SETBACK_TOPOLOGY_COMMIT.load(Ordering::Relaxed);
    }
    true
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
/// caller's structured diagnostics; the vertex-blend path then proceeds.
/// `snapshot` is the caller-held pre-surgery state (the rung's commit
/// snapshot): the surgery restores from it on any inner error via a
/// clone-on-restore, so the success path pays no snapshot clone of
/// its own (PK-278 — the second per-corner whole-model clone was pure
/// duplication of the rung's).
#[allow(clippy::too_many_arguments)]
pub(crate) fn try_setback_boundary_surgery(
    brep: &mut BRep,
    snapshot: &BRep,
    vertex: VertexHandle,
    fillet_surfaces: &[FilletSurfaceResult],
    tolerance: &Tolerance,
    gate_seam_g2: bool,
    cap_presolve: Option<&crate::construct::setback_cap_presolve::SetbackCapPresolveEntry>,
) -> Result<SetbackTopologyOutcome, String> {
    match surgery_inner(
        brep,
        vertex,
        fillet_surfaces,
        tolerance,
        gate_seam_g2,
        cap_presolve,
    ) {
        Ok(outcome) => Ok(outcome),
        Err(e) => {
            *brep = snapshot.clone();
            Err(e)
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn surgery_inner(
    brep: &mut BRep,
    vertex: VertexHandle,
    fillet_surfaces: &[FilletSurfaceResult],
    tolerance: &Tolerance,
    gate_seam_g2: bool,
    cap_presolve: Option<&crate::construct::setback_cap_presolve::SetbackCapPresolveEntry>,
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
    // candidate assembles the cap on). PK-278: consume the parallel
    // presolve when its recorded inputs match this call's inputs
    // exactly; otherwise solve inline (byte-identical either way —
    // the extractor is a pure function of these two arguments).
    let presolved: Option<&crate::construct::setback_cap_presolve::SetbackCapSolve> =
        match cap_presolve.and_then(|entry| entry.matching(fillet_surfaces, vertex_position)) {
            Some(Ok(solve)) => Some(solve),
            Some(Err(e)) => return Err(e.clone()),
            None => None,
        };
    let extract = match presolved {
        Some(solve) => solve.extract.clone(),
        None => kp20_paper_border_from_fillet_setback(fillet_surfaces, vertex_position)?.1,
    };
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
    if popkern_core::trace::trace_enabled("setback") {
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
            let ff = setback_arc_binding::find_fillet_face(
                brep,
                &fillet_surfaces[k],
                &arcs[k].cubic,
                vertex_position,
            )?;
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
        .map(|k| {
            setback_arc_binding::find_fillet_face(
                brep,
                &fillet_surfaces[k],
                &arcs[k].cubic,
                vertex_position,
            )
        })
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
        // than split. Reuse is restricted to vertices ON THIS FACE'S
        // loop: a global nearest-vertex search grabbed unrelated
        // vertices on hexagon6_mixed, whose skew corners (trim gaps
        // to 1.2e-3) pack distinct crossing points inside the reuse
        // band and the stolen vertex then fails downstream wiring
        // ("vertex not on face loop").
        let split_or_reuse = |brep: &mut BRep, edge, point: Point3, which: &str| {
            let mut nearest_on_loop: Option<(f64, VertexHandle)> = None;
            if let Ok(hes) = face_loop_half_edges(brep, fillet_faces[k]) {
                for he in hes {
                    let v = brep.half_edges[he].origin;
                    if let Some(data) = brep.vertices.get(v) {
                        let d = data.position.distance_to(&point);
                        if nearest_on_loop.map(|(bd, _)| d < bd).unwrap_or(true) {
                            nearest_on_loop = Some((d, v));
                        }
                    }
                }
            }
            if let Some((d, v)) = nearest_on_loop
                && d < 1.0e-3
            {
                return Ok(v);
            }
            brep.split_edge_at_point(edge, point, tolerance)
                .map(|r| r.vertex)
                .map_err(|e| format!("split {which} of fillet {k}: {e}"))
        };
        let v1 = split_or_reuse(brep, e1, p1, "crossing 1")?;
        let v2 = split_or_reuse(brep, e2, p2, "crossing 2")?;
        if popkern_core::trace::trace_enabled("setback") {
            let pv1 = brep.vertices[v1].position;
            let pv2 = brep.vertices[v2].position;
            eprintln!(
                "  [crossings] arc {k} face {:?} t1={t1:.4} t2={t2:.4} p1=({:+.3},{:+.3},{:+.3}) p2=({:+.3},{:+.3},{:+.3}) v1-snap {:.2e} v2-snap {:.2e}",
                fillet_faces[k],
                p1.x,
                p1.y,
                p1.z,
                p2.x,
                p2.y,
                p2.z,
                pv1.distance_to(&p1),
                pv2.distance_to(&p2),
            );
        }
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
    let mut corner_ring_edges: Vec<EdgeHandle> = Vec::new();
    for k in 0..n {
        let k_next = (k + 1) % n;
        let x2 = crossing_vertices[k][1].ok_or("missing X2")?;
        let x1n = crossing_vertices[k_next][0].ok_or("missing X1'")?;
        if x2 == x1n {
            // Adjacent arcs' crossings reused the same pre-carve
            // anchor vertex — the ring already closes through it and
            // a corner spur would only add a degenerate flap that
            // breaks the ring's simple-cycle separation (measured on
            // hexagon6_variable: 33v1 -> 39v1 -> 33v1 excursions).
            corner_vertices.push(x2);
            continue;
        }
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
        corner_ring_edges.push(mev.edge);
        corner_ring_edges.push(mef.edge);
        corner_vertices.push(mev.vertex);
    }
    let arc_edges = arc_middle_edges;

    // ---- surgery D: merge the surrendered region ------------------------
    let mut surrendered: Vec<FaceHandle> = Vec::new();
    surrendered.extend(apex_slivers.iter().copied());
    surrendered.extend(wedge_slivers.iter().copied());
    let mut ring_edges: std::collections::BTreeSet<EdgeHandle> =
        arc_edges.iter().copied().collect();
    ring_edges.extend(corner_ring_edges.iter().copied());
    let blend_face = merge_faces(brep, vertex, &surrendered, &ring_edges)
        .map_err(|e| format!("blend region merge: {e}"))?;

    // ---- surgery E: absorb perpendicular-split remnant bands -----------
    // The earlier pre-carve split each fillet face near the apex; the
    // arc mef leaves a narrow same-surface band between the arc and
    // the old split line. The blended fillet face must extend to the
    // arc (design: [arc, trim, base-end, trim]), so scaffold edges
    // between same-surface fillet pieces are dissolved.
    absorb_same_surface_bands(brep, &arc_edges, blend_face)?;
    // ---- cap face fan ---------------------------------------------------
    let cap_faces = {
        let faces = build_cap_face_fan(
            brep,
            blend_face,
            &arc_edges,
            fillet_surfaces,
            vertex_position,
            presolved,
        )?;
        // Rim refit: exact cap-row edge curves (PK-096).
        refit_cap_rim(brep, &faces)?;
        faces
    };

    // ---- pcurves for the new topology -----------------------------------
    // Plane faces: rebuild loop pcurves from the 3D curves (the
    // corner path edges live on the lateral planes). Cap pieces and
    // sweep faces: the new edges follow iso/straight uv paths, so
    // endpoint-projected Linear pcurves are exact.
    author_setback_pcurves(brep, &cap_faces)?;

    // Commit postcondition: the fan's live seams must be commercial.
    // The all-or-nothing contract belongs HERE, on measured seam
    // quality, not only on structural validity — a blunt apex (18°
    // half-angle, pentagon5_blunt) builds a structurally valid fan
    // whose wide flat sectors crease 2.9e-1 at one seam corner, and
    // declining lets the sphere-fit rung serve the gentle
    // corner instead. Sharp-apex fixtures measure g0 ≤ 2e-10 and
    // g1 ≤ 2e-3 here — orders under the gates.
    let (seam_g0, seam_g1, seam_g2) =
        setback_fan_seam_refresh::worst_fan_seam_quality(brep, &cap_faces);
    if popkern_core::trace::trace_enabled("setback") {
        eprintln!("  [seam-postcondition] g0={seam_g0:.3e} g1={seam_g1:.3e} g2={seam_g2:.3e}");
    }
    if seam_g0 > 1.0e-5 || seam_g1 > 5.0e-2 {
        return Err(format!(
            "fan seam postcondition failed: g0={seam_g0:.3e} (gate 1e-5) \
             g1={seam_g1:.3e} (gate 5e-2)"
        ));
    }
    // Fold postcondition: the fast structural validate() the commit
    // runs excludes self-intersection, so a fan whose pieces FOLD
    // would serve and only die later at the Part serving gate
    // (measured on the N=3 cube corner: 3 near-pole pieces fold with
    // interior normal reversals to dot -1.0 while the assembly's
    // unsigned min-regularity 8.7e-5 slid past a gate inflated by the
    // creased-border corner mismatch). Decline HERE with the class
    // named; the vertex-blend path serves the corner instead.
    // Read-only per-piece measurement; ordered collect keeps the
    // first-folded-face report identical to the serial walk (PK-302).
    let fold_of = |fh: &FaceHandle| -> bool {
        brep.faces
            .get(*fh)
            .and_then(|f| f.surface.as_ref())
            .is_some_and(|sof| match &sof.geometry {
                popkern_core::geometry::SurfaceGeometry::Gregory(patch) => {
                    patch.self_intersects_with_tolerance(brep.tolerance.linear().max(1e-6))
                }
                _ => false,
            })
    };
    #[cfg(feature = "parallel")]
    let folds: Vec<bool> = {
        use rayon::prelude::*;
        cap_faces.par_iter().map(fold_of).collect()
    };
    #[cfg(not(feature = "parallel"))]
    let folds: Vec<bool> = cap_faces.iter().map(fold_of).collect();
    for (&fh, folded) in cap_faces.iter().zip(folds) {
        if folded {
            return Err(format!(
                "fan fold postcondition failed: cap piece {fh:?} self-intersects \
                 (assembly-born fold — creased-border wide-sector class)"
            ));
        }
    }
    // Strip-conflict postcondition (PK-225): the setback extends the
    // strips to the arcs, and flat conic sections at small radii
    // genuinely interpenetrate near the shared contact corner. The
    // committed arcs are the strips' new corner-side boundaries —
    // a mid-arc contact with a NEIGHBOUR strip means the kept
    // domains overlap (see setback_commit_gate module docs for the
    // falsified alternatives).
    if let Some(conflict) = super::setback_commit_gate::flat_section_class(fillet_surfaces) {
        return Err(format!("strip-conflict postcondition failed: {conflict}"));
    }
    // A G2 REQUEST additionally gates the seam's normal-curvature
    // residual (PK-122): a blunt apex (18° half-angle, pentagon5_blunt)
    // builds a G0/G1-commercial fan whose wide flat sectors still jump
    // 7.0e-2 in curvature across seams; committing it would report an
    // achieved G2 the geometry does not have. Declining hands the
    // gentle corner to the sphere-fit rung. Sharp-apex fixtures
    // measure well inside the gate (the kp_polar commercial gates
    // assert the same 5e-2 bound on their live seams).
    if gate_seam_g2 && seam_g2 > 5.0e-2 {
        return Err(format!(
            "fan seam G2 postcondition failed under a G2 request: \
             g2={seam_g2:.3e} (gate 5e-2, g0={seam_g0:.3e} g1={seam_g1:.3e})"
        ));
    }

    // Record the fit residual the rim reconciliation legitimately
    // left as per-edge tolerances so the serving curve-surface
    // validation (check 9) can accept the fan — see
    // `record_boundary_tolerances` for the contract and for why this
    // must never feed the tessellator.
    setback_boundary_tolerance::record_boundary_tolerances(brep, &cap_faces);

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
    presolved: Option<&crate::construct::setback_cap_presolve::SetbackCapSolve>,
) -> Result<Vec<FaceHandle>, String> {
    use popkern_core::geometry::gregory::Kp20PieceId;
    let n = arc_edges.len();
    // Assemble the cap (the same call the stage-1 candidate makes).
    // PK-278: the border/assembly pair is pure in (fillet_surfaces,
    // vertex_position); when the caller verified a presolve match the
    // precomputed values are bitwise what the inline calls produce.
    let (border, report) = match presolved {
        Some(solve) => (
            solve.border.clone(),
            solve.assembly.as_ref().map_err(Clone::clone)?.clone(),
        ),
        None => {
            let (border, _extract) =
                kp20_paper_border_from_fillet_setback(fillet_surfaces, vertex_position)?;
            let star =
                crate::construct::setback_cap_presolve::kp20_cap_star(&border, vertex_position);
            let report = popkern_core::geometry::gregory::assemble_kp20_paper_cap_conditioned(
                &border,
                star,
                popkern_core::geometry::gregory::Kp20StripEnforcement::HardG2,
                0.0,
            )
            .map_err(|e| format!("cap assembly: {e}"))?;
            (border, report)
        }
    };
    let star = crate::construct::setback_cap_presolve::kp20_cap_star(&border, vertex_position);
    // Acceptance gate: the paper residuals must close at the relative
    // tolerance (the same criterion as the stage-1 candidate) before
    // any topology is committed. The floor is the INPUT corner
    // mismatch: skew adjacent trims (octagon8_mixed: no exact trim
    // intersection exists, corners are closest-approach midpoints)
    // bound how far any derivation can close, so demanding better
    // than the input rejects an assembly that is exact up to its
    // input (measured: g2 2.1e-7 against a 4.0e-8 pure-scale gate).
    // Fold, regularity, and the live-seam commercial gates still
    // protect the committed quality downstream.
    {
        let scale = border
            .sectors
            .iter()
            .flat_map(|sector| sector.iter())
            .flat_map(|piece| piece.rows.iter().flatten())
            .map(|cp| cp.distance_to(&star))
            .fold(0.0_f64, f64::max);
        let gate = (1.0e-6 * scale.max(f64::MIN_POSITIVE)).max(10.0 * report.input_corner_mismatch);
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
    if popkern_core::trace::trace_enabled("setback") {
        eprintln!(
            "  [assembly] min_reg={:.3e} g1={:.3e} g2={:.3e} c0={:.3e} border={:.3e} strip_move={:.3e} input_corner={:.3e}",
            report.min_regularity,
            report.max_g1_residual,
            report.max_g2_residual,
            report.max_c0_residual,
            report.max_border_residual,
            report.max_strip_movement,
            report.input_corner_mismatch,
        );
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
    if popkern_core::trace::trace_enabled("setback") {
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
        for (s, sec) in sectors.iter().enumerate().take(1) {
            use popkern_core::geometry::gregory::Kp20PieceId as Pid;
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
        // [rim] probe: per sector, sample the four landed rim legs and
        // measure (a) distance to the committed topology ring, (b)
        // distance to the nearest fillet strip surface, (c) distance
        // to the INPUT border rim (the arc cubics). Discriminates
        // "conditioner slid the rim along the strip" (benign — ring
        // can follow) from "conditioner left the strip" (G0 wall).
        {
            use popkern_core::geometry::gregory::Kp20PieceId as Pid;
            let bez7 = |cps: &[Point3; 7], t: f64| -> Point3 {
                let mut p = *cps;
                for level in (1..7).rev() {
                    for i in 0..level {
                        p[i] = Point3::new(
                            p[i].x + (p[i + 1].x - p[i].x) * t,
                            p[i].y + (p[i + 1].y - p[i].y) * t,
                            p[i].z + (p[i + 1].z - p[i].z) * t,
                        );
                    }
                }
                p[0]
            };
            let bez4 = |cps: &[Point3; 4], t: f64| -> Point3 {
                let mut p = *cps;
                for level in (1..4).rev() {
                    for i in 0..level {
                        p[i] = Point3::new(
                            p[i].x + (p[i + 1].x - p[i].x) * t,
                            p[i].y + (p[i + 1].y - p[i].y) * t,
                            p[i].z + (p[i + 1].z - p[i].z) * t,
                        );
                    }
                }
                p[0]
            };
            let strip_dist = |p: Point3| -> f64 {
                let mut best = f64::INFINITY;
                for fillet in fillet_surfaces {
                    if let Some((u, v)) = popkern_core::geometry::uv_projection::project_point_to_uv(
                        &fillet.surface,
                        &p,
                    ) {
                        best = best.min(fillet.surface.point_at(u, v).distance_to(&p));
                    }
                }
                best
            };
            let input_rim_dist = |p: Point3| -> f64 {
                let mut best = f64::INFINITY;
                for sector in &border.sectors {
                    for piece in sector {
                        for i in 0..=128 {
                            let t = i as f64 / 128.0;
                            best = best.min(bez4(&piece.rows[0], t).distance_to(&p));
                        }
                    }
                }
                best
            };
            for (s, sec) in sectors.iter().enumerate() {
                for (label, cps) in [
                    ("A1 P4j6", piece_col(sec.piece(Pid::NearMidpointLeft), 6)),
                    ("A2 P3j6", piece_col(sec.piece(Pid::NearMidpointRight), 6)),
                    ("B1 P3i6", piece_row(sec.piece(Pid::NearMidpointRight), 6)),
                    ("B2 P2i6", piece_row(sec.piece(Pid::NearOriginRight), 6)),
                ] {
                    let mut worst = (0.0_f64, 0.0_f64);
                    for i in 0..=8 {
                        let t = i as f64 / 8.0;
                        let d = dist_to_loop(bez7(&cps, t));
                        if d > worst.0 {
                            worst = (d, t);
                        }
                    }
                    let pw = bez7(&cps, worst.1);
                    eprintln!(
                        "  [rim] s{s} {label}: worst-loop {:.2e} @t={:.2} strip {:.2e} input-rim {:.2e}",
                        worst.0,
                        worst.1,
                        strip_dist(pw),
                        input_rim_dist(pw),
                    );
                    // For badly-off rim points, dump the per-edge ring
                    // distances to localize the missing image span.
                    if worst.0 > 5.0e-3
                        && let Ok(hes) = face_loop_half_edges(brep, blend_face)
                    {
                        for he in hes {
                            let edge = brep.half_edges[he].edge;
                            let v0 = brep.half_edges[he].origin;
                            let p0 = brep.vertices[v0].position;
                            let d = brep.edges[edge].curve.as_ref().map_or(f64::INFINITY, |c| {
                                curve_distance(&c.geometry, c.t_start, c.t_end, pw)
                            });
                            eprintln!(
                                "    [ring-edge] {edge:?} from ({:+.3},{:+.3},{:+.3}) dist {d:.2e}",
                                p0.x, p0.y, p0.z
                            );
                        }
                    }
                }
            }
        }
    }
    // Split at every sector's rim subdivision positions directly —
    // the extractor's sector order is angular and generally differs
    // from the fillet index order, so no arc↔sector pairing is
    // assumed; the splitter locates the right boundary edge globally.
    // Keep every split's vertex handle — re-finding them later by
    // position fails for the G1 arc class (~2e-3 extraction noise vs
    // the G2 class's ~1e-4; PK-097 the consumed-merge stage).
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
    // SEMANTIC rim identification: the rim sides are the rows that
    // interpolate the setback border loop. Sweep proximity alone
    // misclassifies cap-internal spoke seams as rim sides where
    // adjacent fillets pinch the cap corner (hexagon6_variable's
    // equal-radius sector boundaries put whole spoke seams 1.8e-3
    // from the sweeps), and reconciling those seams onto the sweeps
    // folds the cap interior (measured cap_outer g1 7.8e-1). A spoke
    // seam touches the border loop only at one endpoint, so a
    // max-distance-to-loop gate separates the classes exactly.
    let border_loop: Vec<Vec<Point3>> = border
        .sectors
        .iter()
        .flat_map(|sector| sector.iter())
        .map(|piece| {
            (0..=32)
                .map(|i| setback_rim_reconcile::eval_cubic(&piece.rows[0], i as f64 / 32.0))
                .collect()
        })
        .collect();
    // Rim moves are per-patch solves; re-unify the shared boundary
    // CPs afterwards or the fan tears at the spoke seams' rim ends
    // (see setback_shared_net_unify docs).
    let shared_groups = setback_shared_net_unify::record_shared_boundary_groups(&patches);
    // The conditioned assembly's strip enforcement may move the
    // boundary rows off the INPUT border loop (measured 5.06e-3 on
    // the wide pentagon r=0.16 while HardG2 strip movement reported
    // ~1e-2). The rim-side classifier must tolerate that measured
    // movement or genuinely-rim sides are skipped, leaving whole
    // patches unreconciled and creasing every seam they share with a
    // reconciled neighbour (PK-094: fan seam g1 1.1e-1 ramp along
    // one seam). Spoke seams stay excluded by orders of magnitude —
    // their off-loop distance is half the loop diameter.
    let rim_off_loop_allowance = 2.0e-3_f64.max(2.0 * report.max_strip_movement);
    setback_rim_reconcile::reconcile_cap_rim_to_sweeps(
        &mut patches,
        fillet_surfaces,
        &lateral_planes,
        &border_loop,
        rim_off_loop_allowance,
    );
    // Rim cross-jet reconciliation (PK-212): meet the strips G1 at
    // the landed border — see setback_rim_jet.
    {
        let candidates: Vec<popkern_core::geometry::SurfaceGeometry> = fillet_surfaces
            .iter()
            .map(|f| f.surface.clone())
            .chain(lateral_planes.iter().cloned())
            .collect();
        setback_rim_jet::reconcile_rim_jets(
            &mut patches,
            &candidates,
            fillet_surfaces.len(),
            &border_loop,
            rim_off_loop_allowance,
        );
    }
    let rim_tear =
        setback_shared_net_unify::reunify_shared_boundary_cps(&mut patches, &shared_groups);
    if popkern_core::trace::trace_enabled("setback") {
        eprintln!("  [rim-unify] worst shared-CP tear closed: {rim_tear:.3e}");
    }
    // Face ↔ patch assignment by VERTEX-TO-GRID-CORNER agreement,
    // globally greedy over sorted scores. Nearest-centroid matching
    // in face-iteration order was measured order-dependent at n = 8:
    // vertex-average centroids of the small skewed quads cross over
    // between adjacent pieces, a face steals its neighbour's patch,
    // and every seam downstream tears by up to the patch size
    // (octagon8_mixed side-vs-side 3e-1). A face's outer-loop
    // vertices sit at its own patch's grid corners to ~1e-4, but at
    // least one vertex is a full patch width from any WRONG patch's
    // corner set — max-over-vertices of min-corner-distance
    // separates the classes by orders of magnitude.
    let mut cap_faces = Vec::new();
    let face_handles: Vec<FaceHandle> = brep.faces.iter().map(|(fh, _)| fh).collect();
    let mut pairs: Vec<(f64, usize, usize)> = Vec::new();
    let mut face_vertex_sets: Vec<Vec<Point3>> = Vec::new();
    for fh in &face_handles {
        let vertices: Vec<Point3> = match face_loop_half_edges(brep, *fh) {
            Ok(hes) => hes
                .iter()
                .filter_map(|&he| {
                    brep.vertices
                        .get(brep.half_edges[he].origin)
                        .map(|v| v.position)
                })
                .collect(),
            Err(_) => Vec::new(),
        };
        face_vertex_sets.push(vertices);
    }
    for (face_idx, vertices) in face_vertex_sets.iter().enumerate() {
        if vertices.is_empty() {
            continue;
        }
        for (patch_idx, patch) in patches.iter().enumerate() {
            let popkern_core::geometry::gregory::GregoryPatch::RectangularBiNonic(bp) = patch
            else {
                continue;
            };
            let grid_corners = [bp.grid[0][0], bp.grid[0][9], bp.grid[9][0], bp.grid[9][9]];
            let score = vertices
                .iter()
                .map(|v| {
                    grid_corners
                        .iter()
                        .map(|c| v.distance_to(c))
                        .fold(f64::INFINITY, f64::min)
                })
                .fold(0.0_f64, f64::max);
            if score < 0.05 {
                pairs.push((score, face_idx, patch_idx));
            }
        }
    }
    pairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    if popkern_core::trace::trace_enabled("setback") {
        let mut by_face: std::collections::BTreeMap<usize, Vec<(f64, usize)>> =
            std::collections::BTreeMap::new();
        for &(score, face_idx, patch_idx) in &pairs {
            by_face
                .entry(face_idx)
                .or_default()
                .push((score, patch_idx));
        }
        for (face_idx, mut entries) in by_face {
            entries.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
            let head: Vec<String> = entries
                .iter()
                .take(3)
                .map(|(s, p)| format!("patch{p}@{s:.3e}"))
                .collect();
            eprintln!(
                "  [assign] face {:?}: {}",
                face_handles[face_idx],
                head.join(" ")
            );
        }
    }
    let mut face_taken = vec![false; face_handles.len()];
    let mut patch_taken = vec![false; patches.len()];
    for (_, face_idx, patch_idx) in pairs {
        if face_taken[face_idx] || patch_taken[patch_idx] {
            continue;
        }
        if let Some(face) = brep.faces.get_mut(face_handles[face_idx]) {
            face.surface = Some(popkern_core::geometry::SurfaceOnFace {
                geometry: popkern_core::geometry::SurfaceGeometry::Gregory(
                    patches[patch_idx].clone(),
                ),
                orientation: true,
            });
            face_taken[face_idx] = true;
            patch_taken[patch_idx] = true;
            cap_faces.push(face_handles[face_idx]);
        }
    }
    if cap_faces.len() != 4 * n {
        return Err(format!(
            "assigned {} cap faces, expected {}",
            cap_faces.len(),
            4 * n
        ));
    }
    // The seam curves were wired from the pre-reconcile assembly rows;
    // refresh them from the committed (reconciled + re-unified) patch
    // boundary rows (see setback_fan_seam_refresh docs).
    let (refreshed, worst_moved) =
        setback_fan_seam_refresh::refresh_fan_seam_curves(brep, &cap_faces);
    if popkern_core::trace::trace_enabled("setback") {
        eprintln!(
            "  [seam-refresh] {refreshed} seam curves refreshed, worst move {worst_moved:.3e}"
        );
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
            // projections of the 3D curves. Cap pieces get the same
            // projected pcurves on their INTERNAL seams: the seam
            // correspondence between two pieces' boundary rows is a
            // Bezier sub-range reparameterization (non-affine), so an
            // endpoint-Linear pcurve skews same-parameter evaluation
            // by up to half a row (measured 1.4e-2 on the
            // hexagon6_variable fan while both sides trace the seam
            // within 2e-4).
            let cap_internal = !is_sweep && {
                let twin = brep.half_edges[he].twin;
                let twin_face = brep
                    .half_edges
                    .get(twin)
                    .and_then(|t| brep.loops.get(t.loop_))
                    .map(|l| l.face);
                twin_face.is_some_and(|f| cap_faces.contains(&f) && f != fh)
            };
            // Cap-internal seams between committed pieces are CHART
            // EDGES of the face's own patch: after the seam refresh
            // the edge curve IS the patch's boundary-row Bezier, whose
            // parameter equals the chart parameter, so the exact
            // pcurve is the analytic Linear along that chart edge —
            // no point inversion. Projection was measured diverging
            // off-domain on the wide N=3 sectors (uv u=1.409 while the
            // endpoints sit on grid corners to 1e-15), and the seam
            // postcondition read the wrong-basin plateau as
            // g0 2.3e-2 / g1 1.0.
            if cap_internal
                && let Some(popkern_core::geometry::SurfaceGeometry::Gregory(
                    popkern_core::geometry::gregory::GregoryPatch::RectangularBiNonic(bp),
                )) = brep
                    .faces
                    .get(fh)
                    .and_then(|f| f.surface.as_ref())
                    .map(|s| &s.geometry)
            {
                let corners = [
                    ((0.0, 0.0), bp.grid[0][0]),
                    ((0.0, 1.0), bp.grid[0][9]),
                    ((1.0, 0.0), bp.grid[9][0]),
                    ((1.0, 1.0), bp.grid[9][9]),
                ];
                let match_corner = |p: popkern_core::geometry::Point3| {
                    corners
                        .iter()
                        .map(|(uv, c)| (*uv, p.distance_to(c)))
                        .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
                        .filter(|(_, d)| *d <= 1.0e-3)
                        .map(|(uv, _)| uv)
                };
                let he_data = &brep.half_edges[he];
                let origin_v = he_data.origin;
                let dest_v = brep.half_edges[he_data.twin].origin;
                let (Some(origin_p), Some(dest_p)) = (
                    brep.vertices.get(origin_v).map(|v| v.position),
                    brep.vertices.get(dest_v).map(|v| v.position),
                ) else {
                    continue;
                };
                let matched = (match_corner(origin_p), match_corner(dest_p));
                if popkern_core::trace::trace_enabled("setback")
                    && (matched.0.is_none() || matched.1.is_none())
                {
                    eprintln!(
                        "  [iso-pcurve] face {fh:?} he {he:?} off-corner: origin d {:.2e}, dest d {:.2e}",
                        corners
                            .iter()
                            .map(|(_, c)| origin_p.distance_to(c))
                            .fold(f64::INFINITY, f64::min),
                        corners
                            .iter()
                            .map(|(_, c)| dest_p.distance_to(c))
                            .fold(f64::INFINITY, f64::min),
                    );
                }
                if let (Some(uv_o), Some(uv_d)) = matched
                    && uv_o != uv_d
                    // Shared chart edge: exactly one coordinate differs.
                    && ((uv_o.0 == uv_d.0) != (uv_o.1 == uv_d.1))
                {
                    if let Some(he_mut) = brep.half_edges.get_mut(he) {
                        he_mut.pcurve = Some(Pcurve::linear(uv_o, uv_d, 0.0, 1.0));
                    }
                    continue;
                }
            }
            if is_sweep || cap_internal {
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
                let mut samples3d = Vec::with_capacity(N);
                let mut ok = true;
                for i in 0..N {
                    let f = i as f64 / (N - 1) as f64;
                    let f = if forward { f } else { 1.0 - f };
                    let t = curve.t_start + (curve.t_end - curve.t_start) * f;
                    let p3 = ceval.point_at(t);
                    match project_point_to_uv(&surface, &p3) {
                        Some((u, v)) => {
                            cps.push(popkern_core::geometry::Point3::new(u, v, 0.0));
                            samples3d.push(p3);
                        }
                        None => {
                            ok = false;
                            break;
                        }
                    }
                }
                if !ok {
                    continue;
                }
                setback_arc_binding::repair_diverged_projections(&mut cps, &samples3d, &surface);
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
    if popkern_core::trace::trace_enabled("setback") {
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

/// Find the face carrying this fillet's sweep surface (sample match).
/// Merge the surrendered faces (apex slivers, wedges, and everything
/// still touching the apex vertex) into a single blend face by
/// kef-ing shared edges.
fn merge_faces(
    brep: &mut BRep,
    vertex: VertexHandle,
    seed_faces: &[FaceHandle],
    ring_edges: &std::collections::BTreeSet<EdgeHandle>,
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
    // Flood-fill the ring interior: pre-carve slivers inside the
    // setback ring can touch neither the apex nor a seed face
    // (measured on hexagon6_variable — split-strip and mitre-fan
    // pieces left the region disconnected, "3 remain"). Any face
    // reachable from a seed without crossing a ring edge is part of
    // the surrendered region.
    loop {
        let mut grew = false;
        let frontier: Vec<FaceHandle> = region.iter().copied().collect();
        for fh in frontier {
            for he in face_loop_half_edges(brep, fh).unwrap_or_default() {
                let he_data = &brep.half_edges[he];
                if ring_edges.contains(&he_data.edge) {
                    continue;
                }
                let twin_face = brep.loops[brep.half_edges[he_data.twin].loop_].face;
                if region.insert(twin_face) {
                    grew = true;
                }
            }
        }
        if !grew {
            break;
        }
    }
    // Leak guard: the ring must separate the region from the rest of
    // the model — a ring edge with the region on BOTH sides means the
    // fill escaped (a ring segment crossed an un-split edge) and the
    // merge would eat live geometry.
    for &eh in ring_edges {
        let Some(edge) = brep.edges.get(eh) else {
            continue;
        };
        let he = edge.half_edge;
        let f1 = brep.loops[brep.half_edges[he].loop_].face;
        let f2 = brep.loops[brep.half_edges[brep.half_edges[he].twin].loop_].face;
        if region.contains(&f1) && region.contains(&f2) {
            if popkern_core::trace::trace_enabled("setback") {
                eprintln!(
                    "  [merge-leak] ring edge {eh:?} has region on both sides: f1={f1:?} f2={f2:?} \
                     (region size {})",
                    region.len()
                );
            }
            return Err("setback ring does not enclose the surrendered region".to_string());
        }
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
    //
    // Scoped to the surrendered region's own loops: a periodic seam
    // anywhere else on the body also has both half-edges in one loop
    // (measured: the boss cylinder's seam on the flagship boss plate,
    // untouched by any fillet, was dissolved here under a G2 request
    // and left the cylinder's loop open by the seam length, PK-389).
    loop {
        let mut dissolved = false;
        let candidates: Vec<HalfEdgeHandle> = brep
            .half_edges
            .iter()
            .filter(|(_, data)| {
                brep.loops
                    .get(data.loop_)
                    .is_some_and(|l| region.contains(&l.face))
                    && brep
                        .half_edges
                        .get(data.twin)
                        .map(|t| t.loop_ == data.loop_)
                        .unwrap_or(false)
            })
            .map(|(h, _)| h)
            .collect();
        if popkern_core::trace::trace_enabled("setback") {
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
                    if popkern_core::trace::trace_enabled("setback") {
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

pub(crate) use super::setback_vertex_loop_walk::face_loop_half_edges;

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
