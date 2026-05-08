use super::var_radius_geometry::sample_pcurve_on_surface;
use super::*;

// ================================================================
// Vertex-split primitive — un-merges a shared rolling-ball
// footprint vertex into N new vertices, one per incident
// cylinder fillet.
// ================================================================

/// Split a shared apex-region vertex `v` so that `v_keep` retains
/// incidence to `cyl_keep`, and a new vertex `v_new` (at
/// `v_new_pos`) inherits the half-edges incident to `cyl_split`.
///
/// Topology before:
///   v is incident to: cap (1 he-pair), cube_face (1 he-pair),
///   cyl_keep (1 he-pair), cyl_split (1 he-pair).
///
/// Topology after:
///   v continues to be incident to: cap, cube_face, cyl_keep.
///   v_new is incident to: cap, cube_face, cyl_split.
///   A new edge connects v and v_new, with cap on one side and
///   cube_face on the other (this is the "FaceEdge bridge"
///   from Phase 6c-A's substituted boundary).
///
/// Algorithm (raw half-edge surgery, snapshot-rollback):
/// 1. Identify the 4 outgoing half-edges from v (one per incident
///    face) and their incoming counterparts (`he.prev` in each
///    face's loop).
/// 2. Create new vertex `v_new` at `v_new_pos`.
/// 3. Re-route cyl_split's outgoing+incoming half-edges' origin
///    pointer to `v_new` (the half-edges stay in cyl_split's
///    loop; only the origin vertex changes).
/// 4. Create the FaceEdge bridge edge:
///    * he_bridge_cap: cap-side, V → V_new, inserted in cap face's
///      loop after the cap-side half-edge currently from v
///      towards cyl_split.
///    * he_bridge_cube: cube-side, V_new → V, twin of he_bridge_cap,
///      inserted in cube face's loop.
/// 5. Update `cap_face.outer_loop` and `cube_face.outer_loop`
///    half_edge pointers if they pointed to a re-routed he.
/// 6. Update vertex `outgoing` pointers for v and v_new.
///
/// Split a shared apex-region vertex `v` so that `v_keep` retains
/// incidence to `cap_face`, `cube_face`, and `cyl_keep`, while
/// a new vertex `v_new` (at `v_new_pos`) inherits the half-edges
/// incident to `cyl_split`.
///
/// A new edge connects `v` and `v_new`, with `cap_face` on one
/// side and `cube_face` on the other (the FaceEdge bridge from
/// Phase 6c-A's substituted boundary).
///
/// On any sub-step failure, snapshot-rollback restores `brep`.
#[allow(dead_code)]
pub(crate) fn split_apex_footprint_vertex(
    brep: &mut BRep,
    v: VertexHandle,
    cap_face: FaceHandle,
    cube_face: FaceHandle,
    cyl_split: FaceHandle,
    v_new_pos: popkern_core::geometry::Point3,
) -> Result<VertexHandle, RerouteError> {
    use popkern_core::topology::{Edge, HalfEdge, Vertex};

    let snap = brep.clone();

    let he_v_cap = match find_outgoing_he_in_face(brep, v, cap_face) {
        Some(h) => h,
        None => {
            *brep = snap;
            return Err(RerouteError::LoopBroken);
        }
    };
    let he_v_cube = match find_outgoing_he_in_face(brep, v, cube_face) {
        Some(h) => h,
        None => {
            *brep = snap;
            return Err(RerouteError::LoopBroken);
        }
    };
    let he_v_split = match find_outgoing_he_in_face(brep, v, cyl_split) {
        Some(h) => h,
        None => {
            *brep = snap;
            return Err(RerouteError::LoopBroken);
        }
    };

    // 1. Create v_new at V's CURRENT position (zero-length bridge
    //    initially). After topology surgery + validate passes, the
    //    caller can move v_new to v_new_pos via a separate vertex-
    //    move step that updates v_new's incident edge curves.
    let v_starting_pos = brep.vertices[v].position;
    let v_new = brep.vertices.insert(Vertex {
        position: v_starting_pos,
        outgoing: Default::default(),
        name: None,
        tolerance: None,
    });
    let _ = v_new_pos; // applied later via vertex move

    // 2. Re-point cyl_split's outgoing he origin to v_new.
    brep.half_edges[he_v_split].origin = v_new;

    // Detect which face's outgoing-from-V edge is the cyl_seam
    // (twin in cyl_split). Exactly one of {cap, cube} should be the
    // cyl_seam side for the bridge twin pair to work correctly.
    //
    // Case A — `cap_outgoing_seam == true`, `cube_outgoing_seam == false`:
    //   cap's loop walks ..., cap_v_prev (cap-tangent_in V_pre→V),
    //     bridge_cap (V→V_new), he_v_cap (V_new→cyl_split-side), ...
    //   cube's loop walks ..., cube_v_prev (cyl_seam V_pre→V_new),
    //     bridge_cube (V_new→V), he_v_cube (V→cap-tangent_out), ...
    //
    // Case B — swapped: `cap_outgoing_seam == false`, `cube_outgoing_seam == true`:
    //   cap's loop walks ..., cap_v_prev (cyl_seam-cap_n V_pre→V_new),
    //     bridge_cap (V_new→V), he_v_cap (V→cap-tangent_out), ...
    //   cube's loop walks ..., cube_v_prev (cap-tangent_in V_pre→V),
    //     bridge_cube (V→V_new), he_v_cube (V_new→cyl_seam_out), ...
    let cap_outgoing_seam =
        brep.loops[brep.half_edges[brep.half_edges[he_v_cap].twin].loop_].face == cyl_split;
    let cube_outgoing_seam =
        brep.loops[brep.half_edges[brep.half_edges[he_v_cube].twin].loop_].face == cyl_split;
    let trace_split = std::env::var("POPKERN_VAR_RADIUS_TRACE").is_ok_and(|v| v == "1");
    if trace_split {
        eprintln!(
            "[split] cap_outgoing_seam={cap_outgoing_seam} cube_outgoing_seam={cube_outgoing_seam}"
        );
    }
    // Case A: cap_outgoing_seam=true, cube_outgoing_seam=false (original design).
    // Case B: cap_outgoing_seam=false, cube_outgoing_seam=true (swapped).
    // Anything else: configuration not handled by current surgery.
    let case_b = match (cap_outgoing_seam, cube_outgoing_seam) {
        (true, false) => false,
        (false, true) => true,
        _ => {
            *brep = snap;
            return Err(RerouteError::LoopBroken);
        }
    };

    // 3. Create the bridge edge.
    let bridge_edge = brep.edges.insert(Edge {
        half_edge: Default::default(),
        curve: None,
        name: None,
        tolerance: None,
        canonical: Default::default(),
    });

    let cap_loop = brep.half_edges[he_v_cap].loop_;
    let cube_loop = brep.half_edges[he_v_cube].loop_;
    let cap_v_prev = brep.half_edges[he_v_cap].prev;
    let cube_v_prev = brep.half_edges[he_v_cube].prev;

    // Case A: bridge_cap (V→V_new) in cap_loop; bridge_cube (V_new→V) in cube_loop.
    // Case B: bridge_cap (V_new→V) in cap_loop; bridge_cube (V→V_new) in cube_loop.
    let bridge_cap_origin = if case_b { v_new } else { v };
    let bridge_cube_origin = if case_b { v } else { v_new };
    let bridge_cap = brep.half_edges.insert(HalfEdge {
        origin: bridge_cap_origin,
        twin: Default::default(),
        next: he_v_cap,
        prev: cap_v_prev,
        edge: bridge_edge,
        loop_: cap_loop,
        pcurve: None,
    });
    let bridge_cube = brep.half_edges.insert(HalfEdge {
        origin: bridge_cube_origin,
        twin: bridge_cap,
        next: he_v_cube,
        prev: cube_v_prev,
        edge: bridge_edge,
        loop_: cube_loop,
        pcurve: None,
    });
    brep.half_edges[bridge_cap].twin = bridge_cube;
    brep.edges[bridge_edge].half_edge = bridge_cap;

    // 4. Splice into loops.
    brep.half_edges[cap_v_prev].next = bridge_cap;
    brep.half_edges[he_v_cap].prev = bridge_cap;
    brep.half_edges[cube_v_prev].next = bridge_cube;
    brep.half_edges[he_v_cube].prev = bridge_cube;

    // 5. Re-point the outgoing-from-V HE that is on the cyl_seam side.
    //    Case A: he_v_cap is the cyl_seam (cap-cyl_split edge) → re-point.
    //    Case B: he_v_cube is the cyl_seam (cube-cyl_split edge) → re-point.
    if case_b {
        brep.half_edges[he_v_cube].origin = v_new;
    } else {
        brep.half_edges[he_v_cap].origin = v_new;
    }

    // 6. Bridge curve: degenerate (zero-length) initially, both
    //    endpoints at V's position. Vertex-move step (separate)
    //    will update CP[1] to v_new_pos.
    let v_pos = brep.vertices[v].position;
    let bridge_curve = popkern_core::geometry::nurbs_curve::NurbsCurve {
        degree: 1,
        control_points: vec![v_pos, v_pos],
        weights: vec![1.0, 1.0],
        knots: vec![0.0, 0.0, 1.0, 1.0],
    };
    brep.edges[bridge_edge].curve = Some(popkern_core::geometry::curve::CurveOnEdge {
        geometry: popkern_core::geometry::CurveGeometry::Nurbs(bridge_curve),
        t_start: 0.0,
        t_end: 1.0,
    });

    // 7. Update vertex outgoing pointers.  v's outgoing must be a
    //    half-edge with origin == v; same for v_new.  In Case A,
    //    bridge_cap.origin = v and bridge_cube.origin = v_new.
    //    In Case B, the roles are swapped.
    let v_outgoing = if bridge_cap_origin == v {
        bridge_cap
    } else {
        bridge_cube
    };
    let v_new_outgoing = if bridge_cap_origin == v_new {
        bridge_cap
    } else {
        bridge_cube
    };
    brep.vertices[v].outgoing = v_outgoing;
    brep.vertices[v_new].outgoing = v_new_outgoing;
    let _ = he_v_split; // also valid as v_new outgoing; bridge variant chosen for consistency

    // 7b. Re-point twin origins on the cyl_seam edge between
    //     cyl_split and the "incoming-V" face. The incoming-V face
    //     is whichever of {cap, cube} has its INCOMING-V edge as the
    //     cyl_seam (twin in cyl_split):
    //
    //     * Case A (cap outgoing IS cyl_seam → cube_v_prev IS cyl_seam):
    //       cube_v_prev_twin lives in cyl_split (= he_v_split.prev or
    //       he_v_split itself depending on direction). Its origin
    //       is V → re-point to V_new.
    //
    //     * Case B (cube outgoing IS cyl_seam → cap_v_prev IS cyl_seam):
    //       cap_v_prev_twin lives in cyl_split. Its origin is V →
    //       re-point to V_new.
    //
    //     The "incoming-V" face's V-prev twin in cyl_split is
    //     he_v_split.prev — already covered if its origin is V (the
    //     classic step 7b second re-point). After step 5, he_v_split
    //     itself has been re-pointed; he_v_split.prev wasn't, so its
    //     twin still has origin V (= correct to re-point now).
    let cube_v_prev_twin = brep.half_edges[cube_v_prev].twin;
    let cube_v_prev_twin_face = brep.loops[brep.half_edges[cube_v_prev_twin].loop_].face;
    let cap_v_prev_twin = brep.half_edges[cap_v_prev].twin;
    let cap_v_prev_twin_face = brep.loops[brep.half_edges[cap_v_prev_twin].loop_].face;
    if trace_split {
        eprintln!(
            "[split] case_b={case_b} he_v_cap={he_v_cap:?} he_v_cube={he_v_cube:?} he_v_split={he_v_split:?} cap_v_prev={cap_v_prev:?} (twin {cap_v_prev_twin:?} face {cap_v_prev_twin_face:?}) cube_v_prev={cube_v_prev:?} (twin {cube_v_prev_twin:?} face {cube_v_prev_twin_face:?}) (cap={cap_face:?} cube={cube_face:?} cyl_split={cyl_split:?})"
        );
    }
    // Case A: cube_v_prev IS the cyl_seam (twin in cyl_split). Re-point.
    // Case B: cap_v_prev IS the cyl_seam. Re-point.
    let prev_seam_twin = if case_b {
        cap_v_prev_twin
    } else {
        cube_v_prev_twin
    };
    if brep.half_edges[prev_seam_twin].origin == v {
        brep.half_edges[prev_seam_twin].origin = v_new;
    }
    let he_v_split_in = brep.half_edges[he_v_split].prev;
    let he_v_split_in_twin = brep.half_edges[he_v_split_in].twin;
    let he_v_split_in_twin_face = brep.loops[brep.half_edges[he_v_split_in_twin].loop_].face;
    if trace_split {
        eprintln!(
            "[split] he_v_split_in={he_v_split_in:?} he_v_split_in_twin={he_v_split_in_twin:?} (face {he_v_split_in_twin_face:?})"
        );
    }
    if brep.half_edges[he_v_split_in_twin].origin == v {
        brep.half_edges[he_v_split_in_twin].origin = v_new;
    }

    // 8. Set bridge edge pcurves (degree-1 chord on cube_face's
    //    plane and cap_face's surface). Skip if either projection
    //    fails — validate will catch issues.
    let cap_surf = brep
        .faces
        .get(cap_face)
        .and_then(|f| f.surface.as_ref())
        .map(|s| s.geometry.clone());
    let cube_surf = brep
        .faces
        .get(cube_face)
        .and_then(|f| f.surface.as_ref())
        .map(|s| s.geometry.clone());
    if let Some(surf) = cap_surf
        && let Some(bc) = brep.edges[bridge_edge].curve.as_ref()
        && let popkern_core::geometry::CurveGeometry::Nurbs(curve_3d) = &bc.geometry
        && let Some(pc) = sample_pcurve_on_surface(curve_3d, &surf)
    {
        brep.half_edges[bridge_cap].pcurve = Some(pc.into());
    }
    if let Some(surf) = cube_surf
        && let Some(bc) = brep.edges[bridge_edge].curve.as_ref()
        && let popkern_core::geometry::CurveGeometry::Nurbs(curve_3d) = &bc.geometry
        && let Some(pc) = sample_pcurve_on_surface(curve_3d, &surf)
    {
        // Bridge_cube goes V_new -> V, so invert direction.
        let mut inv = pc.clone();
        inv.control_points.reverse();
        inv.weights.reverse();
        brep.half_edges[bridge_cube].pcurve = Some(inv.into());
    }

    // 9. Validate.
    if let Err(e) = brep.validate() {
        if std::env::var("POPKERN_VAR_RADIUS_TRACE").is_ok_and(|v| v == "1") {
            eprintln!("[split] validate failed: {e:?}");
        }
        *brep = snap;
        return Err(RerouteError::ValidationFailed);
    }

    Ok(v_new)
}

/// Find the (unique) Plane-surface face incident to `v`. Returns
/// `None` if zero or multiple plane faces are incident.
pub(super) fn find_incident_plane_face(brep: &BRep, v: VertexHandle) -> Option<FaceHandle> {
    use popkern_core::geometry::surface::SurfaceGeometry;
    let mut found: Option<FaceHandle> = None;
    for he in brep.half_edges.values() {
        if he.origin != v {
            continue;
        }
        let f = brep.loops.get(he.loop_)?.face;
        let surf = brep
            .faces
            .get(f)
            .and_then(|fd| fd.surface.as_ref())
            .map(|s| &s.geometry);
        if matches!(surf, Some(SurfaceGeometry::Plane { .. })) {
            if found.is_some() && found != Some(f) {
                return None;
            }
            found = Some(f);
        }
    }
    found
}

/// Snap whichever of `edge`'s curve's first/last CP is at
/// `from_pos` to `to_pos`.  No-op if neither endpoint matches.
#[allow(dead_code)]
fn snap_edge_curve_endpoint(
    brep: &mut BRep,
    edge: EdgeHandle,
    from_pos: popkern_core::geometry::Point3,
    to_pos: popkern_core::geometry::Point3,
) {
    use popkern_core::geometry::CurveGeometry;
    let tol = 1e-6;
    if let Some(coe) = brep.edges.get_mut(edge).and_then(|e| e.curve.as_mut())
        && let CurveGeometry::Nurbs(nc) = &mut coe.geometry
        && !nc.control_points.is_empty()
    {
        let n = nc.control_points.len();
        let first = nc.control_points[0];
        let last = nc.control_points[n - 1];
        if first.distance_to(&from_pos) < tol {
            nc.control_points[0] = to_pos;
        } else if last.distance_to(&from_pos) < tol {
            nc.control_points[n - 1] = to_pos;
        }
    }
}

/// Find the unique outgoing half-edge from `v` that lies in
/// `face`'s outer loop. Returns `None` if not found or if `v`
/// is multiply incident to `face`'s loop (non-manifold).
fn find_outgoing_he_in_face(
    brep: &BRep,
    v: VertexHandle,
    face: FaceHandle,
) -> Option<HalfEdgeHandle> {
    let face_data = brep.faces.get(face)?;
    let outer_loop = face_data.outer_loop;
    let start = brep.loops.get(outer_loop)?.half_edge;
    let mut cursor = start;
    let mut result = None;
    for _ in 0..32 {
        let he = brep.half_edges.get(cursor)?;
        if he.origin == v {
            if result.is_some() {
                return None; // ambiguous
            }
            result = Some(cursor);
        }
        cursor = he.next;
        if cursor == start {
            break;
        }
    }
    result
}
