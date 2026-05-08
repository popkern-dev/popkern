use super::var_radius_geometry::{project_to_sphere_cylinder_intersection, set_cut_edge_geometry};
use super::var_radius_path_b_vertex_move::update_incident_edges_for_vertex_move;
use super::var_radius_vertex_split::find_incident_plane_face;
use super::*;

// ================================================================
// Path B — `mef`-split + face-absorb commercial-grade surgery
// ================================================================

/// Plan for one `mef`-split + absorb surgery on a NeighbourCylinder
/// edge.  Built before any mutation; `apply_path_b_plan` consumes
/// it transactionally with snapshot-rollback.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub(crate) struct PathBPlan {
    /// The NeighbourCylinder edge in cap_k's loop (cap-side
    /// half-edge).
    pub(crate) cap_k_he: HalfEdgeHandle,
    /// The cap face that owns the NeighCyl edge.
    pub(crate) cap_k: FaceHandle,
    /// The neighbour cylinder face whose loop contains the twin.
    pub(crate) cyl_m: FaceHandle,
    /// The destination cap (whose primary cylinder is cyl_m).
    pub(crate) cap_n: FaceHandle,
    /// V_x = origin of cap_k_he. To be moved to `v_x_new_pos`.
    pub(crate) v_x: VertexHandle,
    /// V_y = destination of cap_k_he (origin of twin). 3-way junction.
    pub(crate) v_y: VertexHandle,
    /// New target position for V_x — on cap_n's sphere AND cube
    /// face AND cyl_m surface.
    pub(crate) v_x_new_pos: popkern_core::geometry::Point3,
    /// The cap_n OwnCyl edge in cyl_m's loop that comes "before"
    /// the apex chip on the cap_n side. Its origin vertex is the
    /// other endpoint of the cap-tangent boundary on cube_face_y0.
    pub(crate) cap_n_own_he_a: HalfEdgeHandle,
    /// New target position for cap_n_own_he_a's origin vertex.
    pub(crate) v_a_new_pos: popkern_core::geometry::Point3,
}

/// Build a Path B plan for one NeighbourCylinder edge.
///
/// Returns `None` when the topology doesn't match the expected
/// pattern (validates the assumptions before mutation).
///
/// `per_fillet_caps` is parallel to `cap_faces` — entry `i`
/// gives the per-fillet sphere `(centre, radius)` for cap face
/// `cap_faces[i]`.  This is read from `BlendBoundary` rather than
/// from the cap face's surface so the surgery works for both the
/// sphere-cap rung (cap surface IS the per-fillet sphere) and the
/// Phase 6c-D Gregory cap rung (cap surface is a Gregory patch
/// but its nominal per-fillet sphere is still meaningful for
/// finding cap-n by primary cylinder match).
#[allow(dead_code)]
pub(crate) fn plan_path_b(
    brep: &BRep,
    cap_k_he: HalfEdgeHandle,
    cap_k: FaceHandle,
    cap_faces: &[FaceHandle],
    per_fillet_caps: &[(popkern_core::geometry::Point3, f64)],
) -> Option<PathBPlan> {
    use popkern_core::geometry::surface::SurfaceGeometry;

    // Extract cap_k's NeighCyl edge: V_x, V_y, twin face cyl_m.
    let tw = brep.half_edges[cap_k_he].twin;
    let v_x = brep.half_edges[cap_k_he].origin;
    let v_y = brep.half_edges[tw].origin;
    let cyl_m = brep.loops[brep.half_edges[tw].loop_].face;

    // cyl_m must be a Cylinder.
    let (cyl_axis, cyl_axis_origin, cyl_radius) = match brep
        .faces
        .get(cyl_m)
        .and_then(|f| f.surface.as_ref())
        .map(|s| &s.geometry)
    {
        Some(SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        }) => (*axis, *origin, *radius),
        _ => return None,
    };

    // dest cap_n: cap whose per-fillet sphere matches `cyl_m` —
    // radius matches `cyl_radius` and centre projects to `cyl_m`'s
    // axis.  Look up via the parallel `per_fillet_caps` array
    // (works for any cap surface type, not just Sphere).
    if cap_faces.len() != per_fillet_caps.len() {
        return None;
    }
    let cap_n_idx = cap_faces.iter().enumerate().find_map(|(i, &fh)| {
        if fh == cap_k {
            return None;
        }
        let (center, radius) = per_fillet_caps[i];
        if (radius - cyl_radius).abs() > cyl_radius * 0.05 {
            return None;
        }
        let to = popkern_core::geometry::Vec3::new(
            center.x - cyl_axis_origin.x,
            center.y - cyl_axis_origin.y,
            center.z - cyl_axis_origin.z,
        );
        let axis_unit = cyl_axis.normalized().unwrap_or(cyl_axis);
        let along = to.dot(&axis_unit);
        let proj = popkern_core::geometry::Point3::new(
            cyl_axis_origin.x + axis_unit.x * along,
            cyl_axis_origin.y + axis_unit.y * along,
            cyl_axis_origin.z + axis_unit.z * along,
        );
        if center.distance_to(&proj) < cyl_radius * 0.05 {
            Some(i)
        } else {
            None
        }
    })?;
    let cap_n = cap_faces[cap_n_idx];
    let (sphere_centre, sphere_radius) = per_fillet_caps[cap_n_idx];

    // Walk cyl_m's outer loop to find the "cap_n OwnCyl edge a" —
    // the cap_n OwnCyl that's NOT shared with V_y (the one whose
    // origin is the cube-face cap-tangent endpoint on the OTHER
    // cube face).
    let cap_set: std::collections::BTreeSet<FaceHandle> = cap_faces.iter().copied().collect();
    let cyl_m_loop = brep.faces.get(cyl_m)?.outer_loop;
    let start = brep.loops.get(cyl_m_loop)?.half_edge;
    let mut cap_n_own_he_a: Option<HalfEdgeHandle> = None;
    let mut cursor = start;
    for _ in 0..32 {
        let he_origin = brep.half_edges[cursor].origin;
        let he_twin = brep.half_edges[cursor].twin;
        let he_next = brep.half_edges[cursor].next;
        let twin_face = brep.loops[brep.half_edges[he_twin].loop_].face;
        // Looking for cap_n's OwnCyl edge whose origin is NOT V_y.
        if twin_face == cap_n && he_origin != v_y && cap_set.contains(&twin_face) {
            cap_n_own_he_a = Some(cursor);
            break;
        }
        cursor = he_next;
        if cursor == start {
            break;
        }
    }
    let cap_n_own_he_a = cap_n_own_he_a?;
    let v_a = brep.half_edges[cap_n_own_he_a].origin;

    // Compute target positions: V_x_new and V_a_new lie at the
    // intersection of (cap_n's sphere) ∩ cyl_m, snapped to the
    // cube face containing the original V_x / V_a.  Mathematically:
    //
    //   V_new = current_V projected onto sphere_n ∩ cyl_m
    //
    // Both surfaces pass through the new position by construction.
    // The cube face plane that V_x or V_a lies on is preserved
    // automatically because the sphere ∩ cyl intersection circle
    // crosses each cube face at exactly one point.
    let v_x_pos = brep.vertices.get(v_x)?.position;
    let v_a_pos = brep.vertices.get(v_a)?.position;
    let v_x_new_pos = project_to_sphere_cylinder_intersection(
        v_x_pos,
        sphere_centre,
        sphere_radius,
        cyl_axis_origin,
        cyl_axis,
    )?;
    let v_a_new_pos = project_to_sphere_cylinder_intersection(
        v_a_pos,
        sphere_centre,
        sphere_radius,
        cyl_axis_origin,
        cyl_axis,
    )?;
    if std::env::var("POPKERN_VAR_RADIUS_TRACE").is_ok_and(|v| v == "1") {
        eprintln!(
            "[plan] cap_k={cap_k:?} cyl_m={cyl_m:?} (axis=({:.3},{:.3},{:.3}) r={cyl_radius:.3}) cap_n={cap_n:?} sphere=({:.3},{:.3},{:.3})",
            cyl_axis.x, cyl_axis.y, cyl_axis.z, sphere_centre.x, sphere_centre.y, sphere_centre.z
        );
        eprintln!(
            "  V_x={v_x:?}=({:.3},{:.3},{:.3}) → ({:.3},{:.3},{:.3})",
            v_x_pos.x, v_x_pos.y, v_x_pos.z, v_x_new_pos.x, v_x_new_pos.y, v_x_new_pos.z
        );
        eprintln!(
            "  V_a={v_a:?}=({:.3},{:.3},{:.3}) → ({:.3},{:.3},{:.3})",
            v_a_pos.x, v_a_pos.y, v_a_pos.z, v_a_new_pos.x, v_a_new_pos.y, v_a_new_pos.z
        );
    }

    Some(PathBPlan {
        cap_k_he,
        cap_k,
        cyl_m,
        cap_n,
        v_x,
        v_y,
        v_x_new_pos,
        cap_n_own_he_a,
        v_a_new_pos,
    })
}

/// Apply a `PathBPlan` to the BRep transactionally.
///
/// Steps:
/// 1. Move V_x and V_a to their new positions on the cap-tangent
///    circle (sphere_n ∩ cyl_m intersection in the plane through
///    sphere_n's centre perpendicular to cyl_m's axis).
/// 2. `mef` cyl_m's outer loop: insert the cut edge from V_x to V_a.
///    cyl_m face slot now holds the apex chip; mef returns a new
///    face for the main cylinder body.
/// 3. Absorb the apex chip into cap_n by splicing cap_n's loop
///    around the chip's 2 shared edges (cap_n's OwnCyl edges with
///    cyl_m), replacing them with the chip's other 2 edges
///    (cap[0]'s NeighCyl edge + the new cut edge).
/// 4. Delete the now-orphaned chip face, lens loop, shared-edge
///    half-edges, and the V_p5 mid-vertex.
/// 5. `validate()`; rollback on any failure.
#[allow(dead_code)]
pub(crate) fn apply_path_b_plan(brep: &mut BRep, plan: &PathBPlan) -> Result<(), RerouteError> {
    let snapshot = brep.clone();
    let v_a_handle = brep.half_edges[plan.cap_n_own_he_a].origin;
    let v_x_old_pos = brep.vertices[plan.v_x].position;
    let v_a_old_pos = brep.vertices[v_a_handle].position;

    // Step 0: split V_x and V_a so each new vertex is incident
    // to ONLY cyl_m (not cyl_keep). This avoids conflicts when
    // multiple surgeries share these footprint vertices.
    let trace = std::env::var("POPKERN_VAR_RADIUS_TRACE").is_ok_and(|v| v == "1");
    let cube_face_x = match find_incident_plane_face(brep, plan.v_x) {
        Some(f) => f,
        None => {
            if trace {
                eprintln!("[apply] no plane face incident to V_x {:?}", plan.v_x);
            }
            *brep = snapshot;
            return Err(RerouteError::LoopBroken);
        }
    };
    let cube_face_a = match find_incident_plane_face(brep, v_a_handle) {
        Some(f) => f,
        None => {
            if trace {
                eprintln!("[apply] no plane face incident to V_a {v_a_handle:?}");
            }
            *brep = snapshot;
            return Err(RerouteError::LoopBroken);
        }
    };
    if trace {
        eprintln!(
            "[apply] split V_x={:?} (cube={cube_face_x:?}) V_a={v_a_handle:?} (cube={cube_face_a:?})",
            plan.v_x
        );
    }
    let _ = (v_x_old_pos, v_a_old_pos);
    let v_x_new = match split_apex_footprint_vertex(
        brep,
        plan.v_x,
        plan.cap_k,
        cube_face_x,
        plan.cyl_m,
        plan.v_x_new_pos,
    ) {
        Ok(v) => v,
        Err(e) => {
            if trace {
                eprintln!("[apply] split V_x failed: {e:?}");
            }
            *brep = snapshot;
            return Err(e);
        }
    };
    let v_a_new = match split_apex_footprint_vertex(
        brep,
        v_a_handle,
        plan.cap_n,
        cube_face_a,
        plan.cyl_m,
        plan.v_a_new_pos,
    ) {
        Ok(v) => v,
        Err(e) => {
            if trace {
                eprintln!("[apply] split V_a failed: {e:?}");
            }
            *brep = snapshot;
            return Err(e);
        }
    };
    let _ = (v_x_new, v_a_new);

    // Step 2: mef cyl_m. After splits, cyl_m's loop has v_x_new
    // (where v_x used to be) and v_a_new (where v_a used to be).
    // he_at_v_x is the half-edge in cyl_m's loop whose origin
    // is now v_x_new — found via `cap_k_tw.next`.
    let cap_k_tw = brep.half_edges[plan.cap_k_he].twin;
    let he_at_v_x = brep.half_edges[cap_k_tw].next;
    if brep.half_edges[he_at_v_x].origin != v_x_new {
        if trace {
            eprintln!(
                "[apply] he_at_v_x origin {:?} != v_x_new {v_x_new:?}",
                brep.half_edges[he_at_v_x].origin
            );
        }
        *brep = snapshot;
        return Err(RerouteError::LoopBroken);
    }

    let mef = match brep.mef(he_at_v_x, plan.cap_n_own_he_a) {
        Ok(r) => r,
        Err(_) => {
            *brep = snapshot;
            return Err(RerouteError::LoopBroken);
        }
    };

    // Step 3: chip face = original cyl_m face (mef keeps original
    // face slot for the loop containing he2). Identify chip's 4
    // half-edges and the 2 shared with cap_n.
    let chip_face = plan.cyl_m;
    let h1 = cap_k_tw; // V_y → V_x (NeighCyl twin, in chip loop now)
    let h2 = mef.half_edge; // V_x → V_a (cut edge, chip side)
    let h3 = plan.cap_n_own_he_a; // V_a → V_p5 (chip side of shared edge 1)
    let h4 = brep.half_edges[h3].next; // V_p5 → V_y (chip side of shared edge 2)

    // Sanity: walk chip loop and confirm 4 edges.
    {
        let chip_loop = brep.faces[chip_face].outer_loop;
        let start = brep.loops[chip_loop].half_edge;
        let mut cnt = 0;
        let mut cursor = start;
        for _ in 0..8 {
            cnt += 1;
            cursor = brep.half_edges[cursor].next;
            if cursor == start {
                break;
            }
        }
        if cnt != 4 {
            *brep = snapshot;
            return Err(RerouteError::LoopBroken);
        }
    }

    // Twin half-edges (in cap_n's loop).
    let t3 = brep.half_edges[h3].twin;
    let t4 = brep.half_edges[h4].twin;

    // Sanity: t3 and t4 must be in cap_n's loop.
    let cap_n_loop = brep.faces[plan.cap_n].outer_loop;
    if brep.half_edges[t3].loop_ != cap_n_loop || brep.half_edges[t4].loop_ != cap_n_loop {
        *brep = snapshot;
        return Err(RerouteError::LoopBroken);
    }

    // V_p5 = origin of h4 ("V_p5 → V_y").
    let v_p5 = brep.half_edges[h4].origin;

    // Find cap_n's neighbours of t3 and t4.
    // Pattern: cap_n's loop walks t3 → s1 → s2 → s3 → t4 → wrap → t3.
    // After absorb, s3.next becomes h1 → h2 → s1 (where s1 = t3.next, s3 = t4.prev).
    let s1 = brep.half_edges[t3].next;
    let s3 = brep.half_edges[t4].prev;

    // Step 4: splice cap_n's loop.
    // s3.next = h1; h1.prev = s3; h1.next = h2 (already); h2.prev = h1 (already);
    // h2.next = s1; s1.prev = h2; h1.loop_ = cap_n_loop; h2.loop_ = cap_n_loop.
    brep.half_edges[s3].next = h1;
    brep.half_edges[h1].prev = s3;
    brep.half_edges[h1].next = h2;
    brep.half_edges[h2].prev = h1;
    brep.half_edges[h2].next = s1;
    brep.half_edges[s1].prev = h2;
    brep.half_edges[h1].loop_ = cap_n_loop;
    brep.half_edges[h2].loop_ = cap_n_loop;

    // If cap_n's loop anchor pointed to t3 or t4, reassign it.
    if brep.loops[cap_n_loop].half_edge == t3 || brep.loops[cap_n_loop].half_edge == t4 {
        brep.loops[cap_n_loop].half_edge = h1;
    }

    // Step 5: delete chip face, chip loop, h3, h4, t3, t4, edges, v_p5.
    let h3_edge = brep.half_edges[h3].edge;
    let h4_edge = brep.half_edges[h4].edge;
    let chip_loop = brep.faces[chip_face].outer_loop;

    brep.half_edges.remove(h3);
    brep.half_edges.remove(h4);
    brep.half_edges.remove(t3);
    brep.half_edges.remove(t4);
    brep.edges.remove(h3_edge);
    brep.edges.remove(h4_edge);
    brep.vertices.remove(v_p5);
    brep.loops.remove(chip_loop);
    brep.faces.remove(chip_face);

    // Step 5b: move v_x_new and v_a_new to their cap-tangent
    // circle positions (they were created at V's old position
    // for zero-length bridge geometry). This snaps their incident
    // edge curves' V_new-side endpoints to the new positions.
    update_incident_edges_for_vertex_move(brep, v_x_new, v_x_old_pos, plan.v_x_new_pos);
    update_incident_edges_for_vertex_move(brep, v_a_new, v_a_old_pos, plan.v_a_new_pos);

    // Step 5c: set the cut edge's 3D curve and pcurves to a
    // cap-tangent arc on cyl_m's surface, from V_x_new to V_a_new.
    // Both endpoints lie on sphere_n ∩ cyl_m by construction;
    // sphere_n is tangent to cyl_m along this circle, so G1
    // between cap_n and main cyl body is perfect along the
    // cut edge.
    let cut_edge = brep.half_edges[h2].edge;
    let main_cyl_face = mef.face;
    set_cut_edge_geometry(
        brep,
        cut_edge,
        h2,
        plan.v_x_new_pos,
        plan.v_a_new_pos,
        plan.cap_n,
        main_cyl_face,
    );

    // Step 6: validate.
    if let Err(e) = brep.validate() {
        if std::env::var("POPKERN_VAR_RADIUS_TRACE").is_ok_and(|v| v == "1") {
            eprintln!("[reroute] validate failed: {e:?}");
            if let popkern_core::error::ValidationError::OpenLoop { loop_, gap: _ } = &e
                && let Some(loop_data) = brep.loops.get(*loop_)
            {
                let face = loop_data.face;
                let kind = match brep
                    .faces
                    .get(face)
                    .and_then(|f| f.surface.as_ref())
                    .map(|s| &s.geometry)
                {
                    Some(popkern_core::geometry::surface::SurfaceGeometry::Plane { .. }) => "Plane",
                    Some(popkern_core::geometry::surface::SurfaceGeometry::Cylinder { .. }) => {
                        "Cylinder"
                    }
                    Some(popkern_core::geometry::surface::SurfaceGeometry::Sphere { .. }) => {
                        "Sphere"
                    }
                    _ => "?",
                };
                eprintln!("  loop {loop_:?} → face {face:?} ({kind})");
                let start = loop_data.half_edge;
                let mut cursor = start;
                for step in 0..16 {
                    let he_origin = brep.half_edges[cursor].origin;
                    let he_edge = brep.half_edges[cursor].edge;
                    let he_next = brep.half_edges[cursor].next;
                    let p = brep.vertices[he_origin].position;
                    eprintln!(
                        "    step {step}: he={cursor:?} origin={he_origin:?} pos=({:.3},{:.3},{:.3}) edge={he_edge:?}",
                        p.x, p.y, p.z
                    );
                    cursor = he_next;
                    if cursor == start {
                        break;
                    }
                }
            }
            if let popkern_core::error::ValidationError::ShellGap { edge, gap: _ } = &e
                && let Some(ed) = brep.edges.get(*edge)
            {
                if let Some(coe) = &ed.curve
                    && let popkern_core::geometry::CurveGeometry::Nurbs(nc) = &coe.geometry
                {
                    eprintln!(
                        "  edge {edge:?}: degree={} cps={:?}",
                        nc.degree, nc.control_points
                    );
                }
                let he = ed.half_edge;
                let twin = brep.half_edges[he].twin;
                let v_he = brep.half_edges[he].origin;
                let v_tw = brep.half_edges[twin].origin;
                let p_he = brep.vertices[v_he].position;
                let p_tw = brep.vertices[v_tw].position;
                eprintln!(
                    "  edge {edge:?} he origin={v_he:?} pos=({:.3},{:.3},{:.3}); twin origin={v_tw:?} pos=({:.3},{:.3},{:.3})",
                    p_he.x, p_he.y, p_he.z, p_tw.x, p_tw.y, p_tw.z
                );
            }
        }
        *brep = snapshot;
        return Err(RerouteError::ValidationFailed);
    }

    Ok(())
}
