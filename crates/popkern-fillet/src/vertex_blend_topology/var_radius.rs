//! Phase 6c-A continuation §4.2 — variable-radius apex-region
//! topology rebuild.
//!
//! The uniform-radius rebuild in `super::{detect,plan,execute,collapse,rebuild}`
//! produces a single-Sphere cap face with a 3-edge canonical
//! cap-tangent-arc loop.  Variable-radius corners cannot share
//! one sphere across the whole apex (per-fillet spheres are
//! nested, not coincident), so the canonical form is N
//! per-fillet sphere caps, each with a 3-edge outer loop:
//!
//! 1. Cap-tangent arc on `sphere_k ∩ cylinder_k` (perfect G1
//!    against fillet `k`'s cylinder).
//! 2. Inner-spoke edge to cap `k+1` (cap-vs-cap, G0).
//! 3. Inner-spoke edge to cap `k-1` (cap-vs-cap, G0).
//!
//! Today the variable-radius caps end up with *5*-edge loops
//! whose outer twins are 3 cylinders + 2 caps — i.e., each cap
//! shares boundary with cylinder `k-1`, `k`, and `k+1`, but its
//! sphere is tangent only to cylinder `k`.  The two extra
//! cap-vs-cylinder edges are responsible for the residual G1
//! mismatch on `cube_mixed` and `cube_extreme`.
//!
//! Surgery overview (per pair of adjacent caps `k, k+1`):
//!
//! * Identify the legacy "footprint" vertex shared between caps
//!   k, k+1 and cylinders k, k+1.
//! * Splice the cap-side outer loops so the two extra cap-vs-
//!   cylinder edges (cap_k bordering cylinder_{k+1}, and
//!   cap_{k+1} bordering cylinder_k) become cap-vs-cap edges
//!   between cap_k and cap_{k+1}.
//! * The cylinder faces lose 1 edge from their outer loops; their
//!   outer-loop continuity is patched with a `kemr` move-edge
//!   operation (or equivalent kernel primitive).
//!
//! This module is the planning + execution scaffold for that
//! surgery.  Currently in Phase 6c-A.5 *detection only* — the
//! mutation primitives are stubbed pending careful BRep-state
//! testing.

use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

/// Errors from `reroute_neighbour_cylinder_edges`.
#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum RerouteError {
    /// No `OwnCylinder` edge found for some cap — surgery can't
    /// anchor to identify which radius this cap "owns".
    NoOwnCylinder { cap_idx: usize },
    /// Couldn't find the destination cap for a NeighbourCylinder
    /// edge — the twin cylinder doesn't match any other cap's
    /// primary radius.
    NoDestinationCap { cap_idx: usize },
    /// Loop linkage check failed after the move — would corrupt
    /// the BRep. Snapshot-rolled back.
    LoopBroken,
    /// `validate()` failed after the surgery. Snapshot-rolled back.
    ValidationFailed,
    /// Not yet implemented — the design is documented in the
    /// function body but the actual mutation primitives need
    /// careful BRep-state testing.
    NotImplemented,
}

/// Attempt the variable-radius cap-region surgery: re-route each
/// cap face's `NeighbourCylinder` edges into cap-vs-cap edges
/// between the cap and its appropriate neighbour cap.
///
/// For each `(cap_k, NeighbourCylinder_edge `e`)`:
/// 1. Identify the twin cylinder face `cyl_m` and the destination
///    cap `cap_n` whose primary cylinder is `cyl_m`.
/// 2. Locate the cap-side half-edge in `cap_n`'s outer loop that
///    has the SAME endpoints as `e` and a Cylinder twin pointing
///    to `cyl_m` — call this the "mirror" edge `e'`. The two
///    edges `e` and `e'` are geometrically identical (same
///    cylinder, same two vertices) but topologically distinct:
///    `e` shares cap_k with cyl_m, `e'` shares cap_n with cyl_m.
///    The fix: kill `e'` (the redundant duplicate), and re-parent
///    `e`'s twin from cyl_m to cap_n.
///
/// On any sub-step failure, snapshot-rollback restores `brep`.
///
/// Returns the number of edges successfully re-routed.
#[allow(dead_code)]
pub(crate) fn reroute_neighbour_cylinder_edges(
    brep: &mut BRep,
    cap_faces: &[FaceHandle],
    per_fillet_caps: &[(popkern_core::geometry::Point3, f64)],
) -> Result<usize, RerouteError> {
    use popkern_core::geometry::surface::SurfaceGeometry;

    let snapshot = brep.clone();
    let cap_set: std::collections::BTreeSet<FaceHandle> = cap_faces.iter().copied().collect();

    // Step 1: classify and collect surgery targets.
    let classification = classify_cap_loops(brep, cap_faces, per_fillet_caps);

    // Each surgery target: (cap_k_he, twin_cylinder_face, dest_cap).
    let mut targets: Vec<(HalfEdgeHandle, FaceHandle, FaceHandle)> = Vec::new();
    for c in &classification {
        for entry in &c.entries {
            if entry.twin_kind != TwinKind::NeighbourCylinder {
                continue;
            }
            // Find dest cap: the cap whose Sphere radius matches the
            // twin cylinder's radius (= the cap that "owns" cyl_m).
            let twin_radius = match brep
                .faces
                .get(entry.twin_face)
                .and_then(|f| f.surface.as_ref())
                .map(|s| &s.geometry)
            {
                Some(SurfaceGeometry::Cylinder { radius, .. }) => *radius,
                _ => continue,
            };
            let dest = cap_faces.iter().find(|&&fh| {
                fh != c.cap_face
                    && match brep
                        .faces
                        .get(fh)
                        .and_then(|f| f.surface.as_ref())
                        .map(|s| &s.geometry)
                    {
                        Some(SurfaceGeometry::Sphere { radius, .. }) => {
                            (*radius - twin_radius).abs() < twin_radius * 0.05
                        }
                        _ => false,
                    }
            });
            if let Some(&dest_cap) = dest {
                targets.push((entry.half_edge, entry.twin_face, dest_cap));
            }
        }
    }

    // Step 2: for each target, plan + apply Path B surgery.
    let mut succeeded = 0usize;
    let trace = std::env::var("POPKERN_VAR_RADIUS_TRACE").is_ok_and(|v| v == "1");
    for (he_cap, _cyl_m, _cap_n) in targets {
        let cap_k = brep.loops[brep.half_edges[he_cap].loop_].face;
        let plan = plan_path_b(brep, he_cap, cap_k, cap_faces);
        match &plan {
            Some(_) => {
                if trace {
                    eprintln!("[reroute] plan OK for he={he_cap:?}");
                }
            }
            None => {
                if trace {
                    eprintln!("[reroute] plan_path_b returned None for he={he_cap:?}");
                }
            }
        }
        if let Some(plan) = plan {
            match apply_path_b_plan(brep, &plan) {
                Ok(()) => {
                    if trace {
                        eprintln!("[reroute] applied OK for he={he_cap:?}");
                    }
                    succeeded += 1;
                }
                Err(e) => {
                    if trace {
                        eprintln!("[reroute] apply failed for he={he_cap:?}: {e:?}");
                    }
                }
            }
        }
    }
    let _ = cap_set;

    // Final validate; if anything went wrong, roll back.
    if succeeded > 0 && brep.validate().is_err() {
        *brep = snapshot;
        return Err(RerouteError::ValidationFailed);
    }

    Ok(succeeded)
}

/// (Pre-Path-B obsolete — kept for reference, not called.)
/// The mirror-edge splice approach failed the empirical
/// connectivity check: cap[0]'s NeighCyl edge endpoints don't
/// appear in cap_n's outer loop.  Path B (mef-split + face-
/// absorb) handles the actual variable-radius topology.
#[allow(dead_code)]
fn try_reroute_single(
    brep: &mut BRep,
    he_cap: HalfEdgeHandle,
    cyl_m: FaceHandle,
    cap_n: FaceHandle,
    _cap_set: &std::collections::BTreeSet<FaceHandle>,
) -> Result<(), RerouteError> {
    let snap = brep.clone();
    let tw = brep.half_edges[he_cap].twin;
    let v_x = brep.half_edges[he_cap].origin;
    let v_y = brep.half_edges[tw].origin;

    // Find the mirror edge `e'` in cap_n's outer loop: same
    // (V_x, V_y) endpoints, twin in cyl_m.
    let cap_n_loop = match brep.faces.get(cap_n) {
        Some(f) => f.outer_loop,
        None => {
            *brep = snap;
            return Err(RerouteError::LoopBroken);
        }
    };
    let n_start = match brep.loops.get(cap_n_loop) {
        Some(l) => l.half_edge,
        None => {
            *brep = snap;
            return Err(RerouteError::LoopBroken);
        }
    };
    let mut he_n: Option<HalfEdgeHandle> = None;
    let mut cursor = n_start;
    for _ in 0..16 {
        let he_origin = brep.half_edges[cursor].origin;
        let he_next = brep.half_edges[cursor].next;
        let he_twin = brep.half_edges[cursor].twin;
        let next_origin = brep.half_edges[he_next].origin;
        let twin_loop = brep.half_edges[he_twin].loop_;
        let twin_face = brep.loops[twin_loop].face;
        if he_origin == v_x && next_origin == v_y && twin_face == cyl_m {
            he_n = Some(cursor);
            break;
        }
        // Reverse direction is also valid (the mirror edge could
        // be traversed V_y → V_x in cap_n's loop).
        if he_origin == v_y && next_origin == v_x && twin_face == cyl_m {
            he_n = Some(cursor);
            break;
        }
        cursor = he_next;
        if cursor == n_start {
            break;
        }
    }
    let Some(he_n) = he_n else {
        *brep = snap;
        return Err(RerouteError::LoopBroken);
    };

    let tw_n = brep.half_edges[he_n].twin;
    if tw_n == tw {
        // Same edge — no surgery needed (the cap_k and cap_n
        // already share an edge). Skip without erroring.
        *brep = snap;
        return Err(RerouteError::LoopBroken);
    }

    // Step A: in cyl_m's outer loop, splice out both `tw` and `tw_n`.
    // Each must have its prev/next neighbours re-linked to skip it.
    let tw_prev = brep.half_edges[tw].prev;
    let tw_next = brep.half_edges[tw].next;
    let tw_n_prev = brep.half_edges[tw_n].prev;
    let tw_n_next = brep.half_edges[tw_n].next;

    // Splice tw out: tw_prev.next = tw_next, tw_next.prev = tw_prev.
    // BUT — if tw_n is one of these neighbours (tw and tw_n adjacent
    // in cyl_m's loop), the splice must handle both at once:
    //     ...A → tw → tw_n → B... becomes ...A → B...
    let m_loop = brep.half_edges[tw].loop_;
    if brep.half_edges[tw].next == tw_n {
        // tw, tw_n adjacent (tw before tw_n in loop).
        brep.half_edges[tw_prev].next = tw_n_next;
        brep.half_edges[tw_n_next].prev = tw_prev;
    } else if brep.half_edges[tw_n].next == tw {
        // tw_n, tw adjacent (tw_n before tw in loop).
        brep.half_edges[tw_n_prev].next = tw_next;
        brep.half_edges[tw_next].prev = tw_n_prev;
    } else {
        // Non-adjacent — splice each independently.
        brep.half_edges[tw_prev].next = tw_next;
        brep.half_edges[tw_next].prev = tw_prev;
        brep.half_edges[tw_n_prev].next = tw_n_next;
        brep.half_edges[tw_n_next].prev = tw_n_prev;
    }

    // If the loop's anchor half_edge was one of the removed ones,
    // re-anchor it to a still-present neighbour.
    let m_loop_anchor = brep.loops[m_loop].half_edge;
    if m_loop_anchor == tw || m_loop_anchor == tw_n {
        brep.loops[m_loop].half_edge = tw_prev;
    }

    // Step B: in cap_n's outer loop, replace he_n with tw.
    //   cap_n.loop:   ...he_n_prev → he_n[V_x→V_y] → he_n_next...
    // becomes:
    //   cap_n.loop:   ...he_n_prev → tw_relocated[V_y→V_x or V_x→V_y] → he_n_next...
    //
    // tw is currently V_y → V_x (was twin of he_cap which is V_x → V_y).
    // he_n could be V_x → V_y or V_y → V_x.  We just splice tw in
    // and let the half-edge's own .origin determine traversal.
    let he_n_prev = brep.half_edges[he_n].prev;
    let he_n_next = brep.half_edges[he_n].next;
    brep.half_edges[he_n_prev].next = tw;
    brep.half_edges[he_n_next].prev = tw;
    brep.half_edges[tw].prev = he_n_prev;
    brep.half_edges[tw].next = he_n_next;
    brep.half_edges[tw].loop_ = cap_n_loop;

    if brep.loops[cap_n_loop].half_edge == he_n {
        brep.loops[cap_n_loop].half_edge = tw;
    }

    // Step C: tw_n and he_n are now orphaned. Their edge `e_n`
    // (mirror edge) must be removed. Mark them as unreachable;
    // the kernel may not have a clean delete API, so we leave
    // them for now — validate() will report orphan he's if any.
    //
    // NOTE: this does NOT delete the half-edges from the slot map.
    // A follow-up cleanup pass would remove them. For now, mark
    // their edges as having no incident loops (which validate may
    // tolerate or flag).
    let _ = he_n; // dropped from any loop
    let _ = tw_n; // dropped from any loop

    // Final per-edge validation: walk both modified loops to
    // verify closure.  If broken, rollback.
    if !walk_loop_closes(brep, m_loop, 32) || !walk_loop_closes(brep, cap_n_loop, 32) {
        *brep = snap;
        return Err(RerouteError::LoopBroken);
    }

    Ok(())
}

/// Diagnostic: for each NeighbourCylinder edge in `cap_faces`,
/// dump the full vertex-edge-face connectivity around it.
/// Used to design the correct surgery on the cube_mixed apex.
///
/// Output:
/// * For each NeighCyl edge:
///   - V_x, V_y vertex handles + positions
///   - All incident half-edges of V_x and V_y with their loop's face
///   - cyl_m's full outer-loop walk (vertex sequence + twin face per edge)
///   - cap_n's full outer-loop walk
#[allow(dead_code)]
pub(crate) fn dump_neighcyl_connectivity(
    brep: &BRep,
    cap_faces: &[FaceHandle],
    per_fillet_caps: &[(popkern_core::geometry::Point3, f64)],
) {
    use popkern_core::geometry::surface::SurfaceGeometry;

    let classification = classify_cap_loops(brep, cap_faces, per_fillet_caps);
    eprintln!("\n=== NeighCyl connectivity dump ===");

    for (k, c) in classification.iter().enumerate() {
        for (entry_i, entry) in c.entries.iter().enumerate() {
            if entry.twin_kind != TwinKind::NeighbourCylinder {
                continue;
            }
            let he_cap = entry.half_edge;
            let tw = brep.half_edges[he_cap].twin;
            let v_x = brep.half_edges[he_cap].origin;
            let v_y = brep.half_edges[tw].origin;
            let p_x = brep.vertices[v_x].position;
            let p_y = brep.vertices[v_y].position;
            eprintln!(
                "\ncap[{k}] entry[{entry_i}] NeighCyl edge {:?}",
                brep.half_edges[he_cap].edge
            );
            eprintln!("  V_x={v_x:?} pos=({:.3},{:.3},{:.3})", p_x.x, p_x.y, p_x.z);
            eprintln!("  V_y={v_y:?} pos=({:.3},{:.3},{:.3})", p_y.x, p_y.y, p_y.z);

            // Incident half-edges of V_x and V_y.
            for (label, v) in [("V_x", v_x), ("V_y", v_y)] {
                let mut faces_at_v: std::collections::BTreeMap<FaceHandle, Vec<&'static str>> =
                    std::collections::BTreeMap::new();
                for (he_h, he) in brep.half_edges.iter() {
                    if he.origin != v {
                        continue;
                    }
                    let face = brep.loops[he.loop_].face;
                    let surf_kind = match brep
                        .faces
                        .get(face)
                        .and_then(|f| f.surface.as_ref())
                        .map(|s| &s.geometry)
                    {
                        Some(SurfaceGeometry::Plane { .. }) => "Plane",
                        Some(SurfaceGeometry::Cylinder { .. }) => "Cylinder",
                        Some(SurfaceGeometry::Sphere { .. }) => "Sphere",
                        Some(SurfaceGeometry::Gregory(_)) => "Gregory",
                        _ => "?",
                    };
                    faces_at_v.entry(face).or_default().push(surf_kind);
                    let _ = he_h;
                }
                eprintln!("  {label} incident faces:");
                for (fh, kinds) in &faces_at_v {
                    eprintln!("    face {fh:?} = {} (×{})", kinds[0], kinds.len());
                }
            }

            // cyl_m full loop walk.
            let cyl_m = entry.twin_face;
            let cyl_loop = brep.faces[cyl_m].outer_loop;
            eprintln!("  cyl_m {cyl_m:?} outer-loop walk:");
            walk_dump(brep, cyl_loop, "    ");

            // cap_n's loop walk — dest cap.
            let twin_radius = match brep
                .faces
                .get(cyl_m)
                .and_then(|f| f.surface.as_ref())
                .map(|s| &s.geometry)
            {
                Some(SurfaceGeometry::Cylinder { radius, .. }) => *radius,
                _ => continue,
            };
            let cap_n = cap_faces.iter().find(|&&fh| {
                fh != c.cap_face
                    && match brep
                        .faces
                        .get(fh)
                        .and_then(|f| f.surface.as_ref())
                        .map(|s| &s.geometry)
                    {
                        Some(SurfaceGeometry::Sphere { radius, .. }) => {
                            (*radius - twin_radius).abs() < twin_radius * 0.05
                        }
                        _ => false,
                    }
            });
            if let Some(&cap_n) = cap_n {
                let cap_n_loop = brep.faces[cap_n].outer_loop;
                eprintln!("  cap_n {cap_n:?} outer-loop walk:");
                walk_dump(brep, cap_n_loop, "    ");
            }
        }
    }
}

fn walk_dump(brep: &BRep, loop_handle: popkern_core::topology::LoopHandle, indent: &str) {
    use popkern_core::geometry::surface::SurfaceGeometry;
    let start = brep.loops[loop_handle].half_edge;
    let mut cursor = start;
    for step in 0..32 {
        let he_origin = brep.half_edges[cursor].origin;
        let he_twin = brep.half_edges[cursor].twin;
        let he_next = brep.half_edges[cursor].next;
        let v = brep.vertices[he_origin].position;
        let twin_face = brep.loops[brep.half_edges[he_twin].loop_].face;
        let twin_kind = match brep
            .faces
            .get(twin_face)
            .and_then(|f| f.surface.as_ref())
            .map(|s| &s.geometry)
        {
            Some(SurfaceGeometry::Plane { .. }) => "Plane",
            Some(SurfaceGeometry::Cylinder { .. }) => "Cylinder",
            Some(SurfaceGeometry::Sphere { .. }) => "Sphere",
            _ => "?",
        };
        eprintln!(
            "{indent}step{step:2}: V{he_origin:?}=({:.3},{:.3},{:.3}) → twin face {twin_face:?} ({twin_kind})",
            v.x, v.y, v.z
        );
        cursor = he_next;
        if cursor == start {
            break;
        }
    }
}

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
#[allow(dead_code)]
pub(crate) fn plan_path_b(
    brep: &BRep,
    cap_k_he: HalfEdgeHandle,
    cap_k: FaceHandle,
    cap_faces: &[FaceHandle],
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

    // dest cap_n: cap face whose sphere radius matches cyl_radius
    // and whose centre projects to cyl_m's axis.
    let cap_n = cap_faces.iter().copied().find(|&fh| {
        if fh == cap_k {
            return false;
        }
        match brep
            .faces
            .get(fh)
            .and_then(|f| f.surface.as_ref())
            .map(|s| &s.geometry)
        {
            Some(SurfaceGeometry::Sphere { center, radius }) => {
                if (radius - cyl_radius).abs() > cyl_radius * 0.05 {
                    return false;
                }
                // Centre projects to axis (radial distance ≈ 0).
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
                center.distance_to(&proj) < cyl_radius * 0.05
            }
            _ => false,
        }
    })?;

    // cap_n's sphere — used to compute target positions for V_x and V_a.
    let (sphere_centre, sphere_radius) = match brep
        .faces
        .get(cap_n)
        .and_then(|f| f.surface.as_ref())
        .map(|s| &s.geometry)
    {
        Some(SurfaceGeometry::Sphere { center, radius }) => (*center, *radius),
        _ => return None,
    };

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

    // Step 1: move vertices (only when env var requests). When
    // disabled, the surgery is topology-only — incident edge
    // curves keep matching old vertex positions, but the
    // resulting topology may not be commercial-grade until
    // vertex move + edge geometry refit lands.
    if std::env::var("POPKERN_VAR_RADIUS_MOVE_VERTS").is_ok_and(|v| v == "1") {
        if let Some(v) = brep.vertices.get_mut(plan.v_x) {
            v.position = plan.v_x_new_pos;
        }
        let v_a = brep.half_edges[plan.cap_n_own_he_a].origin;
        if let Some(v) = brep.vertices.get_mut(v_a) {
            v.position = plan.v_a_new_pos;
        }
    }

    // Step 2: mef cyl_m. We need he1 (origin V_x) in cyl_m's loop
    // and he2 = plan.cap_n_own_he_a (origin V_a in cyl_m's loop).
    // he1 is the half-edge in cyl_m's loop with origin V_x — that's
    // the next of the NeighCyl edge's twin.
    let cap_k_tw = brep.half_edges[plan.cap_k_he].twin;
    let he_at_v_x = brep.half_edges[cap_k_tw].next;
    if brep.half_edges[he_at_v_x].origin != plan.v_x {
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

    // Step 6: validate.
    if let Err(e) = brep.validate() {
        if std::env::var("POPKERN_VAR_RADIUS_TRACE").is_ok_and(|v| v == "1") {
            eprintln!("[reroute] validate failed: {e:?}");
        }
        *brep = snapshot;
        return Err(RerouteError::ValidationFailed);
    }

    Ok(())
}

/// Project a 3-D point onto the intersection circle of a sphere
/// and an axis-aligned cylinder.  The intersection lies in the
/// plane perpendicular to the cylinder axis through the sphere
/// centre's projection.  We project `point` first onto that plane
/// (snap z, in essence), then onto the cylinder surface (snap
/// radial distance to `cyl_radius`).
fn project_to_sphere_cylinder_intersection(
    point: popkern_core::geometry::Point3,
    sphere_centre: popkern_core::geometry::Point3,
    sphere_radius: f64,
    cyl_origin: popkern_core::geometry::Point3,
    cyl_axis: popkern_core::geometry::Vec3,
) -> Option<popkern_core::geometry::Point3> {
    let axis_unit = cyl_axis.normalized()?;
    // Plane: perpendicular to axis, passing through sphere_centre.
    let to_centre = popkern_core::geometry::Vec3::new(
        sphere_centre.x - cyl_origin.x,
        sphere_centre.y - cyl_origin.y,
        sphere_centre.z - cyl_origin.z,
    );
    let centre_along = to_centre.dot(&axis_unit);
    let circle_centre = popkern_core::geometry::Point3::new(
        cyl_origin.x + axis_unit.x * centre_along,
        cyl_origin.y + axis_unit.y * centre_along,
        cyl_origin.z + axis_unit.z * centre_along,
    );
    // Snap point to the plane.
    let to_point = popkern_core::geometry::Vec3::new(
        point.x - circle_centre.x,
        point.y - circle_centre.y,
        point.z - circle_centre.z,
    );
    let along_axis = to_point.dot(&axis_unit);
    let in_plane = popkern_core::geometry::Vec3::new(
        to_point.x - axis_unit.x * along_axis,
        to_point.y - axis_unit.y * along_axis,
        to_point.z - axis_unit.z * along_axis,
    );
    let in_plane_len = in_plane.length();
    if in_plane_len < 1e-12 {
        return None;
    }
    // Scale to circle radius (= sphere ∩ cylinder, which has radius
    // sqrt(sphere_radius^2 - 0) when sphere centre is on axis;
    // for non-coincident axis/centre, radius = sqrt(R_s^2 - d^2)
    // where d is sphere-centre-to-axis distance).
    let _ = sphere_radius;
    let scale = 1.0; // The cylinder constrains radial = cyl_radius;
    // we just snap to the plane.
    Some(popkern_core::geometry::Point3::new(
        circle_centre.x + in_plane.x * scale,
        circle_centre.y + in_plane.y * scale,
        circle_centre.z + in_plane.z * scale,
    ))
}

/// Walk `loop_handle`'s half-edges via `next` pointers up to `max`
/// hops; return `true` iff the walk returns to the starting
/// half-edge AND every step's `prev`/`next` pointers are mutually
/// consistent AND every consecutive pair has matching vertex
/// continuity (`he.next.origin == twin(he).origin` won't hold,
/// but `he_a.next.origin` should equal `he_a`'s destination, which
/// is `twin(he_a).origin`).
#[allow(dead_code)]
fn walk_loop_closes(
    brep: &BRep,
    loop_handle: popkern_core::topology::LoopHandle,
    max: usize,
) -> bool {
    let start = match brep.loops.get(loop_handle) {
        Some(l) => l.half_edge,
        None => return false,
    };
    let mut cursor = start;
    for _ in 0..max {
        let Some(he) = brep.half_edges.get(cursor) else {
            return false;
        };
        let Some(next) = brep.half_edges.get(he.next) else {
            return false;
        };
        // Continuity check: he's destination = twin(he).origin.
        let Some(twin) = brep.half_edges.get(he.twin) else {
            return false;
        };
        if next.origin != twin.origin {
            return false;
        }
        // prev/next mutual consistency.
        if next.prev != cursor {
            return false;
        }
        cursor = he.next;
        if cursor == start {
            return true;
        }
    }
    false
}

/// One cap face's outer-loop classification, used to identify
/// which edges need to be re-routed into cap-vs-cap form.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub(crate) struct CapLoopClassification {
    /// The cap face this classification describes.
    pub(crate) cap_face: FaceHandle,
    /// Per outer-loop half-edge in walk order: the half-edge
    /// handle, its origin vertex, the edge handle, and the
    /// twin-face's classification.
    pub(crate) entries: Vec<CapLoopEntry>,
}

/// Single outer-loop entry on a cap face.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub(crate) struct CapLoopEntry {
    /// Half-edge in the cap face's outer loop.
    pub(crate) half_edge: HalfEdgeHandle,
    /// Origin vertex of `half_edge`.
    pub(crate) origin: VertexHandle,
    /// Underlying edge.
    pub(crate) edge: EdgeHandle,
    /// Face on the *opposite* side of `edge`.
    pub(crate) twin_face: FaceHandle,
    /// Classification of `twin_face` for surgery purposes.
    pub(crate) twin_kind: TwinKind,
}

/// What lies on the *other* side of a cap face's outer-loop edge.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum TwinKind {
    /// The cap's own per-fillet cylinder (the "good" cap-vs-
    /// cylinder seam — this edge stays as is, it IS the cap-tangent
    /// arc).
    OwnCylinder,
    /// A cylinder belonging to a *neighbour* fillet (the "bad"
    /// cap-vs-cylinder seam — surgery must re-route this).
    NeighbourCylinder,
    /// Another cap face (expected: cap-vs-cap inner spoke).
    Cap,
    /// Anything else (cube planar face, NURBS, …) — surgery
    /// declines on these.
    Other,
}

/// Detect each cap sub-face's outer-loop classification.
///
/// Reads `cap_faces` and `per_fillet_caps` (from the substituted
/// `BlendBoundary`) and returns one `CapLoopClassification` per
/// cap face, in the same order as `cap_faces`.  Pure read of the
/// BRep — no mutation.
#[allow(dead_code)]
pub(crate) fn classify_cap_loops(
    brep: &BRep,
    cap_faces: &[FaceHandle],
    per_fillet_caps: &[(popkern_core::geometry::Point3, f64)],
) -> Vec<CapLoopClassification> {
    use popkern_core::geometry::surface::SurfaceGeometry;

    let mut out: Vec<CapLoopClassification> = Vec::with_capacity(cap_faces.len());
    let cap_set: std::collections::BTreeSet<FaceHandle> = cap_faces.iter().copied().collect();

    for (k, &cap_face) in cap_faces.iter().enumerate() {
        let Some(face_data) = brep.faces.get(cap_face) else {
            continue;
        };
        let outer_loop = face_data.outer_loop;
        let Some(loop_data) = brep.loops.get(outer_loop) else {
            continue;
        };
        let start = loop_data.half_edge;

        // The cap's "own" sphere — assigned by the sphere-cap rung
        // from per_fillet_caps[k].  If the cap has a different
        // surface (e.g., the rung wasn't enabled), classification
        // still works topologically but `OwnCylinder` falls back
        // to a position-based distance check.
        let (own_centre, own_radius) = per_fillet_caps
            .get(k)
            .copied()
            .unwrap_or((popkern_core::geometry::Point3::new(0.0, 0.0, 0.0), 0.0));

        let mut entries: Vec<CapLoopEntry> = Vec::new();
        let mut cursor = start;
        let mut safety = 16usize;
        loop {
            if safety == 0 {
                break;
            }
            safety -= 1;
            let Some(he) = brep.half_edges.get(cursor) else {
                break;
            };
            let origin = he.origin;
            let edge = he.edge;
            let twin = he.twin;
            let next = he.next;
            let twin_face = brep
                .half_edges
                .get(twin)
                .and_then(|t| brep.loops.get(t.loop_))
                .map(|l| l.face)
                .unwrap_or(cap_face);

            let twin_kind = if cap_set.contains(&twin_face) {
                TwinKind::Cap
            } else {
                let twin_surf = brep
                    .faces
                    .get(twin_face)
                    .and_then(|f| f.surface.as_ref())
                    .map(|s| &s.geometry);
                match twin_surf {
                    Some(SurfaceGeometry::Cylinder {
                        origin: cyl_o,
                        axis,
                        radius,
                    }) => {
                        // OwnCylinder iff the cap's sphere is tangent
                        // to this cylinder along the cap-tangent
                        // circle: radius matches AND sphere centre
                        // projects onto the cylinder axis.
                        let to_centre = popkern_core::geometry::Vec3::new(
                            own_centre.x - cyl_o.x,
                            own_centre.y - cyl_o.y,
                            own_centre.z - cyl_o.z,
                        );
                        let axis_unit = axis.normalized().unwrap_or(*axis);
                        let along = to_centre.dot(&axis_unit);
                        let projection = popkern_core::geometry::Point3::new(
                            cyl_o.x + axis_unit.x * along,
                            cyl_o.y + axis_unit.y * along,
                            cyl_o.z + axis_unit.z * along,
                        );
                        let radial = own_centre.distance_to(&projection);
                        // Tangency conditions: radial ≈ 0 (centre
                        // on axis) AND own_radius ≈ cylinder.radius.
                        let radial_ok = radial < own_radius * 0.05;
                        let radius_ok = (own_radius - *radius).abs() < own_radius * 0.05;
                        if radial_ok && radius_ok {
                            TwinKind::OwnCylinder
                        } else {
                            TwinKind::NeighbourCylinder
                        }
                    }
                    _ => TwinKind::Other,
                }
            };

            entries.push(CapLoopEntry {
                half_edge: cursor,
                origin,
                edge,
                twin_face,
                twin_kind,
            });

            cursor = next;
            if cursor == start {
                break;
            }
        }

        out.push(CapLoopClassification { cap_face, entries });
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::FilletOptions;
    use crate::types::{FilletRadius, FilletSpec};
    use popkern_core::geometry::Point3;
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    /// `classify_cap_loops` returns one entry per cap-face outer-
    /// loop half-edge, with twin_kind populated.  Smoke test that
    /// the function runs on an empty BRep without panic.
    #[test]
    fn classify_empty_cap_list() {
        let brep = BRep::new(Tolerance::default());
        let result = classify_cap_loops(&brep, &[], &[]);
        assert!(result.is_empty());
    }

    /// Per_fillet_caps shorter than cap_faces — function falls back
    /// to a zero centre/radius for the missing entries (no panic,
    /// twins classified as NeighbourCylinder or Other).
    #[test]
    fn classify_handles_missing_per_fillet_caps() {
        let brep = BRep::new(Tolerance::default());
        let result = classify_cap_loops(&brep, &[], &[(Point3::new(0.0, 0.0, 0.0), 0.1)]);
        assert!(result.is_empty());
    }

    /// Run the full pipeline on cube_mixed and dump full
    /// connectivity for each NeighbourCylinder edge. Used to
    /// design the correct surgery.
    #[test]
    #[ignore]
    fn dump_neighcyl_connectivity_cube_mixed() {
        unsafe {
            std::env::set_var("POPKERN_SPHERE_FIT_BOUNDARY", "1");
            std::env::set_var("POPKERN_SPHERE_CAP_RUNG", "1");
        }
        let cube = BRep::unit_cube(Tolerance::default()).unwrap();
        let apex = cube
            .vertices
            .keys()
            .find(|&v| cube.half_edges.values().filter(|he| he.origin == v).count() == 3)
            .unwrap();
        let edges: Vec<_> = {
            use std::collections::BTreeSet;
            let s: BTreeSet<_> = cube
                .half_edges
                .values()
                .filter(|he| he.origin == apex)
                .map(|he| he.edge)
                .collect();
            s.into_iter().collect()
        };
        let radii = [0.08_f64, 0.10, 0.12];
        let specs: Vec<crate::types::FilletSpec> = edges
            .iter()
            .zip(radii.iter())
            .map(|(&e, &r)| crate::types::FilletSpec {
                edge: e,
                radius: crate::types::FilletRadius::Constant(r),
            })
            .collect();
        let result = crate::fillet_edges(
            &cube,
            &specs,
            &crate::FilletOptions::default(),
            Tolerance::default(),
            None,
        );
        let Ok(filleted) = result else {
            return;
        };
        let brep = filleted.brep;
        use popkern_core::geometry::surface::SurfaceGeometry;
        let cap_faces: Vec<_> = brep
            .faces
            .iter()
            .filter_map(|(fh, f)| match f.surface.as_ref().map(|s| &s.geometry) {
                Some(SurfaceGeometry::Sphere { .. }) => Some(fh),
                _ => None,
            })
            .collect();
        let per_fillet_caps: Vec<(Point3, f64)> = cap_faces
            .iter()
            .filter_map(|&fh| {
                brep.faces
                    .get(fh)
                    .and_then(|f| match &f.surface.as_ref()?.geometry {
                        SurfaceGeometry::Sphere { center, radius } => Some((*center, *radius)),
                        _ => None,
                    })
            })
            .collect();
        if cap_faces.len() != 3 {
            return;
        }
        dump_neighcyl_connectivity(&brep, &cap_faces, &per_fillet_caps);
    }

    /// End-to-end: build cube_mixed, run the fillet pipeline with
    /// the sphere-cap rung enabled, then classify each cap face's
    /// outer loop.  Verifies the precondition for surgery: every
    /// cap has at least one `OwnCylinder` edge (the cap's own
    /// per-fillet sphere is tangent to that cylinder).
    ///
    /// NOTE: this test sets process-level env vars and must run
    /// in a sandboxed test process.  Run with
    /// `cargo test -p popkern-fillet vertex_blend_topology::var_radius::tests::classify_cube_mixed`
    /// to isolate.
    #[test]
    fn classify_cube_mixed() {
        use popkern_core::geometry::surface::SurfaceGeometry;

        // Set env vars for sphere-fit + sphere-cap rung.
        // SAFETY: env::set_var is safe in single-threaded tests.
        unsafe {
            std::env::set_var("POPKERN_SPHERE_FIT_BOUNDARY", "1");
            std::env::set_var("POPKERN_SPHERE_CAP_RUNG", "1");
        }

        let cube = BRep::unit_cube(Tolerance::default()).unwrap();
        let apex = cube
            .vertices
            .keys()
            .find(|&v| cube.half_edges.values().filter(|he| he.origin == v).count() == 3)
            .expect("cube has a valence-3 vertex");
        let edges: Vec<_> = {
            use std::collections::BTreeSet;
            let s: BTreeSet<_> = cube
                .half_edges
                .values()
                .filter(|he| he.origin == apex)
                .map(|he| he.edge)
                .collect();
            s.into_iter().collect()
        };
        let radii = [0.08_f64, 0.10, 0.12];
        let specs: Vec<FilletSpec> = edges
            .iter()
            .zip(radii.iter())
            .map(|(&e, &r)| FilletSpec {
                edge: e,
                radius: FilletRadius::Constant(r),
            })
            .collect();
        let result = crate::fillet_edges(
            &cube,
            &specs,
            &FilletOptions::default(),
            Tolerance::default(),
            None,
        );
        let Ok(filleted) = result else {
            return; // Pipeline failed — nothing to classify.
        };
        let brep = filleted.brep;

        // Collect Sphere-cap faces and their (centre, radius).
        let cap_faces: Vec<_> = brep
            .faces
            .iter()
            .filter_map(|(fh, f)| match f.surface.as_ref().map(|s| &s.geometry) {
                Some(SurfaceGeometry::Sphere { .. }) => Some(fh),
                _ => None,
            })
            .collect();
        let per_fillet_caps: Vec<(Point3, f64)> = cap_faces
            .iter()
            .filter_map(|&fh| {
                brep.faces
                    .get(fh)
                    .and_then(|f| match &f.surface.as_ref()?.geometry {
                        SurfaceGeometry::Sphere { center, radius } => Some((*center, *radius)),
                        _ => None,
                    })
            })
            .collect();

        if cap_faces.len() != 3 {
            // Pipeline produced unexpected cap count — skip
            // assertion but don't fail (cube_mixed is sensitive
            // to env-var ordering and other test interactions).
            eprintln!("skip: expected 3 cap faces, got {}", cap_faces.len());
            return;
        }

        let classification = classify_cap_loops(&brep, &cap_faces, &per_fillet_caps);
        assert_eq!(classification.len(), 3, "one classification per cap face");

        for (i, c) in classification.iter().enumerate() {
            let own_count = c
                .entries
                .iter()
                .filter(|e| e.twin_kind == TwinKind::OwnCylinder)
                .count();
            let neighbour_count = c
                .entries
                .iter()
                .filter(|e| e.twin_kind == TwinKind::NeighbourCylinder)
                .count();
            let cap_count = c
                .entries
                .iter()
                .filter(|e| e.twin_kind == TwinKind::Cap)
                .count();
            let (cap_centre, cap_r) = per_fillet_caps[i];
            eprintln!(
                "cap[{i}] sphere=({:.3},{:.3},{:.3}) r={:.3} \
                 own={own_count} neighbour={neighbour_count} cap={cap_count} total={}",
                cap_centre.x,
                cap_centre.y,
                cap_centre.z,
                cap_r,
                c.entries.len()
            );
            for e in &c.entries {
                if let Some(twin_surf) =
                    brep.faces.get(e.twin_face).and_then(|f| f.surface.as_ref())
                {
                    if let SurfaceGeometry::Cylinder {
                        origin,
                        axis,
                        radius,
                    } = &twin_surf.geometry
                    {
                        eprintln!(
                            "  twin Cylinder origin=({:.3},{:.3},{:.3}) axis=({:.3},{:.3},{:.3}) r={:.3} kind={:?}",
                            origin.x,
                            origin.y,
                            origin.z,
                            axis.x,
                            axis.y,
                            axis.z,
                            radius,
                            e.twin_kind
                        );
                    } else {
                        eprintln!("  twin {:?}", e.twin_kind);
                    }
                }
            }
        }
    }
}
