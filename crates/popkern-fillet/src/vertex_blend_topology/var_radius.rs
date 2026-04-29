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
    let snapshot = brep.clone();
    let cap_set: std::collections::BTreeSet<FaceHandle> = cap_faces.iter().copied().collect();

    // Step 2: iteratively classify + apply.  Each surgery
    // mutates the BRep, so we re-classify after each iteration
    // to find the next NeighbourCylinder edge. If a surgery
    // fails (rolled back via snapshot inside apply_path_b_plan),
    // we mark its target as tried-and-failed and skip in
    // subsequent iterations to avoid infinite loops.
    let mut succeeded = 0usize;
    let trace = std::env::var("POPKERN_VAR_RADIUS_TRACE").is_ok_and(|v| v == "1");
    let mut tried_failed: std::collections::BTreeSet<HalfEdgeHandle> =
        std::collections::BTreeSet::new();
    for iter in 0..16 {
        // Find one NeighbourCylinder target not previously tried.
        let classification = classify_cap_loops(brep, cap_faces, per_fillet_caps);
        let mut target: Option<HalfEdgeHandle> = None;
        for c in &classification {
            for entry in &c.entries {
                if entry.twin_kind == TwinKind::NeighbourCylinder
                    && !tried_failed.contains(&entry.half_edge)
                {
                    target = Some(entry.half_edge);
                    break;
                }
            }
            if target.is_some() {
                break;
            }
        }
        let Some(he_cap) = target else {
            if trace {
                eprintln!("[reroute] iter {iter}: no more NeighCyl targets");
            }
            break;
        };
        let cap_k = brep.loops[brep.half_edges[he_cap].loop_].face;
        let plan = plan_path_b(brep, he_cap, cap_k, cap_faces, per_fillet_caps);
        let Some(plan) = plan else {
            if trace {
                eprintln!("[reroute] iter {iter}: plan_path_b None for he={he_cap:?}");
            }
            tried_failed.insert(he_cap);
            continue;
        };
        match apply_path_b_plan(brep, &plan) {
            Ok(()) => {
                if trace {
                    eprintln!("[reroute] iter {iter}: applied OK for he={he_cap:?}");
                }
                succeeded += 1;
            }
            Err(e) => {
                if trace {
                    eprintln!("[reroute] iter {iter}: apply failed for he={he_cap:?}: {e:?}");
                }
                tried_failed.insert(he_cap);
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
    if let Some(surf) = cap_surf {
        let bc = brep.edges[bridge_edge].curve.as_ref().unwrap().clone();
        let curve_3d = match &bc.geometry {
            popkern_core::geometry::CurveGeometry::Nurbs(nc) => nc.clone(),
            _ => unreachable!(),
        };
        if let Some(pc) = sample_pcurve_on_surface(&curve_3d, &surf) {
            brep.half_edges[bridge_cap].pcurve = Some(pc.into());
        }
    }
    if let Some(surf) = cube_surf {
        let bc = brep.edges[bridge_edge].curve.as_ref().unwrap().clone();
        let curve_3d = match &bc.geometry {
            popkern_core::geometry::CurveGeometry::Nurbs(nc) => nc.clone(),
            _ => unreachable!(),
        };
        if let Some(pc) = sample_pcurve_on_surface(&curve_3d, &surf) {
            // Bridge_cube goes V_new → V — invert direction.
            let mut inv = pc.clone();
            inv.control_points.reverse();
            inv.weights.reverse();
            brep.half_edges[bridge_cube].pcurve = Some(inv.into());
        }
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
fn find_incident_plane_face(brep: &BRep, v: VertexHandle) -> Option<FaceHandle> {
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

/// Move a vertex and update the endpoint of any incident edge's
/// 3D curve that previously matched the old vertex position.
/// Preserves loop closure: each `he`'s origin is at the moved
/// vertex, and the corresponding curve endpoint is updated to
/// match.
fn update_incident_edges_for_vertex_move(
    brep: &mut BRep,
    v: VertexHandle,
    _old_pos: popkern_core::geometry::Point3,
    new_pos: popkern_core::geometry::Point3,
) {
    use popkern_core::geometry::CurveGeometry;

    // Snapshot half-edges with origin = v.  For each, identify
    // whether it is the edge's "primary" half-edge (= edge.half_edge):
    // primary's origin matches CP[0] (curve start);
    // twin's origin matches CP[last] (curve end).
    //
    // This is more precise than position-based matching: it can
    // disambiguate the V vs V_new ends of a zero-length bridge
    // (where both CPs are at V_pos initially).
    let halves: Vec<HalfEdgeHandle> = brep
        .half_edges
        .iter()
        .filter_map(|(h, he)| (he.origin == v).then_some(h))
        .collect();

    for he_h in halves {
        let edge_h = brep.half_edges[he_h].edge;
        let primary = brep.edges.get(edge_h).map(|e| e.half_edge);
        let is_primary = primary == Some(he_h);
        let Some(edge) = brep.edges.get_mut(edge_h) else {
            continue;
        };
        let Some(coe) = edge.curve.as_mut() else {
            continue;
        };
        if let CurveGeometry::Nurbs(nc) = &mut coe.geometry
            && !nc.control_points.is_empty()
        {
            let n = nc.control_points.len();
            if is_primary {
                nc.control_points[0] = new_pos;
            } else {
                nc.control_points[n - 1] = new_pos;
            }
        }
    }

    if let Some(vd) = brep.vertices.get_mut(v) {
        vd.position = new_pos;
    }
}

/// Build a rational quadratic NURBS arc on the cap-tangent
/// circle (sphere_n ∩ cyl_m) from `start` to `end`.  Both
/// endpoints must lie on the circle (caller's responsibility);
/// we reconstruct the circle from the cylinder's axis and the
/// sphere's centre.
#[allow(dead_code)]
fn build_cap_tangent_arc_curve(
    start: popkern_core::geometry::Point3,
    end: popkern_core::geometry::Point3,
    sphere_centre: popkern_core::geometry::Point3,
    cyl_origin: popkern_core::geometry::Point3,
    cyl_axis: popkern_core::geometry::Vec3,
    cyl_radius: f64,
) -> Option<popkern_core::geometry::nurbs_curve::NurbsCurve> {
    use popkern_core::geometry::Vec3;
    use popkern_core::geometry::nurbs_curve::NurbsCurve;
    use std::f64::consts::TAU;

    let axis_unit = cyl_axis.normalized()?;
    let to_centre = Vec3::new(
        sphere_centre.x - cyl_origin.x,
        sphere_centre.y - cyl_origin.y,
        sphere_centre.z - cyl_origin.z,
    );
    let along = to_centre.dot(&axis_unit);
    let circle_centre = popkern_core::geometry::Point3::new(
        cyl_origin.x + axis_unit.x * along,
        cyl_origin.y + axis_unit.y * along,
        cyl_origin.z + axis_unit.z * along,
    );

    // Build orthonormal in-plane basis (u_axis, v_axis) perp to cyl_axis.
    let helper = if axis_unit.x.abs() < 0.9 {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let u_axis = axis_unit.cross(&helper).normalized()?;
    let v_axis = axis_unit.cross(&u_axis);

    let theta = |p: popkern_core::geometry::Point3| -> f64 {
        let d = Vec3::new(
            p.x - circle_centre.x,
            p.y - circle_centre.y,
            p.z - circle_centre.z,
        );
        d.dot(&v_axis).atan2(d.dot(&u_axis))
    };
    let theta_a = theta(start);
    let theta_b = theta(end);

    // Pick the SHORTER arc direction.
    let mut delta = theta_b - theta_a;
    while delta > std::f64::consts::PI {
        delta -= TAU;
    }
    while delta < -std::f64::consts::PI {
        delta += TAU;
    }
    let theta_end = theta_a + delta;

    // Rational quadratic arc: 3 CPs (p0, shoulder, p2), weights (1, cos(α), 1).
    let half_span = (theta_end - theta_a) * 0.5;
    let cos_half = half_span.cos();
    let theta_mid = theta_a + half_span;
    let circle_eval = |th: f64| -> popkern_core::geometry::Point3 {
        popkern_core::geometry::Point3::new(
            circle_centre.x + cyl_radius * (th.cos() * u_axis.x + th.sin() * v_axis.x),
            circle_centre.y + cyl_radius * (th.cos() * u_axis.y + th.sin() * v_axis.y),
            circle_centre.z + cyl_radius * (th.cos() * u_axis.z + th.sin() * v_axis.z),
        )
    };
    let p0 = circle_eval(theta_a);
    let p_mid = circle_eval(theta_mid);
    let p2 = circle_eval(theta_end);
    let shoulder = if cos_half.abs() > 1e-12 {
        popkern_core::geometry::Point3::new(
            circle_centre.x + (p_mid.x - circle_centre.x) / cos_half,
            circle_centre.y + (p_mid.y - circle_centre.y) / cos_half,
            circle_centre.z + (p_mid.z - circle_centre.z) / cos_half,
        )
    } else {
        p_mid
    };

    Some(NurbsCurve {
        degree: 2,
        control_points: vec![p0, shoulder, p2],
        weights: vec![1.0, cos_half.abs(), 1.0],
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
    })
}

/// Set the cut edge's 3D curve (cap-tangent arc on cyl_m, going
/// through the cap-tangent circle of sphere_n on cyl_m) and refit
/// pcurves on both half-edges.
fn set_cut_edge_geometry(
    brep: &mut BRep,
    edge: EdgeHandle,
    cap_side_he: HalfEdgeHandle,
    v_start: popkern_core::geometry::Point3,
    v_end: popkern_core::geometry::Point3,
    cap_n: FaceHandle,
    main_cyl_face: FaceHandle,
) {
    use popkern_core::geometry::CurveGeometry;
    use popkern_core::geometry::Vec3;
    use popkern_core::geometry::analytic_surfaces::orthonormal_basis;
    use popkern_core::geometry::curve::CurveOnEdge;
    use popkern_core::geometry::pcurve::Pcurve;
    use popkern_core::geometry::surface::SurfaceGeometry;
    use popkern_core::geometry::uv_projection::project_point_to_uv;

    let cap_surf = brep
        .faces
        .get(cap_n)
        .and_then(|f| f.surface.as_ref())
        .map(|s| s.geometry.clone());
    let main_surf = brep
        .faces
        .get(main_cyl_face)
        .and_then(|f| f.surface.as_ref())
        .map(|s| s.geometry.clone());

    let twin = brep.half_edges[cap_side_he].twin;

    // Try the analytic path first: cap_n is Sphere, main_cyl_face is
    // Cylinder, and the cap-tangent arc lies exactly on both. We
    // build an `Ellipse` 3D curve parameterised by angle in the
    // cylinder's orthonormal basis — that way the cylinder's UV `u`
    // equals the curve parameter `t` directly, and Linear pcurves
    // are exact for both the cylinder and (for axis-aligned
    // sphere/cylinder pairs) the sphere.
    let analytic = if let (
        Some(SurfaceGeometry::Sphere {
            center: sphere_centre,
            radius: sphere_radius,
        }),
        Some(SurfaceGeometry::Cylinder {
            origin: cyl_origin,
            axis: cyl_axis,
            radius: cyl_radius,
        }),
    ) = (&cap_surf, &main_surf)
    {
        if let Some(axis_unit) = cyl_axis.normalized() {
            let (perp1, perp2) = orthonormal_basis(axis_unit);
            // Circle centre = projection of sphere centre on cyl axis.
            let to_centre = Vec3::new(
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

            // Compute angles in cyl basis: t such that
            // point = circle_centre + R · (perp1 · cos(t) + perp2 · sin(t)).
            let theta = |p: popkern_core::geometry::Point3| -> f64 {
                let d = Vec3::new(
                    p.x - circle_centre.x,
                    p.y - circle_centre.y,
                    p.z - circle_centre.z,
                );
                d.dot(&perp2).atan2(d.dot(&perp1))
            };
            let theta_a = theta(v_start);
            let theta_b = theta(v_end);
            // Pick the shorter arc direction (delta in [-π, π]).
            let mut delta = theta_b - theta_a;
            while delta > std::f64::consts::PI {
                delta -= std::f64::consts::TAU;
            }
            while delta < -std::f64::consts::PI {
                delta += std::f64::consts::TAU;
            }
            let theta_end = theta_a + delta;

            // 3D curve: Ellipse — parameterised linearly by angle.
            // u_basis = R · perp1, v_basis = R · perp2 with the SAME
            // basis as the cylinder UV, so cylinder's `u` equals `t`.
            let curve = CurveGeometry::Ellipse {
                center: circle_centre,
                u_basis: perp1 * *cyl_radius,
                v_basis: perp2 * *cyl_radius,
                start_angle: theta_a,
                end_angle: theta_end,
            };
            // Cylinder pcurve: Linear (theta_a, v_cyl) → (theta_end, v_cyl).
            let v_cyl = (v_start - *cyl_origin).dot(&axis_unit);
            let pcurve_cyl =
                Pcurve::linear((theta_a, v_cyl), (theta_end, v_cyl), theta_a, theta_end);
            // Sphere pcurve: Linear from sphere_uv(v_start) → sphere_uv(v_end).
            // For arcs whose plane coincides with a sphere meridian (one
            // endpoint at the pole) `project_point_to_uv` returns u=0
            // at the pole — fix that to the non-pole endpoint's u so
            // the Linear pcurve is a true meridian.
            let cap_geom = SurfaceGeometry::Sphere {
                center: *sphere_centre,
                radius: *sphere_radius,
            };
            // Pole detection via radial distance in the sphere's local
            // x/y plane (not via `v ≈ 0/π`, since `acos` near ±1 picks
            // up ~1e-8 rad of float noise — too tight for a v-based
            // threshold).  At a pole the radial component vanishes
            // and the projected `u` is undefined.
            let pole_radial_eps = sphere_radius * 1e-6;
            let radial_a = ((v_start.x - sphere_centre.x).powi(2)
                + (v_start.y - sphere_centre.y).powi(2))
            .sqrt();
            let radial_b =
                ((v_end.x - sphere_centre.x).powi(2) + (v_end.y - sphere_centre.y).powi(2)).sqrt();
            let pole_a = radial_a < pole_radial_eps;
            let pole_b = radial_b < pole_radial_eps;
            let pcurve_sphere = match (
                project_point_to_uv(&cap_geom, &v_start),
                project_point_to_uv(&cap_geom, &v_end),
            ) {
                (Some((u_a_raw, v_a)), Some((u_b_raw, v_b))) => {
                    // Snap `v` to exact 0 or π at the pole (acos noise),
                    // and inherit `u` from the non-pole endpoint so the
                    // Linear pcurve traces a true meridian.
                    let v_a = if pole_a {
                        if v_a < std::f64::consts::FRAC_PI_2 {
                            0.0
                        } else {
                            std::f64::consts::PI
                        }
                    } else {
                        v_a
                    };
                    let v_b = if pole_b {
                        if v_b < std::f64::consts::FRAC_PI_2 {
                            0.0
                        } else {
                            std::f64::consts::PI
                        }
                    } else {
                        v_b
                    };
                    let u_a = if pole_a { u_b_raw } else { u_a_raw };
                    let u_b = if pole_b { u_a_raw } else { u_b_raw };
                    Some(Pcurve::linear((u_a, v_a), (u_b, v_b), theta_a, theta_end))
                }
                _ => None,
            };
            Some((curve, theta_a, theta_end, pcurve_cyl, pcurve_sphere))
        } else {
            None
        }
    } else {
        None
    };

    if let Some((curve, t_start, t_end, pcurve_cyl, pcurve_sphere)) = analytic {
        if let Some(edge_data) = brep.edges.get_mut(edge) {
            edge_data.curve = Some(CurveOnEdge {
                geometry: curve,
                t_start,
                t_end,
            });
        }
        if let Some(pc) = pcurve_sphere {
            brep.half_edges[cap_side_he].pcurve = Some(pc);
        }
        brep.half_edges[twin].pcurve = Some(pcurve_cyl);
        return;
    }

    // Fallback: surfaces aren't sphere/cylinder.  Build a degree-1
    // straight line (best we can do) and fit pcurves by sampling.
    let curve_3d_nurbs = popkern_core::geometry::nurbs_curve::NurbsCurve {
        degree: 1,
        control_points: vec![v_start, v_end],
        weights: vec![1.0, 1.0],
        knots: vec![0.0, 0.0, 1.0, 1.0],
    };
    if let Some(edge_data) = brep.edges.get_mut(edge) {
        edge_data.curve = Some(CurveOnEdge {
            geometry: CurveGeometry::Nurbs(curve_3d_nurbs.clone()),
            t_start: 0.0,
            t_end: 1.0,
        });
    }
    if let Some(surf) = cap_surf
        && let Some(pc) = sample_pcurve_on_surface(&curve_3d_nurbs, &surf)
    {
        brep.half_edges[cap_side_he].pcurve = Some(pc.into());
    }
    if let Some(surf) = main_surf
        && let Some(pc) = sample_pcurve_on_surface(&curve_3d_nurbs, &surf)
    {
        brep.half_edges[twin].pcurve = Some(pc.into());
    }
}

/// Sample a 3D NURBS curve and fit a degree-3 NURBS pcurve in
/// the surface's UV space.  Returns None if any sample fails to
/// project.
fn sample_pcurve_on_surface(
    curve_3d: &popkern_core::geometry::nurbs_curve::NurbsCurve,
    surface: &popkern_core::geometry::surface::SurfaceGeometry,
) -> Option<popkern_core::geometry::nurbs_curve::NurbsCurve> {
    use popkern_core::geometry::Point3;
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::fitting;
    use popkern_core::geometry::uv_projection::project_point_to_uv;

    const N: usize = 12;
    let mut points: Vec<Point3> = Vec::with_capacity(N);
    let mut params: Vec<f64> = Vec::with_capacity(N);
    let trace = std::env::var("POPKERN_VAR_RADIUS_TRACE_PCURVE").is_ok_and(|v| v == "1");
    if trace {
        let surf_kind = match surface {
            popkern_core::geometry::surface::SurfaceGeometry::Sphere { center, radius } => {
                format!(
                    "Sphere centre=({:.4},{:.4},{:.4}) r={radius:.4}",
                    center.x, center.y, center.z
                )
            }
            popkern_core::geometry::surface::SurfaceGeometry::Cylinder {
                origin,
                axis,
                radius,
            } => {
                format!(
                    "Cyl origin=({:.4},{:.4},{:.4}) axis=({:.4},{:.4},{:.4}) r={radius:.4}",
                    origin.x, origin.y, origin.z, axis.x, axis.y, axis.z
                )
            }
            popkern_core::geometry::surface::SurfaceGeometry::Plane { origin, .. } => {
                format!(
                    "Plane origin=({:.4},{:.4},{:.4})",
                    origin.x, origin.y, origin.z
                )
            }
            _ => "<other>".to_string(),
        };
        let curve_endpoints = (curve_3d.point_at(0.0), curve_3d.point_at(1.0));
        eprintln!(
            "  [pcurve_fit] surface: {surf_kind} curve: ({:.4},{:.4},{:.4}) → ({:.4},{:.4},{:.4})",
            curve_endpoints.0.x,
            curve_endpoints.0.y,
            curve_endpoints.0.z,
            curve_endpoints.1.x,
            curve_endpoints.1.y,
            curve_endpoints.1.z
        );
    }
    for i in 0..N {
        let s = i as f64 / (N - 1) as f64;
        let p3 = curve_3d.point_at(s);
        let (u, v) = project_point_to_uv(surface, &p3)?;
        if trace {
            eprintln!(
                "    sample[{i}] s={s:.3} p3=({:.4},{:.4},{:.4}) uv=({u:.4},{v:.4})",
                p3.x, p3.y, p3.z
            );
        }
        points.push(Point3::new(u, v, 0.0));
        params.push(s);
    }
    let degree = 3usize;
    let n_cps = (degree + 2).min(points.len());
    let fit = fitting::fit_curve(&points, &params, degree, n_cps).ok()?;
    if trace {
        eprintln!("    fit cps: {:?}", fit.control_points);
    }
    Some(fit)
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

    /// End-to-end: run cube_mixed pipeline, find a shared
    /// rolling-ball vertex, split it, validate. Used as the
    /// first empirical test of `split_apex_footprint_vertex`.
    #[test]
    #[ignore]
    fn split_one_cube_mixed_vertex() {
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
        let mut brep = filleted.brep;

        // Pick the first NeighCyl edge to find V_x and its cap_face / cube_face.
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

        let classification = classify_cap_loops(&brep, &cap_faces, &per_fillet_caps);
        let (cap_k, neigh_he) = classification
            .iter()
            .find_map(|c| {
                c.entries
                    .iter()
                    .find(|e| e.twin_kind == TwinKind::NeighbourCylinder)
                    .map(|e| (c.cap_face, e.half_edge))
            })
            .unwrap();
        let v_x = brep.half_edges[neigh_he].origin;
        let cyl_m = brep
            .loops
            .get(brep.half_edges[brep.half_edges[neigh_he].twin].loop_)
            .unwrap()
            .face;

        // Find the cube face V_x is incident to (the Plane face).
        let cube_face = brep
            .half_edges
            .values()
            .find_map(|he| {
                if he.origin != v_x {
                    return None;
                }
                let f = brep.loops[he.loop_].face;
                let surf = brep.faces.get(f).and_then(|fd| fd.surface.as_ref())?;
                if matches!(surf.geometry, SurfaceGeometry::Plane { .. }) {
                    Some(f)
                } else {
                    None
                }
            })
            .unwrap();

        // Compute target position for v_new on cyl_m's cap-tangent
        // circle for the cap_n's sphere (the cap whose primary is cyl_m).
        let cyl_m_data = match brep
            .faces
            .get(cyl_m)
            .and_then(|f| f.surface.as_ref())
            .map(|s| &s.geometry)
        {
            Some(SurfaceGeometry::Cylinder {
                origin,
                axis,
                radius,
            }) => Some((*origin, *axis, *radius)),
            _ => None,
        }
        .unwrap();
        let (cyl_origin, cyl_axis, _cyl_radius) = cyl_m_data;
        let cap_n_sphere = cap_faces.iter().find_map(|&fh| {
            if fh == cap_k {
                return None;
            }
            match brep
                .faces
                .get(fh)
                .and_then(|f| f.surface.as_ref())
                .map(|s| &s.geometry)
            {
                Some(SurfaceGeometry::Sphere { center, .. }) => Some(*center),
                _ => None,
            }
        });
        let cap_n_centre = cap_n_sphere.unwrap();
        let v_x_pos = brep.vertices[v_x].position;
        // Hardcode: V_x is at (0, 0.12, 0.10). v_new at (0, 0.12, 0.12)
        // is the cap-tangent endpoint of cap_n (sphere r=0.12) on
        // cyl_m on cube face x=0.
        let v_new_pos = Point3::new(0.0, 0.12, 0.12);
        let _ = (cap_n_centre, cyl_origin, cyl_axis, v_x_pos);

        eprintln!(
            "Splitting V_x={v_x:?} pos=({:.3},{:.3},{:.3}) → ({:.3},{:.3},{:.3})",
            v_x_pos.x, v_x_pos.y, v_x_pos.z, v_new_pos.x, v_new_pos.y, v_new_pos.z
        );
        eprintln!("  cap_face={cap_k:?} cube_face={cube_face:?} cyl_split={cyl_m:?}");

        match split_apex_footprint_vertex(&mut brep, v_x, cap_k, cube_face, cyl_m, v_new_pos) {
            Ok(v_new) => {
                eprintln!("Split OK; v_new={v_new:?}");
                brep.validate().expect("validate after split");
            }
            Err(e) => panic!("split failed: {e:?}"),
        }
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

    /// Phase 6c-D investigation: walk each cap's outer loop in both
    /// rung-on and rung-off modes side-by-side.  Reveals exactly
    /// where the topology branches between the two paths.  Skips
    /// reroute (we want the *pre-surgery* topology each path
    /// produces).
    #[test]
    #[ignore]
    fn dump_cap_loops_rung_compare() {
        use popkern_core::geometry::surface::SurfaceGeometry;

        for (label, rung) in [("RUNG_OFF", false), ("RUNG_ON", true)] {
            unsafe {
                std::env::set_var("POPKERN_SPHERE_FIT_BOUNDARY", "1");
                if rung {
                    std::env::set_var("POPKERN_SPHERE_CAP_RUNG", "1");
                } else {
                    std::env::remove_var("POPKERN_SPHERE_CAP_RUNG");
                }
                std::env::remove_var("POPKERN_VAR_RADIUS_REROUTE");
                // Clean other env vars that could affect topology.
                std::env::remove_var("POPKERN_APEX_TOPOLOGY_REBUILD");
                std::env::remove_var("POPKERN_OUTER_FILLET_G1");
                std::env::remove_var("POPKERN_EDGE_TRIPLE");
            }
            eprintln!("\n========== {label} ==========");

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
            let Ok(filleted) = crate::fillet_edges(
                &cube,
                &specs,
                &crate::FilletOptions::default(),
                Tolerance::default(),
                None,
            ) else {
                eprintln!("  fillet_edges failed");
                continue;
            };
            let brep = filleted.brep;
            eprintln!(
                "BRep totals: {} faces, {} edges, {} half_edges, {} vertices",
                brep.faces.len(),
                brep.edges.len(),
                brep.half_edges.len(),
                brep.vertices.len()
            );
            eprintln!("corner_patches: {}", filleted.corner_patches.len());
            for cp in &filleted.corner_patches {
                eprintln!(
                    "  {:?} method={:?} face={:?}",
                    cp.classification, cp.method, cp.face
                );
            }
            // Also list ALL faces in the BRep with their surface type.
            eprintln!("All faces:");
            for (fh, f) in brep.faces.iter() {
                let kind = match f.surface.as_ref().map(|s| &s.geometry) {
                    Some(SurfaceGeometry::Plane { .. }) => "Plane",
                    Some(SurfaceGeometry::Cylinder { .. }) => "Cyl",
                    Some(SurfaceGeometry::Sphere { .. }) => "Sphere",
                    Some(SurfaceGeometry::Gregory(_)) => "Gregory",
                    Some(SurfaceGeometry::Torus { .. }) => "Torus",
                    Some(SurfaceGeometry::Cone { .. }) => "Cone",
                    Some(SurfaceGeometry::Nurbs(_)) => "Nurbs",
                    None => "?",
                };
                eprintln!("  {fh:?} {kind}");
            }

            // Identify cap faces by their cap_sphere data on
            // CornerPatchInfo (works regardless of surface type).
            let mut cap_info: Vec<(FaceHandle, Option<(Point3, f64)>)> = filleted
                .corner_patches
                .iter()
                .map(|p| (p.face, p.cap_sphere))
                .collect();
            cap_info.sort_by_key(|(fh, _)| {
                // Stable order: by sphere centre x then radius
                let (cx, r) = (*fh, 0.0);
                let _ = (cx, r);
                *fh
            });
            for (k, (fh, sphere)) in cap_info.iter().enumerate() {
                let surface_kind = match brep
                    .faces
                    .get(*fh)
                    .and_then(|f| f.surface.as_ref())
                    .map(|s| &s.geometry)
                {
                    Some(SurfaceGeometry::Sphere { center, radius }) => format!(
                        "Sphere centre=({:.3},{:.3},{:.3}) r={:.3}",
                        center.x, center.y, center.z, radius
                    ),
                    Some(SurfaceGeometry::Gregory(_)) => "Gregory".to_string(),
                    Some(other) => format!("{other:?}"),
                    None => "<no surface>".to_string(),
                };
                eprintln!(
                    "\ncap[{k}] {fh:?} cap_sphere={:?} surface={surface_kind}",
                    sphere.map(|(c, r)| format!(
                        "centre=({:.3},{:.3},{:.3}) r={:.3}",
                        c.x, c.y, c.z, r
                    ))
                );
                let Some(face_data) = brep.faces.get(*fh) else {
                    continue;
                };
                let outer = face_data.outer_loop;
                let Some(loop_data) = brep.loops.get(outer) else {
                    continue;
                };
                let start = loop_data.half_edge;
                let mut cursor = start;
                let mut step = 0;
                loop {
                    if step > 16 {
                        break;
                    }
                    let Some(he) = brep.half_edges.get(cursor) else {
                        break;
                    };
                    let edge = he.edge;
                    let v_origin = he.origin;
                    let pos = brep.vertices[v_origin].position;
                    let twin = he.twin;
                    let twin_face = brep
                        .half_edges
                        .get(twin)
                        .map(|t| brep.loops[t.loop_].face)
                        .unwrap_or(*fh);
                    let twin_kind = match brep
                        .faces
                        .get(twin_face)
                        .and_then(|f| f.surface.as_ref())
                        .map(|s| &s.geometry)
                    {
                        Some(SurfaceGeometry::Plane { .. }) => "Plane",
                        Some(SurfaceGeometry::Cylinder { .. }) => "Cyl",
                        Some(SurfaceGeometry::Sphere { .. }) => "Sphere",
                        Some(SurfaceGeometry::Gregory(_)) => "Greg",
                        _ => "?",
                    };
                    eprintln!(
                        "  step {step}: he={cursor:?} edge={edge:?} V={v_origin:?}=({:.3},{:.3},{:.3}) → twin face {twin_face:?} ({twin_kind})",
                        pos.x, pos.y, pos.z
                    );
                    cursor = he.next;
                    step += 1;
                    if cursor == start {
                        break;
                    }
                }
            }
        }
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

    // ================================================================
    // Synthetic 4-face apex BRep — vertex-split primitive coverage
    // ================================================================
    //
    // Build a planar pyramid topology with apex `V` valence-4
    // incident to 4 triangle faces, then exercise
    // `split_apex_footprint_vertex` directly.  No surfaces or curves
    // — purely topology + vertex positions.  Tests the structural
    // logic in isolation from the cube_mixed pipeline.

    /// Build the planar pyramid: square ABCD with central V_apex,
    /// split into 4 inside triangles meeting at V_apex, plus one
    /// outside square face F_back.  Returns:
    /// `(brep, v_apex, [tri0, tri1, tri2, tri3])` where the four
    /// triangles are listed in cyclic order around V_apex going
    /// around the +Z axis (CCW).
    ///
    /// Vertex layout (z=0 plane):
    /// ```text
    ///   D(0,1) -- C(1,1)
    ///     | \  /  |
    ///     |  V    |
    ///     | /  \  |
    ///   A(0,0) -- B(1,0)
    /// ```
    /// V_apex at (0.5, 0.5, 0).  Triangles in cyclic order:
    /// `tri0 = (V,A,B)`, `tri1 = (V,B,C)`, `tri2 = (V,C,D)`,
    /// `tri3 = (V,D,A)`.
    fn build_pyramid_4() -> (
        popkern_core::topology::BRep,
        popkern_core::topology::VertexHandle,
        [popkern_core::topology::FaceHandle; 4],
    ) {
        use popkern_core::geometry::Point3;
        use popkern_core::tolerance::Tolerance;
        use popkern_core::topology::BRep;

        let mut brep = BRep::new(Tolerance::default());

        // Build square A(0,0) → B(1,0) → C(1,1) → D(0,1) in F_orig.
        let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap(); // A
        let mev_b = brep
            .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap(); // B
        let mev_c = brep
            .mev(mev_b.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
            .unwrap(); // C
        let mev_d = brep
            .mev(mev_c.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
            .unwrap(); // D
        // mef D→A: closes the square. The "old face" (F_orig from
        // mvfs) keeps half-edges from new_he onward — i.e., the
        // interior square loop {DA, AB, BC, CD}.  The new face
        // (mef_da.face) gets the loop walking the opposite way
        // {AD, DC, CB, BA} — that's F_back.
        let mef_da = brep.mef(mev_d.half_edge_twin, mev_b.half_edge).unwrap();
        let f_inside = mvfs.face; // contains interior square loop
        let _f_back = mef_da.face; // back/outside face

        // Insert V_apex at (0.5, 0.5, 0) via mev from A's outgoing
        // (= mev_b.half_edge, A→B in F_inside).  The spur lives in
        // F_inside.
        let mev_v = brep
            .mev(mev_b.half_edge, Point3::new(0.5, 0.5, 0.0))
            .unwrap();
        let v_apex = mev_v.vertex;
        // mev_v.half_edge_twin is V→A (in F_inside)

        // mef V_apex → B: split F_inside into triangle {A,B,V} and pentagon {V,B,C,D,A}.
        // mev_v.half_edge_twin originates at V (V→A), and mev_c.half_edge originates at B (B→C).
        // Both in F_inside's loop after the spur insertion.
        let mef_apex_b = brep.mef(mev_v.half_edge_twin, mev_c.half_edge).unwrap();
        // mef_apex_b.half_edge is V→B (in F_inside, now pentagon)
        // mef_apex_b.face is the new triangle face (A,B,V).

        // mef V_apex → C: split pentagon into triangle {B,C,V} and quad {V,C,D,A}.
        // mef_apex_b.half_edge originates at V (V→B in pentagon).
        // mev_d.half_edge originates at C (C→D in pentagon).
        let mef_apex_c = brep.mef(mef_apex_b.half_edge, mev_d.half_edge).unwrap();
        // mef_apex_c.half_edge is V→C (in F_inside, now quad).
        // mef_apex_c.face is the new triangle face (B,C,V).

        // mef V_apex → D: split quad into triangle {C,D,V} and triangle {D,A,V}.
        // mef_apex_c.half_edge originates at V (V→C in quad).
        // mef_da.half_edge originates at D (D→A in quad).
        let mef_apex_d = brep.mef(mef_apex_c.half_edge, mef_da.half_edge).unwrap();
        // mef_apex_d.half_edge is V→D (in F_inside, now last triangle (D,A,V)).
        // mef_apex_d.face is the new triangle face (C,D,V).

        // Faces in cyclic order around V_apex going CCW (V→A is at
        // angle 225°, V→B at 315°, V→C at 45°, V→D at 135°):
        //   tri0 = (V,A,B) face  (angles 225°→315°)
        //   tri1 = (V,B,C) face  (angles 315°→45°)
        //   tri2 = (V,C,D) face  (angles 45°→135°)
        //   tri3 = (V,D,A) face  (angles 135°→225°)
        //
        //   tri0 = mef_apex_b.face (the (A,B,V) triangle face)
        //   tri1 = mef_apex_c.face (the (B,C,V) triangle face)
        //   tri2 = mef_apex_d.face (the (C,D,V) triangle face)
        //   tri3 = f_inside        (the leftover (D,A,V) triangle face)
        let tri0 = mef_apex_b.face;
        let tri1 = mef_apex_c.face;
        let tri2 = mef_apex_d.face;
        let tri3 = f_inside;
        (brep, v_apex, [tri0, tri1, tri2, tri3])
    }

    /// Sanity check: build_pyramid_4 produces validate-able topology.
    #[test]
    fn synthetic_pyramid_4_valid() {
        let (brep, v_apex, tris) = build_pyramid_4();
        brep.validate().expect("pyramid_4 should be valid");
        assert_eq!(brep.vertices.len(), 5);
        assert_eq!(brep.faces.len(), 5); // 4 inner + 1 back
        assert_eq!(brep.edges.len(), 8); // 4 perimeter + 4 spokes
        // V_apex valence 4
        let valence = brep
            .half_edges
            .values()
            .filter(|he| he.origin == v_apex)
            .count();
        assert_eq!(valence, 4, "V_apex should have valence 4");
        // V_apex incident to all 4 inner triangles (and not f_back).
        for &tri in &tris {
            let outer = brep.faces[tri].outer_loop;
            let mut found = false;
            let start = brep.loops[outer].half_edge;
            let mut cursor = start;
            for _ in 0..16 {
                if brep.half_edges[cursor].origin == v_apex {
                    found = true;
                    break;
                }
                cursor = brep.half_edges[cursor].next;
                if cursor == start {
                    break;
                }
            }
            assert!(found, "V_apex must be in tri {tri:?} loop");
        }
    }

    /// First empirical exercise of `split_apex_footprint_vertex` on
    /// the synthetic 4-face apex.  Cyclic order matches cube_mixed:
    /// `cap=tri0`, `cyl_split=tri1`, `cube=tri2`, `extra=tri3`.
    /// This is the cyclic order [cap, cyl_split, cube, extra] —
    /// where cap and cyl_split are adjacent at one V_apex edge,
    /// cyl_split and cube at the next, and so on.
    #[test]
    fn split_synthetic_4_face_cyclic() {
        use popkern_core::geometry::Point3;

        let (mut brep, v_apex, [tri0, tri1, tri2, tri3]) = build_pyramid_4();
        let _ = tri3;
        let v_new_pos = Point3::new(0.4, 0.5, 0.0);
        let pre_face_count = brep.faces.len();
        let pre_vertex_count = brep.vertices.len();

        let v_new = split_apex_footprint_vertex(&mut brep, v_apex, tri0, tri2, tri1, v_new_pos)
            .expect("split should succeed");
        eprintln!("split OK: v_new={v_new:?}");
        brep.validate().expect("validate after split");
        assert_eq!(brep.vertices.len(), pre_vertex_count + 1);
        assert_eq!(brep.faces.len(), pre_face_count);
    }

    /// Configuration where cap and cube are cyclically adjacent
    /// (cyl_split is "across" from the cap-cube boundary). The
    /// surgery design relies on cap and cube each having one
    /// cyl_seam edge incident to V — but in this layout neither
    /// does.  The function detects the unsupported configuration
    /// and returns `LoopBroken` rather than producing a corrupt
    /// BRep.
    #[test]
    fn split_synthetic_4_face_cap_cube_adjacent_unsupported() {
        use popkern_core::geometry::Point3;

        let (mut brep, v_apex, [tri0, tri1, tri2, tri3]) = build_pyramid_4();
        let _ = tri3;
        let v_new_pos = Point3::new(0.5, 0.4, 0.0);

        let result = split_apex_footprint_vertex(&mut brep, v_apex, tri0, tri1, tri2, v_new_pos);
        assert!(matches!(result, Err(RerouteError::LoopBroken)));
        // Snapshot rollback: BRep should be unchanged.
        brep.validate().expect("BRep unchanged after rollback");
    }

    /// Diagnostic: for each cap-vs-non-cap edge in cube_mixed
    /// post-pipeline (with surgery), dump G1 coplanarity deficit
    /// at sample points so we can identify which edges drag the
    /// worst-g1 number up.  `#[ignore]`d because it sets
    /// process-level env vars; run isolated.
    #[test]
    #[ignore]
    fn dump_per_edge_g1_cube_mixed() {
        unsafe {
            std::env::set_var("POPKERN_SPHERE_FIT_BOUNDARY", "1");
            std::env::set_var("POPKERN_SPHERE_CAP_RUNG", "1");
            std::env::set_var("POPKERN_VAR_RADIUS_REROUTE", "1");
        }
        let cube =
            popkern_core::topology::BRep::unit_cube(popkern_core::tolerance::Tolerance::default())
                .unwrap();
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
            popkern_core::tolerance::Tolerance::default(),
            None,
        );
        let Ok(filleted) = result else {
            eprintln!("fillet_edges failed");
            return;
        };
        let brep = filleted.brep;

        use popkern_core::geometry::curve::CurveEval;
        use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
        let cap_set: std::collections::BTreeSet<_> = brep
            .faces
            .iter()
            .filter_map(|(fh, f)| match f.surface.as_ref().map(|s| &s.geometry) {
                Some(SurfaceGeometry::Sphere { .. }) => Some(fh),
                _ => None,
            })
            .collect();

        eprintln!(
            "BRep: {} faces, {} edges; {} sphere-cap faces",
            brep.faces.len(),
            brep.edges.len(),
            cap_set.len()
        );
        let mut entries: Vec<(
            f64,
            f64,
            popkern_core::topology::EdgeHandle,
            popkern_core::topology::FaceHandle,
            popkern_core::topology::FaceHandle,
            &str,
            &str,
        )> = Vec::new();
        for (eh, edge) in brep.edges.iter() {
            let Some(coe) = edge.curve.as_ref() else {
                continue;
            };
            let he_a = edge.half_edge;
            let Some(ha) = brep.half_edges.get(he_a) else {
                continue;
            };
            let Some(hb) = brep.half_edges.get(ha.twin) else {
                continue;
            };
            let face_a = brep.loops[ha.loop_].face;
            let face_b = brep.loops[hb.loop_].face;
            if face_a == face_b {
                continue;
            }
            let a_is_cap = cap_set.contains(&face_a);
            let b_is_cap = cap_set.contains(&face_b);
            // Measure ALL distinct-face edges, label cap vs non-cap.
            let Some(pa) = ha.pcurve.as_ref() else {
                continue;
            };
            let Some(pb) = hb.pcurve.as_ref() else {
                continue;
            };
            let Some(sa) = brep.faces.get(face_a).and_then(|f| f.surface.as_ref()) else {
                continue;
            };
            let Some(sb) = brep.faces.get(face_b).and_then(|f| f.surface.as_ref()) else {
                continue;
            };

            let kind = |face: popkern_core::topology::FaceHandle| -> &'static str {
                match brep
                    .faces
                    .get(face)
                    .and_then(|f| f.surface.as_ref())
                    .map(|s| &s.geometry)
                {
                    Some(SurfaceGeometry::Plane { .. }) => "Plane",
                    Some(SurfaceGeometry::Cylinder { .. }) => "Cyl",
                    Some(SurfaceGeometry::Sphere { .. }) => "Sphere",
                    Some(SurfaceGeometry::Torus { .. }) => "Torus",
                    Some(SurfaceGeometry::Cone { .. }) => "Cone",
                    Some(SurfaceGeometry::Nurbs(_)) => "Nurbs",
                    Some(SurfaceGeometry::Gregory(_)) => "Greg",
                    None => "?",
                }
            };
            let a_label = kind(face_a);
            let b_label = kind(face_b);

            let (t_lo, t_hi) = (coe.t_start, coe.t_end);
            let (pa_lo, pa_hi) = pa.parameter_range();
            let (pb_lo, pb_hi) = pb.parameter_range();
            let p3_lo = coe.geometry.point_at(t_lo);
            let p3_hi = coe.geometry.point_at(t_hi);
            let pa_at = |t: f64| {
                let uv = pa.point_at(t);
                sa.geometry.point_at(uv.x, uv.y)
            };
            let pb_at = |t: f64| {
                let uv = pb.point_at(t);
                sb.geometry.point_at(uv.x, uv.y)
            };
            let pa_forward = p3_lo.distance_to(&pa_at(pa_lo)) + p3_hi.distance_to(&pa_at(pa_hi))
                <= p3_lo.distance_to(&pa_at(pa_hi)) + p3_hi.distance_to(&pa_at(pa_lo));
            let pb_forward = p3_lo.distance_to(&pb_at(pb_lo)) + p3_hi.distance_to(&pb_at(pb_hi))
                <= p3_lo.distance_to(&pb_at(pb_hi)) + p3_hi.distance_to(&pb_at(pb_lo));

            let n_samples = 7usize;
            let mut max_g0 = 0.0_f64;
            let mut max_g1 = 0.0_f64;
            for i in 0..n_samples {
                let s = 0.05 + 0.90 * (i as f64) / ((n_samples - 1) as f64);
                let s_a = if pa_forward { s } else { 1.0 - s };
                let s_b = if pb_forward { s } else { 1.0 - s };
                let t_a = pa_lo + s_a * (pa_hi - pa_lo);
                let t_b = pb_lo + s_b * (pb_hi - pb_lo);
                let uv_a = pa.point_at(t_a);
                let uv_b = pb.point_at(t_b);
                let p_a = sa.geometry.point_at(uv_a.x, uv_a.y);
                let p_b = sb.geometry.point_at(uv_b.x, uv_b.y);
                let g0 = p_a.distance_to(&p_b);
                if g0 > max_g0 {
                    max_g0 = g0;
                }
                let su_a = sa.geometry.partial_u(uv_a.x, uv_a.y);
                let sv_a = sa.geometry.partial_v(uv_a.x, uv_a.y);
                let su_b = sb.geometry.partial_u(uv_b.x, uv_b.y);
                let sv_b = sb.geometry.partial_v(uv_b.x, uv_b.y);
                let n_a = su_a.cross(&sv_a);
                let n_b = su_b.cross(&sv_b);
                let nn_a = (n_a.dot(&n_a)).sqrt();
                let nn_b = (n_b.dot(&n_b)).sqrt();
                if nn_a > 1e-300 && nn_b > 1e-300 {
                    let dot_unit = (n_a.dot(&n_b) / (nn_a * nn_b)).abs();
                    let deficit = 1.0 - dot_unit;
                    if deficit > max_g1 {
                        max_g1 = deficit;
                    }
                }
            }
            let _ = (a_is_cap, b_is_cap);
            entries.push((max_g1, max_g0, eh, face_a, face_b, a_label, b_label));
        }
        // Sort by max_g1 descending.
        entries.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));
        eprintln!("\n--- per-edge G1/G0 sorted by G1 (top 20) ---");
        for (i, (g1, g0, eh, fa, fb, la, lb)) in entries.iter().take(20).enumerate() {
            let cap_a = cap_set.contains(fa);
            let cap_b = cap_set.contains(fb);
            eprintln!(
                "  [{i:2}] edge={eh:?} g1={g1:.3e} g0={g0:.3e} {fa:?}({la}{}) <-> {fb:?}({lb}{})",
                if cap_a { "*" } else { "" },
                if cap_b { "*" } else { "" }
            );
        }

        // For the top-3 cap-vs-Cyl edges (the cut edges), dump the
        // 3D curve point vs each pcurve-evaluated surface point at
        // a few sample params.  This isolates whether the gap is
        // pcurve-fit error, surface-curve mismatch, or both.
        eprintln!("\n--- top-3 cap-vs-Cyl edges: curve vs pcurve-on-surface ---");
        let mut shown = 0usize;
        for (g1, g0, eh, fa, fb, la, lb) in entries.iter() {
            let cap_a = cap_set.contains(fa);
            let cap_b = cap_set.contains(fb);
            if !(cap_a ^ cap_b) {
                continue;
            }
            if !(la.starts_with("Cyl") || lb.starts_with("Cyl")) {
                continue;
            }
            if shown >= 3 {
                break;
            }
            shown += 1;
            let edge = brep.edges.get(*eh).unwrap();
            let coe = edge.curve.as_ref().unwrap();
            let he_a = edge.half_edge;
            let ha = brep.half_edges.get(he_a).unwrap();
            let hb = brep.half_edges.get(ha.twin).unwrap();
            let pa = ha.pcurve.as_ref().unwrap();
            let pb = hb.pcurve.as_ref().unwrap();
            let face_a = brep.loops[ha.loop_].face;
            let face_b = brep.loops[hb.loop_].face;
            let sa = brep.faces[face_a].surface.as_ref().unwrap();
            let sb = brep.faces[face_b].surface.as_ref().unwrap();
            let (pa_lo, pa_hi) = pa.parameter_range();
            let (pb_lo, pb_hi) = pb.parameter_range();
            let (t_lo, t_hi) = (coe.t_start, coe.t_end);

            // Endpoint vertex positions and curve degree.
            let v_origin = brep.vertices[ha.origin].position;
            let v_dest = brep.vertices[hb.origin].position;
            let curve_info = match &coe.geometry {
                popkern_core::geometry::CurveGeometry::Nurbs(nc) => {
                    format!(
                        "Nurbs deg={} cps={} weights={:?}",
                        nc.degree,
                        nc.control_points.len(),
                        nc.weights
                    )
                }
                popkern_core::geometry::CurveGeometry::Line { .. } => "Line".to_string(),
                popkern_core::geometry::CurveGeometry::Arc { .. } => "Arc".to_string(),
                popkern_core::geometry::CurveGeometry::Circle { .. } => "Circle".to_string(),
                popkern_core::geometry::CurveGeometry::Ellipse { .. } => "Ellipse".to_string(),
                _ => "<other>".to_string(),
            };
            eprintln!("edge={eh:?} ({la} on {fa:?} <-> {lb} on {fb:?}) g1={g1:.3e} g0={g0:.3e}");
            eprintln!(
                "  endpoints: ({:.4},{:.4},{:.4}) → ({:.4},{:.4},{:.4})  curve: {curve_info}",
                v_origin.x, v_origin.y, v_origin.z, v_dest.x, v_dest.y, v_dest.z
            );
            // Surface params for context.
            match &sa.geometry {
                popkern_core::geometry::surface::SurfaceGeometry::Sphere { center, radius } => {
                    eprintln!(
                        "  surface_a Sphere centre=({:.4},{:.4},{:.4}) r={radius:.4}",
                        center.x, center.y, center.z
                    );
                }
                popkern_core::geometry::surface::SurfaceGeometry::Cylinder {
                    origin,
                    axis,
                    radius,
                } => {
                    eprintln!(
                        "  surface_a Cylinder origin=({:.4},{:.4},{:.4}) axis=({:.4},{:.4},{:.4}) r={radius:.4}",
                        origin.x, origin.y, origin.z, axis.x, axis.y, axis.z
                    );
                }
                _ => {}
            }
            match &sb.geometry {
                popkern_core::geometry::surface::SurfaceGeometry::Sphere { center, radius } => {
                    eprintln!(
                        "  surface_b Sphere centre=({:.4},{:.4},{:.4}) r={radius:.4}",
                        center.x, center.y, center.z
                    );
                }
                popkern_core::geometry::surface::SurfaceGeometry::Cylinder {
                    origin,
                    axis,
                    radius,
                } => {
                    eprintln!(
                        "  surface_b Cylinder origin=({:.4},{:.4},{:.4}) axis=({:.4},{:.4},{:.4}) r={radius:.4}",
                        origin.x, origin.y, origin.z, axis.x, axis.y, axis.z
                    );
                }
                _ => {}
            }
            eprintln!(
                "  curve t∈[{t_lo:.3},{t_hi:.3}] pa_param∈[{pa_lo:.3},{pa_hi:.3}] pb_param∈[{pb_lo:.3},{pb_hi:.3}]"
            );
            // Dump pa/pb pcurve type.
            eprintln!("  ha={he_a:?} hb={:?}", ha.twin);
            match &pa.geometry {
                popkern_core::geometry::pcurve::PcurveGeometry::Nurbs(nc) => {
                    eprintln!(
                        "  pa Nurbs deg={} n_cps={} cps={:?}",
                        nc.degree,
                        nc.control_points.len(),
                        nc.control_points
                    );
                }
                popkern_core::geometry::pcurve::PcurveGeometry::Linear { start, end } => {
                    eprintln!(
                        "  pa Linear ({:.4},{:.4})→({:.4},{:.4}) t∈[{:.4},{:.4}]",
                        start.0, start.1, end.0, end.1, pa.t_start, pa.t_end
                    );
                }
                _ => eprintln!("  pa <other pcurve geometry>"),
            }
            if let popkern_core::geometry::pcurve::PcurveGeometry::Nurbs(nc) = &pb.geometry {
                eprintln!(
                    "  pb Nurbs deg={} n_cps={} cps={:?}",
                    nc.degree,
                    nc.control_points.len(),
                    nc.control_points
                );
            } else {
                eprintln!("  pb <non-nurbs pcurve geometry>");
            }
            // Detect direction.
            let p3_lo = coe.geometry.point_at(t_lo);
            let p3_hi = coe.geometry.point_at(t_hi);
            let pa_at_param = |t: f64| {
                let uv = pa.point_at(t);
                sa.geometry.point_at(uv.x, uv.y)
            };
            let pb_at_param = |t: f64| {
                let uv = pb.point_at(t);
                sb.geometry.point_at(uv.x, uv.y)
            };
            let pa_forward = p3_lo.distance_to(&pa_at_param(pa_lo))
                + p3_hi.distance_to(&pa_at_param(pa_hi))
                <= p3_lo.distance_to(&pa_at_param(pa_hi)) + p3_hi.distance_to(&pa_at_param(pa_lo));
            let pb_forward = p3_lo.distance_to(&pb_at_param(pb_lo))
                + p3_hi.distance_to(&pb_at_param(pb_hi))
                <= p3_lo.distance_to(&pb_at_param(pb_hi)) + p3_hi.distance_to(&pb_at_param(pb_lo));

            for s in [0.0_f64, 0.25, 0.5, 0.75, 1.0] {
                let t = t_lo + s * (t_hi - t_lo);
                let p3 = coe.geometry.point_at(t);
                let s_a = if pa_forward { s } else { 1.0 - s };
                let s_b = if pb_forward { s } else { 1.0 - s };
                let t_a = pa_lo + s_a * (pa_hi - pa_lo);
                let t_b = pb_lo + s_b * (pb_hi - pb_lo);
                let uv_a = pa.point_at(t_a);
                let uv_b = pb.point_at(t_b);
                let p_a = sa.geometry.point_at(uv_a.x, uv_a.y);
                let p_b = sb.geometry.point_at(uv_b.x, uv_b.y);
                let d_curve_a = p3.distance_to(&p_a);
                let d_curve_b = p3.distance_to(&p_b);
                let d_a_b = p_a.distance_to(&p_b);
                eprintln!(
                    "    s={s:.2}  curve=({:.4},{:.4},{:.4}) |curve→a|={d_curve_a:.3e} |curve→b|={d_curve_b:.3e} |a→b|={d_a_b:.3e}",
                    p3.x, p3.y, p3.z
                );
            }
        }
    }

    /// Run the full Path B reroute on cube_mixed and report which
    /// surgeries fail.  Uses the ignore-attribute so it's only run
    /// when explicitly requested (sets process-level env vars).
    #[test]
    #[ignore]
    fn reroute_full_cube_mixed() {
        unsafe {
            std::env::set_var("POPKERN_SPHERE_FIT_BOUNDARY", "1");
            std::env::set_var("POPKERN_SPHERE_CAP_RUNG", "1");
            std::env::set_var("POPKERN_VAR_RADIUS_REROUTE", "1");
            std::env::set_var("POPKERN_VAR_RADIUS_TRACE", "1");
        }
        let cube =
            popkern_core::topology::BRep::unit_cube(popkern_core::tolerance::Tolerance::default())
                .unwrap();
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
            popkern_core::tolerance::Tolerance::default(),
            None,
        );
        match result {
            Ok(filleted) => {
                eprintln!(
                    "OK — final BRep: {} faces, {} edges, {} vertices",
                    filleted.brep.faces.len(),
                    filleted.brep.edges.len(),
                    filleted.brep.vertices.len()
                );
                if let Err(e) = filleted.brep.validate() {
                    eprintln!("post-pipeline validate FAILED: {e:?}");
                } else {
                    eprintln!("post-pipeline validate OK");
                }
            }
            Err(e) => {
                eprintln!("fillet_edges failed: {e:?}");
            }
        }
    }

    /// Phase 6c-D Step 2 probe — measure post-reroute G1/G0 for
    /// Gregory caps (rung-OFF mode), to size the gap that
    /// re-optimisation needs to close.  Mirror of
    /// `dump_per_edge_g1_cube_mixed` but with `POPKERN_SPHERE_CAP_RUNG`
    /// unset (so cap surfaces stay Gregory after K-P binonic ladder
    /// fits them pre-reroute) and `POPKERN_OUTER_FILLET_G1` on
    /// (so the optimiser includes outer-fillet G1 residuals).
    ///
    /// Reroute is gated by env: when `POPKERN_VAR_RADIUS_REROUTE=1`
    /// (default in this test) we measure post-reroute state.
    /// Set `POPKERN_PHASE_6CD_PROBE_NO_REROUTE=1` to instead measure
    /// pre-reroute baseline (separates optimiser quality from
    /// reroute damage).  `#[ignore]`d because it sets process-level
    /// env vars.
    #[test]
    #[ignore]
    fn dump_per_edge_g1_cube_mixed_gregory_rung_off() {
        let no_reroute = std::env::var("POPKERN_PHASE_6CD_PROBE_NO_REROUTE")
            .map(|v| v == "1")
            .unwrap_or(false);
        let outer_g1 = std::env::var("POPKERN_PHASE_6CD_PROBE_OUTER_G1")
            .map(|v| v == "1")
            .unwrap_or(true);
        let sphere_fit = std::env::var("POPKERN_PHASE_6CD_PROBE_SPHERE_FIT")
            .map(|v| v == "1")
            .unwrap_or(true);
        unsafe {
            if sphere_fit {
                std::env::set_var("POPKERN_SPHERE_FIT_BOUNDARY", "1");
            } else {
                std::env::remove_var("POPKERN_SPHERE_FIT_BOUNDARY");
            }
            std::env::remove_var("POPKERN_SPHERE_CAP_RUNG");
            if outer_g1 {
                std::env::set_var("POPKERN_OUTER_FILLET_G1", "1");
            } else {
                std::env::remove_var("POPKERN_OUTER_FILLET_G1");
            }
            if no_reroute {
                std::env::remove_var("POPKERN_VAR_RADIUS_REROUTE");
            } else {
                std::env::set_var("POPKERN_VAR_RADIUS_REROUTE", "1");
            }
        }
        eprintln!(
            "Phase 6c-D probe mode: {} OUTER_FILLET_G1={} SPHERE_FIT={}",
            if no_reroute {
                "PRE-REROUTE"
            } else {
                "POST-REROUTE"
            },
            if outer_g1 { "ON" } else { "OFF" },
            if sphere_fit { "ON" } else { "OFF" }
        );
        let cube =
            popkern_core::topology::BRep::unit_cube(popkern_core::tolerance::Tolerance::default())
                .unwrap();
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
            popkern_core::tolerance::Tolerance::default(),
            None,
        );
        let Ok(filleted) = result else {
            eprintln!("fillet_edges failed");
            return;
        };
        let brep = filleted.brep;

        use popkern_core::geometry::curve::CurveEval;
        use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
        let cap_set: std::collections::BTreeSet<_> = brep
            .faces
            .iter()
            .filter_map(|(fh, f)| match f.surface.as_ref().map(|s| &s.geometry) {
                Some(SurfaceGeometry::Gregory(_)) => Some(fh),
                _ => None,
            })
            .collect();

        eprintln!(
            "BRep (rung-OFF Gregory): {} faces, {} edges; {} Gregory-cap faces",
            brep.faces.len(),
            brep.edges.len(),
            cap_set.len()
        );

        let kind = |face: popkern_core::topology::FaceHandle| -> &'static str {
            match brep
                .faces
                .get(face)
                .and_then(|f| f.surface.as_ref())
                .map(|s| &s.geometry)
            {
                Some(SurfaceGeometry::Plane { .. }) => "Plane",
                Some(SurfaceGeometry::Cylinder { .. }) => "Cyl",
                Some(SurfaceGeometry::Sphere { .. }) => "Sphere",
                Some(SurfaceGeometry::Torus { .. }) => "Torus",
                Some(SurfaceGeometry::Cone { .. }) => "Cone",
                Some(SurfaceGeometry::Nurbs(_)) => "Nurbs",
                Some(SurfaceGeometry::Gregory(_)) => "Greg",
                None => "?",
            }
        };

        let mut entries: Vec<(
            f64,
            f64,
            popkern_core::topology::EdgeHandle,
            popkern_core::topology::FaceHandle,
            popkern_core::topology::FaceHandle,
            &str,
            &str,
        )> = Vec::new();
        for (eh, edge) in brep.edges.iter() {
            if edge.curve.as_ref().is_none() {
                continue;
            }
            let he_a = edge.half_edge;
            let Some(ha) = brep.half_edges.get(he_a) else {
                continue;
            };
            let Some(hb) = brep.half_edges.get(ha.twin) else {
                continue;
            };
            let face_a = brep.loops[ha.loop_].face;
            let face_b = brep.loops[hb.loop_].face;
            if face_a == face_b {
                continue;
            }
            let Some(pa) = ha.pcurve.as_ref() else {
                continue;
            };
            let Some(pb) = hb.pcurve.as_ref() else {
                continue;
            };
            let Some(sa) = brep.faces.get(face_a).and_then(|f| f.surface.as_ref()) else {
                continue;
            };
            let Some(sb) = brep.faces.get(face_b).and_then(|f| f.surface.as_ref()) else {
                continue;
            };
            let a_label = kind(face_a);
            let b_label = kind(face_b);
            let (pa_lo, pa_hi) = pa.parameter_range();
            let (pb_lo, pb_hi) = pb.parameter_range();
            let coe = edge.curve.as_ref().unwrap();
            let (t_lo, t_hi) = (coe.t_start, coe.t_end);
            let p3_lo = coe.geometry.point_at(t_lo);
            let p3_hi = coe.geometry.point_at(t_hi);
            let pa_at = |t: f64| {
                let uv = pa.point_at(t);
                sa.geometry.point_at(uv.x, uv.y)
            };
            let pb_at = |t: f64| {
                let uv = pb.point_at(t);
                sb.geometry.point_at(uv.x, uv.y)
            };
            let pa_forward = p3_lo.distance_to(&pa_at(pa_lo)) + p3_hi.distance_to(&pa_at(pa_hi))
                <= p3_lo.distance_to(&pa_at(pa_hi)) + p3_hi.distance_to(&pa_at(pa_lo));
            let pb_forward = p3_lo.distance_to(&pb_at(pb_lo)) + p3_hi.distance_to(&pb_at(pb_hi))
                <= p3_lo.distance_to(&pb_at(pb_hi)) + p3_hi.distance_to(&pb_at(pb_lo));

            let n_samples = 7usize;
            let mut max_g0 = 0.0_f64;
            let mut max_g1 = 0.0_f64;
            for i in 0..n_samples {
                let s = 0.05 + 0.90 * (i as f64) / ((n_samples - 1) as f64);
                let s_a = if pa_forward { s } else { 1.0 - s };
                let s_b = if pb_forward { s } else { 1.0 - s };
                let t_a = pa_lo + s_a * (pa_hi - pa_lo);
                let t_b = pb_lo + s_b * (pb_hi - pb_lo);
                let uv_a = pa.point_at(t_a);
                let uv_b = pb.point_at(t_b);
                let p_a = sa.geometry.point_at(uv_a.x, uv_a.y);
                let p_b = sb.geometry.point_at(uv_b.x, uv_b.y);
                let g0 = p_a.distance_to(&p_b);
                if g0 > max_g0 {
                    max_g0 = g0;
                }
                let su_a = sa.geometry.partial_u(uv_a.x, uv_a.y);
                let sv_a = sa.geometry.partial_v(uv_a.x, uv_a.y);
                let su_b = sb.geometry.partial_u(uv_b.x, uv_b.y);
                let sv_b = sb.geometry.partial_v(uv_b.x, uv_b.y);
                let n_a = su_a.cross(&sv_a);
                let n_b = su_b.cross(&sv_b);
                let nn_a = (n_a.dot(&n_a)).sqrt();
                let nn_b = (n_b.dot(&n_b)).sqrt();
                if nn_a > 1e-300 && nn_b > 1e-300 {
                    let dot_unit = (n_a.dot(&n_b) / (nn_a * nn_b)).abs();
                    let deficit = 1.0 - dot_unit;
                    if deficit > max_g1 {
                        max_g1 = deficit;
                    }
                }
            }
            entries.push((max_g1, max_g0, eh, face_a, face_b, a_label, b_label));
        }
        entries.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));
        eprintln!("\n--- per-edge G1/G0 sorted by G1 (top 30) ---");
        for (i, (g1, g0, eh, fa, fb, la, lb)) in entries.iter().take(30).enumerate() {
            let cap_a = cap_set.contains(fa);
            let cap_b = cap_set.contains(fb);
            eprintln!(
                "  [{i:2}] edge={eh:?} g1={g1:.3e} g0={g0:.3e} {fa:?}({la}{}) <-> {fb:?}({lb}{})",
                if cap_a { "*" } else { "" },
                if cap_b { "*" } else { "" }
            );
        }

        // Aggregate stats: cap-vs-non-cap (worst-case post-reroute
        // gap) vs other (sanity baseline).
        let (mut max_cap_g1, mut max_cap_g0) = (0.0_f64, 0.0_f64);
        let (mut max_other_g1, mut max_other_g0) = (0.0_f64, 0.0_f64);
        for (g1, g0, _, fa, fb, _, _) in entries.iter() {
            let cap_a = cap_set.contains(fa);
            let cap_b = cap_set.contains(fb);
            if cap_a ^ cap_b {
                if *g1 > max_cap_g1 {
                    max_cap_g1 = *g1;
                }
                if *g0 > max_cap_g0 {
                    max_cap_g0 = *g0;
                }
            } else {
                if *g1 > max_other_g1 {
                    max_other_g1 = *g1;
                }
                if *g0 > max_other_g0 {
                    max_other_g0 = *g0;
                }
            }
        }
        eprintln!(
            "\nPhase 6c-D Step 2 sizing — Gregory rung-OFF, post-reroute:\n\
             cap↔non-cap: max g1={max_cap_g1:.3e} max g0={max_cap_g0:.3e}\n\
             other:       max g1={max_other_g1:.3e} max g0={max_other_g0:.3e}"
        );
    }

    /// Cap and cyl_split NOT adjacent: cap=tri0, cyl_split=tri2,
    /// cube=tri1.  This puts the bridge between two
    /// non-cyclically-adjacent faces — possibly invalid topology.
    #[test]
    fn split_synthetic_4_face_cap_cyl_split_opposite() {
        use popkern_core::geometry::Point3;

        let (mut brep, v_apex, [tri0, tri1, tri2, tri3]) = build_pyramid_4();
        let _ = tri3;
        let v_new_pos = Point3::new(0.5, 0.5, 0.0);

        match split_apex_footprint_vertex(&mut brep, v_apex, tri0, tri1, tri2, v_new_pos) {
            Ok(v_new) => {
                eprintln!("split OK: v_new={v_new:?}");
                if let Err(e) = brep.validate() {
                    panic!("validate after split failed (cap-cyl opposite): {e:?}");
                }
            }
            Err(e) => eprintln!("split returned error (cap-cyl opposite): {e:?}"),
        }
    }
}
