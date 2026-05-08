use super::var_radius_classify::{TwinKind, classify_cap_loops};
use super::var_radius_geometry::walk_loop_closes;
use super::*;

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
