//! Phase 6c-D Stage 2 — per-vertex snap-only surgery for pyramid
//! and other apex topologies where Path B reroute can't fire.
//!
//! Path B's existing `plan_path_b` requires each per-fillet sphere
//! centre to project onto a *neighbour* cylinder's axis (the cube
//! body-diagonal coincidence). Pyramid violates this — the 4 slant
//! cylinder axes converge above the apex, and no per-fillet sphere
//! centre lies on a neighbour's axis. `plan_path_b` returns `None`
//! for all NeighbourCylinder targets on pyramid (verified empirically).
//!
//! This module's `apply_vertex_snap_for_caps` is a generalised
//! drop-in: it walks each cap face's outer loop, snaps each vertex
//! to its closest-point projection onto the appropriate
//! `sphere_k ∩ cylinder_m` (or `sphere_k ∩ plane`) circle, and
//! uses the same `update_incident_edges` machinery Path B uses
//! to keep loop closure intact.
//!
//! **Geometric basis:**
//! - For a cap-loop vertex `v` incident to an edge whose twin face
//!   is a Cylinder: project `v` onto `sphere_k ∩ cylinder_m`. This
//!   is a circle in 3-space (the perpendicular cross-section of the
//!   cylinder at the projection of the sphere centre onto the
//!   cylinder axis), provided the sphere centre lies *on* the
//!   cylinder axis (Phase 6c-E rolling-ball construction).
//! - For a vertex incident only to plane edges (cap-cap inner spoke
//!   or absorbed Plane fragment vertex): project onto
//!   `sphere_k ∩ plane`. This is a circle in 3-space at the plane's
//!   intersection with the sphere — provided the plane's distance
//!   from the sphere centre is < radius.
//!
//! Both projections are closed-form. After snap, the cap-tangent
//! arc on the cylinder side is (sphere ∩ cylinder), the spoke
//! arc on the plane side is (sphere ∩ plane), and the cap surface
//! (Sphere or Gregory built to interpolate the substituted boundary)
//! interpolates these arcs by construction.
//!
//! Activation: gated on the same env var as Phase 6c-A
//! (`POPKERN_VAR_RADIUS_REROUTE=1`) — the existing
//! `reroute_neighbour_cylinder_edges` call site dispatches to this
//! when cap loops have no NeighbourCylinder edges (pyramid case).

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

use super::var_radius_snap_geometry::{
    collect_loop_half_edges, project_onto_sphere_cylinder, project_onto_sphere_plane,
    re_derive_cap_sphere, update_incident_edges_for_vertex_move,
};

/// Result of `apply_vertex_snap_for_caps` — count of vertices
/// whose position was actually moved (delta > 1e-9).
#[derive(Debug, Default)]
pub(crate) struct SnapResult {
    pub vertices_moved: usize,
    pub edges_with_curves_updated: usize,
}

/// Phase 6c-D Stage 2 entry point. For each cap face, walk its
/// outer loop, project each vertex onto its incident analytic
/// circle (sphere ∩ cylinder or sphere ∩ plane), and apply the
/// snap. Cap-cap inner-spoke vertices (where neither incident
/// twin face is a Cylinder *and* both incident twins lead to
/// Sphere/Gregory caps) are projected onto the *closer* sphere's
/// natural intersection with the spoke direction — handled in a
/// follow-up; this stage's first pass leaves them alone.
///
/// `cap_faces[i]` and `per_fillet_caps[i]` are parallel arrays —
/// `per_fillet_caps[i] = (sphere_centre_i, sphere_radius_i)` for
/// cap face `i`.
pub(crate) fn apply_vertex_snap_for_caps(
    brep: &mut BRep,
    cap_faces: &[FaceHandle],
    per_fillet_caps: &[(Point3, f64)],
) -> SnapResult {
    if cap_faces.len() != per_fillet_caps.len() {
        return SnapResult::default();
    }
    let snap_tol = 1e-9_f64;
    let mut moved = 0usize;
    let mut edges_updated = 0usize;

    // Phase A — collect per-vertex snap targets across all caps.
    // A vertex may be in multiple cap loops; we record one target
    // per cap and average if the cap-side ambiguity matters. In
    // practice the cap-tangent arc projection is the same for
    // both incident caps (they share the cylinder), so averaging
    // reduces to a single point.
    use std::collections::BTreeMap;
    let mut targets: BTreeMap<VertexHandle, (Point3, usize)> = BTreeMap::new();

    for &cap_face in cap_faces.iter() {
        // Re-derive the right per-fillet sphere for THIS cap face by
        // finding which per_fillet_caps[i] has a centre on one of
        // the cap loop's cylinder-twin axes. This is robust against
        // sub_face/per_fillet_caps index misalignment introduced
        // by upstream pick logic that silently falls back to the
        // index-based default when matching fails.
        let Some((sphere_centre, sphere_radius)) =
            re_derive_cap_sphere(brep, cap_face, per_fillet_caps)
        else {
            continue;
        };
        let Some(face_data) = brep.faces.get(cap_face) else {
            continue;
        };
        let Some(loop_data) = brep.loops.get(face_data.outer_loop) else {
            continue;
        };
        let start = loop_data.half_edge;
        let hes = collect_loop_half_edges(brep, start);

        for &he_handle in &hes {
            let Some(he) = brep.half_edges.get(he_handle) else {
                continue;
            };
            let v_origin = he.origin;
            let next_he = he.next;
            let v_dest = brep.half_edges.get(next_he).map(|n| n.origin);
            let twin_face = brep
                .half_edges
                .get(he.twin)
                .and_then(|t| brep.loops.get(t.loop_))
                .map(|l| l.face);
            let twin_surf = twin_face
                .and_then(|f| brep.faces.get(f).and_then(|fd| fd.surface.as_ref()))
                .map(|s| s.geometry.clone());

            // Project BOTH endpoints (origin + destination) of each
            // cap-tangent arc edge so the entire arc lies on
            // sphere ∩ cylinder. The destination vertex would also
            // be processed when iteration reaches it as an origin,
            // but its OWN-side incident edge may have a different
            // (non-Cylinder) twin (e.g., cap-cap inner spoke), so
            // explicitly include it here when we know the cylinder.
            let p_origin = brep.vertices[v_origin].position;
            let p_dest = v_dest.map(|v| brep.vertices[v].position);
            let endpoints: Vec<(VertexHandle, Point3)> = match v_dest {
                Some(vd) if matches!(twin_surf, Some(SurfaceGeometry::Cylinder { .. })) => {
                    if let Some(p_dest) = p_dest {
                        vec![(v_origin, p_origin), (vd, p_dest)]
                    } else {
                        vec![(v_origin, p_origin)]
                    }
                }
                _ => vec![(v_origin, p_origin)],
            };

            for (v, p) in endpoints {
                let target = match &twin_surf {
                    Some(SurfaceGeometry::Cylinder {
                        origin: cyl_origin,
                        axis: cyl_axis,
                        radius: cyl_radius,
                    }) => project_onto_sphere_cylinder(
                        p,
                        sphere_centre,
                        sphere_radius,
                        *cyl_origin,
                        *cyl_axis,
                        *cyl_radius,
                    ),
                    Some(SurfaceGeometry::Plane {
                        origin: plane_origin,
                        normal: plane_normal,
                        ..
                    }) => project_onto_sphere_plane(
                        p,
                        sphere_centre,
                        sphere_radius,
                        *plane_origin,
                        *plane_normal,
                    ),
                    _ => None,
                };
                let Some(target_pos) = target else {
                    continue;
                };
                let move_dist = p.distance_to(&target_pos);
                if move_dist > 2.0 * sphere_radius {
                    continue;
                }
                let entry = targets.entry(v).or_insert((Point3::new(0.0, 0.0, 0.0), 0));
                entry.0 = Point3::new(
                    entry.0.x + target_pos.x,
                    entry.0.y + target_pos.y,
                    entry.0.z + target_pos.z,
                );
                entry.1 += 1;
            }
        }
    }

    // Phase B — apply snaps once per vertex (averaging if multiple
    // caps requested moves). update_incident_edges_for_vertex_move
    // syncs every incident edge's CP[0]/CP[last] to keep loop
    // closure intact; that primitive lives in `var_radius.rs`.
    for (v, (sum, n)) in targets {
        if n == 0 {
            continue;
        }
        let avg = Point3::new(sum.x / n as f64, sum.y / n as f64, sum.z / n as f64);
        let Some(vd) = brep.vertices.get(v) else {
            continue;
        };
        let delta = vd.position.distance_to(&avg);
        if delta < snap_tol {
            continue;
        }
        let edge_count = update_incident_edges_for_vertex_move(brep, v, avg);
        edges_updated += edge_count;
        moved += 1;
    }

    SnapResult {
        vertices_moved: moved,
        edges_with_curves_updated: edges_updated,
    }
}

/// Phase 6c-D Stage 2: detect whether a list of cap faces has
/// **no NeighbourCylinder edges** in their outer loops — the
/// pyramid-style topology where Path B can't fire and vertex-snap
/// is the right surgery. Each cap touches exactly one cylinder
/// face (its OwnCylinder); cube-style caps touch ≥ 2 cylinders.
///
/// Currently kept as a diagnostic helper; the production dispatch
/// in `reroute_neighbour_cylinder_edges` uses Path B's success
/// count as the fallback trigger (more robust to absorbed-Plane
/// edges that grow neighbour-cylinder counts post-sliver-absorb).
#[allow(dead_code)]
pub(crate) fn caps_have_no_neighbour_cylinders(brep: &BRep, cap_faces: &[FaceHandle]) -> bool {
    use std::collections::BTreeSet;
    for &cap_face in cap_faces {
        let mut cyls: BTreeSet<FaceHandle> = BTreeSet::new();
        let Some(face_data) = brep.faces.get(cap_face) else {
            continue;
        };
        let Some(loop_data) = brep.loops.get(face_data.outer_loop) else {
            continue;
        };
        let start = loop_data.half_edge;
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
            if let Some(twin_loop) = brep.half_edges.get(he.twin).map(|t| t.loop_)
                && let Some(twin_face) = brep.loops.get(twin_loop).map(|l| l.face)
                && twin_face != cap_face
                && let Some(twin_sof) = brep.faces.get(twin_face).and_then(|f| f.surface.as_ref())
                && matches!(twin_sof.geometry, SurfaceGeometry::Cylinder { .. })
            {
                cyls.insert(twin_face);
            }
            cursor = he.next;
            if cursor == start {
                break;
            }
        }
        if cyls.len() >= 2 {
            return false;
        }
    }
    true
}
