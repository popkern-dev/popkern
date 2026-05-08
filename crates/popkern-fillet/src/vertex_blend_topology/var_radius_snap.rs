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

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, VertexHandle};

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

/// Re-derive the per-fillet sphere for `cap_face` by walking its
/// outer loop, collecting the (origin, axis, radius) of every
/// cylinder twin face, and finding the `per_fillet_caps[i]` whose
/// centre lies on at least one of those axes (within tolerance).
///
/// Robust against sub_face ↔ per_fillet_caps index misalignment
/// (the upstream `pick_per_fillet_cap_for_face` can silently fall
/// back to its index-based default when its primary-cylinder
/// matching fails). Returns `None` if no cylinder edge in the loop
/// has any per_fillet_caps centre on its axis — typical when the
/// cap's outer loop has been re-routed away from cylinder edges.
fn re_derive_cap_sphere(
    brep: &BRep,
    cap_face: FaceHandle,
    per_fillet_caps: &[(Point3, f64)],
) -> Option<(Point3, f64)> {
    let face_data = brep.faces.get(cap_face)?;
    let loop_data = brep.loops.get(face_data.outer_loop)?;
    let start = loop_data.half_edge;
    let hes = collect_loop_half_edges(brep, start);

    for &he_h in &hes {
        let Some(he) = brep.half_edges.get(he_h) else {
            continue;
        };
        let twin_face = brep
            .half_edges
            .get(he.twin)
            .and_then(|t| brep.loops.get(t.loop_))
            .map(|l| l.face);
        let Some(twin_face) = twin_face else { continue };
        let Some(SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        }) = brep
            .faces
            .get(twin_face)
            .and_then(|f| f.surface.as_ref())
            .map(|s| s.geometry.clone())
        else {
            continue;
        };
        let Some(axis_unit) = axis.normalized() else {
            continue;
        };
        for &(centre, sphere_radius) in per_fillet_caps {
            if (sphere_radius - radius).abs() > radius * 0.05 {
                continue;
            }
            let v = Vec3::new(
                centre.x - origin.x,
                centre.y - origin.y,
                centre.z - origin.z,
            );
            let along = v.dot(&axis_unit);
            let proj = Point3::new(
                origin.x + axis_unit.x * along,
                origin.y + axis_unit.y * along,
                origin.z + axis_unit.z * along,
            );
            if centre.distance_to(&proj) < radius * 0.05 {
                return Some((centre, sphere_radius));
            }
        }
    }
    None
}

/// Project `p` onto `sphere(c, r) ∩ cylinder(origin, axis, cyl_r)`.
/// Returns `None` if the sphere doesn't intersect the cylinder
/// (sphere centre off-axis by more than `r`, or the intersection
/// is empty after restricting to the cylinder's "side" of `p`).
///
/// Geometry: when the sphere centre lies *on* the cylinder axis,
/// the intersection is a circle in the cylinder's perpendicular
/// cross-section through the sphere centre. The circle's centre
/// is the sphere centre, radius = sphere radius (sphere and
/// cylinder share the great-circle equator). Project `p` onto the
/// cross-section plane, then onto the circle.
///
/// When sphere centre is off-axis, the intersection is still a
/// circle but in a plane tilted relative to both surfaces. This
/// path returns `None` for that case (the per-fillet rolling-ball
/// construction guarantees on-axis centres for the cap-tangent
/// arc; off-axis centres indicate Phase 6c-E hasn't fired).
fn project_onto_sphere_cylinder(
    p: Point3,
    sphere_centre: Point3,
    sphere_radius: f64,
    cyl_origin: Point3,
    cyl_axis: Vec3,
    cyl_radius: f64,
) -> Option<Point3> {
    let axis_unit = cyl_axis.normalized()?;
    // Sphere centre projected onto cylinder axis.
    let to_centre = Vec3::new(
        sphere_centre.x - cyl_origin.x,
        sphere_centre.y - cyl_origin.y,
        sphere_centre.z - cyl_origin.z,
    );
    let along = to_centre.dot(&axis_unit);
    let sphere_centre_axial = Point3::new(
        cyl_origin.x + axis_unit.x * along,
        cyl_origin.y + axis_unit.y * along,
        cyl_origin.z + axis_unit.z * along,
    );
    // Off-axis distance — sphere centre must lie on the cylinder
    // axis for sphere ∩ cylinder to be a perpendicular great
    // circle (the rolling-ball cap-tangent geometry). Off-axis
    // configurations (neighbour cylinder, etc.) yield an
    // off-axial intersection that's not the cap-tangent locus —
    // bail.
    let perp = Vec3::new(
        sphere_centre.x - sphere_centre_axial.x,
        sphere_centre.y - sphere_centre_axial.y,
        sphere_centre.z - sphere_centre_axial.z,
    );
    let off_axis = perp.length();
    if off_axis > sphere_radius * 1e-3 {
        return None;
    }
    if (cyl_radius - sphere_radius).abs() > sphere_radius * 1e-2 {
        return None;
    }

    // Project p onto the perpendicular plane through sphere_centre.
    let to_p = Vec3::new(
        p.x - sphere_centre.x,
        p.y - sphere_centre.y,
        p.z - sphere_centre.z,
    );
    let along_axis = to_p.dot(&axis_unit);
    let p_in_plane = Point3::new(
        p.x - axis_unit.x * along_axis,
        p.y - axis_unit.y * along_axis,
        p.z - axis_unit.z * along_axis,
    );
    let radial = Vec3::new(
        p_in_plane.x - sphere_centre.x,
        p_in_plane.y - sphere_centre.y,
        p_in_plane.z - sphere_centre.z,
    );
    let radial_len = radial.length();
    if radial_len < 1e-12 {
        // p is on the cylinder axis; ambiguous direction. Bail.
        return None;
    }
    let scale = cyl_radius / radial_len;
    Some(Point3::new(
        sphere_centre.x + radial.x * scale,
        sphere_centre.y + radial.y * scale,
        sphere_centre.z + radial.z * scale,
    ))
}

/// Project `p` onto `sphere(c, r) ∩ plane(origin, normal)`.
/// Returns `None` if the plane doesn't intersect the sphere
/// (perpendicular distance from sphere centre > radius), or the
/// intersection is degenerate.
///
/// Geometry: sphere ∩ plane = circle in the plane, centred at
/// the sphere centre's projection onto the plane, with radius
/// `sqrt(r² - d²)` where `d` is the perpendicular distance from
/// sphere centre to the plane. Project `p` onto the plane, then
/// onto the circle.
fn project_onto_sphere_plane(
    p: Point3,
    sphere_centre: Point3,
    sphere_radius: f64,
    plane_origin: Point3,
    plane_normal: Vec3,
) -> Option<Point3> {
    let n = plane_normal.normalized()?;
    let to_centre = Vec3::new(
        sphere_centre.x - plane_origin.x,
        sphere_centre.y - plane_origin.y,
        sphere_centre.z - plane_origin.z,
    );
    let d_signed = to_centre.dot(&n);
    let d = d_signed.abs();
    if d > sphere_radius {
        return None;
    }
    let r_circle_sq = sphere_radius * sphere_radius - d * d;
    if r_circle_sq < 1e-18 {
        return None;
    }
    let r_circle = r_circle_sq.sqrt();
    // Circle centre = sphere centre projected onto plane.
    let circle_centre = Point3::new(
        sphere_centre.x - n.x * d_signed,
        sphere_centre.y - n.y * d_signed,
        sphere_centre.z - n.z * d_signed,
    );
    // Project p onto plane.
    let to_p = Vec3::new(
        p.x - plane_origin.x,
        p.y - plane_origin.y,
        p.z - plane_origin.z,
    );
    let p_signed = to_p.dot(&n);
    let p_in_plane = Point3::new(
        p.x - n.x * p_signed,
        p.y - n.y * p_signed,
        p.z - n.z * p_signed,
    );
    // Project p_in_plane onto circle (centre, radius r_circle in plane).
    let radial = Vec3::new(
        p_in_plane.x - circle_centre.x,
        p_in_plane.y - circle_centre.y,
        p_in_plane.z - circle_centre.z,
    );
    let radial_len = radial.length();
    if radial_len < 1e-12 {
        return None;
    }
    let scale = r_circle / radial_len;
    Some(Point3::new(
        circle_centre.x + radial.x * scale,
        circle_centre.y + radial.y * scale,
        circle_centre.z + radial.z * scale,
    ))
}

/// Walk the loop starting at `start` and collect all half-edges,
/// capped at 64 entries to guard against malformed loops.
fn collect_loop_half_edges(brep: &BRep, start: HalfEdgeHandle) -> Vec<HalfEdgeHandle> {
    let mut hes = Vec::new();
    let mut cursor = start;
    let mut safety = 64usize;
    loop {
        if safety == 0 {
            break;
        }
        safety -= 1;
        hes.push(cursor);
        let Some(he) = brep.half_edges.get(cursor) else {
            break;
        };
        cursor = he.next;
        if cursor == start {
            break;
        }
    }
    hes
}

/// Phase 6c-D Stage 2: snap a vertex to `new_pos`, syncing every
/// incident edge's CP[0] (primary side) or CP[last] (twin side) so
/// loop closure stays intact. Returns the count of edge curves
/// updated. Equivalent to the same-named helper in `var_radius.rs`
/// (Path B); duplicated here to avoid pub-promoting that one.
fn update_incident_edges_for_vertex_move(
    brep: &mut BRep,
    v: VertexHandle,
    new_pos: Point3,
) -> usize {
    use popkern_core::geometry::CurveGeometry;

    let halves: Vec<HalfEdgeHandle> = brep
        .half_edges
        .iter()
        .filter_map(|(h, he)| (he.origin == v).then_some(h))
        .collect();

    let mut count = 0usize;
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
            count += 1;
        }
    }

    if let Some(vd) = brep.vertices.get_mut(v) {
        vd.position = new_pos;
    }
    count
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
