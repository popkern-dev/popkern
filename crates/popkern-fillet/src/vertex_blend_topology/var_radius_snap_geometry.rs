use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, VertexHandle};

/// Re-derive the per-fillet sphere for `cap_face` by walking its
/// outer loop, collecting the cylinder twin faces, and finding the
/// per-fillet cap whose centre lies on one of those cylinder axes.
pub(super) fn re_derive_cap_sphere(
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

/// Project `p` onto `sphere(c, r) / cylinder(origin, axis, cyl_r)`.
/// Returns `None` when the rolling-ball on-axis precondition does
/// not hold.
pub(super) fn project_onto_sphere_cylinder(
    p: Point3,
    sphere_centre: Point3,
    sphere_radius: f64,
    cyl_origin: Point3,
    cyl_axis: Vec3,
    cyl_radius: f64,
) -> Option<Point3> {
    let axis_unit = cyl_axis.normalized()?;
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
        return None;
    }
    let scale = cyl_radius / radial_len;
    Some(Point3::new(
        sphere_centre.x + radial.x * scale,
        sphere_centre.y + radial.y * scale,
        sphere_centre.z + radial.z * scale,
    ))
}

/// Project `p` onto `sphere(c, r) / plane(origin, normal)`.
/// Returns `None` if the intersection circle is empty or
/// degenerate.
pub(super) fn project_onto_sphere_plane(
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
    let circle_centre = Point3::new(
        sphere_centre.x - n.x * d_signed,
        sphere_centre.y - n.y * d_signed,
        sphere_centre.z - n.z * d_signed,
    );
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
pub(super) fn collect_loop_half_edges(brep: &BRep, start: HalfEdgeHandle) -> Vec<HalfEdgeHandle> {
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

/// Snap a vertex to `new_pos`, syncing every incident edge endpoint
/// so loop closure stays intact.
pub(super) fn update_incident_edges_for_vertex_move(
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
