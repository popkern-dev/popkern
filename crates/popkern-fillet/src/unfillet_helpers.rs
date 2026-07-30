//! Internal helpers for [`crate::unfillet`].
//!
//! Holds the geometry routines (plane intersection, line-plane
//! intersection, support-face detection, blend-loop walking) that
//! `compute_unfillet_plan` depends on. Extracted from `unfillet.rs`
//! to keep each source file under the 500-line structural cap.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle};

use crate::error::BlendError;

#[derive(Debug, Clone, Copy)]
pub(crate) struct PlaneData {
    pub origin: Point3,
    pub normal: Vec3,
}

/// Find the two support faces of a Torus- or Cylinder-surface
/// blend by "plane neighbours whose normals are perpendicular to
/// the fillet axis". Returns `None` when the blend isn't a Torus
/// or Cylinder, or fewer than 2 such neighbours exist.
///
/// Straight-edge plane-plane fillets now produce Cylinder surfaces
/// instead of Torus, so both surface types are accepted.
///
/// More robust than reusing identify_blend's tangency sampling:
/// the dihedral angle at a boundary midpoint can be non-zero on a
/// torus even when the boundary is geometrically G1 (because the
/// midpoint may not lie exactly on the tangent contact curve).
/// The "perpendicular to axis" criterion is purely
/// surface-type-driven and doesn't depend on sample-point
/// precision.
pub(crate) fn find_torus_support_faces(
    brep: &BRep,
    blend_face: FaceHandle,
) -> Option<(FaceHandle, FaceHandle)> {
    let face_data = brep.faces.get(blend_face)?;
    let surface_on_face = face_data.surface.as_ref()?;
    let torus_axis = match &surface_on_face.geometry {
        SurfaceGeometry::Torus { axis, .. } | SurfaceGeometry::Cylinder { axis, .. } => *axis,
        _ => return None,
    };

    let loop_data = brep.loops.get(face_data.outer_loop)?;
    let start = loop_data.half_edge;

    let mut candidates: Vec<FaceHandle> = Vec::new();
    let mut he = start;
    let mut safety = 0;
    loop {
        safety += 1;
        if safety > 10_000 {
            break;
        }
        let he_data = brep.half_edges.get(he)?;
        let twin_data = brep.half_edges.get(he_data.twin)?;
        let neighbour_loop = brep.loops.get(twin_data.loop_)?;
        let neighbour = neighbour_loop.face;

        if neighbour != blend_face
            && !candidates.contains(&neighbour)
            && let Some(neighbour_face) = brep.faces.get(neighbour)
            && let Some(surface_on_neighbour) = neighbour_face.surface.as_ref()
            && let SurfaceGeometry::Plane { normal, .. } = &surface_on_neighbour.geometry
        {
            let oriented_normal = if surface_on_neighbour.orientation {
                *normal
            } else {
                -*normal
            };
            let dot = oriented_normal.dot(&torus_axis).abs();
            // 0.05 ≈ 3°
            if dot < 0.05 {
                candidates.push(neighbour);
            }
        }

        he = he_data.next;
        if he == start {
            break;
        }
    }

    if candidates.len() >= 2 {
        Some((candidates[0], candidates[1]))
    } else {
        None
    }
}

/// Pull the plane geometry out of a face. Returns
/// `Err(BlendError::InvalidParameter)` when the face surface is
/// not a [`SurfaceGeometry::Plane`].
pub(crate) fn extract_plane(brep: &BRep, face: FaceHandle) -> Result<PlaneData, BlendError> {
    let face_data = brep
        .faces
        .get(face)
        .ok_or_else(|| BlendError::InvalidParameter {
            parameter: "support_face".to_string(),
            value: format!("{face:?}"),
            reason: "face handle does not exist in the BRep".to_string(),
        })?;
    let surface = face_data
        .surface
        .as_ref()
        .ok_or_else(|| BlendError::InvalidParameter {
            parameter: "support_face".to_string(),
            value: format!("{face:?}"),
            reason: "face has no attached surface geometry".to_string(),
        })?;
    match &surface.geometry {
        SurfaceGeometry::Plane { origin, normal, .. } => {
            let oriented_normal = if surface.orientation {
                *normal
            } else {
                -*normal
            };
            Ok(PlaneData {
                origin: *origin,
                normal: oriented_normal,
            })
        }
        other => Err(BlendError::InvalidParameter {
            parameter: "support_face.surface".to_string(),
            value: match other {
                SurfaceGeometry::Cylinder { .. } => "Cylinder".to_string(),
                SurfaceGeometry::Sphere { .. } => "Sphere".to_string(),
                SurfaceGeometry::Torus { .. } => "Torus".to_string(),
                SurfaceGeometry::Cone { .. } => "Cone".to_string(),
                SurfaceGeometry::Nurbs(_) => "Nurbs".to_string(),
                _ => "non-Plane".to_string(),
            },
            reason: "open-blend unfillet plan supports planar supports only; closed \
                     plane-cylinder rims plan through the annulus route, and open \
                     curved-support strips restore via the operation record \
                     (unfillet_via_record / modify_fillet_radius_via_record)"
                .to_string(),
        }),
    }
}

/// Compute the intersection LINE of two planes. Returns
/// `(line_origin, line_direction)` where `line_direction` is unit
/// length. Returns `Err` when the planes are parallel.
pub(crate) fn intersect_planes(a: &PlaneData, b: &PlaneData) -> Result<(Point3, Vec3), BlendError> {
    let n1 = a.normal;
    let n2 = b.normal;

    let direction = n1.cross(&n2);
    let dir_len = direction.length();
    if dir_len < 1e-9 {
        return Err(BlendError::InvalidParameter {
            parameter: "support_planes".to_string(),
            value: "parallel".to_string(),
            reason: "support planes are parallel — no unique intersection line".to_string(),
        });
    }
    let direction = direction * (1.0 / dir_len);

    // Find a point on the intersection line. Solve a 2×2 system
    // in the plane spanned by n1 and n2:
    //   P = α·n1 + β·n2
    //   n1·P = α·|n1|² + β·(n1·n2) = d1
    //   n2·P = α·(n1·n2) + β·|n2|² = d2
    let d1 = n1.dot(&Vec3::new(a.origin.x, a.origin.y, a.origin.z));
    let d2 = n2.dot(&Vec3::new(b.origin.x, b.origin.y, b.origin.z));
    let n1_dot = n1.dot(&n1);
    let n2_dot = n2.dot(&n2);
    let n12 = n1.dot(&n2);
    let det = n1_dot * n2_dot - n12 * n12;
    if det.abs() < 1e-18 {
        return Err(BlendError::InvalidParameter {
            parameter: "support_planes".to_string(),
            value: "degenerate".to_string(),
            reason: "support plane normals are linearly dependent".to_string(),
        });
    }
    let alpha = (d1 * n2_dot - d2 * n12) / det;
    let beta = (d2 * n1_dot - d1 * n12) / det;
    let line_origin = Point3::new(
        alpha * n1.x + beta * n2.x,
        alpha * n1.y + beta * n2.y,
        alpha * n1.z + beta * n2.z,
    );

    Ok((line_origin, direction))
}

/// Walk the blend face's outer loop and collect the planes of every
/// non-support neighbour face. These are the "end face" planes that
/// bound the fillet's running direction.
pub(crate) fn collect_end_face_planes(
    brep: &BRep,
    blend_face: FaceHandle,
    support_a: FaceHandle,
    support_b: FaceHandle,
) -> Result<Vec<PlaneData>, BlendError> {
    let face_data = brep
        .faces
        .get(blend_face)
        .ok_or_else(|| BlendError::InvalidParameter {
            parameter: "blend_face".to_string(),
            value: format!("{blend_face:?}"),
            reason: "face handle does not exist".to_string(),
        })?;
    let loop_data =
        brep.loops
            .get(face_data.outer_loop)
            .ok_or_else(|| BlendError::InvalidParameter {
                parameter: "blend_face.outer_loop".to_string(),
                value: format!("{:?}", face_data.outer_loop),
                reason: "loop handle does not exist".to_string(),
            })?;

    let mut seen_neighbours = Vec::new();
    let mut planes = Vec::new();
    let start = loop_data.half_edge;
    let mut he = start;
    let mut safety = 0;
    loop {
        safety += 1;
        if safety > 10_000 {
            break;
        }
        let Some(he_data) = brep.half_edges.get(he) else {
            break;
        };
        let twin = brep.half_edges.get(he_data.twin);
        if let Some(twin_data) = twin
            && let Some(neighbour_loop) = brep.loops.get(twin_data.loop_)
        {
            let neighbour = neighbour_loop.face;
            if neighbour != support_a
                && neighbour != support_b
                && neighbour != blend_face
                && !seen_neighbours.contains(&neighbour)
                && let Ok(plane) = extract_plane(brep, neighbour)
            {
                seen_neighbours.push(neighbour);
                planes.push(plane);
            }
        }
        he = he_data.next;
        if he == start {
            break;
        }
    }
    Ok(planes)
}

/// Intersect a parametric line with a plane. Returns `None` for
/// parallel cases.
pub(crate) fn intersect_line_plane(
    line_origin: &Point3,
    line_direction: &Vec3,
    plane: &PlaneData,
) -> Option<Point3> {
    let denom = line_direction.dot(&plane.normal);
    if denom.abs() < 1e-12 {
        return None;
    }
    let diff = Vec3::new(
        plane.origin.x - line_origin.x,
        plane.origin.y - line_origin.y,
        plane.origin.z - line_origin.z,
    );
    let t = diff.dot(&plane.normal) / denom;
    Some(Point3::new(
        line_origin.x + t * line_direction.x,
        line_origin.y + t * line_direction.y,
        line_origin.z + t * line_direction.z,
    ))
}

/// Pick the two endpoints from `endpoints` that are furthest apart
/// when projected onto the line direction. They define the extent
/// of the restored sharp edge.
pub(crate) fn pick_extreme_endpoints(endpoints: &[Point3], direction: &Vec3) -> (Point3, Point3) {
    debug_assert!(endpoints.len() >= 2);

    // Project each endpoint onto the line direction (taking
    // endpoints[0] as the parameter origin).
    let origin = endpoints[0];
    let mut min_t = 0.0_f64;
    let mut max_t = 0.0_f64;
    let mut min_p = origin;
    let mut max_p = origin;
    for p in endpoints {
        let d = Vec3::new(p.x - origin.x, p.y - origin.y, p.z - origin.z);
        let t = d.dot(direction);
        if t < min_t {
            min_t = t;
            min_p = *p;
        }
        if t > max_t {
            max_t = t;
            max_p = *p;
        }
    }
    (min_p, max_p)
}

/// Walk the blend's outer loop, identify long G1-tangent boundary
/// half-edges (one to each support), and collect every edge that
/// the surgery would need to delete.
pub(crate) fn walk_blend_boundary(
    brep: &BRep,
    blend_face: FaceHandle,
    support_a: FaceHandle,
    support_b: FaceHandle,
) -> (
    Option<HalfEdgeHandle>,
    Option<HalfEdgeHandle>,
    Vec<EdgeHandle>,
) {
    let mut long_he_to_a = None;
    let mut long_he_to_b = None;
    let mut edges_to_remove = Vec::new();

    let Some(face_data) = brep.faces.get(blend_face) else {
        return (None, None, edges_to_remove);
    };
    let Some(loop_data) = brep.loops.get(face_data.outer_loop) else {
        return (None, None, edges_to_remove);
    };
    let start = loop_data.half_edge;
    let mut he = start;
    let mut safety = 0;
    loop {
        safety += 1;
        if safety > 10_000 {
            break;
        }
        let Some(he_data) = brep.half_edges.get(he) else {
            break;
        };
        let twin = brep.half_edges.get(he_data.twin);
        if let Some(twin_data) = twin
            && let Some(neighbour_loop) = brep.loops.get(twin_data.loop_)
        {
            let neighbour = neighbour_loop.face;
            if neighbour == support_a && long_he_to_a.is_none() {
                long_he_to_a = Some(he);
            } else if neighbour == support_b && long_he_to_b.is_none() {
                long_he_to_b = Some(he);
            }
            if !edges_to_remove.contains(&he_data.edge) {
                edges_to_remove.push(he_data.edge);
            }
        }
        he = he_data.next;
        if he == start {
            break;
        }
    }

    (long_he_to_a, long_he_to_b, edges_to_remove)
}

// Surgery helpers live in `unfillet_surgery` — kept separate so
// each source file stays under the 500-line structural cap.

/// True when `face` is the annulus produced by the Round 5c-3
/// closed-fillet construction path.
///
/// Detection is topological, not geometric: a closed-fillet annulus
/// has exactly one outer loop AND exactly one inner loop, and BOTH
/// loops are single-element (one self-looping half-edge). Non-
/// closed-fillet faces with holes (e.g. a cube top with a cylinder
/// cutout) have a multi-element outer loop and so are excluded.
/// Checking topology rather than surface type is robust to the
/// surface-fix post-pass, which can reassign a closed-fillet face's
/// NURBS sweep to a nearby cylinder when the vertices happen to lie
/// on it.
pub(crate) fn is_closed_fillet_annulus(brep: &BRep, face: FaceHandle) -> bool {
    let Some(face_data) = brep.faces.get(face) else {
        return false;
    };
    if face_data.inner_loops.len() != 1 {
        return false;
    }
    let outer_single = loop_is_single_element(brep, face_data.outer_loop);
    let inner_single = loop_is_single_element(brep, face_data.inner_loops[0]);
    outer_single && inner_single
}

/// True when `loop_` is a closed loop made of a single self-looping
/// half-edge (next == prev == self).
fn loop_is_single_element(brep: &BRep, loop_: popkern_core::topology::LoopHandle) -> bool {
    let Some(l) = brep.loops.get(loop_) else {
        return false;
    };
    let start = l.half_edge;
    let Some(he) = brep.half_edges.get(start) else {
        return false;
    };
    he.next == start && he.prev == start
}
