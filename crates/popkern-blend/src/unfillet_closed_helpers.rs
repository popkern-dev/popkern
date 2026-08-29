//! Helper structs and functions for
//! [`crate::unfillet_closed::unfillet_closed_annulus`].
//!
//! Split out of `unfillet_closed.rs` to keep each source file under
//! the 500-line structural limit. All items here are crate-private.

use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

use crate::error::BlendError;

/// Which side of the annulus owns the planar support.
pub(crate) enum Side {
    /// The outer-loop side is the planar support.
    Outer,
    /// The inner-loop side is the planar support.
    Inner,
    /// Couldn't determine which.
    Ambiguous,
}

pub(crate) fn classify_support_pair(brep: &BRep, outer: FaceHandle, inner: FaceHandle) -> Side {
    let outer_geom = brep
        .faces
        .get(outer)
        .and_then(|f| f.surface.as_ref())
        .map(|s| &s.geometry);
    let inner_geom = brep
        .faces
        .get(inner)
        .and_then(|f| f.surface.as_ref())
        .map(|s| &s.geometry);
    match (outer_geom, inner_geom) {
        (Some(SurfaceGeometry::Plane { .. }), Some(SurfaceGeometry::Cylinder { .. })) => {
            Side::Outer
        }
        (Some(SurfaceGeometry::Cylinder { .. }), Some(SurfaceGeometry::Plane { .. })) => {
            Side::Inner
        }
        _ => Side::Ambiguous,
    }
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct PlaneGeometry {
    pub(crate) origin: Point3,
    pub(crate) normal: Vec3,
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct CylinderGeometry {
    pub(crate) origin: Point3,
    pub(crate) axis: Vec3,
    pub(crate) radius: f64,
}

pub(crate) fn extract_plane_geometry(
    brep: &BRep,
    face: FaceHandle,
) -> Result<PlaneGeometry, BlendError> {
    let face_data = brep
        .faces
        .get(face)
        .ok_or_else(|| BlendError::InvalidParameter {
            parameter: "plane_support".to_string(),
            value: format!("{face:?}"),
            reason: "plane support face missing".to_string(),
        })?;
    let surface = face_data
        .surface
        .as_ref()
        .ok_or_else(|| BlendError::InvalidParameter {
            parameter: "plane_support.surface".to_string(),
            value: format!("{face:?}"),
            reason: "plane support has no attached surface".to_string(),
        })?;
    if let SurfaceGeometry::Plane { origin, normal, .. } = &surface.geometry {
        Ok(PlaneGeometry {
            origin: *origin,
            normal: *normal,
        })
    } else {
        Err(BlendError::InvalidParameter {
            parameter: "plane_support.surface".to_string(),
            value: "non-Plane".to_string(),
            reason: "expected Plane support".to_string(),
        })
    }
}

pub(crate) fn extract_cylinder_geometry(
    brep: &BRep,
    face: FaceHandle,
) -> Result<CylinderGeometry, BlendError> {
    let face_data = brep
        .faces
        .get(face)
        .ok_or_else(|| BlendError::InvalidParameter {
            parameter: "cyl_support".to_string(),
            value: format!("{face:?}"),
            reason: "cylinder support face missing".to_string(),
        })?;
    let surface = face_data
        .surface
        .as_ref()
        .ok_or_else(|| BlendError::InvalidParameter {
            parameter: "cyl_support.surface".to_string(),
            value: format!("{face:?}"),
            reason: "cylinder support has no attached surface".to_string(),
        })?;
    if let SurfaceGeometry::Cylinder {
        origin,
        axis,
        radius,
    } = &surface.geometry
    {
        Ok(CylinderGeometry {
            origin: *origin,
            axis: *axis,
            radius: *radius,
        })
    } else {
        Err(BlendError::InvalidParameter {
            parameter: "cyl_support.surface".to_string(),
            value: "non-Cylinder".to_string(),
            reason: "expected Cylinder support".to_string(),
        })
    }
}

/// Compute the centre and radius of the circle where a plane
/// intersects a cylinder. Assumes the plane is perpendicular to
/// the cylinder axis (the configuration the closed-fillet path
/// produces).
pub(crate) fn plane_cylinder_intersection_circle(
    plane: &PlaneGeometry,
    cyl: &CylinderGeometry,
) -> Result<(Point3, f64), BlendError> {
    let denom = cyl.axis.dot(&plane.normal);
    if denom.abs() < 1e-9 {
        return Err(BlendError::InvalidParameter {
            parameter: "plane,cylinder".to_string(),
            value: "parallel".to_string(),
            reason: "cylinder axis is parallel to the plane — no rim circle".to_string(),
        });
    }
    let diff = Vec3::new(
        plane.origin.x - cyl.origin.x,
        plane.origin.y - cyl.origin.y,
        plane.origin.z - cyl.origin.z,
    );
    let t = diff.dot(&plane.normal) / denom;
    let centre = Point3::new(
        cyl.origin.x + t * cyl.axis.x,
        cyl.origin.y + t * cyl.axis.y,
        cyl.origin.z + t * cyl.axis.z,
    );
    Ok((centre, cyl.radius))
}

/// Compute the new vertex position on the rim circle that
/// preserves the existing vertex's angular position around the
/// cylinder axis. Project the old position onto the plane
/// containing the rim circle, then renormalise to the rim radius.
pub(crate) fn compute_rim_vertex_position(
    old_pos: &Point3,
    rim_centre: &Point3,
    cyl: &CylinderGeometry,
    rim_radius: f64,
) -> Point3 {
    let v = Vec3::new(
        old_pos.x - rim_centre.x,
        old_pos.y - rim_centre.y,
        old_pos.z - rim_centre.z,
    );
    let along_axis = v.dot(&cyl.axis);
    let perp = Vec3::new(
        v.x - along_axis * cyl.axis.x,
        v.y - along_axis * cyl.axis.y,
        v.z - along_axis * cyl.axis.z,
    );
    let perp_len = perp.length();
    if perp_len < 1e-12 {
        return Point3::new(rim_centre.x + rim_radius, rim_centre.y, rim_centre.z);
    }
    let scale = rim_radius / perp_len;
    Point3::new(
        rim_centre.x + perp.x * scale,
        rim_centre.y + perp.y * scale,
        rim_centre.z + perp.z * scale,
    )
}

/// Walk a loop looking for a half-edge whose origin is `vertex`
/// but skipping `skip` (the rim half-edge). Returns the edge of
/// the matching half-edge — typically the seam edge.
pub(crate) fn find_seam_edge_at_vertex_excluding(
    brep: &BRep,
    loop_start: HalfEdgeHandle,
    skip: HalfEdgeHandle,
    vertex: VertexHandle,
) -> Option<EdgeHandle> {
    let mut cursor = loop_start;
    let mut iters = 0;
    loop {
        iters += 1;
        if iters > 10_000 {
            return None;
        }
        if cursor != skip {
            let h = brep.half_edges.get(cursor)?;
            if h.origin == vertex {
                return Some(h.edge);
            }
        }
        let next = brep.half_edges.get(cursor)?.next;
        cursor = next;
        if cursor == loop_start {
            break;
        }
    }
    None
}

/// Replace `edge`'s curve with a fresh straight line spanning its
/// two current endpoint vertex positions.
pub(crate) fn rebuild_line_curve(brep: &mut BRep, edge: EdgeHandle) {
    let Some(edata) = brep.edges.get(edge) else {
        return;
    };
    let he_h = edata.half_edge;
    let Some(he) = brep.half_edges.get(he_h).cloned() else {
        return;
    };
    let Some(twin) = brep.half_edges.get(he.twin).cloned() else {
        return;
    };
    let v0 = brep.vertices[he.origin].position;
    let v1 = brep.vertices[twin.origin].position;
    let dir = Vec3::new(v1.x - v0.x, v1.y - v0.y, v1.z - v0.z);
    if dir.length() < 1e-12 {
        return;
    }
    brep.edges[edge].curve = Some(CurveOnEdge {
        geometry: CurveGeometry::Line {
            origin: v0,
            direction: dir,
        },
        t_start: 0.0,
        t_end: 1.0,
    });
}

/// Build a NURBS approximation of a circle in 3D given the
/// centre, radius, and normal direction. Uses the rational
/// quadratic 9-control-point representation that popkern's other
/// circle helpers use.
pub(crate) fn build_rim_circle_nurbs(centre: &Point3, radius: f64, normal: &Vec3) -> NurbsCurve {
    let n = if normal.length() > 1e-12 {
        let s = 1.0 / normal.length();
        Vec3::new(normal.x * s, normal.y * s, normal.z * s)
    } else {
        Vec3::new(0.0, 0.0, 1.0)
    };
    let u_seed = if n.x.abs() < 0.9 {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let dot = u_seed.dot(&n);
    let u_raw = Vec3::new(
        u_seed.x - dot * n.x,
        u_seed.y - dot * n.y,
        u_seed.z - dot * n.z,
    );
    let u_len = u_raw.length();
    let u = if u_len > 1e-12 {
        Vec3::new(u_raw.x / u_len, u_raw.y / u_len, u_raw.z / u_len)
    } else {
        Vec3::new(1.0, 0.0, 0.0)
    };
    let v = n.cross(&u);

    let r = radius;
    let sqrt2 = std::f64::consts::FRAC_1_SQRT_2;
    let pt = |a: f64, b: f64| -> Point3 {
        Point3::new(
            centre.x + a * r * u.x + b * r * v.x,
            centre.y + a * r * u.y + b * r * v.y,
            centre.z + a * r * u.z + b * r * v.z,
        )
    };
    let cps = vec![
        pt(1.0, 0.0),
        pt(1.0, 1.0),
        pt(0.0, 1.0),
        pt(-1.0, 1.0),
        pt(-1.0, 0.0),
        pt(-1.0, -1.0),
        pt(0.0, -1.0),
        pt(1.0, -1.0),
        pt(1.0, 0.0),
    ];
    let weights = vec![1.0, sqrt2, 1.0, sqrt2, 1.0, sqrt2, 1.0, sqrt2, 1.0];
    let knots = vec![
        0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
    ];
    NurbsCurve {
        degree: 2,
        knots,
        control_points: cps,
        weights,
    }
}
