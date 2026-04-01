//! Resolution of symbolic entity references at replay time.
//!
//! These functions take [`EntityName`] references from [`Feature`] variants
//! and resolve them against the current BRep to get concrete topology handles.
//!
//! [`Feature`]: crate::feature::Feature

use popkern_core::geometry::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{CurveGeometry, Point3, Vec3};
use popkern_core::naming::EntityName;
use popkern_core::naming_resolve::{Resolution, ResolvedEntity, resolve_name};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};
use popkern_sketch::SketchPlane;

use crate::error::ReplayError;
use crate::feature::{AxisDesc, ExtrudeEnd, SketchPlaneDesc, Vec3Desc};
use crate::parameter::{ParameterId, ParameterTable};

/// Resolve an [`EntityName`] to a [`FaceHandle`].
///
/// For split resolutions, returns the first face found in the successor list.
pub(crate) fn resolve_face(name: &EntityName, brep: &BRep) -> Result<FaceHandle, ReplayError> {
    match resolve_name(name, brep) {
        Resolution::Exact(ResolvedEntity::Face(h)) => Ok(h),
        Resolution::Split(entities) => {
            for e in entities {
                if let ResolvedEntity::Face(h) = e {
                    return Ok(h);
                }
            }
            Err(ReplayError::ResolutionFailed(name.clone()))
        }
        _ => Err(ReplayError::ResolutionFailed(name.clone())),
    }
}

/// Resolve an [`EntityName`] to an [`EdgeHandle`].
///
/// For split resolutions, returns the first edge found in the successor list.
pub(crate) fn resolve_edge(name: &EntityName, brep: &BRep) -> Result<EdgeHandle, ReplayError> {
    match resolve_name(name, brep) {
        Resolution::Exact(ResolvedEntity::Edge(h)) => Ok(h),
        Resolution::Split(entities) => {
            for e in entities {
                if let ResolvedEntity::Edge(h) = e {
                    return Ok(h);
                }
            }
            Err(ReplayError::ResolutionFailed(name.clone()))
        }
        _ => Err(ReplayError::ResolutionFailed(name.clone())),
    }
}

/// Resolve an [`EntityName`] to a [`VertexHandle`].
///
/// For split resolutions, returns the first vertex found in the successor list.
pub(crate) fn resolve_vertex(name: &EntityName, brep: &BRep) -> Result<VertexHandle, ReplayError> {
    match resolve_name(name, brep) {
        Resolution::Exact(ResolvedEntity::Vertex(h)) => Ok(h),
        Resolution::Split(entities) => {
            for e in entities {
                if let ResolvedEntity::Vertex(h) = e {
                    return Ok(h);
                }
            }
            Err(ReplayError::ResolutionFailed(name.clone()))
        }
        _ => Err(ReplayError::ResolutionFailed(name.clone())),
    }
}

/// Resolve multiple [`EntityName`]s to [`EdgeHandle`]s.
///
/// For split resolutions, collects ALL edge pieces (e.g. fillet both halves
/// of a split edge).
pub(crate) fn resolve_edges(
    names: &[EntityName],
    brep: &BRep,
) -> Result<Vec<EdgeHandle>, ReplayError> {
    let mut edges = Vec::new();
    for name in names {
        match resolve_name(name, brep) {
            Resolution::Exact(ResolvedEntity::Edge(h)) => edges.push(h),
            Resolution::Split(entities) => {
                for e in entities {
                    if let ResolvedEntity::Edge(h) = e {
                        edges.push(h);
                    }
                }
            }
            _ => return Err(ReplayError::ResolutionFailed(name.clone())),
        }
    }
    Ok(edges)
}

/// Resolve multiple [`EntityName`]s to [`FaceHandle`]s.
///
/// For split resolutions, collects ALL face pieces.
pub(crate) fn resolve_faces(
    names: &[EntityName],
    brep: &BRep,
) -> Result<Vec<FaceHandle>, ReplayError> {
    let mut faces = Vec::new();
    for name in names {
        match resolve_name(name, brep) {
            Resolution::Exact(ResolvedEntity::Face(h)) => faces.push(h),
            Resolution::Split(entities) => {
                for e in entities {
                    if let ResolvedEntity::Face(h) = e {
                        faces.push(h);
                    }
                }
            }
            _ => return Err(ReplayError::ResolutionFailed(name.clone())),
        }
    }
    Ok(faces)
}

/// Resolve a [`ParameterId`] to its `f64` value via the parameter table.
fn resolve_param(table: &ParameterTable, id: ParameterId) -> Result<f64, ReplayError> {
    table
        .get(id)
        .map(|p| p.value)
        .ok_or_else(|| ReplayError::FeatureFailed {
            label: "parameter".into(),
            error: format!("parameter {id:?} not found"),
        })
}

/// Resolve a [`SketchPlaneDesc`] to a concrete [`SketchPlane`].
pub(crate) fn resolve_sketch_plane(
    desc: &SketchPlaneDesc,
    brep: &BRep,
    param_table: &ParameterTable,
) -> Result<SketchPlane, ReplayError> {
    match desc {
        SketchPlaneDesc::XY => Ok(SketchPlane::xy()),
        SketchPlaneDesc::XZ => Ok(SketchPlane::xz()),
        SketchPlaneDesc::YZ => Ok(SketchPlane::yz()),
        SketchPlaneDesc::OnFace(name) => {
            let face = resolve_face(name, brep)?;
            SketchPlane::from_face(brep, face).map_err(|e| ReplayError::FeatureFailed {
                label: "sketch plane".into(),
                error: e.to_string(),
            })
        }
        SketchPlaneDesc::OnFaceOffset(name, offset_param) => {
            let face = resolve_face(name, brep)?;
            let offset = resolve_param(param_table, *offset_param)?;
            SketchPlane::from_face_offset(brep, face, offset).map_err(|e| {
                ReplayError::FeatureFailed {
                    label: "sketch plane".into(),
                    error: e.to_string(),
                }
            })
        }
        SketchPlaneDesc::Custom { origin, normal } => {
            let o = Point3::new(origin[0], origin[1], origin[2]);
            let n = Vec3::new(normal[0], normal[1], normal[2]);
            SketchPlane::from_normal(o, n).map_err(|e| ReplayError::FeatureFailed {
                label: "sketch plane".into(),
                error: e.to_string(),
            })
        }
    }
}

/// Resolve a [`Vec3Desc`] to a concrete [`Vec3`].
///
/// For [`Vec3Desc::FaceNormal`], the normal is computed from the face's
/// surface geometry. Planar surfaces return a constant normal; other surface
/// types evaluate the normal at the centroid of the face's outer loop.
pub(crate) fn resolve_direction(desc: &Vec3Desc, brep: &BRep) -> Result<Vec3, ReplayError> {
    match desc {
        Vec3Desc::Literal(arr) => Ok(Vec3::new(arr[0], arr[1], arr[2])),
        Vec3Desc::FaceNormal(name) => {
            let face_handle = resolve_face(name, brep)?;
            let sof = brep.faces[face_handle].surface.as_ref().ok_or_else(|| {
                ReplayError::FeatureFailed {
                    label: "direction".into(),
                    error: "face has no surface geometry".into(),
                }
            })?;
            // For planes the normal is constant — evaluate at (0, 0).
            // For other surfaces, compute the centroid and use bounded UV.
            let normal = match &sof.geometry {
                popkern_core::geometry::SurfaceGeometry::Plane { normal, .. } => *normal,
                other => {
                    let positions = brep.face_outer_loop_positions(face_handle);
                    if positions.is_empty() {
                        return Err(ReplayError::FeatureFailed {
                            label: "direction".into(),
                            error: "face has no vertices".into(),
                        });
                    }
                    let (u_range, v_range) =
                        popkern_core::geometry::uv_projection::tight_bounded_parameter_range(
                            &positions, other,
                        );
                    let u_mid = (u_range.0 + u_range.1) / 2.0;
                    let v_mid = (v_range.0 + v_range.1) / 2.0;
                    other.normal_at(u_mid, v_mid)
                }
            };
            // Flip normal if the surface orientation is reversed.
            let normal = if sof.orientation { normal } else { -normal };
            Ok(normal)
        }
        Vec3Desc::EdgeDirection(name) => {
            let edge_handle = resolve_edge(name, brep)?;
            let curve_on_edge = brep.edges[edge_handle].curve.as_ref().ok_or_else(|| {
                ReplayError::FeatureFailed {
                    label: "direction".into(),
                    error: "edge has no curve geometry".into(),
                }
            })?;
            let t = (curve_on_edge.t_start + curve_on_edge.t_end) / 2.0;
            Ok(curve_on_edge.geometry.tangent_at(t))
        }
        Vec3Desc::Negated(inner) => {
            let v = resolve_direction(inner, brep)?;
            Ok(-v)
        }
    }
}

/// Compute the AABB diagonal of a BRep by iterating all vertex positions.
///
/// Returns 0.0 for an empty BRep (no vertices).
fn brep_aabb_diagonal(brep: &BRep) -> f64 {
    let mut min = [f64::MAX; 3];
    let mut max = [f64::MIN; 3];
    let mut found = false;
    for (_, v) in &brep.vertices {
        let p = v.position;
        min[0] = min[0].min(p.x);
        min[1] = min[1].min(p.y);
        min[2] = min[2].min(p.z);
        max[0] = max[0].max(p.x);
        max[1] = max[1].max(p.y);
        max[2] = max[2].max(p.z);
        found = true;
    }
    if !found {
        return 0.0;
    }
    let dx = max[0] - min[0];
    let dy = max[1] - min[1];
    let dz = max[2] - min[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// Resolve an [`AxisDesc`] to a concrete axis (origin point + unit direction).
///
/// - [`AxisDesc::Edge`]: resolves the edge, reads start and end vertex
///   positions from the half-edge topology, and computes the unit direction.
///   Returns an error if the edge is degenerate (too short).
/// - [`AxisDesc::Literal`]: converts the stored arrays to [`Point3`]/[`Vec3`]
///   and normalises the direction.  Returns an error if the direction is zero.
pub(crate) fn resolve_axis(desc: &AxisDesc, brep: &BRep) -> Result<(Point3, Vec3), ReplayError> {
    match desc {
        AxisDesc::Edge(name) => {
            let edge_handle = resolve_edge(name, brep)?;
            let edge = &brep.edges[edge_handle];
            let he = &brep.half_edges[edge.half_edge];
            let twin = &brep.half_edges[he.twin];
            let start = brep.vertices[he.origin].position;
            let end = brep.vertices[twin.origin].position;
            let dir = end - start;
            let unit = dir.normalized().ok_or_else(|| ReplayError::FeatureFailed {
                label: "axis".into(),
                error: "axis edge is degenerate (too short)".into(),
            })?;
            Ok((start, unit))
        }
        AxisDesc::Literal { origin, direction } => {
            let o = Point3::new(origin[0], origin[1], origin[2]);
            let d = Vec3::new(direction[0], direction[1], direction[2]);
            let unit = d.normalized().ok_or_else(|| ReplayError::FeatureFailed {
                label: "axis".into(),
                error: "axis direction is zero-length".into(),
            })?;
            Ok((o, unit))
        }
    }
}

/// Resolve an [`ExtrudeEnd`] to a concrete extrude distance in world units.
///
/// - [`ExtrudeEnd::Blind`]: looks up the distance parameter value.
/// - [`ExtrudeEnd::ThroughAll`]: uses the host body AABB diagonal x 2.0 as a
///   safe overshoot; the boolean will clip the result.
/// - [`ExtrudeEnd::UpToSurface`]: intersects a ray (from `extrude_origin`
///   along `extrude_dir`) with the target face's surface and returns the
///   distance to the nearest forward intersection.
/// - [`ExtrudeEnd::OffsetFromSurface`]: delegates to `UpToSurface` and adds
///   the signed offset distance.
/// - [`ExtrudeEnd::UpToVertex`]: projects the target vertex position onto the
///   extrude direction to compute the distance.
/// - [`ExtrudeEnd::MidPlane`]: looks up the half_distance parameter value;
///   callers are responsible for the bidirectional extrude logic.
pub(crate) fn resolve_extrude_distance(
    end: &ExtrudeEnd,
    extrude_origin: Point3,
    extrude_dir: Vec3,
    host_brep: Option<&BRep>,
    brep: &BRep,
    param_table: &ParameterTable,
    tolerance: &popkern_core::tolerance::Tolerance,
) -> Result<f64, ReplayError> {
    match end {
        ExtrudeEnd::Blind { distance } => resolve_param(param_table, *distance),
        ExtrudeEnd::ThroughAll => {
            let diag = host_brep.map_or_else(|| brep_aabb_diagonal(brep), brep_aabb_diagonal);
            Ok(diag * 2.0)
        }
        ExtrudeEnd::UpToSurface { face } => {
            let target_brep = host_brep.unwrap_or(brep);
            let face_handle = resolve_face(face, target_brep)?;
            let sof = target_brep.faces[face_handle]
                .surface
                .as_ref()
                .ok_or_else(|| ReplayError::FeatureFailed {
                    label: String::new(),
                    error: "target face has no surface geometry".into(),
                })?;
            let surface = &sof.geometry;

            // Build ray as a Line curve geometry.
            let ray_curve = CurveGeometry::Line {
                origin: extrude_origin,
                direction: extrude_dir,
            };

            // Get bounded UV range for the target surface.
            let positions = target_brep.face_outer_loop_positions(face_handle);
            let (u_range, v_range) =
                popkern_core::geometry::uv_projection::bounded_parameter_range(&positions, surface);

            // Intersect ray with surface.
            let result = popkern_core::geometry::intersection::intersect_surface_curve(
                surface,
                u_range,
                v_range,
                &ray_curve,
                (0.0, 1e6),
                tolerance,
            )?;

            // Find the nearest forward intersection.
            let mut min_t = f64::MAX;
            for isect in &result.intersections {
                if let popkern_core::geometry::intersection::SurfaceCurveIntersection::Point {
                    t,
                    ..
                } = isect
                    && *t > tolerance.linear
                    && *t < min_t
                {
                    min_t = *t;
                }
            }

            if min_t == f64::MAX {
                return Err(ReplayError::FeatureFailed {
                    label: String::new(),
                    error: "extrude ray does not intersect target surface".into(),
                });
            }

            Ok(min_t)
        }
        ExtrudeEnd::OffsetFromSurface { face, offset } => {
            let base_dist = resolve_extrude_distance(
                &ExtrudeEnd::UpToSurface { face: face.clone() },
                extrude_origin,
                extrude_dir,
                host_brep,
                brep,
                param_table,
                tolerance,
            )?;
            let off = resolve_param(param_table, *offset)?;
            Ok(base_dist + off)
        }
        ExtrudeEnd::UpToVertex { vertex } => {
            let target_brep = host_brep.unwrap_or(brep);
            let vh = resolve_vertex(vertex, target_brep)?;
            let pos = target_brep.vertices[vh].position;
            let diff = pos - extrude_origin;
            let dist = diff.dot(&extrude_dir);
            if dist <= tolerance.linear {
                return Err(ReplayError::FeatureFailed {
                    label: String::new(),
                    error: "target vertex is behind or at extrude origin".into(),
                });
            }
            Ok(dist)
        }
        ExtrudeEnd::MidPlane { half_distance } => resolve_param(param_table, *half_distance),
    }
}

#[path = "resolve_revolve.rs"]
mod resolve_revolve;
pub(crate) use resolve_revolve::{RevolveContext, resolve_revolve_angle};
