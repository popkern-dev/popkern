//! Resolution of symbolic entity references at replay time.
//!
//! These functions take [`EntityName`] references from [`Feature`] variants
//! and resolve them against the current BRep to get concrete topology handles.
//!
//! [`Feature`]: crate::feature::Feature

use popkern_core::geometry::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::naming::EntityName;
use popkern_core::naming_resolve::{Resolution, ResolvedEntity, resolve_name};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};
use popkern_sketch::SketchPlane;

use crate::error::ReplayError;
use crate::feature::{SketchPlaneDesc, Vec3Desc};
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
            let t = curve_on_edge.t_start;
            Ok(curve_on_edge.geometry.tangent_at(t))
        }
        Vec3Desc::Negated(inner) => {
            let v = resolve_direction(inner, brep)?;
            Ok(-v)
        }
    }
}
