//! Dispatch bridge between [`MeasureQuery`] (entity names) and
//! `popkern_measure` (typed handles).

use popkern_core::naming::EntityName;
use popkern_core::topology::BRep;

use crate::error::ReplayError;
use crate::measure_query::{BoundExtent, Component, MeasureContext, MeasureQuery};
use crate::resolve::{resolve_edge, resolve_face, resolve_vertex};

/// Resolve entity names in a [`MeasureQuery`], dispatch to `popkern_measure`,
/// and return a scalar result.
pub(crate) fn resolve_and_evaluate(
    query: &MeasureQuery,
    ctx: &MeasureContext<'_>,
) -> Result<f64, ReplayError> {
    let brep = ctx
        .bodies
        .get(&ctx.active_body)
        .ok_or_else(|| ReplayError::FeatureFailed {
            label: "measure".into(),
            error: "no BRep for active body".into(),
        })?;

    match query {
        MeasureQuery::EdgeLength(name) => {
            let h = resolve_edge(name, brep)?;
            popkern_measure::edge_length(brep, h).map_err(measure_err)
        }
        MeasureQuery::EdgeTangent(name, comp) => {
            let h = resolve_edge(name, brep)?;
            let v = popkern_measure::edge_tangent(brep, h).map_err(measure_err)?;
            Ok(extract_component(v.x, v.y, v.z, comp))
        }
        MeasureQuery::FaceArea(name) => {
            let h = resolve_face(name, brep)?;
            popkern_measure::face_area(brep, h, ctx.tolerance).map_err(measure_err)
        }
        MeasureQuery::FacePerimeter(name) => {
            let h = resolve_face(name, brep)?;
            popkern_measure::face_perimeter(brep, h).map_err(measure_err)
        }
        MeasureQuery::FaceNormal(name, comp) => {
            let h = resolve_face(name, brep)?;
            let n = popkern_measure::face_normal(brep, h).map_err(measure_err)?;
            Ok(extract_component(n.x, n.y, n.z, comp))
        }
        MeasureQuery::FaceCentroid(name, comp) => {
            let h = resolve_face(name, brep)?;
            let c = popkern_measure::face_centroid(brep, h, ctx.tolerance).map_err(measure_err)?;
            Ok(extract_component(c.x, c.y, c.z, comp))
        }
        MeasureQuery::FaceRadius(name) => {
            let h = resolve_face(name, brep)?;
            popkern_measure::face_radius(brep, h).map_err(measure_err)
        }
        MeasureQuery::FaceMinCurvature(name) => {
            let h = resolve_face(name, brep)?;
            let (k_min, _) = popkern_measure::face_curvature(brep, h).map_err(measure_err)?;
            Ok(k_min)
        }
        MeasureQuery::FaceMaxCurvature(name) => {
            let h = resolve_face(name, brep)?;
            let (_, k_max) = popkern_measure::face_curvature(brep, h).map_err(measure_err)?;
            Ok(k_max)
        }
        MeasureQuery::VertexPosition(name, comp) => {
            let h = resolve_vertex(name, brep)?;
            let p = popkern_measure::vertex_position(brep, h).map_err(measure_err)?;
            Ok(extract_component(p.x, p.y, p.z, comp))
        }
        MeasureQuery::MinDistance(name_a, name_b) => {
            let ha = resolve_to_topo_handle(name_a, brep)?;
            let hb = resolve_to_topo_handle(name_b, brep)?;
            popkern_measure::min_distance(brep, ha, brep, hb, ctx.tolerance).map_err(measure_err)
        }
        MeasureQuery::DihedralAngle(name_a, name_b) => {
            let fa = resolve_face(name_a, brep)?;
            let fb = resolve_face(name_b, brep)?;
            popkern_measure::dihedral_angle(brep, fa, fb).map_err(measure_err)
        }
        MeasureQuery::EdgeAngle(name_a, name_b) => {
            let ea = resolve_edge(name_a, brep)?;
            let eb = resolve_edge(name_b, brep)?;
            popkern_measure::edge_angle(brep, ea, eb).map_err(measure_err)
        }
        MeasureQuery::BodyVolume => {
            popkern_measure::body_volume(brep, ctx.tolerance).map_err(measure_err)
        }
        MeasureQuery::BodySurfaceArea => {
            popkern_measure::body_surface_area(brep, ctx.tolerance).map_err(measure_err)
        }
        MeasureQuery::BodyCentroid(comp) => {
            let c = popkern_measure::body_centroid(brep, ctx.tolerance).map_err(measure_err)?;
            Ok(extract_component(c.x, c.y, c.z, comp))
        }
        MeasureQuery::BodyBound(extent) => {
            let bb = popkern_measure::body_bounding_box(brep).map_err(measure_err)?;
            Ok(match extent {
                BoundExtent::MinX => bb.min.x,
                BoundExtent::MinY => bb.min.y,
                BoundExtent::MinZ => bb.min.z,
                BoundExtent::MaxX => bb.max.x,
                BoundExtent::MaxY => bb.max.y,
                BoundExtent::MaxZ => bb.max.z,
                BoundExtent::SizeX => bb.size_x(),
                BoundExtent::SizeY => bb.size_y(),
                BoundExtent::SizeZ => bb.size_z(),
            })
        }
        MeasureQuery::BodyMomentOfInertia(axis) => {
            popkern_measure::body_moment_of_inertia(brep, *axis, ctx.tolerance).map_err(measure_err)
        }
    }
}

/// Extract the requested scalar component from (x, y, z).
fn extract_component(x: f64, y: f64, z: f64, comp: &Component) -> f64 {
    match comp {
        Component::X => x,
        Component::Y => y,
        Component::Z => z,
    }
}

/// Try to resolve an entity name to any topology handle type.
fn resolve_to_topo_handle(
    name: &EntityName,
    brep: &BRep,
) -> Result<popkern_measure::TopoHandle, ReplayError> {
    use popkern_core::naming_resolve::{Resolution, ResolvedEntity, resolve_name};
    match resolve_name(name, brep) {
        Resolution::Exact(ResolvedEntity::Vertex(h)) => Ok(popkern_measure::TopoHandle::Vertex(h)),
        Resolution::Exact(ResolvedEntity::Edge(h)) => Ok(popkern_measure::TopoHandle::Edge(h)),
        Resolution::Exact(ResolvedEntity::Face(h)) => Ok(popkern_measure::TopoHandle::Face(h)),
        Resolution::Split(entities) => {
            if let Some(e) = entities.first() {
                match e {
                    ResolvedEntity::Vertex(h) => Ok(popkern_measure::TopoHandle::Vertex(*h)),
                    ResolvedEntity::Edge(h) => Ok(popkern_measure::TopoHandle::Edge(*h)),
                    ResolvedEntity::Face(h) => Ok(popkern_measure::TopoHandle::Face(*h)),
                }
            } else {
                Err(ReplayError::ResolutionFailed(name.clone()))
            }
        }
        Resolution::NotFound => Err(ReplayError::ResolutionFailed(name.clone())),
    }
}

/// Convert a [`popkern_measure::MeasureError`] into a [`ReplayError`].
fn measure_err(e: popkern_measure::MeasureError) -> ReplayError {
    ReplayError::FeatureFailed {
        label: "measure".into(),
        error: e.to_string(),
    }
}
