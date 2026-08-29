use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{BRep, FaceHandle};

use super::super::constructive_topology;
use super::partition_cylinder_geometry::cylinder_radial_vector;

pub(super) fn generator_anchor_blocker(brep: &BRep, face: FaceHandle, endpoint: Point3) -> String {
    let Some(surface) = brep.faces.get(face).and_then(|face| face.surface.as_ref()) else {
        return "target face has no surface".to_string();
    };
    let SurfaceGeometry::Cylinder {
        origin,
        axis,
        radius,
    } = surface.geometry
    else {
        return "target face is not a cylinder".to_string();
    };
    let Some(endpoint_radial) = cylinder_radial_vector(endpoint, origin, axis) else {
        return "endpoint has degenerate radial vector to cylinder axis".to_string();
    };
    let endpoint_surface_error = (endpoint_radial.length() - radius).abs();
    let Some(endpoint_unit) = endpoint_radial.normalized() else {
        return "endpoint radial vector cannot be normalized".to_string();
    };
    let Some(vertices) = constructive_topology::collect_outer_loop_vertices(brep, face) else {
        return "target cylinder face has no walkable outer loop".to_string();
    };
    let mut reports = vertices
        .into_iter()
        .filter_map(|vertex| {
            let point = brep.vertices.get(vertex)?.position;
            let radial = cylinder_radial_vector(point, origin, axis)?;
            let radial_distance = radial.length();
            let radial_unit = radial.normalized()?;
            let surface_error = (radial_distance - radius).abs();
            let radial_alignment = radial_unit.dot(&endpoint_unit);
            Some((vertex, surface_error, radial_alignment, point))
        })
        .collect::<Vec<_>>();
    reports.sort_by(|a, b| {
        let a_score = a.1 + (1.0 - a.2).abs() * radius.abs().max(1.0e-9);
        let b_score = b.1 + (1.0 - b.2).abs() * radius.abs().max(1.0e-9);
        a_score
            .partial_cmp(&b_score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let nearest = reports
        .iter()
        .take(4)
        .map(|(vertex, surface_error, radial_alignment, point)| {
            format!(
                "{vertex:?}@({:.6},{:.6},{:.6}) surface_error={surface_error:.3e}, \
                 radial_alignment={radial_alignment:.6}",
                point.x, point.y, point.z
            )
        })
        .collect::<Vec<_>>()
        .join("; ");
    format!(
        "endpoint cylinder surface_error={endpoint_surface_error:.3e}; nearest radial loop \
         candidates [{}]. A commercial rebuild must materialize anchor vertices on the stored \
         cylinder surface before wiring generator trims.",
        nearest
    )
}
