//! A B-spline surface closed round u whose two v boundaries collapse to
//! poles (the standard biquadratic sphere, or that net scaled into an
//! ellipsoid): the solid it bounds, on the core sphere's topology (one
//! face, one seam from pole to pole, the two pole vertices, no edge at a
//! pole). The face carries the surface; the seam its exact u = 0
//! iso-curve; the seam's two uses the chart lines u = u1 (running up)
//! and u = u0 (running down), so the loop winds round the chart as the
//! constructor's sphere does; and the body faces out of the volume it
//! bounds.

use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::geometry::nurbs_surface::NurbsSurface;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry, SurfaceOnFace};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use super::failed;
use crate::error::ReplayError;

/// The ball `surface` bounds.
pub(super) fn ball(
    shape: &str,
    surface: &NurbsSurface,
    tolerance: Tolerance,
) -> Result<BRep, ReplayError> {
    let (u0, u1) = (
        surface.knots_u[0],
        surface.knots_u[surface.knots_u.len() - 1],
    );
    let (v0, v1) = (
        surface.knots_v[0],
        surface.knots_v[surface.knots_v.len() - 1],
    );
    let mut brep = BRep::sphere(1.0, tolerance).map_err(|e| failed(shape, "ball", e))?;
    let (south, north) = (surface.point_at(u0, v0), surface.point_at(u0, v1));
    let Some((face, _)) = brep.faces.iter().next() else {
        return Err(failed(shape, "ball", "the sphere has no face"));
    };
    let Some((seam, _)) = brep.edges.iter().next() else {
        return Err(failed(shape, "ball", "the sphere has no seam"));
    };
    // The constructor's poles: south below, north above.
    let vertices: Vec<_> = brep.vertices.keys().collect();
    for v in vertices {
        let z = brep.vertices[v].position.z;
        brep.vertices[v].position = if z < 0.0 { south } else { north };
    }
    brep.faces[face].surface = Some(SurfaceOnFace {
        geometry: SurfaceGeometry::Nurbs(surface.clone()),
        orientation: true,
    });
    let meridian = super::boundary_curves(surface)[0].clone();
    let up = brep.edges[seam].half_edge;
    let starts_south = brep
        .half_edges
        .get(up)
        .and_then(|h| brep.vertices.get(h.origin))
        .is_some_and(|v| v.position.distance_to(&south) <= tolerance.linear());
    let (t_start, t_end) = if starts_south { (v0, v1) } else { (v1, v0) };
    brep.edges[seam].curve = Some(CurveOnEdge {
        geometry: CurveGeometry::Nurbs(meridian),
        t_start,
        t_end,
    });
    let uses: Vec<_> = brep
        .half_edges
        .iter()
        .filter(|(_, h)| h.edge == seam)
        .map(|(handle, h)| (handle, h.origin))
        .collect();
    for (he, origin) in uses {
        let from_south = brep.vertices[origin].position.distance_to(&south) <= tolerance.linear();
        brep.half_edges[he].pcurve = Some(if from_south {
            Pcurve::linear((u1, v0), (u1, v1), v0, v1)
        } else {
            Pcurve::linear((u0, v1), (u0, v0), v1, v0)
        });
    }
    let faces: Vec<_> = brep.faces.keys().collect();
    if brep.folded_wedge_volume(&faces).is_some_and(|v| v < 0.0) {
        brep.reverse_body_orientation();
    }
    let errors = brep.validate_all();
    if !errors.is_empty() {
        return Err(failed(shape, "ball", format!("{errors:?}")));
    }
    Ok(brep)
}
