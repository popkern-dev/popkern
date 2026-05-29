//! Support surface authority selection for exact-boundary ribbons.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::healing::deviation::edge_surface_deviation;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

use super::super::boundary::PreparedExactBoundaryCurve;
use super::super::sector::ExactBoundaryOuterSide;
use super::ribbon_surface_identity::geometrically_unique_surface;

pub(super) fn select_support_surface(
    brep: &BRep,
    cap_face: FaceHandle,
    support_face: FaceHandle,
    side: ExactBoundaryOuterSide,
    edge: EdgeHandle,
    boundary: &PreparedExactBoundaryCurve,
    tolerance: f64,
) -> Option<SurfaceGeometry> {
    let mut surface = face_support_surface(brep, support_face)
        .or_else(|| unique_analytic_surface_for_edge(brep, cap_face, support_face, side, edge))?;
    if !boundary_fits_surface(boundary, &surface, tolerance)
        && let Some(authority) =
            unique_analytic_surface_for_boundary(brep, cap_face, side, boundary)
    {
        surface = authority;
    }
    Some(surface)
}

fn face_support_surface(brep: &BRep, support_face: FaceHandle) -> Option<SurfaceGeometry> {
    brep.faces
        .get(support_face)
        .and_then(|face_data| face_data.surface.as_ref())
        .map(|surface| surface.geometry.clone())
}

fn unique_analytic_surface_for_edge(
    brep: &BRep,
    cap_face: FaceHandle,
    support_face: FaceHandle,
    side: ExactBoundaryOuterSide,
    edge: EdgeHandle,
) -> Option<SurfaceGeometry> {
    let edge = brep.edges.get(edge)?;
    let curve = edge.curve.as_ref()?;
    let tolerance = brep.tolerance.linear().max(1.0e-6);
    let candidates = brep
        .faces
        .iter()
        .filter_map(|(face, face_data)| {
            if face == cap_face || face == support_face {
                return None;
            }
            let surface = face_data.surface.as_ref()?;
            if matches!(
                surface.geometry,
                SurfaceGeometry::Gregory(_) | SurfaceGeometry::Nurbs(_)
            ) {
                return None;
            }
            let residual = edge_surface_deviation(
                &curve.geometry,
                (curve.t_start, curve.t_end),
                &surface.geometry,
                17,
            )
            .max_deviation;
            (residual <= tolerance).then_some(surface.geometry.clone())
        })
        .collect::<Vec<_>>();
    preferred_or_unique_surface(candidates, side, tolerance)
}

fn unique_analytic_surface_for_boundary(
    brep: &BRep,
    cap_face: FaceHandle,
    side: ExactBoundaryOuterSide,
    boundary: &PreparedExactBoundaryCurve,
) -> Option<SurfaceGeometry> {
    let tolerance = brep.tolerance.linear().max(1.0e-6);
    let range = boundary.curve.parameter_range();
    let candidates = brep
        .faces
        .iter()
        .filter_map(|(face, face_data)| {
            if face == cap_face {
                return None;
            }
            let surface = face_data.surface.as_ref()?;
            if matches!(
                surface.geometry,
                SurfaceGeometry::Gregory(_) | SurfaceGeometry::Nurbs(_)
            ) {
                return None;
            }
            let residual =
                edge_surface_deviation(&boundary.curve, range, &surface.geometry, 17).max_deviation;
            (residual <= tolerance).then_some(surface.geometry.clone())
        })
        .collect::<Vec<_>>();
    preferred_or_unique_surface(candidates, side, tolerance)
}

fn preferred_or_unique_surface(
    candidates: Vec<SurfaceGeometry>,
    side: ExactBoundaryOuterSide,
    tolerance: f64,
) -> Option<SurfaceGeometry> {
    let mut preferred = candidates
        .iter()
        .filter(|surface| match side {
            ExactBoundaryOuterSide::FilletEnd => {
                matches!(surface, SurfaceGeometry::Cylinder { .. })
            }
            ExactBoundaryOuterSide::FaceEdge => matches!(surface, SurfaceGeometry::Plane { .. }),
        })
        .cloned()
        .collect::<Vec<_>>();
    if preferred.len() == 1 {
        return preferred.pop();
    }
    if let Some(surface) = geometrically_unique_surface(&preferred, tolerance) {
        return Some(surface);
    }
    if candidates.len() == 1 {
        candidates.into_iter().next()
    } else {
        geometrically_unique_surface(&candidates, tolerance)
    }
}

fn boundary_fits_surface(
    boundary: &PreparedExactBoundaryCurve,
    surface: &SurfaceGeometry,
    tolerance: f64,
) -> bool {
    let range = boundary.curve.parameter_range();
    edge_surface_deviation(&boundary.curve, range, surface, 17).max_deviation <= tolerance
}
