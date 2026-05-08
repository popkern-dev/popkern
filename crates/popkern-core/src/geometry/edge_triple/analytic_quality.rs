use crate::geometry::curve::CurveEval;
use crate::geometry::edge_triple::error::EdgeTripleError;
use crate::geometry::edge_triple::types::{EdgeTriple, EdgeTripleQuality};
use crate::geometry::surface::SurfaceGeometry;
use crate::tolerance::Tolerance;

pub(crate) fn finalize_quality_pub(
    triple: EdgeTriple,
    surf_a: &SurfaceGeometry,
    surf_b: &SurfaceGeometry,
    tol: &Tolerance,
) -> Result<EdgeTriple, EdgeTripleError> {
    finalize_quality(triple, surf_a, surf_b, tol)
}

pub(super) fn finalize_quality(
    mut triple: EdgeTriple,
    surf_a: &SurfaceGeometry,
    surf_b: &SurfaceGeometry,
    tol: &Tolerance,
) -> Result<EdgeTriple, EdgeTripleError> {
    let q = triple.verify(surf_a, surf_b, tol, 32).map_err(|_| {
        EdgeTripleError::DegenerateSurfacePair {
            reason: "post-construction verify failed",
        }
    })?;
    let endpoint_err = endpoint_err_3d(&triple);
    triple.quality = EdgeTripleQuality {
        max_endpoint_error: endpoint_err,
        ..q
    };
    Ok(triple)
}

fn endpoint_err_3d(triple: &EdgeTriple) -> f64 {
    use crate::geometry::curve::CurveGeometry;
    let p_start = triple.curve_3d.point_at(triple.t_start);
    let p_end = triple.curve_3d.point_at(triple.t_end);
    match &triple.curve_3d {
        CurveGeometry::Nurbs(n) => {
            if let (Some(first), Some(last)) = (n.control_points.first(), n.control_points.last()) {
                let d_start = first.distance_to(&p_start);
                let d_end = last.distance_to(&p_end);
                d_start.max(d_end)
            } else {
                f64::INFINITY
            }
        }
        _ => 0.0,
    }
}
