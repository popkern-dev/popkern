//! Edge measurement functions.

use popkern_core::geometry::Vec3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::topology::{BRep, EdgeHandle};

use crate::MeasureError;

/// Compute the arc length of an edge.
///
/// Uses `CurveGeometry::length_between` on the edge's curve.
pub fn edge_length(brep: &BRep, edge: EdgeHandle) -> Result<f64, MeasureError> {
    let curve_on_edge = brep
        .edges
        .get(edge)
        .ok_or_else(|| MeasureError::InvalidHandle("edge".into()))?
        .curve
        .as_ref()
        .ok_or(MeasureError::MissingCurve)?;
    let len = curve_on_edge
        .geometry
        .length_between(curve_on_edge.t_start, curve_on_edge.t_end);
    if len < f64::EPSILON {
        return Err(MeasureError::DegenerateGeometry("zero-length edge".into()));
    }
    Ok(len)
}

/// Compute the tangent direction at the midpoint of an edge.
///
/// Returns a unit tangent vector.
pub fn edge_tangent(brep: &BRep, edge: EdgeHandle) -> Result<Vec3, MeasureError> {
    let curve_on_edge = brep
        .edges
        .get(edge)
        .ok_or_else(|| MeasureError::InvalidHandle("edge".into()))?
        .curve
        .as_ref()
        .ok_or(MeasureError::MissingCurve)?;
    let t_mid = (curve_on_edge.t_start + curve_on_edge.t_end) / 2.0;
    let tangent = curve_on_edge.geometry.tangent_at(t_mid);
    tangent
        .normalized()
        .ok_or_else(|| MeasureError::DegenerateGeometry("zero tangent at midpoint".into()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    fn make_cube() -> BRep {
        BRep::unit_cube(Tolerance::default()).expect("unit cube construction failed")
    }

    #[test]
    fn unit_cube_edge_length_is_one() {
        let cube = make_cube();
        let first_edge = cube.edges.keys().next().expect("cube has edges");
        let len = edge_length(&cube, first_edge).unwrap();
        assert!((len - 1.0).abs() < 1e-10, "expected 1.0, got {len}");
    }

    #[test]
    fn unit_cube_edge_tangent_is_unit() {
        let cube = make_cube();
        let first_edge = cube.edges.keys().next().expect("cube has edges");
        let t = edge_tangent(&cube, first_edge).unwrap();
        assert!(
            (t.length() - 1.0).abs() < 1e-10,
            "tangent should be unit length"
        );
    }

    #[test]
    fn invalid_handle_returns_error() {
        let cube = make_cube();
        let bad = EdgeHandle::default();
        assert!(edge_length(&cube, bad).is_err());
    }
}
