//! Vertex measurement functions.

use popkern_core::geometry::Point3;
use popkern_core::topology::{BRep, VertexHandle};

use crate::MeasureError;

/// Return the position of a vertex.
pub fn vertex_position(brep: &BRep, vertex: VertexHandle) -> Result<Point3, MeasureError> {
    let v = brep
        .vertices
        .get(vertex)
        .ok_or_else(|| MeasureError::InvalidHandle("vertex".into()))?;
    Ok(v.position)
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::tolerance::Tolerance;

    fn make_cube() -> BRep {
        BRep::unit_cube(Tolerance::default()).expect("unit cube construction failed")
    }

    #[test]
    fn cube_vertex_positions_within_bounds() {
        let cube = make_cube();
        for (vh, _) in &cube.vertices {
            let pos = vertex_position(&cube, vh).unwrap();
            assert!(pos.x >= -1e-10 && pos.x <= 1.0 + 1e-10);
            assert!(pos.y >= -1e-10 && pos.y <= 1.0 + 1e-10);
            assert!(pos.z >= -1e-10 && pos.z <= 1.0 + 1e-10);
        }
    }

    #[test]
    fn invalid_vertex_handle() {
        let cube = make_cube();
        assert!(vertex_position(&cube, VertexHandle::default()).is_err());
    }
}
