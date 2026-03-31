//! Angle measurement functions.

use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

use crate::MeasureError;

/// Compute the dihedral angle between two faces.
///
/// Returns the angle in radians between the outward normals of the two faces.
/// Always in [0, pi].
pub fn dihedral_angle(
    brep: &BRep,
    face_a: FaceHandle,
    face_b: FaceHandle,
) -> Result<f64, MeasureError> {
    let na = crate::face::face_normal(brep, face_a)?;
    let nb = crate::face::face_normal(brep, face_b)?;
    Ok(na.angle_to(&nb))
}

/// Compute the angle between two edges.
///
/// Returns the angle in radians between the tangent directions at their
/// midpoints. Always in [0, pi].
pub fn edge_angle(
    brep: &BRep,
    edge_a: EdgeHandle,
    edge_b: EdgeHandle,
) -> Result<f64, MeasureError> {
    let ta = crate::edge::edge_tangent(brep, edge_a)?;
    let tb = crate::edge::edge_tangent(brep, edge_b)?;
    Ok(ta.angle_to(&tb))
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    fn make_cube() -> BRep {
        BRep::unit_cube(Tolerance::default()).expect("cube")
    }

    #[test]
    fn cube_adjacent_faces_are_perpendicular() {
        let cube = make_cube();
        let faces: Vec<_> = cube.faces.keys().collect();
        assert!(faces.len() >= 2);

        // Find two faces that share an edge (adjacent) — should be 90 degrees
        let mut found = false;
        for i in 0..faces.len() {
            for j in (i + 1)..faces.len() {
                if let Ok(angle) = dihedral_angle(&cube, faces[i], faces[j]) {
                    let angle_deg = angle.to_degrees();
                    if (angle_deg - 90.0).abs() < 1.0 {
                        found = true;
                        break;
                    }
                }
            }
            if found {
                break;
            }
        }
        assert!(found, "should find two perpendicular faces on a cube");
    }

    #[test]
    fn cube_edge_angle_between_parallel_edges() {
        let cube = make_cube();
        let edges: Vec<_> = cube.edges.keys().collect();
        // Any two edges on a cube are either parallel (0 or pi) or perpendicular (pi/2)
        let angle = edge_angle(&cube, edges[0], edges[1]).unwrap();
        let deg = angle.to_degrees();
        assert!(
            (deg - 0.0).abs() < 1.0 || (deg - 90.0).abs() < 1.0 || (deg - 180.0).abs() < 1.0,
            "cube edges should be 0, 90, or 180 degrees apart, got {deg}"
        );
    }
}
