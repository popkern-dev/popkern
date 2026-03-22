//! Geometric and topological validation utilities.

use crate::error::CoreError;
use crate::topology_v0::BRepModel;
use std::collections::HashMap;

/// Compute the Euler characteristic: V - E + F.
pub fn euler_characteristic(model: &BRepModel) -> i64 {
    let v = model.vertices.len() as i64;
    let e = model.edges.len() as i64;
    let f = model.faces.len() as i64;
    v - e + f
}

/// Check the Euler characteristic matches the expected value for the given genus.
///
/// For a closed solid of genus g: V - E + F = 2 - 2g.
/// A genus-0 solid (cube, sphere) should have Euler characteristic 2.
pub fn check_euler(model: &BRepModel, genus: u32) -> Result<(), CoreError> {
    let expected = 2 - 2 * genus as i64;
    let actual = euler_characteristic(model);
    if actual != expected {
        return Err(CoreError::ValidationError(format!(
            "Euler characteristic: expected {} (genus {}), got {}. \
             Check V={}, E={}, F={}. Ensure no duplicate or missing entities.",
            expected,
            genus,
            actual,
            model.vertices.len(),
            model.edges.len(),
            model.faces.len()
        )));
    }
    Ok(())
}

/// Check that every edge is shared by exactly two faces (manifold condition).
pub fn check_manifold_edges(model: &BRepModel) -> Result<(), CoreError> {
    let mut edge_face_count: HashMap<usize, usize> = HashMap::new();
    for face in &model.faces {
        for &edge_id in &face.edge_ids {
            *edge_face_count.entry(edge_id).or_insert(0) += 1;
        }
    }
    for (edge_id, count) in &edge_face_count {
        if *count != 2 {
            return Err(CoreError::ValidationError(format!(
                "Edge {} is shared by {} faces (expected exactly 2). \
                 Non-manifold edge. Check face definitions for missing or duplicate edge references.",
                edge_id, count
            )));
        }
    }
    Ok(())
}

/// Check that every vertex bounds at least two edges.
pub fn check_vertex_valence(model: &BRepModel) -> Result<(), CoreError> {
    let mut vertex_edge_count: HashMap<usize, usize> = HashMap::new();
    for edge in &model.edges {
        *vertex_edge_count.entry(edge.start_vertex).or_insert(0) += 1;
        *vertex_edge_count.entry(edge.end_vertex).or_insert(0) += 1;
    }
    for (vertex_id, count) in &vertex_edge_count {
        if *count < 2 {
            return Err(CoreError::ValidationError(format!(
                "Vertex {} has valence {} (expected >= 2). \
                 Isolated or dangling vertex. Remove it or connect it to at least 2 edges.",
                vertex_id, count
            )));
        }
    }
    Ok(())
}

// TODO: Accept genus as parameter for non-genus-0 models (e.g., torus).
/// Run all validation checks on a BRep model.
///
/// The `tolerance` parameter is reserved for future geometric checks
/// (loop closure, curve-surface consistency). Currently used to establish
/// the explicit-tolerance pattern from day one.
pub fn validate_model(model: &BRepModel, _tolerance: f64) -> Result<(), CoreError> {
    check_euler(model, 0)?;
    check_manifold_edges(model)?;
    check_vertex_valence(model)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::Point3;
    use crate::topology_v0::*;

    fn make_unit_cube() -> BRepModel {
        let vertices = vec![
            Vertex::new(0, Point3::new(0.0, 0.0, 0.0)),
            Vertex::new(1, Point3::new(1.0, 0.0, 0.0)),
            Vertex::new(2, Point3::new(1.0, 1.0, 0.0)),
            Vertex::new(3, Point3::new(0.0, 1.0, 0.0)),
            Vertex::new(4, Point3::new(0.0, 0.0, 1.0)),
            Vertex::new(5, Point3::new(1.0, 0.0, 1.0)),
            Vertex::new(6, Point3::new(1.0, 1.0, 1.0)),
            Vertex::new(7, Point3::new(0.0, 1.0, 1.0)),
        ];
        let edges = vec![
            Edge::new(0, 0, 1),
            Edge::new(1, 1, 2),
            Edge::new(2, 2, 3),
            Edge::new(3, 3, 0),
            Edge::new(4, 4, 5),
            Edge::new(5, 5, 6),
            Edge::new(6, 6, 7),
            Edge::new(7, 7, 4),
            Edge::new(8, 0, 4),
            Edge::new(9, 1, 5),
            Edge::new(10, 2, 6),
            Edge::new(11, 3, 7),
        ];
        let faces = vec![
            Face::new(0, vec![0, 1, 2, 3]),
            Face::new(1, vec![4, 5, 6, 7]),
            Face::new(2, vec![0, 9, 4, 8]),
            Face::new(3, vec![1, 10, 5, 9]),
            Face::new(4, vec![2, 11, 6, 10]),
            Face::new(5, vec![3, 8, 7, 11]),
        ];
        let shell = Shell::new(0, vec![0, 1, 2, 3, 4, 5]);
        let solid = Solid::new(0, 0);
        BRepModel {
            vertices,
            edges,
            faces,
            shells: vec![shell],
            solids: vec![solid],
        }
    }

    #[test]
    fn euler_characteristic_cube_is_2() {
        let model = make_unit_cube();
        assert_eq!(euler_characteristic(&model), 2);
    }

    #[test]
    fn check_euler_passes_for_cube() {
        let model = make_unit_cube();
        assert!(check_euler(&model, 0).is_ok());
    }

    #[test]
    fn check_manifold_edges_passes_for_cube() {
        let model = make_unit_cube();
        assert!(check_manifold_edges(&model).is_ok());
    }

    #[test]
    fn check_vertex_valence_passes_for_cube() {
        let model = make_unit_cube();
        assert!(check_vertex_valence(&model).is_ok());
    }

    #[test]
    fn validate_model_passes_for_cube() {
        let model = make_unit_cube();
        assert!(validate_model(&model, 1e-6).is_ok());
    }

    #[test]
    fn check_euler_fails_for_wrong_genus() {
        let model = make_unit_cube();
        // Cube is genus 0 (Euler = 2), asking for genus 1 (Euler = 0) should fail
        let err = check_euler(&model, 1).unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("Euler characteristic"));
        assert!(msg.contains("expected 0"));
        assert!(msg.contains("got 2"));
    }

    #[test]
    fn check_manifold_edges_fails_for_non_manifold() {
        let mut model = make_unit_cube();
        // Add a third face referencing edge 0, making it non-manifold (3 faces share edge 0)
        model.faces.push(Face::new(6, vec![0]));
        let err = check_manifold_edges(&model).unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("Non-manifold edge"));
    }

    #[test]
    fn check_vertex_valence_fails_for_dangling_vertex() {
        let mut model = make_unit_cube();
        // Add two new vertices with only one edge between them (valence 1 each)
        model
            .vertices
            .push(Vertex::new(98, Point3::new(5.0, 5.0, 5.0)));
        model
            .vertices
            .push(Vertex::new(99, Point3::new(6.0, 6.0, 6.0)));
        model.edges.push(Edge::new(12, 98, 99));
        let err = check_vertex_valence(&model).unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("valence"));
    }
}
