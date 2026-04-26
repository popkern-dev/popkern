//! Reference geometry import from BRep face edges.
//!
//! Projects boundary edges of a BRep face onto the sketch plane,
//! creating construction geometry that can be used as reference
//! for constraining new sketch elements.

use std::collections::HashMap;

use popkern_constraint::{Constraint, ElementHandle, PointHandle};
use popkern_core::geometry::curve::CurveGeometry;
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

use crate::error::SketchError;
use crate::SketchOnPlane;

impl SketchOnPlane {
    /// Import edges from a BRep face as fixed construction geometry.
    ///
    /// Projects each boundary edge of the face onto the sketch plane,
    /// creating construction Lines/Arcs with Fixed constraints.
    /// Returns handles to the imported elements for constraining against.
    ///
    /// # Errors
    ///
    /// Returns [`SketchError::InvalidTopology`] if the face or its loop
    /// data is missing or malformed.
    pub fn import_face_edges(
        &mut self,
        brep: &BRep,
        face: FaceHandle,
    ) -> Result<Vec<ElementHandle>, SketchError> {
        let face_data = brep
            .faces
            .get(face)
            .ok_or_else(|| SketchError::InvalidTopology("face not found".into()))?;

        let outer_loop = brep
            .loops
            .get(face_data.outer_loop)
            .ok_or_else(|| SketchError::InvalidTopology("outer loop not found".into()))?;

        // Walk the half-edge chain around the outer loop.
        let start_he = outer_loop.half_edge;
        let mut half_edges = Vec::new();
        let mut cursor = start_he;
        loop {
            half_edges.push(cursor);
            cursor = brep.half_edges[cursor].next;
            if cursor == start_he {
                break;
            }
            if half_edges.len() > 10_000 {
                return Err(SketchError::InvalidTopology(
                    "half-edge chain did not close".into(),
                ));
            }
        }

        // Track already-projected vertices to avoid duplicate points.
        let mut vertex_map: HashMap<VertexHandle, PointHandle> = HashMap::new();
        let mut imported: Vec<ElementHandle> = Vec::new();

        for he_handle in &half_edges {
            let he = &brep.half_edges[*he_handle];
            let start_vertex = he.origin;
            let end_vertex = brep.half_edges[he.twin].origin;

            let edge = &brep.edges[he.edge];

            // Get or create 2D point for start vertex.
            let start_pt = get_or_create_point(
                &mut self.sketch,
                &self.plane,
                brep,
                start_vertex,
                &mut vertex_map,
                &mut imported,
            );

            // Get or create 2D point for end vertex.
            let end_pt = get_or_create_point(
                &mut self.sketch,
                &self.plane,
                brep,
                end_vertex,
                &mut vertex_map,
                &mut imported,
            );

            // Create the sketch element based on the edge curve type.
            match &edge.curve {
                Some(curve_on_edge) => match &curve_on_edge.geometry {
                    CurveGeometry::Line { .. } => {
                        let lh = self.sketch.add_line(start_pt, end_pt);
                        self.sketch.set_construction(ElementHandle::Line(lh), true);
                        imported.push(ElementHandle::Line(lh));
                    }
                    CurveGeometry::Arc {
                        center,
                        radius,
                        start_angle,
                        end_angle,
                        ..
                    } => {
                        let (cx, cy) = self.plane.to_2d(*center);
                        let ah = self
                            .sketch
                            .add_arc(cx, cy, *radius, *start_angle, *end_angle);
                        self.sketch.set_construction(ElementHandle::Arc(ah), true);
                        imported.push(ElementHandle::Arc(ah));
                    }
                    CurveGeometry::Circle { center, radius, .. } => {
                        let (cx, cy) = self.plane.to_2d(*center);
                        // Full circle: 0 to 2*PI.
                        let ah = self
                            .sketch
                            .add_arc(cx, cy, *radius, 0.0, std::f64::consts::TAU);
                        self.sketch.set_construction(ElementHandle::Arc(ah), true);
                        imported.push(ElementHandle::Arc(ah));
                    }
                    CurveGeometry::Ellipse { .. }
                    | CurveGeometry::Parabola { .. }
                    | CurveGeometry::Hyperbola { .. }
                    | CurveGeometry::Nurbs(_) => {
                        // Fall back to a line between the projected endpoints.
                        // Sketches don't have a native ellipse element today.
                        let lh = self.sketch.add_line(start_pt, end_pt);
                        self.sketch.set_construction(ElementHandle::Line(lh), true);
                        imported.push(ElementHandle::Line(lh));
                    }
                },
                None => {
                    // No curve geometry — create a line between endpoints.
                    let lh = self.sketch.add_line(start_pt, end_pt);
                    self.sketch.set_construction(ElementHandle::Line(lh), true);
                    imported.push(ElementHandle::Line(lh));
                }
            }
        }

        Ok(imported)
    }
}

/// Get or create a 2D sketch point for a BRep vertex.
///
/// If the vertex has already been projected, returns the existing
/// [`PointHandle`]. Otherwise projects it to 2D, creates a new point
/// with a `Fixed` constraint, marks it as construction, and records it
/// in the map.
fn get_or_create_point(
    sketch: &mut popkern_constraint::Sketch,
    plane: &crate::plane::SketchPlane,
    brep: &BRep,
    vertex: VertexHandle,
    vertex_map: &mut HashMap<VertexHandle, PointHandle>,
    imported: &mut Vec<ElementHandle>,
) -> PointHandle {
    if let Some(&ph) = vertex_map.get(&vertex) {
        return ph;
    }

    let pos_3d = brep.vertices[vertex].position;
    let (x, y) = plane.to_2d(pos_3d);
    let ph = sketch.add_point(x, y);
    sketch.set_construction(ElementHandle::Point(ph), true);
    sketch.add_constraint(Constraint::Fixed(ph));
    imported.push(ElementHandle::Point(ph));
    vertex_map.insert(vertex, ph);
    ph
}

#[cfg(test)]
#[path = "reference_import_tests.rs"]
mod tests;
