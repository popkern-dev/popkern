//! Export popkern BRep topology to STEP entity instances.
//!
//! Walks a [`BRep`] model and emits VERTEX_POINT, EDGE_CURVE, ORIENTED_EDGE,
//! EDGE_LOOP, FACE_BOUND, ADVANCED_FACE, CLOSED_SHELL, and
//! MANIFOLD_SOLID_BREP entities.

use std::collections::HashMap;

use popkern_core::geometry::Vec3;
use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, LoopHandle, VertexHandle};

use crate::error::IoError;
use crate::step::export::{ExportContext, bool_, list, ref_, str_param};
use crate::step::parser::types::EntityId;

impl ExportContext {
    /// Export a complete [`BRep`] solid to STEP entities.
    ///
    /// Returns the entity ID of the MANIFOLD_SOLID_BREP.
    pub(crate) fn export_brep(&mut self, brep: &BRep) -> Result<EntityId, IoError> {
        // Maps from popkern handles to STEP entity IDs
        let mut vertex_map: HashMap<VertexHandle, EntityId> = HashMap::new();
        let mut edge_map: HashMap<EdgeHandle, EntityId> = HashMap::new();

        // 1. Export all vertices → VERTEX_POINT
        for (vh, vertex) in &brep.vertices {
            let pt_id = self.export_point(&vertex.position);
            let vp_id = self.emit("VERTEX_POINT", vec![str_param(""), ref_(pt_id)]);
            vertex_map.insert(vh, vp_id);
        }

        // 2. Export all edges → EDGE_CURVE
        for (eh, edge) in &brep.edges {
            let he = brep
                .half_edges
                .get(edge.half_edge)
                .ok_or_else(|| IoError::Export("Edge references missing half-edge".into()))?;
            let twin = brep
                .half_edges
                .get(he.twin)
                .ok_or_else(|| IoError::Export("Half-edge references missing twin".into()))?;

            let start_vp = vertex_map
                .get(&he.origin)
                .ok_or_else(|| IoError::Export("Half-edge references unmapped vertex".into()))?;
            let end_vp = vertex_map.get(&twin.origin).ok_or_else(|| {
                IoError::Export("Twin half-edge references unmapped vertex".into())
            })?;

            // Export the curve geometry, or create a default line if absent
            let curve_id = if let Some(curve) = &edge.curve {
                self.export_curve(curve)
            } else {
                // Synthesize a line from start to end vertex
                let start_v = brep
                    .vertices
                    .get(he.origin)
                    .ok_or_else(|| IoError::Export("Vertex not found".into()))?;
                let end_v = brep
                    .vertices
                    .get(twin.origin)
                    .ok_or_else(|| IoError::Export("Vertex not found".into()))?;
                let dir = end_v.position - start_v.position;
                let fallback = CurveOnEdge {
                    geometry: CurveGeometry::Line {
                        origin: start_v.position,
                        direction: dir,
                    },
                    t_start: 0.0,
                    t_end: 1.0,
                };
                self.export_curve(&fallback)
            };

            let ec_id = self.emit(
                "EDGE_CURVE",
                vec![
                    str_param(""),
                    ref_(*start_vp),
                    ref_(*end_vp),
                    ref_(curve_id),
                    bool_(true), // same_sense
                ],
            );
            edge_map.insert(eh, ec_id);
        }

        // 3. Export faces
        let mut face_ids: Vec<EntityId> = Vec::new();
        for (fh, face) in &brep.faces {
            let face_id = self.export_face(brep, fh, face, &vertex_map, &edge_map)?;
            face_ids.push(face_id);
        }

        // 4. Emit CLOSED_SHELL
        let face_refs: Vec<_> = face_ids.iter().map(|&id| ref_(id)).collect();
        let shell_id = self.emit("CLOSED_SHELL", vec![str_param(""), list(face_refs)]);

        // 5. Emit MANIFOLD_SOLID_BREP
        let solid_id = self.emit("MANIFOLD_SOLID_BREP", vec![str_param(""), ref_(shell_id)]);

        Ok(solid_id)
    }

    /// Export a single face with its loops.
    fn export_face(
        &mut self,
        brep: &BRep,
        _fh: FaceHandle,
        face: &popkern_core::topology::Face,
        vertex_map: &HashMap<VertexHandle, EntityId>,
        edge_map: &HashMap<EdgeHandle, EntityId>,
    ) -> Result<EntityId, IoError> {
        // Export outer loop
        let outer_bound_id = self.export_loop_as_bound(
            brep,
            face.outer_loop,
            true, // is outer bound
            vertex_map,
            edge_map,
        )?;

        let mut bounds = vec![ref_(outer_bound_id)];

        // Export inner loops (holes)
        for &inner_lh in &face.inner_loops {
            let inner_bound_id =
                self.export_loop_as_bound(brep, inner_lh, false, vertex_map, edge_map)?;
            bounds.push(ref_(inner_bound_id));
        }

        // Export surface geometry
        let surface_id = if let Some(surface) = &face.surface {
            self.export_surface(surface)
        } else {
            // Synthesize a default plane from face vertices
            let default_surface = synthesize_plane(brep, face.outer_loop)?;
            self.export_surface(&default_surface)
        };

        // Determine face orientation
        let orientation = face.surface.as_ref().map(|s| s.orientation).unwrap_or(true);

        let face_id = self.emit(
            "ADVANCED_FACE",
            vec![
                str_param(""),
                list(bounds),
                ref_(surface_id),
                bool_(orientation),
            ],
        );

        Ok(face_id)
    }

    /// Export a loop as FACE_OUTER_BOUND or FACE_BOUND.
    fn export_loop_as_bound(
        &mut self,
        brep: &BRep,
        lh: LoopHandle,
        is_outer: bool,
        _vertex_map: &HashMap<VertexHandle, EntityId>,
        edge_map: &HashMap<EdgeHandle, EntityId>,
    ) -> Result<EntityId, IoError> {
        let loop_ = brep
            .loops
            .get(lh)
            .ok_or_else(|| IoError::Export("Loop not found".into()))?;

        // Walk the half-edge loop to collect oriented edges
        let mut oriented_edges = Vec::new();
        let start_he = loop_.half_edge;
        let mut current_he = start_he;

        loop {
            let he = brep
                .half_edges
                .get(current_he)
                .ok_or_else(|| IoError::Export("Half-edge not found during loop walk".into()))?;

            let edge_curve_id = edge_map
                .get(&he.edge)
                .ok_or_else(|| IoError::Export("Edge not in export map".into()))?;

            // Determine orientation: forward if this half-edge is the one
            // referenced by the parent edge (i.e., the "forward" half-edge)
            let edge = brep
                .edges
                .get(he.edge)
                .ok_or_else(|| IoError::Export("Edge not found".into()))?;
            let forward = edge.half_edge == current_he;

            let oe_id = self.emit(
                "ORIENTED_EDGE",
                vec![
                    str_param(""),
                    Parameter::Derived, // edge_start
                    Parameter::Derived, // edge_end
                    ref_(*edge_curve_id),
                    bool_(forward),
                ],
            );
            oriented_edges.push(ref_(oe_id));

            current_he = he.next;
            if current_he == start_he {
                break;
            }
        }

        // Emit EDGE_LOOP
        let el_id = self.emit("EDGE_LOOP", vec![str_param(""), list(oriented_edges)]);

        // Emit FACE_OUTER_BOUND or FACE_BOUND
        let bound_type = if is_outer {
            "FACE_OUTER_BOUND"
        } else {
            "FACE_BOUND"
        };
        let bound_id = self.emit(bound_type, vec![str_param(""), ref_(el_id), bool_(true)]);

        Ok(bound_id)
    }
}

use crate::step::parser::types::Parameter;

/// Synthesize a plane surface from the first three vertices of a loop.
fn synthesize_plane(brep: &BRep, lh: LoopHandle) -> Result<SurfaceOnFace, IoError> {
    let loop_ = brep
        .loops
        .get(lh)
        .ok_or_else(|| IoError::Export("Loop not found for plane synthesis".into()))?;

    // Collect up to 3 vertex positions from the loop
    let mut positions = Vec::new();
    let start_he = loop_.half_edge;
    let mut current = start_he;
    loop {
        let he = brep
            .half_edges
            .get(current)
            .ok_or_else(|| IoError::Export("Half-edge not found during plane synthesis".into()))?;
        if let Some(v) = brep.vertices.get(he.origin) {
            positions.push(v.position);
        }
        if positions.len() >= 3 {
            break;
        }
        current = he.next;
        if current == start_he {
            break;
        }
    }

    if positions.len() < 3 {
        return Err(IoError::Export(
            "Cannot synthesize plane: loop has fewer than 3 vertices".into(),
        ));
    }

    let origin = positions[0];
    let u_vec = positions[1] - positions[0];
    let v_vec = positions[2] - positions[0];
    let normal = u_vec.cross(&v_vec);
    let normal = normal.normalized().unwrap_or(Vec3::new(0.0, 0.0, 1.0));
    let u_axis = u_vec.normalized().unwrap_or(Vec3::new(1.0, 0.0, 0.0));
    let v_axis = normal.cross(&u_axis);

    Ok(SurfaceOnFace {
        geometry: SurfaceGeometry::Plane {
            origin,
            normal,
            u_axis,
            v_axis,
        },
        orientation: true,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::tolerance::Tolerance;

    #[test]
    fn export_unit_cube_entity_counts() {
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let mut ctx = ExportContext::new();
        let solid_id = ctx.export_brep(&brep).unwrap();

        // Count entity types
        let count =
            |name: &str| -> usize { ctx.entities.iter().filter(|e| e.type_name == name).count() };

        // A cube has 8 vertices, 12 edges, 6 faces
        assert_eq!(count("VERTEX_POINT"), 8, "Expected 8 VERTEX_POINTs");
        assert_eq!(count("EDGE_CURVE"), 12, "Expected 12 EDGE_CURVEs");
        assert_eq!(count("ADVANCED_FACE"), 6, "Expected 6 ADVANCED_FACEs");
        assert_eq!(count("CLOSED_SHELL"), 1, "Expected 1 CLOSED_SHELL");
        assert_eq!(
            count("MANIFOLD_SOLID_BREP"),
            1,
            "Expected 1 MANIFOLD_SOLID_BREP"
        );

        // Each face has one outer loop
        assert_eq!(count("FACE_OUTER_BOUND"), 6, "Expected 6 FACE_OUTER_BOUNDs");
        assert_eq!(count("EDGE_LOOP"), 6, "Expected 6 EDGE_LOOPs");

        // Each edge appears in two faces → 24 oriented edges
        assert_eq!(count("ORIENTED_EDGE"), 24, "Expected 24 ORIENTED_EDGEs");

        // The solid ID should refer to a real entity
        let solid = ctx.entities.iter().find(|e| e.id == solid_id);
        assert!(solid.is_some());
        assert_eq!(solid.unwrap().type_name, "MANIFOLD_SOLID_BREP");
    }

    #[test]
    fn export_cube_surfaces_are_planes() {
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let mut ctx = ExportContext::new();
        ctx.export_brep(&brep).unwrap();

        let plane_count = ctx
            .entities
            .iter()
            .filter(|e| e.type_name == "PLANE")
            .count();
        // A cube has 6 planar faces; some planes may be shared if
        // the deduplication catches them, but typically 6 distinct planes
        assert!(
            plane_count >= 3 && plane_count <= 6,
            "Expected 3-6 PLANEs for a cube, got {plane_count}"
        );
    }

    #[test]
    fn export_brep_vertex_points_reference_cartesian_points() {
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let mut ctx = ExportContext::new();
        ctx.export_brep(&brep).unwrap();

        for entity in &ctx.entities {
            if entity.type_name == "VERTEX_POINT" {
                // VERTEX_POINT('', #ref_to_cartesian_point)
                assert_eq!(entity.params.len(), 2);
                assert!(matches!(entity.params[1], Parameter::EntityRef(_)));
            }
        }
    }

    #[test]
    fn export_brep_with_inner_loop() {
        use popkern_core::constructors::builder::BRepBuilder;
        use popkern_core::geometry::Point3;

        let mut b = BRepBuilder::new(Tolerance::default());

        // Outer quad: v0-v1-v2-v3
        let v0 = b.add_vertex(Point3::new(0.0, 0.0, 0.0));
        let v1 = b.add_vertex(Point3::new(4.0, 0.0, 0.0));
        let v2 = b.add_vertex(Point3::new(4.0, 4.0, 0.0));
        let v3 = b.add_vertex(Point3::new(0.0, 4.0, 0.0));

        // Inner quad (hole): v4-v5-v6-v7
        let v4 = b.add_vertex(Point3::new(1.0, 1.0, 0.0));
        let v5 = b.add_vertex(Point3::new(3.0, 1.0, 0.0));
        let v6 = b.add_vertex(Point3::new(3.0, 3.0, 0.0));
        let v7 = b.add_vertex(Point3::new(1.0, 3.0, 0.0));

        // Outer loop edges
        let (_, he01, _) = b.add_edge(v0, v1);
        let (_, he12, _) = b.add_edge(v1, v2);
        let (_, he23, _) = b.add_edge(v2, v3);
        let (_, he30, _) = b.add_edge(v3, v0);

        // Inner loop edges (CW from outside = hole)
        let (_, he47, _) = b.add_edge(v4, v7);
        let (_, he76, _) = b.add_edge(v7, v6);
        let (_, he65, _) = b.add_edge(v6, v5);
        let (_, he54, _) = b.add_edge(v5, v4);

        let solid = b.add_solid();
        let shell = b.add_shell(solid);
        let face = b.add_face(shell);
        b.add_loop(&[he01, he12, he23, he30], face).unwrap();
        b.add_inner_loop(&[he47, he76, he65, he54], face).unwrap();

        let brep = b.build_unchecked();
        let mut ctx = ExportContext::new();
        let solid_id = ctx.export_brep(&brep).unwrap();

        let count =
            |name: &str| -> usize { ctx.entities.iter().filter(|e| e.type_name == name).count() };

        // Should have 1 FACE_OUTER_BOUND and 1 FACE_BOUND
        assert_eq!(count("FACE_OUTER_BOUND"), 1);
        assert_eq!(count("FACE_BOUND"), 1);
        // 2 edge loops (outer + inner)
        assert_eq!(count("EDGE_LOOP"), 2);
        assert_eq!(count("ADVANCED_FACE"), 1);

        let solid = ctx.entities.iter().find(|e| e.id == solid_id);
        assert!(solid.is_some());
    }

    #[test]
    fn export_brep_edges_without_curves_synthesize_lines() {
        // unit_cube doesn't set curves, so the exporter should synthesize them
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let mut ctx = ExportContext::new();
        ctx.export_brep(&brep).unwrap();

        let line_count = ctx
            .entities
            .iter()
            .filter(|e| e.type_name == "LINE")
            .count();
        // All 12 edges should have synthesized LINE curves
        assert!(
            line_count >= 12,
            "Expected at least 12 LINEs, got {line_count}"
        );
    }

    #[test]
    fn export_oriented_edge_has_correct_structure() {
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let mut ctx = ExportContext::new();
        ctx.export_brep(&brep).unwrap();

        for entity in &ctx.entities {
            if entity.type_name == "ORIENTED_EDGE" {
                // ORIENTED_EDGE('', *, *, #edge_curve, .T./.F.)
                assert_eq!(entity.params.len(), 5);
                assert_eq!(entity.params[0], Parameter::String("".into()));
                assert_eq!(entity.params[1], Parameter::Derived);
                assert_eq!(entity.params[2], Parameter::Derived);
                assert!(matches!(entity.params[3], Parameter::EntityRef(_)));
                assert!(
                    entity.params[4] == Parameter::Enum("T".into())
                        || entity.params[4] == Parameter::Enum("F".into())
                );
            }
        }
    }
}
