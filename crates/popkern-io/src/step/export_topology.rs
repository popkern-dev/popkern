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

        // 1b. Pre-export face surfaces so we have surface IDs for pcurves
        let surface_map = self.pre_export_surfaces(brep)?;

        // 2. Export all edges → EDGE_CURVE (with SURFACE_CURVE wrappers)
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

            // Export the 3-D curve geometry, or synthesize a line if absent
            let curve_3d_id = if let Some(curve) = &edge.curve {
                self.export_curve(curve)
            } else {
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

            // Build pcurve list from both half-edges
            let curve_id = self.maybe_wrap_surface_curve(
                brep,
                curve_3d_id,
                edge.half_edge,
                he.twin,
                &surface_map,
            );

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

        // 3. Export faces (surfaces already exported; reuses cached IDs)
        let mut face_ids: Vec<EntityId> = Vec::new();
        for (fh, face) in &brep.faces {
            let surface_id = surface_map.get(&fh).copied();
            let face_id = self.export_face(brep, fh, face, &vertex_map, &edge_map, surface_id)?;
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
    ///
    /// If `pre_surface_id` is `Some`, the surface was already exported during
    /// the pcurve pre-pass and we reuse that entity ID.
    fn export_face(
        &mut self,
        brep: &BRep,
        _fh: FaceHandle,
        face: &popkern_core::topology::Face,
        vertex_map: &HashMap<VertexHandle, EntityId>,
        edge_map: &HashMap<EdgeHandle, EntityId>,
        pre_surface_id: Option<EntityId>,
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

        // Use the pre-exported surface, or export now if none was cached
        let surface_id = if let Some(id) = pre_surface_id {
            id
        } else if let Some(surface) = &face.surface {
            self.export_surface(surface)
        } else {
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

    /// Pre-export face surfaces so their entity IDs are available for pcurves.
    fn pre_export_surfaces(
        &mut self,
        brep: &BRep,
    ) -> Result<HashMap<FaceHandle, EntityId>, IoError> {
        let mut surface_map: HashMap<FaceHandle, EntityId> = HashMap::new();
        for (fh, face) in &brep.faces {
            if let Some(surface) = &face.surface {
                let sid = self.export_surface(surface);
                surface_map.insert(fh, sid);
            }
        }
        Ok(surface_map)
    }

    /// If both half-edges have pcurves, wrap the 3-D curve in a SURFACE_CURVE.
    ///
    /// Returns either the original `curve_3d_id` (no pcurves) or a new
    /// SURFACE_CURVE entity ID.
    fn maybe_wrap_surface_curve(
        &mut self,
        brep: &BRep,
        curve_3d_id: EntityId,
        he_handle: popkern_core::topology::HalfEdgeHandle,
        twin_handle: popkern_core::topology::HalfEdgeHandle,
        surface_map: &HashMap<FaceHandle, EntityId>,
    ) -> EntityId {
        let mut pcurve_ids: Vec<EntityId> = Vec::new();

        // Try to export pcurve for the forward half-edge
        if let Some(pcurve_id) = self.try_export_he_pcurve(brep, he_handle, surface_map) {
            pcurve_ids.push(pcurve_id);
        }

        // Try to export pcurve for the twin half-edge
        if let Some(pcurve_id) = self.try_export_he_pcurve(brep, twin_handle, surface_map) {
            pcurve_ids.push(pcurve_id);
        }

        if pcurve_ids.is_empty() {
            curve_3d_id
        } else {
            self.export_surface_curve(curve_3d_id, &pcurve_ids)
        }
    }

    /// Try to export the pcurve for a single half-edge.
    ///
    /// Returns `None` if the half-edge has no pcurve or its face surface was
    /// not pre-exported.
    fn try_export_he_pcurve(
        &mut self,
        brep: &BRep,
        he_handle: popkern_core::topology::HalfEdgeHandle,
        surface_map: &HashMap<FaceHandle, EntityId>,
    ) -> Option<EntityId> {
        let he = brep.half_edges.get(he_handle)?;
        let pcurve = he.pcurve.as_ref()?;
        let loop_ = brep.loops.get(he.loop_)?;
        let surface_id = surface_map.get(&loop_.face)?;
        Some(self.export_pcurve(pcurve, *surface_id))
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
#[path = "export_topology_tests.rs"]
mod tests;
