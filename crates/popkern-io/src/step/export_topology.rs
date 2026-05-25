//! Export popkern BRep topology to STEP entity instances.
//!
//! Walks a [`BRep`] model and emits VERTEX_POINT, EDGE_CURVE, ORIENTED_EDGE,
//! EDGE_LOOP, FACE_BOUND, ADVANCED_FACE, CLOSED_SHELL, and
//! MANIFOLD_SOLID_BREP entities.

use std::collections::HashMap;

use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::geometry::surface::SurfaceOnFace;
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, LoopHandle, VertexHandle};

use crate::error::IoError;
use crate::step::export::{ExportContext, bool_, list, ref_, str_param};
use crate::step::export_pcurve_synthesis::{
    export_pcurve_trace_matches_edge, step_export_pcurve_association_tolerance,
    synthesize_export_half_edge_pcurve_on_surface,
};
use crate::step::export_topology_face_edges::pre_export_face_use_edges;
use crate::step::export_topology_face_uses::{EdgeFaceUse, edge_face_uses, unique_edge_half_edges};
use crate::step::export_topology_helpers::{
    collapsed_trim_missing_curve_edges, edge_curve_same_sense,
    edge_curve_surface_trim_disagreement, edge_end_vertex, pcurve_endpoints_match_edge_segment,
    synthesize_missing_edge_curve_from_surface_trims, synthesize_plane,
};
use crate::step::parser::types::EntityId;

#[derive(Clone)]
pub(super) struct ExportedFaceSurface {
    pub id: EntityId,
    pub surface: SurfaceOnFace,
}

struct TopologyExportMaps<'a> {
    edge: &'a HashMap<EdgeHandle, EntityId>,
    local_edge: &'a HashMap<EdgeFaceUse, EntityId>,
}

impl ExportContext {
    /// Export a complete [`BRep`] to STEP entities.
    ///
    /// Closed BReps become `MANIFOLD_SOLID_BREP`; BReps with explicit boundary
    /// edges become `SHELL_BASED_SURFACE_MODEL`.
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
        let collapsed_trim_edges = collapsed_trim_missing_curve_edges(brep);
        let mut local_edge_map: HashMap<EdgeFaceUse, EntityId> = HashMap::new();

        // 2. Export all edges → EDGE_CURVE (with SURFACE_CURVE wrappers)
        for (eh, edge) in &brep.edges {
            let he = brep
                .half_edges
                .get(edge.half_edge)
                .ok_or_else(|| IoError::Export("Edge references missing half-edge".into()))?;
            let start_vp = vertex_map
                .get(&he.origin)
                .ok_or_else(|| IoError::Export("Half-edge references unmapped vertex".into()))?;
            let end_vertex = edge_end_vertex(brep, eh, he)?;
            let end_vp = vertex_map
                .get(&end_vertex)
                .ok_or_else(|| IoError::Export("Edge references unmapped end vertex".into()))?;

            // Export the 3-D curve geometry. Valid collapsed surface trims get a
            // degenerate B-spline curve; other missing curves retain the legacy
            // line fallback for non-degenerate vertex-bounded edges.
            let surface_trim_curve = synthesize_missing_edge_curve_from_surface_trims(brep, eh);
            let curve_3d_id = if let Some(curve) = &edge.curve {
                let export_curve = surface_trim_curve
                    .as_ref()
                    .filter(|_| {
                        edge_curve_surface_trim_disagreement(brep, eh, curve)
                            .is_some_and(|gap| gap > step_export_pcurve_association_tolerance(brep))
                    })
                    .unwrap_or(curve);
                self.export_curve(export_curve)
            } else if let Some(curve) = surface_trim_curve {
                self.export_curve(&curve)
            } else if collapsed_trim_edges.contains(&eh) {
                let start_v = brep
                    .vertices
                    .get(he.origin)
                    .ok_or_else(|| IoError::Export("Vertex not found".into()))?;
                self.export_degenerate_curve_at(start_v.position)
            } else {
                let start_v = brep
                    .vertices
                    .get(he.origin)
                    .ok_or_else(|| IoError::Export("Vertex not found".into()))?;
                let end_v = brep
                    .vertices
                    .get(end_vertex)
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
            let curve_id = self.maybe_wrap_surface_curve(brep, curve_3d_id, eh, &surface_map);

            let same_sense = edge
                .curve
                .as_ref()
                .map(edge_curve_same_sense)
                .unwrap_or(true);
            let ec_id = self.emit(
                "EDGE_CURVE",
                vec![
                    str_param(""),
                    ref_(*start_vp),
                    ref_(*end_vp),
                    ref_(curve_id),
                    bool_(same_sense),
                ],
            );
            edge_map.insert(eh, ec_id);
        }
        pre_export_face_use_edges(
            self,
            brep,
            &vertex_map,
            &surface_map,
            &edge_map,
            &mut local_edge_map,
        )?;

        // 3. Export faces (surfaces already exported; reuses cached IDs)
        let mut face_ids: Vec<EntityId> = Vec::new();
        for (fh, face) in &brep.faces {
            let surface_id = surface_map.get(&fh).map(|surface| surface.id);
            let maps = TopologyExportMaps {
                edge: &edge_map,
                local_edge: &local_edge_map,
            };
            let face_id = self.export_face(brep, fh, face, &maps, surface_id)?;
            face_ids.push(face_id);
        }

        // 4. Emit shell item. Face-local STEP edge curves preserve geometric
        // fidelity for incoherent face uses, but they also make STEP topology
        // single-sided at those uses. Do not declare that representation as a
        // manifold solid.
        let exports_open_shell = brep.has_boundary_edges() || !local_edge_map.is_empty();
        let face_refs: Vec<_> = face_ids.iter().map(|&id| ref_(id)).collect();
        let shell_type = if exports_open_shell {
            "OPEN_SHELL"
        } else {
            "CLOSED_SHELL"
        };
        let shell_id = self.emit(shell_type, vec![str_param(""), list(face_refs)]);

        // 5. Emit the top-level representation item.
        let item_id = if exports_open_shell {
            self.emit(
                "SHELL_BASED_SURFACE_MODEL",
                vec![str_param(""), list(vec![ref_(shell_id)])],
            )
        } else {
            self.emit("MANIFOLD_SOLID_BREP", vec![str_param(""), ref_(shell_id)])
        };

        Ok(item_id)
    }

    /// Export a single face with its loops.
    ///
    /// If `pre_surface_id` is `Some`, the surface was already exported during
    /// the pcurve pre-pass and we reuse that entity ID.
    fn export_face(
        &mut self,
        brep: &BRep,
        face_handle: FaceHandle,
        face: &popkern_core::topology::Face,
        maps: &TopologyExportMaps<'_>,
        pre_surface_id: Option<EntityId>,
    ) -> Result<EntityId, IoError> {
        // Export outer loop
        let outer_bound_id = self.export_loop_as_bound(
            brep,
            face_handle,
            face.outer_loop,
            true, // is outer bound
            maps,
        )?;

        let mut bounds = vec![ref_(outer_bound_id)];

        // Export inner loops (holes)
        for &inner_lh in &face.inner_loops {
            let inner_bound_id =
                self.export_loop_as_bound(brep, face_handle, inner_lh, false, maps)?;
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
        face_handle: FaceHandle,
        lh: LoopHandle,
        is_outer: bool,
        maps: &TopologyExportMaps<'_>,
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

            let local_use = EdgeFaceUse {
                face: face_handle,
                half_edge: current_he,
            };
            let local_edge_curve_id = maps.local_edge.get(&local_use);
            let edge_curve_id = local_edge_curve_id
                .or_else(|| maps.edge.get(&he.edge))
                .ok_or_else(|| IoError::Export("Edge not in export map".into()))?;

            // Determine orientation: forward if this half-edge is the one
            // referenced by the parent edge (i.e., the "forward" half-edge)
            let edge = brep
                .edges
                .get(he.edge)
                .ok_or_else(|| IoError::Export("Edge not found".into()))?;
            let forward = local_edge_curve_id.is_some() || edge.half_edge == current_he;

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
    ) -> Result<HashMap<FaceHandle, ExportedFaceSurface>, IoError> {
        let mut surface_map: HashMap<FaceHandle, ExportedFaceSurface> = HashMap::new();
        for (fh, face) in &brep.faces {
            let surface = face
                .surface
                .clone()
                .map(Ok)
                .unwrap_or_else(|| synthesize_plane(brep, face.outer_loop))?;
            let id = self.export_surface(&surface);
            surface_map.insert(fh, ExportedFaceSurface { id, surface });
        }
        Ok(surface_map)
    }

    /// If both half-edges have pcurves, wrap the 3-D curve in a SURFACE_CURVE
    /// (or SEAM_CURVE if both half-edges belong to the same face).
    ///
    /// Returns either the original `curve_3d_id` (no pcurves) or a new
    /// SURFACE_CURVE / SEAM_CURVE entity ID.
    fn maybe_wrap_surface_curve(
        &mut self,
        brep: &BRep,
        curve_3d_id: EntityId,
        edge_handle: EdgeHandle,
        surface_map: &HashMap<FaceHandle, ExportedFaceSurface>,
    ) -> EntityId {
        let mut pcurve_ids: Vec<EntityId> = Vec::new();
        let face_uses = edge_face_uses(brep, edge_handle);
        for face_use in &face_uses {
            if let Some(pcurve_id) = self.try_export_he_pcurve(
                brep,
                edge_handle,
                face_use.half_edge,
                face_use.face,
                surface_map,
            ) {
                pcurve_ids.push(pcurve_id);
            }
        }
        if pcurve_ids.is_empty() {
            for he_handle in unique_edge_half_edges(brep, edge_handle) {
                if let Some(pcurve_id) = self.try_export_he_pcurve_via_stored_loop(
                    brep,
                    edge_handle,
                    he_handle,
                    surface_map,
                ) {
                    pcurve_ids.push(pcurve_id);
                }
            }
        }

        if pcurve_ids.is_empty() {
            curve_3d_id
        } else {
            // Detect seam edges: both half-edges belong to the same face
            let mut faces = face_uses.iter().map(|face_use| face_use.face);
            let first_face = faces.next();
            let is_seam = first_face.is_some() && faces.all(|face| Some(face) == first_face);

            if is_seam {
                self.export_seam_curve(curve_3d_id, &pcurve_ids)
            } else {
                self.export_surface_curve(curve_3d_id, &pcurve_ids)
            }
        }
    }

    /// Try to export the pcurve for a single half-edge.
    ///
    /// Returns `None` if the half-edge has no pcurve or its face surface was
    /// not pre-exported.
    fn try_export_he_pcurve(
        &mut self,
        brep: &BRep,
        edge_handle: EdgeHandle,
        he_handle: popkern_core::topology::HalfEdgeHandle,
        face_handle: FaceHandle,
        surface_map: &HashMap<FaceHandle, ExportedFaceSurface>,
    ) -> Option<EntityId> {
        let he = brep.half_edges.get(he_handle)?;
        let surface = surface_map.get(&face_handle)?;
        let pcurve = brep
            .face_use_pcurve(face_handle, he_handle)
            .or(he.pcurve.as_ref())
            .cloned()
            .filter(|pcurve| {
                pcurve_endpoints_match_edge_segment(
                    brep,
                    edge_handle,
                    pcurve,
                    &surface.surface.geometry,
                    step_export_pcurve_association_tolerance(brep),
                ) || export_pcurve_trace_matches_edge(
                    brep,
                    edge_handle,
                    he_handle,
                    pcurve,
                    &surface.surface.geometry,
                    step_export_pcurve_association_tolerance(brep),
                )
            })
            .or_else(|| {
                synthesize_export_half_edge_pcurve_on_surface(
                    brep,
                    edge_handle,
                    he_handle,
                    &surface.surface.geometry,
                )
            })?;
        if !pcurve_endpoints_match_edge_segment(
            brep,
            edge_handle,
            &pcurve,
            &surface.surface.geometry,
            step_export_pcurve_association_tolerance(brep),
        ) && !export_pcurve_trace_matches_edge(
            brep,
            edge_handle,
            he_handle,
            &pcurve,
            &surface.surface.geometry,
            step_export_pcurve_association_tolerance(brep),
        ) {
            return None;
        }
        // STEP doesn't have a native analytic-pcurve type — convert
        // any analytic Pcurve variant (Linear, etc.) to its NURBS
        // form via Pcurve::to_nurbs() before export. NURBS variants
        // are zero-cost.
        let nurbs = pcurve.to_nurbs();
        Some(self.export_pcurve(&nurbs, surface.id))
    }

    fn try_export_he_pcurve_via_stored_loop(
        &mut self,
        brep: &BRep,
        edge_handle: EdgeHandle,
        he_handle: popkern_core::topology::HalfEdgeHandle,
        surface_map: &HashMap<FaceHandle, ExportedFaceSurface>,
    ) -> Option<EntityId> {
        let he = brep.half_edges.get(he_handle)?;
        let loop_ = brep.loops.get(he.loop_)?;
        self.try_export_he_pcurve(brep, edge_handle, he_handle, loop_.face, surface_map)
    }

    fn export_degenerate_curve_at(&mut self, point: Point3) -> EntityId {
        self.export_nurbs_curve(&NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![point, point],
            weights: vec![1.0, 1.0],
        })
    }
}

use crate::step::parser::types::Parameter;

#[cfg(test)]
#[path = "export_topology_tests.rs"]
mod tests;
