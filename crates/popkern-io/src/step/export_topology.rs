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
use crate::step::export_pcurve_synthesis::step_export_pcurve_association_tolerance;
use crate::step::export_topology_face_edges::pre_export_face_use_edges;
use crate::step::export_topology_face_uses::EdgeFaceUse;
use crate::step::export_topology_helpers::{
    collapsed_trim_missing_curve_edges, edge_curve_surface_trim_disagreement,
    edge_curve_surface_trim_geometric_gap, edge_end_vertex,
    synthesize_missing_edge_curve_from_surface_trims, synthesize_plane,
};
use crate::step::parser::types::EntityId;

#[derive(Clone)]
pub(crate) struct ExportedFaceSurface {
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
    /// Export a [`BRep`] to one representation item PER SHELL.
    ///
    /// A multi-shell body (a disjoint union: a bearing's races and
    /// balls in one document export) is NOT one manifold solid — one
    /// `CLOSED_SHELL` swallowing every disjoint shell fails
    /// `BRepCheck` downstream. Each closed shell becomes its own
    /// `MANIFOLD_SOLID_BREP`. Single-shell bodies and open-shell
    /// (sheet) bodies keep their historical single-item forms.
    pub(crate) fn export_brep_items(&mut self, brep: &BRep) -> Result<Vec<EntityId>, IoError> {
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
        // Periodic UV-loop branch normalization: per-face-use shifts
        // that chain each pcurve into its face's closed trim contour
        // (seam pcurve pairs become branch-separated; rim windows
        // align). See export_periodic_branches module docs.
        let branch_shifts =
            super::export_periodic_branches::compute_periodic_branch_shifts(brep, &surface_map);
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
            // degenerate B-spline curve; other missing curves retain the earlier
            // line fallback for non-degenerate vertex-bounded edges.
            let surface_trim_curve = synthesize_missing_edge_curve_from_surface_trims(brep, eh);
            let curve_3d_id = if let Some(curve) = &edge.curve {
                let export_curve = surface_trim_curve
                    .as_ref()
                    .filter(|_| {
                        let gap = edge_curve_surface_trim_disagreement(brep, eh, curve);
                        if popkern_core::trace::trace_enabled("step-io") {
                            eprintln!(
                                "  [step-export] edge {eh:?} trim disagreement {gap:?} (tol {})",
                                step_export_pcurve_association_tolerance(brep)
                            );
                        }
                        // Replace the stored curve only when it is off its
                        // trims as a point set; a curve merely parameterized
                        // differently is kept exact and its pcurves are
                        // rebuilt on it (PK-554).
                        let tolerance = step_export_pcurve_association_tolerance(brep);
                        gap.is_some_and(|g| g > tolerance)
                            && edge_curve_surface_trim_geometric_gap(brep, eh, curve)
                                .is_none_or(|g| g > tolerance)
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
            let curve_id =
                self.maybe_wrap_surface_curve(brep, curve_3d_id, eh, &surface_map, &branch_shifts);

            let same_sense = super::export_edge_sense::edge_same_sense(brep, eh, &branch_shifts);
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

        // 3. Export faces (surfaces already exported; reuses cached
        // IDs), grouped by their owning brep SHELL.
        let mut face_ids_by_shell: Vec<Vec<EntityId>> = Vec::new();
        let mut shell_handles: Vec<_> = brep.shells.keys().collect();
        shell_handles.sort();
        for &sh in &shell_handles {
            let mut ids = Vec::new();
            for (fh, face) in &brep.faces {
                if face.shell != sh {
                    continue;
                }
                let surface_id = surface_map.get(&fh).map(|surface| surface.id);
                let maps = TopologyExportMaps {
                    edge: &edge_map,
                    local_edge: &local_edge_map,
                };
                let face_id = self.export_face(brep, fh, face, &maps, surface_id)?;
                ids.push(face_id);
            }
            if !ids.is_empty() {
                face_ids_by_shell.push(ids);
            }
        }
        // Faces with no live shell (defensive) travel together.
        {
            let mut orphans = Vec::new();
            for (fh, face) in &brep.faces {
                if brep.shells.contains_key(face.shell) {
                    continue;
                }
                let surface_id = surface_map.get(&fh).map(|surface| surface.id);
                let maps = TopologyExportMaps {
                    edge: &edge_map,
                    local_edge: &local_edge_map,
                };
                orphans.push(self.export_face(brep, fh, face, &maps, surface_id)?);
            }
            if !orphans.is_empty() {
                face_ids_by_shell.push(orphans);
            }
        }

        // 4. Emit shell items. Face-local STEP edge curves preserve
        // geometric fidelity for incoherent face uses, but they also
        // make STEP topology single-sided at those uses. Do not
        // declare that representation as a manifold solid.
        let exports_open_shell = brep.has_boundary_edges() || !local_edge_map.is_empty();
        if exports_open_shell {
            // Historical sheet form: one SHELL_BASED_SURFACE_MODEL
            // over every shell.
            let shell_refs: Vec<_> = face_ids_by_shell
                .into_iter()
                .map(|ids| {
                    let refs: Vec<_> = ids.iter().map(|&id| ref_(id)).collect();
                    ref_(self.emit("OPEN_SHELL", vec![str_param(""), list(refs)]))
                })
                .collect();
            return Ok(vec![self.emit(
                "SHELL_BASED_SURFACE_MODEL",
                vec![str_param(""), list(shell_refs)],
            )]);
        }
        let mut items = Vec::with_capacity(face_ids_by_shell.len());
        for ids in face_ids_by_shell {
            let refs: Vec<_> = ids.iter().map(|&id| ref_(id)).collect();
            let shell_id = self.emit("CLOSED_SHELL", vec![str_param(""), list(refs)]);
            items.push(self.emit("MANIFOLD_SOLID_BREP", vec![str_param(""), ref_(shell_id)]));
        }
        Ok(items)
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

        // Determine face orientation. Planes export with a
        // chart-true ẑ (= u_axis × v_axis); when that opposes the
        // stored normal the sense flag flips to compensate (see
        // `step_sense_flip`).
        let orientation = face.surface.as_ref().map(|s| s.orientation).unwrap_or(true)
            ^ face
                .surface
                .as_ref()
                .map(|s| crate::step::export_geometry::step_sense_flip(&s.geometry))
                .unwrap_or(false);

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

    // `maybe_wrap_surface_curve` lives in
    // `export_surface_curve_wrap.rs` (file-size ratchet).

    /// Try to export the pcurve for a single half-edge.
    ///
    /// Returns `None` if the half-edge has no pcurve or its face surface was
    /// not pre-exported.
    pub(crate) fn try_export_he_pcurve(
        &mut self,
        brep: &BRep,
        edge_handle: EdgeHandle,
        he_handle: popkern_core::topology::HalfEdgeHandle,
        face_handle: FaceHandle,
        surface_map: &HashMap<FaceHandle, ExportedFaceSurface>,
        branch_shifts: &super::export_periodic_branches::PeriodicNormalization,
    ) -> Option<EntityId> {
        let he = brep.half_edges.get(he_handle)?;
        if he.edge != edge_handle {
            return None;
        }
        let surface = surface_map.get(&face_handle)?;
        let pcurve = super::export_periodic_branches::export_pcurve_for_use(
            brep,
            face_handle,
            he_handle,
            &surface.surface.geometry,
        )?;
        // STEP has no analytic-pcurve type for most of popkern's
        // variants — convert to NURBS, exact where possible and
        // otherwise co-parameterized within a fraction of the model
        // tolerance (export_pcurve_nurbs, PK-550). NURBS variants
        // are zero-cost. Then, in order: co-parameterize with the
        // edge's basis curve (popkern pcurves run in HALF-EDGE
        // direction; STEP associated pcurves share the 3D curve's
        // parameterization — a pcurve parameterized differently from
        // an exact curve is rebuilt on the curve's parameter, PK-554),
        // apply the periodic branch shift that
        // chains this use into ITS face's UV loop contour, and map
        // popkern's chart convention to the written surface's STEP
        // parameterization (sphere latitude, `pcurve_to_step_chart`).
        let mut nurbs = super::export_pcurve_same_parameter::export_pcurve_on_edge_curve(
            brep,
            he_handle,
            &surface.surface,
            &pcurve,
        );
        if let Some(&(du, dv)) = branch_shifts.shifts.get(&(face_handle, he_handle)) {
            for p in &mut nurbs.control_points {
                p.x += du;
                p.y += dv;
            }
        }
        let nurbs =
            crate::step::export_geometry::pcurve_to_step_chart(&surface.surface.geometry, nurbs);
        Some(self.export_pcurve(&nurbs, surface.id))
    }

    pub(crate) fn try_export_he_pcurve_via_stored_loop(
        &mut self,
        brep: &BRep,
        edge_handle: EdgeHandle,
        he_handle: popkern_core::topology::HalfEdgeHandle,
        surface_map: &HashMap<FaceHandle, ExportedFaceSurface>,
        branch_shifts: &super::export_periodic_branches::PeriodicNormalization,
    ) -> Option<EntityId> {
        let he = brep.half_edges.get(he_handle)?;
        let loop_ = brep.loops.get(he.loop_)?;
        self.try_export_he_pcurve(
            brep,
            edge_handle,
            he_handle,
            loop_.face,
            surface_map,
            branch_shifts,
        )
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

#[cfg(test)]
#[path = "export_topology_fixture_tests.rs"]
mod fixture_tests;
