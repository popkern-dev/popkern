//! SURFACE_CURVE / SEAM_CURVE wrapping for edge export (extracted
//! from `export_topology.rs` under the file-size ratchet).

use std::collections::HashMap;

use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

use super::export::ExportContext;
use super::export_topology::ExportedFaceSurface;
use super::export_topology_face_uses::{edge_face_uses, unique_edge_half_edges};
use super::parser::types::EntityId;

impl ExportContext {
    /// If both half-edges have pcurves, wrap the 3-D curve in a SURFACE_CURVE
    /// (or SEAM_CURVE if both half-edges belong to the same face).
    ///
    /// Returns either the original `curve_3d_id` (no pcurves) or a new
    /// SURFACE_CURVE / SEAM_CURVE entity ID.
    pub(crate) fn maybe_wrap_surface_curve(
        &mut self,
        brep: &BRep,
        curve_3d_id: EntityId,
        edge_handle: EdgeHandle,
        surface_map: &HashMap<FaceHandle, ExportedFaceSurface>,
        branch_shifts: &super::export_periodic_branches::PeriodicNormalization,
    ) -> EntityId {
        let mut pcurve_ids: Vec<EntityId> = Vec::new();
        let mut face_uses = edge_face_uses(brep, edge_handle);
        // Canonical half-edge's pcurve FIRST: OpenCASCADE assigns the
        // first associated pcurve to the edge's FORWARD use in a
        // loop, and `export_loop_as_bound` writes forward exactly for
        // the canonical half-edge. For SEAM_CURVEs (both uses on one
        // face, branch-separated pcurves) the wrong order swaps the
        // seam's two branches and the trim contour breaks.
        if let Some(edge) = brep.edges.get(edge_handle) {
            let canonical = edge.half_edge;
            face_uses.sort_by_key(|face_use| face_use.half_edge != canonical);
        }
        for face_use in &face_uses {
            if let Some(pcurve_id) = self.try_export_he_pcurve(
                brep,
                edge_handle,
                face_use.half_edge,
                face_use.face,
                surface_map,
                branch_shifts,
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
                    branch_shifts,
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
            let same_face = first_face.is_some() && faces.all(|face| Some(face) == first_face);
            // A SEAM_CURVE promises OpenCASCADE that the edge appears
            // TWICE in the face's wire (two branch-separated pcurves,
            // two oriented uses). An edge whose OTHER half was masked
            // out of the loops (boolean strand residue) is used ONCE —
            // written as SEAM_CURVE it raised an exception inside
            // OCC's solid translation and the WHOLE body transferred
            // empty (measured: rotated torus × torus STEP, 0 faces
            // out of 11). Count live-loop uses and demote a once-used
            // edge to a plain SURFACE_CURVE.
            let live_uses = face_uses
                .iter()
                .filter(|face_use| {
                    brep.half_edges
                        .get(face_use.half_edge)
                        .is_some_and(|he| brep.loops.contains_key(he.loop_))
                })
                .count();
            let is_seam = same_face && live_uses >= 2;

            if is_seam {
                self.export_seam_curve(curve_3d_id, &pcurve_ids)
            } else {
                self.export_surface_curve(curve_3d_id, &pcurve_ids)
            }
        }
    }
}
