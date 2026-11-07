//! A body's faces grouped into shells, and the representation items
//! over them (`export_brep_items` steps 3 and 4).

use std::collections::HashMap;

use popkern_core::topology::{BRep, BodyKind, EdgeHandle, FaceHandle};

use super::{ExportedFaceSurface, TopologyExportMaps};
use crate::error::IoError;
use crate::step::export::{ExportContext, list, ref_, str_param};
use crate::step::export_topology_face_uses::EdgeFaceUse;
use crate::step::parser::types::EntityId;

impl ExportContext {
    /// Export every face, grouped by the shell behind it, and emit the
    /// items over them: a general, wire or acorn body as
    /// `export_general` writes it, a body with boundary edges (or
    /// face-local edges) as one surface model, else one solid per solid
    /// region with its cavities.
    pub(super) fn emit_shell_items(
        &mut self,
        brep: &BRep,
        surface_map: &HashMap<FaceHandle, ExportedFaceSurface>,
        edge_map: &HashMap<EdgeHandle, EntityId>,
        local_edge_map: &HashMap<EdgeFaceUse, EntityId>,
    ) -> Result<Vec<EntityId>, IoError> {
        let maps = TopologyExportMaps {
            edge: edge_map,
            local_edge: local_edge_map,
        };
        let mut face_ids_by_shell: Vec<Vec<EntityId>> = Vec::new();
        let mut handles: Vec<(FaceHandle, EntityId)> = Vec::new();
        let mut exported_shells = Vec::new();
        let mut shell_handles = brep.back_shells();
        shell_handles.sort();
        for &sh in &shell_handles {
            self.reverse_face_uses = super::super::export_voids::is_cavity(brep, sh);
            let mut ids = Vec::new();
            for (fh, face) in &brep.faces {
                if face.back != sh {
                    continue;
                }
                let surface_id = surface_map.get(&fh).map(|surface| surface.id);
                let face_id = self.export_face(brep, fh, face, &maps, surface_id)?;
                ids.push(face_id);
                handles.push((fh, face_id));
            }
            if !ids.is_empty() {
                face_ids_by_shell.push(ids);
                exported_shells.push(Some(sh));
            }
        }
        self.reverse_face_uses = false;
        // Faces with no live shell (defensive) travel together.
        let mut orphans = Vec::new();
        for (fh, face) in &brep.faces {
            if brep.face_back_shell(fh).is_some() {
                continue;
            }
            let surface_id = surface_map.get(&fh).map(|surface| surface.id);
            orphans.push(self.export_face(brep, fh, face, &maps, surface_id)?);
        }
        if !orphans.is_empty() {
            face_ids_by_shell.push(orphans);
            exported_shells.push(None);
        }

        // General, wire and acorn bodies, and a closed sheet (bounding a
        // void, no solid), take the general forms. Face-local STEP edge
        // curves preserve geometric fidelity for incoherent face uses,
        // but they also make STEP topology single-sided at those uses:
        // that representation is never declared a solid.
        // A solid opened mid-edit (a free edge on a face between two
        // regions: its open shell still in the solid's region) reads
        // General but is the sheet it will settle into, written in the
        // historical sheet form, never as a solid. A lamina's free edges
        // lie on faces with one region on both sides.
        let open = brep.has_boundary_edges();
        let opened_solid = open && has_free_edge_between_regions(brep);
        let kind = brep.body_kind();
        let closed_sheet = kind == BodyKind::Sheet && !open;
        let general = !opened_solid
            && (closed_sheet
                || matches!(kind, BodyKind::General | BodyKind::Wire | BodyKind::Acorn));
        if general && local_edge_map.is_empty() {
            return super::super::export_general::emit_general(self, brep, &handles);
        }
        if brep.has_boundary_edges() || !local_edge_map.is_empty() {
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
        Ok(super::super::export_voids::emit_solids(
            self,
            brep,
            &exported_shells,
            face_ids_by_shell,
        ))
    }
}

/// Whether some free edge of `brep` lies on a face with different
/// regions behind and in front (an opened solid's boundary).
fn has_free_edge_between_regions(brep: &BRep) -> bool {
    brep.edges
        .keys()
        .filter(|&e| brep.is_free_edge(e))
        .any(|e| {
            brep.radial_fins(e).into_iter().any(|fin| {
                brep.half_edges
                    .get(fin)
                    .and_then(|h| brep.loops.get(h.loop_))
                    .and_then(|l| brep.faces.get(l.face))
                    .is_some_and(|f| brep.shell_region(f.back) != brep.shell_region(f.front))
            })
        })
}
