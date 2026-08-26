//! Region-removal surgery: delete a connected (or disconnected) set of
//! faces, leaving the surviving mesh's rim registered as boundary
//! edges (null-sentinel twins) — the multi-face generalization of the
//! single-face open-shell surgery behind [`BRep::open_box`].
//!
//! Used by global self-intersection trimming (PK-327): the unrealized
//! region of an offset skin dies as a unit and the survivors are
//! re-closed by welding the rim pairs
//! ([`BRep::weld_boundary_edge_pair`]).

use std::collections::BTreeSet;

use crate::error::CoreError;
use crate::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle};

/// Result of [`BRep::remove_face_region_to_boundary`].
#[derive(Debug, Clone, Default)]
pub struct RegionRemovalResult {
    /// Edges whose surviving side became a boundary (the region rim).
    pub rim_edges: Vec<EdgeHandle>,
    /// Edges deleted entirely (interior to the region, or boundary
    /// edges that died with their only face).
    pub removed_edges: usize,
    /// Vertices deleted (no incident half-edges survived).
    pub removed_vertices: usize,
}

impl BRep {
    /// Remove every face in `region`, registering the rim (edges shared
    /// with surviving faces) as boundary edges. Interior edges and
    /// vertices of the region are deleted.
    ///
    /// # Errors
    ///
    /// - [`CoreError::PreconditionFailed`] if the region is empty,
    ///   names a missing face, or covers the whole body (nothing would
    ///   survive).
    /// - [`CoreError::Validation`] if the surviving body fails
    ///   structural validation (the model is left untouched).
    pub fn remove_face_region_to_boundary(
        &mut self,
        region: &BTreeSet<FaceHandle>,
    ) -> Result<RegionRemovalResult, CoreError> {
        if region.is_empty() {
            return Err(CoreError::PreconditionFailed(
                "remove_face_region_to_boundary: empty region".into(),
            ));
        }
        for &fh in region {
            if !self.faces.contains_key(fh) {
                return Err(CoreError::PreconditionFailed(format!(
                    "remove_face_region_to_boundary: face {fh:?} is not in the body"
                )));
            }
        }
        if region.len() >= self.faces.len() {
            return Err(CoreError::PreconditionFailed(
                "remove_face_region_to_boundary: region covers the whole body".into(),
            ));
        }

        let mut candidate = self.clone();
        let mut result = RegionRemovalResult::default();

        // Collect the doomed loops and half-edges.
        let mut doomed_loops = Vec::new();
        let mut doomed_hes: Vec<HalfEdgeHandle> = Vec::new();
        for &fh in region {
            let face = &candidate.faces[fh];
            for lh in std::iter::once(face.outer_loop).chain(face.inner_loops.iter().copied()) {
                doomed_loops.push(lh);
                let Some(loop_data) = candidate.loops.get(lh) else {
                    continue;
                };
                let start = loop_data.half_edge;
                if !candidate.half_edges.contains_key(start) {
                    continue; // Empty ring (kvr-shaped) — loop dies with the face.
                }
                let mut cursor = start;
                let bound = candidate.half_edges.len().max(1);
                for _ in 0..=bound {
                    doomed_hes.push(cursor);
                    cursor = candidate.half_edges[cursor].next;
                    if cursor == start {
                        break;
                    }
                }
            }
        }
        let doomed_set: BTreeSet<HalfEdgeHandle> = doomed_hes.iter().copied().collect();

        // Classify each doomed half-edge's parent edge: rim (surviving
        // twin becomes a boundary) or interior (edge dies).
        let mut dying_edges: BTreeSet<EdgeHandle> = BTreeSet::new();
        for &heh in &doomed_hes {
            let he = &candidate.half_edges[heh];
            let edge_h = he.edge;
            let twin_h = he.twin;
            let twin_survives =
                candidate.half_edges.contains_key(twin_h) && !doomed_set.contains(&twin_h);
            if twin_survives {
                // Rim: repoint the edge at the surviving twin, keeping
                // the canonical-forward convention true (the twin
                // traverses the curve backward — PK-274).
                if candidate.edges[edge_h].half_edge == heh
                    && let Some(curve) = candidate.edges[edge_h].curve.as_mut()
                {
                    std::mem::swap(&mut curve.t_start, &mut curve.t_end);
                }
                candidate.edges[edge_h].half_edge = twin_h;
                candidate.half_edges[twin_h].twin = HalfEdgeHandle::default();
                candidate.boundary_edges.insert(edge_h);
                result.rim_edges.push(edge_h);
            } else {
                dying_edges.insert(edge_h);
            }
        }

        // Remove the doomed half-edges (anchor-aware), then the dying
        // edges, loops, and faces.
        for &heh in &doomed_hes {
            candidate.remove_half_edge_reanchored(heh);
        }
        for &eh in &dying_edges {
            candidate.edges.remove(eh);
            candidate.boundary_edges.remove(&eh);
        }
        result.removed_edges = dying_edges.len();
        for lh in doomed_loops {
            candidate.loops.remove(lh);
        }
        for &fh in region {
            candidate.faces.remove(fh);
        }

        // Vertices with no surviving incident half-edge die with the
        // region.
        let incident: BTreeSet<_> = candidate.half_edges.values().map(|he| he.origin).collect();
        let orphaned: Vec<_> = candidate
            .vertices
            .keys()
            .filter(|vh| !incident.contains(vh))
            .collect();
        for vh in orphaned {
            candidate.vertices.remove(vh);
            result.removed_vertices += 1;
        }

        // Repoint shell anchors that referenced removed faces.
        let shells: Vec<_> = candidate.shells.keys().collect();
        for sh in shells {
            if !candidate.faces.contains_key(candidate.shells[sh].face)
                && let Some((replacement, _)) = candidate.faces.iter().next()
            {
                candidate.shells[sh].face = replacement;
            }
        }

        candidate.validate().map_err(CoreError::Validation)?;
        *self = candidate;
        Ok(result)
    }
}
