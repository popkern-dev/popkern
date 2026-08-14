//! Region-scoped fast validation (PK-273).
//!
//! Whole-model fast validation is priced by model size, not edit
//! size: the per-face Gregory self-intersection sweep was measured at
//! ~25 ms × patch count (2.4 s on a 96-patch body), so acceptance
//! gates that validate once per candidate pay a cost that grows as
//! corners complete. The region API prices the sweep by the edit's
//! neighbourhood instead: callers pass the faces they dirtied, the
//! region dilates by one vertex-sharing ring, and only that region's
//! Gregory patches are checked. Global structural checks stay global
//! — they total well under a millisecond and scoping them would
//! change semantics for no measurable gain.
//!
//! **Containment.** A region-scoped result can miss defects wholly
//! outside the region. This is sound at candidate gates because they
//! diff against a baseline computed over the same scoping rule, and
//! because `validate_all` remains the explicit full serving gate at
//! Part boundaries and in corpus CI. See
//! docs/architecture/region-scoped-validation.md.

use std::collections::BTreeSet;

use crate::error::ValidationError;
use crate::topology::{BRep, FaceHandle, VertexHandle};

impl BRep {
    /// Dilate `dirty` by one face ring: the returned region is every
    /// live dirty face plus every face sharing at least one
    /// topological vertex with one (vertex sharing strictly contains
    /// edge sharing, so corner-straddling neighbours are included —
    /// the boundary-straddling defect test pins this). Unknown or
    /// dead handles are ignored; the result is deduplicated and
    /// deterministically ordered.
    pub fn validation_region(&self, dirty: &[FaceHandle]) -> Vec<FaceHandle> {
        let mut region: BTreeSet<FaceHandle> = dirty
            .iter()
            .copied()
            .filter(|face| self.faces.contains_key(*face))
            .collect();
        if region.is_empty() {
            return Vec::new();
        }
        let mut boundary_vertices: BTreeSet<VertexHandle> = BTreeSet::new();
        for &face in &region {
            self.collect_face_vertices(face, &mut boundary_vertices);
        }
        let mut candidate_vertices = BTreeSet::new();
        for face in self.faces.keys() {
            if region.contains(&face) {
                continue;
            }
            candidate_vertices.clear();
            self.collect_face_vertices(face, &mut candidate_vertices);
            if candidate_vertices
                .iter()
                .any(|vertex| boundary_vertices.contains(vertex))
            {
                region.insert(face);
            }
        }
        region.into_iter().collect()
    }

    /// As [`BRep::validate_collect_fast`], with the Gregory
    /// self-intersection sweep restricted to `region` (typically
    /// [`BRep::validation_region`] of the faces an operator dirtied).
    /// Structural checks (1–9 + precision) still run globally.
    /// Defects on Gregory faces outside `region` are NOT reported —
    /// see the module docs for why that containment is sound at
    /// candidate gates.
    pub fn validate_region_collect_fast(&self, region: &[FaceHandle]) -> Vec<ValidationError> {
        let mut errors = self.validate_collect_structural();
        if let Err(e) = self.validate_gregory_self_intersection_among(region) {
            errors.push(e);
        }
        errors
    }

    /// As [`BRep::validate_gregory_self_intersection`], restricted to
    /// the Gregory-surface faces in `faces` (same tolerance rule and
    /// first-failure semantics). Unknown, dead, duplicate, or
    /// non-Gregory handles are ignored.
    pub fn validate_gregory_self_intersection_among(
        &self,
        faces: &[FaceHandle],
    ) -> Result<(), ValidationError> {
        use crate::geometry::surface::SurfaceGeometry;

        let tol_linear = self.tolerance.linear().max(1e-6);
        let mut seen = BTreeSet::new();
        for &fh in faces {
            if !seen.insert(fh) {
                continue;
            }
            let Some(face) = self.faces.get(fh) else {
                continue;
            };
            let Some(sof) = face.surface.as_ref() else {
                continue;
            };
            let SurfaceGeometry::Gregory(ref patch) = sof.geometry else {
                continue;
            };
            if patch.self_intersects_with_tolerance(tol_linear) {
                return Err(ValidationError::TwistedSurface { face: fh });
            }
        }
        Ok(())
    }

    /// Push every vertex on `face`'s outer and inner loops into
    /// `out`. Walks are cycle-guarded; broken loops contribute the
    /// vertices reached before the break.
    fn collect_face_vertices(&self, face: FaceHandle, out: &mut BTreeSet<VertexHandle>) {
        let Some(face_data) = self.faces.get(face) else {
            return;
        };
        let loops =
            std::iter::once(face_data.outer_loop).chain(face_data.inner_loops.iter().copied());
        for loop_handle in loops {
            let Some(loop_data) = self.loops.get(loop_handle) else {
                continue;
            };
            let start = loop_data.half_edge;
            let mut cursor = start;
            for _ in 0..=self.half_edges.len() {
                let Some(half_edge) = self.half_edges.get(cursor) else {
                    break;
                };
                out.insert(half_edge.origin);
                cursor = half_edge.next;
                if cursor == start {
                    break;
                }
            }
        }
    }
}

#[cfg(test)]
#[path = "region_tests.rs"]
mod tests;
