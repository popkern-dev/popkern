//! Phase 6c-C.1 — apex-region detection.
//!
//! Given an apex `VertexHandle` in a post-fillet BRep (after the
//! cylinder fillets are inserted but before `insert_vertex_blend`
//! creates the cap), find the topological sub-graph that 6c-C
//! must replace.
//!
//! Per design supplement §10.1, detection is **topological** (not
//! distance-based) so the rule is invariant under scale and
//! O(valence(V_apex)) rather than O(BRep size):
//!
//! * The apex vertex itself.
//! * All edges incident to the apex (every half-edge whose
//!   `origin == apex`).
//! * All faces in those edges' loops that *also* have the apex
//!   on their outer loop. These are the cube faces with their
//!   apex corners truncated (cube N=3: three faces) plus any
//!   small "corner" sub-faces left behind by Phase A/B trimming.
//! * The post-`kef` blend face if it already exists (Phase 6c-C
//!   runs before `insert_vertex_blend`, so this is normally
//!   `None`; the field is here for future call-site flexibility).
//!
//! This phase produces the data structure used by 6c-C.2
//! (deletion) and 6c-C.3 (stitching). Pure read of the BRep —
//! no mutation.

use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

/// Sub-graph around an apex vertex that Phase 6c-C will delete
/// and rebuild.  All handles are guaranteed to exist in `brep`
/// at the time `detect_apex_region` returns; callers must not
/// hold the report across topology mutations of unrelated
/// regions either (handles can become stale in slot-map).
///
/// Currently consumed only by tests (Phase 6c-C.1 is the
/// detection deliverable; .2 and .3 land in subsequent commits
/// and will use these handles for deletion + stitching).
#[allow(dead_code)]
#[derive(Debug, Clone, Default)]
pub(crate) struct ApexRegion {
    /// The apex vertex itself.
    pub(crate) apex: Option<VertexHandle>,
    /// All edges with `origin == apex` for at least one of their
    /// half-edges. Deduplicated; deterministic order (sorted by
    /// EdgeHandle key).
    pub(crate) edges: Vec<EdgeHandle>,
    /// Faces whose outer loop touches the apex vertex.
    /// Deduplicated; deterministic order.
    pub(crate) faces: Vec<FaceHandle>,
    /// Vertices in the apex's 1-hop neighborhood via `edges`.
    /// Includes the apex itself plus the "frontier" vertices on
    /// the far side of each incident edge. Deterministic order.
    pub(crate) vertices: Vec<VertexHandle>,
    /// Post-kef triangular blend face if Phase 6c-C runs *after*
    /// `insert_vertex_blend`. Normally `None` for the canonical
    /// "before insert_vertex_blend" call site.
    pub(crate) blend_face: Option<FaceHandle>,
}

impl ApexRegion {
    /// True when the detection found a non-empty apex region —
    /// i.e. the apex vertex is in the BRep and has incident
    /// edges. False for an unsupported apex (e.g. apex was
    /// deleted by an earlier pass, or the BRep is malformed).
    #[allow(dead_code)]
    pub(crate) fn is_valid(&self) -> bool {
        self.apex.is_some() && !self.edges.is_empty() && !self.faces.is_empty()
    }
}

/// Phase 6c-C.1 entry point.  Topological detection — see module
/// docstring for the rule.
#[allow(dead_code)]
pub(crate) fn detect_apex_region(brep: &BRep, apex: VertexHandle) -> ApexRegion {
    if brep.vertices.get(apex).is_none() {
        return ApexRegion::default();
    }

    let mut edges: std::collections::BTreeSet<EdgeHandle> = std::collections::BTreeSet::new();
    let mut faces: std::collections::BTreeSet<FaceHandle> = std::collections::BTreeSet::new();
    let mut vertices: std::collections::BTreeSet<VertexHandle> = std::collections::BTreeSet::new();
    vertices.insert(apex);

    for he in brep.half_edges.values() {
        if he.origin != apex {
            continue;
        }
        edges.insert(he.edge);
        if let Some(loop_data) = brep.loops.get(he.loop_) {
            faces.insert(loop_data.face);
        }
        // The frontier vertex is the destination of this
        // half-edge (the other end of the incident edge).
        if let Some(twin) = brep.half_edges.get(he.twin) {
            vertices.insert(twin.origin);
        }
    }

    ApexRegion {
        apex: Some(apex),
        edges: edges.into_iter().collect(),
        faces: faces.into_iter().collect(),
        vertices: vertices.into_iter().collect(),
        blend_face: None,
    }
}

/// Outer-loop valence of `face` (number of outer-loop
/// half-edges).  Used by post-`insert_vertex_blend` callers to
/// classify the apex-region faces.  Returns 0 if the face or
/// loop is missing.
#[allow(dead_code)]
pub(crate) fn outer_loop_valence(brep: &BRep, face: FaceHandle) -> usize {
    let Some(face_data) = brep.faces.get(face) else {
        return 0;
    };
    let Some(loop_data) = brep.loops.get(face_data.outer_loop) else {
        return 0;
    };
    let start = loop_data.half_edge;
    let mut count = 0usize;
    let mut cursor = start;
    for _ in 0..32 {
        count += 1;
        let Some(he) = brep.half_edges.get(cursor) else {
            return count;
        };
        cursor = he.next;
        if cursor == start {
            break;
        }
    }
    count
}
