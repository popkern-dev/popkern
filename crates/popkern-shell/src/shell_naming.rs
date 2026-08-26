//! Stable-name assignment for shell results via EXACT face provenance
//! (PK-328): the `SkinMap`/`ShellAssemblyMaps` handle chain replaces
//! PK-325's centroid-proximity matching, which could not match a skin
//! displaced by any real thickness (measured: the moved skin of every
//! default-tolerance shell fell through to ancestor-less `Generated`
//! because the match radius was `100·tol` while the displacement was
//! the wall thickness).
//!
//! Roles: faces that CONTINUE a source face (the un-moved outer skin
//! for classic shells; the outward-moved skin for midplane) are
//! `Modified` with the source face as ancestor; the moved void skin's
//! faces are `Generated` with the source face as ancestor (a new face
//! derived from it); join faces are `Generated` with the gap edge /
//! corner vertex as ancestor; side walls and trim-created faces are
//! plain `Generated`. Result handles that died in the trim stage are
//! simply absent from the final body and skipped.

use std::collections::HashMap;

use popkern_core::history::TrackingContext;
use popkern_core::naming::{EntityName, EntityRole, OperationId, RoleCategory};
use popkern_core::topology::{BRep, FaceHandle};

use crate::skin_map::JoinKind;

/// Exact face provenance of a served shell/offset body, filled by
/// `shell_impl`/`midplane_shell_impl` after assembly (handles are
/// final from assembly on; trim-stage deletions are checked at
/// assignment time).
#[derive(Debug, Default, Clone)]
pub(crate) struct ShellProvenance {
    /// Result faces that continue a source face → that source face.
    pub modified: Vec<(FaceHandle, FaceHandle)>,
    /// Result faces offset from a source face (the moved skin) → that
    /// source face.
    pub offset: Vec<(FaceHandle, FaceHandle)>,
    /// Join faces → the source entity's name (gap edge for strips,
    /// corner vertex for sphere caps), resolved at capture; `None`
    /// when the source entity is unnamed or synthesized (rim-split
    /// midpoints).
    pub joins: Vec<(FaceHandle, Option<EntityName>)>,
}

impl ShellProvenance {
    /// Sort every list by result handle: the fill sites iterate
    /// `HashMap`s, and downstream consumers deserve a deterministic
    /// record.
    pub(crate) fn sort(&mut self) {
        self.modified.sort_by_key(|&(res, _)| res);
        self.offset.sort_by_key(|&(res, _)| res);
        self.joins.sort_by_key(|&(res, _)| res);
    }
}

/// Resolve the source entity name a join face descends from: the gap
/// edge for pipe strips, the corner vertex for sphere fills.
pub(crate) fn join_source_name(source: &BRep, kind: JoinKind) -> Option<EntityName> {
    match kind {
        JoinKind::PipeFill { source_edge } => source
            .edges
            .get(source_edge)
            .and_then(|e| e.name)
            .and_then(|id| source.name_pool.get(id))
            .cloned(),
        JoinKind::SphereFill { source_vertex } => source_vertex
            .and_then(|vh| source.vertices.get(vh))
            .and_then(|v| v.name)
            .and_then(|id| source.name_pool.get(id))
            .cloned(),
    }
}

/// Ancestry role of one result face, resolved from the provenance.
enum FaceRole<'a> {
    Modified(FaceHandle),
    OffsetFrom(FaceHandle),
    Join(&'a Option<EntityName>),
}

/// Assign stable names to the faces of a shell/offset result from the
/// exact provenance record. Faces without a provenance entry (side
/// walls, corner sectors, trim-created caps) become plain `Generated`.
pub(crate) fn assign_shell_names(
    result: &mut BRep,
    ctx: &mut TrackingContext,
    source: &BRep,
    provenance: &ShellProvenance,
) {
    use popkern_core::naming_criteria::{GeometricCriteria, capture_face};

    // Source-face names: the input's own name where one exists, else a
    // stable positional placeholder (pre-history bodies) — the same
    // fallback the heuristic used, keyed by source iteration index.
    let mut source_names: HashMap<FaceHandle, EntityName> = HashMap::new();
    for (i, (fh, face)) in source.faces.iter().enumerate() {
        let name = face
            .name
            .and_then(|id| source.name_pool.get(id))
            .cloned()
            .unwrap_or(EntityName {
                operation: OperationId(0),
                role: EntityRole {
                    category: RoleCategory::Generated,
                    index: i as u32,
                },
                ancestors: Vec::new(),
            });
        source_names.insert(fh, name);
    }

    let mut roles: HashMap<FaceHandle, FaceRole<'_>> = HashMap::new();
    for &(res, src) in &provenance.modified {
        roles.insert(res, FaceRole::Modified(src));
    }
    for &(res, src) in &provenance.offset {
        roles.insert(res, FaceRole::OffsetFrom(src));
    }
    for (res, name) in &provenance.joins {
        roles.insert(*res, FaceRole::Join(name));
    }

    let post_faces: Vec<FaceHandle> = result.faces.keys().collect();
    let mut gen_idx: u32 = 0;
    let mut mod_idx: u32 = 0;
    for fh in post_faces {
        let Some(fc) = capture_face(result, fh) else {
            continue;
        };
        // (category, ancestor name, record_modified's "old" name)
        let (category, ancestor, modified_from) = match roles.get(&fh) {
            Some(FaceRole::Modified(src)) => {
                let n = source_names.get(src).cloned();
                (RoleCategory::Modified, n.clone(), n)
            }
            Some(FaceRole::OffsetFrom(src)) => (
                RoleCategory::Generated,
                source_names.get(src).cloned(),
                None,
            ),
            Some(FaceRole::Join(name)) => (RoleCategory::Generated, (*name).clone(), None),
            None => (RoleCategory::Generated, None, None),
        };
        let ancestors = match ancestor {
            Some(a) => vec![result.name_pool.intern(a, None)],
            None => Vec::new(),
        };
        let index = if category == RoleCategory::Modified {
            let i = mod_idx;
            mod_idx += 1;
            i
        } else {
            let i = gen_idx;
            gen_idx += 1;
            i
        };
        let new_name = EntityName {
            operation: ctx.operation,
            role: EntityRole { category, index },
            ancestors,
        };
        let new_id = result
            .name_pool
            .intern(new_name.clone(), Some(GeometricCriteria::Face(fc)));
        result.faces[fh].name = Some(new_id);
        match modified_from {
            Some(old) => ctx.record_modified(old, new_name),
            None => ctx.record_generated(new_name),
        }
    }
}
