//! Resolution of selections against a live body (spec §6.2): names
//! post seam-filter, describe text, labels used.

use popkern_core::naming::EntityName;
use popkern_core::topology::BRep;

use super::{EdgeSelector, EntityRef, FaceSelector};
use crate::error::{ScriptError, SelectionKind};

/// What the document stores for a selection (D3).
#[derive(Debug, Clone, PartialEq)]
pub struct SelectionSnapshot {
    /// Resolved stable names, post seam-filter.
    pub names: Vec<EntityName>,
    /// Canonical spelling of the query (diagnostics only).
    pub describe: String,
    /// Labels the query cited.
    pub labels_used: Vec<String>,
    /// Seam edges dropped from the match.
    pub seams_skipped: usize,
}

fn name_of(brep: &BRep, id: Option<popkern_core::naming::NameId>) -> Option<EntityName> {
    id.and_then(|nid| brep.name_pool.get(nid).cloned())
}

/// Resolve an edge selection against `brep`.
pub(crate) fn resolve_edges(
    brep: &BRep,
    selector: &EdgeSelector,
    resolve: &dyn Fn(EntityRef<'_>) -> Result<EntityName, ScriptError>,
) -> Result<SelectionSnapshot, ScriptError> {
    let describe = selector.describe();
    let query = selector.to_query(brep, resolve)?;
    let infos = popkern_part::edge_infos(brep);
    let handles = query.resolve(&infos);
    let seams: std::collections::HashSet<_> = infos
        .iter()
        .filter(|i| i.is_seam)
        .map(|i| i.handle)
        .collect();
    let seams_skipped = handles.iter().filter(|h| seams.contains(h)).count();
    let kept: Vec<_> = handles.iter().filter(|h| !seams.contains(h)).collect();
    let names: Vec<EntityName> = kept
        .iter()
        .filter_map(|h| brep.edges.get(**h).and_then(|e| name_of(brep, e.name)))
        .collect();
    if names.len() != kept.len() {
        return Err(ScriptError::Argument {
            argument: "edges".into(),
            reason: format!(
                "{} of {} matched edges carry no stable name (selection '{describe}')",
                kept.len() - names.len(),
                kept.len()
            ),
        });
    }
    if names.is_empty() && !selector.allow_empty_flag() {
        return Err(ScriptError::EmptySelection {
            selection: describe,
            kind: SelectionKind::Edges,
            available: infos.len(),
        });
    }
    Ok(SelectionSnapshot {
        names,
        describe,
        labels_used: selector.labels_used(),
        seams_skipped,
    })
}

/// Resolve a face selection against `brep`.
pub(crate) fn resolve_faces(
    brep: &BRep,
    selector: &FaceSelector,
    resolve: &dyn Fn(EntityRef<'_>) -> Result<EntityName, ScriptError>,
) -> Result<SelectionSnapshot, ScriptError> {
    let describe = selector.describe();
    let query = selector.to_query(brep, resolve)?;
    let infos = popkern_part::face_infos(brep);
    let handles = query.resolve(&infos);
    let names: Vec<EntityName> = handles
        .iter()
        .filter_map(|h| brep.faces.get(*h).and_then(|f| name_of(brep, f.name)))
        .collect();
    if names.len() != handles.len() {
        return Err(ScriptError::Argument {
            argument: "faces".into(),
            reason: format!(
                "{} of {} matched faces carry no stable name (selection '{describe}')",
                handles.len() - names.len(),
                handles.len()
            ),
        });
    }
    if names.is_empty() && !selector.allow_empty_flag() {
        return Err(ScriptError::EmptySelection {
            selection: describe,
            kind: SelectionKind::Faces,
            available: infos.len(),
        });
    }
    Ok(SelectionSnapshot {
        names,
        describe,
        labels_used: selector.labels_used(),
        seams_skipped: 0,
    })
}

/// Require exactly one entity.
pub(crate) fn exactly_one(
    snapshot: SelectionSnapshot,
    kind: SelectionKind,
) -> Result<EntityName, ScriptError> {
    if snapshot.names.len() != 1 {
        return Err(ScriptError::SelectionCardinality {
            selection: snapshot.describe,
            kind,
            expected: 1,
            got: snapshot.names.len(),
        });
    }
    snapshot
        .names
        .into_iter()
        .next()
        .ok_or(ScriptError::SelectionCardinality {
            selection: String::new(),
            kind,
            expected: 1,
            got: 0,
        })
}
