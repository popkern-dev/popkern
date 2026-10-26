//! Stable names on the design body.
//!
//! The body keeps the design space's name pool, so its faces can name
//! their ancestors. Every face lying on a design space face (a mounting
//! face, a held bore; the design can meet it in several separate
//! patches) is that face modified: `Modified`, its old name the
//! ancestor, recorded `old → new`, so a reference to the old face
//! resolves to all its patches. The largest patch is indexed by the
//! source face's place in the design space (stable across designs: a
//! setup naming it survives a re-solve); the others by their rank
//! after it. The design's new surface is `Generated`, and edges and
//! vertices are minted as a snapshot mints them (`Lateral`,
//! `Junction`): the design's own geometry has no recipe a later edit
//! could replay onto, only the kept faces do.

use std::collections::BTreeMap;

use popkern_core::history::TrackingContext;
use popkern_core::naming::{EntityName, EntityRole, NameId, RoleCategory};
use popkern_core::naming_criteria::{
    GeometricCriteria, capture_edge, capture_face, capture_vertex,
};
use popkern_core::topology::{BRep, FaceHandle};

/// Twice the area of `face`'s outer loop (Newell), for ranking a source
/// face's fragments.
fn loop_area(brep: &BRep, face: FaceHandle) -> f64 {
    let p = brep.face_outer_loop_positions(face);
    let mut n = [0.0; 3];
    for i in 0..p.len() {
        let (a, b) = (p[i], p[(i + 1) % p.len()]);
        n[0] += (a.y - b.y) * (a.z + b.z);
        n[1] += (a.z - b.z) * (a.x + b.x);
        n[2] += (a.x - b.x) * (a.y + b.y);
    }
    (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt()
}

/// Names every entity of `body` (faces in key order, `sources` by face
/// in that order) under `tracking`'s operation, against the design space
/// `space` it replaces.
pub fn name_design(
    space: &BRep,
    body: &mut BRep,
    sources: &[Option<FaceHandle>],
    tracking: &mut TrackingContext,
) {
    body.name_pool = space.name_pool.clone();
    let operation = tracking.operation;
    let name = |category: RoleCategory, index: u32, ancestors: Vec<NameId>| EntityName {
        operation,
        role: EntityRole { category, index },
        ancestors,
    };
    let ordinal: BTreeMap<FaceHandle, u32> = space
        .faces
        .keys()
        .enumerate()
        .map(|(i, h)| (h, i as u32))
        .collect();
    let faces: Vec<FaceHandle> = body.faces.keys().collect();
    // Each source face's patches, largest first (ties in face order).
    let mut patches: BTreeMap<FaceHandle, Vec<(FaceHandle, f64)>> = BTreeMap::new();
    for (&f, source) in faces.iter().zip(sources) {
        if let Some(s) = source {
            patches.entry(*s).or_default().push((f, loop_area(body, f)));
        }
    }
    let count = space.faces.len() as u32;
    let mut index_of: BTreeMap<FaceHandle, u32> = BTreeMap::new();
    for (s, list) in patches.iter_mut() {
        list.sort_by(|a, b| b.1.total_cmp(&a.1));
        let base = ordinal.get(s).copied().unwrap_or(0);
        for (rank, (f, _)) in list.iter().enumerate() {
            // The largest patch at the source's ordinal; the rest after
            // every source's (unique per name: rank × face count).
            index_of.insert(*f, base + count * rank as u32);
        }
    }
    let mut generated = 0u32;
    for (&f, source) in faces.iter().zip(sources) {
        let old = source.and_then(|s| space.faces.get(s)?.name);
        let entity = match (index_of.get(&f), old) {
            (Some(&index), Some(old_id)) => name(RoleCategory::Modified, index, vec![old_id]),
            _ => {
                generated += 1;
                name(RoleCategory::Generated, generated - 1, Vec::new())
            }
        };
        let criteria = capture_face(body, f).map(GeometricCriteria::Face);
        let id = body.name_pool.intern(entity.clone(), criteria);
        body.faces[f].name = Some(id);
        match old.and_then(|o| space.name_pool.get(o).cloned()) {
            Some(old_name) => tracking.record_modified(old_name, entity),
            None => tracking.record_generated(entity),
        }
    }
    let edges: Vec<_> = body.edges.keys().collect();
    for (i, e) in edges.into_iter().enumerate() {
        let entity = name(RoleCategory::Lateral, i as u32, Vec::new());
        let criteria = capture_edge(body, e).map(GeometricCriteria::Edge);
        body.edges[e].name = Some(body.name_pool.intern(entity.clone(), criteria));
        tracking.record_generated(entity);
    }
    let vertices: Vec<_> = body.vertices.keys().collect();
    for (i, v) in vertices.into_iter().enumerate() {
        let entity = name(RoleCategory::Junction, i as u32, Vec::new());
        let criteria = capture_vertex(body, v).map(GeometricCriteria::Vertex);
        body.vertices[v].name = Some(body.name_pool.intern(entity.clone(), criteria));
        tracking.record_generated(entity);
    }
}
