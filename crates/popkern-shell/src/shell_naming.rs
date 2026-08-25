//! Stable-name assignment for shell results (PK-325 row 474):
//! centroid/normal matching of output faces to input faces (Modified)
//! and of join faces to their source entities (Generated with the gap
//! edge or corner vertex as ancestor).

use popkern_core::geometry::{Point3, Vec3};
use popkern_core::history::TrackingContext;
use popkern_core::naming::EntityName;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use crate::skin_map;

/// Collected face data for ancestry matching during shell naming.
pub(crate) struct InputFaceData {
    /// The entity name (or placeholder).
    name: EntityName,
    /// Boundary centroid.
    centroid: Point3,
    /// Surface normal at the centroid.
    normal: Vec3,
}

/// Naming data for one join face (PK-325): captured on the skin at
/// construction time (handles do not survive assembly, centroids do)
/// so `assign_shell_names` can match the face in the final body and
/// record `Generated` with the SOURCE entity as ancestor (row 474).
pub(crate) struct JoinFaceData {
    /// The source entity's name, when the input carries one: the gap
    /// edge for strips, the corner vertex for sphere caps.
    source_name: Option<EntityName>,
    /// Join-face boundary centroid on the skin.
    centroid: Point3,
    /// Join-face surface normal at the centroid.
    normal: Vec3,
}

/// Capture naming data for freshly constructed join faces.
pub(crate) fn collect_join_face_data(
    skin: &BRep,
    outer: &BRep,
    records: &[skin_map::JoinRecord],
) -> Vec<JoinFaceData> {
    use popkern_core::naming_criteria::capture_face;
    let mut out = Vec::new();
    for record in records {
        let Some(fc) = capture_face(skin, record.face) else {
            continue;
        };
        let source_name = match record.kind {
            skin_map::JoinKind::PipeFill { source_edge } => outer
                .edges
                .get(source_edge)
                .and_then(|e| e.name)
                .and_then(|id| outer.name_pool.get(id))
                .cloned(),
            skin_map::JoinKind::SphereFill { source_vertex } => source_vertex
                .and_then(|vh| outer.vertices.get(vh))
                .and_then(|v| v.name)
                .and_then(|id| outer.name_pool.get(id))
                .cloned(),
        };
        out.push(JoinFaceData {
            source_name,
            centroid: fc.centroid,
            normal: fc.normal,
        });
    }
    out
}

/// Collect face names, centroids, and normals from a BRep for ancestry matching.
pub(crate) fn collect_shell_face_data(brep: &BRep) -> Vec<InputFaceData> {
    use popkern_core::naming::{EntityRole, OperationId, RoleCategory};
    use popkern_core::naming_criteria::capture_face;

    let mut result = Vec::new();
    for (i, (fh, _face)) in brep.faces.iter().enumerate() {
        let Some(fc) = capture_face(brep, fh) else {
            continue;
        };
        let name = if let Some(name_id) = brep.faces[fh].name {
            if let Some(n) = brep.name_pool.get(name_id) {
                n.clone()
            } else {
                EntityName {
                    operation: OperationId(0),
                    role: EntityRole {
                        category: RoleCategory::Generated,
                        index: i as u32,
                    },
                    ancestors: Vec::new(),
                }
            }
        } else {
            EntityName {
                operation: OperationId(0),
                role: EntityRole {
                    category: RoleCategory::Generated,
                    index: i as u32,
                },
                ancestors: Vec::new(),
            }
        };
        result.push(InputFaceData {
            name,
            centroid: fc.centroid,
            normal: fc.normal,
        });
    }
    result
}

/// Assign stable names to entities produced by a shell operation.
///
/// Matches output faces to input faces by centroid proximity and normal
/// alignment. Outer skin faces that match an input face (same-direction
/// normal) get `Modified` with the input face as ancestor. Inner skin
/// faces that match an input face (opposite normal) also get `Modified`
/// with the input face as ancestor. Side wall faces (no match) get
/// `Generated`.
pub(crate) fn assign_shell_names(
    brep: &mut BRep,
    ctx: &mut TrackingContext,
    input_faces: &[InputFaceData],
    join_faces: &[JoinFaceData],
    tol: &Tolerance,
) {
    use popkern_core::naming::{EntityRole, RoleCategory};
    use popkern_core::naming_criteria::{GeometricCriteria, capture_face};

    let post_faces: Vec<FaceHandle> = brep.faces.keys().collect();
    let mut gen_idx: u32 = 0;
    let mut mod_idx: u32 = 0;
    let generous_dist = tol.linear() * 100.0;

    for fh in &post_faces {
        let Some(fc) = capture_face(brep, *fh) else {
            continue;
        };
        let centroid = fc.centroid;
        let normal = fc.normal;

        // Match against input faces: allow same-direction (outer) or
        // opposite-direction (inner/offset) normals.
        let best_match = input_faces
            .iter()
            .filter(|inp| {
                centroid.distance_to(&inp.centroid) < generous_dist
                    && normal.dot(&inp.normal).abs() > 0.9
            })
            .min_by(|a, b| {
                centroid
                    .distance_to(&a.centroid)
                    .partial_cmp(&centroid.distance_to(&b.centroid))
                    .unwrap_or(std::cmp::Ordering::Equal)
            });

        if let Some(input) = best_match {
            let ancestor_id = brep.name_pool.intern(input.name.clone(), None);
            let new_name = EntityName {
                operation: ctx.operation,
                role: EntityRole {
                    category: RoleCategory::Modified,
                    index: mod_idx,
                },
                ancestors: vec![ancestor_id],
            };
            let criteria = Some(GeometricCriteria::Face(fc));
            let new_id = brep.name_pool.intern(new_name.clone(), criteria);
            brep.faces[*fh].name = Some(new_id);
            ctx.record_modified(input.name.clone(), new_name);
            mod_idx += 1;
        } else {
            // Join faces (PK-325 row 474): Generated with the SOURCE
            // entity (the gap edge / corner vertex) as ancestor,
            // matched by the centroid captured at construction time.
            let join_match = join_faces
                .iter()
                .filter(|j| {
                    centroid.distance_to(&j.centroid) < generous_dist
                        && normal.dot(&j.normal).abs() > 0.9
                })
                .min_by(|a, b| {
                    centroid
                        .distance_to(&a.centroid)
                        .partial_cmp(&centroid.distance_to(&b.centroid))
                        .unwrap_or(std::cmp::Ordering::Equal)
                });
            let ancestors = match join_match.and_then(|j| j.source_name.as_ref()) {
                Some(source) => vec![brep.name_pool.intern(source.clone(), None)],
                None => Vec::new(),
            };
            let new_name = EntityName {
                operation: ctx.operation,
                role: EntityRole {
                    category: RoleCategory::Generated,
                    index: gen_idx,
                },
                ancestors,
            };
            let criteria = Some(GeometricCriteria::Face(fc));
            let new_id = brep.name_pool.intern(new_name.clone(), criteria);
            brep.faces[*fh].name = Some(new_id);
            ctx.record_generated(new_name);
            gen_idx += 1;
        }
    }
}
