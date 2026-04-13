//! Stable-name assignment for entities produced by fillet/chamfer
//! operations.
//!
//! Extracted from `lib.rs` to keep each source file under the
//! 500-line structural limit.

use popkern_core::history::TrackingContext;
use popkern_core::topology::{BRep, FaceHandle};

/// Assign stable names to entities produced by a fillet or chamfer
/// operation.
///
/// Surviving faces from the original set receive `Modified` names
/// with the original face name as ancestor (tracked via
/// `record_modified`). New fillet/chamfer faces receive `Generated`
/// names.
pub(crate) fn assign_fillet_names(
    brep: &mut BRep,
    ctx: &mut TrackingContext,
    pre_faces: &std::collections::HashSet<FaceHandle>,
) {
    use popkern_core::naming::{EntityName, EntityRole, RoleCategory};
    use popkern_core::naming_criteria::{GeometricCriteria, capture_face};

    let post_faces: Vec<FaceHandle> = brep.faces.keys().collect();

    let mut gen_idx: u32 = 0;
    let mut mod_idx: u32 = 0;

    for fh in &post_faces {
        if pre_faces.contains(fh) {
            // Surviving face — Modified (trimmed by the fillet).
            let old_name = brep.faces[*fh]
                .name
                .and_then(|nid| brep.name_pool.get(nid).cloned());

            let ancestors: Vec<_> = brep.faces[*fh].name.into_iter().collect();
            let new_name = EntityName {
                operation: ctx.operation,
                role: EntityRole {
                    category: RoleCategory::Modified,
                    index: mod_idx,
                },
                ancestors,
            };
            let criteria = capture_face(brep, *fh).map(GeometricCriteria::Face);
            let name_id = brep.name_pool.intern(new_name.clone(), criteria);
            brep.faces[*fh].name = Some(name_id);

            if let Some(old) = old_name {
                ctx.record_modified(old, new_name);
            } else {
                ctx.record_generated(new_name);
            }
            mod_idx += 1;
        } else {
            // New fillet/chamfer face — Generated.
            let new_name = EntityName {
                operation: ctx.operation,
                role: EntityRole {
                    category: RoleCategory::Generated,
                    index: gen_idx,
                },
                ancestors: Vec::new(),
            };
            let criteria = capture_face(brep, *fh).map(GeometricCriteria::Face);
            let name_id = brep.name_pool.intern(new_name.clone(), criteria);
            brep.faces[*fh].name = Some(name_id);
            ctx.record_generated(new_name);
            gen_idx += 1;
        }
    }
}
