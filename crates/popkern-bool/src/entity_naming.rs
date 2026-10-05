//! Entity naming helpers for boolean operations.
//!
//! `collect_face_data` snapshots input faces before the boolean
//! consumes them; `assign_boolean_names` matches output faces back
//! against those inputs and assigns `Modified`/`Generated` names.
//!
//! Split out of `lib.rs` to keep the parent under the 500-line cap.

use popkern_core::geometry::{Point3, Vec3};
use popkern_core::history::TrackingContext;
use popkern_core::naming::{EntityName, NamePool};
use popkern_core::{tolerance::Tolerance, topology::BRep};

/// Collect face names, centroids, and normals from a BRep for ancestry matching.
///
/// For faces without an existing name, creates a placeholder name so that all
/// faces can participate in the matching step.
pub(crate) fn collect_face_data(brep: &BRep) -> Vec<(EntityName, Point3, Vec3)> {
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
                // Fallback placeholder for faces with an invalid name id.
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
            // Placeholder for unnamed faces.
            EntityName {
                operation: OperationId(0),
                role: EntityRole {
                    category: RoleCategory::Generated,
                    index: i as u32,
                },
                ancestors: Vec::new(),
            }
        };
        result.push((name, fc.centroid, fc.normal));
    }
    result
}

/// Naming scheme 2 (PK-677): replace `name` (a name of `source`'s pool)
/// with its copy imported by value into `target`, so its ancestor ids
/// are `target`'s. A placeholder for an unnamed entity is not in the
/// pool and stays as it is.
pub(crate) fn import_name(name: &mut EntityName, source: &BRep, target: &mut NamePool) {
    let Some(id) = source.name_pool.find(name) else {
        return;
    };
    if let Some(copied) = target.import(&source.name_pool, id)
        && let Some(value) = target.get(copied)
    {
        *name = value.clone();
    }
}

/// Assign stable names to faces in a boolean result by matching against inputs.
///
/// Output faces that geometrically match an input face (by centroid proximity
/// and normal alignment) receive `Modified` names with the input face as
/// ancestor. Genuinely new intersection faces receive `Generated` names.
pub(crate) fn assign_boolean_names(
    brep: &mut BRep,
    input_a_faces: &[(EntityName, Point3, Vec3)],
    input_b_faces: &[(EntityName, Point3, Vec3)],
    ctx: &mut TrackingContext,
    tol: &Tolerance,
) {
    use popkern_core::naming::{EntityRole, RoleCategory};
    use popkern_core::naming_criteria::{GeometricCriteria, capture_face};

    let all_inputs: Vec<&(EntityName, Point3, Vec3)> =
        input_a_faces.iter().chain(input_b_faces.iter()).collect();

    let face_handles: Vec<_> = brep.faces.keys().collect();
    let mut gen_idx: u32 = 0;
    let mut mod_idx: u32 = 0;

    for fh in &face_handles {
        let Some(fc) = capture_face(brep, *fh) else {
            continue;
        };
        let centroid = fc.centroid;
        let normal = fc.normal;

        // Find the closest input face by centroid + normal alignment.
        let generous_dist = tol.linear() * 100.0;
        let best_match = all_inputs
            .iter()
            .filter(|(_, ic, in_)| {
                centroid.distance_to(ic) < generous_dist && normal.dot(in_).abs() > 0.9
            })
            .min_by(|a, b| {
                centroid
                    .distance_to(&a.1)
                    .partial_cmp(&centroid.distance_to(&b.1))
                    .unwrap_or(std::cmp::Ordering::Equal)
            });

        if let Some((input_name, _, _)) = best_match {
            // Modified: this face descended from an input face.
            let ancestor_id = brep.name_pool.intern(input_name.clone(), None);
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
            ctx.record_modified(input_name.clone(), new_name);
            mod_idx += 1;
        } else {
            // Generated: genuinely new intersection face.
            let new_name = EntityName {
                operation: ctx.operation,
                role: EntityRole {
                    category: RoleCategory::Generated,
                    index: gen_idx,
                },
                ancestors: Vec::new(),
            };
            let criteria = Some(GeometricCriteria::Face(fc));
            let new_id = brep.name_pool.intern(new_name.clone(), criteria);
            brep.faces[*fh].name = Some(new_id);
            ctx.record_generated(new_name);
            gen_idx += 1;
        }
    }
}
