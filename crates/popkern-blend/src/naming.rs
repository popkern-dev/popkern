//! Stable-name assignment for entities produced by fillet/chamfer
//! operations (faces since the naming substrate; edges and vertices
//! since PK-407).
//!
//! Faces: a surviving face is `Modified` with its old name as
//! ancestor, a new blend face is `Generated`. Edges and vertices
//! follow the boolean rule (PK-386) with one refinement: an output
//! edge that lies on an input edge's carrier — collinear for straight
//! edges, coincident midpoint otherwise — is `Modified` with that
//! ancestor even when the blend shortened it (a fillet trims every
//! edge that meets the blended one), and the strip rims, rails and
//! cap edges are `Lateral` edges / `Junction` vertices of this
//! operation. Every edge and vertex of the served body carries a
//! name, so scripts can select after a blend and documents can spell
//! the result.

use popkern_core::geometry::{Point3, Vec3};
use popkern_core::history::TrackingContext;
use popkern_core::naming::{EntityName, EntityRole, RoleCategory};
use popkern_core::naming_criteria::{
    CurveType, GeometricCriteria, capture_edge, capture_face, capture_vertex,
};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

/// An input edge: name (or a placeholder for unnamed input), midpoint,
/// tangent, length, and whether it is straight.
pub struct EdgeSnapshot {
    name: EntityName,
    midpoint: Point3,
    tangent: Vec3,
    length: f64,
    straight: bool,
}

/// An input vertex: name (or placeholder) and position.
pub struct VertexSnapshot {
    name: EntityName,
    position: Point3,
}

/// Everything the naming pass needs from the body before a blend.
pub struct PreOperationSnapshot {
    /// Face handles present before the operation.
    pub(crate) faces: std::collections::HashSet<FaceHandle>,
    edges: Vec<EdgeSnapshot>,
    vertices: Vec<VertexSnapshot>,
}

fn placeholder(category: RoleCategory, index: usize) -> EntityName {
    EntityName {
        operation: popkern_core::naming::OperationId(0),
        role: EntityRole {
            category,
            index: index as u32,
        },
        ancestors: Vec::new(),
    }
}

impl PreOperationSnapshot {
    /// Snapshot faces, edges, and vertices of the input body.
    pub fn capture(brep: &BRep) -> PreOperationSnapshot {
        let faces = brep.faces.keys().collect();
        let edges = brep
            .edges
            .iter()
            .enumerate()
            .filter_map(|(i, (eh, edge))| {
                let criteria = capture_edge(brep, eh)?;
                let name = edge
                    .name
                    .and_then(|nid| brep.name_pool.get(nid).cloned())
                    .unwrap_or_else(|| placeholder(RoleCategory::Lateral, i));
                Some(EdgeSnapshot {
                    name,
                    midpoint: criteria.midpoint,
                    tangent: criteria.tangent,
                    length: criteria.length,
                    straight: matches!(criteria.curve_type, CurveType::Line),
                })
            })
            .collect();
        let vertices = brep
            .vertices
            .iter()
            .enumerate()
            .filter_map(|(i, (vh, vertex))| {
                let criteria = capture_vertex(brep, vh)?;
                let name = vertex
                    .name
                    .and_then(|nid| brep.name_pool.get(nid).cloned())
                    .unwrap_or_else(|| placeholder(RoleCategory::Junction, i));
                Some(VertexSnapshot {
                    name,
                    position: criteria.position,
                })
            })
            .collect();
        PreOperationSnapshot {
            faces,
            edges,
            vertices,
        }
    }
}

/// Assign stable names to the faces produced by a fillet or chamfer
/// operation (surviving faces `Modified`, new faces `Generated`).
pub(crate) fn assign_fillet_names(
    brep: &mut BRep,
    ctx: &mut TrackingContext,
    pre_faces: &std::collections::HashSet<FaceHandle>,
) {
    let post_faces: Vec<FaceHandle> = brep.faces.keys().collect();
    let mut gen_idx: u32 = 0;
    let mut mod_idx: u32 = 0;
    for fh in &post_faces {
        if pre_faces.contains(fh) {
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

/// Distance from `point` to the infinite line through `origin` along
/// the unit `direction`.
fn line_distance(point: &Point3, origin: &Point3, direction: &Vec3) -> f64 {
    let d = *point - *origin;
    let along = d.dot(direction);
    (d - *direction * along).length()
}

/// The input edge an output edge continues, if any: same carrier and
/// within the input's extent.
fn matching_input<'a>(
    inputs: &'a [EdgeSnapshot],
    midpoint: &Point3,
    tangent: &Vec3,
    straight: bool,
    band: f64,
) -> Option<&'a EdgeSnapshot> {
    inputs
        .iter()
        .filter(|input| {
            let parallel = input.tangent.dot(tangent).abs() > 0.999;
            if !parallel {
                return false;
            }
            let within_extent = midpoint.distance_to(&input.midpoint) <= input.length * 0.5 + band;
            if straight && input.straight {
                line_distance(midpoint, &input.midpoint, &input.tangent) < band && within_extent
            } else {
                midpoint.distance_to(&input.midpoint) < band
            }
        })
        .min_by(|a, b| {
            midpoint
                .distance_to(&a.midpoint)
                .partial_cmp(&midpoint.distance_to(&b.midpoint))
                .unwrap_or(std::cmp::Ordering::Equal)
        })
}

/// Name every edge and vertex of a served body (fillet, chamfer, or
/// — through `popkern-shell` — shell/offset): survivors `Modified`
/// with their ancestor, the rest `Lateral` / `Junction`.
pub fn assign_output_edge_vertex_names(
    brep: &mut BRep,
    ctx: &mut TrackingContext,
    pre: &PreOperationSnapshot,
    tol: &Tolerance,
) {
    let band = tol.linear() * 100.0;
    let (mut mod_e, mut gen_e, mut mod_v, mut gen_v) = (0u32, 0u32, 0u32, 0u32);
    let edge_handles: Vec<_> = brep.edges.keys().collect();
    for eh in edge_handles {
        let Some(criteria) = capture_edge(brep, eh) else {
            continue;
        };
        let straight = matches!(criteria.curve_type, CurveType::Line);
        let matched = matching_input(
            &pre.edges,
            &criteria.midpoint,
            &criteria.tangent,
            straight,
            band,
        )
        .map(|input| input.name.clone());
        let new_name = match matched {
            Some(input_name) => {
                let ancestor = brep
                    .name_pool
                    .find(&input_name)
                    .unwrap_or_else(|| brep.name_pool.intern(input_name.clone(), None));
                let name = EntityName {
                    operation: ctx.operation,
                    role: EntityRole {
                        category: RoleCategory::Modified,
                        index: mod_e,
                    },
                    ancestors: vec![ancestor],
                };
                mod_e += 1;
                ctx.record_modified(input_name, name.clone());
                name
            }
            None => {
                let name = EntityName {
                    operation: ctx.operation,
                    role: EntityRole {
                        category: RoleCategory::Lateral,
                        index: gen_e,
                    },
                    ancestors: Vec::new(),
                };
                gen_e += 1;
                ctx.record_generated(name.clone());
                name
            }
        };
        let id = brep
            .name_pool
            .intern(new_name, Some(GeometricCriteria::Edge(criteria)));
        brep.edges[eh].name = Some(id);
    }
    let vertex_handles: Vec<_> = brep.vertices.keys().collect();
    for vh in vertex_handles {
        let Some(criteria) = capture_vertex(brep, vh) else {
            continue;
        };
        let matched = pre
            .vertices
            .iter()
            .filter(|input| criteria.position.distance_to(&input.position) < band)
            .min_by(|a, b| {
                criteria
                    .position
                    .distance_to(&a.position)
                    .partial_cmp(&criteria.position.distance_to(&b.position))
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|input| input.name.clone());
        let new_name = match matched {
            Some(input_name) => {
                let ancestor = brep
                    .name_pool
                    .find(&input_name)
                    .unwrap_or_else(|| brep.name_pool.intern(input_name.clone(), None));
                let name = EntityName {
                    operation: ctx.operation,
                    role: EntityRole {
                        category: RoleCategory::Modified,
                        index: mod_v,
                    },
                    ancestors: vec![ancestor],
                };
                mod_v += 1;
                ctx.record_modified(input_name, name.clone());
                name
            }
            None => {
                let name = EntityName {
                    operation: ctx.operation,
                    role: EntityRole {
                        category: RoleCategory::Junction,
                        index: gen_v,
                    },
                    ancestors: Vec::new(),
                };
                gen_v += 1;
                ctx.record_generated(name.clone());
                name
            }
        };
        let id = brep
            .name_pool
            .intern(new_name, Some(GeometricCriteria::Vertex(criteria)));
        brep.vertices[vh].name = Some(id);
    }
}
