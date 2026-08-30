//! Stable names for boolean output EDGES and VERTICES (PK-386).
//!
//! `entity_naming` matches output faces back to operand faces. Edges
//! and vertices were left unnamed, so no dress-up after a boolean
//! could be recorded by name in a document (D3 selection snapshots;
//! measured: plate ∪ boss → 9/9 faces named, 0/18 edges). Same rule
//! as faces: an output entity that coincides with an operand entity
//! (edge: midpoint within the match band and parallel tangent; vertex:
//! position within the band) is `Modified` with that ancestor; the
//! rest — intersection curves and their endpoints — are `Intersection`
//! edges / `Junction` vertices of this operation.

use popkern_core::geometry::{Point3, Vec3};
use popkern_core::history::TrackingContext;
use popkern_core::naming::{EntityName, EntityRole, RoleCategory};
use popkern_core::naming_criteria::{GeometricCriteria, capture_edge, capture_vertex};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

/// An operand edge: its name (or a placeholder for unnamed input),
/// midpoint and tangent.
pub(crate) type EdgeSnapshot = (EntityName, Point3, Vec3);
/// An operand vertex: its name (or a placeholder) and position.
pub(crate) type VertexSnapshot = (EntityName, Point3);

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

/// Snapshot every edge of an operand before the boolean consumes it.
pub(crate) fn collect_edge_data(brep: &BRep) -> Vec<EdgeSnapshot> {
    let mut out = Vec::new();
    for (i, (eh, edge)) in brep.edges.iter().enumerate() {
        let Some(ec) = capture_edge(brep, eh) else {
            continue;
        };
        let name = edge
            .name
            .and_then(|nid| brep.name_pool.get(nid).cloned())
            .unwrap_or_else(|| placeholder(RoleCategory::Lateral, i));
        out.push((name, ec.midpoint, ec.tangent));
    }
    out
}

/// Snapshot every vertex of an operand.
pub(crate) fn collect_vertex_data(brep: &BRep) -> Vec<VertexSnapshot> {
    let mut out = Vec::new();
    for (i, (vh, vertex)) in brep.vertices.iter().enumerate() {
        let Some(vc) = capture_vertex(brep, vh) else {
            continue;
        };
        let name = vertex
            .name
            .and_then(|nid| brep.name_pool.get(nid).cloned())
            .unwrap_or_else(|| placeholder(RoleCategory::Junction, i));
        out.push((name, vc.position));
    }
    out
}

/// Name every output edge and vertex of `brep`.
pub(crate) fn assign_boolean_edge_vertex_names(
    brep: &mut BRep,
    input_edges: &[EdgeSnapshot],
    input_vertices: &[VertexSnapshot],
    ctx: &mut TrackingContext,
    tol: &Tolerance,
) {
    let band = tol.linear() * 100.0;
    let (mut mod_e, mut gen_e, mut mod_v, mut gen_v) = (0u32, 0u32, 0u32, 0u32);

    let edge_handles: Vec<_> = brep.edges.keys().collect();
    for eh in edge_handles {
        let Some(ec) = capture_edge(brep, eh) else {
            continue;
        };
        let best = input_edges
            .iter()
            .filter(|(_, mid, tan)| {
                ec.midpoint.distance_to(mid) < band && ec.tangent.dot(tan).abs() > 0.9
            })
            .min_by(|a, b| {
                ec.midpoint
                    .distance_to(&a.1)
                    .partial_cmp(&ec.midpoint.distance_to(&b.1))
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
        let new_name = match best {
            Some((input_name, _, _)) => {
                let ancestor = brep.name_pool.intern(input_name.clone(), None);
                let name = EntityName {
                    operation: ctx.operation,
                    role: EntityRole {
                        category: RoleCategory::Modified,
                        index: mod_e,
                    },
                    ancestors: vec![ancestor],
                };
                mod_e += 1;
                ctx.record_modified(input_name.clone(), name.clone());
                name
            }
            None => {
                let name = EntityName {
                    operation: ctx.operation,
                    role: EntityRole {
                        category: RoleCategory::Intersection,
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
            .intern(new_name, Some(GeometricCriteria::Edge(ec)));
        brep.edges[eh].name = Some(id);
    }

    let vertex_handles: Vec<_> = brep.vertices.keys().collect();
    for vh in vertex_handles {
        let Some(vc) = capture_vertex(brep, vh) else {
            continue;
        };
        let best = input_vertices
            .iter()
            .filter(|(_, p)| vc.position.distance_to(p) < band)
            .min_by(|a, b| {
                vc.position
                    .distance_to(&a.1)
                    .partial_cmp(&vc.position.distance_to(&b.1))
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
        let new_name = match best {
            Some((input_name, _)) => {
                let ancestor = brep.name_pool.intern(input_name.clone(), None);
                let name = EntityName {
                    operation: ctx.operation,
                    role: EntityRole {
                        category: RoleCategory::Modified,
                        index: mod_v,
                    },
                    ancestors: vec![ancestor],
                };
                mod_v += 1;
                ctx.record_modified(input_name.clone(), name.clone());
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
            .intern(new_name, Some(GeometricCriteria::Vertex(vc)));
        brep.vertices[vh].name = Some(id);
    }
}
