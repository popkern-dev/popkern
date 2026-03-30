//! Name resolution — map persistent names back to live topology handles.
//!
//! After a feature tree is re-executed, entities may have moved, split, or
//! disappeared. [`resolve_name`] looks up a name structurally (via the
//! [`NamePool`]) then falls back to geometric matching when needed.
//! [`resolve_with_criteria`] performs geometry-only matching, useful when
//! the pool has been rebuilt and structural ids are gone.

use crate::naming::{EntityName, NameId};
use crate::naming_criteria::*;
use crate::tolerance::Tolerance;
use crate::topology::*;

// ---------------------------------------------------------------------------
// Resolution types
// ---------------------------------------------------------------------------

/// Result of resolving a persistent name against a BRep.
#[derive(Debug)]
pub enum Resolution {
    /// Exactly one entity matched.
    Exact(ResolvedEntity),
    /// The entity was split into multiple successors.
    Split(Vec<ResolvedEntity>),
    /// No matching entity was found.
    NotFound,
}

/// A resolved topological entity handle.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ResolvedEntity {
    /// A matched vertex.
    Vertex(VertexHandle),
    /// A matched edge.
    Edge(EdgeHandle),
    /// A matched face.
    Face(FaceHandle),
}

// ---------------------------------------------------------------------------
// Structural resolution
// ---------------------------------------------------------------------------

/// Resolve a persistent name against a BRep.
///
/// Algorithm:
/// 1. Find the name in `brep.name_pool` to get a [`NameId`].
/// 2. Scan all vertices, edges, and faces for entities whose `name` field
///    equals that id.
/// 3. If 0 matches and the pool has geometric criteria for the name, fall back
///    to [`resolve_with_criteria`].
/// 4. If 1 match, return [`Resolution::Exact`].
/// 5. If 2+ matches, return [`Resolution::Split`].
pub fn resolve_name(name: &EntityName, brep: &BRep) -> Resolution {
    let Some(name_id) = brep.name_pool.find(name) else {
        return Resolution::NotFound;
    };

    let matches = collect_matches(name_id, brep);

    match matches.len() {
        0 => {
            // Try geometric fallback if criteria are available.
            if let Some(criteria) = brep.name_pool.get_criteria(name_id) {
                resolve_with_criteria(criteria, brep, &brep.tolerance)
            } else {
                Resolution::NotFound
            }
        }
        1 => Resolution::Exact(matches[0]),
        _ => Resolution::Split(matches),
    }
}

/// Collect all entities in a BRep whose `name` field matches `name_id`.
fn collect_matches(name_id: NameId, brep: &BRep) -> Vec<ResolvedEntity> {
    let mut matches = Vec::new();
    for (handle, vertex) in &brep.vertices {
        if vertex.name == Some(name_id) {
            matches.push(ResolvedEntity::Vertex(handle));
        }
    }
    for (handle, edge) in &brep.edges {
        if edge.name == Some(name_id) {
            matches.push(ResolvedEntity::Edge(handle));
        }
    }
    for (handle, face) in &brep.faces {
        if face.name == Some(name_id) {
            matches.push(ResolvedEntity::Face(handle));
        }
    }
    matches
}

// ---------------------------------------------------------------------------
// Geometric fallback resolution
// ---------------------------------------------------------------------------

/// Resolve by geometric criteria alone (no pool lookup).
///
/// Iterates all entities of the appropriate kind and picks the best match.
/// For faces: surface type, centroid, normal, and area must all match within
/// tolerance. For edges: curve type, midpoint, tangent, and length. For
/// vertices: position only.
///
/// Returns [`Resolution::Exact`] for a single best match,
/// [`Resolution::Split`] if multiple entities match equally well,
/// or [`Resolution::NotFound`] if nothing qualifies.
pub fn resolve_with_criteria(
    criteria: &GeometricCriteria,
    brep: &BRep,
    tolerance: &Tolerance,
) -> Resolution {
    match criteria {
        GeometricCriteria::Face(fc) => resolve_face_by_criteria(fc, brep, tolerance),
        GeometricCriteria::Edge(ec) => resolve_edge_by_criteria(ec, brep, tolerance),
        GeometricCriteria::Vertex(vc) => resolve_vertex_by_criteria(vc, brep, tolerance),
    }
}

/// Resolve a face by geometric criteria.
fn resolve_face_by_criteria(
    criteria: &FaceCriteria,
    brep: &BRep,
    tolerance: &Tolerance,
) -> Resolution {
    let generous_linear = tolerance.linear * 10.0;
    let mut matches = Vec::new();

    for (handle, _face) in &brep.faces {
        let Some(fc) = capture_face(brep, handle) else {
            continue;
        };
        // Surface type must match exactly.
        if fc.surface_type != criteria.surface_type {
            continue;
        }
        // Centroid within generous linear tolerance.
        if fc.centroid.distance_to(&criteria.centroid) > generous_linear {
            continue;
        }
        // Normal within angular tolerance (compare via dot product).
        let dot = fc.normal.dot(&criteria.normal);
        // dot ≈ 1.0 means same direction; allow tolerance.angular deviation.
        // cos(theta) >= cos(angular) ≈ 1 - angular^2/2 for small angles.
        if dot < 1.0 - tolerance.angular {
            continue;
        }
        // Area within 20%.
        if criteria.area > 1e-12 {
            let ratio = fc.area / criteria.area;
            if !(0.8..=1.2).contains(&ratio) {
                continue;
            }
        }
        matches.push(ResolvedEntity::Face(handle));
    }

    to_resolution(matches)
}

/// Resolve an edge by geometric criteria.
fn resolve_edge_by_criteria(
    criteria: &EdgeCriteria,
    brep: &BRep,
    tolerance: &Tolerance,
) -> Resolution {
    let generous_linear = tolerance.linear * 10.0;
    let mut matches = Vec::new();

    for (handle, _edge) in &brep.edges {
        let Some(ec) = capture_edge(brep, handle) else {
            continue;
        };
        // Curve type must match exactly.
        if ec.curve_type != criteria.curve_type {
            continue;
        }
        // Midpoint within generous linear tolerance.
        if ec.midpoint.distance_to(&criteria.midpoint) > generous_linear {
            continue;
        }
        // Tangent within angular tolerance (compare absolute dot for direction).
        let dot = ec.tangent.dot(&criteria.tangent).abs();
        if dot < 1.0 - tolerance.angular {
            continue;
        }
        // Length within 20%.
        if criteria.length > 1e-12 {
            let ratio = ec.length / criteria.length;
            if !(0.8..=1.2).contains(&ratio) {
                continue;
            }
        }
        matches.push(ResolvedEntity::Edge(handle));
    }

    to_resolution(matches)
}

/// Resolve a vertex by geometric criteria.
fn resolve_vertex_by_criteria(
    criteria: &VertexCriteria,
    brep: &BRep,
    tolerance: &Tolerance,
) -> Resolution {
    let mut matches = Vec::new();

    for (handle, vertex) in &brep.vertices {
        if vertex.position.distance_to(&criteria.position) <= tolerance.linear {
            matches.push(ResolvedEntity::Vertex(handle));
        }
    }

    to_resolution(matches)
}

/// Convert a match list into a Resolution.
fn to_resolution(matches: Vec<ResolvedEntity>) -> Resolution {
    match matches.len() {
        0 => Resolution::NotFound,
        1 => Resolution::Exact(matches[0]),
        _ => Resolution::Split(matches),
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::naming::{EntityName, EntityRole, NamePool, OperationId, RoleCategory};
    use crate::tolerance::Tolerance;

    fn make_name(op: u64, cat: RoleCategory, idx: u32) -> EntityName {
        EntityName {
            operation: OperationId(op),
            role: EntityRole {
                category: cat,
                index: idx,
            },
            ancestors: Vec::new(),
        }
    }

    #[test]
    fn resolve_named_vertex_exact() {
        let mut brep = BRep::unit_cube(Tolerance::default()).expect("cube");
        let vh = brep.vertices.keys().next().expect("vertex");
        let name = make_name(1, RoleCategory::Generated, 0);
        let criteria = capture_vertex(&brep, vh).map(GeometricCriteria::Vertex);
        let nid = brep.name_pool.intern(name.clone(), criteria);
        brep.vertices[vh].name = Some(nid);

        let resolution = resolve_name(&name, &brep);
        assert!(matches!(resolution, Resolution::Exact(ResolvedEntity::Vertex(h)) if h == vh));
    }

    #[test]
    fn resolve_named_face_exact() {
        let mut brep = BRep::unit_cube(Tolerance::default()).expect("cube");
        let fh = brep.faces.keys().next().expect("face");
        let name = make_name(1, RoleCategory::Generated, 1);
        let criteria = capture_face(&brep, fh).map(GeometricCriteria::Face);
        let nid = brep.name_pool.intern(name.clone(), criteria);
        brep.faces[fh].name = Some(nid);

        let resolution = resolve_name(&name, &brep);
        assert!(matches!(resolution, Resolution::Exact(ResolvedEntity::Face(h)) if h == fh));
    }

    #[test]
    fn resolve_unknown_name_not_found() {
        let brep = BRep::unit_cube(Tolerance::default()).expect("cube");
        let name = make_name(99, RoleCategory::Cap, 0);
        let resolution = resolve_name(&name, &brep);
        assert!(matches!(resolution, Resolution::NotFound));
    }

    #[test]
    fn resolve_split_two_entities_same_name() {
        let mut brep = BRep::unit_cube(Tolerance::default()).expect("cube");
        let name = make_name(1, RoleCategory::Generated, 0);
        let nid = brep.name_pool.intern(name.clone(), None);

        // Assign same name to two faces (simulates a split).
        let mut face_handles: Vec<FaceHandle> = brep.faces.keys().collect();
        face_handles.truncate(2);
        brep.faces[face_handles[0]].name = Some(nid);
        brep.faces[face_handles[1]].name = Some(nid);

        let resolution = resolve_name(&name, &brep);
        match resolution {
            Resolution::Split(entities) => {
                assert_eq!(entities.len(), 2);
            }
            other => panic!("expected Split, got {:?}", other),
        }
    }

    #[test]
    fn resolve_face_geometric_fallback() {
        let mut brep = BRep::unit_cube(Tolerance::default()).expect("cube");
        let fh = brep.faces.keys().next().expect("face");

        // Capture criteria, then name the face.
        let fc = capture_face(&brep, fh).expect("face criteria");
        let name = make_name(1, RoleCategory::Generated, 0);
        let nid = brep
            .name_pool
            .intern(name.clone(), Some(GeometricCriteria::Face(fc.clone())));
        brep.faces[fh].name = Some(nid);

        // Now clear the name from the entity (simulate rebuild).
        brep.faces[fh].name = None;

        // Structural match finds nothing, falls back to geometry.
        let resolution = resolve_name(&name, &brep);
        assert!(matches!(
            resolution,
            Resolution::Exact(ResolvedEntity::Face(_))
        ));
    }

    #[test]
    fn resolve_vertex_by_criteria_finds_by_position() {
        let brep = BRep::unit_cube(Tolerance::default()).expect("cube");
        let vh = brep.vertices.keys().next().expect("vertex");
        let vc = capture_vertex(&brep, vh).expect("vertex criteria");

        let resolution =
            resolve_with_criteria(&GeometricCriteria::Vertex(vc), &brep, &brep.tolerance);
        assert!(matches!(
            resolution,
            Resolution::Exact(ResolvedEntity::Vertex(_))
        ));
    }

    #[test]
    fn resolve_edge_by_criteria_finds_edge() {
        let brep = BRep::unit_cube(Tolerance::default()).expect("cube");
        let eh = brep
            .edges
            .iter()
            .find(|(_, e)| e.curve.is_some())
            .map(|(h, _)| h)
            .expect("edge with curve");
        let ec = capture_edge(&brep, eh).expect("edge criteria");

        let resolution =
            resolve_with_criteria(&GeometricCriteria::Edge(ec), &brep, &brep.tolerance);
        // Should find at least one matching edge (several cube edges
        // share the same geometric profile, so Split is acceptable too).
        assert!(!matches!(resolution, Resolution::NotFound));
    }
}
