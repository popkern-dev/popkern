//! Naming scheme 5 (PK-1121): a boolean result's entities named by the
//! operand entities they are pieces of and by their neighbours, never
//! by position.
//!
//! Scheme 4 matched an output face to an operand face by centroid and
//! an output edge by midpoint, within a hundred tolerances: a face the
//! boolean trimmed, or an edge it cut, moved its centroid or midpoint
//! and was minted a `Generated` or `Intersection` name numbered in slot
//! order. When an upstream edit moved a tool, the same face came out
//! with another number, and a fillet recorded on its edges no longer
//! resolved. Here:
//!
//! - a carried entity (unchanged and clear of every tool) keeps the
//!   target's name, as in scheme 4;
//! - a face is `Modified` with the operand faces it is a piece of as
//!   ancestors (`entity_sources`), `Generated` when it is a piece of
//!   none;
//! - an edge is `Modified` with the operand edges it is a piece of as
//!   ancestors, else `Intersection` with its faces' names;
//! - a vertex is `Modified` with the operand vertex it lies on, else
//!   `Junction` with its edges' names.
//!
//! A piece of an operand entity is named by that entity alone, as a
//! face is, so an upstream edit renaming a face it lies on (a chamfer
//! inserted on the target) leaves its ancestry whole and a reference
//! to it rebinds by provenance; its neighbours tell its siblings apart.
//!
//! Entities sharing a role and ancestors (pieces of one face with the
//! same faces around) are told apart by their index: their rank by the
//! sense against their source (of two coincident faces one faces each
//! way), then the names their neighbours are pieces of, and only where
//! those tie too,
//! by position (the lowest x, then y, then z of a sample point), which
//! is the one place a moving tool can still renumber.

use std::collections::{BTreeMap, BTreeSet, HashMap};

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::history::TrackingContext;
use popkern_core::naming::{EntityName, EntityRole, NameId, RoleCategory};
use popkern_core::naming_criteria::{
    GeometricCriteria, capture_edge, capture_face, capture_vertex,
};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

use crate::entity_lineage_carry::CarryIndex;
use crate::entity_naming_edges::{EdgeSnapshot, VertexSnapshot};
use crate::entity_naming_index::NearIndex;
use crate::entity_sources::{EdgeSources, FaceSources, Source, lineage_of};

/// The operands a result is named after: the target `a` with its carry
/// index, the tools, and their snapshots.
pub(crate) struct Operands<'a> {
    /// The target operand.
    pub a: &'a BRep,
    /// The tool operands.
    pub tools: &'a [&'a BRep],
    /// Every operand face, `a`'s first: name, centroid, normal.
    pub faces: Vec<&'a (EntityName, Point3, popkern_core::geometry::Vec3)>,
    /// Every operand edge.
    pub edges: &'a [EdgeSnapshot],
    /// Every operand vertex.
    pub vertices: &'a [VertexSnapshot],
    /// What the target carries unchanged (scheme 4 and later).
    pub carry: Option<&'a CarryIndex>,
}

/// One result entity awaiting its name: role, ancestors, its tie key
/// and the operand name it is recorded as modified from.
struct Pending<H> {
    handle: H,
    category: RoleCategory,
    ancestors: Vec<NameId>,
    signature: Vec<u32>,
    at: Point3,
    from: Option<EntityName>,
    /// The lineage it records (the sources' names and ancestry).
    lineage: Vec<EntityName>,
}

/// Name every face, edge and vertex of `result` under scheme 5, with
/// each named face's and edge's lineage (what it is a piece of, the
/// relation its name was read from).
pub(crate) fn assign_topological_names(
    result: &mut BRep,
    ops: &Operands<'_>,
    ctx: &mut TrackingContext,
    tol: &Tolerance,
) {
    let band = tol.linear() * 100.0;
    // Faces.
    let sources = FaceSources::new(ops.a, ops.tools, tol, ops.carry);
    let mut piece_of: HashMap<FaceHandle, Vec<NameId>> = HashMap::new();
    let mut pending = Vec::new();
    let faces: Vec<FaceHandle> = result.faces.keys().collect();
    for f in faces {
        let Some(fc) = capture_face(result, f) else {
            continue;
        };
        if let Some(id) = carried_face(result, f, fc.centroid, fc.normal, ops, band) {
            result.faces[f].name = Some(id);
            piece_of.insert(f, vec![id]);
            continue;
        }
        let sensed = sources.sensed(result, f);
        let opposed = sensed.first().is_some_and(|(_, aligned)| !aligned);
        let from: Vec<Source<'_>> = sensed.into_iter().map(|(s, _)| s).collect();
        let ancestors = imported(result, &from);
        piece_of.insert(f, ancestors.clone());
        pending.push(Pending {
            handle: f,
            category: if ancestors.is_empty() {
                RoleCategory::Generated
            } else {
                RoleCategory::Modified
            },
            ancestors,
            // Its sense first: of two coincident faces, pieces of one
            // operand face, one faces the way it does and one opposite.
            signature: vec![u32::from(opposed)],
            at: fc.centroid,
            from: first_name(&from),
            lineage: lineage_of(&from),
        });
    }
    let neighbours = face_neighbours(result);
    for p in &mut pending {
        let around: BTreeSet<u32> = neighbours
            .get(&p.handle)
            .into_iter()
            .flatten()
            .flat_map(|g| piece_of.get(g).into_iter().flatten().map(|id| id.0))
            .collect();
        p.signature.extend(around);
    }
    for (f, name, lineage) in ranked(pending, ctx) {
        let criteria = capture_face(result, f).map(GeometricCriteria::Face);
        let id = result.name_pool.intern(name, criteria);
        result.faces[f].name = Some(id);
        result.name_pool.add_lineage(id, lineage);
    }

    // Edges.
    let sources = EdgeSources::new(ops.a, ops.tools, tol);
    let edge_index = NearIndex::new(ops.edges.iter().map(|e| e.1), band);
    let mut pending = Vec::new();
    let edges: Vec<EdgeHandle> = result.edges.keys().collect();
    for e in edges {
        let Some((mid, tangent)) = mid_tangent(result, e) else {
            continue;
        };
        if let Some(id) = carried_edge(result, e, mid, tangent, ops, &edge_index, band) {
            result.edges[e].name = Some(id);
            continue;
        }
        let from = sources.of(result, e);
        let mut ancestors = imported(result, &from);
        let faces: BTreeSet<u32> = result
            .radial_fins(e)
            .into_iter()
            .filter_map(|h| result.fin_face(h))
            .filter_map(|f| result.faces.get(f).and_then(|f| f.name))
            .map(|id| id.0)
            .collect();
        // A piece of an operand edge is named by that edge, its faces
        // telling its siblings apart: an upstream edit that renames a
        // face it lies on leaves its ancestry whole. An intersection
        // edge has only its faces to be named by.
        let (category, signature) = if ancestors.is_empty() {
            ancestors.extend(faces.into_iter().map(NameId));
            (RoleCategory::Intersection, Vec::new())
        } else {
            (RoleCategory::Modified, faces.into_iter().collect())
        };
        pending.push(Pending {
            handle: e,
            category,
            ancestors,
            signature,
            at: mid,
            from: first_name(&from),
            lineage: lineage_of(&from),
        });
    }
    for (e, name, lineage) in ranked(pending, ctx) {
        let criteria = capture_edge(result, e).map(GeometricCriteria::Edge);
        let id = result.name_pool.intern(name, criteria);
        result.edges[e].name = Some(id);
        result.name_pool.add_lineage(id, lineage);
    }

    // Vertices.
    let vertex_index = NearIndex::new(ops.vertices.iter().map(|v| v.1), band);
    let mut incident: HashMap<VertexHandle, BTreeSet<u32>> = HashMap::new();
    for he in result.half_edges.values() {
        if let Some(id) = result.edges.get(he.edge).and_then(|e| e.name) {
            incident.entry(he.origin).or_default().insert(id.0);
        }
    }
    let mut pending = Vec::new();
    let vertices: Vec<VertexHandle> = result.vertices.keys().collect();
    for v in vertices {
        let Some(vc) = capture_vertex(result, v) else {
            continue;
        };
        let near = vertex_index
            .near(&vc.position)
            .into_iter()
            .map(|i| &ops.vertices[i])
            .filter(|(_, p)| vc.position.distance_to(p) < band)
            .min_by(|x, y| {
                vc.position
                    .distance_to(&x.1)
                    .total_cmp(&vc.position.distance_to(&y.1))
            });
        let mut ancestors = Vec::new();
        let mut from = None;
        if let Some((input, _)) = near {
            let id = result.name_pool.intern(input.clone(), None);
            if let Some(c) = ops.carry
                && CarryIndex::keeps_a_name(result, ops.a, id)
                && c.a_vertex(id)
                    .is_some_and(|av| c.vertex_is_carried(result, v, ops.a, av))
            {
                result.vertices[v].name = Some(id);
                continue;
            }
            ancestors.push(id);
            from = Some(input.clone());
        }
        let edges = incident.get(&v).into_iter().flatten().copied();
        let (category, signature) = if ancestors.is_empty() {
            ancestors.extend(edges.map(NameId));
            (RoleCategory::Junction, Vec::new())
        } else {
            (RoleCategory::Modified, edges.collect())
        };
        pending.push(Pending {
            handle: v,
            category,
            ancestors,
            signature,
            at: vc.position,
            from,
            lineage: Vec::new(),
        });
    }
    for (v, name, _) in ranked(pending, ctx) {
        let criteria = capture_vertex(result, v).map(GeometricCriteria::Vertex);
        let id = result.name_pool.intern(name, criteria);
        result.vertices[v].name = Some(id);
    }
}

/// Each pending entity's name: its index the rank among entities of
/// its role and ancestors, by neighbour signature then position; each
/// recorded in `ctx` as modified from its source or generated.
fn ranked<H: Copy>(
    pending: Vec<Pending<H>>,
    ctx: &mut TrackingContext,
) -> Vec<(H, EntityName, Vec<EntityName>)> {
    let mut groups: BTreeMap<(u8, Vec<u32>), Vec<Pending<H>>> = BTreeMap::new();
    for p in pending {
        let key = (
            category_rank(p.category),
            p.ancestors.iter().map(|a| a.0).collect(),
        );
        groups.entry(key).or_default().push(p);
    }
    let mut out = Vec::new();
    for (_, mut group) in groups {
        group.sort_by(|x, y| {
            x.signature.cmp(&y.signature).then_with(|| {
                (x.at.x.total_cmp(&y.at.x))
                    .then(x.at.y.total_cmp(&y.at.y))
                    .then(x.at.z.total_cmp(&y.at.z))
            })
        });
        for (index, p) in group.into_iter().enumerate() {
            let name = EntityName {
                operation: ctx.operation,
                role: EntityRole {
                    category: p.category,
                    index: index as u32,
                },
                ancestors: p.ancestors,
            };
            match p.from {
                Some(from) => ctx.record_modified(from, name.clone()),
                None => ctx.record_generated(name.clone()),
            }
            out.push((p.handle, name, p.lineage));
        }
    }
    out
}

fn category_rank(c: RoleCategory) -> u8 {
    match c {
        RoleCategory::Modified => 0,
        RoleCategory::Generated => 1,
        RoleCategory::Lateral => 2,
        RoleCategory::Cap => 3,
        RoleCategory::Intersection => 4,
        RoleCategory::Junction => 5,
    }
}

/// `sources`' names imported into `result`'s pool, in order, once each.
fn imported(result: &mut BRep, sources: &[Source<'_>]) -> Vec<NameId> {
    let mut out = Vec::new();
    for (input, id) in sources {
        let copied = if std::ptr::eq(&input.name_pool, &result.name_pool) {
            Some(*id)
        } else {
            result.name_pool.import(&input.name_pool, *id)
        };
        if let Some(c) = copied
            && !out.contains(&c)
        {
            out.push(c);
        }
    }
    out
}

fn first_name(sources: &[Source<'_>]) -> Option<EntityName> {
    sources
        .first()
        .and_then(|(input, id)| input.name_pool.get(*id).cloned())
}

/// Each face's neighbours across its edges.
fn face_neighbours(result: &BRep) -> HashMap<FaceHandle, BTreeSet<FaceHandle>> {
    let mut out: HashMap<FaceHandle, BTreeSet<FaceHandle>> = HashMap::new();
    for e in result.edges.keys() {
        let faces: BTreeSet<FaceHandle> = result
            .radial_fins(e)
            .into_iter()
            .filter_map(|h| result.fin_face(h))
            .collect();
        for &f in &faces {
            out.entry(f)
                .or_default()
                .extend(faces.iter().filter(|&&g| g != f));
        }
    }
    out
}

fn mid_tangent(result: &BRep, e: EdgeHandle) -> Option<(Point3, popkern_core::geometry::Vec3)> {
    let c = result.edges.get(e)?.curve.as_ref()?;
    let t = 0.5 * (c.t_start + c.t_end);
    Some((c.geometry.point_at(t), c.geometry.tangent_at(t)))
}

/// The target's name when face `f` is the target's face, unchanged and
/// clear of every tool (scheme 4's carry, matched as scheme 4 does).
fn carried_face(
    result: &mut BRep,
    f: FaceHandle,
    centroid: Point3,
    normal: popkern_core::geometry::Vec3,
    ops: &Operands<'_>,
    band: f64,
) -> Option<NameId> {
    let carry = ops.carry?;
    let (name, _, _) = ops
        .faces
        .iter()
        .filter(|(_, c, n)| centroid.distance_to(c) < band && normal.dot(n).abs() > 0.9)
        .min_by(|x, y| {
            centroid
                .distance_to(&x.1)
                .total_cmp(&centroid.distance_to(&y.1))
        })?;
    let id = result.name_pool.find(name)?;
    (CarryIndex::keeps_a_name(result, ops.a, id)
        && carry
            .a_face(id)
            .is_some_and(|af| carry.face_is_carried(result, f, ops.a, af)))
    .then_some(id)
}

/// The target's name when edge `e` is the target's edge, unchanged and
/// clear of every tool.
fn carried_edge(
    result: &BRep,
    e: EdgeHandle,
    mid: Point3,
    tangent: popkern_core::geometry::Vec3,
    ops: &Operands<'_>,
    index: &NearIndex,
    band: f64,
) -> Option<NameId> {
    let carry = ops.carry?;
    let (name, _, _) = index
        .near(&mid)
        .into_iter()
        .map(|i| &ops.edges[i])
        .filter(|(_, m, t)| mid.distance_to(m) < band && tangent.dot(t).abs() > 0.9)
        .min_by(|x, y| mid.distance_to(&x.1).total_cmp(&mid.distance_to(&y.1)))?;
    let id = result.name_pool.find(name)?;
    (CarryIndex::keeps_a_name(result, ops.a, id)
        && carry
            .a_edge(id)
            .is_some_and(|ae| carry.edge_is_carried(result, e, ops.a, ae)))
    .then_some(id)
}

#[cfg(test)]
#[path = "entity_naming_topo_tests.rs"]
mod tests;
