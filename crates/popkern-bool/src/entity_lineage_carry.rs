//! Carried entities keep their names and skip lineage (naming scheme 4,
//! PK-803).
//!
//! Through scheme 3 a boolean minted a `Modified` name for EVERY result
//! entity (its input as ancestor) and recorded, by value, the whole
//! ancestry of the input it lies on: per boolean the pool grew by the
//! body's entity count and every entity's lineage was re-flattened —
//! the cost of the 64th cut of a slotted disc doubled with naming on.
//! Scheme 4 treats an entity the boolean CARRIES — operand `a`'s entity
//! unchanged (an edge's end and mid points, a face's loop structure and
//! vertex box, a vertex's position) and clear of the tool's box — as
//! that same entity: it keeps `a`'s name id (the scheme-3+ result pool
//! starts as `a`'s) and records neither a `Modified` event nor lineage.
//! Clear of the tool, it lies on no input entity but itself (two faces
//! or edges of a valid body share no interior point; every tool entity
//! sits in the tool's box), so the lineage it skips is its own, already
//! in the pool. A label or stable name bound to it resolves to it
//! directly, as an unmodified shape does in commercial kernels' history.

use std::collections::HashMap;

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::naming::NameId;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

/// Which entities a role index counts over.
#[derive(Debug, Clone, Copy)]
pub(crate) enum EntityKind {
    /// Faces.
    Face,
    /// Edges.
    Edge,
    /// Vertices.
    Vertex,
}

/// Operand `a`'s named entities by name id, and the tool's box.
pub(crate) struct CarryIndex {
    edges: HashMap<NameId, EdgeHandle>,
    faces: HashMap<NameId, FaceHandle>,
    vertices: HashMap<NameId, popkern_core::topology::VertexHandle>,
    face_boxes: HashMap<FaceHandle, Aabb>,
    a_face_boxes: Vec<(FaceHandle, Aabb)>,
    tool: Option<Aabb>,
    on: f64,
}

impl CarryIndex {
    /// Index `a`'s names; the tool box is `b`'s face boxes' union grown
    /// by `on`.
    pub(crate) fn new(a: &BRep, b: &BRep, on: f64) -> Self {
        let edges = a
            .edges
            .iter()
            .filter_map(|(h, e)| Some((e.name?, h)))
            .collect();
        let faces = a
            .faces
            .iter()
            .filter_map(|(h, f)| Some((f.name?, h)))
            .collect();
        let vertices = a
            .vertices
            .iter()
            .filter_map(|(h, v)| Some((v.name?, h)))
            .collect();
        let a_face_boxes = crate::intersect::compute_face_aabbs(a);
        let face_boxes = a_face_boxes.iter().cloned().collect();
        let corners: Vec<Point3> = crate::intersect::compute_face_aabbs(b)
            .iter()
            .flat_map(|(_, bb)| [bb.min, bb.max])
            .collect();
        Self {
            edges,
            faces,
            vertices,
            face_boxes,
            a_face_boxes,
            tool: Aabb::from_points(&corners),
            on,
        }
    }

    /// `box_` is clear of the tool (grown by `on`).
    fn clear_of_tool(&self, box_: &Aabb) -> bool {
        self.tool
            .as_ref()
            .is_some_and(|tool| !tool.overlaps(box_, 2.0 * self.on))
    }

    /// The name `result` gives `name` (an operand-`a` name value) is
    /// one of `a`'s, with that value — the entity keeps `a`'s name.
    pub(crate) fn keeps_a_name(result: &BRep, a: &BRep, id: NameId) -> bool {
        result.name_pool.get(id).is_some() && result.name_pool.get(id) == a.name_pool.get(id)
    }

    /// The first role index of `category` free under `operation` among
    /// `a`'s live `kind` entities. A carried entity keeps its name, so a
    /// later boolean under the SAME operation (a pattern's instances)
    /// must mint past it — restarting at 0 duplicated a kept face's
    /// spelling and the pattern's instances failed to replay.
    pub(crate) fn first_free_index(
        a: &BRep,
        kind: EntityKind,
        operation: popkern_core::naming::OperationId,
        category: popkern_core::naming::RoleCategory,
    ) -> u32 {
        let ids: Box<dyn Iterator<Item = Option<NameId>>> = match kind {
            EntityKind::Face => Box::new(a.faces.values().map(|f| f.name)),
            EntityKind::Edge => Box::new(a.edges.values().map(|e| e.name)),
            EntityKind::Vertex => Box::new(a.vertices.values().map(|v| v.name)),
        };
        ids.flatten()
            .filter_map(|id| a.name_pool.get(id))
            .filter(|n| n.operation == operation && n.role.category == category)
            .map(|n| n.role.index + 1)
            .max()
            .unwrap_or(0)
    }

    /// `a`'s face boxes, in `a`'s face order (`compute_face_aabbs`).
    pub(crate) fn a_face_boxes(&self) -> &[(FaceHandle, Aabb)] {
        &self.a_face_boxes
    }

    /// `a`'s edge named `id`, if any.
    pub(crate) fn a_edge(&self, id: NameId) -> Option<EdgeHandle> {
        self.edges.get(&id).copied()
    }

    /// `a`'s face named `id`, if any.
    pub(crate) fn a_face(&self, id: NameId) -> Option<FaceHandle> {
        self.faces.get(&id).copied()
    }

    /// `a`'s vertex named `id`, if any.
    pub(crate) fn a_vertex(&self, id: NameId) -> Option<popkern_core::topology::VertexHandle> {
        self.vertices.get(&id).copied()
    }

    /// `result`'s vertex `vertex` is `a`'s vertex `av`, clear of the tool.
    pub(crate) fn vertex_is_carried(
        &self,
        result: &BRep,
        vertex: popkern_core::topology::VertexHandle,
        a: &BRep,
        av: popkern_core::topology::VertexHandle,
    ) -> bool {
        let (Some(rv), Some(avd)) = (result.vertices.get(vertex), a.vertices.get(av)) else {
            return false;
        };
        rv.position.distance_to(&avd.position) <= self.on
            && self.clear_of_tool(&Aabb {
                min: avd.position,
                max: avd.position,
            })
    }

    /// `result`'s edge `edge` is `a`'s edge `ae`, carried clear of the
    /// tool (see the module docs).
    pub(crate) fn edge_is_carried(
        &self,
        result: &BRep,
        edge: EdgeHandle,
        a: &BRep,
        ae: EdgeHandle,
    ) -> bool {
        let (Some(re), Some(ae)) = (result.edges.get(edge), a.edges.get(ae)) else {
            return false;
        };
        let (Some(rc), Some(ac)) = (re.curve.as_ref(), ae.curve.as_ref()) else {
            return false;
        };
        let ends = |c: &popkern_core::geometry::curve::CurveOnEdge| {
            let mid = 0.5 * (c.t_start + c.t_end);
            (
                c.geometry.point_at(c.t_start),
                c.geometry.point_at(c.t_end),
                c.geometry.point_at(mid),
            )
        };
        let ((r0, r1, rm), (a0, a1, am)) = (ends(rc), ends(ac));
        let near = |p: &Point3, q: &Point3| p.distance_to(q) <= self.on;
        let same = near(&rm, &am)
            && ((near(&r0, &a0) && near(&r1, &a1)) || (near(&r0, &a1) && near(&r1, &a0)));
        let lo = ac.t_start.min(ac.t_end);
        let hi = ac.t_start.max(ac.t_end);
        same && self.clear_of_tool(&ac.geometry.curve_bounding_box((lo, hi)))
    }

    /// `result`'s face `face` is `a`'s face `af`, carried clear of the
    /// tool (judged on `a`'s conservative face box).
    pub(crate) fn face_is_carried(
        &self,
        result: &BRep,
        face: FaceHandle,
        a: &BRep,
        af: FaceHandle,
    ) -> bool {
        let shape = |brep: &BRep, f: FaceHandle| -> Option<(Vec<usize>, Aabb)> {
            let data = brep.faces.get(f)?;
            let mut lens = Vec::new();
            let mut points = Vec::new();
            for lh in std::iter::once(data.outer_loop).chain(data.inner_loops.iter().copied()) {
                let start = brep.loops.get(lh)?.half_edge;
                let (mut he, mut n) = (start, 0usize);
                loop {
                    let h = brep.half_edges.get(he)?;
                    points.push(brep.vertices.get(h.origin)?.position);
                    n += 1;
                    he = h.next;
                    if he == start || n > brep.half_edges.len() {
                        break;
                    }
                }
                lens.push(n);
            }
            // Inner-loop order is not identity (the walk lists holes in
            // its own order).
            lens[1..].sort_unstable();
            Some((lens, Aabb::from_points(&points)?))
        };
        let (Some((rl, rb)), Some((al, ab))) = (shape(result, face), shape(a, af)) else {
            return false;
        };
        let near = |p: &Point3, q: &Point3| p.distance_to(q) <= self.on;
        rl == al
            && near(&rb.min, &ab.min)
            && near(&rb.max, &ab.max)
            && self
                .face_boxes
                .get(&af)
                .is_some_and(|bb| self.clear_of_tool(bb))
    }
}
