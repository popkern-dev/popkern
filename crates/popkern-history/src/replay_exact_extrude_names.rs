//! Stable names for an extrude served through the exact periodic wire
//! sweep — the SAME contract `BRep::extrude` mints, so a feature keeps
//! its names whichever route serves it (PK-387: the RawBody-style
//! names the exact route used to mint left `Cap 0` unresolvable for an
//! up-to-surface reference onto an analytic plate).
//!
//! Under the extrude's operation: bottom face `Modified 0` (ancestor:
//! the sketch face), cap face `Cap 0`, lateral face `i` `Lateral i`,
//! vertical edge `i` `Lateral i`, top edge `i` `Cap i`, top vertex `i`
//! `Generated i`. Bottom edges and vertices carry the SKETCH
//! operation's `Lateral i` / `Junction i`, as the sketch face's own
//! boundary does on the generic route. `i` is the profile order: the
//! wire parameter of each entity's footprint.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::history::TrackingContext;
use popkern_core::naming::{EntityName, EntityRole, NameId, NamePool, OperationId, RoleCategory};
use popkern_core::naming_criteria::{
    GeometricCriteria, capture_edge, capture_face, capture_vertex,
};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

/// Which end of the sweep an entity's footprint sits at.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Station {
    Bottom,
    Top,
    Between,
}

struct Frame<'a> {
    origin: Point3,
    normal: Vec3,
    lo: f64,
    hi: f64,
    wire: &'a NurbsCurve,
    wire_samples: Vec<(f64, Point3)>,
    band: f64,
}

impl Frame<'_> {
    fn height(&self, p: &Point3) -> f64 {
        (*p - self.origin).dot(&self.normal)
    }
    fn station(&self, p: &Point3) -> Station {
        let h = self.height(p);
        if (h - self.lo).abs() <= self.band {
            Station::Bottom
        } else if (h - self.hi).abs() <= self.band {
            Station::Top
        } else {
            Station::Between
        }
    }
    /// Wire parameter nearest the point's footprint on the sketch plane.
    fn wire_param(&self, p: &Point3) -> f64 {
        let foot = *p - self.normal * self.height(p);
        let foot = Point3::new(foot.x, foot.y, foot.z);
        self.wire_samples
            .iter()
            .map(|(t, q)| (*t, q.distance_to(&foot)))
            .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(t, _)| t)
            .unwrap_or(0.0)
    }
}

fn sorted_by_param<H: Copy>(mut items: Vec<(f64, H)>) -> Vec<H> {
    items.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    items.into_iter().map(|(_, h)| h).collect()
}

fn face_centroid(brep: &BRep, face: FaceHandle) -> Point3 {
    let pts = brep.face_boundary_sample_points(face, 8);
    if pts.is_empty() {
        return Point3::new(0.0, 0.0, 0.0);
    }
    let n = pts.len() as f64;
    let acc = pts.iter().fold(Vec3::new(0.0, 0.0, 0.0), |a, p| {
        a + (*p - Point3::new(0.0, 0.0, 0.0))
    });
    Point3::new(0.0, 0.0, 0.0) + acc * (1.0 / n)
}

fn edge_midpoint(brep: &BRep, edge: EdgeHandle) -> Point3 {
    let e = &brep.edges[edge];
    match e.curve.as_ref() {
        Some(c) => c.geometry.point_at(0.5 * (c.t_start + c.t_end)),
        None => brep.vertices[brep.half_edges[e.half_edge].origin].position,
    }
}

/// Name every entity of the swept body under the extrude contract.
/// `sketch_op` is the sketch feature's operation; `wire` and the
/// sketch frame come from the exact profile; `span` is the signed
/// extent along the normal.
pub(crate) fn name_exact_extrude(
    brep: &mut BRep,
    ctx: &mut TrackingContext,
    sketch_op: OperationId,
    wire: &NurbsCurve,
    wire_origin: Point3,
    normal: Vec3,
    span: (f64, f64),
) {
    brep.name_pool = NamePool::new();
    let (t0, t1) = wire.parameter_range();
    let wire_samples: Vec<(f64, Point3)> = (0..512)
        .map(|i| {
            let t = t0 + (t1 - t0) * (i as f64 / 512.0);
            (t, wire.point_at(t))
        })
        .collect();
    let frame = Frame {
        origin: wire_origin,
        normal,
        lo: span.0.min(span.1),
        hi: span.0.max(span.1),
        wire,
        wire_samples,
        band: 1.0e-6 * (span.1 - span.0).abs().max(1.0),
    };
    let _ = frame.wire;
    let op = ctx.operation;
    let name_of =
        |operation: OperationId, category: RoleCategory, index: usize, ancestors: Vec<NameId>| {
            EntityName {
                operation,
                role: EntityRole {
                    category,
                    index: index as u32,
                },
                ancestors,
            }
        };

    // --- sketch-side entities: bottom vertices and edges, sketch face ---
    let sketch_face = brep.name_pool.intern(
        name_of(sketch_op, RoleCategory::Generated, 0, Vec::new()),
        None,
    );
    let (mut bottom_v, mut top_v) = (Vec::new(), Vec::new());
    for (vh, v) in brep.vertices.iter() {
        let t = frame.wire_param(&v.position);
        match frame.station(&v.position) {
            Station::Bottom => bottom_v.push((t, vh)),
            _ => top_v.push((t, vh)),
        }
    }
    let bottom_v: Vec<VertexHandle> = sorted_by_param(bottom_v);
    let top_v: Vec<VertexHandle> = sorted_by_param(top_v);
    let mut bottom_vertex_names: Vec<NameId> = Vec::new();
    for (i, &vh) in bottom_v.iter().enumerate() {
        let name = name_of(sketch_op, RoleCategory::Junction, i, Vec::new());
        let criteria = capture_vertex(brep, vh).map(GeometricCriteria::Vertex);
        let id = brep.name_pool.intern(name, criteria);
        brep.vertices[vh].name = Some(id);
        bottom_vertex_names.push(id);
    }
    let (mut bottom_e, mut top_e, mut vertical_e) = (Vec::new(), Vec::new(), Vec::new());
    for (eh, _) in brep.edges.iter() {
        let mid = edge_midpoint(brep, eh);
        let t = frame.wire_param(&mid);
        match frame.station(&mid) {
            Station::Bottom => bottom_e.push((t, eh)),
            Station::Top => top_e.push((t, eh)),
            Station::Between => vertical_e.push((t, eh)),
        }
    }
    let bottom_e: Vec<EdgeHandle> = sorted_by_param(bottom_e);
    let top_e: Vec<EdgeHandle> = sorted_by_param(top_e);
    let vertical_e: Vec<EdgeHandle> = sorted_by_param(vertical_e);
    let mut bottom_edge_names: Vec<NameId> = Vec::new();
    for (i, &eh) in bottom_e.iter().enumerate() {
        let name = name_of(sketch_op, RoleCategory::Lateral, i, Vec::new());
        let criteria = capture_edge(brep, eh).map(GeometricCriteria::Edge);
        let id = brep.name_pool.intern(name, criteria);
        brep.edges[eh].name = Some(id);
        bottom_edge_names.push(id);
    }

    // --- extrude-side entities ---
    let (mut bottom_f, mut top_f, mut lateral_f) = (Vec::new(), Vec::new(), Vec::new());
    for (fh, f) in brep.faces.iter() {
        let c = face_centroid(brep, fh);
        let planar_along = f
            .surface
            .as_ref()
            .map(|s| {
                matches!(
                    &s.geometry,
                    popkern_core::geometry::surface::SurfaceGeometry::Plane { normal: n, .. }
                        if n.normalized().map(|n| n.dot(&frame.normal).abs() > 1.0 - 1e-9).unwrap_or(false)
                )
            })
            .unwrap_or(false);
        match (planar_along, frame.station(&c)) {
            (true, Station::Bottom) => bottom_f.push((0.0, fh)),
            (true, Station::Top) => top_f.push((0.0, fh)),
            _ => lateral_f.push((frame.wire_param(&c), fh)),
        }
    }
    let lateral_f: Vec<FaceHandle> = sorted_by_param(lateral_f);
    for (_, fh) in bottom_f {
        let criteria = capture_face(brep, fh).map(GeometricCriteria::Face);
        let name = name_of(op, RoleCategory::Modified, 0, vec![sketch_face]);
        let id = brep.name_pool.intern(name.clone(), criteria);
        brep.faces[fh].name = Some(id);
        ctx.record_modified(
            name_of(sketch_op, RoleCategory::Generated, 0, Vec::new()),
            name,
        );
    }
    for (_, fh) in top_f {
        let criteria = capture_face(brep, fh).map(GeometricCriteria::Face);
        let name = name_of(op, RoleCategory::Cap, 0, Vec::new());
        let id = brep.name_pool.intern(name.clone(), criteria);
        brep.faces[fh].name = Some(id);
        ctx.record_generated(name);
    }
    for (i, &fh) in lateral_f.iter().enumerate() {
        let ancestors = bottom_edge_names.get(i).copied().into_iter().collect();
        let criteria = capture_face(brep, fh).map(GeometricCriteria::Face);
        let name = name_of(op, RoleCategory::Lateral, i, ancestors);
        let id = brep.name_pool.intern(name.clone(), criteria);
        brep.faces[fh].name = Some(id);
        ctx.record_generated(name);
    }
    for (i, &eh) in vertical_e.iter().enumerate() {
        let ancestors = bottom_vertex_names.get(i).copied().into_iter().collect();
        let criteria = capture_edge(brep, eh).map(GeometricCriteria::Edge);
        let name = name_of(op, RoleCategory::Lateral, i, ancestors);
        let id = brep.name_pool.intern(name.clone(), criteria);
        brep.edges[eh].name = Some(id);
        ctx.record_generated(name);
    }
    for (i, &eh) in top_e.iter().enumerate() {
        let ancestors = bottom_edge_names.get(i).copied().into_iter().collect();
        let criteria = capture_edge(brep, eh).map(GeometricCriteria::Edge);
        let name = name_of(op, RoleCategory::Cap, i, ancestors);
        let id = brep.name_pool.intern(name.clone(), criteria);
        brep.edges[eh].name = Some(id);
        ctx.record_generated(name);
    }
    for (i, &vh) in top_v.iter().enumerate() {
        let ancestors = bottom_vertex_names.get(i).copied().into_iter().collect();
        let criteria = capture_vertex(brep, vh).map(GeometricCriteria::Vertex);
        let name = name_of(op, RoleCategory::Generated, i, ancestors);
        let id = brep.name_pool.intern(name.clone(), criteria);
        brep.vertices[vh].name = Some(id);
        ctx.record_generated(name);
    }
}
