//! Shared helpers for the PK-387 circle-extrude gates.
#![allow(dead_code)]

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::naming::{EntityName, EntityRole, OperationId, RoleCategory};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_history::*;
use popkern_sketch::{AnalyticProfile, ProfileShape, Sketch};

pub fn name(op: u64, category: RoleCategory, index: u32) -> EntityName {
    EntityName {
        operation: OperationId(op),
        role: EntityRole { category, index },
        ancestors: vec![],
    }
}

pub fn add_param(tree: &mut FeatureTree, label: &str, value: f64) -> ParameterId {
    tree.add_global_parameter(label, value, Unit::Meters)
        .expect("param")
}

pub fn sorted_vertices(brep: &BRep) -> Vec<[f64; 3]> {
    let mut v: Vec<[f64; 3]> = brep
        .vertices
        .values()
        .map(|x| [x.position.x, x.position.y, x.position.z])
        .collect();
    v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    v
}

pub fn volume(brep: &BRep) -> f64 {
    let mesh = popkern_mesh::tessellate(brep, &popkern_mesh::TessellationOptions::default())
        .expect("mesh");
    let mut v6 = 0.0;
    for t in &mesh.triangles {
        let (a, b, c) = (
            mesh.vertices[t[0]].position,
            mesh.vertices[t[1]].position,
            mesh.vertices[t[2]].position,
        );
        v6 += a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
            + a.z * (b.x * c.y - b.y * c.x);
    }
    v6 / 6.0
}

/// Same census, same vertices, same mesh volume — and the same
/// triangle count, which only agrees when the two pipelines chart the
/// same faces the same way (one seamed wall, one-edge caps).
pub fn assert_same_body(label: &str, a: &BRep, b: &BRep) {
    assert_eq!(a.faces.len(), b.faces.len(), "{label}: face count");
    assert_eq!(a.edges.len(), b.edges.len(), "{label}: edge count");
    assert_eq!(a.vertices.len(), b.vertices.len(), "{label}: vertex count");
    for (p, q) in sorted_vertices(a).iter().zip(sorted_vertices(b).iter()) {
        for k in 0..3 {
            assert!((p[k] - q[k]).abs() < 1e-9, "{label}: vertex {p:?} vs {q:?}");
        }
    }
    let opts = popkern_mesh::TessellationOptions::default();
    let ma = popkern_mesh::tessellate(a, &opts).expect("mesh a");
    let mb = popkern_mesh::tessellate(b, &opts).expect("mesh b");
    assert_eq!(
        ma.triangles.len(),
        mb.triangles.len(),
        "{label}: triangle count"
    );
    let (va, vb) = (volume(a), volume(b));
    assert!(
        (va - vb).abs() < 1e-9 * va.abs().max(1.0),
        "{label}: volume {va} vs {vb}"
    );
}

pub fn surface_census(brep: &BRep) -> (usize, usize, usize) {
    let (mut planes, mut cylinders, mut other) = (0, 0, 0);
    for f in brep.faces.values() {
        match f.surface.as_ref().map(|s| &s.geometry) {
            Some(SurfaceGeometry::Plane { .. }) => planes += 1,
            Some(SurfaceGeometry::Cylinder { .. }) => cylinders += 1,
            _ => other += 1,
        }
    }
    (planes, cylinders, other)
}

/// Planar faces whose loop is a single closed edge (canonical caps).
pub fn one_edge_planar_loops(brep: &BRep) -> usize {
    brep.faces
        .values()
        .filter(|f| {
            matches!(
                f.surface.as_ref().map(|s| &s.geometry),
                Some(SurfaceGeometry::Plane { .. })
            )
        })
        .filter(|f| {
            let start = brep.loops[f.outer_loop].half_edge;
            brep.half_edges[start].next == start
        })
        .count()
}

pub fn circle_sketch_json(cx: f64, cy: f64, r: f64) -> String {
    let mut sk = Sketch::new();
    sk.add_arc(cx, cy, r, 0.0, std::f64::consts::TAU);
    sk.solve().expect("circle solves");
    serde_json::to_string(&sk).expect("serialize")
}

pub fn analytic_circle(cx: f64, cy: f64, r: f64) -> SketchPayload {
    SketchPayload::Analytic(AnalyticProfile {
        shape: ProfileShape::Circle { radius: r.into() },
        center: [cx.into(), cy.into()],
    })
}

pub fn served(tree: &mut FeatureTree) -> BRep {
    let result = tree.replay(Tolerance::default()).expect("replay");
    assert!(
        result.report.failures().next().is_none(),
        "{:?}",
        result.report
    );
    let body = tree.current_brep().expect("body");
    assert!(body.validate_all().is_empty(), "{:?}", body.validate_all());
    body.clone()
}

/// The served name of the first face under `op` with `category` on the
/// current body (names carry ancestors, so references read the pool).
pub fn find_face_name(tree: &FeatureTree, op: u64, category: RoleCategory) -> EntityName {
    let brep = tree.current_brep().expect("body");
    brep.faces
        .values()
        .filter_map(|f| f.name.and_then(|id| brep.name_pool.get(id).cloned()))
        .find(|n| n.operation == OperationId(op) && n.role.category == category)
        .expect("named face")
}
