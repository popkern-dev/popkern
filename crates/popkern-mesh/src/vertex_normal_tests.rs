use super::*;
use crate::mesh::MeshVertex;
use popkern_core::geometry::Point3;
use popkern_core::topology::FaceHandle;

fn sphere() -> SurfaceGeometry {
    SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.05),
        radius: 0.05,
        axis: popkern_core::geometry::Vec3::Z,
    }
}

fn unit_normal(s: &SurfaceGeometry, u: f64, v: f64) -> Vec3 {
    let n = s.normal_at(u, v);
    let l = n.length();
    Vec3::new(n.x / l, n.y / l, n.z / l)
}

fn chart(s: &SurfaceGeometry) -> std::rc::Rc<ChartDescriptor> {
    ChartDescriptor::of(s, &popkern_core::tolerance::Tolerance::default())
}

/// A sphere's pole is a chart singularity with one normal; a cone's apex
/// is a geometric one with none (the service's rows, in space).
#[test]
fn a_pole_has_a_unique_normal_and_an_apex_none() {
    let s = sphere();
    let pole = s.point_at(3.9, std::f64::consts::PI);
    assert!(!singular_without_unique_normal(&chart(&s), &pole));
    let cone = SurfaceGeometry::Cone {
        apex: Point3::new(0.0, 0.0, 0.0),
        axis: popkern_core::geometry::Vec3::Z,
        half_angle: 0.4,
    };
    let c = chart(&cone);
    let apex = Point3::new(0.0, 0.0, 0.0);
    assert!(
        c.singular_rows
            .iter()
            .any(|r| r.point.distance_to(&apex) < 1e-12),
        "the apex is a singular row: {c:?}"
    );
    assert!(singular_without_unique_normal(&c, &apex));
}

#[test]
fn limit_normal_at_a_pole_matches_the_chart_sense_nearby() {
    // The raw cross product flips just past the pole; the limit normal
    // must agree with the regular chart normal a little inside.
    let s = sphere();
    let pi = std::f64::consts::PI;
    let inside = unit_normal(&s, 3.9, pi - 0.05);
    for v in [pi, pi + 4e-8, pi - 4e-8] {
        let n = limit_normal_at(&chart(&s), &s, (3.9, v)).expect("limit");
        assert!(n.dot(&inside) > 0.99, "v = {v}: {n:?} vs {inside:?}");
    }
}

#[test]
fn repair_fixes_only_the_flipped_pole_vertex() {
    let s = sphere();
    let pi = std::f64::consts::PI;
    // Outward-facing face: stored normals follow the chart (sign +1).
    let uvs = [(3.9, pi + 4e-8), (3.5, 2.75), (3.1, 3.04), (4.7, 2.9)];
    let vertices: Vec<MeshVertex> = uvs
        .iter()
        .map(|&(u, v)| MeshVertex {
            position: s.point_at(u, v),
            normal: unit_normal(&s, u, v),
            uv: (u, v),
        })
        .collect();
    let healthy_before: Vec<Vec3> = vertices[1..].iter().map(|v| v.normal).collect();
    let mut mesh = FaceMesh {
        face: FaceHandle::default(),
        vertices,
        triangles: vec![[0, 1, 2], [0, 2, 3]],
    };
    // The raw normal past the pole is flipped against the chart sense
    // just inside it.
    let inside = unit_normal(&s, 3.9, pi - 0.05);
    assert!(mesh.vertices[0].normal.dot(&inside) < -0.99);
    repair_singular_normals(&mut mesh, &s, &chart(&s));
    assert!(mesh.vertices[0].normal.dot(&inside) > 0.99);
    let healthy_after: Vec<Vec3> = mesh.vertices[1..].iter().map(|v| v.normal).collect();
    assert_eq!(healthy_before, healthy_after, "healthy normals untouched");
}

#[test]
fn repair_keeps_an_inward_face_inward() {
    let s = sphere();
    let pi = std::f64::consts::PI;
    let uvs = [(3.9, pi + 4e-8), (3.5, 2.75), (3.1, 3.04)];
    let vertices: Vec<MeshVertex> = uvs
        .iter()
        .map(|&(u, v)| {
            let n = unit_normal(&s, u, v);
            MeshVertex {
                position: s.point_at(u, v),
                normal: Vec3::new(-n.x, -n.y, -n.z),
                uv: (u, v),
            }
        })
        .collect();
    let mut mesh = FaceMesh {
        face: FaceHandle::default(),
        vertices,
        triangles: vec![[0, 1, 2]],
    };
    repair_singular_normals(&mut mesh, &s, &chart(&s));
    let inside = unit_normal(&s, 3.9, pi - 0.05);
    assert!(mesh.vertices[0].normal.dot(&inside) < -0.99);
}
