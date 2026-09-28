use super::*;
use popkern_core::topology::FaceHandle;

fn unit(v: Vec3) -> Vec3 {
    let l = v.length();
    Vec3::new(v.x / l, v.y / l, v.z / l)
}

fn cylinder() -> SurfaceGeometry {
    SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 0.05,
    }
}

/// A quarter-cylinder chart cell meshed as two triangles plus a
/// boundary-only fan (the evacuated-strip shape): every triangle
/// spans the full 90°.
fn coarse_quarter(surface: &SurfaceGeometry, flip: bool) -> FaceMesh {
    let corners = [
        (0.0, 0.0),
        (std::f64::consts::FRAC_PI_2, 0.0),
        (std::f64::consts::FRAC_PI_2, 1.0),
        (0.0, 1.0),
    ];
    let vertices = corners
        .iter()
        .map(|&(u, v)| {
            let n = unit(surface.normal_at(u, v));
            MeshVertex {
                position: surface.point_at(u, v),
                normal: if flip { -n } else { n },
                uv: (u, v),
            }
        })
        .collect();
    FaceMesh {
        face: FaceHandle::default(),
        vertices,
        triangles: vec![[0, 1, 2], [0, 2, 3]],
    }
}

fn boundary_edges(mesh: &FaceMesh) -> Vec<(usize, usize)> {
    let mut count: HashMap<(usize, usize), usize> = HashMap::new();
    for t in &mesh.triangles {
        for k in 0..3 {
            let (a, b) = (t[k], t[(k + 1) % 3]);
            *count.entry((a.min(b), a.max(b))).or_default() += 1;
        }
    }
    let mut v: Vec<_> = count
        .into_iter()
        .filter(|&(_, c)| c == 1)
        .map(|(e, _)| e)
        .collect();
    v.sort();
    v
}

#[test]
fn boundary_limited_triangles_are_left_unsplit() {
    // Every triangle touches a 90° boundary edge: no interior split can
    // bring it within the bound, so the pass must leave the mesh as is
    // rather than fan slivers onto the boundary.
    let surface = cylinder();
    let mut mesh = coarse_quarter(&surface, false);
    let before = boundary_edges(&mesh);
    refine_face_mesh(&mut mesh, &surface, 12f64.to_radians());
    assert_eq!(boundary_edges(&mesh), before);
    assert_eq!(mesh.triangles.len(), 2);
    assert_eq!(mesh.vertices.len(), 4);
}

#[test]
fn refinement_reaches_the_bound_when_the_boundary_is_fine() {
    // Boundary resampled at the bound (as the edge sampler does), a
    // boundary-only fan across the chart: interior splits must bring
    // every triangle within the bound.
    let surface = cylinder();
    let n = 10usize;
    let du = std::f64::consts::FRAC_PI_2 / n as f64;
    let mut uvs: Vec<(f64, f64)> = (0..=n).map(|i| (i as f64 * du, 0.0)).collect();
    uvs.extend((0..=n).rev().map(|i| (i as f64 * du, 1.0)));
    let vertices: Vec<MeshVertex> = uvs
        .iter()
        .map(|&(u, v)| MeshVertex {
            position: surface.point_at(u, v),
            normal: unit(surface.normal_at(u, v)),
            uv: (u, v),
        })
        .collect();
    // Fan from vertex 0: long rail-to-rail chords, the evacuated shape.
    let triangles: Vec<[usize; 3]> = (1..uvs.len() - 1).map(|k| [0, k, k + 1]).collect();
    let mut mesh = FaceMesh {
        face: FaceHandle::default(),
        vertices,
        triangles,
    };
    let bound = 12f64.to_radians();
    refine_face_mesh(&mut mesh, &surface, bound);
    let worst = mesh
        .triangles
        .iter()
        .map(|t| {
            (0..3)
                .map(|i| {
                    angle_between(
                        mesh.vertices[t[i]].normal,
                        mesh.vertices[t[(i + 1) % 3]].normal,
                    )
                })
                .fold(0.0, f64::max)
        })
        .fold(0.0, f64::max);
    assert!(
        worst <= bound + 1e-9,
        "worst spread {}°",
        worst.to_degrees()
    );
}

#[test]
fn inserted_normals_follow_a_flipped_face() {
    let surface = cylinder();
    let mut mesh = coarse_quarter(&surface, true);
    refine_face_mesh(&mut mesh, &surface, 12f64.to_radians());
    for v in &mesh.vertices {
        assert!(v.normal.dot(&unit(surface.normal_at(v.uv.0, v.uv.1))) < -0.99);
    }
}

#[test]
fn winding_is_preserved_by_splits() {
    let surface = cylinder();
    let mut mesh = coarse_quarter(&surface, false);
    let uv_area = |m: &FaceMesh, t: &[usize; 3]| {
        let (a, b, c) = (
            m.vertices[t[0]].uv,
            m.vertices[t[1]].uv,
            m.vertices[t[2]].uv,
        );
        (b.0 - a.0) * (c.1 - a.1) - (c.0 - a.0) * (b.1 - a.1)
    };
    let sign = uv_area(&mesh, &mesh.triangles[0]).signum();
    refine_face_mesh(&mut mesh, &surface, 12f64.to_radians());
    for t in &mesh.triangles {
        assert_eq!(uv_area(&mesh, t).signum(), sign);
    }
}

/// PK-587: an interior chord along a thin torus's TOP circle has equal
/// vertex normals, so spread alone never splits it, but its midpoint
/// lies over the tube 5.7° around from the top. The chord-foot vote
/// bisects it until every interior chord is within the foot share.
#[test]
fn a_torus_top_chord_is_split_by_its_foot_error() {
    let torus = SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: 1.0,
        minor_radius: 0.05,
    };
    let top = 0.5 * std::f64::consts::PI;
    let corners = [(0.0, top), (0.2, top), (0.1, top + 0.02), (0.1, top - 0.02)];
    let vertex = |&(u, v): &(f64, f64)| MeshVertex {
        position: torus.point_at(u, v),
        normal: unit(torus.normal_at(u, v)),
        uv: (u, v),
    };
    let mut mesh = FaceMesh {
        face: FaceHandle::default(),
        vertices: corners.iter().map(vertex).collect(),
        triangles: vec![[0, 1, 2], [1, 0, 3]],
    };
    let bound = 3.0_f64.to_radians();
    let foot = |m: &FaceMesh, a: usize, b: usize| {
        let (va, vb) = (&m.vertices[a], &m.vertices[b]);
        let mid = (0.5 * (va.uv.0 + vb.uv.0), 0.5 * (va.uv.1 + vb.uv.1));
        crate::facet_angle::chord_foot_error(
            &torus,
            mid,
            (va.position, vb.position),
            (va.normal, vb.normal),
        )
        .unwrap_or(0.0)
    };
    assert!(foot(&mesh, 0, 1) > 5.0_f64.to_radians());
    let boundary_before = boundary_edges(&mesh);
    refine_face_mesh(&mut mesh, &torus, bound);
    assert!(mesh.triangles.len() > 2, "the top chord must split");
    assert_eq!(boundary_edges(&mesh).len(), boundary_before.len());
    let boundary = boundary_edges(&mesh);
    for t in &mesh.triangles {
        for k in 0..3 {
            let (a, b) = (t[k], t[(k + 1) % 3]);
            if boundary.contains(&(a.min(b), a.max(b))) {
                continue;
            }
            let e = foot(&mesh, a, b);
            assert!(
                e <= crate::facet_angle::FOOT_SHARE * bound + 1e-12,
                "interior chord {a}-{b}: foot error {:.3}°",
                e.to_degrees()
            );
        }
    }
}
