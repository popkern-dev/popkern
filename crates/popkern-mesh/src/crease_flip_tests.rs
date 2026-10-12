use super::*;
use crate::mesh::MeshVertex;
use popkern_core::topology::FaceHandle;

fn torus() -> SurfaceGeometry {
    SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: 1.0,
        minor_radius: 0.2,
    }
}

/// An `nu × nv` lattice over `[0, 0.6] × [0, 1.2]` of the torus chart
/// whose cell diagonals alternate by row: the staggered pattern that
/// creases an anisotropic chart (long cells along the tube's length).
fn staggered_lattice(surface: &SurfaceGeometry, nu: usize, nv: usize) -> FaceMesh {
    let mut vertices = Vec::new();
    for j in 0..=nv {
        for i in 0..=nu {
            let (u, v) = (0.6 * i as f64 / nu as f64, 1.2 * j as f64 / nv as f64);
            let n = surface.normal_at(u, v);
            vertices.push(MeshVertex {
                position: surface.point_at(u, v),
                normal: n * (1.0 / n.length()),
                uv: (u, v),
            });
        }
    }
    let at = |i: usize, j: usize| j * (nu + 1) + i;
    let mut triangles = Vec::new();
    for j in 0..nv {
        for i in 0..nu {
            let (a, b, c, d) = (at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1));
            if (i + j) % 2 == 0 {
                triangles.extend([[a, b, c], [a, c, d]]);
            } else {
                triangles.extend([[a, b, d], [b, c, d]]);
            }
        }
    }
    FaceMesh {
        face: FaceHandle::default(),
        vertices,
        triangles,
    }
}

fn edge_census(mesh: &FaceMesh) -> HashMap<EdgeKey, Vec<usize>> {
    let mut edges: HashMap<EdgeKey, Vec<usize>> = HashMap::new();
    for (ti, t) in mesh.triangles.iter().enumerate() {
        for k in 0..3 {
            edges.entry(key(t[k], t[(k + 1) % 3])).or_default().push(ti);
        }
    }
    edges
}

fn max_crease(mesh: &FaceMesh) -> f64 {
    let edges = edge_census(mesh);
    edges
        .keys()
        .map(|&e| crease(mesh, &edges, e))
        .fold(0.0, f64::max)
}

fn boundary(mesh: &FaceMesh) -> Vec<EdgeKey> {
    let mut b: Vec<EdgeKey> = edge_census(mesh)
        .into_iter()
        .filter(|(_, ts)| ts.len() == 1)
        .map(|(e, _)| e)
        .collect();
    b.sort_unstable();
    b
}

fn vertex(surface: &SurfaceGeometry, u: f64, v: f64) -> MeshVertex {
    let n = surface.normal_at(u, v);
    MeshVertex {
        position: surface.point_at(u, v),
        normal: n * (1.0 / n.length()),
        uv: (u, v),
    }
}

#[test]
fn a_chord_sliver_across_the_curvature_flips_to_the_ruling() {
    // Quad P0 (0,0), P1 (0.4,-0.02), P2 (0.8,0), P3 (0.4,1) on a unit
    // cylinder, split along P0-P2: the triangle P2 P0 P1 is a chord
    // sliver across 0.8 rad of the circle, tilted far off the
    // surface. The diagonal P1-P3 follows the ruling.
    let surface = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    let vertices = [(0.0, 0.0), (0.4, -0.02), (0.8, 0.0), (0.4, 1.0)]
        .iter()
        .map(|&(u, v)| vertex(&surface, u, v))
        .collect();
    let mut mesh = FaceMesh {
        face: FaceHandle::default(),
        vertices,
        triangles: vec![[0, 2, 3], [2, 0, 1]],
    };
    let before = max_crease(&mesh);
    flip_creases(&mut mesh, &surface, 0.005, 0.01, 0.21);
    let after = max_crease(&mesh);
    assert!(
        edge_census(&mesh).contains_key(&key(1, 3)),
        "diagonal flipped to the ruling: {:?}",
        mesh.triangles
    );
    assert!(
        after < 0.5 * before,
        "crease {:.1}° -> {:.1}°",
        before.to_degrees(),
        after.to_degrees()
    );
}

#[test]
fn flips_keep_count_boundary_and_chart_orientation() {
    let surface = torus();
    let mut mesh = staggered_lattice(&surface, 6, 3);
    let (tris, edges) = (mesh.triangles.len(), boundary(&mesh));
    let sense = chart_area(&mesh, mesh.triangles[0]).signum();
    let before = max_crease(&mesh);
    flip_creases(&mut mesh, &surface, 0.005, 0.01, 0.21);
    assert!(
        max_crease(&mesh) <= before + 1e-12,
        "worst crease never grows"
    );
    assert_eq!(mesh.triangles.len(), tris, "flips keep the triangle count");
    assert_eq!(boundary(&mesh), edges, "boundary edges never move");
    for &t in &mesh.triangles {
        let area = chart_area(&mesh, t);
        assert!(
            area.signum() == sense && area.abs() > 1e-12,
            "chart fold at {t:?}"
        );
    }
}

#[test]
fn a_plane_is_left_alone() {
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let mut mesh = staggered_lattice(&plane, 3, 3);
    let original = mesh.triangles.clone();
    flip_creases(&mut mesh, &plane, 0.005, 0.01, 0.21);
    assert_eq!(mesh.triangles, original);
}
