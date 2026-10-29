//! The orientation rules on meshes whose winding is known.

use super::*;

/// A unit tetrahedron, wound outward or inward.
fn tetrahedron(outward: bool) -> TriMesh {
    let mut mesh = TriMesh::default();
    for p in [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(0.0, 0.0, 1.0),
    ] {
        mesh.vertices.push(popkern_mesh::MeshVertex {
            position: p,
            normal: Vec3::new(0.0, 0.0, 0.0),
            uv: (0.0, 0.0),
        });
    }
    let out = [[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]];
    mesh.triangles = out
        .iter()
        .map(|&[a, b, c]| if outward { [a, b, c] } else { [a, c, b] })
        .collect();
    mesh
}

const IDENTITY: [[f64; 4]; 4] = [
    [1.0, 0.0, 0.0, 0.0],
    [0.0, 1.0, 0.0, 0.0],
    [0.0, 0.0, 1.0, 0.0],
    [0.0, 0.0, 0.0, 1.0],
];
const MIRROR_X: [[f64; 4]; 4] = [
    [-1.0, 0.0, 0.0, 2.0],
    [0.0, 1.0, 0.0, 0.0],
    [0.0, 0.0, 1.0, 0.0],
    [0.0, 0.0, 0.0, 1.0],
];

#[test]
fn every_winding_and_placement_writes_outward_facets() {
    for outward in [true, false] {
        for (name, placement) in [("identity", IDENTITY), ("mirrored", MIRROR_X)] {
            let mesh = tetrahedron(outward);
            let f = facets(&[Placed {
                mesh: &mesh,
                placement,
            }]);
            let (open, volume) = census(&f);
            assert_eq!(open, 0, "{name}");
            assert!(
                (volume - 1.0 / 6.0).abs() < 1e-12,
                "outward {outward}, {name}: volume {volume}"
            );
        }
    }
}

#[test]
fn the_encoding_carries_unit_normals_along_the_winding() {
    let mesh = tetrahedron(false);
    let bytes = encode(&facets(&[Placed {
        mesh: &mesh,
        placement: IDENTITY,
    }]));
    assert_eq!(bytes.len(), 84 + 50 * 4);
    // The first facet is the base (z = 0): its outward normal is −z.
    let f = |o: usize| f32::from_le_bytes(bytes[o..o + 4].try_into().expect("f32"));
    assert_eq!([f(84), f(88), f(92)], [0.0, 0.0, -1.0]);
}
