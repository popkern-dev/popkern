use popkern_core::geometry::Vec3;
use popkern_math::PortableMath;

use super::*;

/// A latitude-longitude unit sphere, outward triangles (none at the
/// poles' degenerate fans).
fn sphere_mesh(n: usize) -> (Vec<[f64; 3]>, Vec<[u32; 3]>) {
    let mut vertices = Vec::new();
    for i in 0..=n {
        let theta = std::f64::consts::PI * i as f64 / n as f64;
        for j in 0..2 * n {
            let phi = std::f64::consts::PI * j as f64 / n as f64;
            vertices.push([
                theta.psin() * phi.pcos(),
                theta.psin() * phi.psin(),
                theta.pcos(),
            ]);
        }
    }
    let at = |i: usize, j: usize| (i * 2 * n + j % (2 * n)) as u32;
    let mut triangles = Vec::new();
    for i in 1..n - 1 {
        for j in 0..2 * n {
            triangles.push([at(i, j), at(i + 1, j), at(i + 1, j + 1)]);
            triangles.push([at(i, j), at(i + 1, j + 1), at(i, j + 1)]);
        }
    }
    (vertices, triangles)
}

fn area(v: &[[f64; 3]], t: [u32; 3]) -> f64 {
    let [a, b, c] = t.map(|i| v[i as usize]);
    let (u, w) = (
        [b[0] - a[0], b[1] - a[1], b[2] - a[2]],
        [c[0] - a[0], c[1] - a[1], c[2] - a[2]],
    );
    let n = [
        u[1] * w[2] - u[2] * w[1],
        u[2] * w[0] - u[0] * w[2],
        u[0] * w[1] - u[1] * w[0],
    ];
    0.5 * (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt()
}

#[test]
fn flipped_triangles_are_folded_and_far_projections_are_not_counted() {
    let (mut vertices, mut triangles) = sphere_mesh(24);
    let geometry: BTreeMap<FaceId, SurfaceGeometry> = [(
        FaceId(0),
        SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
            axis: Vec3::new(0.0, 0.0, 1.0),
        },
    )]
    .into();
    let faces = vec![FaceId(0); triangles.len()];
    // Chords stay within about 0.01 of a unit sphere at 24 rings.
    let reach = 0.05;
    let sound = folded_areas(&geometry, &vertices, &triangles, &faces, reach);
    assert_eq!(sound.get(&FaceId(0)).copied(), Some(0.0));
    // Fold five triangles back on themselves.
    let mut flipped = 0.0;
    for t in triangles.iter_mut().skip(100).step_by(37).take(5) {
        t.swap(1, 2);
        flipped += area(&vertices, *t);
    }
    let folded = folded_areas(&geometry, &vertices, &triangles, &faces, reach);
    let got = folded.get(&FaceId(0)).copied().unwrap_or(0.0);
    assert!((got - flipped).abs() <= 1e-12, "{got} vs {flipped}");
    // A flipped triangle pushed off the surface projects too far away
    // and is not counted (negative control for the distance filter).
    let far = triangles[100];
    for i in far {
        let p = &mut vertices[i as usize];
        *p = p.map(|x| 1.5 * x);
    }
    let pushed = folded_areas(&geometry, &vertices, &triangles, &faces, reach);
    let less = pushed.get(&FaceId(0)).copied().unwrap_or(0.0);
    assert!(less < got, "{less} vs {got}");
}
