use super::*;
use crate::mesh_test_helpers::{box_mesh, icosphere, tube};

/// (enclosed volume, |∮ n dA|, triangles' total area) of a cut mesh.
fn closure(c: &CutMesh) -> (f64, f64, f64) {
    let (mut vol, mut n, mut area) = (0.0, [0.0; 3], 0.0);
    for t in &c.triangles {
        let [a, b, d] = t.map(|i| c.vertices[i as usize]);
        let (u, w) = (
            [b[0] - a[0], b[1] - a[1], b[2] - a[2]],
            [d[0] - a[0], d[1] - a[1], d[2] - a[2]],
        );
        let x = [
            u[1] * w[2] - u[2] * w[1],
            u[2] * w[0] - u[0] * w[2],
            u[0] * w[1] - u[1] * w[0],
        ];
        vol += (a[0] * (b[1] * d[2] - b[2] * d[1]) - a[1] * (b[0] * d[2] - b[2] * d[0])
            + a[2] * (b[0] * d[1] - b[1] * d[0]))
            / 6.0;
        for k in 0..3 {
            n[k] += 0.5 * x[k];
        }
        area += 0.5 * (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt();
    }
    (vol, (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt(), area)
}

#[test]
fn a_box_cut_to_a_smaller_box_is_that_box() {
    let mesh = box_mesh([0.0; 3], [1.0; 3]);
    let c =
        cut_to_box(&mesh, [0.25, -1.0, 0.3], [0.75, 2.0, 0.6]).unwrap_or_else(|e| panic!("{e}"));
    let (vol, n, _) = closure(&c);
    assert!((vol - 0.5 * 0.3).abs() < 1e-12, "{vol}");
    assert!(n < 1e-12, "{n}");
    // Caps on the four faces that cut; none on the y faces beyond it.
    let caps: std::collections::BTreeSet<(usize, bool)> = c
        .origin
        .iter()
        .filter_map(|o| match o {
            Origin::Cap { axis, upper } => Some((*axis, *upper)),
            Origin::Face(_) => None,
        })
        .collect();
    let want: std::collections::BTreeSet<(usize, bool)> =
        [(0, false), (0, true), (2, false), (2, true)].into();
    assert_eq!(caps, want);
}

#[test]
fn a_tube_cut_across_its_axis_gets_annular_caps() {
    let (ri, ro, length) = (1.0, 1.5, 4.0);
    let mesh = tube(ri, ro, length, 48);
    let full = {
        let whole = CutMesh {
            vertices: mesh.vertices.clone(),
            triangles: mesh.triangles.clone(),
            origin: mesh
                .triangle_face
                .iter()
                .map(|f| Origin::Face(*f))
                .collect(),
            lo: [0.0; 3],
            hi: [0.0; 3],
        };
        closure(&whole).0
    };
    let (lo, hi) = bounds_of(&mesh);
    // Keep the middle half along the axis (z), all of it radially.
    let c = cut_to_box(
        &mesh,
        [lo[0] - 1.0, lo[1] - 1.0, lo[2] + 0.25 * (hi[2] - lo[2])],
        [hi[0] + 1.0, hi[1] + 1.0, lo[2] + 0.75 * (hi[2] - lo[2])],
    )
    .unwrap_or_else(|e| panic!("{e}"));
    let (vol, n, _) = closure(&c);
    assert!(
        (vol - 0.5 * full).abs() < 1e-9 * full,
        "{vol} vs {}",
        0.5 * full
    );
    assert!(n < 1e-9, "{n}");
}

fn bounds_of(m: &SolidMesh) -> ([f64; 3], [f64; 3]) {
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for v in &m.vertices {
        for a in 0..3 {
            lo[a] = lo[a].min(v[a]);
            hi[a] = hi[a].max(v[a]);
        }
    }
    (lo, hi)
}

#[test]
fn a_sphere_cut_through_its_centre_closes() {
    let mesh = icosphere([0.0; 3], 1.0, 3);
    let c =
        cut_to_box(&mesh, [-2.0, -2.0, -2.0], [2.0, 0.0, 0.3]).unwrap_or_else(|e| panic!("{e}"));
    let (vol, n, _) = closure(&c);
    assert!(n < 1e-9, "{n}");
    // Between half the ball and the slab, a spherical-segment volume:
    // y < 0 and z < 0.3 of the unit ball (polyhedral, so a little less).
    assert!(
        vol > 0.5 * 0.55 * 4.18879 && vol < 0.5 * 0.75 * 4.18879,
        "{vol}"
    );
}

#[test]
fn degenerate_fan_triangles_weld_their_crossings() {
    // The quarter plate's fan from its corner leaves zero-area triangles
    // along an edge; their overlapping edges cross a cut at one point.
    let (mesh, _, _) =
        crate::mesh_test_helpers::quarter_plate_with_hole(18.0, 6.0, 2.0, 2.0, 48, 18);
    let c =
        cut_to_box(&mesh, [-1.0, -1.0, -1.0], [6.0, 4.0, 3.0]).unwrap_or_else(|e| panic!("{e}"));
    let (vol, n, _) = closure(&c);
    // (6 · 4 − π · 2² / 4) · 2, the hole a 48-gon.
    let exact = (24.0 - std::f64::consts::PI) * 2.0;
    assert!((vol - exact).abs() < 3e-3, "{vol} vs {exact}");
    assert!(n < 1e-12, "{n}");
}
