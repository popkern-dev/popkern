use super::*;
use crate::mesh_test_helpers::{bounds, box_mesh, face_areas, icosphere, rotated};

fn grid_for(mesh: &SolidMesh, h: f64, shift: [f64; 3]) -> Grid {
    let (lo, hi) = bounds(mesh);
    let lo = [0, 1, 2].map(|a| lo[a] - shift[a]);
    Grid::covering(lo, hi, h).unwrap_or_else(|e| panic!("{e}"))
}

fn clipped(mesh: &SolidMesh, grid: &Grid) -> ClippedSolid {
    clip_solid(grid, mesh).unwrap_or_else(|e| panic!("{e}"))
}

fn piece_areas(c: &ClippedSolid, faces: usize) -> Vec<f64> {
    let mut a = vec![0.0; faces];
    for p in &c.pieces {
        a[p.face.0 as usize] += p.area;
    }
    a
}

fn meshes() -> Vec<(&'static str, SolidMesh, usize)> {
    vec![
        ("box", box_mesh([0.3, -1.1, 2.0], [4.1, 2.2, 3.7]), 6),
        (
            "rotated cube",
            rotated(box_mesh([0.0; 3], [3.0; 3]), 17.0, [1.5; 3]),
            6,
        ),
        ("sphere", icosphere([0.2, 0.1, -0.3], 2.0, 3), 1),
    ]
}

fn assert_close(a: f64, b: f64, rel: f64, what: &str) {
    assert!(
        (a - b).abs() <= rel * b.abs().max(1e-300),
        "{what}: {a} vs {b}"
    );
}

#[test]
fn piece_areas_sum_to_face_areas() {
    for (name, mesh, faces) in meshes() {
        for h in [0.37, 0.8, 1.3] {
            let c = clipped(&mesh, &grid_for(&mesh, h, [0.0; 3]));
            assert_eq!(c.degenerate_triangles, 0);
            let want = face_areas(&mesh, faces);
            for (f, (got, want)) in piece_areas(&c, faces).iter().zip(&want).enumerate() {
                assert_close(*got, *want, 1e-12, &format!("{name} h {h} face {f}"));
            }
        }
    }
}

#[test]
fn pieces_stay_in_their_cells_and_keep_the_normal() {
    for (name, mesh, _) in meshes() {
        let g = grid_for(&mesh, 0.6, [0.0; 3]);
        let c = clipped(&mesh, &g);
        let h = g.cell_size();
        for p in &c.pieces {
            let ijk = g.cell_ijk(p.cell).unwrap_or_else(|| panic!("cell"));
            let lo = g.cell_min(ijk);
            for v in p.polygon.vertices() {
                for a in 0..3 {
                    assert!(
                        v[a] >= lo[a] - 1e-12 * h && v[a] <= lo[a] + h + 1e-12 * h,
                        "{name}: piece outside cell"
                    );
                }
            }
            let [a, b, cc] = mesh.triangles[p.triangle as usize].map(|i| mesh.vertices[i as usize]);
            let n = cross(sub(b, a), sub(cc, a));
            let l = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
            assert_eq!(p.normal, n.map(|x| x / l), "{name}");
        }
        assert!(
            c.pieces
                .windows(2)
                .all(|w| (w[0].cell, w[0].triangle) < (w[1].cell, w[1].triangle))
        );
    }
}

#[test]
fn grid_translation_leaves_face_areas_unchanged() {
    for (name, mesh, faces) in meshes() {
        let reference = piece_areas(&clipped(&mesh, &grid_for(&mesh, 0.7, [0.0; 3])), faces);
        for s in [
            [0.1, 0.0, 0.0],
            [0.0, 0.33, 0.0],
            [0.0, 0.0, 0.5],
            [0.21, 0.42, 0.63],
            [0.69; 3],
        ] {
            let shifted = piece_areas(&clipped(&mesh, &grid_for(&mesh, 0.7, s)), faces);
            for (got, want) in shifted.iter().zip(&reference) {
                assert_close(*got, *want, 1e-12, &format!("{name} shift {s:?}"));
            }
        }
    }
}

#[test]
fn faces_on_grid_planes_are_counted_once() {
    // Box faces at 0 and 2 lie exactly on the planes of a unit grid.
    let mesh = box_mesh([0.0; 3], [2.0; 3]);
    let g = Grid::covering([0.0; 3], [2.0; 3], 1.0).unwrap_or_else(|e| panic!("{e}"));
    let c = clipped(&mesh, &g);
    for (f, a) in piece_areas(&c, 6).iter().enumerate() {
        assert_close(*a, 4.0, 1e-15, &format!("face {f}"));
    }
    // Negative control: the closed-box clip without the tie rule counts
    // every face twice (once in each cell sharing its plane).
    let mut naive = [0.0; 6];
    for (t, corners) in mesh.triangles.iter().enumerate() {
        let [a, b, cc] = corners.map(|i| mesh.vertices[i as usize]);
        for ijk in g.cells() {
            let lo = g.cell_min(ijk);
            let mut p = Polygon::triangle(a, b, cc);
            for (axis, &l) in lo.iter().enumerate() {
                p = p.clip(axis, l, true).clip(axis, l + 1.0, false);
            }
            if p.len() >= 3 {
                naive[mesh.triangle_face[t].0 as usize] += p.area();
            }
        }
    }
    assert!(naive.iter().all(|&a| (a - 8.0).abs() < 1e-12), "{naive:?}");
}

#[test]
fn a_triangle_outside_the_grid_is_an_error() {
    let mesh = box_mesh([0.0; 3], [1.0; 3]);
    let g = Grid::covering([5.0; 3], [6.0; 3], 0.5).unwrap_or_else(|e| panic!("{e}"));
    assert!(matches!(
        clip_solid(&g, &mesh),
        Err(ClipError::OutsideGrid { .. })
    ));
}

/// A block along z with a round bore, its end faces fanned from the
/// bore's ring out to the square's boundary (as a display tessellation
/// fans a holed face): `segments` around, the square's corners among
/// them (a multiple of 8).
fn bored_block(half: f64, r: f64, length: f64, segments: usize) -> SolidMesh {
    let mut vertices: Vec<[f64; 3]> = Vec::new();
    let angle = |i: usize| 2.0 * std::f64::consts::PI * i as f64 / segments as f64;
    use popkern_math::PortableMath;
    // Rings: bore at z 0 and length, square at z 0 and length.
    for (ring, z) in [(0, 0.0), (0, length), (1, 0.0), (1, length)] {
        for i in 0..segments {
            let (c, s) = (angle(i).pcos(), angle(i).psin());
            let p = if ring == 0 {
                [r * c, r * s]
            } else {
                let m = c.abs().max(s.abs());
                [half * c / m, half * s / m]
            };
            vertices.push([p[0], p[1], z]);
        }
    }
    let n = segments as u32;
    let at = |ring: u32, i: usize| ring * n + (i % segments) as u32;
    let mut triangles = Vec::new();
    let mut triangle_face = Vec::new();
    for i in 0..segments {
        let quads = [
            // Bore wall (faces the axis), end at z = 0 (faces −z), end
            // at z = length (+z), outer wall (outward).
            ([at(0, i), at(1, i), at(1, i + 1), at(0, i + 1)], 0),
            ([at(0, i), at(0, i + 1), at(2, i + 1), at(2, i)], 1),
            ([at(1, i), at(3, i), at(3, i + 1), at(1, i + 1)], 1),
            ([at(2, i), at(2, i + 1), at(3, i + 1), at(3, i)], 2),
        ];
        for (q, f) in quads {
            triangles.push([q[0], q[1], q[2]]);
            triangles.push([q[0], q[2], q[3]]);
            triangle_face.extend([FaceId(f), FaceId(f)]);
        }
    }
    SolidMesh {
        id: crate::model::SolidId(0),
        name: None,
        vertices,
        triangles,
        triangle_face,
        exact_volume: None,
    }
}

/// PK-963 reproduction: a bored block whose far end face fans from the
/// bore's ring, one ulp below the grid plane z = 10, out to the square,
/// exactly on it (the bracket's end face as the boolean left it,
/// PK-956). The fan's pieces near its on-plane corners have every
/// clipped vertex on the plane, and the half-open tie rule gives them to
/// the next cell, which clipped the triangle whole against that plane
/// first and holds only a sliver: a third of their area is lost and the
/// rasterised volume reads low. Ignored until the clipping is
/// order-independent under rounding.
#[test]
#[ignore = "known gap: PK-963 facets an ulp off a grid plane lose pieces to the tie rule"]
fn a_facet_one_ulp_off_a_grid_plane_keeps_every_piece() {
    let mut mesh = bored_block(4.0, 1.5, 10.0, 48);
    let below = f64::from_bits(10.0f64.to_bits() - 1);
    for v in mesh.vertices[48..96].iter_mut() {
        v[2] = below;
    }
    let exact = mesh.enclosed_volume();
    assert!(exact > 0.0, "{exact}");
    let grid =
        Grid::covering([-5.0, -5.0, -1.0], [5.0, 5.0, 11.0], 1.0).unwrap_or_else(|e| panic!("{e}"));
    let clipped = clip_solid(&grid, &mesh).unwrap_or_else(|e| panic!("{e}"));
    let m = crate::moments::rasterise_solid(&grid, mesh.id, &clipped, 3)
        .unwrap_or_else(|e| panic!("{e}"));
    let rel = (m.volume - exact) / exact;
    assert!(rel.abs() < 1e-12, "rasterised volume off by {rel:e}");
}
