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
