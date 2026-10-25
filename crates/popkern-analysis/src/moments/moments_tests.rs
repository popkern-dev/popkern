use super::*;
use crate::clip::clip_solid;
use crate::linalg::spd_test_helpers::bits_hash;
use crate::mesh_test_helpers::{bounds, box_mesh, icosphere, rotated};
use crate::model::SolidMesh;

fn grid_for(mesh: &SolidMesh, h: f64, shift: [f64; 3]) -> Grid {
    let (lo, hi) = bounds(mesh);
    Grid::covering([0, 1, 2].map(|a| lo[a] - shift[a]), hi, h).unwrap_or_else(|e| panic!("{e}"))
}

fn rasterise(mesh: &SolidMesh, grid: &Grid) -> SolidMoments {
    let clipped = clip_solid(grid, mesh).unwrap_or_else(|e| panic!("{e}"));
    rasterise_solid(grid, SolidId(0), &clipped, 2).unwrap_or_else(|e| panic!("{e}"))
}

/// ∫ over [lo, hi] (physical) of ξ^a dx, ξ local to a cell starting at
/// `cmin` with size h.
fn axis_integral(a: usize, lo: f64, hi: f64, cmin: f64, h: f64) -> f64 {
    let xi = |x: f64| 2.0 * (x - cmin) / h - 1.0;
    let p = (a + 1) as i32;
    0.5 * h * (xi(hi).powi(p) - xi(lo).powi(p)) / (a + 1) as f64
}

#[test]
fn box_cells_match_closed_forms() {
    let (min, max) = ([0.3, -1.1, 2.0], [4.1, 2.2, 3.7]);
    let mesh = box_mesh(min, max);
    for h in [0.7, 0.45] {
        let g = grid_for(&mesh, h, [0.0; 3]);
        let r = rasterise(&mesh, &g);
        let mut expected_cells = 0;
        for ijk in g.cells() {
            let cmin = g.cell_min(ijk);
            let lo = [0, 1, 2].map(|a| cmin[a].max(min[a]));
            let hi = [0, 1, 2].map(|a| (cmin[a] + h).min(max[a]));
            if (0..3).any(|a| hi[a] <= lo[a]) {
                continue;
            }
            expected_cells += 1;
            let cell = g.cell_index(ijk).unwrap_or_default();
            let got = r
                .cells
                .iter()
                .find(|c| c.cell == cell)
                .unwrap_or_else(|| panic!("cell {ijk:?} missing"));
            for c in 0..3 {
                for b in 0..3 {
                    for a in 0..3 {
                        let want = axis_integral(a, lo[0], hi[0], cmin[0], h)
                            * axis_integral(b, lo[1], hi[1], cmin[1], h)
                            * axis_integral(c, lo[2], hi[2], cmin[2], h);
                        let m = got.m[moment_index(a, b, c)];
                        assert!(
                            (m - want).abs() <= 1e-12 * h * h * h,
                            "h {h} cell {ijk:?} m{a}{b}{c}: {m} vs {want}"
                        );
                    }
                }
            }
        }
        assert_eq!(r.cells.len(), expected_cells, "h {h}");
        assert!(r.full > 0 && r.cut > 0);
    }
}

/// Volume, first and second moments (10 values) of a closed mesh by
/// signed tetrahedra from the origin: an independent reference.
fn tetra_moments(mesh: &SolidMesh) -> [f64; 10] {
    let mut out = [0.0; 10];
    for t in &mesh.triangles {
        let p = t.map(|i| mesh.vertices[i as usize]);
        let det = p[0][0] * (p[1][1] * p[2][2] - p[1][2] * p[2][1])
            - p[0][1] * (p[1][0] * p[2][2] - p[1][2] * p[2][0])
            + p[0][2] * (p[1][0] * p[2][1] - p[1][1] * p[2][0]);
        let v = det / 6.0;
        let s = [0, 1, 2].map(|i| p[0][i] + p[1][i] + p[2][i]);
        out[0] += v;
        for i in 0..3 {
            out[1 + i] += v * s[i] / 4.0;
        }
        let mut k = 4;
        for i in 0..3 {
            for j in i..3 {
                let sq: f64 = (0..3).map(|n| p[n][i] * p[n][j]).sum();
                out[k] += v / 20.0 * (sq + s[i] * s[j]);
                k += 1;
            }
        }
    }
    out
}

/// The same 10 global moments assembled from cell moments.
fn global_moments(r: &SolidMoments, g: &Grid) -> [f64; 10] {
    let h = g.cell_size();
    let e = [
        moment_index(1, 0, 0),
        moment_index(0, 1, 0),
        moment_index(0, 0, 1),
    ];
    let mut out = [0.0; 10];
    for cell in &r.cells {
        let ijk = g.cell_ijk(cell.cell).unwrap_or_default();
        let c = g.cell_min(ijk).map(|x| x + 0.5 * h);
        let m = &cell.m;
        out[0] += m[0];
        for i in 0..3 {
            out[1 + i] += c[i] * m[0] + 0.5 * h * m[e[i]];
        }
        let mut k = 4;
        for i in 0..3 {
            for j in i..3 {
                let mut unit = [0usize; 3];
                unit[i] += 1;
                unit[j] += 1;
                let second = m[moment_index(unit[0], unit[1], unit[2])];
                out[k] += c[i] * c[j] * m[0]
                    + 0.5 * h * (c[i] * m[e[j]] + c[j] * m[e[i]])
                    + 0.25 * h * h * second;
                k += 1;
            }
        }
    }
    out
}

fn assert_moments_close(got: &[f64; 10], want: &[f64; 10], scale: [f64; 3], what: &str) {
    // Compare at the magnitude of each moment family: V, V·L, V·L².
    let [v, l1, l2] = scale;
    for (k, (g, w)) in got.iter().zip(want).enumerate() {
        let s = match k {
            0 => v,
            1..=3 => l1,
            _ => l2,
        };
        assert!((g - w).abs() <= 1e-12 * s, "{what} moment {k}: {g} vs {w}");
    }
}

fn shapes() -> Vec<(&'static str, SolidMesh)> {
    vec![
        ("box", box_mesh([0.3, -1.1, 2.0], [4.1, 2.2, 3.7])),
        (
            "rotated cube",
            rotated(box_mesh([0.0; 3], [3.0; 3]), 17.0, [1.5; 3]),
        ),
        ("sphere", icosphere([0.2, 0.1, -0.3], 2.0, 3)),
    ]
}

fn scale_of(mesh: &SolidMesh, reference: &[f64; 10]) -> [f64; 3] {
    let (lo, hi) = bounds(mesh);
    let l = (0..3)
        .map(|a| lo[a].abs().max(hi[a].abs()))
        .fold(0.0, f64::max);
    [reference[0], reference[0] * l, reference[0] * l * l]
}

#[test]
fn global_moments_match_the_polyhedron() {
    for (name, mesh) in shapes() {
        let want = tetra_moments(&mesh);
        for h in [0.37, 0.8] {
            let g = grid_for(&mesh, h, [0.0; 3]);
            let got = global_moments(&rasterise(&mesh, &g), &g);
            assert_moments_close(
                &got,
                &want,
                scale_of(&mesh, &want),
                &format!("{name} h {h}"),
            );
        }
    }
}

#[test]
fn grid_translation_leaves_global_moments_unchanged() {
    for (name, mesh) in shapes() {
        let want = tetra_moments(&mesh);
        for s in [
            [0.1, 0.0, 0.0],
            [0.0, 0.33, 0.0],
            [0.0, 0.0, 0.5],
            [0.21, 0.42, 0.63],
            [0.69; 3],
        ] {
            let g = grid_for(&mesh, 0.7, s);
            let got = global_moments(&rasterise(&mesh, &g), &g);
            assert_moments_close(
                &got,
                &want,
                scale_of(&mesh, &want),
                &format!("{name} shift {s:?}"),
            );
        }
    }
}

#[test]
fn x5_volume_parity_holds_and_mismatch_is_typed() {
    for (name, mesh) in shapes() {
        let g = grid_for(&mesh, 0.6, [0.0; 3]);
        let r = rasterise(&mesh, &g);
        assert_eq!(
            volume_parity(&r, mesh.enclosed_volume(), 1e-12),
            Ok(()),
            "{name}"
        );
        assert_eq!(
            volume_parity(&r, mesh.enclosed_volume() * 1.001, 1e-12),
            Err(AnalysisError::DegenerateGeometry {
                solid: Some(SolidId(0)),
                detail: Degeneracy::VolumeMismatch
            }),
            "{name}"
        );
    }
}

#[test]
fn io1_an_inside_out_shell_declines() {
    let mut mesh = box_mesh([0.0; 3], [2.0, 1.0, 1.5]);
    for t in &mut mesh.triangles {
        t.swap(1, 2);
    }
    let g = grid_for(&mesh, 0.4, [0.1; 3]);
    let clipped = clip_solid(&g, &mesh).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(
        rasterise_solid(&g, SolidId(0), &clipped, 2),
        Err(AnalysisError::DegenerateGeometry {
            solid: Some(SolidId(0)),
            detail: Degeneracy::InsideOut
        })
    );
}

/// FNV-1a of every moment's bits for the icosphere at h = 0.6. A change
/// means the rasteriser's arithmetic order changed.
const PINNED_MOMENTS_HASH: u64 = 17801455811740660721;

#[test]
fn the_moment_bits_are_pinned() {
    let mesh = icosphere([0.2, 0.1, -0.3], 2.0, 3);
    let g = grid_for(&mesh, 0.6, [0.0; 3]);
    let r = rasterise(&mesh, &g);
    let flat: Vec<f64> = r.cells.iter().flat_map(|c| c.m).collect();
    assert_eq!(
        bits_hash(&flat),
        bits_hash(&{
            let r2 = rasterise(&mesh, &g);
            r2.cells.iter().flat_map(|c| c.m).collect::<Vec<f64>>()
        })
    );
    assert_eq!(bits_hash(&flat), PINNED_MOMENTS_HASH);
}
