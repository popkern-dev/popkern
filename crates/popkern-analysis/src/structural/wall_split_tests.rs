//! Exact fragment moments where no cell wall has a single owner
//! (PH-045): walls no fragment reaches (two fins inside one cell) and
//! walls every fragment's section must be split in the wall's plane
//! (the centre of a cross-shaped slot).

use super::assemble::FaceMechanics;
use super::slit_tests::box_moments;
use super::*;
use crate::element::vector_boundary::PieceSupport;
use crate::mesh_test_helpers::{extrude, model_of};
use crate::moments::Splitting;

pub(super) const H: f64 = 0.25;
const T: f64 = 1.0;

pub(super) const SHIFTS: [[f64; 3]; 5] = [
    [0.0; 3],
    [0.0, 0.125, 0.0],
    [0.07, 0.11, 0.03],
    [0.0, 0.05, 0.0],
    [0.13, 0.2, 0.09],
];

/// Triangles of the simple counter-clockwise polygon `p` by ear
/// clipping (a vertex on a candidate ear blocks it).
fn triangulate(p: &[[f64; 2]]) -> Vec<[u32; 3]> {
    let cross = |o: [f64; 2], a: [f64; 2], b: [f64; 2]| {
        (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
    };
    let mut left: Vec<usize> = (0..p.len()).collect();
    let mut out = Vec::new();
    while left.len() > 3 {
        let n = left.len();
        let ear = (0..n).find(|&i| {
            let (a, b, c) = (left[(i + n - 1) % n], left[i], left[(i + 1) % n]);
            cross(p[a], p[b], p[c]) > 0.0
                && left.iter().all(|&k| {
                    k == a
                        || k == b
                        || k == c
                        || cross(p[a], p[b], p[k]) < 0.0
                        || cross(p[b], p[c], p[k]) < 0.0
                        || cross(p[c], p[a], p[k]) < 0.0
                })
        });
        let i = ear.unwrap_or_else(|| panic!("no ear"));
        out.push([left[(i + n - 1) % n], left[i], left[(i + 1) % n]].map(|k| k as u32));
        left.remove(i);
    }
    out.push([left[0], left[1], left[2]].map(|k| k as u32));
    out
}

/// The prism of outline `p` (z from 0 to T), faces 0 bottom, 1 top,
/// 2 sides.
fn prism(p: &[[f64; 2]], area: f64) -> crate::model::AnalysisModel {
    let mut mesh = extrude(p, &triangulate(p), 0.0, T);
    mesh.exact_volume = Some(area * T);
    model_of(mesh, 3)
}

/// Cross slot: width, centre, and the reach of its east, north and
/// south arms (the west arm runs out of the 4 × 4 square).
const W: f64 = 0.05;
const CENTRE: [f64; 2] = [2.19, 1.9125];
const ARM: f64 = 0.6;

/// A 4 × 4 square with a cross-shaped slot whose west arm runs out of
/// it: around the crossing, four quadrants a slot apart.
fn crossed() -> crate::model::AnalysisModel {
    let ([cx, cy], d) = (CENTRE, 0.5 * W);
    let outline = [
        [0.0, 0.0],
        [4.0, 0.0],
        [4.0, 4.0],
        [0.0, 4.0],
        [0.0, cy + d],
        [cx - d, cy + d],
        [cx - d, cy + ARM],
        [cx + d, cy + ARM],
        [cx + d, cy + d],
        [cx + ARM, cy + d],
        [cx + ARM, cy - d],
        [cx + d, cy - d],
        [cx + d, cy - ARM],
        [cx - d, cy - ARM],
        [cx - d, cy - d],
        [0.0, cy - d],
    ];
    prism(&outline, 16.0 - W * (cx + 3.0 * ARM - W))
}

/// Two fins (0.025 and 0.045 thick, 0.035 apart, from x = 1.14)
/// standing on a base strip: inside one cell column at every shift, so
/// the cells holding both fins have no material on their x walls.
fn finned() -> crate::model::AnalysisModel {
    let (t, u, g, a) = (0.025, 0.045, 0.035, 1.14);
    let outline = [
        [0.0, 0.0],
        [3.0, 0.0],
        [3.0, 0.5],
        [a + t + g + u, 0.5],
        [a + t + g + u, 2.0],
        [a + t + g, 2.0],
        [a + t + g, 0.5],
        [a + t, 0.5],
        [a + t, 2.0],
        [a, 2.0],
        [a, 0.5],
        [0.0, 0.5],
    ];
    prism(&outline, 1.5 + (t + u) * 1.5)
}

/// The fragments of split cells with `fragments` entries inside
/// `region` (lower, upper corner) against the nearest of the boxes
/// `forms` (closed-form moments): the worst error over h³, the
/// fragments checked, and the cells left on lattice moments, over the
/// shifts.
pub(super) fn fragment_errors(
    model: &crate::model::AnalysisModel,
    hi: [f64; 3],
    forms: &[([f64; 3], [f64; 3])],
    fragments: usize,
    region: ([f64; 3], [f64; 3]),
    split: Splitting,
) -> (f64, usize, usize) {
    use crate::clip::clip_solid;
    use crate::moments::rasterise_with;
    let mesh = &model.solids[0];
    let (mut worst, mut checked, mut lattice) = (0.0f64, 0, 0);
    for shift in SHIFTS {
        let options = AnalysisOptions {
            target_cell_size: Some(H),
            grid_shift: shift,
            ..AnalysisOptions::default()
        };
        let grid = Grid::for_model([0.0; 3], hi, &options).unwrap_or_else(|e| panic!("{e}"));
        let clipped = clip_solid(&grid, mesh).unwrap_or_else(|e| panic!("{e}"));
        let m = rasterise_with(&grid, mesh.id, &clipped, 2, &|_| false, split)
            .unwrap_or_else(|e| panic!("{e}"));
        lattice += m.fragments.lattice_moment_cells;
        for sc in &m.fragments.split {
            let entries = m.entries(sc.cell);
            let ijk = grid.cell_ijk(sc.cell).unwrap_or_else(|| panic!("cell"));
            let lo = grid.cell_min(ijk);
            let inside = (0..3).all(|a| lo[a] >= region.0[a] && lo[a] + H <= region.1[a]);
            if entries.len() != fragments || !inside {
                continue;
            }
            let forms: Vec<[f64; 27]> = forms
                .iter()
                .map(|(p, q)| box_moments(*p, *q, lo, H))
                .collect();
            for e in entries {
                let got = &m.cells[e].m;
                let best = forms
                    .iter()
                    .map(|f| {
                        got.iter()
                            .zip(f)
                            .map(|(x, y)| (x - y).abs())
                            .fold(0.0, f64::max)
                    })
                    .fold(f64::INFINITY, f64::min);
                worst = worst.max(best / (H * H * H));
                checked += 1;
            }
        }
    }
    (worst, checked, lattice)
}

#[test]
fn the_quadrants_at_a_slot_crossing_get_exact_moments() {
    // Every wall of the crossing's cells is reached by two or more
    // quadrants: the sections split in the wall's plane. Measured
    // 2026-10-10: exact 5.6e-16 h³ over 80 fragments, lattice moments
    // 1.4e-2 h³.
    let ([cx, cy], d) = (CENTRE, 0.5 * W);
    let quadrants = [
        ([-9.0, -9.0, 0.0], [cx - d, cy - d, T]),
        ([cx + d, -9.0, 0.0], [9.0, cy - d, T]),
        ([-9.0, cy + d, 0.0], [cx - d, 9.0, T]),
        ([cx + d, cy + d, 0.0], [9.0, 9.0, T]),
    ];
    let everywhere = ([-9.0; 3], [9.0; 3]);
    let model = crossed();
    let hi = [4.0, 4.0, T];
    let (exact, checked, lattice) =
        fragment_errors(&model, hi, &quadrants, 4, everywhere, Splitting::Exact);
    assert_eq!(lattice, 0);
    assert!(exact <= 1e-12, "{exact}");
    // Four z layers of four quadrants at each of five shifts.
    assert!(checked >= 16 * 5, "{checked} quadrant fragments");
    let control = fragment_errors(
        &model,
        hi,
        &quadrants,
        4,
        everywhere,
        Splitting::LatticeMoments,
    );
    assert!(control.0 > 1e-6, "control {}", control.0);
}

#[test]
fn fins_inside_one_cell_column_get_exact_moments() {
    // The x walls of the cells holding both fins are reached by neither.
    // Measured 2026-10-10: exact 2.2e-16 h³ over 210 fragments, lattice
    // moments 1.2e-2 h³ (LF1 cannot tell them apart here: 3e-11 either
    // way, the fins being clamped on every face).
    let (t, u, g, a) = (0.025, 0.045, 0.035, 1.14);
    let fins = [
        ([a, 0.5, 0.0], [a + t, 2.0, T]),
        ([a + t + g, 0.5, 0.0], [a + t + g + u, 2.0, T]),
    ];
    let above_base = ([-9.0, 0.5, -9.0], [9.0, 2.0, 9.0]);
    let model = finned();
    let hi = [3.0, 2.0, T];
    let (exact, checked, lattice) =
        fragment_errors(&model, hi, &fins, 2, above_base, Splitting::Exact);
    assert_eq!(lattice, 0);
    assert!(exact <= 1e-12, "{exact}");
    assert!(checked >= 40, "{checked} fin fragments");
    let control = fragment_errors(&model, hi, &fins, 2, above_base, Splitting::LatticeMoments);
    assert!(control.0 > 1e-6, "control {}", control.0);
}

/// The worst vertex error of `model` with u = G x prescribed on every
/// face (relative to |G x|) and its cells left on lattice moments.
fn linear_field_error(
    model: &crate::model::AnalysisModel,
    shift: [f64; 3],
    split: Splitting,
) -> (f64, usize, usize) {
    let g = [[1e-3, 2e-3, 0.0], [-1e-3, 5e-4, 1e-3], [0.0, 1e-3, -2e-3]];
    let fixed = FaceMechanics::Support(PieceSupport::Fixed {
        u0: [0.0; 3],
        gradient: g,
    });
    let problem = ElasticProblem {
        youngs: 1.0,
        poisson: 0.3,
        faces: vec![fixed; 3],
    };
    let options = AnalysisOptions {
        target_cell_size: Some(H),
        grid_shift: shift,
        cg_rel_tol: 1e-12,
        cg_max_iters: 400,
        ..AnalysisOptions::default()
    };
    let r = solve_elastic_with(
        model,
        &problem,
        &options,
        crate::recover::Fallback::InCell,
        None,
        split,
    )
    .unwrap_or_else(|e| panic!("{e}"));
    let (mut worst, mut scale) = (0.0f64, 0.0f64);
    for (v, d) in model.solids[0]
        .vertices
        .iter()
        .zip(&r.vertex_displacement[0])
    {
        for i in 0..3 {
            let exact: f64 = (0..3).map(|j| g[i][j] * v[j]).sum();
            worst = worst.max((d[i] - exact).abs());
            scale = scale.max(exact.abs());
        }
    }
    (
        worst / scale,
        r.report.split_cells,
        r.report.split_cells_lattice_moments,
    )
}

#[test]
fn lf1_is_exact_at_a_slot_crossing_and_lattice_moments_are_the_control() {
    // Measured 2026-10-10: exact 1.2e-12 to 3.9e-12 (65-70 split cells
    // at four shifts; none at the fifth), lattice moments 1.5e-4 to
    // 4.3e-4.
    let model = crossed();
    let (mut control_worst, mut split) = (0.0f64, 0);
    for shift in SHIFTS {
        let (exact, cells, lattice) = linear_field_error(&model, shift, Splitting::Exact);
        assert!(exact <= 1e-8, "{shift:?}: {exact}");
        assert_eq!(lattice, 0, "{shift:?}");
        split += cells;
        let control = linear_field_error(&model, shift, Splitting::LatticeMoments).0;
        control_worst = control_worst.max(control);
    }
    assert!(split > 0, "no split cell");
    assert!(control_worst > 1e-6, "control {control_worst}");
}
