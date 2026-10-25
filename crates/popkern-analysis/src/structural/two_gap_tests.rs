//! PH-066: a fin thinner than the lattice spacing between a slot and a
//! second, narrower gap inside one cell row.

use super::assemble::FaceMechanics;
use super::*;
use crate::element::vector_boundary::PieceSupport;
use crate::mesh_test_helpers::{extrude_faces, model_of};
use crate::moments::Splitting;

const G: f64 = 0.075;
const FIN: f64 = 0.005;
/// A's top: 0.05 into the cell row [2.25, 2.5] at h = 0.25.
const Y1: f64 = 2.3;
/// The second gap: 0.3 lattice spacings (split cells are rasterised at
/// depth 4: h / 16) at h = 0.25.
const GAP: f64 = 0.3 * 0.25 / 16.0;

/// A 6 × 4 × 1 plate (`g2` 0: [`plain`]). Block A (y < `Y1`) hangs from the root x < 1;
/// above it a pocket of width `G` from x = 1 to 5.75, then a fin
/// `FIN` thick from x = 1 + `G` joined to A by the bridge x > 5.75,
/// then a gap of width `g2` from x = 1 to the edge, then B up to the
/// top. Faces: 0/1 z sides, 2 A's bottom (y = 0), 3 the pocket, 4 the
/// fin's top and end and B's bottom, 5 other sides, 6 the top (y = 4).
fn finned(g2: f64) -> crate::model::AnalysisModel {
    let (y2, e) = (Y1 + G, 0.25);
    let (y3, y4) = (y2 + FIN, y2 + FIN + g2);
    if g2 == 0.0 {
        return plain();
    }
    let outline = [
        [0.0, 0.0],
        [6.0, 0.0],
        [6.0, y3],
        [1.0 + G, y3],
        [1.0 + G, y2],
        [6.0 - e, y2],
        [6.0 - e, Y1],
        [1.0, Y1],
        [1.0, y4],
        [6.0, y4],
        [6.0, 4.0],
        [0.0, 4.0],
    ];
    let tris = [
        [0, 1, 6],
        [0, 6, 7],
        [1, 2, 5],
        [1, 5, 6],
        [5, 2, 3],
        [5, 3, 4],
        [0, 7, 8],
        [0, 8, 11],
        [8, 9, 10],
        [8, 10, 11],
    ];
    let mut mesh = extrude_faces(
        &outline,
        &tris,
        0.0,
        1.0,
        &[2, 5, 4, 4, 3, 3, 3, 3, 4, 5, 6, 5],
    );
    let empty = (5.0 - e) * G + G * FIN + 5.0 * g2;
    mesh.exact_volume = Some(24.0 - empty);
    model_of(mesh, 7)
}

/// The same plate without the fin and bridge: one slot from A's top to
/// B's bottom, so the same A and B (the fin adds next to nothing to A).
fn plain() -> crate::model::AnalysisModel {
    let y4 = Y1 + G + FIN + GAP;
    let outline = [
        [0.0, 0.0],
        [6.0, 0.0],
        [6.0, Y1],
        [1.0, Y1],
        [1.0, y4],
        [6.0, y4],
        [6.0, 4.0],
        [0.0, 4.0],
    ];
    let tris = [
        [0, 1, 2],
        [0, 2, 3],
        [0, 3, 4],
        [0, 4, 7],
        [4, 5, 6],
        [4, 6, 7],
    ];
    let mut mesh = extrude_faces(&outline, &tris, 0.0, 1.0, &[2, 5, 3, 3, 4, 5, 6, 5]);
    mesh.exact_volume = Some(24.0 - 5.0 * (y4 - Y1));
    model_of(mesh, 7)
}

/// A's mean downward displacement at its bottom when pulled down with
/// the top clamped, at cell size `h`.
fn sag(g2: f64, h: f64, split: Splitting) -> (f64, StructuralResult) {
    let model = finned(g2);
    let mut faces = vec![FaceMechanics::Free; 7];
    faces[6] = FaceMechanics::Support(PieceSupport::Fixed {
        u0: [0.0; 3],
        gradient: [[0.0; 3]; 3],
    });
    faces[2] = FaceMechanics::Traction([0.0, -1e-3, 0.0]);
    let problem = ElasticProblem {
        youngs: 1.0,
        poisson: 0.3,
        faces,
    };
    let options = AnalysisOptions {
        target_cell_size: Some(h),
        max_active_cells: 400_000,
        coarsest_dofs: 600,
        cg_max_iters: 600,
        ..AnalysisOptions::default()
    };
    let r = solve_elastic_with(
        &model,
        &problem,
        &options,
        crate::recover::Fallback::InCell,
        None,
        split,
    )
    .unwrap_or_else(|e| panic!("{e}"));
    let points = [1.5, 2.5, 3.5, 4.5];
    let sum: f64 = points
        .iter()
        .map(|&x| {
            r.field
                .displacement_at([x, 0.05, 0.5])
                .unwrap_or_else(|e| panic!("{e}"))
                .unwrap_or_else(|| panic!("outside"))[1]
        })
        .sum();
    (-sum / points.len() as f64, r)
}

/// Measured 2026-10-11, A's sag over the plain plate's at h = 0.1: the
/// fin plate at h = 0.25 reads 0.907 (0.037 before, when both the fin's
/// faces counted B's lattice points across the gap as support, and 0.057
/// with only that fixed, when the fin's material boxes stepped back
/// through the gap into B); the plain plate 0.907; welded 0.007.
#[test]
fn a_fin_between_a_slot_and_a_second_gap_does_not_weld() {
    let reference = sag(0.0, 0.1, Splitting::Exact).0;
    let (fin, r) = sag(GAP, 0.25, Splitting::Exact);
    let fin = fin / reference;
    let plain = sag(0.0, 0.25, Splitting::Exact).0 / reference;
    assert!(r.report.split_cells > 0);
    // Every split cell on exact moments (72 on lattice moments before).
    assert_eq!(r.report.split_cells_lattice_moments, 0);
    assert!((fin - plain).abs() <= 0.05, "{fin} vs {plain}");
    // Negative control: welded, A hangs from B too.
    let welded = sag(GAP, 0.25, Splitting::Welded).0 / reference;
    assert!(welded < 0.1, "welded {welded}");
}

/// Along the fin each cell holds A, the fin and B as three fragments with
/// exact volumes. Before, with a gap of 1.2 spacings, the fin's two faces
/// were two sheets meeting only on the cell walls: two slivers with no
/// volume, and the cell fell back to lattice moments.
#[test]
fn the_fin_is_its_own_fragment_with_exact_volumes() {
    use crate::clip::clip_solid;
    use crate::grid::Grid;
    use crate::moments::rasterise_with;
    for gap in [1.2 * 0.25 / 16.0, GAP] {
        let model = finned(gap);
        let options = AnalysisOptions {
            target_cell_size: Some(0.25),
            ..AnalysisOptions::default()
        };
        let grid =
            Grid::for_model([0.0; 3], [6.0, 4.0, 1.0], &options).unwrap_or_else(|e| panic!("{e}"));
        let mesh = &model.solids[0];
        let clipped = clip_solid(&grid, mesh).unwrap_or_else(|e| panic!("{e}"));
        let m = rasterise_with(&grid, mesh.id, &clipped, 2, &|_| false, Splitting::Exact)
            .unwrap_or_else(|e| panic!("{e}"));
        let b = 0.25 - (Y1 - 2.25) - G - FIN - gap;
        // The fin's end cell (it starts at x = 1 + G) and two along it, in
        // the bottom layer (the fin's z = 0 face on the cell wall) and
        // inside.
        let along = [(1.1, 0.25 - G), (3.1, 0.25), (5.6, 0.25)];
        for (x, fin_length, z) in [0.1, 0.4, 0.6]
            .map(|z| along.map(|(x, l)| (x, l, z)))
            .concat()
        {
            let cell = grid
                .cell_of_point([x, 2.3, z])
                .and_then(|ijk| grid.cell_index(ijk))
                .unwrap_or_else(|| panic!("fin cell"));
            let sc = m
                .fragments
                .split
                .iter()
                .find(|s| s.cell == cell)
                .unwrap_or_else(|| panic!("fin cell not split"));
            assert_eq!(sc.slivers.len(), 1, "x {x} z {z}");
            let volumes: Vec<f64> = m.entries(cell).map(|e| m.cells[e].m[0]).collect();
            let exact = [
                (Y1 - 2.25) * 0.25 * 0.25,
                b * 0.25 * 0.25,
                FIN * fin_length * 0.25,
            ];
            assert_eq!(volumes.len(), 3, "x {x}: {volumes:?}");
            for (v, e) in volumes.iter().zip(exact) {
                assert!((v - e).abs() <= 1e-12, "x {x}: {volumes:?} vs {exact:?}");
            }
        }
    }
}

/// PH-067: a linear temperature held on every face is exact on the fin
/// plate (to 1e-9 of its largest value). The fin's fragments are small
/// all along and extend a root up to 19 cells away; with β = factor·k/h
/// alone the system was indefinite (the solve declined: a negative
/// pivot) at both gaps. With β raised there (about 40×), the pocket's end
/// on the grid plane x = 5.75 then read 9.3e-5: A beside the pocket and
/// the bridge beyond took different dof copies at that face.
#[test]
fn a_linear_temperature_is_exact_on_the_fin() {
    use crate::thermal::{ConductionProblem, FaceCondition, solve_conduction};
    let g = [0.7, -0.3, 0.45];
    let problem = ConductionProblem {
        conductivity: 20.0,
        faces: vec![
            FaceCondition::Fixed {
                t0_c: 1.0,
                gradient: g,
            };
            7
        ],
        source_w: 0.0,
        body_sources: Vec::new(),
        hdt_c: None,
    };
    let options = AnalysisOptions {
        target_cell_size: Some(0.25),
        cg_rel_tol: 1e-12,
        coarsest_dofs: 600,
        ..AnalysisOptions::default()
    };
    for gap in [GAP, 1.2 * 0.25 / 16.0] {
        let model = finned(gap);
        let r = solve_conduction(&model, &problem, &options)
            .unwrap_or_else(|e| panic!("gap {gap}: {e}"));
        assert_eq!(r.report.split_cells_lattice_moments, 0);
        let exact = |p: &[f64; 3]| 1.0 + (0..3).map(|a| g[a] * p[a]).sum::<f64>();
        let vertices = &model.solids[0].vertices;
        let scale = vertices.iter().map(|p| exact(p).abs()).fold(0.0, f64::max);
        let worst = vertices
            .iter()
            .zip(&r.vertex_temperature[0])
            .map(|(p, t)| (t - exact(p)).abs())
            .fold(0.0, f64::max);
        assert!(worst <= 1e-9 * scale, "gap {gap}: {worst} of {scale}");
        let heat = &r.heat;
        assert!(heat.imbalance().abs() <= 1e-9 * heat.throughput());
    }
}

/// PH-068: LF1 on the fin plate, a linear displacement held on every
/// face, exact to 1e-8 of its largest value. With β = γ_D(λ + 2μ)/h
/// alone the system was indefinite at both gaps (pivot −8.4e-3 of
/// diagonal 756).
#[test]
fn a_linear_displacement_is_exact_on_the_fin() {
    let g = [[1e-3, 2e-3, 0.0], [-1e-3, 5e-4, 1e-3], [0.0, 1e-3, -2e-3]];
    let exact = |p: &[f64; 3]| -> [f64; 3] {
        std::array::from_fn(|i| (0..3).map(|j| g[i][j] * p[j]).sum::<f64>())
    };
    let options = AnalysisOptions {
        target_cell_size: Some(0.25),
        cg_rel_tol: 1e-12,
        coarsest_dofs: 600,
        ..AnalysisOptions::default()
    };
    for gap in [GAP, 1.2 * 0.25 / 16.0] {
        let model = finned(gap);
        let support = FaceMechanics::Support(PieceSupport::Fixed {
            u0: [0.0; 3],
            gradient: g,
        });
        let problem = ElasticProblem {
            youngs: 1.0,
            poisson: 0.3,
            faces: vec![support; 7],
        };
        let r = solve_elastic_with(
            &model,
            &problem,
            &options,
            crate::recover::Fallback::InCell,
            None,
            Splitting::Exact,
        )
        .unwrap_or_else(|e| panic!("gap {gap}: {e}"));
        let vertices = &model.solids[0].vertices;
        let norm = |u: [f64; 3]| (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt();
        let scale = vertices.iter().map(|p| norm(exact(p))).fold(0.0, f64::max);
        let worst = vertices
            .iter()
            .map(|p| {
                let u = r
                    .field
                    .displacement_at(*p)
                    .ok()
                    .flatten()
                    .unwrap_or_else(|| panic!("no displacement at {p:?}"));
                let e = exact(p);
                norm(std::array::from_fn(|i| u[i] - e[i]))
            })
            .fold(0.0, f64::max);
        assert!(worst <= 1e-8 * scale, "gap {gap}: {worst} of {scale}");
    }
}
