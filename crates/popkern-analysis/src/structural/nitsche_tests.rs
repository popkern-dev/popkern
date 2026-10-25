//! Nitsche supports near incompressibility (PH-034): the flux from the
//! enhanced stress keeps the system definite at ν = 0.49.

use super::assemble::{Elastic, FaceMechanics, SolidInput, assemble_solid};
use super::thin_wall_tests::{rotate_vector, tip};
use super::*;
use crate::aggregate::aggregate;
use crate::clip::clip_solid;
use crate::element::VectorTables;
use crate::element::vector::lame;
use crate::element::vector_boundary::{BoundaryTables, PieceSupport};
use crate::linalg::sparse::Csr;
use crate::mesh_test_helpers::{box_mesh, model_of, rotated};

/// Number of negative Ritz values after `steps` Lanczos steps (full
/// reorthogonalisation) from a fixed start: any is proof that `a` is
/// indefinite (Ritz values lie inside the spectrum).
fn negative_ritz(a: &Csr, steps: usize) -> usize {
    let n = a.rows();
    let mut q: Vec<Vec<f64>> = Vec::new();
    let mut v: Vec<f64> = (0..n)
        .map(|i| 1.0 + ((i * 7919) % 13) as f64 / 13.0)
        .collect();
    let norm = |x: &[f64]| x.iter().map(|y| y * y).sum::<f64>().sqrt();
    let nv = norm(&v);
    v.iter_mut().for_each(|x| *x /= nv);
    let (mut alpha, mut beta) = (Vec::new(), Vec::new());
    for _ in 0..steps.min(n) {
        let mut w = vec![0.0; n];
        if a.mul_vec(&v, &mut w).is_err() {
            break;
        }
        let al: f64 = w.iter().zip(&v).map(|(x, y)| x * y).sum();
        q.push(v.clone());
        for _ in 0..2 {
            for qk in &q {
                let c: f64 = w.iter().zip(qk).map(|(x, y)| x * y).sum();
                w.iter_mut().zip(qk).for_each(|(x, y)| *x -= c * y);
            }
        }
        alpha.push(al);
        let b = norm(&w);
        if b < 1e-12 * al.abs().max(1e-300) {
            break;
        }
        beta.push(b);
        v = w.iter().map(|x| x / b).collect();
    }
    // Sturm count of eigenvalues of T below 0.
    let mut count = 0;
    let mut d = 1.0f64;
    for k in 0..alpha.len() {
        let off = if k == 0 {
            0.0
        } else {
            beta[k - 1] * beta[k - 1]
        };
        d = alpha[k] - off / if d == 0.0 { 1e-300 } else { d };
        if d < 0.0 {
            count += 1;
        }
    }
    count
}

/// The reduced system of the 10 × 1 × 1 cantilever at `cells` through
/// the depth, rotated by `degrees`, with or without the enhanced flux.
fn cantilever_system(cells: f64, degrees: f64, nu: f64, enhanced_flux: bool) -> Csr {
    let centre = [0.0, 0.5, 0.5];
    let model = model_of(
        rotated(box_mesh([0.0; 3], [10.0, 1.0, 1.0]), degrees, centre),
        6,
    );
    let mut faces = vec![FaceMechanics::Free; 6];
    faces[0] = FaceMechanics::Support(PieceSupport::Fixed {
        u0: [0.0; 3],
        gradient: [[0.0; 3]; 3],
    });
    let options = AnalysisOptions {
        target_cell_size: Some(1.0 / cells),
        ..AnalysisOptions::default()
    };
    let solid = &model.solids[0];
    let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    for v in &solid.vertices {
        for a in 0..3 {
            lo[a] = lo[a].min(v[a]);
            hi[a] = hi[a].max(v[a]);
        }
    }
    let grid = Grid::for_model(lo, hi, &options).unwrap_or_else(|e| panic!("{e}"));
    let h = grid.cell_size();
    let clipped = clip_solid(&grid, solid).unwrap_or_else(|e| panic!("{e}"));
    let moments = crate::moments::rasterise_solid(&grid, solid.id, &clipped, 2)
        .unwrap_or_else(|e| panic!("{e}"));
    let area = clipped.pieces.iter().map(|p| p.area).sum();
    let aggregation = aggregate(&grid, &moments, 0.25, area).unwrap_or_else(|e| panic!("{e:?}"));
    let (lambda, mu) = lame(1.0, nu);
    let elastic = Elastic {
        lambda,
        mu,
        beta: options.nitsche_factor * (lambda + 2.0 * mu) / h,
        enhanced_flux,
        core: None,
    };
    let input = SolidInput {
        grid: &grid,
        clipped: &clipped,
        moments: &moments,
        aggregation: &aggregation,
        faces: &faces,
        prescribed: None,
        thermal: None,
        cores: None,
        scales: None,
    };
    assemble_solid(
        &VectorTables::new(),
        &BoundaryTables::new(),
        &input,
        &elastic,
    )
    .unwrap_or_else(|e| panic!("{e}"))
    .a
}

#[test]
fn nu_049_nitsche_stays_definite_and_accurate_with_the_enhanced_flux() {
    // Definiteness: no negative Ritz value with the enhanced flux; the
    // plain Q1 flux (negative control) has some in 5 of these 8 cases
    // (measured: aligned 1 cell; rotated 1, 1.5, 2, 3 cells).
    let mut q1_indefinite = 0;
    for degrees in [0.0, 17.0] {
        for cells in [1.0, 1.5, 2.0, 3.0] {
            let enhanced = negative_ritz(&cantilever_system(cells, degrees, 0.49, true), 300);
            assert_eq!(enhanced, 0, "{degrees}° at {cells} cells");
            if negative_ritz(&cantilever_system(cells, degrees, 0.49, false), 300) > 0 {
                q1_indefinite += 1;
            }
        }
    }
    assert!(
        q1_indefinite >= 1,
        "the Q1 flux control never went indefinite"
    );
    // Accuracy at 2 cells against the 8-cell solution of the same
    // problem (measured −2.77% aligned, −6.53% rotated).
    for (degrees, tol) in [(0.0, 0.05), (17.0, 0.08)] {
        let reference = tip(8.0, 0.49, degrees);
        let e = tip(2.0, 0.49, degrees) / reference - 1.0;
        assert!(e.abs() <= tol, "{degrees}°: {e}");
    }
    let _ = rotate_vector;
}
