//! The cached re-assembly against `assemble_solid` itself.

use super::*;
use crate::element::vector_boundary::PieceSupport;
use crate::mesh_test_helpers::{box_mesh, model_of};
use crate::options::AnalysisOptions;
use crate::structural::assemble::{FaceMechanics, assemble_solid};
use crate::structural::generative::{Evaluator, GenerativeOptions};

/// A clamped, loaded block on a shifted grid (cut cells and fragments
/// on every face), scaled entry by entry: the cached system equals
/// `assemble_solid`'s to the bit, matrices and load.
#[test]
fn the_cached_system_is_assemble_solids_to_the_bit() {
    let model = model_of(box_mesh([0.0; 3], [4.0, 2.0, 1.5]), 6);
    let mut faces = vec![FaceMechanics::Free; 6];
    faces[0] = FaceMechanics::Support(PieceSupport::Fixed {
        u0: [0.0; 3],
        gradient: [[0.0; 3]; 3],
    });
    faces[2] = FaceMechanics::Support(PieceSupport::Roller);
    faces[1] = FaceMechanics::Traction([0.0, -1.0, 0.3]);
    faces[5] = FaceMechanics::Pressure(0.2);
    let problem = ElasticProblem {
        youngs: 2000.0,
        poisson: 0.3,
        faces,
    };
    let options = AnalysisOptions {
        target_cell_size: Some(0.5),
        grid_shift: [0.13, 0.29, 0.07],
        ..AnalysisOptions::default()
    };
    let ev = Evaluator::new(
        &model,
        &problem,
        GenerativeOptions::new(0.5, 0.75),
        &[],
        &options,
    )
    .unwrap_or_else(|e| panic!("{e}"));
    let prepared = &ev.prepared;
    let scales: Vec<Vec<f64>> = prepared
        .solids
        .iter()
        .map(|(_, _, m, _)| {
            (0..m.cells.len())
                .map(|i| 1e-9 + ((13 * i) % 17) as f64 / 16.0)
                .collect()
        })
        .collect();
    let cached = ScaledSystem::new(prepared, &problem).unwrap_or_else(|e| panic!("{e}"));
    let (matrices, b) = cached
        .assemble(&ev.blocks, &scales)
        .unwrap_or_else(|e| panic!("{e}"));
    let mut expected_b = Vec::new();
    for (s, (_, clipped, moments, aggregation)) in prepared.solids.iter().enumerate() {
        let input = SolidInput {
            grid: &prepared.grid,
            clipped,
            moments,
            aggregation,
            faces: &problem.faces,
            prescribed: None,
            thermal: None,
            cores: None,
            scales: Some(&scales[s]),
        };
        let system = assemble_solid(&prepared.vt, &prepared.bt, &input, &prepared.elastic)
            .unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(system.a.rows(), matrices[s].rows());
        assert!(system.a.nnz() == matrices[s].nnz());
        for r in 0..system.a.rows() {
            let same = system
                .a
                .row(r)
                .zip(matrices[s].row(r))
                .all(|(x, y)| x.0 == y.0 && x.1.to_bits() == y.1.to_bits());
            assert!(same, "solid {s} row {r}");
        }
        expected_b.extend(system.b);
    }
    assert_eq!(b.len(), expected_b.len());
    assert!(
        b.iter()
            .zip(&expected_b)
            .all(|(x, y)| x.to_bits() == y.to_bits())
    );
}
