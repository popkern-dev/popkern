//! Gate LF2 on the frozen STEP fixtures: a linear temperature imposed
//! (Nitsche) on every face is reproduced at every vertex, at five grid
//! offsets, with heat balance to 1e-9.

use std::path::Path;

use popkern_analysis::AnalysisOptions;
use popkern_analysis::thermal::{ConductionProblem, FaceCondition, solve_conduction};
use popkern_analysis_cli::step_model::model_from_step_file;

#[test]
fn lf2_on_both_fixtures_at_five_offsets() {
    let gradient = [0.4, -0.25, 0.9];
    for name in ["block_two_bores", "block_dimple"] {
        let path = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join(format!("../../tests/fixtures/analysis/{name}.step"));
        let base = model_from_step_file(&path).unwrap_or_else(|e| panic!("{name}: {e}"));
        let mut dofs = Vec::new();
        for shift in [
            [0.0; 3],
            [0.21, 0.0, 0.0],
            [0.0, 0.37, 0.11],
            [0.43, 0.19, 0.31],
            [0.61; 3],
        ] {
            let model = base.clone();
            let problem = ConductionProblem {
                conductivity: 100.0,
                faces: vec![
                    FaceCondition::Fixed {
                        t0_c: 20.0,
                        gradient
                    };
                    model.faces.len()
                ],
                source_w: 0.0,
                hdt_c: None,
                body_sources: Vec::new(),
            };
            let options = AnalysisOptions {
                target_cell_size: Some(0.8),
                grid_shift: shift,
                cg_rel_tol: 1e-12,
                coarsest_dofs: 500,
                ..AnalysisOptions::default()
            };
            let r = solve_conduction(&model, &problem, &options)
                .unwrap_or_else(|e| panic!("{name}: {e}"));
            let mut worst = 0.0f64;
            for (v, t) in model.solids[0]
                .vertices
                .iter()
                .zip(&r.vertex_temperature[0])
            {
                let exact = 20.0 + gradient[0] * v[0] + gradient[1] * v[1] + gradient[2] * v[2];
                worst = worst.max((t - exact).abs());
            }
            // The field spans about 15 °C over the 10 mm block.
            assert!(worst <= 1e-8 * 15.0, "{name} shift {shift:?}: {worst}");
            assert!(r.heat.imbalance().abs() <= 1e-9 * r.heat.throughput());
            // Identical discretisations solve bit-identically, so the final
            // residual's bits fingerprint the cut pattern.
            dofs.push((
                r.report.dof,
                r.report.cells_cut,
                r.report.residual_history.last().map(|v| v.to_bits()),
            ));
        }
        dofs.sort_unstable();
        dofs.dedup();
        assert!(
            dofs.len() >= 4,
            "{name}: the shifts must change the grid: {dofs:?}"
        );
    }
}
