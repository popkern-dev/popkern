//! Vector gates on STEP fixtures (PH-025): LF1 (a linear displacement
//! field on every face is reproduced at five grid shifts) and E5 (the
//! print-in-place hinge's leaves stay separate bodies).

use std::path::Path;

use popkern_analysis::model::{AnalysisModel, FaceInfo, SolidId, Surface};
use popkern_analysis::structural::assemble::FaceMechanics;
use popkern_analysis::structural::{ElasticProblem, StructuralResult, solve_elastic};
use popkern_analysis::{AnalysisOptions, element::vector_boundary::PieceSupport};
use popkern_analysis_cli::step_model::model_from_step_file;

fn fixture(name: &str) -> AnalysisModel {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join(format!("../../tests/fixtures/analysis/{name}.step"));
    model_from_step_file(&path).unwrap_or_else(|e| panic!("{name}: {e}"))
}

#[test]
fn lf1_a_linear_displacement_field_is_reproduced_at_five_shifts() {
    let g = [
        [1e-3, 2e-4, -3e-4],
        [5e-4, -8e-4, 1e-4],
        [-2e-4, 6e-4, 4e-4],
    ];
    let u0 = [1e-3, -2e-3, 5e-4];
    let field = |x: &[f64; 3]| -> [f64; 3] {
        std::array::from_fn(|i| u0[i] + (0..3).map(|j| g[i][j] * x[j]).sum::<f64>())
    };
    for name in ["block_two_bores", "block_dimple"] {
        let model = fixture(name);
        let problem = ElasticProblem {
            youngs: 2000.0,
            poisson: 0.35,
            faces: vec![
                FaceMechanics::Support(PieceSupport::Fixed { u0, gradient: g });
                model.faces.len()
            ],
        };
        let mut fingerprints = Vec::new();
        for shift in [
            [0.0; 3],
            [0.21, 0.0, 0.0],
            [0.0, 0.37, 0.11],
            [0.43, 0.19, 0.31],
            [0.61; 3],
        ] {
            let options = AnalysisOptions {
                target_cell_size: Some(0.8),
                grid_shift: shift,
                cg_rel_tol: 1e-12,
                coarsest_dofs: 600,
                ..AnalysisOptions::default()
            };
            let r =
                solve_elastic(&model, &problem, &options).unwrap_or_else(|e| panic!("{name}: {e}"));
            let mut worst = 0.0f64;
            for (v, d) in model.solids[0]
                .vertices
                .iter()
                .zip(&r.vertex_displacement[0])
            {
                let exact = field(v);
                for i in 0..3 {
                    worst = worst.max((d[i] - exact[i]).abs());
                }
            }
            // |u| reaches about 1.5e-2 mm over the 10 mm block.
            assert!(worst <= 1e-8 * 1.5e-2, "{name} shift {shift:?}: {worst}");
            fingerprints.push((
                r.report.dof,
                r.report.cells_cut,
                r.report.residual_history.last().map(|v| v.to_bits()),
            ));
        }
        fingerprints.sort_unstable();
        fingerprints.dedup();
        assert!(
            fingerprints.len() >= 4,
            "{name}: the shifts must change the grid"
        );
    }
}

/// Plane faces of `solid` whose outward normal is close to `n`.
fn planes(m: &AnalysisModel, solid: u32, n: [f64; 3]) -> Vec<&FaceInfo> {
    m.faces
        .iter()
        .filter(|f| f.solid == SolidId(solid))
        .filter(|f| matches!(f.surface, Surface::Plane { normal, .. } if (0..3).map(|k| normal[k] * n[k]).sum::<f64>() > 0.9999))
        .collect()
}

/// The far-end face of a leaf: the plane with normal `n` farthest out.
fn far_end(m: &AnalysisModel, solid: u32, n: [f64; 3]) -> usize {
    planes(m, solid, n)
        .into_iter()
        .max_by(|a, b| {
            let off = |f: &FaceInfo| {
                if let Surface::Plane { offset, .. } = f.surface {
                    offset
                } else {
                    f64::NEG_INFINITY
                }
            };
            off(a).total_cmp(&off(b))
        })
        .map(|f| f.id.0 as usize)
        .unwrap_or_else(|| panic!("no far end on solid {solid}"))
}

/// Leaf A clamped at its far end and pushed on its +y faces; leaf B
/// (when present) clamped at its far end.
fn hinge_problem(m: &AnalysisModel) -> ElasticProblem {
    let mut faces = vec![FaceMechanics::Free; m.faces.len()];
    let clamp = FaceMechanics::Support(PieceSupport::Fixed {
        u0: [0.0; 3],
        gradient: [[0.0; 3]; 3],
    });
    faces[far_end(m, 0, [-1.0, 0.0, 0.0])] = clamp;
    for f in planes(m, 0, [0.0, 1.0, 0.0]) {
        faces[f.id.0 as usize] = FaceMechanics::Traction([0.0, -1e-3, 0.0]);
    }
    if m.solids.len() > 1 {
        faces[far_end(m, 1, [1.0, 0.0, 0.0])] = clamp;
    }
    ElasticProblem {
        youngs: 1.0,
        poisson: 0.3,
        faces,
    }
}

fn solve(m: &AnalysisModel) -> StructuralResult {
    solve_with(m, &hinge_problem(m))
}

fn solve_with(m: &AnalysisModel, problem: &ElasticProblem) -> StructuralResult {
    let options = AnalysisOptions {
        target_cell_size: Some(1.0),
        cg_rel_tol: 1e-11,
        cg_max_iters: 400,
        coarsest_dofs: 600,
        ..AnalysisOptions::default()
    };
    solve_elastic(m, problem, &options).unwrap_or_else(|e| panic!("{e}"))
}

#[test]
fn e5_the_hinge_leaves_stay_separate_bodies() {
    let hinge = fixture("hinge_print_in_place_60");
    assert_eq!(hinge.solids.len(), 2);
    let both = solve(&hinge);
    // Leaf A alone: its solid and its faces (numbered first).
    let mut alone = hinge.clone();
    alone.solids.truncate(1);
    alone.faces.retain(|f| f.solid == SolidId(0));
    let single = solve(&alone);
    let scale = single.max_displacement;
    let mut worst = 0.0f64;
    for (a, b) in both.vertex_displacement[0]
        .iter()
        .zip(&single.vertex_displacement[0])
    {
        for i in 0..3 {
            worst = worst.max((a[i] - b[i]).abs());
        }
    }
    assert!(
        worst <= 1e-6 * scale,
        "leaf A differs from leaf A alone by {worst} of {scale}"
    );
    // Stresses too (scope moved from PH-025 to PH-026).
    let peak = single.stresses[0]
        .vertices
        .iter()
        .map(|v| v.von_mises)
        .fold(0.0, f64::max);
    let stress_diff = both.stresses[0]
        .vertices
        .iter()
        .zip(&single.stresses[0].vertices)
        .map(|(a, b)| (a.von_mises - b.von_mises).abs())
        .fold(0.0, f64::max);
    assert!(
        stress_diff <= 1e-6 * peak,
        "leaf A stress differs by {stress_diff} of {peak}"
    );
    assert!(
        both.report.cg_iterations as f64 <= 1.5 * single.report.cg_iterations as f64,
        "{} vs {}",
        both.report.cg_iterations,
        single.report.cg_iterations
    );
    // Negative control: the two leaves as one solid share the cells
    // across the 0.3 mm clearance and weld. Same face conditions (face
    // ids are unchanged by the merge); leaf A's vertices come first.
    let problem = hinge_problem(&hinge);
    let leaf_a_vertices = hinge.solids[0].vertices.len();
    let mut welded = hinge.clone();
    let offset = welded.solids[0].vertices.len() as u32;
    let b = welded.solids.remove(1);
    let a = &mut welded.solids[0];
    a.vertices.extend(b.vertices);
    a.triangles
        .extend(b.triangles.iter().map(|t| t.map(|i| i + offset)));
    a.triangle_face.extend(b.triangle_face);
    a.exact_volume = None;
    for f in &mut welded.faces {
        f.solid = SolidId(0);
    }
    let glued = solve_with(&welded, &problem);
    let tip_a = |r: &StructuralResult| {
        r.vertex_displacement[0][..leaf_a_vertices]
            .iter()
            .map(|d| d[1].abs())
            .fold(0.0, f64::max)
    };
    let (free, stuck) = (tip_a(&single), tip_a(&glued));
    assert!(
        stuck < 0.9 * free,
        "welding should stiffen leaf A: {stuck} vs {free}"
    );
}
