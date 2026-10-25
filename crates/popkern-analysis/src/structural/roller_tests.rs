//! Directional rollers (PH-039): u·d = 0 for a direction d other than
//! the face normal.

use super::assemble::FaceMechanics;
use super::*;
use crate::element::vector_boundary::PieceSupport;
use crate::mesh_test_helpers::{box_mesh, model_of};

/// A 4 × 1 × 1 bar under the uniform stress σ_xx = a, σ_xy = a·tan 30°,
/// whose traction on the end x = 0 is parallel to d = (cos 30°, sin 30°,
/// 0): that end sits on a roller along d, the base z = 0 is held at the
/// exact field, and the other sides carry the exact tractions. Returns
/// the worst nodal error over |u|.
fn bar_on_a_tilted_roller(end: PieceSupport, shift: [f64; 3]) -> f64 {
    let (e, nu, a, l) = (1.0, 0.3, 1e-3, 4.0);
    let (sin, cos) = (0.5, 0.75f64.sqrt());
    let tan = sin / cos;
    let exy = (1.0 + nu) / e * a * tan;
    let eyy = -nu * a / e;
    // A rotation w about z keeps u·d = 0 on x = 0:
    // cos (ε_xy − w) + sin ε_yy = 0.
    let w = exy + tan * eyy;
    let g = [[a / e, exy - w, 0.0], [exy + w, eyy, 0.0], [0.0, 0.0, eyy]];
    let model = model_of(box_mesh([0.0; 3], [l, 1.0, 1.0]), 6);
    let mut faces = vec![FaceMechanics::Free; 6];
    faces[0] = FaceMechanics::Support(end);
    faces[1] = FaceMechanics::Traction([a, a * tan, 0.0]);
    faces[2] = FaceMechanics::Traction([-a * tan, 0.0, 0.0]);
    faces[3] = FaceMechanics::Traction([a * tan, 0.0, 0.0]);
    faces[4] = FaceMechanics::Support(PieceSupport::Fixed {
        u0: [0.0; 3],
        gradient: g,
    });
    let problem = ElasticProblem {
        youngs: e,
        poisson: nu,
        faces,
    };
    let options = AnalysisOptions {
        target_cell_size: Some(0.25),
        grid_shift: shift,
        cg_rel_tol: 1e-12,
        cg_max_iters: 400,
        coarsest_dofs: 600,
        ..AnalysisOptions::default()
    };
    let r = solve_elastic(&model, &problem, &options).unwrap_or_else(|e| panic!("{e}"));
    let mut worst = 0.0f64;
    for (v, d) in model.solids[0]
        .vertices
        .iter()
        .zip(&r.vertex_displacement[0])
    {
        for i in 0..3 {
            let exact: f64 = (0..3).map(|j| g[i][j] * v[j]).sum();
            worst = worst.max((d[i] - exact).abs());
        }
    }
    worst / (a / e * l)
}

#[test]
fn a_roller_30_degrees_off_the_normal_reproduces_the_exact_field() {
    let (sin, cos) = (0.5, 0.75f64.sqrt());
    let tilted = PieceSupport::Directional {
        direction: [cos, sin, 0.0],
    };
    // Measured: 1.3e-12, 1.4e-12, 3.7e-12; the control 0.51-0.52.
    for shift in [[0.0; 3], [0.07, 0.11, 0.03], [0.13, 0.2, 0.09]] {
        let err = bar_on_a_tilted_roller(tilted, shift);
        assert!(err <= 1e-8, "{shift:?}: {err}");
        // Negative control: the face-normal roller forbids the field's
        // u_x on x = 0.
        let control = bar_on_a_tilted_roller(PieceSupport::Roller, shift);
        assert!(control > 1e-3, "{shift:?}: control {control}");
    }
}
