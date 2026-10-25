use super::assemble::FaceMechanics;
use super::rigid::{RigidMotion, symmetric_eigen};
use super::*;
use crate::element::vector_boundary::PieceSupport;
use crate::mesh_test_helpers::{box_mesh, icosphere, model_of};
use crate::model::{FaceId, FaceInfo, RoleSource, SolidId, Surface};

fn options(h: f64) -> AnalysisOptions {
    AnalysisOptions {
        target_cell_size: Some(h),
        cg_rel_tol: 1e-12,
        cg_max_iters: 400,
        coarsest_dofs: 600,
        ..AnalysisOptions::default()
    }
}

const ROLLER: FaceMechanics = FaceMechanics::Support(PieceSupport::Roller);

fn clamp() -> FaceMechanics {
    FaceMechanics::Support(PieceSupport::Fixed {
        u0: [0.0; 3],
        gradient: [[0.0; 3]; 3],
    })
}

fn assert_balanced(r: &StructuralResult, scale: f64) {
    let fb = r
        .report
        .force_balance
        .unwrap_or_else(|| panic!("no force balance"));
    for i in 0..3 {
        assert!(
            (fb.applied[i] + fb.reaction[i]).abs() <= 1e-9 * scale,
            "X6 {fb:?}"
        );
    }
}

#[test]
fn c1_a_bar_in_tension_on_symmetry_rollers_is_exact() {
    let (l, e, nu, p) = (10.0, 1000.0, 0.3, 100.0);
    let model = model_of(box_mesh([0.0; 3], [l, 2.0, 2.0]), 6);
    let mut faces = vec![FaceMechanics::Free; 6];
    faces[0] = ROLLER;
    faces[2] = ROLLER;
    faces[4] = ROLLER;
    faces[1] = FaceMechanics::Traction([p / 4.0, 0.0, 0.0]);
    let problem = ElasticProblem {
        youngs: e,
        poisson: nu,
        faces,
    };
    let r = solve_elastic(&model, &problem, &options(0.5)).unwrap_or_else(|e| panic!("{e}"));
    let s = p / 4.0 / e;
    let mut worst = 0.0f64;
    for (v, d) in model.solids[0]
        .vertices
        .iter()
        .zip(&r.vertex_displacement[0])
    {
        let exact = [s * v[0], -nu * s * v[1], -nu * s * v[2]];
        for i in 0..3 {
            worst = worst.max((d[i] - exact[i]).abs());
        }
    }
    assert!(worst <= 1e-8 * s * l, "worst {worst}");
    assert!((r.max_displacement - (s * l).hypot(nu * s * 2.0 * 2f64.sqrt())).abs() < 1e-8 * s * l);
    assert_balanced(&r, p);
    // Load postcondition: the applied total is the requested force.
    let fb = r
        .report
        .force_balance
        .unwrap_or_else(|| panic!("no force balance"));
    assert!(
        (fb.applied[0] - p).abs() <= 1e-12 * p,
        "applied {:?}",
        fb.applied
    );
}

#[test]
fn x7_a_sphere_under_pressure_compresses_uniformly() {
    let (e, nu, p, radius) = (1.0, 0.3, 1e-3, 5.0);
    let centre = [0.3, -0.2, 0.1];
    let mut mesh = icosphere(centre, radius, 3);
    // The bottom cap carries the exact displacement; the rest the pressure.
    for (t, f) in mesh.triangles.iter().zip(mesh.triangle_face.iter_mut()) {
        let z = t.iter().map(|&i| mesh.vertices[i as usize][2]).sum::<f64>() / 3.0;
        if z < centre[2] - 0.9 * radius {
            *f = FaceId(1);
        }
    }
    let model = model_of(mesh, 2);
    let k = -(1.0 - 2.0 * nu) * p / e;
    let support = PieceSupport::Fixed {
        u0: centre.map(|c| -k * c),
        gradient: std::array::from_fn(|i| std::array::from_fn(|j| if i == j { k } else { 0.0 })),
    };
    let problem = ElasticProblem {
        youngs: e,
        poisson: nu,
        faces: vec![FaceMechanics::Pressure(p), FaceMechanics::Support(support)],
    };
    let r = solve_elastic(&model, &problem, &options(0.9)).unwrap_or_else(|e| panic!("{e}"));
    let mut worst = 0.0f64;
    for (v, d) in model.solids[0]
        .vertices
        .iter()
        .zip(&r.vertex_displacement[0])
    {
        for i in 0..3 {
            worst = worst.max((d[i] - k * (v[i] - centre[i])).abs());
        }
    }
    assert!(worst <= 1e-8 * (k * radius).abs(), "worst {worst}");
    assert_balanced(&r, p * 4.0 * std::f64::consts::PI * radius * radius);
}

/// Two boxes as two solids; faces 0..5 on the first, 6..11 on the second.
fn two_boxes() -> crate::model::AnalysisModel {
    let mut m = model_of(box_mesh([0.0; 3], [2.0, 1.0, 1.0]), 6);
    let mut second = box_mesh([3.0, 0.0, 0.0], [5.0, 1.0, 1.0]);
    second.id = SolidId(1);
    for f in &mut second.triangle_face {
        f.0 += 6;
    }
    for f in 0..6u32 {
        let mut info: FaceInfo = m.faces[f as usize].clone();
        info.id = FaceId(6 + f);
        info.solid = SolidId(1);
        info.surface = Surface::Freeform;
        info.role_source = RoleSource::Inferred;
        m.faces.push(info);
    }
    m.solids.push(second);
    m
}

#[test]
fn rb1_free_motions_are_declined_by_name() {
    // Roller on one face: free to slide in y and z and spin about x.
    let model = model_of(box_mesh([0.0; 3], [4.0, 1.0, 1.0]), 6);
    let mut faces = vec![FaceMechanics::Free; 6];
    faces[0] = ROLLER;
    let problem = ElasticProblem {
        youngs: 1.0,
        poisson: 0.3,
        faces,
    };
    match solve_elastic(&model, &problem, &options(0.5)) {
        Err(AnalysisError::FreeMotion { solid, motions }) => {
            assert_eq!(solid, SolidId(0));
            assert_eq!(motions.len(), 3, "{motions:?}");
            let spins: Vec<_> = motions
                .iter()
                .filter(|m| matches!(m, RigidMotion::Rotation { .. }))
                .collect();
            assert_eq!(spins.len(), 1);
            if let RigidMotion::Rotation { direction, .. } = spins[0] {
                assert!(direction[0].abs() > 1.0 - 1e-6, "{direction:?}");
            }
        }
        other => panic!("expected FreeMotion, got {other:?}"),
    }
    // A second body with no support: declined, naming that body.
    let model = two_boxes();
    let mut faces = vec![FaceMechanics::Free; 12];
    faces[0] = clamp();
    let problem = ElasticProblem {
        youngs: 1.0,
        poisson: 0.3,
        faces,
    };
    match solve_elastic(&model, &problem, &options(0.5)) {
        Err(AnalysisError::FreeMotion { solid, motions }) => {
            assert_eq!(solid, SolidId(1));
            assert_eq!(motions.len(), 6);
        }
        other => panic!("expected FreeMotion, got {other:?}"),
    }
    // Controls: supporting both bodies, and no support at all.
    let mut faces = vec![FaceMechanics::Free; 12];
    faces[0] = clamp();
    faces[6] = clamp();
    faces[1] = FaceMechanics::Traction([0.0, 0.0, -1e-3]);
    let problem = ElasticProblem {
        youngs: 1.0,
        poisson: 0.3,
        faces,
    };
    assert!(solve_elastic(&model, &problem, &options(0.5)).is_ok());
    let problem = ElasticProblem {
        youngs: 1.0,
        poisson: 0.3,
        faces: vec![FaceMechanics::Free; 12],
    };
    assert_eq!(
        solve_elastic(&model, &problem, &options(0.5)).err(),
        Some(AnalysisError::NoSupport)
    );
}

#[test]
fn e6_the_nitsche_factor_barely_moves_the_answer() {
    let model = model_of(box_mesh([0.0; 3], [10.0, 1.0, 1.0]), 6);
    let mut faces = vec![FaceMechanics::Free; 6];
    faces[0] = clamp();
    faces[1] = FaceMechanics::Traction([0.0, 0.0, -1e-3]);
    let mut tips = Vec::new();
    for gamma in [5.0, 10.0, 20.0] {
        let problem = ElasticProblem {
            youngs: 1.0,
            poisson: 0.3,
            faces: faces.clone(),
        };
        let o = AnalysisOptions {
            nitsche_factor: gamma,
            ..options(0.5)
        };
        let r =
            solve_elastic(&model, &problem, &o).unwrap_or_else(|e| panic!("gamma {gamma}: {e}"));
        let tip = r
            .field
            .displacement_at([10.0, 0.5, 0.5])
            .unwrap_or_else(|e| panic!("{e}"))
            .unwrap_or([f64::NAN; 3]);
        tips.push(tip[2]);
        assert_balanced(&r, 1e-3);
    }
    let (lo, hi) = tips
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(l, h), &t| {
            (l.min(t), h.max(t))
        });
    assert!((hi - lo).abs() <= 0.02 * lo.abs(), "{tips:?}");
}

#[test]
fn the_jacobi_eigensolver_diagonalises() {
    let a: [[f64; 6]; 6] =
        std::array::from_fn(|i| std::array::from_fn(|j| 1.0 / (1.0 + i as f64 + j as f64)));
    let (values, vectors) = symmetric_eigen(a);
    for k in 0..6 {
        for i in 0..6 {
            let av: f64 = (0..6).map(|j| a[i][j] * vectors[j][k]).sum();
            assert!((av - values[k] * vectors[i][k]).abs() < 1e-14, "pair {k}");
        }
    }
    assert!(values.windows(2).all(|w| w[0] <= w[1]));
}
