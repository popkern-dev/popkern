//! Point-like supports (PH-033): how firmly a support on a small patch
//! holds the rigid motions.

use super::assemble::FaceMechanics;
use super::*;
use crate::element::vector_boundary::PieceSupport;
use crate::mesh_test_helpers::{box_mesh, model_of};
use crate::model::{FaceId, SolidMesh};

/// The 10 × 1 × 1 box with its x = 0 face split into a central square
/// patch of `fraction` of the face area (face 6) and the ring around it
/// (face 0).
fn patched_cantilever(fraction: f64) -> crate::model::AnalysisModel {
    let mut mesh: SolidMesh = box_mesh([0.0; 3], [10.0, 1.0, 1.0]);
    let keep: Vec<usize> = (0..mesh.triangles.len())
        .filter(|&t| mesh.triangle_face[t] != FaceId(0))
        .collect();
    mesh.triangles = keep.iter().map(|&t| mesh.triangles[t]).collect();
    mesh.triangle_face = keep.iter().map(|&t| mesh.triangle_face[t]).collect();
    let a = 0.5 * fraction.sqrt();
    let outer: [u32; 4] = [0, 4, 6, 2];
    let mut inner = [0u32; 4];
    for (k, &o) in outer.iter().enumerate() {
        let v = mesh.vertices[o as usize];
        let p = [
            0.0,
            0.5 + 2.0 * a * (v[1] - 0.5),
            0.5 + 2.0 * a * (v[2] - 0.5),
        ];
        inner[k] = mesh.vertices.len() as u32;
        mesh.vertices.push(p);
    }
    if fraction >= 1.0 {
        mesh.triangles.push([outer[0], outer[1], outer[2]]);
        mesh.triangles.push([outer[0], outer[2], outer[3]]);
        mesh.triangle_face.extend([FaceId(6); 2]);
        return model_of(mesh, 7);
    }
    mesh.triangles.push([inner[0], inner[1], inner[2]]);
    mesh.triangles.push([inner[0], inner[2], inner[3]]);
    mesh.triangle_face.extend([FaceId(6); 2]);
    for k in 0..4 {
        let (o0, o1, i0, i1) = (outer[k], outer[(k + 1) % 4], inner[k], inner[(k + 1) % 4]);
        mesh.triangles.push([o0, o1, i1]);
        mesh.triangles.push([o0, i1, i0]);
        mesh.triangle_face.extend([FaceId(0); 2]);
    }
    model_of(mesh, 7)
}

fn clamp_patch(fraction: f64, cell: f64) -> Result<StructuralResult, AnalysisError> {
    let model = patched_cantilever(fraction);
    let mut faces = vec![FaceMechanics::Free; 7];
    faces[6] = FaceMechanics::Support(PieceSupport::Fixed {
        u0: [0.0; 3],
        gradient: [[0.0; 3]; 3],
    });
    faces[1] = FaceMechanics::Traction([0.0, 0.0, -1e-3]);
    let problem = ElasticProblem {
        youngs: 1.0,
        poisson: 0.3,
        faces,
    };
    let options = AnalysisOptions {
        target_cell_size: Some(cell),
        coarsest_dofs: 600,
        cg_max_iters: 2000,
        ..AnalysisOptions::default()
    };
    solve_elastic(&model, &problem, &options)
}

fn spans(r: &StructuralResult) -> Vec<f64> {
    r.report
        .warnings
        .iter()
        .filter_map(|w| match w {
            crate::report::Warning::PointLikeSupport { span, .. } => Some(*span),
            _ => None,
        })
        .collect()
}

#[test]
fn rb1_point_like_patch_supports_warn_and_full_clamps_do_not() {
    // Measured (h = 0.25 / 0.125 / 0.0625, cantilever 10 × 1 × 1,
    // E = 1, tip under 1e-3 end traction; beam theory 4.0):
    //   patch area   λ_min/scale (h=.25, .125)   tip
    //   1            2.6e-5, 1.4e-5              4.007, 4.012
    //   0.3          2.7e-6, 1.4e-6, 7.1e-7      5.224, 5.357, 5.335
    //   0.1          3.2e-7, 1.6e-7, 8.0e-8      7.958, 12.31, 11.73
    //   0.01         3.2e-9, 1.6e-9              231, 138
    //   0.001        below 1e-10: FreeMotion     declined
    // The patch side in cells decides the error (1.3 cells: −32%;
    // 2.2: −2%; 2.5: +5%).
    for (fraction, cell, warns) in [
        (1.0, 0.25, false),
        (0.3, 0.125, false),
        (0.1, 0.25, true),
        (0.01, 0.125, true),
    ] {
        let r = clamp_patch(fraction, cell).unwrap_or_else(|e| panic!("{fraction}: {e}"));
        assert_eq!(
            !spans(&r).is_empty(),
            warns,
            "{fraction} at {cell}: {:?}",
            spans(&r)
        );
        assert_eq!(r.report.rigid_mode_ratios.len(), 1);
    }
    // Exact span: a 0.3 patch is √0.3 = 0.548 across.
    let r = clamp_patch(0.3, 0.5).unwrap_or_else(|e| panic!("{e}"));
    assert!(
        (spans(&r)[0] - 0.3f64.sqrt()).abs() < 1e-9,
        "{:?}",
        spans(&r)
    );
    // Below the rigid-mode tolerance the solve declines, naming the
    // two bending rotations through the patch.
    match clamp_patch(0.001, 0.25) {
        Err(AnalysisError::FreeMotion { motions, .. }) => {
            assert_eq!(motions.len(), 2, "{motions:?}");
            assert!(
                motions
                    .iter()
                    .all(|m| matches!(m, rigid::RigidMotion::Rotation { .. }))
            );
        }
        other => panic!("expected FreeMotion, got {other:?}"),
    }
}

#[test]
fn a_single_facet_support_on_a_sphere_names_its_rotations_and_three_rollers_do_not() {
    use crate::mesh_test_helpers::icosphere;
    // One facet of an icosphere clamped: the facet is smaller than two
    // cells, so all three rotations about its centroid are named (or
    // the solve declines outright); never silently solved.
    let mut mesh = icosphere([0.0; 3], 5.0, 2);
    let nt = mesh.triangles.len();
    mesh.triangle_face = (0..nt).map(|t| FaceId(u32::from(t == 0))).collect();
    let model = model_of(mesh, 2);
    let mut faces = vec![FaceMechanics::Free; 2];
    faces[1] = FaceMechanics::Support(PieceSupport::Fixed {
        u0: [0.0; 3],
        gradient: [[0.0; 3]; 3],
    });
    faces[0] = FaceMechanics::Pressure(1e-3);
    let problem = ElasticProblem {
        youngs: 1.0,
        poisson: 0.3,
        faces,
    };
    let options = AnalysisOptions {
        target_cell_size: Some(1.0),
        coarsest_dofs: 600,
        ..AnalysisOptions::default()
    };
    match solve_elastic(&model, &problem, &options) {
        Ok(r) => {
            let named: Vec<_> = r
                .report
                .warnings
                .iter()
                .filter_map(|w| match w {
                    crate::report::Warning::PointLikeSupport { rotations, .. } => {
                        Some(rotations.len())
                    }
                    _ => None,
                })
                .collect();
            assert!(
                !named.is_empty() && named[0] >= 2,
                "{:?}",
                r.report.warnings
            );
        }
        Err(AnalysisError::FreeMotion { motions, .. }) => assert!(!motions.is_empty()),
        Err(e) => panic!("{e}"),
    }
    // Control: rollers on three orthogonal faces of a box hold it
    // without a warning.
    let model = model_of(box_mesh([0.0; 3], [4.0, 2.0, 2.0]), 6);
    let roller = FaceMechanics::Support(PieceSupport::Roller);
    let mut faces = vec![FaceMechanics::Free; 6];
    faces[0] = roller;
    faces[2] = roller;
    faces[4] = roller;
    faces[1] = FaceMechanics::Traction([1e-3, 0.0, 0.0]);
    let problem = ElasticProblem {
        youngs: 1.0,
        poisson: 0.3,
        faces,
    };
    let options = AnalysisOptions {
        target_cell_size: Some(0.5),
        coarsest_dofs: 600,
        ..AnalysisOptions::default()
    };
    let r = solve_elastic(&model, &problem, &options).unwrap_or_else(|e| panic!("{e}"));
    assert!(spans(&r).is_empty(), "{:?}", r.report.warnings);
}
