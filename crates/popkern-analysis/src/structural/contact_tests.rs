//! Sliding contacts (PK-1073): two cubes stacked in z, pressed on top,
//! against the exact field: uniform compression in each cube and a
//! normal jump of σ/β at the contact.

use super::assemble::FaceMechanics;
use super::contact_coupling::SLIDING_PENALTY;
use super::*;
use crate::contact::ContactPatch;
use crate::element::vector_boundary::PieceSupport;
use crate::mesh_test_helpers::box_mesh;
use crate::mesh_test_models::model_of_many;
use crate::model::SolidId;

const E: f64 = 1.0;
const NU: f64 = 0.3;
const P: f64 = 1e-3;

/// Cubes [0, 10]³ and [0, 10]² × [10, 20]; faces x-min, x-max, y-min,
/// y-max, z-min, z-max for each (6..12 the upper cube's).
fn stack() -> AnalysisModel {
    model_of_many(
        vec![
            box_mesh([0.0; 3], [10.0; 3]),
            box_mesh([0.0, 0.0, 10.0], [10.0, 10.0, 20.0]),
        ],
        6,
    )
}

/// The z = 10 square, normal +z (out of the lower cube).
fn interface() -> ContactPatch {
    let p = [
        [0.0, 0.0, 10.0],
        [10.0, 0.0, 10.0],
        [10.0, 10.0, 10.0],
        [0.0, 10.0, 10.0],
    ];
    ContactPatch {
        solids: (SolidId(0), SolidId(1)),
        triangles: vec![[p[0], p[1], p[2]], [p[0], p[2], p[3]]],
    }
}

/// Rollers on z = 0 and on both cubes' x = 0 and y = 0 sides (or the
/// lower cube's alone), pressure on the top.
fn problem(upper_sides_held: bool) -> ElasticProblem {
    let mut faces = vec![FaceMechanics::Free; 12];
    for f in [0, 2, 4] {
        faces[f] = FaceMechanics::Support(PieceSupport::Roller);
    }
    if upper_sides_held {
        for f in [6, 8] {
            faces[f] = FaceMechanics::Support(PieceSupport::Roller);
        }
    }
    faces[11] = FaceMechanics::Traction([0.0, 0.0, -P]);
    ElasticProblem {
        youngs: E,
        poisson: NU,
        faces,
    }
}

fn options(h: f64) -> AnalysisOptions {
    AnalysisOptions {
        target_cell_size: Some(h),
        cg_rel_tol: 1e-12,
        cg_max_iters: 2000,
        coarsest_dofs: 600,
        ..AnalysisOptions::default()
    }
}

#[test]
fn a_pressed_stack_compresses_uniformly_with_the_contact_jump() {
    let model = stack();
    for h in [1.3, 1.0] {
        let r = solve_with_contacts(
            &model,
            &problem(true),
            None,
            None,
            &[interface()],
            &options(h),
        )
        .unwrap_or_else(|e| panic!("h {h}: {e}"));
        let lambda = E * NU / ((1.0 + NU) * (1.0 - 2.0 * NU));
        let mu = E / (2.0 * (1.0 + NU));
        let beta = SLIDING_PENALTY * (lambda + 2.0 * mu) / h;
        let exact = |s: usize, x: [f64; 3]| {
            let jump = if s == 1 { -P / beta } else { 0.0 };
            [NU * P / E * x[0], NU * P / E * x[1], -P / E * x[2] + jump]
        };
        let mut worst = 0.0f64;
        for (s, us) in r.vertex_displacement.iter().enumerate() {
            for (v, u) in model.solids[s].vertices.iter().zip(us) {
                let e = exact(s, *v);
                for i in 0..3 {
                    worst = worst.max((u[i] - e[i]).abs());
                }
            }
        }
        let scale = P / E * 20.0;
        assert!(worst <= 1e-6 * scale, "h {h}: worst {worst} of {scale}");
    }
}

/// Control: without the contact the upper cube is free along z.
#[test]
fn without_the_contact_the_upper_cube_is_free() {
    let e = solve_with_contacts(&stack(), &problem(true), None, None, &[], &options(1.3)).err();
    assert!(
        matches!(
            e,
            Some(AnalysisError::FreeMotion {
                solid: SolidId(1),
                ..
            })
        ),
        "{e:?}"
    );
}

/// A sliding contact holds only the normal: an upper cube held by its
/// contact alone slides, and the pair is free.
#[test]
fn a_cube_held_only_by_a_sliding_contact_slides() {
    let e = solve_with_contacts(
        &stack(),
        &problem(false),
        None,
        None,
        &[interface()],
        &options(1.3),
    )
    .err();
    assert!(matches!(e, Some(AnalysisError::FreeMotion { .. })), "{e:?}");
}

/// Two side-by-side cubes tied along three non-parallel normals (the
/// shared x = 10 face and two patches straddling it in y = 5 and
/// z = 5): each holds the other in every direction, but nothing holds
/// the pair, which declines as free (a third cube, apart, carries the
/// model's one support).
#[test]
fn a_pair_no_support_holds_is_free_whatever_its_contacts() {
    let model = model_of_many(
        vec![
            box_mesh([0.0; 3], [10.0; 3]),
            box_mesh([10.0, 0.0, 0.0], [20.0, 10.0, 10.0]),
            box_mesh([30.0, 0.0, 0.0], [40.0, 10.0, 10.0]),
        ],
        6,
    );
    let square = |a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]| [[a, b, c], [a, c, d]];
    let mut triangles = Vec::new();
    triangles.extend(square(
        [10.0, 0.0, 0.0],
        [10.0, 10.0, 0.0],
        [10.0, 10.0, 10.0],
        [10.0, 0.0, 10.0],
    ));
    let x_face = ContactPatch {
        solids: (SolidId(0), SolidId(1)),
        triangles,
    };
    let in_y = ContactPatch {
        solids: (SolidId(0), SolidId(1)),
        triangles: square(
            [9.5, 5.0, 2.0],
            [9.5, 5.0, 8.0],
            [10.5, 5.0, 8.0],
            [10.5, 5.0, 2.0],
        )
        .to_vec(),
    };
    let in_z = ContactPatch {
        solids: (SolidId(0), SolidId(1)),
        triangles: square(
            [9.5, 2.0, 5.0],
            [10.5, 2.0, 5.0],
            [10.5, 8.0, 5.0],
            [9.5, 8.0, 5.0],
        )
        .to_vec(),
    };
    let mut faces = vec![FaceMechanics::Free; 18];
    faces[12] = FaceMechanics::Support(PieceSupport::Fixed {
        u0: [0.0; 3],
        gradient: [[0.0; 3]; 3],
    });
    faces[7] = FaceMechanics::Traction([0.0, 0.0, -P]);
    let problem = ElasticProblem {
        youngs: E,
        poisson: NU,
        faces,
    };
    let e = solve_with_contacts(
        &model,
        &problem,
        None,
        None,
        &[x_face, in_y, in_z],
        &options(1.3),
    )
    .err();
    assert!(
        matches!(e, Some(AnalysisError::FreeMotion { ref motions, .. }) if motions.is_empty()),
        "{e:?}"
    );
}
