//! Gates promoted from probes E1 and E4 (map, exploration items): EAS
//! cantilever accuracy at 2 cells through the thickness, full and
//! rotated, against a converged 3D solution of the same problem.

use super::assemble::FaceMechanics;
use super::*;
use crate::element::vector_boundary::PieceSupport;
use crate::mesh_test_helpers::{box_mesh, model_of, rotated};

/// Rotation used for the cut-cell beam: 17° about (1, 2, 3)/√14
/// through the root centre (the same Rodrigues map as `rotated`).
pub(crate) fn rotate_vector(v: [f64; 3], degrees: f64) -> [f64; 3] {
    use popkern_math::PortableMath;
    let n = 14.0f64.sqrt();
    let k = [1.0 / n, 2.0 / n, 3.0 / n];
    let theta = degrees * std::f64::consts::PI / 180.0;
    let (s, c) = (theta.psin(), theta.pcos());
    let kxv = [
        k[1] * v[2] - k[2] * v[1],
        k[2] * v[0] - k[0] * v[2],
        k[0] * v[1] - k[1] * v[0],
    ];
    let kv = k[0] * v[0] + k[1] * v[1] + k[2] * v[2];
    std::array::from_fn(|a| v[a] * c + kxv[a] * s + k[a] * kv * (1.0 - c))
}

/// Tip deflection along the load of a 10 × 1 × 1 cantilever, clamped at
/// x = 0 (Nitsche), end shear 1e-3 N/mm² in −z, E = 1, with
/// `cells` cells through the thickness; `degrees` rotates the whole
/// problem about the root centre.
pub(crate) fn tip(cells: f64, nu: f64, degrees: f64) -> f64 {
    let centre = [0.0, 0.5, 0.5];
    let mesh = rotated(box_mesh([0.0; 3], [10.0, 1.0, 1.0]), degrees, centre);
    let model = model_of(mesh, 6);
    let down = rotate_vector([0.0, 0.0, -1.0], degrees);
    let mut faces = vec![FaceMechanics::Free; 6];
    faces[0] = FaceMechanics::Support(PieceSupport::Fixed {
        u0: [0.0; 3],
        gradient: [[0.0; 3]; 3],
    });
    faces[1] = FaceMechanics::Traction(down.map(|d| 1e-3 * d));
    let problem = ElasticProblem {
        youngs: 1.0,
        poisson: nu,
        faces,
    };
    let options = AnalysisOptions {
        target_cell_size: Some(1.0 / cells),
        cg_rel_tol: 1e-10,
        cg_max_iters: 400,
        coarsest_dofs: 600,
        ..AnalysisOptions::default()
    };
    let r =
        solve_elastic(&model, &problem, &options).unwrap_or_else(|e| panic!("{cells} cells: {e}"));
    let tip_point = rotate_vector([10.0, 0.0, 0.0], degrees);
    let p = [0, 1, 2].map(|a| centre[a] + tip_point[a]);
    let d = r
        .field
        .displacement_at(p)
        .unwrap_or_else(|e| panic!("{e}"))
        .unwrap_or_else(|| panic!("tip outside the part"));
    d[0] * down[0] + d[1] * down[1] + d[2] * down[2]
}

/// Mean tip deflection along the load of a cantilevered plate
/// `len × width × thick`, rotated by `degrees`, on cells of size `h`.
fn plate_tip(len: f64, width: f64, thick: f64, h: f64, degrees: f64) -> f64 {
    let centre = [0.0, 0.5 * width, 0.5 * thick];
    let model = model_of(
        rotated(box_mesh([0.0; 3], [len, width, thick]), degrees, centre),
        6,
    );
    let down = rotate_vector([0.0, 0.0, -1.0], degrees);
    let mut faces = vec![FaceMechanics::Free; 6];
    faces[0] = FaceMechanics::Support(PieceSupport::Fixed {
        u0: [0.0; 3],
        gradient: [[0.0; 3]; 3],
    });
    faces[1] = FaceMechanics::Traction(down.map(|d| 1e-4 * d));
    let problem = ElasticProblem {
        youngs: 1.0,
        poisson: 0.38,
        faces,
    };
    let options = AnalysisOptions {
        target_cell_size: Some(h),
        cg_rel_tol: 1e-10,
        cg_max_iters: 400,
        coarsest_dofs: 600,
        ..AnalysisOptions::default()
    };
    let r = solve_elastic(&model, &problem, &options).unwrap_or_else(|e| panic!("h {h}: {e}"));
    let (mut sum, mut n) = (0.0, 0.0);
    for (v, d) in model.solids[0]
        .vertices
        .iter()
        .zip(&r.vertex_displacement[0])
    {
        let local = rotate_vector(
            [v[0] - centre[0], v[1] - centre[1], v[2] - centre[2]],
            -degrees,
        );
        if local[0] > len - 1e-6 {
            sum += d[0] * down[0] + d[1] * down[1] + d[2] * down[2];
            n += 1.0;
        }
    }
    sum / n
}

/// E1, promoted: for the preset range of ν, two cells through the
/// thickness are within 5% (axis-aligned) and 8% (17° rotated, all
/// surface cells cut) of the converged solution.
#[test]
fn e1_two_cells_through_the_thickness_are_within_the_gate() {
    for nu in [0.3, 0.38] {
        let reference = tip(8.0, nu, 0.0);
        let full = tip(2.0, nu, 0.0);
        assert!(
            (full - reference).abs() <= 0.05 * reference,
            "nu {nu} full: {full} vs {reference}"
        );
        let reference = tip(8.0, nu, 17.0);
        let turned = tip(2.0, nu, 17.0);
        assert!(
            (turned - reference).abs() <= 0.08 * reference,
            "nu {nu} rotated: {turned} vs {reference}"
        );
    }
}

/// E4, promoted: deflection within the map's 10% for a 1.2 mm wall at
/// h = 0.5 (2.4 cells) and for the enclosure-class 2 mm wall at h = 1.4
/// (1.43 cells), axis-aligned and rotated. Measured −0.6% to −5.4% for
/// the 2 mm wall depending on orientation and aspect (map, E4).
#[test]
fn e4_thin_walls_at_the_default_resolution() {
    for degrees in [0.0, 17.0] {
        let reference = plate_tip(10.0, 5.0, 1.2, 0.2, degrees);
        let coarse = plate_tip(10.0, 5.0, 1.2, 0.5, degrees);
        assert!(
            (coarse - reference).abs() <= 0.10 * reference,
            "1.2 mm at {degrees}: {coarse} vs {reference}"
        );
        let reference = plate_tip(10.0, 5.0, 2.0, 1.0 / 3.0, degrees);
        let coarse = plate_tip(10.0, 5.0, 2.0, 1.4, degrees);
        assert!(
            (coarse - reference).abs() <= 0.10 * reference,
            "2 mm at {degrees}: {coarse} vs {reference}"
        );
    }
}

#[test]
fn nearly_incompressible_materials_are_declined() {
    let model = model_of(box_mesh([0.0; 3], [4.0, 1.0, 1.0]), 6);
    let mut faces = vec![FaceMechanics::Free; 6];
    faces[0] = FaceMechanics::Support(PieceSupport::Roller);
    let problem = ElasticProblem {
        youngs: 1.0,
        poisson: 0.495,
        faces,
    };
    assert_eq!(
        solve_elastic(&model, &problem, &AnalysisOptions::default()).err(),
        Some(AnalysisError::NearlyIncompressible {
            poisson: 0.495,
            limit: POISSON_LIMIT
        })
    );
}
