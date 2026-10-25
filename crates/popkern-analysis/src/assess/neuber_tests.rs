use super::*;
use crate::element::vector_boundary::PieceSupport;
use crate::material::preset;
use crate::mesh_test_helpers::{box_mesh, model_of, quarter_plate_with_hole};
use crate::options::AnalysisOptions;
use crate::structural::assemble::FaceMechanics;
use crate::structural::{ElasticProblem, solve_elastic};

const ROLLER: FaceMechanics = FaceMechanics::Support(PieceSupport::Roller);

fn options() -> AnalysisOptions {
    AnalysisOptions {
        target_cell_size: Some(0.5),
        cg_rel_tol: 1e-10,
        ..AnalysisOptions::default()
    }
}

/// The quarter plate with a hole (r 2, W 6 × 18) under remote tension
/// `sigma`, MPa, in material `name`.
fn plate(name: &str, sigma: f64) -> Option<LocalYielding> {
    let m = preset(name).unwrap_or_else(|e| panic!("{e}"));
    let e = m.youngs_modulus.value().unwrap_or(f64::NAN);
    let (mesh, nfaces, segments) = quarter_plate_with_hole(18.0, 6.0, 2.0, 2.0, 48, 18);
    let model = model_of(mesh, nfaces);
    let mut faces = vec![FaceMechanics::Free; nfaces];
    for f in [0, 2, 3] {
        faces[f] = ROLLER;
    }
    for (f, _, n) in &segments {
        if n[0] > 0.5 {
            faces[*f as usize] = FaceMechanics::Traction([sigma, 0.0, 0.0]);
        }
    }
    let problem = ElasticProblem {
        youngs: e,
        poisson: m.poisson_ratio.value().unwrap_or(0.33),
        faces,
    };
    let r = solve_elastic(&model, &problem, &options()).unwrap_or_else(|e| panic!("{e}"));
    local_yielding(&r, m, e)
}

/// A hole's peak past yield (K_t ≈ 3.5 on the remote stress): local,
/// with Neuber's elastic–perfectly-plastic strain σ_e²/(Eσ_y).
#[test]
fn a_ductile_peak_past_yield_gets_neubers_strain() {
    let m = preset("Aluminium 6061-T6").unwrap_or_else(|e| panic!("{e}"));
    let (sy, e) = (
        m.strength.inplane.value().unwrap_or(f64::NAN),
        m.youngs_modulus.value().unwrap_or(f64::NAN),
    );
    let y = plate("Aluminium 6061-T6", 0.4 * sy).unwrap_or_else(|| panic!("no local yielding"));
    assert!(y.local, "{y:?}");
    let want = y.elastic_peak_mpa * y.elastic_peak_mpa / (e * sy);
    assert_eq!(y.local_stress_mpa, Some(sy));
    assert!(
        (y.local_strain.unwrap_or(f64::NAN) - want).abs() <= 1e-15 * want,
        "{y:?}"
    );
    assert!((y.plastic_strain.unwrap_or(f64::NAN) - (want - sy / e)).abs() <= 1e-15 * want);
    // Below yield: nothing.
    assert_eq!(plate("Aluminium 6061-T6", 0.1 * sy), None);
}

/// A brittle material's peak is fracture, not local yielding.
#[test]
fn a_brittle_material_never_gets_neuber() {
    let m = preset("PETG (FDM)").unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(
        m.strength.inplane_behaviour,
        crate::material::Behaviour::Brittle
    );
    let s = m.strength.inplane.value().unwrap_or(f64::NAN);
    assert_eq!(plate("PETG (FDM)", s), None);
}

/// A bar stressed past yield everywhere is gross yielding: no Neuber
/// value, the linear solution does not describe it.
#[test]
fn yield_everywhere_is_not_local() {
    let m = preset("Aluminium 6061-T6").unwrap_or_else(|e| panic!("{e}"));
    let (sy, e) = (
        m.strength.inplane.value().unwrap_or(f64::NAN),
        m.youngs_modulus.value().unwrap_or(f64::NAN),
    );
    let model = model_of(box_mesh([0.0; 3], [10.0, 2.0, 2.0]), 6);
    let mut faces = vec![FaceMechanics::Free; 6];
    for f in [0, 2, 4] {
        faces[f] = ROLLER;
    }
    faces[1] = FaceMechanics::Traction([1.2 * sy, 0.0, 0.0]);
    let problem = ElasticProblem {
        youngs: e,
        poisson: 0.33,
        faces,
    };
    let r = solve_elastic(&model, &problem, &options()).unwrap_or_else(|e| panic!("{e}"));
    let y = local_yielding(&r, m, e).unwrap_or_else(|| panic!("no yielding"));
    assert!(!y.local && y.yielded_fraction > 0.9, "{y:?}");
    assert_eq!((y.local_strain, y.local_stress_mpa), (None, None));
}
