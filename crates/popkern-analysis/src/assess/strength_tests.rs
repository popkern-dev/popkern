//! The strength checks with sourced shear, compressive and gyroid
//! strengths, by hand: inter-layer shear (PH-080), Coulomb–Mohr
//! (PH-084), the paraboloidal yield (PH-109) and the infill core
//! (PH-079). Every load is a held linear
//! displacement u = G·x on all faces, so the field is exact (patch
//! test) and every reading is round-off from the hand value.

use super::safety::{Check, SafetyAssessment, Unchecked};
use super::*;
use crate::element::vector_boundary::PieceSupport;
use crate::material::{Material, Property, preset};
use crate::mesh_test_helpers::{box_mesh, model_of};
use crate::options::AnalysisOptions;
use crate::structural::assemble::FaceMechanics;
use crate::structural::infill::{Infill, core_yield};
use crate::structural::{ElasticProblem, solve_with};

fn material(name: &str) -> Material {
    *preset(name).unwrap_or_else(|e| panic!("{e}"))
}

fn props(m: &Material) -> (f64, f64) {
    (
        m.youngs_modulus.value().unwrap_or_else(|| panic!("E")),
        m.poisson_ratio.value().unwrap_or_else(|| panic!("ν")),
    )
}

/// The safety assessment of a 10 mm cube with every face held at
/// u = G·x, optionally infilled, along `direction`.
fn held(
    m: &Material,
    gradient: [[f64; 3]; 3],
    infill: Option<&Infill>,
    direction: Option<[f64; 3]>,
) -> (SafetyAssessment, f64) {
    let model = model_of(box_mesh([0.0; 3], [10.0, 10.0, 10.0]), 6);
    let (e, nu) = props(m);
    let support = FaceMechanics::Support(PieceSupport::Fixed {
        u0: [0.0; 3],
        gradient,
    });
    let problem = ElasticProblem {
        youngs: e,
        poisson: nu,
        faces: vec![support; 6],
    };
    let options = AnalysisOptions {
        target_cell_size: Some(1.0),
        cg_rel_tol: 1e-12,
        coarsest_dofs: 600,
        ..AnalysisOptions::default()
    };
    let r = solve_with(&model, &problem, None, infill, &options).unwrap_or_else(|e| panic!("{e}"));
    let o = AssessOptions {
        build_direction: direction,
        ..AssessOptions::default()
    };
    let a = assess(&model, &r, m, &o).unwrap_or_else(|e| panic!("{e}"));
    (a.safety, r.report.infill_core_volume)
}

fn factor(s: &SafetyAssessment, check: Check) -> Option<f64> {
    s.checks
        .iter()
        .find(|c| c.check == check)
        .and_then(|c| c.factor)
}

fn close(a: Option<f64>, b: f64, tol: f64) -> bool {
    a.is_some_and(|a| (a - b).abs() <= tol * b.abs())
}

/// PH-080: pure shear τ in the xy plane. Layers stacked along x (or y)
/// carry τ on their planes: S_layer/τ. Stacked along z, nothing shears
/// the layers. At the worst-case orientation the largest shear on any
/// plane is τ. PETG has no S_layer and lists the mode as unchecked.
#[test]
fn inter_layer_shear_reads_s_layer_over_the_layer_shear() {
    let tau = 5.0;
    for name in ["PLA (FDM)", "ABS (FDM)", "PC (FDM)", "PETG (FDM)"] {
        let m = material(name);
        let (e, nu) = props(&m);
        let gamma = tau * 2.0 * (1.0 + nu) / e;
        let g = [[0.0, 0.5 * gamma, 0.0], [0.5 * gamma, 0.0, 0.0], [0.0; 3]];
        let shear = m.strength.interlayer_shear.value();
        for direction in [None, Some([1.0, 0.0, 0.0]), Some([0.0, 0.0, 1.0])] {
            let (s, _) = held(&m, g, None, direction);
            match shear {
                Some(sl) => {
                    assert!(s.unchecked.is_empty(), "{name}");
                    let f = factor(&s, Check::InterLayerShear);
                    if direction == Some([0.0, 0.0, 1.0]) {
                        assert!(f.is_none_or(|f| f > 1e6 * sl / tau), "{name}: {f:?}");
                    } else {
                        assert!(close(f, sl / tau, 1e-6), "{name} {direction:?}: {f:?}");
                    }
                }
                None => {
                    assert_eq!(s.unchecked, [Unchecked::InterLayerShear], "{name}");
                    assert!(factor(&s, Check::InterLayerShear).is_none());
                }
            }
        }
    }
}

/// PH-084: uniaxial compression σ = −p along x (held lateral expansion
/// νp/E). Coulomb–Mohr reads σ_c/p: half σ_c passes at 2, 1.25 σ_c
/// fails at 0.8. The negative control: the same PLA without a
/// compressive strength falls back to the tension-only check, which
/// sees only round-off tension: compression goes unchecked.
#[test]
fn coulomb_mohr_reads_the_compressive_strength_in_compression() {
    for name in ["PLA (FDM)", "PETG (FDM)"] {
        let m = material(name);
        let (e, nu) = props(&m);
        let sc = m
            .strength
            .ultimate_compressive
            .value()
            .unwrap_or_else(|| panic!("{name}: σ_c"));
        for (p, expected) in [(0.5 * sc, 2.0), (1.25 * sc, 0.8)] {
            let eps = p / e;
            let g = [[-eps, 0.0, 0.0], [0.0, nu * eps, 0.0], [0.0, 0.0, nu * eps]];
            let (s, _) = held(&m, g, None, None);
            let f = factor(&s, Check::InPlaneCoulombMohr);
            assert!(close(f, expected, 1e-6), "{name} at {p}: {f:?}");
            assert!(factor(&s, Check::InPlaneBrittle).is_none());
            // In-layer, with layers stacked along z: the same.
            let (s, _) = held(&m, g, None, Some([0.0, 0.0, 1.0]));
            assert!(close(factor(&s, Check::InPlaneCoulombMohr), expected, 1e-6));
        }
    }
    let mut bare = material("PLA (FDM)");
    bare.strength.ultimate_compressive = Property::Unsourced { note: "control" };
    let (e, nu) = props(&bare);
    let eps = 30.0 / e;
    let g = [[-eps, 0.0, 0.0], [0.0, nu * eps, 0.0], [0.0, 0.0, nu * eps]];
    let (s, _) = held(&bare, g, None, None);
    assert!(factor(&s, Check::InPlaneCoulombMohr).is_none());
    // The tension it sees is round-off.
    let c = s.checks.iter().find(|c| c.check == Check::InPlaneBrittle);
    assert!(c.is_some_and(|c| c.stress <= 1e-6 * 30.0), "{:?}", s.checks);
}

/// PH-109: ABS yields in-plane and has a compressive yield, so it is
/// checked on Raghava's paraboloid: uniaxial compression p reads σ_c/p,
/// uniaxial tension σ reads σ_t/σ, and pure shear τ (I₁ = 0) reads
/// √(σ_c σ_t)/(√3 τ). The negative control: a metal without σ_c keeps
/// the von Mises check, σ_y/p in compression.
#[test]
fn a_yielding_material_with_a_compressive_yield_reads_the_paraboloid() {
    let m = material("ABS (FDM)");
    let (e, nu) = props(&m);
    let st = &m.strength;
    let t = st.inplane.value().unwrap_or(0.0);
    let c = st.ultimate_compressive.value().unwrap_or(0.0);
    let p = 20.0;
    let eps = p / e;
    let compression = [[-eps, 0.0, 0.0], [0.0, nu * eps, 0.0], [0.0, 0.0, nu * eps]];
    let tension = compression.map(|r| r.map(|x| -x));
    let (s, _) = held(&m, compression, None, None);
    let f = factor(&s, Check::InPlaneParaboloidalYield);
    assert!(close(f, c / p, 1e-6), "{f:?}");
    assert!(factor(&s, Check::InPlaneYield).is_none());
    let (s, _) = held(&m, tension, None, None);
    assert!(close(
        factor(&s, Check::InPlaneParaboloidalYield),
        t / p,
        1e-6
    ));
    let gamma = p * 2.0 * (1.0 + nu) / e;
    let shear = [[0.0, 0.5 * gamma, 0.0], [0.5 * gamma, 0.0, 0.0], [0.0; 3]];
    let (s, _) = held(&m, shear, None, None);
    let expected = (c * t).sqrt() / (3.0f64.sqrt() * p);
    assert!(close(
        factor(&s, Check::InPlaneParaboloidalYield),
        expected,
        1e-6
    ));
    let metal = material("Aluminium 6061-T6");
    let (e, nu) = props(&metal);
    let eps = p / e;
    let g = [[-eps, 0.0, 0.0], [0.0, nu * eps, 0.0], [0.0, 0.0, nu * eps]];
    let (s, _) = held(&metal, g, None, None);
    let y = metal.strength.inplane.value().unwrap_or(0.0);
    assert!(close(factor(&s, Check::InPlaneYield), y / p, 1e-6));
    assert!(factor(&s, Check::InPlaneParaboloidalYield).is_none());
}

/// PH-079, the check's arithmetic: at D = 1 the core is its parent
/// (Bean's fits give 1 in every mode), so the cube held at uniaxial
/// strain is uniform and exact, and the core check reads S_m/σ: 1 at
/// σ = S_m, 1/2 at twice that.
#[test]
fn a_full_density_core_at_the_parent_yield_reads_one() {
    let m = material("PLA (FDM)");
    let (e, nu) = props(&m);
    let infill = Infill {
        density: 1.0,
        walls_mm: 2.0,
        line_width_mm: 0.4,
    };
    let parent = m.strength.inplane.value().unwrap_or(0.0);
    assert!((core_yield(1.0, parent) - parent).abs() <= 1e-12 * parent);
    for (scale, expected) in [(1.0, 1.0), (2.0, 0.5)] {
        let eps = scale * parent / e;
        let g = [
            [eps, 0.0, 0.0],
            [0.0, -nu * eps, 0.0],
            [0.0, 0.0, -nu * eps],
        ];
        let (s, core) = held(&m, g, Some(&infill), None);
        assert!(core > 0.0);
        assert!(s.unchecked.is_empty(), "{:?}", s.unchecked);
        let f = factor(&s, Check::InfillCore);
        assert!(close(f, expected, 1e-6), "{f:?}");
    }
}

/// Bean 2025 eqn 8 by hand: at D = 0.25 the compressive ratio is
/// 0.25 − 0.75^0.67·0.25^1.2 = 0.09375 and the shear ratio
/// 0.25 − 0.75^1.3·0.25^1.3 = 0.13652; at D = 0.5, 0.22643 and 0.33506.
/// Compression governs the core's allowable, and both ratios rise from
/// 0 to 1.
#[test]
fn the_gyroid_yield_fit_matches_bean_by_hand() {
    use crate::structural::infill::gyroid_yield;
    for (d, compressive, shear) in [(0.25, 0.09375, 0.13652), (0.5, 0.22643, 0.33506)] {
        let (c, s) = gyroid_yield(d);
        assert!((c - compressive).abs() < 1e-5, "{d}: {c}");
        assert!((s - shear).abs() < 1e-5, "{d}: {s}");
        assert_eq!(core_yield(d, 10.0), 10.0 * c);
    }
    assert_eq!(gyroid_yield(0.0), (0.0, 0.0));
    assert_eq!(gyroid_yield(1.0), (1.0, 1.0));
    let mut last = 0.0;
    for i in 1..=20 {
        let (c, _) = gyroid_yield(f64::from(i) / 20.0);
        assert!(c > last, "{i}");
        last = c;
    }
}

/// An infilled bar pulled along its length (D = 0.25): the core check
/// is Bean's allowable over the recovered core peak, with the parent's
/// scatter scaled and the gyroid tests' 3.6% combined. The peak sits
/// where the end wall feeds load into the core (measured 2.5× the
/// core's parallel share at 1 mm cells, 2.9× at 0.5 mm), so a pulled
/// bar does not read the core uniformly.
#[test]
fn an_infilled_bar_checks_its_core_against_the_gyroid_yield() {
    use crate::structural::infill::{GYROID_YIELD_SCATTER, gyroid_yield};
    let m = material("PLA (FDM)");
    let (e, nu) = props(&m);
    let infill = Infill {
        density: 0.25,
        walls_mm: 2.0,
        line_width_mm: 0.4,
    };
    let model = model_of(box_mesh([0.0; 3], [60.0, 10.0, 10.0]), 6);
    let roller = FaceMechanics::Support(PieceSupport::Roller);
    let mut faces = vec![FaceMechanics::Free; 6];
    faces[0] = roller;
    faces[2] = roller;
    faces[4] = roller;
    faces[1] = FaceMechanics::Traction([10.0, 0.0, 0.0]);
    let problem = ElasticProblem {
        youngs: e,
        poisson: nu,
        faces,
    };
    let options = AnalysisOptions {
        target_cell_size: Some(1.0),
        cg_rel_tol: 1e-12,
        coarsest_dofs: 600,
        ..AnalysisOptions::default()
    };
    let r = solve_with(&model, &problem, None, Some(&infill), &options)
        .unwrap_or_else(|e| panic!("{e}"));
    let a = assess(&model, &r, &m, &AssessOptions::default()).unwrap_or_else(|e| panic!("{e}"));
    let peak = r.stresses[0]
        .core_peak
        .unwrap_or_else(|| panic!("no core peak"));
    let c = a
        .safety
        .checks
        .iter()
        .find(|c| c.check == Check::InfillCore)
        .unwrap_or_else(|| panic!("no core check"));
    let parent = m.strength.inplane.value().unwrap_or(0.0);
    let parent_sd = m.strength.inplane_scatter.value().unwrap_or(0.0);
    let mean = core_yield(0.25, parent);
    let sd = (gyroid_yield(0.25).0 * parent_sd).hypot(GYROID_YIELD_SCATTER * mean);
    assert_eq!(c.allowable, mean);
    assert!((c.allowable_low - (mean - 2.0 * sd)).abs() <= 1e-12 * mean);
    assert_eq!(c.stress, peak.equivalent);
    assert!(close(c.factor, mean / peak.equivalent, 1e-12));
    assert_eq!(c.at.1, peak.position);
}

/// The negative controls: a core thinner than half a cell is never the
/// majority of a cell, so no sample reads it and the core is listed as
/// unchecked; without infill there is neither a core check nor an
/// unchecked core.
#[test]
fn an_unread_core_is_unchecked_and_a_solid_part_has_no_core_check() {
    let m = material("PLA (FDM)");
    let g = [[1e-3, 0.0, 0.0], [0.0; 3], [0.0; 3]];
    // Walls 4.8 mm leave a 0.4 mm core in 1 mm cells.
    let thin = Infill {
        density: 0.25,
        walls_mm: 4.8,
        line_width_mm: 0.4,
    };
    let (s, core) = held(&m, g, Some(&thin), None);
    assert!(core > 0.0, "no core at all");
    assert_eq!(s.unchecked, [Unchecked::InfillCore]);
    assert!(factor(&s, Check::InfillCore).is_none());
    let (s, core) = held(&m, g, None, None);
    assert_eq!(core, 0.0);
    assert!(s.unchecked.is_empty());
    assert!(s.checks.iter().all(|c| c.check != Check::InfillCore));
}

/// The core peak of a 60 × 10 × 10 PLA bar (D = 0.25, 2 mm walls)
/// pulled at 10 MPa, at cell size `h`, for a gyroid of line width `line`.
fn pulled_core_peak(h: f64, line: f64) -> f64 {
    let m = material("PLA (FDM)");
    let (e, nu) = props(&m);
    let infill = Infill {
        density: 0.25,
        walls_mm: 2.0,
        line_width_mm: line,
    };
    let model = model_of(box_mesh([0.0; 3], [60.0, 10.0, 10.0]), 6);
    let roller = FaceMechanics::Support(PieceSupport::Roller);
    let mut faces = vec![FaceMechanics::Free; 6];
    faces[0] = roller;
    faces[2] = roller;
    faces[4] = roller;
    faces[1] = FaceMechanics::Traction([10.0, 0.0, 0.0]);
    let problem = ElasticProblem {
        youngs: e,
        poisson: nu,
        faces,
    };
    let options = AnalysisOptions {
        target_cell_size: Some(h),
        cg_rel_tol: 1e-10,
        ..AnalysisOptions::default()
    };
    let r = solve_with(&model, &problem, None, Some(&infill), &options)
        .unwrap_or_else(|e| panic!("{e}"));
    r.stresses[0].core_peak.map_or(f64::NAN, |p| p.equivalent)
}

/// PH-110: averaged over the gyroid's cell (λ = 3.09·0.4/0.25 = 4.9
/// mm), the pulled bar's core peak moves 8.9% from 1 mm to 0.5 mm cells
/// (measured 2026-10-13: 1.86 and 2.02 times the core's parallel share,
/// 2.02 again at 0.25 mm). The negative control: a gyroid so fine its
/// cell is below a grid cell is read cell by cell, and its peak, at the
/// corner where the end wall meets the core, grows by more than 10%.
#[test]
fn the_core_peak_averaged_over_the_gyroid_cell_settles_under_refinement() {
    let (coarse, fine) = (pulled_core_peak(1.0, 0.4), pulled_core_peak(0.5, 0.4));
    assert!((fine / coarse - 1.0).abs() < 0.1, "{coarse} -> {fine}");
    let (coarse, fine) = (pulled_core_peak(1.0, 0.01), pulled_core_peak(0.5, 0.01));
    assert!(fine / coarse > 1.1, "{coarse} -> {fine}");
}

/// PH-113: in the pulled infilled bar the walls' tension peaks inside
/// the part, at the end wall's inner corner against the core (measured
/// σ₁ 27.4 against 23.1 on the surface at 1 mm cells), and the
/// inter-layer tension check reads it there. The negative control: the
/// same bar solid has no interior samples, and the check reads the
/// surface.
#[test]
fn an_interior_wall_peak_governs_the_wall_checks_of_an_infilled_bar() {
    let m = material("PLA (FDM)");
    let (e, nu) = props(&m);
    let model = model_of(box_mesh([0.0; 3], [60.0, 10.0, 10.0]), 6);
    let roller = FaceMechanics::Support(PieceSupport::Roller);
    let mut faces = vec![FaceMechanics::Free; 6];
    faces[0] = roller;
    faces[2] = roller;
    faces[4] = roller;
    faces[1] = FaceMechanics::Traction([10.0, 0.0, 0.0]);
    let problem = ElasticProblem {
        youngs: e,
        poisson: nu,
        faces,
    };
    let options = AnalysisOptions {
        target_cell_size: Some(1.0),
        cg_rel_tol: 1e-10,
        ..AnalysisOptions::default()
    };
    let infill = Infill {
        density: 0.25,
        walls_mm: 2.0,
        line_width_mm: 0.4,
    };
    for filled in [true, false] {
        let r = solve_with(&model, &problem, None, filled.then_some(&infill), &options)
            .unwrap_or_else(|e| panic!("{e}"));
        let a = assess(&model, &r, &m, &AssessOptions::default()).unwrap_or_else(|e| panic!("{e}"));
        let c = a
            .safety
            .checks
            .iter()
            .find(|c| c.check == Check::InterLayer)
            .unwrap_or_else(|| panic!("no layer check"));
        let solid = &r.stresses[0];
        let surface = solid
            .surface
            .iter()
            .map(|x| x.stress.max_principal)
            .fold(f64::NEG_INFINITY, f64::max);
        let interior = solid
            .interior
            .iter()
            .map(|x| x.stress.max_principal)
            .fold(f64::NEG_INFINITY, f64::max);
        if filled {
            assert!(interior > 1.1 * surface, "{interior} vs {surface}");
            assert_eq!(c.stress, interior);
            assert!(solid.interior.iter().any(|x| x.position == c.at.1));
        } else {
            assert!(solid.interior.is_empty());
            assert_eq!(c.stress, surface);
        }
    }
}
