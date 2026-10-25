//! Assessment tests: the quantile and edge helpers, the fillet table,
//! and SF1 (hand-computed bar and beam per material preset).

use super::fillet::{concave_at, expected_under_read};
use super::safety::{BandGap, Check, SafetyAssessment, Unchecked};
use super::*;
use crate::element::vector_boundary::PieceSupport;
use crate::material::{Material, PRESETS};
use crate::mesh_test_helpers::{box_mesh, model_of};
use crate::model::Surface;
use crate::options::AnalysisOptions;
use crate::structural::assemble::FaceMechanics;
use crate::structural::{ElasticProblem, solve_elastic};

#[test]
fn the_weighted_quantile_and_polyline_distance() {
    let v = [5.0, 1.0, 3.0, 2.0, 4.0];
    let w = [1.0, 1.0, 1.0, 1.0, 6.0];
    // Cumulative weight 1, 2, 3, 4, 10 in value order 1, 2, 3, 4, 5.
    assert_eq!(weighted_quantile(&v, &w, 0.35), 4.0);
    assert_eq!(weighted_quantile(&v, &w, 0.25), 3.0);
    assert_eq!(weighted_quantile(&v, &w, 0.995), 5.0);
    let line = [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]];
    assert!((distance_to_polyline([5.0, 3.0, 4.0], &line) - 5.0).abs() < 1e-12);
    assert!((distance_to_polyline([-3.0, 4.0, 0.0], &line) - 5.0).abs() < 1e-12);
}

#[test]
fn the_under_read_table_interpolates_and_concavity_follows_the_axis() {
    assert_eq!(expected_under_read(0.5), 0.464);
    assert!((expected_under_read(2.5) - 0.1905).abs() < 1e-12);
    assert_eq!(expected_under_read(4.0), 0.067);
    let cyl = Surface::Cylinder {
        origin: [0.0; 3],
        axis: [0.0, 0.0, 1.0],
        radius: 2.0,
    };
    // Outward normal towards the axis: material outside, a concave blend.
    assert_eq!(
        concave_at(&cyl, [2.0, 0.0, 1.0], [-1.0, 0.0, 0.0]),
        Some(true)
    );
    assert_eq!(
        concave_at(&cyl, [2.0, 0.0, 1.0], [1.0, 0.0, 0.0]),
        Some(false)
    );
    let torus = Surface::Torus {
        centre: [0.0; 3],
        axis: [0.0, 0.0, 1.0],
        major_radius: 10.0,
        minor_radius: 1.0,
    };
    assert_eq!(
        concave_at(&torus, [11.0, 0.0, 0.0], [-1.0, 0.0, 0.0]),
        Some(true)
    );
    assert_eq!(
        concave_at(&Surface::Freeform, [0.0; 3], [1.0, 0.0, 0.0]),
        None
    );
}

fn props(m: &Material) -> (f64, f64) {
    (
        m.youngs_modulus.value().unwrap_or_else(|| panic!("E")),
        m.poisson_ratio.value().unwrap_or_else(|| panic!("ν")),
    )
}

/// Hand rule: (governing check, factor, band low, band high) for a
/// uniaxial peak stress `s` (σ_vM = σ₁ = s, the largest shear s/2), a
/// compressive peak of magnitude `compression` elsewhere, and local
/// error `eta`.
fn hand(m: &Material, s: f64, compression: f64, eta: f64) -> (Check, f64, f64, f64) {
    let st = &m.strength;
    let inplane = st.inplane.value().unwrap_or(0.0);
    let k = AssessOptions::default().scatter_sigmas;
    let low = |mean: f64, sd: &crate::material::Property| mean - k * sd.value().unwrap_or(0.0);
    let simple = |check: Check, mean: f64, sd: &crate::material::Property, s: f64| {
        (
            check,
            mean / s,
            low(mean, sd) / (s * (1.0 + eta)),
            mean / (s * (1.0 - eta)),
        )
    };
    let mut best = match (st.inplane_behaviour, st.ultimate_compressive.value()) {
        (crate::material::Behaviour::Yielding, None) => {
            simple(Check::InPlaneYield, inplane, &st.inplane_scatter, s)
        }
        (behaviour, Some(c)) => {
            // Coulomb–Mohr and the paraboloid agree in uniaxial stress:
            // the tension peak, or the compression peak scaled by
            // σ_t/σ_c, whichever is larger.
            let t_low = low(inplane, &st.inplane_scatter);
            let c_low = low(c, &st.ultimate_compressive_scatter);
            let eq = s.max(inplane / c * compression);
            let eq_low = s.max(t_low / c_low * compression);
            (
                match behaviour {
                    crate::material::Behaviour::Yielding => Check::InPlaneParaboloidalYield,
                    crate::material::Behaviour::Brittle => Check::InPlaneCoulombMohr,
                },
                inplane / eq,
                t_low / (eq_low * (1.0 + eta)),
                inplane / (eq * (1.0 - eta)),
            )
        }
        (crate::material::Behaviour::Brittle, None) => {
            simple(Check::InPlaneBrittle, inplane, &st.inplane_scatter, s)
        }
    };
    let mut consider = |candidate: (Check, f64, f64, f64)| {
        if candidate.1 < best.1 {
            best = candidate;
        }
    };
    if let Some(layer) = st.interlayer.value() {
        consider(simple(Check::InterLayer, layer, &st.interlayer_scatter, s));
    }
    if let Some(shear) = st.interlayer_shear.value() {
        let tau = 0.5 * s.max(compression);
        consider(simple(
            Check::InterLayerShear,
            shear,
            &st.interlayer_shear_scatter,
            tau,
        ));
    }
    best
}

#[test]
fn sf1_bar_in_tension_matches_the_hand_rule_for_every_preset() {
    // C1: rollers on three faces, 25 MPa traction; σ is uniform and
    // exact, so σ_vM = σ₁ = 25 and the local error is round-off.
    let model = model_of(box_mesh([0.0; 3], [10.0, 2.0, 2.0]), 6);
    let roller = FaceMechanics::Support(PieceSupport::Roller);
    for m in &PRESETS {
        let (e, nu) = props(m);
        let mut faces = vec![FaceMechanics::Free; 6];
        faces[0] = roller;
        faces[2] = roller;
        faces[4] = roller;
        faces[1] = FaceMechanics::Traction([25.0, 0.0, 0.0]);
        let problem = ElasticProblem {
            youngs: e,
            poisson: nu,
            faces,
        };
        let options = AnalysisOptions {
            target_cell_size: Some(0.5),
            cg_rel_tol: 1e-12,
            coarsest_dofs: 600,
            ..AnalysisOptions::default()
        };
        let r = solve_elastic(&model, &problem, &options).unwrap_or_else(|e| panic!("{e}"));
        let a = assess(&model, &r, m, &AssessOptions::default()).unwrap_or_else(|e| panic!("{e}"));
        let g = a
            .safety
            .governing
            .unwrap_or_else(|| panic!("{}: none", m.name));
        let c = &a.safety.checks[g];
        let eta = c.error.unwrap_or(0.0);
        assert!(eta < 1e-6, "{}: η {eta}", m.name);
        let (check, factor, low, high) = hand(m, 25.0, 0.0, eta);
        assert_eq!(c.check, check, "{}", m.name);
        let f = c.factor.unwrap_or_else(|| panic!("factor"));
        assert!(
            (f - factor).abs() <= 1e-6 * factor,
            "{}: {f} vs {factor}",
            m.name
        );
        let (lo, hi) = c.band.unwrap_or_else(|| panic!("band"));
        let hi = hi.unwrap_or_else(|| panic!("high"));
        assert!((lo - low).abs() <= 1e-6 * low && (hi - high).abs() <= 1e-6 * high);
        assert!(a.safety.estimated);
        // Scatter gaps are reported exactly where the datasheet has none.
        let unsourced = c.gaps.contains(&BandGap::ScatterUnsourced);
        let scatter = match c.check {
            Check::InterLayer => m.strength.interlayer_scatter.value(),
            Check::InterLayerShear => m.strength.interlayer_shear_scatter.value(),
            Check::InPlaneCoulombMohr | Check::InPlaneParaboloidalYield => m
                .strength
                .inplane_scatter
                .value()
                .and(m.strength.ultimate_compressive_scatter.value()),
            _ => m.strength.inplane_scatter.value(),
        };
        assert_eq!(unsourced, scatter.is_none(), "{}", m.name);
        // Uniform stress: no singularity, and the only fillet check is
        // empty.
        assert!(a.peaks[0].singularity.is_none(), "{}", m.name);
        assert!(a.fillets.is_empty());
    }
}

#[test]
fn sf1_cantilever_matches_the_hand_rule_for_every_preset() {
    // 10 × 1 × 1 cantilever, clamped at x = 0, end shear P = 1/6 N:
    // root σ = 6PL/(bh²) = 10 MPa (beam theory). At 4 cells through the
    // depth the recovered σ₁ peak (at the clamped root) reads 10.30
    // (measured): the factor is within 10% of the hand value, the hand
    // factor lies inside the band, and the governing check is the
    // hand rule's.
    let model = model_of(box_mesh([0.0; 3], [10.0, 1.0, 1.0]), 6);
    for m in &PRESETS {
        let (e, nu) = props(m);
        let mut faces = vec![FaceMechanics::Free; 6];
        faces[0] = FaceMechanics::Support(PieceSupport::Fixed {
            u0: [0.0; 3],
            gradient: [[0.0; 3]; 3],
        });
        faces[1] = FaceMechanics::Traction([0.0, 0.0, -1.0 / 6.0]);
        let problem = ElasticProblem {
            youngs: e,
            poisson: nu,
            faces,
        };
        let options = AnalysisOptions {
            target_cell_size: Some(0.25),
            coarsest_dofs: 600,
            ..AnalysisOptions::default()
        };
        let r = solve_elastic(&model, &problem, &options).unwrap_or_else(|e| panic!("{e}"));
        let a = assess(&model, &r, m, &AssessOptions::default()).unwrap_or_else(|e| panic!("{e}"));
        let g = a
            .safety
            .governing
            .unwrap_or_else(|| panic!("{}: none", m.name));
        let c = &a.safety.checks[g];
        let (check, factor, _, _) = hand(m, 10.0, 10.0, 0.0);
        assert_eq!(c.check, check, "{}", m.name);
        let f = c.factor.unwrap_or_else(|| panic!("factor"));
        assert!(
            (f - factor).abs() <= 0.1 * factor,
            "{}: {f} vs {factor}",
            m.name
        );
        let (lo, hi) = c.band.unwrap_or_else(|| panic!("band"));
        let hi = hi.unwrap_or(f64::INFINITY);
        if c.error.is_some() {
            assert!(
                lo <= factor && factor <= hi,
                "{}: {factor} not in [{lo}, {hi}]",
                m.name
            );
        } else {
            // In-cell fallback at the clamped corner: no error bar is
            // invented, and the band says what it leaves out.
            assert!(c.gaps.contains(&BandGap::NoErrorEstimate), "{}", m.name);
        }
        // A box has no concave edge: no singularity label.
        assert!(a.peaks[0].singularity.is_none());
    }
}

#[test]
fn a_build_direction_picks_the_layer_or_in_plane_check_by_hand() {
    // PH-040: a 25 MPa bar along x. Layers stacked along x take the
    // whole stress across them: the inter-layer check governs at
    // σ_layer/25 and the in-layer block is empty. Layers stacked along
    // z hold the bar in-plane: the in-plane check governs at
    // σ_inplane/25 and nothing crosses the layers.
    let model = model_of(box_mesh([0.0; 3], [10.0, 2.0, 2.0]), 6);
    let roller = FaceMechanics::Support(PieceSupport::Roller);
    let mut layered = 0;
    for m in PRESETS
        .iter()
        .filter(|m| m.strength.interlayer.value().is_some())
    {
        layered += 1;
        let (e, nu) = props(m);
        let mut faces = vec![FaceMechanics::Free; 6];
        faces[0] = roller;
        faces[2] = roller;
        faces[4] = roller;
        faces[1] = FaceMechanics::Traction([25.0, 0.0, 0.0]);
        let problem = ElasticProblem {
            youngs: e,
            poisson: nu,
            faces,
        };
        let options = AnalysisOptions {
            target_cell_size: Some(0.5),
            cg_rel_tol: 1e-12,
            coarsest_dofs: 600,
            ..AnalysisOptions::default()
        };
        let r = solve_elastic(&model, &problem, &options).unwrap_or_else(|e| panic!("{e}"));
        let with = |d: Option<[f64; 3]>| {
            let o = AssessOptions {
                build_direction: d,
                ..AssessOptions::default()
            };
            assess(&model, &r, m, &o)
                .unwrap_or_else(|e| panic!("{e}"))
                .safety
        };
        let factor = |s: &SafetyAssessment, check: Check| {
            s.checks
                .iter()
                .find(|c| c.check == check)
                .and_then(|c| c.factor)
        };
        let inplane = m.strength.inplane.value().unwrap_or(0.0);
        let layer = m.strength.interlayer.value().unwrap_or(0.0);
        let in_plane_check = match (
            m.strength.inplane_behaviour,
            m.strength.ultimate_compressive.value(),
        ) {
            (crate::material::Behaviour::Yielding, Some(_)) => Check::InPlaneParaboloidalYield,
            (crate::material::Behaviour::Yielding, None) => Check::InPlaneYield,
            (crate::material::Behaviour::Brittle, Some(_)) => Check::InPlaneCoulombMohr,
            (crate::material::Behaviour::Brittle, None) => Check::InPlaneBrittle,
        };
        // Only a material without S_layer leaves layer shear unchecked.
        let unsourced_shear: Vec<Unchecked> = m
            .strength
            .interlayer_shear
            .value()
            .map_or(vec![Unchecked::InterLayerShear], |_| Vec::new());
        let close = |a: Option<f64>, b: f64| a.is_some_and(|a| (a - b).abs() <= 1e-6 * b);
        // A stress that is round-off only (the empty block).
        let empty = |s: &SafetyAssessment, check: Check| {
            s.checks
                .iter()
                .find(|c| c.check == check)
                .is_some_and(|c| c.stress <= 1e-6 * 25.0)
        };
        // Across the layers (direction given with any length).
        let s = with(Some([2.0, 0.0, 0.0]));
        assert_eq!(s.build_direction, Some([1.0, 0.0, 0.0]), "{}", m.name);
        assert!(
            close(factor(&s, Check::InterLayer), layer / 25.0),
            "{}",
            m.name
        );
        assert!(empty(&s, in_plane_check), "{}", m.name);
        assert_eq!(s.checks[s.governing.unwrap_or(9)].check, Check::InterLayer);
        assert_eq!(s.unchecked, unsourced_shear, "{}", m.name);
        // Along the layers.
        let s = with(Some([0.0, 0.0, 1.0]));
        assert!(
            close(factor(&s, in_plane_check), inplane / 25.0),
            "{}",
            m.name
        );
        assert!(empty(&s, Check::InterLayer), "{}", m.name);
        assert_eq!(s.checks[s.governing.unwrap_or(9)].check, in_plane_check);
        // No direction: the worst-case orientation.
        let worst = with(None);
        assert_eq!(worst.build_direction, None);
        assert_eq!(worst.unchecked, unsourced_shear, "{}", m.name);
        assert_eq!(worst, with(AssessOptions::default().build_direction));
        assert!(
            close(factor(&worst, Check::InterLayer), layer / 25.0),
            "{}",
            m.name
        );
    }
    assert!(layered >= 4, "{layered} layered presets");
}
