//! Generative design gates GD1–GD3, GD5 and GD6
//! (docs/architecture/analysis-generative-design.md).

use super::*;
use crate::element::vector_boundary::PieceSupport;
use crate::mesh_test_helpers::{box_mesh, model_of};
use crate::model::FaceId;

const E: f64 = 2000.0;
const NU: f64 = 0.3;

fn options(h: f64) -> AnalysisOptions {
    AnalysisOptions {
        target_cell_size: Some(h),
        cg_rel_tol: 1e-12,
        coarsest_dofs: 600,
        ..AnalysisOptions::default()
    }
}

/// A cantilever L × 2 × w clamped at x = 0 and pulled down by a
/// traction on x = L.
fn cantilever(length: f64, width: f64) -> (AnalysisModel, ElasticProblem) {
    let model = model_of(box_mesh([0.0; 3], [length, 2.0, width]), 6);
    let mut faces = vec![FaceMechanics::Free; 6];
    faces[0] = FaceMechanics::Support(PieceSupport::Fixed {
        u0: [0.0; 3],
        gradient: [[0.0; 3]; 3],
    });
    faces[1] = FaceMechanics::Traction([0.0, -1.0, 0.0]);
    (
        model,
        ElasticProblem {
            youngs: E,
            poisson: NU,
            faces,
        },
    )
}

/// GD1: ∂c/∂x against central differences of c, through the filter, on
/// a varied design (every free entry its own density). Measured worst
/// 7.3e-9.
#[test]
fn sensitivities_match_central_differences() {
    let (model, problem) = cantilever(4.0, 1.0);
    let design = GenerativeOptions::new(0.5, 0.75);
    let mut ev = Evaluator::new(&model, &problem, design, &[], &options(0.5))
        .unwrap_or_else(|e| panic!("{e}"));
    let n = ev.entries().len();
    let x: Vec<f64> = (0..n)
        .map(|i| 0.25 + 0.5 * ((7 * i) % 11) as f64 / 10.0)
        .collect();
    let base = ev.evaluate(&x).unwrap_or_else(|e| panic!("{e}"));
    // Every entry is free, the clamped ones (whose Nitsche terms scale
    // with them) among them.
    let free: Vec<usize> = (0..n).filter(|&i| ev.floors()[i] == 0.0).collect();
    assert_eq!(free.len(), n);
    let delta = 1e-5;
    let mut worst = 0.0f64;
    for &j in free.iter().step_by(free.len() / 6).chain(free.first()) {
        let mut up = x.clone();
        up[j] += delta;
        let mut down = x.clone();
        down[j] -= delta;
        let cu = ev
            .evaluate(&up)
            .unwrap_or_else(|e| panic!("{e}"))
            .compliance;
        let cd = ev
            .evaluate(&down)
            .unwrap_or_else(|e| panic!("{e}"))
            .compliance;
        let fd = (cu - cd) / (2.0 * delta);
        let rel = (fd - base.gradient[j]).abs() / fd.abs();
        worst = worst.max(rel);
    }
    assert!(worst < 1e-5, "{worst}");
}

/// The design's fraction of face-neighbour pairs (both full interior
/// cells) whose physical densities differ by more than ½.
fn jumps(d: &Design, h: f64) -> f64 {
    let key = |c: [f64; 3]| c.map(|x| (x / h).floor() as i64);
    let mut cells: Vec<([i64; 3], f64)> = d
        .entries
        .iter()
        .zip(&d.densities)
        .filter(|(e, _)| (e.volume - h * h * h).abs() < 1e-9 * h * h * h)
        .map(|(e, &x)| (key(e.centre), x))
        .collect();
    cells.sort_by_key(|c| c.0);
    let find = |k: [i64; 3]| {
        cells
            .binary_search_by(|c| c.0.cmp(&k))
            .ok()
            .map(|i| cells[i].1)
    };
    let (mut pairs, mut big) = (0usize, 0usize);
    for &(k, x) in &cells {
        for a in 0..3 {
            let mut n = k;
            n[a] += 1;
            if let Some(y) = find(n) {
                pairs += 1;
                big += usize::from((x - y).abs() > 0.5);
            }
        }
    }
    big as f64 / pairs.max(1) as f64
}

/// GD2: with p = 1 compliance is convex in the densities, and for a bar
/// in tension on symmetry rollers the optimum is known: every
/// full-length fibre carries the same strain, so a design uniform along
/// the bar is optimal, c* = F²L/(E f A).
#[test]
fn a_convex_bar_in_tension_reaches_its_analytic_optimum() {
    let (l, side, h) = (12.0, 3.0, 0.5);
    let model = model_of(box_mesh([0.0; 3], [l, side, side]), 6);
    let roller = FaceMechanics::Support(PieceSupport::Roller);
    let mut faces = vec![FaceMechanics::Free; 6];
    faces[0] = roller;
    faces[2] = roller;
    faces[4] = roller;
    faces[1] = FaceMechanics::Traction([1.0, 0.0, 0.0]);
    let problem = ElasticProblem {
        youngs: E,
        poisson: NU,
        faces,
    };
    let area = side * side;
    for f in [0.5, 0.7] {
        let mut g = GenerativeOptions::new(f, 0.75);
        g.penalty = 1.0;
        let d = optimise(&model, &problem, &g, &[], &options(h)).unwrap_or_else(|e| panic!("{e}"));
        let optimum = area * area * l / (E * f * area);
        assert!(d.converged, "f {f}: {} updates", d.history.len());
        assert!(
            (d.volume_fraction - f).abs() < 1e-9,
            "{}",
            d.volume_fraction
        );
        let ratio = d.compliance / optimum;
        assert!((ratio - 1.0).abs() < 0.01, "f {f}: {ratio}");
    }
}

/// A cantilever 12 × 2 × 1 at 0.5 mm cells, volume 0.4.
fn short_cantilever(radius: f64) -> Design {
    let (model, problem) = cantilever(12.0, 1.0);
    let mut g = GenerativeOptions::new(0.4, radius);
    g.max_iterations = 60;
    // The loaded end is kept (GD5).
    optimise(&model, &problem, &g, &[FaceId(1)], &options(0.5)).unwrap_or_else(|e| panic!("{e}"))
}

/// GD3, GD5 and GD6: with the filter (r = 1.5 h) the cantilever's design
/// is smooth, keeps the material by its kept (loaded) end solid, stiffens far
/// beyond its first (uniform) design and is the same bits twice. The
/// negative control: without it (r below a cell) the design is ragged,
/// with neighbouring cells jumping, and reads an artificially low
/// compliance (top88's checkerboard; measured: 28.5% of pairs jump
/// against 4.7%, compliance 3.3 against 11.3, from 27.2 uniform).
#[test]
fn the_filter_smooths_the_cantilever_and_the_design_is_repeatable() {
    let smooth = short_cantilever(0.75);
    let ragged = short_cantilever(0.1);
    let (js, jr) = (jumps(&smooth, 0.5), jumps(&ragged, 0.5));
    eprintln!(
        "GD3 {js} {jr} c {} {} first {}",
        smooth.compliance, ragged.compliance, smooth.history[0].compliance
    );
    assert!(jr > 3.0 * js, "{jr} vs {js}");
    assert!(ragged.compliance < smooth.compliance);
    assert!(smooth.compliance < 0.6 * smooth.history[0].compliance);
    assert!((smooth.volume_fraction - 0.4).abs() < 1e-9);
    // Every entry at least its floor, and the loaded end's cells,
    // wholly within the keep thickness, solid.
    assert!(
        smooth
            .floors
            .iter()
            .zip(&smooth.densities)
            .all(|(&f, &x)| x >= f)
    );
    assert!(smooth.floors.iter().filter(|&&f| f == 1.0).count() >= 4);
    let again = short_cantilever(0.75);
    assert_eq!(again.densities, smooth.densities);
    assert_eq!(again.compliance.to_bits(), smooth.compliance.to_bits());
}

/// Out-of-range options and kept material over the target decline.
#[test]
fn bad_options_and_an_overfull_keep_decline() {
    let (model, problem) = cantilever(4.0, 1.0);
    for g in [
        GenerativeOptions::new(0.0, 0.75),
        GenerativeOptions::new(1.5, 0.75),
        GenerativeOptions::new(0.5, 0.0),
        GenerativeOptions {
            penalty: 0.5,
            ..GenerativeOptions::new(0.5, 0.75)
        },
    ] {
        assert!(matches!(
            optimise(&model, &problem, &g, &[], &options(0.5)),
            Err(AnalysisError::InvalidGenerative { .. })
        ));
    }
    // Both ends of a 4 mm cantilever kept 1.5 mm deep: 3 mm of 4.
    let g = GenerativeOptions::new(0.5, 1.5);
    assert!(matches!(
        optimise(&model, &problem, &g, &[FaceId(0), FaceId(1)], &options(1.0)),
        Err(AnalysisError::KeptVolumeExceedsTarget { .. })
    ));
    // A support at a non-zero displacement moves its load with the
    // densities.
    let mut moved = problem.clone();
    moved.faces[0] = FaceMechanics::Support(PieceSupport::Fixed {
        u0: [0.01, 0.0, 0.0],
        gradient: [[0.0; 3]; 3],
    });
    let g = GenerativeOptions::new(0.5, 0.75);
    assert!(matches!(
        optimise(&model, &moved, &g, &[], &options(0.5)),
        Err(AnalysisError::Unsupported { .. })
    ));
}

/// GD4 (offsets): a 12 × 4 × 2 cantilever (4 cells through its
/// thickness), loaded end kept, at 5 grid offsets, r in mm. A part only
/// 2 cells thick is not offset-stable: its compliance spread 17% (11%
/// without the keep), each offset moving its cut cells and SIMP finding
/// a different local optimum. Refinement is recorded in the plan: the
/// solid part's compliance moves 0.27% from 0.5 to 0.25 mm cells, the
/// design's 5.2%.
#[test]
fn the_cantilever_design_holds_its_compliance_across_grid_offsets() {
    let model = model_of(box_mesh([0.0; 3], [12.0, 4.0, 2.0]), 6);
    let mut faces = vec![FaceMechanics::Free; 6];
    faces[0] = FaceMechanics::Support(PieceSupport::Fixed {
        u0: [0.0; 3],
        gradient: [[0.0; 3]; 3],
    });
    faces[1] = FaceMechanics::Traction([0.0, -1.0, 0.0]);
    let problem = ElasticProblem {
        youngs: E,
        poisson: NU,
        faces,
    };
    let shifts = [
        [0.0; 3],
        [0.21, 0.0, 0.0],
        [0.0, 0.37, 0.11],
        [0.43, 0.19, 0.31],
        [0.61; 3],
    ];
    let mut g = GenerativeOptions::new(0.4, 0.75);
    g.max_iterations = 100;
    let mut c = Vec::new();
    for shift in shifts {
        let mut o = options(0.5);
        o.cg_rel_tol = 1e-8;
        o.grid_shift = shift.map(|f| f * 0.5);
        let d = optimise(&model, &problem, &g, &[FaceId(1)], &o).unwrap_or_else(|e| panic!("{e}"));
        c.push(d.compliance);
    }
    for x in &c {
        assert!((x / c[0] - 1.0).abs() < 0.05, "{c:?}");
    }
}

/// GD1 with the projection (β = 8): the chain rule through
/// x̄(x̃) matches central differences. Entries where the projection
/// has saturated read 0 on both sides.
#[test]
fn projected_sensitivities_match_central_differences() {
    let (model, problem) = cantilever(4.0, 1.0);
    let design = GenerativeOptions::new(0.5, 0.75);
    let mut ev = Evaluator::new(&model, &problem, design, &[], &options(0.5))
        .unwrap_or_else(|e| panic!("{e}"));
    ev.set_projection(Some(projection::Projection {
        beta: 8.0,
        kind: projection::ProjectionKind::Guest,
    }));
    let n = ev.entries().len();
    let x: Vec<f64> = (0..n)
        .map(|i| 0.02 + 0.2 * ((7 * i) % 11) as f64 / 10.0)
        .collect();
    let base = ev.evaluate(&x).unwrap_or_else(|e| panic!("{e}"));
    let delta = 1e-6;
    let mut checked = 0;
    for j in (0..n).step_by(n / 8) {
        let mut up = x.clone();
        up[j] += delta;
        let mut down = x.clone();
        down[j] -= delta;
        let cu = ev
            .evaluate(&up)
            .unwrap_or_else(|e| panic!("{e}"))
            .compliance;
        let cd = ev
            .evaluate(&down)
            .unwrap_or_else(|e| panic!("{e}"))
            .compliance;
        let fd = (cu - cd) / (2.0 * delta);
        if fd.abs() > 1e-6 * base.compliance {
            checked += 1;
            let rel = (fd - base.gradient[j]).abs() / fd.abs();
            assert!(rel < 1e-5, "entry {j}: {fd} against {}", base.gradient[j]);
        }
    }
    assert!(checked >= 4, "{checked}");
}

/// GD1 through the whole chain: density filter, threshold projection
/// (β = 4) and the print filter (built along +y), against central
/// differences of the compliance.
#[test]
fn printed_sensitivities_match_central_differences() {
    let (model, problem) = cantilever(4.0, 1.0);
    let mut design = GenerativeOptions::new(0.5, 0.75);
    design.build_direction = Some(print::BuildDirection {
        axis: 1,
        negative: false,
    });
    let mut ev = Evaluator::new(&model, &problem, design, &[], &options(0.5))
        .unwrap_or_else(|e| panic!("{e}"));
    ev.set_projection(Some(projection::Projection {
        beta: 4.0,
        kind: projection::ProjectionKind::Threshold(0.5),
    }));
    let n = ev.entries().len();
    let x: Vec<f64> = (0..n)
        .map(|i| 0.3 + 0.6 * ((7 * i) % 11) as f64 / 10.0)
        .collect();
    let base = ev.evaluate(&x).unwrap_or_else(|e| panic!("{e}"));
    let delta = 1e-6;
    let mut checked = 0;
    for j in (0..n).step_by(n / 10) {
        let mut up = x.clone();
        up[j] += delta;
        let mut down = x.clone();
        down[j] -= delta;
        let cu = ev
            .evaluate(&up)
            .unwrap_or_else(|e| panic!("{e}"))
            .compliance;
        let cd = ev
            .evaluate(&down)
            .unwrap_or_else(|e| panic!("{e}"))
            .compliance;
        let fd = (cu - cd) / (2.0 * delta);
        if fd.abs() > 1e-6 * base.compliance {
            checked += 1;
            let rel = (fd - base.gradient[j]).abs() / fd.abs();
            assert!(rel < 1e-5, "entry {j}: {fd} against {}", base.gradient[j]);
        }
    }
    assert!(checked >= 5, "{checked}");
}

/// PH-125: the 12 × 2 × 1 cantilever under r = 0.75 mm (members of
/// 1.5 mm) is warned of as thinner than the filter's member, at its
/// 1 mm; the 12 × 4 × 2 block is not (its corners are 1% of it).
#[test]
fn a_part_thinner_than_the_filter_member_is_warned_of() {
    let thin_warning = |size: [f64; 3]| {
        let (_, problem) = cantilever(size[0], size[2]);
        let model = model_of(box_mesh([0.0; 3], size), 6);
        let design = GenerativeOptions::new(0.5, 0.75);
        let ev = Evaluator::new(&model, &problem, design, &[], &options(0.5))
            .unwrap_or_else(|e| panic!("{e}"));
        ev.report().warnings.iter().find_map(|w| match w {
            Warning::ThinDesignSpace {
                fraction,
                thickness,
                member,
            } => Some((*fraction, *thickness, *member)),
            _ => None,
        })
    };
    let plate = thin_warning([12.0, 2.0, 1.0]).unwrap_or_else(|| panic!("no warning"));
    assert!(plate.0 > 0.99 && (plate.1 - 1.0).abs() < 1e-9, "{plate:?}");
    assert_eq!(plate.2, 1.5);
    assert_eq!(thin_warning([12.0, 4.0, 2.0]), None);
}

/// PH-124: a filter radius under 3 cells is warned of (r = 0.75 mm at
/// 0.5 mm cells); at 0.25 mm cells it is not.
#[test]
fn a_filter_under_three_cells_is_warned_of() {
    let (model, problem) = cantilever(4.0, 2.0);
    let warned = |h: f64| {
        let ev = Evaluator::new(
            &model,
            &problem,
            GenerativeOptions::new(0.5, 0.75),
            &[],
            &options(h),
        )
        .unwrap_or_else(|e| panic!("{e}"));
        ev.report()
            .warnings
            .iter()
            .any(|w| matches!(w, Warning::FilterUnderResolved { .. }))
    };
    assert!(warned(0.5));
    assert!(!warned(0.25));
}
