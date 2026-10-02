//! The solve on one level: an energy-only start, a G1 phase, then a
//! damped Gauss–Newton (Levenberg–Marquardt) phase on G1 + G2.
//!
//! Undamped Gauss–Newton on the G2 rows was measured unstable on real
//! fillet corners (the pentagon flagship): the first full step from the
//! energy-only surface flipped the cap normal near a rim (g1 = 1), the
//! next level converged, the level after flipped again. The G1 phase
//! gives the G2 linearisation correct tangent planes to start from, and
//! each G2 step is accepted only when the nonlinear merit (rim residuals
//! over their weights plus the energy) decreases and no rim normal
//! flips; otherwise the damping grows.

use super::chart::CharacteristicMap;
use super::energy::{Anchor, Integrand};
use super::energy_mcv::Reference;
use super::factor::{EnergyFactor, energy_factor};
use super::level::{Level, assemble};
use super::measure::cap_jet;
use super::rim_rows::{Reduction, rim_rows, transversal_rows, transversal_targets};
use super::solve::solve;
use super::{CapEnergy, VariationalCapError, VariationalCapOptions, VariationalCapTargets};

/// Smallest normal-angle sine used to weight the G1 rows. The rim rows
/// are weighted by 1/target against energy rows of scale √λ (0.1 at the
/// default λ); past a weight ratio of ~1e7 the stacked QR cannot resolve
/// the directions only the energy fixes and the solve reads singular
/// (measured with 1e-14 targets). A tighter target still governs the
/// stopping test.
const WEIGHT_FLOOR_G1: f64 = 1.0e-6;

/// Smallest normalised curvature jump used to weight the G2 rows (the
/// same bound as [`WEIGHT_FLOOR_G1`] in the G2 row scale).
const WEIGHT_FLOOR_G2: f64 = 1.0e-5;

/// Step halvings per Gauss–Newton step before the run stops.
const STEP_HALVINGS: usize = 12;

/// Solve one level; returns the coefficients and the iterations used.
pub(super) fn solve_level(
    level: &Level,
    targets: &VariationalCapTargets,
    options: &VariationalCapOptions,
    phi: &CharacteristicMap,
    start: Option<&[Vec<f64>; 3]>,
) -> Result<([Vec<f64>; 3], usize), VariationalCapError> {
    let characteristic = |_s: usize, u: f64, v: f64| phi.metric(u, v);
    let parametric = Integrand::Parametric {
        metric: &characteristic,
        eps: options.membrane_weight,
    };
    let mut energy = energy_factor(level, &parametric, None)?;
    let red = Reduction {
        hs: &level.hs,
        xp: &level.xp,
        fixed: &level.fixed,
    };
    let w1 = targets.g1.max(WEIGHT_FLOOR_G1);
    let w2 = targets.g2.map(|g| g.max(WEIGHT_FLOOR_G2));
    let lambda = options.fairness_weight;
    let mut iterations = 0;
    // A refined level starts from the previous level's surface (exactly
    // representable after knot insertion, and feasible: the spoke rows
    // and the rim rows hold for it unchanged), so the merit is monotone
    // across levels; a first level starts from the energy alone.
    let (mut y, g1_phase) = match start {
        Some(x) => (project(level, x), 0),
        None => (energy.minimiser()?, 1),
    };
    let mut coords = assemble(level, &y);
    // The tangential anchor: the level's start surface, well
    // parameterised by construction (energy-only, or the previous
    // level's accepted result).
    let prior = Prior {
        targets: transversal_targets(&level.space, &coords, &level.samples),
        weight: options.transversal_weight,
    };
    // G1 phase: one solve with the G1 rows against the characteristic
    // energy, which holds the parameterisation; only a first level has
    // it (a warm start is already G1). The step is shortened by halving
    // until the surface keeps the served regularity floor (PK-614): on
    // a low mixed-radius trihedral apex with circular strips (pyramid-3
    // h0.4 radii 0.05/0.15/0.08 at G1) the full step folded the
    // energy-only start (regularity 0.39 -> 0), and every guarded step
    // after it was refused, so the fit declined folded. The damped runs
    // below carry the G1 rows the rest of the way.
    for _ in 0..g1_phase {
        let (rows, rhs) = all_rows(level, &red, &coords, w1, None, &prior);
        let full = solve(&energy, lambda, &rows, &rhs)?;
        let mut alpha = 1.0;
        for _ in 0..STEP_HALVINGS {
            let trial_y: [Vec<f64>; 3] = std::array::from_fn(|c| {
                y[c].iter()
                    .zip(&full[c])
                    .map(|(a, b)| a + alpha * (b - a))
                    .collect()
            });
            let trial = assemble(level, &trial_y);
            if super::measure::served_regularity(&level.space, &trial).ratio
                >= super::measure::SERVED_MIN_REGULARITY
            {
                y = trial_y;
                coords = trial;
                break;
            }
            alpha *= 0.5;
        }
        iterations += 1;
    }
    // The tangential anchor: the G1 phase's surface (or the warm start),
    // G1 and transversal by construction.
    let anchor_coords = coords.clone();
    // Then damped, merit-guarded runs against the chosen energy, its
    // reference frozen per run so the merit compares like with like
    // (`reference_passes` re-freezes it; measured on spherical n-gons,
    // extra passes change the sag by < 1e-5). The G1 rows are
    // homogeneous in the cross derivative, so an unguarded G1 step can
    // shrink it towards a fold once an energy that ignores tangential
    // motion takes over (measured: pyramid4 h=3 r=0.05 regularity
    // 0.29 -> 0 in one undamped step); every step after the first is
    // therefore guarded, G1-only fits included.
    // The sixth-order mean-curvature variation needs the G2 rows for its
    // boundary conditions; a G1-only fit bends instead.
    let chosen = match (options.energy, w2) {
        (CapEnergy::MeanCurvatureVariation, None) => CapEnergy::Bending,
        (energy, _) => energy,
    };
    let run_options = VariationalCapOptions {
        energy: chosen,
        ..*options
    };
    for _ in 0..options.reference_passes.max(1) {
        energy = refresh_energy(level, &run_options, &coords, &anchor_coords, energy)?;
        let (c, yy, it) = damped_run(
            level,
            &red,
            &energy,
            coords,
            y,
            lambda,
            (w1, w2),
            &prior,
            options,
        )?;
        coords = c;
        y = yy;
        iterations += it;
    }
    Ok((coords, iterations))
}

/// Coefficients, null-space coordinates and iterations of one run.
type RunResult = ([Vec<f64>; 3], [Vec<f64>; 3], usize);

/// One damped Gauss–Newton run against a frozen energy: the full step
/// of the linearised problem, shortened by halving until the nonlinear
/// merit decreases with no rim normal flipped and no fold.
///
/// Levenberg–Marquardt damping was measured wrong here: scaled by the
/// rim rows' column norms (~1e8 per column in null-space coordinates,
/// where every column mixes everything) even its weakest setting
/// dwarfed the energy term (λG ~ 1e-2), so every step dropped the
/// energy's direction and the fit's interior never left its G1-phase
/// shape (the fairness weights had no measurable effect). Shortening
/// the undamped step keeps its direction, energy included.
#[allow(clippy::too_many_arguments)]
fn damped_run(
    level: &Level,
    red: &Reduction<'_>,
    energy: &EnergyFactor,
    mut coords: [Vec<f64>; 3],
    mut y: [Vec<f64>; 3],
    lambda: f64,
    (w1, w2): (f64, Option<f64>),
    prior: &Prior,
    options: &VariationalCapOptions,
) -> Result<RunResult, VariationalCapError> {
    let mut iterations = 0;
    for _ in 0..options.max_gauss_newton {
        let current = merit(level, &coords, &y, energy, lambda, w1, w2, prior);
        let (rows, rhs) = all_rows(level, red, &coords, w1, w2, prior);
        let full = solve(energy, lambda, &rows, &rhs)?;
        let mut alpha = 1.0;
        let mut accepted = None;
        for _ in 0..STEP_HALVINGS {
            let trial_y: [Vec<f64>; 3] = std::array::from_fn(|c| {
                y[c].iter()
                    .zip(&full[c])
                    .map(|(a, b)| a + alpha * (b - a))
                    .collect()
            });
            let trial = assemble(level, &trial_y);
            let m = merit(level, &trial, &trial_y, energy, lambda, w1, w2, prior);
            if m < current {
                accepted = Some((trial, trial_y, m, alpha));
                break;
            }
            alpha *= 0.5;
        }
        let Some((trial, trial_y, m, alpha)) = accepted else {
            break;
        };
        coords = trial;
        y = trial_y;
        iterations += 1;
        if alpha == 1.0 && (current - m) <= 1e-4 * current {
            break;
        }
    }
    Ok((coords, y, iterations))
}

/// The bending energy against the current iterate (data-dependent), or
/// the unchanged energy (characteristic).
fn refresh_energy(
    level: &Level,
    options: &VariationalCapOptions,
    coords: &[Vec<f64>; 3],
    anchor_coords: &[Vec<f64>; 3],
    previous: EnergyFactor,
) -> Result<EnergyFactor, VariationalCapError> {
    let space = &level.space;
    let anchor_reference = |s: usize, u: f64, v: f64| {
        let j = cap_jet(space, anchor_coords, s, u, v);
        let e1 = j.su.normalized()?;
        let n = j.su.cross(&j.sv).normalized()?;
        Some((j.p, [e1, n.cross(&e1)]))
    };
    let anchor = Anchor {
        reference: &anchor_reference,
        weight: options.anchor_weight,
    };
    match options.energy {
        CapEnergy::CharacteristicThinPlate => Ok(previous),
        CapEnergy::Bending => {
            let reference = |s: usize, u: f64, v: f64| {
                let j = cap_jet(space, coords, s, u, v);
                let normal = j.su.cross(&j.sv).normalized()?;
                Some((normal, [j.su.dot(&j.su), j.su.dot(&j.sv), j.sv.dot(&j.sv)]))
            };
            energy_factor(
                level,
                &Integrand::Bending {
                    reference: &reference,
                    regular: options.parameter_regularisation,
                    eps: options.membrane_weight,
                },
                Some(&anchor),
            )
        }
        CapEnergy::MeanCurvatureVariation => {
            let reference = |s: usize, u: f64, v: f64| {
                let j = cap_jet(space, coords, s, u, v);
                Reference::from_partials(j.su, j.sv, j.suu, j.suv, j.svv)
            };
            energy_factor(
                level,
                &Integrand::MeanCurvatureVariation {
                    reference: &reference,
                    regular: options.parameter_regularisation,
                    eps: options.membrane_weight,
                },
                Some(&anchor),
            )
        }
    }
}

/// Nonlinear merit: weighted rim G1/G2 residuals plus the energy, or
/// infinity when any rim normal flips against its neighbour.
#[allow(clippy::too_many_arguments)]
fn merit(
    level: &Level,
    coords: &[Vec<f64>; 3],
    y: &[Vec<f64>; 3],
    energy: &EnergyFactor,
    lambda: f64,
    w1: f64,
    w2: Option<f64>,
    prior: &Prior,
) -> f64 {
    let mut total = 0.0;
    for (sample, wt) in level.samples.iter().zip(&prior.targets) {
        let ((u, v), cross_u) = super::measure::rim_uv(sample);
        let j = super::measure::cap_jet_owned(
            &level.space,
            coords,
            sample.sector,
            u,
            v,
            super::measure::rim_owner(sample),
        );
        let Some(n_cap) = j.su.cross(&j.sv).normalized() else {
            return f64::INFINITY;
        };
        let n = sample.jet.normal;
        if n_cap.dot(&n) <= 0.0 {
            return f64::INFINITY;
        }
        let (w, w2v) = if cross_u {
            (j.su, j.suu)
        } else {
            (j.sv, j.svv)
        };
        let len2 = w.dot(&w).max(1e-24);
        let r1 = n.dot(&w) / len2.sqrt() / w1;
        let rp = (w - *wt).length() / wt.length().max(1e-12) * prior.weight;
        total += r1 * r1 + rp * rp;
        if let Some(w2) = w2 {
            let r2 = (n.dot(&w2v) - sample.jet.ii(w, w)) / len2 / w2;
            let tangent = if cross_u { j.sv } else { j.su };
            let r3 = (n.dot(&j.suv) - sample.jet.ii(tangent, w))
                / (len2.sqrt() * tangent.length().max(1e-12))
                / w2;
            total += r2 * r2 + r3 * r3;
        }
    }
    // A step that folds the interior is never an improvement, however
    // much it lowers the rim residuals.
    // A step that flips a normal anywhere on the served grid is never
    // accepted (PK-606: see `measure::FOLD_FLOOR`).
    if super::measure::served_regularity(&level.space, coords).ratio < super::measure::FOLD_FLOOR {
        return f64::INFINITY;
    }
    let e = energy.value(y);
    total + lambda * e
}

/// Null-space coordinates of full coefficients `x` (least squares:
/// `Z` has orthonormal columns).
fn project(level: &Level, x: &[Vec<f64>; 3]) -> [Vec<f64>; 3] {
    let hs = &level.hs;
    std::array::from_fn(|c| {
        let d: Vec<f64> = hs
            .free
            .iter()
            .enumerate()
            .map(|(r, &i)| x[c][i] - level.xp[c][r])
            .collect();
        hs.project(&d)
    })
}

/// The transversality prior (see [`transversal_rows`]).
struct Prior {
    targets: Vec<crate::geometry::Vec3>,
    weight: f64,
}

/// The rim rows stacked over the prior rows.
fn all_rows(
    level: &Level,
    red: &Reduction<'_>,
    coords: &[Vec<f64>; 3],
    w1: f64,
    w2: Option<f64>,
    prior: &Prior,
) -> (faer::Mat<f64>, Vec<f64>) {
    let (a, mut ra) = rim_rows(&level.space, red, coords, &level.samples, w1, w2);
    let (b, rb) = transversal_rows(
        &level.space,
        red,
        &level.samples,
        &prior.targets,
        prior.weight,
    );
    let rows = faer::Mat::from_fn(a.nrows() + b.nrows(), a.ncols(), |i, j| {
        if i < a.nrows() {
            a.read(i, j)
        } else {
            b.read(i - a.nrows(), j)
        }
    });
    ra.extend(rb);
    (rows, ra)
}
