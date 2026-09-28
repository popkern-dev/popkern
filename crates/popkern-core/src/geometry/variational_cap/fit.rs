//! The driver: validation, normalisation, levels, Gauss–Newton,
//! refinement, output.

use super::chart::CharacteristicMap;
use super::energy::Integrand;
use super::level::{Frame, Level, assemble, build_level, energy_factor, solve};
use super::measure::{cap_jet, regularity, rim_report, spoke_report};
use super::rim_rows::{Reduction, rim_rows};
use super::space::SectorSpace;
use super::{
    CapEnergy, RimResidual, VariationalCap, VariationalCapError, VariationalCapOptions,
    VariationalCapReport, VariationalCapSide, VariationalCapTargets, excess,
};
use crate::geometry::curve::CurveEval;
use crate::geometry::{Point3, Vec3};

/// Smallest regularity (min/mean of `|S_u × S_v|`) a served sector may
/// have (E6: revised from the measured fixtures).
const MIN_REGULARITY: f64 = 1.0e-3;

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

/// Fit a cap over the opening bounded by `sides` (see the module docs
/// for the layout and the construction). `tolerance` is the model's
/// linear tolerance, used for the corner-closure check.
///
/// # Errors
///
/// [`VariationalCapError`]: every decline names its side, stage or
/// worst sample and the measured quantity.
pub fn fit_variational_cap(
    sides: &[VariationalCapSide],
    targets: &VariationalCapTargets,
    options: &VariationalCapOptions,
    tolerance: f64,
) -> Result<VariationalCap, VariationalCapError> {
    let n = sides.len();
    if !(3..=12).contains(&n) {
        return Err(VariationalCapError::UnsupportedSideCount { got: n });
    }
    if !(5..=9).contains(&options.degree) {
        return Err(VariationalCapError::UnsupportedDegree {
            degree: options.degree,
        });
    }
    let frame = check_sides(sides, tolerance)?;
    let local_targets = VariationalCapTargets {
        g0: targets.g0 / frame.scale,
        g1: targets.g1,
        g2: targets.g2.map(|g| g * frame.scale),
    };
    let phi = CharacteristicMap::new(n);
    let mut space = SectorSpace::bezier(n, options.degree);
    let mut gn_total = 0;
    let mut levels = 0;
    loop {
        let level = build_level(space.clone(), sides, &frame)?;
        let unknowns = 3 * level.hs.dim();
        let (coords, iterations) = solve_level(&level, &local_targets, options, &phi)?;
        gn_total += iterations;
        let rim = rim_report(&level.space, &coords, &level.samples, n);
        let worst = rim
            .iter()
            .copied()
            .max_by(|a, b| excess(a, &local_targets).total_cmp(&excess(b, &local_targets)));
        let met = rim.iter().all(|r| excess(r, &local_targets) <= 1.0);
        if met {
            return finish(level, coords, rim, levels, unknowns, gn_total, &frame);
        }
        let refined = refine(&level, &coords, &local_targets);
        let next = SectorSpace::new(n, options.degree, refined);
        let next_unknowns_estimate = unknowns * next.m * next.m / (level.space.m * level.space.m);
        if levels >= options.max_refinement_levels
            || next.knots == level.space.knots
            || next_unknowns_estimate > options.max_unknowns
        {
            let worst = worst
                .map(|w| model_residual(w, &frame))
                .unwrap_or(RimResidual {
                    side: 0,
                    half: super::RimHalf::Leading,
                    parameter: f64::NAN,
                    g0: f64::INFINITY,
                    g1: f64::INFINITY,
                    g2: f64::INFINITY,
                });
            return Err(VariationalCapError::BudgetExhausted {
                worst,
                levels,
                unknowns,
            });
        }
        space = next;
        levels += 1;
    }
}

/// Closure, degeneracy, and the normalised frame (corner centroid,
/// mean side chord).
fn check_sides(sides: &[VariationalCapSide], tolerance: f64) -> Result<Frame, VariationalCapError> {
    let n = sides.len();
    let ends: Vec<(Point3, Point3)> = sides
        .iter()
        .map(|s| {
            let (a, b) = s.curve.parameter_range();
            (s.curve.point_at(a), s.curve.point_at(b))
        })
        .collect();
    let mut chord = 0.0;
    for (s, &(a, b)) in ends.iter().enumerate() {
        let length = a.distance_to(&b);
        if length.is_nan() || length <= 1e-9 {
            return Err(VariationalCapError::DegenerateSide { side: s, length });
        }
        chord += length;
        let next = ends[(s + 1) % n].0;
        let gap = b.distance_to(&next);
        if gap > tolerance {
            return Err(VariationalCapError::BoundaryNotClosed {
                side: s,
                gap,
                tolerance,
            });
        }
    }
    let scale = chord / n as f64;
    let mut c = Vec3::new(0.0, 0.0, 0.0);
    for &(a, _) in &ends {
        c = c + Vec3::new(a.x, a.y, a.z) * (1.0 / n as f64);
    }
    Ok(Frame {
        origin: Point3::new(c.x, c.y, c.z),
        scale,
    })
}

/// Energy-only start, then Gauss–Newton on the rim terms. With
/// [`CapEnergy::DataDependent`] every iteration after the first bends
/// against the previous iterate's normal and metric.
fn solve_level(
    level: &Level,
    targets: &VariationalCapTargets,
    options: &VariationalCapOptions,
    phi: &CharacteristicMap,
) -> Result<([Vec<f64>; 3], usize), VariationalCapError> {
    let characteristic = |_s: usize, u: f64, v: f64| phi.metric(u, v);
    let parametric = Integrand::Parametric {
        metric: &characteristic,
        eps: options.membrane_weight,
    };
    let mut energy = energy_factor(level, &parametric)?;
    let red = Reduction {
        hs: &level.hs,
        xp: &level.xp,
        fixed: &level.fixed,
    };
    let mut y = solve(&energy, 1.0, None)?;
    let mut coords = assemble(level, &y);
    let mut previous = f64::INFINITY;
    let mut iterations = 0;
    for _ in 0..options.max_gauss_newton {
        if options.energy == CapEnergy::DataDependent {
            let space = &level.space;
            let current = &coords;
            let reference = |s: usize, u: f64, v: f64| {
                let j = cap_jet(space, current, s, u, v);
                let normal = j.su.cross(&j.sv).normalized()?;
                Some((normal, [j.su.dot(&j.su), j.su.dot(&j.sv), j.sv.dot(&j.sv)]))
            };
            let bending = Integrand::Bending {
                reference: &reference,
                regular: options.parameter_regularisation,
                eps: options.membrane_weight,
            };
            energy = energy_factor(level, &bending)?;
        }
        // Row weights come from the targets with a floor: a target below
        // what double precision can resolve against the energy rows would
        // only make the stacked system numerically singular, while the
        // stopping test still uses the true targets.
        let (rows, rhs) = rim_rows(
            &level.space,
            &red,
            &coords,
            &level.samples,
            targets.g1.max(WEIGHT_FLOOR_G1),
            targets.g2.map(|g| g.max(WEIGHT_FLOOR_G2)),
        );
        y = solve(&energy, options.fairness_weight, Some((&rows, &rhs)))?;
        coords = assemble(level, &y);
        iterations += 1;
        let rim = rim_report(&level.space, &coords, &level.samples, level.space.n);
        let worst = rim.iter().map(|r| excess(r, targets)).fold(0.0, f64::max);
        if !worst.is_finite() {
            continue;
        }
        if (previous - worst).abs() <= 1e-2 * worst.max(1e-300) || worst <= 1e-2 {
            break;
        }
        previous = worst;
    }
    Ok((coords, iterations))
}

/// Insert the midpoint of every span holding a sample over its target.
fn refine(level: &Level, coords: &[Vec<f64>; 3], targets: &VariationalCapTargets) -> Vec<f64> {
    let space = &level.space;
    let mut marked = vec![false; space.spans().len()];
    let spans = space.spans();
    for sample in &level.samples {
        let (g0, g1, g2) = super::measure::rim_residual(space, coords, sample);
        let r = RimResidual {
            side: sample.side,
            half: sample.half,
            parameter: sample.t,
            g0,
            g1,
            g2,
        };
        if excess(&r, targets) > 1.0
            && let Some(i) = spans
                .iter()
                .position(|&(a, b)| sample.x >= a && (sample.x < b || (b == 1.0 && sample.x <= b)))
        {
            marked[i] = true;
        }
    }
    let mut knots = space.knots.clone();
    for (i, &(a, b)) in spans.iter().enumerate() {
        if marked[i] {
            knots.push(0.5 * (a + b));
        }
    }
    knots.sort_by(f64::total_cmp);
    knots
}

fn model_residual(r: RimResidual, frame: &Frame) -> RimResidual {
    RimResidual {
        g0: r.g0 * frame.scale,
        g2: r.g2 / frame.scale,
        ..r
    }
}

fn finish(
    level: Level,
    coords: [Vec<f64>; 3],
    rim: Vec<RimResidual>,
    levels: usize,
    unknowns: usize,
    gauss_newton_iterations: usize,
    frame: &Frame,
) -> Result<VariationalCap, VariationalCapError> {
    let space = &level.space;
    let (min_regularity, sector) = regularity(space, &coords, 24);
    if min_regularity < MIN_REGULARITY {
        return Err(VariationalCapError::Folded {
            sector,
            min_regularity,
        });
    }
    let (s0, s1, s2) = spoke_report(space, &coords, 64);
    let model: [Vec<f64>; 3] = std::array::from_fn(|c| {
        let origin = [frame.origin.x, frame.origin.y, frame.origin.z][c];
        coords[c].iter().map(|v| origin + frame.scale * v).collect()
    });
    let sectors = (0..space.n)
        .map(|s| space.sector_surface(&model, s))
        .collect::<Vec<_>>();
    let centre = sectors[0].control_points[0];
    let m = space.m;
    let spoke_points = (0..space.n)
        .map(|s| sectors[s].control_points[(m - 1) * m])
        .collect();
    let energy = energy_value(&level, &coords);
    Ok(VariationalCap {
        sectors,
        centre,
        spoke_points,
        report: VariationalCapReport {
            knots: space.knots.clone(),
            refinement_levels: levels,
            unknowns,
            gauss_newton_iterations,
            rim: rim.into_iter().map(|r| model_residual(r, frame)).collect(),
            worst_spoke: (s0 * frame.scale, s1, s2 / frame.scale),
            min_regularity,
            energy,
        },
    })
}

/// Thin-plate energy of the result under the characteristic metric
/// (normalised units), for comparison across fits.
fn energy_value(level: &Level, coords: &[Vec<f64>; 3]) -> f64 {
    let phi = CharacteristicMap::new(level.space.n);
    let metric = |_: usize, u: f64, v: f64| phi.metric(u, v);
    let gram = super::energy::energy_gram(
        &level.space,
        &Integrand::Parametric {
            metric: &metric,
            eps: 0.0,
        },
    );
    let per = level.space.m * level.space.m;
    let size = gram.size;
    let mut total = 0.0;
    for (s, block) in gram.blocks.iter().enumerate() {
        let x: Vec<f64> = (0..3)
            .flat_map(|c| coords[c][s * per..(s + 1) * per].iter().copied())
            .collect();
        for i in 0..size {
            let row: f64 = (0..size).map(|j| block[i * size + j] * x[j]).sum();
            total += x[i] * row;
        }
    }
    total
}
