//! The driver: validation, normalisation, levels, Gauss–Newton,
//! refinement, output.

use super::chart::CharacteristicMap;
use super::energy::Integrand;
use super::gauss_newton::solve_level;
use super::layout;
use super::level::{Frame, Level, build_level};
use super::measure::{
    SERVED_MIN_REGULARITY, THIN_SECTOR_MIN_REGULARITY, rim_report, served_regularity, spoke_report,
};
use super::rim_rows::CrossDerivativeTarget;
use super::space::SectorSpace;
use super::{
    RimResidual, VariationalCap, VariationalCapError, VariationalCapOptions, VariationalCapReport,
    VariationalCapSide, VariationalCapTargets, excess,
};
use crate::geometry::curve::CurveEval;
use crate::geometry::{Point3, Vec3};

/// Spans within this factor of the worst span excess are split.
const REFINE_BAND: f64 = 4.0;

/// The most sides the fit takes. The characteristic chart and the
/// spoke conditions hold for every n >= 3; 24 covers the bridged
/// setback border (two sides per strip, PK-611) of a twelve-strip apex.
pub const MAX_CAP_SIDES: usize = 24;

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
    if !(3..=MAX_CAP_SIDES).contains(&n) {
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
    let midpoints = layout::midpoints(sides);
    let first =
        |cross| fit_with_splits(sides, options, &frame, &local_targets, (&midpoints, cross));
    let Some(limit) = options.layout_refit_curvature else {
        return first(CrossDerivativeTarget::InPlaneLength);
    };
    let read = |mut cap: VariationalCap, layout, cross| {
        cap.report.interior_curvature = layout::interior_curvature(&cap);
        cap.report.layout = layout;
        cap.report.cross_target = cross;
        cap
    };
    // A cap bending past the caller's limit is refit, and the fairest
    // kept (the earlier on a tie; a refit that declines is dropped); a
    // cap inside the limit keeps its fit unchanged.
    // - PK-625 (when the caller asks): the transversality prior at the
    //   start's full length, for a cap that turns into a neighbour inside
    //   a thin rim layer. It is fitted alongside the first fit, and is
    //   wasted work only when the first fit is fair.
    // - PK-618: the spokes' feet nearest the first fit's centre, while
    //   the fairest is still past the limit (measured over the 56
    //   corner-patched pieces of the 160-row census: the start-length
    //   refit was fairest 38 times, the centre-feet refit once, by 1%).
    let candidates = if options.cross_target_refit { 2 } else { 1 };
    let mut fits = super::linalg::par_map(candidates, |t| {
        first(if t == 0 {
            CrossDerivativeTarget::InPlaneLength
        } else {
            CrossDerivativeTarget::StartLength
        })
    })
    .into_iter();
    let first_fit = fits.next().ok_or(VariationalCapError::Singular {
        stage: "no level solved",
    })??;
    let first_fit = read(
        first_fit,
        layout::CapLayout::Midpoints,
        CrossDerivativeTarget::InPlaneLength,
    );
    if first_fit.report.interior_curvature <= limit {
        return Ok(first_fit);
    }
    let feet = layout::centre_feet(sides, &first_fit);
    let mut cap = first_fit;
    if let Some(Ok(refit)) = fits.next() {
        let refit = read(
            refit,
            layout::CapLayout::Midpoints,
            CrossDerivativeTarget::StartLength,
        );
        if refit.report.interior_curvature < cap.report.interior_curvature {
            cap = refit;
        }
    }
    if cap.report.interior_curvature <= limit {
        return Ok(cap);
    }
    if let Ok(refit) = fit_with_splits(
        sides,
        options,
        &frame,
        &local_targets,
        (&feet, CrossDerivativeTarget::InPlaneLength),
    ) {
        let refit = read(
            refit,
            layout::CapLayout::CentreFeet,
            CrossDerivativeTarget::InPlaneLength,
        );
        if refit.report.interior_curvature < cap.report.interior_curvature {
            cap = refit;
        }
    }
    Ok(cap)
}

/// One fit (levels, Gauss–Newton, refinement) with side `s` split at
/// `splits[s]`.
fn fit_with_splits(
    sides: &[VariationalCapSide],
    options: &VariationalCapOptions,
    frame: &Frame,
    local_targets: &VariationalCapTargets,
    (splits, cross): (&[f64], CrossDerivativeTarget),
) -> Result<VariationalCap, VariationalCapError> {
    let n = sides.len();
    let frame = *frame;
    let local_targets = *local_targets;
    let phi = CharacteristicMap::new(n);
    let mut space = SectorSpace::bezier(n, options.degree);
    if options.corner_cells {
        space = space.with_corner_cells();
    }
    let mut gn_total = 0;
    let mut levels = 0;
    // The best level so far by worst residual-to-target ratio: a level
    // that refines past the targets' reach is not better than its
    // predecessor, so the result is the best one seen.
    let mut best: Option<(f64, VariationalCap)> = None;
    let mut previous: Option<(SectorSpace, [Vec<f64>; 3])> = None;
    loop {
        let level = build_level(space.clone(), sides, &frame, splits)?;
        let unknowns = 3 * level.hs.dim();
        let start = previous
            .as_ref()
            .and_then(|(old_space, old)| transfer(old_space, old, &level.space));
        let (coords, iterations) =
            solve_level(&level, &local_targets, options, &phi, start.as_ref(), cross)?;
        gn_total += iterations;
        let rim = rim_report(&level.space, &coords, &level.samples, n);
        let worst = rim
            .iter()
            .map(|r| excess(r, &local_targets))
            .fold(0.0, f64::max);
        let met = worst <= 1.0;
        let refined = (!met).then(|| refine(&level, &coords, &local_targets));
        let space_now = level.space.clone();
        previous = Some((space_now.clone(), coords.clone()));
        let cap = finish(
            level,
            coords,
            rim,
            levels,
            unknowns,
            gn_total,
            &frame,
            met,
            options.thin_sector_curvature,
        )?;
        if best.as_ref().is_none_or(|(w, _)| worst < *w) {
            best = Some((worst, cap));
        }
        let Some(refined) = refined else {
            break;
        };
        let next = SectorSpace::new(n, options.degree, refined);
        let next_unknowns_estimate = unknowns * next.m * next.m / (space_now.m * space_now.m);
        if levels >= options.max_refinement_levels
            || space_now.corner_cells
            || next.knots == space_now.knots
            || next_unknowns_estimate > options.max_unknowns
        {
            break;
        }
        space = next;
        levels += 1;
    }
    best.map(|(_, cap)| cap)
        .ok_or(VariationalCapError::Singular {
            stage: "no level solved",
        })
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

/// Insert the midpoint of every span holding a sample over its target.
fn refine(level: &Level, coords: &[Vec<f64>; 3], targets: &VariationalCapTargets) -> Vec<f64> {
    // Split the spans holding the worst samples: those within a factor
    // REFINE_BAND of the level's worst excess (and over target). The
    // residual that survives a level concentrates where the data
    // conflict (the setback corners, where two G2 strips' curvature
    // slopes are third-order incompatible for any polynomial cap), so
    // grading the knots there converges; splitting every over-target
    // span spends the budget uniformly.
    let space = &level.space;
    let spans = space.spans();
    let mut span_excess = vec![0.0f64; spans.len()];
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
        if let Some(i) = spans
            .iter()
            .position(|&(a, b)| sample.x >= a && (sample.x < b || (b == 1.0 && sample.x <= b)))
        {
            span_excess[i] = span_excess[i].max(excess(&r, targets));
        }
    }
    let worst = span_excess.iter().copied().fold(0.0, f64::max);
    let mut knots = space.knots.clone();
    for (i, &(a, b)) in spans.iter().enumerate() {
        if span_excess[i] > 1.0 && span_excess[i] >= worst / REFINE_BAND {
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

#[allow(clippy::too_many_arguments)]
fn finish(
    level: Level,
    coords: [Vec<f64>; 3],
    rim: Vec<RimResidual>,
    levels: usize,
    unknowns: usize,
    gauss_newton_iterations: usize,
    frame: &Frame,
    targets_met: bool,
    thin_limit: Option<f64>,
) -> Result<VariationalCap, VariationalCapError> {
    let space = &level.space;
    let reg = served_regularity(space, &coords);
    let min_regularity = reg.ratio;
    let thin = min_regularity < SERVED_MIN_REGULARITY;
    if thin && (min_regularity < THIN_SECTOR_MIN_REGULARITY || thin_limit.is_none()) {
        return Err(VariationalCapError::Folded {
            sector: reg.sector,
            at: reg.uv,
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
    // Clamped corners are control points of either form (the rational
    // corner-cell sectors have their own net size, PK-595).
    let faces = (0..space.n)
        .map(|s| space.face_surface(&model, s))
        .collect::<Vec<_>>();
    let centre = sectors[0].control_points[0];
    let spoke_points = (0..space.n)
        .map(|s| sectors[s].control_points[(sectors[s].num_u - 1) * sectors[s].num_v])
        .collect();
    let energy = energy_value(&level, &coords);
    let cap = VariationalCap {
        sectors,
        faces,
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
            targets_met,
            layout: layout::CapLayout::Midpoints,
            cross_target: CrossDerivativeTarget::InPlaneLength,
            interior_curvature: 0.0,
        },
    };
    // A thin sector serves only when the cap stays fair.
    if let (true, Some(limit)) = (thin, thin_limit) {
        let curvature = layout::interior_curvature(&cap);
        if curvature > limit {
            return Err(VariationalCapError::ThinSectorBends {
                sector: reg.sector,
                at: reg.uv,
                min_regularity,
                curvature,
                limit,
            });
        }
    }
    Ok(cap)
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
        None,
    );
    let per = level.space.per();
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

/// The previous level's coefficients in the refined space `new` (exact
/// knot insertion per sector); `None` if the insertion fails.
pub(super) fn transfer(
    old_space: &SectorSpace,
    old: &[Vec<f64>; 3],
    new: &SectorSpace,
) -> Option<[Vec<f64>; 3]> {
    if old_space.corner_cells || new.corner_cells {
        return None;
    }
    let mut inserted = Vec::new();
    let mut remaining = old_space.knots.clone();
    for &k in &new.knots {
        match remaining.iter().position(|&r| r == k) {
            Some(i) => {
                remaining.remove(i);
            }
            None => inserted.push(k),
        }
    }
    let mut out: [Vec<f64>; 3] = std::array::from_fn(|_| vec![0.0; new.total()]);
    for s in 0..old_space.n {
        let surface = old_space
            .sector_surface(old, s)
            .insert_knots_u(&inserted)
            .ok()?
            .insert_knots_v(&inserted)
            .ok()?;
        if surface.control_points.len() != new.m * new.m {
            return None;
        }
        for (l, p) in surface.control_points.iter().enumerate() {
            let at = s * new.per() + l;
            out[0][at] = p.x;
            out[1][at] = p.y;
            out[2][at] = p.z;
        }
    }
    Some(out)
}
