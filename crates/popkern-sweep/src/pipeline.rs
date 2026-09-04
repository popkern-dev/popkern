//! Shared, stage-instrumented pipeline bodies for the sweep-family
//! entry points (PK-276). Each public function in `lib.rs` (and each
//! `*_with_control` wrapper in `interactive.rs`) delegates here; the
//! `tracker` parameter threads cooperative cancellation polls and
//! per-stage progress events through every ladder.

use popkern_core::geometry::NurbsCurve;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::history::TrackingContext;
use popkern_core::tolerance::Tolerance;

use crate::error::SweepError;
use crate::error::SweepWarning;
use crate::types::{LoftOptions, LoftSection, PipeRadius, SweepOptions, SweepProfile, SweepResult};
use crate::{arc_length, frame, pipe, profile, self_intersect, surface, topology};

/// Honest per-run stage count of [`sweep_impl`] for progress totals
/// (PK-276).
pub(crate) const SWEEP_STAGE_COUNT: usize = 6;

/// Shared implementation for [`sweep`] (and the helical delegate).
/// `tracker` threads PK-276 interactive control: a cancellation poll
/// before every stage (`SweepError::Cancelled`, nothing produced) and
/// one progress event per stage; the self-intersection stage also
/// polls per surface.
pub(crate) fn sweep_impl(
    profile: &SweepProfile,
    path: &NurbsCurve,
    options: &SweepOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
    tracker: &mut popkern_core::progress::StageTracker<'_, '_>,
) -> Result<SweepResult, SweepError> {
    macro_rules! stage_gate {
        () => {
            if tracker.cancelled() {
                return Err(SweepError::Cancelled);
            }
        };
    }
    // 1. Analyze profile.
    stage_gate!();
    let t0 = popkern_core::time::Stopwatch::start();
    let analyzed = profile::analyze_profile(profile, &tolerance)?;
    tracker.finish_stage("sweep.analyze_profile", t0, analyzed.spans.len());

    // 2. Build arc-length table.
    stage_gate!();
    let t0 = popkern_core::time::Stopwatch::start();
    let table = arc_length::ArcLengthTable::build(path, ARC_LENGTH_SAMPLES);
    tracker.finish_stage("sweep.arc_length_table", t0, ARC_LENGTH_SAMPLES);

    // 3. Compute frames.
    stage_gate!();
    let t0 = popkern_core::time::Stopwatch::start();
    let n_samples = compute_sample_count(&table);
    let is_closed_path = check_path_closed(path, &tolerance);
    let (frames, mut warnings) = frame::compute_frames(
        path,
        &table,
        &options.frame,
        n_samples,
        options.twist,
        is_closed_path,
    )?;
    tracker.finish_stage("sweep.frames", t0, frames.frames.len());

    // 4. Generate surfaces.
    stage_gate!();
    let t0 = popkern_core::time::Stopwatch::start();
    let surface_result =
        surface::generate_sweep_surface(&analyzed, &frames, &options.scale, &tolerance)?;
    tracker.finish_stage("sweep.surfaces", t0, surface_result.surfaces.len());

    // 5. Self-intersection check (per-surface poll: the check is the
    // costliest per-entity step of the pipeline).
    let t0 = popkern_core::time::Stopwatch::start();
    for surf in &surface_result.surfaces {
        stage_gate!();
        if let Some(t) = self_intersect::check_self_intersection(surf, &tolerance) {
            if options.allow_self_intersection {
                warnings.push(SweepWarning::SelfIntersectionAllowed { t });
            } else {
                return Err(SweepError::SelfIntersection { t });
            }
        }
    }
    tracker.finish_stage("sweep.self_intersection", t0, surface_result.surfaces.len());

    // 6. Build topology.
    stage_gate!();
    let t0 = popkern_core::time::Stopwatch::start();
    let brep = topology::side_faces::build_sweep_topology(
        &surface_result.surfaces,
        &surface_result.start_cap,
        &surface_result.end_cap,
        analyzed.is_closed,
        tolerance,
        tracking,
    )?;
    tracker.finish_stage("sweep.topology", t0, brep.faces.len());

    Ok(SweepResult { brep, warnings })
}

/// Honest per-run stage count of [`loft_impl`] for progress totals
/// (PK-276). The guide-validation stage is skipped (its event never
/// emitted) when no guides are supplied.
pub(crate) const LOFT_STAGE_COUNT: usize = 3;

/// Shared implementation for [`loft`] and [`guided_loft`] (`guides`
/// empty = plain loft). `tracker` threads PK-276 interactive control.
pub(crate) fn loft_impl(
    sections: &[LoftSection],
    guides: &[NurbsCurve],
    options: &LoftOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
    tracker: &mut popkern_core::progress::StageTracker<'_, '_>,
) -> Result<SweepResult, SweepError> {
    macro_rules! stage_gate {
        () => {
            if tracker.cancelled() {
                return Err(SweepError::Cancelled);
            }
        };
    }
    stage_gate!();
    let t0 = popkern_core::time::Stopwatch::start();
    // Span-wise lofting: when every section splits into the SAME
    // number of smooth spans (at tangent-discontinuity knots — a
    // stadium's line/arc joints, a polygon's corners), each span
    // family lofts to its own surface, exactly as a swept profile's
    // spans do. Arcs stay arcs per span, and the wall becomes N
    // ordinary bounded faces instead of ONE u-periodic face — a
    // periodic wall cut by a non-contractible intersection loop is a
    // fragmentation class the boolean pipelines do not serve yet
    // (PK-440: the propeller blade's girdle). Sections with
    // mismatched span structures keep the single-surface path.
    let span_sets: Vec<Vec<popkern_core::geometry::NurbsCurve>> = sections
        .iter()
        .map(|s| crate::profile::detect_and_split(&s.curve))
        .collect::<Result<Vec<_>, _>>()?;
    let span_wise =
        span_sets[0].len() > 1 && span_sets.iter().all(|v| v.len() == span_sets[0].len());
    let (surfaces, start_caps, end_caps) = if span_wise {
        let n_spans = span_sets[0].len();
        let mut surfaces = Vec::with_capacity(n_spans);
        let mut start_caps = Vec::with_capacity(n_spans);
        let mut end_caps = Vec::with_capacity(n_spans);
        for j in 0..n_spans {
            let span_sections: Vec<LoftSection> = span_sets
                .iter()
                .map(|set| LoftSection {
                    curve: set[j].clone(),
                    position: None,
                })
                .collect();
            let r = surface::generate_loft_surface(&span_sections, options, &tolerance)?;
            surfaces.push(r.surface);
            start_caps.push(r.start_cap);
            end_caps.push(r.end_cap);
        }
        (surfaces, start_caps, end_caps)
    } else {
        let result = surface::generate_loft_surface(sections, options, &tolerance)?;
        (
            vec![result.surface],
            vec![result.start_cap],
            vec![result.end_cap],
        )
    };
    tracker.finish_stage("loft.surface", t0, sections.len());

    if !guides.is_empty() {
        stage_gate!();
        let t0 = popkern_core::time::Stopwatch::start();
        for surface in &surfaces {
            surface::validate_guide_rails(surface, guides, &tolerance)?;
        }
        tracker.finish_stage("loft.guide_rails", t0, guides.len());
    }

    // Closure must be consistent: the topology builder caps (or not)
    // the whole loft from one answer, so a mix of open and closed
    // sections would silently build garbage.
    let section_closed = |s: &LoftSection| {
        let (t_min, t_max) = s.curve.parameter_range();
        let start = s.curve.point_at(t_min);
        let end = s.curve.point_at(t_max);
        tolerance.points_coincident(&start, &end)
    };
    let is_closed = sections.first().map(section_closed).unwrap_or(false);
    if let Some(mixed) = sections.iter().position(|s| section_closed(s) != is_closed) {
        return Err(SweepError::InvalidProfile(format!(
            "loft sections mix open and closed profiles (section 0 is \
             {}, section {mixed} is not)",
            if is_closed { "closed" } else { "open" }
        )));
    }

    stage_gate!();
    let t0 = popkern_core::time::Stopwatch::start();
    let brep = topology::side_faces::build_sweep_topology(
        &surfaces,
        &start_caps,
        &end_caps,
        is_closed,
        tolerance,
        tracking,
    )?;
    tracker.finish_stage("loft.topology", t0, brep.faces.len());

    Ok(SweepResult {
        brep,
        warnings: vec![],
    })
}

/// Honest per-run stage count of [`pipe_impl`] for progress totals
/// (PK-276).
pub(crate) const PIPE_STAGE_COUNT: usize = 2;

/// Shared implementation for [`pipe`]. `tracker` threads PK-276
/// interactive control.
pub(crate) fn pipe_impl(
    path: &NurbsCurve,
    radius: &PipeRadius,
    options: &SweepOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
    tracker: &mut popkern_core::progress::StageTracker<'_, '_>,
) -> Result<SweepResult, SweepError> {
    if tracker.cancelled() {
        return Err(SweepError::Cancelled);
    }
    let t0 = popkern_core::time::Stopwatch::start();
    let surface_result = pipe::generate_pipe_surface(path, radius, options, &tolerance)?;
    tracker.finish_stage("pipe.surface", t0, surface_result.surfaces.len());

    if tracker.cancelled() {
        return Err(SweepError::Cancelled);
    }
    let t0 = popkern_core::time::Stopwatch::start();
    let brep = topology::side_faces::build_sweep_topology(
        &surface_result.surfaces,
        &surface_result.start_cap,
        &surface_result.end_cap,
        true, // pipes are always closed profiles
        tolerance,
        tracking,
    )?;
    tracker.finish_stage("pipe.topology", t0, brep.faces.len());

    Ok(SweepResult {
        brep,
        warnings: vec![],
    })
}

/// Honest per-run stage count of [`guided_sweep_impl`] for progress
/// totals (PK-276).
pub(crate) const GUIDED_SWEEP_STAGE_COUNT: usize = 6;

/// Shared implementation for [`guided_sweep`]. `tracker` threads
/// PK-276 interactive control; the guide-validation stage polls per
/// surface.
pub(crate) fn guided_sweep_impl(
    profile: &SweepProfile,
    path: &NurbsCurve,
    guides: &[NurbsCurve],
    options: &SweepOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
    tracker: &mut popkern_core::progress::StageTracker<'_, '_>,
) -> Result<SweepResult, SweepError> {
    macro_rules! stage_gate {
        () => {
            if tracker.cancelled() {
                return Err(SweepError::Cancelled);
            }
        };
    }
    // Generate the base sweep.
    stage_gate!();
    let t0 = popkern_core::time::Stopwatch::start();
    let analyzed = profile::analyze_profile(profile, &tolerance)?;
    tracker.finish_stage("guided_sweep.analyze_profile", t0, analyzed.spans.len());

    stage_gate!();
    let t0 = popkern_core::time::Stopwatch::start();
    let table = arc_length::ArcLengthTable::build(path, ARC_LENGTH_SAMPLES);
    tracker.finish_stage("guided_sweep.arc_length_table", t0, ARC_LENGTH_SAMPLES);

    stage_gate!();
    let t0 = popkern_core::time::Stopwatch::start();
    let n_samples = compute_sample_count(&table);
    let is_closed_path = check_path_closed(path, &tolerance);
    let (frames, warnings) = frame::compute_frames(
        path,
        &table,
        &options.frame,
        n_samples,
        options.twist,
        is_closed_path,
    )?;
    tracker.finish_stage("guided_sweep.frames", t0, frames.frames.len());

    stage_gate!();
    let t0 = popkern_core::time::Stopwatch::start();
    let surface_result =
        surface::generate_sweep_surface(&analyzed, &frames, &options.scale, &tolerance)?;
    tracker.finish_stage("guided_sweep.surfaces", t0, surface_result.surfaces.len());

    // Validate guides against each surface (per-surface poll).
    let t0 = popkern_core::time::Stopwatch::start();
    for surf in &surface_result.surfaces {
        stage_gate!();
        surface::validate_guide_rails(surf, guides, &tolerance)?;
    }
    tracker.finish_stage("guided_sweep.guide_rails", t0, guides.len());

    stage_gate!();
    let t0 = popkern_core::time::Stopwatch::start();
    let brep = topology::side_faces::build_sweep_topology(
        &surface_result.surfaces,
        &surface_result.start_cap,
        &surface_result.end_cap,
        analyzed.is_closed,
        tolerance,
        tracking,
    )?;
    tracker.finish_stage("guided_sweep.topology", t0, brep.faces.len());

    Ok(SweepResult { brep, warnings })
}

/// Honest per-run stage count of [`swept_blend_impl`] for progress
/// totals (PK-276).
pub(crate) const SWEPT_BLEND_STAGE_COUNT: usize = 2;

/// Shared implementation for [`swept_blend`]. `tracker` threads
/// PK-276 interactive control.
pub(crate) fn swept_blend_impl(
    profiles: &[(f64, SweepProfile)],
    path: &NurbsCurve,
    options: &SweepOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
    tracker: &mut popkern_core::progress::StageTracker<'_, '_>,
) -> Result<SweepResult, SweepError> {
    if tracker.cancelled() {
        return Err(SweepError::Cancelled);
    }
    let t0 = popkern_core::time::Stopwatch::start();
    let result =
        surface::swept_blend::generate_swept_blend_surface(profiles, path, options, &tolerance)?;
    tracker.finish_stage("swept_blend.surface", t0, profiles.len());

    if tracker.cancelled() {
        return Err(SweepError::Cancelled);
    }
    let t0 = popkern_core::time::Stopwatch::start();
    let brep = topology::side_faces::build_sweep_topology(
        &[result.surface],
        std::slice::from_ref(&result.start_cap),
        std::slice::from_ref(&result.end_cap),
        false, // swept blend sections may not be closed
        tolerance,
        tracking,
    )?;
    tracker.finish_stage("swept_blend.topology", t0, brep.faces.len());

    Ok(SweepResult {
        brep,
        warnings: vec![],
    })
}

/// Default number of arc-length samples for path parameterization.
pub(crate) const ARC_LENGTH_SAMPLES: usize = 100;

/// Minimum number of frames for sweep operations.
pub(crate) const MIN_FRAMES: usize = 10;

/// Compute sample count based on path arc length.
pub(crate) fn compute_sample_count(table: &arc_length::ArcLengthTable) -> usize {
    // Roughly 10 samples per unit length, minimum MIN_FRAMES.
    let count = (table.total_length * 10.0).ceil() as usize;
    count.clamp(MIN_FRAMES, 200)
}

/// Check if a path curve is closed (start point equals end point).
pub(crate) fn check_path_closed(path: &NurbsCurve, tolerance: &Tolerance) -> bool {
    let (t_min, t_max) = path.parameter_range();
    let start = path.point_at(t_min);
    let end = path.point_at(t_max);
    tolerance.points_coincident(&start, &end)
}
