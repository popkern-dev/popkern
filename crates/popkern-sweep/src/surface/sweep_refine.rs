//! Sweep and pipe surfaces skinned densely enough to hold the swept
//! profile within the linear tolerance.
//!
//! A swept surface interpolates the profile placed at each frame; between
//! frames it is a fit. The frame count came from a fixed rule (17 for a
//! pipe, about ten per unit length for a sweep), and a pipe of radius 20
//! along a cubic of length 52 read its volume 1.2e-6 short of πr²L, its
//! wall bulging between frames. The frames are now doubled (n → 2n − 1,
//! the old frames kept and one placed between each pair) until the
//! profile placed at every midpoint lies on the surface skinned through
//! the others within the tolerance, up to `MAX_FRAMES`.

use popkern_core::geometry::NurbsCurve;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::tolerance::Tolerance;

use super::sweep_surface::{
    SweepSurfaceResult, generate_sweep_surface, transform_profile_at_frame,
};
use crate::arc_length::ArcLengthTable;
use crate::error::{SweepError, SweepWarning};
use crate::frame::{FrameSequence, compute_frames};
use crate::profile::AnalyzedProfile;
use crate::types::FrameMethod;
use crate::types::ScaleMethod;

/// The most frames a refined sweep skins through.
const MAX_FRAMES: usize = 1025;
/// Profile samples per span checked at each midpoint frame.
const PROFILE_SAMPLES: usize = 5;

/// Whether a sweep's frames are refined until its skin holds the profile.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum FramePolicy {
    /// Doubled until every midpoint profile lies on the skin.
    Refined,
    /// The frame budget's count as it stands. Helical sweeps (threads,
    /// coils) keep it: their frame counts are tuned per knot span with
    /// the downstream clip and phase bands (`frame_budget.rs`), and
    /// refining them to the tolerance would run to `MAX_FRAMES` on
    /// every multi-turn thread (PK-1208).
    Budgeted,
}

/// Frames and surfaces for sweeping `profile` along `path`, starting from
/// `n_frames` and doubling until the surfaces hold the profile at every
/// midpoint within `tolerance.linear()`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn refined_sweep_surface(
    profile: &AnalyzedProfile,
    path: &NurbsCurve,
    table: &ArcLengthTable,
    method: &FrameMethod,
    n_frames: usize,
    twist: f64,
    is_closed: bool,
    scale: &ScaleMethod,
    tolerance: &Tolerance,
    policy: FramePolicy,
) -> Result<(FrameSequence, Vec<SweepWarning>, SweepSurfaceResult), SweepError> {
    if policy == FramePolicy::Budgeted || prism(path, method, twist, scale) {
        let (frames, warnings) = compute_frames(path, table, method, n_frames, twist, is_closed)?;
        let surface = generate_sweep_surface(profile, &frames, scale, tolerance)?;
        return Ok((frames, warnings, surface));
    }
    let mut n = n_frames.max(2);
    loop {
        let fine_count = 2 * n - 1;
        if fine_count > MAX_FRAMES {
            let (frames, warnings) = compute_frames(path, table, method, n, twist, is_closed)?;
            let surface = generate_sweep_surface(profile, &frames, scale, tolerance)?;
            return Ok((frames, warnings, surface));
        }
        let (fine, warnings) = compute_frames(path, table, method, fine_count, twist, is_closed)?;
        let coarse = FrameSequence {
            frames: fine.frames.iter().step_by(2).cloned().collect(),
            s_values: fine.s_values.iter().step_by(2).copied().collect(),
        };
        let surface = generate_sweep_surface(profile, &coarse, scale, tolerance)?;
        if midpoint_deviation(profile, &fine, scale, &surface) <= tolerance.linear() {
            return Ok((coarse, warnings, surface));
        }
        n = fine_count;
    }
}

/// Whether the sweep is a prism: one straight span, no twist, no scale
/// and a frame that cannot turn. Its skin is exact through any frames,
/// and extrudes built this way keep their frames bit for bit (parity
/// fixtures, and a fillet on a 1.6 mm block read its rounds 0.35% low
/// once the block's walls were reskinned).
fn prism(path: &NurbsCurve, method: &FrameMethod, twist: f64, scale: &ScaleMethod) -> bool {
    path.degree == 1
        && path.control_points.len() == 2
        && twist == 0.0
        && matches!(scale, ScaleMethod::None)
        && !matches!(method, FrameMethod::Auxiliary(_))
}

/// The largest distance from the profile placed at an odd (midpoint)
/// frame of `fine` to the surface skinned through the even ones.
fn midpoint_deviation(
    profile: &AnalyzedProfile,
    fine: &FrameSequence,
    scale: &ScaleMethod,
    surface: &SweepSurfaceResult,
) -> f64 {
    let total = fine.s_values.last().copied().unwrap_or(1.0).max(1e-30);
    let mut worst = 0.0_f64;
    for k in (1..fine.frames.len()).step_by(2) {
        let v = fine.s_values[k] / total;
        let placed = transform_profile_at_frame(&profile.spans, &fine.frames[k], scale, v);
        // The skin keeps each span's own u and runs v by arc length, so
        // the profile point at u lies at (u, v) on an exact skin.
        for (curve, skin) in placed.iter().zip(&surface.surfaces) {
            let (t0, t1) = curve.parameter_range();
            for i in 0..=PROFILE_SAMPLES {
                let u = t0 + (t1 - t0) * i as f64 / PROFILE_SAMPLES as f64;
                worst = worst.max(skin.point_at(u, v).distance_to(&curve.point_at(u)));
            }
        }
    }
    worst
}
