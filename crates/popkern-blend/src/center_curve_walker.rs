//! Walker-based center-curve builder. Walker cone path: routes
//! around the offset-intersect short-cut (which has no analytic
//! offset for cones) by using the rolling-ball contact solver in
//! continuation along the spine.
//!
//! Split out of `center_curve.rs` to keep the parent under the
//! 500-line structural cap.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;

use crate::center_curve::CenterCurveResult;
use crate::center_curve_helpers::{clip_footprint_to_bounds, fit_nurbs_through_points};
use crate::walker::{WalkOptions, WalkRequest, walk_spine};

/// Build a center curve via the spine walker, given an explicit
/// spine curve and parameter range.
///
/// The solver's first sample seeds via `closest_point`; every
/// subsequent sample warm-starts from the previous solution, so the
/// inner-loop cost on cones drops from "20 × multi-start Gauss-Newton"
/// to "20 × ~2 Newton iterations."
///
/// Returns the same [`CenterCurveResult`] shape as the vertex-blend paths,
/// with `center_curve` fitted through the sample ball-centers and
/// `footprint_*` populated from each sample's contact points.
#[allow(clippy::too_many_arguments)]
pub fn walker_based_center_curve(
    surface1: &SurfaceGeometry,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    surface2: &SurfaceGeometry,
    u2_range: (f64, f64),
    v2_range: (f64, f64),
    radius_at: impl Fn(f64) -> f64,
    spine: &dyn CurveEval,
    t_range: (f64, f64),
    tolerance: &Tolerance,
    orientation_1: bool,
    orientation_2: bool,
    convexity: crate::analysis::Convexity,
    pcurve_seed: Option<((f64, f64), (f64, f64))>,
    n_samples: usize,
    refine_stations: bool,
) -> Result<CenterCurveResult, String> {
    // Ball side depends on the edge class (PK-091): inscribed in the
    // material for convex edges, in the void for concave. The old
    // orientation-only mapping here hardcoded the concave case and
    // made the walker converge to the exterior mirror ball on every
    // convex edge.
    let (side_0, side_1) =
        crate::contact::ball_sides_for_edge(convexity, orientation_1, orientation_2);

    // The walker takes `radius_at` as a closure of the spine
    // parameter `t`, but our `radius_at` argument receives the
    // normalized fillet parameter [0, 1]. Adapt by mapping
    // `t_spine ∈ [t_start, t_end]` → `t_norm ∈ [0, 1]`.
    let (t_start, t_end) = t_range;
    let t_span = t_end - t_start;
    let radius_at_spine_t = |t_spine: f64| -> f64 {
        let t_norm = if t_span.abs() < 1.0e-30 {
            0.0
        } else {
            (t_spine - t_start) / t_span
        };
        radius_at(t_norm.clamp(0.0, 1.0))
    };

    let request = WalkRequest {
        surface_0: surface1,
        surface_1: surface2,
        side_0,
        side_1,
        radius_at: &radius_at_spine_t,
        spine,
        t_range,
        uv0_range: (u1_range, v1_range),
        uv1_range: (u2_range, v2_range),
        // Use the pcurve-derived seed when available — it gives
        // the exact UV at the spine start on each face, avoiding
        // the wrapped-UV problem from closest_point on periodic
        // surfaces (§12.3 blend-on-blend).
        initial_seed: pcurve_seed,
        surface_0_periodicity: surface1.periodicity(),
        surface_1_periodicity: surface2.periodicity(),
    };
    // §12.3: when a pcurve seed is available, enable UV clamping
    // as a safety net for periodic surfaces. LM damping helps
    // convergence from pcurve-derived starts that are on the
    // face boundary rather than at the exact contact point.
    let mut walk_opts = WalkOptions {
        mode: crate::walker::WalkMode::Fixed {
            n_samples: n_samples.max(2),
        },
        ..WalkOptions::default()
    };
    let has_curved_support = !matches!(surface1, SurfaceGeometry::Plane { .. })
        || !matches!(surface2, SurfaceGeometry::Plane { .. });
    if pcurve_seed.is_some() || has_curved_support {
        walk_opts.contact_options.clamp_initial_uv = true;
        walk_opts.contact_options.initial_damping = 0.1;
        walk_opts.contact_options.max_iterations = 40;
    }
    let mut walk = walk_spine(&request, &walk_opts, tolerance)
        .map_err(|e| format!("spine walker failed: {e}"))?;
    if refine_stations {
        walk = refine_station_density(&request, walk_opts, walk, n_samples.max(2), tolerance);
    }

    if popkern_core::trace::trace_enabled("walker")
        && let (Some(first), Some(last)) = (walk.samples.first(), walk.samples.last())
    {
        let (a, b) = (spine.point_at(t_start), spine.point_at(t_end));
        eprintln!(
            "[walker] {} samples; spine ends ({:.4},{:.4},{:.4})->({:.4},{:.4},{:.4}); first t {:.4} c ({:.4},{:.4},{:.4}) last t {:.4} c ({:.4},{:.4},{:.4})",
            walk.samples.len(),
            a.x,
            a.y,
            a.z,
            b.x,
            b.y,
            b.z,
            first.t_spine,
            first.solution.center.x,
            first.solution.center.y,
            first.solution.center.z,
            last.t_spine,
            last.solution.center.x,
            last.solution.center.y,
            last.solution.center.z
        );
    }
    let mut center_points: Vec<Point3> = Vec::with_capacity(walk.samples.len());
    let mut footprint_1: Vec<Point3> = Vec::with_capacity(walk.samples.len());
    let mut footprint_2: Vec<Point3> = Vec::with_capacity(walk.samples.len());
    let mut contact_solutions = Vec::with_capacity(walk.samples.len());
    for sample in &walk.samples {
        center_points.push(sample.solution.center);
        footprint_1.push(sample.solution.contact0);
        footprint_2.push(sample.solution.contact1);
        contact_solutions.push(sample.solution.clone());
    }

    let center_curve = fit_nurbs_through_points(&center_points)?;
    let footprint_1 = clip_footprint_to_bounds(surface1, u1_range, v1_range, &footprint_1);
    let footprint_2 = clip_footprint_to_bounds(surface2, u2_range, v2_range, &footprint_2);

    Ok(CenterCurveResult {
        center_curve,
        sweep_center_curve: None,
        sweep_footprint_1: None,
        sweep_footprint_2: None,
        footprint_1,
        footprint_2,
        contact_solutions: Some(contact_solutions),
    })
}

/// Stations to tolerance for a G1 rational-row strip (PK-734): the
/// strip interpolates its exact station rows, so its rails sit off the
/// supports between stations — a cylinder tee's saddle at the fixed 20
/// stations left them 1.3e-3 off. Walk again with the nested midpoints
/// (2n − 1 stations); while the midpoint contacts stray from a
/// quadratic through the coarser stations by more than 10 × the linear
/// tolerance, adopt the finer walk and repeat (to 1281 stations). The
/// proxy has the strip's own rail degree (its rows loft at degree 2):
/// a cubic proxy read 4.6e-6 where the strip's rails parted 1.3e-5
/// from the supports, and the offset tee's support seams read G0
/// (PK-894). A
/// walk already accurate is returned unchanged. G2 sweeps are not
/// refined here: their cubic loft is a different construction, and
/// denser G2 strips broke a later boolean on a filleted cylinder
/// (measured, see the PK-734 notes).
fn refine_station_density(
    request: &WalkRequest<'_>,
    walk_opts: WalkOptions,
    mut walk: crate::walker::WalkResult,
    mut n: usize,
    tolerance: &Tolerance,
) -> crate::walker::WalkResult {
    let target = tolerance.linear().max(1e-12) * 10.0;
    while n < 1281 {
        let finer_n = 2 * n - 1;
        let mut finer_opts = walk_opts;
        finer_opts.mode = crate::walker::WalkMode::Fixed { n_samples: finer_n };
        let Ok(finer) = walk_spine(request, &finer_opts, tolerance) else {
            break;
        };
        if finer.samples.len() != finer_n || walk.samples.len() != n {
            break;
        }
        let rail = |pick: fn(&crate::walker::WalkSample) -> Point3| {
            let pts: Vec<Point3> = walk.samples.iter().map(pick).collect();
            let mut params = Vec::with_capacity(pts.len());
            let mut acc = 0.0;
            params.push(acc);
            for w in pts.windows(2) {
                acc += w[0].distance_to(&w[1]);
                params.push(acc);
            }
            if params.windows(2).any(|w| w[1] <= w[0]) {
                return None;
            }
            popkern_core::geometry::fitting::interpolate_curve(&pts, &params, 2.min(pts.len() - 1))
                .ok()
        };
        let (Some(r0), Some(r1)) = (rail(|s| s.solution.contact0), rail(|s| s.solution.contact1))
        else {
            break;
        };
        let off = |r: &popkern_core::geometry::NurbsCurve, p: Point3| {
            r.closest_point(&p)
                .map(|s| r.point_at(s).distance_to(&p))
                .unwrap_or(f64::INFINITY)
        };
        let worst = finer
            .samples
            .iter()
            .skip(1)
            .step_by(2)
            .map(|s| off(&r0, s.solution.contact0).max(off(&r1, s.solution.contact1)))
            .fold(0.0, f64::max);
        if popkern_core::trace::trace_enabled("fillet") {
            eprintln!(
                "  [stations] walker n={n} worst rail midpoint {worst:.3e} target {target:.1e}"
            );
        }
        if worst <= target {
            break;
        }
        walk = finer;
        n = finer_n;
    }
    walk
}
