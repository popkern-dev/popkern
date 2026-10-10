//! Walker-based center-curve routing. Pulls the edge curve from the
//! BRep, builds a `WalkRequest` against it, and routes through
//! `walker_based_center_curve` which handles the full radius-type
//! matrix (Constant, Variable, ByPoints).
//!
//! **When this path runs:**
//! - Walker trigger: either support surface is a Cone (offset-intersect
//!   path can't run because cones have no analytic offset).
//! - The ComputeGeometry G2 step: the caller requested G2 cross-sections via
//!   `FilletOptions::g2_cross_section`. The walker populates
//!   `CenterCurveResult::contact_solutions`, which lets the downstream
//!   `build_fillet_surface` route through `g2_sweep_from_contact_solutions`
//!   on any surface pair (not just cones).
//!
//! Split out of `geometry_compute.rs` to keep the parent under the
//! 500-line structural cap.

use popkern_core::geometry::Point3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::analysis;
use crate::center_curve;
use crate::center_curve::CenterCurveResult;
use crate::error::BlendFailure;
use crate::types::FilletRadius;

#[allow(clippy::too_many_arguments)]
pub(crate) fn compute_center_via_walker(
    analysis: &analysis::EdgeAnalysis,
    radius: &FilletRadius,
    brep: &BRep,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    u2_range: (f64, f64),
    v2_range: (f64, f64),
    tolerance: &Tolerance,
    edge_midpoint: Point3,
    extend_ends: bool,
) -> Result<center_curve::CenterCurveResult, BlendFailure> {
    let edge = brep
        .edges
        .get(analysis.edge)
        .ok_or_else(|| BlendFailure::InvalidTopology {
            reason: format!("edge {:?} not found", analysis.edge),
            location: None,
        })?;
    let curve_on_edge = edge
        .curve
        .as_ref()
        .ok_or_else(|| BlendFailure::InvalidTopology {
            reason: "cone fillet requires the edge to carry curve geometry".to_string(),
            location: Some(edge_midpoint),
        })?;

    // Map the radius profile to a closure over the [0, 1] parameter
    // that `walker_based_center_curve` expects.
    let radius_at_normalized: Box<dyn Fn(f64) -> f64> = match radius {
        FilletRadius::Constant(r) => {
            let r = *r;
            Box::new(move |_t| r)
        }
        FilletRadius::Variable { start, end } => {
            let start = *start;
            let end = *end;
            Box::new(move |t| start + (end - start) * t)
        }
        FilletRadius::ByPoints(pts) => {
            if pts.is_empty() {
                return Err(BlendFailure::InvalidTopology {
                    reason: "ByPoints radius with empty point list".to_string(),
                    location: Some(edge_midpoint),
                });
            }
            let pts = pts.clone();
            Box::new(move |t| center_curve::bypoints_local_radius(&pts, t))
        }
    };

    // Project a query point onto both surfaces to get an initial
    // UV seed, nudged toward each surface normal to break atan2
    // degeneracies on periodic surfaces.
    //
    // The seed point is the spine START when it lies on both
    // supports: the walk's first station is there, and the contact
    // solver's station constraint (the plane through the spine point
    // normal to its tangent) does not by itself pick the ball at that
    // point — on a closed circular rim the antipode satisfies the
    // same plane, and a midpoint seed converged to THAT ball. The
    // strip then originated half a turn from the rail vertex where
    // the closed-fillet surgery seams the ring (measured: the boss
    // plate's union rim as an exact arc, 70 open mesh edges; the
    // NURBS spelling only escaped by evaluating a few ULPs off the
    // supports, which made an endpoint "closer" than the midpoint —
    // PK-419).
    //
    // For blend-on-blend (§12.3) an endpoint may NOT lie on both
    // surfaces (edge from a cube corner to a cylinder rim): then the
    // other endpoint, else the midpoint, whichever is nearest both.
    let (initial_seed, spine_t_range) = {
        use popkern_core::geometry::SurfaceEval;
        use popkern_core::geometry::curve::CurveEval;
        let ts = curve_on_edge.t_start;
        let te = curve_on_edge.t_end;
        let t_mid = 0.5 * (ts + te);
        let p_mid = curve_on_edge.geometry.point_at(t_mid);
        let p_start = curve_on_edge.geometry.point_at(ts);
        let p_end = curve_on_edge.geometry.point_at(te);

        let dist_to_surfs = |p: &Point3| -> f64 {
            let d = |s: &popkern_core::geometry::surface::SurfaceGeometry| {
                s.closest_point(p)
                    .map(|(u, v)| (s.point_at(u, v) - *p).length())
                    .unwrap_or(f64::INFINITY)
            };
            d(&analysis.surface_1).max(d(&analysis.surface_2))
        };
        let on_supports = tolerance.linear() * 10.0;
        // A TANGENT contact (the supports' normals parallel: the crease
        // vanishes, as at the singular points of an equal-radius tee)
        // is no seed. The ball there touches both supports at the
        // spine point itself, and the next station warm-started from
        // that solution met a singular contact Jacobian before its
        // first iteration (PK-764). The walk starts at the regular end
        // and reaches the tangent end last, converging onto it.
        let tangent_contact = |p: &Point3| -> bool {
            let normal = |s: &popkern_core::geometry::surface::SurfaceGeometry| {
                s.closest_point(p)
                    .map(|(u, v)| s.normal_at(u, v))
                    .and_then(|n| n.normalized())
            };
            match (normal(&analysis.surface_1), normal(&analysis.surface_2)) {
                (Some(a), Some(b)) => a.cross(&b).length() < 1e-6,
                _ => false,
            }
        };
        let best_pt = [p_start, p_end, p_mid]
            .into_iter()
            .find(|p| dist_to_surfs(p) <= on_supports && !tangent_contact(p))
            .unwrap_or_else(|| {
                let mut best = p_mid;
                for candidate in [p_start, p_end] {
                    if dist_to_surfs(&candidate) < dist_to_surfs(&best) {
                        best = candidate;
                    }
                }
                best
            });

        let nudge = 1.0e-4;
        let proj = |s: &popkern_core::geometry::surface::SurfaceGeometry,
                    ur: (f64, f64),
                    vr: (f64, f64)|
         -> Option<(f64, f64)> {
            // Nudge along the normal AT THE POINT'S OWN chart image.
            // The old chart-CENTER normal points anywhere on periodic
            // surfaces — for a boss edge near the cylinder seam the
            // face's u-range midpoint is the OPPOSITE side of the
            // drum, the nudge inverted, and the contact Newton
            // diverged from the mirrored seed on every concave wall
            // config (PK-158: residual 3e-1, exact-pair fallback).
            let (u0, v0) = s.closest_point(&best_pt)?;
            let n = s.normal_at(u0, v0);
            let (u, v) = s.closest_point(&(best_pt + n * nudge))?;
            // A periodic coordinate is carried into the face's window
            // by whole periods BEFORE the clamp: `closest_point`
            // answers in the principal branch (u = −0.25 for a station
            // just behind the seam of a [0, 2π] chart), and clamping
            // that branch parked the seed ON the seam, a quarter
            // radian from the station — the contact Newton diverged
            // from there (measured: the seam-crossing arc piece of a
            // boss-on-wall chain, residual 2.6e-1; the same station
            // seeded from the other side converged at once).
            let (period_u, period_v) = s.parameter_periodicity();
            Some((
                crate::center_curve_projection::wrap_then_clamp(u, ur, period_u),
                crate::center_curve_projection::wrap_then_clamp(v, vr, period_v),
            ))
        };
        let uv0 = proj(&analysis.surface_1, u1_range, v1_range);
        let uv1 = proj(&analysis.surface_2, u2_range, v2_range);
        let seed = match (uv0, uv1) {
            (Some(a), Some(b)) => Some((a, b)),
            _ => None,
        };
        // Walk from the end nearest to the seed point so the
        // initial contact has the best chance of converging.
        let ds = (best_pt - p_start).length();
        let de = (best_pt - p_end).length();
        let t_range = if de < ds { (te, ts) } else { (ts, te) };
        (seed, t_range)
    };
    // PK-531: a fitted (G2) strip must reach past the edge's end
    // stations. An oblique end face leans past the station, and a
    // blend-on-blend end meets the neighbour strips above it. A strip
    // swept only over the edge's own range stops short: the end
    // section lay beyond the strip (2.5e-2 off its clamped border),
    // and the rail extrapolated from the sampled trims landed 1.26e-2
    // off. Straight edges between planar supports extend the walk by
    // up to 4 r past each end (at most half the edge) over padded
    // plane ranges. The face's own loop still bounds the strip.
    let extension = extend_ends
        .then(|| {
            let radius = match radius {
                FilletRadius::Constant(r) => *r,
                _ => return None,
            };
            let is_line = matches!(
                curve_on_edge.geometry,
                popkern_core::geometry::CurveGeometry::Line { .. }
            );
            let planar = matches!(
                (&analysis.surface_1, &analysis.surface_2),
                (
                    popkern_core::geometry::surface::SurfaceGeometry::Plane { .. },
                    popkern_core::geometry::surface::SurfaceGeometry::Plane { .. },
                )
            );
            use popkern_core::geometry::curve::CurveEval;
            let length = (curve_on_edge.geometry.point_at(curve_on_edge.t_end)
                - curve_on_edge.geometry.point_at(curve_on_edge.t_start))
            .length();
            let reach = (4.0 * radius).min(0.5 * length);
            (is_line && planar && length > 0.0)
                .then(|| (reach / length, 2.0 * reach + 4.0 * radius))
        })
        .flatten();
    // The walk's station count scales with the extended span, so the
    // section density over the edge itself is unchanged: the fixed 20
    // stations spread over 1.8x the span measured a 20x larger seam gap
    // on the rim+verticals box (5e-11 -> 1.02e-9).
    const EDGE_STATIONS: f64 = 20.0;
    let (spine_t_range, uv_pad, stations) = match extension {
        Some((fraction, pad)) => {
            // Plane charts are metric, so the pad is a length: the
            // extension twice over plus the contact offset.
            let delta = fraction * (spine_t_range.1 - spine_t_range.0);
            (
                (spine_t_range.0 - delta, spine_t_range.1 + delta),
                pad,
                (EDGE_STATIONS * (1.0 + 2.0 * fraction)).ceil() as usize,
            )
        }
        None => (spine_t_range, 0.0, EDGE_STATIONS as usize),
    };
    let widen = |r: (f64, f64)| (r.0 - uv_pad, r.1 + uv_pad);
    let (u1_range, v1_range, u2_range, v2_range) = (
        widen(u1_range),
        widen(v1_range),
        widen(u2_range),
        widen(v2_range),
    );

    center_curve::walker_based_center_curve(
        &analysis.surface_1,
        u1_range,
        v1_range,
        &analysis.surface_2,
        u2_range,
        v2_range,
        radius_at_normalized,
        &curve_on_edge.geometry,
        spine_t_range,
        tolerance,
        analysis.orientation_1 == analysis.shell_outward,
        analysis.orientation_2 == analysis.shell_outward,
        analysis.convexity,
        initial_seed,
        stations,
        // Station refinement (PK-734), closed spines included: a closed
        // rim's refined strip once broke a later subtraction ('Loop not
        // geometrically closed', gap 0.87). Its fitted rim circle carried
        // centre roundoff, so boundary points sat either side of the
        // torus's v = 0 seam; the contact circles are now exact
        // (`trim::exact_torus_contact_circles`, PK-738).
        !extend_ends,
    )
    .map_err(|e| BlendFailure::IntersectionFailed {
        reason: e,
        location: Some(edge_midpoint),
    })
}

/// Walker-or-offset-intersect center curve dispatch with §7 continuity
/// downgrade tracking. Returns `(center_result, downgrade_reason)` —
/// the reason string is `Some` when the walker failed on a non-cone
/// G2 request and the pipeline degraded to the vertex-blend path.
#[allow(clippy::too_many_arguments)]
pub(crate) fn compute_center_with_fallback(
    analysis: &analysis::EdgeAnalysis,
    radius: &FilletRadius,
    brep: &BRep,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    u2_range: (f64, f64),
    v2_range: (f64, f64),
    tolerance: &Tolerance,
    edge_midpoint: Point3,
    either_needs_walker: bool,
    wants_g2: bool,
) -> Result<(CenterCurveResult, Option<String>), BlendFailure> {
    let walker_wanted = (either_needs_walker || wants_g2)
        && crate::dispatch_overrides::WALKER_DISABLE.get() != Some(true);
    if walker_wanted {
        let w = compute_center_via_walker(
            analysis,
            radius,
            brep,
            u1_range,
            v1_range,
            u2_range,
            v2_range,
            tolerance,
            edge_midpoint,
            wants_g2,
        );
        // Walker-failure fallback boundary: the offset-
        // intersect path is EXACT for plane×cylinder pairs (offset
        // plane ∩ offset cylinder is an analytic line/circle), so a
        // diverged walker on that pair degrades to offset-intersect with the
        // reason recorded (measured: 2 of 4 blocked seam-crossing
        // boss-wall arc fillets complete this way). Cone/torus pairs
        // keep the hard error — analytic offsets are not trustworthy here
        // there (the walker route exists because of them), and a
        // plausible-but-wrong center curve would ship silent wrong
        // shape.
        use popkern_core::geometry::surface::SurfaceGeometry;
        let offset_intersect_exact_pair = matches!(
            (&analysis.surface_1, &analysis.surface_2),
            (
                SurfaceGeometry::Plane { .. },
                SurfaceGeometry::Cylinder { .. }
            ) | (
                SurfaceGeometry::Cylinder { .. },
                SurfaceGeometry::Plane { .. }
            ) | (SurfaceGeometry::Plane { .. }, SurfaceGeometry::Plane { .. })
        );
        match w {
            Ok(r) => Ok((r, None)),
            Err(walker_err)
                if offset_intersect_exact_pair || (!either_needs_walker && wants_g2) =>
            {
                if popkern_core::trace::trace_enabled("fillet") {
                    eprintln!("  [walker-fallback] {walker_err:?}");
                }
                let reason = format!("walker failed: {walker_err:?}");
                let offset_intersect = compute_offset_intersect_center_curve(
                    analysis,
                    radius,
                    brep,
                    u1_range,
                    v1_range,
                    u2_range,
                    v2_range,
                    tolerance,
                    edge_midpoint,
                )?;
                Ok((offset_intersect, Some(reason)))
            }
            Err(e) => Err(e),
        }
    } else {
        compute_offset_intersect_center_curve(
            analysis,
            radius,
            brep,
            u1_range,
            v1_range,
            u2_range,
            v2_range,
            tolerance,
            edge_midpoint,
        )
        .map(|r| (r, None))
    }
}

/// Offset-intersect center curve: the pre-walker path that
/// analytically offsets both surfaces, intersects the offsets, and
/// picks the curve closest to the edge midpoint. Used as the primary
/// route when neither the cone trigger nor the G2 opt-in fires, and
/// as a fallback when the walker fails on a G2 opt-in request
/// (silently degrading the caller's G2 request to G1 rather than
/// propagating a walker-specific error).
#[allow(clippy::too_many_arguments)]
pub(crate) fn compute_offset_intersect_center_curve(
    analysis: &analysis::EdgeAnalysis,
    radius: &FilletRadius,
    brep: &BRep,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    u2_range: (f64, f64),
    v2_range: (f64, f64),
    tolerance: &Tolerance,
    edge_midpoint: Point3,
) -> Result<CenterCurveResult, BlendFailure> {
    let wrap_err = |e: String| BlendFailure::IntersectionFailed {
        reason: e,
        location: Some(edge_midpoint),
    };
    // Per-station anchors along the sharp edge for the variable-law
    // center marches; midpoint fallback when the edge has no curve.
    let edge_curve = brep
        .edges
        .get(analysis.edge)
        .and_then(|edge| edge.curve.clone());
    let edge_anchor_at = move |t: f64| -> Point3 {
        match &edge_curve {
            Some(c) => {
                use popkern_core::geometry::curve::CurveEval;
                c.geometry.point_at(c.t_start + (c.t_end - c.t_start) * t)
            }
            None => edge_midpoint,
        }
    };
    match radius {
        FilletRadius::Constant(r) => center_curve::constant_radius_center_curve(
            &analysis.surface_1,
            u1_range,
            v1_range,
            &analysis.surface_2,
            u2_range,
            v2_range,
            *r,
            &edge_midpoint,
            &edge_anchor_at,
            tolerance,
            analysis.orientation_1 == analysis.shell_outward,
            analysis.orientation_2 == analysis.shell_outward,
            analysis.convexity,
        )
        .map_err(wrap_err),
        FilletRadius::Variable { start, end } => center_curve::variable_radius_center_curve(
            &analysis.surface_1,
            u1_range,
            v1_range,
            &analysis.surface_2,
            u2_range,
            v2_range,
            *start,
            *end,
            &edge_anchor_at,
            tolerance,
            analysis.orientation_1 == analysis.shell_outward,
            analysis.orientation_2 == analysis.shell_outward,
            analysis.convexity,
        )
        .map_err(wrap_err),
        FilletRadius::ByPoints(pts) => {
            if pts.is_empty() {
                return Err(BlendFailure::InvalidTopology {
                    reason: "ByPoints radius with empty point list".to_string(),
                    location: Some(edge_midpoint),
                });
            }
            center_curve::parametric_radius_center_curve(
                &analysis.surface_1,
                u1_range,
                v1_range,
                &analysis.surface_2,
                u2_range,
                v2_range,
                |t| center_curve::bypoints_local_radius(pts, t),
                &edge_anchor_at,
                tolerance,
                analysis.orientation_1 == analysis.shell_outward,
                analysis.orientation_2 == analysis.shell_outward,
                analysis.convexity,
            )
            .map_err(wrap_err)
        }
    }
}
