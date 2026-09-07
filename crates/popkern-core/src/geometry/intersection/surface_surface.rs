//! Surface-surface intersection types and API.

#[path = "surface_surface_tangential.rs"]
mod surface_surface_tangential;
use surface_surface_tangential::tangential_crossing_certificate;

#[path = "surface_surface_fitting.rs"]
mod surface_surface_fitting;
pub(crate) use surface_surface_fitting::{fit_march_samples, fit_marched_trace};
#[path = "surface_surface_torus_seed.rs"]
mod surface_surface_torus_seed;

#[path = "surface_surface_trim.rs"]
mod surface_surface_trim;
use surface_surface_trim::trim_to_windows;

use crate::error::CoreError;
use crate::geometry::Point3;
use crate::geometry::intersection::coincident::{
    closest_point_distance_pub, detect_coincident_regions,
};
use crate::geometry::intersection::coincident_analytic::are_surfaces_coincident_analytic;
use crate::geometry::intersection::degeneracy::{DegeneracyHandler, RobustDegeneracyHandler};
use crate::geometry::intersection::marching::{MarchSample, march};
use crate::geometry::intersection::plane_cylinder;
use crate::geometry::intersection::plane_plane::try_plane_plane;
use crate::geometry::intersection::plane_sphere;
use crate::geometry::intersection::starting_points::{find_branch_points, find_starting_points};
use crate::geometry::intersection::surface_surface_consumed::point_near_polylines;
use crate::geometry::intersection::{
    cylinder_cone, cylinder_cylinder, cylinder_sphere, plane_cone, plane_torus, torus_coaxial,
};
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::tolerance::Tolerance;

#[path = "surface_surface_types.rs"]
mod surface_surface_types;
pub use surface_surface_types::{SurfaceSurfaceIntersection, SurfaceSurfaceIntersectionResult};

/// Intersect two surfaces with default degeneracy handling.
pub fn intersect_surfaces(
    surface1: &SurfaceGeometry,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    surface2: &SurfaceGeometry,
    u2_range: (f64, f64),
    v2_range: (f64, f64),
    tolerance: &Tolerance,
) -> Result<SurfaceSurfaceIntersectionResult, CoreError> {
    intersect_surfaces_with_handler(
        surface1,
        u1_range,
        v1_range,
        surface2,
        u2_range,
        v2_range,
        tolerance,
        &RobustDegeneracyHandler,
    )
}

/// Intersect two surfaces with a custom degeneracy handler.
///
/// The `ssi` trace topic logs one line per call: the surface types,
/// windows, and what came back — the first stop when a boolean's
/// records go missing (an analytic arm, the coincident detector, and
/// the marching path all return through here).
#[allow(clippy::too_many_arguments)]
pub fn intersect_surfaces_with_handler(
    surface1: &SurfaceGeometry,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    surface2: &SurfaceGeometry,
    u2_range: (f64, f64),
    v2_range: (f64, f64),
    tolerance: &Tolerance,
    handler: &dyn DegeneracyHandler,
) -> Result<SurfaceSurfaceIntersectionResult, CoreError> {
    let result = intersect_surfaces_inner(
        surface1, u1_range, v1_range, surface2, u2_range, v2_range, tolerance, handler,
    );
    if crate::trace::trace_enabled("ssi") {
        let kind = |s: &SurfaceGeometry| -> &'static str {
            match s {
                SurfaceGeometry::Plane { .. } => "Plane",
                SurfaceGeometry::Cylinder { .. } => "Cylinder",
                SurfaceGeometry::Sphere { .. } => "Sphere",
                SurfaceGeometry::Cone { .. } => "Cone",
                SurfaceGeometry::Torus { .. } => "Torus",
                SurfaceGeometry::Nurbs(_) => "Nurbs",
                _ => "Other",
            }
        };
        match &result {
            Ok(r) => {
                let curves = r
                    .intersections
                    .iter()
                    .filter(|i| matches!(i, SurfaceSurfaceIntersection::Curve(_)))
                    .count();
                let regions = r.intersections.len() - curves;
                eprintln!(
                    "  [ssi] {}x{} u1[{:.3},{:.3}] v1[{:.3},{:.3}] u2[{:.3},{:.3}] v2[{:.3},{:.3}] -> {curves} curve(s), {regions} region(s)",
                    kind(surface1),
                    kind(surface2),
                    u1_range.0,
                    u1_range.1,
                    v1_range.0,
                    v1_range.1,
                    u2_range.0,
                    u2_range.1,
                    v2_range.0,
                    v2_range.1,
                );
            }
            Err(e) => eprintln!("  [ssi] {}x{} -> ERROR {e}", kind(surface1), kind(surface2)),
        }
    }
    result
}

#[allow(clippy::too_many_arguments)]
fn intersect_surfaces_inner(
    surface1: &SurfaceGeometry,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    surface2: &SurfaceGeometry,
    u2_range: (f64, f64),
    v2_range: (f64, f64),
    tolerance: &Tolerance,
    handler: &dyn DegeneracyHandler,
) -> Result<SurfaceSurfaceIntersectionResult, CoreError> {
    // Fast path: both surfaces are planes → analytic solver.
    if let (
        SurfaceGeometry::Plane {
            origin: o1,
            normal: n1,
            u_axis: ua1,
            v_axis: va1,
        },
        SurfaceGeometry::Plane {
            origin: o2,
            normal: n2,
            u_axis: ua2,
            v_axis: va2,
        },
    ) = (surface1, surface2)
    {
        let intersections = try_plane_plane(
            *o1, *n1, *ua1, *va1, u1_range, v1_range, *o2, *n2, *ua2, *va2, u2_range, v2_range,
            tolerance,
        )?;
        return Ok(SurfaceSurfaceIntersectionResult { intersections });
    }

    // Fast path: plane-sphere pair → analytic solver.
    if let Some(intersections) = plane_sphere::try_plane_sphere(
        surface1, u1_range, v1_range, surface2, u2_range, v2_range, tolerance,
    ) {
        return Ok(SurfaceSurfaceIntersectionResult { intersections });
    }

    // Fast path: plane-cylinder pair → analytic solver.
    if let Some(intersections) = plane_cylinder::try_plane_cylinder(
        surface1, u1_range, v1_range, surface2, u2_range, v2_range, tolerance,
    ) {
        return Ok(SurfaceSurfaceIntersectionResult { intersections });
    }

    // Fast path: plane-cone pair → analytic solver.
    if let Some(intersections) = plane_cone::try_plane_cone(
        surface1, u1_range, v1_range, surface2, u2_range, v2_range, tolerance,
    ) {
        return Ok(SurfaceSurfaceIntersectionResult { intersections });
    }

    // Fast path: COAXIAL cylinder-cone pair → exact circle(s); the
    // marching path fragmented this everyday class (PK-147).
    if let Some(intersections) = cylinder_cone::try_cylinder_cone(
        surface1, u1_range, v1_range, surface2, u2_range, v2_range, tolerance,
    ) {
        return Ok(SurfaceSurfaceIntersectionResult { intersections });
    }

    // Fast path: COAXIAL cylinder-sphere pair → exact circles, or an
    // exact empty answer for the tangent / enclosed regimes. The
    // marching path re-marched a fillet strip × corner sphere tangent
    // circle from every one of ~2100 tangent seeds (PK-390).
    if let Some(intersections) =
        cylinder_sphere::try_cylinder_sphere(surface1, v1_range, surface2, v2_range, tolerance)
    {
        return Ok(SurfaceSurfaceIntersectionResult { intersections });
    }

    // Fast path: cylinder-cylinder with non-parallel axes → exact
    // ruling-swept branches. The marching path fragments the
    // cross-drill class into overlapping partial traces and loses
    // whole branches (PK-170).
    if let Some(intersections) =
        cylinder_cylinder::try_cylinder_cylinder(surface1, v1_range, surface2, v2_range, tolerance)
    {
        return Ok(SurfaceSurfaceIntersectionResult { intersections });
    }

    // PK-279 torus analytic arms: plane×torus in any orientation
    // (exact circles / sampled-analytic spiric sections) and coaxial
    // {cylinder, cone, sphere, torus}×torus (exact circles). A `None`
    // (non-torus pair, skew/offset axes, numeric-guard decline) keeps
    // the marching fallback; coincident coaxial tori fall through to
    // the coincident detector. POPKERN_ANALYTIC_TORUS_SSI=0 is the
    // diagnostic kill-switch (PK-239 precedent) pinning the marching
    // path for A/B measurement and containment.
    // Deliberately a per-call env read (NOT OnceLock'd like the trace
    // gates): the torus-SSI acceptance ratio test toggles this
    // kill-switch mid-process to A/B the marching path, and a cached
    // read froze its second measurement on the analytic arm. One read
    // per intersect call is negligible next to the solve itself.
    let torus_arms_enabled = std::env::var_os("POPKERN_ANALYTIC_TORUS_SSI")
        .map(|v| v != "0")
        .unwrap_or(true);
    if torus_arms_enabled {
        if let Some(intersections) = plane_torus::try_plane_torus(
            surface1, u1_range, v1_range, surface2, u2_range, v2_range, tolerance,
        ) {
            return Ok(SurfaceSurfaceIntersectionResult { intersections });
        }
        if let Some(intersections) = torus_coaxial::try_torus_coaxial(
            surface1, u1_range, v1_range, surface2, u2_range, v2_range, tolerance,
        ) {
            return Ok(SurfaceSurfaceIntersectionResult { intersections });
        }
    }

    // Fast path: coincident surfaces → return CoincidentRegion immediately.
    if let Some(result) = try_coincident_fast(
        surface1, u1_range, v1_range, surface2, u2_range, v2_range, tolerance,
    ) {
        return Ok(result);
    }

    let mut results = Vec::new();

    // 1. Find starting points
    let mut starting_points = find_starting_points(
        surface1, u1_range, v1_range, surface2, u2_range, v2_range, tolerance,
    );
    // CLOSEST-APPROACH seeding fallback for torus pairs (PK-310) —
    // mechanics and measured rationale in
    // `surface_surface_torus_seed.rs`.
    surface_surface_torus_seed::push_torus_closest_approach_seeds(
        &mut starting_points,
        surface1,
        surface2,
        tolerance,
    );

    // 2. Handle tangent starting points -> coincident region detection
    let has_tangent = starting_points.iter().any(|sp| sp.is_tangent);
    let coincident_regions: Vec<SurfaceSurfaceIntersection> = if has_tangent {
        detect_coincident_regions(
            surface1, u1_range, v1_range, surface2, u2_range, v2_range, tolerance,
        )
    } else {
        Vec::new()
    };
    results.extend(coincident_regions.iter().cloned());

    // A tangent-classified seed is only accounted for when a detected
    // coincident region covers it. The tangent gate (|cos| within the
    // angular tolerance of 1) sits far above the march's own tangency
    // epsilon, so grazing CROSSINGS — fillet strip×cap contact — land
    // in this band; dropping them emitted nothing for a pair with a
    // real intersection curve (measured on the served pyramid4 body:
    // free-form-first seeding classified its only hit tangent and the
    // crossing vanished). An uncovered tangent seed is demoted to a
    // march seed when the normals still measurably diverge — sin of
    // the normal angle above the angular tolerance — because there the
    // march direction n1×n2 is numerically meaningful and the trace
    // advances. Seeds whose normals are angularly indistinguishable
    // (e.g. two sub-patches of one smooth surface touching at a split
    // corner, sin ~ 1e-9) stay dropped: their march direction is pure
    // noise and the robust handler creeps at minimum step without an
    // intersection curve to trace. Curvature comparison was tried and
    // falsified as this discriminator: at two approximately-matched
    // contact points the second-derivative readings differ above any
    // fixed gate even for genuine touches.
    let region_covers = |sample: &MarchSample| {
        coincident_regions.iter().any(|region| match region {
            SurfaceSurfaceIntersection::CoincidentRegion {
                u1_range: ur1,
                v1_range: vr1,
                u2_range: ur2,
                v2_range: vr2,
            } => {
                let within = |x: f64, r: &(f64, f64)| x >= r.0.min(r.1) && x <= r.0.max(r.1);
                within(sample.u1, ur1)
                    && within(sample.v1, vr1)
                    && within(sample.u2, ur2)
                    && within(sample.v2, vr2)
            }
            SurfaceSurfaceIntersection::Curve(_) => false,
        })
    };
    let tangent_seed_marchable = |sample: &MarchSample| {
        let n1 = surface1.normal_at(sample.u1, sample.v1);
        let n2 = surface2.normal_at(sample.u2, sample.v2);
        let norms = n1.length() * n2.length();
        if norms < tolerance.numerical_floor() {
            return false;
        }
        let sin_angle = n1.cross(&n2).length() / norms;
        sin_angle > tolerance.angular()
    };

    // 3. March from non-tangent, unconsumed starting points.
    //    Branch detection may add extra seeds; we process them all via a queue.
    const MAX_CURVES: usize = 20;

    // Traced polylines. Seed consumption tests distance to the polyline
    // SEGMENTS, not to the discrete sample points: the adaptive march
    // steps up to 1000× linear tolerance between samples while its
    // curvature-based cap keeps each chord within ~linear tolerance of
    // the true curve, so a converged seed landing mid-segment lies on
    // the traced curve but far from every sample point. The former
    // point-proximity check re-marched the identical curve once per
    // such seed (measured: 18–20 copies of one cylinder×NURBS arc,
    // count varying with operand order).
    let mut consumed: Vec<Vec<Point3>> = Vec::new();
    let consumed_radius = tolerance.linear() * 10.0;

    // Seed the queue with every transversal starting point — plus,
    // ONLY when the pair produced no transversal seed at all, every
    // tangent-classified seed that no coincident region accounts for
    // and whose normals measurably diverge. Demoted seeds carry a
    // marker: their traces must additionally present a signed-distance
    // crossing certificate before emission (below). Both scopes exist
    // because unconditional demotion flooded callers built on the
    // tangency-silent contract — the N>=5 setback pipeline's dispatch
    // bailed on shared-boundary tangency curves, and the Gregory
    // subdivision self-intersection check declined every KP cap when
    // its diagonal sub-patches' corner-touch blob fitted into a
    // "curve" — while the grazing-only class (the sole evidence of a
    // real crossing is a tangent-classified seed) is exactly where
    // demotion is the difference between finding the curve and
    // silence.
    let no_transversal_seeds = starting_points.iter().all(|sp| sp.is_tangent);
    let mut branch_queue: Vec<(MarchSample, bool)> = starting_points
        .iter()
        .filter_map(|sp| {
            if !sp.is_tangent {
                Some((sp.sample.clone(), false))
            } else if no_transversal_seeds
                && !region_covers(&sp.sample)
                && tangent_seed_marchable(&sp.sample)
            {
                Some((sp.sample.clone(), true))
            } else {
                None
            }
        })
        .collect();

    // March windows: pad non-periodic ranges by 1% of extent. A seed
    // converges from a patch midpoint and may land slightly OUTSIDE
    // the query window; in a periodic direction the march wraps and
    // recovers, but in a non-periodic direction the domain-exit check
    // killed the trace at step one and the curve vanished (measured on
    // the served pyramid4 strip×cap pair: the free-form-first seed sat
    // 2.5e-3 outside the cap's knot window). The pad only widens march
    // termination; seeding, coincidence detection, and reported ranges
    // keep the caller's windows, and trim classification downstream
    // clips curve extent.
    let pad_window = |range: (f64, f64), periodic: bool| {
        if periodic {
            range
        } else {
            let entry_band = 0.01 * (range.1 - range.0);
            (range.0 - entry_band, range.1 + entry_band)
        }
    };
    // Periodicity here includes geometrically-closed NURBS charts —
    // the march wraps their seam instead of terminating on it (see
    // `geometry::chart_periods`).
    let (periodic_u1, periodic_v1) = crate::geometry::chart_periods::marching_periodicity(surface1);
    let (periodic_u2, periodic_v2) = crate::geometry::chart_periods::marching_periodicity(surface2);
    let march_u1 = pad_window(u1_range, periodic_u1);
    let march_v1 = pad_window(v1_range, periodic_v1);
    let march_u2 = pad_window(u2_range, periodic_u2);
    let march_v2 = pad_window(v2_range, periodic_v2);

    while let Some((seed, demoted)) = branch_queue.pop() {
        // Cap total curves to avoid runaway branching.
        if results.len() >= MAX_CURVES {
            break;
        }

        // Skip if already traced through this region.
        if point_near_polylines(&seed.point, &consumed, consumed_radius) {
            continue;
        }

        let (samples, is_closed) = march(
            surface1, march_u1, march_v1, surface2, march_u2, march_v2, &seed, tolerance, handler,
        );

        // For EMITTING traces, consumption covers the FULL trace, not
        // the emitted span: seeds converge unclamped and legitimately
        // land on curve stretches beyond the caller's windows (the
        // march windows are padded), and consuming only the trimmed
        // polyline re-marched the identical curve once per
        // out-of-window seed — a thread-groove flank × shank pair
        // re-traced its 28-unit helix 19 times and flooded
        // consolidation with duplicate records (PK-460). Scoped to
        // emitting traces only: a trace TRIMMED TO NOTHING or a
        // failed-certificate blob keeps its tuned consumption (the
        // grazing strip×cap pin measured that full consumption there
        // starves the honest re-seeds — tangent-chart uv registration
        // is sloppy enough that a blob's samples can all read
        // out-of-window while later seeds trace the real crossing).
        let full_trace: Vec<Point3> = samples.iter().map(|s| s.point).collect();

        // The entry band exists so an out-of-window seed can walk into
        // the caller's windows — it must not widen the EMITTED curve:
        // callers key curve endpoints against their own trim windows
        // (measured: 1% overhang declined the variable-radius setback
        // serve). Trim leading/trailing samples outside the caller's
        // windows, refining each trimmed end onto the exact window
        // boundary — mechanism and measured rationale in
        // `surface_surface_trim`.
        let samples = trim_to_windows(
            &samples, surface1, u1_range, v1_range, surface2, u2_range, v2_range, tolerance,
        );

        if samples.len() < 2 {
            continue; // too few samples to fit
        }
        // A demoted-tangent trace must certify an actual CROSSING
        // before emission: probes offset to either side of the traced
        // curve must flip the sign of the analytic operand's implicit
        // distance. Near-contact blobs (surfaces within tolerance
        // around a touch, no sign change) and pairs with no analytic
        // side to certify against stay on the tangency-silent
        // contract.
        if demoted && !tangential_crossing_certificate(surface1, surface2, &samples, tolerance) {
            // The trace is a near-contact blob, not a crossing — but it
            // IS the answer for every other seed lying on it. Mark it
            // consumed so the remaining seeds of a tangent band do not
            // each repeat the identical full-cap march (measured: ~2100
            // seeds × 20 001 steps on a fillet strip × corner sphere
            // pair before the analytic arm covered that class, PK-390).
            // Trimmed span only — see the scoping note on
            // `full_trace` above.
            consumed.push(samples.iter().map(|s| s.point).collect());
            continue;
        }

        // Mark the traced polyline as consumed (full trace — see above).
        consumed.push(full_trace);

        // 4. Segment at UV period jumps, weld marched closure, fit
        //    (PK-304) — mechanism and measured rationale in
        //    `surface_surface_fitting::fit_marched_trace`.
        fit_marched_trace(
            &samples,
            is_closed,
            surface1,
            surface2,
            tolerance,
            &mut results,
        );

        // 5. Branch detection: look for additional intersection
        //    branches near this seed.
        if results.len() < MAX_CURVES {
            let new_branches = find_branch_points(
                surface1, u1_range, v1_range, surface2, u2_range, v2_range, &samples, &seed,
                tolerance,
            );
            for branch in new_branches {
                if !point_near_polylines(&branch.point, &consumed, consumed_radius) {
                    branch_queue.push((branch, demoted));
                }
            }
        }
    }
    Ok(SurfaceSurfaceIntersectionResult {
        intersections: results,
    })
}

/// Attempt to detect coincident surfaces before the expensive marching solver.
///
/// **Tier A** — analytic O(1) comparison for same-type analytic pairs.
/// **Tier B** — sample-based: 4-corner quick rejection followed by the
/// existing 8×8 grid check for NURBS / mixed-type pairs.
#[allow(clippy::too_many_arguments)]
fn try_coincident_fast(
    surface1: &SurfaceGeometry,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    surface2: &SurfaceGeometry,
    u2_range: (f64, f64),
    v2_range: (f64, f64),
    tolerance: &Tolerance,
) -> Option<SurfaceSurfaceIntersectionResult> {
    // Tier A: analytic comparison for same-type analytic surfaces.
    if are_surfaces_coincident_analytic(surface1, surface2, tolerance) {
        return Some(SurfaceSurfaceIntersectionResult {
            intersections: vec![SurfaceSurfaceIntersection::CoincidentRegion {
                u1_range,
                v1_range,
                u2_range,
                v2_range,
            }],
        });
    }

    // Tier B: sample-based detection for NURBS / mixed types.
    // Quick rejection: sample the 4 corners of surface1's parameter domain and
    // project each onto surface2. If any corner is far away, bail out early.
    let corners = [
        (u1_range.0, v1_range.0),
        (u1_range.0, v1_range.1),
        (u1_range.1, v1_range.0),
        (u1_range.1, v1_range.1),
    ];
    let rejection_threshold = tolerance.linear() * 10.0;
    for (u, v) in &corners {
        let p = surface1.point_at(*u, *v);
        let dist = closest_point_distance_pub(p, surface2, u2_range, v2_range);
        if dist > rejection_threshold {
            return None;
        }
    }

    // All corners are close — run the full 8×8 grid check.
    let regions = detect_coincident_regions(
        surface1, u1_range, v1_range, surface2, u2_range, v2_range, tolerance,
    );
    if regions.is_empty() {
        None
    } else {
        Some(SurfaceSurfaceIntersectionResult {
            intersections: regions,
        })
    }
}

#[cfg(test)]
#[path = "surface_surface_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "surface_surface_degenerate_tests.rs"]
mod degenerate_tests;
