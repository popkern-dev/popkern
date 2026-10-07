//! Surface-surface intersection types and API.

#[path = "surface_surface_tangential.rs"]
mod surface_surface_tangential;
use surface_surface_tangential::tangential_crossing_certificate;

#[path = "surface_surface_fitting.rs"]
mod surface_surface_fitting;
pub(crate) use surface_surface_fitting::{fit_march_samples, fit_marched_trace};
#[path = "surface_surface_implicit_seed.rs"]
mod surface_surface_implicit_seed;
#[path = "surface_surface_torus_seed.rs"]
mod surface_surface_torus_seed;

#[path = "surface_surface_trim.rs"]
mod surface_surface_trim;
use surface_surface_trim::trim_to_windows;

use crate::error::CoreError;
use crate::geometry::Point3;
use crate::geometry::intersection::coincident::detect_coincident_regions;
use crate::geometry::intersection::degeneracy::{DegeneracyHandler, RobustDegeneracyHandler};
use crate::geometry::intersection::marching::{MarchSample, march};
use crate::geometry::intersection::starting_points::{find_branch_points, find_starting_points};
use crate::geometry::intersection::surface_surface_consumed::{
    ConsumedTrace, point_near_polylines,
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
        surface1,
        u1_range,
        v1_range,
        surface2,
        u2_range,
        v2_range,
        tolerance,
        handler,
        &[],
    );
    trace_ssi(
        &result, surface1, surface2, u1_range, v1_range, u2_range, v2_range,
    );
    result
}

/// Intersect two surfaces with caller-supplied junction seed points.
///
/// Each seed is a 3D point known to lie ON the intersection — typically
/// an edge-face pierce point from a 1D edge×surface solve, which stays
/// robust exactly where the 2D patch-midpoint seeding goes blind: a
/// narrow near-tangent lens. Measured on the angled tap (PK-465): the
/// thread-runout crest band crosses the tilted bore wall in
/// 0.05–0.15-long lens segments; uniform patch sampling produced zero
/// seeds and the whole face pair came back empty while every lens
/// endpoint was already known as an edge pierce. Seeds are projected
/// onto both surfaces, converged by the same Newton the interior
/// seeding uses, and marched like any transversal seed; a seed landing
/// on an already-traced curve is consumed by the standard segment
/// dedup, so over-seeding costs only the checks.
#[allow(clippy::too_many_arguments)]
pub fn intersect_surfaces_with_seeds(
    surface1: &SurfaceGeometry,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    surface2: &SurfaceGeometry,
    u2_range: (f64, f64),
    v2_range: (f64, f64),
    tolerance: &Tolerance,
    seed_points: &[Point3],
) -> Result<SurfaceSurfaceIntersectionResult, CoreError> {
    let result = intersect_surfaces_inner(
        surface1,
        u1_range,
        v1_range,
        surface2,
        u2_range,
        v2_range,
        tolerance,
        &RobustDegeneracyHandler,
        seed_points,
    );
    trace_ssi(
        &result, surface1, surface2, u1_range, v1_range, u2_range, v2_range,
    );
    result
}

#[path = "surface_surface_trace.rs"]
mod surface_surface_trace;
use surface_surface_trace::{trace_ssi, trace_ssi_seed};
#[path = "surface_surface_fastpath.rs"]
mod surface_surface_fastpath;

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
    seed_points: &[Point3],
) -> Result<SurfaceSurfaceIntersectionResult, CoreError> {
    // Analytic fast paths (exact solvers + coincidence detection);
    // `None` falls through to the marching solver.
    if let Some(result) = surface_surface_fastpath::try_analytic_pair(
        surface1, u1_range, v1_range, surface2, u2_range, v2_range, tolerance,
    )? {
        return Ok(result);
    }

    let mut results = Vec::new();

    // 1. Find starting points
    let mut starting_points = find_starting_points(
        surface1, u1_range, v1_range, surface2, u2_range, v2_range, tolerance,
    );
    // Seeds too are accepted unclamped and can carry a wrong-preimage
    // registration on self-overlapping extended bases — re-branch
    // them into the query windows (same rationale as the march's
    // `rebranch_out_of_window`). A seed the re-branch REJECTS is an
    // extrapolation phantom (far outside a window, no in-window
    // preimage) and is dropped (PK-513): every such march failed on
    // its first step after re-paying the rejected projections
    // (measured: ~18k phantom seeds at u2 in [-7, 16] for a [0, 1]
    // thread-wall window, 36k dead march directions, in one M10 stud
    // tip-chamfer build — the bulk of the cut's cost).
    let mut seed_phantom_cache = Vec::new();
    starting_points.retain_mut(|sp| {
        crate::geometry::intersection::marching::rebranch_sample_pub(
            surface1,
            u1_range,
            v1_range,
            surface2,
            u2_range,
            v2_range,
            &mut sp.sample,
            tolerance,
            &mut seed_phantom_cache,
        )
    });
    // CLOSEST-APPROACH seeding fallback for torus pairs (PK-310) —
    // mechanics and measured rationale in
    // `surface_surface_torus_seed.rs`.
    surface_surface_torus_seed::push_torus_closest_approach_seeds(
        &mut starting_points,
        surface1,
        surface2,
        tolerance,
    );

    // Caller-supplied junction seeds (see
    // [`intersect_surfaces_with_seeds`]): project each onto both
    // surfaces and converge with the same Newton the interior seeding
    // uses. A seed that fails to project or converge is dropped
    // silently — it costs the caller nothing, and the pierce points
    // callers supply can legitimately sit on face pairs that do not
    // intersect near them (an edge bounds two faces; only one may
    // carry the continuing curve).
    crate::geometry::intersection::starting_points::push_junction_seeds(
        &mut starting_points,
        surface1,
        u1_range,
        v1_range,
        surface2,
        u2_range,
        v2_range,
        tolerance,
        seed_points,
        &mut seed_phantom_cache,
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
    let mut consumed: Vec<ConsumedTrace> = Vec::new();
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

    // RESCUE loop: when ordinary seeding-and-marching emits nothing
    // (an extrapolation-phantom basin can swallow every
    // patch-midpoint solve — or leave only phantom seeds that
    // domain-kill at step one), seed once more from exact
    // implicit-distance sign changes on a knot-span grid — a
    // crossing certificate, not a guess (see
    // surface_surface_implicit_seed) — and march those.
    let mut rescue_tried = false;
    'rescue: loop {
        while let Some((seed, demoted)) = branch_queue.pop() {
            // Cap total curves to avoid runaway branching.
            if results.len() >= MAX_CURVES {
                break;
            }

            // Skip if already traced through this region.
            if point_near_polylines(&seed.point, &consumed, consumed_radius) {
                trace_ssi_seed(true, &seed.point, demoted, consumed_radius, consumed.len());
                continue;
            }
            trace_ssi_seed(false, &seed.point, demoted, consumed_radius, consumed.len());

            let (samples, is_closed) = march(
                surface1, march_u1, march_v1, surface2, march_u2, march_v2, &seed, tolerance,
                handler,
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
            if demoted && !tangential_crossing_certificate(surface1, surface2, &samples, tolerance)
            {
                // The trace is a near-contact blob, not a crossing — but it
                // IS the answer for every other seed lying on it. Mark it
                // consumed so the remaining seeds of a tangent band do not
                // each repeat the identical full-cap march (measured: ~2100
                // seeds × 20 001 steps on a fillet strip × corner sphere
                // pair before the analytic arm covered that class, PK-390).
                // Trimmed span only — see the scoping note on
                // `full_trace` above.
                consumed.push(ConsumedTrace::open(
                    samples.iter().map(|s| s.point).collect(),
                ));
                continue;
            }

            // Mark the traced polyline as consumed (full trace — see above).
            // A closed trace consumes as a RING: the closure gap the
            // step loop never walked is claimed at the gap's own
            // radius (PK-484 — see `ConsumedTrace`).
            if is_closed {
                consumed.push(ConsumedTrace::closed(full_trace));
            } else {
                consumed.push(ConsumedTrace::open(full_trace));
            }

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
        if results.is_empty() && !rescue_tried {
            rescue_tried = true;
            surface_surface_implicit_seed::push_implicit_grid_seeds(
                &mut branch_queue,
                surface1,
                u1_range,
                v1_range,
                surface2,
                u2_range,
                v2_range,
                tolerance,
            );
            if !branch_queue.is_empty() {
                continue 'rescue;
            }
        }
        break;
    }
    Ok(SurfaceSurfaceIntersectionResult {
        intersections: results,
    })
}

#[cfg(test)]
#[path = "surface_surface_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "surface_surface_degenerate_tests.rs"]
mod degenerate_tests;
