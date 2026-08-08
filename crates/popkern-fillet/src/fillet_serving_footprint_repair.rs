//! PK-251 serving-boundary footprint repair for direct-apex caps.
//!
//! # Why the repair lives at the serving boundary
//!
//! The construction-time clearance corrections
//! (`exact_boundary_cap::strip_footprint_clearance`) measure the cap
//! against an analytic model of each strip's footprint, but the model
//! cannot be exact at candidate time: the strip trims the candidate
//! sees are pre-rewrite (the cylinder trim rewrite runs after the cap
//! commit), and the served trim classification — tilted rim arcs,
//! angular borders, shared-vertex exemptions — belongs to the
//! validator. Measured on the PK-251 fixtures, every model
//! refinement still disagreed with the served body in one direction
//! or the other. The serving boundary is where the truth lives, so
//! the final repair runs there, driven by the REAL validator:
//! [`popkern_core::topology::BRep::interior_overlap_samples`] returns
//! the Interior-classified overlap points for the offending pair, and
//! the repair pushes the cap's interior control points out of the
//! counterpart face's material (along the analytic surface's outward
//! signed-distance gradient, oriented by the counterpart's honest
//! orientation flag) at exactly those points.
//!
//! # Contract
//!
//! - Only the cap sector face's interior control points move — its
//!   boundary rows carry the served edge curves and never move, so
//!   edge/pcurve exactness (`CurveSurfaceDeviation`) is preserved by
//!   construction.
//! - The counterpart face is untouched.
//! - The repair is accepted only when the scoped region validator
//!   comes back clean; otherwise the caller falls through to the
//!   PK-250 exclusion retry (legacy serves), exactly as before the
//!   repair existed. A repair can therefore only upgrade a serve.
//! - Every repair is recorded on the repair report
//!   (`serving_footprint_repairs`), cleared or not.

use popkern_core::error::ValidationError;
use popkern_core::geometry::analytic_surfaces::analytic_signed_distance;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{NurbsSurface, Point3, Vec3};
use popkern_core::topology::{BRep, FaceHandle};

use crate::apex_construction::exact_boundary_cap::cap_quality_postcondition::candidate_region;
use crate::types::ServingFootprintRepair;

/// Sample rounds: each runs the pair intersection + trim
/// classification (the expensive step). Pushes near the locked rim
/// rows act through small basis weights, so convergence takes tens
/// of push applications (measured: the [0.02, 0.03, 0.04] cube
/// corner cleared after 38 single-push rounds); each sample round
/// therefore drives several push applications whose control-point
/// depths re-measure cheaply from the updated surface. The loop
/// exits as soon as the region validates, so the budget is only
/// spent on bodies that are still defective.
const MAX_REPAIR_ROUNDS: usize = 64;
/// Push applications per sample round.
const PUSHES_PER_SAMPLE_ROUND: usize = 4;
/// Consecutive rounds the SAME face pair may persist before the
/// repair declares a stall. Repairable interior overlaps clear in a
/// handful of rounds (measured 1..5 on the PK-251 fixtures); a pair
/// that survives this many is boundary participation — the crossing
/// runs through a locked boundary curve interior pushes cannot move
/// (measured: the cube [0.15, 0.1, 0.1] spoke dive held its pair for
/// 96 rounds and 972 pushes without clearing).
const MAX_STALLED_ROUNDS: usize = 24;
/// Push overshoot over the measured penetration depth.
const PUSH_OVERSHOOT: f64 = 1.5;
/// Ceiling on the response-scaled push multiplier (see the
/// construction-time corrections for the measured rationale).
const MAX_PUSH_MULTIPLIER: f64 = 3.0;

/// Attempt an in-place footprint repair of one served direct-apex
/// corner region. Returns the repair record; `cleared` is true when
/// the scoped region validator is clean afterwards.
pub(crate) fn repair_served_direct_apex_footprint(
    brep: &mut BRep,
    vertex: popkern_core::topology::VertexHandle,
    built_faces: &[FaceHandle],
) -> ServingFootprintRepair {
    let region = candidate_region(brep, built_faces);
    let mut rounds = 0;
    let mut corrected_control_points = 0;
    let mut detail = String::new();
    // Full-region validation intersects every non-adjacent pair —
    // including NURBS×NURBS cap-cap pairs, the expensive class — so
    // it runs only to DISCOVER the next offending pair. The push
    // iteration re-validates just the offending pair (one
    // cylinder×NURBS intersection per round; measured ~10× cheaper
    // than re-validating the region every round).
    let mut cleared_pairs: Vec<(FaceHandle, FaceHandle)> = Vec::new();
    let cleared = 'region: loop {
        let Err(error) = brep.validate_self_intersection_among(&region) else {
            break true;
        };
        if rounds >= MAX_REPAIR_ROUNDS {
            detail = format!("exhausted: {error}");
            break false;
        }
        let ValidationError::SelfIntersection { face_a, face_b } = error else {
            detail = format!("non-pair validation error: {error}");
            break false;
        };
        // A pair the pair-scoped check declared clean must not come
        // back from the region check — that is a pair-vs-region
        // verdict disagreement, and looping on it would spin without
        // progress. Stop with the evidence instead.
        if cleared_pairs.contains(&(face_a, face_b)) {
            detail = format!("pair/region verdict disagreement on {face_a:?}/{face_b:?}");
            break false;
        }
        cleared_pairs.push((face_a, face_b));
        // The movable face is the built cap sector (NURBS); the other
        // is expected to be an analytic support (strip cylinder,
        // bridge plane, corner sphere, ...).
        let (cap, counterpart) = if built_faces.contains(&face_a) {
            (face_a, face_b)
        } else if built_faces.contains(&face_b) {
            (face_b, face_a)
        } else {
            detail = format!("pair outside built faces: {face_a:?}/{face_b:?}");
            break false;
        };
        let Some(counterpart_surface) = brep
            .faces
            .get(counterpart)
            .and_then(|face| face.surface.as_ref())
            .map(|sof| (sof.geometry.clone(), sof.orientation))
        else {
            detail = format!("stale counterpart: {counterpart:?}");
            break false;
        };
        if analytic_signed_distance(&counterpart_surface.0, &Point3::new(0.0, 0.0, 0.0)).is_none() {
            detail = format!("non-analytic counterpart: {counterpart:?}");
            break false;
        }
        let mut pair_rounds = 0;
        loop {
            if brep
                .validate_self_intersection_among(&[cap, counterpart])
                .is_ok()
            {
                // This pair cleared; rediscover the next one (if any).
                continue 'region;
            }
            if rounds >= MAX_REPAIR_ROUNDS {
                detail = format!("exhausted on {cap:?}/{counterpart:?}");
                break 'region false;
            }
            if pair_rounds >= MAX_STALLED_ROUNDS {
                detail = format!(
                    "stalled (boundary participation): {cap:?}/{counterpart:?} after {pair_rounds} rounds"
                );
                break 'region false;
            }
            let samples = brep.interior_overlap_samples(cap, counterpart);
            if samples.is_empty() {
                detail = format!("no interior samples for {cap:?}/{counterpart:?}");
                break 'region false;
            }
            // Several push applications per sample round: the sample
            // UVs stay valid as attribution anchors while the
            // control-point depths re-measure cheaply from the
            // updated surface, so the expensive
            // intersection/classification amortizes.
            let mut moved_this_round = 0;
            for _ in 0..PUSHES_PER_SAMPLE_ROUND {
                let moved = push_cap_out_of_counterpart(
                    brep,
                    cap,
                    &samples,
                    &counterpart_surface.0,
                    counterpart_surface.1,
                );
                if moved == 0 {
                    break;
                }
                moved_this_round += moved;
            }
            if moved_this_round == 0 {
                detail = format!("no movable control points for {cap:?}/{counterpart:?}");
                break 'region false;
            }
            corrected_control_points += moved_this_round;
            rounds += 1;
            pair_rounds += 1;
        }
    };
    if cleared {
        detail = "cleared".to_string();
    }
    ServingFootprintRepair {
        vertex,
        rounds,
        corrected_control_points,
        cleared,
        detail,
    }
}

/// Push the cap face's interior control points out of the analytic
/// counterpart's material at the given overlap samples. Boundary rows
/// and columns never move. Returns the number of control points moved.
fn push_cap_out_of_counterpart(
    brep: &mut BRep,
    cap: FaceHandle,
    samples: &[popkern_core::validation::InteriorOverlapSample],
    counterpart_surface: &SurfaceGeometry,
    counterpart_orientation: bool,
) -> usize {
    let Some(SurfaceGeometry::Nurbs(nurbs)) = brep
        .faces
        .get(cap)
        .and_then(|face| face.surface.as_ref())
        .map(|sof| sof.geometry.clone())
    else {
        return 0;
    };
    let mut surface = nurbs;
    let num_u = surface.num_u;
    let num_v = surface.num_v;
    if num_u < 4 || num_v < 4 {
        return 0;
    }
    let greville_u = greville(&surface.knots_u, surface.degree_u, num_u);
    let greville_v = greville(&surface.knots_v, surface.degree_v, num_v);
    // Deepest penetration per interior control point, attributed by
    // the sample's exact cap-side UV from the intersection curve —
    // NOT by iterative 3D→UV projection, which misattributes on
    // free-form surfaces (measured: 48 scattered-push rounds without
    // clearing).
    let mut pushes: std::collections::BTreeMap<(usize, usize), (f64, Vec3, f64, f64)> =
        std::collections::BTreeMap::new();
    for sample in samples {
        // Caller order is (cap, counterpart): the cap's UV is `uv_a`.
        let (u, v) = sample.uv_a;
        if !u.is_finite() || !v.is_finite() {
            continue;
        }
        let ci = nearest(&greville_u, u).clamp(1, num_u - 2);
        let cj = nearest(&greville_v, v).clamp(1, num_v - 2);
        // Samples lie ON the intersection curve — depth zero by
        // definition. The dip to remove lives BETWEEN the crossings,
        // so measure it on the cap surface at the attributed control
        // point's Greville parameters, with a tolerance-scaled floor
        // so grazing overlaps still make progress.
        use popkern_core::geometry::SurfaceEval;
        let cap_point =
            SurfaceGeometry::Nurbs(surface.clone()).point_at(greville_u[ci], greville_v[cj]);
        let Some((sd_out, out_direction)) =
            counterpart_outward(counterpart_surface, counterpart_orientation, &cap_point)
        else {
            continue;
        };
        let depth = (-sd_out).max(50.0 * brep.tolerance.linear());
        let entry = pushes
            .entry((ci, cj))
            .or_insert((depth, out_direction, u, v));
        if depth > entry.0 {
            *entry = (depth, out_direction, u, v);
        }
    }
    if pushes.is_empty() {
        return 0;
    }
    let reference = surface.clone();
    let mut moved = 0;
    for ((ci, cj), (depth, radial_out, u, v)) in pushes {
        let index = ci * num_v + cj;
        let response = response_at(&reference, index, u, v, &radial_out);
        let multiplier = (PUSH_OVERSHOOT / response.max(1.0e-3)).min(MAX_PUSH_MULTIPLIER);
        surface.control_points[index] =
            surface.control_points[index] + radial_out * (depth * multiplier);
        moved += 1;
    }
    if moved > 0
        && let Some(sof) = brep
            .faces
            .get_mut(cap)
            .and_then(|face| face.surface.as_mut())
    {
        sof.geometry = SurfaceGeometry::Nurbs(surface);
    }
    moved
}

/// Outward-of-material signed distance and direction at `p` relative
/// to the counterpart's analytic surface: the implicit gradient
/// (finite differences — generic across analytic types), oriented so
/// it points OUT of the counterpart face's material using its honest
/// orientation flag (PK-244/PK-245). Positive `sd_out` means `p` is
/// already outside the material.
fn counterpart_outward(
    surface: &SurfaceGeometry,
    orientation: bool,
    p: &Point3,
) -> Option<(f64, Vec3)> {
    use popkern_core::geometry::SurfaceEval;
    use popkern_core::geometry::analytic_surfaces::analytic_signed_distance;
    let sd = analytic_signed_distance(surface, p)?;
    const H: f64 = 1.0e-6;
    let gradient = Vec3::new(
        (analytic_signed_distance(surface, &Point3::new(p.x + H, p.y, p.z))? - sd) / H,
        (analytic_signed_distance(surface, &Point3::new(p.x, p.y + H, p.z))? - sd) / H,
        (analytic_signed_distance(surface, &Point3::new(p.x, p.y, p.z + H))? - sd) / H,
    );
    let gradient = gradient.normalized()?;
    // Orient by the face's outward normal at the projection of `p`.
    let (u, v) = popkern_core::geometry::uv_projection::project_point_to_uv(surface, p)?;
    let normal = surface.normal_at(u, v);
    let outward_normal = if orientation { normal } else { normal * -1.0 };
    let sign = if outward_normal.dot(&gradient) >= 0.0 {
        1.0
    } else {
        -1.0
    };
    Some((sd * sign, gradient * sign))
}

fn response_at(surface: &NurbsSurface, index: usize, u: f64, v: f64, direction: &Vec3) -> f64 {
    use popkern_core::geometry::SurfaceEval;
    let mut probe = surface.clone();
    probe.control_points[index] = probe.control_points[index] + *direction;
    let base = SurfaceGeometry::Nurbs(surface.clone());
    let probed = SurfaceGeometry::Nurbs(probe);
    (probed.point_at(u, v) - base.point_at(u, v)).dot(direction)
}

fn greville(knots: &[f64], degree: usize, count: usize) -> Vec<f64> {
    (0..count)
        .map(|i| {
            if degree == 0 {
                knots.get(i + 1).copied().unwrap_or(0.0)
            } else {
                knots[i + 1..=i + degree].iter().sum::<f64>() / degree as f64
            }
        })
        .collect()
}

fn nearest(greville: &[f64], parameter: f64) -> usize {
    let mut best = 0;
    let mut best_distance = f64::INFINITY;
    for (index, &abscissa) in greville.iter().enumerate() {
        let distance = (abscissa - parameter).abs();
        if distance < best_distance {
            best_distance = distance;
            best = index;
        }
    }
    best
}
