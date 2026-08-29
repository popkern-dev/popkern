//! Starting point detection for surface-surface intersection.

use crate::geometry::Vec3;
use crate::geometry::intersection::marching::{MarchSample, march_project};
use crate::geometry::intersection::surface_subdivision::{
    SurfacePatch, decompose_to_patches, refine_patch,
};
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::tolerance::Tolerance;

/// Maximum candidate patch pairs during AABB subdivision. Prevents
/// memory explosion for coincident or near-coincident surfaces.
const MAX_CANDIDATE_PAIRS: usize = 50_000;

/// A classified starting point for surface-surface intersection marching.
#[allow(dead_code)]
pub(crate) struct StartingPoint {
    pub sample: MarchSample,
    pub is_tangent: bool,
}

/// Find starting points for surface-surface intersection.
///
/// Uses AABB-based subdivision to narrow down candidate region pairs, then applies
/// Newton 4×4 iteration (via `march_project`) to refine intersection points.
/// Results are deduplicated and classified as transversal or tangent.
#[allow(dead_code)]
pub(crate) fn find_starting_points(
    surface1: &SurfaceGeometry,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    surface2: &SurfaceGeometry,
    u2_range: (f64, f64),
    v2_range: (f64, f64),
    tolerance: &Tolerance,
) -> Vec<StartingPoint> {
    // Step 1: Promote both surfaces to NURBS over the given parameter ranges.
    let nurbs1 = surface1.to_nurbs_in_range(u1_range, v1_range);
    let nurbs2 = surface2.to_nurbs_in_range(u2_range, v2_range);

    // Step 2: Decompose both into patches at internal knots.
    let patches1 = decompose_to_patches(&nurbs1, u1_range, v1_range);
    let patches2 = decompose_to_patches(&nurbs2, u2_range, v2_range);

    // Step 3 & 4: AABB culling with adaptive refinement.
    // We process a queue of (patch1, patch2, depth) pairs.
    let mut candidate_pairs: Vec<(SurfacePatch, SurfacePatch)> = Vec::new();

    let initial_pairs: Vec<(SurfacePatch, SurfacePatch, u32)> = patches1
        .iter()
        .flat_map(|p1| {
            patches2.iter().filter_map(move |p2| {
                if p1.bbox.overlaps(&p2.bbox, tolerance.linear()) {
                    Some((p1.clone(), p2.clone(), 0u32))
                } else {
                    None
                }
            })
        })
        .collect();

    if seed_trace_enabled() {
        let hull = |ps: &[SurfacePatch]| -> (f64, f64) {
            let mut lo = f64::MAX;
            let mut hi = f64::MIN;
            for p in ps {
                lo = lo.min(p.bbox.min.x);
                hi = hi.max(p.bbox.max.x);
            }
            (lo, hi)
        };
        let (a0, a1) = hull(&patches1);
        let (b0, b1) = hull(&patches2);
        eprintln!(
            "  [seed] patches {}x{} initial overlapping pairs {} | x-hulls a[{a0:.3},{a1:.3}] b[{b0:.3},{b1:.3}]",
            patches1.len(),
            patches2.len(),
            initial_pairs.len()
        );
    }
    let mut queue = initial_pairs;

    // Stage 2b.56: depth cap reduced from 6 → 4 and diag threshold
    // raised from 100·tol to 1000·tol. The previous values produced
    // ~1000 candidate pairs on small NURBS fillet face inputs (each
    // requiring a Newton iteration via `march_project`), pushing
    // total `intersect_surfaces` time into the 600+ ms regime —
    // diagnosed in stage 2b.55. The new values keep enough
    // candidate density to find every CAD intersection branch
    // while putting slow calls in the ~50 ms regime.
    while let Some((patch1, patch2, depth)) = queue.pop() {
        // Defense-in-depth: cap total candidates to prevent OOM.
        if candidate_pairs.len() >= MAX_CANDIDATE_PAIRS {
            break;
        }
        let diag = patch1.bbox.diagonal().max(patch2.bbox.diagonal());
        if depth < 4 && diag > 1000.0 * tolerance.linear() {
            // Refine both patches and re-cull sub-pairs.
            let sub1 = refine_patch(&patch1);
            let sub2 = refine_patch(&patch2);
            for sp1 in &sub1 {
                for sp2 in &sub2 {
                    if sp1.bbox.overlaps(&sp2.bbox, tolerance.linear()) {
                        queue.push((sp1.clone(), sp2.clone(), depth + 1));
                    }
                }
            }
        } else {
            candidate_pairs.push((patch1, patch2));
        }
    }

    // Step 5: For each small overlapping pair, try Newton convergence.
    //
    // Newton runs unclamped and converged seeds are accepted wherever
    // they land: clamping iterates stalls the constrained system, and
    // rejecting out-of-window convergences loses real curves — the
    // march itself clips to the query windows, so a seed slightly
    // outside them still traces the in-window stretch (both measured
    // on the served pyramid4 strip×cap pair).
    let mut raw_samples: Vec<MarchSample> = Vec::new();

    if seed_trace_enabled() {
        eprintln!("  [seed] candidates {}", candidate_pairs.len());
    }
    for (patch1, patch2) in &candidate_pairs {
        // Use the midpoint of each patch as the initial guess.
        let u1_mid = (patch1.u_range.0 + patch1.u_range.1) / 2.0;
        let v1_mid = (patch1.v_range.0 + patch1.v_range.1) / 2.0;
        let u2_mid = (patch2.u_range.0 + patch2.u_range.1) / 2.0;
        let v2_mid = (patch2.v_range.0 + patch2.v_range.1) / 2.0;

        // Compute a march direction: cross-product of the two normals gives the
        // tangent direction of the intersection curve. If normals are nearly
        // parallel (tangent case), use an arbitrary non-degenerate direction.
        let n1 = surface1.normal_at(u1_mid, v1_mid);
        let n2 = surface2.normal_at(u2_mid, v2_mid);

        // Newton's advancing plane needs an anchor point, and the two
        // natural choices — either patch midpoint — pose DIFFERENT
        // constrained systems that can converge to different seeds or
        // diverge outright. Anchoring only on operand 1 made the seed
        // set depend on operand order (measured: adjacent strip×cap
        // pairs found their tangential contact in one order and
        // nothing in the other). Both anchors always run; the union is
        // order-symmetric by construction and deduplication collapses
        // agreeing hits.
        let anchor_1 = surface1.point_at(u1_mid, v1_mid);
        let dir_1 = pick_march_dir(n1, n2, surface1, u1_mid, v1_mid);
        if let Some(sample) = march_project(
            surface1, u1_mid, v1_mid, surface2, u2_mid, v2_mid, &anchor_1, &dir_1, tolerance,
        ) {
            raw_samples.push(sample);
        }
        let anchor_2 = surface2.point_at(u2_mid, v2_mid);
        let dir_2 = pick_march_dir(n2, n1, surface2, u2_mid, v2_mid);
        if let Some(sample) = march_project(
            surface1, u1_mid, v1_mid, surface2, u2_mid, v2_mid, &anchor_2, &dir_2, tolerance,
        ) {
            raw_samples.push(sample);
        }
    }

    // Step 6: Deduplicate points within linear tolerance.
    if seed_trace_enabled() {
        eprintln!("  [seed] raw samples {}", raw_samples.len());
    }
    let unique_samples = deduplicate(raw_samples, tolerance);

    // Step 7: Classify each point as transversal or tangent.
    unique_samples
        .into_iter()
        .map(|sample| {
            let is_tangent = classify_tangent(surface1, surface2, &sample, tolerance);
            StartingPoint { sample, is_tangent }
        })
        .collect()
}

/// Choose a march direction for Newton convergence.
///
/// The natural march direction is N1 × N2 (intersection curve tangent).
/// If the normals are nearly parallel (tangent intersection or coincident surfaces),
/// fall back to the first non-degenerate partial derivative of surface 1.
#[allow(dead_code)]
fn pick_march_dir(n1: Vec3, n2: Vec3, surface1: &SurfaceGeometry, u1: f64, v1: f64) -> Vec3 {
    let cross = n1.cross(&n2);
    if cross.length() > 1e-8 {
        return cross;
    }

    // Fallback: use a surface partial (guaranteed non-degenerate for regular surfaces).
    let du = surface1.partial_u(u1, v1);
    if du.length() > 1e-8 {
        return du;
    }
    let dv = surface1.partial_v(u1, v1);
    if dv.length() > 1e-8 {
        return dv;
    }

    // Last resort: cardinal X axis.
    Vec3::new(1.0, 0.0, 0.0)
}

/// Remove duplicate samples: keep first occurrence, skip any within `tolerance.linear()`.
#[allow(dead_code)]
fn deduplicate(samples: Vec<MarchSample>, tolerance: &Tolerance) -> Vec<MarchSample> {
    let mut unique: Vec<MarchSample> = Vec::new();
    'outer: for s in samples {
        for u in &unique {
            if s.point.distance_to(&u.point) < tolerance.linear() {
                continue 'outer;
            }
        }
        unique.push(s);
    }
    unique
}

/// Check for additional intersection branches near a starting point.
///
/// Perturbs the starting point in four cardinal directions in surface 1's parameter space
/// and applies Newton projection to discover intersection points that lie on branches
/// separate from the already-traced curve. Returns candidates not already covered by
/// `traced_samples`.
#[allow(dead_code, clippy::too_many_arguments)]
pub(crate) fn find_branch_points(
    surface1: &SurfaceGeometry,
    _u1_range: (f64, f64),
    _v1_range: (f64, f64),
    surface2: &SurfaceGeometry,
    _u2_range: (f64, f64),
    _v2_range: (f64, f64),
    traced_samples: &[MarchSample],
    start: &MarchSample,
    tolerance: &Tolerance,
) -> Vec<MarchSample> {
    let delta = tolerance.linear() * 50.0;
    let mut branches = Vec::new();

    // Perturb in 4 directions in surface 1's parameter space.
    let perturbations = [(delta, 0.0), (-delta, 0.0), (0.0, delta), (0.0, -delta)];
    for (du, dv) in &perturbations {
        let u1 = start.u1 + du;
        let v1 = start.v1 + dv;

        // Use a march direction derived from the surface normals at this perturbed point.
        let n1 = surface1.normal_at(u1, v1);
        let n2 = surface2.normal_at(start.u2, start.v2);
        let march_dir = pick_march_dir(n1, n2, surface1, u1, v1);

        // Newton project to find a nearby intersection.
        if let Some(candidate) = march_project(
            surface1,
            u1,
            v1,
            surface2,
            start.u2,
            start.v2,
            &start.point,
            &march_dir,
            tolerance,
        ) {
            // Only keep it if it is not already on the existing traced curve.
            let on_existing = traced_samples
                .iter()
                .any(|s| s.point.distance_to(&candidate.point) < tolerance.linear() * 10.0);
            if !on_existing {
                branches.push(candidate);
            }
        }
    }
    branches
}

/// Returns true if the two surface normals at the sample point are nearly parallel,
/// indicating a tangent intersection.
#[allow(dead_code)]
fn classify_tangent(
    surface1: &SurfaceGeometry,
    surface2: &SurfaceGeometry,
    sample: &MarchSample,
    tolerance: &Tolerance,
) -> bool {
    let n1 = surface1.normal_at(sample.u1, sample.v1);
    let n2 = surface2.normal_at(sample.u2, sample.v2);

    let len1 = n1.length();
    let len2 = n2.length();

    // Degenerate normals can't be classified; treat as non-tangent.
    if len1 < tolerance.numerical_floor() || len2 < tolerance.numerical_floor() {
        return false;
    }

    let cos_angle = n1.dot(&n2).abs() / (len1 * len2);
    cos_angle > 1.0 - tolerance.angular()
}

#[cfg(test)]
#[path = "starting_points_tests.rs"]
mod tests;

/// OnceLock'd the `seed` trace topic gate: a per-call env read takes the
/// process environment lock and serializes parallel SSI dispatch
/// (PK-329).
fn seed_trace_enabled() -> bool {
    static TRACE: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *TRACE.get_or_init(|| crate::trace::trace_enabled("seed"))
}
