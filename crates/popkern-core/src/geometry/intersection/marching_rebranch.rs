//! Wrong-preimage re-branching for unclamped Newton landings.
//!
//! Newton runs unclamped by design, and on a self-overlapping
//! extended basis it can converge onto a wrong preimage of the true
//! intersection point (see [`rebranch_out_of_window`]).

use super::MarchSample;
use crate::geometry::Point3;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::tolerance::Tolerance;

/// Re-branch a converged sample whose surface parameters landed FAR
/// outside the query windows onto the in-window preimage of the same
/// 3D point, when one exists.
///
/// Newton runs unclamped by design, and on a self-overlapping
/// extended basis it can converge onto a wrong preimage: a thread
/// root band's linear profile extension covers the same helical
/// cylinder as the band itself, and a cone x band solve converged at
/// u = -16 profile-widths — the 3D point was ON the true
/// intersection, but the domain check killed the march and the pair
/// emitted nothing (measured; whether it happened flipped on a
/// 5e-17 cone-axis tilt — basin luck, not correctness). Re-project
/// the converged POINT into the window (seeded at the window
/// centre); adopt the preimage only when it reproduces the point
/// within tolerance. Landings within 10% of a span outside the
/// window are the tuned slightly-out class and are left alone.
#[allow(clippy::too_many_arguments)]
pub(crate) fn rebranch_out_of_window(
    surface1: &SurfaceGeometry,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    surface2: &SurfaceGeometry,
    u2_range: (f64, f64),
    v2_range: (f64, f64),
    sample: &mut MarchSample,
    tolerance: &Tolerance,
    phantom_cache: &mut Vec<Point3>,
) -> bool {
    // PHANTOM CACHE: consecutive Newton solves fall into the SAME
    // extrapolation basin (1428 seed solves and whole march retry
    // ladders converged onto near-identical off-surface points —
    // ~800k multisheet projections on one two-start rod cut, 5× wall
    // time). A point within the fit-slop band of the last measured
    // phantom is the same phantom; reject it without projecting.
    if phantom_cache
        .iter()
        .any(|prev| sample.point.distance_to(prev) < tolerance.linear() * 100.0)
    {
        return false;
    }
    // Only NON-PERIODIC axes can host an extrapolation phantom: a
    // periodic (or geometrically closed) axis wraps — every branch is
    // the same surface point, the march's wrap_parameter owns it, and
    // treating pre-wrap coordinates as "far" made every step of a
    // multi-turn helical march pay a multisheet projection (measured:
    // the two-start rod serve went 24s → 121s).
    let (p1u, p1v) = crate::geometry::chart_periods::marching_periodicity(surface1);
    let (p2u, p2v) = crate::geometry::chart_periods::marching_periodicity(surface2);
    let far = |x: f64, r: (f64, f64), periodic: bool| -> bool {
        if periodic {
            return false;
        }
        let span = (r.1 - r.0).abs().max(1e-12);
        x < r.0 - 0.1 * span || x > r.1 + 0.1 * span
    };
    let fix = |surface: &SurfaceGeometry,
               ur: (f64, f64),
               vr: (f64, f64),
               periodic: (bool, bool),
               u: &mut f64,
               v: &mut f64,
               p: &Point3|
     -> bool {
        if !far(*u, ur, periodic.0) && !far(*v, vr, periodic.1) {
            return true;
        }
        // HULL PROOF (PK-513): a NURBS surface lies in the union of
        // its span control-net boxes, and adoption below requires a
        // preimage within 10x linear. A point farther than that from
        // every span box can never re-branch — reject it without the
        // multi-start projection (phantom seeds on a thread wall's
        // extended basis paid one each, ~18k per chamfer cut).
        if let SurfaceGeometry::Nurbs(ns) = surface
            && let Some(field) = ns.projection_seed_field()
            && !field.may_reach(p, tolerance.linear() * 10.0)
        {
            return false;
        }
        // The multi-sheet projector (PK-446/PK-452): a window-centre
        // seeded local walk lands on the adjacent coil turn of a
        // helical basis; the multi-start projector picks the true
        // sheet.
        let seed = (0.5 * (ur.0 + ur.1), 0.5 * (vr.0 + vr.1));
        let candidate = crate::geometry::uv_projection::project_point_to_uv_multisheet(surface, p)
            .or_else(|| {
                crate::geometry::uv_projection::project_point_to_uv_near(surface, p, seed.0, seed.1)
            });
        if crate::trace::trace_enabled("march") {
            eprintln!(
                "      [rebranch] far ({:.3},{:.3}) w {ur:?} {vr:?} -> {candidate:?} d={:?}",
                *u,
                *v,
                candidate.map(|(nu, nv)| surface.point_at(nu, nv).distance_to(p))
            );
        }
        if let Some((nu, nv)) = candidate
            && !far(nu, ur, periodic.0)
            && !far(nv, vr, periodic.1)
            && surface.point_at(nu, nv).distance_to(p) < tolerance.linear() * 10.0
        {
            *u = nu;
            *v = nv;
            return true;
        }
        // Far outside AND the point does not re-project onto the
        // in-window surface: the solve converged on the EXTRAPOLATED
        // basis — a phantom, not an intersection of the real faces.
        false
    };
    let point = sample.point;
    let ok1 = fix(
        surface1,
        u1_range,
        v1_range,
        (p1u, p1v),
        &mut sample.u1,
        &mut sample.v1,
        &point,
    );
    let ok2 = ok1
        && fix(
            surface2,
            u2_range,
            v2_range,
            (p2u, p2v),
            &mut sample.u2,
            &mut sample.v2,
            &point,
        );
    let ok = ok1 && ok2;
    if !ok {
        if phantom_cache.len() >= 16 {
            phantom_cache.remove(0);
        }
        phantom_cache.push(point);
    }
    ok
}

/// Public-in-crate wrapper over [`rebranch_out_of_window`] for the
/// seed phase (see the call in `surface_surface.rs`).
#[allow(clippy::too_many_arguments)]
pub(crate) fn rebranch_sample_pub(
    surface1: &SurfaceGeometry,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    surface2: &SurfaceGeometry,
    u2_range: (f64, f64),
    v2_range: (f64, f64),
    sample: &mut MarchSample,
    tolerance: &Tolerance,
    phantom_cache: &mut Vec<Point3>,
) -> bool {
    rebranch_out_of_window(
        surface1,
        u1_range,
        v1_range,
        surface2,
        u2_range,
        v2_range,
        sample,
        tolerance,
        phantom_cache,
    )
}
