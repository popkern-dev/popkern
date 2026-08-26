//! Certified-separation pre-gate for expensive self-intersection
//! pairs (PK-327).
//!
//! Free-form (NURBS/Gregory) pairs route to the marching SSI solver,
//! which degrades badly on near-parallel non-intersecting sheets —
//! exactly what a healthy thickened free-form face produces (two
//! fitted skins a thickness apart; measured 77 s for one wavy-plate
//! validate_all before this gate). Subdivision with sampled,
//! gap-inflated AABBs certifies "minimum pair distance exceeds the
//! clearance" quickly for such pairs; anything the certificate cannot
//! prove falls through to the full solver, so the gate can only
//! REMOVE work, never verdicts. The AABB inflation uses the full
//! largest adjacent-sample gap (not half), keeping the enclosure at
//! least as conservative as the sampled machinery it short-circuits.

use crate::geometry::Point3;
use crate::geometry::intersection::surface_subdivision::Aabb;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};

/// Subdivision budget: pairs examined before giving up (→ `false`).
const MAX_BOXES: usize = 50_000;

/// Maximum subdivision depth per pair.
const MAX_DEPTH: usize = 7;

/// True when the minimum distance between the two windowed patches is
/// PROVABLY greater than `clearance`. `false` means "could not
/// certify" — never "they intersect".
pub(super) fn pair_clearly_separated(
    surface_a: &SurfaceGeometry,
    ua: (f64, f64),
    va: (f64, f64),
    surface_b: &SurfaceGeometry,
    ub: (f64, f64),
    vb: (f64, f64),
    clearance: f64,
) -> bool {
    type Box2 = ((f64, f64), (f64, f64));
    let mut stack: Vec<(Box2, Box2, usize)> = vec![((ua, va), (ub, vb), 0)];
    let mut examined = 0usize;
    while let Some((ba, bb, depth)) = stack.pop() {
        examined += 1;
        if examined >= MAX_BOXES {
            return false;
        }
        let (Some(aab), Some(bab)) = (sample_aabb(surface_a, ba), sample_aabb(surface_b, bb))
        else {
            return false;
        };
        if aabb_gap(&aab, &bab) > clearance {
            continue; // this pair is certified apart
        }
        if depth >= MAX_DEPTH {
            return false; // cannot certify at the depth budget
        }
        // Split the geometrically larger box.
        let a_size = aabb_diag(&aab);
        let b_size = aabb_diag(&bab);
        if a_size >= b_size {
            for half in split_box(ba) {
                stack.push((half, bb, depth + 1));
            }
        } else {
            for half in split_box(bb) {
                stack.push((ba, half, depth + 1));
            }
        }
    }
    true
}

fn split_box(((u0, u1), (v0, v1)): ((f64, f64), (f64, f64))) -> [((f64, f64), (f64, f64)); 2] {
    if (u1 - u0) >= (v1 - v0) {
        let um = (u0 + u1) * 0.5;
        [((u0, um), (v0, v1)), ((um, u1), (v0, v1))]
    } else {
        let vm = (v0 + v1) * 0.5;
        [((u0, u1), (v0, vm)), ((u0, u1), (vm, v1))]
    }
}

/// Sampled AABB inflated by the largest adjacent-sample gap.
fn sample_aabb(
    surface: &SurfaceGeometry,
    ((u0, u1), (v0, v1)): ((f64, f64), (f64, f64)),
) -> Option<Aabb> {
    const N: usize = 3;
    let mut points = Vec::with_capacity((N + 1) * (N + 1));
    for i in 0..=N {
        for j in 0..=N {
            let u = u0 + (u1 - u0) * (i as f64) / N as f64;
            let v = v0 + (v1 - v0) * (j as f64) / N as f64;
            points.push(surface.point_at(u, v));
        }
    }
    let mut max_gap: f64 = 0.0;
    for i in 0..=N {
        for j in 0..=N {
            let p = points[i * (N + 1) + j];
            if i < N {
                max_gap = max_gap.max(p.distance_to(&points[(i + 1) * (N + 1) + j]));
            }
            if j < N {
                max_gap = max_gap.max(p.distance_to(&points[i * (N + 1) + j + 1]));
            }
        }
    }
    let aabb = Aabb::from_points(&points)?;
    Some(Aabb {
        min: Point3::new(
            aabb.min.x - max_gap,
            aabb.min.y - max_gap,
            aabb.min.z - max_gap,
        ),
        max: Point3::new(
            aabb.max.x + max_gap,
            aabb.max.y + max_gap,
            aabb.max.z + max_gap,
        ),
    })
}

fn aabb_gap(a: &Aabb, b: &Aabb) -> f64 {
    let axis_gap = |a0: f64, a1: f64, b0: f64, b1: f64| (b0 - a1).max(a0 - b1).max(0.0);
    let dx = axis_gap(a.min.x, a.max.x, b.min.x, b.max.x);
    let dy = axis_gap(a.min.y, a.max.y, b.min.y, b.max.y);
    let dz = axis_gap(a.min.z, a.max.z, b.min.z, b.max.z);
    (dx * dx + dy * dy + dz * dz).sqrt()
}

fn aabb_diag(a: &Aabb) -> f64 {
    let dx = a.max.x - a.min.x;
    let dy = a.max.y - a.min.y;
    let dz = a.max.z - a.min.z;
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// True when the free-form surface's signed distance to the analytic
/// surface never crosses the tolerance band two-sidedly — tangential
/// contact with no material interpenetration. Applies only when one
/// side has an analytic signed distance and the other is free-form;
/// anything else returns `false` (no verdict).
pub(super) fn one_sided_contact_only(
    surface_a: &SurfaceGeometry,
    surface_b: &SurfaceGeometry,
    b_u_range: (f64, f64),
    b_v_range: (f64, f64),
    linear_tol: f64,
) -> bool {
    use crate::geometry::analytic_surfaces::analytic_signed_distance;
    let origin = Point3::new(0.0, 0.0, 0.0);
    let (analytic, free_form, (ur, vr)) = if analytic_signed_distance(surface_a, &origin).is_some()
        && matches!(
            surface_b,
            SurfaceGeometry::Nurbs(_) | SurfaceGeometry::Gregory(_)
        ) {
        (surface_a, surface_b, (b_u_range, b_v_range))
    } else {
        return false;
    };
    const N: usize = 48;
    let mut min_sd = f64::MAX;
    let mut max_sd = f64::MIN;
    for i in 0..N {
        for j in 0..N {
            let u = ur.0 + (ur.1 - ur.0) * (i as f64) / (N - 1) as f64;
            let v = vr.0 + (vr.1 - vr.0) * (j as f64) / (N - 1) as f64;
            let p = free_form.point_at(u, v);
            let Some(sd) = analytic_signed_distance(analytic, &p) else {
                return false;
            };
            min_sd = min_sd.min(sd);
            max_sd = max_sd.max(sd);
        }
    }
    min_sd > -linear_tol || max_sd < linear_tol
}
