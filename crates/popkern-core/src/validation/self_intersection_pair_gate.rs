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
use popkern_math::PortableMath;

/// Subdivision budget: pairs examined before giving up (→ `false`).
const MAX_BOXES: usize = 50_000;

/// Maximum subdivision depth per pair: halvings of EITHER patch along
/// one path, so a pair of patches each needing `k` halvings per
/// parameter to resolve the gap needs about `4k`. At 7 a fitted corner
/// sector 3e-2 from a non-adjacent fillet strip (sampled-gap inflation
/// needs boxes of about the gap) was never certified and fell to the
/// marching solver, which proved "no curve" in 0.5-14 s per pair
/// (PK-617, the dual-jet sectors' rational form the slowest). At 24
/// those pairs certify in milliseconds. Uncertifiable pairs (touching
/// at a vertex) end along a few paths whose length grows linearly with
/// the depth, and MAX_BOXES still bounds the whole search.
const MAX_DEPTH: usize = 24;

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

/// Exact lower bound on the distance between two UNTRIMMED analytic
/// surfaces, when a closed form exists (PK-329). `None` means "no
/// closed form for this class pair" and the pair falls through to the
/// sampled machinery unchanged. A full-surface minimum never exceeds
/// any trimmed patch's minimum, so `Some(lb)` with `lb > clearance`
/// certifies the trimmed pair apart in O(1) — the offset pipeline's
/// dominant verification class (inner/outer skins of one source face:
/// concentric spheres, parallel cylinders, coaxial cones a thickness
/// apart) burned the subdivision seeder's whole budget per pair
/// (measured 250 ms on the sphere-capped rod's cap pair) proving what
/// these formulas state directly.
///
/// Safety rules: closed forms over bounded geometry (sphere pairs,
/// point-to-cone) are exact in every configuration. Certificates that
/// rely on parallelism of INFINITE surfaces (plane×plane,
/// plane×cylinder, cylinder×cylinder, cone×cone, torus×torus) demand
/// EXACT parallelism — cross products that are identically 0.0, which
/// bitwise-copied axes from the offset pipeline satisfy — because an
/// epsilon-parallel infinite pair diverges without bound and any
/// tolerance gate would eventually mis-certify. Translation residuals
/// off the shared axis are handled by Lipschitz penalty (min distance
/// moves by at most the perpendicular offset).
pub(super) fn analytic_separation_lower_bound(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
) -> Option<f64> {
    use SurfaceGeometry as SG;
    // Distance between two circles/spheres of radii r1, r2 whose
    // centers sit `d` apart (also the cross-section form for parallel
    // cylinders and concentric tori): separated, nested, or touching.
    let radial_gap = |d: f64, r1: f64, r2: f64| -> f64 {
        if d >= r1 + r2 {
            d - r1 - r2
        } else if d + r1 <= r2 {
            r2 - d - r1
        } else if d + r2 <= r1 {
            r1 - d - r2
        } else {
            0.0
        }
    };
    // Min of |x - rc| over x in [lo, hi] (point-to-cylinder radial gap).
    let interval_abs_gap = |lo: f64, hi: f64, rc: f64| -> f64 {
        if rc < lo {
            lo - rc
        } else if rc > hi {
            rc - hi
        } else {
            0.0
        }
    };
    let point_line_distance = |p: &Point3, o: &Point3, axis: &crate::geometry::Vec3| -> f64 {
        let d = *p - *o;
        (d - *axis * d.dot(axis)).length()
    };
    match (a, b) {
        (
            SG::Sphere {
                center: c1,
                radius: r1,
                ..
            },
            SG::Sphere {
                center: c2,
                radius: r2,
                ..
            },
        ) => Some(radial_gap(c1.distance_to(c2), *r1, *r2)),
        (SG::Sphere { center, radius, .. }, SG::Plane { origin, normal, .. })
        | (SG::Plane { origin, normal, .. }, SG::Sphere { center, radius, .. }) => {
            Some((*center - *origin).dot(normal).abs() - radius)
        }
        (
            SG::Sphere { center, radius, .. },
            SG::Cylinder {
                origin,
                axis,
                radius: rc,
            },
        )
        | (
            SG::Cylinder {
                origin,
                axis,
                radius: rc,
            },
            SG::Sphere { center, radius, .. },
        ) => {
            let d_axis = point_line_distance(center, origin, axis);
            Some(interval_abs_gap(
                (d_axis - radius).max(0.0),
                d_axis + radius,
                *rc,
            ))
        }
        (
            SG::Sphere { center, radius, .. },
            SG::Cone {
                apex,
                axis,
                half_angle,
            },
        )
        | (
            SG::Cone {
                apex,
                axis,
                half_angle,
            },
            SG::Sphere { center, radius, .. },
        ) => {
            // The kernel cone evaluates BOTH nappes (negative v mirrors
            // through the apex), so bound against the union of the two
            // full profile generator lines — a superset of the double
            // cone's profile, hence a safe lower bound.
            let dv = *center - *apex;
            let zeta = dv.dot(axis);
            let rho = (dv - *axis * zeta).length();
            let (sin_a, cos_a) = half_angle.psin_cos();
            let d1 = (rho * cos_a - zeta * sin_a).abs();
            let d2 = (rho * cos_a + zeta * sin_a).abs();
            Some(d1.min(d2) - radius)
        }
        (
            SG::Plane {
                origin: o1,
                normal: n1,
                ..
            },
            SG::Plane {
                origin: o2,
                normal: n2,
                ..
            },
        ) => {
            if n1.cross(n2).length() == 0.0 {
                Some((*o2 - *o1).dot(n1).abs())
            } else {
                Some(0.0)
            }
        }
        (
            SG::Plane {
                origin: po, normal, ..
            },
            SG::Cylinder {
                origin,
                axis,
                radius,
            },
        )
        | (
            SG::Cylinder {
                origin,
                axis,
                radius,
            },
            SG::Plane {
                origin: po, normal, ..
            },
        ) => {
            if axis.dot(normal) == 0.0 {
                Some((*origin - *po).dot(normal).abs() - radius)
            } else {
                Some(0.0)
            }
        }
        (
            SG::Cylinder {
                origin: o1,
                axis: a1,
                radius: r1,
            },
            SG::Cylinder {
                origin: o2,
                axis: a2,
                radius: r2,
            },
        ) => {
            if a1.cross(a2).length() == 0.0 {
                Some(radial_gap(point_line_distance(o2, o1, a1), *r1, *r2))
            } else {
                // Skew or intersecting axes: min surface distance is at
                // least the line-line distance minus both radii.
                let n = a1.cross(a2);
                let d_lines = (*o2 - *o1).dot(&n).abs() / n.length();
                Some(d_lines - r1 - r2)
            }
        }
        // Cone×cone has NO untrimmed certificate: the kernel cone is a
        // double cone, and two coaxial double cones offset along their
        // axis genuinely intersect (opposite nappes cross in a circle).
        // The trim-window-aware bound lives in
        // `coaxial_cone_same_nappe_lower_bound`.
        (
            SG::Torus {
                center: c1,
                axis: x1,
                major_radius: mr1,
                minor_radius: tr1,
            },
            SG::Torus {
                center: c2,
                axis: x2,
                major_radius: mr2,
                minor_radius: tr2,
            },
        ) => {
            // Exactly coaxial tori: the profile cross-section is two
            // circles centered (R, ±Δaxial/…) — for centers offset only
            // along the shared axis the tube-center circles sit
            // sqrt(ΔR² + Δz²) apart in profile; radial center residual
            // subtracts via the Lipschitz bound.
            if x1.cross(x2).length() == 0.0 {
                let dv = *c2 - *c1;
                let along = dv.dot(x1);
                let radial = (dv - *x1 * along).length();
                let profile_d = ((mr1 - mr2) * (mr1 - mr2) + along * along).sqrt();
                Some(radial_gap(profile_d, *tr1, *tr2) - radial)
            } else {
                None
            }
        }
        _ => None,
    }
}

/// Window-aware coaxial cone×cone separation lower bound (PK-329).
///
/// The kernel cone is a double cone (negative v mirrors through the
/// apex), so two coaxial same-angle cones offset along their axis have
/// untrimmed distance 0 — opposite nappes cross. But when both TRIMMED
/// v-windows sit on matching nappes (both non-negative or both
/// non-positive), the trimmed patches lie on parallel half-lines in
/// profile and their distance is at least `|Δaxial|·sin α`, minus the
/// radial apex residual (translation Lipschitz bound). Exactness rules
/// match [`analytic_separation_lower_bound`]: axis parallelism and
/// angle equality must be EXACT (bitwise-copied offset-pipeline
/// geometry), anything else returns `None`.
pub(super) fn coaxial_cone_same_nappe_lower_bound(
    a: &SurfaceGeometry,
    va: (f64, f64),
    b: &SurfaceGeometry,
    vb: (f64, f64),
) -> Option<f64> {
    use SurfaceGeometry as SG;
    let (
        SG::Cone {
            apex: p1,
            axis: x1,
            half_angle: h1,
        },
        SG::Cone {
            apex: p2,
            axis: x2,
            half_angle: h2,
        },
    ) = (a, b)
    else {
        return None;
    };
    if !(x1.cross(x2).length() == 0.0 && x1.dot(x2) > 0.0 && h1 == h2) {
        return None;
    }
    let same_positive_nappe = va.0 >= 0.0 && vb.0 >= 0.0;
    let same_negative_nappe = va.1 <= 0.0 && vb.1 <= 0.0;
    if !(same_positive_nappe || same_negative_nappe) {
        return None;
    }
    let dv = *p2 - *p1;
    let along = dv.dot(x1);
    let radial = (dv - *x1 * along).length();
    Some(along.abs() * h1.psin() - radial)
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
