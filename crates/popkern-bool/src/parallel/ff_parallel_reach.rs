//! Whether a face pair can meet at all, before it is solved, and
//! whether its intersection curve can reach a face, before it is
//! clipped (PK-1166).
//!
//! A record lies on both faces, and a trim residue within the clip's
//! boundary band (ten tolerances) of a trim, so a curve farther than
//! that from either face's box yields nothing. Its clip is not free:
//! the boss of a disc union, its bottom plane against the disc's wall,
//! carried the wall's whole circle 27 units from the boss, and the
//! clip that kept none of it, with the retry an empty clip makes,
//! took 2.1 ms (the four boss walls' lines another 0.45 ms each),
//! more than the records the stage keeps.
//!
//! The test is a proof, never a sampling: a NURBS curve of positive
//! weights lies in the hull of its control points, so a piece whose
//! control box is farther than the pad from the face's box is clear,
//! and a piece that is not is split in half and tested again, to a
//! fixed depth (a curve passing near stays undecided and is clipped).
//! Only a face whose box contains it is trusted: a plane, cylinder,
//! cone, sphere or torus (`face_box::face_aabb` boxes their boundary
//! curves exactly and their surfaces' extremes within the trim); a
//! NURBS face's box may be its grid, which can fall short of it.
//!
//! Before the solve: a face lies in its box, so a box farther than the
//! pad from the other face's whole surface leaves the pair nothing to
//! meet in. The same boss's bottom face, within 2.9 of the disc's axis,
//! against the disc's wall of radius 30: its solve, the circle and the
//! circle's co-parametrized pcurve fitted to a fraction of the
//! tolerance, took 1.3 ms. The bound is closed form for a plane (the
//! box's corners on one side), a sphere (the box inside it or apart
//! from it) and a cylinder (its corners' largest distance from the
//! axis under the radius, or its middle's less half its diagonal over
//! it; the distance from an axis is convex and changes no faster than
//! the point).

use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::surface::SurfaceGeometry;

/// The pad beyond a face's box, in linear tolerances: ten times the
/// clip's boundary band.
pub(super) const REACH_TOLERANCES: f64 = 100.0;

/// Halvings before a piece near the box is left undecided.
const DEPTH: usize = 6;

/// Whether a face on `surface` is known to lie in its box.
fn box_contains_face(surface: &SurfaceGeometry) -> bool {
    matches!(
        surface,
        SurfaceGeometry::Plane { .. }
            | SurfaceGeometry::Cylinder { .. }
            | SurfaceGeometry::Cone { .. }
            | SurfaceGeometry::Sphere { .. }
            | SurfaceGeometry::Torus { .. }
    )
}

/// Whether `curve` lies farther than `pad` from the box `face_box` of a
/// face on `surface`: true only when proved (false for a surface whose
/// box is not known to contain its face, or a weight not positive).
pub(crate) fn curve_clear_of_face(
    curve: &NurbsCurve,
    surface: &SurfaceGeometry,
    face_box: &Aabb,
    pad: f64,
) -> bool {
    box_contains_face(surface)
        && curve.weights.iter().all(|&w| w > 0.0)
        && !touches(curve, face_box, pad)
        && clear(curve, face_box, pad, DEPTH)
}

/// Whether a point of `curve` (its ends, its middle) lies within `pad`
/// of the box: then it is not clear, and the halving is not run (every
/// kept record's curve lies on its faces).
fn touches(curve: &NurbsCurve, face_box: &Aabb, pad: f64) -> bool {
    use popkern_core::geometry::curve::CurveEval;
    let (t0, t1) = curve.parameter_range();
    [t0, 0.5 * (t0 + t1), t1].into_iter().any(|t| {
        let p = curve.point_at(t);
        let point = Aabb { min: p, max: p };
        box_gap(&point, face_box) <= pad
    })
}

/// Whether two faces, each on its surface with its box, provably lie
/// farther than `pad` apart: one face's box, known to contain it, is
/// farther than that from the other's whole surface.
pub(crate) fn faces_apart(
    (surface_a, box_a): (&SurfaceGeometry, &Aabb),
    (surface_b, box_b): (&SurfaceGeometry, &Aabb),
    pad: f64,
) -> bool {
    (box_contains_face(surface_a) && box_clear_of_surface(box_a, surface_b, pad))
        || (box_contains_face(surface_b) && box_clear_of_surface(box_b, surface_a, pad))
}

/// Whether box `b` lies farther than `pad` from a face on `surface`
/// whose box is `face_box`: from the whole surface in closed form, or
/// from the face's box when that box is known to hold the face.
pub(crate) fn box_clear_of_face(
    b: &Aabb,
    surface: &SurfaceGeometry,
    face_box: &Aabb,
    pad: f64,
) -> bool {
    box_clear_of_surface(b, surface, pad)
        || (box_contains_face(surface) && box_gap(b, face_box) > pad)
}

/// Whether box `b` lies farther than `pad` from the whole of `surface`:
/// true only when proved (closed form for a plane, sphere, cylinder,
/// cone or torus).
fn box_clear_of_surface(b: &Aabb, surface: &SurfaceGeometry, pad: f64) -> bool {
    use popkern_core::geometry::Point3;
    let corners: [Point3; 8] = std::array::from_fn(|i| {
        Point3::new(
            if i & 1 == 0 { b.min.x } else { b.max.x },
            if i & 2 == 0 { b.min.y } else { b.max.y },
            if i & 4 == 0 { b.min.z } else { b.max.z },
        )
    });
    let middle = Point3::new(
        0.5 * (b.min.x + b.max.x),
        0.5 * (b.min.y + b.max.y),
        0.5 * (b.min.z + b.max.z),
    );
    let half_diagonal = 0.5 * b.min.distance_to(&b.max);
    match surface {
        SurfaceGeometry::Plane { origin, normal, .. } => {
            let Some(n) = normal.normalized() else {
                return false;
            };
            let side = corners.map(|c| (c - *origin).dot(&n));
            side.iter().all(|&d| d > pad) || side.iter().all(|&d| d < -pad)
        }
        SurfaceGeometry::Sphere { center, radius, .. } => {
            let gap = |lo: f64, hi: f64, v: f64| (lo - v).max(v - hi).max(0.0);
            let (dx, dy, dz) = (
                gap(b.min.x, b.max.x, center.x),
                gap(b.min.y, b.max.y, center.y),
                gap(b.min.z, b.max.z, center.z),
            );
            let nearest = (dx * dx + dy * dy + dz * dz).sqrt();
            let farthest = corners
                .iter()
                .map(|c| c.distance_to(center))
                .fold(0.0, f64::max);
            nearest > radius + pad || farthest < radius - pad
        }
        SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        } => {
            let Some(a) = axis.normalized() else {
                return false;
            };
            let off_axis = |p: &Point3| {
                let d = *p - *origin;
                (d - a * d.dot(&a)).length()
            };
            let farthest = corners.iter().map(off_axis).fold(0.0, f64::max);
            farthest < radius - pad || off_axis(&middle) - half_diagonal > radius + pad
        }
        // The distance to a surface is 1-Lipschitz: the box's middle
        // farther than its half-diagonal and the pad clears it.
        SurfaceGeometry::Cone {
            apex,
            axis,
            half_angle,
        } => {
            use popkern_core::geometry::PortableMath;
            let Some(a) = axis.normalized() else {
                return false;
            };
            // In the half-plane through the axis, (|h|, rho) lies
            // |rho cos - |h| sin| from the nearer nappe's generator, its
            // foot on the ray (rho sin + |h| cos >= 0).
            let v = middle - *apex;
            let h = v.dot(&a);
            let rho = (v - a * h).length();
            let d = (rho * half_angle.pcos() - h.abs() * half_angle.psin()).abs();
            d - half_diagonal > pad
        }
        SurfaceGeometry::Torus {
            center,
            axis,
            major_radius,
            minor_radius,
        } => {
            let Some(a) = axis.normalized() else {
                return false;
            };
            let v = middle - *center;
            let h = v.dot(&a);
            let rho = (v - a * h).length();
            let tube = ((rho - major_radius.abs()).powi(2) + h * h).sqrt();
            (tube - minor_radius.abs()).abs() - half_diagonal > pad
        }
        _ => false,
    }
}

fn clear(curve: &NurbsCurve, face_box: &Aabb, pad: f64, depth: usize) -> bool {
    clear_by(curve, depth, &|hull| box_gap(hull, face_box) > pad)
}

/// Whether `curve`, of positive weights, lies farther than `pad` from a
/// face on `surface` with box `face_box`: each piece's control box,
/// halved up to `DEPTH` times, is ([`box_clear_of_face`]: the whole
/// surface in closed form, or the box when it is known to hold the
/// face).
pub(crate) fn curve_clear_of_surface_face(
    curve: &NurbsCurve,
    surface: &SurfaceGeometry,
    face_box: &Aabb,
    pad: f64,
) -> bool {
    curve.weights.iter().all(|&w| w > 0.0)
        && clear_by(curve, DEPTH + 2, &|hull| {
            box_clear_of_face(hull, surface, face_box, pad)
        })
}

/// Whether every piece of `curve`, halved up to `depth` times, has a
/// control box `apart` accepts.
fn clear_by(curve: &NurbsCurve, depth: usize, apart: &dyn Fn(&Aabb) -> bool) -> bool {
    let Some(hull) = Aabb::from_points(&curve.control_points) else {
        return false;
    };
    if apart(&hull) {
        return true;
    }
    if depth == 0 {
        return false;
    }
    let p = curve.degree;
    let (Some(&t0), Some(&t1)) = (
        curve.knots.get(p),
        curve
            .knots
            .len()
            .checked_sub(p + 1)
            .and_then(|i| curve.knots.get(i)),
    ) else {
        return false;
    };
    match curve.split_at(0.5 * (t0 + t1)) {
        Ok((left, right)) => {
            clear_by(&left, depth - 1, apart) && clear_by(&right, depth - 1, apart)
        }
        Err(_) => false,
    }
}

/// The distance between two boxes (zero when they overlap).
pub(crate) fn box_gap(a: &Aabb, b: &Aabb) -> f64 {
    let gap = |alo: f64, ahi: f64, blo: f64, bhi: f64| (blo - ahi).max(alo - bhi).max(0.0);
    let (dx, dy, dz) = (
        gap(a.min.x, a.max.x, b.min.x, b.max.x),
        gap(a.min.y, a.max.y, b.min.y, b.max.y),
        gap(a.min.z, a.max.z, b.min.z, b.max.z),
    );
    (dx * dx + dy * dy + dz * dz).sqrt()
}

#[cfg(test)]
#[path = "ff_parallel_reach_tests.rs"]
mod tests;
