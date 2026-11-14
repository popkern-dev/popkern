//! The plane x torus clear of the self-contact search (PK-1175): a
//! planar face and a toroidal one whose normals are parallel nowhere
//! within the band, and whose boundaries keep apart from the other
//! face, provably lie farther than the band apart.

use popkern_core::geometry::Point3;
use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{BRep, FaceHandle};

use super::clear::{Side, edges_clear_of_plane, plane_of};

/// Whether a planar face (`plane`, its box and edges) and a toroidal
/// face (`other`) provably lie farther than `tol` apart (PK-1175). The
/// torus's normal `cos v e(u) + sin v a` is parallel to the plane's
/// unit normal `p` only where `e(u)` is `±` the direction of `p`'s part
/// across the axis (of length `ρ`, the sine of `p`'s angle to the
/// axis), at heights `h ± R ρ ± r` above the plane, `h` the centre's
/// height, `R` and `r` the radii (at `ρ = 0`, the circles at `h ± r`).
/// When none is within `tol` the faces' nearest points lie on a
/// boundary, and no edge of either may come within `tol` of the other.
pub(super) fn plane_torus_clear(
    brep: &BRep,
    (plane_face, plane, plane_box, plane_edges): Side,
    (torus_face, other, other_box, other_edges): Side,
    tol: f64,
) -> bool {
    let (
        Some((origin, p)),
        SurfaceGeometry::Torus {
            center,
            axis,
            major_radius,
            minor_radius,
        },
    ) = (plane_of(plane), other)
    else {
        return false;
    };
    let Some(a) = axis.normalized() else {
        return false;
    };
    let along = p.dot(&a);
    let rho = (p - a * along).length();
    let h = (*center - origin).dot(&p);
    let (big, small) = (major_radius.abs() * rho, minor_radius.abs());
    let normals_apart = [
        h + big + small,
        h + big - small,
        h - big + small,
        h - big - small,
    ]
    .iter()
    .all(|x| x.abs() > tol);
    if !normals_apart
        || !edges_clear_of_plane(brep, other_edges, (origin, p), (plane_face, plane_box), tol)
    {
        return false;
    }
    // The plane's edges against the torus face's own box: an edge lying
    // on the torus (the plane's cut along it) is near the surface
    // everywhere, and apart only from the face.
    let face_box = torus_face_box(brep, torus_face, other)
        .and_then(|b| meet(&b, other_box))
        .unwrap_or_else(|| other_box.clone());
    plane_edges.iter().all(|&e| {
        brep.edges
            .get(e)
            .and_then(|x| x.curve.as_ref())
            .is_some_and(|c| {
                let sub = c.geometry.to_nurbs_in_range((c.t_start, c.t_end));
                crate::parallel::ff_parallel::ff_parallel_reach::curve_clear_of_surface_face(
                    &sub, other, &face_box, tol,
                )
            })
    })
}

/// The common part of two boxes, `None` when they are apart.
fn meet(a: &Aabb, b: &Aabb) -> Option<Aabb> {
    let min = Point3::new(
        a.min.x.max(b.min.x),
        a.min.y.max(b.min.y),
        a.min.z.max(b.min.z),
    );
    let max = Point3::new(
        a.max.x.min(b.max.x),
        a.max.y.min(b.max.y),
        a.max.z.min(b.max.z),
    );
    (min.x <= max.x && min.y <= max.y && min.z <= max.z).then_some(Aabb { min, max })
}

/// A box holding the toroidal `face` (on `torus`), in closed form over
/// its boundary's uv box: each coordinate of
/// `C + (R + r cos v)(cos u p1 + sin u p2) + r sin v a` bounded by
/// interval products of the ranges of `cos v`, `sin v` and `cos u p1 +
/// sin u p2`'s component (an amplitude times a cosine) over the box.
fn torus_face_box(brep: &BRep, face: FaceHandle, torus: &SurfaceGeometry) -> Option<Aabb> {
    use popkern_core::geometry::PortableMath;
    let SurfaceGeometry::Torus {
        center,
        axis,
        major_radius,
        minor_radius,
    } = torus
    else {
        return None;
    };
    // The evaluation uses the axis as given: read only a unit one.
    if (axis.length() - 1.0).abs() > 1e-12 {
        return None;
    }
    let boxes = super::patches::face_pcurve_boxes(brep, face)?;
    let tau = std::f64::consts::TAU;
    let (u0, u1) = covering_arc(boxes.iter().map(|b| (b.0.0, b.1.0)), tau)?;
    let (v0, v1) = covering_arc(boxes.iter().map(|b| (b.0.1, b.1.1)), tau)?;
    let (p1, p2) = popkern_core::geometry::analytic_surfaces::orthonormal_basis(*axis);
    let half_pi = std::f64::consts::FRAC_PI_2;
    let cos_v = cos_range(v0, v1);
    let sin_v = cos_range(v0 - half_pi, v1 - half_pi);
    let ring = (
        major_radius + minor_radius * cos_v.0,
        major_radius + minor_radius * cos_v.1,
    );
    let ring = (ring.0.min(ring.1), ring.0.max(ring.1));
    let coordinate = |c: f64, e1: f64, e2: f64, ak: f64| {
        // cos u e1 + sin u e2 = A cos(u - phi).
        let (amp, phi) = (e1.phypot(e2), e2.patan2(e1));
        let (lo, hi) = cos_range(u0 - phi, u1 - phi);
        let e = (amp * lo, amp * hi);
        let products = [ring.0 * e.0, ring.0 * e.1, ring.1 * e.0, ring.1 * e.1];
        let tube = (minor_radius * ak * sin_v.0, minor_radius * ak * sin_v.1);
        let (plo, phi_) = products
            .iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(l, h), &x| {
                (l.min(x), h.max(x))
            });
        (c + plo + tube.0.min(tube.1), c + phi_ + tube.0.max(tube.1))
    };
    let x = coordinate(center.x, p1.x, p2.x, axis.x);
    let y = coordinate(center.y, p1.y, p2.y, axis.y);
    let z = coordinate(center.z, p1.z, p2.z, axis.z);
    // A rounding allowance on the closed form's evaluations.
    let pad = 1e-12
        * (center.x.abs()
            + center.y.abs()
            + center.z.abs()
            + major_radius.abs()
            + minor_radius.abs()
            + 1.0);
    Some(Aabb {
        min: Point3::new(x.0 - pad, y.0 - pad, z.0 - pad),
        max: Point3::new(x.1 + pad, y.1 + pad, z.1 + pad),
    })
}

/// The least arc of the circle of `period` covering every interval:
/// the circle less the largest gap between them, as `(start, end)` with
/// `end - start` under the period (the whole period when they cover it).
/// A face that is a disc in its chart spans along each direction what
/// its boundary does, a connected arc, and lies in this one; read off
/// its pcurves' boxes, a face beside the seam would otherwise span the
/// period (one edge at 0, another at the period).
fn covering_arc(intervals: impl Iterator<Item = (f64, f64)>, period: f64) -> Option<(f64, f64)> {
    let mut arcs: Vec<(f64, f64)> = Vec::new();
    for (a, b) in intervals {
        if !(a.is_finite() && b.is_finite()) {
            return None;
        }
        if b - a >= period {
            return Some((0.0, period));
        }
        let s = a.rem_euclid(period);
        arcs.push((s, s + (b - a)));
    }
    arcs.sort_by(|x, y| x.0.total_cmp(&y.0));
    let first = *arcs.first()?;
    // Merge in order; an arc reaching past the period wraps onto the
    // first ones.
    let mut merged: Vec<(f64, f64)> = vec![first];
    for &(s, e) in &arcs[1..] {
        let last = merged.last_mut()?;
        if s <= last.1 {
            last.1 = last.1.max(e);
        } else {
            merged.push((s, e));
        }
    }
    // The gaps between consecutive merged arcs, the last wrapping to
    // the first.
    let n = merged.len();
    let (mut best, mut at) = (f64::NEG_INFINITY, 0);
    for k in 0..n {
        let end = merged[k].1;
        let next = if k + 1 < n {
            merged[k + 1].0
        } else {
            merged[0].0 + period
        };
        if next - end > best {
            (best, at) = (next - end, k);
        }
    }
    if best <= 0.0 {
        return Some((0.0, period));
    }
    // The covering arc starts after the largest gap.
    let start = if at + 1 < n {
        merged[at + 1].0
    } else {
        merged[0].0
    };
    Some((start, start + period - best))
}

/// The least and greatest cosine over `[t0, t1]`: its ends, and 1 and
/// -1 where a multiple of the period or a half period lies inside.
fn cos_range(t0: f64, t1: f64) -> (f64, f64) {
    use popkern_core::geometry::PortableMath;
    let tau = std::f64::consts::TAU;
    let pi = std::f64::consts::PI;
    if !(t0.is_finite() && t1.is_finite()) || t1 - t0 >= tau {
        return (-1.0, 1.0);
    }
    let (a, b) = (t0.pcos(), t1.pcos());
    let (mut lo, mut hi) = (a.min(b), a.max(b));
    // The first multiple of 2 pi (and of 2 pi plus pi) at or past t0.
    if (t0 / tau).ceil() * tau <= t1 {
        hi = 1.0;
    }
    if ((t0 - pi) / tau).ceil() * tau + pi <= t1 {
        lo = -1.0;
    }
    (lo, hi)
}

#[cfg(test)]
#[path = "self_contact_torus_tests.rs"]
mod tests;
