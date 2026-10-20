//! A body's half-space cutter: the solid `{ x : n · (x - o) <= 0 }`
//! clipped to a box around the body, one face exactly on the plane. A
//! prism of a square in the plane, extruded along `-n`, so its cap IS
//! the plane for any normal. Sized from the faces' conservative boxes
//! (a cylinder's rim has one vertex: a vertex box misses most of it).

use popkern_core::geometry::analytic_surfaces::orthonormal_basis;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::error::BoolError;

/// The part of the half-space below the plane through `origin` with
/// normal `normal` (any length) that covers `body`: `None` when the
/// body lies wholly above the plane (nothing of it is in the half), or
/// for a degenerate normal or an empty body.
///
/// # Errors
///
/// [`BoolError::CoreError`] when the prism cannot be built.
pub(crate) fn half_space(
    body: &BRep,
    origin: Point3,
    normal: Vec3,
    tolerance: Tolerance,
) -> Result<Option<BRep>, BoolError> {
    let Some(n) = normal.normalized() else {
        return Ok(None);
    };
    let boxes = crate::intersect::compute_face_aabbs(body);
    if boxes.is_empty() {
        return Ok(None);
    }
    let (mut lo, mut hi) = (
        Point3::new(f64::MAX, f64::MAX, f64::MAX),
        Point3::new(f64::MIN, f64::MIN, f64::MIN),
    );
    for (_, b) in &boxes {
        lo = Point3::new(lo.x.min(b.min.x), lo.y.min(b.min.y), lo.z.min(b.min.z));
        hi = Point3::new(hi.x.max(b.max.x), hi.y.max(b.max.y), hi.z.max(b.max.z));
    }
    // The prism's profile is counter-clockwise about its extrusion
    // direction -n: u x v = -n.
    let d = n * -1.0;
    let (u, _) = orthonormal_basis(d);
    let v = d.cross(&u);
    let corners = [
        (lo.x, lo.y, lo.z),
        (hi.x, lo.y, lo.z),
        (lo.x, hi.y, lo.z),
        (hi.x, hi.y, lo.z),
        (lo.x, lo.y, hi.z),
        (hi.x, lo.y, hi.z),
        (lo.x, hi.y, hi.z),
        (hi.x, hi.y, hi.z),
    ];
    let (mut ua, mut ub, mut va, mut vb, mut below) =
        (f64::MAX, f64::MIN, f64::MAX, f64::MIN, f64::MIN);
    for (x, y, z) in corners {
        let r = Point3::new(x, y, z) - origin;
        ua = ua.min(r.dot(&u));
        ub = ub.max(r.dot(&u));
        va = va.min(r.dot(&v));
        vb = vb.max(r.dot(&v));
        below = below.max(r.dot(&d));
    }
    if below <= tolerance.linear() {
        return Ok(None);
    }
    // A margin around the body (the box faces must stay clear of it).
    let size = (ub - ua).max(vb - va).max(below);
    let margin = 0.05 * size + 10.0 * tolerance.linear();
    let at = |a: f64, b: f64| origin + u * a + v * b;
    let profile = [
        at(ua - margin, va - margin),
        at(ub + margin, va - margin),
        at(ub + margin, vb + margin),
        at(ua - margin, vb + margin),
    ];
    BRep::prism(&profile, d, below + margin, tolerance)
        .map(Some)
        .map_err(BoolError::CoreError)
}
