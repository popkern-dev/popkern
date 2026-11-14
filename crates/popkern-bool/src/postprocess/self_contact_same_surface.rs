//! Two faces of one shell on one surface, cleared by their boundaries
//! (PK-1155).
//!
//! Faces of a valid shell share no interior point, so two faces lying
//! on the same surface meet only where their boundaries do: if a point
//! of one lies within `d` of the other, the shortest path on the surface
//! between them leaves the first face and enters the second, crossing
//! both boundaries, and on a surface of reach `ρ` that path is no longer
//! than `2ρ asin(d / 2ρ)`, within a part in 10⁵ of `d` when `ρ` exceeds
//! a hundred times `d`. So pieces of one torus or one cylinder sharing
//! no vertex, which no bound on a box clears (each box lies on the
//! other's surface), are cleared when every edge of one is farther than
//! the band, so stretched, from every edge of the other: each pair of
//! edges' control boxes apart, the larger halved while they are not.
//! Measured: a torus cut into quadrants and turned through thirty
//! orientations spent 4.1 s searching 200 such pairs.

use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{BRep, EdgeHandle};

/// Halvings allowed across a pair of edges before it is kept.
const DEPTH: usize = 24;

/// The reach of `s`, the least distance from it to its medial axis:
/// infinite for a plane, the radius of a cylinder or sphere, the lesser
/// of a ring torus's tube radius and hole radius; none for a surface
/// with a point of zero reach (a cone's apex, a spindle torus) or none
/// known in closed form.
fn reach(s: &SurfaceGeometry) -> Option<f64> {
    match s {
        SurfaceGeometry::Plane { .. } => Some(f64::INFINITY),
        SurfaceGeometry::Cylinder { radius, .. } | SurfaceGeometry::Sphere { radius, .. } => {
            Some(radius.abs())
        }
        SurfaceGeometry::Torus {
            major_radius,
            minor_radius,
            ..
        } => {
            let (big, small) = (major_radius.abs(), minor_radius.abs());
            (big > small).then(|| small.min(big - small))
        }
        _ => None,
    }
}

/// Whether two faces lying on the identical surface `surface`, with
/// boundary edges `a` and `b`, provably lie farther than `tol` apart:
/// the surface's reach is over a hundred times `tol`, and every edge of
/// one lies farther than `tol`, stretched by a part in a thousand, from
/// every edge of the other.
pub(super) fn same_surface_clear(
    brep: &BRep,
    surface: &SurfaceGeometry,
    (a, b): (&[EdgeHandle], &[EdgeHandle]),
    tol: f64,
) -> bool {
    if !reach(surface).is_some_and(|r| r > 100.0 * tol) {
        return false;
    }
    let curves = |edges: &[EdgeHandle]| -> Option<Vec<NurbsCurve>> {
        edges
            .iter()
            .map(|&e| {
                let c = brep.edges.get(e)?.curve.as_ref()?;
                let n = c.geometry.to_nurbs_in_range((c.t_start, c.t_end));
                n.weights.iter().all(|&w| w > 0.0).then_some(n)
            })
            .collect()
    };
    let (Some(ca), Some(cb)) = (curves(a), curves(b)) else {
        return false;
    };
    let pad = tol * (1.0 + 1e-3);
    ca.iter()
        .all(|p| cb.iter().all(|q| curves_apart(p, q, pad, DEPTH)))
}

/// Whether the curves `p` and `q`, of positive weights (each in its
/// control hull), lie farther than `pad` apart: their control boxes are,
/// or, halving the one of larger box up to `depth` times in all, every
/// pair of pieces is.
pub(super) fn curves_apart(p: &NurbsCurve, q: &NurbsCurve, pad: f64, depth: usize) -> bool {
    let (Some(bp), Some(bq)) = (
        Aabb::from_points(&p.control_points),
        Aabb::from_points(&q.control_points),
    ) else {
        return false;
    };
    if crate::parallel::ff_parallel::ff_parallel_reach::box_gap(&bp, &bq) > pad {
        return true;
    }
    if depth == 0 {
        return false;
    }
    let size = |b: &Aabb| b.min.distance_to(&b.max);
    let (larger, other, swap) = if size(&bp) >= size(&bq) {
        (p, q, false)
    } else {
        (q, p, true)
    };
    let Some((left, right)) = halves(larger) else {
        return false;
    };
    [left, right].iter().all(|half| {
        if swap {
            curves_apart(other, half, pad, depth - 1)
        } else {
            curves_apart(half, other, pad, depth - 1)
        }
    })
}

/// The two halves of `c` at the middle of its parameter range.
fn halves(c: &NurbsCurve) -> Option<(NurbsCurve, NurbsCurve)> {
    let k = c.degree;
    let t0 = *c.knots.get(k)?;
    let t1 = *c.knots.get(c.knots.len().checked_sub(k + 1)?)?;
    c.split_at(0.5 * (t0 + t1)).ok()
}

#[cfg(test)]
#[path = "self_contact_same_surface_tests.rs"]
mod tests;
