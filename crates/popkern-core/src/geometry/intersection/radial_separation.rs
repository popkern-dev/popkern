//! No-crossing certificate for a NURBS wall against a cylinder
//! (PK-513).
//!
//! The thread tool's crest and root walls run coaxial with the bore,
//! and most of them never meet it: the crest sits outside the bore,
//! the root inside it. The marcher could only discover that by hunting
//! seeds over the walls' multi-turn charts (measured: 4.8 s per
//! non-crossing pair on an M6 tap, 0 curves found, two such pairs per
//! tap and nut). A rational surface with positive weights lies in the
//! convex hull of the control points of each knot span, which bounds
//! where it can be:
//!
//! - **Inside**: distance-to-axis is convex, so the whole surface is no
//!   farther from the axis than its farthest control point. A net
//!   strictly inside the radius cannot reach the cylinder.
//! - **Outside**: no global bound exists, because a coil's control net
//!   wraps the axis and its hull contains it. Per knot span, though,
//!   the local net occupies a narrow sector. The half-space beyond the
//!   cylinder's tangent plane in that sector's direction misses the
//!   solid cylinder entirely, so a span net inside that half-space
//!   cannot reach the cylinder.
//!
//! A control polygon bulges off its surface (a helix's cubic net sits
//! up to 1.3% outside it — measured: a tap root wall at true radius
//! 2.4415 carried control points to 2.4734 against a 2.4587 bore, so
//! neither test fired and the pair marched for 4.8 s). Knot refinement
//! pulls the net onto the surface quadratically, so the tests are
//! retried on up to [`REFINE_ROUNDS`] span-midpoint refinements.
//!
//! Both tests are exact sufficient conditions (no sampling). A wall
//! that passes neither falls through to the solvers unchanged.

use crate::geometry::nurbs_surface::NurbsSurface;
use crate::geometry::{Point3, Vec3};

/// Span-midpoint refinement rounds tried before falling through. Each
/// round quarters the polygon's bulge; three take a 1.3% bulge to
/// ~0.02%.
const REFINE_ROUNDS: usize = 3;

/// Whether the NURBS `wall` provably stays clear of the cylinder
/// (`origin`, unit `axis`, `radius`) by more than `gap` everywhere on
/// its knot domain.
pub(crate) fn wall_clear_of_cylinder(
    wall: &NurbsSurface,
    origin: Point3,
    axis: Vec3,
    radius: f64,
    gap: f64,
) -> bool {
    if !wall.weights.iter().all(|w| *w > 0.0)
        || wall.control_points.len() != wall.num_u * wall.num_v
        || wall.control_points.is_empty()
    {
        return false;
    }
    let mut net = std::borrow::Cow::Borrowed(wall);
    for round in 0..=REFINE_ROUNDS {
        if net_clear(&net, origin, axis, radius, gap) {
            return true;
        }
        if round == REFINE_ROUNDS {
            break;
        }
        match refined(&net) {
            Some(r) => net = std::borrow::Cow::Owned(r),
            None => break,
        }
    }
    false
}

/// Insert every live span's midpoint in each direction of degree > 1
/// (a degree-1 direction's net already lies on the surface).
fn refined(wall: &NurbsSurface) -> Option<NurbsSurface> {
    let mids = |knots: &[f64]| -> Vec<f64> {
        knots
            .windows(2)
            .filter(|w| w[1] > w[0])
            .map(|w| 0.5 * (w[0] + w[1]))
            .collect()
    };
    let mut out = wall.clone();
    if wall.degree_u > 1 {
        out = out.insert_knots_u(&mids(&wall.knots_u)).ok()?;
    }
    if wall.degree_v > 1 {
        out = out.insert_knots_v(&mids(&wall.knots_v)).ok()?;
    }
    Some(out)
}

/// The inside and outside certificates on one control net.
fn net_clear(wall: &NurbsSurface, origin: Point3, axis: Vec3, radius: f64, gap: f64) -> bool {
    let radial = |p: &Point3| -> Vec3 {
        let d = *p - origin;
        d - axis * d.dot(&axis)
    };
    let hull_max = wall
        .control_points
        .iter()
        .map(|p| radial(p).length())
        .fold(0.0_f64, f64::max);
    if hull_max < radius - gap {
        return true;
    }
    spans_outside(wall, &radial, radius + gap)
}

/// Every knot span's local control net lies beyond the cylinder's
/// tangent plane along that span's mean radial direction.
fn spans_outside(wall: &NurbsSurface, radial: &dyn Fn(&Point3) -> Vec3, bound: f64) -> bool {
    let (p, q) = (wall.degree_u, wall.degree_v);
    if wall.num_u <= p || wall.num_v <= q {
        return false;
    }
    let live = |knots: &[f64], i: usize| {
        knots
            .get(i + 1)
            .zip(knots.get(i))
            .is_some_and(|(b, a)| b > a)
    };
    for j in q..wall.num_v {
        if !live(&wall.knots_v, j) {
            continue;
        }
        for i in p..wall.num_u {
            if !live(&wall.knots_u, i) {
                continue;
            }
            let net = (j - q..=j)
                .flat_map(|jj| (i - p..=i).map(move |ii| ii * wall.num_v + jj))
                .map(|k| radial(&wall.control_points[k]));
            let vectors: Vec<Vec3> = net.collect();
            let mean = vectors
                .iter()
                .fold(Vec3::new(0.0, 0.0, 0.0), |acc, v| acc + *v)
                * (1.0 / vectors.len() as f64);
            let Some(dir) = mean.normalized() else {
                return false;
            };
            if vectors.iter().any(|v| v.dot(&dir) <= bound) {
                return false;
            }
        }
    }
    true
}

#[cfg(test)]
#[path = "radial_separation_tests.rs"]
mod tests;
