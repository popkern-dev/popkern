//! The rim arc's path on its fillet strip between two setback corners
//! (PK-401).
//!
//! The setback corners are exact trim crossings, and on mixed-radius
//! vertices the two corners of one strip sit at DIFFERENT stations
//! along its edge (the corner with a larger neighbour is farther from
//! the vertex: on a face with apex angle `α` the station of the corner
//! is `(d_next + d_this·cos α) / sin α`, `d` the contact distances).
//! A straight uv line between such corners is skewed along the edge,
//! and where two adjacent arcs are skewed towards the vertex the cap
//! boundary turns FLAT or REFLEX at the corner between them (measured:
//! hexagon6_variable, radii 0.02/0.04/0.06/0.06/0.04/0.02, the corner
//! between the two 0.02 strips turned −25.9°; a 6:1 ratio −103°). A
//! tensor-product cap corner cannot be built on such a turn (PK-399
//! served it through the χ̃-calibrated containment).
//!
//! The path keeps the exact corners and makes both ENDS pure
//! cross-section cuts: the axial uv coordinate follows a C² quintic
//! smoothstep between the two corner stations while the cross
//! coordinate stays linear, so the rim leaves and arrives
//! perpendicular to the edge with the pure section's tangent AND
//! curvature. Every cap corner then turns by the face's apex angle,
//! exactly as on a uniform layout, whatever the radii (measured: +46.0°
//! at all six hexagon corners, +54.6° at all five pentagon corners).
//! A layout with equal stations reproduces the straight line bit for
//! bit, so uniform fixtures are unchanged.

use popkern_core::geometry::curve::{CurveEval, closest_point_on_curve};
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{Point3, SurfaceGeometry, Vec3};

use crate::fillet_surface::FilletSurfaceResult;

/// The strip's parameter axis that runs along the edge.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Axis {
    /// `u` runs along the edge.
    U,
    /// `v` runs along the edge.
    V,
}

/// Minimum `|cos|` between a partial and the contact tangent for that
/// partial to count as axial (18°).
const AXIAL_COS_MIN: f64 = 0.95;
/// Maximum `|cos|` between the OTHER partial and the contact tangent
/// (the cross partial must be at least 72° from the edge).
const CROSS_COS_MAX: f64 = 0.3;

/// One rim arc's uv path from `corner_prev` (t = 0) to `corner_next`
/// (t = 1) on fillet `fillet_index`.
#[derive(Clone, Copy, Debug)]
pub(super) struct SectionPath {
    /// Index into the fillet surfaces.
    pub fillet_index: usize,
    /// uv of `corner_prev`.
    pub start: (f64, f64),
    /// uv of `corner_next`.
    pub end: (f64, f64),
    /// The axial axis when it could be resolved at both corners; the
    /// path is the straight uv line otherwise (recorded fallback).
    pub axial: Option<Axis>,
}

impl SectionPath {
    /// Resolve the path between two projected corners: the axial axis
    /// is the surface partial parallel to the contact tangent at BOTH
    /// corners (the other partial at least 72° from it).
    pub(super) fn resolve(
        fillet: &FilletSurfaceResult,
        fillet_index: usize,
        start: ((f64, f64), Point3),
        end: ((f64, f64), Point3),
    ) -> Self {
        let axis_at = |(u, v): (f64, f64), corner: Point3| -> Option<Axis> {
            let tangent = contact_tangent(fillet, corner)?;
            axial_axis(&fillet.surface, u, v, tangent)
        };
        let axial = match (axis_at(start.0, start.1), axis_at(end.0, end.1)) {
            (Some(a), Some(b)) if a == b => Some(a),
            _ => None,
        };
        Self {
            fillet_index,
            start: start.0,
            end: end.0,
            axial,
        }
    }

    /// uv at path parameter `t ∈ [0, 1]`.
    pub(super) fn uv(&self, t: f64) -> (f64, f64) {
        let (gu, gv) = match self.axial {
            Some(Axis::U) => (smoothstep_quintic(t), t),
            Some(Axis::V) => (t, smoothstep_quintic(t)),
            None => (t, t),
        };
        (
            self.start.0 + (self.end.0 - self.start.0) * gu,
            self.start.1 + (self.end.1 - self.start.1) * gv,
        )
    }

    /// Axial displacement between the two corners in parameter units
    /// (zero on a uniform layout; `None` when the axis is unresolved).
    pub(super) fn axial_displacement(&self) -> Option<f64> {
        match self.axial {
            Some(Axis::U) => Some(self.end.0 - self.start.0),
            Some(Axis::V) => Some(self.end.1 - self.start.1),
            None => None,
        }
    }
}

/// `6t⁵ − 15t⁴ + 10t³`: zero first AND second derivative at both ends.
fn smoothstep_quintic(t: f64) -> f64 {
    t * t * t * (t * (6.0 * t - 15.0) + 10.0)
}

/// Unit tangent of the contact curve nearest `corner` (whichever of
/// the strip's two trims passes closer).
fn contact_tangent(fillet: &FilletSurfaceResult, corner: Point3) -> Option<Vec3> {
    let mut best: Option<(f64, Vec3)> = None;
    for trim in [&fillet.trim_curve_1, &fillet.trim_curve_2] {
        let Some(t) = closest_point_on_curve(trim, &corner) else {
            continue;
        };
        let d = trim.point_at(t).distance_to(&corner);
        if best.is_none_or(|(bd, _)| d < bd) {
            best = Some((d, trim.derivative_at(t, 1)));
        }
    }
    best.and_then(|(_, tangent)| tangent.normalized())
}

/// Which partial at `(u, v)` runs along `tangent`.
fn axial_axis(surface: &SurfaceGeometry, u: f64, v: f64, tangent: Vec3) -> Option<Axis> {
    let cos_u = surface.partial_u(u, v).normalized()?.dot(&tangent).abs();
    let cos_v = surface.partial_v(u, v).normalized()?.dot(&tangent).abs();
    if cos_u >= AXIAL_COS_MIN && cos_v <= CROSS_COS_MAX {
        Some(Axis::U)
    } else if cos_v >= AXIAL_COS_MIN && cos_u <= CROSS_COS_MAX {
        Some(Axis::V)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn path(axial: Option<Axis>) -> SectionPath {
        SectionPath {
            fillet_index: 0,
            start: (0.2, 1.0),
            end: (1.4, 3.0),
            axial,
        }
    }

    #[test]
    fn smoothstep_is_c2_flat_at_both_ends() {
        let h = 1.0e-4;
        for t in [0.0, 1.0] {
            let tt = |s: f64| smoothstep_quintic(s.clamp(0.0, 1.0));
            let d1 = (tt(t + h) - tt(t - h)) / (2.0 * h);
            let d2 = (tt(t + h) - 2.0 * tt(t) + tt(t - h)) / (h * h);
            // One-sided at the ends: a flat quintic still reads ~0.
            assert!(d1.abs() < 1.0e-7, "d1 {d1:e} at {t}");
            assert!(d2.abs() < 1.0e-3, "d2 {d2:e} at {t}");
        }
        assert_eq!(smoothstep_quintic(0.0), 0.0);
        assert_eq!(smoothstep_quintic(1.0), 1.0);
        assert_eq!(smoothstep_quintic(0.5), 0.5);
    }

    #[test]
    fn ends_are_cross_sectional_and_uniform_layouts_are_the_straight_line() {
        let p = path(Some(Axis::V));
        let h = 1.0e-5;
        let (u0, v0) = p.uv(0.0);
        let (u1, v1) = p.uv(h);
        assert_eq!((u0, v0), p.start);
        assert_eq!(p.uv(1.0), p.end);
        // Axial (v) velocity vanishes at the end; cross (u) is linear.
        assert!((v1 - v0).abs() < 1.0e-12, "axial creep {:e}", v1 - v0);
        assert!(((u1 - u0) / h - (p.end.0 - p.start.0)).abs() < 1.0e-9);
        assert_eq!(p.axial_displacement(), Some(2.0));
        // Equal stations: bit-identical to the straight line.
        let uniform = SectionPath {
            end: (1.4, 1.0),
            ..p
        };
        let straight = SectionPath {
            axial: None,
            ..uniform
        };
        for i in 0..=16 {
            let t = i as f64 / 16.0;
            assert_eq!(uniform.uv(t), straight.uv(t));
        }
        assert_eq!(uniform.axial_displacement(), Some(0.0));
        assert_eq!(straight.axial_displacement(), None);
    }

    #[test]
    fn axial_axis_reads_a_cylinder_and_rejects_a_diagonal() {
        let cyl = SurfaceGeometry::Cylinder {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 0.5,
        };
        let z = Vec3::new(0.0, 0.0, 1.0);
        assert_eq!(axial_axis(&cyl, 0.3, 0.2, z), Some(Axis::V));
        let diagonal = (z + cyl.partial_u(0.3, 0.2).normalized().expect("du"))
            .normalized()
            .expect("diag");
        assert_eq!(axial_axis(&cyl, 0.3, 0.2, diagonal), None);
    }
}
