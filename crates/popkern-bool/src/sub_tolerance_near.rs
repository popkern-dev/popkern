//! Which faces of a body can lie within the tolerance of a point: the
//! prefilter of the thinness test (PK-1161).
//!
//! A face lies on its surface, so the distance from a point to the
//! untrimmed surface is a lower bound of its distance to the face. For
//! a plane, cylinder, cone, sphere or torus that distance is closed
//! form, exact, and cheaper than the face's box; a face on any other
//! surface is tested against its box, built when first asked for (a
//! result's boxes cost as much as the rest of its thinness test: 0.07
//! of 0.1 ms on a filleted cube less a pocket, and its test ends at its
//! first face's interior point as often as not).

use std::cell::RefCell;
use std::collections::HashMap;

use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, FaceHandle};

/// The faces of `brep` that may lie near a point, by their surfaces
/// and, for a surface with no closed form, their boxes.
pub(super) struct NearFaces<'a> {
    brep: &'a BRep,
    boxes: RefCell<HashMap<FaceHandle, Option<Aabb>>>,
}

impl<'a> NearFaces<'a> {
    /// The prefilter over `brep`, with whatever face boxes are already
    /// built (`compute_face_aabbs`).
    pub(super) fn new(brep: &'a BRep, known: &[(FaceHandle, Aabb)]) -> Self {
        let boxes = known.iter().map(|(f, b)| (*f, Some(b.clone()))).collect();
        Self {
            brep,
            boxes: RefCell::new(boxes),
        }
    }

    /// Whether face `g` may lie within `reach` of `p`: false only when
    /// its surface, or its box, is farther. A face with no surface is
    /// never near.
    pub(super) fn may_reach(&self, g: FaceHandle, p: &Point3, reach: f64) -> bool {
        let Some(surface) = self.brep.faces.get(g).and_then(|f| f.surface.as_ref()) else {
            return false;
        };
        if let Some(d) = surface_gap(&surface.geometry, p) {
            return d <= reach;
        }
        let mut boxes = self.boxes.borrow_mut();
        let b = boxes
            .entry(g)
            .or_insert_with(|| crate::face_box::face_aabb(self.brep, g));
        b.as_ref().is_some_and(|b| box_gap(b, p) <= reach)
    }
}

/// The distance from `p` to box `b` (zero inside it).
pub(super) fn box_gap(b: &Aabb, p: &Point3) -> f64 {
    let gap = |lo: f64, hi: f64, v: f64| (lo - v).max(v - hi).max(0.0);
    let (dx, dy, dz) = (
        gap(b.min.x, b.max.x, p.x),
        gap(b.min.y, b.max.y, p.y),
        gap(b.min.z, b.max.z, p.z),
    );
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// The distance from `p` to the whole of an analytic surface, or `None`
/// for a surface with no closed form here. A surface of revolution's
/// nearest point lies in the meridian plane through `p`, so each is a
/// distance in that plane: `(ρ, z)` off the axis and along it.
pub(super) fn surface_gap(surface: &SurfaceGeometry, p: &Point3) -> Option<f64> {
    // `(ρ, z)` of `p` about the axis through `o` along `axis`.
    let meridian = |o: &Point3, axis: &Vec3| {
        let a = axis.normalized()?;
        let d = *p - *o;
        let z = d.dot(&a);
        Some(((d - a * z).length(), z))
    };
    match surface {
        SurfaceGeometry::Plane { origin, normal, .. } => {
            Some((*p - *origin).dot(&normal.normalized()?).abs())
        }
        SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        } => {
            let (rho, _) = meridian(origin, axis)?;
            Some((rho - radius).abs())
        }
        SurfaceGeometry::Sphere { center, radius, .. } => {
            Some((p.distance_to(center) - radius).abs())
        }
        SurfaceGeometry::Torus {
            center,
            axis,
            major_radius,
            minor_radius,
        } => {
            // The tube's circles in the meridian plane, either side of
            // the axis (a spindle torus's far circle can be nearer).
            use popkern_core::geometry::PortableMath;
            let (rho, z) = meridian(center, axis)?;
            let to = |c: f64| ((rho - c).phypot(z) - minor_radius).abs();
            Some(to(*major_radius).min(to(-major_radius)))
        }
        SurfaceGeometry::Cone {
            apex,
            axis,
            half_angle,
        } => {
            // Both generator lines through the apex in the meridian
            // plane (both nappes: a lower bound for either).
            use popkern_core::geometry::PortableMath;
            let (rho, z) = meridian(apex, axis)?;
            let (s, c) = (half_angle.psin(), half_angle.pcos());
            Some((rho * c - z * s).abs().min((rho * c + z * s).abs()))
        }
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The closed forms against a dense sampling of each surface: never
    /// over the sampled distance, and within its resolution of it.
    #[test]
    fn surface_gaps_are_the_distance_to_the_surface() {
        use popkern_core::geometry::surface::SurfaceEval;
        let axis = Vec3::new(0.3, -0.5, 0.8).normalized().expect("axis");
        let o = Point3::new(0.4, -1.0, 2.0);
        let surfaces = [
            SurfaceGeometry::Cylinder {
                origin: o,
                axis,
                radius: 1.5,
            },
            SurfaceGeometry::Sphere {
                center: o,
                radius: 1.5,
                axis,
            },
            SurfaceGeometry::Torus {
                center: o,
                axis,
                major_radius: 2.0,
                minor_radius: 0.5,
            },
            SurfaceGeometry::Torus {
                center: o,
                axis,
                major_radius: 0.5,
                minor_radius: 1.0,
            },
            SurfaceGeometry::Cone {
                apex: o,
                axis,
                half_angle: 0.4,
            },
        ];
        let points = [
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(2.0, 1.0, 3.0),
            Point3::new(0.4, -1.0, 2.0),
            Point3::new(-1.0, -3.0, 4.0),
            Point3::new(1.1, -0.2, 1.4),
        ];
        for s in &surfaces {
            let ((u0, u1), (v0, v1)) = s.parameter_range();
            // Both nappes of the cone (its chart runs on past the apex).
            let (v0, v1) = if matches!(s, SurfaceGeometry::Cone { .. }) {
                (-6.0, 6.0)
            } else {
                (v0.max(-6.0), v1.min(6.0))
            };
            let n = 400;
            for p in &points {
                let gap = surface_gap(s, p).expect("closed form");
                let mut sampled = f64::INFINITY;
                for i in 0..=n {
                    for j in 0..=n {
                        let u = u0 + (u1 - u0) * i as f64 / n as f64;
                        let v = v0 + (v1 - v0) * j as f64 / n as f64;
                        sampled = sampled.min(s.point_at(u, v).distance_to(p));
                    }
                }
                assert!(gap <= sampled + 1e-9, "{s:?} {p:?}: {gap} over {sampled}");
                assert!(sampled - gap < 0.05, "{s:?} {p:?}: {gap} under {sampled}");
            }
        }
    }
}
