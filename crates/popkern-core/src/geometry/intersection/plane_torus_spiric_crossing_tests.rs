//! A plane tangent to a torus at its inner equator (a hyperbolic point),
//! opened from the tangency within the tolerance either way, crosses
//! itself there: its loops meet at the point halfway to the tangency.

use crate::geometry::curve::CurveEval;
use crate::geometry::intersection::plane_torus::try_plane_torus;
use crate::geometry::intersection::surface_surface::SurfaceSurfaceIntersection;
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use std::f64::consts::TAU;

fn torus_residual(p: &Point3) -> f64 {
    let rho = (p.x * p.x + p.y * p.y).sqrt();
    (((rho - 2.0) * (rho - 2.0) + p.z * p.z).sqrt() - 0.5).abs()
}

/// The loops of the plane `x = 1.5 + offset` (the inner equator of a
/// torus (2, 0.5) lies at x = 1.5) against the torus.
fn sections(offset: f64) -> Vec<crate::geometry::curve_on_surface::CurveOnSurface> {
    let torus = SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: 2.0,
        minor_radius: 0.5,
    };
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(1.5 + offset, 0.0, 0.0),
        normal: Vec3::new(1.0, 0.0, 0.0),
        u_axis: Vec3::new(0.0, 1.0, 0.0),
        v_axis: Vec3::new(0.0, 0.0, 1.0),
    };
    let tol = Tolerance::default();
    let out = try_plane_torus(
        &plane,
        (-6.0, 6.0),
        (-6.0, 6.0),
        &torus,
        (0.0, TAU),
        (0.0, TAU),
        &tol,
    )
    .expect("claimed");
    out.into_iter()
        .filter_map(|ix| match ix {
            SurfaceSurfaceIntersection::Curve(c) => Some(c),
            SurfaceSurfaceIntersection::CoincidentRegion { .. } => None,
        })
        .collect()
}

/// Failed before the crossing was spelled: opened outward (a waist) the
/// plane cut one loop, opened inward (a gap) two ovals passing beside
/// the point, neither through it.
#[test]
fn a_plane_opened_from_a_saddle_tangency_within_the_tolerance_crosses_at_it() {
    let linear = Tolerance::default().linear();
    for k in [0.5, -0.5] {
        let offset = k * linear;
        let curves = sections(offset);
        assert_eq!(
            curves.len(),
            2,
            "offset {k}: two loops meeting at the crossing"
        );
        let crossing = Point3::new(1.5 + 0.5 * offset, 0.0, 0.0);
        for c in &curves {
            let t = c.curve_3d.closest_point(&crossing).expect("a foot");
            let gap = c.curve_3d.point_at(t).distance_to(&crossing);
            assert!(
                gap < 0.1 * linear,
                "offset {k}: a loop passes the crossing {gap:.3e} away"
            );
            for i in 0..=512 {
                let p = c.curve_3d.point_at(i as f64 / 512.0);
                assert!(
                    torus_residual(&p) <= 0.5 * linear,
                    "offset {k}: off the torus"
                );
                assert!(
                    (p.x - 1.5 - offset).abs() <= 0.5 * linear,
                    "offset {k}: off the plane"
                );
            }
        }
    }
}

/// Opened past the tolerance the plane keeps its exact section: one loop
/// with a waist outward, two ovals inward.
#[test]
fn a_plane_opened_past_the_tolerance_keeps_its_exact_section() {
    let linear = Tolerance::default().linear();
    assert_eq!(sections(3.0 * linear).len(), 1);
    assert_eq!(sections(-3.0 * linear).len(), 2);
    for c in sections(3.0 * linear)
        .iter()
        .chain(sections(-3.0 * linear).iter())
    {
        let crossing = Point3::new(1.5, 0.0, 0.0);
        let t = c.curve_3d.closest_point(&crossing).expect("a foot");
        assert!(c.curve_3d.point_at(t).distance_to(&crossing) > 10.0 * linear);
    }
}
