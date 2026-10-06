use crate::geometry::curve::CurveEval;
use crate::geometry::intersection::surface_surface::SurfaceSurfaceIntersection;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

/// Worst `|cylinder(pcurve(t)) - curve(t)|` over a dense sweep of the
/// record's own parameter: the co-parameterization contract.
fn worst_coparameter_gap(cylinder: &SurfaceGeometry, plane: &SurfaceGeometry) -> f64 {
    let tol = Tolerance::default();
    let found = super::super::plane_cylinder::try_plane_cylinder(
        plane,
        (-5.0, 5.0),
        (-5.0, 5.0),
        cylinder,
        (0.0, std::f64::consts::TAU),
        (-5.0, 5.0),
        &tol,
    );
    let Some(found) = found else {
        panic!("plane-cylinder pair not recognised");
    };
    let Some(SurfaceSurfaceIntersection::Curve(c)) = found.first() else {
        panic!("no conic produced");
    };
    let (t0, t1) = c.curve_3d.parameter_range();
    let mut worst = 0.0f64;
    for i in 0..=2000 {
        let t = t0 + (t1 - t0) * i as f64 / 2000.0;
        let uv = c.curve_uv2.point_at(t);
        let on_cylinder = cylinder.point_at(uv.x, uv.y);
        worst = worst.max(on_cylinder.distance_to(&c.curve_3d.point_at(t)));
    }
    worst
}

fn cylinder(radius: f64) -> SurfaceGeometry {
    SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius,
    }
}

fn plane(normal: Vec3) -> SurfaceGeometry {
    let normal = normal.normalized().unwrap_or(Vec3::new(0.0, 0.0, 1.0));
    let (u_axis, v_axis) = crate::geometry::analytic_surfaces::orthonormal_basis(normal);
    SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 2.0),
        normal,
        u_axis,
        v_axis,
    }
}

#[test]
fn cap_circle_pcurve_follows_the_rational_circle_parameter() {
    // The control-net pcurve u = 2πt sat 0.012 rad off the rational
    // circle mid-quarter (6e-3 in model space at r = 0.5).
    let gap = worst_coparameter_gap(&cylinder(0.5), &plane(Vec3::new(0.0, 0.0, 1.0)));
    assert!(gap < 1.0e-8, "cap circle co-parameter gap {gap:e}");
}

#[test]
fn oblique_ellipse_pcurve_follows_the_rational_ellipse_parameter() {
    let gap = worst_coparameter_gap(&cylinder(3.0), &plane(Vec3::new(0.3, -0.2, 1.0)));
    assert!(gap < 1.0e-8, "oblique ellipse co-parameter gap {gap:e}");
}

#[test]
fn tiny_cylinder_pcurve_meets_the_gate_at_millimetre_scale() {
    // SI-stored millimetre parts (PK-462): a 0.5 mm radius in metres.
    let gap = worst_coparameter_gap(&cylinder(5.0e-4), &plane(Vec3::new(0.1, 0.0, 1.0)));
    assert!(gap < 1.0e-8, "millimetre-scale co-parameter gap {gap:e}");
}
