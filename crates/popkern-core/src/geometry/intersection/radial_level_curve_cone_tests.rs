//! Cone side of the radial level-curve arm (PK-799): a thread flank
//! against a chamfer cone is the same radial level-set model with the
//! level radius growing along the axis.
//!
//! Fixture: the exact helical band of `radial_level_curve_tests`
//! (`ρ = 2 + 2u`, z ∈ [0, 3] over three turns) against a coaxial cone
//! whose radius `(z + 7.5)/3` runs 2.5 → 3.5 across the band, so every
//! station crosses once and the true contact is the curve where the
//! band's radius meets the cone's.

use std::f64::consts::TAU;

use super::radial_level_curve_tests::{curves, env_lock, helical_band};
use crate::geometry::curve::CurveEval;
use crate::geometry::intersection::radial_level_curve::try_radial_level_curve;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

/// Apex 7.5 below the band, slope 1/3, unit axis.
fn chamfer_cone() -> SurfaceGeometry {
    SurfaceGeometry::Cone {
        apex: Point3::new(0.0, 0.0, -7.5),
        axis: Vec3::new(0.0, 0.0, 1.0),
        half_angle: (1.0_f64 / 3.0).atan(),
    }
}

#[test]
fn chamfer_cone_crossing_serves_on_the_cone() {
    let _guard = env_lock();
    let wall = helical_band(2.0, 4.0, 3.0, 1.0);
    let cone = chamfer_cone();
    let tol = Tolerance::default();
    let v_range = (7.0, 11.0);
    let records = curves(
        try_radial_level_curve(
            &wall,
            (0.0, 1.0),
            (0.0, 1.0),
            &cone,
            (0.0, TAU),
            v_range,
            &tol,
        )
        .expect("arm claims the cone x band pair"),
    );
    assert!(
        records.len() >= 3,
        "three wraps, got {} records",
        records.len()
    );
    let band = 10.0 * tol.linear();
    let mut z_lo = f64::INFINITY;
    let mut z_hi = f64::NEG_INFINITY;
    for record in &records {
        let (t0, t1) = record.curve_3d.parameter_range();
        for i in 0..=64 {
            let t = t0 + (t1 - t0) * i as f64 / 64.0;
            let p = record.curve_3d.point_at(t);
            let rho = (p.x * p.x + p.y * p.y).sqrt();
            let cone_r = (p.z + 7.5) / 3.0;
            assert!(
                (rho - cone_r).abs() <= band,
                "record point off the cone by {:.3e}",
                (rho - cone_r).abs()
            );
            // The cone pcurve lands on the same 3D point.
            let cuv = record.curve_uv2.point_at(t);
            let q = cone.point_at(cuv.x, cuv.y);
            assert!(
                q.distance_to(&p) <= band,
                "cone pcurve misses the record point by {:.3e}",
                q.distance_to(&p)
            );
            z_lo = z_lo.min(p.z);
            z_hi = z_hi.max(p.z);
        }
    }
    assert!(
        z_lo < 0.05 && z_hi > 2.95,
        "records cover the whole band, z in [{z_lo}, {z_hi}]"
    );
}

#[test]
fn non_unit_cone_axis_declines() {
    let _guard = env_lock();
    // A Cone's chart scales its second radial direction by |axis|:
    // with axis length 2 the section is an ellipse, not a circle, so
    // the axisymmetric level field must not claim the pair.
    let wall = helical_band(2.0, 4.0, 3.0, 1.0);
    let cone = SurfaceGeometry::Cone {
        apex: Point3::new(0.0, 0.0, -7.5),
        axis: Vec3::new(0.0, 0.0, 2.0),
        half_angle: (2.0_f64 / 3.0).atan(),
    };
    let tol = Tolerance::default();
    let result = try_radial_level_curve(
        &wall,
        (0.0, 1.0),
        (0.0, 1.0),
        &cone,
        (0.0, TAU),
        (3.5, 5.5),
        &tol,
    );
    assert!(result.is_none(), "an elliptic cone declines to marching");
}

#[test]
fn cone_face_reaching_the_apex_declines() {
    let _guard = env_lock();
    // The level field's zero set is the v > 0 nappe only: a face
    // through the apex must fall back to marching, never serve a
    // one-nappe answer.
    let wall = helical_band(2.0, 4.0, 3.0, 1.0);
    let tol = Tolerance::default();
    let result = try_radial_level_curve(
        &wall,
        (0.0, 1.0),
        (0.0, 1.0),
        &chamfer_cone(),
        (0.0, TAU),
        (-1.0, 11.0),
        &tol,
    );
    assert!(result.is_none(), "a face through the apex declines");
}
