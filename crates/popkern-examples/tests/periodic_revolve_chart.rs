//! A full revolution whose profile is not a recognised conic span is a
//! NURBS wall, written in periodic form in its angle direction (PK-893):
//! a ring like a cylinder's lateral face, not a clamped chart holding
//! its seam twice. The body validates, and its volume is Pappus's.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_sweep::SweepProfile;
use popkern_sweep::revolve;

/// A closed profile in the XZ plane: a quadratic bulge from (1, 0) out
/// to x = 1.5 at z = 0.5 and back to (1, 1), closed by the line x = 1.
/// The bulge is a parabola in x(z), not a circle, so no span is a
/// recognised conic.
fn bulge_profile() -> SweepProfile {
    SweepProfile::Wire(NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.5),
            Point3::new(1.0, 0.0, 1.0),
            Point3::new(1.0, 0.0, 0.5),
            Point3::new(1.0, 0.0, 0.0),
        ],
        weights: vec![1.0; 5],
    })
}

#[test]
fn a_full_nurbs_revolution_is_a_periodic_ring() {
    let tol = Tolerance::default();
    let result = revolve(
        &bulge_profile(),
        &Point3::new(0.0, 0.0, 0.0),
        &Vec3::new(0.0, 0.0, 1.0),
        std::f64::consts::TAU,
        tol,
        None,
    )
    .expect("revolve");
    let issues = result.brep.validate_all();
    assert!(issues.is_empty(), "{issues:?}");
    let walls: Vec<_> = result
        .brep
        .faces
        .values()
        .filter_map(|f| match f.surface.as_ref().map(|s| &s.geometry) {
            Some(SurfaceGeometry::Nurbs(n)) => Some(n),
            _ => None,
        })
        .collect();
    assert!(!walls.is_empty(), "the bulge revolves to a NURBS wall");
    for wall in &walls {
        assert!(
            wall.u_period().is_some(),
            "angle direction in periodic form"
        );
    }
    // Pappus: V = 2π ∫ x dA over the profile region. With w(z) =
    // 2z(1 − z) the bulge's offset from x = 1, ∫ x dA = ∫ (w + w²/2) dz
    // = 1/3 + 1/15 = 0.4.
    let volume = popkern_measure::body_volume(&result.brep, &tol).expect("volume");
    let pappus = std::f64::consts::TAU * 0.4;
    assert!(
        (volume - pappus).abs() / pappus < 1e-6,
        "volume {volume} vs Pappus {pappus}"
    );
}
