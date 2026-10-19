//! Projection onto a coiled NURBS wall written in periodic form
//! (PK-892): a 6.5-turn wire coil's tube wall, its closed profile
//! direction converted to periodic form. The 5-start grid projection
//! settled at stationary points on other turns and other sides of the
//! wire (0.18 and 0.23 off for points ON the wall); the clamped chart's
//! u border had caught such starts. Every on-surface point must project
//! back onto itself.

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_sweep::{HelixPitch, HelixSpec, SweepOptions, SweepProfile, helical_sweep};

/// The exact rational circle of radius `r` in the path frame's plane.
fn wire(r: f64) -> NurbsCurve {
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let corners = [
        (1.0, 0.0),
        (1.0, 1.0),
        (0.0, 1.0),
        (-1.0, 1.0),
        (-1.0, 0.0),
        (-1.0, -1.0),
        (0.0, -1.0),
        (1.0, -1.0),
        (1.0, 0.0),
    ];
    NurbsCurve {
        degree: 2,
        knots: vec![
            0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
        ],
        control_points: corners
            .iter()
            .map(|&(x, y)| Point3::new(r * x, r * y, 0.0))
            .collect(),
        weights: (0..9).map(|k| if k % 2 == 1 { w } else { 1.0 }).collect(),
    }
}

#[test]
fn points_on_a_periodic_coil_wall_project_onto_themselves() {
    let spec = HelixSpec {
        axis_origin: Point3::new(0.0, 0.0, 0.0),
        axis_dir: Vec3::new(0.0, 0.0, 1.0),
        start_radius: 0.3,
        end_radius: 0.3,
        pitch: HelixPitch::Constant(0.15),
        turns: 6.5,
        start_angle: 0.0,
        left_handed: false,
    };
    let tol = Tolerance::default();
    let coil = helical_sweep(
        &SweepProfile::Wire(wire(0.04)),
        &spec,
        &SweepOptions::default(),
        tol,
        None,
    )
    .expect("coil");
    let clamped = coil
        .brep
        .faces
        .values()
        .find_map(|f| match f.surface.as_ref().map(|s| &s.geometry) {
            Some(SurfaceGeometry::Nurbs(n)) => Some(n.clone()),
            _ => None,
        })
        .expect("the tube wall");
    let ring = clamped
        .periodic_u_from_closed_clamped(tol.linear())
        .expect("the wire closes in u");
    assert!(ring.u_period().is_some());
    let surface = SurfaceGeometry::Nurbs(ring);
    let mut worst: f64 = 0.0;
    for i in 0..20 {
        for j in 0..=40 {
            let (u, v) = (i as f64 / 20.0, j as f64 / 40.0);
            let p = surface.point_at(u, v);
            let (pu, pv) = project_point_to_uv(&surface, &p).expect("projects");
            worst = worst.max(surface.point_at(pu, pv).distance_to(&p));
        }
    }
    assert!(worst < 1e-9, "an on-wall point projected {worst} off");
}
