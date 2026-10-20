//! A plane section of a coiled NURBS wall written in periodic form
//! (PK-892): a 6.5-turn wire coil's tube wall cut by the plane z = 0.85,
//! one closed loop round the wire. The section's pcurve on the wall must
//! stay on one branch of the u period: read as an open chart (its
//! control net wrapped, not duplicated), the wall took the mixed-pair
//! wrapped-record path and the loop came back as one pcurve jumping a
//! whole period across the seam (a grind's cut edge 0.049 off its face).

use popkern_core::geometry::intersection::{SurfaceSurfaceIntersection, intersect_surfaces};
use popkern_core::geometry::surface::SurfaceGeometry;
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

fn periodic_coil_wall(tol: Tolerance) -> SurfaceGeometry {
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
    SurfaceGeometry::Nurbs(
        clamped
            .periodic_u_from_closed_clamped(tol.linear())
            .expect("the wire closes in u"),
    )
}

#[test]
fn a_plane_section_round_a_periodic_coil_wire_stays_on_one_u_branch() {
    let tol = Tolerance::default();
    let wall = periodic_coil_wall(tol);
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.85),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let result = intersect_surfaces(
        &wall,
        (0.0, 1.0),
        (0.0, 1.0),
        &plane,
        (-1.5, 1.5),
        (-1.5, 1.5),
        &tol,
    )
    .expect("section");
    let curves: Vec<&NurbsCurve> = result
        .intersections
        .iter()
        .filter_map(|i| match i {
            SurfaceSurfaceIntersection::Curve(c) => Some(&c.curve_uv1),
            _ => None,
        })
        .collect();
    assert!(!curves.is_empty(), "the plane cuts the wire");
    // The section winds once round the wire: on one branch it either
    // splits at the seam or runs a full period, never jumping one.
    let mut u_extent = 0.0;
    for uv in &curves {
        let jump = uv
            .control_points
            .windows(2)
            .map(|w| (w[1].x - w[0].x).abs())
            .fold(0.0_f64, f64::max);
        assert!(jump < 0.5, "a pcurve jumps {jump:.3} of the u period");
        let (lo, hi) = uv
            .control_points
            .iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), p| {
                (lo.min(p.x), hi.max(p.x))
            });
        u_extent += hi - lo;
    }
    // Together the pieces go once round the wire.
    assert!((u_extent - 1.0).abs() < 0.05, "u extent {u_extent:.3}");
}
