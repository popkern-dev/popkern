// Structural lint marker: this file is test code, included via #[path] under #[cfg(test)].

//! Partial revolutions of profiles with an edge on the axis (sectors).

use super::*;
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use std::f64::consts::{FRAC_PI_2, PI};

fn polyline(points: &[[f64; 2]]) -> SweepProfile {
    let n = points.len();
    let mut knots = vec![0.0, 0.0];
    knots.extend((1..n - 1).map(|i| i as f64 / (n - 1) as f64));
    knots.extend([1.0, 1.0]);
    SweepProfile::Wire(NurbsCurve {
        degree: 1,
        knots,
        control_points: points
            .iter()
            .map(|p| Point3::new(p[0], 0.0, p[1]))
            .collect(),
        weights: vec![1.0; n],
    })
}

fn counts(brep: &BRep) -> (usize, usize, usize) {
    (brep.vertices.len(), brep.edges.len(), brep.faces.len())
}

fn turn(profile: &SweepProfile, angle: f64) -> Result<SweepResult, SweepError> {
    revolve(
        profile,
        &Point3::new(0.0, 0.0, 0.0),
        &Vec3::new(0.0, 0.0, 1.0),
        angle,
        Tolerance::default(),
        None,
    )
}

/// A quarter of a cylinder: three side faces (bottom fan, wall, top
/// fan) and two rectangular caps meeting on the axis edge.
#[test]
fn a_quarter_cylinder_shares_its_axis_edge_between_the_caps() {
    let rect = polyline(&[[0.0, 0.0], [2.0, 0.0], [2.0, 3.0], [0.0, 3.0], [0.0, 0.0]]);
    let brep = turn(&rect, FRAC_PI_2).expect("sector").brep;
    assert_eq!(counts(&brep), (6, 9, 5));
    brep.validate().expect("sector validates");
}

/// A cone's sector closes its wall on the apex: a triangle face.
#[test]
fn a_cone_sector_closes_on_its_apex() {
    let triangle = polyline(&[[0.0, 0.0], [2.0, 0.0], [0.0, 3.0], [0.0, 0.0]]);
    let brep = turn(&triangle, PI).expect("sector").brep;
    // Poles at the base centre and apex, the rim at both angles.
    assert_eq!(counts(&brep), (4, 6, 4));
    brep.validate().expect("cone sector validates");
}

/// A profile meeting the axis at a vertex between its poles would pinch
/// the sector there: it declines naming the point.
#[test]
fn a_profile_touching_the_axis_between_its_poles_declines() {
    let bowtie = polyline(&[
        [0.0, 0.0],
        [1.0, 0.0],
        [0.0, 1.0],
        [1.0, 2.0],
        [0.0, 2.0],
        [0.0, 0.0],
    ]);
    let err = turn(&bowtie, FRAC_PI_2).expect_err("declines");
    assert!(err.to_string().contains("between its poles"), "{err}");
}

/// Two separate runs on the axis decline.
#[test]
fn a_profile_meeting_the_axis_twice_declines() {
    let notched = polyline(&[
        [0.0, 0.0],
        [1.0, 0.0],
        [1.0, 1.0],
        [0.0, 1.0],
        [0.0, 2.0],
        [1.0, 2.0],
        [1.0, 3.0],
        [0.0, 3.0],
        [0.0, 0.0],
    ]);
    // The axis spans (0,1)-(0,2) and (0,3)-(0,0) are two runs.
    let err = turn(&notched, FRAC_PI_2).expect_err("declines");
    assert!(err.to_string().contains("separate runs"), "{err}");
}

/// A spline profile is no canonical surface: the sector's wall stays
/// NURBS, written angle-in-u like a full revolution's, and still closes
/// on both poles.
#[test]
fn a_spline_sector_closes_on_its_poles() {
    // A quadratic dome from the bottom pole to the top one, closed
    // down the axis by a collinear quadratic run.
    let closed = SweepProfile::Wire(NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(1.5, 0.0, 2.0),
            Point3::new(0.0, 0.0, 2.5),
            Point3::new(0.0, 0.0, 1.25),
            Point3::new(0.0, 0.0, 0.0),
        ],
        weights: vec![1.0; 6],
    });
    let brep = turn(&closed, FRAC_PI_2).expect("sector").brep;
    brep.validate().expect("spline sector validates");
    let nurbs_walls = brep
        .faces
        .values()
        .filter(|f| {
            matches!(
                f.surface.as_ref().map(|s| &s.geometry),
                Some(popkern_core::geometry::surface::SurfaceGeometry::Nurbs(_))
            )
        })
        .count();
    assert!(nurbs_walls >= 1, "{nurbs_walls}");
}
