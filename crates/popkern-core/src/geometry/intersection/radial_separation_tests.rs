use super::wall_clear_of_cylinder;
use crate::geometry::nurbs_surface::NurbsSurface;
use crate::geometry::{Point3, Vec3};

/// A rational degree-2 x degree-1 tube about the z axis: radius
/// `r_bottom` at z = 0 and `r_top` at z = 1 (a cone frustum when they
/// differ), four quarter-circle knot spans around.
fn tube(r_bottom: f64, r_top: f64) -> NurbsSurface {
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let ring = [
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
    // u-major storage: flat index = i_u * num_v + i_v.
    let mut control_points = Vec::new();
    let mut weights = Vec::new();
    for (k, (x, y)) in ring.iter().enumerate() {
        for (z, r) in [(0.0, r_bottom), (1.0, r_top)] {
            control_points.push(Point3::new(x * r, y * r, z));
            weights.push(if k % 2 == 1 { w } else { 1.0 });
        }
    }
    NurbsSurface {
        degree_u: 2,
        degree_v: 1,
        knots_u: vec![
            0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
        ],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points,
        num_u: 9,
        num_v: 2,
        weights,
    }
}

fn clear(wall: &NurbsSurface, radius: f64) -> bool {
    wall_clear_of_cylinder(
        wall,
        Point3::new(0.0, 0.0, -5.0),
        Vec3::new(0.0, 0.0, 1.0),
        radius,
        1e-5,
    )
}

#[test]
fn a_wall_inside_the_radius_is_certified_clear() {
    assert!(clear(&tube(1.0, 1.0), 1.5));
}

#[test]
fn a_wall_outside_the_radius_is_certified_clear_span_by_span() {
    // The global hull of this net contains the axis — only the
    // per-span separating planes can certify it.
    assert!(clear(&tube(2.0, 2.0), 1.0));
}

#[test]
fn a_wall_on_the_radius_is_not_certified() {
    assert!(!clear(&tube(1.0, 1.0), 1.0));
}

#[test]
fn a_wall_straddling_the_radius_is_not_certified() {
    assert!(!clear(&tube(0.8, 1.2), 1.0));
}

#[test]
fn the_gap_keeps_grazing_walls_uncertified() {
    // Inside by 1e-6, below the 1e-5 gap: the solvers decide.
    assert!(!clear(&tube(1.0 - 1e-6, 1.0 - 1e-6), 1.0));
}

#[test]
fn a_bulging_net_is_certified_after_refinement() {
    // The rational circle's corner control points sit at r·√2: the
    // raw net of a radius-1 wall reaches 1.414, past a 1.2 cylinder,
    // though the wall itself stays at 1. Span refinement pulls the net
    // onto the wall and the inside certificate holds.
    assert!(clear(&tube(1.0, 1.0), 1.2));
}
