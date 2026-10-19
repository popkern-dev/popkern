// This file is only compiled under #[cfg(test)] (see periodic_form.rs).
use super::*;
use crate::geometry::{Point3, SurfaceEval};

/// A uniform quadratic ring in periodic form around a hexagon of
/// control points (radius 2), extruded linearly in `v`: six distinct
/// rows, the first two repeated at the end, knots uniform past both
/// ends of the domain [0, 1].
fn periodic_ring() -> NurbsSurface {
    let distinct = 6;
    let mut control_points = Vec::new();
    for k in 0..distinct + 2 {
        let a = std::f64::consts::TAU * (k % distinct) as f64 / distinct as f64;
        for z in [0.0, 3.0] {
            control_points.push(Point3::new(2.0 * a.cos(), 2.0 * a.sin(), z));
        }
    }
    let count = distinct + 2;
    NurbsSurface {
        degree_u: 2,
        degree_v: 1,
        knots_u: (0..count + 3)
            .map(|j| (j as f64 - 2.0) / distinct as f64)
            .collect(),
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        num_u: count,
        num_v: 2,
        weights: vec![1.0; count * 2],
        control_points,
    }
}

#[test]
fn a_periodic_form_ring_reports_its_domain_as_period() {
    let ring = periodic_ring();
    assert_eq!(ring.parameter_range().0, (0.0, 1.0));
    assert!((ring.u_period().unwrap_or(0.0) - 1.0).abs() < 1e-15);
    assert_eq!(ring.v_period(), None);
    assert!(ring.transposed().v_period().is_some());
}

#[test]
fn evaluation_wraps_and_is_smooth_across_the_seam() {
    let ring = periodic_ring();
    for &u in &[0.1, 0.37, 0.5, 0.93] {
        for &shift in &[1.0, -1.0, 3.0] {
            let a = ring.point_at(u, 0.4);
            let b = ring.point_at(u + shift, 0.4);
            assert!(a.distance_to(&b) < 1e-12, "S({u}) vs S({u}{shift:+})");
            let (_, du_a, dv_a) = ring.point_and_first_partials(u, 0.4);
            let (_, du_b, dv_b) = ring.point_and_first_partials(u + shift, 0.4);
            assert!((du_a - du_b).length() < 1e-9 && (dv_a - dv_b).length() < 1e-9);
            let a2 = ring.point_and_two_partials(u, 0.4);
            let b2 = ring.point_and_two_partials(u + shift, 0.4);
            assert!((a2.3 - b2.3).length() < 1e-8, "S_uu wraps");
        }
    }
    // The seam is a knot of a uniform quadratic: C1 across it.
    let (p0, d0, _) = ring.point_and_first_partials(0.0, 0.0);
    let (p1, d1, _) = ring.point_and_first_partials(1.0, 0.0);
    assert!(p0.distance_to(&p1) < 1e-12);
    assert!((d0 - d1).length() < 1e-9 * d0.length());
}

#[test]
fn an_open_direction_does_not_wrap() {
    let ring = periodic_ring();
    // Open in v: past the domain the end span extrapolates (z = 3v),
    // as before; it does not wrap back to v = 0.5.
    assert!((ring.point_at(0.3, 1.5).z - 4.5).abs() < 1e-12);
}

/// Negative controls: a clamped ring whose end rows coincide is closed
/// but not periodic (the exact rational circle wall, like a revolved
/// NURBS or a coil-spring tube), and neither is a periodic-form net
/// whose wrapped row or knot spacing is off.
#[test]
fn closed_clamped_or_broken_wraps_are_not_periodic() {
    let s = std::f64::consts::FRAC_1_SQRT_2;
    let circle = [
        (2.0, 0.0, 1.0),
        (2.0, 2.0, s),
        (0.0, 2.0, 1.0),
        (-2.0, 2.0, s),
        (-2.0, 0.0, 1.0),
        (-2.0, -2.0, s),
        (0.0, -2.0, 1.0),
        (2.0, -2.0, s),
        (2.0, 0.0, 1.0),
    ];
    let mut control_points = Vec::new();
    let mut weights = Vec::new();
    for &(x, y, w) in &circle {
        for z in [0.0, 3.0] {
            control_points.push(Point3::new(x, y, z));
            weights.push(w);
        }
    }
    let clamped = NurbsSurface {
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
    };
    assert_eq!(clamped.u_period(), None);
    // Past its domain a clamped ring extrapolates, as before.
    let past = clamped.point_at(1.1, 0.0);
    assert!(past.distance_to(&clamped.point_at(0.1, 0.0)) > 1e-3);

    let mut moved_row = periodic_ring();
    let last = moved_row.control_points.len() - 1;
    moved_row.control_points[last].x += 1e-9;
    assert_eq!(moved_row.u_period(), None);

    let mut reweighted = periodic_ring();
    let last = reweighted.weights.len() - 1;
    reweighted.weights[last] = 1.0 + 1e-9;
    assert_eq!(reweighted.u_period(), None);

    let mut uneven = periodic_ring();
    let last = uneven.knots_u.len() - 1;
    uneven.knots_u[last] += 1e-6;
    assert_eq!(uneven.u_period(), None);
}

/// The hint-seeded projection runs across a periodic seam instead of
/// stalling on the chart edge (PK-891): seeded just past u = 0, a
/// point a little before the seam projects onto itself. (The clamped
/// projector stalled on the chart edge 0.21 off a closed rail; the
/// angled-hole rim sweep is the end-to-end control.)
#[test]
fn hinted_projection_crosses_a_periodic_seam() {
    use crate::geometry::surface::SurfaceGeometry;
    use crate::geometry::uv_projection::project_point_to_uv_near;
    let ring = periodic_ring();
    let target = ring.point_at(-0.01, 0.5);
    let surface = SurfaceGeometry::Nurbs(ring.clone());
    let (u, v) = project_point_to_uv_near(&surface, &target, 0.02, 0.5).expect("projects");
    assert!(
        ring.point_at(u, v).distance_to(&target) < 1e-9,
        "u {u} v {v}"
    );
    assert!((u + 0.01).abs() < 1e-6, "stays on the hint's branch: u {u}");
}
