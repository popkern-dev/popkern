//! PK-646: rings on a NURBS surface closed in u that arrive a whole
//! period off the knot domain move back into it (a NURBS clamps when
//! evaluated outside its domain: jack body 56's band meshed 84.2 of
//! 127.5 with every interior point on the domain-end iso-line).

use super::*;
use popkern_core::geometry::Point3;
use popkern_core::geometry::nurbs_surface::NurbsSurface;

/// A radius-1 cylinder of height 1 as a rational NURBS: u the angle
/// over [0, 2π] (closed when `closed`), v the height.
fn tube(closed: bool) -> SurfaceGeometry {
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let square = [
        (1.0, 0.0),
        (1.0, 1.0),
        (0.0, 1.0),
        (-1.0, 1.0),
        (-1.0, 0.0),
        (-1.0, -1.0),
        (0.0, -1.0),
        (1.0, -1.0),
        (1.0, if closed { 0.0 } else { -0.5 }),
    ];
    let mut control_points = Vec::new();
    let mut weights = Vec::new();
    for (k, (x, y)) in square.iter().enumerate() {
        for z in [0.0, 1.0] {
            control_points.push(Point3::new(*x, *y, z));
            weights.push(if k % 2 == 0 { 1.0 } else { w });
        }
    }
    use std::f64::consts::{FRAC_PI_2, PI, TAU};
    SurfaceGeometry::Nurbs(NurbsSurface {
        degree_u: 2,
        degree_v: 1,
        knots_u: vec![
            0.0,
            0.0,
            0.0,
            FRAC_PI_2,
            FRAC_PI_2,
            PI,
            PI,
            PI + FRAC_PI_2,
            PI + FRAC_PI_2,
            TAU,
            TAU,
            TAU,
        ],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points,
        num_u: 9,
        num_v: 2,
        weights,
    })
}

fn ring(u_offset: f64) -> BoundaryRing {
    let uv: Vec<(f64, f64)> = (0..8)
        .map(|k| (u_offset + std::f64::consts::TAU * k as f64 / 7.0, 0.25))
        .collect();
    BoundaryRing {
        pos: uv.iter().map(|_| Point3::new(0.0, 0.0, 0.0)).collect(),
        sources: uv.iter().map(|_| (Default::default(), 0.0, 0.0)).collect(),
        uv,
    }
}

#[test]
fn rings_a_period_off_a_closed_nurbs_domain_move_into_it() {
    let tau = std::f64::consts::TAU;
    let mut rings = [ring(tau), ring(tau + 0.5)];
    shift_into_closed_nurbs_domain(&mut rings, &tube(true));
    assert!(
        rings[0]
            .uv
            .iter()
            .all(|q| (-1e-9..=tau + 1e-9).contains(&q.0)),
        "outer ring in the domain: {:?}",
        rings[0].uv
    );
    assert!(
        (rings[1].uv[0].0 - 0.5).abs() < 1e-12,
        "holes move with the outer ring"
    );
}

#[test]
fn rings_already_in_the_domain_or_on_open_surfaces_stay() {
    let tau = std::f64::consts::TAU;
    let mut inside = [ring(0.0)];
    shift_into_closed_nurbs_domain(&mut inside, &tube(true));
    assert_eq!(inside[0].uv[0].0, 0.0);
    // Not closed: no period to shift by.
    let mut open = [ring(tau)];
    shift_into_closed_nurbs_domain(&mut open, &tube(false));
    assert_eq!(open[0].uv[0].0, tau);
}

/// PK-1188: a ring passing beside a pole (a dip in its half-period
/// shift off the pole) widens the gate to the chord height; an edge
/// running into the pole (its samples halved down to the pole sample,
/// a die's pip-cap seam) does not, nor does a pass past a tenth of the
/// chord height.
#[test]
fn only_a_pass_beside_the_pole_widens_its_gate() {
    let chord = 5.0e-3;
    let gate = 2.7e-5;
    let pass = [1.6e-1, 6.0e-5, 0.0, 6.0e-5, 2.2e-5, 4.3e-5, 1.4e-4, 3.2e-4];
    assert_eq!(super::pole_pass_gate(&pass, gate, chord), chord);
    let into_pole = [
        1.3e-5, 6.4e-6, 3.2e-6, 4.3e-19, 3.2e-6, 6.4e-6, 1.3e-5, 2.5e-5,
    ];
    assert_eq!(super::pole_pass_gate(&into_pole, 2.9e-8, chord), 0.0);
    let far = [1.6e-1, 2.0e-3, 1.0e-3, 2.0e-3, 1.6e-1];
    assert_eq!(super::pole_pass_gate(&far, gate, chord), 0.0);
}
