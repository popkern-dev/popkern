use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI, TAU};

use super::super::self_intersection_uv_containment::{
    RawLoopSample, TrimContainment, UvPoint, point_in_uv_loop,
};
use super::ruled_pole_polygon;

fn raw(points: &[(f64, f64)]) -> Vec<RawLoopSample> {
    points
        .iter()
        .map(|&(u, v)| RawLoopSample {
            uv: UvPoint { u, v },
            pole_jump_before: false,
        })
        .collect()
}

/// A sphere chart's poles: v = 0 and v = π.
fn on_row(v: f64) -> bool {
    v.abs() < 1e-9 || (v - PI).abs() < 1e-9
}

fn inside(polygon: &[UvPoint], u: f64, v: f64) -> bool {
    matches!(
        point_in_uv_loop(UvPoint { u, v }, polygon, 1e-9),
        TrimContainment::Interior
    )
}

/// A whole sphere's loop: down one seam side and up the other, both
/// stored at one u. The runs along the poles are a whole period, and
/// the polygon is the whole chart.
#[test]
fn a_seam_loop_on_one_branch_is_the_whole_chart() {
    let samples = raw(&[
        (0.0, PI),
        (0.0, FRAC_PI_2),
        (0.0, 0.0),
        (0.0, 0.0),
        (0.0, FRAC_PI_2),
        (0.0, PI),
    ]);
    let polygon = ruled_pole_polygon(&samples, &on_row, false).expect("closes");
    let width = polygon.iter().map(|p| p.u).fold(f64::MIN, f64::max)
        - polygon.iter().map(|p| p.u).fold(f64::MAX, f64::min);
    assert!((width - TAU).abs() < 1e-9, "{polygon:?}");
}

/// Viviani's cut with the ball's pole (v = π) at the double point: the
/// loop runs down a seam side, along the far pole, up the other side,
/// and through the pole three times round two lobe tongues hanging
/// from it. The tongues are outside the face, the band below them and
/// the pole's sectors between them inside.
#[test]
fn a_loop_visiting_a_pole_three_times_bounds_the_face_once() {
    let s = 3.0 * FRAC_PI_2;
    let samples = raw(&[
        (s, PI),
        (s, FRAC_PI_2),
        (s, 0.0),
        (s + TAU, 0.0),
        (s + TAU, FRAC_PI_2),
        (s + TAU, PI),
        (s + FRAC_PI_4, PI),
        (s + FRAC_PI_2, FRAC_PI_2),
        (s + 3.0 * FRAC_PI_4, PI),
        (s + 5.0 * FRAC_PI_4, PI),
        (s + 3.0 * FRAC_PI_2, FRAC_PI_2),
        (s + 7.0 * FRAC_PI_4, PI),
    ]);
    let polygon = ruled_pole_polygon(&samples, &on_row, false).expect("closes");
    let at = |du: f64, v: f64| inside(&polygon, s + du, v) || inside(&polygon, s + du - TAU, v);
    assert!(!at(FRAC_PI_2, 2.9), "inside the first tongue");
    assert!(!at(3.0 * FRAC_PI_2, 2.9), "inside the second tongue");
    assert!(at(FRAC_PI_2, 1.0), "below the tongues");
    assert!(at(PI, 3.0), "between the tongues at the pole");
    assert!(at(0.1, 3.1), "beside the seam at the pole");
}

/// A loop that never runs along a pole row is left to the folded
/// assembly.
#[test]
fn a_loop_off_the_poles_is_not_ruled() {
    let samples = raw(&[(0.0, 1.0), (1.0, 1.0), (1.0, 2.0), (0.0, 2.0)]);
    assert!(ruled_pole_polygon(&samples, &on_row, true).is_none());
}
