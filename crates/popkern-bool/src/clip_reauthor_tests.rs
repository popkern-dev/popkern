use super::*;
use popkern_math::PortableMath;
use std::f64::consts::TAU;

/// A wave of four periods along x, amplitude 10: 16 cubic spans miss
/// it by far more than a tolerance between their samples.
fn wave(t: f64) -> Point3 {
    Point3::new(40.0 * t, 10.0 * (4.0 * TAU * t).psin(), 0.0)
}

#[test]
fn a_dense_piece_is_reauthored_at_the_density_it_needs() {
    let budget = 1e-5;
    let moved_end = Point3::new(40.0 + 1e-3, 0.0, 0.0);
    let fitted = reauthor_pinned_run(
        (0.0, 1.0),
        &|ts: &[f64]| ts.iter().map(|&t| wave(t)).collect(),
        &|pts: &mut [Point3]| {
            if let Some(last) = pts.last_mut() {
                *last = moved_end;
            }
        },
        &|p: &Point3| *p,
        budget,
    )
    .expect("fits");
    assert!(
        fitted.point_at(1.0).distance_to(&moved_end) < 1e-12,
        "the pin holds"
    );
    assert!(fitted.point_at(0.0).distance_to(&wave(0.0)) < 1e-12);
    // Away from the pinned end segment the refit is the wave.
    let worst = (0..=900)
        .map(|k| {
            let t = 0.9 * k as f64 / 900.0;
            fitted.point_at(t).distance_to(&wave(t))
        })
        .fold(0.0_f64, f64::max);
    assert!(worst <= 2.0 * budget, "worst {worst:.2e}");
    assert!(fitted.control_points.len() > 17, "densified");
}

#[test]
fn a_piece_seventeen_samples_reproduce_keeps_seventeen() {
    let line = |t: f64| Point3::new(t, 2.0 * t, 0.0);
    let fitted = reauthor_pinned_run(
        (0.0, 1.0),
        &|ts: &[f64]| ts.iter().map(|&t| line(t)).collect(),
        &|_: &mut [Point3]| {},
        &|p: &Point3| *p,
        1e-7,
    )
    .expect("fits");
    assert_eq!(fitted.control_points.len(), 17);
}
