//! The spoke normal sine is the true normal-angle sine, invariant to
//! chart skew, and its analytic gradient matches finite differences.

use super::*;

fn normal_sine_direct(u: Vec3, v: Vec3, w: Vec3) -> f64 {
    let (na, nb) = (w.cross(&u), w.cross(&v));
    na.cross(&nb).length() / (na.length() * nb.length())
}

#[test]
fn it_is_the_normal_angle_sine_and_blind_spots_of_coplanarity_are_gone() {
    // A 4° crease between two tangent planes sharing the spoke z.
    let w = Vec3::new(0.0, 0.0, 1.0);
    let theta = 4.0_f64.to_radians();
    let u = Vec3::new(1.0, 0.0, 0.0);
    let v = Vec3::new(theta.cos(), theta.sin(), 0.0);
    let straight = spoke_normal_sine(u, v, w);
    assert!((straight - theta.sin()).abs() < 1e-12);
    // Skew both cross tangents toward the spoke (1.7° away): the
    // normalized triple product collapses, the sine does not move.
    let lean = 1.7_f64.to_radians();
    let (us, vs) = (
        u * lean.sin() + w * lean.cos(),
        v * lean.sin() + w * lean.cos(),
    );
    let triple = us.cross(&vs).dot(&w).abs() / (us.length() * vs.length() * w.length());
    assert!(triple < 1e-4, "coplanarity reads {triple:e}");
    let skewed = spoke_normal_sine(us, vs, w);
    assert!((skewed - normal_sine_direct(us, vs, w)).abs() < 1e-9);
    assert!(
        (skewed - theta.sin()).abs() < 1e-9,
        "skew moved the sine: {skewed}"
    );
}

#[test]
fn a_fold_across_the_spoke_is_not_g1_but_a_mirrored_chart_is() {
    let w = Vec3::new(0.0, 0.0, 1.0);
    let u = Vec3::new(1.0, 0.0, 0.0);
    let aligned = Vec3::new(1.0, 1e-3, 0.0);
    let folded = Vec3::new(-1.0, 1e-3, 0.0);
    // One folded sample among aligned ones reads 1.
    let mut samples = vec![(u, aligned, w); 9];
    samples.push((u, folded, w));
    assert_eq!(worst_spoke_measure(&samples), 1.0);
    // A spoke mirrored throughout is a consistent chart, not a fold.
    let mirrored = vec![(u, folded, w); 10];
    assert!(worst_spoke_measure(&mirrored) < 1e-2);
}

#[test]
fn the_gradient_matches_central_differences() {
    let u = Vec3::new(0.9, 0.3, 0.2);
    let v = Vec3::new(0.7, -0.4, 0.35);
    let w = Vec3::new(0.1, 0.2, 1.1);
    let Some((du, dv, dw)) = spoke_normal_sine_gradient(u, v, w) else {
        panic!("gradient defined");
    };
    let h = 1e-6;
    let axes = [
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
    ];
    for (k, e) in axes.iter().enumerate() {
        let fd = |g: &dyn Fn(Vec3) -> f64| (g(*e * h) - g(*e * -h)) / (2.0 * h);
        let pick = |d: Vec3| [d.x, d.y, d.z][k];
        let gu = fd(&|d| spoke_normal_sine(u + d, v, w));
        let gv = fd(&|d| spoke_normal_sine(u, v + d, w));
        let gw = fd(&|d| spoke_normal_sine(u, v, w + d));
        assert!((gu - pick(du)).abs() < 1e-7, "du[{k}] {gu} vs {}", pick(du));
        assert!((gv - pick(dv)).abs() < 1e-7, "dv[{k}] {gv} vs {}", pick(dv));
        assert!((gw - pick(dw)).abs() < 1e-7, "dw[{k}] {gw} vs {}", pick(dw));
    }
}
