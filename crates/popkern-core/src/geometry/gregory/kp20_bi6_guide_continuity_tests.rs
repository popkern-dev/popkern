use super::*;

const EPS: f64 = 1e-12;

#[test]
fn guide_join_weights_match_kp18_eq1() {
    let (w0, w1, w2) = kp18_guide_join_weights(4).unwrap();
    assert!((w0 + 1.0).abs() < EPS);
    assert!(w1.abs() < EPS);
    assert!((w2 - 2.0).abs() < EPS);
    assert!(kp18_guide_join_weights(2).is_err());
}
